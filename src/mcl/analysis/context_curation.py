from __future__ import annotations

import re

import pandas as pd


_REVIEW_ORDER = {
    "DRIVER_HOTSPOT_FIRST_REVIEW": 0,
    "LOF_SECOND_REVIEW": 1,
    "OTHER_REVIEW": 2,
}


def _as_int(value: object) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y", "t"}


def _alteration_class(value: object) -> str:
    alteration = str(value or "").strip()
    if not alteration:
        return "unknown"
    if "fs" in alteration.lower():
        return "frameshift"
    if "Ter" in alteration or "*" in alteration:
        return "truncating"
    if re.fullmatch(r"p\.[A-Za-z][0-9]+[A-Za-z]", alteration):
        return "missense"
    return "other_exact"


def curate_registry_candidates(frame: pd.DataFrame) -> pd.DataFrame:
    """Create a transparent biological-review queue from READY discovery candidates.

    This function does not promote contexts into config/cancer_contexts.yaml and does
    not calculate an opaque global score. It only derives explicit review flags from
    candidate-level evidence already present in the discovery table.
    """
    required = {
        "candidate_id",
        "event_type",
        "alteration_gene",
        "alteration",
        "comparison_mode",
        "context_models_n",
        "comparator_models_n",
        "driver_context_models_n",
        "hotspot_context_models_n",
        "high_impact_context_models_n",
        "likely_lof_context_models_n",
        "readiness_status",
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"candidate table is missing columns: {sorted(missing)}")

    data = frame.copy()
    if "already_configured" not in data.columns:
        data["already_configured"] = False

    data = data[
        data["readiness_status"].astype(str).eq("READY")
        & data["event_type"].astype(str).eq("exact_protein_change")
        & ~data["already_configured"].map(_as_bool)
    ].copy()
    if data.empty:
        return data

    numeric_cols = [
        "context_models_n",
        "comparator_models_n",
        "driver_context_models_n",
        "hotspot_context_models_n",
        "high_impact_context_models_n",
        "likely_lof_context_models_n",
    ]
    for col in numeric_cols:
        data[col] = data[col].map(_as_int)

    denominator = data["context_models_n"].where(data["context_models_n"] > 0, 1)
    data["driver_fraction"] = data["driver_context_models_n"] / denominator
    data["hotspot_fraction"] = data["hotspot_context_models_n"] / denominator
    data["high_impact_fraction"] = data["high_impact_context_models_n"] / denominator
    data["likely_lof_fraction"] = data["likely_lof_context_models_n"] / denominator
    data["alteration_class"] = data["alteration"].map(_alteration_class)

    data["context_key"] = (
        data.get("oncotree_code", pd.Series("", index=data.index)).fillna("").astype(str)
        + "|"
        + data["alteration_gene"].fillna("").astype(str).str.upper()
        + "|"
        + data["alteration"].fillna("").astype(str)
    )

    comparison_roles: list[str] = []
    for row in data.itertuples(index=False):
        context_key = str(getattr(row, "context_key"))
        mode = str(getattr(row, "comparison_mode"))
        group = data[data["context_key"] == context_key]
        has_wildtype = bool((group["comparison_mode"].astype(str) == "vs_gene_wildtype").any())
        if mode == "vs_gene_wildtype":
            role = "primary_comparator"
        elif has_wildtype:
            role = "secondary_allele_specific_comparator"
        else:
            role = "primary_available_comparator"
        comparison_roles.append(role)
    data["comparison_role"] = comparison_roles

    review_classes: list[str] = []
    for row in data.itertuples(index=False):
        alteration_class = str(getattr(row, "alteration_class"))
        driver_fraction = float(getattr(row, "driver_fraction"))
        hotspot_fraction = float(getattr(row, "hotspot_fraction"))
        high_impact_fraction = float(getattr(row, "high_impact_fraction"))
        likely_lof_fraction = float(getattr(row, "likely_lof_fraction"))

        if alteration_class in {"missense", "other_exact"} and max(driver_fraction, hotspot_fraction) >= 0.8:
            review_class = "DRIVER_HOTSPOT_FIRST_REVIEW"
        elif alteration_class in {"frameshift", "truncating"} and max(
            high_impact_fraction,
            likely_lof_fraction,
            driver_fraction,
            hotspot_fraction,
        ) >= 0.8:
            review_class = "LOF_SECOND_REVIEW"
        else:
            review_class = "OTHER_REVIEW"
        review_classes.append(review_class)
    data["review_class"] = review_classes

    data["co_mutation_burden_review_required"] = data["alteration_class"].isin({"frameshift", "truncating"})
    data["automatic_promotion_allowed"] = False
    data["curation_note"] = (
        "Discovery evidence only. Requires manual review of disease definition, molecular mechanism, "
        "co-mutation structure and comparator before promotion."
    )

    data["_review_order"] = data["review_class"].map(_REVIEW_ORDER).fillna(99)
    data["_comparison_order"] = data["comparison_role"].map(
        {
            "primary_comparator": 0,
            "primary_available_comparator": 1,
            "secondary_allele_specific_comparator": 2,
        }
    ).fillna(9)
    data = data.sort_values(
        [
            "_review_order",
            "_comparison_order",
            "hotspot_fraction",
            "driver_fraction",
            "context_models_n",
            "comparator_models_n",
            "context_key",
        ],
        ascending=[True, True, False, False, False, False, True],
    ).drop(columns=["_review_order", "_comparison_order"])
    return data.reset_index(drop=True)
