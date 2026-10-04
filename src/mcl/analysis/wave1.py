from __future__ import annotations

from collections.abc import Iterable

import pandas as pd


VALID_COMPARISON_MODES = {"vs_gene_wildtype", "vs_other_gene_variants"}


def _clean(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _truthy(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    return _clean(series).str.lower().isin({"true", "1", "yes", "y", "t"})


def prepare_qualifying_mutations(mutations: pd.DataFrame) -> pd.DataFrame:
    """Normalize registry mutation evidence and retain qualifying alterations only."""
    required = {"model_id", "gene"}
    missing = required - set(mutations.columns)
    if missing:
        raise ValueError(f"mutations is missing columns: {sorted(missing)}")

    frame = mutations.copy()
    frame["model_id"] = _clean(frame["model_id"])
    frame["gene"] = _clean(frame["gene"]).str.upper()
    frame["protein_change"] = _clean(
        frame.get("protein_change", pd.Series("", index=frame.index))
    )

    for column in [
        "driver",
        "hotspot",
        "oncogene_high_impact",
        "tumor_suppressor_high_impact",
        "likely_lof",
    ]:
        if column not in frame.columns:
            frame[column] = False
        else:
            frame[column] = _truthy(frame[column])

    impact = _clean(frame.get("vep_impact", pd.Series("", index=frame.index))).str.upper()
    frame["high_impact"] = (
        impact.eq("HIGH")
        | frame["oncogene_high_impact"]
        | frame["tumor_suppressor_high_impact"]
    )
    qualifying = (
        frame["driver"]
        | frame["hotspot"]
        | frame["high_impact"]
        | frame["likely_lof"]
    )
    return frame[qualifying & frame["model_id"].ne("") & frame["gene"].ne("")].copy()


def materialize_wave1_cohorts(
    models: pd.DataFrame,
    mutations: pd.DataFrame,
    contexts: list[dict],
    *,
    crispr_model_ids: Iterable[str],
    sequenced_model_ids: Iterable[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Materialize explicit Wave 1 context/comparator cohorts.

    The grouping contract intentionally mirrors Cancer Context Registry discovery:
    exact qualifying protein-change models form the context; wildtype comparators
    have no qualifying alteration in the defining gene; allele-specific comparators
    carry another qualifying alteration in the same gene.
    """
    required_models = {"model_id", "oncotree_code"}
    missing = required_models - set(models.columns)
    if missing:
        raise ValueError(f"models is missing columns: {sorted(missing)}")

    crispr_ids = {str(x).strip() for x in crispr_model_ids if str(x).strip()}
    sequenced_ids = {str(x).strip() for x in sequenced_model_ids if str(x).strip()}

    model_frame = models.copy()
    model_frame["model_id"] = _clean(model_frame["model_id"])
    model_frame["oncotree_code"] = _clean(model_frame["oncotree_code"])
    model_frame = model_frame[model_frame["model_id"].isin(crispr_ids)].copy()
    model_frame = model_frame.drop_duplicates("model_id")

    qualifying = prepare_qualifying_mutations(mutations)
    qualifying = qualifying[qualifying["model_id"].isin(crispr_ids)].copy()

    cohort_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    seen_ids: set[str] = set()
    for cfg in contexts:
        wave1_id = str(cfg.get("wave1_id") or "").strip()
        if not wave1_id:
            raise ValueError("Every Wave 1 context requires wave1_id")
        if wave1_id in seen_ids:
            raise ValueError(f"Duplicate wave1_id: {wave1_id}")
        seen_ids.add(wave1_id)

        code = str(cfg.get("oncotree_code") or "").strip()
        gene = str(cfg.get("gene") or "").strip().upper()
        protein_change = str(cfg.get("protein_change") or "").strip()
        comparison_mode = str(cfg.get("comparison_mode") or "").strip()
        if comparison_mode not in VALID_COMPARISON_MODES:
            raise ValueError(
                f"{wave1_id}: unsupported comparison_mode {comparison_mode!r}; "
                f"choose one of {sorted(VALID_COMPARISON_MODES)}"
            )
        if not code or not gene or not protein_change:
            raise ValueError(f"{wave1_id}: oncotree_code, gene and protein_change are required")

        disease_models = model_frame[model_frame["oncotree_code"] == code].copy()
        disease_ids = set(disease_models["model_id"]) & sequenced_ids
        disease_mut = qualifying[
            qualifying["model_id"].isin(disease_ids) & qualifying["gene"].eq(gene)
        ].copy()
        altered_ids = set(disease_mut["model_id"])
        exact_ids = set(
            disease_mut.loc[disease_mut["protein_change"].eq(protein_change), "model_id"]
        )

        if comparison_mode == "vs_gene_wildtype":
            comparator_ids = disease_ids - altered_ids
            comparator_reason = f"no qualifying {gene} alteration detected under registry rules"
        else:
            comparator_ids = altered_ids - exact_ids
            comparator_reason = f"another qualifying {gene} alteration, excluding {protein_change}"

        excluded_ids = disease_ids - exact_ids - comparator_ids
        model_lookup = disease_models.set_index("model_id", drop=False)

        for model_id in sorted(disease_ids):
            if model_id in exact_ids:
                group = "context"
                reason = f"qualifying {gene} {protein_change}"
            elif model_id in comparator_ids:
                group = "comparator"
                reason = comparator_reason
            else:
                group = "excluded"
                reason = "same disease but does not satisfy context/comparator rule"

            meta = model_lookup.loc[model_id] if model_id in model_lookup.index else None
            row: dict[str, object] = {
                "wave1_id": wave1_id,
                "label": cfg.get("label") or wave1_id,
                "oncotree_code": code,
                "gene": gene,
                "protein_change": protein_change,
                "comparison_mode": comparison_mode,
                "model_id": model_id,
                "analysis_group": group,
                "assignment_reason": reason,
                "sequencing_available": True,
            }
            for column in [
                "cell_line_name",
                "oncotree_subtype",
                "oncotree_primary_disease",
                "oncotree_lineage",
                "depmap_model_type",
            ]:
                row[column] = meta.get(column) if meta is not None and column in meta.index else None
            cohort_rows.append(row)

        expected_context = int(cfg.get("expected_context_n") or 0)
        expected_comparator = int(cfg.get("expected_comparator_n") or 0)
        actual_context = len(exact_ids)
        actual_comparator = len(comparator_ids)
        summary_rows.append(
            {
                "wave1_id": wave1_id,
                "label": cfg.get("label") or wave1_id,
                "oncotree_code": code,
                "gene": gene,
                "protein_change": protein_change,
                "comparison_mode": comparison_mode,
                "comparison_role": cfg.get("comparison_role"),
                "context_models_n": actual_context,
                "comparator_models_n": actual_comparator,
                "excluded_models_n": len(excluded_ids),
                "expected_context_n": expected_context,
                "expected_comparator_n": expected_comparator,
                "context_count_matches": actual_context == expected_context,
                "comparator_count_matches": actual_comparator == expected_comparator,
                "counts_match_discovery": (
                    actual_context == expected_context
                    and actual_comparator == expected_comparator
                ),
                "comparator_is_operational_proxy": comparison_mode == "vs_gene_wildtype",
            }
        )

    cohort = pd.DataFrame(cohort_rows)
    summary = pd.DataFrame(summary_rows)
    return cohort, summary
