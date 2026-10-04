from __future__ import annotations

import re
from collections.abc import Iterable

import pandas as pd


_READY_ORDER = {
    "READY": 0,
    "EXPLORATORY_LOW_N": 1,
    "EXPLORATORY_LOW_COMPARATOR_N": 2,
    "NO_COMPARATOR": 3,
}


def _truthy(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1", "yes", "y", "t"})
    )


def _clean_text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _slug(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9]+", "-", value.strip()).strip("-")
    return text.upper() or "UNKNOWN"


def _readiness(context_n: int, comparator_n: int, min_context_n: int, min_comparator_n: int) -> str:
    if comparator_n <= 0:
        return "NO_COMPARATOR"
    if context_n < min_context_n:
        return "EXPLORATORY_LOW_N"
    if comparator_n < min_comparator_n:
        return "EXPLORATORY_LOW_COMPARATOR_N"
    return "READY"


def _priority_tier(status: str, event_type: str) -> str:
    if status == "READY" and event_type == "exact_protein_change":
        return "A"
    if status == "READY":
        return "B"
    return "C"


def _disease_key(frame: pd.DataFrame) -> pd.Series:
    code = _clean_text(frame.get("oncotree_code", pd.Series("", index=frame.index)))
    subtype = _clean_text(frame.get("oncotree_subtype", pd.Series("", index=frame.index)))
    primary = _clean_text(frame.get("oncotree_primary_disease", pd.Series("", index=frame.index)))
    fallback = subtype.where(subtype.ne(""), primary)
    return code.where(code.ne(""), fallback)


def _representative(group: pd.DataFrame, column: str) -> str | None:
    if column not in group:
        return None
    values = _clean_text(group[column])
    values = values[values.ne("")]
    return values.iloc[0] if not values.empty else None


def _candidate_row(
    *,
    disease_group: pd.DataFrame,
    disease_key: str,
    gene: str,
    alteration: str,
    event_type: str,
    comparison_mode: str,
    context_ids: set[str],
    comparator_ids: set[str],
    sequenced_crispr_ids: set[str],
    evidence: pd.DataFrame,
    min_context_n: int,
    min_comparator_n: int,
) -> dict[str, object]:
    context_n = len(context_ids)
    comparator_n = len(comparator_ids)
    status = _readiness(context_n, comparator_n, min_context_n, min_comparator_n)

    evidence_context = evidence[evidence["model_id"].isin(context_ids)].copy()
    driver_n = int(evidence_context.loc[evidence_context["driver"], "model_id"].nunique())
    hotspot_n = int(evidence_context.loc[evidence_context["hotspot"], "model_id"].nunique())
    high_impact_n = int(evidence_context.loc[evidence_context["high_impact"], "model_id"].nunique())
    likely_lof_n = int(evidence_context.loc[evidence_context["likely_lof"], "model_id"].nunique())

    comparator_label = {
        "vs_gene_wildtype": f"same disease with mutation profiling and no qualifying {gene} alteration",
        "vs_other_gene_variants": f"same disease with another qualifying {gene} alteration but not {alteration}",
    }[comparison_mode]

    candidate_id = "CTX-" + "-".join(
        [
            _slug(disease_key),
            _slug(gene),
            _slug(alteration),
            _slug(comparison_mode),
        ]
    )

    return {
        "candidate_id": candidate_id,
        "disease_key": disease_key,
        "oncotree_code": _representative(disease_group, "oncotree_code"),
        "oncotree_subtype": _representative(disease_group, "oncotree_subtype"),
        "oncotree_primary_disease": _representative(disease_group, "oncotree_primary_disease"),
        "oncotree_lineage": _representative(disease_group, "oncotree_lineage"),
        "depmap_model_type": _representative(disease_group, "depmap_model_type"),
        "event_type": event_type,
        "alteration_gene": gene,
        "alteration": alteration,
        "comparison_mode": comparison_mode,
        "context_definition": f"{disease_key} models with {gene} {alteration}",
        "comparator_definition": comparator_label,
        "sequenced_crispr_models_n": len(sequenced_crispr_ids),
        "context_models_n": context_n,
        "comparator_models_n": comparator_n,
        "driver_context_models_n": driver_n,
        "hotspot_context_models_n": hotspot_n,
        "high_impact_context_models_n": high_impact_n,
        "likely_lof_context_models_n": likely_lof_n,
        "minimum_context_n": min_context_n,
        "minimum_comparator_n": min_comparator_n,
        "readiness_status": status,
        "priority_tier": _priority_tier(status, event_type),
    }


def discover_context_candidates(
    models: pd.DataFrame,
    mutations: pd.DataFrame,
    *,
    crispr_model_ids: Iterable[str],
    sequenced_model_ids: Iterable[str],
    min_discovery_n: int = 3,
    min_context_n: int = 5,
    min_comparator_n: int = 5,
) -> pd.DataFrame:
    """Discover recurrent molecular contexts suitable for DepMap comparison design.

    This is a discovery layer, not an automatic promotion step. Candidate contexts
    are restricted to models with CRISPR Gene Effect coverage and to mutation-defined
    events supported by driver/hotspot/high-impact/likely-LoF annotations.
    """
    required_models = {"model_id"}
    required_mutations = {"model_id", "gene"}
    missing_models = required_models - set(models.columns)
    missing_mutations = required_mutations - set(mutations.columns)
    if missing_models:
        raise ValueError(f"models is missing columns: {sorted(missing_models)}")
    if missing_mutations:
        raise ValueError(f"mutations is missing columns: {sorted(missing_mutations)}")

    crispr_ids = {str(x).strip() for x in crispr_model_ids if str(x).strip()}
    sequenced_ids = {str(x).strip() for x in sequenced_model_ids if str(x).strip()}

    model_frame = models.copy()
    model_frame["model_id"] = _clean_text(model_frame["model_id"])
    model_frame = model_frame[model_frame["model_id"].isin(crispr_ids)].copy()
    model_frame["disease_key"] = _disease_key(model_frame)
    model_frame = model_frame[model_frame["disease_key"].ne("")].drop_duplicates("model_id")

    mutation_frame = mutations.copy()
    mutation_frame["model_id"] = _clean_text(mutation_frame["model_id"])
    mutation_frame["gene"] = _clean_text(mutation_frame["gene"]).str.upper()
    mutation_frame["protein_change"] = _clean_text(
        mutation_frame.get("protein_change", pd.Series("", index=mutation_frame.index))
    )
    mutation_frame = mutation_frame[
        mutation_frame["model_id"].isin(crispr_ids) & mutation_frame["gene"].ne("")
    ].copy()

    for column in ["driver", "hotspot", "oncogene_high_impact", "tumor_suppressor_high_impact", "likely_lof"]:
        if column not in mutation_frame:
            mutation_frame[column] = False
        else:
            mutation_frame[column] = _truthy(mutation_frame[column])

    impact = _clean_text(mutation_frame.get("vep_impact", pd.Series("", index=mutation_frame.index))).str.upper()
    mutation_frame["high_impact"] = (
        impact.eq("HIGH")
        | mutation_frame["oncogene_high_impact"]
        | mutation_frame["tumor_suppressor_high_impact"]
    )
    qualifying = (
        mutation_frame["driver"]
        | mutation_frame["hotspot"]
        | mutation_frame["high_impact"]
        | mutation_frame["likely_lof"]
    )
    mutation_frame = mutation_frame[qualifying].copy()
    if mutation_frame.empty or model_frame.empty:
        return pd.DataFrame()

    mutation_frame = mutation_frame.merge(
        model_frame[["model_id", "disease_key"]],
        on="model_id",
        how="inner",
    )

    rows: list[dict[str, object]] = []
    for disease_key, disease_models in model_frame.groupby("disease_key", sort=True):
        disease_ids = set(disease_models["model_id"])
        sequenced_crispr_ids = disease_ids & sequenced_ids
        if len(sequenced_crispr_ids) < min_discovery_n:
            continue

        disease_mut = mutation_frame[mutation_frame["disease_key"] == disease_key].copy()
        if disease_mut.empty:
            continue

        # Gene-level altered-vs-wildtype candidates.
        for gene, gene_rows in disease_mut.groupby("gene", sort=True):
            altered_ids = set(gene_rows["model_id"]) & sequenced_crispr_ids
            if len(altered_ids) < min_discovery_n:
                continue
            wt_ids = sequenced_crispr_ids - altered_ids
            rows.append(
                _candidate_row(
                    disease_group=disease_models,
                    disease_key=str(disease_key),
                    gene=str(gene),
                    alteration="qualifying mutation",
                    event_type="gene_altered",
                    comparison_mode="vs_gene_wildtype",
                    context_ids=altered_ids,
                    comparator_ids=wt_ids,
                    sequenced_crispr_ids=sequenced_crispr_ids,
                    evidence=gene_rows,
                    min_context_n=min_context_n,
                    min_comparator_n=min_comparator_n,
                )
            )

            # Exact recurrent protein-change candidates.
            exact_rows = gene_rows[gene_rows["protein_change"].str.startswith("p.")].copy()
            for protein_change, variant_rows in exact_rows.groupby("protein_change", sort=True):
                exact_ids = set(variant_rows["model_id"]) & sequenced_crispr_ids
                if len(exact_ids) < min_discovery_n:
                    continue
                other_gene_ids = altered_ids - exact_ids
                gene_wt_ids = sequenced_crispr_ids - altered_ids

                rows.append(
                    _candidate_row(
                        disease_group=disease_models,
                        disease_key=str(disease_key),
                        gene=str(gene),
                        alteration=str(protein_change),
                        event_type="exact_protein_change",
                        comparison_mode="vs_gene_wildtype",
                        context_ids=exact_ids,
                        comparator_ids=gene_wt_ids,
                        sequenced_crispr_ids=sequenced_crispr_ids,
                        evidence=variant_rows,
                        min_context_n=min_context_n,
                        min_comparator_n=min_comparator_n,
                    )
                )
                if other_gene_ids:
                    rows.append(
                        _candidate_row(
                            disease_group=disease_models,
                            disease_key=str(disease_key),
                            gene=str(gene),
                            alteration=str(protein_change),
                            event_type="exact_protein_change",
                            comparison_mode="vs_other_gene_variants",
                            context_ids=exact_ids,
                            comparator_ids=other_gene_ids,
                            sequenced_crispr_ids=sequenced_crispr_ids,
                            evidence=variant_rows,
                            min_context_n=min_context_n,
                            min_comparator_n=min_comparator_n,
                        )
                    )

    if not rows:
        return pd.DataFrame()

    frame = pd.DataFrame(rows).drop_duplicates("candidate_id")
    frame["_ready_order"] = frame["readiness_status"].map(_READY_ORDER).fillna(99)
    frame["_event_order"] = frame["event_type"].map({"exact_protein_change": 0, "gene_altered": 1}).fillna(9)
    frame = frame.sort_values(
        ["_ready_order", "_event_order", "context_models_n", "comparator_models_n", "disease_key", "alteration_gene", "alteration"],
        ascending=[True, True, False, False, True, True, True],
    ).drop(columns=["_ready_order", "_event_order"])
    return frame.reset_index(drop=True)
