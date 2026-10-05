from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

try:
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
except ImportError as exc:
    raise SystemExit(
        'RDKit is required. Install chemistry extras: '
        '.\\.venv\\Scripts\\python.exe -m pip install -e ".[chemistry]"'
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
TLS = ROOT / "data" / "runtime" / "target_ligand_space"
PHARM = ROOT / "data" / "runtime" / "pharmacology"
OWN = ROOT / "data" / "runtime" / "own_compounds"
PROCESSED = ROOT / "data" / "processed"
RUNTIME = ROOT / "data" / "runtime" / "cellular_evidence"
QC_DIR = ROOT / "outputs" / "qc"

LIGAND_CATALOG = TLS / "ligand_catalog.parquet"
PHARM_COMPOUNDS = PHARM / "compounds.parquet"
PHARM_RESPONSES = PHARM / "responses.parquet"
TARGET_CONTEXT = PROCESSED / "depmap_model_target_context.parquet"
ATLAS = PROCESSED / "depmap_crispr_model_atlas.parquet"
CANDIDATES = PHARM / "candidate_hypotheses_v2.parquet"
PYZ_TOP_HITS = OWN / "target_ligand_space_similarity_top_hits.parquet"

STRUCTURE_MAP_OUT = RUNTIME / "target_ligand_prism_structure_map.parquet"
STRUCTURE_MAP_TSV_OUT = RUNTIME / "target_ligand_prism_structure_map.tsv"
MODEL_OUT = RUNTIME / "target_ligand_cellular_models.parquet"
SUMMARY_OUT = RUNTIME / "target_ligand_cellular_summary.parquet"
SUMMARY_TSV_OUT = RUNTIME / "target_ligand_cellular_summary.tsv"
PYZ_HITS_OUT = RUNTIME / "pyz_cellular_evidence_top_hits.parquet"
PYZ_HITS_TSV_OUT = RUNTIME / "pyz_cellular_evidence_top_hits.tsv"
PYZ_SUMMARY_OUT = RUNTIME / "pyz_cellular_evidence_summary.parquet"
PYZ_SUMMARY_TSV_OUT = RUNTIME / "pyz_cellular_evidence_summary.tsv"
MANIFEST_OUT = RUNTIME / "cellular_evidence_manifest.json"
QC_OUT = QC_DIR / "cellular_evidence_qc.tsv"

PRISM_ACTIVE_LFC_THRESHOLD = -1.0
SELECTIVE_PERCENTILE = 0.75
NONRESPONSIVE_PERCENTILE = 0.25
DEPENDENCY_PROBABILITY_THRESHOLD = 0.5
MIN_PROFILE_MODELS = 20
MIN_LINEAGE_MODELS = 3

EVIDENCE_ORDER = {
    "strong_cross_modal_support": 0,
    "supportive_cross_modal": 1,
    "cell_activity_but_dependency_discordant": 2,
    "cell_activity_only": 3,
    "dependency_without_cell_activity": 4,
    "no_cellular_support": 5,
    "sparse": 6,
}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _json_unique(values: Iterable[Any]) -> str:
    out: list[str] = []
    for value in values:
        text = _text(value)
        if text and text not in out:
            out.append(text)
    return json.dumps(out, ensure_ascii=False)


def _standardize_smiles(smiles: Any) -> tuple[str | None, str | None, str | None]:
    text = _text(smiles)
    if not text:
        return None, None, "missing_smiles"
    try:
        mol = Chem.MolFromSmiles(text, sanitize=True)
        if mol is None:
            return None, None, "rdkit_parse_failed"
        cleaned = rdMolStandardize.Cleanup(mol)
        parent = rdMolStandardize.FragmentParent(cleaned)
        canonical = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
        repaired = Chem.MolFromSmiles(canonical, sanitize=True)
        if repaired is None:
            return None, None, "post_standardization_reparse_failed"
        repaired.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(repaired)
        Chem.GetSymmSSSR(repaired)
        canonical = Chem.MolToSmiles(repaired, canonical=True, isomericSmiles=True)
        inchikey = Chem.MolToInchiKey(repaired)
        if not inchikey:
            return None, None, "missing_inchikey"
        return canonical, inchikey, None
    except Exception as exc:
        return None, None, f"standardization_failed:{type(exc).__name__}"


def _arrow_filtered(path: Path, columns: list[str], expression: Any | None = None) -> pd.DataFrame:
    dataset = ds.dataset(path, format="parquet")
    available = set(dataset.schema.names)
    keep = [c for c in columns if c in available]
    if not keep:
        return pd.DataFrame()
    table = dataset.to_table(columns=keep, filter=expression)
    return table.to_pandas()


def _relative_sensitivity(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    n = int(numeric.notna().sum())
    if n <= 1:
        return pd.Series(1.0, index=values.index, dtype=float)
    ranks = numeric.rank(method="min", ascending=True)
    return 1.0 - ((ranks - 1.0) / float(n - 1))


def _spearman(x: pd.Series, y: pd.Series) -> float | None:
    pair = pd.DataFrame(
        {"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")}
    ).dropna()
    if len(pair) < 3 or pair["x"].nunique() < 2 or pair["y"].nunique() < 2:
        return None
    rho = pair["x"].rank(method="average").corr(
        pair["y"].rank(method="average"), method="pearson"
    )
    return None if pd.isna(rho) else float(rho)


def _lineage_adjusted_spearman(frame: pd.DataFrame) -> tuple[float | None, int, int]:
    needed = frame[["gene_effect", "response_value", "oncotree_lineage"]].copy()
    needed["gene_effect"] = pd.to_numeric(needed["gene_effect"], errors="coerce")
    needed["response_value"] = pd.to_numeric(needed["response_value"], errors="coerce")
    needed["oncotree_lineage"] = needed["oncotree_lineage"].fillna("").astype(str)
    needed = needed.dropna(subset=["gene_effect", "response_value"])
    counts = needed["oncotree_lineage"].value_counts()
    valid = set(counts[counts >= MIN_LINEAGE_MODELS].index) - {""}
    needed = needed[needed["oncotree_lineage"].isin(valid)].copy()
    if len(needed) < MIN_PROFILE_MODELS or len(valid) < 2:
        return None, int(len(needed)), int(len(valid))

    needed["rank_ge"] = needed["gene_effect"].rank(method="average")
    needed["rank_response"] = needed["response_value"].rank(method="average")
    needed["ge_residual"] = (
        needed["rank_ge"]
        - needed.groupby("oncotree_lineage")["rank_ge"].transform("mean")
    )
    needed["response_residual"] = (
        needed["rank_response"]
        - needed.groupby("oncotree_lineage")["rank_response"].transform("mean")
    )
    if needed["ge_residual"].nunique() < 2 or needed["response_residual"].nunique() < 2:
        return None, int(len(needed)), int(len(valid))
    rho = needed["ge_residual"].corr(needed["response_residual"], method="pearson")
    return None if pd.isna(rho) else float(rho), int(len(needed)), int(len(valid))


def _direct_potency_band(value: Any) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "unknown"
    if not math.isfinite(number) or number <= 0:
        return "unknown"
    if number <= 100:
        return "very_strong_le_100nM"
    if number <= 1000:
        return "strong_100_1000nM"
    if number <= 10000:
        return "moderate_1_10uM"
    return "outside_tls_cutoff"


def _cellular_label(row: dict[str, Any]) -> str:
    models_n = int(row.get("models_n") or 0)
    sensitive_n = int(row.get("priority_sensitive_models_n") or 0)
    dep_measured_n = int(row.get("dependency_measured_models_n") or 0)
    dependency_n = int(row.get("dependency_models_n") or 0)
    joint_n = int(row.get("joint_support_models_n") or 0)
    rho = row.get("primary_rho")
    delta = row.get("median_response_delta_dependent_minus_other")

    if models_n < 5:
        return "sparse"
    if (
        models_n >= MIN_PROFILE_MODELS
        and sensitive_n >= 2
        and joint_n >= 2
        and rho is not None
        and not pd.isna(rho)
        and float(rho) >= 0.20
        and delta is not None
        and not pd.isna(delta)
        and float(delta) < 0
    ):
        return "strong_cross_modal_support"
    if (
        sensitive_n >= 1
        and joint_n >= 1
        and delta is not None
        and not pd.isna(delta)
        and float(delta) < 0
    ):
        return "supportive_cross_modal"
    if sensitive_n >= 1 and dep_measured_n >= 5:
        discordant_rho = rho is not None and not pd.isna(rho) and float(rho) <= -0.20
        discordant_delta = delta is not None and not pd.isna(delta) and float(delta) >= 0
        if discordant_rho or discordant_delta:
            return "cell_activity_but_dependency_discordant"
    if sensitive_n >= 1:
        return "cell_activity_only"
    if dependency_n >= 1:
        return "dependency_without_cell_activity"
    return "no_cellular_support"


def _build_structure_map(
    ligand_catalog: pd.DataFrame,
    compounds: pd.DataFrame,
    qc: list[dict[str, Any]],
) -> pd.DataFrame:
    prism_rows: list[dict[str, Any]] = []
    for row in compounds.to_dict("records"):
        compound_id = _text(row.get("compound_id"))
        canonical, inchikey, error = _standardize_smiles(row.get("canonical_smiles"))
        if error:
            qc.append(
                {
                    "stage": "prism_structure_standardization",
                    "entity_id": compound_id,
                    "status": "excluded",
                    "detail": error,
                }
            )
            continue
        prism_rows.append(
            {
                "prism_compound_id": compound_id,
                "prism_preferred_name": _text(row.get("preferred_name")),
                "prism_broad_id": _text(row.get("broad_id")),
                "prism_standardized_smiles": canonical,
                "inchikey": inchikey,
            }
        )
    prism = pd.DataFrame(prism_rows)
    if prism.empty:
        return pd.DataFrame()

    ligands = ligand_catalog.copy()
    ligands["inchikey"] = ligands["inchikey"].map(_text)
    ligands = ligands[ligands["inchikey"] != ""].copy()
    mapped = ligands.merge(prism, on="inchikey", how="inner")
    if mapped.empty:
        return mapped
    mapped["match_basis"] = "standardized_inchikey_exact"
    keep = [
        "target_gene",
        "ligand_id",
        "inchikey",
        "canonical_smiles",
        "sources_json",
        "source_ligand_ids_json",
        "endpoint_types_json",
        "measurements_n",
        "lowest_reported_value_nm_across_endpoints",
        "lowest_value_endpoint_type",
        "lowest_value_relation",
        "lowest_value_source",
        "prism_compound_id",
        "prism_preferred_name",
        "prism_broad_id",
        "prism_standardized_smiles",
        "match_basis",
    ]
    return mapped[[c for c in keep if c in mapped.columns]].drop_duplicates().reset_index(drop=True)


def _priority_cancer_pairs() -> set[tuple[str, str]]:
    if not CANDIDATES.exists():
        return set()
    candidates = pd.read_parquet(CANDIDATES)
    required = {"priority_status", "target_gene", "mcl_cancer_id"}
    if not required.issubset(candidates.columns):
        return set()
    priority = candidates[candidates["priority_status"].astype(str).eq("priority_for_in_vitro")].copy()
    return {
        (str(gene).upper().strip(), str(cancer).strip())
        for gene, cancer in zip(priority["target_gene"], priority["mcl_cancer_id"])
        if _text(gene) and _text(cancer)
    }


def _build_model_layer(structure_map: pd.DataFrame) -> pd.DataFrame:
    prism_ids = sorted(set(structure_map["prism_compound_id"].astype(str)))
    target_genes = sorted(set(structure_map["target_gene"].astype(str).str.upper()))

    response_filter = (
        (ds.field("source") == "PRISM")
        & (ds.field("endpoint") == "LFC")
        & ds.field("compound_id").isin(prism_ids)
    )
    responses = _arrow_filtered(
        PHARM_RESPONSES,
        [
            "model_id",
            "compound_id",
            "value",
            "source_release",
            "dose",
            "dose_unit",
            "exposure_time_h",
        ],
        response_filter,
    )
    if responses.empty:
        return pd.DataFrame()
    responses["model_id"] = responses["model_id"].astype(str)
    responses["compound_id"] = responses["compound_id"].astype(str)
    responses["value"] = pd.to_numeric(responses["value"], errors="coerce")
    responses = responses.dropna(subset=["value"])

    mapping = structure_map[
        [
            "target_gene",
            "ligand_id",
            "prism_compound_id",
            "prism_preferred_name",
            "prism_broad_id",
        ]
    ].copy()
    mapping["target_gene"] = mapping["target_gene"].astype(str).str.upper()
    joined = responses.merge(
        mapping,
        left_on="compound_id",
        right_on="prism_compound_id",
        how="inner",
    )
    if joined.empty:
        return pd.DataFrame()

    model_rows = joined.groupby(
        ["target_gene", "ligand_id", "model_id"], as_index=False
    ).agg(
        response_value=("value", "mean"),
        prism_profiles_n=("compound_id", "nunique"),
        prism_compound_ids_json=("compound_id", _json_unique),
        prism_names_json=("prism_preferred_name", _json_unique),
        prism_broad_ids_json=("prism_broad_id", _json_unique),
        source_release=("source_release", lambda x: next((_text(v) for v in x if _text(v)), "")),
        dose=("dose", "median"),
        exposure_time_h=("exposure_time_h", "median"),
    )

    model_rows["relative_sensitivity"] = model_rows.groupby(
        ["target_gene", "ligand_id"], group_keys=False
    )["response_value"].transform(_relative_sensitivity)
    model_rows["absolute_active"] = model_rows["response_value"] <= PRISM_ACTIVE_LFC_THRESHOLD
    model_rows["relatively_sensitive"] = model_rows["relative_sensitivity"] >= SELECTIVE_PERCENTILE
    model_rows["priority_sensitive"] = (
        model_rows["absolute_active"] & model_rows["relatively_sensitive"]
    )
    model_rows["clearly_nonresponsive"] = (
        (model_rows["response_value"] > PRISM_ACTIVE_LFC_THRESHOLD)
        & (model_rows["relative_sensitivity"] <= NONRESPONSIVE_PERCENTILE)
    )

    context_filter = ds.field("target_gene").isin(target_genes)
    context = _arrow_filtered(
        TARGET_CONTEXT,
        [
            "model_id",
            "target_gene",
            "gene_effect",
            "dependency_probability",
            "dependency_call",
            "expression_log2_tpm1",
            "copy_number_relative",
            "has_hotspot",
            "has_likely_lof",
            "protein_changes_json",
        ],
        context_filter,
    )
    if not context.empty:
        context["model_id"] = context["model_id"].astype(str)
        context["target_gene"] = context["target_gene"].astype(str).str.upper()
        model_rows = model_rows.merge(context, on=["model_id", "target_gene"], how="left")
    else:
        model_rows["gene_effect"] = np.nan
        model_rows["dependency_probability"] = np.nan
        model_rows["dependency_call"] = pd.NA

    atlas = _arrow_filtered(
        ATLAS,
        [
            "model_id",
            "cell_line_name",
            "mcl_cancer_id",
            "mcl_cancer_name",
            "mcl_organ_ru",
            "mcl_system_ru",
            "oncotree_lineage",
            "oncotree_subtype",
        ],
    )
    if not atlas.empty:
        atlas["model_id"] = atlas["model_id"].astype(str)
        atlas = atlas.drop_duplicates("model_id")
        model_rows = model_rows.merge(atlas, on="model_id", how="left")

    model_rows["dependency_probability"] = pd.to_numeric(
        model_rows.get("dependency_probability"), errors="coerce"
    )
    model_rows["dependency_call"] = model_rows["dependency_probability"].gt(
        DEPENDENCY_PROBABILITY_THRESHOLD
    ).where(model_rows["dependency_probability"].notna(), pd.NA)
    model_rows["joint_support"] = (
        model_rows["priority_sensitive"]
        & model_rows["dependency_call"].fillna(False).astype(bool)
    )

    priority_pairs = _priority_cancer_pairs()
    model_rows["is_priority_cancer_for_target"] = [
        (str(gene).upper(), _text(cancer)) in priority_pairs
        for gene, cancer in zip(
            model_rows["target_gene"],
            model_rows.get("mcl_cancer_id", pd.Series("", index=model_rows.index)),
        )
    ]
    return model_rows.sort_values(
        ["target_gene", "ligand_id", "response_value", "model_id"]
    ).reset_index(drop=True)


def _build_summary(model_rows: pd.DataFrame, structure_map: pd.DataFrame) -> pd.DataFrame:
    meta_rows: list[dict[str, Any]] = []
    for (gene, ligand_id), group in structure_map.groupby(["target_gene", "ligand_id"], sort=True):
        first = group.iloc[0]
        meta_rows.append(
            {
                "target_gene": str(gene),
                "ligand_id": str(ligand_id),
                "inchikey": _text(first.get("inchikey")),
                "canonical_smiles": _text(first.get("canonical_smiles")),
                "direct_sources_json": _text(first.get("sources_json")),
                "direct_endpoint_types_json": _text(first.get("endpoint_types_json")),
                "direct_measurements_n": int(first.get("measurements_n") or 0),
                "lowest_reported_value_nm_across_endpoints": first.get(
                    "lowest_reported_value_nm_across_endpoints"
                ),
                "lowest_value_endpoint_type": _text(first.get("lowest_value_endpoint_type")),
                "lowest_value_source": _text(first.get("lowest_value_source")),
                "direct_potency_band": _direct_potency_band(
                    first.get("lowest_reported_value_nm_across_endpoints")
                ),
                "prism_profiles_n": int(group["prism_compound_id"].nunique()),
                "prism_compound_ids_json": _json_unique(group["prism_compound_id"]),
                "prism_names_json": _json_unique(group["prism_preferred_name"]),
            }
        )
    metadata = pd.DataFrame(meta_rows)

    rows: list[dict[str, Any]] = []
    for (gene, ligand_id), group in model_rows.groupby(["target_gene", "ligand_id"], sort=True):
        response = pd.to_numeric(group["response_value"], errors="coerce")
        dep_prob = pd.to_numeric(group.get("dependency_probability"), errors="coerce")
        gene_effect = pd.to_numeric(group.get("gene_effect"), errors="coerce")
        dep_measured = dep_prob.notna()
        dependent = dep_prob > DEPENDENCY_PROBABILITY_THRESHOLD
        other = dep_prob <= DEPENDENCY_PROBABILITY_THRESHOLD
        priority_sensitive = group["priority_sensitive"].fillna(False).astype(bool)
        absolute_active = group["absolute_active"].fillna(False).astype(bool)
        joint = group["joint_support"].fillna(False).astype(bool)
        priority_cancer = group["is_priority_cancer_for_target"].fillna(False).astype(bool)

        dep_median = float(response[dependent].median()) if dependent.any() else None
        other_median = float(response[other].median()) if other.any() else None
        delta = None if dep_median is None or other_median is None else dep_median - other_median
        raw_rho = _spearman(gene_effect, response)
        adjusted_rho, adjusted_n, lineages_n = _lineage_adjusted_spearman(group)
        primary_rho = adjusted_rho if adjusted_rho is not None else raw_rho
        primary_method = (
            "lineage_fixed_effect_rank_residual"
            if adjusted_rho is not None
            else "raw_spearman"
        )

        record = {
            "target_gene": str(gene),
            "ligand_id": str(ligand_id),
            "models_n": int(len(group)),
            "absolute_active_models_n": int(absolute_active.sum()),
            "priority_sensitive_models_n": int(priority_sensitive.sum()),
            "clearly_nonresponsive_models_n": int(
                group["clearly_nonresponsive"].fillna(False).astype(bool).sum()
            ),
            "dependency_measured_models_n": int(dep_measured.sum()),
            "dependency_models_n": int(dependent.fillna(False).sum()),
            "joint_support_models_n": int(joint.sum()),
            "priority_cancer_models_n": int(priority_cancer.sum()),
            "priority_cancer_sensitive_models_n": int((priority_cancer & priority_sensitive).sum()),
            "priority_cancer_joint_models_n": int((priority_cancer & joint).sum()),
            "median_response_all": float(response.median()) if response.notna().any() else None,
            "median_response_dependent": dep_median,
            "median_response_other": other_median,
            "median_response_delta_dependent_minus_other": delta,
            "spearman_raw_rho": raw_rho,
            "spearman_lineage_adjusted_rho": adjusted_rho,
            "lineage_adjusted_models_n": adjusted_n,
            "lineages_n": lineages_n,
            "primary_rho": primary_rho,
            "primary_rho_method": primary_method,
            "built_at": _now(),
        }
        record["cellular_evidence_label"] = _cellular_label(record)
        rows.append(record)

    summary = pd.DataFrame(rows)
    if summary.empty:
        return metadata
    summary = metadata.merge(summary, on=["target_gene", "ligand_id"], how="left")
    summary["evidence_order"] = summary["cellular_evidence_label"].map(EVIDENCE_ORDER).fillna(99)
    return summary.sort_values(
        ["evidence_order", "joint_support_models_n", "priority_sensitive_models_n", "target_gene"],
        ascending=[True, False, False, True],
    ).drop(columns=["evidence_order"]).reset_index(drop=True)


def _build_pyz_outputs(cellular: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not PYZ_TOP_HITS.exists():
        return pd.DataFrame(), pd.DataFrame()
    hits = pd.read_parquet(PYZ_TOP_HITS)
    hits["target_gene"] = hits["target_gene"].astype(str).str.upper()
    cellular = cellular.copy()
    cellular["target_gene"] = cellular["target_gene"].astype(str).str.upper()

    annotated = hits.merge(cellular, on=["target_gene", "ligand_id"], how="left")
    annotated["has_cellular_evidence"] = annotated["models_n"].notna()
    annotated["cellular_evidence_order"] = (
        annotated["cellular_evidence_label"].map(EVIDENCE_ORDER).fillna(99).astype(int)
    )
    annotated = annotated.sort_values(
        [
            "own_compound_id",
            "target_gene",
            "cellular_evidence_order",
            "tanimoto_morgan_r2_2048",
        ],
        ascending=[True, True, True, False],
        na_position="last",
    ).reset_index(drop=True)

    summary_rows: list[dict[str, Any]] = []
    for (own_id, gene), group in annotated.groupby(["own_compound_id", "target_gene"], sort=True):
        cellular_group = group[group["has_cellular_evidence"]].copy()
        best = cellular_group.iloc[0] if not cellular_group.empty else None
        counts = group["cellular_evidence_label"].value_counts(dropna=True).to_dict()
        record: dict[str, Any] = {
            "own_compound_id": str(own_id),
            "target_gene": str(gene),
            "saved_structural_neighbors_n": int(len(group)),
            "neighbors_with_cellular_evidence_n": int(len(cellular_group)),
            "strong_cross_modal_neighbors_n": int(counts.get("strong_cross_modal_support", 0)),
            "supportive_cross_modal_neighbors_n": int(counts.get("supportive_cross_modal", 0)),
            "cell_activity_only_neighbors_n": int(counts.get("cell_activity_only", 0)),
            "discordant_neighbors_n": int(
                counts.get("cell_activity_but_dependency_discordant", 0)
            ),
            "best_cellular_ligand_id": "",
            "best_cellular_tanimoto": None,
            "best_cellular_evidence_label": "",
            "best_cellular_direct_value_nm": None,
            "best_cellular_direct_endpoint": "",
            "best_cellular_priority_sensitive_models_n": 0,
            "best_cellular_joint_support_models_n": 0,
            "best_cellular_primary_rho": None,
            "best_cellular_prism_names_json": "[]",
            "built_at": _now(),
        }
        if best is not None:
            record.update(
                {
                    "best_cellular_ligand_id": _text(best.get("ligand_id")),
                    "best_cellular_tanimoto": best.get("tanimoto_morgan_r2_2048"),
                    "best_cellular_evidence_label": _text(best.get("cellular_evidence_label")),
                    "best_cellular_direct_value_nm": best.get(
                        "lowest_reported_value_nm_across_endpoints_y"
                    )
                    if "lowest_reported_value_nm_across_endpoints_y" in best.index
                    else best.get("lowest_reported_value_nm_across_endpoints"),
                    "best_cellular_direct_endpoint": _text(best.get("lowest_value_endpoint_type_y"))
                    if "lowest_value_endpoint_type_y" in best.index
                    else _text(best.get("lowest_value_endpoint_type")),
                    "best_cellular_priority_sensitive_models_n": int(
                        best.get("priority_sensitive_models_n") or 0
                    ),
                    "best_cellular_joint_support_models_n": int(
                        best.get("joint_support_models_n") or 0
                    ),
                    "best_cellular_primary_rho": best.get("primary_rho"),
                    "best_cellular_prism_names_json": _text(best.get("prism_names_json")),
                }
            )
        summary_rows.append(record)

    summary = pd.DataFrame(summary_rows).sort_values(
        ["own_compound_id", "neighbors_with_cellular_evidence_n", "target_gene"],
        ascending=[True, False, True],
    ).reset_index(drop=True)
    return annotated, summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build MCL Cellular Evidence Layer v1 by exact standardized structure mapping "
            "between Target Ligand Space and PRISM, then join DepMap target dependency."
        )
    )
    _ = parser.parse_args()

    required = [LIGAND_CATALOG, PHARM_COMPOUNDS, PHARM_RESPONSES, TARGET_CONTEXT, ATLAS]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    if missing:
        raise SystemExit(
            "Missing Cellular Evidence Layer inputs:\n"
            + "\n".join(missing)
            + "\nBuild Target Ligand Space, PRISM pharmacology and Candidate v2 first."
        )

    RUNTIME.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    qc: list[dict[str, Any]] = []

    ligand_catalog = pd.read_parquet(LIGAND_CATALOG)
    compounds = pd.read_parquet(PHARM_COMPOUNDS)
    structure_map = _build_structure_map(ligand_catalog, compounds, qc)
    if structure_map.empty:
        raise SystemExit(
            "No exact standardized structures overlap between Target Ligand Space and PRISM compounds."
        )

    structure_map.to_parquet(STRUCTURE_MAP_OUT, index=False, compression="zstd")
    structure_map.to_csv(STRUCTURE_MAP_TSV_OUT, sep="\t", index=False)

    model_rows = _build_model_layer(structure_map)
    if model_rows.empty:
        raise SystemExit("Exact structure matches were found, but no PRISM LFC observations mapped to them.")
    model_rows.to_parquet(MODEL_OUT, index=False, compression="zstd")

    cellular = _build_summary(model_rows, structure_map)
    cellular.to_parquet(SUMMARY_OUT, index=False, compression="zstd")
    cellular.to_csv(SUMMARY_TSV_OUT, sep="\t", index=False)

    pyz_hits, pyz_summary = _build_pyz_outputs(cellular)
    if not pyz_hits.empty:
        pyz_hits.to_parquet(PYZ_HITS_OUT, index=False, compression="zstd")
        pyz_hits.to_csv(PYZ_HITS_TSV_OUT, sep="\t", index=False)
    if not pyz_summary.empty:
        pyz_summary.to_parquet(PYZ_SUMMARY_OUT, index=False, compression="zstd")
        pyz_summary.to_csv(PYZ_SUMMARY_TSV_OUT, sep="\t", index=False)

    pd.DataFrame(qc, columns=["stage", "entity_id", "status", "detail"]).to_csv(
        QC_OUT, sep="\t", index=False
    )

    label_counts = cellular["cellular_evidence_label"].value_counts(dropna=False).to_dict()
    manifest = {
        "contract": "mcl-cellular-evidence-v1",
        "built_at": _now(),
        "target_ligand_structures_n": int(len(ligand_catalog)),
        "exact_target_ligand_prism_map_rows_n": int(len(structure_map)),
        "target_ligands_with_prism_structure_n": int(structure_map["ligand_id"].nunique()),
        "prism_profiles_mapped_n": int(structure_map["prism_compound_id"].nunique()),
        "model_target_ligand_rows_n": int(len(model_rows)),
        "target_ligand_cellular_summaries_n": int(len(cellular)),
        "cellular_evidence_label_counts": {str(k): int(v) for k, v in label_counts.items()},
        "prism_contract": {
            "release": "PRISM Repurposing Public 24Q2",
            "endpoint": "LFC",
            "dose_uM": 2.5,
            "exposure_time_h": 120,
            "absolute_active": f"LFC <= {PRISM_ACTIVE_LFC_THRESHOLD:g}",
            "relative_selectivity": f"relative_sensitivity >= {SELECTIVE_PERCENTILE:g}",
            "priority_sensitive": "absolute_active AND relative_selectivity",
        },
        "dependency_contract": {
            "binary_call": f"Probability of Dependency > {DEPENDENCY_PROBABILITY_THRESHOLD:g}",
            "gene_effect": "continuous Chronos Gene Effect retained for profile correlation",
        },
        "structure_mapping": (
            "RDKit Cleanup -> FragmentParent -> canonical isomeric SMILES -> InChIKey; "
            "Target Ligand Space and PRISM are linked only by exact standardized InChIKey."
        ),
        "scientific_guardrails_ru": [
            "Совпадение структуры известного лиганда с профилем PRISM не доказывает, что клеточный эффект вызван именно заявленной мишенью.",
            "PRISM primary — одноточечный клеточный скрининг; он не заменяет dose-response IC50/GI50 и target-engagement эксперимент.",
            "CRISPR knockout не эквивалентен фармакологическому ингибированию.",
            "Структурное сходство PYZ с известным лигандом не доказывает связывание PYZ с той же мишенью.",
            "Ki, Kd и IC50 прямого белкового опыта сохраняются раздельно и не считаются взаимозаменяемыми.",
        ],
    }
    MANIFEST_OUT.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print("MCL Cellular Evidence Layer v1")
    print(f"Target Ligand Space structures: {len(ligand_catalog)}")
    print(
        "Exact Target Ligand Space -> PRISM structure map: "
        f"{structure_map['ligand_id'].nunique()} target-ligands; "
        f"{structure_map['prism_compound_id'].nunique()} PRISM profiles"
    )
    print(f"Model x target-ligand rows: {len(model_rows)}")
    print(f"Target-ligand cellular summaries: {len(cellular)}")
    print("Cellular evidence labels:")
    for label, count in sorted(label_counts.items(), key=lambda x: EVIDENCE_ORDER.get(str(x[0]), 99)):
        print(f"  {label}: {count}")
    if not pyz_summary.empty:
        covered = int((pyz_summary["neighbors_with_cellular_evidence_n"] > 0).sum())
        print(f"PYZ x target pairs with at least one saved cellular neighbor: {covered}/{len(pyz_summary)}")
        best = pyz_hits[pyz_hits["has_cellular_evidence"]].copy()
        if not best.empty:
            best = best.sort_values(
                ["cellular_evidence_order", "tanimoto_morgan_r2_2048"],
                ascending=[True, False],
            ).head(25)
            print("Best PYZ structural neighbors with cellular evidence:")
            for row in best.to_dict("records"):
                print(
                    f"  {row['own_compound_id']} -> {row['target_gene']} -> {row['ligand_id']}: "
                    f"Tanimoto={float(row['tanimoto_morgan_r2_2048']):.3f}; "
                    f"cellular={row['cellular_evidence_label']}; "
                    f"PRISM-sensitive={int(row.get('priority_sensitive_models_n') or 0)}; "
                    f"joint CRISPR+PRISM={int(row.get('joint_support_models_n') or 0)}"
                )
    print(f"Primary ligand summary: {SUMMARY_TSV_OUT.relative_to(ROOT)}")
    print(f"Model-level evidence: {MODEL_OUT.relative_to(ROOT)}")
    if PYZ_SUMMARY_OUT.exists():
        print(f"PYZ summary: {PYZ_SUMMARY_TSV_OUT.relative_to(ROOT)}")
    if PYZ_HITS_OUT.exists():
        print(f"PYZ annotated neighbors: {PYZ_HITS_TSV_OUT.relative_to(ROOT)}")
    print(f"QC: {QC_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
