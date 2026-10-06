from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OWN = ROOT / "data" / "runtime" / "own_compounds"
CELL = ROOT / "data" / "runtime" / "cellular_evidence"
PHARM = ROOT / "data" / "runtime" / "pharmacology"
PRECLIN = ROOT / "data" / "runtime" / "preclinical"
OUT = ROOT / "data" / "runtime" / "pyz_target_evidence"
QC_DIR = ROOT / "outputs" / "qc"

TLS_SUMMARY = OWN / "target_ligand_space_similarity_summary.parquet"
TLS_TOP = OWN / "target_ligand_space_similarity_top_hits.parquet"
CELL_LIGAND = CELL / "target_ligand_cellular_summary.parquet"
CELL_PYZ = CELL / "pyz_cellular_evidence_summary.parquet"
STRICT_ONCOLOGY = OWN / "oncology_reference_target_similarity_summary.parquet"
CANDIDATES = PHARM / "candidate_hypotheses_v2.parquet"
PRECLIN_TARGET = PRECLIN / "preclinical_target_summary.parquet"

MATRIX_PARQUET = OUT / "pyz_target_evidence_matrix.parquet"
MATRIX_TSV = OUT / "pyz_target_evidence_matrix.tsv"
SHORTLIST_PARQUET = OUT / "pyz_target_validation_shortlist.parquet"
SHORTLIST_TSV = OUT / "pyz_target_validation_shortlist.tsv"
MANIFEST = OUT / "pyz_target_evidence_manifest.json"
QC_OUT = QC_DIR / "pyz_target_evidence_matrix_qc.tsv"

EXPECTED_PYZ = 65
EXPECTED_TARGETS = 13
EXPECTED_ROWS = EXPECTED_PYZ * EXPECTED_TARGETS


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


def _num(value: Any) -> float | None:
    try:
        if value is None or pd.isna(value):
            return None
        x = float(value)
        return x if np.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def _json_unique(values: pd.Series) -> str:
    out: list[str] = []
    for value in values:
        text = _text(value)
        if text and text not in out:
            out.append(text)
    return json.dumps(out, ensure_ascii=False)


def _normalise_keys(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    if "own_compound_id" in out.columns:
        out["own_compound_id"] = out["own_compound_id"].astype(str).str.strip()
    if "target_gene" in out.columns:
        out["target_gene"] = out["target_gene"].astype(str).str.upper().str.strip()
    return out


def _transfer_class(value: Any) -> str:
    x = _num(value)
    if x is None:
        return "no_similarity_result"
    if x >= 0.70:
        return "high_structural_transferability"
    if x >= 0.50:
        return "moderate_structural_transferability"
    if x >= 0.35:
        return "exploratory_low_structural_similarity"
    return "mechanism_not_transferable_by_2d_similarity"


def _transfer_ru(value: Any) -> str:
    cls = _transfer_class(value)
    return {
        "high_structural_transferability": "Высокое 2D-сходство; перенос механистической гипотезы допустим только как проверяемая гипотеза.",
        "moderate_structural_transferability": "Умеренное 2D-сходство; мишень заслуживает направленной проверки, но механизм не доказан.",
        "exploratory_low_structural_similarity": "Низкое, но заметное 2D-сходство; только исследовательская гипотеза для проверки.",
        "mechanism_not_transferable_by_2d_similarity": "Сходство <0,35: механизм известного лиганда нельзя переносить на PYZ по 2D-структуре.",
        "no_similarity_result": "Нет сопоставимого структурного результата.",
    }[cls]


def _next_step(row: pd.Series) -> str:
    sim = _num(row.get("best_tanimoto_morgan_r2_2048"))
    if sim is None:
        return "Не ранжировать по химическому сходству; проверить входные структуры."
    if sim >= 0.50:
        return (
            "Направленная проверка мишени: ортогональный ligand-based метод, docking при пригодной структуре белка и прямой биохимический/биофизический assay."
        )
    if sim >= 0.35:
        return (
            "Исследовательская проверка мишени: docking/ligand-based prediction допустимы для приоритизации, но до прямого assay нельзя утверждать связывание или механизм."
        )
    return (
        "Не переносить механизм ближайшего лиганда на PYZ. Для поиска мишени использовать независимые методы: target prediction/ML, docking по биологически обоснованным мишеням и фенотипические эксперименты."
    )


def _candidate_target_context() -> pd.DataFrame:
    if not CANDIDATES.exists():
        return pd.DataFrame(columns=["target_gene"])
    frame = pd.read_parquet(CANDIDATES)
    required = {"priority_status", "target_gene"}
    if not required.issubset(frame.columns):
        return pd.DataFrame(columns=["target_gene"])
    frame = frame[frame["priority_status"].astype(str).eq("priority_for_in_vitro")].copy()
    if frame.empty:
        return pd.DataFrame(columns=["target_gene"])
    frame["target_gene"] = frame["target_gene"].astype(str).str.upper().str.strip()
    rows: list[dict[str, Any]] = []
    for gene, group in frame.groupby("target_gene", sort=True):
        cancer_id_col = "mcl_cancer_id" if "mcl_cancer_id" in group.columns else None
        cancer_name_col = "mcl_cancer_name" if "mcl_cancer_name" in group.columns else None
        rows.append(
            {
                "target_gene": gene,
                "mcl_priority_hypotheses_n": int(len(group)),
                "mcl_priority_cancers_n": int(group[cancer_id_col].replace("", pd.NA).nunique()) if cancer_id_col else 0,
                "mcl_priority_cancer_ids_json": _json_unique(group[cancer_id_col]) if cancer_id_col else "[]",
                "mcl_priority_cancer_names_json": _json_unique(group[cancer_name_col]) if cancer_name_col else "[]",
            }
        )
    return pd.DataFrame(rows)


def _best_threshold_neighbors(top_hits: pd.DataFrame) -> pd.DataFrame:
    if top_hits.empty:
        return pd.DataFrame(columns=["own_compound_id", "target_gene"])
    top_hits = _normalise_keys(top_hits)
    potency_col = "lowest_reported_value_nm_across_endpoints"
    if potency_col not in top_hits.columns:
        return pd.DataFrame(columns=["own_compound_id", "target_gene"])
    top_hits[potency_col] = pd.to_numeric(top_hits[potency_col], errors="coerce")
    top_hits["tanimoto_morgan_r2_2048"] = pd.to_numeric(
        top_hits["tanimoto_morgan_r2_2048"], errors="coerce"
    )

    base_keys = top_hits[["own_compound_id", "target_gene"]].drop_duplicates().copy()
    result = base_keys
    for threshold, label in ((100.0, "le_100nM"), (1000.0, "le_1000nM")):
        eligible = top_hits[
            top_hits[potency_col].notna() & top_hits[potency_col].le(threshold)
        ].copy()
        if eligible.empty:
            continue
        eligible = eligible.sort_values(
            ["own_compound_id", "target_gene", "tanimoto_morgan_r2_2048"],
            ascending=[True, True, False],
        ).drop_duplicates(["own_compound_id", "target_gene"], keep="first")
        cols = [
            "own_compound_id",
            "target_gene",
            "ligand_id",
            "tanimoto_morgan_r2_2048",
            potency_col,
            "lowest_value_endpoint_type",
            "lowest_value_source",
        ]
        eligible = eligible[[c for c in cols if c in eligible.columns]].copy()
        eligible = eligible.rename(
            columns={
                "ligand_id": f"top20_best_neighbor_{label}_ligand_id",
                "tanimoto_morgan_r2_2048": f"top20_best_neighbor_{label}_tanimoto",
                potency_col: f"top20_best_neighbor_{label}_reported_value_nm",
                "lowest_value_endpoint_type": f"top20_best_neighbor_{label}_endpoint",
                "lowest_value_source": f"top20_best_neighbor_{label}_source",
            }
        )
        result = result.merge(eligible, on=["own_compound_id", "target_gene"], how="left")
    return result


def main() -> None:
    required = [TLS_SUMMARY, TLS_TOP, CELL_LIGAND, CELL_PYZ]
    missing = [str(p.relative_to(ROOT)) for p in required if not p.exists()]
    if missing:
        raise SystemExit(
            "Missing PYZ Target Evidence Matrix inputs:\n"
            + "\n".join(missing)
            + "\nRun Target Ligand Space and Cellular Evidence Layer v1.1 first."
        )

    OUT.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)

    tls = _normalise_keys(pd.read_parquet(TLS_SUMMARY))
    top = _normalise_keys(pd.read_parquet(TLS_TOP))
    cellular_ligand = _normalise_keys(pd.read_parquet(CELL_LIGAND))
    cellular_pyz = _normalise_keys(pd.read_parquet(CELL_PYZ))

    matrix = tls.copy()
    matrix["best_tanimoto_morgan_r2_2048"] = pd.to_numeric(
        matrix["best_tanimoto_morgan_r2_2048"], errors="coerce"
    )

    # Cellular evidence for the SAME ligand that is the best structural neighbor overall.
    same_ligand_cols = [
        "target_gene",
        "ligand_id",
        "cellular_evidence_label",
        "models_n",
        "priority_sensitive_models_n",
        "joint_support_models_n",
        "primary_rho",
        "median_response_delta_dependent_minus_other",
        "prism_names_json",
    ]
    same = cellular_ligand[[c for c in same_ligand_cols if c in cellular_ligand.columns]].copy()
    same = same.rename(
        columns={
            "ligand_id": "best_ligand_id",
            "cellular_evidence_label": "best_structural_ligand_cellular_label",
            "models_n": "best_structural_ligand_prism_models_n",
            "priority_sensitive_models_n": "best_structural_ligand_prism_sensitive_models_n",
            "joint_support_models_n": "best_structural_ligand_joint_crispr_prism_models_n",
            "primary_rho": "best_structural_ligand_crispr_prism_rho",
            "median_response_delta_dependent_minus_other": "best_structural_ligand_dep_vs_other_response_delta",
            "prism_names_json": "best_structural_ligand_prism_names_json",
        }
    )
    same = same.drop_duplicates(["target_gene", "best_ligand_id"], keep="first")
    matrix = matrix.merge(same, on=["target_gene", "best_ligand_id"], how="left")
    matrix["best_structural_ligand_has_cellular_evidence"] = matrix[
        "best_structural_ligand_cellular_label"
    ].notna()

    # Best cellularly characterized ligand among the saved top-20 structural neighbors.
    cell_keep = [
        "own_compound_id",
        "target_gene",
        "neighbors_with_cellular_evidence_n",
        "strong_cross_modal_neighbors_n",
        "supportive_cross_modal_neighbors_n",
        "cell_activity_only_neighbors_n",
        "discordant_neighbors_n",
        "best_cellular_ligand_id",
        "best_cellular_tanimoto",
        "best_cellular_evidence_label",
        "best_cellular_direct_value_nm",
        "best_cellular_direct_endpoint",
        "best_cellular_priority_sensitive_models_n",
        "best_cellular_joint_support_models_n",
        "best_cellular_primary_rho",
        "best_cellular_prism_names_json",
    ]
    cell_join = cellular_pyz[[c for c in cell_keep if c in cellular_pyz.columns]].copy()
    matrix = matrix.merge(cell_join, on=["own_compound_id", "target_gene"], how="left")

    threshold_neighbors = _best_threshold_neighbors(top)
    matrix = matrix.merge(threshold_neighbors, on=["own_compound_id", "target_gene"], how="left")

    # Strict oncology reference comparison, if already built.
    strict_available = STRICT_ONCOLOGY.exists()
    if strict_available:
        strict = _normalise_keys(pd.read_parquet(STRICT_ONCOLOGY))
        strict_keep = [
            "own_compound_id",
            "target_gene",
            "known_comparators_with_structure_n",
            "best_comparator_compound_id",
            "best_comparator_name",
            "best_tanimoto_morgan_r2_2048",
            "best_similarity_band",
        ]
        strict = strict[[c for c in strict_keep if c in strict.columns]].copy()
        strict = strict.rename(
            columns={
                "known_comparators_with_structure_n": "strict_oncology_references_n",
                "best_comparator_compound_id": "strict_oncology_best_reference_id",
                "best_comparator_name": "strict_oncology_best_reference_name",
                "best_tanimoto_morgan_r2_2048": "strict_oncology_best_tanimoto",
                "best_similarity_band": "strict_oncology_best_similarity_band",
            }
        )
        matrix = matrix.merge(strict, on=["own_compound_id", "target_gene"], how="left")

    target_context = _candidate_target_context()
    matrix = matrix.merge(target_context, on="target_gene", how="left")

    preclinical_available = PRECLIN_TARGET.exists()
    if preclinical_available:
        pre = _normalise_keys(pd.read_parquet(PRECLIN_TARGET))
        pre_keep = [
            "target_gene",
            "target_linked_compounds_n",
            "compounds_with_pcdb_evidence_n",
            "compounds_with_oncology_pcdb_evidence_n",
            "compounds_with_exact_priority_cancer_pcdb_evidence_n",
            "priority_cancers_json",
        ]
        pre = pre[[c for c in pre_keep if c in pre.columns]].copy().drop_duplicates("target_gene")
        pre = pre.rename(
            columns={
                "target_linked_compounds_n": "pcdb_target_linked_compounds_n",
                "compounds_with_pcdb_evidence_n": "pcdb_compounds_with_any_evidence_n",
                "compounds_with_oncology_pcdb_evidence_n": "pcdb_compounds_with_oncology_context_n",
                "compounds_with_exact_priority_cancer_pcdb_evidence_n": "pcdb_compounds_with_exact_priority_cancer_context_n",
                "priority_cancers_json": "pcdb_priority_cancers_json",
            }
        )
        matrix = matrix.merge(pre, on="target_gene", how="left")

    matrix["structural_transfer_class"] = matrix["best_tanimoto_morgan_r2_2048"].map(_transfer_class)
    matrix["structural_transfer_interpretation_ru"] = matrix["best_tanimoto_morgan_r2_2048"].map(_transfer_ru)
    matrix["recommended_next_step_ru"] = matrix.apply(_next_step, axis=1)

    best_cell_sim = pd.to_numeric(matrix.get("best_cellular_tanimoto"), errors="coerce")
    matrix["cellular_neighbor_transfer_warning"] = np.where(
        matrix.get("neighbors_with_cellular_evidence_n", 0).fillna(0).astype(int).gt(0)
        & (best_cell_sim.isna() | best_cell_sim.lt(0.35)),
        "Known ligand has cellular/CRISPR evidence, but its 2D similarity to PYZ is <0.35; do not transfer mechanism to PYZ.",
        "",
    )
    matrix["direct_pyz_target_evidence"] = "absent_not_measured"
    matrix["built_at"] = _now()

    # Stable ordering: chemical evidence first, then biological context; no composite 0-100 score.
    matrix = matrix.sort_values(
        ["structural_transfer_class", "best_tanimoto_morgan_r2_2048", "own_compound_id", "target_gene"],
        ascending=[True, False, True, True],
        na_position="last",
    ).reset_index(drop=True)

    shortlist = matrix[matrix["best_tanimoto_morgan_r2_2048"].ge(0.35)].copy()
    shortlist = shortlist.sort_values(
        ["best_tanimoto_morgan_r2_2048", "best_structural_ligand_has_cellular_evidence"],
        ascending=[False, False],
    ).reset_index(drop=True)
    shortlist["shortlist_scope"] = "Tanimoto >= 0.35 only; exploratory target-validation candidates, not validated PYZ targets"

    matrix.to_parquet(MATRIX_PARQUET, index=False, compression="zstd")
    matrix.to_csv(MATRIX_TSV, sep="\t", index=False)
    shortlist.to_parquet(SHORTLIST_PARQUET, index=False, compression="zstd")
    shortlist.to_csv(SHORTLIST_TSV, sep="\t", index=False)

    unique_pyz = int(matrix["own_compound_id"].nunique())
    unique_targets = int(matrix["target_gene"].nunique())
    duplicate_keys = int(matrix.duplicated(["own_compound_id", "target_gene"]).sum())
    qc_rows = [
        {"check": "rows_expected_65x13", "value": len(matrix), "expected": EXPECTED_ROWS, "status": "PASS" if len(matrix) == EXPECTED_ROWS else "WARNING"},
        {"check": "pyz_compounds_n", "value": unique_pyz, "expected": EXPECTED_PYZ, "status": "PASS" if unique_pyz == EXPECTED_PYZ else "WARNING"},
        {"check": "targets_n", "value": unique_targets, "expected": EXPECTED_TARGETS, "status": "PASS" if unique_targets == EXPECTED_TARGETS else "WARNING"},
        {"check": "duplicate_pyz_target_keys", "value": duplicate_keys, "expected": 0, "status": "PASS" if duplicate_keys == 0 else "ERROR"},
        {"check": "direct_pyz_target_measurements_present", "value": 0, "expected": 0, "status": "INFO"},
    ]
    pd.DataFrame(qc_rows).to_csv(QC_OUT, sep="\t", index=False)

    class_counts = matrix["structural_transfer_class"].value_counts(dropna=False).to_dict()
    manifest = {
        "contract": "mcl-pyz-target-evidence-matrix-v1",
        "built_at": _now(),
        "rows_n": int(len(matrix)),
        "pyz_compounds_n": unique_pyz,
        "targets_n": unique_targets,
        "shortlist_rows_n": int(len(shortlist)),
        "strict_oncology_reference_layer_available": strict_available,
        "preclinical_target_context_available": preclinical_available,
        "structural_transfer_class_counts": {str(k): int(v) for k, v in class_counts.items()},
        "no_composite_score": True,
        "transferability_rules": {
            ">=0.70": "high_structural_transferability",
            "0.50-0.70": "moderate_structural_transferability",
            "0.35-0.50": "exploratory_low_structural_similarity",
            "<0.35": "mechanism_not_transferable_by_2d_similarity",
        },
        "scientific_guardrails_ru": [
            "Матрица не содержит прямых измерений связывания PYZ с мишенью; поле direct_pyz_target_evidence поэтому равно absent_not_measured.",
            "Клеточные и CRISPR-данные характеризуют известные экспериментальные лиганды, а не PYZ. Их нельзя переносить на PYZ при слабом структурном сходстве.",
            "Tanimoto является 2D-мерой структурного сходства и не доказывает одинаковую мишень, позу связывания, активность, селективность или противоопухолевый эффект.",
            "Ki, Kd и IC50 не считаются взаимозаменяемыми; числовой минимум сохраняется только как навигационный атрибут вместе с типом endpoint.",
            "PRISM primary — одноточечный клеточный скрининг; он не заменяет dose-response и target engagement.",
            "CRISPR knockout не эквивалентен фармакологическому ингибированию.",
            "PCDB здесь используется только как контекст мишени/известных соединений и не является доказательством эффективности PYZ.",
        ],
        "outputs": [
            str(MATRIX_PARQUET.relative_to(ROOT)),
            str(MATRIX_TSV.relative_to(ROOT)),
            str(SHORTLIST_PARQUET.relative_to(ROOT)),
            str(SHORTLIST_TSV.relative_to(ROOT)),
            str(QC_OUT.relative_to(ROOT)),
        ],
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print("MCL PYZ Target Evidence Matrix v1")
    print(f"Rows: {len(matrix)} = {unique_pyz} PYZ x {unique_targets} targets")
    print("Structural transfer classes:")
    for label, count in class_counts.items():
        print(f"  {label}: {count}")
    print(f"Exploratory shortlist (Tanimoto >= 0.35): {len(shortlist)}")
    if not shortlist.empty:
        print("Shortlist:")
        for row in shortlist.head(30).to_dict("records"):
            direct_value = row.get("best_lowest_reported_value_nm_across_endpoints")
            endpoint = _text(row.get("best_lowest_value_endpoint_type"))
            direct = (
                f"; known ligand {endpoint}={float(direct_value):g} nM"
                if endpoint and _num(direct_value) is not None
                else ""
            )
            same_cell = _text(row.get("best_structural_ligand_cellular_label")) or "no_PRISM_match_for_same_best_ligand"
            print(
                f"  {row['own_compound_id']} -> {row['target_gene']}: "
                f"Tanimoto={float(row['best_tanimoto_morgan_r2_2048']):.3f}; "
                f"best ligand={row['best_ligand_id']}{direct}; same-ligand cellular={same_cell}"
            )
    print(f"Matrix: {MATRIX_TSV.relative_to(ROOT)}")
    print(f"Shortlist: {SHORTLIST_TSV.relative_to(ROOT)}")
    print(f"QC: {QC_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
