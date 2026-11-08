from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "oncology_reference_compounds.tsv"
PHARM = ROOT / "data" / "runtime" / "pharmacology"
PRECLIN = ROOT / "data" / "runtime" / "preclinical"
QC_DIR = ROOT / "outputs" / "qc"

CANDIDATES = PHARM / "candidate_hypotheses_v2.parquet"
TARGET_COMPOUNDS = PHARM / "target_compound_catalog.parquet"
COMPOUNDS = PHARM / "compound_catalog.parquet"
PCDB_SUMMARY = PRECLIN / "preclinical_compound_summary.tsv"

OUTPUT = PHARM / "oncology_reference_compounds.parquet"
OUTPUT_TSV = PHARM / "oncology_reference_compounds.tsv"
TARGET_SUMMARY = PHARM / "oncology_reference_target_summary.tsv"
QC_OUTPUT = QC_DIR / "oncology_reference_compounds_qc.tsv"
MANIFEST = PHARM / "oncology_reference_compounds_manifest.json"

PRIORITY_STATUS = "priority_for_in_vitro"
ALLOWED_LEVELS = {"A", "B", "C", "D"}


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


def _norm(value: Any) -> str:
    return re.sub(r"[^A-Z0-9]+", "", _text(value).upper())


def _keys(row: pd.Series) -> set[str]:
    values = [_text(row.get("compound_name"))]
    values += [x.strip() for x in _text(row.get("aliases")).split("|") if x.strip()]
    return {_norm(x) for x in values if _norm(x)}


def main() -> None:
    required = [CONFIG, CANDIDATES, TARGET_COMPOUNDS, COMPOUNDS]
    missing = [str(p.relative_to(ROOT)) for p in required if not p.exists()]
    if missing:
        raise SystemExit("Missing oncology-reference inputs:\n" + "\n".join(missing))

    refs = pd.read_csv(CONFIG, sep="\t", dtype=str, keep_default_na=False)
    candidates = pd.read_parquet(CANDIDATES)
    target_compounds = pd.read_parquet(TARGET_COMPOUNDS)
    compounds = pd.read_parquet(COMPOUNDS)

    required_cols = {
        "reference_id", "target_gene", "compound_name", "aliases", "reference_level",
        "modality", "target_relation", "action", "include_similarity", "oncology_rationale_ru",
        "source_kind", "source_url",
    }
    missing_cols = sorted(required_cols - set(refs.columns))
    if missing_cols:
        raise SystemExit("oncology_reference_compounds.tsv missing columns: " + ", ".join(missing_cols))
    if refs["reference_id"].duplicated().any():
        dup = refs.loc[refs["reference_id"].duplicated(), "reference_id"].tolist()
        raise SystemExit(f"Duplicate reference_id values: {dup}")
    bad_levels = sorted(set(refs["reference_level"]) - ALLOWED_LEVELS)
    if bad_levels:
        raise SystemExit(f"Unsupported reference levels: {bad_levels}")

    refs = refs.copy()
    refs["target_gene"] = refs["target_gene"].str.upper().str.strip()
    refs["include_similarity"] = refs["include_similarity"].astype(str).isin({"1", "true", "True", "yes"})
    refs["curation_keys"] = refs.apply(lambda r: sorted(_keys(r)), axis=1)

    candidates = candidates.copy()
    candidates["target_gene"] = candidates["target_gene"].astype(str).str.upper()
    priority = candidates[candidates["priority_status"].astype(str).eq(PRIORITY_STATUS)].copy()
    priority_targets = sorted(set(priority["target_gene"].dropna().astype(str)))
    uncovered = sorted(set(priority_targets) - set(refs["target_gene"]))
    if uncovered:
        raise SystemExit("Priority targets missing from oncology reference curation: " + ", ".join(uncovered))

    target_compounds = target_compounds.copy()
    target_compounds["target_gene"] = target_compounds["target_gene"].astype(str).str.upper()
    if "preferred_name" not in target_compounds.columns:
        metadata = compounds[[c for c in ["compound_id", "preferred_name"] if c in compounds.columns]].drop_duplicates("compound_id")
        target_compounds = target_compounds.merge(metadata, on="compound_id", how="left")
    elif "preferred_name" in compounds.columns:
        fallback = compounds[["compound_id", "preferred_name"]].drop_duplicates("compound_id")
        target_compounds = target_compounds.merge(fallback, on="compound_id", how="left", suffixes=("", "_catalog"))
        current = target_compounds["preferred_name"].fillna("").astype(str).str.strip()
        target_compounds.loc[current.eq(""), "preferred_name"] = target_compounds.loc[current.eq(""), "preferred_name_catalog"]
        target_compounds = target_compounds.drop(columns=["preferred_name_catalog"])
    target_compounds["name_key"] = target_compounds.get("preferred_name", "").map(_norm)

    pcdb = pd.DataFrame()
    if PCDB_SUMMARY.exists():
        pcdb = pd.read_csv(PCDB_SUMMARY, sep="\t", dtype=str, keep_default_na=False)
        if not pcdb.empty:
            pcdb["target_gene"] = pcdb["target_gene"].str.upper()
            pcdb["name_key"] = pcdb["compound_name"].map(_norm)

    out_rows: list[dict[str, Any]] = []
    for _, row in refs.iterrows():
        keys = set(row["curation_keys"])
        target = _text(row["target_gene"])
        matches = target_compounds[(target_compounds["target_gene"] == target) & target_compounds["name_key"].isin(keys)].copy()
        matched_ids = sorted(set(matches["compound_id"].dropna().astype(str))) if not matches.empty else []
        matched_names = sorted(set(matches["preferred_name"].dropna().astype(str))) if not matches.empty else []

        pcdb_matches = pcdb[(pcdb["target_gene"] == target) & pcdb["name_key"].isin(keys)].copy() if not pcdb.empty else pd.DataFrame()
        pcdb_records = 0
        pcdb_oncology = 0
        pcdb_exact = 0
        if not pcdb_matches.empty:
            for col, dest in [
                ("pcdb_records_n", "records"),
                ("oncology_records_n", "oncology"),
                ("exact_priority_cancer_records_n", "exact"),
            ]:
                values = pd.to_numeric(pcdb_matches[col], errors="coerce").fillna(0) if col in pcdb_matches.columns else pd.Series([0])
                if dest == "records": pcdb_records = int(values.max())
                elif dest == "oncology": pcdb_oncology = int(values.max())
                else: pcdb_exact = int(values.max())

        record = {k: row[k] for k in refs.columns if k != "curation_keys"}
        record.update({
            "matched_mcl_compounds_n": len(matched_ids),
            "matched_mcl_compound_ids_json": json.dumps(matched_ids, ensure_ascii=False),
            "matched_mcl_names_json": json.dumps(matched_names, ensure_ascii=False),
            "pcdb_records_n": pcdb_records,
            "pcdb_oncology_records_n": pcdb_oncology,
            "pcdb_exact_priority_cancer_records_n": pcdb_exact,
            "similarity_eligible": bool(row["include_similarity"] and row["reference_level"] in {"A", "B"} and row["modality"] == "small_molecule"),
            "built_at": _now(),
        })
        out_rows.append(record)

    out = pd.DataFrame(out_rows).sort_values(["target_gene", "reference_level", "compound_name"]).reset_index(drop=True)

    summary_rows = []
    for target in priority_targets:
        g = out[out["target_gene"] == target]
        summary_rows.append({
            "target_gene": target,
            "priority_hypotheses_n": int((priority["target_gene"] == target).sum()),
            "reference_compounds_n": int(len(g)),
            "level_A_n": int((g["reference_level"] == "A").sum()),
            "level_B_n": int((g["reference_level"] == "B").sum()),
            "level_C_n": int((g["reference_level"] == "C").sum()),
            "level_D_n": int((g["reference_level"] == "D").sum()),
            "similarity_eligible_n": int(g["similarity_eligible"].sum()),
            "references_matched_to_current_mcl_n": int((g["matched_mcl_compounds_n"] > 0).sum()),
            "references_with_exact_pcdb_context_n": int((g["pcdb_exact_priority_cancer_records_n"] > 0).sum()),
        })
    target_summary = pd.DataFrame(summary_rows)

    qc_rows = []
    for _, row in out.iterrows():
        if row["similarity_eligible"] and int(row["matched_mcl_compounds_n"]) == 0:
            qc_rows.append({
                "reference_id": row["reference_id"],
                "target_gene": row["target_gene"],
                "compound_name": row["compound_name"],
                "qc_status": "not_in_current_mcl_catalog_structure_resolution_required",
            })
    qc = pd.DataFrame(qc_rows, columns=["reference_id", "target_gene", "compound_name", "qc_status"])

    PHARM.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUTPUT, index=False, compression="zstd")
    out.to_csv(OUTPUT_TSV, sep="\t", index=False)
    target_summary.to_csv(TARGET_SUMMARY, sep="\t", index=False)
    qc.to_csv(QC_OUTPUT, sep="\t", index=False)
    MANIFEST.write_text(json.dumps({
        "contract": "mcl-oncology-reference-compounds-v1",
        "built_at": _now(),
        "priority_status": PRIORITY_STATUS,
        "priority_targets_n": len(priority_targets),
        "reference_rows_n": int(len(out)),
        "similarity_eligible_levels": ["A", "B"],
        "similarity_modality": "small_molecule",
        "scientific_guardrails": [
            "Reference level is manually curated and does not prove efficacy in every MCL cancer context.",
            "Level A denotes approved oncology use with a direct/strong mechanistic relation to the target; level B denotes oncology clinical/investigational direct pharmacology.",
            "Level C is a pharmacology/preclinical tool and is excluded from the primary similarity universe.",
            "Level D is an unsuitable/uncertain target-specific oncology reference and is excluded.",
            "For FNTA/FNTB, farnesyltransferase inhibitors support the heterodimeric enzyme complex rather than isolated-subunit binding.",
            "Biologics and ADCs remain in the evidence registry but are excluded from small-molecule Tanimoto comparison.",
        ],
    }, ensure_ascii=False, indent=2), encoding="utf-8")

    print("MCL Oncology Reference Compounds v1")
    print(f"Priority targets: {len(priority_targets)}")
    print(f"Curated reference rows: {len(out)}")
    print(f"  A: {(out['reference_level'] == 'A').sum()}")
    print(f"  B: {(out['reference_level'] == 'B').sum()}")
    print(f"  C: {(out['reference_level'] == 'C').sum()}")
    print(f"  D: {(out['reference_level'] == 'D').sum()}")
    print(f"Small-molecule A+B similarity references: {int(out['similarity_eligible'].sum())}")
    print(f"Already matched to current MCL catalog: {int((out['matched_mcl_compounds_n'] > 0).sum())}/{len(out)}")
    if PCDB_SUMMARY.exists():
        print(f"References with exact priority-cancer PCDB evidence: {int((out['pcdb_exact_priority_cancer_records_n'] > 0).sum())}")
    print(f"Wrote {OUTPUT_TSV.relative_to(ROOT)}")
    print(f"Target summary: {TARGET_SUMMARY.relative_to(ROOT)}")
    print(f"QC: {QC_OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
