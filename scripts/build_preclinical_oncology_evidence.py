from __future__ import annotations

import argparse
import ast
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PHARM = ROOT / "data" / "runtime" / "pharmacology"
PCDB_BASE = ROOT / "data" / "external" / "preclinical_db"
OUT = ROOT / "data" / "runtime" / "preclinical"
QC = ROOT / "outputs" / "qc"

CANDIDATES = PHARM / "candidate_hypotheses_v2.parquet"
TARGET_COMPOUNDS = PHARM / "target_compound_catalog.parquet"
COMPOUND_CATALOG = PHARM / "compound_catalog.parquet"
ENRICHED_COMPOUND_CATALOG = PHARM / "compound_catalog_structures_enriched.parquet"

PRIORITY_STATUS = "priority_for_in_vitro"
ONCOLOGY_PATTERN = re.compile(
    r"\b(cancer|carcinoma|adenocarcinoma|leukemia|leukaemia|lymphoma|melanoma|sarcoma|glioma|glioblastoma|"
    r"neoplasm|tumou?r|myeloma|mesothelioma|neuroblastoma|retinoblastoma|medulloblastoma)\b",
    re.I,
)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _norm(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    text = str(value).strip().lower()
    text = text.replace("β", "beta").replace("α", "alpha")
    return re.sub(r"[^a-z0-9]+", "", text)


def _text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _split_entities(value: Any) -> list[str]:
    text = _text(value)
    if not text:
        return []
    if text[0] in "[{(" and text[-1] in "]})":
        for loader in (json.loads, ast.literal_eval):
            try:
                parsed = loader(text)
                if isinstance(parsed, dict):
                    parsed = list(parsed.values())
                if isinstance(parsed, (list, tuple, set)):
                    return [str(x).strip() for x in parsed if str(x).strip()]
            except Exception:
                continue
    for sep in ("|", ";"):
        if sep in text:
            return [part.strip() for part in text.split(sep) if part.strip()]
    return [text]


def _find_col(columns: Iterable[str], aliases: tuple[str, ...], contains: tuple[str, ...] = ()) -> str | None:
    lookup = {str(c).strip().lower(): str(c) for c in columns}
    for alias in aliases:
        if alias.lower() in lookup:
            return lookup[alias.lower()]
    for column in columns:
        low = str(column).lower()
        if contains and all(token in low for token in contains):
            return str(column)
    return None


def _find_pcdb_complete(release: str) -> Path:
    release_dir = PCDB_BASE / release / "extracted"
    matches = list(release_dir.rglob("pcdb_complete.tsv"))
    if len(matches) != 1:
        raise SystemExit(
            f"Expected exactly one pcdb_complete.tsv under {release_dir}; found {len(matches)}. "
            "Run scripts/fetch_preclinical_db.py first."
        )
    return matches[0]


def _load_mcl_universe() -> tuple[pd.DataFrame, pd.DataFrame]:
    for path in (CANDIDATES, TARGET_COMPOUNDS):
        if not path.exists():
            raise SystemExit(f"Missing required MCL input: {path.relative_to(ROOT)}")

    candidates = pd.read_parquet(CANDIDATES)
    if "priority_status" not in candidates.columns or "target_gene" not in candidates.columns:
        raise SystemExit("candidate_hypotheses_v2.parquet lacks priority_status/target_gene")
    priority = candidates[candidates["priority_status"].astype(str).eq(PRIORITY_STATUS)].copy()
    priority["target_gene"] = priority["target_gene"].astype(str).str.upper()

    tc = pd.read_parquet(TARGET_COMPOUNDS).copy()
    tc["target_gene"] = tc["target_gene"].astype(str).str.upper()
    tc = tc.merge(priority[["target_gene"]].drop_duplicates(), on="target_gene", how="inner")

    catalog_path = ENRICHED_COMPOUND_CATALOG if ENRICHED_COMPOUND_CATALOG.exists() else COMPOUND_CATALOG
    if catalog_path.exists():
        catalog = pd.read_parquet(catalog_path)
        keep = [c for c in ("compound_id", "preferred_name", "pubchem_cid", "chembl_id", "canonical_smiles") if c in catalog.columns]
        if "compound_id" in keep:
            tc = tc.merge(catalog[keep].drop_duplicates("compound_id"), on="compound_id", how="left", suffixes=("", "_catalog"))
            for col in ("preferred_name", "pubchem_cid", "chembl_id", "canonical_smiles"):
                alt = f"{col}_catalog"
                if alt in tc.columns:
                    if col not in tc.columns:
                        tc[col] = tc[alt]
                    else:
                        current = tc[col].fillna("").astype(str).str.strip()
                        tc[col] = tc[col].astype("object")
                        tc.loc[current.eq(""), col] = tc.loc[current.eq(""), alt]
                    tc = tc.drop(columns=alt)

    tc = tc.drop_duplicates(["target_gene", "compound_id"]).copy()
    return priority, tc


def _priority_cancers(priority: pd.DataFrame) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    name_col = "mcl_cancer_name" if "mcl_cancer_name" in priority.columns else None
    for target, group in priority.groupby("target_gene"):
        if not name_col:
            result[str(target)] = []
            continue
        vals = []
        for value in group[name_col]:
            text = _text(value)
            if text and text not in vals:
                vals.append(text)
        result[str(target)] = vals
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Link PCDB in-vivo study records to compounds around MCL priority targets.")
    parser.add_argument("--release", default="v1.0")
    args = parser.parse_args()

    pcdb_path = _find_pcdb_complete(args.release)
    priority, target_compounds = _load_mcl_universe()
    target_cancers = _priority_cancers(priority)

    pcdb = pd.read_csv(pcdb_path, sep="\t", dtype=str, keep_default_na=False, low_memory=False)
    cols = list(pcdb.columns)
    drug_col = _find_col(cols, ("drug", "drugs", "drug_name", "drug_names"), ("drug",))
    disease_col = _find_col(cols, ("disease", "diseases", "disease_name", "disease_names"), ("disease",))
    pmcid_col = _find_col(cols, ("pmcid", "pmc_id"), ("pmc",))
    title_col = _find_col(cols, ("title", "publication_title", "article_title"), ("title",))
    animal_col = _find_col(cols, ("animal", "animals", "animal_name", "animal_names"), ("animal",))
    confidence_col = _find_col(cols, ("confidence_score", "confidence", "classification_confidence"), ("confidence",))
    subject_col = _find_col(cols, ("total_subject_size", "subject_size", "n_animals"), ("subject", "size"))
    date_col = _find_col(cols, ("date", "publication_date", "year"), ())
    link_col = _find_col(cols, ("link", "url", "pmc_link"), ())

    if not drug_col or not disease_col:
        QC.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"column": cols}).to_csv(QC / "preclinical_db_detected_columns.tsv", sep="\t", index=False)
        raise SystemExit(
            "Could not identify drug and disease columns in pcdb_complete.tsv. "
            "Detected columns were written to outputs/qc/preclinical_db_detected_columns.tsv."
        )

    # Conservative compound matching: normalized exact name only. We never use fuzzy matching automatically.
    compound_rows: list[dict[str, Any]] = []
    name_to_mcl: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for _, row in target_compounds.iterrows():
        name = _text(row.get("preferred_name"))
        if not name:
            continue
        item = {
            "target_gene": _text(row.get("target_gene")).upper(),
            "compound_id": _text(row.get("compound_id")),
            "preferred_name": name,
            "pubchem_cid": _text(row.get("pubchem_cid")),
            "chembl_id": _text(row.get("chembl_id")),
        }
        name_to_mcl[_norm(name)].append(item)
        compound_rows.append(item)

    linked: list[dict[str, Any]] = []
    unmatched_pcdb_drugs: set[str] = set()
    for idx, row in pcdb.iterrows():
        drugs = _split_entities(row.get(drug_col))
        diseases = _split_entities(row.get(disease_col))
        animals = _split_entities(row.get(animal_col)) if animal_col else []
        for drug_name in drugs:
            matches = name_to_mcl.get(_norm(drug_name), [])
            if not matches:
                if drug_name:
                    unmatched_pcdb_drugs.add(drug_name)
                continue
            for match in matches:
                candidate_names = target_cancers.get(match["target_gene"], [])
                candidate_norms = {_norm(x): x for x in candidate_names}
                if not diseases:
                    diseases_iter = [""]
                else:
                    diseases_iter = diseases
                for disease in diseases_iter:
                    exact_candidate = candidate_norms.get(_norm(disease), "") if disease else ""
                    looks_oncologic = bool(ONCOLOGY_PATTERN.search(disease)) if disease else False
                    linked.append(
                        {
                            "target_gene": match["target_gene"],
                            "compound_id": match["compound_id"],
                            "mcl_compound_name": match["preferred_name"],
                            "pcdb_drug_name": drug_name,
                            "compound_match_method": "normalized_exact_name",
                            "mcl_pubchem_cid": match["pubchem_cid"],
                            "mcl_chembl_id": match["chembl_id"],
                            "pcdb_disease_name": disease,
                            "disease_looks_oncologic": looks_oncologic,
                            "exact_priority_cancer_match": bool(exact_candidate),
                            "matched_priority_cancer_name": exact_candidate,
                            "priority_cancers_json": json.dumps(candidate_names, ensure_ascii=False),
                            "pcdb_animals_json": json.dumps(animals, ensure_ascii=False),
                            "pmcid": _text(row.get(pmcid_col)) if pmcid_col else "",
                            "title": _text(row.get(title_col)) if title_col else "",
                            "confidence_score": _text(row.get(confidence_col)) if confidence_col else "",
                            "total_subject_size": _text(row.get(subject_col)) if subject_col else "",
                            "publication_date": _text(row.get(date_col)) if date_col else "",
                            "link": _text(row.get(link_col)) if link_col else "",
                            "pcdb_source_row": int(idx),
                        }
                    )

    OUT.mkdir(parents=True, exist_ok=True)
    QC.mkdir(parents=True, exist_ok=True)
    records = pd.DataFrame(linked)
    if records.empty:
        records = pd.DataFrame(
            columns=[
                "target_gene", "compound_id", "mcl_compound_name", "pcdb_drug_name", "compound_match_method",
                "pcdb_disease_name", "disease_looks_oncologic", "exact_priority_cancer_match", "pmcid", "title",
            ]
        )
    records.to_parquet(OUT / "preclinical_compound_records.parquet", index=False, compression="zstd")
    records.to_csv(OUT / "preclinical_compound_records.tsv", sep="\t", index=False)

    summary_rows: list[dict[str, Any]] = []
    if not records.empty:
        for (target, compound_id, name), group in records.groupby(["target_gene", "compound_id", "mcl_compound_name"], dropna=False):
            onc = group[group["disease_looks_oncologic"].astype(bool)]
            exact = group[group["exact_priority_cancer_match"].astype(bool)]
            summary_rows.append(
                {
                    "target_gene": target,
                    "compound_id": compound_id,
                    "compound_name": name,
                    "pcdb_records_n": int(len(group)),
                    "pcdb_publications_n": int(group["pmcid"].replace("", pd.NA).nunique()),
                    "oncology_records_n": int(len(onc)),
                    "oncology_publications_n": int(onc["pmcid"].replace("", pd.NA).nunique()),
                    "exact_priority_cancer_records_n": int(len(exact)),
                    "exact_priority_cancers_json": json.dumps(
                        sorted({x for x in exact["matched_priority_cancer_name"].astype(str) if x}), ensure_ascii=False
                    ),
                    "pcdb_diseases_json": json.dumps(
                        sorted({x for x in group["pcdb_disease_name"].astype(str) if x}), ensure_ascii=False
                    ),
                }
            )
    summary = pd.DataFrame(summary_rows)
    if summary.empty:
        summary = pd.DataFrame(columns=["target_gene", "compound_id", "compound_name", "pcdb_records_n"])
    summary.to_parquet(OUT / "preclinical_compound_summary.parquet", index=False, compression="zstd")
    summary.to_csv(OUT / "preclinical_compound_summary.tsv", sep="\t", index=False)

    target_rows: list[dict[str, Any]] = []
    for target in sorted(set(priority["target_gene"].astype(str))):
        group = summary[summary["target_gene"].astype(str).eq(target)] if not summary.empty else summary
        target_rows.append(
            {
                "target_gene": target,
                "priority_hypotheses_n": int((priority["target_gene"].astype(str) == target).sum()),
                "target_linked_compounds_n": int(target_compounds[target_compounds["target_gene"].astype(str).eq(target)]["compound_id"].nunique()),
                "compounds_with_pcdb_evidence_n": int(group["compound_id"].nunique()) if not group.empty else 0,
                "compounds_with_oncology_pcdb_evidence_n": int(group.loc[group.get("oncology_records_n", 0).astype(int).gt(0), "compound_id"].nunique()) if not group.empty and "oncology_records_n" in group else 0,
                "compounds_with_exact_priority_cancer_pcdb_evidence_n": int(group.loc[group.get("exact_priority_cancer_records_n", 0).astype(int).gt(0), "compound_id"].nunique()) if not group.empty and "exact_priority_cancer_records_n" in group else 0,
                "priority_cancers_json": json.dumps(target_cancers.get(target, []), ensure_ascii=False),
            }
        )
    target_summary = pd.DataFrame(target_rows)
    target_summary.to_parquet(OUT / "preclinical_target_summary.parquet", index=False, compression="zstd")
    target_summary.to_csv(OUT / "preclinical_target_summary.tsv", sep="\t", index=False)

    pd.DataFrame({"pcdb_drug_name_unmatched": sorted(unmatched_pcdb_drugs)}).to_csv(
        QC / "preclinical_unmatched_pcdb_drug_names.tsv", sep="\t", index=False
    )

    manifest = {
        "contract": "mcl-preclinical-evidence-v1",
        "built_at": _now(),
        "pcdb_release": args.release,
        "pcdb_complete": str(pcdb_path.relative_to(ROOT)),
        "priority_targets_n": int(priority["target_gene"].nunique()),
        "priority_hypotheses_n": int(len(priority)),
        "target_linked_compounds_n": int(target_compounds["compound_id"].nunique()),
        "linked_pcdb_rows_n": int(len(records)),
        "compound_match_policy": "normalized exact preferred_name only; no fuzzy matching",
        "same_cancer_policy": "only normalized exact PCDB disease name == MCL priority cancer name is called an exact priority-cancer match",
        "oncology_flag_policy": "descriptive keyword flag only; it does not prove antitumor efficacy",
        "scientific_guardrail": (
            "A PCDB record shows that a compound was reported in an in-vivo drug study for a disease/animal context. "
            "PCDB v1.0 does not provide outcome values, effect sizes, direct target engagement, or proof that the observed effect is mediated by the MCL target."
        ),
        "outputs": [
            "data/runtime/preclinical/preclinical_compound_records.parquet",
            "data/runtime/preclinical/preclinical_compound_summary.tsv",
            "data/runtime/preclinical/preclinical_target_summary.tsv",
        ],
    }
    (OUT / "preclinical_evidence_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print("MCL Preclinical Evidence v1")
    print(f"PCDB release: {args.release}")
    print(f"Priority targets: {priority['target_gene'].nunique()}")
    print(f"Priority hypotheses: {len(priority)}")
    print(f"MCL target-linked compounds: {target_compounds['compound_id'].nunique()}")
    print(f"Linked PCDB evidence rows: {len(records)}")
    if not summary.empty and "oncology_records_n" in summary.columns:
        print(f"Compounds with oncology-like PCDB records: {(summary['oncology_records_n'].astype(int) > 0).sum()}")
        print(f"Compounds with exact priority-cancer records: {(summary['exact_priority_cancer_records_n'].astype(int) > 0).sum()}")
    print(f"Wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
