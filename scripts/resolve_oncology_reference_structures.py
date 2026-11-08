from __future__ import annotations

import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import httpx
import pandas as pd

from enrich_priority_compound_structures_pubchem import _request_properties

ROOT = Path(__file__).resolve().parents[1]
PHARM = ROOT / "data" / "runtime" / "pharmacology"
QC_DIR = ROOT / "outputs" / "qc"

REFERENCES = PHARM / "oncology_reference_compounds.parquet"
MCL_STRUCTURES = PHARM / "compound_structure_registry.parquet"
OUTPUT = PHARM / "oncology_reference_structures.parquet"
OUTPUT_TSV = PHARM / "oncology_reference_structures.tsv"
QC_OUTPUT = QC_DIR / "oncology_reference_structures_qc.tsv"
MANIFEST = PHARM / "oncology_reference_structures_manifest.json"

USER_AGENT = "Master-Cancer-Landscape/0.6 oncology-reference-structure-resolution"


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


def _alias_queries(row: pd.Series) -> list[str]:
    values = [_text(row.get("compound_name"))]
    values += [x.strip() for x in _text(row.get("aliases")).split("|") if x.strip()]
    result: list[str] = []
    for value in values:
        if value and value not in result:
            result.append(value)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Resolve structures for strict MCL oncology reference compounds.")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--pause", type=float, default=0.18)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=20.0)
    args = parser.parse_args()

    if not REFERENCES.exists():
        raise SystemExit("Missing oncology_reference_compounds.parquet; run build_oncology_reference_compounds.py first.")

    refs = pd.read_parquet(REFERENCES)
    refs = refs[refs["similarity_eligible"].astype(bool)].copy()
    if refs.empty:
        raise SystemExit("No A/B small-molecule oncology references are eligible for similarity.")

    source_map: dict[str, dict[str, Any]] = {}
    if MCL_STRUCTURES.exists():
        registry = pd.read_parquet(MCL_STRUCTURES)
        for rec in registry.to_dict("records"):
            cid = _text(rec.get("compound_id"))
            if cid and _text(rec.get("canonical_smiles")) and _text(rec.get("resolution_status")) == "resolved":
                source_map[cid] = rec

    existing: dict[str, dict[str, Any]] = {}
    if OUTPUT.exists() and not args.force:
        old = pd.read_parquet(OUTPUT)
        for rec in old.to_dict("records"):
            rid = _text(rec.get("reference_id"))
            if rid and _text(rec.get("canonical_smiles")) and _text(rec.get("resolution_status")) == "resolved":
                existing[rid] = rec

    rows: list[dict[str, Any]] = []
    name_cache: dict[str, dict[str, Any]] = {}
    requests_n = 0

    with httpx.Client(timeout=args.timeout, headers={"User-Agent": USER_AGENT, "Accept": "application/json"}, follow_redirects=True) as client:
        for _, row in refs.iterrows():
            rid = _text(row.get("reference_id"))
            if rid in existing:
                cached = dict(existing[rid])
                rows.append(cached)
                continue

            base = row.to_dict()
            matched_ids: list[str] = []
            try:
                matched_ids = json.loads(_text(row.get("matched_mcl_compound_ids_json")) or "[]")
            except json.JSONDecodeError:
                matched_ids = []

            resolved = None
            for compound_id in matched_ids:
                record = source_map.get(str(compound_id))
                if record:
                    resolved = {
                        "canonical_smiles": _text(record.get("canonical_smiles")),
                        "isomeric_smiles": _text(record.get("isomeric_smiles")),
                        "inchikey": _text(record.get("inchikey")),
                        "pubchem_cid_resolved": _text(record.get("pubchem_cid_resolved")),
                        "resolution_status": "resolved",
                        "resolution_source": "mcl_compound_structure_registry",
                        "resolution_query": str(compound_id),
                        "resolution_note": "Reused structure already resolved for a matched MCL pharmacology compound.",
                    }
                    break

            notes: list[str] = []
            if resolved is None:
                for query in _alias_queries(row):
                    cache_key = query.upper()
                    if cache_key in name_cache:
                        resolved = dict(name_cache[cache_key])
                        resolved["resolution_source"] = "pubchem_name_cache"
                        resolved["resolution_query"] = query
                        break
                    status, properties, note = _request_properties(client, "preferred_name", query, args.retries, args.pause)
                    requests_n += 1
                    if status == "resolved" and properties:
                        resolved = {
                            **properties,
                            "resolution_status": "resolved",
                            "resolution_source": "pubchem_name",
                            "resolution_query": query,
                            "resolution_note": "Resolved conservatively to one PubChem CID by reference name/alias.",
                        }
                        name_cache[cache_key] = dict(resolved)
                        break
                    notes.append(f"{query}:{status}:{note}")
                    time.sleep(max(args.pause, 0.0))

            if resolved is None:
                resolved = {
                    "canonical_smiles": "",
                    "isomeric_smiles": "",
                    "inchikey": "",
                    "pubchem_cid_resolved": "",
                    "resolution_status": "unresolved",
                    "resolution_source": "",
                    "resolution_query": "",
                    "resolution_note": "; ".join(notes),
                }

            base.update(resolved)
            base["retrieved_at"] = _now()
            rows.append(base)

    out = pd.DataFrame(rows).sort_values(["target_gene", "reference_level", "compound_name"]).reset_index(drop=True)
    qc = out[out["resolution_status"] != "resolved"][["reference_id", "target_gene", "compound_name", "resolution_status", "resolution_note"]].copy()

    PHARM.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUTPUT, index=False, compression="zstd")
    out.to_csv(OUTPUT_TSV, sep="\t", index=False)
    qc.to_csv(QC_OUTPUT, sep="\t", index=False)
    MANIFEST.write_text(json.dumps({
        "contract": "mcl-oncology-reference-structures-v1",
        "built_at": _now(),
        "eligible_reference_rows_n": int(len(out)),
        "resolved_n": int((out["resolution_status"] == "resolved").sum()),
        "unresolved_n": int((out["resolution_status"] != "resolved").sum()),
        "pubchem_requests_n": requests_n,
        "guardrail": "Structure resolution does not strengthen target/mechanism evidence; it only enables chemical comparison.",
    }, ensure_ascii=False, indent=2), encoding="utf-8")

    print("MCL Oncology Reference Structures v1")
    print(f"Eligible A/B small-molecule references: {len(out)}")
    print(f"Resolved structures: {(out['resolution_status'] == 'resolved').sum()}/{len(out)}")
    print(f"PubChem requests: {requests_n}")
    print(f"Wrote {OUTPUT_TSV.relative_to(ROOT)}")
    print(f"QC: {QC_OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
