from __future__ import annotations

import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import quote

import httpx
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PHARM = ROOT / "data" / "runtime" / "pharmacology"
RUNTIME = PHARM
QC_DIR = ROOT / "outputs" / "qc"

COMPOUNDS = PHARM / "compound_catalog.parquet"
TARGET_COMPOUNDS = PHARM / "target_compound_catalog.parquet"
CANDIDATES = PHARM / "candidate_hypotheses_v2.parquet"
OUTPUT = RUNTIME / "compound_structure_registry.parquet"
OUTPUT_TSV = RUNTIME / "compound_structure_registry.tsv"
QC_OUTPUT = QC_DIR / "compound_structure_registry_qc.tsv"
MANIFEST = RUNTIME / "compound_structure_registry_manifest.json"

PRIORITY_STATUS = "priority_for_in_vitro"
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PROPERTIES = "CanonicalSMILES,IsomericSMILES,InChIKey"
USER_AGENT = "Master-Cancer-Landscape/0.5 structure-enrichment (research; PubChem PUG REST)"


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


def _existing_structure(row: pd.Series) -> str:
    return _text(row.get("canonical_smiles"))


def _candidate_queries(row: pd.Series) -> list[tuple[str, str]]:
    queries: list[tuple[str, str]] = []
    cid = _text(row.get("pubchem_cid"))
    if cid:
        # Some imported CIDs arrive as integer-looking floats.
        if cid.endswith(".0") and cid[:-2].isdigit():
            cid = cid[:-2]
        if cid.isdigit():
            queries.append(("pubchem_cid", cid))

    name = _text(row.get("preferred_name"))
    if name:
        queries.append(("preferred_name", name))

    chembl = _text(row.get("chembl_id"))
    if chembl:
        queries.append(("chembl_id_as_pubchem_name", chembl))

    broad = _text(row.get("broad_id"))
    if broad:
        queries.append(("broad_id_as_pubchem_name", broad))

    # Preserve order while dropping identical values.
    seen: set[tuple[str, str]] = set()
    result: list[tuple[str, str]] = []
    for item in queries:
        if item not in seen:
            result.append(item)
            seen.add(item)
    return result


def _property_url(query_type: str, value: str) -> str:
    if query_type == "pubchem_cid":
        return f"{PUBCHEM_BASE}/compound/cid/{quote(value, safe='')}/property/{PROPERTIES}/JSON"
    return f"{PUBCHEM_BASE}/compound/name/{quote(value, safe='')}/property/{PROPERTIES}/JSON"


def _request_properties(
    client: httpx.Client,
    query_type: str,
    value: str,
    retries: int,
    pause_s: float,
) -> tuple[str, dict[str, Any] | None, str]:
    url = _property_url(query_type, value)
    last_error = ""
    for attempt in range(1, retries + 1):
        try:
            response = client.get(url)
        except httpx.HTTPError as exc:
            last_error = f"http_error:{type(exc).__name__}"
            if attempt < retries:
                time.sleep(max(pause_s, 0.25) * attempt)
                continue
            return "request_failed", None, last_error

        if response.status_code == 404:
            return "not_found", None, "pubchem_404"
        if response.status_code == 429 or response.status_code >= 500:
            last_error = f"pubchem_http_{response.status_code}"
            if attempt < retries:
                time.sleep(max(pause_s, 0.5) * attempt)
                continue
            return "request_failed", None, last_error
        if response.status_code != 200:
            return "not_found", None, f"pubchem_http_{response.status_code}"

        try:
            payload = response.json()
        except ValueError:
            return "request_failed", None, "invalid_json"

        properties = payload.get("PropertyTable", {}).get("Properties", [])
        if not isinstance(properties, list) or not properties:
            return "not_found", None, "empty_property_table"

        by_cid: dict[str, dict[str, Any]] = {}
        for record in properties:
            cid = _text(record.get("CID"))
            if cid:
                by_cid[cid] = record
        if len(by_cid) != 1:
            return "ambiguous", None, f"unique_cids={len(by_cid)}"

        record = next(iter(by_cid.values()))
        canonical = _text(record.get("ConnectivitySMILES")) or _text(record.get("CanonicalSMILES"))
        isomeric = _text(record.get("SMILES")) or _text(record.get("IsomericSMILES"))
        smiles = isomeric or canonical
        inchikey = _text(record.get("InChIKey"))
        if not smiles:
            return "not_found", None, "pubchem_record_without_smiles"

        return (
            "resolved",
            {
                "pubchem_cid_resolved": _text(record.get("CID")),
                "canonical_smiles": canonical or smiles,
                "isomeric_smiles": isomeric,
                "inchikey": inchikey,
            },
            "",
        )
    return "request_failed", None, last_error or "unknown_request_failure"


def _load_universe() -> pd.DataFrame:
    required = [COMPOUNDS, TARGET_COMPOUNDS, CANDIDATES]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    if missing:
        raise SystemExit("Missing structure-enrichment inputs:\n" + "\n".join(missing))

    compounds = pd.read_parquet(COMPOUNDS)
    target_compounds = pd.read_parquet(TARGET_COMPOUNDS)
    candidates = pd.read_parquet(CANDIDATES)

    if "priority_status" not in candidates.columns:
        raise SystemExit("candidate_hypotheses_v2.parquet has no priority_status column")
    if "target_gene" not in candidates.columns or "target_gene" not in target_compounds.columns:
        raise SystemExit("Expected target_gene in Candidate v2 and target_compound_catalog")

    priority_targets = set(
        candidates.loc[
            candidates["priority_status"].astype(str).eq(PRIORITY_STATUS), "target_gene"
        ].dropna().astype(str).str.upper()
    )
    target_compounds = target_compounds.copy()
    target_compounds["target_gene"] = target_compounds["target_gene"].astype(str).str.upper()
    linked_ids = set(
        target_compounds.loc[
            target_compounds["target_gene"].isin(priority_targets), "compound_id"
        ].dropna().astype(str)
    )

    universe = compounds[compounds["compound_id"].astype(str).isin(linked_ids)].copy()
    if universe.empty:
        raise SystemExit("No known compounds are linked to Candidate v2 priority targets.")
    universe["compound_id"] = universe["compound_id"].astype(str)
    return universe.drop_duplicates("compound_id", keep="first").reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve missing structures for known compounds linked to Candidate v2 priority targets "
            "using conservative PubChem PUG REST lookups."
        )
    )
    parser.add_argument("--force", action="store_true", help="Retry compounds already present in the registry.")
    parser.add_argument("--pause", type=float, default=0.18, help="Pause between PubChem requests in seconds.")
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=20.0)
    args = parser.parse_args()

    universe = _load_universe()
    existing_registry = pd.DataFrame()
    if OUTPUT.exists() and not args.force:
        existing_registry = pd.read_parquet(OUTPUT)
        if "compound_id" in existing_registry.columns:
            existing_registry["compound_id"] = existing_registry["compound_id"].astype(str)

    existing_map: dict[str, dict[str, Any]] = {}
    if not existing_registry.empty:
        for record in existing_registry.to_dict("records"):
            compound_id = _text(record.get("compound_id"))
            if compound_id:
                existing_map[compound_id] = record

    rows: list[dict[str, Any]] = []
    requests_n = 0
    with httpx.Client(
        timeout=args.timeout,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
        follow_redirects=True,
    ) as client:
        for idx, row in universe.iterrows():
            compound_id = _text(row.get("compound_id"))
            preferred_name = _text(row.get("preferred_name"))
            source_smiles = _existing_structure(row)

            if source_smiles:
                rows.append(
                    {
                        "compound_id": compound_id,
                        "preferred_name": preferred_name,
                        "canonical_smiles": source_smiles,
                        "isomeric_smiles": "",
                        "inchikey": _text(row.get("inchikey")),
                        "pubchem_cid_resolved": _text(row.get("pubchem_cid")),
                        "resolution_status": "resolved",
                        "resolution_source": "mcl_source_catalog",
                        "resolution_method": "existing_canonical_smiles",
                        "resolution_query": "",
                        "resolution_note": "Structure already present in the normalized pharmacology source.",
                        "retrieved_at": _now(),
                    }
                )
                continue

            cached = existing_map.get(compound_id)
            if cached and _text(cached.get("canonical_smiles")) and _text(cached.get("resolution_status")) == "resolved":
                rows.append(cached)
                continue

            resolved: dict[str, Any] | None = None
            final_status = "not_found"
            notes: list[str] = []
            used_type = ""
            used_value = ""

            queries = _candidate_queries(row)
            if not queries:
                notes.append("no_pubchem_query_identifier")

            for query_type, value in queries:
                requests_n += 1
                status, record, note = _request_properties(
                    client, query_type, value, max(1, args.retries), max(0.0, args.pause)
                )
                notes.append(f"{query_type}:{status}:{note}".rstrip(":"))
                if status == "resolved" and record is not None:
                    resolved = record
                    final_status = "resolved"
                    used_type = query_type
                    used_value = value
                    break
                if status == "ambiguous":
                    # Do not fall through to weaker identifiers after an ambiguous preferred identifier.
                    final_status = "ambiguous"
                    used_type = query_type
                    used_value = value
                    break
                if status == "request_failed":
                    final_status = "request_failed"
                time.sleep(max(0.0, args.pause))

            if resolved:
                rows.append(
                    {
                        "compound_id": compound_id,
                        "preferred_name": preferred_name,
                        "canonical_smiles": resolved["canonical_smiles"],
                        "isomeric_smiles": resolved["isomeric_smiles"],
                        "inchikey": resolved["inchikey"],
                        "pubchem_cid_resolved": resolved["pubchem_cid_resolved"],
                        "resolution_status": "resolved",
                        "resolution_source": "PubChem PUG REST",
                        "resolution_method": used_type,
                        "resolution_query": used_value,
                        "resolution_note": " | ".join(notes),
                        "retrieved_at": _now(),
                    }
                )
            else:
                rows.append(
                    {
                        "compound_id": compound_id,
                        "preferred_name": preferred_name,
                        "canonical_smiles": "",
                        "isomeric_smiles": "",
                        "inchikey": "",
                        "pubchem_cid_resolved": "",
                        "resolution_status": final_status,
                        "resolution_source": "PubChem PUG REST",
                        "resolution_method": used_type,
                        "resolution_query": used_value,
                        "resolution_note": " | ".join(notes),
                        "retrieved_at": _now(),
                    }
                )

            if (idx + 1) % 25 == 0:
                print(f"Processed {idx + 1}/{len(universe)} known compounds ...")

    registry = pd.DataFrame(rows).drop_duplicates("compound_id", keep="last")
    registry = registry.sort_values(["resolution_status", "preferred_name", "compound_id"]).reset_index(drop=True)

    RUNTIME.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    registry.to_parquet(OUTPUT, index=False, compression="zstd")
    registry.to_csv(OUTPUT_TSV, sep="\t", index=False)

    qc = registry[registry["resolution_status"].astype(str).ne("resolved")].copy()
    qc.to_csv(QC_OUTPUT, sep="\t", index=False)

    resolved_n = int(registry["resolution_status"].astype(str).eq("resolved").sum())
    source_catalog_n = int(registry["resolution_source"].astype(str).eq("mcl_source_catalog").sum())
    pubchem_n = int(registry["resolution_source"].astype(str).eq("PubChem PUG REST").mul(
        registry["resolution_status"].astype(str).eq("resolved")
    ).sum())
    unresolved_n = int(len(registry) - resolved_n)
    coverage = resolved_n / len(registry) if len(registry) else 0.0

    manifest = {
        "contract": "mcl-priority-compound-structure-registry-v1",
        "built_at": _now(),
        "priority_target_linked_compounds_n": int(len(registry)),
        "resolved_structures_n": resolved_n,
        "resolved_from_source_catalog_n": source_catalog_n,
        "resolved_from_pubchem_n": pubchem_n,
        "unresolved_or_ambiguous_n": unresolved_n,
        "structure_coverage_fraction": coverage,
        "pubchem_requests_n": requests_n,
        "resolution_policy": (
            "Use existing normalized MCL structure first. Otherwise query PubChem in order: explicit PubChem CID, "
            "preferred compound name, ChEMBL ID as PubChem name, Broad ID as PubChem name. Accept only a query "
            "that resolves to exactly one unique PubChem CID. Ambiguous results are not auto-selected."
        ),
        "outputs": [
            str(OUTPUT.relative_to(ROOT)),
            str(OUTPUT_TSV.relative_to(ROOT)),
            str(QC_OUTPUT.relative_to(ROOT)),
        ],
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print("MCL Priority Compound Structure Registry v1")
    print(f"Known compounds linked to priority targets: {len(registry)}")
    print(f"Resolved structures: {resolved_n} ({coverage:.1%})")
    print(f"  from MCL source catalog: {source_catalog_n}")
    print(f"  from PubChem: {pubchem_n}")
    print(f"Unresolved / ambiguous: {unresolved_n}")
    print(f"PubChem requests: {requests_n}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")
    print(f"QC: {QC_OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
