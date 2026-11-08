from __future__ import annotations

import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import httpx
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

from enrich_priority_compound_structures_pubchem import _request_properties


ROOT = Path(__file__).resolve().parents[1]
PHARM = ROOT / "data" / "runtime" / "pharmacology"
QC_DIR = ROOT / "outputs" / "qc"

COMPOUNDS = PHARM / "compounds.parquet"
RESPONSES = PHARM / "responses.parquet"
GENERIC_REGISTRY = PHARM / "compound_structure_registry.parquet"
OUTPUT = PHARM / "prism_structure_registry.parquet"
OUTPUT_TSV = PHARM / "prism_structure_registry.tsv"
QC_OUTPUT = QC_DIR / "prism_structure_registry_qc.tsv"
MANIFEST = PHARM / "prism_structure_registry_manifest.json"

USER_AGENT = "Master-Cancer-Landscape/0.5 PRISM-structure-resolution (research; PubChem PUG REST)"


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


def _standardize(smiles: Any) -> tuple[str | None, str | None, str | None]:
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


def _prism_compound_ids() -> set[str]:
    dataset = ds.dataset(RESPONSES, format="parquet")
    expression = (ds.field("source") == "PRISM") & (ds.field("endpoint") == "LFC")
    table = dataset.to_table(columns=["compound_id"], filter=expression)
    frame = table.to_pandas()
    return set(frame["compound_id"].dropna().astype(str))


def _generic_map() -> dict[str, dict[str, Any]]:
    if not GENERIC_REGISTRY.exists():
        return {}
    frame = pd.read_parquet(GENERIC_REGISTRY)
    if "compound_id" not in frame.columns:
        return {}
    out: dict[str, dict[str, Any]] = {}
    for row in frame.to_dict("records"):
        compound_id = _text(row.get("compound_id"))
        smiles = _text(row.get("canonical_smiles"))
        status = _text(row.get("resolution_status"))
        if compound_id and smiles and (not status or status == "resolved"):
            out[compound_id] = row
    return out


def _candidate_queries(row: dict[str, Any]) -> list[tuple[str, str]]:
    queries: list[tuple[str, str]] = []
    cid = _text(row.get("pubchem_cid"))
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

    seen: set[tuple[str, str]] = set()
    result: list[tuple[str, str]] = []
    for item in queries:
        if item not in seen:
            seen.add(item)
            result.append(item)
    return result


def _write_checkpoint(rows: list[dict[str, Any]]) -> None:
    frame = pd.DataFrame(rows)
    if frame.empty:
        return
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUTPUT.with_suffix(".checkpoint.parquet")
    frame.to_parquet(tmp, index=False, compression="zstd")
    tmp.replace(OUTPUT)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve structures for PRISM 24Q2 compounds using existing MCL structures first, "
            "then conservative PubChem name/ID lookups. Results are cached and resumable."
        )
    )
    parser.add_argument("--force", action="store_true", help="Retry all PRISM compounds, including cached rows.")
    parser.add_argument(
        "--retry-unresolved",
        action="store_true",
        help="Retry cached unresolved rows while retaining cached resolved structures.",
    )
    parser.add_argument("--pause", type=float, default=0.15)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=20.0)
    parser.add_argument("--checkpoint-every", type=int, default=100)
    args = parser.parse_args()

    for path in (COMPOUNDS, RESPONSES):
        if not path.exists():
            raise SystemExit(f"Missing {path.relative_to(ROOT)}")

    prism_ids = _prism_compound_ids()
    compounds = pd.read_parquet(COMPOUNDS)
    compounds["compound_id"] = compounds["compound_id"].astype(str)
    universe = compounds[compounds["compound_id"].isin(prism_ids)].copy()
    universe = universe.drop_duplicates("compound_id", keep="first").reset_index(drop=True)
    if universe.empty:
        raise SystemExit("No PRISM compounds were found in the pharmacology compound catalog.")

    generic_map = _generic_map()
    cached_map: dict[str, dict[str, Any]] = {}
    if OUTPUT.exists() and not args.force:
        cached = pd.read_parquet(OUTPUT)
        if "compound_id" in cached.columns:
            for row in cached.to_dict("records"):
                compound_id = _text(row.get("compound_id"))
                if compound_id:
                    cached_map[compound_id] = row

    rows: list[dict[str, Any]] = []
    requests_n = 0
    reused_n = 0
    source_n = 0
    generic_n = 0
    pubchem_n = 0
    unresolved_n = 0

    print("MCL PRISM structure registry v1.1")
    print(f"PRISM profiles in responses: {len(prism_ids)}")
    print(f"PRISM compounds in catalog: {len(universe)}")
    print(f"Cached registry rows available: {len(cached_map)}")
    print(f"Reusable generic MCL structures: {len(generic_map)}")

    with httpx.Client(
        timeout=args.timeout,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
        follow_redirects=True,
    ) as client:
        for ordinal, row in enumerate(universe.to_dict("records"), start=1):
            compound_id = _text(row.get("compound_id"))
            preferred_name = _text(row.get("preferred_name"))
            broad_id = _text(row.get("broad_id"))

            cached = cached_map.get(compound_id)
            if cached and not args.force:
                cached_status = _text(cached.get("resolution_status"))
                if cached_status == "resolved" or not args.retry_unresolved:
                    rows.append(cached)
                    reused_n += 1
                    if ordinal % max(1, args.checkpoint_every) == 0:
                        _write_checkpoint(rows)
                        print(
                            f"  {ordinal}/{len(universe)} processed; resolved="
                            f"{sum(_text(x.get('resolution_status')) == 'resolved' for x in rows)}; requests={requests_n}"
                        )
                    continue

            source_smiles = _text(row.get("canonical_smiles"))
            canonical, inchikey, error = _standardize(source_smiles)
            if canonical and inchikey:
                rows.append(
                    {
                        "compound_id": compound_id,
                        "preferred_name": preferred_name,
                        "broad_id": broad_id,
                        "canonical_smiles": canonical,
                        "inchikey": inchikey,
                        "pubchem_cid_resolved": _text(row.get("pubchem_cid")),
                        "resolution_status": "resolved",
                        "resolution_source": "prism_source_catalog",
                        "resolution_method": "existing_canonical_smiles",
                        "resolution_query": "",
                        "resolution_note": "Structure present in normalized PRISM source metadata.",
                        "retrieved_at": _now(),
                    }
                )
                source_n += 1
            else:
                generic = generic_map.get(compound_id)
                generic_smiles = _text(generic.get("canonical_smiles")) if generic else ""
                gcanonical, ginchikey, gerror = _standardize(generic_smiles)
                if gcanonical and ginchikey:
                    rows.append(
                        {
                            "compound_id": compound_id,
                            "preferred_name": preferred_name,
                            "broad_id": broad_id,
                            "canonical_smiles": gcanonical,
                            "inchikey": ginchikey,
                            "pubchem_cid_resolved": _text(generic.get("pubchem_cid_resolved")),
                            "resolution_status": "resolved",
                            "resolution_source": "mcl_compound_structure_registry",
                            "resolution_method": _text(generic.get("resolution_method")) or "cached_structure",
                            "resolution_query": _text(generic.get("resolution_query")),
                            "resolution_note": "Reused previously resolved MCL pharmacology structure.",
                            "retrieved_at": _now(),
                        }
                    )
                    generic_n += 1
                else:
                    resolved_payload: dict[str, Any] | None = None
                    used_type = ""
                    used_value = ""
                    final_status = "not_found"
                    notes: list[str] = []
                    for query_type, value in _candidate_queries(row):
                        requests_n += 1
                        status, payload, detail = _request_properties(
                            client,
                            query_type,
                            value,
                            retries=max(1, args.retries),
                            pause_s=max(0.0, args.pause),
                        )
                        notes.append(f"{query_type}:{status}:{detail}")
                        if status == "resolved" and payload:
                            qcanonical, qinchikey, qerror = _standardize(payload.get("canonical_smiles"))
                            if qcanonical and qinchikey:
                                resolved_payload = {
                                    **payload,
                                    "canonical_smiles": qcanonical,
                                    "inchikey": qinchikey,
                                }
                                used_type = query_type
                                used_value = value
                                final_status = "resolved"
                                break
                            notes.append(f"{query_type}:rdkit:{qerror}")
                        elif status in {"ambiguous", "request_failed"}:
                            final_status = status
                        time.sleep(max(0.0, args.pause))

                    if resolved_payload:
                        rows.append(
                            {
                                "compound_id": compound_id,
                                "preferred_name": preferred_name,
                                "broad_id": broad_id,
                                "canonical_smiles": _text(resolved_payload.get("canonical_smiles")),
                                "inchikey": _text(resolved_payload.get("inchikey")),
                                "pubchem_cid_resolved": _text(resolved_payload.get("pubchem_cid_resolved")),
                                "resolution_status": "resolved",
                                "resolution_source": "pubchem_pug_rest",
                                "resolution_method": used_type,
                                "resolution_query": used_value,
                                "resolution_note": " | ".join(notes),
                                "retrieved_at": _now(),
                            }
                        )
                        pubchem_n += 1
                    else:
                        rows.append(
                            {
                                "compound_id": compound_id,
                                "preferred_name": preferred_name,
                                "broad_id": broad_id,
                                "canonical_smiles": "",
                                "inchikey": "",
                                "pubchem_cid_resolved": "",
                                "resolution_status": final_status,
                                "resolution_source": "pubchem_pug_rest",
                                "resolution_method": "",
                                "resolution_query": "",
                                "resolution_note": " | ".join(notes) or error or gerror or "no_resolvable_identifier",
                                "retrieved_at": _now(),
                            }
                        )
                        unresolved_n += 1

            if ordinal % max(1, args.checkpoint_every) == 0:
                _write_checkpoint(rows)
                print(
                    f"  {ordinal}/{len(universe)} processed; resolved="
                    f"{sum(_text(x.get('resolution_status')) == 'resolved' for x in rows)}; requests={requests_n}"
                )

    registry = pd.DataFrame(rows).drop_duplicates("compound_id", keep="last").reset_index(drop=True)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    registry.to_parquet(OUTPUT, index=False, compression="zstd")
    registry.to_csv(OUTPUT_TSV, sep="\t", index=False)

    unresolved = registry[registry["resolution_status"].astype(str).ne("resolved")].copy()
    unresolved.to_csv(QC_OUTPUT, sep="\t", index=False)

    resolved = registry[registry["resolution_status"].astype(str).eq("resolved")].copy()
    manifest = {
        "contract": "mcl-prism-structure-registry-v1.1",
        "built_at": _now(),
        "prism_compounds_n": int(len(universe)),
        "registry_rows_n": int(len(registry)),
        "resolved_n": int(len(resolved)),
        "unresolved_n": int(len(unresolved)),
        "unique_resolved_inchikeys_n": int(resolved["inchikey"].replace("", pd.NA).dropna().nunique()),
        "requests_n_this_run": int(requests_n),
        "reused_cached_n": int(reused_n),
        "resolved_from_prism_source_n": int(source_n),
        "resolved_from_existing_mcl_registry_n": int(generic_n),
        "resolved_from_pubchem_n": int(pubchem_n),
        "standardization": "RDKit Cleanup -> FragmentParent -> canonical isomeric SMILES -> InChIKey",
        "scientific_guardrail_ru": (
            "PubChem используется только для восстановления химической идентичности PRISM-профиля. "
            "Совпадение по структуре не является доказательством механизма клеточного ответа."
        ),
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print("\nPRISM structure resolution completed")
    print(f"Resolved: {len(resolved)}/{len(registry)}")
    print(f"Unique resolved structures: {manifest['unique_resolved_inchikeys_n']}")
    print(f"  source metadata: {source_n}")
    print(f"  existing MCL registry: {generic_n}")
    print(f"  PubChem this run: {pubchem_n}")
    print(f"  cached rows reused: {reused_n}")
    print(f"Unresolved: {len(unresolved)}")
    print(f"PubChem requests this run: {requests_n}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")
    print(f"QC: {QC_OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
