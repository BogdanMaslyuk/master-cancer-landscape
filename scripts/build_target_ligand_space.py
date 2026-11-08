from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import httpx
import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
except ImportError as exc:
    raise SystemExit(
        'RDKit is required. Install chemistry extras: .\\.venv\\Scripts\\python.exe -m pip install -e ".[chemistry]"'
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
PHARM = ROOT / "data" / "runtime" / "pharmacology"
CANDIDATES = PHARM / "candidate_hypotheses_v2.parquet"
TARGET_REGISTRY = ROOT / "data" / "processed" / "target_registry" / "protein_targets.parquet"
CACHE = ROOT / "data" / "external" / "target_ligand_space" / "v1"
RUNTIME = ROOT / "data" / "runtime" / "target_ligand_space"
QC_DIR = ROOT / "outputs" / "qc"

TARGET_MAP_OUT = RUNTIME / "target_map.tsv"
MEASUREMENTS_OUT = RUNTIME / "ligand_measurements.parquet"
CATALOG_OUT = RUNTIME / "ligand_catalog.parquet"
CATALOG_TSV_OUT = RUNTIME / "ligand_catalog.tsv"
TARGET_SUMMARY_OUT = RUNTIME / "target_summary.tsv"
MANIFEST_OUT = RUNTIME / "target_ligand_space_manifest.json"
QC_OUT = QC_DIR / "target_ligand_space_qc.tsv"

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
BINDINGDB_URL = "https://www.bindingdb.org/rest/getLigandsByUniprots"
PRIORITY_STATUS = "priority_for_in_vitro"
ENDPOINTS = {"KI": "Ki", "KD": "Kd", "IC50": "IC50"}
USER_AGENT = "Master-Cancer-Landscape/0.5 Target-Ligand-Space-v1 (research)"


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


def _float(value: Any) -> float | None:
    try:
        if value is None or pd.isna(value):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    tmp.replace(path)


def _request_json(
    client: httpx.Client,
    url: str,
    *,
    params: dict[str, Any] | None = None,
    retries: int = 5,
) -> Any:
    last_error = ""
    for attempt in range(1, retries + 1):
        try:
            response = client.get(url, params=params)
        except httpx.HTTPError as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            if attempt < retries:
                time.sleep(min(2**attempt, 20))
                continue
            raise RuntimeError(last_error) from exc

        if response.status_code == 429 or response.status_code >= 500:
            last_error = f"HTTP {response.status_code}"
            retry_after = _float(response.headers.get("Retry-After"))
            if attempt < retries:
                time.sleep(retry_after if retry_after is not None else min(2**attempt, 30))
                continue
        response.raise_for_status()
        if not response.text.strip():
            return {}
        return response.json()
    raise RuntimeError(last_error or "request_failed")


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
            return None, None, "missing_inchikey_after_standardization"
        return canonical, inchikey, None
    except Exception as exc:
        return None, None, f"standardization_failed:{type(exc).__name__}"


def _load_target_map() -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    missing = [p for p in (CANDIDATES, TARGET_REGISTRY) if not p.exists()]
    if missing:
        raise SystemExit("Missing Target Ligand Space inputs:\n" + "\n".join(str(p.relative_to(ROOT)) for p in missing))

    candidates = pd.read_parquet(CANDIDATES)
    required = {"priority_status", "target_gene"}
    if not required.issubset(candidates.columns):
        raise SystemExit(f"Candidate v2 must contain {sorted(required)}")
    priority = candidates[candidates["priority_status"].astype(str).eq(PRIORITY_STATUS)].copy()
    priority["target_gene"] = priority["target_gene"].astype(str).str.upper()

    registry = pd.read_parquet(TARGET_REGISTRY)
    needed_registry = {"target_gene", "protein_mapping_status", "uniprot_primary_accession"}
    if not needed_registry.issubset(registry.columns):
        raise SystemExit(f"Protein target registry must contain {sorted(needed_registry)}")
    registry = registry.copy()
    registry["target_gene"] = registry["target_gene"].astype(str).str.upper()

    qc: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []
    for gene, group in priority.groupby("target_gene", sort=True):
        matches = registry[registry["target_gene"].eq(gene)].drop_duplicates("target_gene")
        accession = ""
        mapping_status = "missing_registry_row"
        protein_name = ""
        if len(matches) == 1:
            row = matches.iloc[0]
            mapping_status = _text(row.get("protein_mapping_status"))
            accession = _text(row.get("uniprot_primary_accession"))
            protein_name = _text(row.get("protein_preferred_name"))
        usable = mapping_status == "unique_swissprot" and bool(accession)
        if not usable:
            qc.append(
                {
                    "stage": "target_mapping",
                    "target_gene": gene,
                    "entity_id": gene,
                    "status": "unresolved",
                    "detail": f"protein_mapping_status={mapping_status}; accession={accession or 'missing'}",
                }
            )
        rows.append(
            {
                "target_gene": gene,
                "uniprot_accession": accession,
                "protein_preferred_name": protein_name,
                "protein_mapping_status": mapping_status,
                "usable_for_ligand_query": usable,
                "priority_hypotheses_n": int(len(group)),
                "priority_cancers_n": int(group["mcl_cancer_id"].nunique()) if "mcl_cancer_id" in group.columns else 0,
                "priority_cancers_json": _json_unique(group["mcl_cancer_name"] if "mcl_cancer_name" in group.columns else []),
            }
        )
    return pd.DataFrame(rows), qc


def _chembl_target(
    client: httpx.Client,
    gene: str,
    accession: str,
    *,
    force: bool,
) -> tuple[str | None, dict[str, Any] | None, str]:
    cache_file = CACHE / "chembl" / f"{gene}__target.json"
    if cache_file.exists() and not force:
        payload = json.loads(cache_file.read_text(encoding="utf-8"))
    else:
        payload = _request_json(
            client,
            f"{CHEMBL_BASE}/target.json",
            params={"target_components__accession": accession, "limit": 100},
        )
        _atomic_json(cache_file, payload)

    candidates = payload.get("targets", []) if isinstance(payload, dict) else []
    exact: list[dict[str, Any]] = []
    for item in candidates:
        if _text(item.get("target_type")).upper() != "SINGLE PROTEIN":
            continue
        if _text(item.get("organism")) != "Homo sapiens":
            continue
        components = item.get("target_components") or []
        accessions = {_text(c.get("accession")) for c in components if isinstance(c, dict)}
        if accession in accessions:
            exact.append(item)
    if len(exact) != 1:
        return None, None, f"exact_human_single_protein_targets={len(exact)}"
    item = exact[0]
    chembl_id = _text(item.get("target_chembl_id"))
    return chembl_id or None, item, "ok" if chembl_id else "missing_target_chembl_id"


def _chembl_activities(
    client: httpx.Client,
    gene: str,
    target_chembl_id: str,
    *,
    cutoff_nm: float,
    force: bool,
) -> list[dict[str, Any]]:
    cache_file = CACHE / "chembl" / f"{gene}__activities_{int(cutoff_nm)}nM.json"
    if cache_file.exists() and not force:
        payload = json.loads(cache_file.read_text(encoding="utf-8"))
        return payload.get("activities", []) if isinstance(payload, dict) else []

    activities: list[dict[str, Any]] = []
    offset = 0
    limit = 1000
    while True:
        payload = _request_json(
            client,
            f"{CHEMBL_BASE}/activity.json",
            params={
                "target_chembl_id": target_chembl_id,
                "standard_type__in": "IC50,Ki,Kd",
                "standard_units": "nM",
                "standard_value__lte": cutoff_nm,
                "limit": limit,
                "offset": offset,
            },
        )
        page = payload.get("activities", []) if isinstance(payload, dict) else []
        if not page:
            break
        activities.extend(page)
        page_meta = payload.get("page_meta", {}) if isinstance(payload, dict) else {}
        total = int(page_meta.get("total_count") or 0)
        offset += len(page)
        if len(page) < limit or (total and offset >= total):
            break
        time.sleep(0.08)
    _atomic_json(
        cache_file,
        {
            "retrieved_at": _now(),
            "target_chembl_id": target_chembl_id,
            "cutoff_nm": cutoff_nm,
            "activities": activities,
        },
    )
    return activities


def _normalize_chembl(
    gene: str,
    accession: str,
    target_chembl_id: str,
    activities: list[dict[str, Any]],
    cutoff_nm: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    out: list[dict[str, Any]] = []
    qc: list[dict[str, Any]] = []
    for activity in activities:
        endpoint = ENDPOINTS.get(_text(activity.get("standard_type")).upper())
        relation = _text(activity.get("standard_relation"))
        units = _text(activity.get("standard_units"))
        value = _float(activity.get("standard_value"))
        smiles = _text(activity.get("canonical_smiles"))
        bao = _text(activity.get("bao_label")).lower()
        assay_type = _text(activity.get("assay_type")).upper()
        validity = _text(activity.get("data_validity_comment"))
        duplicate = _text(activity.get("potential_duplicate"))

        reason = ""
        if endpoint is None:
            reason = "endpoint_not_allowed"
        elif relation != "=":
            reason = f"non_exact_relation:{relation or 'missing'}"
        elif units.lower() != "nm":
            reason = f"units_not_nm:{units or 'missing'}"
        elif value is None or value <= 0 or value > cutoff_nm:
            reason = "value_outside_cutoff"
        elif not smiles:
            reason = "missing_smiles"
        elif validity:
            reason = f"data_validity_comment:{validity}"
        elif duplicate in {"1", "True", "true"}:
            reason = "potential_duplicate"
        else:
            excluded_format = any(token in bao for token in ("cell-based", "organism-based", "tissue-based"))
            direct_format = "single protein format" in bao
            if excluded_format:
                reason = f"non_direct_bao:{bao}"
            elif endpoint in {"Ki", "Kd"}:
                if not (direct_format or assay_type == "B"):
                    reason = f"direct_context_not_resolved:bao={bao or 'missing'};assay_type={assay_type or 'missing'}"
            elif endpoint == "IC50":
                if not (direct_format or assay_type == "B"):
                    reason = f"ic50_not_direct_binding_context:bao={bao or 'missing'};assay_type={assay_type or 'missing'}"

        if reason:
            qc.append(
                {
                    "stage": "chembl_filter",
                    "target_gene": gene,
                    "entity_id": _text(activity.get("activity_id")),
                    "status": "excluded",
                    "detail": reason,
                }
            )
            continue

        out.append(
            {
                "source": "ChEMBL",
                "target_gene": gene,
                "uniprot_accession": accession,
                "source_target_id": target_chembl_id,
                "source_ligand_id": _text(activity.get("molecule_chembl_id")),
                "raw_smiles": smiles,
                "endpoint_type": endpoint,
                "relation": relation,
                "value_nm": float(value),
                "evidence_class": "direct_single_protein_measurement",
                "assay_id": _text(activity.get("assay_chembl_id")),
                "assay_type": assay_type,
                "bao_label": _text(activity.get("bao_label")),
                "assay_description": _text(activity.get("assay_description")),
                "document_id": _text(activity.get("document_chembl_id")),
                "pchembl_value": _float(activity.get("pchembl_value")),
                "source_record_id": _text(activity.get("activity_id")),
                "source_url": f"https://www.ebi.ac.uk/chembl/explore/target/{target_chembl_id}",
            }
        )
    return out, qc


def _normalized_key(key: Any) -> str:
    return re.sub(r"[^a-z0-9]", "", str(key).lower())


def _record_field(record: dict[str, Any], *suffixes: str) -> Any:
    wanted = tuple(_normalized_key(x) for x in suffixes)
    for key, value in record.items():
        normalized = _normalized_key(key)
        if any(normalized == suffix or normalized.endswith(suffix) for suffix in wanted):
            return value
    return None


def _bindingdb_candidate_records(payload: Any) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []

    def walk(node: Any) -> None:
        if isinstance(node, dict):
            has_smiles = _record_field(node, "smile", "smiles") is not None
            has_affinity = _record_field(node, "affinity") is not None
            has_type = _record_field(node, "affinitytype") is not None
            if has_smiles and has_affinity and has_type:
                records.append(node)
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            for value in node:
                walk(value)

    walk(payload)
    return records


def _parse_affinity(value: Any) -> tuple[str, float | None]:
    text = _text(value).replace(",", "").strip()
    if not text:
        return "", None
    match = re.match(r"^\s*(<=|>=|<|>|=)?\s*([0-9]*\.?[0-9]+)", text)
    if not match:
        return "", None
    relation = match.group(1) or "="
    number = _float(match.group(2))
    return relation, number


def _bindingdb_measurements(
    client: httpx.Client,
    gene: str,
    accession: str,
    *,
    cutoff_nm: float,
    force: bool,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    cache_file = CACHE / "bindingdb" / f"{gene}__{int(cutoff_nm)}nM.json"
    if cache_file.exists() and not force:
        payload = json.loads(cache_file.read_text(encoding="utf-8"))
    else:
        payload = _request_json(
            client,
            BINDINGDB_URL,
            params={"uniprot": accession, "cutoff": int(cutoff_nm), "response": "application/json"},
        )
        _atomic_json(cache_file, payload)
        time.sleep(0.8)

    out: list[dict[str, Any]] = []
    qc: list[dict[str, Any]] = []
    records = _bindingdb_candidate_records(payload)
    uniprot_pattern = re.compile(r"\b[OPQ][0-9][A-Z0-9]{3}[0-9]\b")
    for idx, record in enumerate(records):
        endpoint_raw = _text(_record_field(record, "affinitytype"))
        endpoint = ENDPOINTS.get(endpoint_raw.upper())
        relation, value = _parse_affinity(_record_field(record, "affinity"))
        smiles = _text(_record_field(record, "smile", "smiles"))
        query = _text(_record_field(record, "query", "uniprot"))
        monomer = _text(_record_field(record, "monomerid", "monomer"))

        reason = ""
        parsed_uniprots = set(uniprot_pattern.findall(query.upper())) if query else set()
        if parsed_uniprots and accession not in parsed_uniprots:
            reason = f"related_target_not_exact:{query}"
        elif endpoint is None:
            reason = f"endpoint_not_allowed:{endpoint_raw or 'missing'}"
        elif relation in {">", ">="}:
            reason = f"lower_bound_not_usable:{relation}"
        elif value is None or value <= 0 or value > cutoff_nm:
            reason = "value_outside_cutoff_or_unreported"
        elif not smiles:
            reason = "missing_smiles"

        if reason:
            qc.append(
                {
                    "stage": "bindingdb_filter",
                    "target_gene": gene,
                    "entity_id": monomer or f"row-{idx}",
                    "status": "excluded",
                    "detail": reason,
                }
            )
            continue

        out.append(
            {
                "source": "BindingDB",
                "target_gene": gene,
                "uniprot_accession": accession,
                "source_target_id": accession,
                "source_ligand_id": monomer,
                "raw_smiles": smiles,
                "endpoint_type": endpoint,
                "relation": relation,
                "value_nm": float(value),
                "evidence_class": "direct_bindingdb_affinity",
                "assay_id": "",
                "assay_type": "binding",
                "bao_label": "",
                "assay_description": "BindingDB getLigandsByUniprots direct affinity record",
                "document_id": "",
                "pchembl_value": None,
                "source_record_id": monomer or f"{gene}:{idx}",
                "source_url": f"https://www.bindingdb.org/rwd/bind/ByUniProtids.jsp?uniprot={accession}",
            }
        )
    return out, qc


def _ligand_id(target_gene: str, inchikey: str) -> str:
    digest = hashlib.sha1(f"{target_gene}|{inchikey}".encode("utf-8")).hexdigest()[:14].upper()
    return f"TLS-{digest}"


def _build_catalog(measurements: pd.DataFrame, qc: list[dict[str, Any]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    standardized_rows: list[dict[str, Any]] = []
    for record in measurements.to_dict("records"):
        canonical, inchikey, error = _standardize_smiles(record.get("raw_smiles"))
        if error:
            qc.append(
                {
                    "stage": "structure_standardization",
                    "target_gene": _text(record.get("target_gene")),
                    "entity_id": _text(record.get("source_ligand_id")),
                    "status": "excluded",
                    "detail": error,
                }
            )
            continue
        record["canonical_smiles"] = canonical
        record["inchikey"] = inchikey
        standardized_rows.append(record)

    standardized = pd.DataFrame(standardized_rows)
    if standardized.empty:
        return standardized, pd.DataFrame()

    catalog_rows: list[dict[str, Any]] = []
    for (gene, inchikey), group in standardized.groupby(["target_gene", "inchikey"], sort=True):
        ordered = group.sort_values(["value_nm", "source", "source_ligand_id"], ascending=[True, True, True])
        best = ordered.iloc[0]
        catalog_rows.append(
            {
                "ligand_id": _ligand_id(str(gene), str(inchikey)),
                "target_gene": str(gene),
                "uniprot_accession": _text(best.get("uniprot_accession")),
                "inchikey": str(inchikey),
                "canonical_smiles": _text(best.get("canonical_smiles")),
                "sources_json": _json_unique(group["source"]),
                "source_ligand_ids_json": _json_unique(group["source_ligand_id"]),
                "endpoint_types_json": _json_unique(group["endpoint_type"]),
                "measurements_n": int(len(group)),
                "chembl_measurements_n": int(group["source"].astype(str).eq("ChEMBL").sum()),
                "bindingdb_measurements_n": int(group["source"].astype(str).eq("BindingDB").sum()),
                "lowest_reported_value_nm_across_endpoints": float(best["value_nm"]),
                "lowest_value_endpoint_type": _text(best.get("endpoint_type")),
                "lowest_value_relation": _text(best.get("relation")),
                "lowest_value_source": _text(best.get("source")),
                "built_at": _now(),
            }
        )
    catalog = pd.DataFrame(catalog_rows).sort_values(["target_gene", "ligand_id"]).reset_index(drop=True)
    return standardized, catalog


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build Target Ligand Space v1 for current MCL priority targets from ChEMBL and BindingDB."
    )
    parser.add_argument("--cutoff-nm", type=float, default=10000.0)
    parser.add_argument("--force-refresh", action="store_true")
    parser.add_argument("--timeout", type=float, default=120.0)
    args = parser.parse_args()
    if args.cutoff_nm <= 0:
        raise SystemExit("--cutoff-nm must be positive")

    RUNTIME.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    target_map, qc = _load_target_map()
    target_map.to_csv(TARGET_MAP_OUT, sep="\t", index=False)

    measurements: list[dict[str, Any]] = []
    with httpx.Client(
        timeout=args.timeout,
        follow_redirects=True,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
    ) as client:
        for row in target_map.to_dict("records"):
            gene = _text(row.get("target_gene"))
            accession = _text(row.get("uniprot_accession"))
            if not bool(row.get("usable_for_ligand_query")):
                continue

            print(f"[{gene}] UniProt {accession}")
            try:
                chembl_id, _, status = _chembl_target(
                    client, gene, accession, force=args.force_refresh
                )
                if chembl_id:
                    activities = _chembl_activities(
                        client,
                        gene,
                        chembl_id,
                        cutoff_nm=args.cutoff_nm,
                        force=args.force_refresh,
                    )
                    accepted, chembl_qc = _normalize_chembl(
                        gene, accession, chembl_id, activities, args.cutoff_nm
                    )
                    measurements.extend(accepted)
                    qc.extend(chembl_qc)
                    print(f"  ChEMBL: fetched={len(activities)}; accepted direct={len(accepted)}; target={chembl_id}")
                else:
                    qc.append(
                        {
                            "stage": "chembl_target_mapping",
                            "target_gene": gene,
                            "entity_id": accession,
                            "status": "unresolved",
                            "detail": status,
                        }
                    )
                    print(f"  ChEMBL: target unresolved ({status})")
            except Exception as exc:
                qc.append(
                    {
                        "stage": "chembl_request",
                        "target_gene": gene,
                        "entity_id": accession,
                        "status": "request_failed",
                        "detail": f"{type(exc).__name__}: {exc}",
                    }
                )
                print(f"  ChEMBL: ERROR {type(exc).__name__}: {exc}")

            try:
                bdb, bdb_qc = _bindingdb_measurements(
                    client,
                    gene,
                    accession,
                    cutoff_nm=args.cutoff_nm,
                    force=args.force_refresh,
                )
                measurements.extend(bdb)
                qc.extend(bdb_qc)
                print(f"  BindingDB: accepted direct={len(bdb)}")
            except Exception as exc:
                qc.append(
                    {
                        "stage": "bindingdb_request",
                        "target_gene": gene,
                        "entity_id": accession,
                        "status": "request_failed",
                        "detail": f"{type(exc).__name__}: {exc}",
                    }
                )
                print(f"  BindingDB: ERROR {type(exc).__name__}: {exc}")

    raw_measurements = pd.DataFrame(measurements)
    if raw_measurements.empty:
        raise SystemExit("No direct ligand measurements were retrieved from ChEMBL or BindingDB.")

    standardized, catalog = _build_catalog(raw_measurements, qc)
    if catalog.empty:
        raise SystemExit("No ligand structures survived standardization.")

    standardized.to_parquet(MEASUREMENTS_OUT, index=False, compression="zstd")
    catalog.to_parquet(CATALOG_OUT, index=False, compression="zstd")
    catalog.to_csv(CATALOG_TSV_OUT, sep="\t", index=False)

    summary_rows: list[dict[str, Any]] = []
    for row in target_map.to_dict("records"):
        gene = _text(row.get("target_gene"))
        m = standardized[standardized["target_gene"].astype(str).eq(gene)] if not standardized.empty else pd.DataFrame()
        c = catalog[catalog["target_gene"].astype(str).eq(gene)]
        summary_rows.append(
            {
                "target_gene": gene,
                "uniprot_accession": _text(row.get("uniprot_accession")),
                "priority_hypotheses_n": int(row.get("priority_hypotheses_n") or 0),
                "measurements_n": int(len(m)),
                "chembl_measurements_n": int(m["source"].astype(str).eq("ChEMBL").sum()) if not m.empty else 0,
                "bindingdb_measurements_n": int(m["source"].astype(str).eq("BindingDB").sum()) if not m.empty else 0,
                "unique_ligands_n": int(len(c)),
                "chembl_unique_ligands_n": int(c["sources_json"].astype(str).str.contains("ChEMBL", regex=False).sum()) if not c.empty else 0,
                "bindingdb_unique_ligands_n": int(c["sources_json"].astype(str).str.contains("BindingDB", regex=False).sum()) if not c.empty else 0,
            }
        )
    target_summary = pd.DataFrame(summary_rows)
    target_summary.to_csv(TARGET_SUMMARY_OUT, sep="\t", index=False)

    pd.DataFrame(qc, columns=["stage", "target_gene", "entity_id", "status", "detail"]).to_csv(
        QC_OUT, sep="\t", index=False
    )

    manifest = {
        "contract": "mcl-target-ligand-space-v1",
        "built_at": _now(),
        "priority_status": PRIORITY_STATUS,
        "targets_n": int(len(target_map)),
        "queryable_targets_n": int(target_map["usable_for_ligand_query"].astype(bool).sum()),
        "affinity_cutoff_nm": args.cutoff_nm,
        "allowed_endpoints": ["Ki", "Kd", "IC50"],
        "chembl_contract": {
            "base": CHEMBL_BASE,
            "target": "exact UniProt component, Homo sapiens, SINGLE PROTEIN",
            "relation": "= only",
            "units": "nM",
            "directness": "single protein BAO or binding assay type; cell/organism/tissue formats excluded",
        },
        "bindingdb_contract": {
            "endpoint": BINDINGDB_URL,
            "query": "one UniProt accession per request",
            "cutoff_nm": args.cutoff_nm,
            "relations": ["=", "<", "<="],
            "note": "Rows explicitly identifying another UniProt accession are excluded; zero/unreported affinities are excluded.",
        },
        "standardization": "RDKit Cleanup -> FragmentParent -> canonical isomeric SMILES -> InChIKey",
        "deduplication": "target_gene + standardized InChIKey",
        "measurements_n": int(len(standardized)),
        "unique_target_ligands_n": int(len(catalog)),
        "scientific_guardrail": (
            "The ligand space is an experimentally annotated chemical neighborhood, not proof that a PYZ compound binds the same target. "
            "Ki, Kd and IC50 are retained as separate endpoint types and must not be treated as interchangeable potency measurements."
        ),
    }
    MANIFEST_OUT.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print("\nMCL Target Ligand Space v1")
    print(f"Priority targets: {len(target_map)}; queryable by unique Swiss-Prot: {int(target_map['usable_for_ligand_query'].astype(bool).sum())}")
    print(f"Accepted direct measurements: {len(standardized)}")
    print(f"Unique target × ligand structures: {len(catalog)}")
    print("Per-target ligand counts:")
    for row in target_summary.to_dict("records"):
        print(
            f"  {row['target_gene']}: ligands={row['unique_ligands_n']}; "
            f"measurements={row['measurements_n']}; ChEMBL={row['chembl_measurements_n']}; BindingDB={row['bindingdb_measurements_n']}"
        )
    print(f"Wrote {CATALOG_OUT.relative_to(ROOT)}")
    print(f"Summary: {TARGET_SUMMARY_OUT.relative_to(ROOT)}")
    print(f"QC: {QC_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
