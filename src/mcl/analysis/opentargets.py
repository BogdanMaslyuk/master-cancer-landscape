from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from mcl.models.opentargets import CancerTargetPairSeed, OpenTargetsAssociation
from mcl.models.provenance import ProvenanceRecord
from mcl.normalize.disease_mapping import load_opentargets_disease_mappings
from mcl.sources.opentargets import OpenTargetsClient, PARSER_VERSION, save_raw_response, utc_now_iso
from mcl.utils.hash import sha256_bytes


def load_wave1_pairs(path: str | Path) -> list[CancerTargetPairSeed]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return [
        CancerTargetPairSeed(
            evidence_id=r["Evidence_ID"],
            cancer_id=r["Cancer_ID"],
            target_id=r["Target_ID"],
            cancer_name=r["Cancer_name"],
            molecular_context=r["Molecular_context"],
            hgnc_symbol=r["HGNC_symbol"],
        )
        for r in rows
    ]


def load_target_identifier_index(path: str | Path) -> dict[str, dict[str, str]]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {r["target_id"]: r for r in rows}


def _score_dict(rows: list[dict[str, Any]] | None) -> dict[str, float]:
    out: dict[str, float] = {}
    for row in rows or []:
        if row.get("id") is not None and row.get("score") is not None:
            out[str(row["id"])] = float(row["score"])
    return out


def _extract_score(scores: dict[str, float], *keys: str) -> float | None:
    for key in keys:
        if key in scores:
            return scores[key]
    return None


def _platform_release_from_meta(payload: dict[str, Any]) -> dict[str, str | None]:
    meta = payload.get("data", {}).get("meta") or {}
    api = meta.get("apiVersion") or {}
    data = meta.get("dataVersion") or {}
    api_version = None
    if api:
        api_version = ".".join(str(api.get(k)) for k in ("x", "y", "z"))
    data_version = None
    if data:
        data_version = f"{data.get('year')}.{str(data.get('month')).zfill(2)}"
    return {"api_version": api_version, "data_version": data_version}


def _exact_disease_hits(payload: dict[str, Any], expected_name: str) -> list[dict[str, Any]]:
    mappings = payload.get("data", {}).get("mapIds", {}).get("mappings") or []
    expected = expected_name.strip().casefold()
    hits: list[dict[str, Any]] = []
    for mapping in mappings:
        for hit in mapping.get("hits") or []:
            if str(hit.get("name") or "").strip().casefold() == expected:
                hits.append(hit)
    unique: dict[str, dict[str, Any]] = {}
    for hit in hits:
        if hit.get("id"):
            unique[str(hit["id"])] = hit
    return list(unique.values())


def resolve_disease_for_release(
    *,
    client: OpenTargetsClient,
    root: Path,
    raw_root: Path,
    cancer_id: str,
    mapping,
) -> dict[str, Any]:
    """Resolve a configured disease mapping against the current Open Targets release.

    The configured ID is tried first. If that identifier no longer exists in the current
    Platform release, `mapIds` is used on the curated expected disease name. A fuzzy first
    hit is never accepted: exactly one case-insensitive exact-name hit is required.
    """
    configured_payload = client.fetch_disease(mapping.disease_id)
    configured_path = raw_root / f"{cancer_id}_{mapping.disease_id}_configured_disease.json"
    save_raw_response(configured_path, configured_payload)
    disease = configured_payload.get("data", {}).get("disease")

    if disease:
        return {
            "configured_id": mapping.disease_id,
            "resolved_id": disease.get("id"),
            "resolved_name": disease.get("name"),
            "resolution_method": "configured_id",
            "id_changed": disease.get("id") != mapping.disease_id,
            "raw_files": [str(configured_path.relative_to(root))],
            "resolution_error": None,
        }

    resolver_payload = client.resolve_disease_names([mapping.expected_disease_name])
    resolver_path = raw_root / f"{cancer_id}_disease_name_resolution.json"
    save_raw_response(resolver_path, resolver_payload)
    exact_hits = _exact_disease_hits(resolver_payload, mapping.expected_disease_name)

    if len(exact_hits) != 1:
        return {
            "configured_id": mapping.disease_id,
            "resolved_id": None,
            "resolved_name": None,
            "resolution_method": "unresolved",
            "id_changed": False,
            "raw_files": [
                str(configured_path.relative_to(root)),
                str(resolver_path.relative_to(root)),
            ],
            "resolution_error": (
                f"Expected exactly one exact-name Open Targets disease hit for "
                f"{mapping.expected_disease_name!r}; found {len(exact_hits)}"
            ),
        }

    candidate = exact_hits[0]
    resolved_id = str(candidate["id"])
    confirm_payload = client.fetch_disease(resolved_id)
    confirm_path = raw_root / f"{cancer_id}_{resolved_id}_resolved_disease.json"
    save_raw_response(confirm_path, confirm_payload)
    confirmed = confirm_payload.get("data", {}).get("disease")
    if not confirmed:
        return {
            "configured_id": mapping.disease_id,
            "resolved_id": None,
            "resolved_name": None,
            "resolution_method": "unresolved",
            "id_changed": False,
            "raw_files": [
                str(configured_path.relative_to(root)),
                str(resolver_path.relative_to(root)),
                str(confirm_path.relative_to(root)),
            ],
            "resolution_error": f"Resolved candidate {resolved_id} did not round-trip through disease(efoId:)",
        }

    return {
        "configured_id": mapping.disease_id,
        "resolved_id": confirmed.get("id"),
        "resolved_name": confirmed.get("name"),
        "resolution_method": "mapIds_exact_name",
        "id_changed": confirmed.get("id") != mapping.disease_id,
        "raw_files": [
            str(configured_path.relative_to(root)),
            str(resolver_path.relative_to(root)),
            str(confirm_path.relative_to(root)),
        ],
        "resolution_error": None,
    }


def collect_opentargets_wave1(
    *,
    root: str | Path,
    source_release: str,
    endpoint: str,
    page_size: int = 500,
) -> tuple[list[OpenTargetsAssociation], list[ProvenanceRecord], dict[str, Any]]:
    root = Path(root)
    pairs = load_wave1_pairs(root / "data/input/cancer_target_pairs_wave1.tsv")
    target_idx = load_target_identifier_index(root / "data/processed/target_identifiers.tsv")
    disease_mappings = load_opentargets_disease_mappings(root / "config/cancer_contexts.yaml")
    client = OpenTargetsClient(endpoint=endpoint)
    retrieved_at = utc_now_iso()

    required_cancers = sorted({p.cancer_id for p in pairs})
    missing_mappings = [x for x in required_cancers if x not in disease_mappings]
    if missing_mappings:
        raise ValueError(f"Missing Open Targets disease mapping(s): {missing_mappings}")

    missing_targets = [p.target_id for p in pairs if p.target_id not in target_idx]
    if missing_targets:
        raise ValueError(
            "Target normalization output is missing target(s): " + ", ".join(sorted(set(missing_targets)))
        )

    raw_root = root / "data/raw/opentargets" / source_release
    raw_root.mkdir(parents=True, exist_ok=True)

    meta_payload = client.fetch_meta()
    meta_path = raw_root / "platform_meta.json"
    save_raw_response(meta_path, meta_payload)
    platform_meta = _platform_release_from_meta(meta_payload)

    # Use the maintained single-target endpoint and cache one raw response per target.
    ensembl_ids = sorted({target_idx[p.target_id]["ensembl_gene_id"] for p in pairs})
    tract_idx: dict[str, list[dict[str, Any]]] = {}
    tract_files: dict[str, str] = {}
    for ensembl_id in ensembl_ids:
        tract_payload = client.fetch_target_tractability(ensembl_id)
        tract_path = raw_root / f"target_{ensembl_id}_tractability.json"
        save_raw_response(tract_path, tract_payload)
        target = tract_payload.get("data", {}).get("target")
        if target:
            tract_idx[ensembl_id] = target.get("tractability") or []
        else:
            tract_idx[ensembl_id] = []
        tract_files[ensembl_id] = str(tract_path.relative_to(root))

    disease_rows: dict[str, dict[str, Any]] = {}
    raw_file_for_cancer: dict[str, str] = {}
    disease_meta: dict[str, dict[str, Any]] = {}

    resolved_mappings: dict[str, dict[str, Any]] = {}
    for cancer_id in required_cancers:
        mapping = disease_mappings[cancer_id]
        resolved = resolve_disease_for_release(
            client=client,
            root=root,
            raw_root=raw_root,
            cancer_id=cancer_id,
            mapping=mapping,
        )
        resolved_mappings[cancer_id] = resolved
        if not resolved.get("resolved_id"):
            raise ValueError(
                f"Open Targets disease mapping unresolved for {cancer_id}: {resolved.get('resolution_error')}"
            )

        resolved_id = str(resolved["resolved_id"])
        page_index = 0
        all_rows: list[dict[str, Any]] = []
        disease_name: str | None = None
        disease_id: str | None = None
        count: int | None = None
        page_files: list[str] = []

        while True:
            payload = client.fetch_associated_targets_page(resolved_id, page_index, page_size)
            page_path = raw_root / f"{cancer_id}_{resolved_id}_associated_targets_page_{page_index:04d}.json"
            save_raw_response(page_path, payload)
            page_files.append(str(page_path.relative_to(root)))
            disease = payload.get("data", {}).get("disease")
            if not disease:
                raise ValueError(f"Open Targets returned no disease for resolved ID {resolved_id}")
            disease_id = disease.get("id")
            disease_name = disease.get("name")
            assoc = disease.get("associatedTargets") or {}
            if count is None:
                count = int(assoc.get("count") or 0)
            page_rows = assoc.get("rows") or []
            all_rows.extend(page_rows)
            if len(all_rows) >= count or not page_rows:
                break
            page_index += 1
            if page_index > 200:
                raise RuntimeError(f"Pagination safety limit exceeded for {resolved_id}")

        disease_meta[cancer_id] = {
            "configured_id": mapping.disease_id,
            "resolved_id": disease_id,
            "returned_id": disease_id,
            "returned_name": disease_name,
            "expected_name": mapping.expected_disease_name,
            "resolution_method": resolved.get("resolution_method"),
            "id_changed": bool(resolved.get("id_changed")),
            "count": str(count or 0),
            "resolution_raw_files": resolved.get("raw_files") or [],
        }
        raw_file_for_cancer[cancer_id] = json.dumps(page_files, ensure_ascii=False)
        disease_rows[cancer_id] = {
            row.get("target", {}).get("id"): row
            for row in all_rows
            if row.get("target", {}).get("id")
        }

    associations: list[OpenTargetsAssociation] = []
    provenance: list[ProvenanceRecord] = []

    for cancer_id in required_cancers:
        mapping = disease_mappings[cancer_id]
        resolved = resolved_mappings[cancer_id]
        resolved_id = str(resolved["resolved_id"])
        for field_name, value, nature in (
            ("Configured_Open_Targets_disease_ID", mapping.disease_id, "Inferred"),
            ("Resolved_Open_Targets_disease_ID", resolved_id, "Observed"),
            ("Open_Targets_resolution_method", resolved["resolution_method"], "Observed"),
            ("Open_Targets_mapping_scope", mapping.mapping_scope.value, "Inferred"),
            ("Open_Targets_mapping_confidence", mapping.mapping_confidence.value, "Inferred"),
            ("Open_Targets_mapping_rationale", mapping.rationale, "Inferred"),
        ):
            provenance.append(
                ProvenanceRecord(
                    entity_type="cancer_context",
                    entity_id=cancer_id,
                    field_name=field_name,
                    value=str(value),
                    source_name="Open Targets disease mapping",
                    source_record_id=resolved_id,
                    source_release=source_release,
                    source_url_or_endpoint=endpoint,
                    retrieved_at=retrieved_at,
                    evidence_nature=nature,
                    parser_version=PARSER_VERSION,
                    raw_file=json.dumps(resolved.get("raw_files") or [], ensure_ascii=False),
                    raw_record_hash=sha256_bytes(json.dumps(resolved, sort_keys=True).encode("utf-8")),
                )
            )

    for pair in pairs:
        mapping = disease_mappings[pair.cancer_id]
        target = target_idx[pair.target_id]
        ensembl_id = target["ensembl_gene_id"]
        resolved_id = str(resolved_mappings[pair.cancer_id]["resolved_id"])
        found = disease_rows[pair.cancer_id].get(ensembl_id)
        datatype_scores = _score_dict(found.get("datatypeScores") if found else None)
        datasource_scores = _score_dict(found.get("datasourceScores") if found else None)
        tractability = tract_idx.get(ensembl_id, [])
        returned_name = str(disease_meta[pair.cancer_id].get("returned_name") or mapping.expected_disease_name)
        raw_assoc_files = raw_file_for_cancer[pair.cancer_id]
        raw_hash = sha256_bytes(json.dumps(found or {}, sort_keys=True).encode("utf-8"))

        row = OpenTargetsAssociation(
            evidence_id=pair.evidence_id,
            cancer_id=pair.cancer_id,
            target_id=pair.target_id,
            target_ensembl_id=ensembl_id,
            target_symbol=target["hgnc_symbol"],
            disease_id=resolved_id,
            disease_name=returned_name,
            expected_disease_name=mapping.expected_disease_name,
            molecular_context=pair.molecular_context,
            mapping_scope=mapping.mapping_scope,
            mapping_confidence=mapping.mapping_confidence,
            molecular_context_encoded_in_ot=mapping.molecular_context_encoded_in_ot,
            association_found=found is not None,
            association_score=float(found["score"]) if found and found.get("score") is not None else None,
            direct_association=True,
            genetic_association_score=_extract_score(datatype_scores, "genetic_association"),
            somatic_mutation_score=_extract_score(datatype_scores, "somatic_mutation"),
            known_drug_score=_extract_score(datatype_scores, "known_drug"),
            literature_score=_extract_score(datatype_scores, "literature"),
            datatype_scores_json=json.dumps(datatype_scores, ensure_ascii=False, sort_keys=True),
            datasource_scores_json=json.dumps(datasource_scores, ensure_ascii=False, sort_keys=True),
            tractability_json=json.dumps(tractability, ensure_ascii=False, sort_keys=True),
            source_release=source_release,
            api_endpoint=endpoint,
            retrieved_at=retrieved_at,
            raw_association_file=raw_assoc_files,
            raw_tractability_file=tract_files.get(ensembl_id),
            raw_record_hash=raw_hash,
            notes=(
                "Open Targets association is disease-level only; molecular context is retained separately "
                "and is not claimed to be encoded by the Open Targets disease ontology. "
                f"Configured disease ID {mapping.disease_id}; release-resolved disease ID {resolved_id}."
            ),
        )
        associations.append(row)

        for field_name in (
            "association_score",
            "genetic_association_score",
            "somatic_mutation_score",
            "known_drug_score",
            "literature_score",
            "tractability_json",
        ):
            value = getattr(row, field_name)
            provenance.append(
                ProvenanceRecord(
                    entity_type="cancer_target_pair",
                    entity_id=pair.evidence_id,
                    field_name=field_name,
                    value=None if value is None else str(value),
                    source_name="Open Targets Platform",
                    source_record_id=f"{ensembl_id}|{resolved_id}",
                    source_release=source_release,
                    source_url_or_endpoint=endpoint,
                    retrieved_at=retrieved_at,
                    evidence_nature="Aggregated evidence",
                    parser_version=PARSER_VERSION,
                    raw_file=raw_assoc_files if field_name != "tractability_json" else str(tract_files.get(ensembl_id) or ""),
                    raw_record_hash=raw_hash,
                )
            )

    return associations, provenance, {
        "platform": platform_meta,
        "platform_meta_file": str(meta_path.relative_to(root)),
        "diseases": disease_meta,
        "tractability_files": tract_files,
    }


def validate_opentargets_diseases(
    *, root: str | Path, source_release: str, endpoint: str
) -> tuple[list[dict[str, Any]], list[ProvenanceRecord], dict[str, Any]]:
    root = Path(root)
    mappings = load_opentargets_disease_mappings(root / "config/cancer_contexts.yaml")
    client = OpenTargetsClient(endpoint=endpoint)
    raw_root = root / "data/raw/opentargets" / source_release
    raw_root.mkdir(parents=True, exist_ok=True)
    retrieved_at = utc_now_iso()
    validations: list[dict[str, Any]] = []
    provenance: list[ProvenanceRecord] = []

    meta_payload = client.fetch_meta()
    meta_path = raw_root / "platform_meta.json"
    save_raw_response(meta_path, meta_payload)
    platform_meta = _platform_release_from_meta(meta_payload)

    for cancer_id, mapping in mappings.items():
        resolved = resolve_disease_for_release(
            client=client,
            root=root,
            raw_root=raw_root,
            cancer_id=cancer_id,
            mapping=mapping,
        )
        returned_id = resolved.get("resolved_id")
        returned_name = resolved.get("resolved_name")
        name_match = (
            str(returned_name or "").strip().casefold()
            == mapping.expected_disease_name.strip().casefold()
        )
        if not returned_id:
            status = "ERROR"
        elif not name_match:
            status = "WARNING"
        elif bool(resolved.get("id_changed")):
            status = "WARNING"
        else:
            status = "PASS"

        validations.append({
            "cancer_id": cancer_id,
            "cancer_name": mapping.cancer_name,
            "molecular_context": mapping.molecular_context,
            "configured_disease_id": mapping.disease_id,
            "expected_disease_name": mapping.expected_disease_name,
            "resolved_disease_id": returned_id,
            "returned_disease_name": returned_name,
            "configured_id_still_current": returned_id == mapping.disease_id if returned_id else False,
            "name_match": name_match,
            "resolution_method": resolved.get("resolution_method"),
            "id_changed": bool(resolved.get("id_changed")),
            "mapping_scope": mapping.mapping_scope.value,
            "mapping_confidence": mapping.mapping_confidence.value,
            "molecular_context_encoded_in_ot": mapping.molecular_context_encoded_in_ot,
            "status": status,
            "resolution_error": resolved.get("resolution_error"),
            "rationale": mapping.rationale,
            "raw_files": json.dumps(resolved.get("raw_files") or [], ensure_ascii=False),
        })
        provenance.append(ProvenanceRecord(
            entity_type="cancer_context",
            entity_id=cancer_id,
            field_name="Open_Targets_disease_validation",
            value=json.dumps({
                "configured_id": mapping.disease_id,
                "resolved_id": returned_id,
                "name": returned_name,
                "resolution_method": resolved.get("resolution_method"),
            }, ensure_ascii=False),
            source_name="Open Targets Platform",
            source_record_id=str(returned_id or mapping.disease_id),
            source_release=source_release,
            source_url_or_endpoint=endpoint,
            retrieved_at=retrieved_at,
            evidence_nature="Observed",
            parser_version=PARSER_VERSION,
            raw_file=json.dumps(resolved.get("raw_files") or [], ensure_ascii=False),
            raw_record_hash=sha256_bytes(json.dumps(resolved, sort_keys=True).encode("utf-8")),
        ))
    return validations, provenance, {
        "platform": platform_meta,
        "platform_meta_file": str(meta_path.relative_to(root)),
    }
