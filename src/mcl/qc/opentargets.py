from __future__ import annotations

import json
from collections import Counter

from mcl.models.opentargets import OpenTargetsAssociation
from mcl.models.qc import QCRecord


def opentargets_qc(rows: list[OpenTargetsAssociation], expected_n: int) -> list[QCRecord]:
    qc: list[QCRecord] = []

    if len(rows) != expected_n:
        qc.append(QCRecord(
            severity="ERROR", check="row_count", entity_id="opentargets",
            field="rows", observed=str(len(rows)), expected=str(expected_n),
            message="Open Targets output row count does not match Wave 1 pair count",
        ))

    keys = [(r.cancer_id, r.target_id) for r in rows]
    duplicates = [k for k, n in Counter(keys).items() if n > 1]
    for cancer_id, target_id in duplicates:
        qc.append(QCRecord(
            severity="ERROR", check="unique_cancer_target", entity_id=f"{cancer_id}|{target_id}",
            field="Cancer_ID × Target_ID", observed="duplicate", expected="unique",
            message="Duplicate Cancer × Target pair in Open Targets output",
        ))

    for r in rows:
        if not r.target_ensembl_id.startswith("ENSG"):
            qc.append(QCRecord(
                severity="ERROR", check="target_identifier", entity_id=r.target_id,
                field="target_ensembl_id", observed=r.target_ensembl_id, expected="ENSG...",
                message="Target is not mapped to a valid Ensembl gene identifier",
            ))
        if r.association_score is not None and not (0.0 <= r.association_score <= 1.0):
            qc.append(QCRecord(
                severity="ERROR", check="association_score_range", entity_id=r.evidence_id,
                field="association_score", observed=str(r.association_score), expected="0..1",
                message="Open Targets association score is outside expected range",
            ))
        if r.molecular_context_encoded_in_ot:
            qc.append(QCRecord(
                severity="WARNING", check="molecular_context_claim", entity_id=r.evidence_id,
                field="molecular_context_encoded_in_ot", observed="True", expected="False for Wave 1 mappings",
                message="Review molecular-context ontology claim before using association as subtype-specific evidence",
            ))
        if not r.association_found:
            qc.append(QCRecord(
                severity="INFO", check="association_absent", entity_id=r.evidence_id,
                field="association_found", observed="False", expected="context-dependent",
                message="No direct disease-level Open Targets association row was found; this is not treated as zero evidence",
            ))
        try:
            json.loads(r.datatype_scores_json)
            json.loads(r.datasource_scores_json)
            json.loads(r.tractability_json)
        except json.JSONDecodeError:
            qc.append(QCRecord(
                severity="ERROR", check="json_serialization", entity_id=r.evidence_id,
                field="OpenTargets JSON fields", observed="invalid JSON", expected="valid JSON",
                message="Serialized Open Targets evidence cannot be parsed",
            ))

    if not [x for x in qc if x.severity in {"ERROR", "WARNING"}]:
        qc.append(QCRecord(
            severity="INFO", check="opentargets", entity_id="", field="", observed="", expected="",
            message="All Open Targets structural QC checks passed",
        ))
    return qc


def opentargets_disease_mapping_qc(disease_meta: dict, mappings: dict) -> list[QCRecord]:
    qc: list[QCRecord] = []
    for cancer_id, mapping in mappings.items():
        meta = (disease_meta or {}).get(cancer_id)
        if not meta:
            qc.append(QCRecord(
                severity="ERROR", check="disease_mapping_response", entity_id=cancer_id,
                field="disease", observed="missing", expected=mapping.expected_disease_name,
                message="No Open Targets disease response metadata was captured",
            ))
            continue

        resolved_id = str(meta.get("resolved_id") or meta.get("returned_id") or "")
        returned_name = str(meta.get("returned_name") or "")
        if not resolved_id:
            qc.append(QCRecord(
                severity="ERROR", check="disease_resolution", entity_id=cancer_id,
                field="disease_id", observed="unresolved", expected=mapping.expected_disease_name,
                message="Disease could not be resolved in the current Open Targets release",
            ))
            continue

        if returned_name.strip().casefold() != mapping.expected_disease_name.strip().casefold():
            qc.append(QCRecord(
                severity="WARNING", check="disease_name_validation", entity_id=cancer_id,
                field="disease_name", observed=returned_name, expected=mapping.expected_disease_name,
                message="Open Targets canonical disease name differs from the curated expected name; manual ontology review is required",
            ))

        if resolved_id != mapping.disease_id:
            qc.append(QCRecord(
                severity="WARNING", check="disease_identifier_drift", entity_id=cancer_id,
                field="disease_id", observed=resolved_id, expected=mapping.disease_id,
                message=(
                    "Configured disease identifier is no longer the current Open Targets identifier. "
                    "The run-time exact-name resolver supplied a release-specific canonical ID; "
                    "the configured historical ID is retained for provenance."
                ),
            ))
    return qc
