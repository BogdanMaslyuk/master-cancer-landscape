"""Conservative normalization of pinned ChEMBL / BindingDB assay extracts.

P0.3 pilot contract. Not a direct-binding prediction, and not a generic parser
for an unmodified ChEMBL SQL dump. The ChEMBL extract must contain joined
target/assay/activity fields; missing target provenance fails closed.
"""
from __future__ import annotations

import hashlib
import json
import re
from decimal import Decimal, InvalidOperation
from typing import Any

MISSING = {"", "nan", "none", "null"}
UNIT_MULTIPLIER = {
    "m": Decimal("1000000000"),
    "mm": Decimal("1000000"),
    "um": Decimal("1000"),
    "nm": Decimal("1"),
    "pm": Decimal("0.001"),
}
VALUE_PATTERN = re.compile(r"^\s*(<=|>=|<|>|=|≤|≥)?\s*(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*$")
PRIMARY_FIELDS = (
    "source", "source_release", "source_record_id", "raw_row_locator",
    "source_target_id", "target_gene", "uniprot_accession", "target_type",
    "organism_id", "molecule_id", "original_smiles", "endpoint", "relation",
    "original_value", "original_units", "numeric_value_nm", "assay_id",
    "assay_type", "assay_format", "confidence_score", "relationship_type",
    "document_id", "doi", "pmid", "source_origin", "original_chembl_activity_id",
    "quality_tier", "measurement_status", "reason_codes", "raw_file_sha256",
)


def clean(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in MISSING else text


def field(row: dict, *keys: str) -> str:
    for key in keys:
        value = clean(row.get(key))
        if value:
            return value
    return ""


def measurement(value: str, units: str, separate_relation: str = "") -> tuple[str, str, str]:
    """Return operator, numeric value in nM, or a machine-readable failure."""
    raw = clean(value).replace("µ", "u").replace("μ", "u")
    match = VALUE_PATTERN.fullmatch(raw)
    if not match:
        return "", "", "invalid_numeric_value"
    embedded = (match.group(1) or "=").replace("≤", "<=").replace("≥", ">=")
    explicit = clean(separate_relation).replace("≤", "<=").replace("≥", ">=")
    if explicit and explicit not in {"=", "<", ">", "<=", ">="}:
        return "", "", "invalid_relation"
    if explicit and embedded != "=" and embedded != explicit:
        return "", "", "conflicting_relations"
    relation = explicit or embedded
    normalized_unit = clean(units).replace("µ", "u").replace("μ", "u").lower()
    if normalized_unit not in UNIT_MULTIPLIER:
        return "", "", "unknown_units"
    try:
        number = Decimal(match.group(2)) * UNIT_MULTIPLIER[normalized_unit]
    except (InvalidOperation, ValueError):
        return "", "", "invalid_numeric_value"
    if not number.is_finite() or number <= 0:
        return "", "", "nonpositive_value"
    plain = format(number, "f")
    if "." in plain:
        plain = plain.rstrip("0").rstrip(".")
    return relation, plain, ""


def cutoff_state(relation: str, value_nm: str, cutoff_nm: str = "10000") -> str:
    number, cutoff = Decimal(value_nm), Decimal(cutoff_nm)
    if relation == "=":
        return "exact_eligible" if number <= cutoff else "outside_cutoff"
    if relation in {"<", "<="}:
        return "censored_eligible" if number <= cutoff else "censored_indeterminate"
    return "censored_indeterminate"


def _base(source: str, release: str, sha: str, locator: str) -> dict:
    return dict.fromkeys(PRIMARY_FIELDS, "") | {
        "source": source, "source_release": release,
        "raw_file_sha256": sha, "raw_row_locator": locator,
        "quality_tier": "D", "measurement_status": "excluded",
    }


def normalize_chembl(
    row: dict, *, release: str, sha: str, locator: str,
    registry: dict[str, str], cutoff_nm: str = "10000",
) -> dict:
    """Normalize a joined target/assay/activity row; never guess missing joins."""
    out = _base("ChEMBL", release, sha, locator)
    out.update({
        "source_record_id": field(row, "activity_id"),
        "source_target_id": field(row, "target_chembl_id"),
        "uniprot_accession": field(row, "uniprot_accession"),
        "target_type": field(row, "target_type"),
        "organism_id": field(row, "organism_id"),
        "molecule_id": field(row, "molecule_chembl_id"),
        "original_smiles": field(row, "canonical_smiles"),
        "endpoint": field(row, "standard_type"),
        "original_value": field(row, "standard_value"),
        "original_units": field(row, "standard_units"),
        "assay_id": field(row, "assay_chembl_id"),
        "assay_type": field(row, "assay_type"),
        "assay_format": field(row, "bao_label"),
        "confidence_score": field(row, "confidence_score"),
        "relationship_type": field(row, "relationship_type"),
        "document_id": field(row, "document_chembl_id"),
        "doi": field(row, "doi"),
        "pmid": field(row, "pubmed_id"),
        "source_origin": "chembl",
    })
    reasons = []
    if not out["source_record_id"] or not out["assay_id"] or not out["source_target_id"]:
        reasons.append("missing_primary_identifier")
    accession = out["uniprot_accession"]
    if accession not in registry:
        reasons.append("target_not_in_pinned_registry")
    else:
        out["target_gene"] = registry[accession]
    if out["organism_id"] != "9606" or out["target_type"].upper() != "SINGLE PROTEIN":
        reasons.append("non_unique_human_single_protein_target")
    if out["endpoint"].upper() not in {"KI", "KD", "IC50"}:
        reasons.append("endpoint_not_supported")
    operator, value_nm, issue = measurement(
        out["original_value"], out["original_units"], field(row, "standard_relation")
    )
    if issue:
        reasons.append(issue)
    else:
        out["relation"], out["numeric_value_nm"] = operator, value_nm
        state = cutoff_state(operator, value_nm, cutoff_nm)
        out["measurement_status"] = state
        if state not in {"exact_eligible", "censored_eligible"}:
            reasons.append(state)
    if field(row, "data_validity_comment"):
        reasons.append("data_validity_comment")
    if field(row, "potential_duplicate").lower() in {"true", "1", "yes"}:
        reasons.append("potential_duplicate_requires_review")
    if out["confidence_score"] != "9" or out["relationship_type"].upper() != "D":
        reasons.append("target_confidence_below_strict")
    assay = out["assay_type"].upper()
    bao = out["assay_format"].lower()
    if assay not in {"B", "F"} or "single protein format" not in bao:
        reasons.append("assay_not_verified_single_protein")
    if assay != "B":
        reasons.append("assay_type_requires_review")
    if not out["original_smiles"]:
        reasons.append("missing_structure")
    if reasons:
        out["quality_tier"] = "C" if out["target_gene"] and out["original_smiles"] else "D"
    elif out["measurement_status"] == "censored_eligible":
        out["quality_tier"] = "C"  # never promote a bound to exact primary data
    else:
        out["quality_tier"] = "candidate_A" if out["endpoint"].upper() in {"KI", "KD"} else "candidate_B"
    out["reason_codes"] = json.dumps(sorted(set(reasons)), ensure_ascii=False)
    return out


def normalize_bindingdb(
    row: dict, *, endpoint: str, release: str, sha: str, locator: str,
    registry: dict[str, str], cutoff_nm: str = "10000",
) -> dict:
    """BindingDB bulk TSV: preserve source evidence; absent assay context => C."""
    out = _base("BindingDB", release, sha, locator)
    out.update({
        "source_record_id": field(row, "BindingDB Reactant_set_id"),
        "source_target_id": field(row, "UniProt (SwissProt) Primary ID of Target Chain"),
        "uniprot_accession": field(row, "UniProt (SwissProt) Primary ID of Target Chain"),
        "target_type": "unverified",
        "organism_id": field(row, "Target Source Organism According to Curator or DataSource"),
        "molecule_id": field(row, "BindingDB MonomerID"),
        "original_smiles": field(row, "Ligand SMILES"),
        "endpoint": endpoint,
        "original_value": field(row, endpoint + " (nM)"),
        "original_units": "nM",
        "assay_id": field(row, "BindingDB Assay ID", "Assay ID"),
        "assay_type": "unverified",
        "assay_format": "unverified",
        "doi": field(row, "Article DOI"),
        "pmid": field(row, "PMID"),
        "document_id": field(row, "BindingDB Entry DOI"),
        "source_origin": field(row, "Curation/DataSource"),
        "original_chembl_activity_id": field(row, "ChEMBL Activity ID"),
    })
    reasons = []
    if not out["source_record_id"] or not out["molecule_id"]:
        reasons.append("missing_primary_identifier")
    acc = out["uniprot_accession"]
    if acc not in registry:
        reasons.append("target_not_in_pinned_registry")
    else:
        out["target_gene"] = registry[acc]
    # BindingDB 'Target Source Organism...' may be blank or human spelled out.
    if out["organism_id"].lower() not in {"homo sapiens", "human", "9606"}:
        reasons.append("organism_not_verified")
    if not out["original_smiles"]:
        reasons.append("missing_structure")
    operator, value_nm, issue = measurement(out["original_value"], "nM")
    if issue:
        reasons.append(issue)
    else:
        out["relation"], out["numeric_value_nm"] = operator, value_nm
        out["measurement_status"] = cutoff_state(operator, value_nm, cutoff_nm)
        if out["measurement_status"] != "exact_eligible":
            reasons.append(out["measurement_status"])
    # Official bulk TSV does not necessarily establish measured construct or
    # direct binding. Require separate assay provenance review to upgrade.
    reasons.append("assay_method_requires_review")
    out["quality_tier"] = "C" if out["target_gene"] and out["original_smiles"] else "D"
    out["reason_codes"] = json.dumps(sorted(set(reasons)), ensure_ascii=False)
    return out


def mirror_group_key(row: dict) -> str:
    """Only explicit original activity IDs prove cross-provider mirroring."""
    chembl_id = clean(row.get("original_chembl_activity_id"))
    if row.get("source") == "ChEMBL":
        chembl_id = clean(row.get("source_record_id"))
    if not chembl_id:
        return ""
    return hashlib.sha256(("chembl-activity:" + chembl_id).encode()).hexdigest()[:16]
