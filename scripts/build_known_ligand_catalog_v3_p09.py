#!/usr/bin/env python3
"""MCL TLS V3/P0.9: build a reviewable known-ligand evidence catalog.

The program reads only the validated P0.8 report ZIP. It does not open the
multi-gigabyte ChEMBL or BindingDB sources, does not infer activity for PYZ,
and never upgrades scientific acceptance automatically.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import zipfile
from collections import Counter, defaultdict
from collections.abc import Iterable
from datetime import UTC, datetime
from decimal import Decimal, InvalidOperation
from pathlib import Path, PurePosixPath
from typing import Any

CONTRACT = "mcl-tls-v3-p09-known-ligand-catalog-v1"
P08_PACKAGE_CONTRACT = "mcl-tls-v3-p08-combined-validated-package-v1"
P08_CONTRACT = "mcl-tls-v3-p08-local-source-validation-v2"
FOCUS6 = {
    "PRKDC": "P78527",
    "BCL2L1": "Q07817",
    "IGF1R": "P08069",
    "ERBB3": "P21860",
    "FNTB": "P49356",
    "DDC": "P20711",
}
ENDPOINTS = ("Ki", "Kd", "IC50")
SHA256_RE = re.compile(r"[0-9a-f]{64}")

P08_REQUIRED_MEMBERS = {
    "validated_package_manifest.json",
    "p08_summary.json",
    "p08_combined_validation.json",
    "p08_source_manifest.json",
    "p08_target_panel.tsv",
    "p08_chembl_observations.tsv",
    "p08_bindingdb_observations.tsv",
    "p08_cross_source_potential_duplicates.tsv",
}

CHEMBL_FIELDS = [
    "gene_symbol",
    "uniprot_accession",
    "chembl_target_id",
    "activity_id",
    "assay_chembl_id",
    "molecule_chembl_id",
    "canonical_smiles",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    "assay_type",
    "confidence_score",
    "potential_duplicate",
    "data_validity_comment",
    "document_chembl_id",
    "doi",
    "pmid",
    "measurement_value_status",
    "evidence_status",
    "primary_experiment_review_status",
    "cross_source_potential_duplicate",
    "scientific_acceptance",
]

BINDINGDB_FIELDS = [
    "gene_symbol",
    "uniprot_accession",
    "source_row",
    "bindingdb_reactant_set_id",
    "bindingdb_monomer_id",
    "source_row_field_count",
    "expected_field_count",
    "source_row_width_status",
    "raw_chain_accessions",
    "all_chain_accessions",
    "populated_chain_count",
    "declared_chain_count",
    "target_name",
    "organism",
    "raw_ki_nm",
    "raw_kd_nm",
    "raw_ic50_nm",
    "endpoint",
    "raw_value_nm",
    "relation",
    "numeric_value_nm",
    "ligand_smiles",
    "ligand_inchikey",
    "article_doi",
    "pmid",
    "curation_source",
    "chain_status",
    "species_status",
    "measurement_value_status",
    "measurement_status",
    "primary_experiment_review_status",
    "cross_source_potential_duplicate",
    "scientific_acceptance",
]

CROSS_FIELDS = [
    "gene_symbol",
    "doi",
    "endpoint",
    "relation",
    "numeric_value_nm",
    "ligand_smiles",
    "chembl_activity_id",
    "bindingdb_source_row",
    "assessment",
]

ASSESSMENT_FIELDS = [
    "observation_id",
    "source",
    "gene_symbol",
    "uniprot_accession",
    "source_target_id",
    "source_measurement_id",
    "source_experiment_id",
    "source_molecule_id",
    "source_record_id",
    "source_document_id",
    "publication_key",
    "publication_url",
    "doi",
    "pmid",
    "endpoint",
    "relation",
    "numeric_value_nm",
    "raw_value_nm",
    "source_units",
    "ligand_smiles",
    "ligand_inchikey",
    "structure_identity_key",
    "protein_match_status",
    "species_status",
    "experiment_type_status",
    "measurement_precision_status",
    "structure_status",
    "publication_status",
    "curator_quality_status",
    "potential_duplicate_status",
    "within_source_duplicate_group_ids",
    "cross_source_duplicate_group_ids",
    "catalog_inclusion_status",
    "evidence_completeness_status",
    "warning_codes",
    "primary_experiment_review_status",
    "original_scientific_acceptance",
    "scientific_acceptance",
    "chembl_target_id",
    "activity_id",
    "assay_chembl_id",
    "molecule_chembl_id",
    "document_chembl_id",
    "assay_type",
    "confidence_score",
    "potential_duplicate",
    "data_validity_comment",
    "chembl_evidence_status",
    "bindingdb_source_row",
    "bindingdb_reactant_set_id",
    "bindingdb_monomer_id",
    "source_row_field_count",
    "expected_field_count",
    "source_row_width_status",
    "raw_chain_accessions",
    "all_chain_accessions",
    "populated_chain_count",
    "declared_chain_count",
    "target_name",
    "organism",
    "raw_ki_nm",
    "raw_kd_nm",
    "raw_ic50_nm",
    "curation_source",
    "bindingdb_chain_status",
    "bindingdb_species_status",
    "bindingdb_measurement_status",
]

CATALOG_FIELDS = [
    "gene_symbol",
    "uniprot_accession",
    "source",
    "source_molecule_id",
    "ligand_smiles",
    "ligand_inchikey",
    "endpoint",
    "relation",
    "numeric_value_nm",
    "raw_value_nm",
    "source_units",
    "doi",
    "pmid",
    "publication_key",
    "publication_url",
    "protein_match_status",
    "species_status",
    "experiment_type_status",
    "measurement_precision_status",
    "structure_status",
    "potential_duplicate_status",
    "catalog_inclusion_status",
    "evidence_completeness_status",
    "warning_codes",
    "observation_id",
    "source_target_id",
    "source_measurement_id",
    "source_experiment_id",
    "source_record_id",
    "source_document_id",
    "primary_experiment_review_status",
    "scientific_acceptance",
]

DUPLICATE_GROUP_FIELDS = [
    "duplicate_group_id",
    "scope",
    "gene_symbol",
    "endpoint",
    "match_basis",
    "member_observations_n",
    "distinct_source_records_n",
    "sources",
    "publication_key",
    "ligand_identity_key",
    "relation",
    "numeric_value_nm",
    "review_status",
]

DUPLICATE_MEMBER_FIELDS = [
    "duplicate_group_id",
    "scope",
    "observation_id",
    "source",
    "gene_symbol",
    "endpoint",
    "source_molecule_id",
    "source_measurement_id",
    "source_experiment_id",
    "source_record_id",
    "publication_key",
]

CROSS_PAIR_FIELDS = [
    "duplicate_pair_id",
    "duplicate_group_id",
    "chembl_observation_id",
    "bindingdb_observation_id",
] + CROSS_FIELDS

LIGAND_SUMMARY_FIELDS = [
    "gene_symbol",
    "uniprot_accession",
    "source",
    "source_molecule_id",
    "structure_identity_key",
    "ligand_smiles",
    "ligand_inchikey",
    "endpoint",
    "observation_records_n",
    "exact_numeric_n",
    "censored_n",
    "invalid_or_missing_n",
    "strict_candidate_n",
    "review_candidate_n",
    "context_only_n",
    "excluded_invalid_n",
    "highest_eligibility_status",
    "reported_measurements",
    "publication_keys",
    "warning_codes",
    "observation_ids",
    "scientific_acceptance",
]

TARGET_SUMMARY_FIELDS = [
    "gene_symbol",
    "uniprot_accession",
    "chembl_observation_records_n",
    "bindingdb_observation_records_n",
    "all_observation_records_n",
    "ki_observation_records_n",
    "kd_observation_records_n",
    "ic50_observation_records_n",
    "exact_numeric_observation_records_n",
    "censored_observation_records_n",
    "invalid_or_missing_observation_records_n",
    "exact_single_protein_human_records_n",
    "multi_component_or_unknown_protein_records_n",
    "species_missing_or_mismatch_records_n",
    "structure_present_records_n",
    "structure_missing_records_n",
    "reference_identifier_present_records_n",
    "reference_identifier_missing_records_n",
    "strict_candidate_records_n",
    "review_candidate_records_n",
    "context_only_records_n",
    "excluded_invalid_records_n",
    "unique_source_ligand_ids_n",
    "strict_candidate_source_ligand_ids_n",
    "within_source_duplicate_groups_n",
    "within_source_duplicate_member_records_n",
    "cross_source_duplicate_groups_n",
    "cross_source_duplicate_member_records_n",
    "priority_publications_n",
    "coverage_interpretation",
    "scientific_acceptance",
]

PRIORITY_PUBLICATION_FIELDS = [
    "review_queue_rank",
    "gene_symbol",
    "publication_key",
    "publication_url",
    "doi",
    "pmid",
    "selection_role",
    "selection_reason_code",
    "selection_rationale_ru",
    "sources",
    "endpoints",
    "observation_records_n",
    "exact_numeric_records_n",
    "censored_records_n",
    "strict_candidate_records_n",
    "review_candidate_records_n",
    "context_only_records_n",
    "potential_duplicate_member_records_n",
    "multi_component_records_n",
    "primary_review_status",
    "scientific_acceptance",
]

PRIORITY_LINK_FIELDS = [
    "review_queue_rank",
    "gene_symbol",
    "publication_key",
    "observation_id",
    "source",
    "endpoint",
    "catalog_inclusion_status",
    "warning_codes",
    "source_measurement_id",
    "source_experiment_id",
    "source_molecule_id",
]

QC_FIELDS = ["severity", "category", "gene_symbol", "affected_records_n", "finding_ru", "required_action_ru"]


def text(value: Any) -> str:
    return "" if value is None else str(value).strip()


def truthy(value: Any) -> bool:
    return text(value).lower() in {"1", "true", "yes"}


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def require_sha256(value: str, label: str) -> str:
    value = text(value)
    if not SHA256_RE.fullmatch(value):
        raise ValueError(f"{label} must be a lowercase SHA-256 digest")
    return value


def normalized_doi(value: Any) -> str:
    doi = text(value).lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        doi = doi.removeprefix(prefix)
    return doi.strip()


def publication_key(doi: Any, pmid: Any) -> str:
    doi_norm = normalized_doi(doi)
    if doi_norm:
        return f"DOI:{doi_norm}"
    pmid_norm = text(pmid)
    return f"PMID:{pmid_norm}" if pmid_norm else ""


def publication_url(doi: Any, pmid: Any) -> str:
    doi_norm = normalized_doi(doi)
    if doi_norm:
        return f"https://doi.org/{doi_norm}"
    pmid_norm = text(pmid)
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmid_norm}/" if pmid_norm else ""


def normalized_relation(value: Any) -> str:
    return {"≤": "<=", "≥": ">="}.get(text(value), text(value))


def normalized_numeric(value: Any) -> str:
    raw = text(value)
    try:
        number = Decimal(raw)
    except InvalidOperation:
        return ""
    if not number.is_finite() or number <= 0:
        return ""
    fixed = format(number.normalize(), "f")
    return fixed.rstrip("0").rstrip(".") if "." in fixed else fixed


class ValidatedZip:
    """Strict, read-only ZIP reader with bounded members and CRC verification."""

    def __init__(self, path: Path):
        self.path = path
        self.archive = zipfile.ZipFile(path, "r")
        infos = [item for item in self.archive.infolist() if not item.is_dir()]
        names = [item.filename for item in infos]
        if len(names) != len(set(names)):
            self.close()
            raise ValueError("P0.8 ZIP contains duplicate member names")
        if len(names) > 100:
            self.close()
            raise ValueError("P0.8 ZIP has too many members")
        for item in infos:
            posix = PurePosixPath(item.filename)
            mode = (item.external_attr >> 16) & 0o170000
            if (
                posix.is_absolute()
                or ".." in posix.parts
                or "\\" in item.filename
                or item.flag_bits & 0x1
                or mode == 0o120000
                or item.file_size > 100_000_000
            ):
                self.close()
                raise ValueError(f"Unsafe P0.8 ZIP member: {item.filename}")
        self.names = set(names)
        bad_member = self.archive.testzip()
        if bad_member:
            self.close()
            raise ValueError(f"P0.8 ZIP CRC failed for {bad_member}")

    def read(self, name: str) -> bytes:
        if name not in self.names:
            raise ValueError(f"P0.8 ZIP missing required member {name}")
        return self.archive.read(name)

    def close(self) -> None:
        self.archive.close()


def read_json(data: bytes, label: str) -> dict[str, Any]:
    value = json.loads(data.decode("utf-8-sig"))
    if not isinstance(value, dict):
        raise TypeError(f"{label} must contain a JSON object")
    return value


def read_tsv(data: bytes, expected_fields: list[str], label: str) -> list[dict[str, str]]:
    with io.StringIO(data.decode("utf-8-sig"), newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != expected_fields:
            raise ValueError(f"{label} header does not match the P0.8 contract")
        rows: list[dict[str, str]] = []
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{label} malformed TSV row {line_number}")
            rows.append({key: value for key, value in row.items()})
        return rows


def verify_p08_report(path: Path, expected_sha256: str) -> dict[str, Any]:
    expected_sha256 = require_sha256(expected_sha256, "Expected P0.8 ZIP SHA-256")
    actual_sha256 = sha_file(path)
    if actual_sha256 != expected_sha256:
        raise ValueError("P0.8 ZIP SHA-256 mismatch")

    archive = ValidatedZip(path)
    try:
        if not P08_REQUIRED_MEMBERS.issubset(archive.names):
            missing = sorted(P08_REQUIRED_MEMBERS - archive.names)
            raise ValueError(f"P0.8 ZIP lacks required members: {missing}")
        package = read_json(archive.read("validated_package_manifest.json"), "P0.8 package manifest")
        if package.get("contract") != P08_PACKAGE_CONTRACT:
            raise ValueError("P0.8 validated package contract mismatch")
        if package.get("validation_status") != "PASSED":
            raise ValueError("P0.8 validated package status is not PASSED")
        if package.get("scientific_acceptance") != "NOT_APPROVED" or package.get("merge_permitted") is not False:
            raise ValueError("P0.8 scientific acceptance safeguards mismatch")
        if package.get("source_releases") != {"chembl": "chembl_37", "bindingdb": "202611"}:
            raise ValueError("P0.8 source release mismatch")
        if package.get("source_input_hashes_unchanged") is not True:
            raise ValueError("P0.8 did not preserve source input hashes")

        manifest_files = package.get("files")
        if not isinstance(manifest_files, dict):
            raise TypeError("P0.8 package file manifest missing")
        if archive.names != set(manifest_files) | {"validated_package_manifest.json"}:
            raise ValueError("P0.8 ZIP members do not exactly match its package manifest")
        file_verification: dict[str, dict[str, Any]] = {}
        for name, recorded in manifest_files.items():
            if not isinstance(recorded, dict):
                raise TypeError(f"P0.8 manifest entry invalid for {name}")
            content = archive.read(name)
            actual = {"size_bytes": len(content), "sha256": sha_bytes(content)}
            if actual != recorded:
                raise ValueError(f"P0.8 member hash or size mismatch: {name}")
            file_verification[name] = actual

        summary = read_json(archive.read("p08_summary.json"), "P0.8 summary")
        if summary.get("contract") != P08_CONTRACT:
            raise ValueError("P0.8 report contract mismatch")
        if summary.get("panel_mode") != "focus6" or summary.get("panel_genes_n") != 6:
            raise ValueError("P0.8 report is not the focus6 panel")
        if summary.get("source_status") != {
            "chembl": "LOCAL_SOURCE_CHECKED",
            "bindingdb": "LOCAL_SOURCE_CHECKED",
        }:
            raise ValueError("P0.8 did not assess both local sources")
        if summary.get("releases") != {"chembl": "chembl_37", "bindingdb": "202611"}:
            raise ValueError("P0.8 summary release mismatch")
        if summary.get("scientific_acceptance") != "NOT_APPROVED" or summary.get("merge_permitted") is not False:
            raise ValueError("P0.8 summary acceptance safeguards mismatch")
        if summary.get("source_input_hashes_unchanged") is not True or summary.get("baseline_65x13_untouched") is not True:
            raise ValueError("P0.8 source or baseline immutability was not established")
        if summary.get("external_requests_n") != 0:
            raise ValueError("P0.8 summary reports external requests")
        if summary.get("bindingdb_source_rows_short_padded_n") != 0 or summary.get("bindingdb_source_rows_long_rejected_n") != 0:
            raise ValueError("P0.8 BindingDB contains structural-width exceptions")

        validation = read_json(archive.read("p08_combined_validation.json"), "P0.8 combined validation")
        required_validation = {
            "base_zip_crc": "PASSED",
            "base_zip_duplicate_members_n": 0,
            "base_zip_unsafe_members_n": 0,
            "base_zip_vs_directory_mismatches_n": 0,
            "sqlite_audit_status": "PASSED",
            "source_hashes_unchanged": True,
            "chembl_sqlite_value_mismatches_n": 0,
            "bindingdb_raw_endpoint_value_mismatches_n": 0,
            "truncated_observations_n": 0,
            "qc_rows_n": 0,
            "all_observations_not_approved": True,
            "technical_validation_status": "PASSED",
        }
        for key, expected in required_validation.items():
            if validation.get(key) != expected:
                raise ValueError(f"P0.8 combined validation mismatch: {key}")
        for source_counts in validation.get("missing_identifier_counts", {}).values():
            if any(value != 0 for value in source_counts.values()):
                raise ValueError("P0.8 contains missing required source identifiers")

        source_manifest = read_json(archive.read("p08_source_manifest.json"), "P0.8 source manifest")
        if source_manifest.get("contract") != P08_CONTRACT:
            raise ValueError("P0.8 source manifest contract mismatch")
        for key in ("no_network", "read_only_scientific_data", "exact_target_uniprot_required", "no_automatic_scientific_acceptance"):
            if source_manifest.get(key) is not True:
                raise ValueError(f"P0.8 source manifest safeguard mismatch: {key}")
        source_inputs = source_manifest.get("input_paths_sha256")
        if not isinstance(source_inputs, dict) or not source_inputs:
            raise ValueError("P0.8 source input hashes missing")
        for source_path, hashes in source_inputs.items():
            before = require_sha256(hashes.get("sha256_before", ""), f"P0.8 source hash before {source_path}")
            after = require_sha256(hashes.get("sha256_after", ""), f"P0.8 source hash after {source_path}")
            if before != after:
                raise ValueError(f"P0.8 source changed during analysis: {source_path}")

        panel = read_tsv(
            archive.read("p08_target_panel.tsv"),
            ["gene_symbol", "uniprot_accession", "mapping_status", "uniprot_hgnc", "sample_stratum", "source_identity_status"],
            "P0.8 target panel",
        )
        if {row["gene_symbol"]: row["uniprot_accession"] for row in panel} != FOCUS6:
            raise ValueError("P0.8 focus6 gene-to-UniProt panel mismatch")
        if any(row["source_identity_status"] != "VERIFIED_UNIPROT_HGNC" for row in panel):
            raise ValueError("P0.8 focus6 contains unresolved protein identity")

        chembl_rows = read_tsv(archive.read("p08_chembl_observations.tsv"), CHEMBL_FIELDS, "P0.8 ChEMBL observations")
        binding_rows = read_tsv(archive.read("p08_bindingdb_observations.tsv"), BINDINGDB_FIELDS, "P0.8 BindingDB observations")
        cross_rows = read_tsv(archive.read("p08_cross_source_potential_duplicates.tsv"), CROSS_FIELDS, "P0.8 cross-source duplicates")
        if len(chembl_rows) != summary.get("chembl_rows_reviewed_n") or len(chembl_rows) != validation.get("chembl_observation_rows_n"):
            raise ValueError("P0.8 ChEMBL observation denominator mismatch")
        if len(binding_rows) != summary.get("bindingdb_measurement_observations_retained_n") or len(binding_rows) != validation.get("bindingdb_observation_rows_n"):
            raise ValueError("P0.8 BindingDB observation denominator mismatch")
        if len(cross_rows) != summary.get("cross_source_potential_duplicate_pairs_n") or len(cross_rows) != validation.get("cross_source_potential_duplicate_pairs_n"):
            raise ValueError("P0.8 cross-source duplicate denominator mismatch")
        if any(row["scientific_acceptance"] != "NOT_APPROVED" for row in chembl_rows + binding_rows):
            raise ValueError("P0.8 observation acceptance state mismatch")
        if any(row["gene_symbol"] not in FOCUS6 or row["uniprot_accession"] != FOCUS6[row["gene_symbol"]] for row in chembl_rows + binding_rows):
            raise ValueError("P0.8 observation target identity mismatch")

        return {
            "actual_zip_sha256": actual_sha256,
            "zip_size_bytes": path.stat().st_size,
            "zip_crc_status": "PASSED",
            "duplicate_members_n": 0,
            "unsafe_members_n": 0,
            "manifest_member_files_n": len(manifest_files),
            "manifest_hash_or_size_mismatches_n": 0,
            "file_verification": file_verification,
            "package_manifest": package,
            "summary": summary,
            "validation": validation,
            "source_manifest": source_manifest,
            "panel": panel,
            "chembl_rows": chembl_rows,
            "binding_rows": binding_rows,
            "cross_rows": cross_rows,
        }
    finally:
        archive.close()


def structure_key(smiles: str, inchikey: str = "") -> str:
    return f"INCHIKEY:{text(inchikey)}" if text(inchikey) else f"SMILES:{text(smiles)}"


def build_observations(chembl_rows: list[dict[str, str]], binding_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    observations: list[dict[str, str]] = []
    for source in chembl_rows:
        observation = {field: "" for field in ASSESSMENT_FIELDS}
        observation.update(
            {
                "observation_id": f"CHEMBL:activity:{source['activity_id']}",
                "source": "ChEMBL",
                "gene_symbol": source["gene_symbol"],
                "uniprot_accession": source["uniprot_accession"],
                "source_target_id": source["chembl_target_id"],
                "source_measurement_id": source["activity_id"],
                "source_experiment_id": source["assay_chembl_id"],
                "source_molecule_id": source["molecule_chembl_id"],
                "source_record_id": source["activity_id"],
                "source_document_id": source["document_chembl_id"],
                "publication_key": publication_key(source["doi"], source["pmid"]),
                "publication_url": publication_url(source["doi"], source["pmid"]),
                "doi": source["doi"],
                "pmid": source["pmid"],
                "endpoint": source["standard_type"],
                "relation": source["standard_relation"],
                "numeric_value_nm": source["standard_value"],
                "raw_value_nm": source["standard_value"],
                "source_units": source["standard_units"],
                "ligand_smiles": source["canonical_smiles"],
                "structure_identity_key": structure_key(source["canonical_smiles"]),
                "protein_match_status": "EXACT_SINGLE_PROTEIN_UNIPROT_MATCH",
                "species_status": "HUMAN_SPECIES_MATCH",
                "experiment_type_status": (
                    "BINDING_ASSAY_METADATA_MATCHED"
                    if source["assay_type"] == "B" and source["confidence_score"] == "9"
                    else "REVIEW_ASSAY_METADATA_MISMATCH"
                ),
                "measurement_precision_status": source["measurement_value_status"],
                "structure_status": "STRUCTURE_PRESENT" if text(source["canonical_smiles"]) else "REVIEW_STRUCTURE_MISSING",
                "publication_status": (
                    "REFERENCE_IDENTIFIER_PRESENT_PRIMARY_STATUS_UNVERIFIED"
                    if publication_key(source["doi"], source["pmid"])
                    else "REVIEW_REFERENCE_IDENTIFIER_MISSING"
                ),
                "curator_quality_status": (
                    "REVIEW_CURATOR_FLAGGED"
                    if truthy(source["potential_duplicate"]) or text(source["data_validity_comment"])
                    else "NO_CURATOR_WARNING_RECORDED"
                ),
                "primary_experiment_review_status": source["primary_experiment_review_status"],
                "original_scientific_acceptance": source["scientific_acceptance"],
                "scientific_acceptance": "NOT_APPROVED",
                "chembl_target_id": source["chembl_target_id"],
                "activity_id": source["activity_id"],
                "assay_chembl_id": source["assay_chembl_id"],
                "molecule_chembl_id": source["molecule_chembl_id"],
                "document_chembl_id": source["document_chembl_id"],
                "assay_type": source["assay_type"],
                "confidence_score": source["confidence_score"],
                "potential_duplicate": source["potential_duplicate"],
                "data_validity_comment": source["data_validity_comment"],
                "chembl_evidence_status": source["evidence_status"],
            }
        )
        observations.append(observation)

    for source in binding_rows:
        observation = {field: "" for field in ASSESSMENT_FIELDS}
        chain_status = source["chain_status"]
        if chain_status == "SWISSPROT_PRIMARY_SINGLE_CHAIN_MATCH":
            protein_status = "EXACT_SINGLE_PROTEIN_UNIPROT_MATCH"
        elif chain_status == "REVIEW_MULTI_CHAIN":
            protein_status = "REVIEW_MULTI_COMPONENT_PROTEIN_OBJECT"
        else:
            protein_status = "REVIEW_CHAIN_COUNT_UNKNOWN"
        observation.update(
            {
                "observation_id": f"BINDINGDB:source_row:{source['source_row']}:endpoint:{source['endpoint']}",
                "source": "BindingDB",
                "gene_symbol": source["gene_symbol"],
                "uniprot_accession": source["uniprot_accession"],
                "source_target_id": source["uniprot_accession"],
                "source_measurement_id": f"{source['source_row']}:{source['endpoint']}",
                "source_experiment_id": source["bindingdb_reactant_set_id"],
                "source_molecule_id": source["bindingdb_monomer_id"],
                "source_record_id": source["source_row"],
                "publication_key": publication_key(source["article_doi"], source["pmid"]),
                "publication_url": publication_url(source["article_doi"], source["pmid"]),
                "doi": source["article_doi"],
                "pmid": source["pmid"],
                "endpoint": source["endpoint"],
                "relation": source["relation"],
                "numeric_value_nm": source["numeric_value_nm"],
                "raw_value_nm": source["raw_value_nm"],
                "source_units": "nM",
                "ligand_smiles": source["ligand_smiles"],
                "ligand_inchikey": source["ligand_inchikey"],
                "structure_identity_key": structure_key(source["ligand_smiles"], source["ligand_inchikey"]),
                "protein_match_status": protein_status,
                "species_status": source["species_status"],
                "experiment_type_status": "REPORTED_AFFINITY_ENDPOINT_ASSAY_UNVERIFIED",
                "measurement_precision_status": source["measurement_value_status"],
                "structure_status": "STRUCTURE_PRESENT" if text(source["ligand_smiles"]) else "REVIEW_STRUCTURE_MISSING",
                "publication_status": (
                    "REFERENCE_IDENTIFIER_PRESENT_PRIMARY_STATUS_UNVERIFIED"
                    if publication_key(source["article_doi"], source["pmid"])
                    else "REVIEW_REFERENCE_IDENTIFIER_MISSING"
                ),
                "curator_quality_status": "CURATOR_QUALITY_FLAG_NOT_PROVIDED",
                "primary_experiment_review_status": source["primary_experiment_review_status"],
                "original_scientific_acceptance": source["scientific_acceptance"],
                "scientific_acceptance": "NOT_APPROVED",
                "bindingdb_source_row": source["source_row"],
                "bindingdb_reactant_set_id": source["bindingdb_reactant_set_id"],
                "bindingdb_monomer_id": source["bindingdb_monomer_id"],
                "source_row_field_count": source["source_row_field_count"],
                "expected_field_count": source["expected_field_count"],
                "source_row_width_status": source["source_row_width_status"],
                "raw_chain_accessions": source["raw_chain_accessions"],
                "all_chain_accessions": source["all_chain_accessions"],
                "populated_chain_count": source["populated_chain_count"],
                "declared_chain_count": source["declared_chain_count"],
                "target_name": source["target_name"],
                "organism": source["organism"],
                "raw_ki_nm": source["raw_ki_nm"],
                "raw_kd_nm": source["raw_kd_nm"],
                "raw_ic50_nm": source["raw_ic50_nm"],
                "curation_source": source["curation_source"],
                "bindingdb_chain_status": source["chain_status"],
                "bindingdb_species_status": source["species_status"],
                "bindingdb_measurement_status": source["measurement_status"],
            }
        )
        observations.append(observation)

    identifiers = [row["observation_id"] for row in observations]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("P0.9 observation identifiers are not unique")
    return observations


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def find(self, item: str) -> str:
        self.parent.setdefault(item, item)
        if self.parent[item] != item:
            self.parent[item] = self.find(self.parent[item])
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self.parent[max(left_root, right_root)] = min(left_root, right_root)


def duplicate_key(row: dict[str, str]) -> tuple[str, ...]:
    return (
        row["gene_symbol"],
        row["structure_identity_key"],
        row["endpoint"],
        normalized_relation(row["relation"]),
        normalized_numeric(row["numeric_value_nm"]),
        row["publication_key"] or "NO_REFERENCE_IDENTIFIER",
    )


def assign_duplicate_groups(
    observations: list[dict[str, str]], cross_rows: list[dict[str, str]]
) -> tuple[list[dict[str, Any]], list[dict[str, str]], list[dict[str, str]]]:
    by_id = {row["observation_id"]: row for row in observations}
    within_memberships: dict[str, list[str]] = defaultdict(list)
    cross_memberships: dict[str, list[str]] = defaultdict(list)
    groups: list[dict[str, Any]] = []
    members: list[dict[str, str]] = []

    for source, prefix in (("ChEMBL", "P09C"), ("BindingDB", "P09B")):
        keyed: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
        for row in observations:
            if row["source"] == source:
                keyed[duplicate_key(row)].append(row)
        duplicate_items = sorted((key, rows) for key, rows in keyed.items() if len(rows) > 1)
        for index, (key, rows) in enumerate(duplicate_items, start=1):
            group_id = f"{prefix}{index:06d}"
            gene, ligand_key, endpoint, relation, numeric, ref_key = key
            groups.append(
                {
                    "duplicate_group_id": group_id,
                    "scope": f"WITHIN_{source.upper()}",
                    "gene_symbol": gene,
                    "endpoint": endpoint,
                    "match_basis": "EXACT_TARGET_LIGAND_ENDPOINT_RELATION_VALUE_REFERENCE_IDENTIFIER",
                    "member_observations_n": len(rows),
                    "distinct_source_records_n": len({row["source_record_id"] for row in rows}),
                    "sources": source,
                    "publication_key": "" if ref_key == "NO_REFERENCE_IDENTIFIER" else ref_key,
                    "ligand_identity_key": ligand_key,
                    "relation": relation,
                    "numeric_value_nm": numeric,
                    "review_status": "POTENTIAL_DUPLICATE_PRIMARY_REVIEW_REQUIRED",
                }
            )
            for row in sorted(rows, key=lambda item: item["observation_id"]):
                within_memberships[row["observation_id"]].append(group_id)
                members.append(duplicate_member(group_id, f"WITHIN_{source.upper()}", row))

    chembl_by_activity = {
        row["activity_id"]: row for row in observations if row["source"] == "ChEMBL"
    }
    binding_by_row_endpoint = {
        (row["bindingdb_source_row"], row["endpoint"]): row
        for row in observations
        if row["source"] == "BindingDB"
    }
    union = UnionFind()
    pair_records: list[dict[str, str]] = []
    for index, source_pair in enumerate(cross_rows, start=1):
        chembl = chembl_by_activity.get(source_pair["chembl_activity_id"])
        binding = binding_by_row_endpoint.get((source_pair["bindingdb_source_row"], source_pair["endpoint"]))
        if chembl is None or binding is None:
            raise ValueError("P0.8 cross-source duplicate link cannot be resolved")
        if chembl["gene_symbol"] != source_pair["gene_symbol"] or binding["gene_symbol"] != source_pair["gene_symbol"]:
            raise ValueError("P0.8 cross-source duplicate target mismatch")
        pair_signature = (
            normalized_doi(source_pair["doi"]),
            source_pair["endpoint"],
            normalized_relation(source_pair["relation"]),
            normalized_numeric(source_pair["numeric_value_nm"]),
            source_pair["ligand_smiles"],
        )
        chembl_signature = (
            normalized_doi(chembl["doi"]),
            chembl["endpoint"],
            normalized_relation(chembl["relation"]),
            normalized_numeric(chembl["numeric_value_nm"]),
            chembl["ligand_smiles"],
        )
        binding_signature = (
            normalized_doi(binding["doi"]),
            binding["endpoint"],
            normalized_relation(binding["relation"]),
            normalized_numeric(binding["numeric_value_nm"]),
            binding["ligand_smiles"],
        )
        if pair_signature != chembl_signature or pair_signature != binding_signature:
            raise ValueError("P0.8 cross-source duplicate exact-field signature mismatch")
        union.union(chembl["observation_id"], binding["observation_id"])
        pair_records.append(
            {
                "duplicate_pair_id": f"P09XP{index:06d}",
                "duplicate_group_id": "",
                "chembl_observation_id": chembl["observation_id"],
                "bindingdb_observation_id": binding["observation_id"],
                **source_pair,
            }
        )
    components: dict[str, list[str]] = defaultdict(list)
    for observation_id in union.parent:
        components[union.find(observation_id)].append(observation_id)
    ordered_components = sorted((sorted(ids) for ids in components.values()), key=lambda ids: ids[0])
    cross_group_by_observation: dict[str, str] = {}
    for index, observation_ids in enumerate(ordered_components, start=1):
        group_id = f"P09X{index:06d}"
        rows = [by_id[item] for item in observation_ids]
        for observation_id in observation_ids:
            cross_group_by_observation[observation_id] = group_id
            cross_memberships[observation_id].append(group_id)
            members.append(duplicate_member(group_id, "CROSS_SOURCE", by_id[observation_id]))
        groups.append(
            {
                "duplicate_group_id": group_id,
                "scope": "CROSS_SOURCE",
                "gene_symbol": rows[0]["gene_symbol"],
                "endpoint": rows[0]["endpoint"],
                "match_basis": "P08_EXACT_DOI_ENDPOINT_RELATION_VALUE_SMILES",
                "member_observations_n": len(rows),
                "distinct_source_records_n": len({(row["source"], row["source_record_id"]) for row in rows}),
                "sources": "|".join(sorted({row["source"] for row in rows})),
                "publication_key": rows[0]["publication_key"],
                "ligand_identity_key": f"SMILES:{rows[0]['ligand_smiles']}",
                "relation": normalized_relation(rows[0]["relation"]),
                "numeric_value_nm": normalized_numeric(rows[0]["numeric_value_nm"]),
                "review_status": "POTENTIAL_DUPLICATE_PRIMARY_REVIEW_REQUIRED",
            }
        )
    for pair in pair_records:
        pair["duplicate_group_id"] = cross_group_by_observation[pair["chembl_observation_id"]]

    for row in observations:
        within = sorted(within_memberships.get(row["observation_id"], []))
        cross = sorted(cross_memberships.get(row["observation_id"], []))
        row["within_source_duplicate_group_ids"] = "|".join(within)
        row["cross_source_duplicate_group_ids"] = "|".join(cross)
        if within and cross:
            row["potential_duplicate_status"] = "WITHIN_AND_CROSS_SOURCE_POTENTIAL_DUPLICATE"
        elif within:
            row["potential_duplicate_status"] = "WITHIN_SOURCE_POTENTIAL_DUPLICATE"
        elif cross:
            row["potential_duplicate_status"] = "CROSS_SOURCE_POTENTIAL_DUPLICATE"
        else:
            row["potential_duplicate_status"] = "NO_EXACT_FIELD_DUPLICATE_DETECTED"

    groups.sort(key=lambda row: row["duplicate_group_id"])
    members.sort(key=lambda row: (row["duplicate_group_id"], row["observation_id"]))
    return groups, members, pair_records


def duplicate_member(group_id: str, scope: str, row: dict[str, str]) -> dict[str, str]:
    return {
        "duplicate_group_id": group_id,
        "scope": scope,
        "observation_id": row["observation_id"],
        "source": row["source"],
        "gene_symbol": row["gene_symbol"],
        "endpoint": row["endpoint"],
        "source_molecule_id": row["source_molecule_id"],
        "source_measurement_id": row["source_measurement_id"],
        "source_experiment_id": row["source_experiment_id"],
        "source_record_id": row["source_record_id"],
        "publication_key": row["publication_key"],
    }


def classify_observation(row: dict[str, str]) -> None:
    warnings: list[str] = ["PRIMARY_EXPERIMENT_REVIEW_REQUIRED"]
    if row["source"] == "BindingDB":
        warnings.extend(["BINDINGDB_ASSAY_CONTEXT_UNVERIFIED", "CELLULAR_VS_BIOCHEMICAL_CONTEXT_NOT_VERIFIED"])
    if row["protein_match_status"] == "REVIEW_MULTI_COMPONENT_PROTEIN_OBJECT":
        warnings.append("MULTI_COMPONENT_PROTEIN_OBJECT")
    elif row["protein_match_status"] == "REVIEW_CHAIN_COUNT_UNKNOWN":
        warnings.append("CHAIN_COUNT_UNKNOWN")
    if row["species_status"] == "REVIEW_SPECIES_MISSING":
        warnings.append("SPECIES_MISSING")
    elif row["species_status"] == "REVIEW_SPECIES_MISMATCH":
        warnings.append("SPECIES_MISMATCH")
    if row["measurement_precision_status"] == "CENSORED":
        warnings.append("CENSORED_VALUE")
    elif row["measurement_precision_status"] not in {"EXACT_NUMERIC", "CENSORED"}:
        warnings.append("INVALID_OR_MISSING_VALUE")
    if row["endpoint"] == "IC50":
        warnings.append("IC50_DO_NOT_EQUATE_TO_KI_OR_KD")
    if row["structure_status"] != "STRUCTURE_PRESENT":
        warnings.append("STRUCTURE_MISSING")
    if row["publication_status"] != "REFERENCE_IDENTIFIER_PRESENT_PRIMARY_STATUS_UNVERIFIED":
        warnings.append("PRIMARY_REFERENCE_IDENTIFIER_MISSING")
    if row["curator_quality_status"] == "REVIEW_CURATOR_FLAGGED":
        warnings.append("CURATOR_FLAGGED")
    if row["within_source_duplicate_group_ids"]:
        warnings.append("WITHIN_SOURCE_POTENTIAL_DUPLICATE")
    if row["cross_source_duplicate_group_ids"]:
        warnings.append("CROSS_SOURCE_POTENTIAL_DUPLICATE")

    fundamental_invalid = (
        row["endpoint"] not in ENDPOINTS
        or row["source_units"] != "nM"
        or row["measurement_precision_status"] not in {"EXACT_NUMERIC", "CENSORED"}
        or not normalized_numeric(row["numeric_value_nm"])
        or row["structure_status"] != "STRUCTURE_PRESENT"
        or (
            row["source"] == "ChEMBL"
            and row["experiment_type_status"] != "BINDING_ASSAY_METADATA_MATCHED"
        )
    )
    target_context_ambiguous = (
        row["protein_match_status"] != "EXACT_SINGLE_PROTEIN_UNIPROT_MATCH"
        or row["species_status"] != "HUMAN_SPECIES_MATCH"
    )
    duplicate_or_curator_warning = (
        row["potential_duplicate_status"] != "NO_EXACT_FIELD_DUPLICATE_DETECTED"
        or row["curator_quality_status"] == "REVIEW_CURATOR_FLAGGED"
    )
    strict = (
        not fundamental_invalid
        and not target_context_ambiguous
        and row["source"] == "ChEMBL"
        and row["measurement_precision_status"] == "EXACT_NUMERIC"
        and row["publication_status"] == "REFERENCE_IDENTIFIER_PRESENT_PRIMARY_STATUS_UNVERIFIED"
        and not duplicate_or_curator_warning
    )

    if fundamental_invalid:
        row["catalog_inclusion_status"] = "EXCLUDED_INVALID"
        row["evidence_completeness_status"] = "NOT_CATALOG_ELIGIBLE"
    elif target_context_ambiguous:
        row["catalog_inclusion_status"] = "CONTEXT_ONLY"
        row["evidence_completeness_status"] = "CONTEXT_ONLY_PRIMARY_REVIEW_REQUIRED"
    elif strict:
        row["catalog_inclusion_status"] = "STRICT_CANDIDATE"
        row["evidence_completeness_status"] = "COMPLETE_METADATA_PRIMARY_REVIEW_PENDING"
    else:
        row["catalog_inclusion_status"] = "REVIEW_CANDIDATE"
        row["evidence_completeness_status"] = (
            "PARTIAL_ASSAY_CONTEXT_PRIMARY_REVIEW_REQUIRED"
            if row["source"] == "BindingDB"
            else "PARTIAL_METADATA_PRIMARY_REVIEW_REQUIRED"
        )
    row["warning_codes"] = "|".join(dict.fromkeys(warnings))


def build_ligand_summary(observations: list[dict[str, str]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in observations:
        grouped[
            (
                row["gene_symbol"],
                row["uniprot_accession"],
                row["source"],
                row["source_molecule_id"],
                row["structure_identity_key"],
                row["endpoint"],
            )
        ].append(row)
    priority = {
        "STRICT_CANDIDATE": 0,
        "REVIEW_CANDIDATE": 1,
        "CONTEXT_ONLY": 2,
        "EXCLUDED_INVALID": 3,
    }
    result: list[dict[str, Any]] = []
    for key, rows in sorted(grouped.items()):
        gene, accession, source, molecule, identity, endpoint = key
        ordered = sorted(rows, key=lambda row: row["observation_id"])
        counts = Counter(row["catalog_inclusion_status"] for row in ordered)
        measurements = [
            f"{row['observation_id']}={row['relation']}{row['raw_value_nm']} {row['source_units']}"
            for row in ordered
        ]
        result.append(
            {
                "gene_symbol": gene,
                "uniprot_accession": accession,
                "source": source,
                "source_molecule_id": molecule,
                "structure_identity_key": identity,
                "ligand_smiles": ordered[0]["ligand_smiles"],
                "ligand_inchikey": ordered[0]["ligand_inchikey"],
                "endpoint": endpoint,
                "observation_records_n": len(ordered),
                "exact_numeric_n": sum(row["measurement_precision_status"] == "EXACT_NUMERIC" for row in ordered),
                "censored_n": sum(row["measurement_precision_status"] == "CENSORED" for row in ordered),
                "invalid_or_missing_n": sum(row["measurement_precision_status"] not in {"EXACT_NUMERIC", "CENSORED"} for row in ordered),
                "strict_candidate_n": counts["STRICT_CANDIDATE"],
                "review_candidate_n": counts["REVIEW_CANDIDATE"],
                "context_only_n": counts["CONTEXT_ONLY"],
                "excluded_invalid_n": counts["EXCLUDED_INVALID"],
                "highest_eligibility_status": min(counts, key=lambda item: priority[item]),
                "reported_measurements": " | ".join(measurements),
                "publication_keys": "|".join(sorted({row["publication_key"] for row in ordered if row["publication_key"]})),
                "warning_codes": "|".join(sorted({code for row in ordered for code in row["warning_codes"].split("|") if code})),
                "observation_ids": "|".join(row["observation_id"] for row in ordered),
                "scientific_acceptance": "NOT_APPROVED",
            }
        )
    return result


def publication_metrics(rows: list[dict[str, str]]) -> dict[str, Any]:
    return {
        "sources": sorted({row["source"] for row in rows}),
        "endpoints": sorted({row["endpoint"] for row in rows}),
        "exact": sum(row["measurement_precision_status"] == "EXACT_NUMERIC" for row in rows),
        "censored": sum(row["measurement_precision_status"] == "CENSORED" for row in rows),
        "strict": sum(row["catalog_inclusion_status"] == "STRICT_CANDIDATE" for row in rows),
        "review": sum(row["catalog_inclusion_status"] == "REVIEW_CANDIDATE" for row in rows),
        "context": sum(row["catalog_inclusion_status"] == "CONTEXT_ONLY" for row in rows),
        "duplicate": sum(row["potential_duplicate_status"] != "NO_EXACT_FIELD_DUPLICATE_DETECTED" for row in rows),
        "multi": sum(row["protein_match_status"] != "EXACT_SINGLE_PROTEIN_UNIPROT_MATCH" for row in rows),
        "curator": sum(row["curator_quality_status"] == "REVIEW_CURATOR_FLAGGED" for row in rows),
    }


def select_priority_publications(
    observations: list[dict[str, str]], limit_per_target: int = 2
) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    by_target_publication: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in observations:
        if row["publication_key"]:
            by_target_publication[(row["gene_symbol"], row["publication_key"])].append(row)

    selected: list[tuple[str, str, str, list[dict[str, str]], dict[str, Any]]] = []
    for gene in FOCUS6:
        candidates = []
        for (row_gene, key), rows in by_target_publication.items():
            if row_gene != gene:
                continue
            metrics = publication_metrics(rows)
            dense_sort = (
                metrics["strict"],
                len(metrics["sources"]),
                len(metrics["endpoints"]),
                metrics["exact"],
                len(rows),
            )
            risk_sort = (
                metrics["multi"] if gene in {"ERBB3", "FNTB"} else 0,
                metrics["context"],
                metrics["duplicate"],
                metrics["censored"],
                metrics["curator"],
                len(rows),
            )
            candidates.append((key, rows, metrics, dense_sort, risk_sort))
        if not candidates:
            continue
        dense = max(candidates, key=lambda item: (item[3], item[0]))
        selected.append((gene, dense[0], "EVIDENCE_DENSE_METADATA", dense[1], dense[2]))
        if limit_per_target > 1:
            remaining = [item for item in candidates if item[0] != dense[0]]
            if remaining:
                risk = max(remaining, key=lambda item: (item[4], item[0]))
                selected.append((gene, risk[0], "QUALITY_RISK_REVIEW", risk[1], risk[2]))

    publication_rows: list[dict[str, Any]] = []
    link_rows: list[dict[str, str]] = []
    for rank, (gene, key, role, rows, metrics) in enumerate(selected, start=1):
        doi_values = sorted({row["doi"] for row in rows if row["doi"]})
        pmid_values = sorted({row["pmid"] for row in rows if row["pmid"]})
        publication_urls = sorted({row["publication_url"] for row in rows if row["publication_url"]})
        if role == "EVIDENCE_DENSE_METADATA":
            reason_code = "MAXIMIZE_EXACT_ENDPOINT_AND_SOURCE_COVERAGE"
            rationale = "Приоритет для проверки насыщенного блока точных измерений, конечных точек и источников."
        elif gene in {"ERBB3", "FNTB"} and metrics["multi"]:
            reason_code = "MULTI_COMPONENT_TARGET_COVERAGE_RISK"
            rationale = "Приоритет для проверки многокомпонентного белкового объекта и расхождения покрытия источников."
        elif metrics["duplicate"] or metrics["censored"] or metrics["curator"]:
            reason_code = "DUPLICATE_CENSORING_OR_CURATOR_RISK"
            rationale = "Приоритет для проверки потенциальных дублей, цензурирования или предупреждений куратора."
        else:
            reason_code = "ASSAY_CONTEXT_OR_METADATA_REVIEW"
            rationale = "Приоритет для проверки типа опыта и полноты экспериментальных метаданных."
        publication_rows.append(
            {
                "review_queue_rank": rank,
                "gene_symbol": gene,
                "publication_key": key,
                "publication_url": "|".join(publication_urls),
                "doi": "|".join(doi_values),
                "pmid": "|".join(pmid_values),
                "selection_role": role,
                "selection_reason_code": reason_code,
                "selection_rationale_ru": rationale,
                "sources": "|".join(metrics["sources"]),
                "endpoints": "|".join(metrics["endpoints"]),
                "observation_records_n": len(rows),
                "exact_numeric_records_n": metrics["exact"],
                "censored_records_n": metrics["censored"],
                "strict_candidate_records_n": metrics["strict"],
                "review_candidate_records_n": metrics["review"],
                "context_only_records_n": metrics["context"],
                "potential_duplicate_member_records_n": metrics["duplicate"],
                "multi_component_records_n": metrics["multi"],
                "primary_review_status": "REVIEW_REQUIRED",
                "scientific_acceptance": "NOT_APPROVED",
            }
        )
        for row in sorted(rows, key=lambda item: item["observation_id"]):
            link_rows.append(
                {
                    "review_queue_rank": rank,
                    "gene_symbol": gene,
                    "publication_key": key,
                    "observation_id": row["observation_id"],
                    "source": row["source"],
                    "endpoint": row["endpoint"],
                    "catalog_inclusion_status": row["catalog_inclusion_status"],
                    "warning_codes": row["warning_codes"],
                    "source_measurement_id": row["source_measurement_id"],
                    "source_experiment_id": row["source_experiment_id"],
                    "source_molecule_id": row["source_molecule_id"],
                }
            )
    return publication_rows, link_rows


def build_target_summary(
    observations: list[dict[str, str]],
    duplicate_groups: list[dict[str, Any]],
    duplicate_members: list[dict[str, str]],
    priority_publications: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for gene, accession in FOCUS6.items():
        rows = [row for row in observations if row["gene_symbol"] == gene]
        inclusion = Counter(row["catalog_inclusion_status"] for row in rows)
        within_groups = {row["duplicate_group_id"] for row in duplicate_groups if row["gene_symbol"] == gene and row["scope"].startswith("WITHIN_")}
        cross_groups = {row["duplicate_group_id"] for row in duplicate_groups if row["gene_symbol"] == gene and row["scope"] == "CROSS_SOURCE"}
        within_member_ids = {row["observation_id"] for row in duplicate_members if row["duplicate_group_id"] in within_groups}
        cross_member_ids = {row["observation_id"] for row in duplicate_members if row["duplicate_group_id"] in cross_groups}
        chembl_n = sum(row["source"] == "ChEMBL" for row in rows)
        binding_n = sum(row["source"] == "BindingDB" for row in rows)
        interpretation = (
            "COVERAGE_DISCORDANCE_MULTI_COMPONENT_REVIEW_REQUIRED"
            if gene in {"ERBB3", "FNTB"}
            else "SOURCE_COVERAGE_REVIEW_REQUIRED"
        )
        result.append(
            {
                "gene_symbol": gene,
                "uniprot_accession": accession,
                "chembl_observation_records_n": chembl_n,
                "bindingdb_observation_records_n": binding_n,
                "all_observation_records_n": len(rows),
                "ki_observation_records_n": sum(row["endpoint"] == "Ki" for row in rows),
                "kd_observation_records_n": sum(row["endpoint"] == "Kd" for row in rows),
                "ic50_observation_records_n": sum(row["endpoint"] == "IC50" for row in rows),
                "exact_numeric_observation_records_n": sum(row["measurement_precision_status"] == "EXACT_NUMERIC" for row in rows),
                "censored_observation_records_n": sum(row["measurement_precision_status"] == "CENSORED" for row in rows),
                "invalid_or_missing_observation_records_n": sum(row["measurement_precision_status"] not in {"EXACT_NUMERIC", "CENSORED"} for row in rows),
                "exact_single_protein_human_records_n": sum(row["protein_match_status"] == "EXACT_SINGLE_PROTEIN_UNIPROT_MATCH" and row["species_status"] == "HUMAN_SPECIES_MATCH" for row in rows),
                "multi_component_or_unknown_protein_records_n": sum(row["protein_match_status"] != "EXACT_SINGLE_PROTEIN_UNIPROT_MATCH" for row in rows),
                "species_missing_or_mismatch_records_n": sum(row["species_status"] != "HUMAN_SPECIES_MATCH" for row in rows),
                "structure_present_records_n": sum(row["structure_status"] == "STRUCTURE_PRESENT" for row in rows),
                "structure_missing_records_n": sum(row["structure_status"] != "STRUCTURE_PRESENT" for row in rows),
                "reference_identifier_present_records_n": sum(bool(row["publication_key"]) for row in rows),
                "reference_identifier_missing_records_n": sum(not row["publication_key"] for row in rows),
                "strict_candidate_records_n": inclusion["STRICT_CANDIDATE"],
                "review_candidate_records_n": inclusion["REVIEW_CANDIDATE"],
                "context_only_records_n": inclusion["CONTEXT_ONLY"],
                "excluded_invalid_records_n": inclusion["EXCLUDED_INVALID"],
                "unique_source_ligand_ids_n": len({(row["source"], row["source_molecule_id"]) for row in rows}),
                "strict_candidate_source_ligand_ids_n": len({(row["source"], row["source_molecule_id"]) for row in rows if row["catalog_inclusion_status"] == "STRICT_CANDIDATE"}),
                "within_source_duplicate_groups_n": len(within_groups),
                "within_source_duplicate_member_records_n": len(within_member_ids),
                "cross_source_duplicate_groups_n": len(cross_groups),
                "cross_source_duplicate_member_records_n": len(cross_member_ids),
                "priority_publications_n": sum(row["gene_symbol"] == gene for row in priority_publications),
                "coverage_interpretation": interpretation,
                "scientific_acceptance": "NOT_APPROVED",
            }
        )
    return result


def build_qc_rows(observations: list[dict[str, str]]) -> list[dict[str, Any]]:
    definitions = [
        ("REVIEW", "BINDINGDB_ASSAY_CONTEXT_UNVERIFIED", lambda row: row["source"] == "BindingDB", "BindingDB P0.8 не содержит подтверждённого типа первичного опыта.", "Проверить методику и биохимический/клеточный контекст по первоисточнику."),
        ("REVIEW", "MULTI_COMPONENT_PROTEIN_OBJECT", lambda row: row["protein_match_status"] != "EXACT_SINGLE_PROTEIN_UNIPROT_MATCH", "Наблюдение относится к многокомпонентному объекту или объекту с неизвестным числом цепей.", "Не переносить измерение на одиночный белок без ручной проверки состава объекта."),
        ("REVIEW", "SPECIES_UNCERTAIN", lambda row: row["species_status"] != "HUMAN_SPECIES_MATCH", "Вид отсутствует или не соответствует Homo sapiens.", "Проверить вид и исходный эксперимент."),
        ("REVIEW", "CENSORED_VALUE", lambda row: row["measurement_precision_status"] == "CENSORED", "Измерение содержит неравенство.", "Хранить отдельно от точных чисел; не использовать как точное значение."),
        ("HIGH", "INVALID_OR_MISSING_VALUE", lambda row: row["measurement_precision_status"] not in {"EXACT_NUMERIC", "CENSORED"}, "Численное измерение отсутствует или непригодно.", "Исключить из качественного каталога, сохранив исходную запись."),
        ("HIGH", "STRUCTURE_MISSING", lambda row: row["structure_status"] != "STRUCTURE_PRESENT", "Структура лиганда отсутствует.", "Исключить из строгого каталога до восстановления химической идентичности."),
        ("REVIEW", "REFERENCE_IDENTIFIER_MISSING", lambda row: not row["publication_key"], "DOI и PMID отсутствуют.", "Найти первичную публикацию или оставить запись требующей проверки."),
        ("REVIEW", "POTENTIAL_DUPLICATE", lambda row: row["potential_duplicate_status"] != "NO_EXACT_FIELD_DUPLICATE_DETECTED", "Запись входит в потенциальную группу дублей.", "Проверить первичный эксперимент; исходные строки не удалять."),
    ]
    result: list[dict[str, Any]] = []
    for severity, category, predicate, finding, action in definitions:
        for gene in FOCUS6:
            affected = sum(row["gene_symbol"] == gene and predicate(row) for row in observations)
            if affected:
                result.append(
                    {
                        "severity": severity,
                        "category": category,
                        "gene_symbol": gene,
                        "affected_records_n": affected,
                        "finding_ru": finding,
                        "required_action_ru": action,
                    }
                )
    return result


def output_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, delimiter="\t", fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def report_ru(
    summary: dict[str, Any], target_rows: list[dict[str, Any]], priority_rows: list[dict[str, Any]]
) -> str:
    totals = summary["quality_classification"]
    duplicate = summary["duplicates"]
    lines = [
        "# MCL — TLS V3 / P0.9. Каталог экспериментальных взаимодействий известных молекул",
        "",
        f"**Контракт:** `{CONTRACT}`  ",
        "**Научная приёмка:** `NOT_APPROVED`  ",
        "**Слияние:** `merge_permitted=false`",
        "",
        "## Резюме",
        "",
        (
            f"Проверенный вход P0.8 содержит {summary['observation_records_n']} построчных наблюдений: "
            f"{summary['source_observation_records']['ChEMBL']} ChEMBL и "
            f"{summary['source_observation_records']['BindingDB']} BindingDB. Это число записей, а не число независимых экспериментов."
        ),
        "",
        (
            f"К строгим кандидатам отнесено {totals['STRICT_CANDIDATE']} наблюдений. "
            f"Ещё {totals['REVIEW_CANDIDATE']} требуют проверки, {totals['CONTEXT_ONLY']} сохранены только как контекст, "
            f"{totals['EXCLUDED_INVALID']} исключены из качественного каталога без удаления исходных строк. "
            "Даже строгие кандидаты не считаются научно утверждёнными до проверки первичной публикации."
        ),
        "",
        "## Проверка входа P0.8",
        "",
        f"- ZIP SHA-256: `{summary['p08_input']['sha256']}`.",
        f"- CRC: `{summary['p08_input']['crc_status']}`; несоответствий размеров/SHA-256 внутренних таблиц: 0.",
        "- Релизы: ChEMBL 37 и BindingDB 202611; усечённых наблюдений: 0.",
        "- Исходные научные базы не открывались и не обрабатывались повторно.",
        "",
        "## Критерии включения",
        "",
        "`STRICT_CANDIDATE` требует точного одиночного белка, подтверждённого Homo sapiens, положительного точного Ki/Kd/IC50 в нМ, структуры, DOI или PMID, ChEMBL binding-assay с confidence score 9 и отсутствия признаков дубля/предупреждений куратора. BindingDB остаётся `REVIEW_CANDIDATE`, поскольку тип и условия первичного опыта не подтверждены в P0.8.",
        "",
        "Цензурированные значения, Ki, Kd и IC50 хранятся раздельно. IC50 не приравнивается к Ki или Kd. Наблюдения многокомпонентных объектов и неопределённого вида не повышаются до доказательства взаимодействия с одиночным человеческим белком.",
        "",
        "## Результаты по мишеням",
        "",
        "| Ген | ChEMBL, записи | BindingDB, записи | Ki | Kd | IC50 | Строгие | Требуют проверки | Контекст | Исключены |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in target_rows:
        lines.append(
            f"| {row['gene_symbol']} | {row['chembl_observation_records_n']} | {row['bindingdb_observation_records_n']} | "
            f"{row['ki_observation_records_n']} | {row['kd_observation_records_n']} | {row['ic50_observation_records_n']} | "
            f"{row['strict_candidate_records_n']} | {row['review_candidate_records_n']} | "
            f"{row['context_only_records_n']} | {row['excluded_invalid_records_n']} |"
        )
    erbb3 = next(row for row in target_rows if row["gene_symbol"] == "ERBB3")
    fntb = next(row for row in target_rows if row["gene_symbol"] == "FNTB")
    lines.extend(
        [
            "",
            "## ERBB3 и FNTB",
            "",
            (
                f"ERBB3: ChEMBL — {erbb3['chembl_observation_records_n']} записей, BindingDB — {erbb3['bindingdb_observation_records_n']}; "
                f"многокомпонентных/неопределённых белковых объектов — {erbb3['multi_component_or_unknown_protein_records_n']}. "
                "Расхождение покрытия нельзя интерпретировать как отсутствие или наличие лигандов без проверки состава объекта."
            ),
            "",
            (
                f"FNTB: ChEMBL — {fntb['chembl_observation_records_n']} записи, BindingDB — {fntb['bindingdb_observation_records_n']}; "
                f"многокомпонентных/неопределённых белковых объектов — {fntb['multi_component_or_unknown_protein_records_n']}. "
                "Для фарнезилтрансферазы особенно важно отличать измерение комплекса от свидетельства для одной субъединицы FNTB."
            ),
            "",
            "## Потенциальные дубли",
            "",
            (
                f"Выявлено {duplicate['within_chembl_groups_n']} групп внутри ChEMBL, "
                f"{duplicate['within_bindingdb_groups_n']} групп внутри BindingDB и "
                f"{duplicate['cross_source_groups_n']} межбазовых групп ({duplicate['cross_source_pairs_n']} исходных пар P0.8). "
                "Дубликаты не удалялись: таблица членов групп сохраняет каждое наблюдение и все исходные идентификаторы."
            ),
            "",
            "## Публикации для ручной проверки",
            "",
            "Очередь формируется детерминированно: до двух публикаций на белок — одна с наиболее насыщенным набором точных конечных точек/источников и одна с максимальным риском дублей, цензурирования или многокомпонентности. Это приоритет ручной работы, а не оценка научной значимости.",
            "",
            "| № | Ген | Публикация | Роль | Источники | Записи | Причина |",
            "|---:|---|---|---|---|---:|---|",
        ]
    )
    for row in priority_rows:
        lines.append(
            f"| {row['review_queue_rank']} | {row['gene_symbol']} | {row['publication_key']} | "
            f"{row['selection_role']} | {row['sources']} | {row['observation_records_n']} | {row['selection_rationale_ru']} |"
        )
    lines.extend(
        [
            "",
            "## Ограничения до P1",
            "",
            "- Наличие DOI/PMID не доказывает, что запись правильно передаёт методику первичного эксперимента.",
            "- BindingDB P0.8 не позволяет автоматически отделить биохимический binding assay от клеточного контекста; все его записи требуют ручной проверки методики.",
            "- Точное совпадение структуры/значения/публикации обозначает потенциальный дубль, но не доказывает идентичность экспериментов.",
            "- Источник с большим числом записей не обязательно имеет более полное или независимое экспериментальное покрытие.",
            "- Результаты относятся к известным лигандам и не являются измерениями 65 PYZ. Общий показатель противоопухолевой эффективности не рассчитывался.",
            "",
            "`scientific_acceptance=NOT_APPROVED`  ",
            "`merge_permitted=false`",
            "",
        ]
    )
    return "\n".join(lines)


def build_catalog(
    *,
    p08_zip: Path,
    expected_p08_sha256: str,
    output_dir: Path,
    output_zip: Path,
    git_commit: str = "UNKNOWN",
    created_utc: str | None = None,
) -> dict[str, Any]:
    if output_dir.resolve() == p08_zip.resolve() or output_zip.resolve() == p08_zip.resolve():
        raise ValueError("P0.9 output cannot overwrite the P0.8 input")
    if output_dir.exists() or output_zip.exists():
        raise FileExistsError("P0.9 output directory or ZIP already exists")
    created_utc = created_utc or datetime.now(UTC).isoformat().replace("+00:00", "Z")
    input_before = sha_file(p08_zip)
    verified = verify_p08_report(p08_zip, expected_p08_sha256)

    observations = build_observations(verified["chembl_rows"], verified["binding_rows"])
    duplicate_groups, duplicate_members, cross_pairs = assign_duplicate_groups(observations, verified["cross_rows"])
    for row in observations:
        classify_observation(row)
    observations.sort(key=lambda row: (list(FOCUS6).index(row["gene_symbol"]), row["source"], row["observation_id"]))

    ligand_summary = build_ligand_summary(observations)
    priority_publications, priority_links = select_priority_publications(observations)
    target_summary = build_target_summary(
        observations, duplicate_groups, duplicate_members, priority_publications
    )
    qc_rows = build_qc_rows(observations)

    inclusion = Counter(row["catalog_inclusion_status"] for row in observations)
    source_counts = Counter(row["source"] for row in observations)
    endpoint_counts = Counter(row["endpoint"] for row in observations)
    measurement_counts = Counter(row["measurement_precision_status"] for row in observations)
    duplicate_scope_counts = Counter(row["scope"] for row in duplicate_groups)

    summary = {
        "contract": CONTRACT,
        "created_utc": created_utc,
        "git_commit": git_commit,
        "technical_status": "REVIEW_REQUIRED",
        "scientific_acceptance": "NOT_APPROVED",
        "merge_permitted": False,
        "external_requests_n": 0,
        "input_scope": "VALIDATED_P08_REPORT_ONLY",
        "p08_input": {
            "path": str(p08_zip.resolve()),
            "sha256": verified["actual_zip_sha256"],
            "size_bytes": verified["zip_size_bytes"],
            "crc_status": verified["zip_crc_status"],
            "manifest_hash_or_size_mismatches_n": verified["manifest_hash_or_size_mismatches_n"],
            "source_releases": verified["package_manifest"]["source_releases"],
            "truncated_observations_n": verified["validation"]["truncated_observations_n"],
        },
        "observation_records_n": len(observations),
        "source_observation_records": dict(source_counts),
        "endpoint_observation_records": {endpoint: endpoint_counts[endpoint] for endpoint in ENDPOINTS},
        "measurement_precision": dict(measurement_counts),
        "quality_classification": {
            status: inclusion[status]
            for status in ("STRICT_CANDIDATE", "REVIEW_CANDIDATE", "CONTEXT_ONLY", "EXCLUDED_INVALID")
        },
        "duplicates": {
            "within_chembl_groups_n": duplicate_scope_counts["WITHIN_CHEMBL"],
            "within_bindingdb_groups_n": duplicate_scope_counts["WITHIN_BINDINGDB"],
            "cross_source_groups_n": duplicate_scope_counts["CROSS_SOURCE"],
            "cross_source_pairs_n": len(cross_pairs),
            "source_records_deleted_n": 0,
        },
        "ligand_endpoint_summary_rows_n": len(ligand_summary),
        "priority_publications_n": len(priority_publications),
        "priority_publication_links_n": len(priority_links),
        "qc_findings_n": len(qc_rows),
        "scientific_verified_binding_to_pyz_n": 0,
        "anticancer_efficacy_score_created": False,
        "baseline_65x13_untouched": True,
        "limitations": [
            "All BindingDB experiment contexts require primary-method review",
            "Potential duplicate groups are exact-field review candidates, not proven duplicate experiments",
            "Ki, Kd, and IC50 remain separate and are not pooled",
            "Known-ligand measurements are not transferred to PYZ",
        ],
    }

    output_dir.mkdir(parents=True, exist_ok=False)
    output_tsv(output_dir / "p09_observation_assessment.tsv", observations, ASSESSMENT_FIELDS)
    output_tsv(output_dir / "p09_known_ligand_catalog.tsv", observations, CATALOG_FIELDS)
    output_tsv(output_dir / "p09_ligand_endpoint_summary.tsv", ligand_summary, LIGAND_SUMMARY_FIELDS)
    output_tsv(output_dir / "p09_duplicate_groups.tsv", duplicate_groups, DUPLICATE_GROUP_FIELDS)
    output_tsv(output_dir / "p09_duplicate_group_members.tsv", duplicate_members, DUPLICATE_MEMBER_FIELDS)
    output_tsv(output_dir / "p09_cross_source_duplicate_pairs.tsv", cross_pairs, CROSS_PAIR_FIELDS)
    output_tsv(output_dir / "p09_priority_publications.tsv", priority_publications, PRIORITY_PUBLICATION_FIELDS)
    output_tsv(output_dir / "p09_priority_publication_observation_links.tsv", priority_links, PRIORITY_LINK_FIELDS)
    output_tsv(output_dir / "p09_target_summary.tsv", target_summary, TARGET_SUMMARY_FIELDS)
    output_tsv(output_dir / "p09_QC.tsv", qc_rows, QC_FIELDS)
    write_json(output_dir / "p09_summary.json", summary)

    input_verification = {
        "contract": CONTRACT,
        "p08_validated_package_contract": P08_PACKAGE_CONTRACT,
        "p08_report_contract": P08_CONTRACT,
        "p08_zip_sha256_expected": require_sha256(expected_p08_sha256, "Expected P0.8 ZIP SHA-256"),
        "p08_zip_sha256_actual": verified["actual_zip_sha256"],
        "p08_zip_crc_status": verified["zip_crc_status"],
        "p08_zip_duplicate_members_n": verified["duplicate_members_n"],
        "p08_zip_unsafe_members_n": verified["unsafe_members_n"],
        "p08_manifest_files_n": verified["manifest_member_files_n"],
        "p08_manifest_hash_or_size_mismatches_n": 0,
        "p08_internal_files": verified["file_verification"],
        "p08_validation_status": verified["package_manifest"]["validation_status"],
        "p08_source_input_hashes_unchanged": verified["package_manifest"]["source_input_hashes_unchanged"],
        "p08_truncated_observations_n": verified["validation"]["truncated_observations_n"],
        "scientific_acceptance": "NOT_APPROVED",
        "merge_permitted": False,
    }
    write_json(output_dir / "p09_input_verification.json", input_verification)

    input_after = sha_file(p08_zip)
    if input_after != input_before:
        raise RuntimeError("P0.8 input ZIP changed during P0.9")
    source_manifest = {
        "contract": CONTRACT,
        "created_utc": created_utc,
        "input_p08_zip": {
            "path": str(p08_zip.resolve()),
            "size_bytes": p08_zip.stat().st_size,
            "sha256_before": input_before,
            "sha256_after": input_after,
        },
        "upstream_source_releases": verified["package_manifest"]["source_releases"],
        "upstream_source_metadata": verified["summary"]["source_metadata"],
        "upstream_source_input_hashes": verified["source_manifest"]["input_paths_sha256"],
        "no_network": True,
        "read_only_p08_input": True,
        "multi_gigabyte_sources_opened": False,
        "baseline_65x13_untouched": True,
        "scientific_acceptance": "NOT_APPROVED",
        "merge_permitted": False,
    }
    write_json(output_dir / "p09_source_manifest.json", source_manifest)
    (output_dir / "P09_REPORT_RU.md").write_text(
        report_ru(summary, target_summary, priority_publications), encoding="utf-8"
    )

    package_files: dict[str, dict[str, Any]] = {}
    for path in sorted(item for item in output_dir.iterdir() if item.is_file()):
        package_files[path.name] = {"size_bytes": path.stat().st_size, "sha256": sha_file(path)}
    package_manifest = {
        "contract": CONTRACT,
        "created_utc": created_utc,
        "input_p08_zip_sha256": input_before,
        "scientific_acceptance": "NOT_APPROVED",
        "merge_permitted": False,
        "files": package_files,
    }
    write_json(output_dir / "p09_package_manifest.json", package_manifest)

    output_zip.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output_zip, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in sorted(item for item in output_dir.iterdir() if item.is_file()):
            info = zipfile.ZipInfo(path.name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, path.read_bytes(), compresslevel=9)
    with zipfile.ZipFile(output_zip) as output_archive:
        if output_archive.testzip() is not None:
            raise RuntimeError("P0.9 output ZIP CRC validation failed")
    summary["output_zip"] = {
        "path": str(output_zip.resolve()),
        "size_bytes": output_zip.stat().st_size,
        "sha256": sha_file(output_zip),
        "crc_status": "PASSED",
    }
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p08-zip", type=Path, required=True)
    parser.add_argument("--expected-p08-sha256", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--output-zip", type=Path, required=True)
    parser.add_argument("--git-commit", default="UNKNOWN")
    parser.add_argument("--created-utc")
    args = parser.parse_args()
    result = build_catalog(
        p08_zip=args.p08_zip,
        expected_p08_sha256=args.expected_p08_sha256,
        output_dir=args.output_dir,
        output_zip=args.output_zip,
        git_commit=args.git_commit,
        created_utc=args.created_utc,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
