#!/usr/bin/env python3
"""MCL TLS V3/P0.10: primary-publication review of selected P0.9 evidence.

The program reads only the validated P0.9 ZIP and a version-controlled curation
journal. It never opens ChEMBL or BindingDB source databases, never infers PYZ
activity and never upgrades scientific acceptance automatically.
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
from typing import Any, Self

CONTRACT = "mcl-tls-v3-p010-primary-publication-review-v1"
CURATION_CONTRACT = "mcl-tls-v3-p010-publication-curation-input-v1"
P09_CONTRACT = "mcl-tls-v3-p09-known-ligand-catalog-v1"
SHA256_RE = re.compile(r"[0-9a-f]{64}")

FOCUS6 = {
    "PRKDC": "P78527",
    "BCL2L1": "Q07817",
    "IGF1R": "P08069",
    "ERBB3": "P21860",
    "FNTB": "P49356",
    "DDC": "P20711",
}

REVIEW_STATUSES = {
    "CONFIRMED_PRIMARY_SOURCE",
    "REQUIRES_REVIEW",
    "EVIDENCE_TYPE_MISMATCH",
}
ACCESS_STATUSES = {
    "FULL_TEXT_PRIMARY_OPEN",
    "FULL_TEXT_AUTHOR_MANUSCRIPT",
    "PUBLISHER_METADATA_AND_ABSTRACT_ONLY",
    "NONAUTHORITATIVE_COPY_NOT_USED_FOR_CONFIRMATION",
}
PRIMARY_ACCESS = {"FULL_TEXT_PRIMARY_OPEN", "FULL_TEXT_AUTHOR_MANUSCRIPT"}

P09_REQUIRED_MEMBERS = {
    "p09_package_manifest.json",
    "p09_summary.json",
    "p09_source_manifest.json",
    "p09_observation_assessment.tsv",
    "p09_priority_publications.tsv",
    "p09_priority_publication_observation_links.tsv",
    "p09_duplicate_groups.tsv",
    "p09_duplicate_group_members.tsv",
    "p09_cross_source_duplicate_pairs.tsv",
}

PUBLICATION_FIELDS = [
    "publication_id",
    "review_queue_rank",
    "gene_symbol",
    "uniprot_accession",
    "doi",
    "pmid",
    "title",
    "publisher_url",
    "pubmed_url",
    "full_text_url",
    "access_status",
    "access_limitation_ru",
    "protein_object_ru",
    "species_ru",
    "construct_ru",
    "assay_class",
    "endpoint_interpretation_ru",
    "method_ru",
    "conditions_ru",
    "controls_ru",
    "primary_source_locator",
    "publication_review_status",
    "conclusion_ru",
    "warning_codes",
    "p09_observation_records_n",
    "p09_exact_numeric_records_n",
    "p09_censored_records_n",
    "p09_sources",
    "p09_endpoints",
    "scientific_acceptance",
]

INTERACTION_FIELDS = [
    "interaction_id",
    "publication_id",
    "gene_symbol",
    "uniprot_accession",
    "ligand_label",
    "protein_object_ru",
    "organism_ru",
    "construct_ru",
    "assay_class",
    "endpoint",
    "relation",
    "numeric_value_nm",
    "reported_value_text",
    "evidence_scope_ru",
    "primary_source_locator",
    "review_status",
    "rationale_ru",
    "warning_codes",
    "representative_observation_id",
    "linked_observations_n",
    "linked_sources",
    "source_molecule_ids",
    "source_protein_ids",
    "source_experiment_ids",
    "source_measurement_ids",
    "source_record_ids",
    "doi",
    "pmid",
    "scientific_acceptance",
]

LINK_FIELDS = [
    "interaction_id",
    "publication_id",
    "review_status",
    "observation_id",
    "source",
    "gene_symbol",
    "uniprot_accession",
    "source_target_id",
    "source_molecule_id",
    "source_experiment_id",
    "source_measurement_id",
    "source_record_id",
    "source_document_id",
    "endpoint",
    "relation",
    "numeric_value_nm",
    "raw_value_nm",
    "doi",
    "pmid",
    "protein_match_status",
    "species_status",
    "experiment_type_status",
    "measurement_precision_status",
    "scientific_acceptance",
]

DUPLICATE_FIELDS = [
    "interaction_id",
    "publication_id",
    "gene_symbol",
    "curator_duplicate_assessment",
    "curator_rationale_ru",
    "within_chembl_group_ids",
    "within_bindingdb_group_ids",
    "cross_source_group_ids",
    "cross_source_pair_ids",
    "linked_observations_n",
    "records_removed_n",
    "review_status",
]

UNCERTAINTY_FIELDS = [
    "uncertainty_id",
    "publication_id",
    "interaction_id",
    "gene_symbol",
    "warning_code",
    "description_ru",
    "resolution_required_ru",
    "review_status",
]

QC_FIELDS = ["check_id", "status", "value", "expected", "details_ru"]


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def require_sha256(value: str, label: str) -> str:
    normalized = value.strip().lower()
    if not SHA256_RE.fullmatch(normalized):
        raise ValueError(f"{label} must be a lowercase SHA-256")
    return normalized


def text(value: Any) -> str:
    return "" if value is None else str(value).strip()


def normalized_doi(value: Any) -> str:
    result = text(value).lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        result = result.removeprefix(prefix)
    return result.strip()


def normalized_numeric(value: Any) -> str:
    raw = text(value)
    if not raw:
        return ""
    try:
        number = Decimal(raw)
    except InvalidOperation as error:
        raise ValueError(f"invalid numeric value: {raw}") from error
    normalized = format(number.normalize(), "f")
    return "0" if normalized in {"-0", ""} else normalized


def split_pipe(value: Any) -> set[str]:
    return {part for part in text(value).split("|") if part}


def join_values(values: Iterable[Any]) -> str:
    return "|".join(sorted({text(value) for value in values if text(value)}))


def read_json(data: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(data.decode("utf-8-sig"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid JSON in {label}: {error}") from error
    if not isinstance(value, dict):
        raise TypeError(f"{label} must contain a JSON object")
    return value


def read_tsv(data: bytes, label: str) -> list[dict[str, str]]:
    try:
        content = data.decode("utf-8-sig")
    except UnicodeDecodeError as error:
        raise ValueError(f"invalid UTF-8 in {label}") from error
    reader = csv.DictReader(io.StringIO(content), delimiter="\t")
    if not reader.fieldnames:
        raise ValueError(f"missing TSV header in {label}")
    rows: list[dict[str, str]] = []
    for index, row in enumerate(reader, start=2):
        if None in row:
            raise ValueError(f"extra fields in {label} row {index}")
        rows.append({key: text(value) for key, value in row.items()})
    return rows


class ValidatedZip:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.archive = zipfile.ZipFile(path, "r")
        infos = self.archive.infolist()
        names = [info.filename for info in infos]
        duplicates = sorted(name for name, count in Counter(names).items() if count > 1)
        if duplicates:
            self.archive.close()
            raise ValueError(f"duplicate ZIP members: {duplicates}")
        unsafe: list[str] = []
        for info in infos:
            name = info.filename
            pure = PurePosixPath(name.replace("\\", "/"))
            if (
                not name
                or name.startswith(("/", "\\"))
                or re.match(r"^[A-Za-z]:", name)
                or ".." in pure.parts
                or info.is_dir()
            ):
                unsafe.append(name)
        if unsafe:
            self.archive.close()
            raise ValueError(f"unsafe ZIP members: {unsafe}")
        bad_crc = self.archive.testzip()
        if bad_crc is not None:
            self.archive.close()
            raise ValueError(f"ZIP CRC failure: {bad_crc}")
        self.names = set(names)

    def read(self, name: str) -> bytes:
        if name not in self.names:
            raise ValueError(f"missing ZIP member: {name}")
        return self.archive.read(name)

    def close(self) -> None:
        self.archive.close()

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *_: object) -> None:
        self.close()


def verify_p09_report(path: Path, expected_sha256: str) -> dict[str, Any]:
    expected = require_sha256(expected_sha256, "expected P0.9 SHA-256")
    actual = sha_file(path)
    if actual != expected:
        raise ValueError(f"P0.9 ZIP SHA-256 mismatch: expected {expected}, actual {actual}")

    with ValidatedZip(path) as archive:
        missing = sorted(P09_REQUIRED_MEMBERS - archive.names)
        if missing:
            raise ValueError(f"P0.9 ZIP missing required members: {missing}")
        manifest = read_json(archive.read("p09_package_manifest.json"), "p09_package_manifest.json")
        if manifest.get("contract") != P09_CONTRACT:
            raise ValueError("P0.9 package contract mismatch")
        if (
            manifest.get("scientific_acceptance") != "NOT_APPROVED"
            or manifest.get("merge_permitted") is not False
        ):
            raise ValueError("P0.9 acceptance safeguards mismatch")
        files = manifest.get("files")
        if not isinstance(files, dict) or not files:
            raise ValueError("P0.9 package manifest has no file entries")
        expected_names = set(files) | {"p09_package_manifest.json"}
        if archive.names != expected_names:
            unexpected = sorted(archive.names - expected_names)
            unlisted_missing = sorted(expected_names - archive.names)
            raise ValueError(
                f"P0.9 ZIP/manifest member set mismatch: unexpected={unexpected}, missing={unlisted_missing}"
            )
        for name, metadata in files.items():
            if name not in archive.names:
                raise ValueError(f"manifested P0.9 member missing: {name}")
            if not isinstance(metadata, dict):
                raise TypeError(f"invalid manifest metadata for {name}")
            payload = archive.read(name)
            if len(payload) != metadata.get("size_bytes"):
                raise ValueError(f"P0.9 member size mismatch: {name}")
            if sha_bytes(payload) != metadata.get("sha256"):
                raise ValueError(f"P0.9 member SHA-256 mismatch: {name}")

        summary = read_json(archive.read("p09_summary.json"), "p09_summary.json")
        source_manifest = read_json(
            archive.read("p09_source_manifest.json"), "p09_source_manifest.json"
        )
        if (
            summary.get("contract") != P09_CONTRACT
            or source_manifest.get("contract") != P09_CONTRACT
        ):
            raise ValueError("P0.9 internal contract mismatch")
        for value, label in ((summary, "summary"), (source_manifest, "source manifest")):
            if (
                value.get("scientific_acceptance") != "NOT_APPROVED"
                or value.get("merge_permitted") is not False
            ):
                raise ValueError(f"P0.9 {label} acceptance safeguards mismatch")
        if source_manifest.get("multi_gigabyte_sources_opened") is not False:
            raise ValueError(
                "P0.9 source manifest does not prove the intended read-only derived-data workflow"
            )

        result = {
            "zip_sha256": actual,
            "zip_size_bytes": path.stat().st_size,
            "manifest": manifest,
            "summary": summary,
            "source_manifest": source_manifest,
            "observations": read_tsv(
                archive.read("p09_observation_assessment.tsv"), "p09_observation_assessment.tsv"
            ),
            "priority_publications": read_tsv(
                archive.read("p09_priority_publications.tsv"), "p09_priority_publications.tsv"
            ),
            "priority_links": read_tsv(
                archive.read("p09_priority_publication_observation_links.tsv"),
                "p09_priority_publication_observation_links.tsv",
            ),
            "duplicate_groups": read_tsv(
                archive.read("p09_duplicate_groups.tsv"), "p09_duplicate_groups.tsv"
            ),
            "duplicate_members": read_tsv(
                archive.read("p09_duplicate_group_members.tsv"), "p09_duplicate_group_members.tsv"
            ),
            "cross_pairs": read_tsv(
                archive.read("p09_cross_source_duplicate_pairs.tsv"),
                "p09_cross_source_duplicate_pairs.tsv",
            ),
            "member_names_n": len(archive.names),
        }
    return result


def load_curation(path: Path) -> dict[str, Any]:
    return read_json(path.read_bytes(), str(path))


def _require_fields(record: dict[str, Any], fields: Iterable[str], label: str) -> None:
    missing = [field for field in fields if not text(record.get(field))]
    if missing:
        raise ValueError(f"{label} missing required fields: {missing}")


def validate_curation(curation: dict[str, Any], p09: dict[str, Any]) -> dict[str, Any]:
    if curation.get("contract") != CURATION_CONTRACT:
        raise ValueError("P0.10 curation contract mismatch")
    if (
        curation.get("scientific_acceptance") != "NOT_APPROVED"
        or curation.get("merge_permitted") is not False
    ):
        raise ValueError("P0.10 curation acceptance safeguards mismatch")
    publications = curation.get("publications")
    interactions = curation.get("interactions")
    uncertainties = curation.get("uncertainties")
    sources = curation.get("source_access_log")
    if (
        not isinstance(publications, list)
        or not isinstance(interactions, list)
        or not isinstance(uncertainties, list)
    ):
        raise TypeError("P0.10 curation lists are missing")
    if not isinstance(sources, list) or not sources:
        raise ValueError("P0.10 source access log is missing")
    if len(publications) != 12:
        raise ValueError(
            f"P0.10 requires exactly 12 publication reviews, found {len(publications)}"
        )
    if len(interactions) != 12:
        raise ValueError(
            f"P0.10 requires exactly 12 representative interactions, found {len(interactions)}"
        )

    p09_priority = {
        (row["gene_symbol"], normalized_doi(row["doi"]), row["pmid"]): row
        for row in p09["priority_publications"]
    }
    publication_by_id: dict[str, dict[str, Any]] = {}
    reviewed_keys: set[tuple[str, str, str]] = set()
    for record in publications:
        if not isinstance(record, dict):
            raise TypeError("publication review must be an object")
        _require_fields(
            record,
            (
                "publication_id",
                "gene_symbol",
                "doi",
                "pmid",
                "title",
                "publisher_url",
                "pubmed_url",
                "access_status",
                "protein_object_ru",
                "species_ru",
                "construct_ru",
                "assay_class",
                "endpoint_interpretation_ru",
                "method_ru",
                "conditions_ru",
                "controls_ru",
                "primary_source_locator",
                "publication_review_status",
                "conclusion_ru",
            ),
            "publication review",
        )
        publication_id = text(record["publication_id"])
        if publication_id in publication_by_id:
            raise ValueError(f"duplicate publication_id: {publication_id}")
        key = (text(record["gene_symbol"]), normalized_doi(record["doi"]), text(record["pmid"]))
        if key not in p09_priority:
            raise ValueError(f"curation publication is not in P0.9 priority queue: {key}")
        if key in reviewed_keys:
            raise ValueError(f"duplicate publication review: {key}")
        reviewed_keys.add(key)
        status = text(record["publication_review_status"])
        access = text(record["access_status"])
        if status not in REVIEW_STATUSES:
            raise ValueError(f"invalid publication review status: {status}")
        if access not in ACCESS_STATUSES:
            raise ValueError(f"invalid access status: {access}")
        if status == "CONFIRMED_PRIMARY_SOURCE" and access not in PRIMARY_ACCESS:
            raise ValueError(f"full primary text is required for confirmation: {publication_id}")
        if access not in PRIMARY_ACCESS and not text(record.get("access_limitation_ru")):
            raise ValueError(f"access limitation must be explicit: {publication_id}")
        publication_by_id[publication_id] = record
    if reviewed_keys != set(p09_priority):
        raise ValueError("P0.10 reviews do not exactly cover the 12 P0.9 priority publications")

    observations = {row["observation_id"]: row for row in p09["observations"]}
    interaction_by_id: dict[str, dict[str, Any]] = {}
    for record in interactions:
        if not isinstance(record, dict):
            raise TypeError("interaction review must be an object")
        _require_fields(
            record,
            (
                "interaction_id",
                "publication_id",
                "gene_symbol",
                "uniprot_accession",
                "ligand_label",
                "protein_object_ru",
                "organism_ru",
                "construct_ru",
                "assay_class",
                "endpoint",
                "relation",
                "numeric_value_nm",
                "reported_value_text",
                "evidence_scope_ru",
                "primary_source_locator",
                "review_status",
                "rationale_ru",
                "representative_observation_id",
                "duplicate_assessment",
                "duplicate_rationale_ru",
            ),
            "interaction review",
        )
        interaction_id = text(record["interaction_id"])
        if interaction_id in interaction_by_id:
            raise ValueError(f"duplicate interaction_id: {interaction_id}")
        publication_id = text(record["publication_id"])
        publication = publication_by_id.get(publication_id)
        if publication is None:
            raise ValueError(f"unknown publication_id in interaction: {publication_id}")
        status = text(record["review_status"])
        if status not in REVIEW_STATUSES:
            raise ValueError(f"invalid interaction review status: {status}")
        if (
            status == "CONFIRMED_PRIMARY_SOURCE"
            and text(publication["access_status"]) not in PRIMARY_ACCESS
        ):
            raise ValueError(
                f"interaction confirmation without primary full text: {interaction_id}"
            )
        if status != text(publication["publication_review_status"]):
            raise ValueError(f"publication/interaction review status mismatch: {interaction_id}")
        warning_codes = record.get("warning_codes")
        if not isinstance(warning_codes, list):
            raise TypeError(f"warning_codes must be a list: {interaction_id}")
        if status == "EVIDENCE_TYPE_MISMATCH" and "EVIDENCE_TYPE_MISMATCH" not in warning_codes:
            raise ValueError(f"evidence mismatch must have warning code: {interaction_id}")
        if text(record["gene_symbol"]) != text(publication["gene_symbol"]):
            raise ValueError(f"interaction/publication target mismatch: {interaction_id}")
        if FOCUS6.get(text(record["gene_symbol"])) != text(record["uniprot_accession"]):
            raise ValueError(f"interaction UniProt mismatch: {interaction_id}")
        if text(record["gene_symbol"]) == "FNTB" and status == "CONFIRMED_PRIMARY_SOURCE":
            combined = " ".join(
                [
                    text(record["protein_object_ru"]),
                    text(record["construct_ru"]),
                    text(record["evidence_scope_ru"]),
                ]
            ).lower()
            if "комплекс" not in combined and "гетеродимер" not in combined:
                raise ValueError(
                    f"confirmed FNTB evidence must name the FTase complex: {interaction_id}"
                )

        linked_ids = record.get("observation_ids")
        if not isinstance(linked_ids, list) or not linked_ids:
            raise ValueError(f"observation_ids must be a non-empty list: {interaction_id}")
        if len(linked_ids) != len(set(linked_ids)):
            raise ValueError(f"duplicate observation links: {interaction_id}")
        representative = text(record["representative_observation_id"])
        if representative not in linked_ids:
            raise ValueError(f"representative observation is not linked: {interaction_id}")
        for observation_id in linked_ids:
            row = observations.get(observation_id)
            if row is None:
                raise ValueError(f"unknown P0.9 observation: {observation_id}")
            if row["gene_symbol"] != text(record["gene_symbol"]):
                raise ValueError(f"observation target mismatch: {observation_id}")
            if normalized_doi(row["doi"]) != normalized_doi(publication["doi"]) or row[
                "pmid"
            ] != text(publication["pmid"]):
                raise ValueError(f"observation publication mismatch: {observation_id}")
            if row["endpoint"] != text(record["endpoint"]) or row["relation"] != text(
                record["relation"]
            ):
                raise ValueError(f"observation endpoint mismatch: {observation_id}")
            if normalized_numeric(row["numeric_value_nm"]) != normalized_numeric(
                record["numeric_value_nm"]
            ):
                raise ValueError(f"observation value mismatch: {observation_id}")
            if row.get("scientific_acceptance") != "NOT_APPROVED":
                raise ValueError(f"upstream observation acceptance mismatch: {observation_id}")
        interaction_by_id[interaction_id] = record

    reviewed_publication_ids = {text(row["publication_id"]) for row in interactions}
    if reviewed_publication_ids != set(publication_by_id):
        raise ValueError("each publication must have exactly one representative interaction")

    source_ids: set[str] = set()
    for source in sources:
        if not isinstance(source, dict):
            raise TypeError("source access log entry must be an object")
        _require_fields(
            source, ("source_id", "url", "accessed_utc", "source_role"), "source access log"
        )
        source_id = text(source["source_id"])
        if source_id in source_ids:
            raise ValueError(f"duplicate source_id: {source_id}")
        source_ids.add(source_id)

    for uncertainty in uncertainties:
        if not isinstance(uncertainty, dict):
            raise TypeError("uncertainty must be an object")
        _require_fields(
            uncertainty,
            (
                "uncertainty_id",
                "publication_id",
                "warning_code",
                "description_ru",
                "resolution_required_ru",
            ),
            "uncertainty",
        )
        if text(uncertainty["publication_id"]) not in publication_by_id:
            raise ValueError(f"unknown publication in uncertainty: {uncertainty['publication_id']}")
        interaction_id = text(uncertainty.get("interaction_id"))
        if interaction_id and interaction_id not in interaction_by_id:
            raise ValueError(f"unknown interaction in uncertainty: {interaction_id}")

    return {
        "publications": publications,
        "publication_by_id": publication_by_id,
        "interactions": interactions,
        "interaction_by_id": interaction_by_id,
        "uncertainties": uncertainties,
        "source_access_log": sources,
        "observations": observations,
        "priority_by_key": p09_priority,
    }


def build_outputs(curation: dict[str, Any], p09: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    validated = validate_curation(curation, p09)
    publication_rows: list[dict[str, Any]] = []
    for publication in validated["publications"]:
        key = (
            text(publication["gene_symbol"]),
            normalized_doi(publication["doi"]),
            text(publication["pmid"]),
        )
        priority = validated["priority_by_key"][key]
        publication_rows.append(
            {
                **{field: text(publication.get(field)) for field in PUBLICATION_FIELDS},
                "review_queue_rank": priority["review_queue_rank"],
                "uniprot_accession": FOCUS6[text(publication["gene_symbol"])],
                "warning_codes": join_values(publication.get("warning_codes", [])),
                "p09_observation_records_n": priority["observation_records_n"],
                "p09_exact_numeric_records_n": priority["exact_numeric_records_n"],
                "p09_censored_records_n": priority["censored_records_n"],
                "p09_sources": priority["sources"],
                "p09_endpoints": priority["endpoints"],
                "scientific_acceptance": "NOT_APPROVED",
            }
        )
    publication_rows.sort(key=lambda row: int(row["review_queue_rank"]))

    link_rows: list[dict[str, Any]] = []
    interaction_rows: list[dict[str, Any]] = []
    duplicate_rows: list[dict[str, Any]] = []
    duplicate_members_by_observation: dict[str, list[dict[str, str]]] = defaultdict(list)
    for member in p09["duplicate_members"]:
        duplicate_members_by_observation[member["observation_id"]].append(member)
    cross_pairs_by_observation: dict[str, list[dict[str, str]]] = defaultdict(list)
    for pair in p09["cross_pairs"]:
        cross_pairs_by_observation[pair["chembl_observation_id"]].append(pair)
        cross_pairs_by_observation[pair["bindingdb_observation_id"]].append(pair)

    for interaction in validated["interactions"]:
        publication = validated["publication_by_id"][text(interaction["publication_id"])]
        linked = [validated["observations"][item] for item in interaction["observation_ids"]]
        for row in linked:
            link_rows.append(
                {
                    "interaction_id": interaction["interaction_id"],
                    "publication_id": interaction["publication_id"],
                    "review_status": interaction["review_status"],
                    **{
                        field: row.get(field, "")
                        for field in LINK_FIELDS
                        if field not in {"interaction_id", "publication_id", "review_status"}
                    },
                    "scientific_acceptance": "NOT_APPROVED",
                }
            )
        interaction_rows.append(
            {
                **{field: text(interaction.get(field)) for field in INTERACTION_FIELDS},
                "warning_codes": join_values(interaction.get("warning_codes", [])),
                "linked_observations_n": len(linked),
                "linked_sources": join_values(row["source"] for row in linked),
                "source_molecule_ids": join_values(row["source_molecule_id"] for row in linked),
                "source_protein_ids": join_values(row["source_target_id"] for row in linked),
                "source_experiment_ids": join_values(row["source_experiment_id"] for row in linked),
                "source_measurement_ids": join_values(
                    row["source_measurement_id"] for row in linked
                ),
                "source_record_ids": join_values(row["source_record_id"] for row in linked),
                "doi": normalized_doi(publication["doi"]),
                "pmid": text(publication["pmid"]),
                "scientific_acceptance": "NOT_APPROVED",
            }
        )
        within_chembl: set[str] = set()
        within_bindingdb: set[str] = set()
        cross_groups: set[str] = set()
        cross_pair_ids: set[str] = set()
        for row in linked:
            for member in duplicate_members_by_observation.get(row["observation_id"], []):
                scope = member.get("scope", "")
                group_id = member.get("duplicate_group_id", "")
                if scope == "WITHIN_CHEMBL":
                    within_chembl.add(group_id)
                elif scope == "WITHIN_BINDINGDB":
                    within_bindingdb.add(group_id)
                elif scope == "CROSS_SOURCE":
                    cross_groups.add(group_id)
            for pair in cross_pairs_by_observation.get(row["observation_id"], []):
                cross_groups.add(pair.get("duplicate_group_id", ""))
                cross_pair_ids.add(pair.get("duplicate_pair_id", ""))
        duplicate_rows.append(
            {
                "interaction_id": interaction["interaction_id"],
                "publication_id": interaction["publication_id"],
                "gene_symbol": interaction["gene_symbol"],
                "curator_duplicate_assessment": interaction["duplicate_assessment"],
                "curator_rationale_ru": interaction["duplicate_rationale_ru"],
                "within_chembl_group_ids": join_values(within_chembl),
                "within_bindingdb_group_ids": join_values(within_bindingdb),
                "cross_source_group_ids": join_values(cross_groups),
                "cross_source_pair_ids": join_values(cross_pair_ids),
                "linked_observations_n": len(linked),
                "records_removed_n": 0,
                "review_status": "PRESERVED_NOT_DEDUPLICATED",
            }
        )

    uncertainty_rows = [
        {
            **{field: text(row.get(field)) for field in UNCERTAINTY_FIELDS},
            "review_status": "OPEN",
        }
        for row in validated["uncertainties"]
    ]
    interaction_rows.sort(key=lambda row: row["interaction_id"])
    link_rows.sort(key=lambda row: (row["interaction_id"], row["source"], row["observation_id"]))
    duplicate_rows.sort(key=lambda row: row["interaction_id"])
    uncertainty_rows.sort(key=lambda row: row["uncertainty_id"])
    return {
        "publications": publication_rows,
        "interactions": interaction_rows,
        "links": link_rows,
        "duplicates": duplicate_rows,
        "uncertainties": uncertainty_rows,
    }


def output_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, delimiter="\t", fieldnames=fields, lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def make_qc_rows(
    outputs: dict[str, list[dict[str, Any]]],
    p09_hash_unchanged: bool,
    curation_hash_unchanged: bool,
) -> list[dict[str, Any]]:
    status_counts = Counter(row["review_status"] for row in outputs["interactions"])
    checks = [
        (
            "P010_QC_001",
            len(outputs["publications"]) == 12,
            len(outputs["publications"]),
            12,
            "Проверены все 12 публикаций P0.9.",
        ),
        (
            "P010_QC_002",
            len(outputs["interactions"]) == 12,
            len(outputs["interactions"]),
            12,
            "Сформировано по одной репрезентативной связи на публикацию.",
        ),
        (
            "P010_QC_003",
            p09_hash_unchanged,
            str(p09_hash_unchanged).upper(),
            "TRUE",
            "ZIP P0.9 не изменён.",
        ),
        (
            "P010_QC_004",
            curation_hash_unchanged,
            str(curation_hash_unchanged).upper(),
            "TRUE",
            "Журнал курации не изменён программой.",
        ),
        (
            "P010_QC_005",
            all(row["scientific_acceptance"] == "NOT_APPROVED" for row in outputs["interactions"]),
            "NOT_APPROVED",
            "NOT_APPROVED",
            "Автоматическая научная приёмка запрещена.",
        ),
        (
            "P010_QC_006",
            all(row["records_removed_n"] == 0 for row in outputs["duplicates"]),
            0,
            0,
            "Потенциальные дубли не удалены.",
        ),
        (
            "P010_QC_007",
            status_counts["CONFIRMED_PRIMARY_SOURCE"] > 0,
            status_counts["CONFIRMED_PRIMARY_SOURCE"],
            ">0",
            "Есть связи, подтверждённые доступным первичным полным текстом.",
        ),
        (
            "P010_QC_008",
            status_counts["REQUIRES_REVIEW"] > 0,
            status_counts["REQUIRES_REVIEW"],
            ">0",
            "Ограничение доступа не скрыто.",
        ),
        (
            "P010_QC_009",
            status_counts["EVIDENCE_TYPE_MISMATCH"] > 0,
            status_counts["EVIDENCE_TYPE_MISMATCH"],
            ">0",
            "Несоответствия типа доказательства выделены отдельно.",
        ),
    ]
    return [
        {
            "check_id": check_id,
            "status": "PASS" if passed else "FAIL",
            "value": value,
            "expected": expected,
            "details_ru": details,
        }
        for check_id, passed, value, expected, details in checks
    ]


def report_ru(summary: dict[str, Any], outputs: dict[str, list[dict[str, Any]]]) -> str:
    status_counts = summary["review_status_counts"]
    target_counts = defaultdict(Counter)
    for row in outputs["interactions"]:
        target_counts[row["gene_symbol"]][row["review_status"]] += 1
    lines = [
        "# MCL TLS V3 / P0.10 — проверка первичных публикаций",
        "",
        "## Статус",
        "",
        f"- Контракт: `{CONTRACT}`.",
        "- Проверены 12 публикаций, выбранных P0.9; исходные ChEMBL и BindingDB не открывались.",
        "- Результаты известных лигандов не перенесены на 65 PYZ и не объединены в показатель эффективности.",
        "- `scientific_acceptance=NOT_APPROVED`; `merge_permitted=false`.",
        "",
        "## Критерии статуса",
        "",
        "- `CONFIRMED_PRIMARY_SOURCE`: доступен оригинальный полный текст или авторская рукопись; объект, опыт, конечная точка и значение сопоставлены в заявленном объёме.",
        "- `REQUIRES_REVIEW`: полный текст/методика недоступны либо существенная деталь конструкции или контроля не подтверждена.",
        "- `EVIDENCE_TYPE_MISMATCH`: первичный источник или описание опыта не подтверждают заявленную индивидуальную мишень либо тип доказательства.",
        "",
        "## Итог",
        "",
        f"- Подтверждено первоисточником: {status_counts.get('CONFIRMED_PRIMARY_SOURCE', 0)}.",
        f"- Требует проверки: {status_counts.get('REQUIRES_REVIEW', 0)}.",
        f"- Не соответствует заявленному типу доказательства: {status_counts.get('EVIDENCE_TYPE_MISMATCH', 0)}.",
        f"- Открытых неопределённостей: {summary['open_uncertainties_n']}.",
        "",
        "## По мишеням",
        "",
        "| Белок | Подтверждено | Требует проверки | Несоответствие |",
        "|---|---:|---:|---:|",
    ]
    for gene in FOCUS6:
        counts = target_counts[gene]
        lines.append(
            f"| {gene} | {counts['CONFIRMED_PRIMARY_SOURCE']} | {counts['REQUIRES_REVIEW']} | {counts['EVIDENCE_TYPE_MISMATCH']} |"
        )
    lines.extend(
        [
            "",
            "## Научно существенные выводы",
            "",
            "- ERBB3: две приоритетные публикации исследуют EGFR или совокупный объект семейства ERBB; индивидуальная активность ERBB3 из этих строк не подтверждается.",
            "- FNTB: подтверждённые результаты относятся к активному комплексу фарнезилтрансферазы FNTA/FNTB, а не к изолированной β-субъединице.",
            "- DDC: значение из работы о NVP-BHG712 получено Kinobead pull-down; это не очищенное прямое связывание DDC.",
            "- IC50 сохранены как контекстно-зависимые показатели ингибирования и не приравнивались к Ki или Kd.",
            "- Предполагаемые дубли сохранены как отдельные исходные записи со всеми идентификаторами.",
            "",
            "## Таблицы",
            "",
            "- `p010_publication_reviews.tsv` — аудит 12 публикаций и доступности методики.",
            "- `p010_reviewed_interactions.tsv` — 12 репрезентативных связей с научным статусом.",
            "- `p010_interaction_source_links.tsv` — исходные идентификаторы ChEMBL/BindingDB без удаления записей.",
            "- `p010_duplicate_assessment.tsv` — внутрибазовые и межбазовые потенциальные дубли.",
            "- `p010_uncertainties.tsv` — открытые вопросы и требуемые действия.",
            "",
            "## Ограничения",
            "",
            "- Для девяти публикаций не был доступен авторитетный полный первичный текст: для шести связей поэтому сохранён статус `REQUIRES_REVIEW`; для трёх других несоответствие типа доказательства установлено по метаданным статьи и опыта без подмены ими полнотекстовой проверки методики.",
            "- Каталог мал и репрезентативен: он не является исчерпывающим перечнем известных лигандов шести белков.",
            "- Сопоставление структуры между ChEMBL и BindingDB остаётся потенциальным дублем, пока не проверены исходные таблицы/схемы соединений публикации.",
            "- Научная приёмка требует ручного рассмотрения отчёта и открытых неопределённостей.",
            "",
        ]
    )
    return "\n".join(lines)


def verify_p010_package(path: Path) -> dict[str, Any]:
    """Validate the final P0.10 ZIP, including exact membership and checksums."""
    with ValidatedZip(path) as archive:
        manifest = read_json(
            archive.read("p010_package_manifest.json"), "p010_package_manifest.json"
        )
        if manifest.get("contract") != CONTRACT:
            raise ValueError("P0.10 package contract mismatch")
        if (
            manifest.get("scientific_acceptance") != "NOT_APPROVED"
            or manifest.get("merge_permitted") is not False
        ):
            raise ValueError("P0.10 package acceptance safeguards mismatch")
        files = manifest.get("files")
        if not isinstance(files, dict) or not files:
            raise ValueError("P0.10 package manifest has no file entries")
        expected_names = set(files) | {"p010_package_manifest.json"}
        if archive.names != expected_names:
            raise ValueError("P0.10 ZIP/manifest member set mismatch")
        for name, metadata in files.items():
            payload = archive.read(name)
            if len(payload) != metadata.get("size_bytes"):
                raise ValueError(f"P0.10 member size mismatch: {name}")
            if sha_bytes(payload) != metadata.get("sha256"):
                raise ValueError(f"P0.10 member SHA-256 mismatch: {name}")
        return {
            "contract": manifest["contract"],
            "files_n": len(files),
            "zip_sha256": sha_file(path),
            "scientific_acceptance": manifest["scientific_acceptance"],
            "merge_permitted": manifest["merge_permitted"],
        }


def build_review(
    p09_zip: Path,
    expected_p09_sha256: str,
    curation_path: Path,
    output_dir: Path,
    output_zip: Path,
    git_commit: str,
) -> dict[str, Any]:
    if output_dir.exists():
        raise FileExistsError(f"output directory already exists: {output_dir}")
    if output_zip.exists():
        raise FileExistsError(f"output ZIP already exists: {output_zip}")
    p09_sha_before = sha_file(p09_zip)
    curation_sha_before = sha_file(curation_path)
    p09 = verify_p09_report(p09_zip, expected_p09_sha256)
    curation = load_curation(curation_path)
    outputs = build_outputs(curation, p09)
    p09_sha_after = sha_file(p09_zip)
    curation_sha_after = sha_file(curation_path)
    if p09_sha_before != p09_sha_after or curation_sha_before != curation_sha_after:
        raise RuntimeError("input changed during P0.10")

    created_utc = datetime.now(UTC).isoformat().replace("+00:00", "Z")
    status_counts = Counter(row["review_status"] for row in outputs["interactions"])
    summary = {
        "contract": CONTRACT,
        "created_utc": created_utc,
        "git_commit": git_commit,
        "input_p09_zip_sha256": p09_sha_before,
        "curation_sha256": curation_sha_before,
        "publications_reviewed_n": len(outputs["publications"]),
        "representative_interactions_n": len(outputs["interactions"]),
        "source_observation_links_n": len(outputs["links"]),
        "review_status_counts": dict(sorted(status_counts.items())),
        "open_uncertainties_n": len(outputs["uncertainties"]),
        "records_removed_as_duplicates_n": 0,
        "multi_gigabyte_sources_opened": False,
        "network_used_by_builder": False,
        "known_ligand_results_transferred_to_pyz": False,
        "baseline_65x13_untouched": True,
        "scientific_acceptance": "NOT_APPROVED",
        "merge_permitted": False,
    }
    verification = {
        "contract": CONTRACT,
        "p09_zip": str(p09_zip.resolve()),
        "p09_zip_size_bytes": p09["zip_size_bytes"],
        "expected_p09_zip_sha256": require_sha256(expected_p09_sha256, "expected P0.9 SHA-256"),
        "actual_p09_zip_sha256_before": p09_sha_before,
        "actual_p09_zip_sha256_after": p09_sha_after,
        "p09_package_contract": p09["manifest"]["contract"],
        "p09_manifested_files_n": len(p09["manifest"]["files"]),
        "p09_zip_members_n": p09["member_names_n"],
        "curation_path": str(curation_path.resolve()),
        "curation_sha256_before": curation_sha_before,
        "curation_sha256_after": curation_sha_after,
        "source_files_unchanged": True,
        "scientific_acceptance": "NOT_APPROVED",
        "merge_permitted": False,
    }
    source_manifest = {
        "contract": CONTRACT,
        "created_utc": created_utc,
        "input_p09_zip": {
            "path": str(p09_zip.resolve()),
            "size_bytes": p09["zip_size_bytes"],
            "sha256_before": p09_sha_before,
            "sha256_after": p09_sha_after,
        },
        "curation_journal": {
            "path": str(curation_path.resolve()),
            "sha256_before": curation_sha_before,
            "sha256_after": curation_sha_after,
        },
        "publication_sources": curation["source_access_log"],
        "upstream_source_releases": p09["source_manifest"].get("upstream_source_releases", {}),
        "upstream_source_metadata": p09["source_manifest"].get("upstream_source_metadata", {}),
        "multi_gigabyte_sources_opened": False,
        "baseline_65x13_untouched": True,
        "scientific_acceptance": "NOT_APPROVED",
        "merge_permitted": False,
    }

    output_dir.mkdir(parents=True)
    output_tsv(
        output_dir / "p010_publication_reviews.tsv", outputs["publications"], PUBLICATION_FIELDS
    )
    output_tsv(
        output_dir / "p010_reviewed_interactions.tsv", outputs["interactions"], INTERACTION_FIELDS
    )
    output_tsv(output_dir / "p010_interaction_source_links.tsv", outputs["links"], LINK_FIELDS)
    output_tsv(
        output_dir / "p010_duplicate_assessment.tsv", outputs["duplicates"], DUPLICATE_FIELDS
    )
    output_tsv(output_dir / "p010_uncertainties.tsv", outputs["uncertainties"], UNCERTAINTY_FIELDS)
    qc_rows = make_qc_rows(
        outputs, p09_sha_before == p09_sha_after, curation_sha_before == curation_sha_after
    )
    output_tsv(output_dir / "p010_QC.tsv", qc_rows, QC_FIELDS)
    if any(row["status"] != "PASS" for row in qc_rows):
        raise RuntimeError("P0.10 QC failed")
    write_json(output_dir / "p010_summary.json", summary)
    write_json(output_dir / "p010_input_verification.json", verification)
    write_json(output_dir / "p010_source_manifest.json", source_manifest)
    (output_dir / "P010_REPORT_RU.md").write_text(report_ru(summary, outputs), encoding="utf-8")

    package_files = sorted(path for path in output_dir.iterdir() if path.is_file())
    package_manifest = {
        "contract": CONTRACT,
        "created_utc": created_utc,
        "input_p09_zip_sha256": p09_sha_before,
        "curation_sha256": curation_sha_before,
        "scientific_acceptance": "NOT_APPROVED",
        "merge_permitted": False,
        "files": {
            path.name: {"size_bytes": path.stat().st_size, "sha256": sha_file(path)}
            for path in package_files
        },
    }
    write_json(output_dir / "p010_package_manifest.json", package_manifest)
    with zipfile.ZipFile(
        output_zip, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=9
    ) as archive:
        for path in sorted(output_dir.iterdir()):
            archive.write(path, arcname=path.name)
    package_verification = verify_p010_package(output_zip)
    return {
        **summary,
        "output_zip": str(output_zip.resolve()),
        "output_zip_sha256": package_verification["zip_sha256"],
        "output_manifested_files_n": package_verification["files_n"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p09-zip", type=Path, required=True)
    parser.add_argument("--expected-p09-sha256", required=True)
    parser.add_argument("--curation", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--output-zip", type=Path, required=True)
    parser.add_argument("--git-commit", required=True)
    args = parser.parse_args()
    result = build_review(
        args.p09_zip,
        args.expected_p09_sha256,
        args.curation,
        args.output_dir,
        args.output_zip,
        args.git_commit,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
