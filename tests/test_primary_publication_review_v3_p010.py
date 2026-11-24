"""Synthetic P0.10 tests; no network or multi-gigabyte source databases."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import sys
import warnings
import zipfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from build_primary_publication_review_v3_p010 import (
    CURATION_CONTRACT,
    FOCUS6,
    P09_CONTRACT,
    build_outputs,
    build_review,
    load_curation,
    validate_curation,
    verify_p010_package,
    verify_p09_report,
)

OBS_FIELDS = [
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
    "doi",
    "pmid",
    "endpoint",
    "relation",
    "numeric_value_nm",
    "raw_value_nm",
    "protein_match_status",
    "species_status",
    "experiment_type_status",
    "measurement_precision_status",
    "scientific_acceptance",
]

PRIORITY_FIELDS = [
    "review_queue_rank",
    "gene_symbol",
    "publication_key",
    "publication_url",
    "doi",
    "pmid",
    "sources",
    "endpoints",
    "observation_records_n",
    "exact_numeric_records_n",
    "censored_records_n",
]

LINK_FIELDS = [
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

DUP_GROUP_FIELDS = [
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

DUP_MEMBER_FIELDS = [
    "duplicate_group_id",
    "scope",
    "observation_id",
    "source",
    "source_record_id",
]

CROSS_FIELDS = [
    "duplicate_pair_id",
    "duplicate_group_id",
    "chembl_observation_id",
    "bindingdb_observation_id",
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


def tsv_bytes(rows: list[dict[str, str]], fields: list[str]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, delimiter="\t", fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def synthetic_records() -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    genes = list(FOCUS6) * 2
    observations: list[dict[str, str]] = []
    priorities: list[dict[str, str]] = []
    links: list[dict[str, str]] = []
    for index, gene in enumerate(genes, start=1):
        doi = f"10.1000/p010.{index}"
        pmid = str(1000 + index)
        observation_id = f"CHEMBL:activity:{index}"
        observations.append(
            {
                "observation_id": observation_id,
                "source": "ChEMBL",
                "gene_symbol": gene,
                "uniprot_accession": FOCUS6[gene],
                "source_target_id": f"CHEMBL_TARGET_{gene}",
                "source_measurement_id": str(index),
                "source_experiment_id": f"CHEMBL_ASSAY_{index}",
                "source_molecule_id": f"CHEMBL_MOL_{index}",
                "source_record_id": str(index),
                "source_document_id": f"CHEMBL_DOC_{index}",
                "doi": doi,
                "pmid": pmid,
                "endpoint": "IC50",
                "relation": "=",
                "numeric_value_nm": str(index),
                "raw_value_nm": str(index),
                "protein_match_status": "EXACT_SINGLE_PROTEIN_UNIPROT_MATCH",
                "species_status": "HUMAN_SPECIES_MATCH",
                "experiment_type_status": "BINDING_ASSAY_METADATA_MATCHED",
                "measurement_precision_status": "EXACT_NUMERIC",
                "scientific_acceptance": "NOT_APPROVED",
            }
        )
        priorities.append(
            {
                "review_queue_rank": str(index),
                "gene_symbol": gene,
                "publication_key": f"DOI:{doi}",
                "publication_url": f"https://doi.org/{doi}",
                "doi": doi,
                "pmid": pmid,
                "sources": "ChEMBL",
                "endpoints": "IC50",
                "observation_records_n": "1",
                "exact_numeric_records_n": "1",
                "censored_records_n": "0",
            }
        )
        links.append(
            {
                "review_queue_rank": str(index),
                "gene_symbol": gene,
                "publication_key": f"DOI:{doi}",
                "observation_id": observation_id,
                "source": "ChEMBL",
                "endpoint": "IC50",
                "catalog_inclusion_status": "REVIEW_CANDIDATE",
                "warning_codes": "PRIMARY_EXPERIMENT_REVIEW_REQUIRED",
                "source_measurement_id": str(index),
                "source_experiment_id": f"CHEMBL_ASSAY_{index}",
                "source_molecule_id": f"CHEMBL_MOL_{index}",
            }
        )
    return observations, priorities, links


def make_p09_zip(
    tmp_path: Path, *, bad_manifest_hash: bool = False, duplicate_member: bool = False
) -> Path:
    observations, priorities, links = synthetic_records()
    files = {
        "p09_summary.json": json.dumps(
            {
                "contract": P09_CONTRACT,
                "scientific_acceptance": "NOT_APPROVED",
                "merge_permitted": False,
            }
        ).encode(),
        "p09_source_manifest.json": json.dumps(
            {
                "contract": P09_CONTRACT,
                "scientific_acceptance": "NOT_APPROVED",
                "merge_permitted": False,
                "multi_gigabyte_sources_opened": False,
                "upstream_source_releases": {"chembl": "chembl_37", "bindingdb": "202611"},
            }
        ).encode(),
        "p09_observation_assessment.tsv": tsv_bytes(observations, OBS_FIELDS),
        "p09_priority_publications.tsv": tsv_bytes(priorities, PRIORITY_FIELDS),
        "p09_priority_publication_observation_links.tsv": tsv_bytes(links, LINK_FIELDS),
        "p09_duplicate_groups.tsv": tsv_bytes([], DUP_GROUP_FIELDS),
        "p09_duplicate_group_members.tsv": tsv_bytes([], DUP_MEMBER_FIELDS),
        "p09_cross_source_duplicate_pairs.tsv": tsv_bytes([], CROSS_FIELDS),
    }
    manifest_files = {
        name: {
            "size_bytes": len(payload),
            "sha256": (
                "0" * 64
                if bad_manifest_hash and name == "p09_summary.json"
                else hashlib.sha256(payload).hexdigest()
            ),
        }
        for name, payload in files.items()
    }
    manifest = {
        "contract": P09_CONTRACT,
        "scientific_acceptance": "NOT_APPROVED",
        "merge_permitted": False,
        "files": manifest_files,
    }
    path = tmp_path / "p09.zip"
    mode = "w"
    with zipfile.ZipFile(path, mode, compression=zipfile.ZIP_DEFLATED) as archive:
        for name, payload in files.items():
            archive.writestr(name, payload)
        archive.writestr("p09_package_manifest.json", json.dumps(manifest).encode())
        if duplicate_member:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                archive.writestr("p09_summary.json", files["p09_summary.json"])
    return path


def make_curation(tmp_path: Path) -> Path:
    genes = list(FOCUS6) * 2
    publications = []
    interactions = []
    for index, gene in enumerate(genes, start=1):
        status = (
            "CONFIRMED_PRIMARY_SOURCE"
            if index in {1, 11}
            else "EVIDENCE_TYPE_MISMATCH"
            if index == 7
            else "REQUIRES_REVIEW"
        )
        access = (
            "FULL_TEXT_PRIMARY_OPEN"
            if status == "CONFIRMED_PRIMARY_SOURCE"
            else "PUBLISHER_METADATA_AND_ABSTRACT_ONLY"
        )
        publication_id = f"P010-PUB-{index:02d}"
        publications.append(
            {
                "publication_id": publication_id,
                "gene_symbol": gene,
                "doi": f"10.1000/p010.{index}",
                "pmid": str(1000 + index),
                "title": f"Synthetic publication {index}",
                "publisher_url": f"https://example.test/publisher/{index}",
                "pubmed_url": f"https://example.test/pubmed/{index}",
                "full_text_url": f"https://example.test/full/{index}"
                if access == "FULL_TEXT_PRIMARY_OPEN"
                else "",
                "access_status": access,
                "access_limitation_ru": "Ограничение доступа."
                if access != "FULL_TEXT_PRIMARY_OPEN"
                else "Полный текст доступен.",
                "protein_object_ru": "Комплекс фарнезилтрансферазы FNTA/FNTB"
                if gene == "FNTB"
                else f"Белок {gene}",
                "species_ru": "Homo sapiens",
                "construct_ru": "Функциональный гетеродимерный комплекс"
                if gene == "FNTB"
                else "Синтетическая конструкция",
                "assay_class": "BIOCHEMICAL_ENZYME_INHIBITION",
                "endpoint_interpretation_ru": "Контекстный IC50.",
                "method_ru": "Синтетический метод.",
                "conditions_ru": "Синтетические условия.",
                "controls_ru": "Синтетический контроль.",
                "primary_source_locator": "Synthetic source locator",
                "publication_review_status": status,
                "conclusion_ru": "Синтетический вывод.",
                "warning_codes": ["EVIDENCE_TYPE_MISMATCH"]
                if status == "EVIDENCE_TYPE_MISMATCH"
                else [],
            }
        )
        interactions.append(
            {
                "interaction_id": f"P010-INT-{index:03d}",
                "publication_id": publication_id,
                "gene_symbol": gene,
                "uniprot_accession": FOCUS6[gene],
                "ligand_label": f"Ligand {index}",
                "protein_object_ru": "Комплекс фарнезилтрансферазы FNTA/FNTB"
                if gene == "FNTB"
                else f"Белок {gene}",
                "organism_ru": "Homo sapiens",
                "construct_ru": "Функциональный гетеродимерный комплекс"
                if gene == "FNTB"
                else "Синтетическая конструкция",
                "assay_class": "BIOCHEMICAL_ENZYME_INHIBITION",
                "endpoint": "IC50",
                "relation": "=",
                "numeric_value_nm": str(index),
                "reported_value_text": f"{index} nM",
                "evidence_scope_ru": "Комплекс FTase"
                if gene == "FNTB"
                else "Синтетическая область",
                "primary_source_locator": "Synthetic source locator",
                "review_status": status,
                "rationale_ru": "Синтетическое обоснование.",
                "warning_codes": ["EVIDENCE_TYPE_MISMATCH"]
                if status == "EVIDENCE_TYPE_MISMATCH"
                else [],
                "representative_observation_id": f"CHEMBL:activity:{index}",
                "observation_ids": [f"CHEMBL:activity:{index}"],
                "duplicate_assessment": "NO_SELECTED_DUPLICATE_SINGLE_SOURCE",
                "duplicate_rationale_ru": "Одна синтетическая запись.",
            }
        )
    value = {
        "contract": CURATION_CONTRACT,
        "scientific_acceptance": "NOT_APPROVED",
        "merge_permitted": False,
        "publications": publications,
        "interactions": interactions,
        "uncertainties": [],
        "source_access_log": [
            {
                "source_id": "SRC-1",
                "url": "https://example.test/source",
                "accessed_utc": "2026-11-24",
                "source_role": "synthetic_fixture",
            }
        ],
    }
    path = tmp_path / "curation.json"
    path.write_text(json.dumps(value, ensure_ascii=False), encoding="utf-8")
    return path


def verified_inputs(tmp_path: Path) -> tuple[Path, str, Path, dict[str, object], dict[str, object]]:
    p09_zip = make_p09_zip(tmp_path)
    p09_sha = hashlib.sha256(p09_zip.read_bytes()).hexdigest()
    curation_path = make_curation(tmp_path)
    p09 = verify_p09_report(p09_zip, p09_sha)
    curation = load_curation(curation_path)
    return p09_zip, p09_sha, curation_path, p09, curation


def test_end_to_end_preserves_sources_and_acceptance(tmp_path: Path) -> None:
    p09_zip, p09_sha, curation_path, _, _ = verified_inputs(tmp_path)
    before_p09 = p09_zip.read_bytes()
    before_curation = curation_path.read_bytes()
    output_dir = tmp_path / "out"
    output_zip = tmp_path / "out.zip"

    result = build_review(p09_zip, p09_sha, curation_path, output_dir, output_zip, "deadbeef")

    assert result["publications_reviewed_n"] == 12
    assert result["representative_interactions_n"] == 12
    assert result["scientific_acceptance"] == "NOT_APPROVED"
    assert result["merge_permitted"] is False
    assert p09_zip.read_bytes() == before_p09
    assert curation_path.read_bytes() == before_curation
    with zipfile.ZipFile(output_zip) as archive:
        names = set(archive.namelist())
        assert "p010_package_manifest.json" in names
        assert "p010_interaction_source_links.tsv" in names
    verified = verify_p010_package(output_zip)
    assert verified["files_n"] == 10


def test_duplicate_zip_member_is_rejected(tmp_path: Path) -> None:
    path = make_p09_zip(tmp_path, duplicate_member=True)
    with pytest.raises(ValueError, match="duplicate ZIP members"):
        verify_p09_report(path, hashlib.sha256(path.read_bytes()).hexdigest())


def test_corrupt_manifest_hash_is_rejected_without_breaking_zip(tmp_path: Path) -> None:
    path = make_p09_zip(tmp_path, bad_manifest_hash=True)
    with zipfile.ZipFile(path) as archive:
        assert archive.testzip() is None
    with pytest.raises(ValueError, match="member SHA-256 mismatch"):
        verify_p09_report(path, hashlib.sha256(path.read_bytes()).hexdigest())


def test_confirmation_without_primary_full_text_is_rejected(tmp_path: Path) -> None:
    _, _, _, p09, curation = verified_inputs(tmp_path)
    curation["publications"][0]["access_status"] = "PUBLISHER_METADATA_AND_ABSTRACT_ONLY"
    with pytest.raises(ValueError, match="full primary text is required"):
        validate_curation(curation, p09)


def test_value_mismatch_with_p09_is_rejected(tmp_path: Path) -> None:
    _, _, _, p09, curation = verified_inputs(tmp_path)
    curation["interactions"][0]["numeric_value_nm"] = "999"
    with pytest.raises(ValueError, match="observation value mismatch"):
        build_outputs(curation, p09)


def test_confirmed_fntb_cannot_be_upgraded_to_isolated_subunit(tmp_path: Path) -> None:
    _, _, _, p09, curation = verified_inputs(tmp_path)
    fntb = next(
        row
        for row in curation["interactions"]
        if row["gene_symbol"] == "FNTB" and row["review_status"] == "CONFIRMED_PRIMARY_SOURCE"
    )
    fntb["protein_object_ru"] = "Изолированная FNTB"
    fntb["construct_ru"] = "Изолированная субъединица"
    fntb["evidence_scope_ru"] = "Изолированная субъединица"
    with pytest.raises(ValueError, match="must name the FTase complex"):
        validate_curation(curation, p09)
