"""Synthetic P0.9 tests; no network or multi-gigabyte scientific sources."""
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
from build_known_ligand_catalog_v3_p09 import (
    ASSESSMENT_FIELDS,
    BINDINGDB_FIELDS,
    CHEMBL_FIELDS,
    CROSS_FIELDS,
    FOCUS6,
    P08_CONTRACT,
    P08_PACKAGE_CONTRACT,
    assign_duplicate_groups,
    build_catalog,
    build_observations,
    sha_file,
    verify_p08_report,
)


def tsv_bytes(rows: list[dict[str, str]], fields: list[str]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, delimiter="\t", fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


def chembl_row(
    activity_id: str,
    gene: str,
    endpoint: str,
    relation: str,
    value: str,
    *,
    molecule: str,
    smiles: str,
    doi: str = "",
    pmid: str = "",
    potential_duplicate: str = "0",
    validity: str = "",
) -> dict[str, str]:
    row = dict.fromkeys(CHEMBL_FIELDS, "")
    row.update(
        {
            "gene_symbol": gene,
            "uniprot_accession": FOCUS6[gene],
            "chembl_target_id": f"CHEMBL_TARGET_{gene}",
            "activity_id": activity_id,
            "assay_chembl_id": f"CHEMBL_ASSAY_{activity_id}",
            "molecule_chembl_id": molecule,
            "canonical_smiles": smiles,
            "standard_type": endpoint,
            "standard_relation": relation,
            "standard_value": value,
            "standard_units": "nM",
            "assay_type": "B",
            "confidence_score": "9",
            "potential_duplicate": potential_duplicate,
            "data_validity_comment": validity,
            "document_chembl_id": f"CHEMBL_DOC_{activity_id}",
            "doi": doi,
            "pmid": pmid,
            "measurement_value_status": (
                "INVALID_OR_MISSING" if not value else ("EXACT_NUMERIC" if relation == "=" else "CENSORED")
            ),
            "evidence_status": "ASSAY_METADATA_MATCHED_PRIMARY_REVIEW_REQUIRED",
            "primary_experiment_review_status": "REVIEW_REQUIRED",
            "cross_source_potential_duplicate": "FALSE",
            "scientific_acceptance": "NOT_APPROVED",
        }
    )
    return row


def binding_row(
    source_row: str,
    gene: str,
    endpoint: str,
    relation: str,
    value: str,
    *,
    molecule: str,
    smiles: str,
    inchikey: str,
    doi: str = "",
    pmid: str = "",
    chain_status: str = "SWISSPROT_PRIMARY_SINGLE_CHAIN_MATCH",
    species_status: str = "HUMAN_SPECIES_MATCH",
) -> dict[str, str]:
    row = dict.fromkeys(BINDINGDB_FIELDS, "")
    raw_field = {"Ki": "raw_ki_nm", "Kd": "raw_kd_nm", "IC50": "raw_ic50_nm"}[endpoint]
    row.update(
        {
            "gene_symbol": gene,
            "uniprot_accession": FOCUS6[gene],
            "source_row": source_row,
            "bindingdb_reactant_set_id": f"RS{source_row}",
            "bindingdb_monomer_id": molecule,
            "source_row_field_count": "640",
            "expected_field_count": "640",
            "source_row_width_status": "EXPECTED_WIDTH",
            "raw_chain_accessions": FOCUS6[gene],
            "all_chain_accessions": FOCUS6[gene] + ("|P99999" if chain_status == "REVIEW_MULTI_CHAIN" else ""),
            "populated_chain_count": "2" if chain_status == "REVIEW_MULTI_CHAIN" else "1",
            "declared_chain_count": "2" if chain_status == "REVIEW_MULTI_CHAIN" else "1",
            "target_name": f"Synthetic {gene}",
            "organism": (
                "" if species_status == "REVIEW_SPECIES_MISSING" else "Mus musculus" if species_status == "REVIEW_SPECIES_MISMATCH" else "Homo sapiens"
            ),
            raw_field: ("" if relation == "=" else relation) + value,
            "endpoint": endpoint,
            "raw_value_nm": ("" if relation == "=" else relation) + value,
            "relation": relation,
            "numeric_value_nm": value,
            "ligand_smiles": smiles,
            "ligand_inchikey": inchikey,
            "article_doi": doi,
            "pmid": pmid,
            "curation_source": "Synthetic fixture",
            "chain_status": chain_status,
            "species_status": species_status,
            "measurement_value_status": "EXACT_NUMERIC" if relation == "=" else "CENSORED",
            "measurement_status": "EXACT_NUMERIC_STRUCTURE_PRESENT" if relation == "=" else "CENSORED_NOT_EXACT",
            "primary_experiment_review_status": "REVIEW_REQUIRED",
            "cross_source_potential_duplicate": "FALSE",
            "scientific_acceptance": "NOT_APPROVED",
        }
    )
    return row


def synthetic_rows() -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    chembl = [
        chembl_row("1", "PRKDC", "Ki", "=", "12.0", molecule="CHEMBL1", smiles="CCO", doi="10.1000/strict", pmid="101"),
        chembl_row("2", "PRKDC", "IC50", "<", "50", molecule="CHEMBL1", smiles="CCO"),
        chembl_row("3", "BCL2L1", "Ki", "=", "25", molecule="CHEMBL2", smiles="CCC", doi="10.1000/dup"),
        chembl_row("4", "BCL2L1", "Ki", "=", "25.0", molecule="CHEMBL2", smiles="CCC", doi="https://doi.org/10.1000/dup"),
        chembl_row("5", "ERBB3", "Kd", "=", "", molecule="CHEMBL3", smiles="", doi="10.1000/invalid"),
        chembl_row("6", "FNTB", "IC50", "=", "70", molecule="CHEMBL4", smiles="CCN", doi="10.1000/cross", pmid="606"),
        chembl_row("7", "IGF1R", "Kd", "=", "8", molecule="CHEMBL5", smiles="CCCl", doi="10.1000/igf"),
        chembl_row("8", "DDC", "Ki", "=", "40", molecule="CHEMBL6", smiles="CCBr", doi="10.1000/ddc"),
    ]
    binding = [
        binding_row("101", "PRKDC", "Ki", "=", "12", molecule="MON1", smiles="CCO", inchikey="KEY1", doi="10.1000/strict", pmid="101"),
        binding_row("102", "ERBB3", "IC50", "=", "90", molecule="MON2", smiles="CCCC", inchikey="KEY2", doi="10.1000/erbb3", chain_status="REVIEW_MULTI_CHAIN"),
        binding_row("103", "FNTB", "IC50", "=", "70.0", molecule="MON3", smiles="CCN", inchikey="KEY3", doi="10.1000/cross", pmid="606"),
        binding_row("104", "DDC", "Ki", "=", "40", molecule="MON4", smiles="CCBr", inchikey="KEY4", doi="10.1000/ddc", species_status="REVIEW_SPECIES_MISSING"),
        binding_row("105", "BCL2L1", "Kd", "=", "15", molecule="MON5", smiles="CCF", inchikey="KEY5", doi="10.1000/bcl", chain_status="REVIEW_MULTI_CHAIN", species_status="REVIEW_SPECIES_MISMATCH"),
        binding_row("106", "IGF1R", "Kd", ">", "500", molecule="MON6", smiles="CCI", inchikey="KEY6", doi="10.1000/igf"),
    ]
    cross = [
        {
            "gene_symbol": "FNTB",
            "doi": "10.1000/cross",
            "endpoint": "IC50",
            "relation": "=",
            "numeric_value_nm": "70",
            "ligand_smiles": "CCN",
            "chembl_activity_id": "6",
            "bindingdb_source_row": "103",
            "assessment": "POTENTIAL_DUPLICATE_EXACT_FIELDS_PRIMARY_REVIEW_REQUIRED",
        }
    ]
    return chembl, binding, cross


def make_p08_zip(tmp_path: Path) -> Path:
    chembl, binding, cross = synthetic_rows()
    panel = [
        {
            "gene_symbol": gene,
            "uniprot_accession": accession,
            "mapping_status": "EXACT_REVIEWED_WITH_HGNC",
            "uniprot_hgnc": f"HGNC:{index}",
            "sample_stratum": "CONTROL",
            "source_identity_status": "VERIFIED_UNIPROT_HGNC",
        }
        for index, (gene, accession) in enumerate(FOCUS6.items(), start=1)
    ]
    source_metadata = {
        "chembl": {"release": "chembl_37", "actual_input_sha256": "1" * 64},
        "bindingdb": {"release": "202611", "actual_input_sha256": "2" * 64},
    }
    summary = {
        "contract": P08_CONTRACT,
        "scientific_acceptance": "NOT_APPROVED",
        "merge_permitted": False,
        "external_requests_n": 0,
        "panel_mode": "focus6",
        "panel_genes_n": 6,
        "source_status": {"chembl": "LOCAL_SOURCE_CHECKED", "bindingdb": "LOCAL_SOURCE_CHECKED"},
        "releases": {"chembl": "chembl_37", "bindingdb": "202611"},
        "chembl_rows_reviewed_n": len(chembl),
        "bindingdb_measurement_observations_retained_n": len(binding),
        "bindingdb_source_rows_short_padded_n": 0,
        "bindingdb_source_rows_long_rejected_n": 0,
        "cross_source_potential_duplicate_pairs_n": len(cross),
        "source_input_hashes_unchanged": True,
        "baseline_65x13_untouched": True,
        "source_metadata": source_metadata,
    }
    validation = {
        "base_zip_crc": "PASSED",
        "base_zip_duplicate_members_n": 0,
        "base_zip_unsafe_members_n": 0,
        "base_zip_vs_directory_mismatches_n": 0,
        "sqlite_audit_status": "PASSED",
        "source_hashes_unchanged": True,
        "chembl_observation_rows_n": len(chembl),
        "bindingdb_observation_rows_n": len(binding),
        "chembl_sqlite_value_mismatches_n": 0,
        "bindingdb_raw_endpoint_value_mismatches_n": 0,
        "missing_identifier_counts": {
            "chembl": {"gene_symbol": 0, "uniprot_accession": 0, "activity_id": 0},
            "bindingdb": {"gene_symbol": 0, "uniprot_accession": 0, "source_row": 0},
        },
        "truncated_observations_n": 0,
        "qc_rows_n": 0,
        "all_observations_not_approved": True,
        "cross_source_potential_duplicate_pairs_n": len(cross),
        "technical_validation_status": "PASSED",
    }
    source_manifest = {
        "contract": P08_CONTRACT,
        "input_paths_sha256": {
            "synthetic_chembl.db": {"sha256_before": "1" * 64, "sha256_after": "1" * 64},
            "synthetic_bindingdb.zip": {"sha256_before": "2" * 64, "sha256_after": "2" * 64},
        },
        "source_metadata": source_metadata,
        "no_network": True,
        "read_only_scientific_data": True,
        "exact_target_uniprot_required": True,
        "no_automatic_scientific_acceptance": True,
    }
    files = {
        "p08_summary.json": (json.dumps(summary) + "\n").encode(),
        "p08_combined_validation.json": (json.dumps(validation) + "\n").encode(),
        "p08_source_manifest.json": (json.dumps(source_manifest) + "\n").encode(),
        "p08_target_panel.tsv": tsv_bytes(
            panel,
            ["gene_symbol", "uniprot_accession", "mapping_status", "uniprot_hgnc", "sample_stratum", "source_identity_status"],
        ),
        "p08_chembl_observations.tsv": tsv_bytes(chembl, CHEMBL_FIELDS),
        "p08_bindingdb_observations.tsv": tsv_bytes(binding, BINDINGDB_FIELDS),
        "p08_cross_source_potential_duplicates.tsv": tsv_bytes(cross, CROSS_FIELDS),
    }
    package = {
        "contract": P08_PACKAGE_CONTRACT,
        "validation_status": "PASSED",
        "source_releases": {"chembl": "chembl_37", "bindingdb": "202611"},
        "source_input_hashes_unchanged": True,
        "scientific_acceptance": "NOT_APPROVED",
        "merge_permitted": False,
        "files": {
            name: {"size_bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}
            for name, content in files.items()
        },
    }
    files["validated_package_manifest.json"] = (json.dumps(package) + "\n").encode()
    output = tmp_path / "p08_validated.zip"
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, content in files.items():
            archive.writestr(name, content)
    return output


@pytest.fixture
def p08_zip(tmp_path: Path) -> Path:
    return make_p08_zip(tmp_path)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def run_catalog(tmp_path: Path, p08_zip: Path):
    output_dir = tmp_path / "p09"
    output_zip = tmp_path / "p09.zip"
    summary = build_catalog(
        p08_zip=p08_zip,
        expected_p08_sha256=sha_file(p08_zip),
        output_dir=output_dir,
        output_zip=output_zip,
        git_commit="synthetic",
        created_utc="2026-11-24T00:00:00Z",
    )
    return summary, output_dir, output_zip


@pytest.fixture(scope="module")
def catalog_output(tmp_path_factory: pytest.TempPathFactory):
    root = tmp_path_factory.mktemp("p09_catalog")
    input_zip = make_p08_zip(root)
    input_hash = sha_file(input_zip)
    summary, output_dir, output_zip = run_catalog(root, input_zip)
    return summary, output_dir, output_zip, input_zip, input_hash


def test_p08_zip_manifest_crc_and_table_hashes_are_strictly_verified(p08_zip: Path):
    verified = verify_p08_report(p08_zip, sha_file(p08_zip))
    assert verified["zip_crc_status"] == "PASSED"
    assert verified["manifest_hash_or_size_mismatches_n"] == 0
    assert len(verified["chembl_rows"]) == 8
    assert len(verified["binding_rows"]) == 6


def test_p08_outer_sha_mismatch_is_rejected(p08_zip: Path):
    with pytest.raises(ValueError, match="P0.8 ZIP SHA-256 mismatch"):
        verify_p08_report(p08_zip, "0" * 64)


def test_p08_internal_table_corruption_is_rejected(tmp_path: Path, p08_zip: Path):
    corrupt = tmp_path / "corrupt.zip"
    with zipfile.ZipFile(p08_zip) as source, zipfile.ZipFile(corrupt, "w") as target:
        for info in source.infolist():
            content = source.read(info.filename)
            if info.filename == "p08_chembl_observations.tsv":
                content += b"corruption"
            target.writestr(info.filename, content)
    with pytest.raises(ValueError, match="hash or size mismatch"):
        verify_p08_report(corrupt, sha_file(corrupt))


def test_p08_duplicate_zip_member_is_rejected(tmp_path: Path, p08_zip: Path):
    duplicate = tmp_path / "duplicate.zip"
    with zipfile.ZipFile(p08_zip) as source, zipfile.ZipFile(duplicate, "w") as target:
        for info in source.infolist():
            target.writestr(info.filename, source.read(info.filename))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            target.writestr("p08_summary.json", source.read("p08_summary.json"))
    with pytest.raises(ValueError, match="duplicate member names"):
        verify_p08_report(duplicate, sha_file(duplicate))


def test_cross_source_duplicate_signature_must_match_observations(p08_zip: Path):
    verified = verify_p08_report(p08_zip, sha_file(p08_zip))
    observations = build_observations(verified["chembl_rows"], verified["binding_rows"])
    verified["cross_rows"][0]["numeric_value_nm"] = "71"
    with pytest.raises(ValueError, match="exact-field signature mismatch"):
        assign_duplicate_groups(observations, verified["cross_rows"])


def test_quality_tiers_keep_bindingdb_assay_context_under_review(catalog_output):
    _summary, output_dir, _output_zip, _input_zip, _input_hash = catalog_output
    rows = {row["observation_id"]: row for row in read_tsv(output_dir / "p09_observation_assessment.tsv")}
    assert rows["CHEMBL:activity:1"]["catalog_inclusion_status"] == "STRICT_CANDIDATE"
    assert rows["CHEMBL:activity:1"]["publication_url"] == "https://doi.org/10.1000/strict"
    assert rows["BINDINGDB:source_row:101:endpoint:Ki"]["catalog_inclusion_status"] == "REVIEW_CANDIDATE"
    assert rows["BINDINGDB:source_row:101:endpoint:Ki"]["experiment_type_status"] == "REPORTED_AFFINITY_ENDPOINT_ASSAY_UNVERIFIED"
    assert "CELLULAR_VS_BIOCHEMICAL_CONTEXT_NOT_VERIFIED" in rows["BINDINGDB:source_row:101:endpoint:Ki"]["warning_codes"]


def test_multicomponent_and_species_uncertain_rows_are_context_only(catalog_output):
    _summary, output_dir, _output_zip, _input_zip, _input_hash = catalog_output
    rows = {row["observation_id"]: row for row in read_tsv(output_dir / "p09_observation_assessment.tsv")}
    assert rows["BINDINGDB:source_row:102:endpoint:IC50"]["catalog_inclusion_status"] == "CONTEXT_ONLY"
    assert rows["BINDINGDB:source_row:104:endpoint:Ki"]["catalog_inclusion_status"] == "CONTEXT_ONLY"
    assert rows["BINDINGDB:source_row:105:endpoint:Kd"]["species_status"] == "REVIEW_SPECIES_MISMATCH"


def test_invalid_structure_and_value_are_retained_but_excluded(catalog_output):
    summary, output_dir, _output_zip, _input_zip, _input_hash = catalog_output
    rows = {row["observation_id"]: row for row in read_tsv(output_dir / "p09_observation_assessment.tsv")}
    invalid = rows["CHEMBL:activity:5"]
    assert invalid["catalog_inclusion_status"] == "EXCLUDED_INVALID"
    assert invalid["raw_value_nm"] == ""
    assert invalid["source_measurement_id"] == "5"
    assert summary["observation_records_n"] == 14


def test_censored_values_and_endpoints_remain_separate(catalog_output):
    _summary, output_dir, _output_zip, _input_zip, _input_hash = catalog_output
    rows = read_tsv(output_dir / "p09_observation_assessment.tsv")
    censored = next(row for row in rows if row["observation_id"] == "CHEMBL:activity:2")
    assert censored["measurement_precision_status"] == "CENSORED"
    assert censored["relation"] == "<"
    assert censored["raw_value_nm"] == "50"
    summaries = read_tsv(output_dir / "p09_ligand_endpoint_summary.tsv")
    prkdc_chembl1 = [row for row in summaries if row["gene_symbol"] == "PRKDC" and row["source_molecule_id"] == "CHEMBL1"]
    assert {row["endpoint"] for row in prkdc_chembl1} == {"Ki", "IC50"}


def test_within_and_cross_source_duplicates_keep_all_members(catalog_output):
    summary, output_dir, _output_zip, _input_zip, _input_hash = catalog_output
    groups = read_tsv(output_dir / "p09_duplicate_groups.tsv")
    members = read_tsv(output_dir / "p09_duplicate_group_members.tsv")
    pairs = read_tsv(output_dir / "p09_cross_source_duplicate_pairs.tsv")
    assert any(row["scope"] == "WITHIN_CHEMBL" and row["member_observations_n"] == "2" for row in groups)
    assert any(row["scope"] == "CROSS_SOURCE" and row["member_observations_n"] == "2" for row in groups)
    assert len(pairs) == 1
    assert {row["observation_id"] for row in members if row["scope"] == "CROSS_SOURCE"} == {
        "CHEMBL:activity:6",
        "BINDINGDB:source_row:103:endpoint:IC50",
    }
    assert summary["duplicates"]["source_records_deleted_n"] == 0


def test_priority_publication_queue_is_bounded_and_linked(catalog_output):
    summary, output_dir, _output_zip, _input_zip, _input_hash = catalog_output
    publications = read_tsv(output_dir / "p09_priority_publications.tsv")
    links = read_tsv(output_dir / "p09_priority_publication_observation_links.tsv")
    assert len(publications) <= 12
    assert summary["priority_publications_n"] == len(publications)
    assert {row["review_queue_rank"] for row in links}.issubset({row["review_queue_rank"] for row in publications})
    erbb3 = next(row for row in publications if row["gene_symbol"] == "ERBB3")
    assert erbb3["multi_component_records_n"] == "1"


def test_output_package_is_self_verifying_and_not_approved(catalog_output):
    summary, output_dir, output_zip, input_zip, input_hash = catalog_output
    assert sha_file(input_zip) == input_hash
    assert zipfile.ZipFile(output_zip).testzip() is None
    with zipfile.ZipFile(output_zip) as archive:
        manifest = json.loads(archive.read("p09_package_manifest.json"))
        for name, expected in manifest["files"].items():
            content = archive.read(name)
            assert len(content) == expected["size_bytes"]
            assert hashlib.sha256(content).hexdigest() == expected["sha256"]
    disk_summary = json.loads((output_dir / "p09_summary.json").read_text(encoding="utf-8"))
    assert disk_summary["scientific_acceptance"] == "NOT_APPROVED"
    assert disk_summary["merge_permitted"] is False
    assert disk_summary["scientific_verified_binding_to_pyz_n"] == 0
    assert disk_summary["anticancer_efficacy_score_created"] is False
    assert summary["baseline_65x13_untouched"] is True


def test_assessment_header_preserves_source_identifiers(catalog_output):
    _summary, output_dir, _output_zip, _input_zip, _input_hash = catalog_output
    with (output_dir / "p09_observation_assessment.tsv").open(encoding="utf-8", newline="") as stream:
        header = next(csv.reader(stream, delimiter="\t"))
    assert header == ASSESSMENT_FIELDS
    for required in (
        "activity_id",
        "assay_chembl_id",
        "molecule_chembl_id",
        "bindingdb_source_row",
        "bindingdb_reactant_set_id",
        "bindingdb_monomer_id",
        "raw_ki_nm",
        "raw_kd_nm",
        "raw_ic50_nm",
    ):
        assert required in header


def test_output_paths_cannot_overwrite_existing_or_input(tmp_path: Path, p08_zip: Path):
    existing = tmp_path / "existing"
    existing.mkdir()
    with pytest.raises(FileExistsError):
        build_catalog(
            p08_zip=p08_zip,
            expected_p08_sha256=sha_file(p08_zip),
            output_dir=existing,
            output_zip=tmp_path / "new.zip",
        )
    with pytest.raises(ValueError, match="cannot overwrite"):
        build_catalog(
            p08_zip=p08_zip,
            expected_p08_sha256=sha_file(p08_zip),
            output_dir=tmp_path / "new",
            output_zip=p08_zip,
        )
