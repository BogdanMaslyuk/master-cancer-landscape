"""P0.3.5: fixed chemical identity, assay verification, and staged FP index."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
import pytest

from mcl.analysis.expanded_ligand_chemistry import characterize_smiles, FP_BYTES

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_expanded_target_index.py"
SPEC = importlib.util.spec_from_file_location("build_expanded_target_index", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def digest(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def item(*, source="ChEMBL", source_id="1", gene="CDK4", accession="P11802",
         smiles="CCO", endpoint="Ki", assay="A1", nm="350", origin="", doi=""):
    return {
        "source": source, "source_release": "fixture-test", "source_record_id": source_id,
        "target_gene": gene, "uniprot_accession": accession, "original_smiles": smiles,
        "endpoint": endpoint, "assay_id": assay,
        "quality_tier": "candidate_A" if endpoint in ("Ki", "Kd") and source=="ChEMBL"
                        else "candidate_B" if source=="ChEMBL" else "C",
        "relation": "=", "numeric_value_nm": nm,
        "measurement_status": "exact_eligible",
        "original_chembl_activity_id": origin,
        "doi": doi,
    }


def prepared(tmp_path: Path, rows=None):
    if rows is None:
        rows = [
            item(source_id="101",gene="CDK4",smiles="CCO",endpoint="Ki",assay="A1"),
            item(source="BindingDB",source_id="222",gene="CDK4",smiles="CCO",
                 endpoint="Ki",assay="",origin="101"),
            item(source_id="102",gene="EGFR",accession="P00533",smiles="CCN",
                 endpoint="IC50",assay="A2",nm="500"),
            item(source_id="103",gene="CDK4",smiles="CCO.Cl",assay="A1"),
        ]
    ingest = tmp_path / "p03" / "staging" / "fixture"
    ingest.mkdir(parents=True)
    normalized = ingest / "measurements_normalized.parquet"
    pd.DataFrame(rows).fillna("").to_parquet(normalized, index=False)
    reg = tmp_path / "protein_target_universe.tsv"
    reg.write_text("uniprot_accession\ttarget_gene\nP11802\tCDK4\nP00533\tEGFR\n")
    source_manifest = {
        "contract": "mcl-tls-v2-normalized-staging-v1",
        "status": "staging_not_published",
        "normalized_sha256": digest(normalized),
        "registry_sha256": digest(reg),
        "measurement_rows": len(rows),
    }
    (ingest/"ingest_manifest.json").write_text(json.dumps(source_manifest))
    review = tmp_path / "assays.tsv"
    review.write_text(
        "source\tsource_release\tassay_id\tmethod_class\treviewer\tevidence_url\tverification_status\n"
        "ChEMBL\tfixture-test\tA1\tbinding_direct\tfixture-reviewer\thttps://example.org/assay-1\tapproved\n"
        "ChEMBL\tfixture-test\tA2\tbiochemical_function\tfixture-reviewer\thttps://example.org/assay-2\tapproved\n"
    )
    return ingest, reg, review


def test_exact_chemistry_keeps_stereoisomers_distinct():
    a, afp = characterize_smiles("C[C@H](O)Cl")
    b, bfp = characterize_smiles("C[C@@H](O)Cl")
    assert a["exact_structure_id"] != b["exact_structure_id"]
    assert a["parent_connectivity_id"] == b["parent_connectivity_id"]
    assert a["chemistry_status"] == b["chemistry_status"] == "verified_pilot_scope"
    assert len(afp) == len(bfp) == FP_BYTES
    assert afp == bfp  # non-chiral fingerprint is deliberately not a stereochemical identity


def test_salt_and_undefined_chirality_are_not_automatically_approved():
    salt, salt_fp = characterize_smiles("CCO.Cl")
    stereo, stereo_fp = characterize_smiles("CC(O)Cl")
    invalid, invalid_fp = characterize_smiles("this is not a molecule")
    assert salt["chemistry_status"] == "manual_review"
    assert stereo["chemistry_status"] == "manual_review"
    assert salt_fp is None and stereo_fp is None
    assert invalid["chemistry_status"] == "excluded" and invalid_fp is None


def test_verified_catalog_and_index_end_to_end(tmp_path: Path):
    ingest, registry, reviews = prepared(tmp_path)
    inputs = {p: digest(p) for p in [
        ingest/"measurements_normalized.parquet", ingest/"ingest_manifest.json",
        registry, reviews,
    ]}
    report = MODULE.build(ingest, registry, reviews, tmp_path / "runtime")
    assert report["status"] == "staging_not_published"
    assert report["measurement_rows_n"] == 4
    assert report["verified_measurements_n"] == 2
    assert report["unique_target_structure_pairs_n"] == 2
    assert report["proteins_with_verified_ligands_n"] == 2
    release = tmp_path / "runtime" / "staging" / report["build_id"]
    catalog = pd.read_parquet(release/"ligand_catalog.parquet")
    assert set(catalog.target_gene) == {"CDK4", "EGFR"}
    assert (release/"fp_index"/"fingerprints.bin").stat().st_size == 2*FP_BYTES
    assert len(pd.read_csv(release/"fp_index"/"entries.tsv",sep="\t")) == 2
    assert len(pd.read_parquet(release/"ligand_measurement_links.parquet")) == 2
    curated = pd.read_parquet(release/"measurements_curated.parquet")
    assert set(curated.curated_quality_tier) == {"A", "B", "C"}
    assert "confirmed_source_mirror" in set(curated.duplicate_status)
    salt = curated[curated.original_smiles=="CCO.Cl"].iloc[0]
    assert salt.curated_quality_tier == "C"
    assert not (tmp_path/"runtime"/"current.json").exists()
    assert not (tmp_path/"runtime"/"releases").exists()
    assert inputs == {p: digest(p) for p in inputs}
    with pytest.raises(FileExistsError):
        MODULE.build(ingest,registry,reviews,tmp_path/"runtime")


def test_no_assay_review_means_no_verified_catalog(tmp_path: Path):
    ingest, registry, reviews = prepared(tmp_path)
    reviews.write_text(
        "source\tsource_release\tassay_id\tmethod_class\treviewer\tevidence_url\tverification_status\n"
    )
    with pytest.raises(ValueError, match="No approved assay methods"):
        MODULE.build(ingest,registry,reviews,tmp_path/"runtime")
    assert not (tmp_path/"runtime"/"staging").exists()


def test_same_doi_different_measurements_remain_independent(tmp_path: Path):
    rows = [
        item(source_id="1",gene="CDK4",smiles="CCO",assay="A1",
             doi="10.1111/example",nm="100"),
        item(source_id="2",gene="CDK4",smiles="CCO",assay="A1",
             doi="10.1111/example",nm="500"),
    ]
    ingest, registry, review = prepared(tmp_path, rows=rows)
    result = MODULE.build(ingest,registry,review,tmp_path/"runtime")
    assert result["verified_measurements_n"] == 2
    assert result["unique_target_structure_pairs_n"] == 1
    release = tmp_path/"runtime"/"staging"/result["build_id"]
    row = pq.read_table(release/"ligand_catalog.parquet").to_pylist()[0]
    assert row["independent_measurements_n"] == "2"


def test_same_doi_same_value_requires_review(tmp_path: Path):
    rows = [
        item(source_id="1",gene="CDK4",smiles="CCO",assay="A1",
             doi="10.1111/example",nm="100"),
        item(source_id="2",gene="CDK4",smiles="CCO",assay="A1",
             doi="10.1111/example",nm="100"),
    ]
    ingest, registry, review = prepared(tmp_path, rows=rows)
    with pytest.raises(ValueError, match="Zero verified"):
        MODULE.build(ingest,registry,review,tmp_path/"runtime")


def test_tampered_normalized_file_blocks_staging(tmp_path: Path):
    ingest, registry, review = prepared(tmp_path)
    (ingest/"measurements_normalized.parquet").write_bytes(b"corrupted")
    with pytest.raises(ValueError, match="Normalized Parquet changed"):
        MODULE.build(ingest,registry,review,tmp_path/"runtime")
    assert not (tmp_path/"runtime"/"staging").exists()
