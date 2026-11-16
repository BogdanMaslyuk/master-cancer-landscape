"""E1-P0-B inventory tests: synthetic links only, never assay approval."""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_pyz_e1_p0b_sources.py"
SPEC = importlib.util.spec_from_file_location("audit_pyz_e1_p0b_sources", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def _fixtures():
    hits = pd.DataFrame([
        {"own_compound_id": "PYZ-035", "target_gene": "CDK4", "rank_in_target": 1,
         "ligand_id": "CDK4-CHEM", "ligand_inchikey": "IK-CHEM",
         "tanimoto_morgan_r2_2048": 0.36},
        {"own_compound_id": "PYZ-055", "target_gene": "BCL2L1", "rank_in_target": 1,
         "ligand_id": "BCL2L1-BIND", "ligand_inchikey": "IK-BIND",
         "tanimoto_morgan_r2_2048": 0.41},
    ])
    catalog = pd.DataFrame([
        {"target_gene": "CDK4", "ligand_id": "CDK4-CHEM", "inchikey": "IK-CHEM"},
        {"target_gene": "BCL2L1", "ligand_id": "BCL2L1-BIND", "inchikey": "IK-BIND"},
    ])
    measurements = pd.DataFrame([
        {"target_gene": "CDK4", "inchikey": "IK-CHEM", "source": "ChEMBL",
         "source_ligand_id": "CHEMBL123", "source_record_id": "456",
         "endpoint_type": "IC50", "relation": "=", "value_nm": 175,
         "assay_id": "CHEMBL_ASSAY", "document_id": "CHEMBL_DOC"},
        {"target_gene": "CDK4", "inchikey": "IK-CHEM", "source": "ChEMBL",
         "source_ligand_id": "CHEMBL123", "source_record_id": "457",
         "endpoint_type": "Kd", "relation": "=", "value_nm": 1000,
         "assay_id": "CHEMBL_ASSAY2", "document_id": "CHEMBL_DOC"},
        {"target_gene": "BCL2L1", "inchikey": "IK-BIND", "source": "BindingDB",
         "source_ligand_id": "987", "source_record_id": "987",
         "endpoint_type": "Ki", "relation": "<", "value_nm": 2500,
         "assay_id": "", "document_id": ""},
        # Correct chemical key but wrong target; MUST NOT join to CDK4.
        {"target_gene": "EGFR", "inchikey": "IK-CHEM", "source": "BindingDB",
         "source_ligand_id": "fake", "source_record_id": "99",
         "endpoint_type": "Ki", "relation": "=", "value_nm": 100,
         "assay_id": "", "document_id": ""},
    ])
    return hits, catalog, measurements


def test_exact_gene_and_inchikey_join_keeps_endpoints_and_gaps():
    hits, catalog, measurements = _fixtures()
    n, e, stats = MODULE.prepare(hits, catalog, measurements)
    assert len(n) == 2
    assert len(e) == 3
    assert stats["assay_id_present_rows_n"] == 2
    assert stats["assay_id_missing_rows_n"] == 1
    assert stats["primary_source_approved_n"] == 0
    by_key = {(x["target_gene"], x["endpoint_type"]): x for x in e}
    assert set(by_key) == {("CDK4", "IC50"), ("CDK4", "Kd"), ("BCL2L1", "Ki")}
    assert by_key[("BCL2L1", "Ki")]["relation"] == "<"
    assert by_key[("BCL2L1", "Ki")]["provenance_status"] == "ASSAY_ID_MISSING"
    assert by_key[("CDK4", "IC50")]["provenance_status"] == "SOURCE_ASSAY_ID_UNREVIEWED"
    assert all(x["measurement_owner"] == "known_reference_ligand_not_pyz" for x in e)
    assert all(x["scientific_review_status"] == "PENDING_PRIMARY_SOURCE_VERIFICATION" for x in e)


def test_structure_mismatch_blocks_mapping():
    hits, catalog, measurements = _fixtures()
    hits.loc[hits["own_compound_id"].eq("PYZ-035"), "ligand_inchikey"] = "WRONG-KEY"
    n, e, _ = MODULE.prepare(hits, catalog, measurements)
    assert n[0]["link_status"] == "INCHIKEY_MISMATCH_NO_LINK"
    assert len(e) == 1
    assert e[0]["target_gene"] == "BCL2L1"


def test_ambiguous_catalog_ligand_does_not_choose_first():
    hits, catalog, measurements = _fixtures()
    catalog = pd.concat([catalog, catalog.iloc[[0]]], ignore_index=True)
    n, e, _ = MODULE.prepare(hits, catalog, measurements)
    assert n[0]["link_status"] == "AMBIGUOUS_CATALOG_KEY"
    assert len(e) == 1


def test_missing_data_fail_closed():
    hits, _, _ = _fixtures()
    n, e, stats = MODULE.prepare(hits, None, None)
    assert len(n) == 2
    assert e == []
    assert stats["measurement_candidates_n"] == 0
    assert {x["link_status"] for x in n} == {"CATALOG_UNAVAILABLE"}


def test_invalid_measurement_schema_rejected():
    hits, catalog, measurements = _fixtures()
    with pytest.raises(ValueError, match="missing required columns"):
        MODULE.prepare(hits, catalog, measurements.drop(columns=["inchikey"]))


def test_inventory_cli_never_touches_original_sources(tmp_path):
    science = tmp_path / "science"
    ligand = tmp_path / "ligands"
    output = tmp_path / "reports" / "audit"
    top_path = science / "data" / "runtime" / "own_compounds" / "target_ligand_space_similarity_top_hits.parquet"
    catalog_path = ligand / "data" / "runtime" / "target_ligand_space" / "ligand_catalog.parquet"
    measurement_path = ligand / "data" / "runtime" / "target_ligand_space" / "ligand_measurements.parquet"
    for path in (top_path, catalog_path, measurement_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    hits, catalog, measurements = _fixtures()
    hits.to_parquet(top_path, index=False)
    catalog.to_parquet(catalog_path, index=False)
    measurements.to_parquet(measurement_path, index=False)
    originals = [(x, x.read_bytes()) for x in (top_path, catalog_path, measurement_path)]
    command = [
        sys.executable, str(SCRIPT),
        "--science-root", str(science),
        "--ligand-root", str(ligand),
        "--out-dir", str(output),
    ]
    result = subprocess.run(command, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    report = json.loads((output / "inventory_report.json").read_text(encoding="utf-8"))
    assert report["cross_worktree_ligand_source"] is True
    assert report["scientific_status"] == "NO_PRIMARY_ASSAYS_APPROVED"
    assert report["stats"]["measurement_candidates_n"] == 3
    assert report["stats"]["assay_id_missing_rows_n"] == 1
    assert "DATA_GAPS" in report["inventory_status"]
    assert (output / "neighbors_audit.tsv").is_file()
    assert (output / "measurement_source_candidates.tsv").is_file()
    assert (output / "local_source_inventory.tsv").is_file()
    assert all(path.read_bytes() == data for path, data in originals)
    again = subprocess.run(command, text=True, capture_output=True)
    assert again.returncode != 0
    assert "already exists" in again.stderr


def test_discovery_avoids_browser_cache_and_preserves_relevant_archive(tmp_path):
    root = tmp_path / "source"
    original = root / "data" / "external" / "bindingdb" / "BindingDB_Assays.zip"
    ignored = root / "outputs" / "ui" / "browser_profiles" / "BindingDB_Assays.zip"
    original.parent.mkdir(parents=True)
    ignored.parent.mkdir(parents=True)
    original.write_bytes(b"synthetic")
    ignored.write_bytes(b"browser-cache")
    found = MODULE.discover_sources([root])
    assert any(x["path"] == str(original.resolve()) for x in found)
    assert all(x["path"] != str(ignored.resolve()) for x in found)
