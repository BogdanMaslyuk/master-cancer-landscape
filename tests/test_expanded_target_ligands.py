"""P0.3 staging regression tests: numeric bounds, provenance, and pinned sources."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

from mcl.analysis.expanded_ligand_measurements import (
    cutoff_state,
    measurement,
    mirror_group_key,
    normalize_bindingdb,
    normalize_chembl,
)

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "ingest_expanded_target_ligands.py"
SPEC = importlib.util.spec_from_file_location("ingest_expanded_target_ligands", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)

REGISTRY = {"P11802": "CDK4", "P00533": "EGFR"}


def chembl_row(**overrides):
    row = {
        "activity_id": "123", "assay_chembl_id": "CHEMBL-A1",
        "target_chembl_id": "CHEMBL-T1", "uniprot_accession": "P11802",
        "organism_id": "9606", "target_type": "SINGLE PROTEIN",
        "confidence_score": "9", "relationship_type": "D",
        "assay_type": "B", "bao_label": "single protein format",
        "standard_type": "Ki", "standard_value": "350",
        "standard_units": "nM", "standard_relation": "=",
        "canonical_smiles": "CCO", "molecule_chembl_id": "CHEMBL-M1",
        "document_chembl_id": "CHEMBL-D1",
        "data_validity_comment": "", "potential_duplicate": "0",
    }
    row.update(overrides)
    return row


def normalize(row):
    return normalize_chembl(row, release="fixture", sha="0" * 64, locator="x:2", registry=REGISTRY)


@pytest.mark.parametrize("raw,unit,expected", [
    ("25", "nM", "25"),
    ("0.04", "µM", "40"),
    ("0.5", "uM", "500"),
    ("<500", "nM", "500"),
    (">20000", "nM", "20000"),
])
def test_numeric_unit_conversion_preserves_bound(raw, unit, expected):
    operator, nm, err = measurement(raw, unit)
    assert not err
    assert nm == expected
    if "<" in raw or ">" in raw:
        assert operator == raw[0]


def test_conflicting_relations_and_zero_are_not_exact():
    assert measurement("<500", "nM", ">")[2] == "conflicting_relations"
    assert measurement("0", "nM")[2] == "nonpositive_value"
    assert measurement("-100", "nM")[2] == "invalid_numeric_value"
    assert cutoff_state("<", "50000") == "censored_indeterminate"
    assert cutoff_state("<", "500") == "censored_eligible"
    assert cutoff_state(">", "20000") == "censored_indeterminate"


def test_chembl_strict_binding_is_provisional_only():
    row = normalize(chembl_row())
    assert row["quality_tier"] == "candidate_A"
    assert row["numeric_value_nm"] == "350"
    assert json.loads(row["reason_codes"]) == []
    assert row["measurement_status"] == "exact_eligible"
    assert mirror_group_key(row)


def test_functional_assay_is_not_claimed_as_physical_binding():
    row = normalize(chembl_row(standard_type="IC50"))
    assert row["quality_tier"] == "candidate_B"


def test_chembl_cell_based_and_other_species_never_get_primary_tier():
    cell = normalize(chembl_row(bao_label="cell-based format"))
    mouse = normalize(chembl_row(organism_id="10090"))
    low = normalize(chembl_row(confidence_score="8"))
    assert all(x["quality_tier"] not in {"candidate_A", "candidate_B"} for x in (cell, mouse, low))


def test_bound_not_upgraded_to_exact():
    row = normalize(chembl_row(standard_value="<500", standard_relation="<"))
    assert row["relation"] == "<"
    assert row["numeric_value_nm"] == "500"
    assert row["quality_tier"] == "C"


def test_bindingdb_method_unverified_and_mirror_key_only_if_explicit():
    row = {
        "BindingDB Reactant_set_id": "456",
        "BindingDB MonomerID": "789",
        "Ligand SMILES": "CCO",
        "UniProt (SwissProt) Primary ID of Target Chain": "P11802",
        "Target Source Organism According to Curator or DataSource": "Homo sapiens",
        "Ki (nM)": "<500",
        "Kd (nM)": "",
        "IC50 (nM)": "",
        "Article DOI": "10.1000/test",
        "ChEMBL Activity ID": "123",
        "Curation/DataSource": "ChEMBL",
    }
    parsed = normalize_bindingdb(row, endpoint="Ki", release="fixture", sha="0" * 64, locator="y:2:Ki", registry=REGISTRY)
    assert parsed["quality_tier"] == "C"
    assert parsed["measurement_status"] == "censored_eligible"
    assert mirror_group_key(parsed) == mirror_group_key(normalize(chembl_row()))
    row.pop("ChEMBL Activity ID")
    independent = normalize_bindingdb(row, endpoint="Ki", release="fixture", sha="0" * 64, locator="y:2:Ki", registry=REGISTRY)
    assert not mirror_group_key(independent)


def write_tsv(path: Path, rows: list[dict]):
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_manifest_hash_validation_and_streamed_staging(tmp_path: Path):
    chembl = tmp_path / "chembl.tsv"
    bdb = tmp_path / "bindingdb.tsv"
    write_tsv(chembl, [chembl_row()])
    write_tsv(bdb, [{
        "BindingDB Reactant_set_id": "456",
        "BindingDB MonomerID": "789",
        "Ligand SMILES": "CCO",
        "UniProt (SwissProt) Primary ID of Target Chain": "P11802",
        "Target Source Organism According to Curator or DataSource": "Homo sapiens",
        "Ki (nM)": "25", "Kd (nM)": "", "IC50 (nM)": "0.5",
    }])
    registry = tmp_path / "protein_target_universe.tsv"
    pd.DataFrame([{
        "target_gene": "CDK4", "uniprot_accession": "P11802",
        "mapping_status": "unique_reviewed_human", "organism_id": "9606",
    }]).to_csv(registry, sep="\t", index=False)
    def source(name, path, format_):
        return {
            "source": name, "path": path.name, "release": "fixture-pinned",
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "format": format_, "url": "https://example.org/fixture",
            "license": "test-only",
        }
    src = {
        "contract": "mcl-tls-v2-pinned-extracts-v1",
        "sources": [
            source("ChEMBL", chembl, "joined_activity_tsv"),
            source("BindingDB", bdb, "bindingdb_bulk_tsv"),
        ],
    }
    manifest = tmp_path / "source_manifest.json"
    manifest.write_text(json.dumps(src), encoding="utf-8")
    result = MODULE.build(manifest, registry, tmp_path / "output")
    assert result["measurement_rows"] == 3
    assert result["candidate_quality_tiers"]["candidate_A"] == 1
    assert result["status"] == "staging_not_published"
    saved = pd.read_parquet(tmp_path / "output" / "staging" / result["build_id"] / "measurements_normalized.parquet")
    assert len(saved) == 3
    assert set(saved["endpoint"]) == {"Ki", "IC50"}
    with pytest.raises(FileExistsError):
        MODULE.build(manifest, registry, tmp_path / "output")
    src["sources"][0]["sha256"] = "0" * 64
    manifest.write_text(json.dumps(src), encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        MODULE.build(manifest, registry, tmp_path / "changed")
    assert not (tmp_path / "changed").exists()


def test_pilot_hash_scope_filters_both_chembl_and_bindingdb(tmp_path: Path):
    """A 13+extra pilot selector cannot silently ingest unrelated proteins."""
    chembl = tmp_path / "chembl.tsv"
    binding = tmp_path / "binding.tsv"
    write_tsv(chembl, [
        chembl_row(activity_id="111", uniprot_accession="P11802"),
        chembl_row(activity_id="222", uniprot_accession="P00533",
                   target_chembl_id="CHEMBL-T2"),
    ])
    bdb_base = {
        "BindingDB Reactant_set_id": "501",
        "BindingDB MonomerID": "509",
        "Ligand SMILES": "CCO",
        "UniProt (SwissProt) Primary ID of Target Chain": "P11802",
        "Target Source Organism According to Curator or DataSource": "Homo sapiens",
        "Ki (nM)": "100", "Kd (nM)": "", "IC50 (nM)": "",
    }
    write_tsv(binding, [
        bdb_base,
        {**bdb_base, "BindingDB Reactant_set_id": "502",
         "UniProt (SwissProt) Primary ID of Target Chain": "P00533"},
    ])
    registry=tmp_path/"protein_target_universe.tsv"
    pd.DataFrame([
        {"target_gene":"CDK4","uniprot_accession":"P11802",
         "mapping_status":"unique_reviewed_human","organism_id":"9606"},
        {"target_gene":"EGFR","uniprot_accession":"P00533",
         "mapping_status":"unique_reviewed_human","organism_id":"9606"},
    ]).to_csv(registry,sep="\t",index=False)
    pilot=tmp_path/"pilot_panel.tsv"
    pd.DataFrame([
        {"target_gene":"CDK4","uniprot_accession":"P11802","scope":"core"}
    ]).to_csv(pilot,sep="\t",index=False)
    digest=lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=tmp_path/"source_manifest.json"
    sources=[
        {"source":name, "path":file.name,"release":"fixture",
         "sha256":digest(file),"format":format_,"url":"https://example.org",
         "license":"synthetic-only"}
        for name,file,format_ in [
            ("ChEMBL",chembl,"joined_activity_tsv"),
            ("BindingDB",binding,"bindingdb_bulk_tsv"),
        ]
    ]
    manifest.write_text(json.dumps({
        "contract":"mcl-tls-v2-pinned-extracts-v1",
        "source_extract_provenance":{"pilot_panel_sha256":digest(pilot)},
        "sources":sources,
    }))
    result=MODULE.build(manifest,registry,tmp_path/"out",pilot_panel=pilot)
    assert result["measurement_rows"]==2
    assert result["selection_mode"]=="reviewed_pilot_panel"
    assert result["pilot_panel_sha256"]==digest(pilot)
    rows=pd.read_parquet(tmp_path/"out"/"staging"/result["build_id"]/
                         "measurements_normalized.parquet")
    assert set(rows.target_gene)=={"CDK4"}
    pilot.write_text("target_gene\tuniprot_accession\tscope\nEGFR\tP00533\texpanded\n")
    with pytest.raises(ValueError,match="does not match pinned ChEMBL extract"):
        MODULE.build(manifest,registry,tmp_path/"mismatch",pilot_panel=pilot)
