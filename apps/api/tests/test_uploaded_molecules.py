from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mcl_api.own_compounds import OwnCompoundStore
from mcl_api.uploaded_molecules import UploadedMoleculeStore


def make_store(tmp_path: Path) -> UploadedMoleculeStore:
    return UploadedMoleculeStore(tmp_path, OwnCompoundStore(tmp_path))


def test_uploaded_molecule_standardization_and_descriptors(tmp_path: Path):
    store = make_store(tmp_path)
    mol, canonical, inchikey = store._standardize("CCO")

    assert canonical == "CCO"
    assert len(inchikey) == 27

    descriptors = store._descriptors(mol)
    assert descriptors["molecular_formula"] == "C2H6O"
    assert 40 < descriptors["molecular_weight"] < 50
    assert descriptors["lipinski_violations"] == 0
    assert descriptors["veber_pass"] is True


def test_uploaded_molecule_analysis_persists_and_reuses_result(tmp_path: Path):
    store = make_store(tmp_path)

    first = store.analyze("CCO", "ethanol")
    second = store.analyze("CCO", "another name")

    assert first["available"] is True
    assert first["analysis_id"].startswith("MCL-UPL-")
    assert first["targets"] == []
    assert first["experiment_panel"]["available"] is False
    assert first["reused_existing_analysis"] is False

    assert second["analysis_id"] == first["analysis_id"]
    assert second["reused_existing_analysis"] is True

    record = tmp_path / "data" / "runtime" / "uploaded_molecules" / f"{first['analysis_id']}.json"
    assert record.exists()
    payload = json.loads(record.read_text(encoding="utf-8"))
    assert payload["standardized_smiles"] == "CCO"


def test_target_ligand_space_exact_match_is_ranked_first(tmp_path: Path):
    target_dir = tmp_path / "data" / "runtime" / "target_ligand_space"
    target_dir.mkdir(parents=True)
    pd.DataFrame(
        [
            {
                "target_gene": "CDK4",
                "ligand_id": "TLS-EXACT",
                "canonical_smiles": "CCO",
                "lowest_reported_value_nm_across_endpoints": 100.0,
                "lowest_value_endpoint_type": "IC50",
                "lowest_value_source": "test",
            },
            {
                "target_gene": "CDK4",
                "ligand_id": "TLS-OTHER",
                "canonical_smiles": "c1ccccc1",
                "lowest_reported_value_nm_across_endpoints": 500.0,
                "lowest_value_endpoint_type": "IC50",
                "lowest_value_source": "test",
            },
        ]
    ).to_parquet(target_dir / "ligand_catalog.parquet", index=False)

    store = make_store(tmp_path)
    mol, _, _ = store._standardize("CCO")
    targets = store._target_search(mol)

    assert len(targets) == 1
    assert targets[0]["target_gene"] == "CDK4"
    assert targets[0]["best_ligand_id"] == "TLS-EXACT"
    assert targets[0]["best_tanimoto"] == 1.0
    assert targets[0]["support_class"] == "high_structural_transferability"


def test_force_recompute_replaces_cached_analysis(tmp_path: Path):
    store = make_store(tmp_path)
    first = store.analyze("CCO", "first")
    refreshed = store.analyze("CCO", "refreshed", force_recompute=True)

    assert refreshed["analysis_id"] == first["analysis_id"]
    assert refreshed["name"] == "refreshed"
    assert refreshed["reused_existing_analysis"] is False
