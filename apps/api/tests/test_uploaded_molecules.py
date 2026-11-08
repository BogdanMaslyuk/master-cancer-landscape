from __future__ import annotations

import json
from pathlib import Path

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
