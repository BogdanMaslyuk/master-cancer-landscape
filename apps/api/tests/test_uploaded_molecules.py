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


def _model(model_id: str, role: str, probability, *, lab=False, cancer="CANCER-001"):
    return {
        "model_id": model_id,
        "cell_line_name": model_id,
        "role": role,
        "dependency_probability": probability,
        "gene_effect": -0.8 if probability is not None else None,
        "available_in_laboratory": lab,
        "laboratory_lines": [model_id] if lab else [],
        "known_reference_compounds": [],
    }


def test_missing_crispr_never_becomes_negative_control(tmp_path: Path):
    store = make_store(tmp_path)
    store.own_compounds._contexts_for_target = lambda gene: [{
        "mcl_cancer_id": "CANCER-001", "mcl_cancer_name": "Test cancer",
        "models": [
            _model("P1", "positive", 0.91, lab=True),
            _model("N_UNMEASURED", "negative_same_cancer", None, lab=True),
            _model("N_MEASURED", "negative_same_cancer", 0.12, lab=False),
        ],
    }]
    store._non_tumor_control = lambda: {"lab_name": "BJ5ta"}
    targets = [{"target_gene": "CDK4", "best_tanimoto": 0.75}]
    cells = store._cell_lines(targets)
    unmeasured = next(c for c in cells if c["model_id"] == "N_UNMEASURED")
    assert unmeasured["role"] == "insufficient_data"
    assert unmeasured["priority_level"] == "insufficient_data"
    panel = store._experiment_panel(targets, cells)
    assert panel["negative_tumor_control"]["model_id"] == "N_MEASURED"
    assert panel["status"] == "partial"
    assert panel["complete"] is False
    assert "independent_confirmatory_model" in panel["missing_components"]


def test_missing_positive_blocks_ready_panel(tmp_path: Path):
    store = make_store(tmp_path)
    store.own_compounds._contexts_for_target = lambda gene: [{
        "mcl_cancer_id": "CANCER-001",
        "models": [_model("N1", "negative_same_cancer", 0.1)],
    }]
    targets = [{"target_gene": "CDK4", "best_tanimoto": 0.6}]
    cells = store._cell_lines(targets)
    panel = store._experiment_panel(targets, cells)
    assert panel["available"] is False
    assert panel["status"] == "not_ready"
    assert "positive_tumour_model" in panel["missing_components"]


def test_weak_similarity_does_not_generate_cell_recommendations(tmp_path: Path):
    store = make_store(tmp_path)
    store.own_compounds._contexts_for_target = lambda gene: [{
        "models": [_model("P1", "positive", 0.9, lab=True)],
    }]
    targets = [{"target_gene": "CDK4", "best_tanimoto": 0.34}]
    assert store._cell_lines(targets) == []
    panel = store._experiment_panel(targets, [])
    assert panel["available"] is False
    assert panel["complete"] is False
    assert panel["status"] == "not_ready"
    assert "supported_target_hypothesis" in panel["missing_components"]


def test_discordant_candidate_role_is_not_promoted(tmp_path: Path):
    store = make_store(tmp_path)
    store.own_compounds._contexts_for_target = lambda gene: [{
        "models": [
            _model("P_CONFLICT", "positive", 0.1, lab=True),
            _model("N_CONFLICT", "negative_same_cancer", 0.95, lab=True),
        ],
    }]
    cells = store._cell_lines([{"target_gene": "CDK4", "best_tanimoto": 0.7}])
    assert {c["role"] for c in cells} == {"discordant_or_unclassified"}
    assert all(c["priority_level"] == "exploratory" for c in cells)


def test_experiment_focus_skips_target_without_positive_model(tmp_path: Path):
    store = make_store(tmp_path)
    def contexts(gene):
        if gene == "CDK4":
            return [{"mcl_cancer_id": "CANCER-001", "models": [
                _model("N1", "negative_same_cancer", 0.1)
            ]}]
        return [{"mcl_cancer_id": "CANCER-002", "models": [
            _model("P1", "positive", 0.92, lab=True)
        ]}]
    store.own_compounds._contexts_for_target = contexts
    targets = [
        {"target_gene": "CDK4", "best_tanimoto": 0.8},
        {"target_gene": "BCL2L1", "best_tanimoto": 0.6},
    ]
    panel = store._experiment_panel(targets, store._cell_lines(targets))
    assert panel["focus_target_gene"] == "BCL2L1"
    assert panel["available"] is True
    assert panel["complete"] is False


def test_in_silico_panel_without_lab_cultures_is_not_complete(tmp_path: Path):
    store = make_store(tmp_path)
    store._non_tumor_control = lambda: {"lab_name": "BJ5ta"}
    target = [{"target_gene": "CDK4", "best_tanimoto": 0.7}]
    cells = [
        {"target_gene": "CDK4", "model_id": "P1", "role": "positive",
         "dependency_probability": 0.9, "available_in_laboratory": False,
         "cancer_id": "CANCER-001"},
        {"target_gene": "CDK4", "model_id": "P2", "role": "positive",
         "dependency_probability": 0.8, "available_in_laboratory": False,
         "cancer_id": "CANCER-001"},
        {"target_gene": "CDK4", "model_id": "N1", "role": "negative_same_cancer",
         "dependency_probability": 0.1, "available_in_laboratory": False,
         "cancer_id": "CANCER-001"},
    ]
    panel = store._experiment_panel(target, cells)
    assert panel["evidence_complete"] is True
    assert panel["laboratory_ready"] is False
    assert panel["complete"] is False
    assert panel["status"] == "partial"
    assert "laboratory_positive_tumour_model" in panel["missing_components"]


def test_mismatched_negative_cancer_is_reported(tmp_path: Path):
    store = make_store(tmp_path)
    store._non_tumor_control = lambda: {"lab_name": "BJ5ta"}
    target = [{"target_gene": "CDK4", "best_tanimoto": 0.7}]
    cells = [
        {"target_gene": "CDK4", "model_id": "P1", "role": "positive",
         "dependency_probability": 0.9, "available_in_laboratory": True,
         "cancer_id": "CANCER-001"},
        {"target_gene": "CDK4", "model_id": "P2", "role": "positive",
         "dependency_probability": 0.8, "available_in_laboratory": True,
         "cancer_id": "CANCER-001"},
        {"target_gene": "CDK4", "model_id": "N1", "role": "negative_same_cancer",
         "dependency_probability": 0.1, "available_in_laboratory": True,
         "cancer_id": "CANCER-002"},
    ]
    panel = store._experiment_panel(target, cells)
    assert panel["laboratory_ready"] is True
    assert panel["complete"] is False
    assert panel["status"] == "partial"
    assert "same_cancer_negative_control" in panel["missing_components"]


def test_full_gene_context_uses_all_depmapped_models_and_never_infers_pyz_response(tmp_path: Path, monkeypatch):
    monkeypatch.delenv("MCL_SCIENTIFIC_DATA_ROOT", raising=False)
    processed = tmp_path / "data" / "processed"
    processed.mkdir(parents=True)
    rows = []
    atlas = []
    for i in range(45):
        mid = f"ACH-TEST-{i:03d}"
        rows.append({
            "model_id": mid, "target_gene": "CDK4",
            "dependency_probability": 0.91 if i % 2 == 0 else 0.12,
            "gene_effect": -1.3 if i % 2 == 0 else -0.2,
            "mcl_cancer_id": "C-BREAST",
            "mcl_cancer_name": "Breast carcinoma",
            "mcl_organ_ru": "Молочная железа",
            "depmap_release": "26Q1",
        })
        atlas.append({"model_id": mid, "cell_line_name": f"Line-{i:03d}"})
    rows.append({
        "model_id": "ACH-UNKNOWN", "target_gene": "CDK4",
        "dependency_probability": None, "gene_effect": None,
        "mcl_cancer_id": "C-BREAST",
        "mcl_cancer_name": "Breast carcinoma",
        "mcl_organ_ru": "Молочная железа",
        "depmap_release": "26Q1",
    })
    atlas.append({"model_id": "ACH-UNKNOWN", "cell_line_name": "Unknown dependency"})
    pd.DataFrame(rows).to_parquet(processed / "depmap_model_target_context.parquet", index=False)
    pd.DataFrame(atlas).to_parquet(processed / "depmap_crispr_model_atlas.parquet", index=False)
    pd.DataFrame([{
        "model_id": "ACH-TEST-044", "lab_name": "Available but formerly outside Top-30"
    }]).to_parquet(processed / "laboratory_panel.parquet", index=False)

    store = make_store(tmp_path)
    targets = [{"target_gene": "CDK4", "best_tanimoto": 0.396}]
    cells = store._cell_lines(targets)

    assert len(cells) == 46
    assert all(c["source_layer"] == "full_processed_depmap_target_context" for c in cells)
    selected = next(c for c in cells if c["model_id"] == "ACH-TEST-044")
    assert selected["available_in_laboratory"] is True
    assert selected["physical_stock_verified"] is False
    assert selected["role"] == "genetically_dependent_candidate"
    assert selected["direct_uploaded_molecule_response"] == "absent_not_measured"
    assert set(c["role"] for c in cells) == {
        "genetically_dependent_candidate",
        "genetically_nondependent_candidate",
        "insufficient_data",
    }
    unknown = next(c for c in cells if c["model_id"] == "ACH-UNKNOWN")
    assert unknown["genetic_dependency_status"] == "not_measured"

    panel = store._experiment_panel(targets, cells)
    assert panel["panel_type"] == "exploratory_genetic_context"
    assert panel["complete"] is False
    assert panel["laboratory_ready"] is False
    assert panel["positive_models"] == []
    assert panel["negative_tumor_control"] is None
    assert panel["genetic_candidate_models_n"] == 23
    assert panel["genetic_comparator_models_n"] == 22


def test_full_gene_context_treats_missing_processed_source_as_unavailable(tmp_path: Path, monkeypatch):
    monkeypatch.delenv("MCL_SCIENTIFIC_DATA_ROOT", raising=False)
    store = make_store(tmp_path)
    context = store.own_compounds.full_target_models("CDK4")
    assert context["available"] is False
    assert context["status"] == "source_not_available"
    assert context["items"] == []
