from __future__ import annotations

from pathlib import Path

import pandas as pd

from mcl_api.scientific_cell_context import FullGeneModelContext
from mcl_api.scientific_reference_context import KnownDrugReferenceContext


def _make_scientific_fixture(root: Path) -> None:
    processed = root / "data" / "processed"
    pharm = root / "data" / "runtime" / "pharmacology"
    processed.mkdir(parents=True)
    pharm.mkdir(parents=True)

    pd.DataFrame([
        {
            "target_gene": "CDK4", "model_id": "M-A", "dependency_probability": 1.0,
            "gene_effect": -1.9, "depmap_release": "26Q1",
            "mcl_cancer_id": "BREAST", "mcl_cancer_name": "Breast",
        },
        {
            "target_gene": "CDK4", "model_id": "M-B", "dependency_probability": 0.2,
            "gene_effect": -0.2, "depmap_release": "26Q1",
            "mcl_cancer_id": "BREAST", "mcl_cancer_name": "Breast",
        },
        {
            "target_gene": "CDK4", "model_id": "M-C", "dependency_probability": 0.8,
            "gene_effect": -0.9, "depmap_release": "26Q1",
            "mcl_cancer_id": "COLON", "mcl_cancer_name": "Colon",
        },
    ]).to_parquet(processed / "depmap_model_target_context.parquet", index=False)
    pd.DataFrame([
        {"model_id": "M-A", "cell_line_name": "MCF-7"},
        {"model_id": "M-B", "cell_line_name": "NOT-IN-LAB"},
        {"model_id": "M-C", "cell_line_name": "HCT 116"},
    ]).to_parquet(processed / "depmap_crispr_model_atlas.parquet", index=False)
    pd.DataFrame([
        {"model_id": "M-A", "lab_name": "MCF-7", "laboratory_role": "human_tumor"},
        {"model_id": "M-C", "lab_name": "HCT 116", "laboratory_role": "human_tumor"},
        # M-X is in the laboratory register but has no measured CDK4 row.
        {"model_id": "M-X", "lab_name": "OUTSIDE-CDK4", "laboratory_role": "human_technical"},
    ]).to_parquet(processed / "laboratory_panel.parquet", index=False)

    pd.DataFrame([
        {"model_id": "M-A", "compound_id": "RIB", "target_gene": "CDK4"},
        {"model_id": "M-A", "compound_id": "RIB", "target_gene": "CDK4"},
        {"model_id": "M-A", "compound_id": "ON", "target_gene": "CDK4"},
        {"model_id": "M-B", "compound_id": "ON", "target_gene": "CDK4"},
        {"model_id": "M-X", "compound_id": "ON", "target_gene": "CDK4"},
        {"model_id": "M-A", "compound_id": "NOT-CDK4", "target_gene": "CDK5"},
    ]).to_parquet(pharm / "model_compound_target_links.parquet", index=False)
    pd.DataFrame([
        {
            "observation_id": "PRISM-RIB-1", "model_id": "M-A", "compound_id": "RIB",
            "source": "PRISM", "source_release": "24Q2", "source_assay_id": "ASSAY-X",
            "endpoint": "LFC", "value": -3.42, "dose": 2.5,
            "dose_unit": "um", "exposure_time_h": 120.0,
        },
        {
            "observation_id": "PRISM-RIB-2", "model_id": "M-A", "compound_id": "RIB",
            "source": "PRISM", "source_release": "24Q2", "source_assay_id": "ASSAY-Y",
            "endpoint": "LFC", "value": -2.7, "dose": 2.5,
            "dose_unit": "um", "exposure_time_h": 120.0,
        },
        {
            "observation_id": "PRISM-ON-1", "model_id": "M-A", "compound_id": "ON",
            "source": "PRISM", "source_release": "24Q2", "source_assay_id": "ASSAY-X",
            "endpoint": "LFC", "value": -0.869, "dose": 2.5,
            "dose_unit": "um", "exposure_time_h": 120.0,
        },
        {
            "observation_id": "PRISM-ON-B", "model_id": "M-B", "compound_id": "ON",
            "source": "PRISM", "source_release": "24Q2", "source_assay_id": "ASSAY-X",
            "endpoint": "LFC", "value": -0.1, "dose": 2.5,
            "dose_unit": "um", "exposure_time_h": 120.0,
        },
        {
            "observation_id": "PRISM-OUTSIDE-GENE", "model_id": "M-X", "compound_id": "ON",
            "source": "PRISM", "source_release": "24Q2", "source_assay_id": "ASSAY-X",
            "endpoint": "LFC", "value": -5.0, "dose": 2.5,
            "dose_unit": "um", "exposure_time_h": 120.0,
        },
        {
            "observation_id": "OTHER", "model_id": "M-A", "compound_id": "ON",
            "source": "OTHER", "source_release": "24Q2", "source_assay_id": "OTHER",
            "endpoint": "LFC", "value": -10.0, "dose": 2.5,
            "dose_unit": "um", "exposure_time_h": 120.0,
        },
    ]).to_parquet(pharm / "responses.parquet", index=False)
    pd.DataFrame([
        {"compound_id": "RIB", "preferred_name": "RIBOCICLIB"},
        {"compound_id": "ON", "preferred_name": "ON123300"},
    ]).to_parquet(pharm / "compound_catalog.parquet", index=False)
    pd.DataFrame([
        {
            "evidence_id": "E-RIB", "compound_id": "RIB", "target_gene": "CDK4",
            "evidence_type": "known_compound_annotation",
            "directness": "uncertain", "action": "inhibitor", "source": "curated",
        },
        {
            "evidence_id": "E-ON", "compound_id": "ON", "target_gene": "CDK4",
            "evidence_type": "known_compound_annotation",
            "directness": "uncertain", "action": "inhibitor", "source": "curated",
        },
    ]).to_parquet(pharm / "target_evidence.parquet", index=False)


def test_prism_known_drug_observations_stay_independent_from_pyz(tmp_path: Path, monkeypatch):
    monkeypatch.delenv("MCL_SCIENTIFIC_DATA_ROOT", raising=False)
    _make_scientific_fixture(tmp_path)
    result = KnownDrugReferenceContext(tmp_path).load("CDK4", ["M-A", "M-C"])

    assert result["available"] is True
    assert result["queried_models_n"] == 2
    assert result["observations_n"] == 3
    records = result["observations_by_model"]["M-A"]
    assert len(records) == 3
    assert {r["compound_id"] for r in records} == {"RIB", "ON"}
    assert {r["observation_id"] for r in records} == {
        "PRISM-RIB-1", "PRISM-RIB-2", "PRISM-ON-1"
    }
    assert all(r["source_release"] == "24Q2" for r in records)
    assert all(r["response_is_for_uploaded_molecule"] is False for r in records)
    assert all(r["reference_model_role"] is None for r in records)
    assert not any("M-B" == r["model_id"] for r in records)
    assert len([r for r in records if r["compound_id"] == "RIB"]) == 2
    assert len(next(r for r in records if r["compound_id"] == "RIB")["target_annotations"]) == 1
    assert "M-C" not in result["observations_by_model"]


def test_full_depmapped_context_keeps_all_models_and_only_attaches_lab_prism(tmp_path: Path, monkeypatch):
    monkeypatch.delenv("MCL_SCIENTIFIC_DATA_ROOT", raising=False)
    _make_scientific_fixture(tmp_path)
    result = FullGeneModelContext(tmp_path).load("CDK4")

    assert result["total_model_rows"] == 3
    assert result["total_measured_models"] == 3
    assert result["laboratory_registered_models_n"] == 2
    assert result["laboratory_registry_mapped_models_n"] == 3
    assert result["laboratory_registry_outside_context_n"] == 1
    assert result["laboratory_registry_outside_context_model_ids"] == ["M-X"]
    assert result["laboratory_registry_outside_context_details"] == [{
        "model_id": "M-X",
        "lab_names": ["OUTSIDE-CDK4"],
        "laboratory_roles": ["human_technical"],
        "in_full_crispr_atlas": False,
        "in_requested_gene_context": False,
        "diagnostic_status": "outside_requested_gene_processed_context",
        "interpretation_ru": (
            "Нет строки для данного гена в текущем обработанном контексте. "
            "Причина и наличие измерений по другим генам требуют отдельной сверки."
        ),
    }]
    assert result["laboratory_genetic_triage"]["available"] is True
    assert result["laboratory_genetic_triage"]["laboratory_human_tumor_model_ids_n"] == 2
    assert result["laboratory_genetic_triage"]["paired_cancer_cohorts_n"] == 0
    assert result["reference_drug_context"]["queried_models_n"] == 2
    assert result["reference_drug_context"]["observations_n"] == 3
    assert result["reference_drug_context"]["query_scope"] == (
        "laboratory_registry_intersection_with_measured_target_context"
    )
    models = {r["model_id"]: r for r in result["items"]}
    assert models["M-A"]["reference_drug_observations_n"] == 3
    assert len(models["M-A"]["known_reference_compounds"]) == 2
    assert models["M-A"]["genetic_dependency_status"] == "dependent"
    assert models["M-A"]["physical_stock_verified"] is False
    assert models["M-B"]["reference_drug_evidence_status"] == "not_queried_non_laboratory_model"
    assert models["M-B"]["reference_drug_observations"] == []
    assert models["M-C"]["reference_drug_evidence_status"] == "no_linked_prism_lfc_observations"
    assert all(r["direct_uploaded_molecule_response"] == "absent_not_measured"
               for r in result["items"])


def test_missing_known_drug_tables_are_not_misread_as_negative(tmp_path: Path, monkeypatch):
    monkeypatch.delenv("MCL_SCIENTIFIC_DATA_ROOT", raising=False)
    _make_scientific_fixture(tmp_path)
    (tmp_path / "data" / "runtime" / "pharmacology" / "target_evidence.parquet").unlink()
    result = FullGeneModelContext(tmp_path).load("CDK4")
    assert result["available"] is True
    assert result["reference_drug_context"]["available"] is False
    models = {r["model_id"]: r for r in result["items"]}
    assert models["M-A"]["reference_drug_evidence_status"] == "normalized_source_unavailable"
    assert models["M-A"]["reference_drug_observations"] == []
    assert models["M-A"]["genetic_dependency_status"] == "dependent"


def test_lab_model_without_target_data_does_not_inflate_prism_measurement_count(
    tmp_path: Path, monkeypatch
):
    monkeypatch.delenv("MCL_SCIENTIFIC_DATA_ROOT", raising=False)
    _make_scientific_fixture(tmp_path)
    raw = KnownDrugReferenceContext(tmp_path).load("CDK4", ["M-A", "M-C", "M-X"])
    assert raw["queried_models_n"] == 3
    assert raw["observations_n"] == 4
    assert len(raw["observations_by_model"]["M-X"]) == 1

    scoped = FullGeneModelContext(tmp_path).load("CDK4")
    assert scoped["reference_drug_context"]["queried_models_n"] == 2
    assert scoped["reference_drug_context"]["observations_n"] == 3
    assert len(scoped["items"]) == 3
    assert set(item["model_id"] for item in scoped["items"]) == {"M-A", "M-B", "M-C"}
    assert sum(item["reference_drug_observations_n"] for item in scoped["items"]) == 3
    assert scoped["laboratory_registry_outside_context_model_ids"] == ["M-X"]
