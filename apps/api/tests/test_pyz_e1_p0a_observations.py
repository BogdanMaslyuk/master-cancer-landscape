"""E1-P0-A: preserve per-reference-drug candidate observations without claiming PYZ activity."""
from __future__ import annotations

import pandas as pd

from mcl_api.own_compounds import OwnCompoundStore
from mcl_api.pyz_profile_v2 import make_profile
from mcl_api.schemas.own_compounds import PyzProfileV2


def _candidate_rows():
    hypotheses = pd.DataFrame([
        {
            "hypothesis_id": "HYP-DRUG-A",
            "compound_id": "REF-A",
            "preferred_name": "Known drug A",
            "target_gene": "BCL2L1",
            "mcl_organ_ru": "Молочная железа",
            "mcl_cancer_id": "BREAST-TEST",
            "mcl_cancer_name": "Synthetic breast model",
            "priority_status": "supported_hypothesis",
        },
        {
            "hypothesis_id": "HYP-DRUG-B",
            "compound_id": "REF-B",
            "preferred_name": "Known drug B",
            "target_gene": "BCL2L1",
            "mcl_organ_ru": "Молочная железа",
            "mcl_cancer_id": "BREAST-TEST",
            "mcl_cancer_name": "Synthetic breast model",
            "priority_status": "exploratory_hypothesis",
        },
    ])
    models = pd.DataFrame([
        {
            "hypothesis_id": "HYP-DRUG-A",
            "model_id": "ACH-TEST",
            "cell_line_name": "SYNTHETIC-LINE",
            "model_role": "positive",
            "response_value": -2.1,
            "dependency_probability": 0.8,
            "gene_effect": -1.1,
        },
        {
            "hypothesis_id": "HYP-DRUG-B",
            "model_id": "ACH-TEST",
            "cell_line_name": "SYNTHETIC-LINE",
            "model_role": "negative_same_cancer",
            "response_value": 0.2,
            "dependency_probability": 0.8,
            "gene_effect": -1.1,
        },
    ])
    return hypotheses, models


def test_candidate_store_keeps_two_drugs_with_opposite_roles(tmp_path):
    hypotheses, models = _candidate_rows()
    store = OwnCompoundStore(tmp_path)
    store.candidate_hypotheses = lambda: hypotheses
    store.candidate_models = lambda: models
    store.lab_panel = lambda: pd.DataFrame()

    contexts = store._contexts_for_target("BCL2L1")
    assert len(contexts) == 1
    assert len(contexts[0]["models"]) == 1
    model = contexts[0]["models"][0]
    assert model["role"] == "positive"  # Legacy aggregate only.
    assert model["legacy_candidate_role"] == "positive"
    assert model["roles"] == ["negative_same_cancer", "positive"]

    observations = model["reference_drug_observations"]
    assert len(observations) == 2
    by_drug = {row["compound_id"]: row for row in observations}
    assert set(by_drug) == {"REF-A", "REF-B"}
    assert by_drug["REF-A"]["reference_model_role"] == "positive"
    assert by_drug["REF-B"]["reference_model_role"] == "negative_same_cancer"
    assert by_drug["REF-A"]["prism_lfc"] == -2.1
    assert by_drug["REF-B"]["prism_lfc"] == 0.2
    for row in observations:
        assert row["hypothesis_id"] in {"HYP-DRUG-A", "HYP-DRUG-B"}
        assert row["source_layer"] == "candidate_v2"
        assert row["source_scope"] == "candidate_v2_top_five_per_role"


def test_profile_keeps_each_drug_and_distinguishes_genetics_from_prism(tmp_path):
    hypotheses, models = _candidate_rows()
    store = OwnCompoundStore(tmp_path)
    store.candidate_hypotheses = lambda: hypotheses
    store.candidate_models = lambda: models
    store.lab_panel = lambda: pd.DataFrame()
    contexts = store._contexts_for_target("BCL2L1")

    payload = {
        "compound": {
            "own_compound_id": "PYZ-055",
            "standardized_smiles": "CCO",
            "qed": 0.53,
        },
        "target_evidence": [{
            "target_gene": "BCL2L1",
            "best_tanimoto_morgan_r2_2048": 0.4118,
        }],
        "relevant_contexts": contexts,
    }
    profile = PyzProfileV2.model_validate(make_profile(payload))
    cell = profile.cell_context.items[0].models[0]

    assert cell.role == "positive"  # Compatibility field.
    assert cell.legacy_aggregated_role == "positive"
    assert cell.dependency_probability == 0.8  # Separate CRISPR context.
    assert cell.gene_effect == -1.1
    assert len(cell.reference_drug_observations) == 2
    by_drug = {x.compound_id: x for x in cell.reference_drug_observations}
    assert by_drug["REF-A"].reference_model_role == "positive"
    assert by_drug["REF-B"].reference_model_role == "negative_same_cancer"
    assert by_drug["REF-A"].prism_lfc == -2.1
    assert by_drug["REF-B"].prism_lfc == 0.2
    assert by_drug["REF-A"].candidate_priority_status == "supported_hypothesis"
    assert by_drug["REF-B"].candidate_priority_status == "exploratory_hypothesis"
    assert by_drug["REF-A"].measurement_owner == "known_reference_compound_not_pyz"
    assert by_drug["REF-B"].source_scope == "candidate_v2_top_five_per_role"
    assert profile.scientific_summary.own_cellular_activity == "not_measured"
    assert profile.scientific_summary.own_target_binding == "not_measured"


def test_missing_observations_do_not_become_measured_activity():
    payload = {
        "compound": {"own_compound_id": "PYZ-001", "standardized_smiles": "CCO"},
        "target_evidence": [],
        "relevant_contexts": [{
            "target_gene": "BCL2L1",
            "models": [{
                "model_id": "ACH-TEST",
                "role": "positive",
                "dependency_probability": 0.8,
                "gene_effect": -1.1,
                "available_in_laboratory": False,
                "known_reference_compounds": ["Drug-without-detailed-rows"],
            }],
        }],
    }
    profile = PyzProfileV2.model_validate(make_profile(payload))
    cell = profile.cell_context.items[0].models[0]
    assert cell.reference_drug_observations == []
    assert cell.legacy_aggregated_role == "positive"
    assert profile.scientific_summary.own_cellular_activity == "not_measured"


def test_nonfinite_prism_value_is_not_displayed_as_valid_measurement():
    payload = {
        "compound": {"own_compound_id": "PYZ-001", "standardized_smiles": "CCO"},
        "target_evidence": [],
        "relevant_contexts": [{
            "target_gene": "CDK4",
            "models": [{
                "model_id": "ACH-TEST",
                "role": "unclassified",
                "available_in_laboratory": False,
                "reference_drug_observations": [{
                    "compound_id": "REF-A",
                    "reference_model_role": "positive",
                    "prism_lfc": float("nan"),
                    "source_layer": "candidate_v2",
                    "source_scope": "candidate_v2_top_five_per_role",
                }],
            }],
        }],
    }
    profile = PyzProfileV2.model_validate(make_profile(payload))
    observation = profile.cell_context.items[0].models[0].reference_drug_observations[0]
    assert observation.prism_lfc is None
    assert observation.measurement_owner == "known_reference_compound_not_pyz"
