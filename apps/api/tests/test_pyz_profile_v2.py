"""Scientific safeguards and strict additive PYZ profile API contract."""
from __future__ import annotations

from copy import deepcopy

import pytest

from mcl_api.pyz_profile_v2 import make_profile
from mcl_api.schemas.own_compounds import PyzDetailResponse, PyzProfileV2


def detail(*, similarity=0.42, with_context=True, with_admet=True):
    return {
        "available": True,
        "compound": {
            "own_compound_id": "PYZ-035",
            "standardized_smiles": "CCO",
            "chemotype_ru": "антипиринат",
            "ring_variant": "Pir",
            "series_range": "PYZ-031…PYZ-035",
            "qed": 0.54 if with_admet else None,
            "integrated_risk": "signal_requires_verification" if with_admet else None,
            "admet_note_ru": "Требует проверки." if with_admet else None,
            "source_report": "SYL-RPT-2026-022 v1.6" if with_admet else None,
        },
        "admet_summary": {
            "source_report": "SYL-RPT-2026-022 v1.6",
            "ames_positive_n": 48,
            "dili_positive_n": 55,
        },
        "target_evidence": [{
            "target_gene": "CDK4",
            "best_tanimoto_morgan_r2_2048": similarity,
            "ligands_compared_n": 12,
            "best_ligand_id": "CHEMBL-REF",
            "best_lowest_reported_value_nm_across_endpoints": 18.0,
            "best_lowest_value_endpoint_type": "IC50",
            "top_neighbors": [{
                "ligand_id": "CHEMBL-REF",
                "tanimoto_morgan_r2_2048": similarity,
                "lowest_reported_value_nm_across_endpoints": 18.0,
                "lowest_value_endpoint_type": "IC50",
            }],
        }],
        "priority_target_genes": ["CDK4"] if similarity is not None and similarity >= .35 else [],
        "relevant_contexts": [{
            "target_gene": "CDK4",
            "mcl_organ_ru": "Молочная железа",
            "mcl_cancer_name": "Пример опухолевого контекста",
            "models": [{
                "model_id": "ACH-TEST",
                "cell_line_name": "TEST-LINE",
                "role": "positive",
                "gene_effect": -0.7,
                "dependency_probability": 0.8,
                "available_in_laboratory": True,
                "known_reference_compounds": ["KnownDrug"],
            }],
        }] if with_context else [],
        "docking_status": "pending_external_results",
    }


def test_profile_preserves_provenance_and_never_claims_pyz_activity():
    payload = detail()
    profile = PyzProfileV2.model_validate(make_profile(payload))
    assert profile.schema_version == "mcl-pyz-profile-v2-p0"
    assert profile.identity.standardized_smiles == "CCO"
    assert profile.target_hypotheses.scope == "existing_13_core_targets"
    target = profile.target_hypotheses.items[0]
    assert target.reference_value_nm == 18
    assert target.reference_measurement_owner == "known_ligand_not_pyz"
    assert target.passes_navigation_cutoff
    assert target.scientific_status == "exploratory_not_directly_validated"
    assert profile.scientific_summary.own_cellular_activity == "not_measured"
    assert profile.scientific_summary.own_target_binding == "not_measured"
    assert profile.provenance.bindingdb_p036_status == "NOT_INGESTED_UNREVIEWED"
    # Nested strict schema and original response fields coexist.
    validated = PyzDetailResponse.model_validate({**payload, "profile_v2": profile})
    assert validated.compound["own_compound_id"] == "PYZ-035"
    assert validated.profile_v2.compound_id == "PYZ-035"


@pytest.mark.parametrize("similarity", [0.0, 0.12, 0.349, None])
def test_low_similarity_is_missing_support_not_inactivity(similarity):
    profile = PyzProfileV2.model_validate(make_profile(detail(similarity=similarity)))
    assert not profile.target_hypotheses.items[0].passes_navigation_cutoff
    assert profile.scientific_summary.hypotheses == []
    assert "P0-NO-STRUCTURAL-SUPPORT" in profile.scientific_summary.rule_ids
    assert not any("неактивно" in x.title.lower() for x in profile.scientific_summary.gaps)


def test_admet_aggregate_does_not_become_individual_prediction():
    profile = PyzProfileV2.model_validate(make_profile(detail(with_admet=False)))
    assert profile.properties.qed is None
    assert profile.admet.integrated_risk is None
    assert profile.admet.individual_ames is None
    assert profile.admet.individual_dili is None
    assert profile.admet.status == "not_available"
    assert "P0-ADMET-MISSING" in profile.scientific_summary.rule_ids


def test_missing_target_and_cell_data_do_not_hide_compound():
    payload = detail(with_context=False)
    payload["target_evidence"] = []
    profile = PyzProfileV2.model_validate(make_profile(payload))
    assert profile.identity.data_available
    assert not profile.target_hypotheses.data_available
    assert profile.cell_context.status == "not_available"
    assert "P0-NO-CONTEXT" in profile.scientific_summary.rule_ids
    assert profile.scientific_summary.next_step.rule_id == "P0-EXPLORE-INDEPENDENTLY"


def test_invalid_confidence_does_not_create_priority():
    payload = detail()
    payload["target_evidence"][0]["best_tanimoto_morgan_r2_2048"] = "nan"
    profile = PyzProfileV2.model_validate(make_profile(payload))
    assert profile.target_hypotheses.items[0].tanimoto is None
    assert not profile.target_hypotheses.items[0].passes_navigation_cutoff


def test_profile_model_forbids_unexpected_new_fields():
    data = make_profile(detail())
    data = deepcopy(data)
    data["invented_probability_of_pyz_binding"] = 0.96
    with pytest.raises(ValueError):
        PyzProfileV2.model_validate(data)


def test_all_65_registered_pyz_have_versioned_card_without_runtime_library():
    """Use the real immutable registry; optional scientific runtimes may be absent."""
    from pathlib import Path
    from mcl_api.own_compounds import OwnCompoundStore

    root = Path(__file__).resolve().parents[3]
    store = OwnCompoundStore(root)
    frame = store.registry()
    assert frame["own_compound_id"].tolist() == [
        "PYZ-" + str(i).zfill(3) for i in range(1, 66)
    ]
    for compound_id in frame["own_compound_id"]:
        payload = store.detail(compound_id)
        result = PyzDetailResponse.model_validate(payload)
        profile = result.profile_v2
        assert profile.compound_id == compound_id
        assert profile.identity.standardized_smiles
        assert profile.scientific_summary.own_target_binding == "not_measured"
        assert profile.provenance.bindingdb_p036_status == "NOT_INGESTED_UNREVIEWED"
