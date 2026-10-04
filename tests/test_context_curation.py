from __future__ import annotations

import pandas as pd

from mcl.analysis.context_curation import curate_registry_candidates


def _row(**overrides):
    row = {
        "candidate_id": "CTX-MEL-BRAF-V600E",
        "oncotree_code": "MEL",
        "oncotree_subtype": "Melanoma",
        "event_type": "exact_protein_change",
        "alteration_gene": "BRAF",
        "alteration": "p.V600E",
        "comparison_mode": "vs_gene_wildtype",
        "context_models_n": 10,
        "comparator_models_n": 8,
        "driver_context_models_n": 10,
        "hotspot_context_models_n": 10,
        "high_impact_context_models_n": 0,
        "likely_lof_context_models_n": 0,
        "readiness_status": "READY",
        "already_configured": False,
    }
    row.update(overrides)
    return row


def test_driver_hotspot_missense_enters_first_review():
    curated = curate_registry_candidates(pd.DataFrame([_row()]))
    row = curated.iloc[0]
    assert row["review_class"] == "DRIVER_HOTSPOT_FIRST_REVIEW"
    assert row["comparison_role"] == "primary_comparator"
    assert row["alteration_class"] == "missense"
    assert row["automatic_promotion_allowed"] is False or row["automatic_promotion_allowed"] == False


def test_truncating_event_is_separated_for_lof_review():
    curated = curate_registry_candidates(
        pd.DataFrame(
            [
                _row(
                    candidate_id="CTX-OCSC-CDKN2A-R58TER",
                    oncotree_code="OCSC",
                    alteration_gene="CDKN2A",
                    alteration="p.R58Ter",
                    driver_context_models_n=5,
                    hotspot_context_models_n=5,
                    high_impact_context_models_n=5,
                    context_models_n=5,
                )
            ]
        )
    )
    row = curated.iloc[0]
    assert row["review_class"] == "LOF_SECOND_REVIEW"
    assert row["alteration_class"] == "truncating"
    assert bool(row["co_mutation_burden_review_required"]) is True


def test_wildtype_comparator_is_preferred_over_other_variant_comparator():
    frame = pd.DataFrame(
        [
            _row(candidate_id="CTX-WT", comparison_mode="vs_gene_wildtype"),
            _row(candidate_id="CTX-OTHER", comparison_mode="vs_other_gene_variants", comparator_models_n=12),
        ]
    )
    curated = curate_registry_candidates(frame)
    roles = dict(zip(curated["comparison_mode"], curated["comparison_role"]))
    assert roles["vs_gene_wildtype"] == "primary_comparator"
    assert roles["vs_other_gene_variants"] == "secondary_allele_specific_comparator"


def test_other_variant_comparator_becomes_primary_available_when_no_wildtype_row_exists():
    curated = curate_registry_candidates(
        pd.DataFrame([_row(comparison_mode="vs_other_gene_variants")])
    )
    assert curated.iloc[0]["comparison_role"] == "primary_available_comparator"


def test_non_ready_or_already_configured_candidates_are_excluded():
    frame = pd.DataFrame(
        [
            _row(candidate_id="READY"),
            _row(candidate_id="LOW", readiness_status="EXPLORATORY_LOW_N"),
            _row(candidate_id="EXISTING", already_configured=True),
        ]
    )
    curated = curate_registry_candidates(frame)
    assert list(curated["candidate_id"]) == ["READY"]
