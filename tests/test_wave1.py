from __future__ import annotations

import pandas as pd

from mcl.analysis.wave1 import materialize_wave1_cohorts


def _models() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "model_id": f"ACH-{idx:03d}",
                "cell_line_name": f"CL{idx}",
                "oncotree_code": "MEL",
                "oncotree_subtype": "Melanoma",
            }
            for idx in range(1, 11)
        ]
    )


def test_wave1_materializes_operational_wildtype_comparator():
    models = _models()
    mutations = pd.DataFrame(
        [
            {
                "model_id": f"ACH-{idx:03d}",
                "gene": "BRAF",
                "protein_change": "p.V600E",
                "driver": True,
                "hotspot": True,
                "vep_impact": "MODERATE",
            }
            for idx in range(1, 5)
        ]
        + [
            {
                "model_id": "ACH-005",
                "gene": "BRAF",
                "protein_change": "p.K601E",
                "driver": True,
                "hotspot": True,
                "vep_impact": "MODERATE",
            },
            {
                "model_id": "ACH-006",
                "gene": "BRAF",
                "protein_change": "p.A1V",
                "driver": False,
                "hotspot": False,
                "vep_impact": "LOW",
            },
        ]
    )
    ids = set(models["model_id"])
    contexts = [
        {
            "wave1_id": "W1-TEST",
            "label": "Melanoma — BRAF p.V600E",
            "oncotree_code": "MEL",
            "gene": "BRAF",
            "protein_change": "p.V600E",
            "comparison_mode": "vs_gene_wildtype",
            "expected_context_n": 4,
            "expected_comparator_n": 5,
        }
    ]

    cohort, summary = materialize_wave1_cohorts(
        models,
        mutations,
        contexts,
        crispr_model_ids=ids,
        sequenced_model_ids=ids,
    )

    assert (cohort["analysis_group"] == "context").sum() == 4
    assert (cohort["analysis_group"] == "comparator").sum() == 5
    assert (cohort["analysis_group"] == "excluded").sum() == 1
    # Low-impact BRAF p.A1V is not a qualifying registry alteration and therefore
    # remains in the operational comparator by design.
    row6 = cohort[cohort["model_id"] == "ACH-006"].iloc[0]
    assert row6["analysis_group"] == "comparator"
    assert bool(summary.iloc[0]["counts_match_discovery"])
    assert bool(summary.iloc[0]["comparator_is_operational_proxy"])


def test_wave1_materializes_other_variant_comparator():
    models = _models()
    mutations = pd.DataFrame(
        [
            {
                "model_id": f"ACH-{idx:03d}",
                "gene": "KRAS",
                "protein_change": "p.G12V",
                "driver": True,
                "hotspot": True,
            }
            for idx in range(1, 4)
        ]
        + [
            {
                "model_id": f"ACH-{idx:03d}",
                "gene": "KRAS",
                "protein_change": "p.G12D",
                "driver": True,
                "hotspot": True,
            }
            for idx in range(4, 7)
        ]
    )
    ids = set(models["model_id"])
    contexts = [
        {
            "wave1_id": "W1-ALLELE",
            "label": "Test",
            "oncotree_code": "MEL",
            "gene": "KRAS",
            "protein_change": "p.G12V",
            "comparison_mode": "vs_other_gene_variants",
            "expected_context_n": 3,
            "expected_comparator_n": 3,
        }
    ]

    cohort, summary = materialize_wave1_cohorts(
        models,
        mutations,
        contexts,
        crispr_model_ids=ids,
        sequenced_model_ids=ids,
    )

    assert (cohort["analysis_group"] == "context").sum() == 3
    assert (cohort["analysis_group"] == "comparator").sum() == 3
    assert (cohort["analysis_group"] == "excluded").sum() == 4
    assert bool(summary.iloc[0]["counts_match_discovery"])
    assert not bool(summary.iloc[0]["comparator_is_operational_proxy"])
