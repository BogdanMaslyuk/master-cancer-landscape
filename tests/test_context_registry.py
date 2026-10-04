from __future__ import annotations

import pandas as pd

from mcl.analysis.context_registry import discover_context_candidates


def _models() -> pd.DataFrame:
    rows = []
    for idx in range(1, 13):
        rows.append(
            {
                "model_id": f"ACH-{idx:03d}",
                "oncotree_code": "MEL",
                "oncotree_subtype": "Melanoma",
                "oncotree_primary_disease": "Melanoma",
                "oncotree_lineage": "Skin",
                "depmap_model_type": "MEL",
            }
        )
    return pd.DataFrame(rows)


def test_discovers_ready_exact_variant_with_wildtype_and_other_variant_comparators():
    mutations = []
    for idx in range(1, 6):
        mutations.append(
            {
                "model_id": f"ACH-{idx:03d}",
                "gene": "BRAF",
                "protein_change": "p.V600E",
                "driver": True,
                "hotspot": True,
                "vep_impact": "MODERATE",
            }
        )
    for idx in range(6, 11):
        mutations.append(
            {
                "model_id": f"ACH-{idx:03d}",
                "gene": "BRAF",
                "protein_change": "p.K601E",
                "driver": True,
                "hotspot": True,
                "vep_impact": "MODERATE",
            }
        )

    models = _models()
    all_ids = set(models["model_id"])
    frame = discover_context_candidates(
        models,
        pd.DataFrame(mutations),
        crispr_model_ids=all_ids,
        sequenced_model_ids=all_ids,
        min_discovery_n=3,
        min_context_n=5,
        min_comparator_n=2,
    )

    exact = frame[
        (frame["event_type"] == "exact_protein_change")
        & (frame["alteration_gene"] == "BRAF")
        & (frame["alteration"] == "p.V600E")
    ]
    assert set(exact["comparison_mode"]) == {"vs_gene_wildtype", "vs_other_gene_variants"}

    wt = exact[exact["comparison_mode"] == "vs_gene_wildtype"].iloc[0]
    assert wt["context_models_n"] == 5
    assert wt["comparator_models_n"] == 2
    assert wt["readiness_status"] == "READY"
    assert wt["priority_tier"] == "A"

    other = exact[exact["comparison_mode"] == "vs_other_gene_variants"].iloc[0]
    assert other["context_models_n"] == 5
    assert other["comparator_models_n"] == 5
    assert other["readiness_status"] == "READY"


def test_gene_level_candidate_uses_qualifying_mutations_only():
    models = _models()
    mutations = pd.DataFrame(
        [
            {
                "model_id": f"ACH-{idx:03d}",
                "gene": "TP53",
                "protein_change": f"p.R{idx}*",
                "driver": False,
                "hotspot": False,
                "likely_lof": True,
                "vep_impact": "HIGH",
            }
            for idx in range(1, 7)
        ]
        + [
            {
                "model_id": "ACH-007",
                "gene": "TTN",
                "protein_change": "p.A1V",
                "driver": False,
                "hotspot": False,
                "likely_lof": False,
                "vep_impact": "LOW",
            }
        ]
    )
    all_ids = set(models["model_id"])

    frame = discover_context_candidates(
        models,
        mutations,
        crispr_model_ids=all_ids,
        sequenced_model_ids=all_ids,
        min_discovery_n=3,
        min_context_n=5,
        min_comparator_n=5,
    )

    tp53 = frame[
        (frame["event_type"] == "gene_altered")
        & (frame["alteration_gene"] == "TP53")
    ].iloc[0]
    assert tp53["context_models_n"] == 6
    assert tp53["comparator_models_n"] == 6
    assert tp53["readiness_status"] == "READY"
    assert tp53["priority_tier"] == "B"
    assert not (frame["alteration_gene"] == "TTN").any()


def test_low_sample_candidate_remains_exploratory():
    models = _models()
    mutations = pd.DataFrame(
        [
            {
                "model_id": f"ACH-{idx:03d}",
                "gene": "NRAS",
                "protein_change": "p.Q61R",
                "driver": True,
                "hotspot": True,
                "vep_impact": "MODERATE",
            }
            for idx in range(1, 4)
        ]
    )
    all_ids = set(models["model_id"])

    frame = discover_context_candidates(
        models,
        mutations,
        crispr_model_ids=all_ids,
        sequenced_model_ids=all_ids,
        min_discovery_n=3,
        min_context_n=5,
        min_comparator_n=5,
    )

    exact = frame[
        (frame["event_type"] == "exact_protein_change")
        & (frame["alteration_gene"] == "NRAS")
        & (frame["alteration"] == "p.Q61R")
        & (frame["comparison_mode"] == "vs_gene_wildtype")
    ].iloc[0]
    assert exact["context_models_n"] == 3
    assert exact["readiness_status"] == "EXPLORATORY_LOW_N"
    assert exact["priority_tier"] == "C"
