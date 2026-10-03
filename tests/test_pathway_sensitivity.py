from __future__ import annotations

import json

import pandas as pd

from mcl.analysis.pathway_sensitivity import (
    build_gene_stability,
    build_term_stability,
)


def test_build_gene_stability_tracks_recurrence_across_thresholds():
    recurrence = {
        50: pd.DataFrame([
            {"gene_symbol": "A", "comparisons_n": 2, "best_rank": 5, "mean_rank": 10.0, "mean_delta_gene_effect": -0.5},
            {"gene_symbol": "B", "comparisons_n": 1, "best_rank": 4, "mean_rank": 4.0, "mean_delta_gene_effect": -0.2},
        ]),
        100: pd.DataFrame([
            {"gene_symbol": "A", "comparisons_n": 2, "best_rank": 5, "mean_rank": 12.0, "mean_delta_gene_effect": -0.5},
            {"gene_symbol": "C", "comparisons_n": 2, "best_rank": 20, "mean_rank": 30.0, "mean_delta_gene_effect": -0.3},
        ]),
        200: pd.DataFrame([
            {"gene_symbol": "A", "comparisons_n": 3, "best_rank": 5, "mean_rank": 15.0, "mean_delta_gene_effect": -0.5},
            {"gene_symbol": "C", "comparisons_n": 2, "best_rank": 20, "mean_rank": 35.0, "mean_delta_gene_effect": -0.3},
        ]),
    }

    result = build_gene_stability(
        recurrence,
        [50, 100, 200],
        minimum_comparisons=2,
    )

    a = result.loc[result["gene_symbol"] == "A"].iloc[0]
    assert bool(a["present_all_thresholds"]) is True
    assert int(a["thresholds_n"]) == 3
    assert a["thresholds_present"] == "50;100;200"

    c = result.loc[result["gene_symbol"] == "C"].iloc[0]
    assert bool(c["recurrent_top50"]) is False
    assert bool(c["recurrent_top100"]) is True
    assert bool(c["recurrent_top200"]) is True
    assert int(c["thresholds_n"]) == 2

    assert "B" not in result["gene_symbol"].tolist()


def test_build_term_stability_tracks_significance_and_driver_genes():
    enrichment = pd.DataFrame([
        {
            "top_n": 50,
            "source": "KEGG",
            "term_id": "KEGG:1",
            "term_name": "Pathway A",
            "p_value_adjusted": 0.01,
            "significant": True,
            "intersection_size": 2,
            "intersecting_gene_symbols_json": json.dumps(["A", "B"]),
        },
        {
            "top_n": 100,
            "source": "KEGG",
            "term_id": "KEGG:1",
            "term_name": "Pathway A",
            "p_value_adjusted": 0.02,
            "significant": True,
            "intersection_size": 3,
            "intersecting_gene_symbols_json": json.dumps(["A", "B", "C"]),
        },
        {
            "top_n": 200,
            "source": "KEGG",
            "term_id": "KEGG:1",
            "term_name": "Pathway A",
            "p_value_adjusted": 0.03,
            "significant": True,
            "intersection_size": 3,
            "intersecting_gene_symbols_json": json.dumps(["A", "C", "B"]),
        },
        {
            "top_n": 100,
            "source": "CORUM",
            "term_id": "CORUM:2",
            "term_name": "Complex B",
            "p_value_adjusted": 0.04,
            "significant": True,
            "intersection_size": 2,
            "intersecting_gene_symbols_json": json.dumps(["X", "Y"]),
        },
    ])

    result = build_term_stability(enrichment, [50, 100, 200])

    stable = result.loc[result["term_id"] == "KEGG:1"].iloc[0]
    assert bool(stable["significant_all_thresholds"]) is True
    assert int(stable["significant_thresholds_n"]) == 3
    assert stable["thresholds_significant"] == "50;100;200"
    assert stable["intersecting_gene_symbols_top200"] == '["A", "B", "C"]'

    partial = result.loc[result["term_id"] == "CORUM:2"].iloc[0]
    assert bool(partial["significant_all_thresholds"]) is False
    assert int(partial["significant_thresholds_n"]) == 1
    assert partial["thresholds_significant"] == "100"
