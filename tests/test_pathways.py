from __future__ import annotations

import pandas as pd

from mcl.analysis.pathways import (
    build_background_sets,
    build_candidate_recurrence,
    build_entrez_identifier_map,
    eligible_gene_set,
    normalize_gprofiler_result,
    resolve_ensembl_identifiers,
    select_top_candidates,
)


def _frame(rows):
    return pd.DataFrame(rows)


def test_eligible_gene_set_excludes_broad_and_low_sample_size():
    frame = _frame([
        {"gene_symbol": "A", "delta_gene_effect": -0.4, "broad_dependency_warning": False, "low_sample_size": False},
        {"gene_symbol": "B", "delta_gene_effect": -0.3, "broad_dependency_warning": True, "low_sample_size": False},
        {"gene_symbol": "C", "delta_gene_effect": -0.2, "broad_dependency_warning": False, "low_sample_size": True},
    ])
    assert eligible_gene_set(frame) == {"A"}


def test_select_top_candidates_ranks_negative_effects():
    frame = _frame([
        {"gene_symbol": "A", "delta_gene_effect": -0.2, "broad_dependency_warning": False, "low_sample_size": False},
        {"gene_symbol": "B", "delta_gene_effect": -0.7, "broad_dependency_warning": False, "low_sample_size": False},
        {"gene_symbol": "C", "delta_gene_effect": 0.1, "broad_dependency_warning": False, "low_sample_size": False},
    ])
    result = select_top_candidates(frame, comparison_label="X", top_n=2)
    assert result["gene_symbol"].tolist() == ["B", "A"]
    assert result["rank"].tolist() == [1, 2]


def test_candidate_recurrence_and_background_thresholds():
    long = pd.DataFrame([
        {"comparison": "A", "rank": 1, "gene_symbol": "KRAS", "delta_gene_effect": -1.0},
        {"comparison": "B", "rank": 2, "gene_symbol": "KRAS", "delta_gene_effect": -0.5},
        {"comparison": "A", "rank": 3, "gene_symbol": "X", "delta_gene_effect": -0.3},
    ])
    recurrence = build_candidate_recurrence(long)
    kras = recurrence.loc[recurrence["gene_symbol"] == "KRAS"].iloc[0]
    assert int(kras["comparisons_n"]) == 2
    assert float(kras["mean_delta_gene_effect"]) == -0.75

    backgrounds = build_background_sets(
        {"A": {"KRAS", "X", "Y"}, "B": {"KRAS", "Y", "Z"}, "C": {"Y", "Z"}},
        {"recurrent": 2, "core": 3},
    )
    assert backgrounds["recurrent"] == {"KRAS", "Y", "Z"}
    assert backgrounds["core"] == {"Y"}


def test_normalize_gprofiler_result():
    response = {
        "result": [
            {
                "source": "REAC",
                "native": "REAC:R-HSA-1",
                "name": "Example pathway",
                "description": "Example pathway",
                "p_value": 0.01,
                "significant": True,
                "query_size": 10,
                "term_size": 20,
                "intersection_size": 3,
                "effective_domain_size": 1000,
                "precision": 0.3,
                "recall": 0.15,
                "parents": [],
            }
        ]
    }
    result = normalize_gprofiler_result(response, "recurrent")
    assert result.loc[0, "analysis_set"] == "recurrent"
    assert result.loc[0, "source"] == "REAC"
    assert bool(result.loc[0, "significant"]) is True


def test_entrez_identifier_map_normalizes_numeric_values():
    frame = pd.DataFrame([
        {"gene_symbol": "KRAS", "entrez_gene_id": 3845},
        {"gene_symbol": "FIS1", "entrez_gene_id": 51024.0},
    ])
    mapping = build_entrez_identifier_map(frame)
    assert mapping == {"KRAS": "3845", "FIS1": "51024"}


def test_normalize_gprofiler_result_excludes_root_and_recovers_intersection_ids():
    response = {
        "meta": {
            "genes_metadata": {
                "query": {
                    "query_1": {
                        "ensgs": ["ENSG1", "ENSG2"],
                        "mapping": {"3845": "ENSG1", "51024": "ENSG2"},
                    }
                }
            }
        },
        "result": [
            {
                "source": "CORUM",
                "native": "CORUM:0000000",
                "name": "CORUM root",
                "p_value": 0.001,
                "significant": True,
                "query_size": 2,
                "term_size": 100,
                "intersection_size": 2,
                "effective_domain_size": 1000,
                "precision": 1.0,
                "recall": 0.02,
                "parents": [],
                "query": "query_1",
                "intersections": [["IPI"], ["IPI"]],
            },
            {
                "source": "CORUM",
                "native": "CORUM:4200",
                "name": "DNM1L-FIS1 complex",
                "p_value": 0.01,
                "significant": True,
                "query_size": 2,
                "term_size": 2,
                "intersection_size": 1,
                "effective_domain_size": 1000,
                "precision": 0.5,
                "recall": 0.5,
                "parents": ["CORUM:0000000"],
                "query": "query_1",
                "intersections": [[], ["IPI"]],
            },
        ],
    }
    result = normalize_gprofiler_result(
        response,
        "recurrent",
        excluded_term_ids={"CORUM:0000000"},
        input_id_to_symbol={"3845": "KRAS", "51024": "FIS1"},
    )
    assert result["term_id"].tolist() == ["CORUM:4200"]
    assert result.loc[0, "intersecting_input_ids_json"] == '["51024"]'
    assert result.loc[0, "intersecting_gene_symbols_json"] == '["FIS1"]'

def test_resolve_ensembl_identifiers_uses_symbol_to_disambiguate_readthrough():
    mapping = resolve_ensembl_identifiers(
        {"SYNCRIP": "10492", "KRAS": "3845"},
        [
            {"incoming": "10492", "converted": "ENSG00000271793", "name": "SYNCRIP-SNX14"},
            {"incoming": "10492", "converted": "ENSG00000135316", "name": "SYNCRIP"},
            {"incoming": "3845", "converted": "ENSG00000133703", "name": "KRAS"},
        ],
    )
    syncrip = mapping.loc[mapping["gene_symbol"] == "SYNCRIP"].iloc[0]
    assert syncrip["ensembl_gene_id"] == "ENSG00000135316"
    assert syncrip["mapping_status"] == "resolved_symbol_match"
    kras = mapping.loc[mapping["gene_symbol"] == "KRAS"].iloc[0]
    assert kras["ensembl_gene_id"] == "ENSG00000133703"


def test_resolve_ensembl_identifiers_keeps_unresolved_ambiguity_explicit():
    mapping = resolve_ensembl_identifiers(
        {"X": "123"},
        [
            {"incoming": "123", "converted": "ENSG1", "name": "A"},
            {"incoming": "123", "converted": "ENSG2", "name": "B"},
        ],
    )
    row = mapping.iloc[0]
    assert pd.isna(row["ensembl_gene_id"])
    assert row["mapping_status"] == "ambiguous"
