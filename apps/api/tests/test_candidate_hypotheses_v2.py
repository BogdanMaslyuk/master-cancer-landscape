from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mcl_api.candidate_hypotheses import CandidateHypothesisStore


def _runtime_v2(root: Path) -> None:
    runtime = root / "data" / "runtime" / "pharmacology"
    runtime.mkdir(parents=True, exist_ok=True)

    pd.DataFrame([
        {
            "hypothesis_id": "HYP2-TEST",
            "compound_id": "CMP-2",
            "preferred_name": "Example v2 inhibitor",
            "target_gene": "BRAF",
            "protein_preferred_name": "B-Raf proto-oncogene, serine/threonine kinase",
            "uniprot_primary_accession": "P15056",
            "mcl_cancer_id": "skin-melanoma",
            "mcl_cancer_name": "Melanoma",
            "mcl_organ_ru": "Кожа",
            "models_n": 5,
            "active_models_n": 3,
            "sensitive_models_n": 2,
            "dependency_measured_models_n": 5,
            "crispr_models_n": 5,
            "dependency_models_n": 3,
            "joint_support_models_n": 2,
            "active_fraction": 0.6,
            "sensitive_fraction": 0.4,
            "dependency_fraction_in_cancer": 0.6,
            "median_dependency_probability": 0.84,
            "median_gene_effect": -0.71,
            "hotspot_models_n": 2,
            "likely_lof_models_n": 0,
            "phenotype_axis": "strong",
            "dependency_axis": "strong",
            "mechanism_axis": "supportive",
            "specificity_axis": "partial_enrichment",
            "molecular_context_axis": "mutation_context_present_direction_unresolved",
            "priority_status": "priority_for_in_vitro",
            "priority_status_ru": "Приоритет для in vitro проверки",
            "priority_reasons_json": json.dumps(["absolute activity", "dependency"], ensure_ascii=False),
            "evidence_gaps_json": json.dumps(["normal tissues not assessed"], ensure_ascii=False),
        }
    ]).to_parquet(runtime / "candidate_hypotheses_v2.parquet", index=False)

    pd.DataFrame([
        {
            "hypothesis_id": "HYP2-TEST",
            "model_id": "ACH-1",
            "cell_line_name": "MODEL1",
            "mcl_cancer_name": "Melanoma",
            "model_role": "positive",
            "response_value": -1.8,
            "relative_sensitivity": 0.95,
            "absolute_active": True,
            "dependency_probability": 0.91,
            "dependency_call": True,
            "gene_effect": -0.9,
            "expression_log2_tpm1": 5.2,
            "copy_number_relative": 2.4,
            "has_hotspot": True,
            "has_likely_lof": False,
            "protein_changes_json": json.dumps(["p.V600E"]),
            "role_reason_ru": "support",
        }
    ]).to_parquet(runtime / "candidate_hypothesis_models_v2.parquet", index=False)

    (runtime / "candidate_hypotheses_v2_manifest.json").write_text(
        json.dumps(
            {
                "contract": "mcl-candidate-hypotheses-v2",
                "hypotheses_n": 1,
                "status_counts": {"priority_for_in_vitro": 1},
                "dependency_probability_threshold": 0.5,
                "prism_active_lfc_threshold": -1.0,
            },
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )


def test_store_prefers_v2_and_exposes_probability_dependency(tmp_path: Path):
    _runtime_v2(tmp_path)
    store = CandidateHypothesisStore(tmp_path)
    payload = store.search(
        search="BRAF",
        status="priority_for_in_vitro",
        target_gene=None,
        cancer=None,
        limit=10,
        offset=0,
    )
    assert payload["candidate_version"] == "v2"
    assert payload["total"] == 1
    item = payload["items"][0]
    assert item["dependency_measured_models_n"] == 5
    assert item["dependency_models_n"] == 3
    assert item["molecular_context_axis"] == "mutation_context_present_direction_unresolved"


def test_v2_detail_keeps_mutation_context_descriptive(tmp_path: Path):
    _runtime_v2(tmp_path)
    store = CandidateHypothesisStore(tmp_path)
    payload = store.detail("HYP2-TEST")
    assert payload["candidate_version"] == "v2"
    model = payload["models_by_role"]["positive"][0]
    assert model["dependency_probability"] == 0.91
    assert model["protein_changes"] == ["p.V600E"]
    assert model["has_hotspot"] is True
    assert payload["falsification_criteria"]
