from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mcl_api.candidate_hypotheses import CandidateHypothesisStore


def _runtime(root: Path) -> None:
    runtime = root / "data" / "runtime" / "pharmacology"
    runtime.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([
        {
            "hypothesis_id": "HYP-TEST",
            "compound_id": "CMP-1",
            "preferred_name": "Example inhibitor",
            "target_gene": "EGFR",
            "protein_preferred_name": "Epidermal growth factor receptor",
            "uniprot_primary_accession": "P00533",
            "mcl_cancer_id": "lung-test",
            "mcl_cancer_name": "Lung cancer",
            "mcl_organ_ru": "Лёгкое",
            "models_n": 4,
            "sensitive_models_n": 2,
            "crispr_models_n": 4,
            "dependency_models_n": 2,
            "joint_support_models_n": 2,
            "sensitive_fraction": 0.5,
            "dependency_fraction_in_cancer": 0.5,
            "phenotype_axis": "strong",
            "dependency_axis": "strong",
            "mechanism_axis": "supportive",
            "specificity_axis": "partial_enrichment",
            "priority_status": "priority_for_in_vitro",
            "priority_status_ru": "Приоритет для in vitro проверки",
            "priority_reasons_json": json.dumps(["reason"], ensure_ascii=False),
            "evidence_gaps_json": json.dumps(["gap"], ensure_ascii=False),
        }
    ]).to_parquet(runtime / "candidate_hypotheses.parquet", index=False)
    pd.DataFrame([
        {
            "hypothesis_id": "HYP-TEST",
            "model_id": "ACH-1",
            "cell_line_name": "MODEL1",
            "model_role": "positive",
            "response_value": -1.0,
            "relative_sensitivity": 0.95,
            "gene_effect": -1.1,
            "role_reason_ru": "support",
        },
        {
            "hypothesis_id": "HYP-TEST",
            "model_id": "ACH-2",
            "cell_line_name": "MODEL2",
            "model_role": "negative_same_cancer",
            "response_value": 0.1,
            "relative_sensitivity": 0.1,
            "gene_effect": -0.1,
            "role_reason_ru": "control",
        },
    ]).to_parquet(runtime / "candidate_hypothesis_models.parquet", index=False)
    (runtime / "candidate_hypotheses_manifest.json").write_text(
        json.dumps(
            {
                "contract": "mcl-candidate-hypotheses-v1",
                "hypotheses_n": 1,
                "status_counts": {"priority_for_in_vitro": 1},
            },
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )


def test_hypothesis_search_keeps_separate_evidence_axes(tmp_path: Path):
    _runtime(tmp_path)
    store = CandidateHypothesisStore(tmp_path)
    payload = store.search(
        search="EGFR",
        status="priority_for_in_vitro",
        target_gene=None,
        cancer=None,
        limit=10,
        offset=0,
    )
    assert payload["total"] == 1
    item = payload["items"][0]
    assert item["phenotype_axis"] == "strong"
    assert item["dependency_axis"] == "strong"
    assert item["priority_reasons"] == ["reason"]


def test_hypothesis_detail_exposes_positive_and_negative_models(tmp_path: Path):
    _runtime(tmp_path)
    store = CandidateHypothesisStore(tmp_path)
    payload = store.detail("HYP-TEST")
    assert payload["hypothesis"]["target_gene"] == "EGFR"
    assert payload["models_by_role"]["positive"][0]["model_id"] == "ACH-1"
    assert payload["models_by_role"]["negative_same_cancer"][0]["model_id"] == "ACH-2"
    assert len(payload["laboratory_route"]) == 3
    assert payload["falsification_criteria"]
