from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mcl_api.laboratory import LaboratoryPanelStore


def _runtime(root: Path) -> None:
    processed = root / "data" / "processed"
    runtime = root / "data" / "runtime" / "laboratory"
    processed.mkdir(parents=True, exist_ok=True)
    runtime.mkdir(parents=True, exist_ok=True)

    pd.DataFrame([
        {
            "lab_id": "LAB-001",
            "lab_name": "A549",
            "species": "human",
            "model_class": "tumor",
            "laboratory_role": "human_tumor",
            "match_status": "matched",
            "model_id": "ACH-1",
            "mcl_cancer_name": "Lung cancer",
            "mcl_organ_ru": "Лёгкое",
            "has_crispr_atlas": True,
            "highlighted_on_source": True,
            "control_role": "",
        },
        {
            "lab_id": "LAB-002",
            "lab_name": "BJ5ta",
            "species": "human",
            "model_class": "non_tumor_immortalized",
            "laboratory_role": "human_non_tumor_control",
            "match_status": "matched",
            "model_id": "ACH-BJ",
            "mcl_cancer_name": None,
            "mcl_organ_ru": None,
            "has_crispr_atlas": False,
            "highlighted_on_source": False,
            "control_role": "general_non_tumor_control",
        },
    ]).to_parquet(processed / "laboratory_panel.parquet", index=False)
    (processed / "laboratory_panel_manifest.json").write_text(
        json.dumps({
            "laboratory_lines_n": 2,
            "human_lines_n": 2,
            "human_tumor_lines_n": 1,
            "human_non_tumor_controls_n": 1,
            "control_interpretation_ru": "BJ5ta is a non-tumor control.",
        }),
        encoding="utf-8",
    )

    pd.DataFrame([
        {
            "hypothesis_id": "HYP2-X",
            "preferred_name": "Example inhibitor",
            "compound_id": "CMP-1",
            "target_gene": "EGFR",
            "protein_preferred_name": "EGFR",
            "mcl_cancer_name": "Lung cancer",
            "mcl_organ_ru": "Лёгкое",
            "priority_status": "priority_for_in_vitro",
            "laboratory_readiness": "ready_with_internal_tumor_control",
            "laboratory_positive_models_n": 1,
            "laboratory_negative_models_n": 1,
            "joint_support_models_n": 3,
            "positive_lab_lines": "A549",
            "negative_lab_lines": "CONTROL-X",
            "bj5ta_available": True,
            "bj5ta_prism_available": False,
        }
    ]).to_parquet(runtime / "laboratory_candidate_hypotheses.parquet", index=False)
    pd.DataFrame([
        {
            "hypothesis_id": "HYP2-X",
            "lab_name": "A549",
            "model_id": "ACH-1",
            "laboratory_model_role": "positive",
        }
    ]).to_parquet(runtime / "laboratory_candidate_models.parquet", index=False)
    (runtime / "laboratory_candidates_manifest.json").write_text(
        json.dumps({
            "candidate_priorities_n": 1,
            "readiness_counts": {"ready_with_internal_tumor_control": 1},
            "bj5ta_available": True,
        }),
        encoding="utf-8",
    )


def test_laboratory_summary_and_control_are_separate(tmp_path: Path):
    _runtime(tmp_path)
    store = LaboratoryPanelStore(tmp_path)
    summary = store.summary()
    assert summary["available"] is True
    lines = store.lines(species="human")
    assert lines["total"] == 2
    control = [row for row in lines["items"] if row["lab_name"] == "BJ5ta"][0]
    assert control["laboratory_role"] == "human_non_tumor_control"
    assert control["has_crispr_atlas"] is False


def test_laboratory_candidates_prioritize_experimental_readiness(tmp_path: Path):
    _runtime(tmp_path)
    store = LaboratoryPanelStore(tmp_path)
    payload = store.candidate_list(readiness="ready_with_internal_tumor_control", limit=10, offset=0)
    assert payload["total"] == 1
    item = payload["items"][0]
    assert item["target_gene"] == "EGFR"
    assert item["laboratory_positive_models_n"] == 1
    assert item["bj5ta_available"] is True

    detail = store.candidate_detail("HYP2-X")
    assert detail["laboratory_models_by_role"]["positive"][0]["lab_name"] == "A549"
