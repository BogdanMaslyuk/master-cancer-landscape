from __future__ import annotations

from pathlib import Path

import pandas as pd

from mcl_api.cohort import MCLModelCohortStore


def _write_tsv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_cohort_summary_tracks_group_and_data_coverage(tmp_path: Path):
    (tmp_path / "config").mkdir(parents=True)
    (tmp_path / "config/cancer_contexts.yaml").write_text(
        """contexts:
  C1:
    minimum_context_n_warning: 2
    depmap:
      context_definition: target models
      comparator_definition: comparator models
""",
        encoding="utf-8",
    )

    _write_tsv(
        tmp_path / "data/processed/depmap_context_audit.tsv",
        [
            {"cancer_id": "C1", "model_id": "ACH-1", "assigned_group": "context", "sequencing_available": True},
            {"cancer_id": "C1", "model_id": "ACH-2", "assigned_group": "comparator", "sequencing_available": True},
            {"cancer_id": "C1", "model_id": "ACH-3", "assigned_group": "excluded", "sequencing_available": False},
        ],
    )
    _write_tsv(
        tmp_path / "data/processed/depmap_model_metadata.tsv",
        [
            {"model_id": "ACH-1", "depmap_release": "test"},
            {"model_id": "ACH-2", "depmap_release": "test"},
        ],
    )
    _write_tsv(
        tmp_path / "data/processed/depmap_model_mutations.tsv",
        [
            {"model_id": "ACH-1", "gene": "KRAS"},
        ],
    )

    cohort = MCLModelCohortStore(tmp_path).summary("C1")

    assert cohort["total"]["models_n"] == 3
    assert cohort["total"]["sequencing_n"] == 2
    assert cohort["total"]["mutation_profile_n"] == 1
    assert cohort["groups"]["context"]["models_n"] == 1
    assert cohort["groups"]["comparator"]["models_n"] == 1
    assert cohort["groups"]["excluded"]["models_n"] == 1
    assert cohort["status"] == "caution"
    codes = {flag["code"] for flag in cohort["flags"]}
    assert "low_context_sample" in codes
    assert "mutation_index_incomplete" in codes
    assert "metadata_incomplete" in codes
    assert cohort["representativeness"]["status"] == "not_assessed"


def test_cohort_is_ready_when_configured_sample_rules_are_met(tmp_path: Path):
    (tmp_path / "config").mkdir(parents=True)
    (tmp_path / "config/cancer_contexts.yaml").write_text(
        """contexts:
  C1:
    minimum_context_n_warning: 1
    depmap: {}
""",
        encoding="utf-8",
    )
    _write_tsv(
        tmp_path / "data/processed/depmap_context_audit.tsv",
        [
            {"cancer_id": "C1", "model_id": "ACH-1", "assigned_group": "context", "sequencing_available": True},
            {"cancer_id": "C1", "model_id": "ACH-2", "assigned_group": "comparator", "sequencing_available": True},
        ],
    )
    _write_tsv(
        tmp_path / "data/processed/depmap_model_metadata.tsv",
        [
            {"model_id": "ACH-1"},
            {"model_id": "ACH-2"},
        ],
    )
    _write_tsv(
        tmp_path / "data/processed/depmap_model_mutations.tsv",
        [
            {"model_id": "ACH-1"},
            {"model_id": "ACH-2"},
        ],
    )

    cohort = MCLModelCohortStore(tmp_path).summary("C1")
    assert cohort["status"] == "ready"
    assert not [flag for flag in cohort["flags"] if flag["level"] == "warning"]
