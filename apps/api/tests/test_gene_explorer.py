from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mcl_api.gene_explorer import GeneExplorerStore
from mcl_api.store import MCLDataStore


def _fixture_root(tmp_path: Path) -> Path:
    (tmp_path / "config").mkdir(parents=True)
    (tmp_path / "data/processed/depmap_genomewide").mkdir(parents=True)
    (tmp_path / "data/processed/pathways_sensitivity").mkdir(parents=True)
    (tmp_path / "outputs/reports").mkdir(parents=True)
    (tmp_path / "outputs/qc").mkdir(parents=True)

    (tmp_path / "config/pathways.yaml").write_text(
        """organism: hsapiens
comparisons:
  - label: TEST_AHR_vs_control
    cancer_id: CANCER-X
    comparison: test
""",
        encoding="utf-8",
    )
    (tmp_path / "config/cancer_contexts.yaml").write_text(
        """contexts:
  CANCER-X:
    name: Test cancer — molecular context
""",
        encoding="utf-8",
    )
    (tmp_path / "outputs/reports/depmap_genomewide_CANCER-X__test_meta.json").write_text(
        json.dumps(
            {
                "context_models_n": 2,
                "comparator_models_n": 2,
                "genes_analyzed_n": 2,
                "depmap_release": "26Q1",
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "outputs/qc/depmap_genomewide_CANCER-X__test_qc.json").write_text("[]", encoding="utf-8")

    pd.DataFrame(
        [
            {
                "gene_symbol": "AHR",
                "context_median_gene_effect": -0.65,
                "comparator_median_gene_effect": -0.20,
                "delta_gene_effect": -0.45,
                "cliffs_delta": -0.60,
                "p_value": 0.001,
                "q_value": 0.01,
                "fdr_0_05": True,
                "context_models_n": 2,
                "comparator_models_n": 2,
                "broad_dependency_warning": False,
                "low_sample_size": False,
            },
            {
                "gene_symbol": "KRAS",
                "context_median_gene_effect": -1.2,
                "comparator_median_gene_effect": -0.4,
                "delta_gene_effect": -0.8,
                "cliffs_delta": -0.8,
                "p_value": 0.0001,
                "q_value": 0.001,
                "fdr_0_05": True,
                "context_models_n": 2,
                "comparator_models_n": 2,
                "broad_dependency_warning": False,
                "low_sample_size": False,
            },
        ]
    ).to_csv(tmp_path / "data/processed/depmap_genomewide/CANCER-X__test_genes.tsv", sep="\t", index=False)

    # AHR is deliberately absent from stability: arbitrary genes must still be searchable.
    pd.DataFrame(
        [
            {
                "gene_symbol": "KRAS",
                "recurrent_top50": True,
                "recurrent_top100": True,
                "recurrent_top200": True,
                "thresholds_n": 3,
                "present_all_thresholds": True,
                "thresholds_present": "50;100;200",
            }
        ]
    ).to_csv(tmp_path / "data/processed/pathways_sensitivity/gene_stability.tsv", sep="\t", index=False)

    for layer, values in {
        "gene_effect": [("ACH-1", -0.9, -1.1), ("ACH-2", -0.3, -0.8)],
        "expression": [("ACH-1", 6.2, 5.1), ("ACH-2", 4.8, 6.0)],
        "copy_number": [("ACH-1", 1.05, 0.98), ("ACH-2", 0.97, 1.02)],
    }.items():
        pd.DataFrame(values, columns=["model_id", "AHR", "KRAS"]).to_parquet(
            tmp_path / f"data/processed/depmap_model_{layer}.parquet", index=False
        )
        (tmp_path / f"data/processed/depmap_model_{layer}_genes.json").write_text(
            json.dumps([{"gene": "AHR"}, {"gene": "KRAS"}]), encoding="utf-8"
        )

    pd.DataFrame(
        [
            {"model_id": "ACH-1", "cell_line_name": "MODEL1", "oncotree_subtype": "Test tumor"},
            {"model_id": "ACH-2", "cell_line_name": "MODEL2", "oncotree_subtype": "Test tumor"},
        ]
    ).to_parquet(tmp_path / "data/processed/depmap_model_metadata.parquet", index=False)
    pd.DataFrame(
        [
            {"model_id": "ACH-1", "cancer_id": "CANCER-X", "assigned_group": "context"},
            {"model_id": "ACH-2", "cancer_id": "CANCER-X", "assigned_group": "comparator"},
        ]
    ).to_parquet(tmp_path / "data/processed/depmap_context_audit.parquet", index=False)

    return tmp_path


def test_arbitrary_gene_is_searchable_outside_stable_topn(tmp_path: Path) -> None:
    root = _fixture_root(tmp_path)
    explorer = GeneExplorerStore(root, MCLDataStore(root))

    suggestions = explorer.suggest("AH")
    assert suggestions[0]["gene_symbol"] == "AHR"

    result = explorer.search(query="AHR", delta_gene_effect_max=-0.2, q_value_max=0.05)
    assert result["total"] == 1
    assert result["items"][0]["gene_symbol"] == "AHR"
    assert result["items"][0]["present_all_thresholds"] is False

    identity = explorer.identity("AHR")
    assert identity["identity"]["gene_symbol"] == "AHR"
    assert identity["summary"]["best_delta_gene_effect"] == -0.45

    contexts = explorer.contexts("AHR")
    assert contexts["total"] == 1
    assert contexts["items"][0]["comparison_label"] == "TEST_AHR_vs_control"


def test_gene_models_join_crispr_expression_copy_number_and_context(tmp_path: Path) -> None:
    root = _fixture_root(tmp_path)
    explorer = GeneExplorerStore(root, MCLDataStore(root))

    payload = explorer.models("AHR", cancer_id="CANCER-X")
    assert payload["available"] is True
    assert payload["available_layers"] == {
        "gene_effect": True,
        "expression": True,
        "copy_number": True,
    }
    assert payload["total"] == 2
    assert payload["items"][0]["model_id"] == "ACH-1"
    assert payload["items"][0]["gene_effect"] == -0.9
    assert payload["items"][0]["expression"] == 6.2
    assert payload["items"][0]["copy_number"] == 1.05
    assert payload["items"][0]["memberships"][0]["assigned_group"] == "context"
