from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mcl_api.matrix_gene_explorer import MatrixGeneExplorerStore
from mcl_api.store import MCLDataStore


def _root(tmp_path: Path) -> Path:
    (tmp_path / "config").mkdir(parents=True)
    (tmp_path / "data/processed/depmap_genomewide").mkdir(parents=True)
    (tmp_path / "data/processed/pathways_sensitivity").mkdir(parents=True)
    (tmp_path / "outputs/reports").mkdir(parents=True)
    (tmp_path / "outputs/qc").mkdir(parents=True)

    (tmp_path / "config/pathways.yaml").write_text(
        """organism: hsapiens
comparisons:
  - label: C1_target_vs_control
    cancer_id: CANCER-1
    comparison: test
  - label: C2_target_vs_control
    cancer_id: CANCER-2
    comparison: test
""",
        encoding="utf-8",
    )
    (tmp_path / "config/cancer_contexts.yaml").write_text(
        """contexts:
  CANCER-1:
    name: Cancer one
  CANCER-2:
    name: Cancer two
""",
        encoding="utf-8",
    )
    (tmp_path / "config/mcl_functional_domains.yaml").write_text(
        """version: 'test'
domains: []
protein_classes: []
compartments: []
hallmarks: []
""",
        encoding="utf-8",
    )

    rows = {
        "CANCER-1__test": [
            {"gene_symbol": "AHR", "context_median_gene_effect": -0.7, "comparator_median_gene_effect": -0.2, "delta_gene_effect": -0.5, "cliffs_delta": -0.6, "q_value": 0.01, "fdr_0_05": True, "broad_dependency_warning": False, "low_sample_size": False},
            {"gene_symbol": "KRAS", "context_median_gene_effect": -1.1, "comparator_median_gene_effect": -0.5, "delta_gene_effect": -0.6, "cliffs_delta": -0.7, "q_value": 0.002, "fdr_0_05": True, "broad_dependency_warning": False, "low_sample_size": False},
        ],
        "CANCER-2__test": [
            {"gene_symbol": "AHR", "context_median_gene_effect": -0.3, "comparator_median_gene_effect": -0.25, "delta_gene_effect": -0.05, "cliffs_delta": -0.1, "q_value": 0.8, "fdr_0_05": False, "broad_dependency_warning": False, "low_sample_size": False},
            {"gene_symbol": "KRAS", "context_median_gene_effect": -0.6, "comparator_median_gene_effect": -0.55, "delta_gene_effect": -0.05, "cliffs_delta": -0.05, "q_value": 0.7, "fdr_0_05": False, "broad_dependency_warning": False, "low_sample_size": False},
        ],
    }
    for comparison_id, data in rows.items():
        pd.DataFrame(data).to_csv(
            tmp_path / f"data/processed/depmap_genomewide/{comparison_id}_genes.tsv",
            sep="\t",
            index=False,
        )
        (tmp_path / f"outputs/reports/depmap_genomewide_{comparison_id}_meta.json").write_text(
            json.dumps({"context_models_n": 4, "comparator_models_n": 5, "genes_analyzed_n": 2, "depmap_release": "26Q1"}),
            encoding="utf-8",
        )
        (tmp_path / f"outputs/qc/depmap_genomewide_{comparison_id}_qc.json").write_text("[]", encoding="utf-8")
    return tmp_path


def test_matrix_returns_same_gene_across_comparisons(tmp_path: Path) -> None:
    root = _root(tmp_path)
    explorer = MatrixGeneExplorerStore(root, MCLDataStore(root))

    payload = explorer.matrix(query="AHR", exclude_broad=True, exclude_low_sample=True)
    assert payload["genes_returned_n"] == 1
    assert payload["comparisons_n"] == 2
    row = payload["rows"][0]
    assert row["gene_symbol"] == "AHR"
    assert row["cells"]["CANCER-1__test"]["delta_gene_effect"] == -0.5
    assert row["cells"]["CANCER-1__test"]["fdr_0_05"] is True
    assert row["cells"]["CANCER-2__test"]["delta_gene_effect"] == -0.05
    assert row["cells"]["CANCER-2__test"]["fdr_0_05"] is False


def test_matrix_cancer_filter_limits_columns_and_selects_by_that_context(tmp_path: Path) -> None:
    root = _root(tmp_path)
    explorer = MatrixGeneExplorerStore(root, MCLDataStore(root))

    payload = explorer.matrix(cancer_id="CANCER-1", q_value_max=0.05)
    assert payload["comparisons_n"] == 1
    assert payload["comparisons"][0]["cancer_id"] == "CANCER-1"
    assert {row["gene_symbol"] for row in payload["rows"]} == {"AHR", "KRAS"}
