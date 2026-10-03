from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mcl_api.store import MCLDataStore


def _write_tsv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_summary_and_gene_stability(tmp_path: Path):
    (tmp_path / "config").mkdir()
    (tmp_path / "config/pathways.yaml").write_text(
        "comparisons:\n  - label: A\n    cancer_id: C1\n    comparison: x\n",
        encoding="utf-8",
    )
    _write_tsv(
        tmp_path / "outputs/reports/M3_3_1_stable_recurrent_genes.tsv",
        [{"gene_symbol": "KRAS"}],
    )
    _write_tsv(
        tmp_path / "outputs/reports/M3_3_1_stable_pathways.tsv",
        [{"source": "CORUM", "term_id": "C1", "term_name": "Complex"}],
    )
    meta = {
        "thresholds": [50, 100, 200],
        "analysis_version": "test",
        "timestamp": "2026-01-01T00:00:00Z",
    }
    (tmp_path / "outputs/reports/M3_3_1_pathway_sensitivity_meta.json").write_text(
        json.dumps(meta), encoding="utf-8"
    )
    depmeta = {
        "genes_analyzed_n": 18531,
        "depmap_release": "26Q1",
        "context_models_n": 5,
        "comparator_models_n": 8,
    }
    (tmp_path / "outputs/reports/depmap_genomewide_C1__x_meta.json").write_text(
        json.dumps(depmeta), encoding="utf-8"
    )

    store = MCLDataStore(tmp_path)
    summary = store.summary()
    assert summary["genes_analyzed_n"] == 18531
    assert summary["stable_recurrent_genes_n"] == 1
    assert summary["stable_pathways_n"] == 1
    assert summary["data_release"] == "26Q1"


def test_network_builds_gene_term_edges(tmp_path: Path):
    _write_tsv(
        tmp_path / "data/processed/pathways_sensitivity/term_stability.tsv",
        [
            {
                "source": "CORUM",
                "term_id": "4200",
                "term_name": "DNM1L-FIS1 complex",
                "significant_all_thresholds": True,
                "intersecting_gene_symbols_top50": '["DNM1L", "FIS1"]',
                "intersecting_gene_symbols_top100": '["DNM1L", "FIS1"]',
                "intersecting_gene_symbols_top200": '["DNM1L", "FIS1"]',
            }
        ],
    )
    store = MCLDataStore(tmp_path)
    graph = store.network()
    assert len(graph["nodes"]) == 3
    assert len(graph["edges"]) == 2
