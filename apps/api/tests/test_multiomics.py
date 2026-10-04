from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mcl_api.multiomics import MCLMultiOmicsStore


def _write_parquet(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_parquet(path, index=False)


def _write_tsv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def _store(tmp_path: Path) -> MCLMultiOmicsStore:
    processed = tmp_path / "data" / "processed"
    processed.mkdir(parents=True)
    (processed / "depmap_model_multiomics_manifest.json").write_text(
        json.dumps(
            {
                "depmap_release": "26Q1",
                "built_at": "2026-10-04T00:00:00+00:00",
                "layers": {
                    "expression": {"models_n": 2, "genes_n": 3, "source_file": "expr.csv", "value_semantics": "log2(TPM + 1)"},
                    "copy_number": {"models_n": 2, "genes_n": 3, "source_file": "cn.csv", "value_semantics": "relative linear"},
                    "gene_effect": {"models_n": 2, "genes_n": 3, "source_file": "ge.csv", "value_semantics": "Chronos"},
                },
            }
        ),
        encoding="utf-8",
    )
    _write_parquet(
        processed / "depmap_model_expression.parquet",
        [
            {"model_id": "ACH-1", "KRAS": 7.0, "DNM1L": 6.4, "FIS1": 5.2},
            {"model_id": "ACH-2", "KRAS": 6.5, "DNM1L": 6.1, "FIS1": 5.0},
        ],
    )
    _write_parquet(
        processed / "depmap_model_copy_number.parquet",
        [
            {"model_id": "ACH-1", "KRAS": 1.2, "DNM1L": 1.0, "FIS1": 0.9},
            {"model_id": "ACH-2", "KRAS": 1.0, "DNM1L": 1.0, "FIS1": 1.1},
        ],
    )
    _write_parquet(
        processed / "depmap_model_gene_effect.parquet",
        [
            {"model_id": "ACH-1", "KRAS": -1.4, "DNM1L": -0.8, "FIS1": -0.5},
            {"model_id": "ACH-2", "KRAS": -0.4, "DNM1L": -0.2, "FIS1": -0.1},
        ],
    )
    _write_tsv(
        processed / "depmap_context_audit.tsv",
        [
            {"cancer_id": "C1", "model_id": "ACH-1", "assigned_group": "context"},
            {"cancer_id": "C1", "model_id": "ACH-2", "assigned_group": "comparator"},
        ],
    )
    reports = tmp_path / "outputs" / "reports"
    _write_tsv(reports / "M3_3_1_stable_recurrent_genes.tsv", [{"gene_symbol": "DNM1L"}, {"gene_symbol": "FIS1"}])
    return MCLMultiOmicsStore(tmp_path)


def test_model_multiomics_combines_three_layers(tmp_path: Path):
    store = _store(tmp_path)
    payload = store.model("ACH-1", genes=["KRAS"], limit=2)

    assert payload["availability"]["complete"] is True
    panel = {row["gene"]: row for row in payload["candidate_panel"]}
    assert panel["KRAS"]["gene_effect"] == -1.4
    assert panel["KRAS"]["expression"] == 7.0
    assert panel["KRAS"]["copy_number"] == 1.2
    assert payload["top_dependencies"][0]["gene"] == "KRAS"


def test_context_multiomics_reports_group_medians(tmp_path: Path):
    store = _store(tmp_path)
    payload = store.context("C1", genes=["KRAS"], limit=3)

    assert payload["coverage"]["gene_effect"]["context_models_n"] == 1
    assert payload["coverage"]["gene_effect"]["comparator_models_n"] == 1
    panel = {row["gene"]: row for row in payload["candidate_panel"]}
    assert panel["KRAS"]["gene_effect_context_median"] == -1.4
    assert panel["KRAS"]["gene_effect_comparator_median"] == -0.4
    assert round(panel["KRAS"]["gene_effect_delta"], 6) == -1.0
