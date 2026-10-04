from __future__ import annotations

import pandas as pd
import pytest

from mcl_api.runtime_gene_explorer import RuntimeGeneExplorerStore
from mcl_api.store import MCLDataError, MCLDataStore


def _runtime_store(tmp_path):
    (tmp_path / "data" / "processed" / "gene_explorer").mkdir(parents=True)
    return RuntimeGeneExplorerStore(tmp_path, MCLDataStore(tmp_path))


def test_runtime_catalog_fails_fast_without_materialized_index(tmp_path):
    runtime = _runtime_store(tmp_path)

    with pytest.raises(MCLDataError, match="build-explorer.ps1"):
        runtime.catalog_frame()


def test_runtime_reads_only_materialized_gene_indexes(tmp_path):
    runtime = _runtime_store(tmp_path)
    index_dir = tmp_path / "data" / "processed" / "gene_explorer"

    pd.DataFrame(
        [
            {
                "gene_symbol": "AHR",
                "gene_name": "aryl hydrocarbon receptor",
                "mcl_domains_json": "[]",
                "mcl_subdomains_json": "[]",
                "protein_classes_json": "[]",
                "compartments_json": "[]",
                "hallmarks_json": "[]",
            }
        ]
    ).to_parquet(index_dir / "gene_catalog.parquet", index=False)

    pd.DataFrame(
        [
            {
                "gene_symbol": "AHR",
                "comparison_id": "TEST",
                "cancer_id": "CANCER-TEST",
                "delta_gene_effect": -0.2,
            }
        ]
    ).to_parquet(index_dir / "gene_context_metrics.parquet", index=False)

    pd.DataFrame(
        [
            {
                "gene_symbol": "AHR",
                "annotation_type": "formal_term",
                "annotation_id": "GO:TEST",
                "annotation_label_ru": "test term",
                "source": "GO",
                "source_id": "GO:TEST",
                "source_term_name": "test term",
                "source_version": "test",
                "source_significant": False,
            }
        ]
    ).to_parquet(index_dir / "gene_annotations.parquet", index=False)

    catalog = runtime.catalog_frame()
    metrics = runtime.context_metrics_frame()
    formal = runtime.formal_annotation_frame()

    assert catalog["gene_symbol"].tolist() == ["AHR"]
    assert metrics["comparison_id"].tolist() == ["TEST"]
    assert formal[["source", "term_id", "term_name"]].to_dict("records") == [
        {"source": "GO", "term_id": "GO:TEST", "term_name": "test term"}
    ]
