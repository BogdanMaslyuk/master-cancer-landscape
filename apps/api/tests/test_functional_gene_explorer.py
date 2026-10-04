from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mcl_api.functional_gene_explorer import FunctionalGeneExplorerStore
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
  - label: TEST_AHR_vs_control
    cancer_id: CANCER-X
    comparison: test
""",
        encoding="utf-8",
    )
    (tmp_path / "config/cancer_contexts.yaml").write_text(
        """contexts:
  CANCER-X:
    name: Test cancer
""",
        encoding="utf-8",
    )
    (tmp_path / "config/mcl_functional_domains.yaml").write_text(
        """version: 'test-1'
status: test
domains:
  - id: signaling
    label_ru: Сигнальная трансдукция
    keywords: [aryl hydrocarbon, signaling]
    subdomains:
      - id: ahr
        label_ru: AhR-сигналинг
        keywords: [aryl hydrocarbon]
  - id: mitochondrial_biology
    label_ru: Митохондриальная биология
    keywords: [mitochond]
""",
        encoding="utf-8",
    )
    (tmp_path / "outputs/reports/depmap_genomewide_CANCER-X__test_meta.json").write_text(
        json.dumps({"context_models_n": 2, "comparator_models_n": 2, "genes_analyzed_n": 2, "depmap_release": "26Q1"}),
        encoding="utf-8",
    )
    (tmp_path / "outputs/qc/depmap_genomewide_CANCER-X__test_qc.json").write_text("[]", encoding="utf-8")

    pd.DataFrame(
        [
            {
                "gene_symbol": "AHR",
                "context_median_gene_effect": -0.7,
                "comparator_median_gene_effect": -0.2,
                "delta_gene_effect": -0.5,
                "cliffs_delta": -0.6,
                "q_value": 0.01,
                "fdr_0_05": True,
                "broad_dependency_warning": False,
                "low_sample_size": False,
            },
            {
                "gene_symbol": "KRAS",
                "context_median_gene_effect": -1.0,
                "comparator_median_gene_effect": -0.4,
                "delta_gene_effect": -0.6,
                "cliffs_delta": -0.7,
                "q_value": 0.005,
                "fdr_0_05": True,
                "broad_dependency_warning": False,
                "low_sample_size": False,
            },
        ]
    ).to_csv(tmp_path / "data/processed/depmap_genomewide/CANCER-X__test_genes.tsv", sep="\t", index=False)

    pd.DataFrame(
        [
            {
                "source": "REAC",
                "term_id": "R-HSA-TEST",
                "term_name": "Aryl hydrocarbon receptor signaling",
                "significant": True,
                "top_n": 100,
                "intersecting_gene_symbols_json": json.dumps(["AHR"]),
            },
            {
                "source": "GO:BP",
                "term_id": "GO:TEST2",
                "term_name": "mitochondrial fission",
                "significant": True,
                "top_n": 100,
                "intersecting_gene_symbols_json": json.dumps(["KRAS"]),
            },
        ]
    ).to_csv(tmp_path / "data/processed/pathways_sensitivity/enrichment_all.tsv", sep="\t", index=False)
    return tmp_path


def test_functional_domains_are_source_backed_and_filterable(tmp_path: Path) -> None:
    root = _root(tmp_path)
    explorer = FunctionalGeneExplorerStore(root, MCLDataStore(root))

    facets = explorer.facets()
    signaling = next(row for row in facets["domains"] if row["id"] == "signaling")
    assert signaling["genes_n"] == 1

    result = explorer.search(domain="signaling")
    assert result["total"] == 1
    assert result["items"][0]["gene_symbol"] == "AHR"

    annotations = explorer.annotations("AHR")
    assert "Сигнальная трансдукция" in annotations["mcl_domains"]
    assert "AhR-сигналинг" in annotations["subdomains"]
    assert annotations["mcl_domain_details"][0]["source"] == "REAC"
    assert annotations["mcl_domain_details"][0]["source_id"] == "R-HSA-TEST"
    assert annotations["mcl_domain_details"][0]["mapping_method"] == "rule_based_from_formal_annotation"


def test_materialized_index_keeps_functional_provenance(tmp_path: Path) -> None:
    root = _root(tmp_path)
    explorer = FunctionalGeneExplorerStore(root, MCLDataStore(root))
    manifest = explorer.materialize_indexes()

    assert manifest["functionally_annotated_genes_n"] == 2
    assert (root / "data/processed/gene_explorer/gene_annotations.parquet").exists()
    catalog = pd.read_parquet(root / "data/processed/gene_explorer/gene_catalog.parquet")
    ahr = catalog[catalog["gene_symbol"] == "AHR"].iloc[0]
    assert "signaling" in json.loads(ahr["mcl_domains_json"])
