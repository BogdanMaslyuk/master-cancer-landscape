from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mcl_api.annotated_gene_explorer import AnnotatedGeneExplorerStore
from mcl_api.store import MCLDataStore


def _root(tmp_path: Path) -> Path:
    (tmp_path / "config").mkdir(parents=True)
    (tmp_path / "data/processed/depmap_genomewide").mkdir(parents=True)
    (tmp_path / "data/processed/pathways_sensitivity").mkdir(parents=True)
    (tmp_path / "data/processed/gene_explorer").mkdir(parents=True)
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
        """version: 'test-2'
status: test
domains:
  - id: signaling
    label_ru: Сигнальная трансдукция
    keywords: [aryl hydrocarbon, signaling]
    subdomains:
      - id: ahr
        label_ru: AhR-сигналинг
        keywords: [aryl hydrocarbon]
protein_classes:
  - id: transcription_factor
    label_ru: Транскрипционный фактор
    sources: [GO:MF]
    keywords: [dna-binding transcription factor]
compartments:
  - id: nucleus
    label_ru: Ядро
    sources: [GO:CC]
    keywords: [nucleus]
hallmarks:
  - id: proliferative_signaling
    label_ru: Поддержание пролиферативного сигналинга
    sources: [GO:BP]
    keywords: [cell proliferation]
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
                "gene_symbol": "AHR",
                "matched_symbol": "AHR",
                "gene_name": "aryl hydrocarbon receptor",
                "aliases_json": json.dumps(["bHLHe76", "Ah receptor"]),
                "entrez_gene_id": 196,
                "ensembl_gene_ids_json": json.dumps(["ENSG00000106546"]),
                "uniprot_swissprot_ids_json": json.dumps(["P35869"]),
                "type_of_gene": "protein-coding",
                "retrieved_at": "2026-10-04T00:00:00+00:00",
                "aggregator": "MyGene.info v3",
            }
        ]
    ).to_parquet(tmp_path / "data/processed/gene_explorer/gene_reference.parquet", index=False)
    pd.DataFrame(
        [
            {"gene_symbol": "AHR", "source": "REAC", "term_id": "R-HSA-AHR", "term_name": "Aryl hydrocarbon receptor signaling", "source_version": None},
            {"gene_symbol": "AHR", "source": "GO:MF", "term_id": "GO:0003700", "term_name": "DNA-binding transcription factor activity", "source_version": None},
            {"gene_symbol": "AHR", "source": "GO:CC", "term_id": "GO:0005634", "term_name": "nucleus", "source_version": None},
            {"gene_symbol": "AHR", "source": "GO:BP", "term_id": "GO:0008283", "term_name": "cell proliferation", "source_version": None},
        ]
    ).to_parquet(tmp_path / "data/processed/gene_explorer/gene_reference_terms.parquet", index=False)
    return tmp_path


def test_reference_snapshot_enriches_identity_and_name_search(tmp_path: Path) -> None:
    root = _root(tmp_path)
    explorer = AnnotatedGeneExplorerStore(root, MCLDataStore(root))

    identity = explorer.identity("AHR")["identity"]
    assert identity["gene_name"] == "aryl hydrocarbon receptor"
    assert identity["ensembl_gene_id"] == "ENSG00000106546"
    assert identity["uniprot_id"] == "P35869"
    assert "bHLHe76" in identity["aliases"]

    suggestions = explorer.suggest("aryl hydrocarbon")
    assert suggestions[0]["gene_symbol"] == "AHR"
    result = explorer.search(query="Ah receptor")
    assert result["total"] == 1
    assert result["items"][0]["gene_symbol"] == "AHR"


def test_reference_terms_drive_protein_class_compartment_and_hallmark(tmp_path: Path) -> None:
    root = _root(tmp_path)
    explorer = AnnotatedGeneExplorerStore(root, MCLDataStore(root))

    annotations = explorer.annotations("AHR")
    assert "Сигнальная трансдукция" in annotations["mcl_domains"]
    assert "Транскрипционный фактор" in annotations["protein_classes"]
    assert "Ядро" in annotations["compartments"]
    assert "Поддержание пролиферативного сигналинга" in annotations["hallmarks"]

    assert explorer.search(protein_class="transcription_factor")["total"] == 1
    assert explorer.search(compartment="nucleus")["total"] == 1
    assert explorer.search(hallmark="proliferative_signaling")["total"] == 1


def test_reference_layer_is_optional_and_fallback_stays_operational(tmp_path: Path) -> None:
    root = _root(tmp_path)
    (root / "data/processed/gene_explorer/gene_reference.parquet").unlink()
    (root / "data/processed/gene_explorer/gene_reference_terms.parquet").unlink()
    explorer = AnnotatedGeneExplorerStore(root, MCLDataStore(root))

    assert explorer.reference_coverage()["available"] is False
    result = explorer.search(query="AHR")
    assert result["total"] == 1
    assert explorer.identity("AHR")["identity"]["gene_symbol"] == "AHR"
