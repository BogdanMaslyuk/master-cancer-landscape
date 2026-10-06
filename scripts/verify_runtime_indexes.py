from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INDEX_DIR = ROOT / "data" / "runtime" / "explorer"
CONTRACT = "mcl-gene-explorer-runtime-v1"
SCHEMA_VERSION = "1.0"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"ERROR: {message}")


def main() -> None:
    required_files = {
        "gene_catalog.parquet",
        "gene_context_metrics.parquet",
        "gene_annotations.parquet",
        "gene_dependency_summary.parquet",
        "gene_dependency_summary.json",
        "manifest.json",
    }
    missing = sorted(name for name in required_files if not (INDEX_DIR / name).exists())
    require(not missing, f"Missing Explorer runtime indexes: {', '.join(missing)}. Run scripts/build-explorer.ps1")

    manifest = json.loads((INDEX_DIR / "manifest.json").read_text(encoding="utf-8"))
    require(manifest.get("index_contract") == CONTRACT, f"Unexpected index contract: {manifest.get('index_contract')!r}")
    require(str(manifest.get("schema_version")) == SCHEMA_VERSION, f"Unexpected schema version: {manifest.get('schema_version')!r}")
    require(manifest.get("runtime_root") == "data/runtime/explorer", f"Unexpected runtime root: {manifest.get('runtime_root')!r}")

    catalog = pd.read_parquet(INDEX_DIR / "gene_catalog.parquet")
    contexts = pd.read_parquet(INDEX_DIR / "gene_context_metrics.parquet")
    annotations = pd.read_parquet(INDEX_DIR / "gene_annotations.parquet")
    dependency = pd.read_parquet(INDEX_DIR / "gene_dependency_summary.parquet")

    catalog_required = {
        "gene_symbol",
        "mcl_domains_json",
        "protein_classes_json",
        "compartments_json",
        "hallmarks_json",
    }
    context_required = {"gene_symbol", "comparison_id", "delta_gene_effect"}
    annotation_required = {"gene_symbol", "annotation_type", "annotation_id"}
    dependency_required = {
        "gene_symbol",
        "dependency_models_n",
        "dependency_models_total_n",
        "dependency_fraction",
        "dependency_type",
        "best_cancer_name",
        "best_cancer_dependency_fraction",
        "specificity_score",
        "specificity_label_ru",
    }

    require(catalog_required.issubset(catalog.columns), f"gene_catalog.parquet missing columns: {sorted(catalog_required - set(catalog.columns))}")
    require(context_required.issubset(contexts.columns), f"gene_context_metrics.parquet missing columns: {sorted(context_required - set(contexts.columns))}")
    require(annotation_required.issubset(annotations.columns), f"gene_annotations.parquet missing columns: {sorted(annotation_required - set(annotations.columns))}")
    require(dependency_required.issubset(dependency.columns), f"gene_dependency_summary.parquet missing columns: {sorted(dependency_required - set(dependency.columns))}")
    require(not catalog.empty, "gene_catalog.parquet is empty")
    require(not dependency.empty, "gene_dependency_summary.parquet is empty")
    require(catalog["gene_symbol"].notna().all(), "gene_catalog.parquet contains null gene_symbol values")
    require(catalog["gene_symbol"].astype(str).str.len().gt(0).all(), "gene_catalog.parquet contains empty gene_symbol values")

    manifest_genes = int(manifest.get("genes_n") or 0)
    actual_genes = int(catalog["gene_symbol"].astype(str).nunique())
    dependency_genes = int(dependency["gene_symbol"].astype(str).nunique())
    require(manifest_genes == actual_genes, f"Manifest genes_n={manifest_genes} but catalog contains {actual_genes} unique genes")
    require(dependency_genes > 0, "Dependency summary contains no genes")

    model_layers = INDEX_DIR / "model_layers"
    for layer in ("gene_effect", "expression", "copy_number"):
        require((model_layers / f"{layer}.npy").exists(), f"Missing runtime model layer: {layer}.npy")
        require((model_layers / f"{layer}.json").exists(), f"Missing runtime model layer metadata: {layer}.json")

    print("Explorer runtime index verification: PASS")
    print(f"Runtime root: {INDEX_DIR.relative_to(ROOT)}")
    print(f"Contract: {CONTRACT}")
    print(f"Schema version: {SCHEMA_VERSION}")
    print(f"Genes: {actual_genes}")
    print(f"Gene dependency summaries: {dependency_genes}")
    print(f"Gene x comparison rows: {len(contexts)}")
    print(f"Annotation rows: {len(annotations)}")
    print("Fast model layers: gene_effect, expression, copy_number")


if __name__ == "__main__":
    main()
