from __future__ import annotations

import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
API_DIR = ROOT / "apps" / "api"
if str(API_DIR) not in sys.path:
    sys.path.insert(0, str(API_DIR))

from mcl_api.annotated_gene_explorer import AnnotatedGeneExplorerStore  # noqa: E402
from mcl_api.store import MCLDataStore  # noqa: E402


INDEX_CONTRACT = "mcl-gene-explorer-runtime-v1"
SCHEMA_VERSION = "1.0"
BUILD_INDEX_DIR = ROOT / "data" / "processed" / "gene_explorer"
RUNTIME_INDEX_DIR = ROOT / "data" / "runtime" / "explorer"
RUNTIME_FILES = (
    "gene_catalog.parquet",
    "gene_context_metrics.parquet",
    "gene_annotations.parquet",
    "manifest.json",
)
OPTIONAL_REFERENCE_FILES = (
    "gene_reference.parquet",
    "gene_reference_terms.parquet",
)


def _publish_runtime() -> None:
    RUNTIME_INDEX_DIR.mkdir(parents=True, exist_ok=True)
    for name in RUNTIME_FILES + OPTIONAL_REFERENCE_FILES:
        source = BUILD_INDEX_DIR / name
        target = RUNTIME_INDEX_DIR / name
        if source.exists():
            shutil.copy2(source, target)
        elif target.exists() and name in OPTIONAL_REFERENCE_FILES:
            target.unlink()


def main() -> None:
    store = MCLDataStore(ROOT)
    explorer = AnnotatedGeneExplorerStore(ROOT, store)
    BUILD_INDEX_DIR.mkdir(parents=True, exist_ok=True)

    # Build-stage artifacts remain under processed/. Runtime consumers never read
    # these files directly; after the build completes a serving snapshot is copied
    # to data/runtime/explorer/. Optional reference snapshots are preserved.
    for name in ("gene_catalog.parquet", "gene_context_metrics.parquet", "gene_annotations.parquet"):
        path = BUILD_INDEX_DIR / name
        if path.exists():
            path.unlink()

    explorer.catalog_frame.cache_clear()
    explorer.context_metrics_frame.cache_clear()
    explorer.formal_annotation_frame.cache_clear()
    explorer.domain_mapping_frame.cache_clear()
    explorer.secondary_mapping_frame.cache_clear()
    manifest = explorer.materialize_indexes()
    manifest["index_contract"] = INDEX_CONTRACT
    manifest["schema_version"] = SCHEMA_VERSION
    manifest["builder"] = "scripts/build_gene_explorer_index.py"
    manifest["built_at"] = datetime.now(timezone.utc).isoformat()
    manifest["source"] = (
        "current processed MCL genome-wide comparisons + local DepMap multi-omics indexes + "
        "provenance-aware functional taxonomy + optional local human-gene reference snapshot"
    )
    manifest["runtime_root"] = "data/runtime/explorer"
    (BUILD_INDEX_DIR / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    _publish_runtime()

    print("Gene Explorer index built")
    print(f"Contract: {INDEX_CONTRACT}")
    print(f"Schema version: {SCHEMA_VERSION}")
    print(f"Genes: {manifest['genes_n']}")
    print(f"Gene × comparison rows: {manifest['context_rows_n']}")
    print(f"Comparisons: {manifest['comparisons_n']}")
    print(f"Formal annotation rows: {manifest.get('formal_annotation_rows_n', 0)}")
    print(f"Functional mapping rows: {manifest.get('functional_mapping_rows_n', 0)}")
    print(f"Secondary mapping rows: {manifest.get('secondary_mapping_rows_n', 0)}")
    print(f"Functionally annotated genes: {manifest.get('functionally_annotated_genes_n', 0)}")
    print(f"Reference genes: {manifest.get('reference_genes_n', 0)}")
    print(f"Reference terms: {manifest.get('reference_term_rows_n', 0)}")
    print("Published runtime snapshot to data/runtime/explorer/")


if __name__ == "__main__":
    main()
