from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
API_DIR = ROOT / "apps" / "api"
if str(API_DIR) not in sys.path:
    sys.path.insert(0, str(API_DIR))

from mcl_api.annotated_gene_explorer import AnnotatedGeneExplorerStore  # noqa: E402
from mcl_api.store import MCLDataStore  # noqa: E402


def main() -> None:
    store = MCLDataStore(ROOT)
    explorer = AnnotatedGeneExplorerStore(ROOT, store)
    index_dir = ROOT / "data" / "processed" / "gene_explorer"
    index_dir.mkdir(parents=True, exist_ok=True)

    # Remove previous materialized copies so the builder always reflects the current
    # processed MCL outputs rather than simply re-reading an older index. Optional
    # reference snapshot files are preserved and reused when present.
    for name in ("gene_catalog.parquet", "gene_context_metrics.parquet", "gene_annotations.parquet"):
        path = index_dir / name
        if path.exists():
            path.unlink()

    explorer.catalog_frame.cache_clear()
    explorer.context_metrics_frame.cache_clear()
    explorer.formal_annotation_frame.cache_clear()
    explorer.domain_mapping_frame.cache_clear()
    explorer.secondary_mapping_frame.cache_clear()
    manifest = explorer.materialize_indexes()
    manifest["built_at"] = datetime.now(timezone.utc).isoformat()
    manifest["source"] = (
        "current processed MCL genome-wide comparisons + local DepMap multi-omics indexes + "
        "provenance-aware functional taxonomy + optional local human-gene reference snapshot"
    )
    (index_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("Gene Explorer index built")
    print(f"Genes: {manifest['genes_n']}")
    print(f"Gene × comparison rows: {manifest['context_rows_n']}")
    print(f"Comparisons: {manifest['comparisons_n']}")
    print(f"Formal annotation rows: {manifest.get('formal_annotation_rows_n', 0)}")
    print(f"Functional mapping rows: {manifest.get('functional_mapping_rows_n', 0)}")
    print(f"Secondary mapping rows: {manifest.get('secondary_mapping_rows_n', 0)}")
    print(f"Functionally annotated genes: {manifest.get('functionally_annotated_genes_n', 0)}")
    print(f"Reference genes: {manifest.get('reference_genes_n', 0)}")
    print(f"Reference terms: {manifest.get('reference_term_rows_n', 0)}")
    print("Wrote data/processed/gene_explorer/gene_catalog.parquet")
    print("Wrote data/processed/gene_explorer/gene_context_metrics.parquet")
    print("Wrote data/processed/gene_explorer/gene_annotations.parquet")


if __name__ == "__main__":
    main()
