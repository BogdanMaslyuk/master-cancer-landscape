from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
API_DIR = ROOT / "apps" / "api"
if str(API_DIR) not in sys.path:
    sys.path.insert(0, str(API_DIR))

from mcl_api.gene_explorer import GeneExplorerStore  # noqa: E402
from mcl_api.store import MCLDataStore  # noqa: E402


def main() -> None:
    store = MCLDataStore(ROOT)
    explorer = GeneExplorerStore(ROOT, store)
    index_dir = ROOT / "data" / "processed" / "gene_explorer"
    index_dir.mkdir(parents=True, exist_ok=True)

    # Remove previous materialized copies so the builder always reflects the current
    # processed MCL outputs rather than simply re-reading an older index.
    for name in ("gene_catalog.parquet", "gene_context_metrics.parquet"):
        path = index_dir / name
        if path.exists():
            path.unlink()

    explorer.catalog_frame.cache_clear()
    explorer.context_metrics_frame.cache_clear()
    manifest = explorer.materialize_indexes()
    manifest["built_at"] = datetime.now(timezone.utc).isoformat()
    manifest["source"] = "current processed MCL genome-wide comparisons + local DepMap multi-omics indexes"
    (index_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("Gene Explorer index built")
    print(f"Genes: {manifest['genes_n']}")
    print(f"Gene × comparison rows: {manifest['context_rows_n']}")
    print(f"Comparisons: {manifest['comparisons_n']}")
    print("Wrote data/processed/gene_explorer/gene_catalog.parquet")
    print("Wrote data/processed/gene_explorer/gene_context_metrics.parquet")


if __name__ == "__main__":
    main()
