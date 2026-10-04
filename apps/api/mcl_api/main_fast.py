from __future__ import annotations

from . import main as base
from .runtime_gene_explorer import RuntimeGeneExplorerStore

# Keep the same FastAPI app/routes but replace only the Gene Explorer store with
# a runtime implementation that reads materialized indexes instead of rebuilding
# ontology mappings during interactive requests.
base.gene_explorer_store = RuntimeGeneExplorerStore(base.MCL_ROOT, base.store)

# Clear route-level caches in case this module is reloaded in the same process.
for cache in (
    base._gene_cached,
    base._gene_suggest_cached,
    base._gene_facets_cached,
    base._gene_search_cached,
    base._gene_matrix_cached,
    base._gene_contexts_cached,
    base._gene_models_cached,
    base._gene_annotations_cached,
):
    cache.cache_clear()

app = base.app
