from __future__ import annotations

from .atlas import MCLAtlas
from .cohort import MCLModelCohortStore
from .fast_runtime_gene_explorer import FastRuntimeGeneExplorerStore
from .multiomics import MCLMultiOmicsStore
from .settings import MCL_ROOT
from .store import MCLDataStore


store = MCLDataStore(MCL_ROOT)
atlas_store = MCLAtlas(MCL_ROOT, store)
cohort_store = MCLModelCohortStore(MCL_ROOT)
multiomics_store = MCLMultiOmicsStore(MCL_ROOT)
gene_explorer_store = FastRuntimeGeneExplorerStore(MCL_ROOT, store)
