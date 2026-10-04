from __future__ import annotations

from .atlas import MCLAtlas
from .cohort import MCLModelCohortStore
from .multiomics import MCLMultiOmicsStore
from .runtime_gene_explorer import RuntimeGeneExplorerStore
from .settings import MCL_ROOT
from .store import MCLDataStore


store = MCLDataStore(MCL_ROOT)
atlas_store = MCLAtlas(MCL_ROOT, store)
cohort_store = MCLModelCohortStore(MCL_ROOT)
multiomics_store = MCLMultiOmicsStore(MCL_ROOT)
gene_explorer_store = RuntimeGeneExplorerStore(MCL_ROOT, store)
