from __future__ import annotations

from .atlas import MCLAtlas
from .cohort import MCLModelCohortStore
from .crispr_catalog import CRISPRModelCatalog
from .fast_runtime_gene_explorer import FastRuntimeGeneExplorerStore
from .model_dependencies import ModelDependencyStore
from .multiomics import MCLMultiOmicsStore
from .repositories.atlas_repository import AtlasRepository
from .repositories.comparison_repository import ComparisonRepository
from .repositories.gene_repository import GeneRepository
from .repositories.model_repository import ModelRepository
from .repositories.overview_repository import OverviewRepository
from .repositories.pathway_repository import PathwayRepository
from .repositories.qc_repository import QCRepository
from .services.atlas_service import AtlasService
from .services.comparison_service import ComparisonService
from .services.gene_service import GeneService
from .services.model_service import ModelService
from .services.overview_service import OverviewService
from .services.pathway_service import PathwayService
from .services.qc_service import QCService
from .settings import MCL_ROOT
from .store import MCLDataStore


# Canonical read-only stores. Scientific/runtime store construction lives here so
# HTTP routers never need to instantiate data-access objects themselves.
store = MCLDataStore(MCL_ROOT)
atlas_store = MCLAtlas(MCL_ROOT, store)
crispr_catalog = CRISPRModelCatalog(MCL_ROOT)
cohort_store = MCLModelCohortStore(MCL_ROOT)
multiomics_store = MCLMultiOmicsStore(MCL_ROOT)
model_dependency_store = ModelDependencyStore(MCL_ROOT)
gene_explorer_store = FastRuntimeGeneExplorerStore(MCL_ROOT, store)

# Repository -> service application boundaries.
gene_repository = GeneRepository(gene_explorer_store, store)
gene_service = GeneService(gene_repository)

atlas_repository = AtlasRepository(atlas_store, cohort_store, multiomics_store)
atlas_service = AtlasService(atlas_repository)

model_repository = ModelRepository(atlas_store, multiomics_store)
model_service = ModelService(model_repository)

comparison_repository = ComparisonRepository(store)
comparison_service = ComparisonService(comparison_repository)

pathway_repository = PathwayRepository(store)
pathway_service = PathwayService(pathway_repository)

qc_repository = QCRepository(store)
qc_service = QCService(qc_repository)

overview_repository = OverviewRepository(store)
overview_service = OverviewService(overview_repository)
