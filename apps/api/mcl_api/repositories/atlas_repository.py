from __future__ import annotations

from typing import Any


class AtlasRepository:
    """Read-only data access for disease/context navigation and context multi-omics."""

    def __init__(self, atlas_store: Any, cohort_store: Any, multiomics_store: Any):
        self.atlas_store = atlas_store
        self.cohort_store = cohort_store
        self.multiomics_store = multiomics_store

    def atlas(self):
        return self.atlas_store.atlas()

    def cohort(self, cancer_id: str):
        return self.cohort_store.summary(cancer_id)

    def context(self, cancer_id: str):
        return self.atlas_store.context(cancer_id)

    def multiomics_availability(self):
        return self.multiomics_store.availability()

    def context_multiomics(self, cancer_id: str, genes: tuple[str, ...], limit: int):
        return self.multiomics_store.context(cancer_id, list(genes) or None, limit=limit)
