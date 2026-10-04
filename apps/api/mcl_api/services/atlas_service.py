from __future__ import annotations

from functools import lru_cache

from ..repositories.atlas_repository import AtlasRepository


class AtlasService:
    """Cached application service for disease/context navigation use-cases."""

    def __init__(self, repository: AtlasRepository):
        self.repository = repository

    @lru_cache(maxsize=1)
    def atlas(self):
        return self.repository.atlas()

    @lru_cache(maxsize=32)
    def cohort(self, cancer_id: str):
        return self.repository.cohort(cancer_id)

    @lru_cache(maxsize=32)
    def context(self, cancer_id: str):
        return self.repository.context(cancer_id)

    @lru_cache(maxsize=1)
    def multiomics_availability(self):
        return self.repository.multiomics_availability()

    @lru_cache(maxsize=128)
    def context_multiomics(self, cancer_id: str, genes: tuple[str, ...], limit: int):
        return self.repository.context_multiomics(cancer_id, genes, limit)
