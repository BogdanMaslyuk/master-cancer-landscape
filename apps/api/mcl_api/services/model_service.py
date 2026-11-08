from __future__ import annotations

from functools import lru_cache

from ..repositories.model_repository import ModelRepository


class ModelService:
    """Cached application service for cell-model use-cases."""

    def __init__(self, repository: ModelRepository):
        self.repository = repository

    @lru_cache(maxsize=256)
    def models(
        self,
        cancer_id: str | None,
        group: str | None,
        search: str | None,
        sequencing_only: bool,
        limit: int,
    ):
        return self.repository.models(
            cancer_id=cancer_id,
            group=group,
            search=search,
            sequencing_only=sequencing_only,
            limit=limit,
        )

    @lru_cache(maxsize=512)
    def model(self, model_id: str):
        return self.repository.model(model_id)

    @lru_cache(maxsize=512)
    def multiomics(self, model_id: str, genes: tuple[str, ...], limit: int):
        return self.repository.multiomics(model_id, genes, limit)

    @lru_cache(maxsize=1024)
    def dependencies(
        self,
        model_id: str,
        search: str | None,
        dependency_type: str | None,
        domain: str | None,
        limit: int,
        offset: int,
    ):
        return self.repository.dependencies(
            model_id,
            search=search,
            dependency_type=dependency_type,
            domain=domain,
            limit=limit,
            offset=offset,
        )
