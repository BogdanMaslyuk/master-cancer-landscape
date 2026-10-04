from __future__ import annotations

from functools import lru_cache

from ..repositories.pathway_repository import PathwayRepository


class PathwayService:
    """Cached application service for pathway and network use-cases."""

    def __init__(self, repository: PathwayRepository):
        self.repository = repository

    @lru_cache(maxsize=256)
    def pathways(
        self,
        top_n: int | None,
        source: str | None,
        stable_only: bool,
        significant_only: bool,
        search: str | None,
        limit: int,
    ):
        return self.repository.pathways(
            top_n=top_n,
            source=source,
            stable_only=stable_only,
            significant_only=significant_only,
            search=search,
            limit=limit,
        )

    @lru_cache(maxsize=32)
    def stability(
        self,
        stable_only: bool = False,
        min_significant_thresholds: int | None = None,
        limit: int = 5000,
    ):
        return self.repository.stability(
            stable_only=stable_only,
            min_significant_thresholds=min_significant_thresholds,
            limit=limit,
        )

    @lru_cache(maxsize=32)
    def network(self, stable_only: bool, limit_terms: int):
        return self.repository.network(stable_only=stable_only, limit_terms=limit_terms)
