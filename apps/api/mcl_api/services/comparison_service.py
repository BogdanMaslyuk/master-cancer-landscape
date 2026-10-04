from __future__ import annotations

from functools import lru_cache

from ..repositories.comparison_repository import ComparisonRepository


class ComparisonService:
    """Cached application service for genome-wide comparison use-cases."""

    def __init__(self, repository: ComparisonRepository):
        self.repository = repository

    @lru_cache(maxsize=1)
    def comparisons(self):
        return self.repository.comparisons()

    @lru_cache(maxsize=64)
    def comparison(self, comparison_id: str):
        return self.repository.comparison(comparison_id)

    @lru_cache(maxsize=512)
    def genes(
        self,
        comparison_id: str,
        page: int,
        page_size: int,
        search: str | None,
        negative_delta_only: bool,
        exclude_broad: bool,
        exclude_low_sample: bool,
        fdr_only: bool,
        stable_only: bool,
        sort_by: str,
        sort_order: str,
    ):
        return self.repository.genes(
            comparison_id,
            page=page,
            page_size=page_size,
            search=search,
            negative_delta_only=negative_delta_only,
            exclude_broad=exclude_broad,
            exclude_low_sample=exclude_low_sample,
            fdr_only=fdr_only,
            stable_only=stable_only,
            sort_by=sort_by,
            sort_order=sort_order,
        )
