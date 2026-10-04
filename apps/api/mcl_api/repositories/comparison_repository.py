from __future__ import annotations

from typing import Any


class ComparisonRepository:
    """Read-only data access for configured genome-wide comparison results."""

    def __init__(self, store: Any):
        self.store = store

    def comparisons(self):
        return self.store.comparison_specs()

    def comparison(self, comparison_id: str):
        return self.store.comparison(comparison_id)

    def genes(
        self,
        comparison_id: str,
        *,
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
        return self.store.comparison_genes(
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
