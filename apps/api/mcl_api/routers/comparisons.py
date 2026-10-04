from __future__ import annotations

from functools import lru_cache

from fastapi import APIRouter, Query

from ..api_utils import guard
from ..schemas.comparisons import ComparisonGenesResponse, ComparisonListResponse, ComparisonSummary
from ..state import store


router = APIRouter()


@lru_cache(maxsize=1)
def _comparisons_cached():
    return store.comparison_specs()


@lru_cache(maxsize=64)
def _comparison_cached(comparison_id: str):
    return store.comparison(comparison_id)


@lru_cache(maxsize=512)
def _comparison_genes_cached(
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
    return store.comparison_genes(
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


@router.get("/api/comparisons", response_model=ComparisonListResponse)
def comparisons():
    return guard(_comparisons_cached)


@router.get("/api/comparisons/{comparison_id}", response_model=ComparisonSummary)
def comparison(comparison_id: str):
    return guard(lambda: _comparison_cached(comparison_id))


@router.get("/api/comparisons/{comparison_id}/genes", response_model=ComparisonGenesResponse)
def comparison_genes(
    comparison_id: str,
    page: int = Query(1, ge=1),
    page_size: int = Query(100, ge=1, le=2000),
    search: str | None = None,
    negative_delta_only: bool = False,
    exclude_broad: bool = False,
    exclude_low_sample: bool = False,
    fdr_only: bool = False,
    stable_only: bool = False,
    sort_by: str = "delta_gene_effect",
    sort_order: str = "asc",
):
    return guard(
        lambda: _comparison_genes_cached(
            comparison_id,
            page,
            page_size,
            search,
            negative_delta_only,
            exclude_broad,
            exclude_low_sample,
            fdr_only,
            stable_only,
            sort_by,
            sort_order,
        )
    )
