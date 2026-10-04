from __future__ import annotations

from fastapi import APIRouter, Query

from ..api_utils import guard
from ..schemas.comparisons import ComparisonGenesResponse, ComparisonListResponse, ComparisonSummary
from ..state import comparison_service


router = APIRouter()


@router.get("/api/comparisons", response_model=ComparisonListResponse)
def comparisons():
    return guard(comparison_service.comparisons)


@router.get("/api/comparisons/{comparison_id}", response_model=ComparisonSummary)
def comparison(comparison_id: str):
    return guard(lambda: comparison_service.comparison(comparison_id))


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
        lambda: comparison_service.genes(
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
