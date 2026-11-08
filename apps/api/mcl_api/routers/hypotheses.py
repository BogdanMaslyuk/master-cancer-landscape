from __future__ import annotations

from fastapi import APIRouter, Query

from ..api_utils import guard
from ..schemas.hypotheses import (
    HypothesisCatalogResponse,
    HypothesisDetailResponse,
    HypothesisSummaryResponse,
)
from ..state import hypothesis_service


router = APIRouter()


@router.get("/api/hypotheses/summary", response_model=HypothesisSummaryResponse)
def hypothesis_summary():
    return guard(hypothesis_service.summary)


@router.get("/api/hypotheses", response_model=HypothesisCatalogResponse)
def hypotheses(
    q: str | None = None,
    status: str | None = None,
    target_gene: str | None = None,
    cancer: str | None = None,
    limit: int = Query(100, ge=1, le=500),
    offset: int = Query(0, ge=0),
):
    return guard(
        lambda: hypothesis_service.search(q, status, target_gene, cancer, limit, offset)
    )


@router.get("/api/hypotheses/{hypothesis_id}", response_model=HypothesisDetailResponse)
def hypothesis_detail(hypothesis_id: str):
    return guard(lambda: hypothesis_service.detail(hypothesis_id))
