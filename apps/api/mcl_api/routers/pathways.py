from __future__ import annotations

from fastapi import APIRouter, Query

from ..api_utils import guard
from ..schemas.pathways import NetworkResponse, PathwayListResponse, PathwayStabilityResponse
from ..state import pathway_service


router = APIRouter()


@router.get("/api/pathways", response_model=PathwayListResponse)
def pathways(
    top_n: int | None = None,
    source: str | None = None,
    stable_only: bool = False,
    significant_only: bool = True,
    search: str | None = None,
    limit: int = Query(1000, ge=1, le=5000),
):
    return guard(lambda: pathway_service.pathways(top_n, source, stable_only, significant_only, search, limit))


@router.get("/api/pathways/stability", response_model=PathwayStabilityResponse)
def pathway_stability():
    return guard(pathway_service.stability)


@router.get("/api/network", response_model=NetworkResponse)
def network(stable_only: bool = True, limit_terms: int = Query(100, ge=1, le=500)):
    return guard(lambda: pathway_service.network(stable_only, limit_terms))
