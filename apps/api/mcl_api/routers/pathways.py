from __future__ import annotations

from functools import lru_cache

from fastapi import APIRouter, Query

from ..api_utils import guard
from ..state import store


router = APIRouter()


@lru_cache(maxsize=256)
def _pathways_cached(
    top_n: int | None,
    source: str | None,
    stable_only: bool,
    significant_only: bool,
    search: str | None,
    limit: int,
):
    return store.pathways(
        top_n=top_n,
        source=source,
        stable_only=stable_only,
        significant_only=significant_only,
        search=search,
        limit=limit,
    )


@lru_cache(maxsize=1)
def _pathway_stability_cached():
    return store.pathway_stability()


@lru_cache(maxsize=32)
def _network_cached(stable_only: bool, limit_terms: int):
    return store.network(stable_only=stable_only, limit_terms=limit_terms)


@router.get("/api/pathways")
def pathways(
    top_n: int | None = None,
    source: str | None = None,
    stable_only: bool = False,
    significant_only: bool = True,
    search: str | None = None,
    limit: int = Query(1000, ge=1, le=5000),
):
    return guard(lambda: _pathways_cached(top_n, source, stable_only, significant_only, search, limit))


@router.get("/api/pathways/stability")
def pathway_stability():
    return guard(_pathway_stability_cached)


@router.get("/api/network")
def network(stable_only: bool = True, limit_terms: int = Query(100, ge=1, le=500)):
    return guard(lambda: _network_cached(stable_only, limit_terms))
