from __future__ import annotations

from functools import lru_cache

from fastapi import APIRouter, Query

from ..api_utils import guard, split_genes
from ..schemas.atlas import AtlasResponse, CancerContextResponse, CohortResponse, MultiomicsResponse
from ..state import atlas_store, cohort_store, multiomics_store


router = APIRouter()


@lru_cache(maxsize=1)
def _atlas_cached():
    return atlas_store.atlas()


@lru_cache(maxsize=32)
def _cohort_cached(cancer_id: str):
    return cohort_store.summary(cancer_id)


@lru_cache(maxsize=32)
def _context_cached(cancer_id: str):
    return atlas_store.context(cancer_id)


@lru_cache(maxsize=1)
def _multiomics_availability_cached():
    return multiomics_store.availability()


@lru_cache(maxsize=128)
def _context_multiomics_cached(cancer_id: str, genes: tuple[str, ...], limit: int):
    return multiomics_store.context(cancer_id, list(genes) or None, limit=limit)


@router.get("/api/atlas", response_model=AtlasResponse)
def atlas():
    return guard(_atlas_cached)


@router.get("/api/multiomics", response_model=MultiomicsResponse)
def multiomics_availability():
    return guard(_multiomics_availability_cached)


@router.get("/api/atlas/{cancer_id}/cohort", response_model=CohortResponse)
def cancer_model_cohort(cancer_id: str):
    return guard(lambda: _cohort_cached(cancer_id))


@router.get("/api/atlas/{cancer_id}/multiomics", response_model=MultiomicsResponse)
def cancer_multiomics(
    cancer_id: str,
    genes: str | None = None,
    limit: int = Query(12, ge=1, le=50),
):
    gene_tuple = split_genes(genes)
    return guard(lambda: _context_multiomics_cached(cancer_id, gene_tuple, limit))


@router.get("/api/atlas/{cancer_id}", response_model=CancerContextResponse)
def cancer_context(cancer_id: str):
    return guard(lambda: _context_cached(cancer_id))
