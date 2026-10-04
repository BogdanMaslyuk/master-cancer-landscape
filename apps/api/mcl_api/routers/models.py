from __future__ import annotations

from functools import lru_cache

from fastapi import APIRouter, Query

from ..api_utils import guard, split_genes
from ..schemas.models import ModelDetailResponse, ModelMultiomicsResponse, ModelsResponse
from ..state import atlas_store, multiomics_store


router = APIRouter()


@lru_cache(maxsize=256)
def _models_cached(
    cancer_id: str | None,
    group: str | None,
    search: str | None,
    sequencing_only: bool,
    limit: int,
):
    return atlas_store.models(
        cancer_id=cancer_id,
        group=group,
        search=search,
        sequencing_only=sequencing_only,
        limit=limit,
    )


@lru_cache(maxsize=512)
def _model_cached(model_id: str):
    return atlas_store.model(model_id)


@lru_cache(maxsize=512)
def _model_multiomics_cached(model_id: str, genes: tuple[str, ...], limit: int):
    return multiomics_store.model(model_id, list(genes) or None, limit=limit)


@router.get("/api/models", response_model=ModelsResponse)
def models(
    cancer_id: str | None = None,
    group: str | None = None,
    search: str | None = None,
    sequencing_only: bool = False,
    limit: int = Query(1000, ge=1, le=5000),
):
    return guard(lambda: _models_cached(cancer_id, group, search, sequencing_only, limit))


@router.get("/api/models/{model_id}/multiomics", response_model=ModelMultiomicsResponse)
def model_multiomics(
    model_id: str,
    genes: str | None = None,
    limit: int = Query(30, ge=1, le=100),
):
    gene_tuple = split_genes(genes)
    return guard(lambda: _model_multiomics_cached(model_id, gene_tuple, limit))


@router.get("/api/models/{model_id}", response_model=ModelDetailResponse)
def model(model_id: str):
    return guard(lambda: _model_cached(model_id))
