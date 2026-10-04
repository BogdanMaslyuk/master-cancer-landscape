from __future__ import annotations

from fastapi import APIRouter, Query

from ..api_utils import guard, split_genes
from ..schemas.models import ModelDetailResponse, ModelMultiomicsResponse, ModelsResponse
from ..state import model_service


router = APIRouter()


@router.get("/api/models", response_model=ModelsResponse)
def models(
    cancer_id: str | None = None,
    group: str | None = None,
    search: str | None = None,
    sequencing_only: bool = False,
    limit: int = Query(1000, ge=1, le=5000),
):
    return guard(lambda: model_service.models(cancer_id, group, search, sequencing_only, limit))


@router.get("/api/models/{model_id}/multiomics", response_model=ModelMultiomicsResponse)
def model_multiomics(
    model_id: str,
    genes: str | None = None,
    limit: int = Query(30, ge=1, le=100),
):
    gene_tuple = split_genes(genes)
    return guard(lambda: model_service.multiomics(model_id, gene_tuple, limit))


@router.get("/api/models/{model_id}", response_model=ModelDetailResponse)
def model(model_id: str):
    return guard(lambda: model_service.model(model_id))
