from __future__ import annotations

from fastapi import APIRouter, Query

from ..api_utils import guard, split_genes
from ..schemas.models import (
    ModelDependenciesResponse,
    ModelDetailResponse,
    ModelMultiomicsResponse,
    ModelsResponse,
)
from ..state import crispr_catalog, model_service


router = APIRouter()


@router.get("/api/models", response_model=ModelsResponse)
def models(
    cancer_id: str | None = None,
    group: str | None = None,
    search: str | None = None,
    sequencing_only: bool = False,
    limit: int = Query(1000, ge=1, le=5000),
):
    # This endpoint remains the curated context-membership API used by Wave 1 pages.
    # The unique CRISPR model explorer uses /api/crispr-models instead.
    return guard(lambda: model_service.models(cancer_id, group, search, sequencing_only, limit))


@router.get("/api/models/{model_id}/multiomics", response_model=ModelMultiomicsResponse)
def model_multiomics(
    model_id: str,
    genes: str | None = None,
    limit: int = Query(30, ge=1, le=100),
):
    gene_tuple = split_genes(genes)

    def load():
        try:
            return model_service.multiomics(model_id, gene_tuple, limit)
        except Exception:
            # A model can already belong to the CRISPR universe before optional
            # RNA/CNV/multi-omics indexes have been materialized locally.
            return {
                "model_id": model_id,
                "available": False,
                "status": "not_indexed",
                "note_ru": "Дополнительные молекулярные слои для этой модели ещё не проиндексированы.",
            }

    return guard(load)


@router.get("/api/models/{model_id}/dependencies", response_model=ModelDependenciesResponse)
def model_dependencies(
    model_id: str,
    search: str | None = None,
    dependency_type: str | None = None,
    domain: str | None = None,
    limit: int = Query(200, ge=1, le=1000),
    offset: int = Query(0, ge=0),
):
    return guard(
        lambda: model_service.dependencies(
            model_id,
            search,
            dependency_type,
            domain,
            limit,
            offset,
        )
    )


@router.get("/api/models/{model_id}", response_model=ModelDetailResponse)
def model(model_id: str):
    def load():
        try:
            return model_service.model(model_id)
        except Exception:
            return crispr_catalog.model(model_id)

    return guard(load)
