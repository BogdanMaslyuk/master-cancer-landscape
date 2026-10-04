from __future__ import annotations

from fastapi import APIRouter, Query

from ..api_utils import guard, split_genes
from ..schemas.atlas import AtlasResponse, CancerContextResponse, CohortResponse, MultiomicsResponse
from ..state import atlas_service


router = APIRouter()


@router.get("/api/atlas", response_model=AtlasResponse)
def atlas():
    return guard(atlas_service.atlas)


@router.get("/api/multiomics", response_model=MultiomicsResponse)
def multiomics_availability():
    return guard(atlas_service.multiomics_availability)


@router.get("/api/atlas/{cancer_id}/cohort", response_model=CohortResponse)
def cancer_model_cohort(cancer_id: str):
    return guard(lambda: atlas_service.cohort(cancer_id))


@router.get("/api/atlas/{cancer_id}/multiomics", response_model=MultiomicsResponse)
def cancer_multiomics(
    cancer_id: str,
    genes: str | None = None,
    limit: int = Query(12, ge=1, le=50),
):
    gene_tuple = split_genes(genes)
    return guard(lambda: atlas_service.context_multiomics(cancer_id, gene_tuple, limit))


@router.get("/api/atlas/{cancer_id}", response_model=CancerContextResponse)
def cancer_context(cancer_id: str):
    return guard(lambda: atlas_service.context(cancer_id))
