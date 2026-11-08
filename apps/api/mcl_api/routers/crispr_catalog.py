from __future__ import annotations

from fastapi import APIRouter, Query

from ..api_utils import guard
from ..schemas.crispr_catalog import CRISPRAtlasResponse, CRISPRModelResponse, CRISPRModelsResponse
from ..state import crispr_catalog


router = APIRouter()


@router.get("/api/crispr-atlas", response_model=CRISPRAtlasResponse)
def crispr_atlas():
    return guard(crispr_catalog.atlas)


@router.get("/api/crispr-models", response_model=CRISPRModelsResponse)
def crispr_models(
    organ_id: str | None = None,
    cancer_id: str | None = None,
    subtype: str | None = None,
    search: str | None = None,
    limit: int = Query(5000, ge=1, le=5000),
):
    return guard(
        lambda: crispr_catalog.models(
            organ_id=organ_id,
            cancer_id=cancer_id,
            subtype=subtype,
            search=search,
            limit=limit,
        )
    )


@router.get("/api/crispr-models/{model_id}", response_model=CRISPRModelResponse)
def crispr_model(model_id: str):
    return guard(lambda: crispr_catalog.model(model_id))
