from __future__ import annotations

from fastapi import APIRouter, Query

from ..api_utils import guard
from ..schemas.pharmacology import (
    CompoundPharmacologyResponse,
    ModelPharmacologyResponse,
    PharmacologyConcordanceResponse,
    PharmacologySummaryResponse,
)
from ..state import pharmacology_service


router = APIRouter()


@router.get("/api/pharmacology/summary", response_model=PharmacologySummaryResponse)
def pharmacology_summary():
    return guard(pharmacology_service.summary)


@router.get("/api/pharmacology/concordance", response_model=PharmacologyConcordanceResponse)
def pharmacology_concordance(
    label: str | None = None,
    target_gene: str | None = None,
    compound_id: str | None = None,
    limit: int = Query(200, ge=1, le=1000),
):
    return guard(lambda: pharmacology_service.concordance(label, target_gene, compound_id, limit))


@router.get("/api/models/{model_id}/pharmacology", response_model=ModelPharmacologyResponse)
def model_pharmacology(
    model_id: str,
    source: str | None = None,
    limit: int = Query(100, ge=1, le=500),
):
    return guard(lambda: pharmacology_service.model(model_id, limit, source))


@router.get("/api/compounds/{compound_id}", response_model=CompoundPharmacologyResponse)
def compound_pharmacology(compound_id: str):
    return guard(lambda: pharmacology_service.compound(compound_id))
