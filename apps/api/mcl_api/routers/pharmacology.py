from __future__ import annotations

from fastapi import APIRouter, Query

from ..api_utils import guard
from ..schemas.pharmacology import (
    CompoundCatalogResponse,
    CompoundPharmacologyResponse,
    ModelPharmacologyResponse,
    PharmacologyConcordanceResponse,
    PharmacologySummaryResponse,
    TargetCatalogResponse,
    TargetPharmacologyResponse,
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


@router.get("/api/compounds", response_model=CompoundCatalogResponse)
def compounds(
    q: str | None = None,
    target_gene: str | None = None,
    has_smiles: bool = False,
    limit: int = Query(100, ge=1, le=500),
    offset: int = Query(0, ge=0),
):
    return guard(lambda: pharmacology_service.compounds(q, target_gene, has_smiles, limit, offset))


@router.get("/api/compounds/{compound_id}", response_model=CompoundPharmacologyResponse)
def compound_pharmacology(
    compound_id: str,
    limit: int = Query(60, ge=1, le=200),
):
    return guard(lambda: pharmacology_service.compound(compound_id, limit))


@router.get("/api/targets", response_model=TargetCatalogResponse)
def targets(
    q: str | None = None,
    limit: int = Query(100, ge=1, le=500),
    offset: int = Query(0, ge=0),
):
    return guard(lambda: pharmacology_service.targets(q, limit, offset))


@router.get("/api/targets/{target_id}", response_model=TargetPharmacologyResponse)
def target_pharmacology(
    target_id: str,
    limit: int = Query(100, ge=1, le=250),
):
    return guard(lambda: pharmacology_service.target(target_id, limit))
