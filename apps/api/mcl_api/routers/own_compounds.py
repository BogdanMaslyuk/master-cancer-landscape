from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query, Response

from ..api_utils import guard
from ..schemas.own_compounds import (
    PyzCatalogResponse,
    PyzDetailResponse,
    PyzMatrixResponse,
    PyzSummaryResponse,
)
from ..state import own_compound_service


router = APIRouter()


@router.get("/api/pyz/summary", response_model=PyzSummaryResponse)
def pyz_summary():
    return guard(own_compound_service.summary)


@router.get("/api/pyz/matrix", response_model=PyzMatrixResponse)
def pyz_matrix(
    q: str | None = None,
    min_tanimoto: float = Query(0.0, ge=0.0, le=1.0),
    shortlist_only: bool = False,
):
    return guard(
        lambda: own_compound_service.matrix(
            q,
            min_tanimoto,
            shortlist_only,
        )
    )


@router.get("/api/pyz", response_model=PyzCatalogResponse)
def pyz_catalog(
    q: str | None = None,
    target_gene: str | None = None,
    shortlist_only: bool = False,
    limit: int = Query(100, ge=1, le=500),
    offset: int = Query(0, ge=0),
):
    return guard(
        lambda: own_compound_service.catalog(
            q,
            target_gene,
            shortlist_only,
            limit,
            offset,
        )
    )


# Media is deliberately outside /api/: the API contract is JSON-only and every
# /api GET endpoint has a named response model. This route serves a visual aid.
@router.get("/media/pyz/{compound_id}/structure.svg", response_class=Response)
def pyz_structure(compound_id: str):
    smiles = guard(lambda: own_compound_service.smiles(compound_id))
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError as exc:
        raise HTTPException(status_code=503, detail="RDKit is required for 2D structure rendering") from exc

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(status_code=422, detail="Stored PYZ structure could not be parsed by RDKit")
    drawer = rdMolDraw2D.MolDraw2DSVG(760, 440)
    options = drawer.drawOptions()
    options.clearBackground = False
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return Response(content=drawer.GetDrawingText(), media_type="image/svg+xml")


@router.get("/api/pyz/{compound_id}", response_model=PyzDetailResponse)
def pyz_detail(compound_id: str, target_gene: str | None = None):
    return guard(lambda: own_compound_service.detail(compound_id, target_gene))
