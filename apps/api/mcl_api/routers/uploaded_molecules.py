from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query, Response

from ..api_utils import guard
from ..schemas.uploaded_molecules import (
    MoleculeAnalysisResponse,
    MoleculeAnalyzeRequest,
    MoleculeCellLinesResponse,
    MoleculeDetailResponse,
    MoleculeExperimentResponse,
    MoleculeTargetsResponse,
)
from ..state import uploaded_molecule_service


router = APIRouter()


@router.post("/api/molecules/analyze", response_model=MoleculeAnalysisResponse)
def analyze_molecule(request: MoleculeAnalyzeRequest):
    return guard(
        lambda: uploaded_molecule_service.analyze(
            request.smiles, request.name, request.force_recompute
        )
    )


@router.get("/api/molecules/{analysis_id}", response_model=MoleculeDetailResponse)
def molecule_detail(analysis_id: str):
    return guard(lambda: uploaded_molecule_service.detail(analysis_id))


@router.get("/api/molecules/{analysis_id}/targets", response_model=MoleculeTargetsResponse)
def molecule_targets(analysis_id: str):
    return guard(lambda: uploaded_molecule_service.targets(analysis_id))


@router.get("/api/molecules/{analysis_id}/cell-lines", response_model=MoleculeCellLinesResponse)
def molecule_cell_lines(
    analysis_id: str,
    target_gene: str | None = None,
    laboratory_only: bool = Query(False),
):
    return guard(lambda: uploaded_molecule_service.cell_lines(analysis_id, target_gene, laboratory_only))


@router.get("/api/molecules/{analysis_id}/experiment", response_model=MoleculeExperimentResponse)
def molecule_experiment(analysis_id: str):
    return guard(lambda: uploaded_molecule_service.experiment(analysis_id))


@router.get("/media/molecules/{analysis_id}/structure.svg", response_class=Response)
def molecule_structure(analysis_id: str):
    smiles = guard(lambda: uploaded_molecule_service.smiles(analysis_id))
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError as exc:
        raise HTTPException(status_code=503, detail="RDKit is required for 2D structure rendering") from exc

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(status_code=422, detail="Stored structure could not be parsed by RDKit")
    drawer = rdMolDraw2D.MolDraw2DSVG(760, 440)
    drawer.drawOptions().clearBackground = False
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return Response(content=drawer.GetDrawingText(), media_type="image/svg+xml")
