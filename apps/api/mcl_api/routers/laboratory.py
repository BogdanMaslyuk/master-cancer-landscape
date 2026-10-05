from __future__ import annotations

from fastapi import APIRouter, Query

from ..api_utils import guard
from ..schemas.laboratory import (
    LaboratoryCandidateDetailResponse,
    LaboratoryCandidatesResponse,
    LaboratoryLinesResponse,
    LaboratoryMechanismPanelsResponse,
    LaboratorySummaryResponse,
)
from ..state import laboratory_service


router = APIRouter()


@router.get("/api/laboratory/summary", response_model=LaboratorySummaryResponse)
def laboratory_summary():
    return guard(laboratory_service.summary)


@router.get("/api/laboratory/lines", response_model=LaboratoryLinesResponse)
def laboratory_lines(
    q: str | None = None,
    species: str | None = None,
    role: str | None = None,
    matched_only: bool = False,
):
    return guard(lambda: laboratory_service.lines(q, species, role, matched_only))


@router.get("/api/laboratory/candidates", response_model=LaboratoryCandidatesResponse)
def laboratory_candidates(
    q: str | None = None,
    readiness: str | None = None,
    cancer: str | None = None,
    target_gene: str | None = None,
    limit: int = Query(100, ge=1, le=500),
    offset: int = Query(0, ge=0),
):
    return guard(
        lambda: laboratory_service.candidates(q, readiness, cancer, target_gene, limit, offset)
    )


@router.get("/api/laboratory/mechanism-panels", response_model=LaboratoryMechanismPanelsResponse)
def laboratory_mechanism_panels():
    return guard(laboratory_service.mechanism_panels)


@router.get("/api/laboratory/candidates/{hypothesis_id}", response_model=LaboratoryCandidateDetailResponse)
def laboratory_candidate_detail(hypothesis_id: str):
    return guard(lambda: laboratory_service.candidate_detail(hypothesis_id))
