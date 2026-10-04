from __future__ import annotations

from fastapi import APIRouter

from ..api_utils import guard
from ..schemas.qc import QCResponse
from ..state import qc_service


router = APIRouter()


@router.get("/api/qc", response_model=QCResponse)
def qc():
    return guard(qc_service.qc)
