from __future__ import annotations

from functools import lru_cache

from fastapi import APIRouter

from ..api_utils import guard
from ..schemas.qc import QCResponse
from ..state import store


router = APIRouter()


@lru_cache(maxsize=1)
def _qc_cached():
    return store.qc()


@router.get("/api/qc", response_model=QCResponse)
def qc():
    return guard(_qc_cached)
