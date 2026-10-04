from __future__ import annotations

from functools import lru_cache

from fastapi import APIRouter

from ..api_utils import guard
from ..state import store


router = APIRouter()


@lru_cache(maxsize=1)
def _qc_cached():
    return store.qc()


@router.get("/api/qc")
def qc():
    return guard(_qc_cached)
