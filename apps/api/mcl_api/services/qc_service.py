from __future__ import annotations

from functools import lru_cache

from ..repositories.qc_repository import QCRepository


class QCService:
    """Cached application service for Explorer quality-control summaries."""

    def __init__(self, repository: QCRepository):
        self.repository = repository

    @lru_cache(maxsize=1)
    def qc(self):
        return self.repository.qc()
