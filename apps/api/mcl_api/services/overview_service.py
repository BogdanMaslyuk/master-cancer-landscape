from __future__ import annotations

from functools import lru_cache

from ..repositories.overview_repository import OverviewRepository


class OverviewService:
    """Cached application service for the Explorer-wide summary."""

    def __init__(self, repository: OverviewRepository):
        self.repository = repository

    @lru_cache(maxsize=1)
    def summary(self):
        return self.repository.summary()
