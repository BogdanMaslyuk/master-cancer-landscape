from __future__ import annotations

from functools import lru_cache

from ..repositories.laboratory_repository import LaboratoryRepository


class LaboratoryService:
    def __init__(self, repository: LaboratoryRepository):
        self.repository = repository

    @lru_cache(maxsize=1)
    def summary(self):
        return self.repository.summary()

    @lru_cache(maxsize=64)
    def lines(self, search: str | None, species: str | None, role: str | None, matched_only: bool):
        return self.repository.lines(search=search, species=species, role=role, matched_only=matched_only)

    @lru_cache(maxsize=256)
    def candidates(
        self,
        search: str | None,
        readiness: str | None,
        cancer: str | None,
        target_gene: str | None,
        limit: int,
        offset: int,
    ):
        return self.repository.candidates(
            search=search,
            readiness=readiness,
            cancer=cancer,
            target_gene=target_gene,
            limit=limit,
            offset=offset,
        )

    @lru_cache(maxsize=512)
    def candidate_detail(self, hypothesis_id: str):
        return self.repository.candidate_detail(hypothesis_id)
