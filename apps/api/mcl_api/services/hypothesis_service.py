from __future__ import annotations

from functools import lru_cache

from ..repositories.hypothesis_repository import HypothesisRepository


class HypothesisService:
    def __init__(self, repository: HypothesisRepository):
        self.repository = repository

    @lru_cache(maxsize=1)
    def summary(self):
        return self.repository.summary()

    @lru_cache(maxsize=256)
    def search(
        self,
        search: str | None,
        status: str | None,
        target_gene: str | None,
        cancer: str | None,
        limit: int,
        offset: int,
    ):
        return self.repository.search(
            search=search,
            status=status,
            target_gene=target_gene,
            cancer=cancer,
            limit=limit,
            offset=offset,
        )

    @lru_cache(maxsize=512)
    def detail(self, hypothesis_id: str):
        return self.repository.detail(hypothesis_id)
