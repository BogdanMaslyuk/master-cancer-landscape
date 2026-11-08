from __future__ import annotations

from functools import lru_cache

from ..repositories.own_compound_repository import OwnCompoundRepository


class OwnCompoundService:
    def __init__(self, repository: OwnCompoundRepository):
        self.repository = repository

    @lru_cache(maxsize=1)
    def summary(self):
        return self.repository.summary()

    @lru_cache(maxsize=256)
    def catalog(
        self,
        search: str | None,
        target_gene: str | None,
        shortlist_only: bool,
        limit: int,
        offset: int,
    ):
        return self.repository.catalog(
            search=search,
            target_gene=target_gene,
            shortlist_only=shortlist_only,
            limit=limit,
            offset=offset,
        )

    @lru_cache(maxsize=256)
    def matrix(
        self,
        search: str | None,
        min_tanimoto: float,
        shortlist_only: bool,
    ):
        return self.repository.matrix(
            search=search,
            min_tanimoto=min_tanimoto,
            shortlist_only=shortlist_only,
        )

    @lru_cache(maxsize=512)
    def detail(self, compound_id: str, target_gene: str | None):
        return self.repository.detail(compound_id, target_gene)

    @lru_cache(maxsize=128)
    def smiles(self, compound_id: str) -> str:
        return self.repository.smiles(compound_id)
