from __future__ import annotations

from functools import lru_cache

from ..repositories.pharmacology_repository import PharmacologyRepository


class PharmacologyService:
    """Cached application service for pharmacology use-cases."""

    def __init__(self, repository: PharmacologyRepository):
        self.repository = repository

    @lru_cache(maxsize=1)
    def summary(self):
        return self.repository.summary()

    @lru_cache(maxsize=512)
    def concordance(
        self,
        label: str | None,
        target_gene: str | None,
        compound_id: str | None,
        limit: int,
    ):
        return self.repository.concordance(
            label=label,
            target_gene=target_gene,
            compound_id=compound_id,
            limit=limit,
        )

    @lru_cache(maxsize=2048)
    def model(self, model_id: str, limit: int, source: str | None):
        return self.repository.model(model_id, limit=limit, source=source)

    @lru_cache(maxsize=2048)
    def compound(self, compound_id: str, limit: int = 60):
        return self.repository.compound(compound_id, limit=limit)

    @lru_cache(maxsize=2048)
    def compounds(
        self,
        search: str | None,
        target_gene: str | None,
        has_smiles: bool,
        limit: int,
        offset: int,
    ):
        return self.repository.compounds(
            search=search,
            target_gene=target_gene,
            has_smiles=has_smiles,
            limit=limit,
            offset=offset,
        )

    @lru_cache(maxsize=1024)
    def targets(self, search: str | None, limit: int, offset: int):
        return self.repository.targets(search=search, limit=limit, offset=offset)

    @lru_cache(maxsize=2048)
    def target(self, target_id: str, limit: int = 100):
        return self.repository.target(target_id, limit=limit)
