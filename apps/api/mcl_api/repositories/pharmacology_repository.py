from __future__ import annotations

from typing import Any


class PharmacologyRepository:
    """Read-only adapter over pharmacology runtime and catalog indexes."""

    def __init__(self, pharmacology_store: Any, catalog_store: Any):
        self.pharmacology_store = pharmacology_store
        self.catalog_store = catalog_store

    def summary(self):
        return self.pharmacology_store.summary()

    def concordance(
        self,
        *,
        label: str | None,
        target_gene: str | None,
        compound_id: str | None,
        limit: int,
    ):
        return self.pharmacology_store.concordance(
            label=label,
            target_gene=target_gene,
            compound_id=compound_id,
            limit=limit,
        )

    def model(self, model_id: str, *, limit: int, source: str | None):
        return self.pharmacology_store.model(model_id, limit=limit, source=source)

    def compound(self, compound_id: str, *, limit: int = 60):
        return self.catalog_store.compound_detail(compound_id, limit=limit)

    def compounds(
        self,
        *,
        search: str | None,
        target_gene: str | None,
        has_smiles: bool,
        limit: int,
        offset: int,
    ):
        return self.catalog_store.compound_search(
            search=search,
            target_gene=target_gene,
            has_smiles=has_smiles,
            limit=limit,
            offset=offset,
        )

    def targets(self, *, search: str | None, limit: int, offset: int):
        return self.catalog_store.target_search(search=search, limit=limit, offset=offset)

    def target(self, target_id: str, *, limit: int = 100):
        return self.catalog_store.target_detail(target_id, limit=limit)
