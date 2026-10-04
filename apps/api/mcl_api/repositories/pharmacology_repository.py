from __future__ import annotations

from typing import Any


class PharmacologyRepository:
    """Read-only adapter over the materialized MCL pharmacology runtime."""

    def __init__(self, pharmacology_store: Any):
        self.pharmacology_store = pharmacology_store

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

    def compound(self, compound_id: str):
        return self.pharmacology_store.compound(compound_id)
