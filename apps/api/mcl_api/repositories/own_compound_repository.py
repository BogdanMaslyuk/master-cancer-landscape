from __future__ import annotations

from ..own_compounds import OwnCompoundStore


class OwnCompoundRepository:
    def __init__(self, store: OwnCompoundStore):
        self.store = store

    def summary(self):
        return self.store.summary()

    def catalog(
        self,
        *,
        search: str | None,
        target_gene: str | None,
        shortlist_only: bool,
        limit: int,
        offset: int,
    ):
        return self.store.catalog(
            search=search,
            target_gene=target_gene,
            shortlist_only=shortlist_only,
            limit=limit,
            offset=offset,
        )

    def matrix(
        self,
        *,
        search: str | None,
        min_tanimoto: float,
        shortlist_only: bool,
    ):
        return self.store.matrix_view(
            search=search,
            min_tanimoto=min_tanimoto,
            shortlist_only=shortlist_only,
        )

    def detail(self, compound_id: str, target_gene: str | None):
        return self.store.detail(compound_id, focus_target=target_gene)

    def smiles(self, compound_id: str) -> str:
        return self.store.smiles(compound_id)
