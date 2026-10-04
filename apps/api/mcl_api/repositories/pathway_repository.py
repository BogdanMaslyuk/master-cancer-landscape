from __future__ import annotations

from typing import Any


class PathwayRepository:
    """Read-only data access for pathway enrichment, stability and networks."""

    def __init__(self, store: Any):
        self.store = store

    def pathways(
        self,
        *,
        top_n: int | None,
        source: str | None,
        stable_only: bool,
        significant_only: bool,
        search: str | None,
        limit: int,
    ):
        return self.store.pathways(
            top_n=top_n,
            source=source,
            stable_only=stable_only,
            significant_only=significant_only,
            search=search,
            limit=limit,
        )

    def stability(self):
        return self.store.pathway_stability()

    def network(self, *, stable_only: bool, limit_terms: int):
        return self.store.network(stable_only=stable_only, limit_terms=limit_terms)
