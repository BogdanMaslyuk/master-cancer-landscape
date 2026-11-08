from __future__ import annotations

from typing import Any


class GeneRepository:
    """Read-only adapter over the legacy scientific store and runtime Gene Explorer."""

    def __init__(self, explorer: Any, legacy_store: Any):
        self.explorer = explorer
        self.legacy_store = legacy_store

    def genes(self, search: str | None, limit: int):
        return self.legacy_store.genes(search=search, limit=limit)

    def stable_genes(self):
        return self.legacy_store.stable_genes()

    def detail(self, gene_symbol: str):
        base = self.explorer.identity(gene_symbol)
        legacy = self.legacy_store.gene(gene_symbol)
        return {
            **base,
            "comparisons": legacy.get("comparisons") or [],
            "pathways": legacy.get("pathways") or [],
        }

    def suggest(self, query: str, limit: int):
        return self.explorer.suggest(query, limit)

    def facets(self):
        return self.explorer.facets()

    def search(self, **filters: Any):
        return self.explorer.search(**filters)

    def matrix(self, **filters: Any):
        return self.explorer.matrix(**filters)

    def contexts(self, gene_symbol: str):
        return self.explorer.contexts(gene_symbol)

    def models(self, gene_symbol: str, **filters: Any):
        return self.explorer.models(gene_symbol, **filters)

    def annotations(self, gene_symbol: str):
        return self.explorer.annotations(gene_symbol)

    def dependency_landscape(self, gene_symbol: str):
        return self.explorer.dependency_landscape(gene_symbol)

    def mutation_associations(self, gene_symbol: str, *, limit: int = 30):
        return self.explorer.mutation_associations(gene_symbol, limit=limit)
