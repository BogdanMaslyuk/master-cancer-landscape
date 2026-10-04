from __future__ import annotations

from functools import lru_cache

from fastapi import APIRouter

from ..api_utils import guard
from ..state import gene_explorer_store


router = APIRouter()


@lru_cache(maxsize=512)
def _mutation_associations_cached(gene_symbol: str):
    return gene_explorer_store.mutation_associations(gene_symbol, limit=30)


@router.get("/api/genes/{gene_symbol}/mutation-associations")
def gene_mutation_associations(gene_symbol: str):
    """Run the expensive mutation-background screen only when explicitly requested."""
    symbol = gene_symbol.strip().upper()
    return guard(lambda: _mutation_associations_cached(symbol))
