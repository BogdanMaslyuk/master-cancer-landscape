from __future__ import annotations

from collections.abc import Callable
from typing import Any

from fastapi import HTTPException

from .store import MCLDataError


def guard(call: Callable[[], Any]) -> Any:
    """Translate data-layer errors into stable HTTP responses."""
    try:
        return call()
    except MCLDataError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except (ValueError, KeyError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


def split_genes(value: str | None) -> tuple[str, ...]:
    """Normalize a comma/semicolon separated gene list for cached API calls."""
    if not value:
        return ()
    genes: list[str] = []
    for part in value.replace(";", ",").split(","):
        gene = part.strip().upper()
        if gene and gene not in genes:
            genes.append(gene)
    return tuple(genes[:100])
