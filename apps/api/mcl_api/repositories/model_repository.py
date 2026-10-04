from __future__ import annotations

from typing import Any


class ModelRepository:
    """Read-only data access for cell-model catalog, detail and multi-omics."""

    def __init__(self, atlas_store: Any, multiomics_store: Any):
        self.atlas_store = atlas_store
        self.multiomics_store = multiomics_store

    def models(
        self,
        *,
        cancer_id: str | None,
        group: str | None,
        search: str | None,
        sequencing_only: bool,
        limit: int,
    ):
        return self.atlas_store.models(
            cancer_id=cancer_id,
            group=group,
            search=search,
            sequencing_only=sequencing_only,
            limit=limit,
        )

    def model(self, model_id: str):
        return self.atlas_store.model(model_id)

    def multiomics(self, model_id: str, genes: tuple[str, ...], limit: int):
        return self.multiomics_store.model(model_id, list(genes) or None, limit=limit)
