from __future__ import annotations

from typing import Any

from .common import ApiResponseModel


class ModelListItem(ApiResponseModel):
    cancer_id: str
    model_id: str
    cell_line_name: str | None
    assigned_group: str
    sequencing_available: bool | None
    variants: list[dict[str, Any]]


class ModelsResponse(ApiResponseModel):
    total: int
    items: list[ModelListItem]


class ModelDetailResponse(ApiResponseModel):
    pass


class ModelMultiomicsResponse(ApiResponseModel):
    pass
