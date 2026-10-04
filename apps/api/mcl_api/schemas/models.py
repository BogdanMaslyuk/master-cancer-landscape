from __future__ import annotations

from typing import Any

from pydantic import Field

from .common import ApiResponseModel


class ModelVariant(ApiResponseModel):
    gene: str | None = None
    protein_change: str | None = None
    dna_change: str | None = None
    hotspot: Any = None
    driver: Any = None
    classification: str | None = None
    kind: str | None = None


class ModelListItem(ApiResponseModel):
    cancer_id: str
    model_id: str
    cell_line_name: str | None = None
    cancer_ru: str | None = None
    organ_ru: str | None = None
    molecular_context: str | None = None
    assigned_group: str
    assigned_group_ru: str | None = None
    sequencing_available: bool | None = None
    alteration_status: str | None = None
    assignment_reason: str | None = None
    variants: list[ModelVariant] = Field(default_factory=list)


class ModelsResponse(ApiResponseModel):
    total: int
    items: list[ModelListItem]


class ModelDetailResponse(ApiResponseModel):
    pass


class ModelMultiomicsResponse(ApiResponseModel):
    pass
