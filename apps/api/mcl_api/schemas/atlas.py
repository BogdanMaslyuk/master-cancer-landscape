from __future__ import annotations

from typing import Any

from .common import ApiResponseModel


class AtlasOrgan(ApiResponseModel):
    id: str
    name_ru: str | None
    name_en: str | None
    icon: str | None
    contexts: list[dict[str, Any]]
    contexts_n: int
    models_n: int
    analyses_n: int


class AtlasResponse(ApiResponseModel):
    organs_n: int
    contexts_n: int
    models_n: int
    analyses_n: int
    organs: list[AtlasOrgan]


class CancerContextResponse(ApiResponseModel):
    id: str
    name: str
    organ_id: str
    organ_ru: str | None
    cancer_ru: str | None
    molecular_context: str | None
    models_n: int
    context_models_n: int
    comparator_models_n: int
    comparisons_n: int


class CohortResponse(ApiResponseModel):
    pass


class MultiomicsResponse(ApiResponseModel):
    pass
