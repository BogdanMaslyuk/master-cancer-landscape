from __future__ import annotations

from typing import Any

from pydantic import RootModel

from .common import ApiResponseModel


class ComparisonSummary(ApiResponseModel):
    id: str
    label: str
    cancer_id: str
    comparison: str
    context_definition: str | None
    comparator_definition: str | None
    context_models_n: int | None
    comparator_models_n: int | None
    genes_analyzed_n: int | None
    depmap_release: str | None
    retrieved_at: str | None
    qc_status: str


class ComparisonListResponse(RootModel[list[ComparisonSummary]]):
    pass


class ComparisonGenesResponse(ApiResponseModel):
    comparison_id: str
    page: int
    page_size: int
    total: int
    items: list[dict[str, Any]]
