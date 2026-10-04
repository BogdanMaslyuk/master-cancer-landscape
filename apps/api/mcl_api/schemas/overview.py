from __future__ import annotations

from typing import Any

from pydantic import Field

from .common import ApiResponseModel


class OverviewResponse(ApiResponseModel):
    genes_analyzed_n: int
    comparisons_n: int
    stable_recurrent_genes_n: int
    stable_pathways_n: int
    thresholds: list[int] = Field(default_factory=list)
    per_threshold: dict[str, Any] = Field(default_factory=dict)
    analysis_version: str | None = None
    generated_at: str | None = None
    data_release: str | None = None
