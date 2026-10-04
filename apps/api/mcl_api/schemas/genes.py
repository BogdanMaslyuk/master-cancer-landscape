from __future__ import annotations

from typing import Any

from pydantic import Field

from .common import ApiResponseModel


class CoverageSummary(ApiResponseModel):
    annotated_genes_n: int | None = None
    gene_universe_n: int | None = None
    coverage_fraction: float | None = None
    status: str | None = None
    reason: str | None = None


class ReferenceCoverageSummary(ApiResponseModel):
    resolved_genes_n: int | None = None
    gene_universe_n: int | None = None
    coverage_fraction: float | None = None
    available: bool | None = None


class GeneSearchResponse(ApiResponseModel):
    page: int
    page_size: int
    total: int
    pages: int
    sort_by: str
    sort_order: str
    items: list[dict[str, Any]]
    functional_coverage: CoverageSummary | None = None
    reference_coverage: ReferenceCoverageSummary | None = None


class GeneFacetResponse(ApiResponseModel):
    taxonomy_version: str | None = None
    status: str | None = None
    domains: list[dict[str, Any]] = Field(default_factory=list)
    protein_classes: list[dict[str, Any]] = Field(default_factory=list)
    compartments: list[dict[str, Any]] = Field(default_factory=list)
    hallmarks: list[dict[str, Any]] = Field(default_factory=list)
    sources: list[dict[str, Any]] = Field(default_factory=list)
    coverage: CoverageSummary | dict[str, Any] | None = None
    reference_coverage: ReferenceCoverageSummary | dict[str, Any] | None = None
    provenance_note: str | None = None


class GeneMatrixResponse(ApiResponseModel):
    genes_total_after_filters: int
    genes_returned_n: int
    comparisons_n: int
    comparisons: list[dict[str, Any]]
    rows: list[dict[str, Any]]
    filters: dict[str, Any]
    interpretation: dict[str, Any]


class GeneDetailResponse(ApiResponseModel):
    identity: dict[str, Any]
    model_insights: dict[str, Any] | None = None
    comparisons: list[dict[str, Any]] = Field(default_factory=list)
    pathways: list[dict[str, Any]] = Field(default_factory=list)
