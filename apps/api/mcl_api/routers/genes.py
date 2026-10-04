from __future__ import annotations

from fastapi import APIRouter, Query

from ..api_utils import guard
from ..schemas.genes import (
    GeneDetailResponse,
    GeneFacetResponse,
    GeneMatrixResponse,
    GeneSearchResponse,
)
from ..state import gene_service


router = APIRouter()


@router.get("/api/genes")
def genes(search: str | None = None, limit: int = Query(200, ge=1, le=2000)):
    return guard(lambda: gene_service.genes(search, limit))


@router.get("/api/genes/suggest")
def gene_suggest(q: str = Query(..., min_length=1), limit: int = Query(12, ge=1, le=30)):
    return guard(lambda: gene_service.suggest(q.strip(), limit))


@router.get("/api/genes/facets", response_model=GeneFacetResponse)
def gene_facets():
    return guard(gene_service.facets)


@router.get("/api/genes/search", response_model=GeneSearchResponse)
def gene_search(
    q: str | None = None,
    domain: str | None = None,
    subdomain: str | None = None,
    pathway: str | None = None,
    annotation_source: str | None = None,
    protein_class: str | None = None,
    compartment: str | None = None,
    hallmark: str | None = None,
    cancer_id: str | None = None,
    comparison_id: str | None = None,
    gene_effect_max: float | None = None,
    delta_gene_effect_max: float | None = None,
    q_value_max: float | None = None,
    cliffs_delta_abs_min: float | None = None,
    stable_only: bool = False,
    exclude_broad: bool = False,
    exclude_low_sample: bool = False,
    page: int = Query(1, ge=1),
    page_size: int = Query(100, ge=1, le=250),
    sort_by: str = "best_delta_gene_effect",
    sort_order: str = "asc",
):
    return guard(
        lambda: gene_service.search(
            q,
            domain,
            subdomain,
            pathway,
            annotation_source,
            protein_class,
            compartment,
            hallmark,
            cancer_id,
            comparison_id,
            gene_effect_max,
            delta_gene_effect_max,
            q_value_max,
            cliffs_delta_abs_min,
            stable_only,
            exclude_broad,
            exclude_low_sample,
            page,
            page_size,
            sort_by,
            sort_order,
        )
    )


@router.get("/api/gene-matrix", response_model=GeneMatrixResponse)
def gene_matrix(
    q: str | None = None,
    domain: str | None = None,
    subdomain: str | None = None,
    protein_class: str | None = None,
    compartment: str | None = None,
    hallmark: str | None = None,
    cancer_id: str | None = None,
    gene_effect_max: float | None = None,
    delta_gene_effect_max: float | None = None,
    q_value_max: float | None = None,
    cliffs_delta_abs_min: float | None = None,
    stable_only: bool = False,
    exclude_broad: bool = True,
    exclude_low_sample: bool = True,
    limit: int = Query(60, ge=1, le=120),
    sort_by: str = "best_delta_gene_effect",
    sort_order: str = "asc",
):
    return guard(
        lambda: gene_service.matrix(
            q,
            domain,
            subdomain,
            protein_class,
            compartment,
            hallmark,
            cancer_id,
            gene_effect_max,
            delta_gene_effect_max,
            q_value_max,
            cliffs_delta_abs_min,
            stable_only,
            exclude_broad,
            exclude_low_sample,
            limit,
            sort_by,
            sort_order,
        )
    )


@router.get("/api/genes/stable")
def stable_genes():
    return guard(gene_service.stable_genes)


@router.get("/api/genes/{gene_symbol}/contexts")
def gene_contexts(gene_symbol: str):
    symbol = gene_symbol.strip().upper()
    return guard(lambda: gene_service.contexts(symbol))


@router.get("/api/genes/{gene_symbol}/models")
def gene_models(
    gene_symbol: str,
    cancer_id: str | None = None,
    gene_effect_max: float | None = None,
    page: int = Query(1, ge=1),
    page_size: int = Query(100, ge=1, le=500),
    sort_by: str = "gene_effect",
    sort_order: str = "asc",
):
    symbol = gene_symbol.strip().upper()
    return guard(
        lambda: gene_service.models(
            symbol,
            cancer_id,
            gene_effect_max,
            page,
            page_size,
            sort_by,
            sort_order,
        )
    )


@router.get("/api/genes/{gene_symbol}/annotations")
def gene_annotations(gene_symbol: str):
    symbol = gene_symbol.strip().upper()
    return guard(lambda: gene_service.annotations(symbol))


@router.get("/api/genes/{gene_symbol}/mutation-associations")
def gene_mutation_associations(
    gene_symbol: str,
    limit: int = Query(30, ge=1, le=100),
):
    symbol = gene_symbol.strip().upper()
    return guard(lambda: gene_service.mutation_associations(symbol, limit))


@router.get("/api/genes/{gene_symbol}", response_model=GeneDetailResponse)
def gene(gene_symbol: str):
    symbol = gene_symbol.strip().upper()
    return guard(lambda: gene_service.detail(symbol))
