from __future__ import annotations

from functools import lru_cache

from ..repositories.gene_repository import GeneRepository


class GeneService:
    """Cached application service for Gene Explorer use-cases.

    Routers delegate orchestration here so HTTP concerns stay separate from data
    access and scientific runtime stores. Cache sizes mirror the previous router
    implementation to preserve behavior while making the boundary explicit.
    """

    def __init__(self, repository: GeneRepository):
        self.repository = repository

    @lru_cache(maxsize=256)
    def genes(self, search: str | None, limit: int):
        return self.repository.genes(search, limit)

    @lru_cache(maxsize=1)
    def stable_genes(self):
        return self.repository.stable_genes()

    @lru_cache(maxsize=512)
    def detail(self, gene_symbol: str):
        return self.repository.detail(gene_symbol)

    @lru_cache(maxsize=512)
    def suggest(self, query: str, limit: int):
        return self.repository.suggest(query, limit)

    @lru_cache(maxsize=1)
    def facets(self):
        return self.repository.facets()

    @lru_cache(maxsize=4096)
    def search(
        self,
        query: str | None,
        domain: str | None,
        subdomain: str | None,
        pathway: str | None,
        annotation_source: str | None,
        protein_class: str | None,
        compartment: str | None,
        hallmark: str | None,
        cancer_id: str | None,
        comparison_id: str | None,
        gene_effect_max: float | None,
        delta_gene_effect_max: float | None,
        q_value_max: float | None,
        cliffs_delta_abs_min: float | None,
        stable_only: bool,
        exclude_broad: bool,
        exclude_low_sample: bool,
        page: int,
        page_size: int,
        sort_by: str,
        sort_order: str,
    ):
        return self.repository.search(
            query=query,
            domain=domain,
            subdomain=subdomain,
            pathway=pathway,
            annotation_source=annotation_source,
            protein_class=protein_class,
            compartment=compartment,
            hallmark=hallmark,
            cancer_id=cancer_id,
            comparison_id=comparison_id,
            gene_effect_max=gene_effect_max,
            delta_gene_effect_max=delta_gene_effect_max,
            q_value_max=q_value_max,
            cliffs_delta_abs_min=cliffs_delta_abs_min,
            stable_only=stable_only,
            exclude_broad=exclude_broad,
            exclude_low_sample=exclude_low_sample,
            page=page,
            page_size=page_size,
            sort_by=sort_by,
            sort_order=sort_order,
        )

    @lru_cache(maxsize=1024)
    def matrix(
        self,
        query: str | None,
        domain: str | None,
        subdomain: str | None,
        protein_class: str | None,
        compartment: str | None,
        hallmark: str | None,
        cancer_id: str | None,
        gene_effect_max: float | None,
        delta_gene_effect_max: float | None,
        q_value_max: float | None,
        cliffs_delta_abs_min: float | None,
        stable_only: bool,
        exclude_broad: bool,
        exclude_low_sample: bool,
        limit: int,
        sort_by: str,
        sort_order: str,
    ):
        return self.repository.matrix(
            query=query,
            domain=domain,
            subdomain=subdomain,
            protein_class=protein_class,
            compartment=compartment,
            hallmark=hallmark,
            cancer_id=cancer_id,
            gene_effect_max=gene_effect_max,
            delta_gene_effect_max=delta_gene_effect_max,
            q_value_max=q_value_max,
            cliffs_delta_abs_min=cliffs_delta_abs_min,
            stable_only=stable_only,
            exclude_broad=exclude_broad,
            exclude_low_sample=exclude_low_sample,
            limit=limit,
            sort_by=sort_by,
            sort_order=sort_order,
        )

    @lru_cache(maxsize=512)
    def contexts(self, gene_symbol: str):
        return self.repository.contexts(gene_symbol)

    @lru_cache(maxsize=2048)
    def models(
        self,
        gene_symbol: str,
        cancer_id: str | None,
        gene_effect_max: float | None,
        page: int,
        page_size: int,
        sort_by: str,
        sort_order: str,
    ):
        return self.repository.models(
            gene_symbol,
            cancer_id=cancer_id,
            gene_effect_max=gene_effect_max,
            page=page,
            page_size=page_size,
            sort_by=sort_by,
            sort_order=sort_order,
        )

    @lru_cache(maxsize=512)
    def annotations(self, gene_symbol: str):
        return self.repository.annotations(gene_symbol)

    @lru_cache(maxsize=512)
    def mutation_associations(self, gene_symbol: str, limit: int):
        return self.repository.mutation_associations(gene_symbol, limit=limit)
