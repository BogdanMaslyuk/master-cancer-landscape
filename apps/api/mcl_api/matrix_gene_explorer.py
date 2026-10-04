from __future__ import annotations

from typing import Any

import pandas as pd

from .annotated_gene_explorer import AnnotatedGeneExplorerStore
from .gene_explorer import _clean, _truthy


class MatrixGeneExplorerStore(AnnotatedGeneExplorerStore):
    """Adds a compact Gene × MCL comparison matrix for cross-context inspection."""

    def matrix(
        self,
        *,
        query: str | None = None,
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
        limit: int = 60,
        sort_by: str = "best_delta_gene_effect",
        sort_order: str = "asc",
    ) -> dict[str, Any]:
        limit = max(1, min(int(limit), 120))
        selected = self.search(
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
            page=1,
            page_size=limit,
            sort_by=sort_by,
            sort_order=sort_order,
        )
        catalog_rows = selected.get("items") or []
        genes = [str(row.get("gene_symbol") or "").upper() for row in catalog_rows if row.get("gene_symbol")]

        specs = self.store.comparison_specs()
        if cancer_id:
            specs = [row for row in specs if str(row.get("cancer_id")) == str(cancer_id)]
        comparison_ids = [str(row.get("id")) for row in specs if row.get("id")]

        metrics = self.context_metrics_frame().copy()
        if metrics.empty or not genes:
            metric_lookup: dict[tuple[str, str], dict[str, Any]] = {}
        else:
            metrics = metrics[
                metrics["gene_symbol"].astype(str).str.upper().isin(set(genes))
                & metrics["comparison_id"].astype(str).isin(set(comparison_ids))
            ].copy()
            metric_lookup = {}
            for _, row in metrics.iterrows():
                key = (str(row.get("gene_symbol") or "").upper(), str(row.get("comparison_id") or ""))
                metric_lookup[key] = {
                    "delta_gene_effect": row.get("delta_gene_effect"),
                    "context_median_gene_effect": row.get("context_median_gene_effect"),
                    "comparator_median_gene_effect": row.get("comparator_median_gene_effect"),
                    "q_value": row.get("q_value"),
                    "cliffs_delta": row.get("cliffs_delta"),
                    "fdr_0_05": _truthy(row.get("fdr_0_05")),
                    "broad_dependency_warning": _truthy(row.get("broad_dependency_warning")),
                    "low_sample_size": _truthy(row.get("low_sample_size")),
                    "context_models_n": row.get("context_models_n"),
                    "comparator_models_n": row.get("comparator_models_n"),
                }

        rows: list[dict[str, Any]] = []
        for catalog_row in catalog_rows:
            symbol = str(catalog_row.get("gene_symbol") or "").upper()
            cells = {
                comparison_id: metric_lookup.get((symbol, comparison_id))
                for comparison_id in comparison_ids
            }
            rows.append(
                {
                    "gene_symbol": symbol,
                    "gene_name": catalog_row.get("gene_name"),
                    "present_all_thresholds": _truthy(catalog_row.get("present_all_thresholds")),
                    "best_delta_gene_effect": catalog_row.get("best_delta_gene_effect"),
                    "best_model_gene_effect": catalog_row.get("best_model_gene_effect"),
                    "mcl_domains_json": catalog_row.get("mcl_domains_json"),
                    "cells": cells,
                }
            )

        comparisons = [
            {
                "id": row.get("id"),
                "label": row.get("label"),
                "cancer_id": row.get("cancer_id"),
                "context_models_n": row.get("context_models_n"),
                "comparator_models_n": row.get("comparator_models_n"),
                "qc_status": row.get("qc_status"),
            }
            for row in specs
        ]
        return _clean(
            {
                "genes_total_after_filters": selected.get("total", 0),
                "genes_returned_n": len(rows),
                "comparisons_n": len(comparisons),
                "comparisons": comparisons,
                "rows": rows,
                "filters": {
                    "query": query,
                    "domain": domain,
                    "subdomain": subdomain,
                    "protein_class": protein_class,
                    "compartment": compartment,
                    "hallmark": hallmark,
                    "cancer_id": cancer_id,
                    "gene_effect_max": gene_effect_max,
                    "delta_gene_effect_max": delta_gene_effect_max,
                    "q_value_max": q_value_max,
                    "cliffs_delta_abs_min": cliffs_delta_abs_min,
                    "stable_only": stable_only,
                    "exclude_broad": exclude_broad,
                    "exclude_low_sample": exclude_low_sample,
                },
                "interpretation": {
                    "cell_value": "delta_gene_effect = median(context) - median(comparator)",
                    "negative_delta": "More negative values indicate stronger dependency in the target/context group relative to comparator.",
                    "guardrail": "The matrix summarizes CRISPR genetic dependency comparisons, not drug response or patient prevalence.",
                },
            }
        )
