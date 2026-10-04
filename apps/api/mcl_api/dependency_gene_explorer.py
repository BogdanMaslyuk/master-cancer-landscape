from __future__ import annotations

from functools import lru_cache
from typing import Any

import pandas as pd

from .fast_runtime_gene_explorer import FastRuntimeGeneExplorerStore
from .gene_explorer import _bool_series, _clean, _number_series, _records


class DependencyAwareGeneExplorerStore(FastRuntimeGeneExplorerStore):
    """Gene Explorer enriched with CRISPR dependency prevalence and specificity.

    The dependency summary is precomputed at build time from the complete CRISPR
    model universe. Interactive requests only read the compact summary and never
    rescan the wide Gene Effect matrix.
    """

    @lru_cache(maxsize=1)
    def _dependency_summary(self) -> pd.DataFrame:
        path = self._runtime_path("gene_dependency_summary.parquet")
        frame = pd.read_parquet(path)
        if "gene_symbol" in frame.columns:
            frame["gene_symbol"] = frame["gene_symbol"].astype(str).str.upper()
        return frame

    @lru_cache(maxsize=1)
    def catalog_frame(self) -> pd.DataFrame:
        base = super().catalog_frame().copy()
        summary = self._dependency_summary()
        if base.empty or summary.empty:
            return base
        overlap = [c for c in summary.columns if c != "gene_symbol" and c in base.columns]
        if overlap:
            base = base.drop(columns=overlap)
        return base.merge(summary.drop_duplicates("gene_symbol"), on="gene_symbol", how="left")

    def search(
        self,
        *,
        query: str | None = None,
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
        dependency_type: str | None = None,
        dependency_fraction_min: float | None = None,
        specificity_score_min: float | None = None,
        page: int = 1,
        page_size: int = 100,
        sort_by: str = "specificity_score",
        sort_order: str = "desc",
    ) -> dict[str, Any]:
        catalog = self.catalog_frame().copy()

        if query:
            q = str(query).strip().casefold()
            mask = catalog["gene_symbol"].astype(str).str.casefold().str.contains(q, regex=False)
            if "gene_name" in catalog.columns:
                mask = mask | catalog["gene_name"].fillna("").astype(str).str.casefold().str.contains(q, regex=False)
            if "aliases_json" in catalog.columns:
                mask = mask | catalog["aliases_json"].fillna("").astype(str).str.casefold().str.contains(q, regex=False)
            catalog = catalog[mask]

        for value, column in (
            (domain, "mcl_domains_json"),
            (subdomain, "mcl_subdomains_json"),
            (protein_class, "protein_classes_json"),
            (compartment, "compartments_json"),
            (hallmark, "hallmarks_json"),
        ):
            if value and column in catalog.columns:
                catalog = catalog[catalog[column].map(lambda item, needle=value: self._json_contains(item, needle))]

        if pathway or annotation_source:
            formal = self.formal_annotation_frame().copy()
            if pathway:
                needle = str(pathway).casefold()
                term_name = formal.get("term_name", pd.Series("", index=formal.index)).astype(str).str.casefold()
                term_id = formal.get("term_id", pd.Series("", index=formal.index)).astype(str).str.casefold()
                formal = formal[term_name.str.contains(needle, regex=False) | term_id.str.contains(needle, regex=False)]
            if annotation_source:
                formal = formal[formal["source"].astype(str).str.upper() == str(annotation_source).upper()]
            allowed = set(formal["gene_symbol"].astype(str).str.upper()) if not formal.empty else set()
            catalog = catalog[catalog["gene_symbol"].astype(str).str.upper().isin(allowed)]

        context_filter_active = any(
            value is not None and value is not False
            for value in (
                cancer_id,
                comparison_id,
                delta_gene_effect_max,
                q_value_max,
                cliffs_delta_abs_min,
                exclude_broad,
                exclude_low_sample,
            )
        )
        if context_filter_active:
            metrics = self.context_metrics_frame().copy()
            if cancer_id:
                metrics = metrics[metrics["cancer_id"].astype(str) == str(cancer_id)]
            if comparison_id:
                metrics = metrics[metrics["comparison_id"].astype(str) == str(comparison_id)]
            if delta_gene_effect_max is not None:
                metrics = metrics[_number_series(metrics, "delta_gene_effect") <= float(delta_gene_effect_max)]
            if q_value_max is not None:
                q_values = _number_series(metrics, "q_value")
                metrics = metrics[q_values.notna() & (q_values <= float(q_value_max))]
            if cliffs_delta_abs_min is not None:
                cliffs = _number_series(metrics, "cliffs_delta").abs()
                metrics = metrics[cliffs.notna() & (cliffs >= float(cliffs_delta_abs_min))]
            if exclude_broad and "broad_dependency_warning" in metrics.columns:
                metrics = metrics[~_bool_series(metrics["broad_dependency_warning"])]
            if exclude_low_sample and "low_sample_size" in metrics.columns:
                metrics = metrics[~_bool_series(metrics["low_sample_size"])]
            allowed = set(metrics["gene_symbol"].astype(str).str.upper()) if not metrics.empty else set()
            catalog = catalog[catalog["gene_symbol"].astype(str).str.upper().isin(allowed)]

        if gene_effect_max is not None and "best_model_gene_effect" in catalog.columns:
            values = _number_series(catalog, "best_model_gene_effect")
            catalog = catalog[values.notna() & (values <= float(gene_effect_max))]
        if stable_only and "present_all_thresholds" in catalog.columns:
            catalog = catalog[_bool_series(catalog["present_all_thresholds"])]

        if dependency_type and "dependency_type" in catalog.columns:
            catalog = catalog[catalog["dependency_type"].astype(str) == str(dependency_type)]
        if dependency_fraction_min is not None and "dependency_fraction" in catalog.columns:
            values = _number_series(catalog, "dependency_fraction")
            catalog = catalog[values.notna() & (values >= float(dependency_fraction_min))]
        if specificity_score_min is not None and "specificity_score" in catalog.columns:
            values = _number_series(catalog, "specificity_score")
            catalog = catalog[values.notna() & (values >= float(specificity_score_min))]

        allowed_sort = {
            "gene_symbol", "gene_name", "dependency_fraction", "dependency_models_n",
            "specificity_score", "best_cancer_dependency_fraction", "global_median_gene_effect",
            "best_delta_gene_effect", "best_model_gene_effect", "best_q_value", "best_cliffs_delta",
            "comparisons_n", "significant_comparisons_n", "functional_annotations_n",
        }
        if sort_by not in allowed_sort or sort_by not in catalog.columns:
            sort_by = "specificity_score" if "specificity_score" in catalog.columns else "gene_symbol"
        ascending = str(sort_order).lower() != "desc"
        catalog = catalog.sort_values(sort_by, ascending=ascending, na_position="last")

        total = int(len(catalog))
        page = max(1, int(page))
        page_size = max(1, min(int(page_size), 250))
        start = (page - 1) * page_size
        page_frame = catalog.iloc[start : start + page_size].copy()
        return _clean(
            {
                "page": page,
                "page_size": page_size,
                "total": total,
                "pages": max(1, (total + page_size - 1) // page_size),
                "sort_by": sort_by,
                "sort_order": "asc" if ascending else "desc",
                "functional_coverage": self.functional_coverage(),
                "reference_coverage": self.reference_coverage(),
                "items": _records(page_frame),
            }
        )
