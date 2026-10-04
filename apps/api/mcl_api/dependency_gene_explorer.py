from __future__ import annotations

from functools import lru_cache
from typing import Any

import numpy as np
import pandas as pd

from .fast_runtime_gene_explorer import FastRuntimeGeneExplorerStore
from .gene_explorer import _bool_series, _clean, _number_series, _records


DEPENDENCY_THRESHOLD = -0.5
MIN_CANCER_MODELS = 5


class DependencyAwareGeneExplorerStore(FastRuntimeGeneExplorerStore):
    """Gene Explorer enriched with interpretable CRISPR dependency summaries."""

    @lru_cache(maxsize=1)
    def _dependency_summary(self) -> pd.DataFrame:
        path = self._runtime_path("gene_dependency_summary.parquet")
        frame = pd.read_parquet(path)
        if "gene_symbol" in frame.columns:
            frame["gene_symbol"] = frame["gene_symbol"].astype(str).str.upper()
        return frame

    @lru_cache(maxsize=1)
    def _crispr_atlas_frame(self) -> pd.DataFrame:
        path = self.processed / "depmap_crispr_model_atlas.parquet"
        if not path.exists():
            return pd.DataFrame()
        frame = pd.read_parquet(path)
        if "model_id" in frame.columns:
            frame["model_id"] = frame["model_id"].astype(str)
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

    @staticmethod
    def _strength_label(rho: float | None) -> str:
        if rho is None or not np.isfinite(rho):
            return "не определена"
        value = abs(float(rho))
        if value < 0.10:
            return "практически отсутствует"
        if value < 0.20:
            return "очень слабая"
        if value < 0.40:
            return "слабая"
        if value < 0.60:
            return "умеренная"
        if value < 0.80:
            return "сильная"
        return "очень сильная"

    @classmethod
    def _relationship_interpretation(cls, layer: str, rho: float | None, p_value: float | None, n: int) -> dict[str, Any]:
        strength = cls._strength_label(rho)
        if rho is None or not np.isfinite(rho):
            direction = "направление не определено"
            sentence = "Недостаточно вариабельности или данных, чтобы оценить связь."
        elif rho < 0:
            direction = "больше показателя → сильнее зависимость"
            subject = "экспрессия" if layer == "expression" else "относительное число копий"
            sentence = f"Чем выше {subject}, тем в среднем более отрицателен Gene Effect, то есть зависимость несколько сильнее."
        elif rho > 0:
            direction = "больше показателя → слабее зависимость"
            subject = "экспрессия" if layer == "expression" else "относительное число копий"
            sentence = f"Чем выше {subject}, тем в среднем менее отрицателен Gene Effect, то есть зависимость несколько слабее."
        else:
            direction = "направленной связи нет"
            sentence = "Монотонная связь практически отсутствует."

        statistically_clear = p_value is not None and np.isfinite(p_value) and float(p_value) < 0.05
        if strength in {"практически отсутствует", "очень слабая", "слабая"}:
            conclusion = sentence + " Связь по величине невелика и сама по себе плохо объясняет различия между моделями."
        else:
            conclusion = sentence + " Связь достаточно выражена, чтобы рассматривать этот показатель как один из кандидатов на объяснение межмодельных различий."
        if statistically_clear and n >= 100:
            conclusion += " Статистическая уверенность высокая, но она не заменяет оценку величины эффекта."
        conclusion += " Корреляция не доказывает причинный механизм."
        return {
            "strength": strength,
            "direction": direction,
            "statistically_clear": statistically_clear,
            "conclusion_ru": conclusion,
        }

    def correlations(self, gene_symbol: str) -> dict[str, Any]:
        payload = super().correlations(gene_symbol)
        relationships = payload.get("relationships") or []
        points = pd.DataFrame(payload.get("points") or [])
        atlas = self._crispr_atlas_frame()

        if not points.empty and not atlas.empty and "model_id" in points.columns:
            meta_cols = [c for c in ("model_id", "mcl_cancer_id", "mcl_cancer_name", "mcl_organ_ru") if c in atlas.columns]
            if meta_cols:
                points["model_id"] = points["model_id"].astype(str)
                points = points.merge(atlas[meta_cols].drop_duplicates("model_id"), on="model_id", how="left", suffixes=("", "_atlas"))

        for relationship in relationships:
            layer = str(relationship.get("layer") or "")
            rho = relationship.get("rho")
            p_value = relationship.get("p_value")
            n = int(relationship.get("n") or 0)
            relationship["interpretation"] = self._relationship_interpretation(layer, rho, p_value, n)
            within: list[dict[str, Any]] = []
            if layer in {"expression", "copy_number"} and not points.empty and "mcl_cancer_id" in points.columns:
                for cancer_id, group in points.dropna(subset=["mcl_cancer_id"]).groupby("mcl_cancer_id", sort=False):
                    if len(group) < 8:
                        continue
                    local_rho, local_p, local_n = self._spearman(group[layer], group["gene_effect"])
                    if local_rho is None or local_n < 8:
                        continue
                    first = group.iloc[0]
                    within.append({
                        "cancer_id": str(cancer_id),
                        "cancer_name": first.get("mcl_cancer_name"),
                        "organ_ru": first.get("mcl_organ_ru"),
                        "rho": local_rho,
                        "p_value": local_p,
                        "n": local_n,
                        "strength": self._strength_label(local_rho),
                    })
                within.sort(key=lambda row: abs(float(row.get("rho") or 0)), reverse=True)
            relationship["within_cancers"] = within[:8]

        payload["relationships"] = relationships
        payload["points"] = _records(points) if not points.empty else payload.get("points", [])
        return _clean(payload)

    def dependency_landscape(self, gene_symbol: str) -> dict[str, Any]:
        symbol, catalog_row = self._gene_or_raise(gene_symbol)
        ge = self._model_gene_layer("gene_effect", symbol)
        atlas = self._crispr_atlas_frame()
        if ge.empty or atlas.empty or "model_id" not in atlas.columns:
            return {"gene_symbol": symbol, "available": False, "summary": {}, "cancers": []}

        ge = ge[["model_id", "gene_effect"]].copy()
        ge["model_id"] = ge["model_id"].astype(str)
        merged = ge.merge(atlas, on="model_id", how="inner")
        merged["gene_effect"] = pd.to_numeric(merged["gene_effect"], errors="coerce")
        merged = merged.dropna(subset=["gene_effect"])
        if merged.empty:
            return {"gene_symbol": symbol, "available": False, "summary": {}, "cancers": []}

        merged["dependent"] = merged["gene_effect"] <= DEPENDENCY_THRESHOLD
        global_n = int(len(merged))
        global_dep_n = int(merged["dependent"].sum())
        global_fraction = float(global_dep_n / global_n) if global_n else 0.0
        global_median = float(merged["gene_effect"].median())

        rows: list[dict[str, Any]] = []
        group_col = "mcl_cancer_id" if "mcl_cancer_id" in merged.columns else "mcl_cancer_name"
        for cancer_id, group in merged.dropna(subset=[group_col]).groupby(group_col, sort=False):
            values = group["gene_effect"].dropna()
            if values.empty:
                continue
            n = int(len(values))
            dependent_n = int((values <= DEPENDENCY_THRESHOLD).sum())
            fraction = float(dependent_n / n)
            first = group.iloc[0]
            rows.append({
                "cancer_id": str(cancer_id),
                "cancer_name": first.get("mcl_cancer_name") or cancer_id,
                "organ_ru": first.get("mcl_organ_ru"),
                "system_ru": first.get("mcl_system_ru"),
                "models_n": n,
                "dependent_models_n": dependent_n,
                "dependent_fraction": fraction,
                "median_gene_effect": float(values.median()),
                "q1_gene_effect": float(values.quantile(0.25)),
                "q3_gene_effect": float(values.quantile(0.75)),
                "specificity_vs_global": float(fraction - global_fraction),
                "eligible_for_label": n >= MIN_CANCER_MODELS,
            })
        rows.sort(key=lambda row: (row["eligible_for_label"], row["dependent_fraction"], -row["median_gene_effect"]), reverse=True)

        best = next((row for row in rows if row["eligible_for_label"]), None)
        summary = {
            "models_n": global_n,
            "dependent_models_n": global_dep_n,
            "dependent_fraction": global_fraction,
            "global_median_gene_effect": global_median,
            "dependency_threshold": DEPENDENCY_THRESHOLD,
            "dependency_type": catalog_row.get("dependency_type"),
            "dependency_type_ru": catalog_row.get("dependency_type_ru"),
            "specificity_score": catalog_row.get("specificity_score"),
            "specificity_label_ru": catalog_row.get("specificity_label_ru"),
            "best_cancer_id": best.get("cancer_id") if best else catalog_row.get("best_cancer_id"),
            "best_cancer_name": best.get("cancer_name") if best else catalog_row.get("best_cancer_name"),
            "best_cancer_organ_ru": best.get("organ_ru") if best else catalog_row.get("best_cancer_organ_ru"),
            "best_cancer_dependency_fraction": best.get("dependent_fraction") if best else catalog_row.get("best_cancer_dependency_fraction"),
        }
        return _clean({
            "gene_symbol": symbol,
            "available": True,
            "summary": summary,
            "cancers": rows,
            "interpretation": {
                "dependency_definition": f"Рабочий навигационный порог MCL: Gene Effect ≤ {DEPENDENCY_THRESHOLD}.",
                "specificity": "Специфичность — описательное превышение доли зависимых моделей в опухолевой группе над общей долей среди CRISPR-моделей.",
                "minimum_group": f"Автоматический ярлык опухоли формируется только для групп с n ≥ {MIN_CANCER_MODELS}.",
                "guardrail": "Это функциональная зависимость после CRISPR-выключения в клеточных моделях, а не доказательство эффективности лекарственного ингибитора у пациентов.",
            },
        })

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

        context_filter_active = any(value is not None and value is not False for value in (
            cancer_id, comparison_id, delta_gene_effect_max, q_value_max, cliffs_delta_abs_min, exclude_broad, exclude_low_sample,
        ))
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
            "gene_symbol", "gene_name", "dependency_fraction", "dependency_models_n", "specificity_score",
            "best_cancer_dependency_fraction", "global_median_gene_effect", "best_delta_gene_effect", "best_model_gene_effect",
            "best_q_value", "best_cliffs_delta", "comparisons_n", "significant_comparisons_n", "functional_annotations_n",
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
        return _clean({
            "page": page,
            "page_size": page_size,
            "total": total,
            "pages": max(1, (total + page_size - 1) // page_size),
            "sort_by": sort_by,
            "sort_order": "asc" if ascending else "desc",
            "functional_coverage": self.functional_coverage(),
            "reference_coverage": self.reference_coverage(),
            "items": _records(page_frame),
        })
