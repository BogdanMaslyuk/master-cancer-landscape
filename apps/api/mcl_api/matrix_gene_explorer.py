from __future__ import annotations

from functools import lru_cache
from typing import Any

import numpy as np
import pandas as pd

from .annotated_gene_explorer import AnnotatedGeneExplorerStore
from .gene_explorer import _clean, _records, _truthy

try:  # scipy is an optional analysis dependency in MCL
    from scipy.stats import mannwhitneyu, spearmanr
except Exception:  # pragma: no cover
    mannwhitneyu = None
    spearmanr = None


class MatrixGeneExplorerStore(AnnotatedGeneExplorerStore):
    """Gene Explorer with cross-context matrix and model-level hypothesis tools."""

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
            cells = {comparison_id: metric_lookup.get((symbol, comparison_id)) for comparison_id in comparison_ids}
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

    @staticmethod
    def _distribution(values: pd.Series) -> dict[str, Any]:
        numeric = pd.to_numeric(values, errors="coerce").dropna()
        if numeric.empty:
            return {"n": 0, "median": None, "q1": None, "q3": None, "min": None, "max": None, "fraction_lt_minus_0_5": None, "fraction_lt_minus_1": None}
        return {
            "n": int(len(numeric)),
            "median": float(numeric.median()),
            "q1": float(numeric.quantile(0.25)),
            "q3": float(numeric.quantile(0.75)),
            "min": float(numeric.min()),
            "max": float(numeric.max()),
            "fraction_lt_minus_0_5": float((numeric < -0.5).mean()),
            "fraction_lt_minus_1": float((numeric < -1.0).mean()),
        }

    def descriptive_contexts(self, gene_symbol: str) -> dict[str, Any]:
        symbol, _ = self._gene_or_raise(gene_symbol)
        ge = self._model_gene_layer("gene_effect", symbol)
        audit = self._audit()
        if ge.empty or audit.empty or "cancer_id" not in audit.columns:
            return {"gene_symbol": symbol, "available": False, "items": []}
        ge = ge[["model_id", "gene_effect"]].copy()
        ge["model_id"] = ge["model_id"].astype(str)
        audit = audit.copy()
        audit["model_id"] = audit["model_id"].astype(str)
        cfgs = self._contexts_config()
        items: list[dict[str, Any]] = []
        for cancer_id, group in audit.groupby("cancer_id", sort=False):
            group = group.drop_duplicates("model_id")
            assigned = group["assigned_group"].astype(str).str.lower() if "assigned_group" in group.columns else pd.Series("context", index=group.index)
            eligible = group[assigned != "excluded"]
            context = group[assigned == "context"]
            comparator = group[assigned == "comparator"]
            eligible_values = ge[ge["model_id"].isin(set(eligible["model_id"]))]["gene_effect"]
            context_values = ge[ge["model_id"].isin(set(context["model_id"]))]["gene_effect"]
            comparator_values = ge[ge["model_id"].isin(set(comparator["model_id"]))]["gene_effect"]
            if eligible_values.dropna().empty and context_values.dropna().empty:
                continue
            cfg = cfgs.get(str(cancer_id)) or {}
            items.append(
                {
                    "cancer_id": str(cancer_id),
                    "cancer_name": cfg.get("name") or str(cancer_id),
                    "context_definition": (cfg.get("depmap") or {}).get("context_definition"),
                    "context": self._distribution(context_values),
                    "comparator": self._distribution(comparator_values),
                    "eligible": self._distribution(eligible_values),
                    "has_comparator": bool(len(comparator)),
                    "interpretation_level": "descriptive_only",
                }
            )
        items.sort(key=lambda row: (row["context"]["median"] is None, row["context"]["median"] if row["context"]["median"] is not None else 999))
        return _clean({
            "gene_symbol": symbol,
            "available": bool(items),
            "items": items,
            "guardrail": "Descriptive summaries show observed CRISPR Gene Effect in MCL cell models. Context selectivity requires a valid comparator and inferential test.",
        })

    @staticmethod
    def _spearman(x: pd.Series, y: pd.Series) -> tuple[float | None, float | None, int]:
        frame = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")}).dropna()
        n = int(len(frame))
        if n < 3 or frame["x"].nunique() < 2 or frame["y"].nunique() < 2:
            return None, None, n
        if spearmanr is not None:
            result = spearmanr(frame["x"], frame["y"], nan_policy="omit")
            return float(result.statistic), float(result.pvalue), n
        rho = frame["x"].rank().corr(frame["y"].rank())
        return (float(rho) if pd.notna(rho) else None), None, n

    def correlations(self, gene_symbol: str) -> dict[str, Any]:
        symbol, _ = self._gene_or_raise(gene_symbol)
        ge = self._model_gene_layer("gene_effect", symbol)
        if ge.empty:
            return {"gene_symbol": symbol, "available": False, "relationships": [], "points": []}
        merged = ge.copy()
        for layer in ("expression", "copy_number"):
            frame = self._model_gene_layer(layer, symbol)
            if not frame.empty:
                merged = merged.merge(frame, on="model_id", how="left")
        metadata = self._model_metadata()
        if not metadata.empty and "model_id" in metadata.columns:
            keep = [c for c in ("model_id", "cell_line_name", "oncotree_lineage", "oncotree_primary_disease", "oncotree_subtype") if c in metadata.columns]
            merged = merged.merge(metadata[keep].drop_duplicates("model_id"), on="model_id", how="left")
        audit = self._audit()
        memberships: dict[str, list[str]] = {}
        if not audit.empty and {"model_id", "cancer_id"}.issubset(audit.columns):
            for model_id, group in audit.groupby("model_id"):
                memberships[str(model_id)] = sorted(set(group["cancer_id"].dropna().astype(str)))
        relationships = []
        for layer, label in (("expression", "RNA expression"), ("copy_number", "Copy number")):
            if layer not in merged.columns:
                relationships.append({"layer": layer, "label": label, "available": False, "rho": None, "p_value": None, "n": 0})
                continue
            rho, p_value, n = self._spearman(merged[layer], merged["gene_effect"])
            relationships.append({"layer": layer, "label": label, "available": n >= 3, "rho": rho, "p_value": p_value, "n": n})
        points = _records(merged.sort_values("gene_effect", ascending=True, na_position="last"))
        for point in points:
            point["cancer_ids"] = memberships.get(str(point.get("model_id")), [])
        return _clean({
            "gene_symbol": symbol,
            "available": True,
            "statistics_engine": "scipy" if spearmanr is not None else "rank-correlation fallback; p-value unavailable",
            "relationships": relationships,
            "points": points,
            "guardrails": {
                "correlation": "Exploratory association; correlation does not establish causality.",
                "expression": "RNA expression is not equivalent to dependency or protein activity.",
                "copy_number": "Relative copy number is continuous and is not automatically called amplification/deletion.",
            },
        })

    @lru_cache(maxsize=1)
    def _mutation_frame(self) -> pd.DataFrame:
        parquet = self.processed / "depmap_model_mutations.parquet"
        tsv = self.processed / "depmap_model_mutations.tsv"
        if parquet.exists():
            frame = pd.read_parquet(parquet)
        elif tsv.exists():
            frame = pd.read_csv(tsv, sep="\t", low_memory=False)
        else:
            return pd.DataFrame()
        if "model_id" in frame.columns:
            frame["model_id"] = frame["model_id"].astype(str)
        if "gene" in frame.columns:
            frame["gene"] = frame["gene"].astype(str).str.upper()
        return frame

    @staticmethod
    def _cliffs_delta(x: np.ndarray, y: np.ndarray) -> float | None:
        if len(x) == 0 or len(y) == 0:
            return None
        y_sorted = np.sort(y)
        greater = 0
        less = 0
        for value in x:
            greater += int(np.searchsorted(y_sorted, value, side="left"))
            less += int(len(y_sorted) - np.searchsorted(y_sorted, value, side="right"))
        return float((greater - less) / (len(x) * len(y)))

    @staticmethod
    def _bh_qvalues(p_values: list[float | None]) -> list[float | None]:
        valid = [(i, p) for i, p in enumerate(p_values) if p is not None and np.isfinite(p)]
        out: list[float | None] = [None] * len(p_values)
        if not valid:
            return out
        ordered = sorted(valid, key=lambda item: item[1])
        m = len(ordered)
        previous = 1.0
        for rank_from_end, (index, p) in enumerate(reversed(ordered), start=1):
            rank = m - rank_from_end + 1
            q = min(previous, float(p) * m / rank)
            previous = q
            out[index] = min(q, 1.0)
        return out

    def mutation_associations(self, gene_symbol: str, *, min_mutated: int = 3, min_wildtype: int = 3, limit: int = 30) -> dict[str, Any]:
        symbol, _ = self._gene_or_raise(gene_symbol)
        ge = self._model_gene_layer("gene_effect", symbol)
        mutations = self._mutation_frame()
        if ge.empty or mutations.empty or not {"model_id", "gene"}.issubset(mutations.columns):
            return {"gene_symbol": symbol, "available": False, "items": [], "note": "Mutation profile index is not available."}
        mutations = mutations.copy()
        if "is_functional" in mutations.columns:
            functional = mutations[mutations["is_functional"].map(_truthy)].copy()
        else:
            flags = [c for c in ("driver", "hotspot", "likely_lof", "oncogene_high_impact", "tumor_suppressor_high_impact") if c in mutations.columns]
            mask = pd.Series(False, index=mutations.index)
            for col in flags:
                mask = mask | mutations[col].map(_truthy)
            functional = mutations[mask].copy() if flags else mutations.copy()
        profiled_ids = set(mutations["model_id"].dropna().astype(str))
        ge = ge[ge["model_id"].astype(str).isin(profiled_ids)].dropna(subset=["gene_effect"]).copy()
        ge_lookup = ge.set_index("model_id")["gene_effect"]
        all_ids = set(ge_lookup.index.astype(str))
        rows: list[dict[str, Any]] = []
        for mutated_gene, group in functional.groupby("gene"):
            mutated_ids = set(group["model_id"].dropna().astype(str)) & all_ids
            wildtype_ids = all_ids - mutated_ids
            if len(mutated_ids) < min_mutated or len(wildtype_ids) < min_wildtype:
                continue
            x = pd.to_numeric(ge_lookup.loc[list(mutated_ids)], errors="coerce").dropna().to_numpy(dtype=float)
            y = pd.to_numeric(ge_lookup.loc[list(wildtype_ids)], errors="coerce").dropna().to_numpy(dtype=float)
            if len(x) < min_mutated or len(y) < min_wildtype:
                continue
            p_value = None
            if mannwhitneyu is not None:
                try:
                    p_value = float(mannwhitneyu(x, y, alternative="two-sided").pvalue)
                except ValueError:
                    p_value = None
            rows.append({
                "mutation_gene": str(mutated_gene),
                "mutated_n": int(len(x)),
                "wildtype_n": int(len(y)),
                "mutated_median_gene_effect": float(np.median(x)),
                "wildtype_median_gene_effect": float(np.median(y)),
                "delta_median_gene_effect": float(np.median(x) - np.median(y)),
                "cliffs_delta": self._cliffs_delta(x, y),
                "p_value": p_value,
                "functional_variant_rows_n": int(len(group)),
            })
        q_values = self._bh_qvalues([row["p_value"] for row in rows])
        for row, q_value in zip(rows, q_values):
            row["q_value"] = q_value
        rows.sort(key=lambda row: (row["q_value"] is None, row["q_value"] if row["q_value"] is not None else 1.0, -abs(row["delta_median_gene_effect"])))
        rows = rows[: max(1, min(int(limit), 100))]
        return _clean({
            "gene_symbol": symbol,
            "available": True,
            "profiled_models_n": int(len(all_ids)),
            "min_mutated": int(min_mutated),
            "min_wildtype": int(min_wildtype),
            "statistics_engine": "scipy" if mannwhitneyu is not None else "descriptive fallback; p/q unavailable",
            "items": rows,
            "interpretation": {
                "delta": "mutated median Gene Effect - mutation-negative median Gene Effect; more negative means stronger dependency in the mutated group",
                "mutation_definition": "functional mutation rows from the local DepMap mutation index",
                "guardrail": "Exploratory cell-model association screen; not proof of causality or patient response.",
            },
        })

    def identity(self, gene_symbol: str) -> dict[str, Any]:
        payload = super().identity(gene_symbol)
        symbol = str(payload["identity"]["gene_symbol"]).upper()
        payload["model_insights"] = {
            "descriptive_contexts": self.descriptive_contexts(symbol),
            "correlations": self.correlations(symbol),
            "mutation_associations": self.mutation_associations(symbol, limit=30),
        }
        return payload
