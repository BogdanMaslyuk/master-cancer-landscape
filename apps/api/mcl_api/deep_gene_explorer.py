from __future__ import annotations

from functools import lru_cache
from typing import Any

import numpy as np
import pandas as pd

from .gene_explorer import _clean, _number_series, _records, _truthy
from .matrix_gene_explorer import MatrixGeneExplorerStore

try:  # scipy is an optional analysis dependency in MCL
    from scipy.stats import mannwhitneyu, spearmanr
except Exception:  # pragma: no cover - graceful fallback when analysis extras are absent
    mannwhitneyu = None
    spearmanr = None


class DeepGeneExplorerStore(MatrixGeneExplorerStore):
    """Adds model-level descriptive and hypothesis-generation analyses to Gene Explorer.

    These methods deliberately separate descriptive cell-model summaries from validated
    target-vs-comparator inference. They never reinterpret cell-line observations as
    patient prevalence or pharmacological response.
    """

    @staticmethod
    def _summary(values: pd.Series) -> dict[str, Any]:
        numeric = pd.to_numeric(values, errors="coerce").dropna()
        if numeric.empty:
            return {
                "n": 0,
                "median": None,
                "q1": None,
                "q3": None,
                "min": None,
                "max": None,
                "fraction_lt_minus_0_5": None,
                "fraction_lt_minus_1": None,
            }
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
        gene_effect = self._model_gene_layer("gene_effect", symbol)
        audit = self._audit()
        if gene_effect.empty or audit.empty or "cancer_id" not in audit.columns:
            return {
                "gene_symbol": symbol,
                "available": False,
                "items": [],
                "note": "Model-level Gene Effect or MCL context audit is not available.",
            }

        ge = gene_effect[["model_id", "gene_effect"]].copy()
        ge["model_id"] = ge["model_id"].astype(str)
        audit = audit.copy()
        audit["model_id"] = audit["model_id"].astype(str)
        rows: list[dict[str, Any]] = []
        contexts = self._contexts_config()

        for cancer_id, group in audit.groupby("cancer_id", sort=False):
            group = group.drop_duplicates("model_id")
            eligible = group[group.get("assigned_group", "").astype(str).str.lower() != "excluded"] if "assigned_group" in group.columns else group
            context = group[group.get("assigned_group", "").astype(str).str.lower() == "context"] if "assigned_group" in group.columns else group
            comparator = group[group.get("assigned_group", "").astype(str).str.lower() == "comparator"] if "assigned_group" in group.columns else group.iloc[0:0]

            eligible_values = ge[ge["model_id"].isin(set(eligible["model_id"]))]["gene_effect"]
            context_values = ge[ge["model_id"].isin(set(context["model_id"]))]["gene_effect"]
            comparator_values = ge[ge["model_id"].isin(set(comparator["model_id"]))]["gene_effect"]
            if eligible_values.dropna().empty and context_values.dropna().empty:
                continue
            cfg = contexts.get(str(cancer_id)) or {}
            rows.append(
                {
                    "cancer_id": str(cancer_id),
                    "cancer_name": cfg.get("name") or str(cancer_id),
                    "context_definition": (cfg.get("depmap") or {}).get("context_definition"),
                    "eligible": self._summary(eligible_values),
                    "context": self._summary(context_values),
                    "comparator": self._summary(comparator_values),
                    "has_comparator": bool(len(comparator) > 0),
                    "interpretation_level": "descriptive_only",
                }
            )

        rows.sort(
            key=lambda row: (
                row["context"]["median"] is None,
                row["context"]["median"] if row["context"]["median"] is not None else 999,
            )
        )
        return _clean(
            {
                "gene_symbol": symbol,
                "available": bool(rows),
                "items": rows,
                "guardrail": "These are descriptive Gene Effect summaries for MCL cell models. They do not establish context selectivity without a pre-specified comparator and inferential test.",
            }
        )

    @staticmethod
    def _spearman(x: pd.Series, y: pd.Series) -> tuple[float | None, float | None, int]:
        frame = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")}).dropna()
        n = int(len(frame))
        if n < 3 or frame["x"].nunique() < 2 or frame["y"].nunique() < 2:
            return None, None, n
        if spearmanr is not None:
            result = spearmanr(frame["x"], frame["y"], nan_policy="omit")
            return float(result.statistic), float(result.pvalue), n
        rho = frame["x"].rank().corr(frame["y"].rank(), method="pearson")
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

        relationships: list[dict[str, Any]] = []
        for layer, label in (("expression", "RNA expression"), ("copy_number", "Copy number")):
            if layer not in merged.columns:
                relationships.append({"layer": layer, "label": label, "available": False, "rho": None, "p_value": None, "n": 0})
                continue
            rho, p_value, n = self._spearman(merged[layer], merged["gene_effect"])
            relationships.append(
                {
                    "layer": layer,
                    "label": label,
                    "available": n >= 3,
                    "rho": rho,
                    "p_value": p_value,
                    "n": n,
                    "x_semantics": "RNA log2(TPM+1)" if layer == "expression" else "Relative copy number",
                    "y_semantics": "CRISPR Gene Effect; more negative = stronger dependency",
                }
            )

        points = _records(merged.sort_values("gene_effect", ascending=True, na_position="last"))
        for point in points:
            point["cancer_ids"] = memberships.get(str(point.get("model_id")), [])

        return _clean(
            {
                "gene_symbol": symbol,
                "available": True,
                "statistics_engine": "scipy" if spearmanr is not None else "rank-correlation fallback; p-value unavailable",
                "relationships": relationships,
                "points": points,
                "guardrails": {
                    "correlation": "Correlation is exploratory and does not establish causality.",
                    "expression": "High RNA expression is not equivalent to dependency or protein activity.",
                    "copy_number": "Relative copy number is shown continuously and is not automatically called amplification/deletion.",
                },
            }
        )

    @lru_cache(maxsize=1)
    def _mutations(self) -> pd.DataFrame:
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

    def mutation_associations(
        self,
        gene_symbol: str,
        *,
        min_mutated: int = 3,
        min_wildtype: int = 3,
        limit: int = 30,
    ) -> dict[str, Any]:
        symbol, _ = self._gene_or_raise(gene_symbol)
        ge = self._model_gene_layer("gene_effect", symbol)
        mutations = self._mutations()
        if ge.empty or mutations.empty or not {"model_id", "gene"}.issubset(mutations.columns):
            return {
                "gene_symbol": symbol,
                "available": False,
                "items": [],
                "note": "Mutation profile index is not available.",
            }

        mutations = mutations.copy()
        functional_mask = pd.Series(True, index=mutations.index)
        if "is_functional" in mutations.columns:
            functional_mask = mutations["is_functional"].map(_truthy)
        else:
            flags = [c for c in ("driver", "hotspot", "likely_lof", "oncogene_high_impact", "tumor_suppressor_high_impact") if c in mutations.columns]
            if flags:
                functional_mask = pd.Series(False, index=mutations.index)
                for col in flags:
                    functional_mask = functional_mask | mutations[col].map(_truthy)
        functional = mutations[functional_mask].copy()
        profiled_ids = set(mutations["model_id"].dropna().astype(str))
        ge = ge[ge["model_id"].astype(str).isin(profiled_ids)].dropna(subset=["gene_effect"]).copy()
        ge_lookup = ge.set_index("model_id")["gene_effect"]
        all_ids = set(ge_lookup.index.astype(str))
        if len(all_ids) < (min_mutated + min_wildtype):
            return {"gene_symbol": symbol, "available": True, "items": [], "profiled_models_n": len(all_ids)}

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
            rows.append(
                {
                    "mutation_gene": str(mutated_gene),
                    "mutated_n": int(len(x)),
                    "wildtype_n": int(len(y)),
                    "mutated_median_gene_effect": float(np.median(x)),
                    "wildtype_median_gene_effect": float(np.median(y)),
                    "delta_median_gene_effect": float(np.median(x) - np.median(y)),
                    "cliffs_delta": self._cliffs_delta(x, y),
                    "p_value": p_value,
                    "functional_variant_rows_n": int(len(group)),
                }
            )

        q_values = self._bh_qvalues([row["p_value"] for row in rows])
        for row, q_value in zip(rows, q_values):
            row["q_value"] = q_value
        rows.sort(
            key=lambda row: (
                row["q_value"] is None,
                row["q_value"] if row["q_value"] is not None else 1.0,
                -abs(row["delta_median_gene_effect"]),
            )
        )
        rows = rows[: max(1, min(int(limit), 100))]

        return _clean(
            {
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
                    "guardrail": "This is exploratory association screening across MCL cell models, not proof that the mutation causes dependency or predicts patient response.",
                },
            }
        )
