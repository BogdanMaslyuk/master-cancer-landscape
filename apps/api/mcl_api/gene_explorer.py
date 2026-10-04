from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import yaml

from .store import MCLDataError, MCLDataStore


_CONTEXT_METRIC_ALIASES: dict[str, tuple[str, ...]] = {
    "gene_symbol": ("gene_symbol", "gene", "GeneSymbol", "HugoSymbol"),
    "context_median_gene_effect": (
        "context_median_gene_effect",
        "target_median_gene_effect",
        "context_median",
    ),
    "comparator_median_gene_effect": (
        "comparator_median_gene_effect",
        "control_median_gene_effect",
        "comparator_median",
    ),
    "delta_gene_effect": ("delta_gene_effect", "delta", "effect_delta"),
    "cliffs_delta": ("cliffs_delta", "cliff_delta", "cliffs_d"),
    "p_value": (
        "p_value",
        "mannwhitney_p_value",
        "mann_whitney_p_value",
        "mw_p_value",
    ),
    "q_value": (
        "q_value",
        "fdr_q_value",
        "mannwhitney_q_value",
        "mann_whitney_q_value",
        "p_value_adjusted",
    ),
    "fdr_0_05": ("fdr_0_05", "fdr05", "significant_fdr_0_05"),
    "context_models_n": ("context_models_n", "target_models_n", "context_n"),
    "comparator_models_n": ("comparator_models_n", "control_models_n", "comparator_n"),
    "broad_dependency_warning": (
        "broad_dependency_warning",
        "broad_dependency",
        "common_essential_warning",
    ),
    "low_sample_size": ("low_sample_size", "low_sample_warning"),
}


_MODEL_META_COLUMNS = (
    "model_id",
    "cell_line_name",
    "stripped_cell_line_name",
    "oncotree_lineage",
    "oncotree_primary_disease",
    "oncotree_subtype",
    "oncotree_code",
    "depmap_model_type",
    "primary_or_metastasis",
    "sample_collection_site",
)


def _clean(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean(v) for v in value]
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(value, "item"):
        try:
            return value.item()
        except (TypeError, ValueError):
            pass
    return value


def _records(frame: pd.DataFrame) -> list[dict[str, Any]]:
    return [_clean(row) for row in frame.to_dict("records")]


def _truthy(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    try:
        if pd.isna(value):
            return False
    except (TypeError, ValueError):
        pass
    return str(value).strip().lower() in {"true", "1", "yes", "y", "t"}


def _bool_series(series: pd.Series) -> pd.Series:
    return series.map(_truthy).astype(bool)


def _first_existing(columns: Iterable[str], aliases: Iterable[str]) -> str | None:
    available = set(str(x) for x in columns)
    return next((name for name in aliases if name in available), None)


def _number_series(frame: pd.DataFrame, name: str) -> pd.Series:
    if name not in frame.columns:
        return pd.Series(float("nan"), index=frame.index, dtype="float64")
    return pd.to_numeric(frame[name], errors="coerce")


class GeneExplorerStore:
    """Target -> cancer -> cell-model navigation over existing MCL/DepMap outputs.

    The store prefers materialized Parquet indexes when present, but can derive the
    same compact indexes from the existing genome-wide comparison outputs. Heavy
    derivations are cached for the lifetime of the API process.
    """

    def __init__(self, root: Path, store: MCLDataStore):
        self.root = Path(root).resolve()
        self.store = store
        self.processed = self.root / "data" / "processed"
        self.index_dir = self.processed / "gene_explorer"

    @lru_cache(maxsize=1)
    def _contexts_config(self) -> dict[str, Any]:
        path = self.root / "config" / "cancer_contexts.yaml"
        if not path.exists():
            return {}
        payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        return payload.get("contexts") or {}

    @lru_cache(maxsize=32)
    def _comparison_frame(self, comparison_id: str) -> pd.DataFrame:
        stem = self.processed / "depmap_genomewide" / f"{comparison_id}_genes"
        parquet = stem.with_suffix(".parquet")
        tsv = stem.with_suffix(".tsv")
        if parquet.exists():
            frame = pd.read_parquet(parquet)
        elif tsv.exists():
            frame = pd.read_csv(tsv, sep="\t", low_memory=False)
        else:
            raise MCLDataError(f"Genome-wide comparison output is missing: {comparison_id}")
        return frame

    def _normalize_comparison(self, comparison_id: str, spec: dict[str, Any]) -> pd.DataFrame:
        raw = self._comparison_frame(comparison_id)
        if raw.empty:
            return pd.DataFrame()
        out = pd.DataFrame(index=raw.index)
        for canonical, aliases in _CONTEXT_METRIC_ALIASES.items():
            source = _first_existing(raw.columns, aliases)
            if source is not None:
                out[canonical] = raw[source]
        if "gene_symbol" not in out.columns:
            return pd.DataFrame()
        out["gene_symbol"] = out["gene_symbol"].astype(str).str.strip().str.upper()
        out = out[out["gene_symbol"] != ""].copy()
        out["comparison_id"] = comparison_id
        out["comparison_label"] = spec.get("label") or comparison_id
        out["cancer_id"] = spec.get("cancer_id")
        cfg = self._contexts_config().get(str(spec.get("cancer_id"))) or {}
        out["cancer_name"] = cfg.get("name") or spec.get("cancer_id")
        out["depmap_release"] = spec.get("depmap_release")
        for col in (
            "context_median_gene_effect",
            "comparator_median_gene_effect",
            "delta_gene_effect",
            "cliffs_delta",
            "p_value",
            "q_value",
            "context_models_n",
            "comparator_models_n",
        ):
            if col in out.columns:
                out[col] = pd.to_numeric(out[col], errors="coerce")
        for col in ("fdr_0_05", "broad_dependency_warning", "low_sample_size"):
            if col in out.columns:
                out[col] = _bool_series(out[col])
        if "fdr_0_05" not in out.columns:
            q = _number_series(out, "q_value")
            out["fdr_0_05"] = q.lt(0.05) & q.notna()
        if "broad_dependency_warning" not in out.columns:
            out["broad_dependency_warning"] = False
        if "low_sample_size" not in out.columns:
            out["low_sample_size"] = False
        return out.reset_index(drop=True)

    @lru_cache(maxsize=1)
    def context_metrics_frame(self) -> pd.DataFrame:
        materialized = self.index_dir / "gene_context_metrics.parquet"
        if materialized.exists():
            return pd.read_parquet(materialized)
        frames: list[pd.DataFrame] = []
        for spec in self.store.comparison_specs():
            comparison_id = str(spec.get("id") or "")
            if not comparison_id:
                continue
            normalized = self._normalize_comparison(comparison_id, spec)
            if not normalized.empty:
                frames.append(normalized)
        if not frames:
            return pd.DataFrame(columns=["gene_symbol", "comparison_id", "cancer_id"])
        return pd.concat(frames, ignore_index=True)

    @lru_cache(maxsize=1)
    def _stability(self) -> pd.DataFrame:
        path = self.processed / "pathways_sensitivity" / "gene_stability.tsv"
        if not path.exists():
            return pd.DataFrame()
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        if "gene_symbol" in frame.columns:
            frame["gene_symbol"] = frame["gene_symbol"].astype(str).str.upper()
        return frame

    @lru_cache(maxsize=3)
    def _layer_gene_map(self, layer: str) -> list[str]:
        path = self.processed / f"depmap_model_{layer}_genes.json"
        if not path.exists():
            return []
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            return []
        genes: list[str] = []
        for row in payload if isinstance(payload, list) else []:
            if not isinstance(row, dict):
                continue
            symbol = str(row.get("gene") or "").strip().upper()
            if symbol and symbol not in genes:
                genes.append(symbol)
        return genes

    @lru_cache(maxsize=1)
    def _target_identifiers(self) -> pd.DataFrame:
        path = self.processed / "target_identifiers.tsv"
        if not path.exists():
            return pd.DataFrame()
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        symbol_col = _first_existing(frame.columns, ("hgnc_symbol", "gene_symbol", "input_symbol"))
        if symbol_col is not None:
            frame = frame.copy()
            frame["_symbol"] = frame[symbol_col].astype(str).str.upper()
        return frame

    @lru_cache(maxsize=1)
    def _best_model_gene_effect(self) -> dict[str, float]:
        path = self.processed / "depmap_model_gene_effect.parquet"
        if not path.exists():
            return {}
        try:
            frame = pd.read_parquet(path)
        except Exception:
            return {}
        if "model_id" in frame.columns:
            frame = frame.drop(columns=["model_id"])
        if frame.empty:
            return {}
        minima = frame.apply(pd.to_numeric, errors="coerce").min(axis=0, skipna=True)
        return {
            str(gene).upper(): float(value)
            for gene, value in minima.items()
            if pd.notna(value)
        }

    def _build_catalog(self) -> pd.DataFrame:
        metrics = self.context_metrics_frame()
        universe: set[str] = set()
        if not metrics.empty:
            universe.update(metrics["gene_symbol"].dropna().astype(str).str.upper())
        for layer in ("gene_effect", "expression", "copy_number"):
            universe.update(self._layer_gene_map(layer))
        stability = self._stability()
        if not stability.empty and "gene_symbol" in stability.columns:
            universe.update(stability["gene_symbol"].dropna().astype(str).str.upper())
        ids = self._target_identifiers()
        if not ids.empty and "_symbol" in ids.columns:
            universe.update(ids["_symbol"].dropna().astype(str).str.upper())
        if not universe:
            return pd.DataFrame(columns=["gene_symbol"])

        catalog = pd.DataFrame({"gene_symbol": sorted(universe)})
        if not metrics.empty:
            rows: list[dict[str, Any]] = []
            for gene, group in metrics.groupby("gene_symbol", sort=False):
                delta = _number_series(group, "delta_gene_effect")
                q = _number_series(group, "q_value")
                cliffs = _number_series(group, "cliffs_delta")
                best_delta_idx = delta.idxmin() if delta.notna().any() else None
                abs_cliffs = cliffs.abs()
                best_cliff_idx = abs_cliffs.idxmax() if abs_cliffs.notna().any() else None
                significant = _bool_series(group["fdr_0_05"]) if "fdr_0_05" in group.columns else q.lt(0.05)
                rows.append(
                    {
                        "gene_symbol": str(gene),
                        "comparisons_n": int(group["comparison_id"].nunique()),
                        "significant_comparisons_n": int(significant.sum()),
                        "best_delta_gene_effect": float(delta.loc[best_delta_idx]) if best_delta_idx is not None else None,
                        "best_context_id": group.loc[best_delta_idx, "comparison_id"] if best_delta_idx is not None else None,
                        "best_context_label": group.loc[best_delta_idx, "comparison_label"] if best_delta_idx is not None else None,
                        "best_q_value": float(q.min()) if q.notna().any() else None,
                        "best_cliffs_delta": float(cliffs.loc[best_cliff_idx]) if best_cliff_idx is not None else None,
                        "broad_dependency_any": bool(_bool_series(group["broad_dependency_warning"]).any()) if "broad_dependency_warning" in group.columns else False,
                        "low_sample_any": bool(_bool_series(group["low_sample_size"]).any()) if "low_sample_size" in group.columns else False,
                    }
                )
            catalog = catalog.merge(pd.DataFrame(rows), on="gene_symbol", how="left")

        if not stability.empty and "gene_symbol" in stability.columns:
            stable_cols = [
                c
                for c in (
                    "gene_symbol",
                    "recurrent_top50",
                    "recurrent_top100",
                    "recurrent_top200",
                    "thresholds_n",
                    "present_all_thresholds",
                    "thresholds_present",
                )
                if c in stability.columns
            ]
            catalog = catalog.merge(stability[stable_cols].drop_duplicates("gene_symbol"), on="gene_symbol", how="left")

        best_model = self._best_model_gene_effect()
        catalog["best_model_gene_effect"] = catalog["gene_symbol"].map(best_model)
        for col in ("comparisons_n", "significant_comparisons_n", "thresholds_n"):
            if col in catalog.columns:
                catalog[col] = pd.to_numeric(catalog[col], errors="coerce").fillna(0).astype(int)
        for col in ("broad_dependency_any", "low_sample_any", "present_all_thresholds"):
            if col not in catalog.columns:
                catalog[col] = False
            else:
                catalog[col] = _bool_series(catalog[col])
        return catalog.sort_values("gene_symbol").reset_index(drop=True)

    @lru_cache(maxsize=1)
    def catalog_frame(self) -> pd.DataFrame:
        materialized = self.index_dir / "gene_catalog.parquet"
        if materialized.exists():
            return pd.read_parquet(materialized)
        return self._build_catalog()

    def materialize_indexes(self) -> dict[str, Any]:
        self.index_dir.mkdir(parents=True, exist_ok=True)
        # Build directly from source outputs rather than re-reading an older materialized copy.
        self.context_metrics_frame.cache_clear()
        metrics = self.context_metrics_frame()
        metrics.to_parquet(self.index_dir / "gene_context_metrics.parquet", index=False, compression="zstd")
        self.catalog_frame.cache_clear()
        catalog = self._build_catalog()
        catalog.to_parquet(self.index_dir / "gene_catalog.parquet", index=False, compression="zstd")
        manifest = {
            "genes_n": int(catalog["gene_symbol"].nunique()) if not catalog.empty else 0,
            "context_rows_n": int(len(metrics)),
            "comparisons_n": int(metrics["comparison_id"].nunique()) if not metrics.empty else 0,
            "files": ["gene_catalog.parquet", "gene_context_metrics.parquet"],
        }
        (self.index_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
        self.catalog_frame.cache_clear()
        self.context_metrics_frame.cache_clear()
        return manifest

    def _gene_or_raise(self, gene_symbol: str) -> tuple[str, pd.Series]:
        symbol = str(gene_symbol).strip().upper()
        if not symbol:
            raise MCLDataError("Gene symbol is empty")
        catalog = self.catalog_frame()
        hit = catalog[catalog["gene_symbol"].astype(str).str.upper() == symbol]
        if hit.empty:
            raise MCLDataError(f"Gene is not present in the current MCL/DepMap gene universe: {symbol}")
        return symbol, hit.iloc[0]

    def suggest(self, query: str, limit: int = 12) -> list[dict[str, Any]]:
        q = str(query or "").strip().upper()
        if not q:
            return []
        frame = self.catalog_frame().copy()
        symbols = frame["gene_symbol"].astype(str).str.upper()
        frame = frame[symbols.str.contains(q, regex=False)].copy()
        if frame.empty:
            return []
        frame["_rank"] = frame["gene_symbol"].astype(str).str.upper().map(
            lambda value: 0 if value == q else (1 if value.startswith(q) else 2)
        )
        frame = frame.sort_values(["_rank", "gene_symbol"], ascending=[True, True]).head(max(1, min(int(limit), 30)))
        fields = [
            c
            for c in (
                "gene_symbol",
                "best_context_label",
                "best_delta_gene_effect",
                "best_model_gene_effect",
                "present_all_thresholds",
            )
            if c in frame.columns
        ]
        return _records(frame[fields])

    def search(
        self,
        *,
        query: str | None = None,
        cancer_id: str | None = None,
        comparison_id: str | None = None,
        gene_effect_max: float | None = None,
        delta_gene_effect_max: float | None = None,
        q_value_max: float | None = None,
        cliffs_delta_abs_min: float | None = None,
        stable_only: bool = False,
        exclude_broad: bool = False,
        exclude_low_sample: bool = False,
        page: int = 1,
        page_size: int = 100,
        sort_by: str = "best_delta_gene_effect",
        sort_order: str = "asc",
    ) -> dict[str, Any]:
        catalog = self.catalog_frame().copy()
        if query:
            q = str(query).strip()
            catalog = catalog[catalog["gene_symbol"].astype(str).str.contains(q, case=False, regex=False)]

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

        allowed_sort = {
            "gene_symbol",
            "best_delta_gene_effect",
            "best_model_gene_effect",
            "best_q_value",
            "best_cliffs_delta",
            "comparisons_n",
            "significant_comparisons_n",
        }
        if sort_by not in allowed_sort or sort_by not in catalog.columns:
            sort_by = "best_delta_gene_effect" if "best_delta_gene_effect" in catalog.columns else "gene_symbol"
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
                "items": _records(page_frame),
            }
        )

    def identity(self, gene_symbol: str) -> dict[str, Any]:
        symbol, row = self._gene_or_raise(gene_symbol)
        identity: dict[str, Any] = {
            "gene_symbol": symbol,
            "gene_name": None,
            "hgnc_id": None,
            "ensembl_gene_id": None,
            "uniprot_id": None,
        }
        ids = self._target_identifiers()
        if not ids.empty and "_symbol" in ids.columns:
            hit = ids[ids["_symbol"] == symbol]
            if not hit.empty:
                source = hit.iloc[0].to_dict()
                identity.update(
                    {
                        "hgnc_id": source.get("hgnc_id") or source.get("HGNC_ID"),
                        "ensembl_gene_id": source.get("ensembl_gene_id") or source.get("existing_ensembl_gene_id"),
                        "uniprot_id": source.get("uniprot_id") or source.get("existing_uniprot_id"),
                    }
                )
        stability = {
            key: _clean(row.get(key))
            for key in ("recurrent_top50", "recurrent_top100", "recurrent_top200", "thresholds_n", "present_all_thresholds", "thresholds_present")
            if key in row.index
        }
        summary = {
            key: _clean(row.get(key))
            for key in (
                "comparisons_n",
                "significant_comparisons_n",
                "best_delta_gene_effect",
                "best_context_id",
                "best_context_label",
                "best_q_value",
                "best_cliffs_delta",
                "best_model_gene_effect",
                "broad_dependency_any",
                "low_sample_any",
            )
            if key in row.index
        }
        return {"identity": identity, "summary": summary, "stability": stability}

    def contexts(self, gene_symbol: str) -> dict[str, Any]:
        symbol, _ = self._gene_or_raise(gene_symbol)
        frame = self.context_metrics_frame()
        if frame.empty:
            items: list[dict[str, Any]] = []
        else:
            sub = frame[frame["gene_symbol"].astype(str).str.upper() == symbol].copy()
            if "delta_gene_effect" in sub.columns:
                sub = sub.sort_values("delta_gene_effect", ascending=True, na_position="last")
            items = _records(sub)
        return {"gene_symbol": symbol, "total": len(items), "items": items}

    @lru_cache(maxsize=1)
    def _model_metadata(self) -> pd.DataFrame:
        parquet = self.processed / "depmap_model_metadata.parquet"
        tsv = self.processed / "depmap_model_metadata.tsv"
        if parquet.exists():
            frame = pd.read_parquet(parquet)
        elif tsv.exists():
            frame = pd.read_csv(tsv, sep="\t", low_memory=False)
        else:
            return pd.DataFrame()
        if "model_id" in frame.columns:
            frame["model_id"] = frame["model_id"].astype(str)
        return frame

    @lru_cache(maxsize=1)
    def _audit(self) -> pd.DataFrame:
        parquet = self.processed / "depmap_context_audit.parquet"
        tsv = self.processed / "depmap_context_audit.tsv"
        if parquet.exists():
            frame = pd.read_parquet(parquet)
        elif tsv.exists():
            frame = pd.read_csv(tsv, sep="\t", low_memory=False)
        else:
            return pd.DataFrame()
        if "model_id" in frame.columns:
            frame["model_id"] = frame["model_id"].astype(str)
        return frame

    @lru_cache(maxsize=1536)
    def _model_gene_layer(self, layer: str, gene_symbol: str) -> pd.DataFrame:
        symbol = gene_symbol.upper()
        path = self.processed / f"depmap_model_{layer}.parquet"
        if not path.exists():
            return pd.DataFrame(columns=["model_id", layer])
        gene_map = set(self._layer_gene_map(layer))
        if gene_map and symbol not in gene_map:
            return pd.DataFrame(columns=["model_id", layer])
        try:
            frame = pd.read_parquet(path, columns=["model_id", symbol])
        except (KeyError, ValueError):
            return pd.DataFrame(columns=["model_id", layer])
        if symbol not in frame.columns:
            return pd.DataFrame(columns=["model_id", layer])
        frame = frame.rename(columns={symbol: layer})
        frame["model_id"] = frame["model_id"].astype(str)
        frame[layer] = pd.to_numeric(frame[layer], errors="coerce")
        return frame

    def models(
        self,
        gene_symbol: str,
        *,
        cancer_id: str | None = None,
        gene_effect_max: float | None = None,
        page: int = 1,
        page_size: int = 100,
        sort_by: str = "gene_effect",
        sort_order: str = "asc",
    ) -> dict[str, Any]:
        symbol, _ = self._gene_or_raise(gene_symbol)
        layers = {
            "gene_effect": self._model_gene_layer("gene_effect", symbol),
            "expression": self._model_gene_layer("expression", symbol),
            "copy_number": self._model_gene_layer("copy_number", symbol),
        }
        available_layers = {key: not value.empty for key, value in layers.items()}
        frames = [frame for frame in layers.values() if not frame.empty]
        if not frames:
            return {
                "gene_symbol": symbol,
                "available": False,
                "available_layers": available_layers,
                "total": 0,
                "page": 1,
                "page_size": page_size,
                "items": [],
                "note": "Model-level multi-omics indexes are not built for this gene yet.",
            }
        merged = frames[0].copy()
        for frame in frames[1:]:
            merged = merged.merge(frame, on="model_id", how="outer")

        metadata = self._model_metadata()
        if not metadata.empty and "model_id" in metadata.columns:
            cols = [c for c in _MODEL_META_COLUMNS if c in metadata.columns]
            merged = merged.merge(metadata[cols].drop_duplicates("model_id"), on="model_id", how="left")

        audit = self._audit()
        memberships: dict[str, list[dict[str, Any]]] = {}
        if not audit.empty and "model_id" in audit.columns:
            audit_sub = audit.copy()
            if cancer_id:
                if "cancer_id" not in audit_sub.columns:
                    audit_sub = audit_sub.iloc[0:0]
                else:
                    audit_sub = audit_sub[audit_sub["cancer_id"].astype(str) == str(cancer_id)]
                allowed = set(audit_sub["model_id"].astype(str))
                merged = merged[merged["model_id"].astype(str).isin(allowed)]
            membership_cols = [c for c in ("cancer_id", "assigned_group") if c in audit_sub.columns]
            if membership_cols:
                for model_id, group in audit_sub.groupby("model_id"):
                    memberships[str(model_id)] = _records(group[["model_id"] + membership_cols].drop_duplicates())

        if gene_effect_max is not None:
            values = _number_series(merged, "gene_effect")
            merged = merged[values.notna() & (values <= float(gene_effect_max))]

        allowed_sort = {"gene_effect", "expression", "copy_number", "model_id", "cell_line_name"}
        if sort_by not in allowed_sort or sort_by not in merged.columns:
            sort_by = "gene_effect" if "gene_effect" in merged.columns else "model_id"
        ascending = str(sort_order).lower() != "desc"
        merged = merged.sort_values(sort_by, ascending=ascending, na_position="last")
        total = int(len(merged))
        page = max(1, int(page))
        page_size = max(1, min(int(page_size), 500))
        start = (page - 1) * page_size
        page_frame = merged.iloc[start : start + page_size].copy()
        items = _records(page_frame)
        for item in items:
            item["memberships"] = memberships.get(str(item.get("model_id")), [])
        return _clean(
            {
                "gene_symbol": symbol,
                "available": True,
                "available_layers": available_layers,
                "total": total,
                "page": page,
                "page_size": page_size,
                "pages": max(1, (total + page_size - 1) // page_size),
                "sort_by": sort_by,
                "sort_order": "asc" if ascending else "desc",
                "items": items,
                "guardrails": {
                    "gene_effect": "Более отрицательный Gene Effect означает более сильную зависимость после CRISPR-выключения; это не эквивалент фармакологического ингибирования.",
                    "expression": "RNA expression не доказывает активность белка.",
                    "copy_number": "Copy number показан как относительное значение DepMap и не переводится автоматически в amplification/deletion call.",
                },
            }
        )

    def annotations(self, gene_symbol: str) -> dict[str, Any]:
        symbol, _ = self._gene_or_raise(gene_symbol)
        existing = self.store.gene(symbol)
        pathways = existing.get("pathways") or []
        significant = [row for row in pathways if _truthy(row.get("significant"))]
        unique: dict[tuple[str, str], dict[str, Any]] = {}
        for row in significant:
            key = (str(row.get("source") or ""), str(row.get("term_id") or ""))
            if key not in unique:
                unique[key] = row
        return {
            "gene_symbol": symbol,
            "status": "partial",
            "mcl_domains": [],
            "subdomains": [],
            "protein_classes": [],
            "compartments": [],
            "hallmarks": [],
            "formal_annotations": list(unique.values()),
            "note": "v1 exposes current statistically supported MCL pathway/complex links. Full multi-label MCL Functional Domains and complete GO/Reactome/KEGG/CORUM mappings are the next annotation layer and will retain source/version provenance.",
        }
