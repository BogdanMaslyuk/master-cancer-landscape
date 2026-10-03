from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


class MCLDataError(RuntimeError):
    pass


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
        except (ValueError, TypeError):
            pass
    return value


def _records(frame: pd.DataFrame) -> list[dict[str, Any]]:
    return [_clean(row) for row in frame.to_dict("records")]


def _bool_mask(series: pd.Series, default: bool = False) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(default).astype(bool)
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .map({"true": True, "1": True, "yes": True, "false": False, "0": False, "no": False})
        .fillna(default)
        .astype(bool)
    )


class MCLDataStore:
    def __init__(self, root: Path):
        self.root = Path(root).resolve()

    def _path(self, relative: str) -> Path:
        return self.root / relative

    @lru_cache(maxsize=64)
    def _read_table(self, relative: str) -> pd.DataFrame:
        path = self._path(relative)
        if not path.exists():
            raise MCLDataError(f"Required MCL output is missing: {relative}")
        if path.suffix == ".parquet":
            return pd.read_parquet(path)
        return pd.read_csv(path, sep="\t", low_memory=False)

    def _read_preferred(self, stem: str) -> pd.DataFrame:
        parquet = self._path(stem + ".parquet")
        if parquet.exists():
            return self._read_table(stem + ".parquet").copy()
        return self._read_table(stem + ".tsv").copy()

    def _read_json(self, relative: str, default: Any = None) -> Any:
        path = self._path(relative)
        if not path.exists():
            return default
        return json.loads(path.read_text(encoding="utf-8-sig"))

    def _pathway_config(self) -> dict[str, Any]:
        path = self._path("config/pathways.yaml")
        if not path.exists():
            return {}
        return yaml.safe_load(path.read_text(encoding="utf-8")) or {}

    def comparison_specs(self) -> list[dict[str, Any]]:
        specs: list[dict[str, Any]] = []
        for item in self._pathway_config().get("comparisons") or []:
            cancer_id = str(item["cancer_id"])
            comparison = str(item["comparison"])
            comparison_id = f"{cancer_id}__{comparison}"
            meta = self._read_json(
                f"outputs/reports/depmap_genomewide_{comparison_id}_meta.json", {}
            ) or {}
            qc_rows = self._read_json(
                f"outputs/qc/depmap_genomewide_{comparison_id}_qc.json", []
            ) or []
            severities = [str(x.get("severity", "")) for x in qc_rows if isinstance(x, dict)]
            qc_status = "ERROR" if "ERROR" in severities else ("WARNING" if "WARNING" in severities else "PASS")
            specs.append(
                {
                    "id": comparison_id,
                    "label": str(item.get("label") or comparison_id),
                    "cancer_id": cancer_id,
                    "comparison": comparison,
                    "context_definition": meta.get("context_definition"),
                    "comparator_definition": meta.get("comparator_definition"),
                    "context_models_n": meta.get("context_models_n"),
                    "comparator_models_n": meta.get("comparator_models_n"),
                    "genes_analyzed_n": meta.get("genes_analyzed_n"),
                    "depmap_release": meta.get("depmap_release"),
                    "retrieved_at": meta.get("retrieved_at"),
                    "qc_status": qc_status,
                }
            )
        return _clean(specs)

    def _comparison_or_raise(self, comparison_id: str) -> dict[str, Any]:
        for spec in self.comparison_specs():
            if spec["id"] == comparison_id:
                return spec
        raise MCLDataError(f"Unknown comparison: {comparison_id}")

    def _gene_stability(self) -> pd.DataFrame:
        path = self._path("data/processed/pathways_sensitivity/gene_stability.tsv")
        if not path.exists():
            return pd.DataFrame()
        return self._read_table("data/processed/pathways_sensitivity/gene_stability.tsv").copy()

    def _recurrence(self) -> pd.DataFrame:
        path = self._path("data/processed/pathways/candidate_recurrence.tsv")
        if not path.exists():
            return pd.DataFrame()
        return self._read_table("data/processed/pathways/candidate_recurrence.tsv").copy()

    def summary(self) -> dict[str, Any]:
        comparisons = self.comparison_specs()
        genes_n = max(
            [int(x["genes_analyzed_n"]) for x in comparisons if x.get("genes_analyzed_n") is not None]
            or [0]
        )
        stable_genes_path = self._path("outputs/reports/M3_3_1_stable_recurrent_genes.tsv")
        stable_pathways_path = self._path("outputs/reports/M3_3_1_stable_pathways.tsv")
        stable_genes_n = len(self._read_table("outputs/reports/M3_3_1_stable_recurrent_genes.tsv")) if stable_genes_path.exists() else 0
        stable_pathways_n = len(self._read_table("outputs/reports/M3_3_1_stable_pathways.tsv")) if stable_pathways_path.exists() else 0
        sensitivity_meta = self._read_json("outputs/reports/M3_3_1_pathway_sensitivity_meta.json", {}) or {}
        per_threshold = sensitivity_meta.get("per_threshold") or {}
        return _clean(
            {
                "genes_analyzed_n": genes_n,
                "comparisons_n": len(comparisons),
                "stable_recurrent_genes_n": stable_genes_n,
                "stable_pathways_n": stable_pathways_n,
                "thresholds": sensitivity_meta.get("thresholds") or [50, 100, 200],
                "per_threshold": per_threshold,
                "analysis_version": sensitivity_meta.get("analysis_version"),
                "generated_at": sensitivity_meta.get("timestamp"),
                "data_release": next((x.get("depmap_release") for x in comparisons if x.get("depmap_release")), None),
            }
        )

    def comparison(self, comparison_id: str) -> dict[str, Any]:
        return self._comparison_or_raise(comparison_id)

    def comparison_genes(
        self,
        comparison_id: str,
        *,
        page: int = 1,
        page_size: int = 100,
        search: str | None = None,
        negative_delta_only: bool = False,
        exclude_broad: bool = False,
        exclude_low_sample: bool = False,
        fdr_only: bool = False,
        stable_only: bool = False,
        sort_by: str = "delta_gene_effect",
        sort_order: str = "asc",
    ) -> dict[str, Any]:
        self._comparison_or_raise(comparison_id)
        frame = self._read_preferred(f"data/processed/depmap_genomewide/{comparison_id}_genes")
        if search:
            frame = frame[frame["gene_symbol"].astype(str).str.contains(search, case=False, regex=False)]
        if negative_delta_only and "delta_gene_effect" in frame:
            frame = frame[pd.to_numeric(frame["delta_gene_effect"], errors="coerce") < 0]
        if exclude_broad and "broad_dependency_warning" in frame:
            frame = frame[~_bool_mask(frame["broad_dependency_warning"])]
        if exclude_low_sample and "low_sample_size" in frame:
            frame = frame[~_bool_mask(frame["low_sample_size"])]
        if fdr_only and "fdr_0_05" in frame:
            frame = frame[_bool_mask(frame["fdr_0_05"])]

        stability = self._gene_stability()
        if not stability.empty:
            stable_cols = [c for c in ["gene_symbol", "thresholds_n", "present_all_thresholds", "thresholds_present"] if c in stability.columns]
            frame = frame.merge(stability[stable_cols], on="gene_symbol", how="left")
        if stable_only:
            if "present_all_thresholds" not in frame:
                frame = frame.iloc[0:0]
            else:
                frame = frame[_bool_mask(frame["present_all_thresholds"])]

        allowed_sort = set(frame.columns)
        if sort_by not in allowed_sort:
            sort_by = "delta_gene_effect" if "delta_gene_effect" in frame else "gene_symbol"
        ascending = sort_order.lower() != "desc"
        frame = frame.sort_values(sort_by, ascending=ascending, na_position="last")
        total = int(len(frame))
        page = max(1, int(page))
        page_size = min(2000, max(1, int(page_size)))
        start = (page - 1) * page_size
        page_frame = frame.iloc[start : start + page_size].copy()
        return {
            "comparison_id": comparison_id,
            "page": page,
            "page_size": page_size,
            "total": total,
            "items": _records(page_frame),
        }

    def stable_genes(self) -> list[dict[str, Any]]:
        path = self._path("outputs/reports/M3_3_1_stable_recurrent_genes.tsv")
        if not path.exists():
            return []
        return _records(self._read_table("outputs/reports/M3_3_1_stable_recurrent_genes.tsv"))

    def genes(self, search: str | None = None, limit: int = 200) -> list[dict[str, Any]]:
        frame = self._gene_stability()
        if frame.empty:
            frame = self._recurrence()
        if frame.empty:
            return []
        if search:
            frame = frame[frame["gene_symbol"].astype(str).str.contains(search, case=False, regex=False)]
        return _records(frame.head(min(max(int(limit), 1), 2000)))

    def gene(self, gene_symbol: str) -> dict[str, Any]:
        symbol = gene_symbol.strip().upper()
        identity: dict[str, Any] = {"gene_symbol": symbol}
        identifiers = self._path("data/processed/target_identifiers.tsv")
        if identifiers.exists():
            ids = self._read_table("data/processed/target_identifiers.tsv")
            symbol_col = "hgnc_symbol" if "hgnc_symbol" in ids.columns else ("gene_symbol" if "gene_symbol" in ids.columns else None)
            if symbol_col:
                hit = ids[ids[symbol_col].astype(str).str.upper() == symbol]
                if not hit.empty:
                    identity.update(_clean(hit.iloc[0].to_dict()))

        stability_frame = self._gene_stability()
        stability = None
        if not stability_frame.empty:
            hit = stability_frame[stability_frame["gene_symbol"].astype(str).str.upper() == symbol]
            if not hit.empty:
                stability = _clean(hit.iloc[0].to_dict())

        comparison_rows: list[dict[str, Any]] = []
        for spec in self.comparison_specs():
            frame = self._read_preferred(f"data/processed/depmap_genomewide/{spec['id']}_genes")
            hit = frame[frame["gene_symbol"].astype(str).str.upper() == symbol]
            if hit.empty:
                continue
            row = _clean(hit.iloc[0].to_dict())
            row["comparison_id"] = spec["id"]
            row["comparison_label"] = spec["label"]
            comparison_rows.append(row)

        pathway_rows: list[dict[str, Any]] = []
        path = self._path("data/processed/pathways_sensitivity/enrichment_all.tsv")
        if path.exists():
            enrichment = self._read_table("data/processed/pathways_sensitivity/enrichment_all.tsv")
            gene_col = "intersecting_gene_symbols_json"
            if gene_col in enrichment.columns:
                def contains_gene(value: Any) -> bool:
                    try:
                        genes = json.loads(str(value))
                    except (json.JSONDecodeError, TypeError):
                        return False
                    return symbol in {str(x).upper() for x in genes if x is not None}
                pathway_rows = _records(enrichment[enrichment[gene_col].map(contains_gene)])

        return _clean(
            {
                "identity": identity,
                "stability": stability,
                "comparisons": comparison_rows,
                "pathways": pathway_rows,
            }
        )

    def pathways(
        self,
        *,
        top_n: int | None = None,
        source: str | None = None,
        stable_only: bool = False,
        significant_only: bool = True,
        search: str | None = None,
        limit: int = 1000,
    ) -> list[dict[str, Any]]:
        path = self._path("data/processed/pathways_sensitivity/enrichment_all.tsv")
        if not path.exists():
            return []
        frame = self._read_table("data/processed/pathways_sensitivity/enrichment_all.tsv").copy()
        if top_n is not None and "top_n" in frame:
            frame = frame[pd.to_numeric(frame["top_n"], errors="coerce") == int(top_n)]
        if source and "source" in frame:
            frame = frame[frame["source"].astype(str).str.upper() == source.upper()]
        if significant_only and "significant" in frame:
            frame = frame[_bool_mask(frame["significant"])]
        if search:
            text = frame.get("term_name", pd.Series("", index=frame.index)).astype(str)
            frame = frame[text.str.contains(search, case=False, regex=False)]
        if stable_only:
            stable_path = self._path("outputs/reports/M3_3_1_stable_pathways.tsv")
            if not stable_path.exists():
                return []
            stable = self._read_table("outputs/reports/M3_3_1_stable_pathways.tsv")
            keys = set(zip(stable["source"].astype(str), stable["term_id"].astype(str)))
            mask = [(str(s), str(t)) in keys for s, t in zip(frame["source"], frame["term_id"])]
            frame = frame[pd.Series(mask, index=frame.index)]
        sort_cols = [c for c in ["top_n", "p_value_adjusted"] if c in frame.columns]
        if sort_cols:
            frame = frame.sort_values(sort_cols, ascending=True, na_position="last")
        return _records(frame.head(min(max(int(limit), 1), 5000)))

    def pathway_stability(self) -> list[dict[str, Any]]:
        path = self._path("data/processed/pathways_sensitivity/term_stability.tsv")
        if not path.exists():
            return []
        return _records(self._read_table("data/processed/pathways_sensitivity/term_stability.tsv"))

    def network(self, stable_only: bool = True, limit_terms: int = 100) -> dict[str, Any]:
        terms = self.pathway_stability()
        if stable_only:
            terms = [x for x in terms if str(x.get("significant_all_thresholds", "")).lower() in {"true", "1"}]
        terms = terms[: max(1, min(int(limit_terms), 500))]
        nodes: dict[str, dict[str, Any]] = {}
        edges: list[dict[str, Any]] = []
        for term in terms:
            term_key = f"term:{term.get('source')}:{term.get('term_id')}"
            nodes[term_key] = {
                "id": term_key,
                "type": "pathway",
                "label": term.get("term_name") or term.get("term_id"),
                "source": term.get("source"),
                "stable": term.get("significant_all_thresholds"),
            }
            gene_union: set[str] = set()
            for key, value in term.items():
                if not str(key).startswith("intersecting_gene_symbols_top") or not value:
                    continue
                try:
                    gene_union.update(str(x) for x in json.loads(str(value)))
                except (json.JSONDecodeError, TypeError):
                    continue
            for gene in sorted(gene_union):
                gene_key = f"gene:{gene}"
                nodes.setdefault(gene_key, {"id": gene_key, "type": "gene", "label": gene})
                edges.append({"id": f"{gene_key}->{term_key}", "source": gene_key, "target": term_key})
        return {"nodes": list(nodes.values()), "edges": edges}

    def qc(self) -> dict[str, Any]:
        qc_dir = self._path("outputs/qc")
        records: list[dict[str, Any]] = []
        if qc_dir.exists():
            for path in sorted(qc_dir.glob("*_qc.json")):
                try:
                    payload = json.loads(path.read_text(encoding="utf-8-sig"))
                except json.JSONDecodeError:
                    continue
                rows = payload if isinstance(payload, list) else [payload]
                for row in rows:
                    if isinstance(row, dict):
                        records.append({"source_file": path.name, **_clean(row)})
        counts = {"ERROR": 0, "WARNING": 0, "INFO": 0}
        for row in records:
            sev = str(row.get("severity", "INFO")).upper()
            counts[sev] = counts.get(sev, 0) + 1
        return {"counts": counts, "records": records}
