from __future__ import annotations

import json
from functools import lru_cache
from typing import Any

import pandas as pd

from .gene_explorer import _clean, _records
from .matrix_gene_explorer import MatrixGeneExplorerStore
from .store import MCLDataError


class RuntimeGeneExplorerStore(MatrixGeneExplorerStore):
    """Fast read-only Gene Explorer over materialized runtime indexes.

    Build-time code owns ontology projection and catalog construction. Interactive
    requests must never fall back to rebuilding those layers from scientific source
    tables. Missing runtime artifacts therefore fail fast with an actionable error.
    """

    def _runtime_path(self, name: str):
        path = self.index_dir / name
        if not path.exists():
            raise MCLDataError(
                f"Gene Explorer runtime index is missing: {name}. "
                "Run scripts/build-explorer.ps1 before starting MCL Explorer."
            )
        return path

    def _annotation_path(self):
        return self._runtime_path("gene_annotations.parquet")

    @lru_cache(maxsize=1)
    def catalog_frame(self) -> pd.DataFrame:
        return pd.read_parquet(self._runtime_path("gene_catalog.parquet"))

    @lru_cache(maxsize=1)
    def context_metrics_frame(self) -> pd.DataFrame:
        return pd.read_parquet(self._runtime_path("gene_context_metrics.parquet"))

    @lru_cache(maxsize=1)
    def _annotation_frame(self) -> pd.DataFrame:
        return pd.read_parquet(self._annotation_path())

    @lru_cache(maxsize=1)
    def formal_annotation_frame(self) -> pd.DataFrame:
        rows = self._annotation_frame()
        columns = [
            "gene_symbol",
            "source",
            "term_id",
            "term_name",
            "significant",
            "top_n",
            "source_version",
        ]
        if rows.empty or "annotation_type" not in rows.columns:
            return pd.DataFrame(columns=columns)

        formal = rows[rows["annotation_type"].astype(str) == "formal_term"].copy()
        if formal.empty:
            return pd.DataFrame(columns=columns)

        out = pd.DataFrame(index=formal.index)
        out["gene_symbol"] = formal["gene_symbol"].astype(str).str.upper()
        out["source"] = formal.get("source")
        out["term_id"] = formal.get("source_id")
        out["term_name"] = formal.get("source_term_name")
        out["significant"] = formal.get("source_significant", False)
        out["top_n"] = None
        out["source_version"] = formal.get("source_version")
        return out[columns].drop_duplicates(
            ["gene_symbol", "source", "term_id", "term_name"], keep="first"
        ).reset_index(drop=True)

    @staticmethod
    def _json_ids(value: Any) -> list[str]:
        if isinstance(value, list):
            return [str(x) for x in value if x is not None]
        try:
            parsed = json.loads(str(value or "[]"))
        except (json.JSONDecodeError, TypeError):
            return []
        return [str(x) for x in parsed if x is not None] if isinstance(parsed, list) else []

    def functional_coverage(self) -> dict[str, Any]:
        catalog = self.catalog_frame()
        total = int(catalog["gene_symbol"].nunique()) if not catalog.empty else 0
        annotated = 0
        if not catalog.empty and "mcl_domains_json" in catalog.columns:
            annotated = int(
                catalog["mcl_domains_json"].map(lambda value: bool(self._json_ids(value))).sum()
            )
        return {
            "annotated_genes_n": annotated,
            "gene_universe_n": total,
            "coverage_fraction": (annotated / total) if total else None,
            "status": "partial",
            "reason": "Coverage is read from the materialized MCL functional-domain index; absence of a category is not interpreted as absence of function.",
        }

    @lru_cache(maxsize=1)
    def _source_facets(self) -> list[dict[str, Any]]:
        frame = self._annotation_frame()
        if frame.empty or "source" not in frame.columns:
            return []
        if "annotation_type" in frame.columns:
            formal = frame[frame["annotation_type"].astype(str) == "formal_term"]
            if not formal.empty:
                frame = formal
        frame = frame.dropna(subset=["source", "gene_symbol"])
        if frame.empty:
            return []
        counts = frame.groupby("source")["gene_symbol"].nunique().sort_values(ascending=False)
        return [{"id": str(source), "genes_n": int(count)} for source, count in counts.items()]

    def facets(self) -> dict[str, Any]:
        catalog = self.catalog_frame()
        config = self._functional_config()

        def counts(column: str) -> dict[str, int]:
            if catalog.empty or column not in catalog.columns:
                return {}
            result: dict[str, int] = {}
            for value in catalog[column]:
                for item in set(self._json_ids(value)):
                    result[item] = result.get(item, 0) + 1
            return result

        domain_counts = counts("mcl_domains_json")
        subdomain_counts = counts("mcl_subdomains_json")
        protein_counts = counts("protein_classes_json")
        compartment_counts = counts("compartments_json")
        hallmark_counts = counts("hallmarks_json")

        domains: list[dict[str, Any]] = []
        for domain in config.get("domains") or []:
            item = {
                "id": domain.get("id"),
                "label_ru": domain.get("label_ru"),
                "label_en": domain.get("label_en"),
                "genes_n": int(domain_counts.get(str(domain.get("id")), 0)),
                "subdomains": [],
            }
            for subdomain in domain.get("subdomains") or []:
                item["subdomains"].append(
                    {
                        "id": subdomain.get("id"),
                        "label_ru": subdomain.get("label_ru"),
                        "genes_n": int(subdomain_counts.get(str(subdomain.get("id")), 0)),
                    }
                )
            domains.append(item)

        def secondary(config_key: str, values: dict[str, int]) -> list[dict[str, Any]]:
            return [
                {
                    "id": rule.get("id"),
                    "label_ru": rule.get("label_ru"),
                    "genes_n": int(values.get(str(rule.get("id")), 0)),
                }
                for rule in config.get(config_key) or []
            ]

        return _clean(
            {
                "taxonomy_version": config.get("version"),
                "status": "materialized_runtime",
                "domains": domains,
                "protein_classes": secondary("protein_classes", protein_counts),
                "compartments": secondary("compartments", compartment_counts),
                "hallmarks": secondary("hallmarks", hallmark_counts),
                "sources": self._source_facets(),
                "coverage": self.functional_coverage(),
                "reference_coverage": self.reference_coverage(),
                "provenance_note": "Interactive facets are read from precomputed Gene Explorer indexes. Ontology projection is performed during index build, not during page loading.",
            }
        )

    def _read_annotation_rows(
        self,
        gene_symbol: str,
        annotation_types: set[str] | None = None,
    ) -> pd.DataFrame:
        rows = self._annotation_frame()
        symbol = str(gene_symbol).strip().upper()
        if rows.empty or "gene_symbol" not in rows.columns:
            return rows.iloc[0:0].copy()
        frame = rows[rows["gene_symbol"].astype(str).str.upper() == symbol].copy()
        if annotation_types and "annotation_type" in frame.columns:
            frame = frame[frame["annotation_type"].astype(str).isin(annotation_types)]
        return frame.reset_index(drop=True)

    def annotations(self, gene_symbol: str) -> dict[str, Any]:
        symbol, _ = self._gene_or_raise(gene_symbol)
        rows = self._read_annotation_rows(symbol)

        def labels(annotation_type: str) -> list[str]:
            if rows.empty or "annotation_type" not in rows.columns:
                return []
            sub = rows[rows["annotation_type"].astype(str) == annotation_type]
            if sub.empty:
                return []
            column = "annotation_label_ru" if "annotation_label_ru" in sub.columns else "annotation_id"
            return sorted(set(sub[column].dropna().astype(str)))

        domains = labels("mcl_domain")
        subdomains = labels("mcl_subdomain")
        protein_classes = labels("protein_class")
        compartments = labels("compartment")
        hallmarks = labels("hallmark")
        formal = (
            rows[rows["annotation_type"].astype(str) == "formal_term"].copy()
            if not rows.empty and "annotation_type" in rows.columns
            else rows.iloc[0:0].copy()
        )
        formal_records: list[dict[str, Any]] = []
        for _, row in formal.iterrows():
            formal_records.append(
                {
                    "gene_symbol": symbol,
                    "source": row.get("source"),
                    "term_id": row.get("source_id"),
                    "term_name": row.get("source_term_name"),
                    "source_version": row.get("source_version"),
                }
            )

        return _clean(
            {
                "gene_symbol": symbol,
                "status": "materialized_runtime",
                "taxonomy_version": self._functional_config().get("version"),
                "mcl_domains": domains,
                "subdomains": subdomains,
                "protein_classes": protein_classes,
                "compartments": compartments,
                "hallmarks": hallmarks,
                "mcl_domain_details": _records(
                    rows[rows["annotation_type"].astype(str).isin(["mcl_domain", "mcl_subdomain"])]
                )
                if not rows.empty and "annotation_type" in rows.columns
                else [],
                "protein_class_details": _records(
                    rows[rows["annotation_type"].astype(str) == "protein_class"]
                )
                if not rows.empty and "annotation_type" in rows.columns
                else [],
                "compartment_details": _records(
                    rows[rows["annotation_type"].astype(str) == "compartment"]
                )
                if not rows.empty and "annotation_type" in rows.columns
                else [],
                "hallmark_details": _records(
                    rows[rows["annotation_type"].astype(str) == "hallmark"]
                )
                if not rows.empty and "annotation_type" in rows.columns
                else [],
                "formal_annotations": formal_records,
                "coverage": self.functional_coverage(),
                "reference_coverage": self.reference_coverage(),
                "reference_available": bool(self.reference_coverage().get("available")),
                "note": "Annotations are served only from the materialized provenance-aware Gene Explorer index.",
            }
        )
