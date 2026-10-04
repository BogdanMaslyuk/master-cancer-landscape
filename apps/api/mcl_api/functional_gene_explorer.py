from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from .gene_explorer import (
    GeneExplorerStore,
    _bool_series,
    _clean,
    _number_series,
    _records,
    _truthy,
)


class FunctionalGeneExplorerStore(GeneExplorerStore):
    """Gene Explorer with provenance-aware MCL Functional Domains.

    The taxonomy is deliberately a projection layer over formal annotations already
    present in MCL. It never assigns a domain merely from a gene symbol. Every domain
    row keeps the formal source term that triggered the rule, so users can distinguish
    source evidence from the human-readable MCL grouping.
    """

    @lru_cache(maxsize=1)
    def _functional_config(self) -> dict[str, Any]:
        path = self.root / "config" / "mcl_functional_domains.yaml"
        if not path.exists():
            return {"version": None, "domains": []}
        return yaml.safe_load(path.read_text(encoding="utf-8")) or {"version": None, "domains": []}

    @staticmethod
    def _contains_keyword(text: str, keywords: list[str]) -> bool:
        value = str(text or "").casefold()
        return any(str(keyword).casefold() in value for keyword in keywords if str(keyword).strip())

    @lru_cache(maxsize=1)
    def formal_annotation_frame(self) -> pd.DataFrame:
        path = self.processed / "pathways_sensitivity" / "enrichment_all.tsv"
        columns = [
            "gene_symbol",
            "source",
            "term_id",
            "term_name",
            "significant",
            "top_n",
            "source_version",
        ]
        if not path.exists():
            return pd.DataFrame(columns=columns)
        raw = pd.read_csv(path, sep="\t", low_memory=False)
        gene_col = "intersecting_gene_symbols_json"
        if raw.empty or gene_col not in raw.columns:
            return pd.DataFrame(columns=columns)

        rows: list[dict[str, Any]] = []
        for _, row in raw.iterrows():
            try:
                genes = json.loads(str(row.get(gene_col) or "[]"))
            except (json.JSONDecodeError, TypeError):
                continue
            if not isinstance(genes, list):
                continue
            source = str(row.get("source") or "").strip()
            term_id = str(row.get("term_id") or "").strip()
            term_name = str(row.get("term_name") or "").strip()
            if not term_name and not term_id:
                continue
            for gene in genes:
                symbol = str(gene or "").strip().upper()
                if not symbol:
                    continue
                rows.append(
                    {
                        "gene_symbol": symbol,
                        "source": source or None,
                        "term_id": term_id or None,
                        "term_name": term_name or term_id or None,
                        "significant": _truthy(row.get("significant")),
                        "top_n": row.get("top_n"),
                        "source_version": row.get("source_version") if "source_version" in raw.columns else None,
                    }
                )
        if not rows:
            return pd.DataFrame(columns=columns)
        frame = pd.DataFrame(rows).drop_duplicates(
            ["gene_symbol", "source", "term_id", "term_name"], keep="first"
        )
        return frame.reset_index(drop=True)

    @lru_cache(maxsize=1)
    def domain_mapping_frame(self) -> pd.DataFrame:
        formal = self.formal_annotation_frame()
        columns = [
            "gene_symbol",
            "annotation_type",
            "annotation_id",
            "annotation_label_ru",
            "domain_id",
            "domain_label_ru",
            "source",
            "source_id",
            "source_term_name",
            "source_version",
            "source_significant",
            "mapping_method",
            "taxonomy_version",
        ]
        if formal.empty:
            return pd.DataFrame(columns=columns)
        config = self._functional_config()
        domains = config.get("domains") or []
        taxonomy_version = config.get("version")
        rows: list[dict[str, Any]] = []
        seen: set[tuple[str, str, str, str | None]] = set()

        for _, term in formal.iterrows():
            term_name = str(term.get("term_name") or "")
            for domain in domains:
                domain_id = str(domain.get("id") or "").strip()
                domain_label = str(domain.get("label_ru") or domain_id).strip()
                if not domain_id or not self._contains_keyword(term_name, list(domain.get("keywords") or [])):
                    continue
                key = (str(term["gene_symbol"]), "mcl_domain", domain_id, str(term.get("term_id") or ""))
                if key not in seen:
                    seen.add(key)
                    rows.append(
                        {
                            "gene_symbol": term["gene_symbol"],
                            "annotation_type": "mcl_domain",
                            "annotation_id": domain_id,
                            "annotation_label_ru": domain_label,
                            "domain_id": domain_id,
                            "domain_label_ru": domain_label,
                            "source": term.get("source"),
                            "source_id": term.get("term_id"),
                            "source_term_name": term.get("term_name"),
                            "source_version": term.get("source_version"),
                            "source_significant": bool(term.get("significant")),
                            "mapping_method": "rule_based_from_formal_annotation",
                            "taxonomy_version": taxonomy_version,
                        }
                    )

                for subdomain in domain.get("subdomains") or []:
                    sub_id = str(subdomain.get("id") or "").strip()
                    sub_label = str(subdomain.get("label_ru") or sub_id).strip()
                    if not sub_id or not self._contains_keyword(term_name, list(subdomain.get("keywords") or [])):
                        continue
                    sub_key = (str(term["gene_symbol"]), "mcl_subdomain", sub_id, str(term.get("term_id") or ""))
                    if sub_key in seen:
                        continue
                    seen.add(sub_key)
                    rows.append(
                        {
                            "gene_symbol": term["gene_symbol"],
                            "annotation_type": "mcl_subdomain",
                            "annotation_id": sub_id,
                            "annotation_label_ru": sub_label,
                            "domain_id": domain_id,
                            "domain_label_ru": domain_label,
                            "source": term.get("source"),
                            "source_id": term.get("term_id"),
                            "source_term_name": term.get("term_name"),
                            "source_version": term.get("source_version"),
                            "source_significant": bool(term.get("significant")),
                            "mapping_method": "rule_based_from_formal_annotation",
                            "taxonomy_version": taxonomy_version,
                        }
                    )
        if not rows:
            return pd.DataFrame(columns=columns)
        return pd.DataFrame(rows, columns=columns).reset_index(drop=True)

    def _catalog_with_taxonomy(self, catalog: pd.DataFrame) -> pd.DataFrame:
        mapping = self.domain_mapping_frame()
        if catalog.empty:
            return catalog
        out = catalog.copy()
        if mapping.empty:
            out["mcl_domains_json"] = "[]"
            out["mcl_subdomains_json"] = "[]"
            out["functional_annotations_n"] = 0
            return out

        domain_rows = mapping[mapping["annotation_type"] == "mcl_domain"]
        sub_rows = mapping[mapping["annotation_type"] == "mcl_subdomain"]

        def group_json(frame: pd.DataFrame) -> dict[str, str]:
            if frame.empty:
                return {}
            return {
                str(gene): json.dumps(sorted(set(group["annotation_id"].astype(str))), ensure_ascii=False)
                for gene, group in frame.groupby("gene_symbol")
            }

        domains = group_json(domain_rows)
        subdomains = group_json(sub_rows)
        counts = mapping.groupby("gene_symbol").size().to_dict()
        symbols = out["gene_symbol"].astype(str).str.upper()
        out["mcl_domains_json"] = symbols.map(domains).fillna("[]")
        out["mcl_subdomains_json"] = symbols.map(subdomains).fillna("[]")
        out["functional_annotations_n"] = symbols.map(counts).fillna(0).astype(int)
        return out

    @lru_cache(maxsize=1)
    def catalog_frame(self) -> pd.DataFrame:
        materialized = self.index_dir / "gene_catalog.parquet"
        if materialized.exists():
            frame = pd.read_parquet(materialized)
            if "mcl_domains_json" in frame.columns:
                return frame
        return self._catalog_with_taxonomy(self._build_catalog())

    def materialize_indexes(self) -> dict[str, Any]:
        self.index_dir.mkdir(parents=True, exist_ok=True)
        self.context_metrics_frame.cache_clear()
        metrics = self.context_metrics_frame()
        metrics.to_parquet(self.index_dir / "gene_context_metrics.parquet", index=False, compression="zstd")

        self.formal_annotation_frame.cache_clear()
        self.domain_mapping_frame.cache_clear()
        formal = self.formal_annotation_frame()
        mapping = self.domain_mapping_frame()
        annotation_rows: list[pd.DataFrame] = []
        if not formal.empty:
            formal_out = formal.rename(
                columns={"term_id": "source_id", "term_name": "source_term_name", "significant": "source_significant"}
            ).copy()
            formal_out["annotation_type"] = "formal_term"
            formal_out["annotation_id"] = formal_out["source"].astype(str) + ":" + formal_out["source_id"].astype(str)
            formal_out["annotation_label_ru"] = formal_out["source_term_name"]
            formal_out["domain_id"] = None
            formal_out["domain_label_ru"] = None
            formal_out["mapping_method"] = "source_term_intersection_from_mcl_enrichment"
            formal_out["taxonomy_version"] = self._functional_config().get("version")
            keep = [
                "gene_symbol", "annotation_type", "annotation_id", "annotation_label_ru",
                "domain_id", "domain_label_ru", "source", "source_id", "source_term_name",
                "source_version", "source_significant", "mapping_method", "taxonomy_version",
            ]
            annotation_rows.append(formal_out[keep])
        if not mapping.empty:
            annotation_rows.append(mapping)
        annotations = pd.concat(annotation_rows, ignore_index=True) if annotation_rows else pd.DataFrame()
        annotations.to_parquet(self.index_dir / "gene_annotations.parquet", index=False, compression="zstd")

        catalog = self._catalog_with_taxonomy(self._build_catalog())
        catalog.to_parquet(self.index_dir / "gene_catalog.parquet", index=False, compression="zstd")
        manifest = {
            "genes_n": int(catalog["gene_symbol"].nunique()) if not catalog.empty else 0,
            "context_rows_n": int(len(metrics)),
            "comparisons_n": int(metrics["comparison_id"].nunique()) if not metrics.empty else 0,
            "formal_annotation_rows_n": int(len(formal)),
            "functional_mapping_rows_n": int(len(mapping)),
            "functionally_annotated_genes_n": int(mapping["gene_symbol"].nunique()) if not mapping.empty else 0,
            "functional_taxonomy_version": self._functional_config().get("version"),
            "files": ["gene_catalog.parquet", "gene_context_metrics.parquet", "gene_annotations.parquet"],
        }
        (self.index_dir / "manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        self.catalog_frame.cache_clear()
        self.context_metrics_frame.cache_clear()
        return manifest

    @staticmethod
    def _json_contains(value: Any, needle: str) -> bool:
        try:
            items = json.loads(str(value or "[]"))
        except (json.JSONDecodeError, TypeError):
            return False
        return str(needle) in {str(item) for item in items if item is not None}

    def search(
        self,
        *,
        query: str | None = None,
        domain: str | None = None,
        subdomain: str | None = None,
        pathway: str | None = None,
        annotation_source: str | None = None,
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
        if domain and "mcl_domains_json" in catalog.columns:
            catalog = catalog[catalog["mcl_domains_json"].map(lambda value: self._json_contains(value, domain))]
        if subdomain and "mcl_subdomains_json" in catalog.columns:
            catalog = catalog[catalog["mcl_subdomains_json"].map(lambda value: self._json_contains(value, subdomain))]
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

        allowed_sort = {
            "gene_symbol", "best_delta_gene_effect", "best_model_gene_effect", "best_q_value",
            "best_cliffs_delta", "comparisons_n", "significant_comparisons_n", "functional_annotations_n",
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
                "functional_coverage": self.functional_coverage(),
                "items": _records(page_frame),
            }
        )

    def functional_coverage(self) -> dict[str, Any]:
        catalog = self.catalog_frame()
        mapping = self.domain_mapping_frame()
        total = int(catalog["gene_symbol"].nunique()) if not catalog.empty else 0
        annotated = int(mapping["gene_symbol"].nunique()) if not mapping.empty else 0
        return {
            "annotated_genes_n": annotated,
            "gene_universe_n": total,
            "coverage_fraction": (annotated / total) if total else None,
            "status": "partial",
            "reason": "Domains are projected only from formal terms already observed in the current MCL enrichment layer.",
        }

    def facets(self) -> dict[str, Any]:
        config = self._functional_config()
        mapping = self.domain_mapping_frame()
        domain_counts = {}
        subdomain_counts = {}
        if not mapping.empty:
            domains = mapping[mapping["annotation_type"] == "mcl_domain"]
            subdomains = mapping[mapping["annotation_type"] == "mcl_subdomain"]
            domain_counts = domains.groupby("annotation_id")["gene_symbol"].nunique().to_dict()
            subdomain_counts = subdomains.groupby("annotation_id")["gene_symbol"].nunique().to_dict()
        output_domains: list[dict[str, Any]] = []
        for domain in config.get("domains") or []:
            item = {
                "id": domain.get("id"),
                "label_ru": domain.get("label_ru"),
                "label_en": domain.get("label_en"),
                "genes_n": int(domain_counts.get(domain.get("id"), 0)),
                "subdomains": [],
            }
            for subdomain in domain.get("subdomains") or []:
                item["subdomains"].append(
                    {
                        "id": subdomain.get("id"),
                        "label_ru": subdomain.get("label_ru"),
                        "genes_n": int(subdomain_counts.get(subdomain.get("id"), 0)),
                    }
                )
            output_domains.append(item)
        formal = self.formal_annotation_frame()
        sources = []
        if not formal.empty:
            for source, group in formal.groupby("source"):
                sources.append({"id": str(source), "genes_n": int(group["gene_symbol"].nunique())})
            sources.sort(key=lambda row: row["id"])
        return {
            "taxonomy_version": config.get("version"),
            "status": config.get("status"),
            "domains": output_domains,
            "sources": sources,
            "coverage": self.functional_coverage(),
            "provenance_note": config.get("provenance_note"),
        }

    def annotations(self, gene_symbol: str) -> dict[str, Any]:
        symbol, _ = self._gene_or_raise(gene_symbol)
        formal = self.formal_annotation_frame()
        formal_sub = formal[formal["gene_symbol"].astype(str).str.upper() == symbol].copy() if not formal.empty else formal
        mapping = self.domain_mapping_frame()
        mapped = mapping[mapping["gene_symbol"].astype(str).str.upper() == symbol].copy() if not mapping.empty else mapping

        domain_rows = mapped[mapped["annotation_type"] == "mcl_domain"] if not mapped.empty else mapped
        sub_rows = mapped[mapped["annotation_type"] == "mcl_subdomain"] if not mapped.empty else mapped
        domains = sorted(set(domain_rows["annotation_label_ru"].dropna().astype(str))) if not domain_rows.empty else []
        subdomains = sorted(set(sub_rows["annotation_label_ru"].dropna().astype(str))) if not sub_rows.empty else []
        domain_details = _records(domain_rows.drop_duplicates(["annotation_id", "source", "source_id"])) if not domain_rows.empty else []

        return {
            "gene_symbol": symbol,
            "status": "partial" if formal_sub.empty else "partial_formal_mcl",
            "mcl_domains": domains,
            "subdomains": subdomains,
            "protein_classes": [],
            "compartments": [],
            "hallmarks": [],
            "mcl_domain_details": domain_details,
            "formal_annotations": _records(formal_sub.sort_values(["source", "term_name"]).head(250)) if not formal_sub.empty else [],
            "coverage": self.functional_coverage(),
            "taxonomy_version": self._functional_config().get("version"),
            "note": (
                "MCL Functional Domains are rule-based projections of formal GO/Reactome/KEGG/CORUM terms observed in the current MCL enrichment layer; every assignment retains its source term. "
                "Coverage is intentionally partial until a complete gene-to-ontology snapshot is ingested. Protein class, compartment and Hallmarks remain unassigned rather than inferred without a dedicated source."
            ),
        }
