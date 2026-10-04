from __future__ import annotations

import json
from functools import lru_cache
from typing import Any

import pandas as pd

from .functional_gene_explorer import FunctionalGeneExplorerStore
from .gene_explorer import _bool_series, _clean, _number_series, _records


class AnnotatedGeneExplorerStore(FunctionalGeneExplorerStore):
    """Functional Gene Explorer enriched by an optional local reference snapshot.

    The reference snapshot is created by scripts/build_gene_reference_snapshot.py.
    Direct MCL-derived annotations remain usable when that snapshot is absent.
    """

    @lru_cache(maxsize=1)
    def _reference_frame(self) -> pd.DataFrame:
        path = self.index_dir / "gene_reference.parquet"
        if not path.exists():
            return pd.DataFrame()
        frame = pd.read_parquet(path)
        if "gene_symbol" in frame.columns:
            frame["gene_symbol"] = frame["gene_symbol"].astype(str).str.upper()
        return frame

    @lru_cache(maxsize=1)
    def _reference_terms(self) -> pd.DataFrame:
        path = self.index_dir / "gene_reference_terms.parquet"
        if not path.exists():
            return pd.DataFrame()
        frame = pd.read_parquet(path)
        if "gene_symbol" in frame.columns:
            frame["gene_symbol"] = frame["gene_symbol"].astype(str).str.upper()
        return frame

    @lru_cache(maxsize=1)
    def formal_annotation_frame(self) -> pd.DataFrame:
        current = super().formal_annotation_frame().copy()
        reference = self._reference_terms().copy()
        if reference.empty:
            return current

        if "significant" not in reference.columns:
            reference["significant"] = False
        if "top_n" not in reference.columns:
            reference["top_n"] = None
        if "source_version" not in reference.columns:
            reference["source_version"] = None

        columns = ["gene_symbol", "source", "term_id", "term_name", "significant", "top_n", "source_version"]
        for col in columns:
            if col not in current.columns:
                current[col] = None
            if col not in reference.columns:
                reference[col] = None
        merged = pd.concat([current[columns], reference[columns]], ignore_index=True)
        merged["gene_symbol"] = merged["gene_symbol"].astype(str).str.upper()
        return merged.drop_duplicates(["gene_symbol", "source", "term_id", "term_name"], keep="first").reset_index(drop=True)

    def _rules_for(self, group_name: str) -> list[dict[str, Any]]:
        return list(self._functional_config().get(group_name) or [])

    @lru_cache(maxsize=1)
    def secondary_mapping_frame(self) -> pd.DataFrame:
        formal = self.formal_annotation_frame()
        columns = [
            "gene_symbol",
            "annotation_type",
            "annotation_id",
            "annotation_label_ru",
            "source",
            "source_id",
            "source_term_name",
            "source_version",
            "mapping_method",
            "taxonomy_version",
        ]
        if formal.empty:
            return pd.DataFrame(columns=columns)

        groups = (
            ("protein_classes", "protein_class"),
            ("compartments", "compartment"),
            ("hallmarks", "hallmark"),
        )
        rows: list[dict[str, Any]] = []
        seen: set[tuple[str, str, str, str, str]] = set()
        version = self._functional_config().get("version")

        for _, term in formal.iterrows():
            source = str(term.get("source") or "").strip()
            term_name = str(term.get("term_name") or "").strip()
            if not term_name:
                continue
            for config_key, annotation_type in groups:
                for rule in self._rules_for(config_key):
                    allowed_sources = {str(x).casefold() for x in (rule.get("sources") or [])}
                    if allowed_sources and source.casefold() not in allowed_sources:
                        continue
                    if not self._contains_keyword(term_name, list(rule.get("keywords") or [])):
                        continue
                    annotation_id = str(rule.get("id") or "").strip()
                    if not annotation_id:
                        continue
                    key = (
                        str(term.get("gene_symbol") or ""),
                        annotation_type,
                        annotation_id,
                        source,
                        str(term.get("term_id") or ""),
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                    rows.append(
                        {
                            "gene_symbol": str(term.get("gene_symbol") or "").upper(),
                            "annotation_type": annotation_type,
                            "annotation_id": annotation_id,
                            "annotation_label_ru": str(rule.get("label_ru") or annotation_id),
                            "source": source or None,
                            "source_id": term.get("term_id"),
                            "source_term_name": term.get("term_name"),
                            "source_version": term.get("source_version"),
                            "mapping_method": f"rule_based_{annotation_type}_projection",
                            "taxonomy_version": version,
                        }
                    )
        return pd.DataFrame(rows, columns=columns) if rows else pd.DataFrame(columns=columns)

    def _catalog_with_taxonomy(self, catalog: pd.DataFrame) -> pd.DataFrame:
        out = super()._catalog_with_taxonomy(catalog)
        if out.empty:
            return out
        mapping = self.secondary_mapping_frame()

        def grouped(annotation_type: str) -> dict[str, str]:
            sub = mapping[mapping["annotation_type"] == annotation_type] if not mapping.empty else mapping
            if sub.empty:
                return {}
            return {
                str(gene): json.dumps(sorted(set(group["annotation_id"].astype(str))), ensure_ascii=False)
                for gene, group in sub.groupby("gene_symbol")
            }

        symbols = out["gene_symbol"].astype(str).str.upper()
        for annotation_type, column in (
            ("protein_class", "protein_classes_json"),
            ("compartment", "compartments_json"),
            ("hallmark", "hallmarks_json"),
        ):
            mapping_dict = grouped(annotation_type)
            out[column] = symbols.map(mapping_dict).fillna("[]")

        reference = self._reference_frame()
        if not reference.empty:
            keep = [
                c
                for c in (
                    "gene_symbol",
                    "gene_name",
                    "aliases_json",
                    "entrez_gene_id",
                    "ensembl_gene_ids_json",
                    "uniprot_swissprot_ids_json",
                    "type_of_gene",
                )
                if c in reference.columns
            ]
            out = out.merge(reference[keep].drop_duplicates("gene_symbol"), on="gene_symbol", how="left")
        return out

    @lru_cache(maxsize=1)
    def catalog_frame(self) -> pd.DataFrame:
        materialized = self.index_dir / "gene_catalog.parquet"
        required = {"mcl_domains_json", "protein_classes_json", "compartments_json", "hallmarks_json"}
        if materialized.exists():
            frame = pd.read_parquet(materialized)
            if required.issubset(frame.columns):
                return frame
        return self._catalog_with_taxonomy(self._build_catalog())

    def materialize_indexes(self) -> dict[str, Any]:
        manifest = super().materialize_indexes()
        annotations_path = self.index_dir / "gene_annotations.parquet"
        secondary = self.secondary_mapping_frame()
        if annotations_path.exists():
            primary = pd.read_parquet(annotations_path)
        else:
            primary = pd.DataFrame()
        if not secondary.empty:
            annotations = pd.concat([primary, secondary], ignore_index=True, sort=False)
            annotations = annotations.drop_duplicates(
                ["gene_symbol", "annotation_type", "annotation_id", "source", "source_id"], keep="first"
            )
            annotations.to_parquet(annotations_path, index=False, compression="zstd")
        manifest.update(
            {
                "secondary_mapping_rows_n": int(len(secondary)),
                "reference_genes_n": int(self._reference_frame()["gene_symbol"].nunique()) if not self._reference_frame().empty else 0,
                "reference_term_rows_n": int(len(self._reference_terms())),
            }
        )
        (self.index_dir / "manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        self.catalog_frame.cache_clear()
        return manifest

    def identity(self, gene_symbol: str) -> dict[str, Any]:
        payload = super().identity(gene_symbol)
        symbol = str(payload["identity"]["gene_symbol"]).upper()
        reference = self._reference_frame()
        if reference.empty:
            return payload
        hit = reference[reference["gene_symbol"].astype(str).str.upper() == symbol]
        if hit.empty:
            return payload
        row = hit.iloc[0]
        identity = dict(payload["identity"])
        identity["gene_name"] = row.get("gene_name") or identity.get("gene_name")
        identity["aliases"] = self._json_list(row.get("aliases_json"))
        identity["entrez_gene_id"] = row.get("entrez_gene_id")
        ensembl = self._json_list(row.get("ensembl_gene_ids_json"))
        uniprot = self._json_list(row.get("uniprot_swissprot_ids_json"))
        if not identity.get("ensembl_gene_id") and ensembl:
            identity["ensembl_gene_id"] = ensembl[0]
        if not identity.get("uniprot_id") and uniprot:
            identity["uniprot_id"] = uniprot[0]
        identity["ensembl_gene_ids"] = ensembl
        identity["uniprot_swissprot_ids"] = uniprot
        identity["type_of_gene"] = row.get("type_of_gene")
        identity["reference_source"] = row.get("aggregator")
        identity["reference_retrieved_at"] = row.get("retrieved_at")
        payload["identity"] = _clean(identity)
        return payload

    @staticmethod
    def _json_list(value: Any) -> list[str]:
        try:
            parsed = json.loads(str(value or "[]"))
        except (json.JSONDecodeError, TypeError):
            return []
        if not isinstance(parsed, list):
            return []
        return [str(x) for x in parsed if x is not None]

    def suggest(self, query: str, limit: int = 12) -> list[dict[str, Any]]:
        q = str(query or "").strip().casefold()
        if not q:
            return []
        frame = self.catalog_frame().copy()
        symbols = frame["gene_symbol"].astype(str)
        mask = symbols.str.casefold().str.contains(q, regex=False)
        if "gene_name" in frame.columns:
            mask = mask | frame["gene_name"].fillna("").astype(str).str.casefold().str.contains(q, regex=False)
        if "aliases_json" in frame.columns:
            mask = mask | frame["aliases_json"].fillna("").astype(str).str.casefold().str.contains(q, regex=False)
        frame = frame[mask].copy()
        if frame.empty:
            return []
        frame["_rank"] = frame.apply(
            lambda row: 0
            if str(row.get("gene_symbol") or "").casefold() == q
            else (1 if str(row.get("gene_symbol") or "").casefold().startswith(q) else 2),
            axis=1,
        )
        frame = frame.sort_values(["_rank", "gene_symbol"]).head(max(1, min(int(limit), 30)))
        fields = [
            c
            for c in (
                "gene_symbol",
                "gene_name",
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
        page: int = 1,
        page_size: int = 100,
        sort_by: str = "best_delta_gene_effect",
        sort_order: str = "asc",
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

        allowed_sort = {
            "gene_symbol", "gene_name", "best_delta_gene_effect", "best_model_gene_effect", "best_q_value",
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
                "reference_coverage": self.reference_coverage(),
                "items": _records(page_frame),
            }
        )

    def reference_coverage(self) -> dict[str, Any]:
        catalog = self.catalog_frame()
        reference = self._reference_frame()
        total = int(catalog["gene_symbol"].nunique()) if not catalog.empty else 0
        resolved = int(reference["gene_symbol"].nunique()) if not reference.empty else 0
        return {
            "resolved_genes_n": resolved,
            "gene_universe_n": total,
            "coverage_fraction": (resolved / total) if total else None,
            "available": not reference.empty,
        }

    def facets(self) -> dict[str, Any]:
        payload = super().facets()
        mapping = self.secondary_mapping_frame()
        config = self._functional_config()

        def facet_group(config_key: str, annotation_type: str) -> list[dict[str, Any]]:
            sub = mapping[mapping["annotation_type"] == annotation_type] if not mapping.empty else mapping
            counts = sub.groupby("annotation_id")["gene_symbol"].nunique().to_dict() if not sub.empty else {}
            return [
                {
                    "id": rule.get("id"),
                    "label_ru": rule.get("label_ru"),
                    "genes_n": int(counts.get(rule.get("id"), 0)),
                }
                for rule in config.get(config_key) or []
            ]

        payload["protein_classes"] = facet_group("protein_classes", "protein_class")
        payload["compartments"] = facet_group("compartments", "compartment")
        payload["hallmarks"] = facet_group("hallmarks", "hallmark")
        payload["reference_coverage"] = self.reference_coverage()
        return payload

    def annotations(self, gene_symbol: str) -> dict[str, Any]:
        payload = super().annotations(gene_symbol)
        symbol = str(payload["gene_symbol"]).upper()
        mapping = self.secondary_mapping_frame()
        sub = mapping[mapping["gene_symbol"].astype(str).str.upper() == symbol].copy() if not mapping.empty else mapping

        for annotation_type, output_key in (
            ("protein_class", "protein_classes"),
            ("compartment", "compartments"),
            ("hallmark", "hallmarks"),
        ):
            rows = sub[sub["annotation_type"] == annotation_type] if not sub.empty else sub
            payload[output_key] = sorted(set(rows["annotation_label_ru"].dropna().astype(str))) if not rows.empty else []
            payload[f"{annotation_type}_details"] = _records(rows) if not rows.empty else []
        payload["reference_coverage"] = self.reference_coverage()
        payload["reference_available"] = not self._reference_frame().empty
        return payload
