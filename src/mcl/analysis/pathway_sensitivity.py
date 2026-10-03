from __future__ import annotations

import argparse
import json
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from mcl.analysis.pathways import (
    build_pathway_inputs,
    load_pathway_config,
    normalize_gprofiler_result,
)
from mcl.sources.gprofiler import DEFAULT_PROFILE_ENDPOINT, run_gprofiler_profile


PATHWAY_SENSITIVITY_VERSION = "m3.3.1-pathway-sensitivity-v0.1"


def _json_list(value: Any) -> list[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    if isinstance(value, list):
        return [str(x) for x in value]
    text = str(value).strip()
    if not text:
        return []
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        return []
    if not isinstance(parsed, list):
        return []
    return [str(x) for x in parsed]


def build_term_stability(
    enrichment: pd.DataFrame,
    thresholds: list[int],
) -> pd.DataFrame:
    """Summarize whether each enriched term is significant across top-N cut-offs."""
    thresholds = sorted({int(x) for x in thresholds})
    base_columns = [
        "source",
        "term_id",
        "term_name",
        "significant_thresholds_n",
        "significant_all_thresholds",
        "thresholds_significant",
        "min_p_value_adjusted",
    ]
    dynamic_columns: list[str] = []
    for top_n in thresholds:
        dynamic_columns.extend(
            [
                f"significant_top{top_n}",
                f"p_value_adjusted_top{top_n}",
                f"intersection_size_top{top_n}",
                f"intersecting_gene_symbols_top{top_n}",
            ]
        )
    if enrichment.empty:
        return pd.DataFrame(columns=base_columns + dynamic_columns)

    records: list[dict[str, Any]] = []
    for (source, term_id), rows in enrichment.groupby(["source", "term_id"], dropna=False):
        first = rows.iloc[0]
        record: dict[str, Any] = {
            "source": source,
            "term_id": term_id,
            "term_name": first.get("term_name"),
        }
        significant_at: list[int] = []
        p_values: list[float] = []
        for top_n in thresholds:
            subset = rows.loc[pd.to_numeric(rows["top_n"], errors="coerce") == top_n]
            if subset.empty:
                record[f"significant_top{top_n}"] = False
                record[f"p_value_adjusted_top{top_n}"] = None
                record[f"intersection_size_top{top_n}"] = 0
                record[f"intersecting_gene_symbols_top{top_n}"] = "[]"
                continue
            row = subset.sort_values("p_value_adjusted", na_position="last").iloc[0]
            significant = bool(row.get("significant", False))
            p_value = pd.to_numeric(
                pd.Series([row.get("p_value_adjusted")]), errors="coerce"
            ).iloc[0]
            genes = sorted(set(_json_list(row.get("intersecting_gene_symbols_json"))))
            record[f"significant_top{top_n}"] = significant
            record[f"p_value_adjusted_top{top_n}"] = (
                float(p_value) if pd.notna(p_value) else None
            )
            record[f"intersection_size_top{top_n}"] = int(
                pd.to_numeric(
                    pd.Series([row.get("intersection_size", 0)]), errors="coerce"
                )
                .fillna(0)
                .iloc[0]
            )
            record[f"intersecting_gene_symbols_top{top_n}"] = json.dumps(
                genes, ensure_ascii=False
            )
            if significant:
                significant_at.append(top_n)
            if pd.notna(p_value):
                p_values.append(float(p_value))

        record["significant_thresholds_n"] = len(significant_at)
        record["significant_all_thresholds"] = len(significant_at) == len(thresholds)
        record["thresholds_significant"] = ";".join(str(x) for x in significant_at)
        record["min_p_value_adjusted"] = min(p_values) if p_values else None
        records.append(record)

    result = pd.DataFrame(records)
    return result.sort_values(
        ["significant_thresholds_n", "min_p_value_adjusted", "source", "term_name"],
        ascending=[False, True, True, True],
        na_position="last",
    ).reset_index(drop=True)


def build_gene_stability(
    recurrence_by_threshold: dict[int, pd.DataFrame],
    thresholds: list[int],
    *,
    minimum_comparisons: int,
) -> pd.DataFrame:
    """Summarize recurrent-gene persistence across top-N cut-offs."""
    thresholds = sorted({int(x) for x in thresholds})
    query_by_threshold: dict[int, pd.DataFrame] = {}
    all_genes: set[str] = set()
    for top_n in thresholds:
        frame = recurrence_by_threshold.get(top_n, pd.DataFrame()).copy()
        if frame.empty:
            query = frame
        else:
            query = frame.loc[
                pd.to_numeric(frame["comparisons_n"], errors="coerce")
                >= int(minimum_comparisons)
            ].copy()
        query_by_threshold[top_n] = query
        if not query.empty:
            all_genes.update(query["gene_symbol"].dropna().astype(str))

    records: list[dict[str, Any]] = []
    for gene in sorted(all_genes):
        record: dict[str, Any] = {"gene_symbol": gene}
        present_at: list[int] = []
        for top_n in thresholds:
            rows = query_by_threshold[top_n]
            hit = (
                rows.loc[rows["gene_symbol"].astype(str) == gene]
                if not rows.empty
                else rows
            )
            present = not hit.empty
            record[f"recurrent_top{top_n}"] = present
            if present:
                present_at.append(top_n)
                row = hit.iloc[0]
                record[f"comparisons_n_top{top_n}"] = int(row["comparisons_n"])
                record[f"best_rank_top{top_n}"] = int(row["best_rank"])
                record[f"mean_rank_top{top_n}"] = float(row["mean_rank"])
                record[f"mean_delta_gene_effect_top{top_n}"] = float(
                    row["mean_delta_gene_effect"]
                )
            else:
                record[f"comparisons_n_top{top_n}"] = None
                record[f"best_rank_top{top_n}"] = None
                record[f"mean_rank_top{top_n}"] = None
                record[f"mean_delta_gene_effect_top{top_n}"] = None

        record["thresholds_n"] = len(present_at)
        record["present_all_thresholds"] = len(present_at) == len(thresholds)
        record["thresholds_present"] = ";".join(str(x) for x in present_at)
        records.append(record)

    if not records:
        return pd.DataFrame(
            columns=["gene_symbol", "thresholds_n", "present_all_thresholds"]
        )
    result = pd.DataFrame(records)
    return result.sort_values(
        ["thresholds_n", "gene_symbol"], ascending=[False, True]
    ).reset_index(drop=True)


def _load_sensitivity_config(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise ValueError("Sensitivity configuration must be a YAML mapping")
    thresholds = sorted({int(x) for x in data.get("thresholds") or []})
    if not thresholds or any(x <= 0 for x in thresholds):
        raise ValueError("Sensitivity thresholds must contain positive integers")
    data["thresholds"] = thresholds
    return data


def _load_profile_endpoint(root: Path) -> str:
    path = root / "config/source_versions.yaml"
    if not path.exists():
        return DEFAULT_PROFILE_ENDPOINT
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    gprofiler = data.get("gprofiler") or {}
    return str(gprofiler.get("profile_endpoint") or DEFAULT_PROFILE_ENDPOINT)


def _qc_record(
    severity: str,
    check: str,
    entity_id: str,
    observed: str,
    expected: str,
    message: str,
) -> dict[str, str]:
    return {
        "severity": severity,
        "check": check,
        "entity_id": entity_id,
        "observed": observed,
        "expected": expected,
        "message": message,
    }


def analyze_pathway_sensitivity(
    root: Path,
    *,
    pathway_config_path: Path | None = None,
    sensitivity_config_path: Path | None = None,
) -> dict[str, Any]:
    root = root.resolve()
    pathway_config_path = pathway_config_path or (root / "config/pathways.yaml")
    sensitivity_config_path = sensitivity_config_path or (
        root / "config/pathway_sensitivity.yaml"
    )
    pathway_cfg = load_pathway_config(pathway_config_path)
    sensitivity_cfg = _load_sensitivity_config(sensitivity_config_path)

    thresholds = sensitivity_cfg["thresholds"]
    analysis_set = str(sensitivity_cfg.get("analysis_set") or "recurrent")
    baseline_threshold = int(sensitivity_cfg.get("baseline_threshold", 100))
    analysis_sets_cfg = pathway_cfg.get("analysis_sets") or {}
    if analysis_set not in analysis_sets_cfg:
        raise ValueError(f"Unknown analysis_set in sensitivity config: {analysis_set}")
    minimum_comparisons = int(
        (analysis_sets_cfg.get(analysis_set) or {}).get("minimum_comparisons", 1)
    )

    resolution_path = root / str(
        sensitivity_cfg.get(
            "identifier_resolution",
            "data/processed/pathways/identifier_resolution.tsv",
        )
    )
    if not resolution_path.exists():
        raise FileNotFoundError(
            f"Missing canonical identifier-resolution table: {resolution_path}. "
            "Run M3.3 analyze-pathways first."
        )
    identifier_resolution = pd.read_csv(
        resolution_path, sep="\t", dtype=str, keep_default_na=False
    )
    required_resolution = {"gene_symbol", "ensembl_gene_id", "mapping_status"}
    missing_resolution = required_resolution - set(identifier_resolution.columns)
    if missing_resolution:
        raise ValueError(
            "Identifier-resolution table is missing columns: "
            + ", ".join(sorted(missing_resolution))
        )
    resolution_by_symbol = identifier_resolution.set_index("gene_symbol", drop=False)

    enrichment_cfg = pathway_cfg.get("enrichment") or {}
    sources = [
        str(x)
        for x in enrichment_cfg.get("sources")
        or ["GO:BP", "REAC", "KEGG", "CORUM"]
    ]
    method = str(enrichment_cfg.get("significance_threshold_method") or "g_SCS")
    user_threshold = float(enrichment_cfg.get("user_threshold", 0.05))
    excluded_term_ids = {
        str(x) for x in enrichment_cfg.get("excluded_term_ids") or []
    }
    organism = str(pathway_cfg.get("organism") or "hsapiens")
    profile_endpoint = _load_profile_endpoint(root)

    recurrence_by_threshold: dict[int, pd.DataFrame] = {}
    candidate_outputs: list[pd.DataFrame] = []
    enrichment_outputs: list[pd.DataFrame] = []
    raw_responses: dict[int, dict[str, Any]] = {}
    per_threshold_meta: dict[str, dict[str, Any]] = {}
    qc: list[dict[str, str]] = []
    canonical_background_symbols: set[str] | None = None
    canonical_background_ensg: list[str] | None = None

    for top_n in thresholds:
        cfg = deepcopy(pathway_cfg)
        selection = dict(cfg.get("candidate_selection") or {})
        selection["top_n_per_comparison"] = int(top_n)
        cfg["candidate_selection"] = selection

        candidate_long, recurrence, candidate_sets, backgrounds, _, build_meta = (
            build_pathway_inputs(root, cfg)
        )
        recurrence_by_threshold[top_n] = recurrence
        if not candidate_long.empty:
            candidate_long = candidate_long.copy()
            candidate_long.insert(0, "top_n", top_n)
            candidate_outputs.append(candidate_long)

        query_genes = set(candidate_sets.get(analysis_set, set()))
        background_symbols = set(backgrounds.get(analysis_set, set()))
        if not query_genes:
            qc.append(
                _qc_record(
                    "ERROR",
                    "sensitivity_query_nonempty",
                    f"top{top_n}",
                    "0",
                    ">0",
                    "No recurrent pathway query genes were produced at this threshold.",
                )
            )
            continue

        if canonical_background_symbols is None:
            canonical_background_symbols = background_symbols
        elif background_symbols != canonical_background_symbols:
            qc.append(
                _qc_record(
                    "ERROR",
                    "sensitivity_background_stability",
                    f"top{top_n}",
                    str(len(background_symbols)),
                    str(len(canonical_background_symbols)),
                    "Eligibility-matched background changed with top-N cut-off; "
                    "background must be threshold-independent.",
                )
            )

        missing_symbols = sorted(
            (query_genes | background_symbols)
            - set(resolution_by_symbol.index.astype(str))
        )
        if missing_symbols:
            raise ValueError(
                "Sensitivity genes are absent from canonical identifier resolution: "
                + ", ".join(missing_symbols[:20])
            )

        query_rows = resolution_by_symbol.loc[sorted(query_genes)].copy()
        background_rows = resolution_by_symbol.loc[sorted(background_symbols)].copy()
        query_unresolved = query_rows.loc[
            query_rows["ensembl_gene_id"].astype(str) == ""
        ]
        background_unresolved = background_rows.loc[
            background_rows["ensembl_gene_id"].astype(str) == ""
        ]
        if not query_unresolved.empty:
            genes = query_unresolved["gene_symbol"].astype(str).tolist()
            qc.append(
                _qc_record(
                    "ERROR",
                    "sensitivity_query_identifier_mapping",
                    f"top{top_n}",
                    ", ".join(genes[:20]),
                    "0 unresolved query genes",
                    "One or more recurrent sensitivity genes lack canonical ENSG mapping.",
                )
            )
            continue

        query_ensg = sorted(set(query_rows["ensembl_gene_id"].astype(str)) - {""})
        background_ensg = sorted(
            set(background_rows["ensembl_gene_id"].astype(str)) - {""}
        )
        if canonical_background_ensg is None:
            canonical_background_ensg = background_ensg
        elif background_ensg != canonical_background_ensg:
            qc.append(
                _qc_record(
                    "ERROR",
                    "sensitivity_resolved_background_stability",
                    f"top{top_n}",
                    str(len(background_ensg)),
                    str(len(canonical_background_ensg)),
                    "Resolved canonical ENSG background changed with top-N cut-off.",
                )
            )

        if not set(query_ensg).issubset(set(background_ensg)):
            raise ValueError(
                f"Top-{top_n} recurrent query is not contained in canonical background"
            )

        ensg_to_symbol = {
            str(row.ensembl_gene_id): str(row.gene_symbol)
            for row in background_rows.itertuples(index=False)
            if str(row.ensembl_gene_id)
        }
        response = run_gprofiler_profile(
            query_ensg,
            background_ensg,
            organism=organism,
            sources=sources,
            significance_threshold_method=method,
            user_threshold=user_threshold,
            all_results=bool(enrichment_cfg.get("all_results", True)),
            domain_scope=str(enrichment_cfg.get("domain_scope") or "custom"),
            no_evidences=bool(enrichment_cfg.get("no_evidences", False)),
            numeric_ns="",
            endpoint=profile_endpoint,
        )
        raw_responses[top_n] = response
        normalized = normalize_gprofiler_result(
            response,
            f"top{top_n}_{analysis_set}",
            excluded_term_ids=excluded_term_ids,
            input_id_to_symbol=ensg_to_symbol,
        )
        if not normalized.empty:
            normalized.insert(0, "top_n", top_n)
            enrichment_outputs.append(normalized)

        sig_n = (
            int(normalized["significant"].fillna(False).astype(bool).sum())
            if not normalized.empty
            else 0
        )
        observed_query_sizes = sorted(
            {
                int(x)
                for x in normalized.get("query_size", pd.Series(dtype=float))
                .dropna()
                .tolist()
            }
        )
        observed_domains = sorted(
            {
                int(x)
                for x in normalized.get(
                    "effective_domain_size", pd.Series(dtype=float)
                )
                .dropna()
                .tolist()
            }
        )
        per_threshold_meta[str(top_n)] = {
            "top_n": top_n,
            "query_genes_n": len(query_genes),
            "query_identifiers_n": len(query_ensg),
            "background_genes_n": len(background_symbols),
            "background_identifiers_n": len(background_ensg),
            "background_unresolved_symbols_n": int(len(background_unresolved)),
            "significant_terms_n": sig_n,
            "gprofiler_query_sizes": observed_query_sizes,
            "gprofiler_effective_domain_sizes": observed_domains,
            "build_meta": build_meta,
        }

        qc.append(
            _qc_record(
                "INFO" if len(query_ensg) == len(query_genes) else "WARNING",
                "sensitivity_query_identifier_mapping",
                f"top{top_n}",
                f"genes={len(query_genes)}; ENSG={len(query_ensg)}",
                "one ENSG per query gene",
                "Canonical identifier mapping for the recurrent sensitivity query.",
            )
        )
        qc.append(
            _qc_record(
                "WARNING" if len(background_unresolved) else "INFO",
                "sensitivity_background_identifier_mapping",
                f"top{top_n}",
                str(len(background_unresolved)),
                "0 unresolved background genes",
                "Unresolved eligible background symbols are excluded before g:Profiler.",
            )
        )
        if observed_query_sizes:
            qc.append(
                _qc_record(
                    "INFO"
                    if observed_query_sizes == [len(query_ensg)]
                    else "WARNING",
                    "sensitivity_gprofiler_query_size_match",
                    f"top{top_n}",
                    str(observed_query_sizes),
                    str(len(query_ensg)),
                    "g:Profiler query size should match the submitted canonical ENSG query.",
                )
            )
        if observed_domains:
            qc.append(
                _qc_record(
                    "INFO"
                    if observed_domains == [len(background_ensg)]
                    else "WARNING",
                    "sensitivity_gprofiler_background_size_match",
                    f"top{top_n}",
                    str(observed_domains),
                    str(len(background_ensg)),
                    "g:Profiler effective domain should match the canonical custom background.",
                )
            )

    candidate_all = (
        pd.concat(candidate_outputs, ignore_index=True)
        if candidate_outputs
        else pd.DataFrame()
    )
    enrichment_all = (
        pd.concat(enrichment_outputs, ignore_index=True)
        if enrichment_outputs
        else pd.DataFrame()
    )
    gene_stability = build_gene_stability(
        recurrence_by_threshold,
        thresholds,
        minimum_comparisons=minimum_comparisons,
    )
    term_stability = build_term_stability(enrichment_all, thresholds)

    expected_thresholds = set(thresholds)
    observed_thresholds = (
        {
            int(x)
            for x in pd.to_numeric(enrichment_all.get("top_n"), errors="coerce")
            .dropna()
            .tolist()
        }
        if not enrichment_all.empty
        else set()
    )
    qc.append(
        _qc_record(
            "INFO" if expected_thresholds == observed_thresholds else "ERROR",
            "sensitivity_thresholds_complete",
            "M3.3.1",
            ",".join(str(x) for x in sorted(observed_thresholds)),
            ",".join(str(x) for x in thresholds),
            "All configured top-N thresholds should complete enrichment.",
        )
    )

    baseline_report = root / str(
        sensitivity_cfg.get(
            "baseline_significant_report",
            "outputs/reports/M3_3_pathway_significant.tsv",
        )
    )
    baseline_match: bool | None = None
    if (
        baseline_threshold in thresholds
        and baseline_report.exists()
        and not enrichment_all.empty
    ):
        baseline = pd.read_csv(baseline_report, sep="\t", dtype=str).fillna("")
        if "analysis_set" in baseline.columns:
            baseline = baseline.loc[baseline["analysis_set"] == analysis_set]
        baseline_pairs = set(
            zip(
                baseline.get("source", pd.Series(dtype=str)),
                baseline.get("term_id", pd.Series(dtype=str)),
            )
        )
        current = enrichment_all.loc[
            (pd.to_numeric(enrichment_all["top_n"], errors="coerce") == baseline_threshold)
            & enrichment_all["significant"].fillna(False).astype(bool)
        ]
        current_pairs = set(
            zip(current["source"].astype(str), current["term_id"].astype(str))
        )
        baseline_match = baseline_pairs == current_pairs
        qc.append(
            _qc_record(
                "INFO" if baseline_match else "WARNING",
                "sensitivity_top100_baseline_reproduction",
                f"top{baseline_threshold}",
                f"baseline={len(baseline_pairs)}; current={len(current_pairs)}",
                "identical significant source/term_id set",
                "Top-100 sensitivity result is compared with the locked M3.3 baseline.",
            )
        )

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    processed = root / "data/processed/pathways_sensitivity"
    reports = root / "outputs/reports"
    qc_dir = root / "outputs/qc"
    raw_dir = root / "data/raw/gprofiler_sensitivity" / timestamp
    for directory in (processed, reports, qc_dir, raw_dir):
        directory.mkdir(parents=True, exist_ok=True)

    candidate_path = processed / "candidate_topn_long.tsv"
    enrichment_path = processed / "enrichment_all.tsv"
    gene_stability_path = processed / "gene_stability.tsv"
    term_stability_path = processed / "term_stability.tsv"
    candidate_all.to_csv(candidate_path, sep="\t", index=False)
    enrichment_all.to_csv(enrichment_path, sep="\t", index=False)
    gene_stability.to_csv(gene_stability_path, sep="\t", index=False)
    term_stability.to_csv(term_stability_path, sep="\t", index=False)

    parquet_written = True
    try:
        candidate_all.to_parquet(
            processed / "candidate_topn_long.parquet", index=False
        )
        enrichment_all.to_parquet(processed / "enrichment_all.parquet", index=False)
        gene_stability.to_parquet(processed / "gene_stability.parquet", index=False)
        term_stability.to_parquet(processed / "term_stability.parquet", index=False)
    except ImportError:
        parquet_written = False

    stable_genes = gene_stability.loc[
        gene_stability.get("present_all_thresholds", pd.Series(dtype=bool))
        .fillna(False)
        .astype(bool)
    ].copy()
    stable_terms = term_stability.loc[
        term_stability.get("significant_all_thresholds", pd.Series(dtype=bool))
        .fillna(False)
        .astype(bool)
    ].copy()
    stable_genes_path = reports / "M3_3_1_stable_recurrent_genes.tsv"
    stable_terms_path = reports / "M3_3_1_stable_pathways.tsv"
    stable_genes.to_csv(stable_genes_path, sep="\t", index=False)
    stable_terms.to_csv(stable_terms_path, sep="\t", index=False)

    for top_n, payload in raw_responses.items():
        (raw_dir / f"top{top_n}_{analysis_set}.json").write_text(
            json.dumps(payload, ensure_ascii=False, indent=2, default=str),
            encoding="utf-8",
        )

    qc_path = qc_dir / "pathway_sensitivity_qc.tsv"
    qc_json_path = qc_dir / "pathway_sensitivity_qc.json"
    pd.DataFrame(qc).to_csv(qc_path, sep="\t", index=False)
    qc_json_path.write_text(
        json.dumps(qc, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    meta = {
        "analysis_version": PATHWAY_SENSITIVITY_VERSION,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "pathway_config": str(pathway_config_path),
        "sensitivity_config": str(sensitivity_config_path),
        "thresholds": thresholds,
        "analysis_set": analysis_set,
        "minimum_comparisons": minimum_comparisons,
        "baseline_threshold": baseline_threshold,
        "baseline_top100_matches_m3_3": baseline_match,
        "profile_endpoint": profile_endpoint,
        "sources": sources,
        "significance_threshold_method": method,
        "user_threshold": user_threshold,
        "per_threshold": per_threshold_meta,
        "stable_recurrent_genes_n": int(len(stable_genes)),
        "stable_significant_terms_n": int(len(stable_terms)),
        "parquet_written": parquet_written,
        "outputs": {
            "candidate_topn_long": str(candidate_path),
            "enrichment_all": str(enrichment_path),
            "gene_stability": str(gene_stability_path),
            "term_stability": str(term_stability_path),
            "stable_recurrent_genes": str(stable_genes_path),
            "stable_pathways": str(stable_terms_path),
            "qc_tsv": str(qc_path),
            "qc_json": str(qc_json_path),
            "raw_dir": str(raw_dir),
        },
    }
    meta_path = reports / "M3_3_1_pathway_sensitivity_meta.json"
    meta_path.write_text(
        json.dumps(meta, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8",
    )

    errors = [row for row in qc if row["severity"] == "ERROR"]
    warnings = [row for row in qc if row["severity"] == "WARNING"]
    return {
        "meta": meta,
        "qc": qc,
        "errors_n": len(errors),
        "warnings_n": len(warnings),
        "stable_genes": stable_genes,
        "stable_terms": stable_terms,
        "gene_stability": gene_stability,
        "term_stability": term_stability,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "M3.3.1 pathway-enrichment sensitivity analysis across top-N cut-offs."
        )
    )
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--pathway-config", type=Path, default=None)
    parser.add_argument("--sensitivity-config", type=Path, default=None)
    args = parser.parse_args()

    result = analyze_pathway_sensitivity(
        args.root,
        pathway_config_path=args.pathway_config,
        sensitivity_config_path=args.sensitivity_config,
    )
    meta = result["meta"]
    print("M3.3.1 pathway sensitivity completed")
    for top_n in meta["thresholds"]:
        details = meta["per_threshold"].get(str(top_n), {})
        print(
            f"Top-{top_n}: recurrent query={details.get('query_genes_n', 0)} | "
            f"significant terms={details.get('significant_terms_n', 0)}"
        )
    print(
        f"Stable recurrent genes: {meta['stable_recurrent_genes_n']} | "
        f"stable significant terms: {meta['stable_significant_terms_n']}"
    )
    print(
        f"QC ERROR: {result['errors_n']} | WARNING: {result['warnings_n']} | "
        f"TOTAL: {len(result['qc'])}"
    )
    if result["errors_n"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
