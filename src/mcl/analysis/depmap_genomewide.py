from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from mcl.analysis.depmap import (
    _target_variant_label,
    build_context_membership,
    build_kras_sensitivity_membership,
    cliffs_delta,
    load_context_config,
)
from mcl.sources.depmap import (
    genomewide_gene_schema,
    load_genomewide_dependency_subset_with_background,
    load_genomewide_matrix_subset,
    load_model_metadata,
    load_relevant_mutations,
    load_sequenced_model_ids,
    parse_gene_label,
    resolve_depmap_files,
)
from mcl.utils.hash import sha256_file


GENOMEWIDE_PARSER_VERSION = "depmap-m3.2-genomewide-v0.1"
VALID_COMPARISONS = {"primary", "kras-wt", "other-kras"}


def _gene_metadata(path: Path) -> pd.DataFrame:
    header = pd.read_csv(path, nrows=0).columns.tolist()
    rows = []
    for raw in header[1:]:
        symbol, entrez = parse_gene_label(raw)
        rows.append({"gene_symbol": symbol, "entrez_gene_id": entrez, "source_column": raw})
    return pd.DataFrame(rows)


def _mutation_genes_for_context(cfg: dict) -> list[str]:
    mcfg = cfg.get("depmap", {}).get("mutation", {})
    genes: set[str] = set()
    if mcfg.get("gene"):
        genes.add(str(mcfg["gene"]))
    for rule in mcfg.get("defining_variants", []) or []:
        if rule.get("gene"):
            genes.add(str(rule["gene"]))
    return sorted(genes)


def _resolve_groups(
    root: Path,
    files: dict[str, Path],
    cancer_id: str,
    comparison: str,
) -> tuple[list[str], list[str], pd.DataFrame, str, str, dict]:
    contexts = load_context_config(root / "config/cancer_contexts.yaml")
    if cancer_id not in contexts:
        raise ValueError(f"Unknown cancer context: {cancer_id}")
    if comparison not in VALID_COMPARISONS:
        raise ValueError(
            f"Unsupported comparison {comparison!r}; choose one of {sorted(VALID_COMPARISONS)}"
        )

    cfg = contexts[cancer_id]
    models = load_model_metadata(files["models"])
    sequenced_models, profile_schema = load_sequenced_model_ids(files["omics_profiles"])

    if comparison == "primary":
        mutation_genes = _mutation_genes_for_context(cfg)
        mutations, mut_schema = load_relevant_mutations(files["mutations"], mutation_genes)
        memberships, audit = build_context_membership(
            models=models,
            mutations=mutations,
            sequenced_models=sequenced_models,
            contexts={cancer_id: cfg},
            mutation_schema=mut_schema,
        )
        groups = memberships[cancer_id]
        context_ids = list(groups.get("context", []))
        comparator_ids = list(groups.get("comparator", []))
        dep_cfg = cfg.get("depmap", {})
        context_definition = str(
            dep_cfg.get("context_definition", cfg.get("name", cancer_id))
        )
        comparator_definition = str(
            dep_cfg.get("comparator_definition", cfg.get("comparator", ""))
        )
        audit_rows = []
        for row in audit:
            payload = row.model_dump(mode="json")
            audit_rows.append(
                {
                    "cancer_id": cancer_id,
                    "comparison_type": comparison,
                    "model_id": payload["model_id"],
                    "cell_line_name": payload.get("cell_line_name"),
                    "group": payload["assigned_group"],
                    "sequencing_available": payload["sequencing_available"],
                    "molecular_status": payload["alteration_status"],
                    "assignment_reason": payload["assignment_reason"],
                }
            )
    else:
        mcfg = cfg.get("depmap", {}).get("mutation", {})
        if (
            str(mcfg.get("mode", "")) != "exact_protein_change_present"
            or str(mcfg.get("gene", "")).upper() != "KRAS"
        ):
            raise ValueError(
                f"{comparison} is only available for exact KRAS-variant contexts; {cancer_id} is not one"
            )
        mutations, mut_schema = load_relevant_mutations(files["mutations"], ["KRAS"])
        memberships, audit = build_kras_sensitivity_membership(
            models=models,
            mutations=mutations,
            sequenced_models=sequenced_models,
            contexts={cancer_id: cfg},
            mutation_schema=mut_schema,
        )
        groups = memberships[cancer_id]
        context_ids = list(groups["target_variant"])
        comparator_key = (
            "kras_wildtype_proxy" if comparison == "kras-wt" else "other_kras_driver_hotspot"
        )
        comparator_ids = list(groups[comparator_key])
        patterns = list(mcfg.get("protein_change_regexes", []) or [])
        target_variant = _target_variant_label(patterns)
        context_definition = f"same-disease models with KRAS {target_variant}"
        if comparison == "kras-wt":
            comparator_definition = (
                "same-disease models with WES/WGS coverage and no detected KRAS variant"
            )
        else:
            comparator_definition = (
                "same-disease models carrying a non-target KRAS variant annotated as driver/hotspot"
            )
        audit_rows = []
        for row in audit:
            payload = row.model_dump(mode="json")
            audit_rows.append(
                {
                    "cancer_id": cancer_id,
                    "comparison_type": comparison,
                    "model_id": payload["model_id"],
                    "cell_line_name": payload.get("cell_line_name"),
                    "group": payload["assigned_group"],
                    "sequencing_available": payload["sequencing_available"],
                    "molecular_status": payload["kras_status"],
                    "assignment_reason": payload["assignment_reason"],
                }
            )

    cohort = pd.DataFrame(audit_rows)
    if not cohort.empty:
        context_set = set(context_ids)
        comparator_set = set(comparator_ids)
        cohort["analysis_group"] = cohort["model_id"].map(
            lambda x: "context" if x in context_set else ("comparator" if x in comparator_set else "excluded")
        )
    return (
        context_ids,
        comparator_ids,
        cohort,
        context_definition,
        comparator_definition,
        {
            "omics_profiles_schema": profile_schema,
            "mutation_schema": {k: v for k, v in mut_schema.items() if k != "all_columns"},
        },
    )


def _safe_fraction(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    return numerator.astype(float).div(denominator.where(denominator > 0).astype(float))


def _subset_rows(matrix: pd.DataFrame, model_ids: list[str]) -> pd.DataFrame:
    available = [x for x in model_ids if x in matrix.index]
    if not available:
        return matrix.iloc[0:0].copy()
    return matrix.loc[available].copy()


def analyze_depmap_genomewide(
    root: Path,
    release: str,
    cancer_id: str,
    comparison: str = "primary",
    depmap_dir: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """Genome-wide DepMap comparison for one molecular cancer context.

    This is a discovery layer, not a target-ranking score.  Each gene retains
    separate observed statistics: dependency in context, dependency in comparator,
    effect size, uncertainty, multiple-testing correction, and broad dependency.
    """
    comparison = str(comparison).strip().lower()
    files = resolve_depmap_files(root, release, depmap_dir)
    (
        context_ids,
        comparator_ids,
        cohort,
        context_definition,
        comparator_definition,
        schema_meta,
    ) = _resolve_groups(root, files, cancer_id, comparison)

    selected_ids = list(dict.fromkeys(context_ids + comparator_ids))
    if not context_ids:
        raise ValueError(f"{cancer_id}: no context models are available for {comparison}")

    thresholds = yaml.safe_load((root / "config/thresholds.yaml").read_text(encoding="utf-8")) or {}
    dcfg = thresholds.get("depmap", {})
    gcfg = dcfg.get("genomewide", {}) or {}
    dep_threshold = float(dcfg.get("dependency_probability_threshold", 0.5))
    broad_threshold = float(dcfg.get("broad_dependency_fraction_warning", 0.8))
    min_context_n = int(gcfg.get("minimum_context_n_warning", dcfg.get("minimum_context_n_warning", 5)))
    min_comparator_n = int(
        gcfg.get("minimum_comparator_n_warning", dcfg.get("minimum_comparator_n_warning", 5))
    )
    min_test_n = int(gcfg.get("minimum_n_for_statistical_test", dcfg.get("minimum_n_for_statistical_test", 2)))

    gene_effect = load_genomewide_matrix_subset(files["gene_effect"], selected_ids)
    dependency, broad_fraction, broad_n = load_genomewide_dependency_subset_with_background(
        files["gene_dependency"],
        selected_ids,
        dependency_threshold=dep_threshold,
    )

    effect_genes = list(gene_effect.columns)
    dependency_genes = set(dependency.columns)
    genes = [g for g in effect_genes if g in dependency_genes]
    if not genes:
        raise ValueError("Gene Effect and Gene Dependency matrices have no common gene symbols")
    gene_effect = gene_effect[genes]
    dependency = dependency[genes]

    c_ge = _subset_rows(gene_effect, context_ids)
    r_ge = _subset_rows(gene_effect, comparator_ids)
    c_dp = _subset_rows(dependency, context_ids)
    r_dp = _subset_rows(dependency, comparator_ids)

    c_ge_n = c_ge.notna().sum(axis=0)
    r_ge_n = r_ge.notna().sum(axis=0)
    c_dp_n = c_dp.notna().sum(axis=0)
    r_dp_n = r_dp.notna().sum(axis=0)

    c_ge_median = c_ge.median(axis=0, skipna=True)
    r_ge_median = r_ge.median(axis=0, skipna=True) if len(r_ge) else pd.Series(np.nan, index=genes)
    c_ge_mean = c_ge.mean(axis=0, skipna=True)
    c_ge_q1 = c_ge.quantile(0.25, axis=0, numeric_only=True)
    c_ge_q3 = c_ge.quantile(0.75, axis=0, numeric_only=True)
    c_dep_numerator = (c_dp > dep_threshold).sum(axis=0)
    r_dep_numerator = (r_dp > dep_threshold).sum(axis=0) if len(r_dp) else pd.Series(0, index=genes)
    c_dep_fraction = _safe_fraction(c_dep_numerator, c_dp_n)
    r_dep_fraction = _safe_fraction(r_dep_numerator, r_dp_n)

    delta = c_ge_median - r_ge_median
    p_values = pd.Series(np.nan, index=genes, dtype=float)
    cliff_values = pd.Series(np.nan, index=genes, dtype=float)
    performed = pd.Series(False, index=genes, dtype=bool)

    if comparator_ids:
        c_values = c_ge.to_numpy(dtype=float)
        r_values = r_ge.to_numpy(dtype=float)
        for idx, gene in enumerate(genes):
            x = c_values[:, idx]
            y = r_values[:, idx]
            x = x[~np.isnan(x)]
            y = y[~np.isnan(y)]
            if len(x) >= min_test_n and len(y) >= min_test_n:
                test = mannwhitneyu(x, y, alternative="two-sided", method="auto")
                p_values.loc[gene] = float(test.pvalue)
                cliff_values.loc[gene] = cliffs_delta(x, y)
                performed.loc[gene] = True

    q_values = pd.Series(np.nan, index=genes, dtype=float)
    valid_p = p_values.dropna()
    if len(valid_p):
        adjusted = multipletests(valid_p.to_numpy(), alpha=0.05, method="fdr_bh")[1]
        q_values.loc[valid_p.index] = adjusted

    low_sample = (c_ge_n < min_context_n) | (
        (len(comparator_ids) > 0) & (r_ge_n < min_comparator_n)
    )

    gene_meta = _gene_metadata(files["gene_effect"]).set_index("gene_symbol")
    results = pd.DataFrame(
        {
            "gene_symbol": genes,
            "entrez_gene_id": [gene_meta.at[g, "entrez_gene_id"] if g in gene_meta.index else None for g in genes],
            "depmap_release": release,
            "cancer_id": cancer_id,
            "comparison_type": comparison,
            "context_definition": context_definition,
            "comparator_definition": comparator_definition,
            "context_models_n": len(context_ids),
            "comparator_models_n": len(comparator_ids),
            "context_gene_effect_n": c_ge_n.reindex(genes).to_numpy(),
            "comparator_gene_effect_n": r_ge_n.reindex(genes).to_numpy(),
            "context_dependency_probability_n": c_dp_n.reindex(genes).to_numpy(),
            "comparator_dependency_probability_n": r_dp_n.reindex(genes).to_numpy(),
            "context_median_gene_effect": c_ge_median.reindex(genes).to_numpy(),
            "context_mean_gene_effect": c_ge_mean.reindex(genes).to_numpy(),
            "context_gene_effect_iqr": (c_ge_q3 - c_ge_q1).reindex(genes).to_numpy(),
            "comparator_median_gene_effect": r_ge_median.reindex(genes).to_numpy(),
            "delta_gene_effect": delta.reindex(genes).to_numpy(),
            "cliffs_delta": cliff_values.reindex(genes).to_numpy(),
            "context_dependent_models_n": c_dep_numerator.reindex(genes).astype(int).to_numpy(),
            "context_dependency_fraction": c_dep_fraction.reindex(genes).to_numpy(),
            "comparator_dependent_models_n": r_dep_numerator.reindex(genes).astype(int).to_numpy(),
            "comparator_dependency_fraction": r_dep_fraction.reindex(genes).to_numpy(),
            "broad_dependency_n": broad_n.reindex(genes).to_numpy(),
            "broad_dependency_fraction": broad_fraction.reindex(genes).to_numpy(),
            "broad_dependency_warning": (
                broad_fraction.reindex(genes) >= broad_threshold
            ).fillna(False).to_numpy(),
            "statistical_test": ["Mann–Whitney U, two-sided" if x else None for x in performed.reindex(genes)],
            "p_value": p_values.reindex(genes).to_numpy(),
            "q_value": q_values.reindex(genes).to_numpy(),
            "fdr_0_05": (q_values.reindex(genes) < 0.05).fillna(False).to_numpy(),
            "low_sample_size": low_sample.reindex(genes).fillna(True).to_numpy(),
            "statistical_test_performed": performed.reindex(genes).to_numpy(),
        }
    )

    # Transparent ordering only: strongest negative context-vs-comparator median difference first.
    if comparator_ids:
        results = results.sort_values(
            ["delta_gene_effect", "q_value", "gene_symbol"],
            ascending=[True, True, True],
            na_position="last",
        ).reset_index(drop=True)
    else:
        results = results.sort_values(
            ["context_median_gene_effect", "gene_symbol"],
            ascending=[True, True],
            na_position="last",
        ).reset_index(drop=True)

    retrieved_at = datetime.now(timezone.utc).isoformat()
    meta = {
        "parser_version": GENOMEWIDE_PARSER_VERSION,
        "retrieved_at": retrieved_at,
        "depmap_release": release,
        "cancer_id": cancer_id,
        "comparison_type": comparison,
        "context_definition": context_definition,
        "comparator_definition": comparator_definition,
        "context_models_n": len(context_ids),
        "comparator_models_n": len(comparator_ids),
        "genes_gene_effect_n": len(effect_genes),
        "genes_dependency_n": len(dependency_genes),
        "genes_analyzed_n": len(genes),
        "dependency_probability_threshold": dep_threshold,
        "broad_dependency_fraction_warning": broad_threshold,
        "minimum_context_n_warning": min_context_n,
        "minimum_comparator_n_warning": min_comparator_n,
        "minimum_n_for_statistical_test": min_test_n,
        "multiple_testing": "Benjamini–Hochberg across all calculable genes within this one context/comparison",
        "interpretation_guardrail": (
            "Genome-wide results are genetic dependency evidence. CRISPR knockout is not equivalent "
            "to pharmacological inhibition, and broad dependency is not a safety verdict."
        ),
        "cohort_schema": schema_meta,
        "files": {
            key: {
                "path": str(path),
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
            for key, path in files.items()
        },
    }
    return results, cohort, gene_effect, dependency, meta
