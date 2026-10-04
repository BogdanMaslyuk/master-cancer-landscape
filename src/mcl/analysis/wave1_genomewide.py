from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from mcl.analysis.depmap import cliffs_delta
from mcl.analysis.depmap_genomewide import _gene_metadata, _safe_fraction, _subset_rows
from mcl.sources.depmap import (
    load_genomewide_dependency_subset_with_background,
    load_genomewide_matrix_subset,
    resolve_depmap_files,
)


WAVE1_GENOMEWIDE_VERSION = "scientific-expansion-wave1-genomewide-v1"


def analyze_wave1_genomewide(
    root: Path,
    release: str,
    cohorts: pd.DataFrame,
    depmap_dir: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Run genome-wide CRISPR dependency comparisons for explicit Wave 1 cohorts.

    Multiple-testing correction is performed independently within each molecular
    context. No opaque global target score is calculated.
    """
    required = {
        "wave1_id",
        "label",
        "gene",
        "protein_change",
        "comparison_mode",
        "model_id",
        "analysis_group",
    }
    missing = required - set(cohorts.columns)
    if missing:
        raise ValueError(f"cohorts is missing columns: {sorted(missing)}")

    analysis_rows = cohorts[cohorts["analysis_group"].isin(["context", "comparator"])].copy()
    if analysis_rows.empty:
        raise ValueError("Wave 1 cohorts contain no context/comparator models")

    selected_ids = sorted(set(analysis_rows["model_id"].dropna().astype(str)))
    files = resolve_depmap_files(root, release, depmap_dir)

    thresholds = yaml.safe_load((root / "config" / "thresholds.yaml").read_text(encoding="utf-8")) or {}
    dcfg = thresholds.get("depmap", {}) or {}
    gcfg = dcfg.get("genomewide", {}) or {}
    dep_threshold = float(dcfg.get("dependency_probability_threshold", 0.5))
    broad_threshold = float(dcfg.get("broad_dependency_fraction_warning", 0.8))
    min_context_n = int(gcfg.get("minimum_context_n_warning", dcfg.get("minimum_context_n_warning", 5)))
    min_comparator_n = int(
        gcfg.get("minimum_comparator_n_warning", dcfg.get("minimum_comparator_n_warning", 5))
    )
    min_test_n = int(
        gcfg.get("minimum_n_for_statistical_test", dcfg.get("minimum_n_for_statistical_test", 2))
    )

    gene_effect = load_genomewide_matrix_subset(files["gene_effect"], selected_ids)
    dependency, broad_fraction, broad_n = load_genomewide_dependency_subset_with_background(
        files["gene_dependency"],
        selected_ids,
        dependency_threshold=dep_threshold,
    )

    dependency_genes = set(dependency.columns)
    genes = [gene for gene in gene_effect.columns if gene in dependency_genes]
    if not genes:
        raise ValueError("Gene Effect and Gene Dependency matrices have no common gene symbols")
    gene_effect = gene_effect[genes]
    dependency = dependency[genes]
    gene_meta = _gene_metadata(files["gene_effect"]).set_index("gene_symbol")

    all_results: list[pd.DataFrame] = []
    summary_rows: list[dict[str, object]] = []

    for wave1_id, group in cohorts.groupby("wave1_id", sort=True):
        context_ids = sorted(
            set(group.loc[group["analysis_group"].eq("context"), "model_id"].astype(str))
        )
        comparator_ids = sorted(
            set(group.loc[group["analysis_group"].eq("comparator"), "model_id"].astype(str))
        )
        if not context_ids or not comparator_ids:
            raise ValueError(
                f"{wave1_id}: explicit cohort requires non-empty context and comparator groups"
            )

        c_ge = _subset_rows(gene_effect, context_ids)
        r_ge = _subset_rows(gene_effect, comparator_ids)
        c_dp = _subset_rows(dependency, context_ids)
        r_dp = _subset_rows(dependency, comparator_ids)

        c_ge_n = c_ge.notna().sum(axis=0)
        r_ge_n = r_ge.notna().sum(axis=0)
        c_dp_n = c_dp.notna().sum(axis=0)
        r_dp_n = r_dp.notna().sum(axis=0)
        c_ge_median = c_ge.median(axis=0, skipna=True)
        r_ge_median = r_ge.median(axis=0, skipna=True)
        c_ge_mean = c_ge.mean(axis=0, skipna=True)
        c_ge_q1 = c_ge.quantile(0.25, axis=0, numeric_only=True)
        c_ge_q3 = c_ge.quantile(0.75, axis=0, numeric_only=True)
        c_dep_num = (c_dp > dep_threshold).sum(axis=0)
        r_dep_num = (r_dp > dep_threshold).sum(axis=0)
        c_dep_fraction = _safe_fraction(c_dep_num, c_dp_n)
        r_dep_fraction = _safe_fraction(r_dep_num, r_dp_n)
        delta = c_ge_median - r_ge_median

        p_values = pd.Series(np.nan, index=genes, dtype=float)
        cliff_values = pd.Series(np.nan, index=genes, dtype=float)
        performed = pd.Series(False, index=genes, dtype=bool)
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

        low_sample = (c_ge_n < min_context_n) | (r_ge_n < min_comparator_n)
        first = group.iloc[0]
        result = pd.DataFrame(
            {
                "wave1_id": wave1_id,
                "label": first["label"],
                "defining_gene": first["gene"],
                "defining_protein_change": first["protein_change"],
                "comparison_mode": first["comparison_mode"],
                "gene_symbol": genes,
                "entrez_gene_id": [
                    gene_meta.at[g, "entrez_gene_id"] if g in gene_meta.index else None for g in genes
                ],
                "depmap_release": release,
                "context_models_n": len(context_ids),
                "comparator_models_n": len(comparator_ids),
                "context_gene_effect_n": c_ge_n.reindex(genes).to_numpy(),
                "comparator_gene_effect_n": r_ge_n.reindex(genes).to_numpy(),
                "context_median_gene_effect": c_ge_median.reindex(genes).to_numpy(),
                "context_mean_gene_effect": c_ge_mean.reindex(genes).to_numpy(),
                "context_gene_effect_iqr": (c_ge_q3 - c_ge_q1).reindex(genes).to_numpy(),
                "comparator_median_gene_effect": r_ge_median.reindex(genes).to_numpy(),
                "delta_gene_effect": delta.reindex(genes).to_numpy(),
                "cliffs_delta": cliff_values.reindex(genes).to_numpy(),
                "context_dependency_probability_n": c_dp_n.reindex(genes).to_numpy(),
                "comparator_dependency_probability_n": r_dp_n.reindex(genes).to_numpy(),
                "context_dependent_models_n": c_dep_num.reindex(genes).astype(int).to_numpy(),
                "context_dependency_fraction": c_dep_fraction.reindex(genes).to_numpy(),
                "comparator_dependent_models_n": r_dep_num.reindex(genes).astype(int).to_numpy(),
                "comparator_dependency_fraction": r_dep_fraction.reindex(genes).to_numpy(),
                "broad_dependency_n": broad_n.reindex(genes).to_numpy(),
                "broad_dependency_fraction": broad_fraction.reindex(genes).to_numpy(),
                "broad_dependency_warning": (
                    broad_fraction.reindex(genes) >= broad_threshold
                ).fillna(False).to_numpy(),
                "statistical_test": [
                    "Mann–Whitney U, two-sided" if x else None for x in performed.reindex(genes)
                ],
                "p_value": p_values.reindex(genes).to_numpy(),
                "q_value": q_values.reindex(genes).to_numpy(),
                "fdr_0_05": (q_values.reindex(genes) < 0.05).fillna(False).to_numpy(),
                "low_sample_size": low_sample.reindex(genes).fillna(True).to_numpy(),
                "statistical_test_performed": performed.reindex(genes).to_numpy(),
            }
        )
        result = result.sort_values(
            ["delta_gene_effect", "q_value", "gene_symbol"],
            ascending=[True, True, True],
            na_position="last",
        ).reset_index(drop=True)
        all_results.append(result)

        summary_rows.append(
            {
                "wave1_id": wave1_id,
                "label": first["label"],
                "defining_gene": first["gene"],
                "defining_protein_change": first["protein_change"],
                "comparison_mode": first["comparison_mode"],
                "context_models_n": len(context_ids),
                "comparator_models_n": len(comparator_ids),
                "genes_tested_n": int(performed.sum()),
                "fdr_0_05_genes_n": int((q_values < 0.05).fillna(False).sum()),
                "negative_delta_genes_n": int((delta < 0).fillna(False).sum()),
                "low_sample_size_genes_n": int(low_sample.fillna(True).sum()),
            }
        )

    results = pd.concat(all_results, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    meta = {
        "contract": WAVE1_GENOMEWIDE_VERSION,
        "depmap_release": release,
        "contexts_n": int(summary["wave1_id"].nunique()),
        "genes_rows_n": int(len(results)),
        "unique_genes_n": int(results["gene_symbol"].nunique()),
        "dependency_probability_threshold": dep_threshold,
        "broad_dependency_fraction_warning": broad_threshold,
        "minimum_context_n_warning": min_context_n,
        "minimum_comparator_n_warning": min_comparator_n,
        "minimum_n_for_statistical_test": min_test_n,
        "multiple_testing": "Benjamini-Hochberg FDR independently within each Wave 1 context",
        "ranking_note": "Rows are ordered by delta_gene_effect; no composite target score is calculated.",
    }
    return results, summary, meta
