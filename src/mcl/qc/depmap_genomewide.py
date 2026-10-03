from __future__ import annotations

import pandas as pd

from mcl.models.qc import QCRecord


def depmap_genomewide_qc(results: pd.DataFrame, meta: dict) -> list[QCRecord]:
    """Compact QC for one genome-wide context analysis.

    Genome-wide analysis can contain ~18k rows, so QC deliberately aggregates
    repeated limitations instead of emitting one warning per gene.
    """
    qc: list[QCRecord] = []
    entity = f"{meta.get('cancer_id','')}|{meta.get('comparison_type','')}"

    if results.empty:
        qc.append(QCRecord(
            severity="ERROR", check="genomewide_nonempty", entity_id=entity,
            observed="0", expected=">0 genes",
            message="Genome-wide analysis produced no gene rows.",
        ))
        return qc

    if results["gene_symbol"].duplicated().any():
        duplicates = results.loc[results["gene_symbol"].duplicated(keep=False), "gene_symbol"].unique()
        qc.append(QCRecord(
            severity="ERROR", check="gene_symbol_unique", entity_id=entity,
            observed=str(list(duplicates[:20])), expected="unique gene symbols",
            message="Genome-wide result contains duplicate normalized gene symbols.",
        ))

    context_n = int(meta.get("context_models_n", 0))
    comparator_n = int(meta.get("comparator_models_n", 0))
    if context_n <= 0:
        qc.append(QCRecord(
            severity="ERROR", check="context_nonempty", entity_id=entity,
            observed=str(context_n), expected=">0",
            message="No context models were available for genome-wide analysis.",
        ))
    if comparator_n == 0:
        qc.append(QCRecord(
            severity="WARNING", check="comparator_nonempty", entity_id=entity,
            observed="0", expected=">0 for context-selective inference",
            message="No comparator cohort is available; results are descriptive dependencies only.",
        ))

    for col in ("context_dependency_fraction", "comparator_dependency_fraction", "broad_dependency_fraction"):
        if col in results:
            values = pd.to_numeric(results[col], errors="coerce").dropna()
            bad = values[(values < 0) | (values > 1)]
            if len(bad):
                qc.append(QCRecord(
                    severity="ERROR", check="dependency_fraction_range", entity_id=entity,
                    field=col, observed=str(len(bad)), expected="all values in 0..1",
                    message=f"{col} contains values outside the valid probability/fraction range.",
                ))

    performed = results.get("statistical_test_performed", pd.Series(False, index=results.index)).fillna(False).astype(bool)
    p = pd.to_numeric(results.get("p_value", pd.Series(index=results.index, dtype=float)), errors="coerce")
    q = pd.to_numeric(results.get("q_value", pd.Series(index=results.index, dtype=float)), errors="coerce")
    incomplete = performed & (p.isna() | q.isna())
    if incomplete.any():
        qc.append(QCRecord(
            severity="ERROR", check="statistics_complete", entity_id=entity,
            observed=str(int(incomplete.sum())), expected="0",
            message="Genes marked as statistically tested must have both p and q values.",
        ))
    if comparator_n == 0 and (p.notna().any() or q.notna().any()):
        qc.append(QCRecord(
            severity="ERROR", check="descriptive_mode_statistics", entity_id=entity,
            observed="p/q values present", expected="no p/q without comparator",
            message="Context-selective statistics were produced despite an absent comparator cohort.",
        ))

    low = results.get("low_sample_size", pd.Series(False, index=results.index)).fillna(False).astype(bool)
    if low.any():
        qc.append(QCRecord(
            severity="WARNING", check="low_sample_size", entity_id=entity,
            observed=f"{int(low.sum())}/{len(results)} genes",
            expected="available groups meet configured minimums",
            message="Some or all genes have fewer available models than the pre-specified warning threshold; inference remains exploratory.",
        ))

    qc.append(QCRecord(
        severity="INFO", check="genes_analyzed", entity_id=entity,
        observed=str(len(results)), expected=">0",
        message="Number of normalized genes evaluated in this genome-wide comparison.",
    ))
    if "fdr_0_05" in results:
        qc.append(QCRecord(
            severity="INFO", check="fdr_hits", entity_id=entity,
            observed=str(int(results["fdr_0_05"].fillna(False).astype(bool).sum())),
            expected="descriptive count",
            message="Genes with Benjamini-Hochberg q < 0.05; this is not an overall target-priority score.",
        ))
    if "broad_dependency_warning" in results:
        qc.append(QCRecord(
            severity="INFO", check="broad_dependency_flags", entity_id=entity,
            observed=str(int(results["broad_dependency_warning"].fillna(False).astype(bool).sum())),
            expected="descriptive count",
            message="Genes crossing the configured broad-dependency screening threshold.",
        ))

    if not [x for x in qc if x.severity in {"ERROR", "WARNING"}]:
        qc.append(QCRecord(
            severity="INFO", check="depmap_genomewide_analysis", entity_id=entity,
            message="Genome-wide DepMap analysis passed QC without errors or warnings.",
        ))
    return qc
