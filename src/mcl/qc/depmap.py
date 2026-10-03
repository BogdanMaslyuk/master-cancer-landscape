from __future__ import annotations

from collections import Counter

from mcl.models.depmap import (
    DepMapContextAuditRow,
    DepMapEvidence,
    DepMapKrasSensitivityAuditRow,
    DepMapKrasSensitivityEvidence,
)
from mcl.models.qc import QCRecord


def depmap_qc(
    rows: list[DepMapEvidence],
    audit: list[DepMapContextAuditRow],
    expected_pairs_n: int,
) -> list[QCRecord]:
    qc: list[QCRecord] = []

    if len(rows) != expected_pairs_n:
        qc.append(
            QCRecord(
                severity="ERROR",
                check="row_count",
                observed=str(len(rows)),
                expected=str(expected_pairs_n),
                message="DepMap evidence row count does not match Wave 1 Cancer×Target pairs.",
            )
        )

    pair_keys = [(x.cancer_id, x.target_id) for x in rows]
    dup_pairs = [k for k, n in Counter(pair_keys).items() if n > 1]
    if dup_pairs:
        qc.append(
            QCRecord(
                severity="ERROR",
                check="unique_cancer_target",
                observed=str(dup_pairs[:10]),
                expected="unique Cancer_ID × Target_ID",
                message="Duplicate Cancer×Target records in DepMap evidence.",
            )
        )

    for row in rows:
        entity = f"{row.cancer_id}|{row.target_id}"
        if row.context_models_n == 0:
            qc.append(
                QCRecord(
                    severity="ERROR",
                    check="context_nonempty",
                    entity_id=entity,
                    observed="0",
                    expected=">0",
                    message="No models were assigned to the molecular context.",
                )
            )
        if row.comparator_models_n == 0:
            qc.append(
                QCRecord(
                    severity="WARNING",
                    check="comparator_nonempty",
                    entity_id=entity,
                    observed="0",
                    expected=">0",
                    message="No same-disease comparator models are available; context-specific statistics cannot be interpreted.",
                )
            )
        if row.context_gene_effect_n == 0:
            qc.append(
                QCRecord(
                    severity="ERROR",
                    check="context_gene_effect_available",
                    entity_id=entity,
                    observed="0",
                    expected=">0",
                    message="Context exists but target Gene Effect is unavailable for all context models.",
                )
            )
        if row.context_dependency_probability_n == 0:
            qc.append(
                QCRecord(
                    severity="ERROR",
                    check="context_dependency_probability_available",
                    entity_id=entity,
                    observed="0",
                    expected=">0",
                    message="Context exists but target dependency probability is unavailable for all context models.",
                )
            )
        if row.low_sample_size:
            # A completely absent comparator is already reported by comparator_nonempty.
            # Do not additionally label a large descriptive context (e.g. GBM IDH-WT)
            # as a "low sample" solely because comparator_n == 0.
            qc.append(
                QCRecord(
                    severity="WARNING",
                    check="low_sample_size",
                    entity_id=entity,
                    observed=f"context={row.context_gene_effect_n}; comparator={row.comparator_gene_effect_n}",
                    expected="available groups meet configured minimums",
                    message="At least one available analysis group is below the pre-specified sample-size warning threshold; categorical interpretation must remain uncertain.",
                )
            )
        if row.dependency_fraction is not None and not (0 <= row.dependency_fraction <= 1):
            qc.append(
                QCRecord(
                    severity="ERROR",
                    check="dependency_fraction_range",
                    entity_id=entity,
                    observed=str(row.dependency_fraction),
                    expected="0..1",
                    message="Dependency fraction is outside the valid range.",
                )
            )
        if row.broad_dependency_fraction is not None and not (0 <= row.broad_dependency_fraction <= 1):
            qc.append(
                QCRecord(
                    severity="ERROR",
                    check="broad_dependency_fraction_range",
                    entity_id=entity,
                    observed=str(row.broad_dependency_fraction),
                    expected="0..1",
                    message="Broad dependency fraction is outside the valid range.",
                )
            )
        if row.statistical_test_performed and (row.p_value is None or row.q_value is None):
            qc.append(
                QCRecord(
                    severity="ERROR",
                    check="statistics_complete",
                    entity_id=entity,
                    observed=f"p={row.p_value}; q={row.q_value}",
                    expected="both p and q",
                    message="A performed statistical comparison must have both raw and adjusted significance values.",
                )
            )
        if not row.statistical_test_performed and row.p_value is not None:
            qc.append(
                QCRecord(
                    severity="ERROR",
                    check="statistics_consistency",
                    entity_id=entity,
                    observed=str(row.p_value),
                    expected="empty",
                    message="p-value exists although the statistical test is marked as not performed.",
                )
            )

    audit_keys = [(x.cancer_id, x.model_id) for x in audit]
    dup_audit = [k for k, n in Counter(audit_keys).items() if n > 1]
    if dup_audit:
        qc.append(
            QCRecord(
                severity="ERROR",
                check="context_audit_unique",
                observed=str(dup_audit[:10]),
                expected="one row per Cancer_ID × ModelID",
                message="Context audit contains duplicate model assignments.",
            )
        )

    groups_by_context: dict[str, Counter] = {}
    for x in audit:
        groups_by_context.setdefault(x.cancer_id, Counter())[x.assigned_group] += 1
    for cancer_id, groups in groups_by_context.items():
        if groups.get("context", 0) == 0:
            qc.append(
                QCRecord(
                    severity="ERROR",
                    check="context_audit_has_context",
                    entity_id=cancer_id,
                    observed=str(dict(groups)),
                    expected="context > 0",
                    message="Audit table contains no context models for this cancer context.",
                )
            )

    if not [x for x in qc if x.severity in {"ERROR", "WARNING"}]:
        qc.append(
            QCRecord(
                severity="INFO",
                check="depmap_analysis",
                message="All DepMap analysis checks passed without errors or warnings.",
            )
        )
    return qc


def depmap_kras_sensitivity_qc(
    rows: list[DepMapKrasSensitivityEvidence],
    audit: list[DepMapKrasSensitivityAuditRow],
    expected_rows_n: int,
) -> list[QCRecord]:
    """QC for M3.1 KRAS comparator sensitivity analysis."""
    qc: list[QCRecord] = []
    if len(rows) != expected_rows_n:
        qc.append(QCRecord(
            severity="ERROR",
            check="row_count",
            observed=str(len(rows)),
            expected=str(expected_rows_n),
            message="KRAS sensitivity row count does not match eligible pairs × comparator types.",
        ))

    keys = [(x.cancer_id, x.target_id, x.comparison_type) for x in rows]
    duplicates = [k for k, n in Counter(keys).items() if n > 1]
    if duplicates:
        qc.append(QCRecord(
            severity="ERROR",
            check="unique_sensitivity_comparison",
            observed=str(duplicates[:10]),
            expected="unique Cancer_ID × Target_ID × comparison_type",
            message="Duplicate KRAS sensitivity evidence rows detected.",
        ))

    for row in rows:
        entity = f"{row.cancer_id}|{row.target_id}|{row.comparison_type}"
        if row.context_models_n == 0 or row.context_gene_effect_n == 0:
            qc.append(QCRecord(
                severity="ERROR",
                check="context_available",
                entity_id=entity,
                observed=f"models={row.context_models_n}; gene_effect={row.context_gene_effect_n}",
                expected=">0",
                message="Exact KRAS-variant context is unavailable for this sensitivity comparison.",
            ))
        if row.comparator_models_n == 0:
            qc.append(QCRecord(
                severity="WARNING",
                check="comparator_nonempty",
                entity_id=entity,
                observed="0",
                expected=">0",
                message="This KRAS sensitivity comparator is absent in the release; no comparison can be inferred.",
            ))
        elif row.comparator_gene_effect_n == 0:
            qc.append(QCRecord(
                severity="WARNING",
                check="comparator_gene_effect_available",
                entity_id=entity,
                observed="0",
                expected=">0",
                message="Comparator models exist but target Gene Effect is unavailable.",
            ))
        if row.low_sample_size:
            qc.append(QCRecord(
                severity="WARNING",
                check="low_sample_size",
                entity_id=entity,
                observed=f"context={row.context_gene_effect_n}; comparator={row.comparator_gene_effect_n}",
                expected="available groups meet configured minimums",
                message="At least one available group is below the pre-specified n threshold; inference remains exploratory.",
            ))
        if row.statistical_test_performed and (row.p_value is None or row.q_value is None):
            qc.append(QCRecord(
                severity="ERROR",
                check="statistics_complete",
                entity_id=entity,
                observed=f"p={row.p_value}; q={row.q_value}",
                expected="both p and q",
                message="Performed M3.1 comparison must have raw and BH-adjusted values.",
            ))
        if not row.statistical_test_performed and row.p_value is not None:
            qc.append(QCRecord(
                severity="ERROR",
                check="statistics_consistency",
                entity_id=entity,
                observed=str(row.p_value),
                expected="empty",
                message="p-value exists although no M3.1 statistical test was performed.",
            ))

    audit_keys = [(x.cancer_id, x.model_id) for x in audit]
    duplicates = [k for k, n in Counter(audit_keys).items() if n > 1]
    if duplicates:
        qc.append(QCRecord(
            severity="ERROR",
            check="audit_unique",
            observed=str(duplicates[:10]),
            expected="one row per Cancer_ID × ModelID",
            message="KRAS sensitivity audit contains duplicate model assignments.",
        ))

    for cancer_id in sorted({x.cancer_id for x in audit}):
        subset = [x for x in audit if x.cancer_id == cancer_id]
        counts = Counter(x.assigned_group for x in subset)
        if counts.get("target_variant", 0) == 0:
            qc.append(QCRecord(
                severity="ERROR",
                check="audit_target_context_nonempty",
                entity_id=cancer_id,
                observed=str(dict(counts)),
                expected="target_variant > 0",
                message="No exact KRAS target-variant models in sensitivity audit.",
            ))
        if counts.get("kras_wildtype_proxy", 0) == 0:
            qc.append(QCRecord(
                severity="WARNING",
                check="audit_wildtype_comparator_nonempty",
                entity_id=cancer_id,
                observed=str(dict(counts)),
                expected="kras_wildtype_proxy > 0",
                message="No sequencing-backed KRAS WT proxy models are available.",
            ))
        if counts.get("other_kras_driver_hotspot", 0) == 0:
            qc.append(QCRecord(
                severity="WARNING",
                check="audit_other_kras_comparator_nonempty",
                entity_id=cancer_id,
                observed=str(dict(counts)),
                expected="other_kras_driver_hotspot > 0",
                message="No other KRAS driver/hotspot comparator models are available.",
            ))

    if not [x for x in qc if x.severity in {"ERROR", "WARNING"}]:
        qc.append(QCRecord(
            severity="INFO",
            check="depmap_kras_sensitivity",
            message="All M3.1 KRAS sensitivity checks passed without errors or warnings.",
        ))
    return qc
