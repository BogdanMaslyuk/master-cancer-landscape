from __future__ import annotations

"""Exploratory CDK4/target laboratory selection with transparent evidence rules.

No prediction of uploaded-compound activity. No synthetic numerical score.
The *same cancer* requirement applies to candidate genetic comparators;
a different lineage never becomes a matched negative control.
"""

from collections import defaultdict
from typing import Any


def _dep(row: dict[str, Any]) -> float | None:
    p = row.get("dependency_probability")
    if p is None:
        return None
    try:
        number = float(p)
    except (ValueError, TypeError):
        return None
    if not 0 <= number <= 1:
        return None
    return number


def lab_genetic_triage(
    full_context: dict[str, Any],
    *,
    eligible_tumor_model_ids: set[str],
) -> dict[str, Any]:
    """Rank *lab-registered human tumor* models only within genetic evidence.

    The caller must supply human-tumor ModelIDs from the laboratory registry.
    For missing registry roles, do not quietly infer a tumor identity.
    """
    gene = str(full_context.get("target_gene") or "").upper()
    if not full_context.get("available"):
        return {
            "available": False,
            "target_gene": gene,
            "status": "gene_context_unavailable",
            "items": [],
            "cohorts": [],
        }
    lab = [
        row for row in full_context.get("items", [])
        if row.get("registered_in_laboratory")
    ]
    lab_ids = {str(row.get("model_id")) for row in lab}
    human_tumor = [
        row for row in lab
        if str(row.get("model_id")) in eligible_tumor_model_ids
    ]
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in human_tumor:
        cancer = row.get("cancer_id")
        # Unknown tumor class is not a valid same-cancer comparator.
        if cancer is not None and str(cancer).strip():
            groups[str(cancer)].append(row)

    cohorts: list[dict[str, Any]] = []
    paired_cancer_ids: set[str] = set()
    for cancer_id, rows in sorted(groups.items()):
        dep = [r for r in rows if _dep(r) is not None and _dep(r) > 0.5]
        indep = [r for r in rows if _dep(r) is not None and _dep(r) <= 0.5]
        if dep and indep:
            paired_cancer_ids.add(cancer_id)
        cohorts.append({
            "cancer_id": cancer_id,
            "cancer_name": next((r.get("cancer_name") for r in rows if r.get("cancer_name")), None),
            "laboratory_human_tumor_models_n": len(rows),
            "dependent_models_n": len(dep),
            "nondependent_measured_models_n": len(indep),
            "missing_dependency_models_n": len(rows) - len(dep) - len(indep),
            "same_cancer_genetic_comparison_possible": bool(dep and indep),
            "dependent_model_ids": sorted(str(r["model_id"]) for r in dep),
            "nondependent_model_ids": sorted(str(r["model_id"]) for r in indep),
            "note_ru": (
                "Общая опухолевая классификация и различие CRISPR не доказывают "
                "сопоставимость остальных характеристик линий и не устанавливают "
                "чувствительность к новой молекуле."
            ),
        })

    items: list[dict[str, Any]] = []
    for row in human_tumor:
        model_id = str(row["model_id"])
        p = _dep(row)
        cancer_id = str(row["cancer_id"]) if row.get("cancer_id") else None
        has_pair = cancer_id in paired_cancer_ids if cancer_id else False
        if p is None:
            role, tier = "missing_dependency", "unranked"
        elif p > 0.5:
            role = "genetically_dependent_hypothesis_model"
            tier = "genetic_pair_available" if has_pair else "comparator_missing"
        elif has_pair:
            role, tier = "same_cancer_genetic_comparator_candidate", "genetic_pair_available"
        else:
            role, tier = "genetically_nondependent_context_only", "dependent_model_missing"
        items.append({
            "target_gene": gene,
            "model_id": model_id,
            "lab_names": row.get("laboratory_lines") or [],
            "cell_line_name": row.get("cell_line_name"),
            "cancer_id": cancer_id,
            "cancer_name": row.get("cancer_name"),
            "dependency_probability": p,
            "gene_effect": row.get("gene_effect"),
            "depmap_release": row.get("depmap_release"),
            "genetic_role": role,
            "genetic_triage_group": tier,
            "same_cancer_comparison_possible": has_pair,
            "reference_drug_observations_n": row.get("reference_drug_observations_n"),
            "reference_drug_evidence_status": row.get("reference_drug_evidence_status"),
            "compound_specific_activity": "not_measured",
            "laboratory_ready": False,
            "laboratory_stock_status": "not_verified",
            "ranking_basis": "measured_gene_dependency_and_same_cancer_lab_registry",
            "interpretation_ru": (
                "Место в списке означает исследовательскую пригодность для сравнения "
                "CRISPR-показателей, а не вероятность противоопухолевой активности PYZ."
            ),
        })
    group_order = {
        "genetic_pair_available": 0,
        "comparator_missing": 1,
        "dependent_model_missing": 2,
        "unranked": 3,
    }
    role_order = {
        "genetically_dependent_hypothesis_model": 0,
        "same_cancer_genetic_comparator_candidate": 1,
        "genetically_nondependent_context_only": 2,
        "missing_dependency": 3,
    }
    items.sort(key=lambda r: (
        group_order.get(r["genetic_triage_group"], 9),
        str(r["cancer_id"] or ""),
        role_order.get(r["genetic_role"], 9),
        -(r["dependency_probability"] if r["dependency_probability"] is not None else -1)
        if r["genetic_role"] == "genetically_dependent_hypothesis_model"
        else (r["dependency_probability"] if r["dependency_probability"] is not None else 2),
        str(r["model_id"]),
    ))
    return {
        "available": True,
        "status": "exploratory_genetic_lab_triage",
        "target_gene": gene,
        "genetic_threshold": "P(dep) > 0.5 is a screening rule, not a drug sensitivity threshold",
        "laboratory_registered_model_ids_n": len(lab_ids),
        "laboratory_human_tumor_model_ids_n": len(human_tumor),
        "lab_roles_not_classified_as_human_tumor_n": len(lab_ids) - len(human_tumor),
        "paired_cancer_cohorts_n": len(paired_cancer_ids),
        "cohorts": cohorts,
        "items": items,
        "limitations_ru": (
            "Только зарегистрированные человеческие опухолевые линии с DepMap-профилем "
            "данной мишени. Возможность сравнения по одному типу опухоли не означает "
            "полную сопоставимость. Принадлежность к лабораторному реестру не подтверждает "
            "фактическое наличие или чистоту культуры. Данные PRISM отдельных известных "
            "препаратов отображаются для справки и не влияют на ранжирование PYZ."
        ),
    }
