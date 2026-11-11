"""Deterministic scientific presentation profile for own PYZ molecules.

This view-model describes the existing 13-target research space only. It never
infers PYZ potency, binding or cellular selectivity from reference ligands,
CRISPR, ADMET or structural similarity. BindingDB P0.3.6 is not ingested.
"""
from __future__ import annotations

from typing import Any

CUTOFF = 0.35
PROFILE_CONTRACT = "mcl-pyz-profile-v2-p0"
TARGET_SOURCE = "existing-mcl-13-core-target-ligand-space"


def _number(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        num = float(value)
    except (ValueError, TypeError):
        return None
    import math
    return num if math.isfinite(num) else None


def _statement(code: str, title: str, meaning: str, *, source: str) -> dict[str, str]:
    return {"rule_id": code, "title": title, "meaning": meaning, "source": source}


def make_profile(detail: dict[str, Any]) -> dict[str, Any]:
    """Build a typed additive API block from legacy detail, with no data writes."""
    c = detail["compound"]
    compound_id = str(c["own_compound_id"])
    evidence = detail.get("target_evidence") or []
    contexts = detail.get("relevant_contexts") or []
    admet = detail.get("admet_summary") or {}

    targets = []
    for raw in evidence:
        similarity = _number(raw.get("best_tanimoto_morgan_r2_2048"))
        target = {
            "target_gene": str(raw.get("target_gene") or ""),
            "uniprot_accession": raw.get("uniprot_accession"),
            "tanimoto": similarity,
            "passes_navigation_cutoff": similarity is not None and similarity >= CUTOFF,
            "ligands_compared_n": raw.get("ligands_compared_n"),
            "best_ligand_id": raw.get("best_ligand_id"),
            "reference_value_nm": _number(raw.get("best_lowest_reported_value_nm_across_endpoints")),
            "reference_endpoint": raw.get("best_lowest_value_endpoint_type"),
            "reference_measurement_owner": "known_ligand_not_pyz",
            "evidence_type": "computational_prediction",
            "scientific_status": "exploratory_not_directly_validated",
            "source": TARGET_SOURCE,
            "neighbors": [
                {
                    "ligand_id": n.get("ligand_id"),
                    "tanimoto": _number(n.get("tanimoto_morgan_r2_2048")),
                    "reference_endpoint": n.get("lowest_value_endpoint_type"),
                    "reference_value_nm": _number(n.get("lowest_reported_value_nm_across_endpoints")),
                }
                for n in raw.get("top_neighbors", [])[:5]
            ],
        }
        targets.append(target)
    targets.sort(key=lambda x: (x["tanimoto"] is None, -(x["tanimoto"] or 0), x["target_gene"]))
    shortlist = [t for t in targets if t["passes_navigation_cutoff"]]
    admet_source = c.get("source_report") or admet.get("source_report")
    admet_available = any(c.get(key) is not None for key in ("qed", "integrated_risk", "admet_note_ru"))
    context_items = []
    for context in contexts:
        models = []
        for model in context.get("models") or []:
            models.append({
                "model_id": model.get("model_id"),
                "cell_line_name": model.get("cell_line_name"),
                "role": model.get("role"),
                "dependency_probability": _number(model.get("dependency_probability")),
                "gene_effect": _number(model.get("gene_effect")),
                "available_in_laboratory_registry": bool(model.get("available_in_laboratory")),
                "known_reference_compounds": list(model.get("known_reference_compounds") or []),
            })
        context_items.append({
            "target_gene": context.get("target_gene"),
            "mcl_organ_ru": context.get("mcl_organ_ru"),
            "mcl_cancer_name": context.get("mcl_cancer_name"),
            "mcl_cancer_id": context.get("mcl_cancer_id"),
            "tanimoto": _number(context.get("pyz_target_similarity")),
            "models": models,
            "context_status": "known_target_or_reference_context_not_pyz",
            "source": "Candidate v2 / DepMap / PRISM / laboratory registry",
        })

    established = [
        _statement("P0-IDENTITY", "Химическая структура зарегистрирована",
                   "Идентификатор и SMILES взяты из собственной библиотеки PYZ.",
                   source="config/own_compounds.tsv")
    ]
    hypotheses = []
    gaps = [
        _statement("P0-NO-PYZ-BINDING", "Прямое связывание с белком не подтверждено",
                   "Данные известных лигандов не являются измерениями для PYZ.",
                   source="MCL: отсутствует подтверждённый прямой опыт PYZ–мишень"),
        _statement("P0-NO-PYZ-CELL", "Жизнеспособность и селективность PYZ не установлены",
                   "Чувствительность клеток к другим веществам нельзя переносить на PYZ.",
                   source="MCL: собственные верифицированные клеточные измерения отсутствуют"),
    ]
    if shortlist:
        main = shortlist[0]
        hypotheses.append(_statement(
            "P0-SIMILARITY-SHORTLIST", f"Есть структурная гипотеза по {main['target_gene']}",
            f"Сходство по Танимото {main['tanimoto']:.3f} с известным лигандом. "
            "Это только основание для независимой проверки, не прогноз связывания.",
            source=TARGET_SOURCE))
        next_step = _statement(
            "P0-TEST-SHORTLIST", "Проверить наиболее обоснованную гипотезу",
            f"Сначала проверить экспериментальное происхождение аналога для {main['target_gene']}; "
            "затем рассмотреть прямой биохимический опыт и подходящие клеточные модели с контролями.",
            source="existing-mcl-target-evidence-and-model-context")
    else:
        gaps.append(_statement("P0-NO-STRUCTURAL-SUPPORT",
                   "Нет пары с Танимото ≥ 0,35 в действующей библиотеке",
                   "Отсутствие сходства не означает биологическую неактивность PYZ.",
                   source=TARGET_SOURCE))
        next_step = _statement(
            "P0-EXPLORE-INDEPENDENTLY", "Проверить независимые основания для гипотез",
            "Расширить поиск экспериментальных аналогов либо использовать отдельную "
            "методику поиска мишени; не объявлять соединение неактивным.",
            source="current-13-target-coverage-gap")
    if not admet_available:
        gaps.append(_statement("P0-ADMET-MISSING", "Индивидуальный ADMET не загружен",
                               "Прогноз не подменяется средними значениями по серии.",
                               source="config/pyz_admet_snapshot.tsv"))
    if not contexts:
        gaps.append(_statement("P0-NO-CONTEXT", "Сведения о клеточных моделях ограничены",
                               "Нет достаточного контекста Candidate v2 для данного маршрута.",
                               source="Candidate v2 / DepMap"))

    rules = [x["rule_id"] for x in [*established, *hypotheses, *gaps, next_step]]
    return {
        "schema_version": PROFILE_CONTRACT,
        "compound_id": compound_id,
        "identity": {
            "status": "available",
            "data_available": True,
            "own_compound_id": compound_id,
            "standardized_smiles": c.get("standardized_smiles"),
            "chemotype_ru": c.get("chemotype_ru"),
            "ring_variant": c.get("ring_variant"),
            "series_range": c.get("series_range"),
            "source": "config/own_compounds.tsv",
        },
        "properties": {
            "status": "available" if c.get("qed") is not None else "not_available",
            "data_available": c.get("qed") is not None,
            "qed": _number(c.get("qed")),
            "source": admet_source,
        },
        "admet": {
            "status": "available" if admet_available else "not_available",
            "data_available": admet_available,
            "integrated_risk": c.get("integrated_risk"),
            "note_ru": c.get("admet_note_ru"),
            "source_report": admet_source,
            "prediction_not_experiment": True,
            "individual_ames": None,
            "individual_dili": None,
        },
        "target_hypotheses": {
            "status": "available" if targets else "not_available",
            "data_available": bool(targets),
            "scope": "existing_13_core_targets",
            "cutoff_navigational": CUTOFF,
            "items": targets,
            "source": TARGET_SOURCE,
        },
        "cell_context": {
            "status": "available" if context_items else "not_available",
            "data_available": bool(context_items),
            "items": context_items,
        },
        "scientific_summary": {
            "established": established,
            "hypotheses": hypotheses,
            "gaps": gaps,
            "next_step": next_step,
            "rule_ids": rules,
            "own_cellular_activity": "not_measured",
            "own_target_binding": "not_measured",
        },
        "provenance": {
            "own_compounds_registry": "config/own_compounds.tsv",
            "admet_report": admet_source,
            "target_library": TARGET_SOURCE,
            "candidate_context": "Candidate v2 / DepMap / PRISM",
            "bindingdb_p036_status": "NOT_INGESTED_UNREVIEWED",
            "docking_status": detail.get("docking_status") or "pending_external_results",
        },
    }
