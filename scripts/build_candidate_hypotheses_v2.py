from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

try:
    import duckdb
except ImportError as exc:
    raise SystemExit(
        'DuckDB is required. Install analysis extras: .\\.venv\\Scripts\\python.exe -m pip install -e ".[analysis]"'
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
PHARM = ROOT / "data" / "runtime" / "pharmacology"
PROCESSED = ROOT / "data" / "processed"
ATLAS = PROCESSED / "depmap_crispr_model_atlas.parquet"
TARGET_CONTEXT = PROCESSED / "depmap_model_target_context.parquet"
TARGET_REGISTRY = PROCESSED / "target_registry" / "protein_targets.parquet"
GENE_SUMMARY = ROOT / "data" / "runtime" / "explorer" / "gene_dependency_summary.parquet"
OUTPUT = PHARM / "candidate_hypotheses_v2.parquet"
MODEL_OUTPUT = PHARM / "candidate_hypothesis_models_v2.parquet"
MANIFEST = PHARM / "candidate_hypotheses_v2_manifest.json"
QC_DIR = ROOT / "outputs" / "qc"

PRISM_ACTIVE_LFC_THRESHOLD = -1.0
SELECTIVE_PERCENTILE = 0.75
NONRESPONSIVE_PERCENTILE = 0.25
DEPENDENCY_PROBABILITY_THRESHOLD = 0.5
MIN_CONTEXT_MODELS = 5

STATUS_LABELS = {
    "priority_for_in_vitro": "Приоритет для in vitro проверки",
    "supported_hypothesis": "Поддержанная гипотеза",
    "exploratory_hypothesis": "Исследовательская гипотеза",
    "insufficient_evidence": "Недостаточно данных",
}
STATUS_ORDER = {
    "priority_for_in_vitro": 0,
    "supported_hypothesis": 1,
    "exploratory_hypothesis": 2,
    "insufficient_evidence": 3,
}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _p(path: Path) -> str:
    return path.resolve().as_posix().replace("'", "''")


def _hypothesis_id(compound_id: str, target_gene: str, cancer_id: str) -> str:
    return "HYP2-" + hashlib.sha1(f"{compound_id}|{target_gene}|{cancer_id}".encode("utf-8")).hexdigest()[:14].upper()


def _float(value: Any) -> float | None:
    try:
        if value is None or pd.isna(value):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _fraction(n: Any, d: Any) -> float | None:
    try:
        denominator = int(d)
        return float(n) / denominator if denominator > 0 else None
    except (TypeError, ValueError):
        return None


def _mechanism_axis(label: Any) -> str:
    value = str(label or "")
    if value == "strong_support":
        return "strong"
    if value == "supportive":
        return "supportive"
    if value == "discordant":
        return "conflicting"
    if value in {"direction_not_resolved", "endpoint_not_resolved"}:
        return "not_assessable"
    if value in {"inconclusive", "insufficient"}:
        return "uncertain"
    return "not_available"


def _phenotype_axis(active_fraction: float | None, sensitive_fraction: float | None, active_n: int, sensitive_n: int, models_n: int) -> str:
    if models_n < 2:
        return "sparse"
    if active_n >= 2 and sensitive_n >= 2 and (active_fraction or 0) >= 0.30 and (sensitive_fraction or 0) >= 0.20:
        return "strong"
    if active_n >= 1 and (active_fraction or 0) >= 0.20:
        return "supportive"
    return "weak"


def _dependency_axis(dep_fraction: float | None, dependent_n: int, measured_n: int) -> str:
    if measured_n < 2:
        return "sparse"
    if dependent_n >= 2 and (dep_fraction or 0) >= 0.50:
        return "strong"
    if dependent_n >= 1 and (dep_fraction or 0) >= 0.25:
        return "supportive"
    return "weak"


def _molecular_context_axis(row: pd.Series) -> str:
    measured = int(row.get("context_models_n") or 0)
    if measured == 0:
        return "not_available"
    hotspot = int(row.get("hotspot_models_n") or 0)
    lof = int(row.get("likely_lof_models_n") or 0)
    expr = int(row.get("expression_detected_models_n") or 0)
    if hotspot or lof:
        return "mutation_context_present_direction_unresolved"
    if expr:
        return "expression_context_available"
    return "omics_context_available"


def _priority_status(row: pd.Series) -> tuple[str, list[str], list[str]]:
    models_n = int(row.get("models_n") or 0)
    active_n = int(row.get("active_models_n") or 0)
    sensitive_n = int(row.get("sensitive_models_n") or 0)
    dependency_measured_n = int(row.get("dependency_measured_models_n") or 0)
    dependency_n = int(row.get("dependency_models_n") or 0)
    joint_n = int(row.get("joint_support_models_n") or 0)
    active_fraction = _float(row.get("active_fraction")) or 0.0
    dependency_fraction = _float(row.get("dependency_fraction_in_cancer")) or 0.0
    mechanism = _mechanism_axis(row.get("concordance_label"))

    reasons: list[str] = []
    gaps: list[str] = [
        "PRISM primary screen — single-dose клеточный фенотип; он не заменяет полноценную dose-response валидацию.",
        "CRISPR knockout не эквивалентен фармакологическому ингибированию.",
        "Мутация, экспрессия и relative copy number пока не интерпретируются автоматически как причинный биомаркер.",
        "Риск для нормальных тканей, терапевтическое окно, прямое target engagement и ADMET пока не включены в статус.",
    ]
    if models_n < MIN_CONTEXT_MODELS:
        gaps.append("Малое число моделей в опухолевой группе; специфичность контекста считается исследовательской.")
    if dependency_measured_n == 0:
        gaps.append("Для моделей гипотезы отсутствует Probability of Dependency; высокий in vitro приоритет запрещён.")
    if mechanism == "conflicting":
        gaps.append("Панельная фармакология расходится с простой CRISPR loss-of-function моделью заявленной мишени.")
    if mechanism == "not_assessable":
        gaps.append("Направление действия вещества нельзя корректно сопоставить с CRISPR loss-of-function.")

    if active_n:
        reasons.append(f"{active_n} модель(и) имеют PRISM LFC ≤ {PRISM_ACTIVE_LFC_THRESHOLD:g}.")
    if sensitive_n:
        reasons.append(f"{sensitive_n} модель(и) одновременно проходят абсолютный activity threshold и входят в наиболее чувствительную четверть.")
    if dependency_n:
        reasons.append(f"{dependency_n} модель(и) имеют Probability of Dependency > {DEPENDENCY_PROBABILITY_THRESHOLD:g} для кодирующего гена мишени.")
    if joint_n:
        reasons.append(f"{joint_n} модель(и) одновременно фармакологически чувствительны и CRISPR-зависимы по Probability of Dependency.")
    if mechanism in {"strong", "supportive"}:
        reasons.append("Профиль чувствительности к веществу согласуется с непрерывным CRISPR Gene Effect заявленной мишени.")

    if (
        mechanism in {"strong", "supportive"}
        and models_n >= 3
        and dependency_measured_n >= 2
        and active_n >= 2
        and sensitive_n >= 2
        and dependency_n >= 2
        and joint_n >= 2
        and active_fraction >= 0.30
        and dependency_fraction >= 0.30
    ):
        return "priority_for_in_vitro", reasons, gaps

    if (
        mechanism != "conflicting"
        and dependency_measured_n >= 1
        and active_n >= 1
        and joint_n >= 1
    ):
        return "supported_hypothesis", reasons, gaps

    ge_depletion_n = int(row.get("gene_effect_depletion_models_n") or 0)
    if active_n >= 1 and (dependency_n >= 1 or ge_depletion_n >= 1):
        return "exploratory_hypothesis", reasons, gaps

    return "insufficient_evidence", reasons, gaps


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build Candidate Prioritization v2 with absolute PRISM activity, DepMap Probability of Dependency, "
            "continuous Gene Effect and model-level mutation/expression/relative-CN context."
        )
    )
    _ = parser.parse_args()

    required = {
        "responses": PHARM / "responses.parquet",
        "links": PHARM / "model_compound_target_links.parquet",
        "compounds": PHARM / "compound_catalog.parquet",
        "atlas": ATLAS,
        "context": TARGET_CONTEXT,
    }
    missing = [str(path.relative_to(ROOT)) for path in required.values() if not path.exists()]
    if missing:
        raise SystemExit(
            "Missing Candidate v2 inputs:\n" + "\n".join(missing) +
            "\nRun scripts/build_depmap_multiomics.py and scripts/build_depmap_target_context.py first."
        )

    con = duckdb.connect(database=":memory:")
    con.execute("PRAGMA threads=4")
    con.execute("PRAGMA preserve_insertion_order=false")

    con.execute(
        f"""
        CREATE VIEW responses_raw AS
        SELECT CAST(model_id AS VARCHAR) AS model_id,
               CAST(compound_id AS VARCHAR) AS compound_id,
               CAST(value AS DOUBLE) AS response_value
        FROM read_parquet('{_p(required['responses'])}')
        WHERE upper(CAST(source AS VARCHAR)) = 'PRISM'
          AND upper(CAST(endpoint AS VARCHAR)) = 'LFC'
          AND value IS NOT NULL;

        CREATE VIEW response_by_model AS
        SELECT model_id, compound_id, avg(response_value) AS response_value
        FROM responses_raw
        GROUP BY model_id, compound_id;

        CREATE VIEW response_profile AS
        SELECT
            model_id,
            compound_id,
            response_value,
            1.0 - percent_rank() OVER (PARTITION BY compound_id ORDER BY response_value ASC) AS relative_sensitivity
        FROM response_by_model;

        CREATE VIEW response_flagged AS
        SELECT *,
            response_value <= {PRISM_ACTIVE_LFC_THRESHOLD} AS absolute_active,
            relative_sensitivity >= {SELECTIVE_PERCENTILE} AS relatively_sensitive,
            response_value <= {PRISM_ACTIVE_LFC_THRESHOLD}
              AND relative_sensitivity >= {SELECTIVE_PERCENTILE} AS priority_sensitive,
            response_value > {PRISM_ACTIVE_LFC_THRESHOLD}
              AND relative_sensitivity <= {NONRESPONSIVE_PERCENTILE} AS clearly_nonresponsive
        FROM response_profile;

        CREATE VIEW compound_global AS
        SELECT
            compound_id,
            count(*) AS global_response_models_n,
            avg(CASE WHEN absolute_active THEN 1.0 ELSE 0.0 END) AS global_active_fraction,
            avg(CASE WHEN priority_sensitive THEN 1.0 ELSE 0.0 END) AS global_sensitive_fraction
        FROM response_flagged
        GROUP BY compound_id;

        CREATE VIEW links_dedup AS
        SELECT
            CAST(model_id AS VARCHAR) AS model_id,
            CAST(compound_id AS VARCHAR) AS compound_id,
            upper(CAST(target_gene AS VARCHAR)) AS target_gene,
            min(CAST(gene_effect AS DOUBLE)) AS link_gene_effect,
            string_agg(DISTINCT CAST(action AS VARCHAR), ' | ') FILTER (
                WHERE action IS NOT NULL AND CAST(action AS VARCHAR) <> ''
            ) AS actions
        FROM read_parquet('{_p(required['links'])}')
        GROUP BY model_id, compound_id, upper(CAST(target_gene AS VARCHAR));

        CREATE VIEW context AS
        SELECT * FROM read_parquet('{_p(required['context'])}');

        CREATE VIEW atlas_small AS
        SELECT
            CAST(model_id AS VARCHAR) AS model_id,
            cell_line_name,
            mcl_cancer_id,
            mcl_cancer_name,
            mcl_organ_ru,
            mcl_system_ru,
            oncotree_lineage,
            oncotree_subtype
        FROM read_parquet('{_p(required['atlas'])}');

        CREATE VIEW candidate_joined AS
        SELECT
            l.model_id,
            l.compound_id,
            l.target_gene,
            coalesce(CAST(c.gene_effect AS DOUBLE), l.link_gene_effect) AS gene_effect,
            CAST(c.dependency_probability AS DOUBLE) AS dependency_probability,
            CASE
                WHEN c.dependency_probability IS NULL THEN NULL
                ELSE CAST(c.dependency_probability AS DOUBLE) > {DEPENDENCY_PROBABILITY_THRESHOLD}
            END AS dependency_call,
            CAST(c.expression_log2_tpm1 AS DOUBLE) AS expression_log2_tpm1,
            CAST(c.expression_percentile_in_cancer AS DOUBLE) AS expression_percentile_in_cancer,
            CAST(c.copy_number_relative AS DOUBLE) AS copy_number_relative,
            CAST(c.copy_number_percentile_in_cancer AS DOUBLE) AS copy_number_percentile_in_cancer,
            coalesce(CAST(c.has_hotspot AS BOOLEAN), false) AS has_hotspot,
            coalesce(CAST(c.has_likely_lof AS BOOLEAN), false) AS has_likely_lof,
            c.mutation_resolution,
            c.protein_changes_json,
            l.actions,
            r.response_value,
            r.relative_sensitivity,
            r.absolute_active,
            r.relatively_sensitive,
            r.priority_sensitive,
            r.clearly_nonresponsive,
            a.cell_line_name,
            a.mcl_cancer_id,
            a.mcl_cancer_name,
            a.mcl_organ_ru,
            a.mcl_system_ru,
            a.oncotree_lineage,
            a.oncotree_subtype
        FROM links_dedup l
        INNER JOIN response_flagged r USING (model_id, compound_id)
        LEFT JOIN context c
          ON l.model_id = c.model_id AND l.target_gene = upper(c.target_gene)
        LEFT JOIN atlas_small a USING (model_id)
        WHERE a.mcl_cancer_id IS NOT NULL AND a.mcl_cancer_id <> '';
        """
    )

    grouped = con.execute(
        f"""
        SELECT
            j.compound_id,
            j.target_gene,
            j.mcl_cancer_id,
            any_value(j.mcl_cancer_name) AS mcl_cancer_name,
            any_value(j.mcl_organ_ru) AS mcl_organ_ru,
            any_value(j.mcl_system_ru) AS mcl_system_ru,
            any_value(j.actions) AS actions,
            count(DISTINCT j.model_id) AS models_n,
            count(DISTINCT j.model_id) FILTER (WHERE j.absolute_active) AS active_models_n,
            count(DISTINCT j.model_id) FILTER (WHERE j.relatively_sensitive) AS relatively_sensitive_models_n,
            count(DISTINCT j.model_id) FILTER (WHERE j.priority_sensitive) AS sensitive_models_n,
            count(DISTINCT j.model_id) FILTER (WHERE j.dependency_probability IS NOT NULL) AS dependency_measured_models_n,
            count(DISTINCT j.model_id) FILTER (WHERE j.dependency_call) AS dependency_models_n,
            count(DISTINCT j.model_id) FILTER (WHERE j.gene_effect IS NOT NULL) AS gene_effect_models_n,
            count(DISTINCT j.model_id) FILTER (WHERE j.gene_effect <= -0.5) AS gene_effect_depletion_models_n,
            count(DISTINCT j.model_id) FILTER (WHERE j.gene_effect <= -1.0) AS strong_gene_effect_models_n,
            count(DISTINCT j.model_id) FILTER (WHERE j.priority_sensitive AND j.dependency_call) AS joint_support_models_n,
            count(DISTINCT j.model_id) FILTER (
                WHERE j.priority_sensitive AND j.dependency_probability IS NOT NULL AND NOT j.dependency_call
            ) AS sensitive_without_dependency_n,
            count(DISTINCT j.model_id) FILTER (
                WHERE NOT j.absolute_active AND j.dependency_call
            ) AS dependency_without_activity_n,
            count(DISTINCT j.model_id) FILTER (WHERE j.has_hotspot) AS hotspot_models_n,
            count(DISTINCT j.model_id) FILTER (WHERE j.has_likely_lof) AS likely_lof_models_n,
            count(DISTINCT j.model_id) FILTER (WHERE j.expression_log2_tpm1 IS NOT NULL) AS context_models_n,
            count(DISTINCT j.model_id) FILTER (WHERE j.expression_log2_tpm1 > 0) AS expression_detected_models_n,
            median(j.response_value) AS median_response,
            min(j.response_value) AS best_response,
            median(j.relative_sensitivity) AS median_relative_sensitivity,
            median(j.gene_effect) FILTER (WHERE j.gene_effect IS NOT NULL) AS median_gene_effect,
            median(j.dependency_probability) FILTER (WHERE j.dependency_probability IS NOT NULL) AS median_dependency_probability,
            median(j.expression_log2_tpm1) FILTER (WHERE j.expression_log2_tpm1 IS NOT NULL) AS median_expression_log2_tpm1,
            median(j.expression_percentile_in_cancer) FILTER (WHERE j.expression_percentile_in_cancer IS NOT NULL) AS median_expression_percentile_in_cancer,
            median(j.copy_number_relative) FILTER (WHERE j.copy_number_relative IS NOT NULL) AS median_copy_number_relative,
            median(j.copy_number_percentile_in_cancer) FILTER (WHERE j.copy_number_percentile_in_cancer IS NOT NULL) AS median_copy_number_percentile_in_cancer,
            g.global_response_models_n,
            g.global_active_fraction,
            g.global_sensitive_fraction
        FROM candidate_joined j
        LEFT JOIN compound_global g USING (compound_id)
        GROUP BY j.compound_id, j.target_gene, j.mcl_cancer_id,
                 g.global_response_models_n, g.global_active_fraction, g.global_sensitive_fraction
        """
    ).df()

    compounds = pd.read_parquet(required["compounds"])
    keep = [c for c in ("compound_id", "preferred_name", "canonical_smiles", "inchikey", "chembl_id", "broad_id") if c in compounds.columns]
    grouped = grouped.merge(compounds[keep].drop_duplicates("compound_id"), on="compound_id", how="left")

    concordance_path = PHARM / "target_concordance.parquet"
    if concordance_path.exists():
        concordance = pd.read_parquet(concordance_path)
        concordance["target_gene"] = concordance["target_gene"].astype(str).str.upper()
        order = {
            "strong_support": 0, "supportive": 1, "inconclusive": 2, "discordant": 3,
            "direction_not_resolved": 4, "endpoint_not_resolved": 5, "insufficient": 6,
        }
        concordance["_order"] = concordance["concordance_label"].map(order).fillna(99)
        concordance = concordance.sort_values(
            ["compound_id", "target_gene", "_order", "q_value", "spearman_rho"],
            ascending=[True, True, True, True, False], na_position="last"
        ).drop_duplicates(["compound_id", "target_gene"], keep="first")
        ckeep = [c for c in (
            "compound_id", "target_gene", "concordance_label", "spearman_rho", "q_value", "models_n"
        ) if c in concordance.columns]
        concordance = concordance[ckeep].rename(columns={"models_n": "concordance_models_n"})
        grouped = grouped.merge(concordance, on=["compound_id", "target_gene"], how="left")

    if GENE_SUMMARY.exists():
        gs = pd.read_parquet(GENE_SUMMARY)
        gkeep = [c for c in (
            "gene_symbol", "dependency_fraction", "dependency_call_method", "dependency_type",
            "dependency_type_ru", "specificity_label_ru"
        ) if c in gs.columns]
        gs = gs[gkeep].rename(columns={
            "gene_symbol": "target_gene",
            "dependency_fraction": "global_dependency_fraction",
            "specificity_label_ru": "gene_global_specificity_label_ru",
        })
        gs["target_gene"] = gs["target_gene"].astype(str).str.upper()
        grouped = grouped.merge(gs, on="target_gene", how="left")
    else:
        grouped["global_dependency_fraction"] = pd.NA

    if TARGET_REGISTRY.exists():
        reg = pd.read_parquet(TARGET_REGISTRY)
        rkeep = [c for c in (
            "target_gene", "protein_preferred_name", "uniprot_primary_accession",
            "protein_mapping_status", "protein_families", "source_resolution"
        ) if c in reg.columns]
        reg = reg[rkeep].drop_duplicates("target_gene")
        reg["target_gene"] = reg["target_gene"].astype(str).str.upper()
        grouped = grouped.merge(reg, on="target_gene", how="left")

    grouped["active_fraction"] = [_fraction(a, b) for a, b in zip(grouped["active_models_n"], grouped["models_n"])]
    grouped["sensitive_fraction"] = [_fraction(a, b) for a, b in zip(grouped["sensitive_models_n"], grouped["models_n"])]
    grouped["dependency_fraction_in_cancer"] = [
        _fraction(a, b) for a, b in zip(grouped["dependency_models_n"], grouped["dependency_measured_models_n"])
    ]
    grouped["joint_support_fraction"] = [_fraction(a, b) for a, b in zip(grouped["joint_support_models_n"], grouped["models_n"])]
    grouped["active_enrichment"] = [
        (_float(a) - _float(b)) if _float(a) is not None and _float(b) is not None else None
        for a, b in zip(grouped["active_fraction"], grouped["global_active_fraction"])
    ]
    grouped["sensitivity_enrichment"] = [
        (_float(a) - _float(b)) if _float(a) is not None and _float(b) is not None else None
        for a, b in zip(grouped["sensitive_fraction"], grouped["global_sensitive_fraction"])
    ]
    grouped["dependency_enrichment"] = [
        (_float(a) - _float(b)) if _float(a) is not None and _float(b) is not None else None
        for a, b in zip(grouped["dependency_fraction_in_cancer"], grouped.get("global_dependency_fraction", pd.Series(pd.NA, index=grouped.index)))
    ]

    grouped["phenotype_axis"] = [
        _phenotype_axis(_float(a), _float(s), int(an), int(sn), int(n))
        for a, s, an, sn, n in zip(
            grouped["active_fraction"], grouped["sensitive_fraction"], grouped["active_models_n"], grouped["sensitive_models_n"], grouped["models_n"]
        )
    ]
    grouped["dependency_axis"] = [
        _dependency_axis(_float(f), int(dn), int(mn))
        for f, dn, mn in zip(grouped["dependency_fraction_in_cancer"], grouped["dependency_models_n"], grouped["dependency_measured_models_n"])
    ]
    grouped["mechanism_axis"] = grouped.get("concordance_label", pd.Series("", index=grouped.index)).map(_mechanism_axis)
    grouped["molecular_context_axis"] = grouped.apply(_molecular_context_axis, axis=1)
    grouped["normal_tissue_axis"] = "not_assessed"

    statuses: list[str] = []
    reasons: list[str] = []
    gaps: list[str] = []
    for _, row in grouped.iterrows():
        status, why, missing_evidence = _priority_status(row)
        statuses.append(status)
        reasons.append(json.dumps(why, ensure_ascii=False))
        gaps.append(json.dumps(missing_evidence, ensure_ascii=False))
    grouped["priority_status"] = statuses
    grouped["priority_status_ru"] = grouped["priority_status"].map(STATUS_LABELS)
    grouped["priority_reasons_json"] = reasons
    grouped["evidence_gaps_json"] = gaps
    grouped["hypothesis_id"] = [
        _hypothesis_id(str(c), str(t), str(k))
        for c, t, k in zip(grouped["compound_id"], grouped["target_gene"], grouped["mcl_cancer_id"])
    ]
    grouped["built_at"] = _now()

    grouped["_status_order"] = grouped["priority_status"].map(STATUS_ORDER).fillna(99)
    grouped = grouped.sort_values(
        ["_status_order", "joint_support_models_n", "sensitive_models_n", "active_models_n", "dependency_models_n", "models_n"],
        ascending=[True, False, False, False, False, False], na_position="last"
    ).drop(columns="_status_order").reset_index(drop=True)

    supported = grouped[grouped["priority_status"].isin(
        ["priority_for_in_vitro", "supported_hypothesis", "exploratory_hypothesis"]
    )][["hypothesis_id", "compound_id", "target_gene", "mcl_cancer_id"]].copy()
    con.register("supported", supported)

    model_rows = con.execute(
        f"""
        WITH labelled AS (
            SELECT
                k.hypothesis_id,
                j.*,
                CASE
                    WHEN j.priority_sensitive AND j.dependency_call THEN 'positive'
                    WHEN j.clearly_nonresponsive AND j.dependency_probability IS NOT NULL AND NOT j.dependency_call THEN 'negative_same_cancer'
                    WHEN j.priority_sensitive AND j.dependency_probability IS NOT NULL AND NOT j.dependency_call THEN 'discordant_sensitive_without_dependency'
                    WHEN NOT j.absolute_active AND j.dependency_call THEN 'discordant_dependency_without_activity'
                    WHEN j.dependency_probability IS NULL THEN 'unclassified_missing_dependency_probability'
                    ELSE 'unclassified'
                END AS model_role
            FROM candidate_joined j
            INNER JOIN supported k
              ON j.compound_id = k.compound_id
             AND j.target_gene = k.target_gene
             AND j.mcl_cancer_id = k.mcl_cancer_id
        ), ranked AS (
            SELECT *, row_number() OVER (
                PARTITION BY hypothesis_id, model_role
                ORDER BY
                    CASE
                        WHEN model_role IN ('positive','discordant_sensitive_without_dependency') THEN response_value
                        WHEN model_role = 'negative_same_cancer' THEN -response_value
                        ELSE coalesce(gene_effect, 0.0)
                    END ASC NULLS LAST,
                    CASE WHEN model_role = 'positive' THEN -dependency_probability ELSE dependency_probability END ASC NULLS LAST,
                    model_id
            ) AS role_rank
            FROM labelled
            WHERE model_role <> 'unclassified'
        )
        SELECT * EXCLUDE(role_rank)
        FROM ranked
        WHERE role_rank <= 5
        ORDER BY hypothesis_id, model_role, model_id
        """
    ).df()

    if not model_rows.empty:
        role_reason = {
            "positive": "PRISM LFC проходит абсолютный activity threshold, модель входит в чувствительную четверть и Probability of Dependency > 0.5.",
            "negative_same_cancer": "Модель того же опухолевого контекста находится среди наименее чувствительных, не проходит activity threshold и имеет измеренную Probability of Dependency ≤ 0.5.",
            "discordant_sensitive_without_dependency": "Фармакологическая чувствительность выражена, но Probability of Dependency для заявленной мишени ≤ 0.5; возможен альтернативный механизм или полифармакология.",
            "discordant_dependency_without_activity": "Probability of Dependency > 0.5, но выраженного PRISM activity signal нет; простая модель target dependency → drug response не подтверждается.",
            "unclassified_missing_dependency_probability": "Probability of Dependency отсутствует; модель нельзя считать ни положительной, ни отрицательной по бинарной CRISPR-зависимости.",
        }
        model_rows["role_reason_ru"] = model_rows["model_role"].map(role_reason)

    PHARM.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    grouped.to_parquet(OUTPUT, index=False, compression="zstd")
    model_rows.to_parquet(MODEL_OUTPUT, index=False, compression="zstd")

    status_counts = grouped["priority_status"].value_counts().to_dict()
    role_counts = model_rows["model_role"].value_counts().to_dict() if not model_rows.empty else {}
    manifest = {
        "contract": "mcl-candidate-hypotheses-v2",
        "built_at": _now(),
        "hypotheses_n": int(len(grouped)),
        "status_counts": {str(k): int(v) for k, v in status_counts.items()},
        "model_role_counts": {str(k): int(v) for k, v in role_counts.items()},
        "hypotheses_with_model_recommendations_n": int(supported["hypothesis_id"].nunique()),
        "model_recommendation_rows_n": int(len(model_rows)),
        "prism_active_lfc_threshold": PRISM_ACTIVE_LFC_THRESHOLD,
        "selective_percentile_threshold": SELECTIVE_PERCENTILE,
        "dependency_probability_threshold": DEPENDENCY_PROBABILITY_THRESHOLD,
        "ranking_contract": (
            "No composite numeric score. High priority requires absolute PRISM activity plus relative selectivity, "
            "measured Probability of Dependency and non-conflicting pharmacology-CRISPR concordance. Molecular omics context is shown separately and does not automatically add priority."
        ),
        "copy_number_contract": "relative linear CN only; no automatic amplification/deletion call",
        "mutation_contract": "Hotspot and LikelyLoF are separate annotations; hotspot is not automatically treated as gain-of-function",
        "not_yet_in_status": [
            "curated direction-aware biomarker rules",
            "normal-tissue expression and human loss-of-function tolerance",
            "direct target engagement/binding validation",
            "dose-response confirmation",
            "ADMET/developability",
        ],
        "outputs": [str(OUTPUT.relative_to(ROOT)), str(MODEL_OUTPUT.relative_to(ROOT))],
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    pd.DataFrame([
        {"priority_status": status, "priority_status_ru": STATUS_LABELS.get(status, status), "hypotheses_n": status_counts.get(status, 0)}
        for status in STATUS_ORDER
    ]).to_csv(QC_DIR / "candidate_hypotheses_v2_status_counts.tsv", sep="\t", index=False)
    pd.DataFrame([
        {"model_role": role, "models_n": count} for role, count in sorted(role_counts.items())
    ]).to_csv(QC_DIR / "candidate_hypothesis_models_v2_roles.tsv", sep="\t", index=False)

    print("MCL Candidate Prioritization v2")
    print(f"Hypotheses: {len(grouped)}")
    for status in STATUS_ORDER:
        print(f"  {status}: {status_counts.get(status, 0)}")
    print(f"Hypotheses with model recommendations: {manifest['hypotheses_with_model_recommendations_n']}")
    print(f"Model recommendation rows: {len(model_rows)}")
    for role, count in sorted(role_counts.items()):
        print(f"  model_role::{role}: {count}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")
    print(f"Wrote {MODEL_OUTPUT.relative_to(ROOT)}")
    print(f"Wrote {MANIFEST.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
