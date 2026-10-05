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
except ImportError as exc:  # pragma: no cover - runtime dependency check
    raise SystemExit(
        'DuckDB is required for the memory-safe hypothesis build. Install project analysis extras: '
        '.\\.venv\\Scripts\\python.exe -m pip install -e ".[analysis]"'
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
PHARM = ROOT / "data" / "runtime" / "pharmacology"
ATLAS = ROOT / "data" / "processed" / "depmap_crispr_model_atlas.parquet"
GENE_SUMMARY = ROOT / "data" / "runtime" / "explorer" / "gene_dependency_summary.parquet"
TARGET_REGISTRY = ROOT / "data" / "processed" / "target_registry" / "protein_targets.parquet"
OUTPUT = PHARM / "candidate_hypotheses.parquet"
MODEL_OUTPUT = PHARM / "candidate_hypothesis_models.parquet"
MANIFEST = PHARM / "candidate_hypotheses_manifest.json"
QC_DIR = ROOT / "outputs" / "qc"

DEPENDENCY_THRESHOLD = -0.5
STRONG_DEPENDENCY_THRESHOLD = -1.0
SENSITIVE_PERCENTILE = 0.75

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


def _sql_path(path: Path) -> str:
    return path.resolve().as_posix().replace("'", "''")


def _hypothesis_id(compound_id: str, target_gene: str, cancer_id: str) -> str:
    payload = f"{compound_id}|{target_gene}|{cancer_id}".encode("utf-8")
    return "HYP-" + hashlib.sha1(payload).hexdigest()[:14].upper()


def _float(value: Any) -> float | None:
    try:
        if value is None or pd.isna(value):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _fraction(numerator: Any, denominator: Any) -> float | None:
    try:
        d = int(denominator)
        return float(numerator) / d if d > 0 else None
    except (TypeError, ValueError):
        return None


def _mechanism_axis(label: str) -> str:
    if label == "strong_support":
        return "strong"
    if label == "supportive":
        return "supportive"
    if label == "discordant":
        return "conflicting"
    if label in {"direction_not_resolved", "endpoint_not_resolved"}:
        return "not_assessable"
    if label in {"inconclusive", "insufficient"}:
        return "uncertain"
    return "not_available"


def _phenotype_axis(sensitive_fraction: float | None, sensitive_n: int, models_n: int) -> str:
    if models_n < 2:
        return "sparse"
    if sensitive_n >= 2 and sensitive_fraction is not None and sensitive_fraction >= 0.40:
        return "strong"
    if sensitive_n >= 1 and sensitive_fraction is not None and sensitive_fraction >= 0.25:
        return "supportive"
    return "weak"


def _dependency_axis(dependency_fraction: float | None, dependency_n: int, models_n: int) -> str:
    if models_n < 2:
        return "sparse"
    if dependency_n >= 2 and dependency_fraction is not None and dependency_fraction >= 0.50:
        return "strong"
    if dependency_n >= 1 and dependency_fraction is not None and dependency_fraction >= 0.25:
        return "supportive"
    return "weak"


def _specificity_axis(
    sensitivity_enrichment: float | None,
    dependency_enrichment: float | None,
    models_n: int,
) -> str:
    if models_n < 5:
        return "limited_sample"
    s = sensitivity_enrichment if sensitivity_enrichment is not None else 0.0
    d = dependency_enrichment if dependency_enrichment is not None else 0.0
    if s >= 0.15 and d >= 0.15:
        return "concordant_enrichment"
    if s >= 0.15 or d >= 0.15:
        return "partial_enrichment"
    if s <= 0 and d <= 0:
        return "not_enriched"
    return "weak_enrichment"


def _priority_status(row: pd.Series) -> tuple[str, list[str], list[str]]:
    models_n = int(row.get("models_n") or 0)
    sensitive_n = int(row.get("sensitive_models_n") or 0)
    dependency_n = int(row.get("dependency_models_n") or 0)
    joint_n = int(row.get("joint_support_models_n") or 0)
    sensitive_fraction = _float(row.get("sensitive_fraction")) or 0.0
    dependency_fraction = _float(row.get("dependency_fraction_in_cancer")) or 0.0
    sens_enrichment = _float(row.get("sensitivity_enrichment")) or 0.0
    dep_enrichment = _float(row.get("dependency_enrichment")) or 0.0
    concordance = str(row.get("concordance_label") or "")
    mechanism = _mechanism_axis(concordance)

    reasons: list[str] = []
    gaps: list[str] = [
        "Молекулярный контекст модели (мутации/амплификации/экспрессия) пока не включён в итоговый статус.",
        "Риск для нормальных тканей и терапевтическое окно пока не оценены.",
        "CRISPR-согласованность не заменяет прямую проверку связывания или target engagement.",
    ]
    if models_n < 5:
        gaps.append("Малое число моделей в опухолевой группе: вывод о специфичности считается исследовательским.")
    if mechanism == "not_assessable":
        gaps.append("Направление фармакологического действия нельзя корректно сопоставить с CRISPR loss-of-function.")
    if mechanism == "conflicting":
        gaps.append("Глобальный профиль чувствительности расходится с CRISPR-зависимостью заявленной мишени.")

    if sensitive_n:
        reasons.append(f"{sensitive_n} модель(и) входят в наиболее чувствительную четверть для этого вещества.")
    if dependency_n:
        reasons.append(f"{dependency_n} модель(и) имеют Gene Effect ≤ {DEPENDENCY_THRESHOLD} для кодирующего гена мишени.")
    if joint_n:
        reasons.append(f"{joint_n} модель(и) одновременно чувствительны к веществу и CRISPR-зависимы от мишени.")
    if concordance in {"strong_support", "supportive"}:
        reasons.append("Глобальный профиль вещество ↔ мишень согласуется с CRISPR-профилем по панели моделей.")
    if sens_enrichment >= 0.15:
        reasons.append("Чувствительность в этой опухолевой группе повышена относительно общей панели данного вещества.")
    if dep_enrichment >= 0.15:
        reasons.append("CRISPR-зависимость в этой опухолевой группе выше глобальной зависимости от гена.")

    if (
        mechanism in {"strong", "supportive"}
        and models_n >= 3
        and joint_n >= 2
        and sensitive_fraction >= 0.30
        and dependency_fraction >= 0.30
    ):
        return "priority_for_in_vitro", reasons, gaps

    if mechanism != "conflicting" and models_n >= 2 and joint_n >= 1 and (
        mechanism in {"strong", "supportive"}
        or sens_enrichment >= 0.10
        or dep_enrichment >= 0.10
    ):
        return "supported_hypothesis", reasons, gaps

    if joint_n >= 1 or (sensitive_n >= 1 and dependency_n >= 1):
        return "exploratory_hypothesis", reasons, gaps

    return "insufficient_evidence", reasons, gaps


def _register_optional_parquet(
    con: "duckdb.DuckDBPyConnection",
    view_name: str,
    path: Path,
    empty_columns: list[str],
) -> None:
    if path.exists():
        con.execute(f"CREATE OR REPLACE VIEW {view_name} AS SELECT * FROM read_parquet('{_sql_path(path)}')")
    else:
        con.register(f"_{view_name}_empty", pd.DataFrame(columns=empty_columns))
        con.execute(f"CREATE OR REPLACE VIEW {view_name} AS SELECT * FROM _{view_name}_empty")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build transparent compound × target × cancer hypotheses from pharmacology, CRISPR, "
            "Cancer Atlas and target-concordance layers. No single composite score is produced."
        )
    )
    parser.add_argument(
        "--include-insufficient",
        action="store_true",
        help="Keep hypotheses with insufficient evidence in the output (default: yes in full table, flag retained for compatibility).",
    )
    args = parser.parse_args()
    _ = args.include_insufficient

    required = {
        "responses": PHARM / "responses.parquet",
        "links": PHARM / "model_compound_target_links.parquet",
        "compounds": PHARM / "compound_catalog.parquet",
        "atlas": ATLAS,
    }
    missing = [str(path.relative_to(ROOT)) for path in required.values() if not path.exists()]
    if missing:
        raise SystemExit("Missing required inputs:\n" + "\n".join(missing))

    con = duckdb.connect(database=":memory:")
    con.execute("PRAGMA threads=4")
    con.execute("PRAGMA preserve_insertion_order=false")

    con.execute(
        f"""
        CREATE VIEW responses_raw AS
        SELECT model_id, compound_id, source, endpoint, CAST(value AS DOUBLE) AS response_value
        FROM read_parquet('{_sql_path(required['responses'])}')
        WHERE source = 'PRISM' AND upper(endpoint) = 'LFC' AND value IS NOT NULL;

        CREATE VIEW response_by_model AS
        SELECT model_id, compound_id, avg(response_value) AS response_value
        FROM responses_raw
        GROUP BY model_id, compound_id;

        CREATE VIEW response_profile AS
        SELECT
            model_id,
            compound_id,
            response_value,
            1.0 - percent_rank() OVER (
                PARTITION BY compound_id ORDER BY response_value ASC
            ) AS relative_sensitivity
        FROM response_by_model;

        CREATE VIEW response_profile_flagged AS
        SELECT *, relative_sensitivity >= {SENSITIVE_PERCENTILE} AS sensitive
        FROM response_profile;

        CREATE VIEW compound_global AS
        SELECT
            compound_id,
            count(*) AS global_response_models_n,
            avg(CASE WHEN sensitive THEN 1.0 ELSE 0.0 END) AS global_sensitive_fraction
        FROM response_profile_flagged
        GROUP BY compound_id;

        CREATE VIEW links_dedup AS
        SELECT
            model_id,
            compound_id,
            upper(target_gene) AS target_gene,
            min(CAST(gene_effect AS DOUBLE)) AS gene_effect,
            string_agg(DISTINCT CAST(action AS VARCHAR), ' | ') FILTER (
                WHERE action IS NOT NULL AND CAST(action AS VARCHAR) <> ''
            ) AS actions
        FROM read_parquet('{_sql_path(required['links'])}')
        GROUP BY model_id, compound_id, upper(target_gene);

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
        FROM read_parquet('{_sql_path(required['atlas'])}');

        CREATE VIEW candidate_joined AS
        SELECT
            l.model_id,
            l.compound_id,
            l.target_gene,
            l.gene_effect,
            l.actions,
            r.response_value,
            r.relative_sensitivity,
            r.sensitive,
            a.cell_line_name,
            a.mcl_cancer_id,
            a.mcl_cancer_name,
            a.mcl_organ_ru,
            a.mcl_system_ru,
            a.oncotree_lineage,
            a.oncotree_subtype
        FROM links_dedup l
        INNER JOIN response_profile_flagged r USING (model_id, compound_id)
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
            count(DISTINCT j.model_id) FILTER (WHERE j.sensitive) AS sensitive_models_n,
            count(DISTINCT j.model_id) FILTER (WHERE j.gene_effect IS NOT NULL) AS crispr_models_n,
            count(DISTINCT j.model_id) FILTER (WHERE j.gene_effect <= {DEPENDENCY_THRESHOLD}) AS dependency_models_n,
            count(DISTINCT j.model_id) FILTER (WHERE j.gene_effect <= {STRONG_DEPENDENCY_THRESHOLD}) AS strong_dependency_models_n,
            count(DISTINCT j.model_id) FILTER (
                WHERE j.sensitive AND j.gene_effect <= {DEPENDENCY_THRESHOLD}
            ) AS joint_support_models_n,
            count(DISTINCT j.model_id) FILTER (
                WHERE j.sensitive AND (j.gene_effect IS NULL OR j.gene_effect > {DEPENDENCY_THRESHOLD})
            ) AS sensitive_without_dependency_n,
            count(DISTINCT j.model_id) FILTER (
                WHERE NOT j.sensitive AND j.gene_effect <= {DEPENDENCY_THRESHOLD}
            ) AS dependency_without_sensitivity_n,
            median(j.response_value) AS median_response,
            min(j.response_value) AS best_response,
            median(j.relative_sensitivity) AS median_relative_sensitivity,
            median(j.gene_effect) FILTER (WHERE j.gene_effect IS NOT NULL) AS median_gene_effect,
            g.global_response_models_n,
            g.global_sensitive_fraction
        FROM candidate_joined j
        LEFT JOIN compound_global g USING (compound_id)
        GROUP BY
            j.compound_id, j.target_gene, j.mcl_cancer_id,
            g.global_response_models_n, g.global_sensitive_fraction
        """
    ).df()

    compounds = pd.read_parquet(required["compounds"])
    compound_keep = [c for c in (
        "compound_id", "preferred_name", "canonical_smiles", "inchikey", "chembl_id", "broad_id"
    ) if c in compounds.columns]
    grouped = grouped.merge(
        compounds[compound_keep].drop_duplicates("compound_id"),
        on="compound_id",
        how="left",
    )

    _register_optional_parquet(
        con,
        "concordance_input",
        PHARM / "target_concordance.parquet",
        ["compound_id", "target_gene", "concordance_label", "spearman_rho", "q_value", "models_n"],
    )
    concordance = con.execute(
        """
        SELECT * EXCLUDE (_rn) FROM (
            SELECT
                compound_id,
                upper(target_gene) AS target_gene,
                concordance_label,
                CAST(spearman_rho AS DOUBLE) AS spearman_rho,
                CAST(q_value AS DOUBLE) AS q_value,
                CAST(models_n AS BIGINT) AS concordance_models_n,
                row_number() OVER (
                    PARTITION BY compound_id, upper(target_gene)
                    ORDER BY
                        CASE concordance_label
                            WHEN 'strong_support' THEN 0
                            WHEN 'supportive' THEN 1
                            WHEN 'inconclusive' THEN 2
                            WHEN 'discordant' THEN 3
                            WHEN 'direction_not_resolved' THEN 4
                            WHEN 'endpoint_not_resolved' THEN 5
                            WHEN 'insufficient' THEN 6
                            ELSE 99
                        END,
                        q_value ASC NULLS LAST,
                        spearman_rho DESC NULLS LAST
                ) AS _rn
            FROM concordance_input
        ) WHERE _rn = 1
        """
    ).df()
    if not concordance.empty:
        grouped = grouped.merge(concordance, on=["compound_id", "target_gene"], how="left")

    if GENE_SUMMARY.exists():
        gene_summary = pd.read_parquet(GENE_SUMMARY)
        keep = [c for c in (
            "gene_symbol", "dependency_fraction", "dependency_type", "dependency_type_ru",
            "best_cancer_id", "best_cancer_name", "specificity_label_ru"
        ) if c in gene_summary.columns]
        gene_summary = gene_summary[keep].rename(columns={
            "gene_symbol": "target_gene",
            "dependency_fraction": "global_dependency_fraction",
            "specificity_label_ru": "gene_global_specificity_label_ru",
        })
        gene_summary["target_gene"] = gene_summary["target_gene"].astype(str).str.upper()
        grouped = grouped.merge(gene_summary, on="target_gene", how="left")
    else:
        grouped["global_dependency_fraction"] = pd.NA

    if TARGET_REGISTRY.exists():
        registry = pd.read_parquet(TARGET_REGISTRY)
        keep = [c for c in (
            "target_gene", "protein_preferred_name", "uniprot_primary_accession",
            "protein_mapping_status", "protein_families", "source_resolution"
        ) if c in registry.columns]
        registry = registry[keep].drop_duplicates("target_gene")
        registry["target_gene"] = registry["target_gene"].astype(str).str.upper()
        grouped = grouped.merge(registry, on="target_gene", how="left")

    grouped["sensitive_fraction"] = [
        _fraction(a, b) for a, b in zip(grouped["sensitive_models_n"], grouped["models_n"])
    ]
    grouped["dependency_fraction_in_cancer"] = [
        _fraction(a, b) for a, b in zip(grouped["dependency_models_n"], grouped["crispr_models_n"])
    ]
    grouped["joint_support_fraction"] = [
        _fraction(a, b) for a, b in zip(grouped["joint_support_models_n"], grouped["models_n"])
    ]
    grouped["sensitivity_enrichment"] = [
        (_float(a) - _float(b)) if _float(a) is not None and _float(b) is not None else None
        for a, b in zip(grouped["sensitive_fraction"], grouped["global_sensitive_fraction"])
    ]
    grouped["dependency_enrichment"] = [
        (_float(a) - _float(b)) if _float(a) is not None and _float(b) is not None else None
        for a, b in zip(grouped["dependency_fraction_in_cancer"], grouped["global_dependency_fraction"])
    ]

    grouped["phenotype_axis"] = [
        _phenotype_axis(_float(f), int(n), int(m))
        for f, n, m in zip(grouped["sensitive_fraction"], grouped["sensitive_models_n"], grouped["models_n"])
    ]
    grouped["dependency_axis"] = [
        _dependency_axis(_float(f), int(n), int(m))
        for f, n, m in zip(grouped["dependency_fraction_in_cancer"], grouped["dependency_models_n"], grouped["models_n"])
    ]
    grouped["mechanism_axis"] = grouped.get(
        "concordance_label", pd.Series("", index=grouped.index)
    ).fillna("").astype(str).map(_mechanism_axis)
    grouped["specificity_axis"] = [
        _specificity_axis(_float(s), _float(d), int(m))
        for s, d, m in zip(grouped["sensitivity_enrichment"], grouped["dependency_enrichment"], grouped["models_n"])
    ]
    grouped["molecular_context_axis"] = "not_assessed"
    grouped["normal_tissue_axis"] = "not_assessed"

    statuses: list[str] = []
    reasons_json: list[str] = []
    gaps_json: list[str] = []
    for _, row in grouped.iterrows():
        status, reasons, gaps = _priority_status(row)
        statuses.append(status)
        reasons_json.append(json.dumps(reasons, ensure_ascii=False))
        gaps_json.append(json.dumps(gaps, ensure_ascii=False))
    grouped["priority_status"] = statuses
    grouped["priority_status_ru"] = grouped["priority_status"].map(STATUS_LABELS)
    grouped["priority_reasons_json"] = reasons_json
    grouped["evidence_gaps_json"] = gaps_json
    grouped["hypothesis_id"] = [
        _hypothesis_id(str(c), str(t), str(k))
        for c, t, k in zip(grouped["compound_id"], grouped["target_gene"], grouped["mcl_cancer_id"])
    ]
    grouped["built_at"] = _now()

    grouped["_status_order"] = grouped["priority_status"].map(STATUS_ORDER).fillna(99)
    grouped = grouped.sort_values(
        ["_status_order", "joint_support_models_n", "sensitive_models_n", "dependency_models_n", "models_n"],
        ascending=[True, False, False, False, False],
        na_position="last",
    ).drop(columns="_status_order").reset_index(drop=True)

    # Select concrete experimental models only for hypotheses that have at least some support.
    supported_keys = grouped[
        grouped["priority_status"].isin(
            ["priority_for_in_vitro", "supported_hypothesis", "exploratory_hypothesis"]
        )
    ][["hypothesis_id", "compound_id", "target_gene", "mcl_cancer_id"]].copy()
    con.register("supported_keys", supported_keys)
    model_rows = con.execute(
        f"""
        WITH labelled AS (
            SELECT
                k.hypothesis_id,
                j.model_id,
                j.cell_line_name,
                j.compound_id,
                j.target_gene,
                j.mcl_cancer_id,
                j.mcl_cancer_name,
                j.mcl_organ_ru,
                j.oncotree_subtype,
                j.response_value,
                j.relative_sensitivity,
                j.gene_effect,
                CASE
                    WHEN j.sensitive AND j.gene_effect <= {DEPENDENCY_THRESHOLD} THEN 'positive'
                    WHEN NOT j.sensitive AND (j.gene_effect IS NULL OR j.gene_effect > {DEPENDENCY_THRESHOLD}) THEN 'negative_same_cancer'
                    WHEN j.sensitive AND (j.gene_effect IS NULL OR j.gene_effect > {DEPENDENCY_THRESHOLD}) THEN 'discordant_sensitive_without_dependency'
                    WHEN NOT j.sensitive AND j.gene_effect <= {DEPENDENCY_THRESHOLD} THEN 'discordant_dependency_without_sensitivity'
                    ELSE 'unclassified'
                END AS model_role
            FROM candidate_joined j
            INNER JOIN supported_keys k
              ON j.compound_id = k.compound_id
             AND j.target_gene = k.target_gene
             AND j.mcl_cancer_id = k.mcl_cancer_id
        ), ranked AS (
            SELECT *,
                row_number() OVER (
                    PARTITION BY hypothesis_id, model_role
                    ORDER BY
                        CASE
                            WHEN model_role IN ('positive', 'discordant_sensitive_without_dependency')
                                THEN -relative_sensitivity
                            WHEN model_role = 'negative_same_cancer'
                                THEN relative_sensitivity
                            ELSE coalesce(gene_effect, 0.0)
                        END ASC,
                        CASE
                            WHEN model_role = 'positive' THEN gene_effect
                            WHEN model_role = 'negative_same_cancer' THEN -coalesce(gene_effect, 0.0)
                            ELSE response_value
                        END ASC NULLS LAST,
                        model_id
                ) AS role_rank
            FROM labelled
            WHERE model_role <> 'unclassified'
        )
        SELECT * EXCLUDE (role_rank)
        FROM ranked
        WHERE role_rank <= 5
        ORDER BY hypothesis_id, model_role, model_id
        """
    ).df()

    if not model_rows.empty:
        role_reason = {
            "positive": "Высокая относительная чувствительность к веществу сочетается с CRISPR-зависимостью от кодирующего гена мишени.",
            "negative_same_cancer": "Модель того же опухолевого контекста менее чувствительна и не показывает выраженной CRISPR-зависимости; кандидат на отрицательный контроль.",
            "discordant_sensitive_without_dependency": "Модель чувствительна к веществу без соответствующей CRISPR-зависимости; возможен альтернативный механизм или полифармакология.",
            "discordant_dependency_without_sensitivity": "Модель зависит от гена, но не входит в чувствительную четверть для вещества; это ограничивает простую механистическую интерпретацию.",
        }
        model_rows["role_reason_ru"] = model_rows["model_role"].map(role_reason)

    PHARM.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    grouped.to_parquet(OUTPUT, index=False, compression="zstd")
    model_rows.to_parquet(MODEL_OUTPUT, index=False, compression="zstd")

    status_counts = grouped["priority_status"].value_counts().to_dict()
    manifest = {
        "contract": "mcl-candidate-hypotheses-v1",
        "built_at": _now(),
        "hypotheses_n": int(len(grouped)),
        "hypotheses_with_model_recommendations_n": int(supported_keys["hypothesis_id"].nunique()),
        "model_recommendation_rows_n": int(len(model_rows)),
        "status_counts": {str(k): int(v) for k, v in status_counts.items()},
        "dependency_threshold": DEPENDENCY_THRESHOLD,
        "strong_dependency_threshold": STRONG_DEPENDENCY_THRESHOLD,
        "sensitive_percentile_threshold": SENSITIVE_PERCENTILE,
        "response_scope": "PRISM LFC only; lower LFC is treated as stronger depletion relative to control.",
        "ranking_contract": (
            "No composite numeric score. Priority status is a transparent rule over separate phenotype, CRISPR, "
            "profile-concordance, context-enrichment and sample-size axes."
        ),
        "not_yet_in_status": [
            "model-level mutation/amplification/expression context",
            "normal-tissue expression and human loss-of-function tolerance",
            "direct target engagement or binding validation",
            "ADMET/developability",
        ],
        "sources": {
            "responses": str(required["responses"].relative_to(ROOT)),
            "links": str(required["links"].relative_to(ROOT)),
            "atlas": str(required["atlas"].relative_to(ROOT)),
            "concordance": str((PHARM / "target_concordance.parquet").relative_to(ROOT)),
            "gene_dependency_summary": str(GENE_SUMMARY.relative_to(ROOT)),
            "protein_target_registry": str(TARGET_REGISTRY.relative_to(ROOT)),
        },
        "outputs": [str(OUTPUT.relative_to(ROOT)), str(MODEL_OUTPUT.relative_to(ROOT))],
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame(
        [{"priority_status": k, "priority_status_ru": STATUS_LABELS.get(k, k), "hypotheses_n": v}
         for k, v in status_counts.items()]
    ).to_csv(QC_DIR / "candidate_hypotheses_status_counts.tsv", sep="\t", index=False)

    print("MCL Candidate Prioritization v1")
    print(f"Hypotheses: {manifest['hypotheses_n']}")
    for status in STATUS_ORDER:
        print(f"  {status}: {status_counts.get(status, 0)}")
    print(f"Hypotheses with model recommendations: {manifest['hypotheses_with_model_recommendations_n']}")
    print(f"Model recommendation rows: {manifest['model_recommendation_rows_n']}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")
    print(f"Wrote {MODEL_OUTPUT.relative_to(ROOT)}")
    print(f"Wrote {MANIFEST.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
