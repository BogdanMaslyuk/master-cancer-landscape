from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

try:
    import duckdb
except ImportError as exc:
    raise SystemExit(
        'DuckDB is required. Install analysis extras: .\\.venv\\Scripts\\python.exe -m pip install -e ".[analysis]"'
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
PHARM = ROOT / "data" / "runtime" / "pharmacology"
RUNTIME = ROOT / "data" / "runtime" / "laboratory"
QC_DIR = ROOT / "outputs" / "qc"

PANEL = PROCESSED / "laboratory_panel.parquet"
ATLAS = PROCESSED / "depmap_crispr_model_atlas.parquet"
CONTEXT = PROCESSED / "depmap_model_target_context.parquet"
HYPOTHESES = PHARM / "candidate_hypotheses_v2.parquet"
RESPONSES = PHARM / "responses.parquet"
LINKS = PHARM / "model_compound_target_links.parquet"
OUTPUT = RUNTIME / "laboratory_candidate_hypotheses.parquet"
MODEL_OUTPUT = RUNTIME / "laboratory_candidate_models.parquet"
MANIFEST = RUNTIME / "laboratory_candidates_manifest.json"
QC = QC_DIR / "laboratory_candidates_qc.tsv"

PRISM_ACTIVE_LFC_THRESHOLD = -1.0
SELECTIVE_PERCENTILE = 0.75
NONRESPONSIVE_PERCENTILE = 0.25
DEPENDENCY_PROBABILITY_THRESHOLD = 0.5

STATUS_ORDER = {
    "priority_for_in_vitro": 0,
    "supported_hypothesis": 1,
    "exploratory_hypothesis": 2,
    "insufficient_evidence": 3,
}
READINESS_ORDER = {
    "ready_with_internal_tumor_control": 0,
    "ready_positive_model": 1,
    "discordant_current_panel": 2,
    "no_support_in_current_panel": 3,
}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _p(path: Path) -> str:
    return path.resolve().as_posix().replace("'", "''")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Intersect Candidate v2 with the department's physically available human cell-line panel."
    )
    parser.parse_args()

    required = [PANEL, ATLAS, CONTEXT, HYPOTHESES, RESPONSES, LINKS]
    missing = [str(p.relative_to(ROOT)) for p in required if not p.exists()]
    if missing:
        raise SystemExit("Missing Laboratory Panel inputs:\n" + "\n".join(missing))

    panel = pd.read_parquet(PANEL)
    human_tumor = panel[
        panel["laboratory_role"].astype(str).eq("human_tumor")
        & panel["model_id"].notna()
        & panel.get("has_crispr_atlas", pd.Series(False, index=panel.index)).fillna(False).astype(bool)
    ].copy()
    if human_tumor.empty:
        raise SystemExit("No human tumor laboratory lines are mapped into the CRISPR Atlas. Review laboratory_panel.tsv first.")

    control = panel[panel["control_role"].astype(str).eq("general_non_tumor_control")].copy()
    control_name = str(control.iloc[0]["lab_name"]) if not control.empty else "BJ5ta"
    control_model_id = (
        str(control.iloc[0]["model_id"])
        if not control.empty and pd.notna(control.iloc[0].get("model_id")) and str(control.iloc[0].get("model_id")).strip()
        else None
    )

    con = duckdb.connect(database=":memory:")
    con.execute("PRAGMA threads=4")
    con.execute("PRAGMA preserve_insertion_order=false")
    con.register(
        "lab_models_df",
        human_tumor[["lab_id", "lab_name", "model_id", "highlighted_on_source", "priority_group"]].copy(),
    )
    con.execute("CREATE VIEW lab_models AS SELECT * FROM lab_models_df")

    con.execute(
        f"""
        CREATE VIEW response_by_model AS
        SELECT CAST(model_id AS VARCHAR) AS model_id,
               CAST(compound_id AS VARCHAR) AS compound_id,
               avg(CAST(value AS DOUBLE)) AS response_value
        FROM read_parquet('{_p(RESPONSES)}')
        WHERE upper(CAST(source AS VARCHAR)) = 'PRISM'
          AND upper(CAST(endpoint AS VARCHAR)) = 'LFC'
          AND value IS NOT NULL
        GROUP BY model_id, compound_id;

        CREATE VIEW response_profile AS
        SELECT *,
               1.0 - percent_rank() OVER (
                   PARTITION BY compound_id ORDER BY response_value ASC
               ) AS relative_sensitivity
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

        CREATE VIEW links_dedup AS
        SELECT CAST(model_id AS VARCHAR) AS model_id,
               CAST(compound_id AS VARCHAR) AS compound_id,
               upper(CAST(target_gene AS VARCHAR)) AS target_gene
        FROM read_parquet('{_p(LINKS)}')
        WHERE target_gene IS NOT NULL AND CAST(target_gene AS VARCHAR) <> ''
        GROUP BY model_id, compound_id, upper(CAST(target_gene AS VARCHAR));

        CREATE VIEW context AS
        SELECT CAST(model_id AS VARCHAR) AS model_id,
               upper(CAST(target_gene AS VARCHAR)) AS target_gene,
               CAST(gene_effect AS DOUBLE) AS gene_effect,
               CAST(dependency_probability AS DOUBLE) AS dependency_probability,
               CAST(expression_log2_tpm1 AS DOUBLE) AS expression_log2_tpm1,
               CAST(expression_percentile_in_cancer AS DOUBLE) AS expression_percentile_in_cancer,
               CAST(copy_number_relative AS DOUBLE) AS copy_number_relative,
               CAST(copy_number_percentile_in_cancer AS DOUBLE) AS copy_number_percentile_in_cancer,
               CAST(has_hotspot AS BOOLEAN) AS has_hotspot,
               CAST(has_likely_lof AS BOOLEAN) AS has_likely_lof,
               protein_changes_json
        FROM read_parquet('{_p(CONTEXT)}');

        CREATE VIEW atlas_small AS
        SELECT CAST(model_id AS VARCHAR) AS model_id,
               mcl_cancer_id, mcl_cancer_name, mcl_organ_ru,
               oncotree_lineage, oncotree_subtype
        FROM read_parquet('{_p(ATLAS)}');

        CREATE VIEW hypotheses AS
        SELECT * FROM read_parquet('{_p(HYPOTHESES)}')
        WHERE priority_status IN ('priority_for_in_vitro','supported_hypothesis','exploratory_hypothesis');

        CREATE VIEW lab_joined AS
        SELECT
            h.hypothesis_id,
            h.priority_status,
            h.priority_status_ru,
            h.preferred_name,
            h.compound_id,
            h.target_gene,
            h.protein_preferred_name,
            h.uniprot_primary_accession,
            h.mcl_cancer_id,
            h.mcl_cancer_name,
            h.mcl_organ_ru,
            h.mechanism_axis,
            h.concordance_label,
            lm.lab_id,
            lm.lab_name,
            lm.model_id,
            lm.highlighted_on_source,
            r.response_value,
            r.relative_sensitivity,
            r.absolute_active,
            r.priority_sensitive,
            r.clearly_nonresponsive,
            c.gene_effect,
            c.dependency_probability,
            CASE WHEN c.dependency_probability IS NULL THEN NULL
                 ELSE c.dependency_probability > {DEPENDENCY_PROBABILITY_THRESHOLD} END AS dependency_call,
            c.expression_log2_tpm1,
            c.expression_percentile_in_cancer,
            c.copy_number_relative,
            c.copy_number_percentile_in_cancer,
            coalesce(c.has_hotspot, false) AS has_hotspot,
            coalesce(c.has_likely_lof, false) AS has_likely_lof,
            c.protein_changes_json,
            a.oncotree_subtype,
            CASE
                WHEN r.priority_sensitive AND c.dependency_probability > {DEPENDENCY_PROBABILITY_THRESHOLD}
                    THEN 'positive'
                WHEN r.clearly_nonresponsive AND c.dependency_probability IS NOT NULL
                     AND c.dependency_probability <= {DEPENDENCY_PROBABILITY_THRESHOLD}
                    THEN 'negative_same_cancer'
                WHEN r.priority_sensitive AND c.dependency_probability IS NOT NULL
                     AND c.dependency_probability <= {DEPENDENCY_PROBABILITY_THRESHOLD}
                    THEN 'discordant_sensitive_without_dependency'
                WHEN NOT r.absolute_active AND c.dependency_probability > {DEPENDENCY_PROBABILITY_THRESHOLD}
                    THEN 'discordant_dependency_without_activity'
                WHEN c.dependency_probability IS NULL THEN 'unclassified_missing_dependency_probability'
                ELSE 'unclassified'
            END AS laboratory_model_role
        FROM hypotheses h
        INNER JOIN lab_models lm ON true
        INNER JOIN atlas_small a
          ON lm.model_id = a.model_id AND h.mcl_cancer_id = a.mcl_cancer_id
        INNER JOIN response_flagged r
          ON lm.model_id = r.model_id AND h.compound_id = r.compound_id
        INNER JOIN links_dedup l
          ON lm.model_id = l.model_id
         AND h.compound_id = l.compound_id
         AND upper(h.target_gene) = l.target_gene
        LEFT JOIN context c
          ON lm.model_id = c.model_id AND upper(h.target_gene) = c.target_gene;
        """
    )

    model_rows = con.execute(
        """
        SELECT * FROM lab_joined
        ORDER BY
          CASE laboratory_model_role
            WHEN 'positive' THEN 0
            WHEN 'negative_same_cancer' THEN 1
            WHEN 'discordant_sensitive_without_dependency' THEN 2
            WHEN 'discordant_dependency_without_activity' THEN 3
            ELSE 4 END,
          response_value ASC NULLS LAST,
          lab_name
        """
    ).df()

    aggregated = con.execute(
        """
        SELECT
          hypothesis_id,
          any_value(priority_status) AS priority_status,
          any_value(priority_status_ru) AS priority_status_ru,
          any_value(preferred_name) AS preferred_name,
          any_value(compound_id) AS compound_id,
          any_value(target_gene) AS target_gene,
          any_value(protein_preferred_name) AS protein_preferred_name,
          any_value(uniprot_primary_accession) AS uniprot_primary_accession,
          any_value(mcl_cancer_id) AS mcl_cancer_id,
          any_value(mcl_cancer_name) AS mcl_cancer_name,
          any_value(mcl_organ_ru) AS mcl_organ_ru,
          any_value(mechanism_axis) AS mechanism_axis,
          any_value(concordance_label) AS concordance_label,
          count(DISTINCT model_id) AS laboratory_models_with_data_n,
          count(DISTINCT model_id) FILTER (WHERE laboratory_model_role='positive') AS laboratory_positive_models_n,
          count(DISTINCT model_id) FILTER (WHERE laboratory_model_role='negative_same_cancer') AS laboratory_negative_models_n,
          count(DISTINCT model_id) FILTER (WHERE laboratory_model_role='discordant_sensitive_without_dependency') AS laboratory_sensitive_without_dependency_n,
          count(DISTINCT model_id) FILTER (WHERE laboratory_model_role='discordant_dependency_without_activity') AS laboratory_dependency_without_activity_n,
          string_agg(DISTINCT lab_name, ' | ') FILTER (WHERE laboratory_model_role='positive') AS positive_lab_lines,
          string_agg(DISTINCT lab_name, ' | ') FILTER (WHERE laboratory_model_role='negative_same_cancer') AS negative_lab_lines,
          string_agg(DISTINCT lab_name, ' | ') FILTER (WHERE laboratory_model_role LIKE 'discordant%') AS discordant_lab_lines
        FROM lab_joined
        GROUP BY hypothesis_id
        """
    ).df()

    hypotheses = pd.read_parquet(HYPOTHESES)
    # Keep every v2 priority candidate in the laboratory triage, including those with no currently matched evidence.
    priority = hypotheses[hypotheses["priority_status"].eq("priority_for_in_vitro")].copy()
    core_cols = [c for c in (
        "hypothesis_id", "priority_status", "priority_status_ru", "preferred_name", "compound_id",
        "target_gene", "protein_preferred_name", "uniprot_primary_accession", "mcl_cancer_id",
        "mcl_cancer_name", "mcl_organ_ru", "mechanism_axis", "concordance_label",
        "active_fraction", "dependency_fraction_in_cancer", "joint_support_models_n",
        "molecular_context_axis",
    ) if c in priority.columns]
    triage = priority[core_cols].drop_duplicates("hypothesis_id")
    if not aggregated.empty:
        agg_drop = [c for c in (
            "priority_status", "priority_status_ru", "preferred_name", "compound_id", "target_gene",
            "protein_preferred_name", "uniprot_primary_accession", "mcl_cancer_id", "mcl_cancer_name",
            "mcl_organ_ru", "mechanism_axis", "concordance_label"
        ) if c in aggregated.columns]
        triage = triage.merge(aggregated.drop(columns=agg_drop), on="hypothesis_id", how="left")

    count_cols = [
        "laboratory_models_with_data_n", "laboratory_positive_models_n", "laboratory_negative_models_n",
        "laboratory_sensitive_without_dependency_n", "laboratory_dependency_without_activity_n",
    ]
    for col in count_cols:
        if col not in triage.columns:
            triage[col] = 0
        triage[col] = pd.to_numeric(triage[col], errors="coerce").fillna(0).astype(int)

    def readiness(row: pd.Series) -> str:
        if int(row["laboratory_positive_models_n"]) >= 1 and int(row["laboratory_negative_models_n"]) >= 1:
            return "ready_with_internal_tumor_control"
        if int(row["laboratory_positive_models_n"]) >= 1:
            return "ready_positive_model"
        if int(row["laboratory_sensitive_without_dependency_n"]) + int(row["laboratory_dependency_without_activity_n"]) >= 1:
            return "discordant_current_panel"
        return "no_support_in_current_panel"

    triage["laboratory_readiness"] = triage.apply(readiness, axis=1)
    triage["bj5ta_available"] = not control.empty
    triage["bj5ta_model_id"] = control_model_id
    triage["bj5ta_role"] = "general_non_tumor_control"
    triage["bj5ta_interpretation_ru"] = (
        f"{control_name} доступна как общий человеческий неопухолевый контроль. Это не органоспецифическая нормальная ткань; "
        "сравнение с ней не следует называть доказанным терапевтическим окном."
    )

    if control_model_id:
        control_response = con.execute(
            "SELECT compound_id, response_value FROM response_by_model WHERE model_id = ?",
            [control_model_id],
        ).df()
        if not control_response.empty:
            control_response = control_response.rename(columns={"response_value": "bj5ta_prism_lfc"})
            triage = triage.merge(control_response, on="compound_id", how="left")
        else:
            triage["bj5ta_prism_lfc"] = pd.NA
    else:
        triage["bj5ta_prism_lfc"] = pd.NA
    triage["bj5ta_prism_available"] = pd.to_numeric(triage["bj5ta_prism_lfc"], errors="coerce").notna()

    triage["_readiness_order"] = triage["laboratory_readiness"].map(READINESS_ORDER).fillna(99)
    triage = triage.sort_values(
        ["_readiness_order", "laboratory_positive_models_n", "laboratory_negative_models_n", "joint_support_models_n"],
        ascending=[True, False, False, False], na_position="last"
    ).drop(columns="_readiness_order").reset_index(drop=True)

    RUNTIME.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    triage.to_parquet(OUTPUT, index=False, compression="zstd")
    model_rows.to_parquet(MODEL_OUTPUT, index=False, compression="zstd")

    readiness_counts = triage["laboratory_readiness"].value_counts().to_dict()
    manifest = {
        "contract": "mcl-laboratory-candidate-triage-v1",
        "built_at": _now(),
        "candidate_scope": "Candidate v2 priority_for_in_vitro only",
        "candidate_priorities_n": int(len(triage)),
        "mapped_human_tumor_lines_n": int(len(human_tumor)),
        "readiness_counts": {str(k): int(v) for k, v in readiness_counts.items()},
        "laboratory_model_rows_n": int(len(model_rows)),
        "bj5ta_available": bool(not control.empty),
        "bj5ta_model_id": control_model_id,
        "role_contract": {
            "positive": "PRISM LFC <= -1 AND top sensitivity quartile AND Probability of Dependency > 0.5",
            "negative_same_cancer": "same MCL cancer context AND PRISM LFC > -1 AND bottom sensitivity quartile AND measured Probability of Dependency <= 0.5",
        },
        "interpretation_ru": (
            "Laboratory readiness означает только то, что вычислительную гипотезу можно проверить на физически доступной панели кафедры. "
            "Она не повышает доказательность механизма и не означает клиническую перспективность вещества."
        ),
        "outputs": [str(OUTPUT.relative_to(ROOT)), str(MODEL_OUTPUT.relative_to(ROOT))],
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame([
        {"laboratory_readiness": key, "hypotheses_n": value}
        for key, value in readiness_counts.items()
    ]).to_csv(QC, sep="\t", index=False)

    print("MCL Laboratory Candidate Triage v1")
    print(f"Candidate v2 priorities: {len(triage)}")
    print(f"Mapped human tumor lines: {len(human_tumor)}")
    for key in READINESS_ORDER:
        print(f"  {key}: {readiness_counts.get(key, 0)}")
    print(f"Laboratory model evidence rows: {len(model_rows)}")
    print(f"BJ5ta available: {not control.empty}; DepMap ID: {control_model_id or 'not mapped'}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")
    print(f"Wrote {MODEL_OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
