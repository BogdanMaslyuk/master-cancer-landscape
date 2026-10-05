from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

try:
    import duckdb
except ImportError as exc:
    raise SystemExit(
        'DuckDB is required for the pharmacology audit. Install analysis extras: '
        '.\\.venv\\Scripts\\python.exe -m pip install -e ".[analysis]"'
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
PHARM = ROOT / "data" / "runtime" / "pharmacology"
QC = ROOT / "outputs" / "qc"
UNRESOLVED = QC / "pharmacology_unresolved_responses.tsv"
OUTPUT_JSON = QC / "pharmacology_layer_audit.json"
OUTPUT_TSV = QC / "pharmacology_layer_audit.tsv"


def _p(path: Path) -> str:
    return path.resolve().as_posix().replace("'", "''")


def main() -> None:
    responses = PHARM / "responses.parquet"
    compounds = PHARM / "compounds.parquet"
    targets = PHARM / "target_evidence.parquet"
    missing = [p for p in (responses, compounds, targets) if not p.exists()]
    if missing:
        raise SystemExit("Missing pharmacology runtime files: " + ", ".join(str(x) for x in missing))

    con = duckdb.connect(database=":memory:")
    response_summary = con.execute(
        f"""
        SELECT
            count(*) AS observations_n,
            count(DISTINCT CAST(model_id AS VARCHAR)) AS models_n,
            count(DISTINCT CAST(compound_id AS VARCHAR)) AS compounds_with_response_n,
            count(DISTINCT CAST(source AS VARCHAR)) AS sources_n
        FROM read_parquet('{_p(responses)}')
        """
    ).df().iloc[0].to_dict()

    assay_rows = con.execute(
        f"""
        SELECT
            source, source_release, assay_type, endpoint, unit,
            dose, dose_unit, exposure_time_h,
            count(*) AS observations_n,
            count(DISTINCT CAST(model_id AS VARCHAR)) AS models_n,
            count(DISTINCT CAST(compound_id AS VARCHAR)) AS compounds_n
        FROM read_parquet('{_p(responses)}')
        GROUP BY source, source_release, assay_type, endpoint, unit, dose, dose_unit, exposure_time_h
        ORDER BY observations_n DESC
        """
    ).df()

    # Some optional identifier/chemistry columns can be inferred by pandas/Parquet as
    # DOUBLE when every source value is missing. Never compare those columns directly
    # with string literals in DuckDB. Cast to VARCHAR first and normalize textual null
    # sentinels so the audit is schema-tolerant without changing the source tables.
    compound_stats = con.execute(
        f"""
        WITH normalized AS (
            SELECT
                trim(CAST(compound_id AS VARCHAR)) AS compound_id_text,
                CASE
                    WHEN broad_id IS NULL THEN NULL
                    WHEN lower(trim(CAST(broad_id AS VARCHAR))) IN ('', 'nan', 'none', 'null') THEN NULL
                    ELSE trim(CAST(broad_id AS VARCHAR))
                END AS broad_id_text,
                CASE
                    WHEN canonical_smiles IS NULL THEN NULL
                    WHEN lower(trim(CAST(canonical_smiles AS VARCHAR))) IN ('', 'nan', 'none', 'null') THEN NULL
                    ELSE trim(CAST(canonical_smiles AS VARCHAR))
                END AS smiles_text
            FROM read_parquet('{_p(compounds)}')
        )
        SELECT
            count(*) AS compound_rows_n,
            count(DISTINCT compound_id_text) AS compound_ids_n,
            count(DISTINCT broad_id_text) AS nonempty_broad_ids_n,
            count(DISTINCT smiles_text) AS nonempty_smiles_n,
            sum(CASE WHEN smiles_text IS NULL THEN 1 ELSE 0 END) AS missing_smiles_n
        FROM normalized
        """
    ).df().iloc[0].to_dict()

    duplicate_identity = con.execute(
        f"""
        WITH normalized AS (
            SELECT
                CASE
                    WHEN broad_id IS NULL THEN NULL
                    WHEN lower(trim(CAST(broad_id AS VARCHAR))) IN ('', 'nan', 'none', 'null') THEN NULL
                    ELSE trim(CAST(broad_id AS VARCHAR))
                END AS broad_id_text,
                CASE
                    WHEN canonical_smiles IS NULL THEN NULL
                    WHEN lower(trim(CAST(canonical_smiles AS VARCHAR))) IN ('', 'nan', 'none', 'null') THEN NULL
                    ELSE trim(CAST(canonical_smiles AS VARCHAR))
                END AS smiles_text
            FROM read_parquet('{_p(compounds)}')
        ), broad AS (
            SELECT broad_id_text, count(*) AS n
            FROM normalized
            WHERE broad_id_text IS NOT NULL
            GROUP BY broad_id_text HAVING count(*) > 1
        ), smiles AS (
            SELECT smiles_text, count(*) AS n
            FROM normalized
            WHERE smiles_text IS NOT NULL
            GROUP BY smiles_text HAVING count(*) > 1
        )
        SELECT
            (SELECT count(*) FROM broad) AS duplicated_broad_ids_n,
            (SELECT coalesce(sum(n),0) FROM broad) AS rows_with_duplicated_broad_id_n,
            (SELECT count(*) FROM smiles) AS duplicated_smiles_n,
            (SELECT coalesce(sum(n),0) FROM smiles) AS rows_with_duplicated_smiles_n
        """
    ).df().iloc[0].to_dict()

    target_stats = con.execute(
        f"""
        WITH normalized AS (
            SELECT
                trim(CAST(compound_id AS VARCHAR)) AS compound_id_text,
                CASE
                    WHEN target_gene IS NULL THEN NULL
                    WHEN lower(trim(CAST(target_gene AS VARCHAR))) IN ('', 'nan', 'none', 'null') THEN NULL
                    ELSE upper(trim(CAST(target_gene AS VARCHAR)))
                END AS target_gene_text
            FROM read_parquet('{_p(targets)}')
        )
        SELECT
            count(*) AS target_evidence_rows_n,
            count(DISTINCT compound_id_text) AS compounds_with_target_n,
            count(DISTINCT target_gene_text) AS target_genes_n,
            sum(CASE WHEN target_gene_text IS NULL THEN 1 ELSE 0 END) AS missing_target_gene_n
        FROM normalized
        """
    ).df().iloc[0].to_dict()

    unresolved_counts: dict[str, int] = {}
    unresolved_n = 0
    if UNRESOLVED.exists():
        unresolved = pd.read_csv(UNRESOLVED, sep="\t", low_memory=False)
        unresolved_n = int(len(unresolved))
        if "qc_reason" in unresolved.columns:
            unresolved_counts = {
                str(k): int(v) for k, v in unresolved["qc_reason"].fillna("unknown").value_counts().to_dict().items()
            }

    resolved_n = int(response_summary.get("observations_n") or 0)
    total_normalized = resolved_n + unresolved_n
    payload = {
        "resolved_response_observations_n": resolved_n,
        "unresolved_response_rows_n": unresolved_n,
        "normalized_response_rows_total_n": total_normalized,
        "resolved_fraction": (resolved_n / total_normalized) if total_normalized else None,
        "unresolved_reason_counts": unresolved_counts,
        "response_summary": {k: (int(v) if isinstance(v, int) or str(v).isdigit() else v) for k, v in response_summary.items()},
        "compound_identity": {**compound_stats, **duplicate_identity},
        "target_evidence": target_stats,
        "assay_conditions": assay_rows.to_dict("records"),
        "interpretation_ru": [
            "Неразрешённые response-строки не следует автоматически считать ошибочными экспериментами: чаще это отсутствие модели в текущем CRISPR-Атласе или проблема идентификаторов.",
            "compound_id в PRISM v1 является treatment/profile ID. Число compound_id не следует автоматически трактовать как число уникальных химических структур.",
            "Дублирование Broad ID или SMILES должно быть учтено отдельным химическим identity-слоем до межисточникового объединения веществ.",
            "Assay/endpoint/dose/time сохраняются отдельно; разные endpoints нельзя смешивать в единый показатель чувствительности без заранее описанной нормализации.",
        ],
    }
    QC.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    flat_rows = [
        {"metric": "resolved_response_observations_n", "value": resolved_n},
        {"metric": "unresolved_response_rows_n", "value": unresolved_n},
        {"metric": "normalized_response_rows_total_n", "value": total_normalized},
        {"metric": "resolved_fraction", "value": payload["resolved_fraction"]},
    ]
    for reason, count in unresolved_counts.items():
        flat_rows.append({"metric": f"unresolved::{reason}", "value": count})
    for key, value in compound_stats.items():
        flat_rows.append({"metric": f"compound::{key}", "value": value})
    for key, value in duplicate_identity.items():
        flat_rows.append({"metric": f"compound::{key}", "value": value})
    for key, value in target_stats.items():
        flat_rows.append({"metric": f"target::{key}", "value": value})
    pd.DataFrame(flat_rows).to_csv(OUTPUT_TSV, sep="\t", index=False)

    print("MCL Pharmacology Layer audit")
    print(f"Resolved responses: {resolved_n}")
    print(f"Unresolved responses: {unresolved_n}")
    if total_normalized:
        print(f"Resolved fraction: {resolved_n / total_normalized:.3%}")
    for reason, count in sorted(unresolved_counts.items(), key=lambda x: (-x[1], x[0])):
        print(f"  {reason}: {count}")
    print(f"Wrote {OUTPUT_JSON.relative_to(ROOT)}")
    print(f"Wrote {OUTPUT_TSV.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
