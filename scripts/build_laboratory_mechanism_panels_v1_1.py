from __future__ import annotations

"""Compatibility entry point for Laboratory Mechanism Context v1.1.

The main implementation remains in build_laboratory_mechanism_panels.py. This
entry point applies two narrow runtime fixes:
1) DuckDB receives fully-qualified GROUP BY expressions after joining the
   selected compound table.
2) zero-argument temporary DataFrames have a stable empty schema, so physical
   controls without a DepMap model (BJ5ta) can pass through the mechanism-panel
   builder without a KeyError on target_gene.

No scientific thresholds or classifications are changed here.
"""

import warnings

import pandas as pd

import build_laboratory_mechanism_panels as impl


warnings.filterwarnings(
    "ignore",
    message=r"Downcasting object dtype arrays on \.fillna.*",
    category=FutureWarning,
)


# Keep a reference to the real pandas constructor before replacing the module
# attribute used by the implementation.
_REAL_DATAFRAME = pd.DataFrame
_EMPTY_SCHEMA = [
    "model_id",
    "target_gene",
    "compound_id",
    "hypothesis_id",
    "lab_name",
    "laboratory_model_role",
    "response_value",
    "relative_sensitivity",
    "dependency_probability",
    "gene_effect",
]


def _schema_safe_dataframe(*args, **kwargs):
    """Preserve normal DataFrame construction; stabilize only DataFrame()."""
    if not args and not kwargs:
        return _REAL_DATAFRAME(columns=_EMPTY_SCHEMA)
    return _REAL_DATAFRAME(*args, **kwargs)


def _load_responses(compound_ids: list[str], model_ids: list[str]) -> pd.DataFrame:
    columns = ["model_id", "compound_id", "response_value", "relative_sensitivity"]
    if not impl.RESPONSES.exists() or not compound_ids or not model_ids:
        return _REAL_DATAFRAME(columns=columns)

    con = impl.duckdb.connect(database=":memory:")
    try:
        con.register(
            "selected_compounds",
            _REAL_DATAFRAME({"compound_id": sorted(set(compound_ids))}),
        )
        con.register(
            "selected_models",
            _REAL_DATAFRAME({"model_id": sorted(set(model_ids))}),
        )
        path = impl.RESPONSES.resolve().as_posix().replace("'", "''")
        frame = con.execute(
            f"""
            WITH response_by_model AS (
              SELECT CAST(r.model_id AS VARCHAR) AS model_id,
                     CAST(r.compound_id AS VARCHAR) AS compound_id,
                     avg(CAST(r.value AS DOUBLE)) AS response_value
              FROM read_parquet('{path}') r
              INNER JOIN selected_compounds c
                ON CAST(r.compound_id AS VARCHAR) = c.compound_id
              WHERE upper(CAST(r.source AS VARCHAR)) = 'PRISM'
                AND upper(CAST(r.endpoint AS VARCHAR)) = 'LFC'
                AND r.value IS NOT NULL
              GROUP BY r.model_id, r.compound_id
            ), ranked AS (
              SELECT *,
                     1.0 - percent_rank() OVER (
                       PARTITION BY compound_id ORDER BY response_value ASC
                     ) AS relative_sensitivity
              FROM response_by_model
            )
            SELECT r.*
            FROM ranked r
            INNER JOIN selected_models m ON r.model_id = m.model_id
            """
        ).df()
    finally:
        con.close()

    return frame if not frame.empty else _REAL_DATAFRAME(columns=columns)


impl._load_responses = _load_responses
impl.pd.DataFrame = _schema_safe_dataframe


if __name__ == "__main__":
    impl.main()
