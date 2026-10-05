from __future__ import annotations

"""Compatibility entry point for Laboratory Mechanism Context v1.1.

The main implementation remains in build_laboratory_mechanism_panels.py. This
entry point overrides only the PRISM response loader so DuckDB receives fully
qualified GROUP BY expressions after joining the selected compound table.
"""

import warnings

import pandas as pd

import build_laboratory_mechanism_panels as impl


# pandas 2.x warns about a future dtype change in two defensive fillna calls in
# the implementation. The values are immediately normalized through _safe_bool,
# so the warning does not affect the scientific result.
warnings.filterwarnings(
    "ignore",
    message=r"Downcasting object dtype arrays on \.fillna.*",
    category=FutureWarning,
)


def _load_responses(compound_ids: list[str], model_ids: list[str]) -> pd.DataFrame:
    columns = ["model_id", "compound_id", "response_value", "relative_sensitivity"]
    if not impl.RESPONSES.exists() or not compound_ids or not model_ids:
        return pd.DataFrame(columns=columns)

    con = impl.duckdb.connect(database=":memory:")
    try:
        con.register(
            "selected_compounds",
            pd.DataFrame({"compound_id": sorted(set(compound_ids))}),
        )
        con.register(
            "selected_models",
            pd.DataFrame({"model_id": sorted(set(model_ids))}),
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

    return frame if not frame.empty else pd.DataFrame(columns=columns)


impl._load_responses = _load_responses


if __name__ == "__main__":
    impl.main()
