from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "runtime" / "pharmacology" / "candidate_hypothesis_models_v2.parquet"
OUTPUT = ROOT / "outputs" / "qc" / "candidate_hypothesis_models_v2_qc.tsv"

ACTIVE_LFC = -1.0
SENSITIVE_PERCENTILE = 0.75
NONRESPONSIVE_PERCENTILE = 0.25
DEPENDENCY_PROBABILITY = 0.5


def main() -> None:
    if not INPUT.exists():
        raise SystemExit(f"Missing {INPUT.relative_to(ROOT)}")
    frame = pd.read_parquet(INPUT)
    if frame.empty:
        raise SystemExit("Candidate v2 model recommendation table is empty.")

    response = pd.to_numeric(frame["response_value"], errors="coerce")
    relative = pd.to_numeric(frame["relative_sensitivity"], errors="coerce")
    probability = pd.to_numeric(frame["dependency_probability"], errors="coerce")
    role = frame["model_role"].astype(str)

    checks: list[dict[str, object]] = []

    def add(name: str, mask: pd.Series, expectation: str) -> None:
        bad = frame[mask].copy()
        checks.append(
            {
                "check": name,
                "violations_n": int(len(bad)),
                "expectation": expectation,
                "example_model_ids": ",".join(bad.get("model_id", pd.Series(dtype=str)).astype(str).head(5)),
            }
        )

    positive = role.eq("positive")
    add(
        "positive_requires_absolute_activity",
        positive & ~(response <= ACTIVE_LFC),
        f"positive response_value <= {ACTIVE_LFC}",
    )
    add(
        "positive_requires_relative_selectivity",
        positive & ~(relative >= SENSITIVE_PERCENTILE),
        f"positive relative_sensitivity >= {SENSITIVE_PERCENTILE}",
    )
    add(
        "positive_requires_dependency_probability",
        positive & ~(probability > DEPENDENCY_PROBABILITY),
        f"positive Probability of Dependency > {DEPENDENCY_PROBABILITY}",
    )

    negative = role.eq("negative_same_cancer")
    add(
        "negative_requires_measured_probability",
        negative & probability.isna(),
        "negative control must have measured Probability of Dependency",
    )
    add(
        "negative_requires_non_dependency",
        negative & ~(probability <= DEPENDENCY_PROBABILITY),
        f"negative Probability of Dependency <= {DEPENDENCY_PROBABILITY}",
    )
    add(
        "negative_requires_no_absolute_activity",
        negative & ~(response > ACTIVE_LFC),
        f"negative response_value > {ACTIVE_LFC}",
    )
    add(
        "negative_requires_low_relative_sensitivity",
        negative & ~(relative <= NONRESPONSIVE_PERCENTILE),
        f"negative relative_sensitivity <= {NONRESPONSIVE_PERCENTILE}",
    )

    missing = role.eq("unclassified_missing_dependency_probability")
    add(
        "missing_probability_role_must_really_be_missing",
        missing & probability.notna(),
        "unclassified_missing_dependency_probability must have missing probability",
    )

    out = pd.DataFrame(checks)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT, sep="\t", index=False)

    total_violations = int(out["violations_n"].sum()) if not out.empty else 0
    print("Candidate v2 model-role QC")
    print(f"Rows: {len(frame)}")
    print(f"Total rule violations: {total_violations}")
    for _, row in out.iterrows():
        print(f"  {row['check']}: {row['violations_n']}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")
    if total_violations:
        raise SystemExit("Candidate v2 model-role QC failed. Do not use the model recommendations until fixed.")


if __name__ == "__main__":
    main()
