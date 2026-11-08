from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "data" / "runtime" / "pharmacology" / "candidate_hypothesis_models.parquet"
QC = ROOT / "outputs" / "qc" / "candidate_hypothesis_model_roles.tsv"


def main() -> None:
    if not PATH.exists():
        raise SystemExit("Missing candidate_hypothesis_models.parquet. Build hypotheses first.")
    frame = pd.read_parquet(PATH)
    if frame.empty:
        print("Candidate model table is empty; nothing to QC.")
        return

    gene_effect = pd.to_numeric(frame.get("gene_effect"), errors="coerce")
    needs_measured_crispr = frame["model_role"].isin(
        ["negative_same_cancer", "discordant_sensitive_without_dependency"]
    )
    invalid = needs_measured_crispr & gene_effect.isna()
    invalid_n = int(invalid.sum())
    if invalid_n:
        frame.loc[invalid, "model_role"] = "unclassified_missing_crispr"
        frame.loc[invalid, "role_reason_ru"] = (
            "CRISPR Gene Effect для мишени отсутствует; отсутствие измерения нельзя трактовать как отсутствие зависимости."
        )
        frame.to_parquet(PATH, index=False, compression="zstd")

    counts = frame["model_role"].value_counts().rename_axis("model_role").reset_index(name="models_n")
    QC.parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(QC, sep="\t", index=False)

    print("Candidate hypothesis model-role QC")
    print(f"Reclassified missing-CRISPR rows: {invalid_n}")
    for row in counts.itertuples(index=False):
        print(f"  {row.model_role}: {row.models_n}")
    print(f"Wrote {QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
