from __future__ import annotations

from pathlib import Path

import pandas as pd

from mcl.analysis.context_curation import curate_registry_candidates


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
REPORTS = ROOT / "outputs" / "reports"


def main() -> None:
    source = PROCESSED / "cancer_context_registry_candidates.tsv"
    if not source.exists():
        raise SystemExit(
            "Missing data/processed/cancer_context_registry_candidates.tsv. "
            "Run scripts/discover_cancer_context_registry.py first."
        )

    frame = pd.read_csv(source, sep="\t", low_memory=False)
    curated = curate_registry_candidates(frame)

    PROCESSED.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    curated.to_csv(PROCESSED / "cancer_context_registry_curated.tsv", sep="\t", index=False)
    curated.to_parquet(PROCESSED / "cancer_context_registry_curated.parquet", index=False)

    unique = curated[
        curated["comparison_role"].isin({"primary_comparator", "primary_available_comparator"})
    ].copy()
    first_wave = unique[unique["review_class"] == "DRIVER_HOTSPOT_FIRST_REVIEW"].copy()
    lof_review = unique[unique["review_class"] == "LOF_SECOND_REVIEW"].copy()

    md = [
        "# Cancer Context Registry v1 — biological curation queue",
        "",
        "This report is a review queue, not an automatically accepted registry.",
        "No candidate is promoted into config/cancer_contexts.yaml by this script.",
        "",
        f"READY exact-variant comparisons reviewed: **{len(curated)}**",
        f"Unique contexts with a preferred comparator: **{len(unique)}**",
        f"Driver/hotspot first-review contexts: **{len(first_wave)}**",
        f"Loss-of-function second-review contexts: **{len(lof_review)}**",
        "",
        "## First-review driver/hotspot contexts",
        "",
    ]

    columns = [
        "oncotree_code",
        "oncotree_subtype",
        "alteration_gene",
        "alteration",
        "comparison_mode",
        "context_models_n",
        "comparator_models_n",
        "driver_context_models_n",
        "hotspot_context_models_n",
        "comparison_role",
    ]
    if first_wave.empty:
        md.append("No driver/hotspot contexts met the current transparent curation rule.")
    else:
        subset = first_wave[columns].head(100)
        md.append("| " + " | ".join(columns) + " |")
        md.append("| " + " | ".join(["---"] * len(columns)) + " |")
        for row in subset.itertuples(index=False):
            md.append("| " + " | ".join(str(x).replace("|", "/") for x in row) + " |")

    md.extend(
        [
            "",
            "## Interpretation rules",
            "",
            "- DRIVER_HOTSPOT_FIRST_REVIEW: exact non-truncating event with driver or hotspot support in at least 80% of context models.",
            "- LOF_SECOND_REVIEW: truncating/frameshift event with strong high-impact/LoF/driver/hotspot support; co-mutation burden must be reviewed before use.",
            "- OTHER_REVIEW: READY exact event that does not meet either first-pass evidence rule.",
            "- Wild-type comparator is preferred when available; otherwise another-variant comparator is retained as the primary available allele-specific comparison.",
            "- No opaque global score is calculated.",
        ]
    )
    (REPORTS / "cancer_context_registry_curation_v1.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    print("Curation classes:")
    if curated.empty:
        print("No READY exact-variant candidates available for curation.")
    else:
        print(curated["review_class"].value_counts().to_string())
        print(f"\nUnique contexts with preferred comparator: {len(unique)}")
        print(f"Driver/hotspot first-review contexts: {len(first_wave)}")
        if not first_wave.empty:
            show = [
                "oncotree_code",
                "oncotree_subtype",
                "alteration_gene",
                "alteration",
                "comparison_mode",
                "context_models_n",
                "comparator_models_n",
                "driver_context_models_n",
                "hotspot_context_models_n",
            ]
            print("\nFIRST REVIEW:")
            print(first_wave[show].head(40).to_string(index=False))

    print("\nWrote:")
    print("  data/processed/cancer_context_registry_curated.tsv")
    print("  data/processed/cancer_context_registry_curated.parquet")
    print("  outputs/reports/cancer_context_registry_curation_v1.md")


if __name__ == "__main__":
    main()
