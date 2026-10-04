from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml

from mcl.analysis.wave1_genomewide import analyze_wave1_genomewide


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
REPORTS = ROOT / "outputs" / "reports"


def main() -> None:
    manifest_path = ROOT / "config" / "scientific_expansion_wave1.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
    release = str(manifest.get("depmap_release") or "").strip()
    if not release:
        raise SystemExit("scientific_expansion_wave1.yaml has no depmap_release")

    wave1_dir = PROCESSED / "scientific_expansion_wave1"
    cohort_path = wave1_dir / "cohorts.parquet"
    summary_path = wave1_dir / "cohort_summary.parquet"
    if not cohort_path.exists() or not summary_path.exists():
        raise SystemExit(
            "Wave 1 cohorts are missing. Run scripts/materialize_scientific_expansion_wave1.py first."
        )

    cohort_summary = pd.read_parquet(summary_path)
    if cohort_summary.empty or not bool(cohort_summary["counts_match_discovery"].all()):
        raise SystemExit(
            "Wave 1 cohort counts do not reproduce registry discovery. Genome-wide analysis is blocked."
        )

    cohorts = pd.read_parquet(cohort_path)
    print(f"DepMap release: {release}")
    print(f"Wave 1 contexts: {cohort_summary['wave1_id'].nunique()}")
    print("Loading genome-wide Gene Effect and Gene Dependency matrices once for the Wave 1 union cohort...")
    results, summary, meta = analyze_wave1_genomewide(ROOT, release, cohorts)

    out_dir = wave1_dir / "genomewide"
    out_dir.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)

    results.to_parquet(out_dir / "genes.parquet", index=False)
    results.to_csv(out_dir / "genes.tsv", sep="\t", index=False)
    summary.to_parquet(out_dir / "summary.parquet", index=False)
    summary.to_csv(out_dir / "summary.tsv", sep="\t", index=False)

    top_rows = []
    for wave1_id, group in results.groupby("wave1_id", sort=True):
        eligible = group[~group["broad_dependency_warning"].fillna(True)].copy()
        eligible = eligible.sort_values(
            ["delta_gene_effect", "q_value", "gene_symbol"],
            ascending=[True, True, True],
            na_position="last",
        ).head(50)
        top_rows.append(eligible)
    top = pd.concat(top_rows, ignore_index=True) if top_rows else results.iloc[0:0].copy()
    top.to_csv(REPORTS / "scientific_expansion_wave1_top50_by_delta.tsv", sep="\t", index=False)

    meta_payload = {
        **meta,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "manifest": "config/scientific_expansion_wave1.yaml",
        "cohort_source": "data/processed/scientific_expansion_wave1/cohorts.parquet",
        "top50_note": (
            "Top50 is a descriptive view ordered by most negative delta_gene_effect after excluding "
            "broad-dependency warnings. It is not a validated target list or composite ranking."
        ),
        "contexts": summary.to_dict("records"),
    }
    (REPORTS / "scientific_expansion_wave1_genomewide_meta.json").write_text(
        json.dumps(meta_payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("\nWave 1 genome-wide summary:")
    show = [
        "wave1_id",
        "label",
        "context_models_n",
        "comparator_models_n",
        "genes_tested_n",
        "fdr_0_05_genes_n",
        "negative_delta_genes_n",
        "low_sample_size_genes_n",
    ]
    print(summary[show].to_string(index=False))

    print("\nWrote:")
    print("  data/processed/scientific_expansion_wave1/genomewide/genes.tsv")
    print("  data/processed/scientific_expansion_wave1/genomewide/genes.parquet")
    print("  data/processed/scientific_expansion_wave1/genomewide/summary.tsv")
    print("  outputs/reports/scientific_expansion_wave1_top50_by_delta.tsv")
    print("  outputs/reports/scientific_expansion_wave1_genomewide_meta.json")


if __name__ == "__main__":
    main()
