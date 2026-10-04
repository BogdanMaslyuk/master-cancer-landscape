from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml

from mcl.analysis.wave1 import materialize_wave1_cohorts


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
RAW_DEPMAP = ROOT / "data" / "raw" / "depmap"
REPORTS = ROOT / "outputs" / "reports"


def _truthy(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y", "t"})


def _release_from_inventory() -> str | None:
    path = PROCESSED / "depmap_input_inventory.tsv"
    if not path.exists():
        return None
    frame = pd.read_csv(path, sep="\t", usecols=["depmap_release"], nrows=1)
    if frame.empty:
        return None
    value = str(frame.iloc[0]["depmap_release"]).strip()
    return value or None


def _load_models(path: Path) -> pd.DataFrame:
    header = pd.read_csv(path, nrows=0)
    wanted = [
        "ModelID",
        "CellLineName",
        "DepmapModelType",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "OncotreeCode",
    ]
    usecols = [c for c in wanted if c in header.columns]
    if "ModelID" not in usecols or "OncotreeCode" not in usecols:
        raise SystemExit("Model.csv must contain ModelID and OncotreeCode")
    frame = pd.read_csv(path, usecols=usecols, low_memory=False)
    return frame.rename(
        columns={
            "ModelID": "model_id",
            "CellLineName": "cell_line_name",
            "DepmapModelType": "depmap_model_type",
            "OncotreeLineage": "oncotree_lineage",
            "OncotreePrimaryDisease": "oncotree_primary_disease",
            "OncotreeSubtype": "oncotree_subtype",
            "OncotreeCode": "oncotree_code",
        }
    )


def _load_crispr_model_ids(path: Path) -> set[str]:
    header = pd.read_csv(path, nrows=0)
    if len(header.columns) == 0:
        raise SystemExit("CRISPRGeneEffect.csv has no columns")
    first = header.columns[0]
    frame = pd.read_csv(path, usecols=[first], low_memory=False)
    return set(frame[first].dropna().astype(str).str.strip())


def _scan_mutations(
    path: Path,
    crispr_model_ids: set[str],
    genes: set[str],
    chunksize: int,
) -> tuple[pd.DataFrame, set[str]]:
    header = pd.read_csv(path, nrows=0)
    wanted = [
        "ModelID",
        "IsDefaultEntryForModel",
        "ProteinChange",
        "HugoSymbol",
        "VepImpact",
        "OncogeneHighImpact",
        "TumorSuppressorHighImpact",
        "LikelyLoF",
        "HessDriver",
        "Hotspot",
    ]
    usecols = [c for c in wanted if c in header.columns]
    if not {"ModelID", "HugoSymbol"}.issubset(usecols):
        raise SystemExit("OmicsSomaticMutations.csv is missing ModelID/HugoSymbol")

    sequenced: set[str] = set()
    retained: list[pd.DataFrame] = []
    chunks_n = 0
    for chunk in pd.read_csv(path, usecols=usecols, chunksize=chunksize, low_memory=False):
        chunks_n += 1
        chunk["ModelID"] = chunk["ModelID"].astype(str).str.strip()
        chunk = chunk[chunk["ModelID"].isin(crispr_model_ids)].copy()
        if chunk.empty:
            continue
        if "IsDefaultEntryForModel" in chunk.columns:
            default = _truthy(chunk["IsDefaultEntryForModel"])
            if default.any():
                chunk = chunk[default].copy()
        if chunk.empty:
            continue

        sequenced.update(chunk["ModelID"].dropna().astype(str).str.strip())
        gene_mask = chunk["HugoSymbol"].fillna("").astype(str).str.upper().isin(genes)
        chunk = chunk[gene_mask].copy()
        if chunk.empty:
            continue

        rename = {
            "ModelID": "model_id",
            "ProteinChange": "protein_change",
            "HugoSymbol": "gene",
            "VepImpact": "vep_impact",
            "OncogeneHighImpact": "oncogene_high_impact",
            "TumorSuppressorHighImpact": "tumor_suppressor_high_impact",
            "LikelyLoF": "likely_lof",
            "HessDriver": "driver",
            "Hotspot": "hotspot",
        }
        chunk = chunk.rename(columns=rename)
        chunk = chunk.drop(columns=["IsDefaultEntryForModel"], errors="ignore")
        retained.append(chunk)
        if chunks_n % 10 == 0:
            print(f"Scanned {chunks_n} mutation chunks")

    mutations = pd.concat(retained, ignore_index=True) if retained else pd.DataFrame(
        columns=["model_id", "gene", "protein_change"]
    )
    return mutations, sequenced


def main() -> None:
    parser = argparse.ArgumentParser(description="Materialize explicit Scientific Expansion Wave 1 cohorts.")
    parser.add_argument("--chunksize", type=int, default=250_000)
    parser.add_argument("--allow-count-mismatch", action="store_true")
    args = parser.parse_args()

    manifest_path = ROOT / "config" / "scientific_expansion_wave1.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
    contexts = list(manifest.get("contexts") or [])
    if not contexts:
        raise SystemExit("Wave 1 manifest has no contexts")

    release = str(manifest.get("depmap_release") or _release_from_inventory() or "").strip()
    if not release:
        raise SystemExit("Cannot resolve DepMap release")
    release_dir = RAW_DEPMAP / release
    model_path = release_dir / "Model.csv"
    gene_effect_path = release_dir / "CRISPRGeneEffect.csv"
    mutation_path = release_dir / "OmicsSomaticMutations.csv"
    for path in [model_path, gene_effect_path, mutation_path]:
        if not path.exists():
            raise SystemExit(f"Missing required DepMap input: {path}")

    genes = {str(x.get("gene") or "").strip().upper() for x in contexts}
    genes.discard("")

    print(f"DepMap release: {release}")
    print("Loading model metadata...")
    models = _load_models(model_path)
    print("Loading CRISPR model coverage...")
    crispr_ids = _load_crispr_model_ids(gene_effect_path)
    print(f"CRISPR-covered models: {len(crispr_ids):,}")
    print(f"Scanning mutations for Wave 1 genes: {', '.join(sorted(genes))}")
    mutations, sequenced_ids = _scan_mutations(mutation_path, crispr_ids, genes, args.chunksize)
    print(f"Mutation-profiled CRISPR models observed: {len(sequenced_ids):,}")

    cohort, summary = materialize_wave1_cohorts(
        models,
        mutations,
        contexts,
        crispr_model_ids=crispr_ids,
        sequenced_model_ids=sequenced_ids,
    )

    out_dir = PROCESSED / "scientific_expansion_wave1"
    out_dir.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    cohort.to_csv(out_dir / "cohorts.tsv", sep="\t", index=False)
    cohort.to_parquet(out_dir / "cohorts.parquet", index=False)
    summary.to_csv(out_dir / "cohort_summary.tsv", sep="\t", index=False)
    summary.to_parquet(out_dir / "cohort_summary.parquet", index=False)

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "contract": "mcl-scientific-expansion-wave1-cohorts-v1",
        "depmap_release": release,
        "manifest": "config/scientific_expansion_wave1.yaml",
        "contexts_n": int(len(summary)),
        "all_counts_match_discovery": bool(summary["counts_match_discovery"].all()),
        "contexts": summary.to_dict("records"),
    }
    (REPORTS / "scientific_expansion_wave1_cohorts.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("\nWave 1 cohorts:")
    show = [
        "wave1_id", "label", "comparison_mode", "context_models_n", "comparator_models_n",
        "expected_context_n", "expected_comparator_n", "counts_match_discovery",
    ]
    print(summary[show].to_string(index=False))
    print("\nWrote:")
    print("  data/processed/scientific_expansion_wave1/cohorts.tsv")
    print("  data/processed/scientific_expansion_wave1/cohorts.parquet")
    print("  data/processed/scientific_expansion_wave1/cohort_summary.tsv")
    print("  outputs/reports/scientific_expansion_wave1_cohorts.json")

    if not bool(summary["counts_match_discovery"].all()) and not args.allow_count_mismatch:
        raise SystemExit(
            "Wave 1 cohort counts do not reproduce registry discovery. Genome-wide analysis is blocked until the mismatch is reviewed."
        )


if __name__ == "__main__":
    main()
