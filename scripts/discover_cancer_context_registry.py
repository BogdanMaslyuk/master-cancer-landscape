from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from mcl.analysis.context_registry import discover_context_candidates


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


def _resolve_release(explicit: str | None) -> str:
    if explicit:
        return explicit
    from_inventory = _release_from_inventory()
    if from_inventory:
        return from_inventory
    candidates = sorted([p.name for p in RAW_DEPMAP.iterdir() if p.is_dir()], reverse=True) if RAW_DEPMAP.exists() else []
    if candidates:
        return candidates[0]
    raise SystemExit("Cannot infer DepMap release. Pass --release, e.g. --release 26Q1")


def _load_models(path: Path) -> pd.DataFrame:
    header = pd.read_csv(path, nrows=0)
    wanted = [
        "ModelID",
        "DepmapModelType",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "OncotreeCode",
    ]
    usecols = [c for c in wanted if c in header.columns]
    if "ModelID" not in usecols:
        raise SystemExit("Model.csv is missing ModelID")
    frame = pd.read_csv(path, usecols=usecols, low_memory=False)
    return frame.rename(
        columns={
            "ModelID": "model_id",
            "DepmapModelType": "depmap_model_type",
            "OncotreeLineage": "oncotree_lineage",
            "OncotreePrimaryDisease": "oncotree_primary_disease",
            "OncotreeSubtype": "oncotree_subtype",
            "OncotreeCode": "oncotree_code",
        }
    )


def _load_crispr_model_ids(path: Path) -> set[str]:
    header = pd.read_csv(path, nrows=0)
    if header.empty:
        raise SystemExit("CRISPRGeneEffect.csv has no columns")
    first = header.columns[0]
    frame = pd.read_csv(path, usecols=[first], low_memory=False)
    return set(frame[first].dropna().astype(str).str.strip())


def _scan_mutations(path: Path, crispr_model_ids: set[str], chunksize: int) -> tuple[pd.DataFrame, set[str]]:
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
    required = {"ModelID", "HugoSymbol"}
    if not required.issubset(usecols):
        raise SystemExit(f"OmicsSomaticMutations.csv is missing columns: {sorted(required - set(usecols))}")

    sequenced_model_ids: set[str] = set()
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

        sequenced_model_ids.update(chunk["ModelID"].dropna().astype(str).str.strip())

        flags: dict[str, pd.Series] = {}
        for source in ["HessDriver", "Hotspot", "OncogeneHighImpact", "TumorSuppressorHighImpact", "LikelyLoF"]:
            flags[source] = _truthy(chunk[source]) if source in chunk.columns else pd.Series(False, index=chunk.index)
        impact = chunk.get("VepImpact", pd.Series("", index=chunk.index)).fillna("").astype(str).str.upper()
        qualifying = (
            flags["HessDriver"]
            | flags["Hotspot"]
            | flags["OncogeneHighImpact"]
            | flags["TumorSuppressorHighImpact"]
            | flags["LikelyLoF"]
            | impact.eq("HIGH")
        )
        chunk = chunk[qualifying].copy()
        if chunk.empty:
            continue

        chunk = chunk.rename(
            columns={
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
        )
        retained.append(chunk)
        if chunks_n % 10 == 0:
            print(f"Scanned {chunks_n} mutation chunks; retained {sum(len(x) for x in retained):,} qualifying rows")

    if retained:
        mutations = pd.concat(retained, ignore_index=True)
    else:
        mutations = pd.DataFrame(columns=["model_id", "gene", "protein_change"])
    return mutations, sequenced_model_ids


def _configured_keys() -> set[tuple[str, str, str]]:
    path = ROOT / "config" / "cancer_contexts.yaml"
    if not path.exists():
        return set()
    import yaml

    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    output: set[tuple[str, str, str]] = set()
    for cfg in (payload.get("contexts") or {}).values():
        cfg = cfg or {}
        inclusion = cfg.get("inclusion") or {}
        depmap = cfg.get("depmap") or {}
        metadata = depmap.get("metadata") or {}
        gene = str(inclusion.get("alteration_gene") or "").strip().upper()
        alteration = str(inclusion.get("alteration") or "").strip()
        codes = metadata.get("oncotree_codes") or []
        for code in codes:
            output.add((str(code).strip(), gene, alteration))
    return output


def _mark_existing(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        frame["already_configured"] = pd.Series(dtype=bool)
        return frame
    configured = _configured_keys()
    result = frame.copy()
    flags = []
    for row in result.itertuples(index=False):
        code = str(getattr(row, "oncotree_code", "") or "").strip()
        gene = str(getattr(row, "alteration_gene", "") or "").strip().upper()
        alteration = str(getattr(row, "alteration", "") or "").strip()
        flags.append((code, gene, alteration) in configured)
    result["already_configured"] = flags
    return result


def _write_summary(frame: pd.DataFrame, release: str, args: argparse.Namespace) -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)

    frame.to_csv(PROCESSED / "cancer_context_registry_candidates.tsv", sep="\t", index=False)
    frame.to_parquet(PROCESSED / "cancer_context_registry_candidates.parquet", index=False)

    status_counts = frame["readiness_status"].value_counts().to_dict() if not frame.empty else {}
    tier_counts = frame["priority_tier"].value_counts().to_dict() if not frame.empty else {}
    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "depmap_release": release,
        "candidate_contract": "mcl-cancer-context-registry-candidates-v1",
        "principle": (
            "Discovery candidates are not automatically promoted into config/cancer_contexts.yaml. "
            "Promotion requires biological review of disease scope, molecular event and comparator design."
        ),
        "thresholds": {
            "min_discovery_n": args.min_discovery_n,
            "min_context_n": args.min_context_n,
            "min_comparator_n": args.min_comparator_n,
        },
        "candidates_n": int(len(frame)),
        "readiness_counts": status_counts,
        "priority_tier_counts": tier_counts,
        "ready_not_configured_n": int(
            ((frame.get("readiness_status") == "READY") & (~frame.get("already_configured", False))).sum()
        ) if not frame.empty else 0,
    }
    (REPORTS / "cancer_context_registry_candidates.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    ready = frame[(frame["readiness_status"] == "READY") & (~frame["already_configured"])].copy() if not frame.empty else frame
    columns = [
        "candidate_id",
        "oncotree_code",
        "oncotree_subtype",
        "alteration_gene",
        "alteration",
        "comparison_mode",
        "context_models_n",
        "comparator_models_n",
        "priority_tier",
    ]
    md = [
        "# Cancer Context Registry v1 — discovery report",
        "",
        f"DepMap release: **{release}**",
        "",
        f"Discovered candidates: **{len(frame)}**",
        f"READY and not yet configured: **{len(ready)}**",
        "",
        "Candidates are discovery hypotheses, not automatically accepted cancer contexts.",
        "",
        "## READY candidates not yet configured",
        "",
    ]
    if ready.empty:
        md.append("No additional READY candidates were found under the current thresholds.")
    else:
        subset = ready[columns].head(args.report_top)
        md.append("| " + " | ".join(columns) + " |")
        md.append("| " + " | ".join(["---"] * len(columns)) + " |")
        for row in subset.itertuples(index=False):
            md.append("| " + " | ".join(str(x).replace("|", "/") for x in row) + " |")
    (REPORTS / "cancer_context_registry_v1.md").write_text("\n".join(md) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Discover scalable DepMap cancer-context comparison candidates.")
    parser.add_argument("--release", default=None)
    parser.add_argument("--chunksize", type=int, default=250_000)
    parser.add_argument("--min-discovery-n", type=int, default=3)
    parser.add_argument("--min-context-n", type=int, default=5)
    parser.add_argument("--min-comparator-n", type=int, default=5)
    parser.add_argument("--report-top", type=int, default=100)
    args = parser.parse_args()

    release = _resolve_release(args.release)
    release_dir = RAW_DEPMAP / release
    model_path = release_dir / "Model.csv"
    mutation_path = release_dir / "OmicsSomaticMutations.csv"
    gene_effect_path = release_dir / "CRISPRGeneEffect.csv"
    for path in [model_path, mutation_path, gene_effect_path]:
        if not path.exists():
            raise SystemExit(f"Missing required DepMap input: {path}")

    print(f"DepMap release: {release}")
    print("Loading model metadata...")
    models = _load_models(model_path)
    print("Loading CRISPR model coverage...")
    crispr_model_ids = _load_crispr_model_ids(gene_effect_path)
    print(f"CRISPR-covered models: {len(crispr_model_ids):,}")
    print("Scanning somatic mutations once for recurrent driver/hotspot/high-impact contexts...")
    mutations, sequenced_model_ids = _scan_mutations(mutation_path, crispr_model_ids, args.chunksize)
    print(f"Mutation-profiled CRISPR models observed: {len(sequenced_model_ids):,}")
    print(f"Qualifying mutation rows retained: {len(mutations):,}")

    candidates = discover_context_candidates(
        models,
        mutations,
        crispr_model_ids=crispr_model_ids,
        sequenced_model_ids=sequenced_model_ids,
        min_discovery_n=args.min_discovery_n,
        min_context_n=args.min_context_n,
        min_comparator_n=args.min_comparator_n,
    )
    candidates = _mark_existing(candidates)
    _write_summary(candidates, release, args)

    if candidates.empty:
        print("No candidates met the discovery threshold.")
    else:
        print("\nReadiness:")
        print(candidates["readiness_status"].value_counts().to_string())
        ready = candidates[(candidates["readiness_status"] == "READY") & (~candidates["already_configured"])].copy()
        print(f"\nREADY and not configured: {len(ready)}")
        if not ready.empty:
            show = [
                "oncotree_code", "oncotree_subtype", "alteration_gene", "alteration",
                "comparison_mode", "context_models_n", "comparator_models_n", "priority_tier",
            ]
            print(ready[show].head(30).to_string(index=False))

    print("\nWrote:")
    print("  data/processed/cancer_context_registry_candidates.tsv")
    print("  data/processed/cancer_context_registry_candidates.parquet")
    print("  outputs/reports/cancer_context_registry_candidates.json")
    print("  outputs/reports/cancer_context_registry_v1.md")


if __name__ == "__main__":
    main()
