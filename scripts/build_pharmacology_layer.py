from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "data" / "raw" / "pharmacology" / "normalized"
RUNTIME = ROOT / "data" / "runtime" / "pharmacology"
QC = ROOT / "outputs" / "qc"
ATLAS = ROOT / "data" / "processed" / "depmap_crispr_model_atlas.parquet"
GENE_EFFECT = ROOT / "data" / "processed" / "depmap_model_gene_effect.parquet"

COMPOUND_COLUMNS = [
    "compound_id", "preferred_name", "canonical_smiles", "inchikey", "pubchem_cid",
    "chembl_id", "broad_id", "gdsc_id", "source_ids_json",
]
RESPONSE_COLUMNS = [
    "observation_id", "model_id", "compound_id", "source", "source_release", "source_assay_id",
    "assay_type", "endpoint", "value", "unit", "auc", "ic50", "ec50", "gi50", "viability",
    "dose", "dose_unit", "exposure_time_h", "replicate_n", "quality_flag", "publication",
]
TARGET_COLUMNS = [
    "evidence_id", "compound_id", "target_gene", "action", "evidence_type", "activity_type",
    "activity_value", "activity_unit", "source", "source_assay_id", "publication", "confidence",
    "directness",
]
VALID_SUFFIXES = {".parquet", ".tsv", ".csv"}


def _read(path: Path, columns: list[str]) -> pd.DataFrame:
    if path.suffix.lower() == ".parquet":
        frame = pd.read_parquet(path)
    else:
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        frame = pd.read_csv(path, sep=sep, low_memory=False)
    for column in columns:
        if column not in frame.columns:
            frame[column] = None
    out = frame[columns].copy()
    out["_normalized_source_dir"] = str(path.parent)
    return out


def _find_all(input_dir: Path, stem: str) -> list[Path]:
    paths: list[Path] = []
    if not input_dir.exists():
        return paths
    for path in input_dir.rglob(f"{stem}.*"):
        if path.is_file() and path.suffix.lower() in VALID_SUFFIXES:
            paths.append(path)
    return sorted(set(paths))


def _read_many(input_dir: Path, stem: str, columns: list[str]) -> pd.DataFrame:
    paths = _find_all(input_dir, stem)
    if not paths:
        return pd.DataFrame(columns=[*columns, "_normalized_source_dir"])
    frames = [_read(path, columns) for path in paths]
    return pd.concat(frames, ignore_index=True, sort=False)


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _normalise_compounds(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["compound_id"] = out["compound_id"].map(_text)
    out = out[out["compound_id"] != ""].copy()
    if out.empty:
        return out.drop(columns=["_normalized_source_dir"], errors="ignore")

    # A stable compound_id is the join key. If multiple sources deliberately use
    # the same ID, retain the first non-empty metadata record rather than multiplying
    # response rows. Cross-source chemical identity by InChIKey is a later layer.
    out["_completeness"] = out[COMPOUND_COLUMNS].notna().sum(axis=1)
    out = out.sort_values("_completeness", ascending=False).drop_duplicates("compound_id", keep="first")
    return out.drop(columns=["_completeness", "_normalized_source_dir"], errors="ignore").reset_index(drop=True)


def _normalise_responses(
    frame: pd.DataFrame,
    valid_models: set[str],
    valid_compounds: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = frame.copy()
    for col in ("observation_id", "model_id", "compound_id", "source", "endpoint", "unit"):
        out[col] = out[col].map(_text)
    missing = out["observation_id"].eq("")
    out.loc[missing, "observation_id"] = [f"PHARM-OBS-{i:08d}" for i in range(1, int(missing.sum()) + 1)]
    for col in ("value", "auc", "ic50", "ec50", "gi50", "viability", "dose", "exposure_time_h", "replicate_n"):
        out[col] = pd.to_numeric(out[col], errors="coerce")

    reason = pd.Series("", index=out.index, dtype=object)
    reason.loc[~out["model_id"].isin(valid_models)] = "unknown_model_id"
    bad_compound = ~out["compound_id"].isin(valid_compounds)
    reason.loc[bad_compound] = reason.loc[bad_compound].map(
        lambda x: f"{x};unknown_compound_id".strip(";")
    )
    unresolved = out[reason.ne("")].copy()
    unresolved["qc_reason"] = reason[reason.ne("")]
    resolved = out[reason.eq("")].copy().drop_duplicates("observation_id", keep="first")
    return (
        resolved.drop(columns=["_normalized_source_dir"], errors="ignore").reset_index(drop=True),
        unresolved.reset_index(drop=True),
    )


def _normalise_targets(
    frame: pd.DataFrame,
    valid_compounds: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = frame.copy()
    for col in ("evidence_id", "compound_id", "target_gene", "action", "evidence_type", "source", "confidence", "directness"):
        out[col] = out[col].map(_text)
    out["target_gene"] = out["target_gene"].str.upper()
    missing = out["evidence_id"].eq("")
    out.loc[missing, "evidence_id"] = [f"PHARM-TGT-{i:08d}" for i in range(1, int(missing.sum()) + 1)]
    out["activity_value"] = pd.to_numeric(out["activity_value"], errors="coerce")

    reason = pd.Series("", index=out.index, dtype=object)
    bad_compound = ~out["compound_id"].isin(valid_compounds)
    reason.loc[bad_compound] = "unknown_compound_id"
    missing_gene = out["target_gene"].eq("")
    reason.loc[missing_gene] = reason.loc[missing_gene].map(
        lambda x: f"{x};missing_target_gene".strip(";")
    )
    unresolved = out[reason.ne("")].copy()
    unresolved["qc_reason"] = reason[reason.ne("")]
    resolved = out[reason.eq("")].copy().drop_duplicates("evidence_id", keep="first")
    return (
        resolved.drop(columns=["_normalized_source_dir"], errors="ignore").reset_index(drop=True),
        unresolved.reset_index(drop=True),
    )


def _gene_effect_lookup(model_ids: set[str], target_genes: set[str]) -> dict[tuple[str, str], float]:
    if not GENE_EFFECT.exists() or not model_ids or not target_genes:
        return {}
    schema_names = set(pq.ParquetFile(GENE_EFFECT).schema.names)
    genes = sorted(g for g in target_genes if g in schema_names)
    if not genes:
        return {}
    frame = pd.read_parquet(GENE_EFFECT, columns=["model_id", *genes])
    frame["model_id"] = frame["model_id"].astype(str)
    frame = frame[frame["model_id"].isin(model_ids)].copy()
    lookup: dict[tuple[str, str], float] = {}
    for gene in genes:
        values = pd.to_numeric(frame[gene], errors="coerce")
        for model_id, value in zip(frame["model_id"], values):
            if pd.notna(value):
                lookup[(str(model_id), gene)] = float(value)
    return lookup


def _crispr_level(value: float | None) -> str:
    if value is None:
        return "not_available"
    if value <= -1.0:
        return "strong_dependency"
    if value <= -0.5:
        return "dependency"
    return "weak_or_none"


def _build_links(responses: pd.DataFrame, targets: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "model_id", "compound_id", "target_gene", "action", "evidence_type", "confidence",
        "directness", "source", "activity_type", "activity_value", "activity_unit", "gene_effect",
        "crispr_support_level",
    ]
    if responses.empty or targets.empty:
        return pd.DataFrame(columns=columns)
    pairs = responses[["model_id", "compound_id"]].drop_duplicates()
    joined = pairs.merge(targets, on="compound_id", how="inner")
    if joined.empty:
        return pd.DataFrame(columns=columns)
    lookup = _gene_effect_lookup(
        set(joined["model_id"].astype(str)),
        set(joined["target_gene"].astype(str)),
    )
    joined["gene_effect"] = [
        lookup.get((str(model_id), str(gene)))
        for model_id, gene in zip(joined["model_id"], joined["target_gene"])
    ]
    joined["crispr_support_level"] = joined["gene_effect"].map(_crispr_level)
    keep = [c for c in columns if c in joined.columns]
    return joined[keep].drop_duplicates().reset_index(drop=True)


def _write_templates(input_dir: Path) -> None:
    template_dir = input_dir / "manual_template"
    template_dir.mkdir(parents=True, exist_ok=True)
    for name, columns in {
        "compounds.tsv": COMPOUND_COLUMNS,
        "responses.tsv": RESPONSE_COLUMNS,
        "target_evidence.tsv": TARGET_COLUMNS,
    }.items():
        path = template_dir / name
        if not path.exists():
            pd.DataFrame(columns=columns).to_csv(path, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build MCL Pharmacology Layer v1 by merging normalized source subdirectories."
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument(
        "--create-templates",
        action="store_true",
        help="Create empty manual templates under normalized/manual_template.",
    )
    args = parser.parse_args()
    input_dir = args.input_dir.resolve()
    if args.create_templates:
        _write_templates(input_dir)

    if not ATLAS.exists():
        raise SystemExit("Missing data/processed/depmap_crispr_model_atlas.parquet")
    atlas = pd.read_parquet(ATLAS)
    valid_models = set(atlas["model_id"].dropna().astype(str))

    raw_compounds = _read_many(input_dir, "compounds", COMPOUND_COLUMNS)
    raw_responses = _read_many(input_dir, "responses", RESPONSE_COLUMNS)
    raw_targets = _read_many(input_dir, "target_evidence", TARGET_COLUMNS)
    compounds = _normalise_compounds(raw_compounds)
    valid_compounds = set(compounds["compound_id"].astype(str)) if not compounds.empty else set()
    responses, response_qc = _normalise_responses(raw_responses, valid_models, valid_compounds)
    targets, target_qc = _normalise_targets(raw_targets, valid_compounds)
    links = _build_links(responses, targets)

    RUNTIME.mkdir(parents=True, exist_ok=True)
    QC.mkdir(parents=True, exist_ok=True)
    compounds.to_parquet(RUNTIME / "compounds.parquet", index=False, compression="zstd")
    responses.to_parquet(RUNTIME / "responses.parquet", index=False, compression="zstd")
    targets.to_parquet(RUNTIME / "target_evidence.parquet", index=False, compression="zstd")
    links.to_parquet(RUNTIME / "model_compound_target_links.parquet", index=False, compression="zstd")
    response_qc.to_csv(QC / "pharmacology_unresolved_responses.tsv", sep="\t", index=False)
    target_qc.to_csv(QC / "pharmacology_unresolved_targets.tsv", sep="\t", index=False)

    sources = sorted(set(responses["source"].dropna().astype(str))) if not responses.empty else []
    source_dirs = sorted(
        set(
            str(path.parent.relative_to(input_dir))
            for stem in ("compounds", "responses", "target_evidence")
            for path in _find_all(input_dir, stem)
        )
    )
    manifest = {
        "contract": "mcl-pharmacology-v1",
        "schema_version": "1.0",
        "status": "available" if not responses.empty else "empty_inputs",
        "built_at": datetime.now(timezone.utc).isoformat(),
        "input_dir": str(input_dir),
        "normalized_source_dirs": source_dirs,
        "models_in_crispr_atlas_n": len(valid_models),
        "compounds_n": int(compounds["compound_id"].nunique()) if not compounds.empty else 0,
        "responses_n": int(len(responses)),
        "models_with_response_n": int(responses["model_id"].nunique()) if not responses.empty else 0,
        "target_evidence_n": int(len(targets)),
        "compound_target_pairs_n": int(
            targets[["compound_id", "target_gene"]].drop_duplicates().shape[0]
        ) if not targets.empty else 0,
        "model_compound_target_links_n": int(len(links)),
        "unresolved_responses_n": int(len(response_qc)),
        "unresolved_targets_n": int(len(target_qc)),
        "sources": sources,
        "build_command": ".\\.venv\\Scripts\\python.exe .\\scripts\\build_pharmacology_layer.py",
        "note_ru": (
            "IC50, AUC, GI50, viability, PRISM LFC и другие endpoints сохраняются раздельно "
            "и не сводятся автоматически к общей шкале."
        ),
    }
    (RUNTIME / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print("MCL Pharmacology Layer v1")
    print(f"Normalized source directories: {', '.join(source_dirs) if source_dirs else 'none'}")
    print(f"CRISPR Atlas models: {len(valid_models)}")
    print(f"Compounds: {manifest['compounds_n']}")
    print(f"Response observations: {manifest['responses_n']}")
    print(f"Models with pharmacology: {manifest['models_with_response_n']}")
    print(f"Target evidence rows: {manifest['target_evidence_n']}")
    print(f"Model × compound × target links: {manifest['model_compound_target_links_n']}")
    print(f"Unresolved responses: {manifest['unresolved_responses_n']}")
    print(f"Unresolved targets: {manifest['unresolved_targets_n']}")
    print(f"Wrote {RUNTIME.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
