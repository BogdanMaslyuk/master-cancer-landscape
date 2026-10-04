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


def _read(path: Path, columns: list[str]) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=columns)
    if path.suffix.lower() == ".parquet":
        frame = pd.read_parquet(path)
    else:
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        frame = pd.read_csv(path, sep=sep, low_memory=False)
    for column in columns:
        if column not in frame.columns:
            frame[column] = None
    return frame[columns].copy()


def _find(input_dir: Path, stem: str) -> Path | None:
    for suffix in (".parquet", ".tsv", ".csv"):
        path = input_dir / f"{stem}{suffix}"
        if path.exists():
            return path
    return None


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _normalise_compounds(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["compound_id"] = out["compound_id"].map(_text)
    out = out[out["compound_id"] != ""].copy()
    if out.empty:
        return out
    out = out.drop_duplicates("compound_id", keep="first")
    return out.reset_index(drop=True)


def _normalise_responses(frame: pd.DataFrame, valid_models: set[str], valid_compounds: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = frame.copy()
    for col in ("observation_id", "model_id", "compound_id", "source", "endpoint", "unit"):
        out[col] = out[col].map(_text)
    if "observation_id" in out.columns:
        missing = out["observation_id"].eq("")
        out.loc[missing, "observation_id"] = [f"PHARM-OBS-{i:08d}" for i in range(1, int(missing.sum()) + 1)]
    for col in ("value", "auc", "ic50", "ec50", "gi50", "viability", "dose", "exposure_time_h", "replicate_n"):
        out[col] = pd.to_numeric(out[col], errors="coerce")
    reason = pd.Series("", index=out.index, dtype=object)
    reason.loc[~out["model_id"].isin(valid_models)] = "unknown_model_id"
    reason.loc[~out["compound_id"].isin(valid_compounds)] = reason.loc[~out["compound_id"].isin(valid_compounds)].map(
        lambda x: f"{x};unknown_compound_id".strip(";")
    )
    unresolved = out[reason.ne("")].copy()
    unresolved["qc_reason"] = reason[reason.ne("")]
    resolved = out[reason.eq("")].copy()
    resolved = resolved.drop_duplicates("observation_id", keep="first")
    return resolved.reset_index(drop=True), unresolved.reset_index(drop=True)


def _normalise_targets(frame: pd.DataFrame, valid_compounds: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = frame.copy()
    for col in ("evidence_id", "compound_id", "target_gene", "action", "evidence_type", "source", "confidence", "directness"):
        out[col] = out[col].map(_text)
    out["target_gene"] = out["target_gene"].str.upper()
    missing = out["evidence_id"].eq("")
    out.loc[missing, "evidence_id"] = [f"PHARM-TGT-{i:08d}" for i in range(1, int(missing.sum()) + 1)]
    out["activity_value"] = pd.to_numeric(out["activity_value"], errors="coerce")
    reason = pd.Series("", index=out.index, dtype=object)
    reason.loc[~out["compound_id"].isin(valid_compounds)] = "unknown_compound_id"
    reason.loc[out["target_gene"].eq("")] = reason.loc[out["target_gene"].eq("")].map(
        lambda x: f"{x};missing_target_gene".strip(";")
    )
    unresolved = out[reason.ne("")].copy()
    unresolved["qc_reason"] = reason[reason.ne("")]
    resolved = out[reason.eq("")].copy().drop_duplicates("evidence_id", keep="first")
    return resolved.reset_index(drop=True), unresolved.reset_index(drop=True)


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
    lookup = _gene_effect_lookup(set(joined["model_id"].astype(str)), set(joined["target_gene"].astype(str)))
    joined["gene_effect"] = [lookup.get((str(m), str(g))) for m, g in zip(joined["model_id"], joined["target_gene"])]
    joined["crispr_support_level"] = joined["gene_effect"].map(_crispr_level)
    keep = [c for c in columns if c in joined.columns]
    return joined[keep].drop_duplicates().reset_index(drop=True)


def _write_templates(input_dir: Path) -> None:
    input_dir.mkdir(parents=True, exist_ok=True)
    templates = {
        "compounds.tsv": COMPOUND_COLUMNS,
        "responses.tsv": RESPONSE_COLUMNS,
        "target_evidence.tsv": TARGET_COLUMNS,
    }
    for name, columns in templates.items():
        path = input_dir / name
        if not path.exists():
            pd.DataFrame(columns=columns).to_csv(path, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build MCL Pharmacology Layer v1 from normalized source tables.")
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--create-templates", action="store_true", help="Create empty normalized TSV templates when files are missing.")
    args = parser.parse_args()
    input_dir = args.input_dir.resolve()
    if args.create_templates:
        _write_templates(input_dir)

    if not ATLAS.exists():
        raise SystemExit("Missing data/processed/depmap_crispr_model_atlas.parquet")
    atlas = pd.read_parquet(ATLAS)
    valid_models = set(atlas["model_id"].dropna().astype(str))

    compound_path = _find(input_dir, "compounds")
    response_path = _find(input_dir, "responses")
    target_path = _find(input_dir, "target_evidence")
    compounds = _normalise_compounds(_read(compound_path, COMPOUND_COLUMNS) if compound_path else pd.DataFrame(columns=COMPOUND_COLUMNS))
    valid_compounds = set(compounds["compound_id"].astype(str)) if not compounds.empty else set()
    responses, response_qc = _normalise_responses(
        _read(response_path, RESPONSE_COLUMNS) if response_path else pd.DataFrame(columns=RESPONSE_COLUMNS),
        valid_models,
        valid_compounds,
    )
    targets, target_qc = _normalise_targets(
        _read(target_path, TARGET_COLUMNS) if target_path else pd.DataFrame(columns=TARGET_COLUMNS),
        valid_compounds,
    )
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
    manifest = {
        "contract": "mcl-pharmacology-v1",
        "schema_version": "1.0",
        "status": "available" if not responses.empty else "empty_inputs",
        "built_at": datetime.now(timezone.utc).isoformat(),
        "input_dir": str(input_dir),
        "models_in_crispr_atlas_n": len(valid_models),
        "compounds_n": int(compounds["compound_id"].nunique()) if not compounds.empty else 0,
        "responses_n": int(len(responses)),
        "models_with_response_n": int(responses["model_id"].nunique()) if not responses.empty else 0,
        "target_evidence_n": int(len(targets)),
        "compound_target_pairs_n": int(targets[["compound_id", "target_gene"]].drop_duplicates().shape[0]) if not targets.empty else 0,
        "model_compound_target_links_n": int(len(links)),
        "unresolved_responses_n": int(len(response_qc)),
        "unresolved_targets_n": int(len(target_qc)),
        "sources": sources,
        "build_command": ".\\.venv\\Scripts\\python.exe .\\scripts\\build_pharmacology_layer.py",
        "note_ru": "IC50, AUC, GI50, viability и другие endpoints сохраняются раздельно и не сводятся автоматически к общей шкале.",
    }
    (RUNTIME / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print("MCL Pharmacology Layer v1")
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
