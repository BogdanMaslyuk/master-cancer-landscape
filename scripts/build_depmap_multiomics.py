from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
RAW_DEPMAP = ROOT / "data" / "raw" / "depmap"

LAYER_SOURCES: dict[str, tuple[str, ...]] = {
    "expression": (
        "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
        "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
    ),
    "copy_number": ("OmicsCNGeneWGS.csv", "OmicsCNGene.csv"),
    "gene_effect": ("CRISPRGeneEffect.csv",),
}

OUTPUTS = {
    "expression": "depmap_model_expression.parquet",
    "copy_number": "depmap_model_copy_number.parquet",
    "gene_effect": "depmap_model_gene_effect.parquet",
}

LAYER_META = {
    "expression": {
        "label": "RNA expression",
        "value_semantics": "log2(TPM + 1)",
        "note": "DepMap model-level protein-coding gene expression matrix.",
    },
    "copy_number": {
        "label": "Relative copy number",
        "value_semantics": "relative copy number, linear scale",
        "note": "MCL stores the DepMap relative WGS copy-number values as supplied; it does not convert them to absolute integer copy number or call amplification/deletion here.",
    },
    "gene_effect": {
        "label": "CRISPR Gene Effect",
        "value_semantics": "Chronos Gene Effect; more negative values indicate stronger loss-of-function dependency",
        "note": "CRISPR knockout dependency is not equivalent to pharmacological inhibition.",
    },
}

MODEL_ID_ALIASES = ("ModelID", "DepMap_ID", "DepMapID")
GENE_LABEL_RE = re.compile(r"^(.*?)\s*\((\d+)\)\s*$")


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _release_from_inventory() -> str | None:
    path = PROCESSED / "depmap_input_inventory.tsv"
    if not path.exists():
        return None
    try:
        frame = pd.read_csv(path, sep="\t", usecols=["depmap_release"], nrows=1)
    except (ValueError, OSError):
        return None
    if frame.empty:
        return None
    value = str(frame.iloc[0]["depmap_release"]).strip()
    return value or None


def _resolve_release(explicit: str | None) -> str:
    if explicit:
        return explicit
    release = _release_from_inventory()
    if release:
        return release
    if RAW_DEPMAP.exists():
        candidates = sorted((p.name for p in RAW_DEPMAP.iterdir() if p.is_dir()), reverse=True)
        if candidates:
            return candidates[0]
    raise SystemExit("Cannot infer DepMap release. Pass --release, e.g. --release 26Q1")


def _load_model_universe() -> tuple[set[str], str]:
    """Prefer the full CRISPR Atlas; retain the older Wave 1 audit as fallback."""
    atlas_parquet = PROCESSED / "depmap_crispr_model_atlas.parquet"
    atlas_tsv = PROCESSED / "depmap_crispr_model_atlas.tsv"
    audit_parquet = PROCESSED / "depmap_context_audit.parquet"
    audit_tsv = PROCESSED / "depmap_context_audit.tsv"

    if atlas_parquet.exists():
        frame = pd.read_parquet(atlas_parquet, columns=["model_id"])
        source = atlas_parquet.name
    elif atlas_tsv.exists():
        frame = pd.read_csv(atlas_tsv, sep="\t", usecols=["model_id"], low_memory=False)
        source = atlas_tsv.name
    elif audit_parquet.exists():
        frame = pd.read_parquet(audit_parquet, columns=["model_id"])
        source = audit_parquet.name
    elif audit_tsv.exists():
        frame = pd.read_csv(audit_tsv, sep="\t", usecols=["model_id"], low_memory=False)
        source = audit_tsv.name
    else:
        raise SystemExit(
            "Missing CRISPR model universe. Build scripts/build_crispr_cancer_atlas.py first "
            "or provide the legacy depmap_context_audit.tsv/parquet."
        )

    ids = set(frame["model_id"].dropna().astype(str).str.strip())
    ids.discard("")
    return ids, source


def _resolve_source(release_dir: Path, candidates: Iterable[str]) -> Path | None:
    return next((release_dir / name for name in candidates if (release_dir / name).exists()), None)


def _normalize_gene_label(label: str) -> tuple[str, str | None]:
    text = str(label).strip()
    match = GENE_LABEL_RE.match(text)
    if match:
        return match.group(1).strip(), match.group(2)
    return text, None


def _column_plan(columns: list[str]) -> tuple[str, list[str], dict[str, str], list[dict[str, str | None]]]:
    if not columns:
        raise ValueError("Matrix has an empty CSV header")
    id_col = next((c for c in MODEL_ID_ALIASES if c in columns), columns[0])
    gene_cols = [c for c in columns if c != id_col]
    rename: dict[str, str] = {id_col: "model_id"}
    gene_map: list[dict[str, str | None]] = []
    seen: set[str] = set()
    keep: list[str] = [id_col]
    duplicate_symbols: set[str] = set()

    for raw in gene_cols:
        symbol, entrez = _normalize_gene_label(raw)
        if not symbol:
            continue
        if symbol in seen:
            duplicate_symbols.add(symbol)
            continue
        seen.add(symbol)
        keep.append(raw)
        rename[raw] = symbol
        gene_map.append({"gene": symbol, "entrez_id": entrez, "source_column": raw})

    if duplicate_symbols:
        print(f"Warning: ignored duplicate normalized gene symbols: {sorted(duplicate_symbols)[:20]}")
    return id_col, keep, rename, gene_map


def _build_matrix(path: Path, model_ids: set[str], chunksize: int) -> tuple[pd.DataFrame, list[dict[str, str | None]]]:
    header = pd.read_csv(path, nrows=0).columns.tolist()
    id_col, keep, rename, gene_map = _column_plan(header)
    chunks: list[pd.DataFrame] = []

    # Chunking by rows keeps peak memory bounded even for ~20k-gene matrices.
    for chunk in pd.read_csv(path, usecols=keep, chunksize=max(8, chunksize), low_memory=False):
        mask = chunk[id_col].astype(str).str.strip().isin(model_ids)
        if not mask.any():
            continue
        chunks.append(chunk.loc[mask].copy())

    if not chunks:
        return pd.DataFrame(columns=["model_id"] + [x["gene"] for x in gene_map]), gene_map

    frame = pd.concat(chunks, ignore_index=True).rename(columns=rename)
    frame["model_id"] = frame["model_id"].astype(str).str.strip()
    frame = frame.drop_duplicates("model_id", keep="first")

    gene_cols = [c for c in frame.columns if c != "model_id"]
    non_numeric = [c for c in gene_cols if not pd.api.types.is_numeric_dtype(frame[c])]
    for col in non_numeric:
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    if gene_cols:
        frame[gene_cols] = frame[gene_cols].astype("float32")

    return frame.sort_values("model_id").reset_index(drop=True), gene_map


def _write_layer(layer: str, frame: pd.DataFrame) -> Path:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    path = PROCESSED / OUTPUTS[layer]
    frame.to_parquet(path, index=False, compression="zstd")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build model-level DepMap expression, relative copy-number and CRISPR Gene Effect indexes for the full MCL CRISPR model universe."
    )
    parser.add_argument("--release", default=None, help="DepMap release directory, e.g. 26Q1. Inferred by default.")
    parser.add_argument("--chunksize", type=int, default=48, help="Rows per wide-matrix CSV chunk. Lower this if memory is limited.")
    parser.add_argument(
        "--allow-partial",
        action="store_true",
        help="Build available layers even when expression or copy-number files are missing.",
    )
    args = parser.parse_args()

    release = _resolve_release(args.release)
    release_dir = RAW_DEPMAP / release
    model_ids, universe_source = _load_model_universe()
    print(f"DepMap release: {release}")
    print(f"MCL model universe: {len(model_ids)} models from {universe_source}")

    manifest: dict[str, object] = {
        "depmap_release": release,
        "built_at": _utc_now(),
        "model_universe_n": len(model_ids),
        "model_universe_source": universe_source,
        "raw_dir": str(release_dir),
        "layers": {},
    }
    missing: list[str] = []

    for layer, candidates in LAYER_SOURCES.items():
        source = _resolve_source(release_dir, candidates)
        if source is None:
            missing.append(" or ".join(candidates))
            manifest["layers"][layer] = {
                **LAYER_META[layer],
                "available": False,
                "source_file": None,
                "output_file": OUTPUTS[layer],
                "models_n": 0,
                "genes_n": 0,
                "missing_models_n": len(model_ids),
            }
            print(f"[{layer}] missing source: {' or '.join(candidates)}")
            continue

        print(f"[{layer}] reading {source.name} ...")
        frame, gene_map = _build_matrix(source, model_ids, args.chunksize)
        output = _write_layer(layer, frame)
        indexed_ids = set(frame["model_id"].astype(str)) if "model_id" in frame.columns else set()
        map_path = PROCESSED / f"depmap_model_{layer}_genes.json"
        map_path.write_text(json.dumps(gene_map, ensure_ascii=False), encoding="utf-8")
        manifest["layers"][layer] = {
            **LAYER_META[layer],
            "available": True,
            "source_file": source.name,
            "output_file": output.name,
            "gene_map_file": map_path.name,
            "models_n": len(indexed_ids),
            "genes_n": max(0, len(frame.columns) - 1),
            "missing_models_n": len(model_ids - indexed_ids),
        }
        print(f"[{layer}] indexed {len(indexed_ids)} models × {max(0, len(frame.columns) - 1)} genes -> {output.name}")

    manifest_path = PROCESSED / "depmap_model_multiomics_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote {manifest_path.relative_to(ROOT)}")

    if missing and not args.allow_partial:
        expected = "\n".join(f"  - {release_dir / name.split(' or ')[0]}" for name in missing)
        raise SystemExit(
            "Multi-omics index is incomplete because source files are missing.\n"
            "Download the missing files from the same pinned DepMap release and rerun:\n"
            f"{expected}\n"
            "CRISPRGeneEffect may already be present. Use --allow-partial only if you intentionally want the available layers."
        )


if __name__ == "__main__":
    main()
