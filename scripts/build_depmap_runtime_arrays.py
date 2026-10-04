from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
RUNTIME_DIR = ROOT / "data" / "runtime" / "explorer" / "model_layers"
LAYERS = ("gene_effect", "expression", "copy_number")


def build_layer(layer: str) -> None:
    source = PROCESSED / f"depmap_model_{layer}.parquet"
    if not source.exists():
        raise SystemExit(f"Missing processed multi-omics matrix: {source}")

    print(f"[{layer}] reading {source.name} ...")
    frame = pd.read_parquet(source)
    if "model_id" not in frame.columns:
        raise SystemExit(f"{source.name} does not contain model_id")

    models = frame["model_id"].astype(str).tolist()
    genes = [str(column).upper() for column in frame.columns if column != "model_id"]
    if not models or not genes:
        raise SystemExit(f"{source.name} has no model/gene data")

    values = frame.drop(columns=["model_id"]).to_numpy(dtype=np.float32, copy=False)
    if values.shape != (len(models), len(genes)):
        raise SystemExit(
            f"Unexpected matrix shape for {layer}: {values.shape}; "
            f"expected {(len(models), len(genes))}"
        )

    # Transpose once at build time. Each runtime gene lookup then reads one
    # contiguous row instead of reopening a ~20k-column Parquet schema.
    gene_by_model = np.ascontiguousarray(values.T, dtype=np.float32)

    RUNTIME_DIR.mkdir(parents=True, exist_ok=True)
    matrix_path = RUNTIME_DIR / f"{layer}.npy"
    metadata_path = RUNTIME_DIR / f"{layer}.json"

    np.save(matrix_path, gene_by_model, allow_pickle=False)
    metadata = {
        "schema_version": "1.0",
        "layer": layer,
        "orientation": "gene_by_model",
        "dtype": "float32",
        "models_n": len(models),
        "genes_n": len(genes),
        "models": models,
        "genes": genes,
        "source_file": source.name,
        "matrix_file": matrix_path.name,
    }
    metadata_path.write_text(json.dumps(metadata, ensure_ascii=False), encoding="utf-8")

    size_mb = matrix_path.stat().st_size / (1024 * 1024)
    print(
        f"[{layer}] wrote {len(genes)} genes x {len(models)} models "
        f"-> {matrix_path.relative_to(ROOT)} ({size_mb:.1f} MB)"
    )


def main() -> None:
    print("Building fast DepMap runtime arrays (gene x model)...")
    for layer in LAYERS:
        build_layer(layer)
    print("Fast DepMap runtime arrays: PASS")


if __name__ == "__main__":
    main()
