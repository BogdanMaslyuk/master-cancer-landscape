from __future__ import annotations

import json
from functools import lru_cache
from typing import Any

import numpy as np
import pandas as pd

from .runtime_gene_explorer import RuntimeGeneExplorerStore


class FastRuntimeGeneExplorerStore(RuntimeGeneExplorerStore):
    """Runtime Gene Explorer with memory-mapped gene-by-model multi-omics arrays.

    The ordinary processed Parquet matrices are intentionally wide (models x genes).
    Opening an arbitrary gene column from a ~20k-column Parquet file can spend tens
    of seconds parsing wide-file metadata on Windows. Explorer therefore prefers a
    build-time transposed NumPy runtime array (genes x models). One gene lookup then
    reads a single contiguous row from a memory-mapped file.

    When the optimized runtime array is absent we retain the parent implementation
    as a compatibility fallback, but normal local builds should always create the
    arrays through scripts/build_depmap_runtime_arrays.py.
    """

    @lru_cache(maxsize=3)
    def _runtime_layer_bundle(
        self, layer: str
    ) -> tuple[np.ndarray, tuple[str, ...], dict[str, int]] | None:
        runtime_dir = self.index_dir / "model_layers"
        metadata_path = runtime_dir / f"{layer}.json"
        matrix_path = runtime_dir / f"{layer}.npy"
        if not metadata_path.exists() or not matrix_path.exists():
            return None

        try:
            metadata: dict[str, Any] = json.loads(metadata_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            return None

        models = tuple(str(value) for value in (metadata.get("models") or []))
        genes = [str(value).upper() for value in (metadata.get("genes") or [])]
        if not models or not genes:
            return None

        try:
            matrix = np.load(matrix_path, mmap_mode="r", allow_pickle=False)
        except (OSError, ValueError):
            return None

        if matrix.ndim != 2 or matrix.shape != (len(genes), len(models)):
            return None

        gene_index = {symbol: index for index, symbol in enumerate(genes)}
        return matrix, models, gene_index

    @lru_cache(maxsize=1536)
    def _model_gene_layer(self, layer: str, gene_symbol: str) -> pd.DataFrame:
        symbol = str(gene_symbol).strip().upper()
        bundle = self._runtime_layer_bundle(layer)
        if bundle is None:
            return super()._model_gene_layer(layer, symbol)

        matrix, models, gene_index = bundle
        index = gene_index.get(symbol)
        if index is None:
            return pd.DataFrame(columns=["model_id", layer])

        # The runtime matrix is gene x model, so this is a contiguous row read.
        values = np.asarray(matrix[index, :], dtype=np.float32)
        return pd.DataFrame(
            {
                "model_id": list(models),
                layer: pd.to_numeric(values, errors="coerce"),
            }
        )
