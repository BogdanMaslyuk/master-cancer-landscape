from __future__ import annotations

import json
from functools import lru_cache
from typing import Any

import numpy as np
import pandas as pd

from .runtime_gene_explorer import RuntimeGeneExplorerStore
from .store import MCLDataError


class FastRuntimeGeneExplorerStore(RuntimeGeneExplorerStore):
    """Runtime Gene Explorer with required memory-mapped gene-by-model arrays.

    The ordinary processed Parquet matrices are intentionally wide (models x genes).
    Opening an arbitrary gene column from a ~20k-column Parquet file can spend tens
    of seconds parsing wide-file metadata on Windows. Explorer therefore requires a
    build-time transposed NumPy runtime array (genes x models). One gene lookup then
    reads a single contiguous row from a memory-mapped file.

    Architecture v1 deliberately fails closed when an optimized runtime layer is
    missing or invalid. Interactive HTTP requests must never fall back to parsing the
    wide processed scientific matrices. Rebuild serving artifacts with
    scripts/build-explorer.ps1 instead.
    """

    @lru_cache(maxsize=3)
    def _runtime_layer_bundle(
        self, layer: str
    ) -> tuple[np.ndarray, tuple[str, ...], dict[str, int]]:
        metadata_path = self._runtime_path(f"model_layers/{layer}.json")
        matrix_path = self._runtime_path(f"model_layers/{layer}.npy")

        try:
            metadata: dict[str, Any] = json.loads(metadata_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError) as exc:
            raise MCLDataError(
                f"Gene Explorer runtime metadata is invalid for layer '{layer}'. "
                "Run scripts/build-explorer.ps1 to rebuild serving artifacts."
            ) from exc

        models = tuple(str(value) for value in (metadata.get("models") or []))
        genes = [str(value).upper() for value in (metadata.get("genes") or [])]
        if not models or not genes:
            raise MCLDataError(
                f"Gene Explorer runtime metadata is incomplete for layer '{layer}'. "
                "Run scripts/build-explorer.ps1 to rebuild serving artifacts."
            )

        try:
            matrix = np.load(matrix_path, mmap_mode="r", allow_pickle=False)
        except (OSError, ValueError) as exc:
            raise MCLDataError(
                f"Gene Explorer runtime matrix is invalid for layer '{layer}'. "
                "Run scripts/build-explorer.ps1 to rebuild serving artifacts."
            ) from exc

        if matrix.ndim != 2 or matrix.shape != (len(genes), len(models)):
            raise MCLDataError(
                f"Gene Explorer runtime matrix shape does not match metadata for layer '{layer}'. "
                "Run scripts/build-explorer.ps1 to rebuild serving artifacts."
            )

        gene_index = {symbol: index for index, symbol in enumerate(genes)}
        return matrix, models, gene_index

    @lru_cache(maxsize=1536)
    def _model_gene_layer(self, layer: str, gene_symbol: str) -> pd.DataFrame:
        symbol = str(gene_symbol).strip().upper()
        matrix, models, gene_index = self._runtime_layer_bundle(layer)
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
