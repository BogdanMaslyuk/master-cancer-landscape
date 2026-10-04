from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .store import MCLDataError


LAYER_FILES = {
    "expression": "depmap_model_expression.parquet",
    "copy_number": "depmap_model_copy_number.parquet",
    "gene_effect": "depmap_model_gene_effect.parquet",
}

LAYER_LABELS = {
    "expression": "RNA expression",
    "copy_number": "Relative copy number",
    "gene_effect": "CRISPR Gene Effect",
}


def _clean(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean(v) for v in value]
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(value, "item"):
        try:
            return value.item()
        except (TypeError, ValueError):
            pass
    return value


class MCLMultiOmicsStore:
    """Read compact model-level DepMap multi-omics indexes built for MCL Explorer."""

    def __init__(self, root: Path):
        self.root = Path(root).resolve()
        self.processed = self.root / "data" / "processed"

    @lru_cache(maxsize=1)
    def _manifest(self) -> dict[str, Any]:
        path = self.processed / "depmap_model_multiomics_manifest.json"
        if not path.exists():
            return {}
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            return {}

    @lru_cache(maxsize=3)
    def _matrix(self, layer: str) -> pd.DataFrame:
        if layer not in LAYER_FILES:
            raise MCLDataError(f"Unknown multi-omics layer: {layer}")
        path = self.processed / LAYER_FILES[layer]
        if not path.exists():
            return pd.DataFrame()
        frame = pd.read_parquet(path)
        if frame.empty or "model_id" not in frame.columns:
            return pd.DataFrame()
        frame["model_id"] = frame["model_id"].astype(str)
        return frame.set_index("model_id", drop=True)

    @lru_cache(maxsize=1)
    def _audit(self) -> pd.DataFrame:
        parquet = self.processed / "depmap_context_audit.parquet"
        tsv = self.processed / "depmap_context_audit.tsv"
        if parquet.exists():
            return pd.read_parquet(parquet)
        if tsv.exists():
            return pd.read_csv(tsv, sep="\t", low_memory=False)
        return pd.DataFrame()

    @lru_cache(maxsize=1)
    def _stable_genes(self) -> list[str]:
        path = self.root / "outputs" / "reports" / "M3_3_1_stable_recurrent_genes.tsv"
        if not path.exists():
            return []
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        col = "gene_symbol" if "gene_symbol" in frame.columns else ("gene" if "gene" in frame.columns else None)
        if col is None:
            return []
        return list(dict.fromkeys(frame[col].dropna().astype(str).str.upper().tolist()))

    def availability(self) -> dict[str, Any]:
        manifest = self._manifest()
        manifest_layers = manifest.get("layers") if isinstance(manifest, dict) else {}
        layers: dict[str, Any] = {}
        for layer, filename in LAYER_FILES.items():
            path = self.processed / filename
            meta = (manifest_layers or {}).get(layer) or {}
            layers[layer] = {
                "id": layer,
                "label": LAYER_LABELS[layer],
                "available": path.exists(),
                "models_n": meta.get("models_n"),
                "genes_n": meta.get("genes_n"),
                "source_file": meta.get("source_file"),
                "value_semantics": meta.get("value_semantics"),
                "note": meta.get("note"),
            }
        return {
            "available": any(x["available"] for x in layers.values()),
            "complete": all(x["available"] for x in layers.values()),
            "depmap_release": manifest.get("depmap_release") if isinstance(manifest, dict) else None,
            "built_at": manifest.get("built_at") if isinstance(manifest, dict) else None,
            "layers": layers,
            "build_command": ".\\.venv\\Scripts\\python.exe .\\scripts\\build_depmap_multiomics.py",
        }

    def _row(self, layer: str, model_id: str) -> pd.Series | None:
        frame = self._matrix(layer)
        if frame.empty:
            return None
        model_id = model_id.strip()
        if model_id not in frame.index:
            matches = [idx for idx in frame.index if str(idx).upper() == model_id.upper()]
            if not matches:
                return None
            model_id = matches[0]
        row = frame.loc[model_id]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        return row

    @staticmethod
    def _value(row: pd.Series | None, gene: str) -> float | None:
        if row is None:
            return None
        symbol = gene.upper()
        if symbol not in row.index:
            return None
        value = pd.to_numeric(pd.Series([row[symbol]]), errors="coerce").iloc[0]
        return None if pd.isna(value) else float(value)

    def _candidate_genes(self, extra: list[str] | None = None) -> list[str]:
        ordered: list[str] = []
        for gene in (extra or []) + self._stable_genes():
            symbol = str(gene).strip().upper()
            if symbol and symbol not in ordered:
                ordered.append(symbol)
        return ordered

    def model(self, model_id: str, genes: list[str] | None = None, limit: int = 30) -> dict[str, Any]:
        availability = self.availability()
        rows = {layer: self._row(layer, model_id) for layer in LAYER_FILES}
        if not any(row is not None for row in rows.values()):
            return {
                "model_id": model_id,
                "availability": availability,
                "candidate_panel": [],
                "top_dependencies": [],
                "note": "Для этой модели multi-omics индекс пока не построен или модель отсутствует в доступных матрицах.",
            }

        candidate_genes = self._candidate_genes(genes)[:100]
        candidate_panel = []
        for gene in candidate_genes:
            values = {
                "gene": gene,
                "gene_effect": self._value(rows["gene_effect"], gene),
                "expression": self._value(rows["expression"], gene),
                "copy_number": self._value(rows["copy_number"], gene),
            }
            if any(values[key] is not None for key in ("gene_effect", "expression", "copy_number")):
                candidate_panel.append(values)

        top_dependencies: list[dict[str, Any]] = []
        ge = rows["gene_effect"]
        if ge is not None:
            numeric = pd.to_numeric(ge, errors="coerce").dropna().sort_values(ascending=True)
            for gene, value in numeric.head(max(1, min(int(limit), 100))).items():
                symbol = str(gene).upper()
                top_dependencies.append(
                    {
                        "gene": symbol,
                        "gene_effect": float(value),
                        "expression": self._value(rows["expression"], symbol),
                        "copy_number": self._value(rows["copy_number"], symbol),
                    }
                )

        layer_model_status = {
            layer: {
                "available": availability["layers"][layer]["available"],
                "model_present": rows[layer] is not None,
            }
            for layer in LAYER_FILES
        }

        return _clean(
            {
                "model_id": model_id,
                "availability": availability,
                "model_layers": layer_model_status,
                "candidate_panel": candidate_panel,
                "top_dependencies": top_dependencies,
                "interpretation": {
                    "gene_effect": "Более отрицательный Gene Effect означает более сильную зависимость клетки от функции гена после CRISPR-выключения; это не равнозначно фармакологическому ингибированию.",
                    "expression": "Expression показан в шкале DepMap log2(TPM + 1) и используется как характеристика экспрессии модели, а не как доказательство активности белка.",
                    "copy_number": "Copy number — относительное линейное значение DepMap WGS. MCL v0.1 не превращает его автоматически в абсолютное число копий или клинический вызов amplification/deletion.",
                },
            }
        )

    @staticmethod
    def _median(frame: pd.DataFrame, model_ids: list[str], gene: str) -> float | None:
        if frame.empty or gene not in frame.columns:
            return None
        ids = [x for x in model_ids if x in frame.index]
        if not ids:
            return None
        values = pd.to_numeric(frame.loc[ids, gene], errors="coerce").dropna()
        return None if values.empty else float(values.median())

    def context(self, cancer_id: str, genes: list[str] | None = None, limit: int = 12) -> dict[str, Any]:
        audit = self._audit()
        if audit.empty or "cancer_id" not in audit.columns or "model_id" not in audit.columns:
            raise MCLDataError("DepMap context audit is not available")
        sub = audit[audit["cancer_id"].astype(str) == cancer_id].copy()
        if sub.empty:
            raise MCLDataError(f"Unknown cancer context: {cancer_id}")

        sub["model_id"] = sub["model_id"].astype(str)
        context_ids = sub.loc[sub.get("assigned_group", "").astype(str) == "context", "model_id"].drop_duplicates().tolist()
        comparator_ids = sub.loc[sub.get("assigned_group", "").astype(str) == "comparator", "model_id"].drop_duplicates().tolist()
        all_ids = sub["model_id"].drop_duplicates().tolist()
        availability = self.availability()

        coverage: dict[str, Any] = {}
        matrices: dict[str, pd.DataFrame] = {}
        for layer in LAYER_FILES:
            frame = self._matrix(layer)
            matrices[layer] = frame
            indexed = set(frame.index.astype(str)) if not frame.empty else set()
            context_n = len(set(context_ids) & indexed)
            comparator_n = len(set(comparator_ids) & indexed)
            total_n = len(set(all_ids) & indexed)
            coverage[layer] = {
                "available": not frame.empty,
                "models_n": total_n,
                "models_total_n": len(all_ids),
                "context_models_n": context_n,
                "context_total_n": len(context_ids),
                "comparator_models_n": comparator_n,
                "comparator_total_n": len(comparator_ids),
            }

        candidate_genes = self._candidate_genes(genes)[: max(1, min(int(limit), 50))]
        panel: list[dict[str, Any]] = []
        for gene in candidate_genes:
            row: dict[str, Any] = {"gene": gene}
            for layer, prefix in (("gene_effect", "gene_effect"), ("expression", "expression"), ("copy_number", "copy_number")):
                frame = matrices[layer]
                a = self._median(frame, context_ids, gene)
                b = self._median(frame, comparator_ids, gene)
                row[f"{prefix}_context_median"] = a
                row[f"{prefix}_comparator_median"] = b
                row[f"{prefix}_delta"] = (a - b) if a is not None and b is not None else None
            if any(value is not None for key, value in row.items() if key != "gene"):
                panel.append(row)

        return _clean(
            {
                "cancer_id": cancer_id,
                "availability": availability,
                "coverage": coverage,
                "candidate_panel": panel,
                "note": "Multi-omics значения относятся к экспериментальным моделям DepMap. Они не являются пациентскими частотами или клиническими биомаркерами.",
            }
        )
