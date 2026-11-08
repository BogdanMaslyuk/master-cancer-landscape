from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from .store import MCLDataError


DEPENDENCY_THRESHOLD = -0.5

DEPENDENCY_TYPES: dict[str, dict[str, str]] = {
    "broad_core": {
        "label_ru": "Широкая базовая зависимость",
        "description_ru": "Сильная зависимость встречается у большинства опухолевых моделей DepMap; это не доказывает необходимость гена для всех нормальных клеток.",
    },
    "broad_tumor": {
        "label_ru": "Широкая опухолевая зависимость",
        "description_ru": "Сильная зависимость встречается во многих, но не в большинстве всех моделей.",
    },
    "cancer_enriched": {
        "label_ru": "Обогащена в этом типе опухоли",
        "description_ru": "Зависимость заметно чаще встречается среди моделей того же типа опухоли, чем во всём наборе DepMap.",
    },
    "selective": {
        "label_ru": "Селективная зависимость",
        "description_ru": "Сильная зависимость встречается только в части моделей DepMap.",
    },
    "model_selective": {
        "label_ru": "Особенно выражена в этой модели",
        "description_ru": "В этой модели Gene Effect существенно отрицательнее общего фона, а сильная зависимость редко встречается в других моделях.",
    },
    "weak_or_none": {
        "label_ru": "Слабая / невыраженная",
        "description_ru": "Gene Effect выше рабочего порога MCL −0.5; выраженная зависимость в этой модели не показана.",
    },
}


def _clean(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean(v) for v in value]
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


def _json_list(value: Any) -> list[str]:
    if isinstance(value, list):
        return [str(x) for x in value]
    try:
        payload = json.loads(str(value or "[]"))
    except (json.JSONDecodeError, TypeError):
        return []
    return [str(x) for x in payload] if isinstance(payload, list) else []


class ModelDependencyStore:
    """Interpret one model's genome-wide CRISPR Gene Effect profile.

    The ranking itself is direct DepMap Chronos Gene Effect. Dependency classes are
    transparent MCL heuristics derived from prevalence across the currently indexed
    DepMap cancer models; they are navigation aids, not clinical target labels.
    """

    def __init__(self, root: Path):
        self.root = Path(root).resolve()
        self.processed = self.root / "data" / "processed"

    @lru_cache(maxsize=1)
    def _matrix(self) -> pd.DataFrame:
        path = self.processed / "depmap_model_gene_effect.parquet"
        if not path.exists():
            return pd.DataFrame()
        frame = pd.read_parquet(path)
        if frame.empty or "model_id" not in frame.columns:
            return pd.DataFrame()
        frame["model_id"] = frame["model_id"].astype(str).str.strip()
        return frame.set_index("model_id", drop=True)

    @lru_cache(maxsize=1)
    def _atlas(self) -> pd.DataFrame:
        parquet = self.processed / "depmap_crispr_model_atlas.parquet"
        tsv = self.processed / "depmap_crispr_model_atlas.tsv"
        if parquet.exists():
            frame = pd.read_parquet(parquet)
        elif tsv.exists():
            frame = pd.read_csv(tsv, sep="\t", low_memory=False)
        else:
            return pd.DataFrame()
        if "model_id" in frame.columns:
            frame["model_id"] = frame["model_id"].astype(str).str.strip()
        return frame

    @lru_cache(maxsize=1)
    def _gene_catalog(self) -> pd.DataFrame:
        path = self.processed / "gene_explorer" / "gene_catalog.parquet"
        if not path.exists():
            return pd.DataFrame()
        frame = pd.read_parquet(path)
        if "gene_symbol" in frame.columns:
            frame["gene_symbol"] = frame["gene_symbol"].astype(str).str.upper()
            frame = frame.drop_duplicates("gene_symbol").set_index("gene_symbol", drop=False)
        return frame

    @lru_cache(maxsize=1)
    def _functional_labels(self) -> dict[str, str]:
        path = self.root / "config" / "mcl_functional_domains.yaml"
        if not path.exists():
            return {}
        payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        labels: dict[str, str] = {}
        for domain in payload.get("domains") or []:
            did = str(domain.get("id") or "").strip()
            if did:
                labels[did] = str(domain.get("label_ru") or did)
            for sub in domain.get("subdomains") or []:
                sid = str(sub.get("id") or "").strip()
                if sid:
                    labels[sid] = str(sub.get("label_ru") or sid)
        return labels

    @lru_cache(maxsize=1)
    def _pan_stats(self) -> pd.DataFrame:
        matrix = self._matrix()
        if matrix.empty:
            return pd.DataFrame()
        numeric = matrix.apply(pd.to_numeric, errors="coerce")
        n = numeric.notna().sum(axis=0)
        return pd.DataFrame(
            {
                "pan_median_gene_effect": numeric.median(axis=0, skipna=True),
                "pan_models_n": n,
                "pan_dependency_fraction": numeric.le(DEPENDENCY_THRESHOLD).sum(axis=0).div(n.where(n > 0)),
            }
        )

    @lru_cache(maxsize=128)
    def _cancer_stats(self, cancer_id: str) -> pd.DataFrame:
        matrix = self._matrix()
        atlas = self._atlas()
        if matrix.empty or atlas.empty or "mcl_cancer_id" not in atlas.columns:
            return pd.DataFrame()
        ids = atlas.loc[atlas["mcl_cancer_id"].astype(str) == cancer_id, "model_id"].astype(str)
        ids = [x for x in ids if x in matrix.index]
        if not ids:
            return pd.DataFrame()
        numeric = matrix.loc[ids].apply(pd.to_numeric, errors="coerce")
        n = numeric.notna().sum(axis=0)
        return pd.DataFrame(
            {
                "cancer_median_gene_effect": numeric.median(axis=0, skipna=True),
                "cancer_models_n": n,
                "cancer_dependency_fraction": numeric.le(DEPENDENCY_THRESHOLD).sum(axis=0).div(n.where(n > 0)),
            }
        )

    def _model_meta(self, model_id: str) -> dict[str, Any]:
        atlas = self._atlas()
        if atlas.empty or "model_id" not in atlas.columns:
            return {}
        hit = atlas[atlas["model_id"].astype(str).str.upper() == model_id.upper()]
        return _clean(hit.iloc[0].to_dict()) if not hit.empty else {}

    @staticmethod
    def _dependency_type(
        gene_effect: float,
        pan_fraction: float | None,
        pan_median: float | None,
        cancer_fraction: float | None,
        cancer_n: int,
    ) -> str:
        if gene_effect > DEPENDENCY_THRESHOLD:
            return "weak_or_none"
        pf = float(pan_fraction) if pan_fraction is not None else 0.0
        cf = float(cancer_fraction) if cancer_fraction is not None else 0.0
        if pf >= 0.80:
            return "broad_core"
        if pf >= 0.40:
            return "broad_tumor"
        if cancer_n >= 5 and cf >= 0.50 and cf >= pf + 0.20:
            return "cancer_enriched"
        if pf <= 0.10 and gene_effect <= -0.80 and pan_median is not None and gene_effect <= float(pan_median) - 0.50:
            return "model_selective"
        return "selective"

    def _annotation(self, gene: str) -> dict[str, Any]:
        catalog = self._gene_catalog()
        if catalog.empty or gene not in catalog.index:
            return {"gene_name": None, "domains": [], "subdomains": [], "protein_classes": []}
        row = catalog.loc[gene]
        labels = self._functional_labels()
        domains = _json_list(row.get("mcl_domains_json"))
        subdomains = _json_list(row.get("mcl_subdomains_json"))
        protein_classes = _json_list(row.get("protein_classes_json"))
        return {
            "gene_name": row.get("gene_name"),
            "domains": [{"id": x, "label_ru": labels.get(x, x)} for x in domains],
            "subdomains": [{"id": x, "label_ru": labels.get(x, x)} for x in subdomains],
            "protein_classes": protein_classes,
        }

    def model(
        self,
        model_id: str,
        *,
        search: str | None = None,
        dependency_type: str | None = None,
        domain: str | None = None,
        limit: int = 200,
        offset: int = 0,
    ) -> dict[str, Any]:
        matrix = self._matrix()
        if matrix.empty:
            return {
                "model_id": model_id,
                "available": False,
                "status": "not_indexed",
                "build_command": ".\\.venv\\Scripts\\python.exe .\\scripts\\build_depmap_multiomics.py --allow-partial",
                "note_ru": "Индекс индивидуального CRISPR Gene Effect ещё не построен.",
            }

        matches = [idx for idx in matrix.index if str(idx).upper() == model_id.strip().upper()]
        if not matches:
            return {
                "model_id": model_id,
                "available": False,
                "status": "model_not_indexed",
                "indexed_models_n": int(len(matrix.index)),
                "build_command": ".\\.venv\\Scripts\\python.exe .\\scripts\\build_depmap_multiomics.py --allow-partial",
                "note_ru": "Эта модель входит в CRISPR-атлас, но текущий multi-omics индекс был построен для более узкого набора моделей. Перестройте индекс после обновления скрипта.",
            }

        resolved = matches[0]
        row = pd.to_numeric(matrix.loc[resolved], errors="coerce").dropna().sort_values(ascending=True)
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        total = int(len(row))
        pan = self._pan_stats()
        meta = self._model_meta(str(resolved))
        cancer_id = str(meta.get("mcl_cancer_id") or "")
        cancer_stats = self._cancer_stats(cancer_id) if cancer_id else pd.DataFrame()

        records: list[dict[str, Any]] = []
        type_counts = {key: 0 for key in DEPENDENCY_TYPES}
        domain_counts: dict[str, dict[str, Any]] = {}

        for rank, (raw_gene, raw_value) in enumerate(row.items(), start=1):
            gene = str(raw_gene).upper()
            ge = float(raw_value)
            pan_row = pan.loc[gene] if not pan.empty and gene in pan.index else None
            cancer_row = cancer_stats.loc[gene] if not cancer_stats.empty and gene in cancer_stats.index else None
            pan_fraction = None if pan_row is None or pd.isna(pan_row.get("pan_dependency_fraction")) else float(pan_row["pan_dependency_fraction"])
            pan_median = None if pan_row is None or pd.isna(pan_row.get("pan_median_gene_effect")) else float(pan_row["pan_median_gene_effect"])
            cancer_fraction = None if cancer_row is None or pd.isna(cancer_row.get("cancer_dependency_fraction")) else float(cancer_row["cancer_dependency_fraction"])
            cancer_median = None if cancer_row is None or pd.isna(cancer_row.get("cancer_median_gene_effect")) else float(cancer_row["cancer_median_gene_effect"])
            cancer_n = 0 if cancer_row is None or pd.isna(cancer_row.get("cancer_models_n")) else int(cancer_row["cancer_models_n"])
            dep_type = self._dependency_type(ge, pan_fraction, pan_median, cancer_fraction, cancer_n)
            type_counts[dep_type] += 1
            annotation = self._annotation(gene)
            for d in annotation["domains"]:
                bucket = domain_counts.setdefault(d["id"], {"id": d["id"], "label_ru": d["label_ru"], "genes_n": 0})
                bucket["genes_n"] += 1
            records.append(
                {
                    "rank": rank,
                    "rank_percentile": 100.0 * (total - rank + 1) / total if total else None,
                    "gene": gene,
                    "gene_name": annotation["gene_name"],
                    "gene_effect": ge,
                    "dependency_type": dep_type,
                    "dependency_type_ru": DEPENDENCY_TYPES[dep_type]["label_ru"],
                    "dependency_description_ru": DEPENDENCY_TYPES[dep_type]["description_ru"],
                    "pan_dependency_fraction": pan_fraction,
                    "pan_median_gene_effect": pan_median,
                    "cancer_dependency_fraction": cancer_fraction,
                    "cancer_median_gene_effect": cancer_median,
                    "cancer_models_n": cancer_n,
                    "delta_vs_pan": ge - pan_median if pan_median is not None else None,
                    "delta_vs_cancer": ge - cancer_median if cancer_median is not None else None,
                    "domains": annotation["domains"],
                    "subdomains": annotation["subdomains"],
                    "protein_classes": annotation["protein_classes"],
                }
            )

        filtered = records
        if search:
            needle = search.strip().casefold()
            filtered = [x for x in filtered if needle in x["gene"].casefold() or needle in str(x.get("gene_name") or "").casefold()]
        if dependency_type and dependency_type in DEPENDENCY_TYPES:
            filtered = [x for x in filtered if x["dependency_type"] == dependency_type]
        if domain:
            filtered = [x for x in filtered if any(d["id"] == domain for d in x.get("domains") or [])]

        filtered_total = len(filtered)
        start = max(0, int(offset))
        stop = start + max(1, min(int(limit), 1000))
        strong_n = sum(1 for x in records if x["gene_effect"] <= DEPENDENCY_THRESHOLD)

        return _clean(
            {
                "model_id": str(resolved),
                "available": True,
                "depmap_release": meta.get("depmap_release"),
                "model": {
                    "cell_line_name": meta.get("cell_line_name"),
                    "organ_ru": meta.get("mcl_organ_ru"),
                    "cancer_name": meta.get("mcl_cancer_name"),
                    "cancer_id": cancer_id or None,
                },
                "genes_measured_n": total,
                "strong_dependencies_n": strong_n,
                "filtered_total": filtered_total,
                "offset": start,
                "limit": min(int(limit), 1000),
                "items": filtered[start:stop],
                "dependency_threshold": DEPENDENCY_THRESHOLD,
                "type_counts": type_counts,
                "facets": {
                    "dependency_types": [
                        {"id": key, "label_ru": value["label_ru"], "genes_n": type_counts[key]}
                        for key, value in DEPENDENCY_TYPES.items()
                    ],
                    "domains": sorted(domain_counts.values(), key=lambda x: (-int(x["genes_n"]), str(x["label_ru"]))),
                },
                "interpretation": {
                    "ranking_ru": "Гены отсортированы по Chronos Gene Effect: чем значение отрицательнее, тем сильнее потеря гена ухудшает рост или выживание этой клеточной модели.",
                    "classification_ru": "Тип зависимости — эвристическая классификация MCL по распространённости сильной зависимости среди проиндексированных опухолевых моделей DepMap. Она не оценивает токсичность для нормальных тканей.",
                    "pharmacology_ru": "CRISPR-нокаут гена не равнозначен фармакологическому ингибированию его белкового продукта.",
                },
            }
        )
