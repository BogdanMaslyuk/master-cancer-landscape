from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from .store import MCLDataError, MCLDataStore


_DISPLAY = {
    "CANCER-001": {
        "organ_ru": "Лёгкое",
        "cancer_ru": "Аденокарцинома лёгкого",
        "molecular_ru": "KRAS G12C",
        "short_ru": "LUAD · KRAS G12C",
        "organ_icon": "lung",
    },
    "CANCER-004": {
        "organ_ru": "Поджелудочная железа",
        "cancer_ru": "Аденокарцинома поджелудочной железы",
        "molecular_ru": "KRAS G12D",
        "short_ru": "PDAC · KRAS G12D",
        "organ_icon": "pancreas",
    },
    "CANCER-011": {
        "organ_ru": "Центральная нервная система",
        "cancer_ru": "Глиобластома",
        "molecular_ru": "IDH-wildtype",
        "short_ru": "Глиобластома · IDH-wildtype",
        "organ_icon": "brain",
    },
    "CANCER-013": {
        "organ_ru": "Кроветворная система",
        "cancer_ru": "Острый миелоидный лейкоз",
        "molecular_ru": "FLT3-mutant",
        "short_ru": "ОМЛ · FLT3-mutant",
        "organ_icon": "blood",
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


def _first(series: pd.Series, default: str | None = None) -> str | None:
    values = series.dropna().astype(str).str.strip()
    values = values[values != ""]
    return values.iloc[0] if not values.empty else default


def _slug(value: str) -> str:
    return "-".join(part for part in "".join(ch.lower() if ch.isalnum() else " " for ch in value).split() if part)


def _parse_variant_json(value: Any, kind: str) -> list[dict[str, Any]]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    try:
        payload = json.loads(str(value))
    except (json.JSONDecodeError, TypeError):
        return []
    if not isinstance(payload, list):
        return []
    output: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str]] = set()
    for row in payload:
        if not isinstance(row, dict):
            continue
        gene = str(row.get("GeneSymbol") or "").strip()
        protein = str(row.get("ProteinChange") or "").strip()
        dna = str(row.get("DNAChange") or "").strip()
        key = (gene, protein, dna)
        if key in seen:
            continue
        seen.add(key)
        output.append(
            {
                "gene": gene or None,
                "protein_change": protein or None,
                "dna_change": dna or None,
                "hotspot": row.get("Hotspot"),
                "driver": row.get("HessDriver"),
                "classification": row.get("VariantClassification"),
                "kind": kind,
            }
        )
    return output


class MCLAtlas:
    """Read-only disease → molecular context → cell model navigation layer."""

    def __init__(self, root: Path, store: MCLDataStore):
        self.root = Path(root).resolve()
        self.store = store

    @lru_cache(maxsize=1)
    def _contexts_config(self) -> dict[str, Any]:
        path = self.root / "config/cancer_contexts.yaml"
        if not path.exists():
            return {}
        payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        return payload.get("contexts") or {}

    @lru_cache(maxsize=1)
    def _audit(self) -> pd.DataFrame:
        parquet = self.root / "data/processed/depmap_context_audit.parquet"
        tsv = self.root / "data/processed/depmap_context_audit.tsv"
        if parquet.exists():
            return pd.read_parquet(parquet)
        if tsv.exists():
            return pd.read_csv(tsv, sep="\t", low_memory=False)
        return pd.DataFrame()

    def _comparison_specs(self, cancer_id: str) -> list[dict[str, Any]]:
        return [x for x in self.store.comparison_specs() if str(x.get("cancer_id")) == cancer_id]

    def _context_summary(self, cancer_id: str, cfg: dict[str, Any]) -> dict[str, Any]:
        audit = self._audit()
        sub = audit[audit["cancer_id"].astype(str) == cancer_id].copy() if not audit.empty and "cancer_id" in audit else pd.DataFrame()
        display = _DISPLAY.get(cancer_id, {})
        inclusion = cfg.get("inclusion") or {}
        ot = cfg.get("opentargets") or {}
        depmap = cfg.get("depmap") or {}

        organ_en = _first(sub.get("oncotree_lineage", pd.Series(dtype=str)), str(inclusion.get("lineage") or "Unknown"))
        primary_disease = _first(sub.get("oncotree_primary_disease", pd.Series(dtype=str)), str(inclusion.get("disease") or cfg.get("name") or cancer_id))
        subtype = _first(sub.get("oncotree_subtype", pd.Series(dtype=str)), primary_disease)
        model_type = _first(sub.get("depmap_model_type", pd.Series(dtype=str)), None)

        groups = sub.get("assigned_group", pd.Series(dtype=str)).astype(str) if not sub.empty else pd.Series(dtype=str)
        sequencing = sub.get("sequencing_available", pd.Series(dtype=bool))
        sequencing_n = int(sequencing.fillna(False).astype(bool).sum()) if not sequencing.empty else 0
        comparisons = self._comparison_specs(cancer_id)

        molecular = str(ot.get("molecular_context") or inclusion.get("alteration") or cfg.get("name") or cancer_id)
        cancer_ru = display.get("cancer_ru") or subtype or primary_disease
        organ_ru = display.get("organ_ru") or organ_en
        molecular_ru = display.get("molecular_ru") or molecular.replace("p.", "")

        return _clean(
            {
                "id": cancer_id,
                "name": cfg.get("name") or cancer_id,
                "organ_en": organ_en,
                "organ_ru": organ_ru,
                "organ_id": _slug(str(organ_en or organ_ru or cancer_id)),
                "organ_icon": display.get("organ_icon") or "tissue",
                "primary_disease": primary_disease,
                "subtype": subtype,
                "cancer_ru": cancer_ru,
                "molecular_context": molecular,
                "molecular_ru": molecular_ru,
                "short_ru": display.get("short_ru") or f"{cancer_ru} · {molecular_ru}",
                "depmap_model_type": model_type,
                "context_definition": depmap.get("context_definition"),
                "comparator_definition": depmap.get("comparator_definition"),
                "project_status": cfg.get("status"),
                "models_n": int(sub["model_id"].nunique()) if not sub.empty and "model_id" in sub else 0,
                "context_models_n": int((groups == "context").sum()),
                "comparator_models_n": int((groups == "comparator").sum()),
                "excluded_models_n": int((groups == "excluded").sum()),
                "sequenced_models_n": sequencing_n,
                "analysis_available": bool(comparisons),
                "comparisons_n": len(comparisons),
                "comparisons": comparisons,
                "minimum_context_n_warning": cfg.get("minimum_context_n_warning"),
            }
        )

    def atlas(self) -> dict[str, Any]:
        contexts = [self._context_summary(cid, cfg or {}) for cid, cfg in self._contexts_config().items()]
        organ_map: dict[str, dict[str, Any]] = {}
        for context in contexts:
            organ_id = str(context["organ_id"])
            organ = organ_map.setdefault(
                organ_id,
                {
                    "id": organ_id,
                    "name_ru": context["organ_ru"],
                    "name_en": context["organ_en"],
                    "icon": context["organ_icon"],
                    "contexts": [],
                },
            )
            organ["contexts"].append(context)

        organs = []
        for organ in organ_map.values():
            contexts_here = organ["contexts"]
            organ["contexts_n"] = len(contexts_here)
            organ["models_n"] = sum(int(x.get("models_n") or 0) for x in contexts_here)
            organ["analyses_n"] = sum(int(x.get("comparisons_n") or 0) for x in contexts_here)
            organs.append(organ)
        organs.sort(key=lambda x: (not bool(x.get("analyses_n")), str(x.get("name_ru"))))

        audit = self._audit()
        unique_models = int(audit["model_id"].nunique()) if not audit.empty and "model_id" in audit else 0
        return _clean(
            {
                "organs_n": len(organs),
                "contexts_n": len(contexts),
                "models_n": unique_models,
                "analyses_n": sum(int(x.get("comparisons_n") or 0) for x in contexts),
                "organs": organs,
            }
        )

    def context(self, cancer_id: str) -> dict[str, Any]:
        cfg = self._contexts_config().get(cancer_id)
        if cfg is None:
            raise MCLDataError(f"Unknown cancer context: {cancer_id}")
        return self._context_summary(cancer_id, cfg or {})

    def models(
        self,
        *,
        cancer_id: str | None = None,
        group: str | None = None,
        search: str | None = None,
        sequencing_only: bool = False,
        limit: int = 1000,
    ) -> dict[str, Any]:
        frame = self._audit().copy()
        if frame.empty:
            return {"total": 0, "items": []}
        if cancer_id:
            if cancer_id not in self._contexts_config():
                raise MCLDataError(f"Unknown cancer context: {cancer_id}")
            frame = frame[frame["cancer_id"].astype(str) == cancer_id]
        if group:
            frame = frame[frame["assigned_group"].astype(str).str.lower() == group.lower()]
        if sequencing_only and "sequencing_available" in frame:
            frame = frame[frame["sequencing_available"].fillna(False).astype(bool)]
        if search:
            needle = search.strip().lower()
            mask = pd.Series(False, index=frame.index)
            for col in ["model_id", "cell_line_name", "depmap_model_type", "oncotree_subtype", "oncotree_primary_disease"]:
                if col in frame:
                    mask = mask | frame[col].astype(str).str.lower().str.contains(needle, regex=False)
            frame = frame[mask]

        total = int(len(frame))
        frame = frame.head(min(max(int(limit), 1), 5000))
        items: list[dict[str, Any]] = []
        for _, row in frame.iterrows():
            cancer = str(row.get("cancer_id") or "")
            cfg = self._contexts_config().get(cancer) or {}
            display = _DISPLAY.get(cancer, {})
            variants = _parse_variant_json(row.get("qualifying_variants_json"), "defining")
            variants += _parse_variant_json(row.get("other_relevant_variants_json"), "other")
            group_value = str(row.get("assigned_group") or "")
            group_ru = {"context": "целевая группа", "comparator": "группа сравнения", "excluded": "исключена"}.get(group_value, group_value)
            items.append(
                _clean(
                    {
                        "cancer_id": cancer,
                        "model_id": row.get("model_id"),
                        "cell_line_name": row.get("cell_line_name"),
                        "depmap_model_type": row.get("depmap_model_type"),
                        "oncotree_lineage": row.get("oncotree_lineage"),
                        "oncotree_primary_disease": row.get("oncotree_primary_disease"),
                        "oncotree_subtype": row.get("oncotree_subtype"),
                        "oncotree_code": row.get("oncotree_code"),
                        "sequencing_available": row.get("sequencing_available"),
                        "alteration_status": row.get("alteration_status"),
                        "assigned_group": group_value,
                        "assigned_group_ru": group_ru,
                        "assignment_reason": row.get("assignment_reason"),
                        "variants": variants,
                        "molecular_context": (cfg.get("opentargets") or {}).get("molecular_context"),
                        "cancer_ru": display.get("cancer_ru") or row.get("oncotree_subtype"),
                        "organ_ru": display.get("organ_ru") or row.get("oncotree_lineage"),
                    }
                )
            )
        return {"total": total, "items": items}

    def model(self, model_id: str) -> dict[str, Any]:
        frame = self._audit()
        if frame.empty or "model_id" not in frame:
            raise MCLDataError(f"Cell model not found: {model_id}")
        hits = frame[frame["model_id"].astype(str).str.upper() == model_id.strip().upper()]
        if hits.empty:
            raise MCLDataError(f"Cell model not found: {model_id}")
        row = hits.iloc[0]
        memberships = self.models(search=str(row.get("model_id")), limit=100)["items"]
        return _clean(
            {
                "model_id": row.get("model_id"),
                "cell_line_name": row.get("cell_line_name"),
                "depmap_model_type": row.get("depmap_model_type"),
                "oncotree_lineage": row.get("oncotree_lineage"),
                "oncotree_primary_disease": row.get("oncotree_primary_disease"),
                "oncotree_subtype": row.get("oncotree_subtype"),
                "oncotree_code": row.get("oncotree_code"),
                "sequencing_available": row.get("sequencing_available"),
                "memberships": memberships,
            }
        )
