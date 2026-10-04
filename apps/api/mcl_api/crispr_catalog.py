from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

from .store import MCLDataError


_BUILD_COMMAND = ".\\.venv\\Scripts\\python.exe .\\scripts\\build_crispr_cancer_atlas.py"


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


def _slug(value: Any) -> str:
    text = str(value or "").strip().lower()
    return "-".join("".join(ch if ch.isalnum() else " " for ch in text).split()) or "unknown"


def _truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value or "").strip().lower() in {"1", "true", "yes", "y", "t"}


def _first_nonempty(row: pd.Series, *columns: str) -> str | None:
    for column in columns:
        value = row.get(column)
        if value is None:
            continue
        try:
            if pd.isna(value):
                continue
        except (TypeError, ValueError):
            pass
        text = str(value).strip()
        if text:
            return text
    return None


class CRISPRModelCatalog:
    """Navigation layer over every DepMap model with a CRISPR Gene Effect profile.

    The canonical input is data/processed/depmap_crispr_model_atlas.*, built by
    scripts/build_crispr_cancer_atlas.py from CRISPRGeneEffect.csv + Model.csv.
    Until that local index exists, the catalog degrades to the already-versioned
    context audit so the Explorer remains usable instead of failing at startup.
    """

    def __init__(self, root: Path):
        self.root = Path(root).resolve()

    @lru_cache(maxsize=1)
    def _frame(self) -> pd.DataFrame:
        parquet = self.root / "data/processed/depmap_crispr_model_atlas.parquet"
        tsv = self.root / "data/processed/depmap_crispr_model_atlas.tsv"
        if parquet.exists():
            frame = pd.read_parquet(parquet)
        elif tsv.exists():
            frame = pd.read_csv(tsv, sep="\t", low_memory=False)
        else:
            frame = self._fallback_frame()
        if frame.empty:
            return frame
        if "model_id" not in frame.columns:
            raise MCLDataError("CRISPR model atlas is missing model_id")
        frame = frame.copy()
        frame["model_id"] = frame["model_id"].astype(str).str.strip()
        return frame.drop_duplicates("model_id")

    @lru_cache(maxsize=1)
    def _audit(self) -> pd.DataFrame:
        parquet = self.root / "data/processed/depmap_context_audit.parquet"
        tsv = self.root / "data/processed/depmap_context_audit.tsv"
        if parquet.exists():
            return pd.read_parquet(parquet)
        if tsv.exists():
            return pd.read_csv(tsv, sep="\t", low_memory=False)
        return pd.DataFrame()

    def _fallback_frame(self) -> pd.DataFrame:
        audit = self._audit().copy()
        if audit.empty:
            return pd.DataFrame()
        rename = {}
        # The audit already uses normalized snake_case names in current MCL builds.
        frame = audit.rename(columns=rename).drop_duplicates("model_id").copy()
        lineage = frame.get("oncotree_lineage", pd.Series("", index=frame.index)).fillna("").astype(str)
        cancer = frame.get("oncotree_primary_disease", pd.Series("", index=frame.index)).fillna("").astype(str)
        subtype = frame.get("oncotree_subtype", pd.Series("", index=frame.index)).fillna("").astype(str)
        frame["mcl_system_ru"] = lineage.where(lineage != "", "Не классифицировано")
        frame["mcl_organ_ru"] = lineage.where(lineage != "", "Не классифицировано")
        frame["mcl_organ_icon"] = "tissue"
        frame["mcl_system_id"] = frame["mcl_system_ru"].map(_slug)
        frame["mcl_organ_id"] = frame["mcl_organ_ru"].map(_slug)
        frame["mcl_cancer_name"] = cancer.where(cancer != "", subtype)
        frame["mcl_cancer_id"] = (lineage + "|" + frame["mcl_cancer_name"]).map(_slug)
        frame["mcl_subtype_name"] = subtype
        frame["mcl_subtype_id"] = (frame["mcl_cancer_id"] + "|" + subtype).map(_slug)
        frame["classification_status"] = "fallback"
        frame["classification_confidence"] = "medium"
        frame["has_crispr"] = True
        if "depmap_release" not in frame.columns:
            frame["depmap_release"] = None
        return frame

    def source_status(self) -> dict[str, Any]:
        canonical = (
            (self.root / "data/processed/depmap_crispr_model_atlas.parquet").exists()
            or (self.root / "data/processed/depmap_crispr_model_atlas.tsv").exists()
        )
        return {
            "status": "ready" if canonical else "fallback_context_audit",
            "canonical_index_available": canonical,
            "build_command": None if canonical else _BUILD_COMMAND,
            "note_ru": (
                "Атлас построен по всем моделям, присутствующим в CRISPRGeneEffect.csv."
                if canonical
                else "Полный индекс CRISPR-моделей ещё не построен; временно показаны модели из существующих контекстов MCL."
            ),
        }

    def _memberships(self, model_id: str) -> list[dict[str, Any]]:
        audit = self._audit()
        if audit.empty or "model_id" not in audit.columns:
            return []
        rows = audit[audit["model_id"].astype(str).str.upper() == model_id.upper()]
        columns = [
            "cancer_id", "assigned_group", "assignment_reason", "alteration_status",
            "sequencing_available", "qualifying_variants_json", "other_relevant_variants_json",
        ]
        keep = [c for c in columns if c in rows.columns]
        return [_clean(x) for x in rows[keep].to_dict(orient="records")]

    def atlas(self) -> dict[str, Any]:
        frame = self._frame().copy()
        status = self.source_status()
        if frame.empty:
            return {
                **status,
                "models_n": 0,
                "organs_n": 0,
                "cancers_n": 0,
                "subtypes_n": 0,
                "requires_review_n": 0,
                "organs": [],
            }

        organs: list[dict[str, Any]] = []
        for organ_id, organ_rows in frame.groupby(frame.get("mcl_organ_id", pd.Series("unknown", index=frame.index)).fillna("unknown").astype(str)):
            first = organ_rows.iloc[0]
            cancers: list[dict[str, Any]] = []
            cancer_key = organ_rows.get("mcl_cancer_id", pd.Series("unknown", index=organ_rows.index)).fillna("unknown").astype(str)
            for cancer_id, cancer_rows in organ_rows.groupby(cancer_key):
                cancer_first = cancer_rows.iloc[0]
                subtypes: list[dict[str, Any]] = []
                subtype_series = cancer_rows.get("mcl_subtype_name", pd.Series("", index=cancer_rows.index)).fillna("").astype(str).str.strip()
                for subtype_name in sorted({x for x in subtype_series.tolist() if x}):
                    subtype_rows = cancer_rows[subtype_series == subtype_name]
                    subtypes.append(
                        {
                            "id": _slug(f"{cancer_id}|{subtype_name}"),
                            "name": subtype_name,
                            "models_n": int(subtype_rows["model_id"].nunique()),
                        }
                    )

                memberships: set[str] = set()
                audit = self._audit()
                if not audit.empty and {"model_id", "cancer_id"}.issubset(audit.columns):
                    model_ids = set(cancer_rows["model_id"].astype(str))
                    memberships = set(
                        audit.loc[audit["model_id"].astype(str).isin(model_ids), "cancer_id"]
                        .dropna().astype(str).tolist()
                    )
                cancers.append(
                    {
                        "id": str(cancer_id),
                        "name": _first_nonempty(cancer_first, "mcl_cancer_name", "oncotree_primary_disease", "oncotree_subtype") or "Не классифицировано",
                        "oncotree_primary_disease": _first_nonempty(cancer_first, "oncotree_primary_disease"),
                        "models_n": int(cancer_rows["model_id"].nunique()),
                        "subtypes_n": len(subtypes),
                        "subtypes": subtypes,
                        "curated_context_ids": sorted(memberships),
                    }
                )
            cancers.sort(key=lambda x: (-int(x["models_n"]), str(x["name"])))
            organs.append(
                {
                    "id": str(organ_id),
                    "system_ru": _first_nonempty(first, "mcl_system_ru") or "Прочие опухоли",
                    "name_ru": _first_nonempty(first, "mcl_organ_ru", "oncotree_lineage") or "Не классифицировано",
                    "name_en": _first_nonempty(first, "oncotree_lineage"),
                    "icon": _first_nonempty(first, "mcl_organ_icon") or "tissue",
                    "models_n": int(organ_rows["model_id"].nunique()),
                    "cancers_n": len(cancers),
                    "cancers": cancers,
                }
            )
        organs.sort(key=lambda x: (-int(x["models_n"]), str(x["name_ru"])))

        classification = frame.get("classification_status", pd.Series("", index=frame.index)).fillna("").astype(str)
        subtype_values = frame.get("mcl_subtype_name", pd.Series("", index=frame.index)).fillna("").astype(str).str.strip()
        return _clean(
            {
                **status,
                "depmap_release": _first_nonempty(frame.iloc[0], "depmap_release"),
                "models_n": int(frame["model_id"].nunique()),
                "organs_n": len(organs),
                "cancers_n": int(frame.get("mcl_cancer_id", pd.Series(dtype=str)).dropna().astype(str).nunique()),
                "subtypes_n": int(frame.loc[subtype_values != "", "mcl_subtype_id"].dropna().astype(str).nunique()) if "mcl_subtype_id" in frame.columns else int(subtype_values[subtype_values != ""].nunique()),
                "requires_review_n": int(classification.eq("requires_review").sum()),
                "organs": organs,
            }
        )

    def models(
        self,
        *,
        organ_id: str | None = None,
        cancer_id: str | None = None,
        subtype: str | None = None,
        search: str | None = None,
        limit: int = 5000,
    ) -> dict[str, Any]:
        frame = self._frame().copy()
        if frame.empty:
            return {"total": 0, "items": [], **self.source_status()}
        if organ_id and "mcl_organ_id" in frame.columns:
            frame = frame[frame["mcl_organ_id"].astype(str) == organ_id]
        if cancer_id and "mcl_cancer_id" in frame.columns:
            frame = frame[frame["mcl_cancer_id"].astype(str) == cancer_id]
        if subtype:
            needle = subtype.strip().lower()
            values = frame.get("mcl_subtype_name", pd.Series("", index=frame.index)).fillna("").astype(str).str.lower()
            frame = frame[values == needle]
        if search:
            needle = search.strip().lower()
            mask = pd.Series(False, index=frame.index)
            for col in [
                "model_id", "cell_line_name", "depmap_model_type", "oncotree_lineage",
                "oncotree_primary_disease", "oncotree_subtype", "oncotree_code",
                "mcl_organ_ru", "mcl_cancer_name", "patient_subtype_features",
            ]:
                if col in frame.columns:
                    mask = mask | frame[col].fillna("").astype(str).str.lower().str.contains(needle, regex=False)
            frame = frame[mask]

        total = int(frame["model_id"].nunique())
        frame = frame.head(min(max(int(limit), 1), 5000))
        items: list[dict[str, Any]] = []
        for _, row in frame.iterrows():
            model_id = str(row.get("model_id") or "")
            memberships = self._memberships(model_id)
            items.append(
                _clean(
                    {
                        "model_id": model_id,
                        "cell_line_name": row.get("cell_line_name"),
                        "depmap_model_type": row.get("depmap_model_type"),
                        "oncotree_lineage": row.get("oncotree_lineage"),
                        "oncotree_primary_disease": row.get("oncotree_primary_disease"),
                        "oncotree_subtype": row.get("oncotree_subtype"),
                        "oncotree_code": row.get("oncotree_code"),
                        "mcl_system_ru": row.get("mcl_system_ru"),
                        "mcl_organ_id": row.get("mcl_organ_id"),
                        "mcl_organ_ru": row.get("mcl_organ_ru"),
                        "mcl_cancer_id": row.get("mcl_cancer_id"),
                        "mcl_cancer_name": row.get("mcl_cancer_name"),
                        "mcl_subtype_name": row.get("mcl_subtype_name"),
                        "primary_or_metastasis": row.get("primary_or_metastasis"),
                        "sample_collection_site": row.get("sample_collection_site"),
                        "sex": row.get("sex"),
                        "age": row.get("age"),
                        "classification_status": row.get("classification_status"),
                        "classification_confidence": row.get("classification_confidence"),
                        "has_crispr": _truthy(row.get("has_crispr", True)),
                        "curated_context_ids": sorted({str(x.get("cancer_id")) for x in memberships if x.get("cancer_id")}),
                        "curated_contexts_n": len({str(x.get("cancer_id")) for x in memberships if x.get("cancer_id")}),
                    }
                )
            )
        return {"total": total, "items": items, **self.source_status()}

    def model(self, model_id: str) -> dict[str, Any]:
        frame = self._frame()
        if frame.empty:
            raise MCLDataError(f"CRISPR cell model not found: {model_id}")
        hits = frame[frame["model_id"].astype(str).str.upper() == model_id.strip().upper()]
        if hits.empty:
            raise MCLDataError(f"CRISPR cell model not found: {model_id}")
        row = hits.iloc[0]
        resolved = str(row.get("model_id"))
        memberships = self._memberships(resolved)
        return _clean(
            {
                **row.to_dict(),
                "memberships": memberships,
                "sequencing_available": bool(memberships and any(_truthy(x.get("sequencing_available")) for x in memberships)),
                "metadata": row.to_dict(),
                "genetics": {
                    "availability": "not_indexed",
                    "availability_ru": "мутационный профиль не загружен на этом уровне",
                    "mutations_n": 0,
                    "mutated_genes_n": 0,
                    "driver_mutations_n": 0,
                    "hotspot_mutations_n": 0,
                    "likely_lof_n": 0,
                    "high_impact_n": 0,
                    "priority_variants": [],
                    "note": "Карточка уже входит в CRISPR-атлас; полный мутационный индекс строится отдельным шагом.",
                },
                "source_status": self.source_status(),
            }
        )
