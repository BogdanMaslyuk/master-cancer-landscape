from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

from .store import MCLDataError


READINESS_ORDER = {
    "ready_with_internal_tumor_control": 0,
    "ready_positive_model": 1,
    "discordant_current_panel": 2,
    "no_support_in_current_panel": 3,
}


class LaboratoryPanelStore:
    """Read-only view of the department's physical cell-line collection and testable Candidate v2 hypotheses."""

    def __init__(self, root: Path):
        self.root = Path(root).resolve()
        self.panel_path = self.root / "data" / "processed" / "laboratory_panel.parquet"
        self.panel_manifest = self.root / "data" / "processed" / "laboratory_panel_manifest.json"
        self.runtime = self.root / "data" / "runtime" / "laboratory"
        self.candidates_path = self.runtime / "laboratory_candidate_hypotheses.parquet"
        self.models_path = self.runtime / "laboratory_candidate_models.parquet"
        self.candidates_manifest = self.runtime / "laboratory_candidates_manifest.json"
        self.mechanism_panels_path = self.runtime / "laboratory_mechanism_panels.parquet"
        self.mechanism_models_path = self.runtime / "laboratory_mechanism_panel_models.parquet"
        self.mechanism_manifest = self.runtime / "laboratory_mechanism_panels_manifest.json"

    @staticmethod
    def _clean(value: Any) -> Any:
        if isinstance(value, dict):
            return {str(k): LaboratoryPanelStore._clean(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [LaboratoryPanelStore._clean(v) for v in value]
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

    @staticmethod
    def _json_value(value: Any, default: Any) -> Any:
        if isinstance(value, (list, dict)):
            return value
        if value is None:
            return default
        try:
            parsed = json.loads(str(value))
        except (json.JSONDecodeError, TypeError, ValueError):
            return default
        return parsed

    @classmethod
    def _records(cls, frame: pd.DataFrame) -> list[dict[str, Any]]:
        return [cls._clean(row) for row in frame.to_dict("records")]

    @lru_cache(maxsize=1)
    def panel(self) -> pd.DataFrame:
        if not self.panel_path.exists():
            return pd.DataFrame()
        return pd.read_parquet(self.panel_path)

    @lru_cache(maxsize=1)
    def candidates(self) -> pd.DataFrame:
        if not self.candidates_path.exists():
            return pd.DataFrame()
        return pd.read_parquet(self.candidates_path)

    def summary(self) -> dict[str, Any]:
        if not self.panel_manifest.exists():
            return {
                "available": False,
                "status": "not_built",
                "build_command": ".\\scripts\\build-laboratory-panel.ps1",
                "note_ru": "Лабораторная панель ещё не сопоставлена с DepMap.",
            }
        panel = json.loads(self.panel_manifest.read_text(encoding="utf-8"))
        candidate = (
            json.loads(self.candidates_manifest.read_text(encoding="utf-8"))
            if self.candidates_manifest.exists() else {}
        )
        mechanism = (
            json.loads(self.mechanism_manifest.read_text(encoding="utf-8"))
            if self.mechanism_manifest.exists() else {}
        )
        return self._clean({
            "available": self.panel_path.exists(),
            "candidate_triage_available": self.candidates_path.exists(),
            "mechanism_panels_available": self.mechanism_panels_path.exists(),
            "panel": panel,
            "candidate_triage": candidate,
            "mechanism_context": mechanism,
            "interpretation_ru": (
                "Лабораторная панель — физическое ограничение эксперимента, а не новый источник биологической доказательности. "
                "BJ5ta используется только как доступный общий человеческий неопухолевый контроль и не считается органоспецифической нормальной тканью."
            ),
        })

    def lines(
        self,
        *,
        search: str | None = None,
        species: str | None = None,
        role: str | None = None,
        matched_only: bool = False,
    ) -> dict[str, Any]:
        frame = self.panel().copy()
        if frame.empty:
            return {"available": False, "total": 0, "items": [], "build_command": ".\\scripts\\build-laboratory-panel.ps1"}
        if search:
            needle = search.strip().casefold()
            mask = pd.Series(False, index=frame.index)
            for col in ("lab_name", "aliases", "disease_group_ru", "model_id", "depmap_cell_line_name", "preferred_candidate_model_id"):
                if col in frame.columns:
                    mask |= frame[col].fillna("").astype(str).str.casefold().str.contains(needle, regex=False)
            frame = frame[mask]
        if species:
            frame = frame[frame["species"].astype(str).str.lower() == species.strip().lower()]
        if role:
            frame = frame[frame["laboratory_role"].astype(str) == role]
        if matched_only:
            frame = frame[frame["match_status"].astype(str) == "matched"]
        role_order = {"human_non_tumor_control": 0, "human_tumor": 1, "human_technical": 2, "other": 3}
        frame["_role"] = frame["laboratory_role"].map(role_order).fillna(9)
        frame = frame.sort_values(["_role", "highlighted_on_source", "lab_name"], ascending=[True, False, True]).drop(columns="_role")
        return {"available": True, "total": int(len(frame)), "items": self._records(frame)}

    def candidate_list(
        self,
        *,
        search: str | None = None,
        readiness: str | None = None,
        cancer: str | None = None,
        target_gene: str | None = None,
        limit: int = 100,
        offset: int = 0,
    ) -> dict[str, Any]:
        frame = self.candidates().copy()
        if frame.empty:
            return {
                "available": False,
                "total": 0,
                "items": [],
                "build_command": ".\\scripts\\build-laboratory-panel.ps1",
            }
        if search:
            needle = search.strip().casefold()
            mask = pd.Series(False, index=frame.index)
            for col in ("preferred_name", "compound_id", "target_gene", "protein_preferred_name", "mcl_cancer_name", "positive_lab_lines", "negative_lab_lines"):
                if col in frame.columns:
                    mask |= frame[col].fillna("").astype(str).str.casefold().str.contains(needle, regex=False)
            frame = frame[mask]
        if readiness:
            frame = frame[frame["laboratory_readiness"].astype(str) == readiness]
        if cancer:
            needle = cancer.strip().casefold()
            frame = frame[
                frame["mcl_cancer_name"].fillna("").astype(str).str.casefold().str.contains(needle, regex=False)
                | frame["mcl_organ_ru"].fillna("").astype(str).str.casefold().str.contains(needle, regex=False)
            ]
        if target_gene:
            frame = frame[frame["target_gene"].astype(str).str.upper() == target_gene.strip().upper()]
        frame["_readiness"] = frame["laboratory_readiness"].map(READINESS_ORDER).fillna(99)
        frame = frame.sort_values(
            ["_readiness", "laboratory_positive_models_n", "laboratory_negative_models_n", "joint_support_models_n"],
            ascending=[True, False, False, False], na_position="last"
        ).drop(columns="_readiness")
        total = int(len(frame))
        start = max(0, int(offset))
        page = frame.iloc[start:start + max(1, min(int(limit), 500))]
        return self._clean({
            "available": True,
            "total": total,
            "offset": start,
            "limit": int(limit),
            "items": self._records(page),
            "interpretation_ru": (
                "В список входят только 222 Candidate v2 приоритета. Статус лабораторной готовности отвечает на вопрос, "
                "есть ли среди реально доступных линий подходящая положительная модель и, отдельно, внутренний опухолевый контроль."
            ),
        })

    def mechanism_panels(self) -> dict[str, Any]:
        if not self.mechanism_panels_path.exists():
            return {
                "available": False,
                "total": 0,
                "items": [],
                "build_command": ".\\scripts\\build-laboratory-panel.ps1",
            }
        panels = pd.read_parquet(self.mechanism_panels_path)
        models = pd.read_parquet(self.mechanism_models_path) if self.mechanism_models_path.exists() else pd.DataFrame()
        items: list[dict[str, Any]] = []
        json_fields = (
            "compound_names_json", "context_genes_json", "focus_lab_lines_json",
            "candidate_names_json", "candidate_cancers_json", "positive_lab_lines_json",
            "negative_lab_lines_json", "discordant_lab_lines_json",
        )
        for row in self._records(panels):
            panel_id = str(row.get("panel_id") or "")
            for field in json_fields:
                if field in row:
                    row[field.removesuffix("_json")] = self._json_value(row.pop(field), [])
            panel_models = models[models["panel_id"].astype(str) == panel_id].copy() if not models.empty else pd.DataFrame()
            model_records = self._records(panel_models)
            for model in model_records:
                if "mechanism_context_json" in model:
                    model["mechanism_context"] = self._json_value(model.pop("mechanism_context_json"), [])
            row["models"] = model_records
            items.append(row)
        manifest = (
            json.loads(self.mechanism_manifest.read_text(encoding="utf-8"))
            if self.mechanism_manifest.exists() else {}
        )
        return self._clean({
            "available": True,
            "total": len(items),
            "items": items,
            "contract": manifest,
            "interpretation_ru": (
                "Механистические панели объединяют доступность клеток, фармакологический ответ, CRISPR и контекст других генов механизма. "
                "Контекст не является дополнительным баллом и не превращает ассоциацию в доказанную причинность."
            ),
        })

    def candidate_detail(self, hypothesis_id: str) -> dict[str, Any]:
        frame = self.candidates()
        hit = frame[frame["hypothesis_id"].astype(str) == hypothesis_id] if not frame.empty else pd.DataFrame()
        if hit.empty:
            raise MCLDataError(f"Laboratory candidate not found: {hypothesis_id}")
        models = pd.DataFrame()
        if self.models_path.exists():
            try:
                models = pd.read_parquet(self.models_path, filters=[("hypothesis_id", "==", hypothesis_id)], engine="pyarrow")
            except Exception:
                all_models = pd.read_parquet(self.models_path)
                models = all_models[all_models["hypothesis_id"].astype(str) == hypothesis_id].copy()
        roles: dict[str, list[dict[str, Any]]] = {}
        if not models.empty:
            for role, group in models.groupby("laboratory_model_role", sort=False):
                roles[str(role)] = self._records(group)
        return self._clean({
            "available": True,
            "candidate": self._records(hit.head(1))[0],
            "laboratory_models_by_role": roles,
        })
