from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

from .store import MCLDataError


STATUS_ORDER = {
    "priority_for_in_vitro": 0,
    "supported_hypothesis": 1,
    "exploratory_hypothesis": 2,
    "insufficient_evidence": 3,
}


class CandidateHypothesisStore:
    """Read-only translational hypothesis layer.

    Candidate v2 is preferred when materialized. The store keeps evidence axes
    separate and never turns a hypothesis into a validated drug mechanism or a
    therapeutic claim. Candidate v1 remains a fallback for older local builds.
    """

    def __init__(self, root: Path):
        self.root = Path(root).resolve()
        self.runtime = self.root / "data" / "runtime" / "pharmacology"

    def _paths(self) -> tuple[Path, Path, Path, str]:
        v2 = (
            self.runtime / "candidate_hypotheses_v2.parquet",
            self.runtime / "candidate_hypothesis_models_v2.parquet",
            self.runtime / "candidate_hypotheses_v2_manifest.json",
        )
        if all(path.exists() for path in v2):
            return (*v2, "v2")
        v1 = (
            self.runtime / "candidate_hypotheses.parquet",
            self.runtime / "candidate_hypothesis_models.parquet",
            self.runtime / "candidate_hypotheses_manifest.json",
        )
        return (*v1, "v1")

    @staticmethod
    def _clean(value: Any) -> Any:
        if isinstance(value, dict):
            return {str(k): CandidateHypothesisStore._clean(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [CandidateHypothesisStore._clean(v) for v in value]
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
    def _json_list(value: Any) -> list[str]:
        if isinstance(value, list):
            return [str(x) for x in value]
        if value is None:
            return []
        try:
            payload = json.loads(str(value))
        except (json.JSONDecodeError, TypeError, ValueError):
            return []
        return [str(x) for x in payload] if isinstance(payload, list) else []

    @classmethod
    def _records(cls, frame: pd.DataFrame) -> list[dict[str, Any]]:
        rows = [cls._clean(row) for row in frame.to_dict("records")]
        for row in rows:
            if "priority_reasons_json" in row:
                row["priority_reasons"] = cls._json_list(row.pop("priority_reasons_json"))
            if "evidence_gaps_json" in row:
                row["evidence_gaps"] = cls._json_list(row.pop("evidence_gaps_json"))
            if "protein_changes_json" in row:
                row["protein_changes"] = cls._json_list(row.get("protein_changes_json"))
        return rows

    @lru_cache(maxsize=1)
    def manifest(self) -> dict[str, Any]:
        hypotheses_path, _, manifest_path, version = self._paths()
        if not manifest_path.exists():
            return {
                "available": False,
                "status": "not_built",
                "hypotheses_n": 0,
                "candidate_version": version,
                "build_command": ".\\scripts\\build-candidate-prioritization-v2.ps1",
                "note_ru": "Слой исследовательских гипотез ещё не построен.",
            }
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        payload["available"] = hypotheses_path.exists()
        payload["candidate_version"] = version
        return payload

    @lru_cache(maxsize=1)
    def hypotheses(self) -> pd.DataFrame:
        hypotheses_path, _, _, _ = self._paths()
        if not hypotheses_path.exists():
            return pd.DataFrame()
        frame = pd.read_parquet(hypotheses_path)
        if "target_gene" in frame.columns:
            frame["target_gene"] = frame["target_gene"].astype(str).str.upper()
        return frame

    def summary(self) -> dict[str, Any]:
        payload = dict(self.manifest())
        version = payload.get("candidate_version")
        if version == "v2":
            payload["contract_note_ru"] = (
                "Candidate v2 разделяет абсолютную активность PRISM, относительную селективность, Probability of Dependency, "
                "непрерывный Gene Effect и молекулярный контекст. Мутации, экспрессия и relative copy number помогают "
                "объяснить модель, но сами по себе не повышают гипотезу до причинного механизма."
            )
        else:
            payload["contract_note_ru"] = (
                "Отображается Candidate v1. Для лабораторного выбора предпочтительно построить v2, где бинарная "
                "CRISPR-зависимость основана на Probability of Dependency и добавлен абсолютный порог активности PRISM."
            )
        return self._clean(payload)

    def search(
        self,
        *,
        search: str | None = None,
        status: str | None = None,
        target_gene: str | None = None,
        cancer: str | None = None,
        limit: int = 100,
        offset: int = 0,
    ) -> dict[str, Any]:
        frame = self.hypotheses().copy()
        if frame.empty:
            return {
                **self.summary(),
                "total": 0,
                "items": [],
                "limit": int(limit),
                "offset": int(offset),
            }

        if search:
            needle = search.strip().casefold()
            mask = pd.Series(False, index=frame.index)
            for column in (
                "preferred_name", "compound_id", "protein_preferred_name", "target_gene",
                "uniprot_primary_accession", "mcl_cancer_name", "mcl_organ_ru",
            ):
                if column in frame.columns:
                    mask = mask | frame[column].fillna("").astype(str).str.casefold().str.contains(
                        needle, regex=False
                    )
            frame = frame[mask]

        if status and "priority_status" in frame.columns:
            frame = frame[frame["priority_status"].astype(str) == status]
        if target_gene and "target_gene" in frame.columns:
            frame = frame[
                frame["target_gene"].astype(str).str.upper() == target_gene.strip().upper()
            ]
        if cancer:
            needle = cancer.strip().casefold()
            mask = pd.Series(False, index=frame.index)
            for column in ("mcl_cancer_name", "mcl_organ_ru", "mcl_system_ru"):
                if column in frame.columns:
                    mask = mask | frame[column].fillna("").astype(str).str.casefold().str.contains(
                        needle, regex=False
                    )
            frame = frame[mask]

        frame["_status_order"] = frame.get(
            "priority_status", pd.Series("", index=frame.index)
        ).map(STATUS_ORDER).fillna(99)
        sort_cols = [
            c for c in (
                "_status_order", "joint_support_models_n", "sensitive_models_n",
                "active_models_n", "dependency_models_n", "models_n"
            ) if c in frame.columns
        ]
        ascending = [True if c == "_status_order" else False for c in sort_cols]
        if sort_cols:
            frame = frame.sort_values(sort_cols, ascending=ascending, na_position="last")
        frame = frame.drop(columns="_status_order", errors="ignore")

        total = int(len(frame))
        start = max(0, int(offset))
        page = frame.iloc[start:start + max(1, min(int(limit), 500))]
        version = self.manifest().get("candidate_version")
        interpretation = (
            "Candidate v2: высокий статус требует измеренного абсолютного PRISM-эффекта, относительной селективности, "
            "Probability of Dependency и непротиворечащей профильной согласованности. Универсальный числовой score не используется."
            if version == "v2"
            else
            "Candidate v1 показан для совместимости; для решений о лабораторной проверке постройте Candidate v2."
        )
        return self._clean(
            {
                "available": True,
                "candidate_version": version,
                "total": total,
                "limit": int(limit),
                "offset": start,
                "items": self._records(page),
                "interpretation_ru": interpretation,
            }
        )

    def detail(self, hypothesis_id: str) -> dict[str, Any]:
        frame = self.hypotheses()
        if frame.empty or "hypothesis_id" not in frame.columns:
            raise MCLDataError("Candidate hypothesis index is not available")
        hit = frame[frame["hypothesis_id"].astype(str) == hypothesis_id]
        if hit.empty:
            raise MCLDataError(f"Unknown candidate hypothesis: {hypothesis_id}")
        hypothesis = self._records(hit.head(1))[0]

        _, models_path, _, version = self._paths()
        models = pd.DataFrame()
        if models_path.exists():
            try:
                models = pd.read_parquet(
                    models_path,
                    filters=[("hypothesis_id", "==", hypothesis_id)],
                    engine="pyarrow",
                )
            except Exception:
                all_models = pd.read_parquet(models_path)
                models = all_models[
                    all_models["hypothesis_id"].astype(str) == hypothesis_id
                ].copy()

        roles: dict[str, list[dict[str, Any]]] = {}
        if not models.empty and "model_role" in models.columns:
            for role, group in models.groupby("model_role", sort=False):
                roles[str(role)] = self._records(group)

        lab_route = [
            {
                "stage": "1. Подтвердить фенотип",
                "goal": (
                    "Построить независимую dose-response кривую на выбранных положительных моделях и сравнить "
                    "её с отрицательными моделями того же опухолевого контекста. Single-dose PRISM используется "
                    "для выбора моделей, а не как финальная оценка мощности препарата."
                ),
            },
            {
                "stage": "2. Проверить механизм",
                "goal": (
                    "Добавить мишень-специфический биохимический или клеточный readout: target engagement, "
                    "активность белка либо ближайший downstream-маркер. Конкретный метод зависит от класса мишени."
                ),
            },
            {
                "stage": "3. Проверить причинность",
                "goal": (
                    "Сопоставить фармакологическое воздействие с генетической perturbation/rescue-валидацией. "
                    "Несовпадение не игнорировать: оно может указывать на полифармакологию или иной механизм."
                ),
            },
        ]

        falsification = [
            "Положительные модели не подтверждают воспроизводимую dose-response чувствительность при независимом эксперименте.",
            "Фармакологический эффект не сопровождается ожидаемым изменением мишени или downstream-маркера.",
            "Отрицательные модели отвечают так же сильно, как предполагаемые положительные, без ожидаемой зависимости от мишени.",
            "Генетическая потеря функции и фармакологическое воздействие дают принципиально несогласованные фенотипы без объяснимого механизма.",
        ]

        return self._clean(
            {
                "available": True,
                "candidate_version": version,
                "hypothesis": hypothesis,
                "models_by_role": roles,
                "laboratory_route": lab_route,
                "falsification_criteria": falsification,
                "interpretation_ru": (
                    "Эта карточка является планом проверки исследовательской гипотезы. Молекулярный контекст "
                    "помогает выбирать и объяснять модели, но не является автоматическим доказательством причинности."
                ),
            }
        )
