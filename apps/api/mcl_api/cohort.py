from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from .store import MCLDataError


def _truthy(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and not pd.isna(value):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes", "y", "t"}


def _pct(numerator: int, denominator: int) -> float | None:
    if denominator <= 0:
        return None
    return round((numerator / denominator) * 100.0, 1)


class MCLModelCohortStore:
    """Read-only audit of the cell-model cohort supporting each cancer context."""

    def __init__(self, root: Path):
        self.root = Path(root).resolve()

    @lru_cache(maxsize=1)
    def _contexts(self) -> dict[str, Any]:
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

    @lru_cache(maxsize=1)
    def _metadata(self) -> pd.DataFrame:
        parquet = self.root / "data/processed/depmap_model_metadata.parquet"
        tsv = self.root / "data/processed/depmap_model_metadata.tsv"
        if parquet.exists():
            return pd.read_parquet(parquet)
        if tsv.exists():
            return pd.read_csv(tsv, sep="\t", low_memory=False)
        return pd.DataFrame()

    @lru_cache(maxsize=1)
    def _mutations(self) -> pd.DataFrame:
        parquet = self.root / "data/processed/depmap_model_mutations.parquet"
        tsv = self.root / "data/processed/depmap_model_mutations.tsv"
        if parquet.exists():
            return pd.read_parquet(parquet, columns=["model_id"])
        if tsv.exists():
            return pd.read_csv(tsv, sep="\t", usecols=["model_id"], low_memory=False)
        return pd.DataFrame(columns=["model_id"])

    @staticmethod
    def _ids(frame: pd.DataFrame) -> set[str]:
        if frame.empty or "model_id" not in frame.columns:
            return set()
        return set(frame["model_id"].dropna().astype(str).str.strip())

    def _group_summary(
        self,
        frame: pd.DataFrame,
        group: str,
        sequenced_ids: set[str],
        profiled_ids: set[str],
        metadata_ids: set[str],
    ) -> dict[str, Any]:
        if frame.empty or "assigned_group" not in frame.columns:
            ids: set[str] = set()
        else:
            rows = frame[frame["assigned_group"].astype(str).str.lower() == group]
            ids = self._ids(rows)
        sequenced = ids & sequenced_ids
        profiled = ids & profiled_ids
        metadata = ids & metadata_ids
        return {
            "id": group,
            "models_n": len(ids),
            "sequencing_n": len(sequenced),
            "sequencing_pct": _pct(len(sequenced), len(ids)),
            "mutation_profile_n": len(profiled),
            "mutation_profile_pct": _pct(len(profiled), len(ids)),
            "metadata_n": len(metadata),
            "metadata_pct": _pct(len(metadata), len(ids)),
        }

    def summary(self, cancer_id: str) -> dict[str, Any]:
        cfg = self._contexts().get(cancer_id)
        if cfg is None:
            raise MCLDataError(f"Unknown cancer context: {cancer_id}")

        audit = self._audit()
        if audit.empty or "cancer_id" not in audit.columns:
            raise MCLDataError("DepMap context audit is not available")
        sub = audit[audit["cancer_id"].astype(str) == cancer_id].copy()

        all_ids = self._ids(sub)
        if "sequencing_available" in sub.columns:
            sequenced_ids = self._ids(sub[sub["sequencing_available"].map(_truthy)])
        else:
            sequenced_ids = set()

        mutations = self._mutations()
        profiled_ids = self._ids(mutations) & all_ids
        metadata_ids = self._ids(self._metadata()) & all_ids

        groups = {
            name: self._group_summary(sub, name, sequenced_ids, profiled_ids, metadata_ids)
            for name in ("context", "comparator", "excluded")
        }

        total = {
            "models_n": len(all_ids),
            "sequencing_n": len(sequenced_ids),
            "sequencing_pct": _pct(len(sequenced_ids), len(all_ids)),
            "mutation_profile_n": len(profiled_ids),
            "mutation_profile_pct": _pct(len(profiled_ids), len(all_ids)),
            "metadata_n": len(metadata_ids),
            "metadata_pct": _pct(len(metadata_ids), len(all_ids)),
        }

        minimum_n = int(cfg.get("minimum_context_n_warning") or 0)
        flags: list[dict[str, str]] = []
        context_n = int(groups["context"]["models_n"])
        comparator_n = int(groups["comparator"]["models_n"])

        if minimum_n and context_n < minimum_n:
            flags.append(
                {
                    "level": "warning",
                    "code": "low_context_sample",
                    "title_ru": "Малая целевая группа",
                    "text_ru": f"В целевой группе {context_n} моделей при настроенном пороге предупреждения n={minimum_n}. Выводы следует считать исследовательскими.",
                }
            )
        if comparator_n == 0:
            flags.append(
                {
                    "level": "warning",
                    "code": "no_comparator",
                    "title_ru": "Нет группы сравнения",
                    "text_ru": "Для этого молекулярного контекста в текущем аудите нет моделей, назначенных в контрольную группу. Контекст-специфический эффект нельзя интерпретировать как сравнительный.",
                }
            )
        if len(profiled_ids) < len(sequenced_ids):
            missing = len(sequenced_ids - profiled_ids)
            flags.append(
                {
                    "level": "info",
                    "code": "mutation_index_incomplete",
                    "title_ru": "Мутационный индекс покрывает не все секвенированные модели",
                    "text_ru": f"Для {missing} секвенированных моделей не найдено мутационных записей в компактном индексе. Это отражается как покрытие данных, а не как отсутствие мутаций.",
                }
            )
        if len(metadata_ids) < len(all_ids):
            missing = len(all_ids - metadata_ids)
            flags.append(
                {
                    "level": "info",
                    "code": "metadata_incomplete",
                    "title_ru": "Метаданные моделей неполны",
                    "text_ru": f"Для {missing} моделей отсутствует строка в локальном индексе Model.csv. Основная принадлежность к контексту при этом остаётся доступна из DepMap-аудита.",
                }
            )

        warning_n = sum(1 for flag in flags if flag["level"] == "warning")
        status = "caution" if warning_n else "ready"
        status_ru = (
            "Есть ограничения состава модельной когорты"
            if status == "caution"
            else "Базовые критерии состава модельной когорты выполнены"
        )

        depmap = cfg.get("depmap") or {}
        release = None
        metadata = self._metadata()
        if not metadata.empty and "depmap_release" in metadata.columns:
            values = metadata["depmap_release"].dropna().astype(str).str.strip()
            if not values.empty:
                release = values.iloc[0]

        return {
            "cancer_id": cancer_id,
            "status": status,
            "status_ru": status_ru,
            "minimum_context_n_warning": minimum_n or None,
            "total": total,
            "groups": groups,
            "flags": flags,
            "definitions": {
                "context": depmap.get("context_definition"),
                "comparator": depmap.get("comparator_definition"),
            },
            "representativeness": {
                "status": "not_assessed",
                "status_ru": "Репрезентативность относительно опухолей пациентов ещё не оценена",
                "reason_ru": "Для такой оценки нужна подключённая пациентская молекулярная когорта. MCL пока не заменяет её статистикой клеточных линий.",
            },
            "provenance": {
                "depmap_release": release,
                "membership_source": "data/processed/depmap_context_audit",
                "metadata_source": "data/processed/depmap_model_metadata",
                "mutation_source": "data/processed/depmap_model_mutations",
            },
        }
