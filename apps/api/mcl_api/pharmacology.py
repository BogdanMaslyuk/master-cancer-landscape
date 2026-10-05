from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

from .store import MCLDataError


class MCLPharmacologyStore:
    """Read-only pharmacology layer over materialized runtime tables.

    Large response/link tables are never cached as whole pandas DataFrames during
    interactive use. Model, compound and target pages read only the required Parquet
    slice. This keeps the rest of MCL responsive after pharmacology pages are opened.

    Drug-response observations and compound-target evidence are intentionally kept
    separate. A cell-line response does not by itself prove the annotated target is
    responsible for the phenotype. CRISPR concordance is therefore exposed as an
    additional evidence layer, never as a mechanistic conclusion.
    """

    def __init__(self, root: Path):
        self.root = Path(root).resolve()
        self.runtime = self.root / "data" / "runtime" / "pharmacology"

    @staticmethod
    def _clean(value: Any) -> Any:
        if isinstance(value, dict):
            return {str(k): MCLPharmacologyStore._clean(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [MCLPharmacologyStore._clean(v) for v in value]
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
    def _records(frame: pd.DataFrame) -> list[dict[str, Any]]:
        return [MCLPharmacologyStore._clean(x) for x in frame.to_dict("records")]

    def _read_slice(
        self,
        filename: str,
        *,
        filters: list[tuple[str, str, Any]] | None = None,
        columns: list[str] | None = None,
    ) -> pd.DataFrame:
        path = self.runtime / filename
        if not path.exists():
            return pd.DataFrame(columns=columns or [])
        return pd.read_parquet(path, columns=columns, filters=filters, engine="pyarrow")

    @lru_cache(maxsize=1)
    def manifest(self) -> dict[str, Any]:
        path = self.runtime / "manifest.json"
        if not path.exists():
            return {
                "available": False,
                "status": "not_built",
                "note_ru": "Фармакологический индекс MCL ещё не построен.",
                "build_command": ".\\.venv\\Scripts\\python.exe .\\scripts\\build_pharmacology_layer.py",
            }
        payload = json.loads(path.read_text(encoding="utf-8"))
        payload["available"] = bool(payload.get("responses_n", 0))
        return payload

    @lru_cache(maxsize=1)
    def concordance_summary(self) -> dict[str, Any]:
        path = self.runtime / "target_concordance_summary.json"
        if not path.exists():
            return {"status": "not_built", "pairs_n": 0}
        return json.loads(path.read_text(encoding="utf-8"))

    @lru_cache(maxsize=1)
    def compounds(self) -> pd.DataFrame:
        path = self.runtime / "compounds.parquet"
        return pd.read_parquet(path) if path.exists() else pd.DataFrame()

    @lru_cache(maxsize=1)
    def target_evidence(self) -> pd.DataFrame:
        path = self.runtime / "target_evidence.parquet"
        frame = pd.read_parquet(path) if path.exists() else pd.DataFrame()
        if "target_gene" in frame.columns:
            frame["target_gene"] = frame["target_gene"].astype(str).str.upper()
        return frame

    @lru_cache(maxsize=1)
    def target_concordance(self) -> pd.DataFrame:
        path = self.runtime / "target_concordance.parquet"
        frame = pd.read_parquet(path) if path.exists() else pd.DataFrame()
        if "target_gene" in frame.columns:
            frame["target_gene"] = frame["target_gene"].astype(str).str.upper()
        return frame

    # Full-table accessors are retained for offline compatibility only. Interactive
    # routes below deliberately use filtered slice methods instead.
    def responses(self) -> pd.DataFrame:
        return self._read_slice("responses.parquet")

    def model_target_links(self) -> pd.DataFrame:
        return self._read_slice("model_compound_target_links.parquet")

    def responses_for_model(self, model_id: str, source: str | None = None) -> pd.DataFrame:
        filters: list[tuple[str, str, Any]] = [("model_id", "==", model_id.strip().upper())]
        if source:
            filters.append(("source", "==", source.strip()))
        frame = self._read_slice("responses.parquet", filters=filters)
        if "model_id" in frame.columns:
            frame["model_id"] = frame["model_id"].astype(str)
        return frame

    def responses_for_compound(self, compound_id: str) -> pd.DataFrame:
        frame = self._read_slice(
            "responses.parquet",
            filters=[("compound_id", "==", str(compound_id))],
        )
        if "model_id" in frame.columns:
            frame["model_id"] = frame["model_id"].astype(str)
        return frame

    def links_for_model(
        self,
        model_id: str,
        compound_ids: list[str] | None = None,
    ) -> pd.DataFrame:
        filters: list[tuple[str, str, Any]] = [("model_id", "==", model_id.strip().upper())]
        if compound_ids:
            filters.append(("compound_id", "in", [str(x) for x in compound_ids]))
        frame = self._read_slice("model_compound_target_links.parquet", filters=filters)
        if "model_id" in frame.columns:
            frame["model_id"] = frame["model_id"].astype(str)
        if "target_gene" in frame.columns:
            frame["target_gene"] = frame["target_gene"].astype(str).str.upper()
        return frame

    def links_for_target(self, target_gene: str) -> pd.DataFrame:
        gene = target_gene.strip().upper()
        frame = self._read_slice(
            "model_compound_target_links.parquet",
            filters=[("target_gene", "==", gene)],
        )
        if "model_id" in frame.columns:
            frame["model_id"] = frame["model_id"].astype(str)
        if "target_gene" in frame.columns:
            frame["target_gene"] = frame["target_gene"].astype(str).str.upper()
        return frame

    def links_for_compound(self, compound_id: str) -> pd.DataFrame:
        frame = self._read_slice(
            "model_compound_target_links.parquet",
            filters=[("compound_id", "==", str(compound_id))],
        )
        if "model_id" in frame.columns:
            frame["model_id"] = frame["model_id"].astype(str)
        if "target_gene" in frame.columns:
            frame["target_gene"] = frame["target_gene"].astype(str).str.upper()
        return frame

    def summary(self) -> dict[str, Any]:
        manifest = dict(self.manifest())
        manifest.setdefault("sources", [])
        manifest["target_concordance"] = self.concordance_summary()
        manifest["contract_note_ru"] = (
            "Ответ клеточной модели на вещество, известная мишень вещества и CRISPR-зависимость "
            "хранятся как независимые типы доказательств. Их согласованность поддерживает гипотезу, "
            "но не доказывает причинный механизм."
        )
        return self._clean(manifest)

    def concordance(
        self,
        *,
        label: str | None = None,
        target_gene: str | None = None,
        compound_id: str | None = None,
        limit: int = 200,
    ) -> dict[str, Any]:
        frame = self.target_concordance().copy()
        if frame.empty:
            return {
                "available": False,
                "status": self.concordance_summary().get("status", "not_built"),
                "total": 0,
                "items": [],
                "note_ru": "Слой согласованности чувствительности препарата с CRISPR-профилем ещё не построен.",
            }
        if label and "concordance_label" in frame.columns:
            frame = frame[frame["concordance_label"].astype(str) == label]
        if target_gene and "target_gene" in frame.columns:
            frame = frame[frame["target_gene"].astype(str).str.upper() == target_gene.strip().upper()]
        if compound_id and "compound_id" in frame.columns:
            frame = frame[frame["compound_id"].astype(str) == str(compound_id)]

        order = {
            "strong_support": 0,
            "supportive": 1,
            "inconclusive": 2,
            "discordant": 3,
            "direction_not_resolved": 4,
            "endpoint_not_resolved": 5,
            "insufficient": 6,
        }
        if "concordance_label" in frame.columns:
            frame["_order"] = frame["concordance_label"].map(order).fillna(99)
            sort_cols = ["_order"]
            ascending = [True]
            if "spearman_rho" in frame.columns:
                sort_cols.append("spearman_rho")
                ascending.append(False)
            if "models_n" in frame.columns:
                sort_cols.append("models_n")
                ascending.append(False)
            frame = frame.sort_values(sort_cols, ascending=ascending, na_position="last").drop(columns=["_order"])
        total = int(len(frame))
        return self._clean(
            {
                "available": True,
                "status": "available",
                "total": total,
                "summary": self.concordance_summary(),
                "items": self._records(frame.head(max(1, min(int(limit), 1000)))),
                "interpretation_ru": (
                    "Для PRISM LFC положительная корреляция с Chronos Gene Effect ожидается для ингибитора, "
                    "если клеточные модели с более сильной CRISPR-зависимостью также сильнее чувствительны к препарату. "
                    "Это поддерживает механизм, но не доказывает прямое target engagement."
                ),
            }
        )

    def model(self, model_id: str, *, limit: int = 100, source: str | None = None) -> dict[str, Any]:
        manifest = self.manifest()
        frame = self.responses_for_model(model_id, source)
        if frame.empty and not manifest.get("available"):
            return self._clean(
                {
                    "model_id": model_id,
                    "available": False,
                    "status": manifest.get("status", "not_built"),
                    "observations_n": 0,
                    "compounds_n": 0,
                    "sources": manifest.get("sources", []),
                    "items": [],
                    "build_command": manifest.get("build_command", ".\\.venv\\Scripts\\python.exe .\\scripts\\build_pharmacology_layer.py"),
                    "note_ru": manifest.get("note_ru", "Фармакологические данные для модели ещё не проиндексированы."),
                }
            )
        if frame.empty:
            return {
                "model_id": model_id,
                "available": True,
                "status": "no_observations_for_model",
                "observations_n": 0,
                "compounds_n": 0,
                "sources": [],
                "items": [],
                "note_ru": "Фармакологический слой построен, но для этой модели пока нет нормализованных наблюдений.",
            }

        total = int(len(frame))
        compounds_n = int(frame["compound_id"].nunique()) if "compound_id" in frame.columns else 0
        sources = sorted(set(frame["source"].dropna().astype(str))) if "source" in frame.columns else []

        compounds = self.compounds()
        if not compounds.empty and "compound_id" in compounds.columns:
            keep = [c for c in (
                "compound_id", "preferred_name", "canonical_smiles", "inchikey", "pubchem_cid",
                "chembl_id", "broad_id", "gdsc_id"
            ) if c in compounds.columns]
            frame = frame.merge(compounds[keep].drop_duplicates("compound_id"), on="compound_id", how="left")

        sort_cols = [c for c in ("source", "compound_id", "endpoint") if c in frame.columns]
        if sort_cols:
            frame = frame.sort_values(sort_cols, na_position="last")
        preview = frame.head(max(1, min(int(limit), 500))).copy()
        preview_compounds = preview["compound_id"].dropna().astype(str).unique().tolist() if "compound_id" in preview.columns else []

        links = self.links_for_model(model_id, preview_compounds)
        link_lookup: dict[str, list[dict[str, Any]]] = {}
        if not links.empty and "compound_id" in links.columns:
            for compound, group in links.groupby("compound_id", sort=False):
                cols = [c for c in (
                    "target_gene", "action", "evidence_type", "confidence", "gene_effect",
                    "crispr_support_level", "source", "activity_type", "activity_value", "activity_unit"
                ) if c in group.columns]
                ranked = group.copy()
                if "gene_effect" in ranked.columns:
                    ranked["_ge"] = pd.to_numeric(ranked["gene_effect"], errors="coerce")
                    ranked = ranked.sort_values("_ge", ascending=True, na_position="last").drop(columns=["_ge"])
                link_lookup[str(compound)] = self._records(ranked[cols].head(20))

        concordance = self.target_concordance()
        concordance_lookup: dict[str, list[dict[str, Any]]] = {}
        if not concordance.empty and "compound_id" in concordance.columns and preview_compounds:
            sub_concordance = concordance[concordance["compound_id"].astype(str).isin(preview_compounds)].copy()
            cols = [c for c in (
                "target_gene", "concordance_label", "spearman_rho", "q_value", "models_n",
                "median_response_delta_dependent_minus_other", "interpretation_ru"
            ) if c in sub_concordance.columns]
            for compound, group in sub_concordance.groupby("compound_id", sort=False):
                concordance_lookup[str(compound)] = self._records(group[cols].head(20))

        items = self._records(preview)
        for item in items:
            compound = str(item.get("compound_id"))
            item["target_hypotheses"] = link_lookup.get(compound, [])
            item["target_concordance"] = concordance_lookup.get(compound, [])

        return self._clean(
            {
                "model_id": model_id,
                "available": True,
                "status": "available",
                "observations_n": total,
                "compounds_n": compounds_n,
                "sources": sources,
                "items": items,
                "interpretation": {
                    "response": "Каждая строка — отдельное экспериментальное наблюдение model × compound; endpoints и единицы не смешиваются в единый рейтинг.",
                    "target": "Мишени происходят из отдельного слоя compound-target evidence.",
                    "crispr": "Отрицательный Gene Effect у аннотированной мишени поддерживает функциональную согласованность, но не доказывает механизм препарата.",
                    "concordance": "Профильная согласованность проверяет связь чувствительности к препарату и CRISPR-зависимости мишени сразу по общей панели моделей.",
                },
            }
        )

    def compound(self, compound_id: str) -> dict[str, Any]:
        compounds = self.compounds()
        if compounds.empty or "compound_id" not in compounds.columns:
            raise MCLDataError("Pharmacology compound index is not available")
        hit = compounds[compounds["compound_id"].astype(str) == str(compound_id)]
        if hit.empty:
            raise MCLDataError(f"Unknown pharmacology compound: {compound_id}")
        identity = self._clean(hit.iloc[0].to_dict())

        targets = self.target_evidence()
        target_rows = targets[targets["compound_id"].astype(str) == str(compound_id)].copy() if not targets.empty else targets
        response_rows = self.responses_for_compound(compound_id)
        concordance = self.target_concordance()
        concordance_rows = concordance[concordance["compound_id"].astype(str) == str(compound_id)].copy() if not concordance.empty and "compound_id" in concordance.columns else pd.DataFrame()
        return self._clean(
            {
                "identity": identity,
                "target_evidence": self._records(target_rows.head(250)) if not target_rows.empty else [],
                "target_concordance": self._records(concordance_rows.head(250)) if not concordance_rows.empty else [],
                "response_summary": {
                    "observations_n": int(len(response_rows)),
                    "models_n": int(response_rows["model_id"].nunique()) if not response_rows.empty and "model_id" in response_rows.columns else 0,
                    "sources": sorted(set(response_rows["source"].dropna().astype(str))) if not response_rows.empty and "source" in response_rows.columns else [],
                },
            }
        )
