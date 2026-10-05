from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

from .store import MCLDataError


class MCLPharmacologyCatalogStore:
    """Fast read-only indexes for compound and gene-mapped protein target pages.

    Catalog pages use compact precomputed indexes. Detail pages request only filtered
    slices from the large response/link Parquet files through MCLPharmacologyStore.
    """

    def __init__(self, root: Path, base_store: Any):
        self.root = Path(root).resolve()
        self.runtime = self.root / "data" / "runtime" / "pharmacology"
        self.processed = self.root / "data" / "processed"
        self.base = base_store

    @staticmethod
    def _clean(value: Any) -> Any:
        if isinstance(value, dict):
            return {str(k): MCLPharmacologyCatalogStore._clean(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [MCLPharmacologyCatalogStore._clean(v) for v in value]
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

    @classmethod
    def _records(cls, frame: pd.DataFrame) -> list[dict[str, Any]]:
        return [cls._clean(x) for x in frame.to_dict("records")]

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

    @lru_cache(maxsize=1)
    def compounds(self) -> pd.DataFrame:
        path = self.runtime / "compound_catalog.parquet"
        return pd.read_parquet(path) if path.exists() else pd.DataFrame()

    @lru_cache(maxsize=1)
    def targets(self) -> pd.DataFrame:
        path = self.runtime / "target_catalog.parquet"
        frame = pd.read_parquet(path) if path.exists() else pd.DataFrame()
        if "target_gene" in frame.columns:
            frame["target_gene"] = frame["target_gene"].astype(str).str.upper()
        return frame

    @lru_cache(maxsize=1)
    def target_compounds(self) -> pd.DataFrame:
        path = self.runtime / "target_compound_catalog.parquet"
        frame = pd.read_parquet(path) if path.exists() else pd.DataFrame()
        if "target_gene" in frame.columns:
            frame["target_gene"] = frame["target_gene"].astype(str).str.upper()
        return frame

    @lru_cache(maxsize=1)
    def atlas(self) -> pd.DataFrame:
        path = self.processed / "depmap_crispr_model_atlas.parquet"
        if not path.exists():
            return pd.DataFrame()
        keep = [
            "model_id", "cell_line_name", "mcl_cancer_name", "mcl_organ_ru", "mcl_system_ru"
        ]
        try:
            frame = pd.read_parquet(path, columns=keep)
        except Exception:
            frame = pd.read_parquet(path)
            frame = frame[[c for c in keep if c in frame.columns]]
        if "model_id" in frame.columns:
            frame["model_id"] = frame["model_id"].astype(str)
        return frame

    def compound_search(
        self,
        *,
        search: str | None = None,
        target_gene: str | None = None,
        has_smiles: bool = False,
        limit: int = 100,
        offset: int = 0,
    ) -> dict[str, Any]:
        frame = self.compounds().copy()
        if frame.empty:
            return {
                "available": False,
                "status": "catalog_not_built",
                "total": 0,
                "items": [],
                "build_command": ".\\.venv\\Scripts\\python.exe .\\scripts\\build_pharmacology_layer.py",
            }

        if search:
            needle = search.strip().casefold()
            mask = pd.Series(False, index=frame.index)
            for column in (
                "preferred_name", "compound_id", "canonical_smiles", "inchikey",
                "chembl_id", "broad_id", "pubchem_cid",
            ):
                if column in frame.columns:
                    mask = mask | frame[column].fillna("").astype(str).str.casefold().str.contains(
                        needle, regex=False
                    )
            frame = frame[mask]

        if target_gene and "target_genes_json" in frame.columns:
            gene = target_gene.strip().upper()
            frame = frame[
                frame["target_genes_json"].map(
                    lambda value: gene in {item.upper() for item in self._json_list(value)}
                )
            ]

        if has_smiles and "canonical_smiles" in frame.columns:
            frame = frame[frame["canonical_smiles"].fillna("").astype(str).str.strip().ne("")]

        sort_cols = [c for c in ("models_n", "observations_n", "preferred_name") if c in frame.columns]
        if sort_cols:
            frame = frame.sort_values(
                sort_cols,
                ascending=[False if c in {"models_n", "observations_n"} else True for c in sort_cols],
                na_position="last",
            )

        total = int(len(frame))
        start = max(0, int(offset))
        page = frame.iloc[start : start + max(1, min(int(limit), 500))].copy()
        items = self._records(page)
        for item in items:
            item["target_genes"] = self._json_list(item.pop("target_genes_json", "[]"))
            item["sources"] = self._json_list(item.pop("sources_json", "[]"))
            item["endpoints"] = self._json_list(item.pop("endpoints_json", "[]"))

        return self._clean(
            {
                "available": True,
                "total": total,
                "limit": int(limit),
                "offset": start,
                "items": items,
                "note_ru": (
                    "Одна строка соответствует внутреннему compound_id. До химической дедупликации "
                    "разные PRISM treatment-profile одного вещества могут оставаться отдельными записями."
                ),
            }
        )

    def target_search(
        self,
        *,
        search: str | None = None,
        limit: int = 100,
        offset: int = 0,
    ) -> dict[str, Any]:
        frame = self.targets().copy()
        if frame.empty:
            return {
                "available": False,
                "status": "catalog_not_built",
                "total": 0,
                "items": [],
                "build_command": ".\\.venv\\Scripts\\python.exe .\\scripts\\build_pharmacology_layer.py",
            }

        if search:
            needle = search.strip().upper()
            frame = frame[frame["target_gene"].astype(str).str.upper().str.contains(needle, regex=False)]

        sort_cols = [c for c in ("compounds_n", "models_n", "target_gene") if c in frame.columns]
        if sort_cols:
            frame = frame.sort_values(
                sort_cols,
                ascending=[False if c in {"compounds_n", "models_n"} else True for c in sort_cols],
                na_position="last",
            )

        total = int(len(frame))
        start = max(0, int(offset))
        page = frame.iloc[start : start + max(1, min(int(limit), 500))].copy()
        items = self._records(page)
        for item in items:
            item["actions"] = self._json_list(item.pop("actions_json", "[]"))
            item["evidence_types"] = self._json_list(item.pop("evidence_types_json", "[]"))
            crispr_n = int(item.get("crispr_models_n") or 0)
            dependent_n = int(item.get("dependent_models_n") or 0)
            item["dependency_fraction"] = dependent_n / crispr_n if crispr_n else None

        return self._clean(
            {
                "available": True,
                "total": total,
                "limit": int(limit),
                "offset": start,
                "items": items,
                "note_ru": (
                    "Текущий каталог мишеней нормализован до кодирующего гена. Конкретная белковая "
                    "изоформа или белковый комплекс не назначаются без отдельного источника."
                ),
            }
        )

    def compound_detail(self, compound_id: str, *, limit: int = 60) -> dict[str, Any]:
        catalog = self.compounds()
        if catalog.empty or "compound_id" not in catalog.columns:
            raise MCLDataError("Pharmacology compound catalog is not available")
        hit = catalog[catalog["compound_id"].astype(str) == str(compound_id)]
        if hit.empty:
            raise MCLDataError(f"Unknown pharmacology compound: {compound_id}")

        identity = self._clean(hit.iloc[0].to_dict())
        identity["target_genes"] = self._json_list(identity.pop("target_genes_json", "[]"))
        identity["sources"] = self._json_list(identity.pop("sources_json", "[]"))
        identity["endpoints"] = self._json_list(identity.pop("endpoints_json", "[]"))

        evidence = self.base.target_evidence()
        target_rows = (
            evidence[evidence["compound_id"].astype(str) == str(compound_id)].copy()
            if not evidence.empty else pd.DataFrame()
        )
        target_summary: list[dict[str, Any]] = []
        if not target_rows.empty:
            for gene, group in target_rows.groupby("target_gene", sort=True):
                target_summary.append(
                    {
                        "target_gene": str(gene),
                        "evidence_rows_n": int(len(group)),
                        "actions": sorted({str(x) for x in group.get("action", pd.Series(dtype=object)).dropna() if str(x).strip()}),
                        "evidence_types": sorted({str(x) for x in group.get("evidence_type", pd.Series(dtype=object)).dropna() if str(x).strip()}),
                        "sources": sorted({str(x) for x in group.get("source", pd.Series(dtype=object)).dropna() if str(x).strip()}),
                    }
                )

        # Critical performance rule: read only this compound from the 3.46M-row table.
        response_rows = self.base.responses_for_compound(str(compound_id))
        if not response_rows.empty:
            endpoints = set(response_rows["endpoint"].dropna().astype(str).str.upper()) if "endpoint" in response_rows.columns else set()
            if endpoints == {"LFC"} and "value" in response_rows.columns:
                response_rows["value"] = pd.to_numeric(response_rows["value"], errors="coerce")
                response_rows = response_rows.sort_values("value", ascending=True, na_position="last")
            atlas = self.atlas()
            if not atlas.empty:
                response_rows = response_rows.merge(atlas.drop_duplicates("model_id"), on="model_id", how="left")

        example_cols = [c for c in (
            "model_id", "cell_line_name", "mcl_cancer_name", "mcl_organ_ru", "source", "endpoint",
            "value", "unit", "dose", "dose_unit", "exposure_time_h", "quality_flag",
        ) if c in response_rows.columns]

        concordance = self.base.concordance(compound_id=str(compound_id), limit=100)
        return self._clean(
            {
                "identity": identity,
                "targets": target_summary,
                "response_summary": {
                    "observations_n": int(len(response_rows)),
                    "models_n": int(response_rows["model_id"].nunique()) if not response_rows.empty and "model_id" in response_rows.columns else 0,
                },
                "response_examples": self._records(
                    response_rows[example_cols].head(max(1, min(int(limit), 200)))
                ) if not response_rows.empty else [],
                "target_concordance": concordance,
                "interpretation_ru": (
                    "Чувствительность модели к веществу и аннотация белковой мишени являются разными "
                    "наблюдениями. Даже совпадение с CRISPR-зависимостью поддерживает механизм, но не доказывает его причинно."
                ),
            }
        )

    def target_detail(self, target_id: str, *, limit: int = 100) -> dict[str, Any]:
        gene = target_id.strip().upper()
        catalog = self.targets()
        if catalog.empty or "target_gene" not in catalog.columns:
            raise MCLDataError("Pharmacology target catalog is not available")
        hit = catalog[catalog["target_gene"].astype(str).str.upper() == gene]
        if hit.empty:
            raise MCLDataError(f"Unknown pharmacology target: {gene}")

        identity = self._clean(hit.iloc[0].to_dict())
        identity["actions"] = self._json_list(identity.pop("actions_json", "[]"))
        identity["evidence_types"] = self._json_list(identity.pop("evidence_types_json", "[]"))
        crispr_n = int(identity.get("crispr_models_n") or 0)
        dependent_n = int(identity.get("dependent_models_n") or 0)
        identity["dependency_fraction"] = dependent_n / crispr_n if crispr_n else None

        compounds = self.target_compounds()
        compound_rows = (
            compounds[compounds["target_gene"].astype(str).str.upper() == gene].copy()
            if not compounds.empty else pd.DataFrame()
        )
        if not compound_rows.empty:
            sort_cols = [c for c in ("models_n", "dependent_models_n", "preferred_name") if c in compound_rows.columns]
            if sort_cols:
                compound_rows = compound_rows.sort_values(
                    sort_cols,
                    ascending=[False if c in {"models_n", "dependent_models_n"} else True for c in sort_cols],
                    na_position="last",
                )
            for column in ("actions_json", "evidence_types_json"):
                if column in compound_rows.columns:
                    compound_rows[column.replace("_json", "")] = compound_rows[column].map(self._json_list)
            compound_rows = compound_rows.drop(columns=["actions_json", "evidence_types_json"], errors="ignore")

        # Critical performance rule: read only this target from the 6.03M-row link table.
        model_rows = self.base.links_for_target(gene)
        if not model_rows.empty:
            model_rows["gene_effect"] = pd.to_numeric(model_rows["gene_effect"], errors="coerce")
            model_rows = model_rows.sort_values("gene_effect", ascending=True, na_position="last")
            model_rows = model_rows.drop_duplicates("model_id", keep="first")
            atlas = self.atlas()
            if not atlas.empty:
                model_rows = model_rows.merge(atlas.drop_duplicates("model_id"), on="model_id", how="left")

        model_cols = [c for c in (
            "model_id", "cell_line_name", "mcl_cancer_name", "mcl_organ_ru", "gene_effect",
            "crispr_support_level",
        ) if c in model_rows.columns]

        concordance = self.base.concordance(target_gene=gene, limit=100)
        return self._clean(
            {
                "identity": identity,
                "coding_gene": gene,
                "resolution": "gene_mapped_protein_target",
                "resolution_note_ru": (
                    "Мишень связана с кодирующим геном по аннотации источника. Конкретная изоформа белка, "
                    "протеоформа или белковый комплекс требуют отдельного подтверждения."
                ),
                "compounds": self._records(
                    compound_rows.head(max(1, min(int(limit), 250)))
                ) if not compound_rows.empty else [],
                "model_examples": self._records(model_rows[model_cols].head(40)) if not model_rows.empty else [],
                "target_concordance": concordance,
            }
        )
