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


def _truthy(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and not pd.isna(value):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes", "y", "t"}


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
        gene = str(row.get("GeneSymbol") or row.get("HugoSymbol") or "").strip()
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
                "classification": row.get("VariantClassification") or row.get("MolecularConsequence"),
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

    @lru_cache(maxsize=1)
    def _model_metadata(self) -> pd.DataFrame:
        parquet = self.root / "data/processed/depmap_model_metadata.parquet"
        tsv = self.root / "data/processed/depmap_model_metadata.tsv"
        if parquet.exists():
            return pd.read_parquet(parquet)
        if tsv.exists():
            return pd.read_csv(tsv, sep="\t", low_memory=False)
        return pd.DataFrame()

    @lru_cache(maxsize=1)
    def _model_mutations(self) -> pd.DataFrame:
        parquet = self.root / "data/processed/depmap_model_mutations.parquet"
        tsv = self.root / "data/processed/depmap_model_mutations.tsv"
        if parquet.exists():
            return pd.read_parquet(parquet)
        if tsv.exists():
            return pd.read_csv(tsv, sep="\t", low_memory=False)
        return pd.DataFrame()

    def _comparison_specs(self, cancer_id: str) -> list[dict[str, Any]]:
        return [x for x in self.store.comparison_specs() if str(x.get("cancer_id")) == cancer_id]

    def _patient_layer(self, cancer_id: str, cfg: dict[str, Any], cancer_ru: str) -> dict[str, Any]:
        ot = cfg.get("opentargets") or {}
        return {
            "status": "not_connected",
            "status_ru": "Пациентская когорта пока не подключена",
            "disease_name": cancer_ru,
            "disease_id": ot.get("disease_id") or ot.get("disease_efo_id"),
            "description": (
                "MCL пока не использует частоты мутаций из пациентских когорт для этого экрана. "
                "Поэтому данные клеточных линий ниже нельзя трактовать как частоту изменений у пациентов."
            ),
            "planned_sources": ["cBioPortal / TCGA", "AACR Project GENIE"],
        }

    def _context_model_genetics(self, cancer_id: str, audit_sub: pd.DataFrame) -> dict[str, Any]:
        mutations = self._model_mutations()
        if mutations.empty or "model_id" not in mutations.columns or "gene" not in mutations.columns:
            return {
                "available": False,
                "status_ru": "Полный генетический профиль моделей ещё не проиндексирован",
                "build_command": ".\\.venv\\Scripts\\python.exe .\\scripts\\build_depmap_model_profiles.py",
                "note": "Сейчас доступны только варианты, использованные при формировании молекулярных групп.",
                "top_genes": [],
            }

        model_ids = set(audit_sub.get("model_id", pd.Series(dtype=str)).dropna().astype(str))
        sub = mutations[mutations["model_id"].astype(str).isin(model_ids)].copy()
        if sub.empty:
            return {
                "available": True,
                "status_ru": "Генетический профиль проиндексирован, но вариантов для этих моделей не найдено",
                "top_genes": [],
                "profiled_models_n": 0,
            }

        if "is_functional" in sub.columns:
            functional_mask = sub["is_functional"].map(_truthy)
            functional = sub[functional_mask].copy()
        else:
            impact = sub.get("vep_impact", pd.Series("", index=sub.index)).astype(str).str.upper()
            driver = sub.get("driver", pd.Series(False, index=sub.index)).map(_truthy)
            hotspot = sub.get("hotspot", pd.Series(False, index=sub.index)).map(_truthy)
            lof = sub.get("likely_lof", pd.Series(False, index=sub.index)).map(_truthy)
            functional = sub[impact.isin(["HIGH", "MODERATE"]) | driver | hotspot | lof].copy()

        if functional.empty:
            functional = sub.copy()

        group_map = {}
        if not audit_sub.empty and {"model_id", "assigned_group"}.issubset(audit_sub.columns):
            group_map = dict(
                audit_sub[["model_id", "assigned_group"]]
                .drop_duplicates("model_id")
                .astype(str)
                .itertuples(index=False, name=None)
            )
        functional["assigned_group"] = functional["model_id"].astype(str).map(group_map)

        genes: list[dict[str, Any]] = []
        for gene, gene_rows in functional.groupby(functional["gene"].astype(str)):
            context_n = int(gene_rows.loc[gene_rows["assigned_group"] == "context", "model_id"].nunique())
            comparator_n = int(gene_rows.loc[gene_rows["assigned_group"] == "comparator", "model_id"].nunique())
            models_n = int(gene_rows["model_id"].nunique())
            genes.append(
                {
                    "gene": gene,
                    "models_n": models_n,
                    "context_models_n": context_n,
                    "comparator_models_n": comparator_n,
                }
            )
        genes.sort(key=lambda x: (-int(x["models_n"]), str(x["gene"])))

        profiled_models_n = int(sub["model_id"].nunique())
        sequenced_models_n = int(audit_sub.get("sequencing_available", pd.Series(dtype=bool)).map(_truthy).sum()) if not audit_sub.empty else 0
        for row in genes[:20]:
            denominator = sequenced_models_n or profiled_models_n
            row["model_fraction"] = (row["models_n"] / denominator) if denominator else None

        return {
            "available": True,
            "status_ru": "Полный мутационный профиль клеточных моделей подключён",
            "source": "DepMap OmicsSomaticMutations",
            "profiled_models_n": profiled_models_n,
            "sequenced_models_n": sequenced_models_n,
            "functional_variants_n": int(len(functional)),
            "top_genes": genes[:20],
            "note": "Частоты рассчитаны только среди доступных клеточных моделей MCL и не являются частотами у пациентов.",
        }

    def _context_summary(
        self,
        cancer_id: str,
        cfg: dict[str, Any],
        *,
        include_model_genetics: bool = True,
    ) -> dict[str, Any]:
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
        sequencing_n = int(sequencing.map(_truthy).sum()) if not sequencing.empty else 0
        comparisons = self._comparison_specs(cancer_id)

        molecular = str(ot.get("molecular_context") or inclusion.get("alteration") or cfg.get("name") or cancer_id)
        cancer_ru = display.get("cancer_ru") or subtype or primary_disease
        organ_ru = display.get("organ_ru") or organ_en
        molecular_ru = display.get("molecular_ru") or molecular.replace("p.", "")

        patient_layer = self._patient_layer(cancer_id, cfg, str(cancer_ru))
        if include_model_genetics:
            model_genetics = self._context_model_genetics(cancer_id, sub)
        else:
            model_genetics = {
                "status": "deferred",
                "status_ru": "Подробный мутационный профиль загружается только на странице выбранного контекста",
                "top_genes": [],
            }

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
                "patient_layer": patient_layer,
                "model_genetics": model_genetics,
                "evidence_layers": [
                    {
                        "id": "patient",
                        "label_ru": "Опухоль у пациентов",
                        "status": patient_layer["status"],
                        "source_ru": "Пациентские молекулярные когорты",
                    },
                    {
                        "id": "molecular_context",
                        "label_ru": "Молекулярный подтип",
                        "status": "available",
                        "source_ru": "Конфигурация MCL + мутационные данные DepMap",
                    },
                    {
                        "id": "models",
                        "label_ru": "Экспериментальные модели",
                        "status": "available",
                        "source_ru": "DepMap / OncoTree",
                    },
                ],
            }
        )

    def atlas(self) -> dict[str, Any]:
        contexts = [
            self._context_summary(cid, cfg or {}, include_model_genetics=False)
            for cid, cfg in self._contexts_config().items()
        ]
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
        return self._context_summary(cancer_id, cfg or {}, include_model_genetics=True)

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
            frame = frame[frame["sequencing_available"].map(_truthy)]
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

    def _model_genetics(self, model_id: str, memberships: list[dict[str, Any]]) -> dict[str, Any]:
        mutations = self._model_mutations()
        if not mutations.empty and "model_id" in mutations.columns:
            rows = mutations[mutations["model_id"].astype(str).str.upper() == model_id.upper()].copy()
            if not rows.empty:
                if "priority_score" in rows.columns:
                    rows = rows.sort_values(["priority_score", "gene"], ascending=[False, True], na_position="last")
                elif "gene" in rows.columns:
                    rows = rows.sort_values("gene")

                driver = rows.get("driver", pd.Series(False, index=rows.index)).map(_truthy)
                hotspot = rows.get("hotspot", pd.Series(False, index=rows.index)).map(_truthy)
                lof = rows.get("likely_lof", pd.Series(False, index=rows.index)).map(_truthy)
                impact = rows.get("vep_impact", pd.Series("", index=rows.index)).astype(str).str.upper()
                priority_mask = driver | hotspot | lof | impact.eq("HIGH")
                priority = rows[priority_mask].copy()
                if priority.empty:
                    priority = rows.head(25).copy()

                columns = [
                    "gene", "protein_change", "dna_change", "variant_type", "molecular_consequence",
                    "vep_impact", "driver", "hotspot", "likely_lof", "allele_fraction", "depth",
                    "clin_sig", "civic_description", "priority_score",
                ]
                keep = [c for c in columns if c in rows.columns]
                return {
                    "availability": "full",
                    "availability_ru": "Полный мутационный профиль подключён",
                    "source": "DepMap OmicsSomaticMutations",
                    "mutations_n": int(len(rows)),
                    "mutated_genes_n": int(rows["gene"].dropna().astype(str).nunique()) if "gene" in rows else 0,
                    "driver_mutations_n": int(driver.sum()),
                    "hotspot_mutations_n": int(hotspot.sum()),
                    "likely_lof_n": int(lof.sum()),
                    "high_impact_n": int(impact.eq("HIGH").sum()),
                    "priority_variants": _clean(priority[keep].head(50).to_dict("records")),
                    "variants": _clean(rows[keep].head(300).to_dict("records")),
                    "variants_returned_n": min(int(len(rows)), 300),
                    "note": "Это генетика конкретной экспериментальной модели, а не частоты изменений в опухолях пациентов.",
                }

        fallback: list[dict[str, Any]] = []
        seen: set[tuple[Any, Any, Any]] = set()
        for membership in memberships:
            for variant in membership.get("variants") or []:
                key = (variant.get("gene"), variant.get("protein_change"), variant.get("dna_change"))
                if key in seen:
                    continue
                seen.add(key)
                fallback.append(variant)
        return {
            "availability": "limited",
            "availability_ru": "Доступны только варианты, использованные при формировании групп",
            "source": "MCL DepMap context audit",
            "mutations_n": len(fallback),
            "mutated_genes_n": len({x.get("gene") for x in fallback if x.get("gene")}),
            "driver_mutations_n": sum(1 for x in fallback if _truthy(x.get("driver"))),
            "hotspot_mutations_n": sum(1 for x in fallback if _truthy(x.get("hotspot"))),
            "likely_lof_n": 0,
            "high_impact_n": 0,
            "priority_variants": fallback,
            "variants": fallback,
            "variants_returned_n": len(fallback),
            "build_command": ".\\.venv\\Scripts\\python.exe .\\scripts\\build_depmap_model_profiles.py",
            "note": "Для полной генетики один раз постройте компактный индекс из локального OmicsSomaticMutations.csv.",
        }

    def model(self, model_id: str) -> dict[str, Any]:
        frame = self._audit()
        if frame.empty or "model_id" not in frame:
            raise MCLDataError(f"Cell model not found: {model_id}")
        hits = frame[frame["model_id"].astype(str).str.upper() == model_id.strip().upper()]
        if hits.empty:
            raise MCLDataError(f"Cell model not found: {model_id}")
        row = hits.iloc[0]
        resolved_model_id = str(row.get("model_id"))
        memberships = self.models(search=resolved_model_id, limit=100)["items"]

        metadata: dict[str, Any] = {}
        metadata_frame = self._model_metadata()
        if not metadata_frame.empty and "model_id" in metadata_frame.columns:
            meta_hit = metadata_frame[metadata_frame["model_id"].astype(str).str.upper() == resolved_model_id.upper()]
            if not meta_hit.empty:
                metadata = _clean(meta_hit.iloc[0].to_dict())

        genetics = self._model_genetics(resolved_model_id, memberships)

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
                "metadata": metadata,
                "genetics": genetics,
            }
        )
