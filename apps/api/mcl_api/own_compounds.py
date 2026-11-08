from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

from .store import MCLDataError


SERIES_RULES = [
    (1, 5, "ацетилацетонат", "PyrNOx"),
    (6, 10, "димедон", "PyrNOx"),
    (11, 15, "барбитурат", "PyrNOx"),
    (16, 20, "антипиринат", "PyrNOx"),
    (21, 25, "ацетилацетонат", "Pir"),
    (26, 30, "димедон", "Pir"),
    (31, 35, "антипиринат", "Pir"),
    (36, 40, "барбитурат", "Pir"),
    (41, 45, "N,N-диметилбарбитурат", "PyrNOx"),
    (46, 50, "циклогександион", "PyrNOx"),
    (51, 55, "фенилметилпиразолонат", "PyrNOx"),
    (56, 60, "4-гидроксикумарин", "PyrNOx"),
    (61, 65, "малонодинитрил", "PyrNOx"),
]

ROLE_ORDER = {
    "positive": 0,
    "negative_same_cancer": 1,
    "negative_panel_comparator": 1,
    "discordant_sensitive_without_dependency": 2,
    "discordant_dependency_without_activity": 2,
}


class OwnCompoundStore:
    """Read-only explorer for the 65 PYZ compounds and their MCL evidence.

    The store deliberately separates evidence about PYZ from evidence about known
    ligands and from target/cancer context. Target-relevant cell lines are proposed
    for testing a PYZ-target hypothesis; they are not presented as cell lines on
    which the PYZ compound has already been shown to work.
    """

    def __init__(self, root: Path):
        self.root = Path(root).resolve()
        self.registry_path = self.root / "config" / "own_compounds.tsv"
        self.admet_path = self.root / "config" / "pyz_admet_snapshot.tsv"
        self.admet_summary_path = self.root / "config" / "pyz_admet_summary.json"
        self.matrix_path = self.root / "data" / "runtime" / "pyz_target_evidence" / "pyz_target_evidence_matrix.parquet"
        self.matrix_tsv = self.root / "data" / "runtime" / "pyz_target_evidence" / "pyz_target_evidence_matrix.tsv"
        self.top_hits_path = self.root / "data" / "runtime" / "own_compounds" / "target_ligand_space_similarity_top_hits.parquet"
        self.top_hits_tsv = self.root / "data" / "runtime" / "own_compounds" / "target_ligand_space_similarity_top_hits.tsv"
        self.candidate_path = self.root / "data" / "runtime" / "pharmacology" / "candidate_hypotheses_v2.parquet"
        self.candidate_models_path = self.root / "data" / "runtime" / "pharmacology" / "candidate_hypothesis_models_v2.parquet"
        self.lab_path = self.root / "data" / "processed" / "laboratory_panel.parquet"

    @staticmethod
    def _clean(value: Any) -> Any:
        if isinstance(value, dict):
            return {str(k): OwnCompoundStore._clean(v) for k, v in value.items()}
        if isinstance(value, (list, tuple, set)):
            return [OwnCompoundStore._clean(v) for v in value]
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
        return [cls._clean(row) for row in frame.to_dict("records")]

    @staticmethod
    def _json_list(value: Any) -> list[str]:
        if isinstance(value, list):
            return [str(v) for v in value]
        if value is None:
            return []
        try:
            parsed = json.loads(str(value))
        except (json.JSONDecodeError, TypeError, ValueError):
            return []
        return [str(v) for v in parsed] if isinstance(parsed, list) else []

    @staticmethod
    def _series_meta(compound_id: str) -> dict[str, Any]:
        try:
            number = int(str(compound_id).split("-")[-1])
        except (TypeError, ValueError):
            return {"tier": None, "chemotype_ru": "—", "ring_variant": "—"}
        for start, end, chemotype, ring in SERIES_RULES:
            if start <= number <= end:
                return {
                    "tier": 1 if number <= 40 else 2,
                    "chemotype_ru": chemotype,
                    "ring_variant": ring,
                    "series_range": f"PYZ-{start:03d}…PYZ-{end:03d}",
                }
        return {"tier": 1 if number <= 40 else 2, "chemotype_ru": "—", "ring_variant": "—"}

    @staticmethod
    def _read_table(parquet: Path, tsv: Path | None = None) -> pd.DataFrame:
        if parquet.exists():
            return pd.read_parquet(parquet)
        if tsv is not None and tsv.exists():
            return pd.read_csv(tsv, sep="\t")
        return pd.DataFrame()

    @lru_cache(maxsize=1)
    def registry(self) -> pd.DataFrame:
        if not self.registry_path.exists():
            return pd.DataFrame()
        frame = pd.read_csv(self.registry_path, sep="\t").copy()
        meta = frame["own_compound_id"].map(self._series_meta).apply(pd.Series)
        frame = pd.concat([frame, meta], axis=1)
        if self.admet_path.exists():
            admet = pd.read_csv(self.admet_path, sep="\t")
            frame = frame.merge(admet, on="own_compound_id", how="left")
        return frame

    @lru_cache(maxsize=1)
    def evidence_matrix(self) -> pd.DataFrame:
        frame = self._read_table(self.matrix_path, self.matrix_tsv)
        if not frame.empty:
            frame = frame.copy()
            frame["own_compound_id"] = frame["own_compound_id"].astype(str)
            frame["target_gene"] = frame["target_gene"].astype(str).str.upper()
            if "best_tanimoto_morgan_r2_2048" in frame.columns:
                frame["best_tanimoto_morgan_r2_2048"] = pd.to_numeric(
                    frame["best_tanimoto_morgan_r2_2048"], errors="coerce"
                )
        return frame

    @lru_cache(maxsize=1)
    def top_hits(self) -> pd.DataFrame:
        frame = self._read_table(self.top_hits_path, self.top_hits_tsv)
        if not frame.empty:
            frame = frame.copy()
            frame["own_compound_id"] = frame["own_compound_id"].astype(str)
            frame["target_gene"] = frame["target_gene"].astype(str).str.upper()
        return frame

    @lru_cache(maxsize=1)
    def candidate_hypotheses(self) -> pd.DataFrame:
        if not self.candidate_path.exists():
            return pd.DataFrame()
        frame = pd.read_parquet(self.candidate_path).copy()
        if "target_gene" in frame.columns:
            frame["target_gene"] = frame["target_gene"].astype(str).str.upper()
        return frame

    @lru_cache(maxsize=1)
    def candidate_models(self) -> pd.DataFrame:
        if not self.candidate_models_path.exists():
            return pd.DataFrame()
        return pd.read_parquet(self.candidate_models_path)

    @lru_cache(maxsize=1)
    def lab_panel(self) -> pd.DataFrame:
        if not self.lab_path.exists():
            return pd.DataFrame()
        return pd.read_parquet(self.lab_path)

    def _admet_summary(self) -> dict[str, Any]:
        if not self.admet_summary_path.exists():
            return {"available": False}
        payload = json.loads(self.admet_summary_path.read_text(encoding="utf-8"))
        payload["available"] = True
        return payload

    def summary(self) -> dict[str, Any]:
        registry = self.registry()
        matrix = self.evidence_matrix()
        shortlist = matrix[
            pd.to_numeric(matrix.get("best_tanimoto_morgan_r2_2048"), errors="coerce").ge(0.35)
        ] if not matrix.empty else pd.DataFrame()
        transfer_counts = (
            matrix["structural_transfer_class"].value_counts().to_dict()
            if not matrix.empty and "structural_transfer_class" in matrix.columns else {}
        )
        targets = sorted(matrix["target_gene"].dropna().astype(str).unique().tolist()) if not matrix.empty else []
        return self._clean({
            "available": not registry.empty,
            "evidence_matrix_available": not matrix.empty,
            "compounds_n": int(len(registry)),
            "targets_n": int(len(targets)),
            "target_genes": targets,
            "matrix_rows_n": int(len(matrix)),
            "shortlist_rows_n": int(len(shortlist)),
            "shortlist_compounds_n": int(shortlist["own_compound_id"].nunique()) if not shortlist.empty else 0,
            "shortlist_targets": sorted(shortlist["target_gene"].dropna().astype(str).unique().tolist()) if not shortlist.empty else [],
            "structural_transfer_class_counts": transfer_counts,
            "admet": self._admet_summary(),
            "build_command": ".\\scripts\\build-pyz-target-evidence-matrix.ps1",
            "interpretation_ru": (
                "PYZ Explorer связывает собственную молекулу с target-based контекстом опухоли и клеточных моделей. "
                "Клеточная релевантность означает 'где рационально проверить гипотезу', а не доказанную активность PYZ на этой культуре."
            ),
        })

    def catalog(
        self,
        *,
        search: str | None = None,
        target_gene: str | None = None,
        shortlist_only: bool = False,
        limit: int = 100,
        offset: int = 0,
    ) -> dict[str, Any]:
        registry = self.registry().copy()
        if registry.empty:
            return {"available": False, "total": 0, "items": []}
        matrix = self.evidence_matrix()
        if not matrix.empty:
            best = matrix.sort_values(
                ["own_compound_id", "best_tanimoto_morgan_r2_2048"], ascending=[True, False], na_position="last"
            ).drop_duplicates("own_compound_id")
            best = best[[c for c in [
                "own_compound_id", "target_gene", "best_tanimoto_morgan_r2_2048", "best_ligand_id",
                "best_lowest_reported_value_nm_across_endpoints", "best_lowest_value_endpoint_type",
                "structural_transfer_class"
            ] if c in best.columns]].rename(columns={
                "target_gene": "best_target_gene",
                "best_tanimoto_morgan_r2_2048": "best_tanimoto",
            })
            registry = registry.merge(best, on="own_compound_id", how="left")
            eligible = matrix[matrix["best_tanimoto_morgan_r2_2048"].ge(0.35)]
            grouped = eligible.groupby("own_compound_id")["target_gene"].agg(list).to_dict()
            registry["shortlist_target_genes"] = registry["own_compound_id"].map(
                lambda x: sorted(set(grouped.get(x, [])))
            )
            registry["shortlist_targets_n"] = registry["shortlist_target_genes"].map(len)
        else:
            registry["shortlist_target_genes"] = [[] for _ in range(len(registry))]
            registry["shortlist_targets_n"] = 0

        if search:
            needle = search.strip().casefold()
            mask = pd.Series(False, index=registry.index)
            for col in ("own_compound_id", "preferred_name", "standardized_smiles", "chemotype_ru", "ring_variant"):
                if col in registry.columns:
                    mask |= registry[col].fillna("").astype(str).str.casefold().str.contains(needle, regex=False)
            registry = registry[mask]
        if target_gene:
            gene = target_gene.strip().upper()
            registry = registry[registry["shortlist_target_genes"].map(lambda xs: gene in xs)]
        if shortlist_only:
            registry = registry[registry["shortlist_targets_n"].gt(0)]

        registry["_priority"] = registry["shortlist_targets_n"].fillna(0).astype(int)
        registry["_best"] = pd.to_numeric(registry.get("best_tanimoto"), errors="coerce")
        registry = registry.sort_values(["_priority", "_best", "own_compound_id"], ascending=[False, False, True])
        registry = registry.drop(columns=["_priority", "_best"], errors="ignore")
        total = int(len(registry))
        start = max(0, int(offset))
        page = registry.iloc[start:start + max(1, min(int(limit), 500))]
        return self._clean({
            "available": True,
            "total": total,
            "limit": int(limit),
            "offset": start,
            "items": self._records(page),
        })

    def matrix_view(
        self,
        *,
        search: str | None = None,
        min_tanimoto: float = 0.0,
        shortlist_only: bool = False,
    ) -> dict[str, Any]:
        matrix = self.evidence_matrix().copy()
        registry = self.registry()
        if matrix.empty:
            return {"available": False, "targets": [], "compounds": [], "items": []}
        matrix = matrix[matrix["best_tanimoto_morgan_r2_2048"].fillna(-1).ge(float(min_tanimoto))]
        if shortlist_only:
            matrix = matrix[matrix["best_tanimoto_morgan_r2_2048"].ge(0.35)]
        if search:
            needle = search.strip().casefold()
            matrix = matrix[
                matrix["own_compound_id"].astype(str).str.casefold().str.contains(needle, regex=False)
                | matrix["target_gene"].astype(str).str.casefold().str.contains(needle, regex=False)
            ]
        targets = sorted(self.evidence_matrix()["target_gene"].dropna().astype(str).unique().tolist())
        compound_meta = registry[[c for c in ["own_compound_id", "chemotype_ru", "ring_variant", "tier", "integrated_risk", "qed"] if c in registry.columns]]
        return self._clean({
            "available": True,
            "targets": targets,
            "compounds": self._records(compound_meta),
            "items": self._records(matrix),
            "interpretation_ru": "Цвет/значение ячейки — только 2D Tanimoto до лучшего экспериментального лиганда этой мишени; это не вероятность связывания PYZ.",
        })

    @staticmethod
    def _first_existing(columns: list[str], frame: pd.DataFrame) -> str | None:
        return next((c for c in columns if c in frame.columns), None)

    @staticmethod
    def _first_value(series: pd.Series) -> Any:
        values = series.dropna()
        return values.iloc[0] if not values.empty else None

    def _contexts_for_target(self, target_gene: str) -> list[dict[str, Any]]:
        hypotheses = self.candidate_hypotheses()
        if hypotheses.empty or "target_gene" not in hypotheses.columns:
            return []
        gene = target_gene.strip().upper()
        target = hypotheses[hypotheses["target_gene"].astype(str).str.upper().eq(gene)].copy()
        if "priority_status" in target.columns:
            target = target[target["priority_status"].astype(str).eq("priority_for_in_vitro")]
        if target.empty:
            return []

        group_cols = [c for c in ["mcl_organ_ru", "mcl_system_ru", "mcl_cancer_id", "mcl_cancer_name"] if c in target.columns]
        if not group_cols:
            return []
        models = self.candidate_models()
        lab = self.lab_panel()
        lab_by_model: dict[str, list[dict[str, Any]]] = {}
        if not lab.empty and "model_id" in lab.columns:
            for model_id, group in lab[lab["model_id"].notna()].groupby("model_id"):
                lab_by_model[str(model_id)] = self._records(group)

        contexts: list[dict[str, Any]] = []
        grouped = target.groupby(group_cols, dropna=False, sort=False)
        for key, group in grouped:
            key_values = key if isinstance(key, tuple) else (key,)
            meta = dict(zip(group_cols, key_values))
            ids = set(group["hypothesis_id"].dropna().astype(str).tolist()) if "hypothesis_id" in group.columns else set()
            context_models = pd.DataFrame()
            if ids and not models.empty and "hypothesis_id" in models.columns:
                context_models = models[models["hypothesis_id"].astype(str).isin(ids)].copy()
                join_cols = [c for c in ["hypothesis_id", "compound_id", "preferred_name"] if c in group.columns]
                if join_cols and "hypothesis_id" in join_cols:
                    mapping = group[join_cols].drop_duplicates("hypothesis_id")
                    extra = [c for c in join_cols if c != "hypothesis_id" and c not in context_models.columns]
                    if extra:
                        context_models = context_models.merge(mapping[["hypothesis_id", *extra]], on="hypothesis_id", how="left")

            model_items: list[dict[str, Any]] = []
            if not context_models.empty:
                model_id_col = self._first_existing(["model_id", "ModelID", "depmap_model_id"], context_models)
                name_col = self._first_existing(["cell_line_name", "model_name", "depmap_cell_line_name", "stripped_cell_line_name"], context_models)
                role_col = self._first_existing(["model_role", "laboratory_model_role"], context_models)
                dep_col = self._first_existing(["dependency_probability", "probability_of_dependency"], context_models)
                ge_col = self._first_existing(["gene_effect", "crispr_gene_effect"], context_models)
                if model_id_col:
                    for model_id, model_group in context_models.groupby(model_id_col, dropna=False, sort=False):
                        model_key = str(model_id) if pd.notna(model_id) else ""
                        roles = sorted(set(model_group[role_col].dropna().astype(str).tolist())) if role_col else []
                        role = min(roles, key=lambda x: ROLE_ORDER.get(x, 9)) if roles else "unclassified"
                        known_compounds_col = self._first_existing(["preferred_name", "compound_id"], model_group)
                        known_compounds = sorted(set(model_group[known_compounds_col].dropna().astype(str).tolist())) if known_compounds_col else []
                        lab_rows = lab_by_model.get(model_key, [])
                        item = {
                            "model_id": model_key or None,
                            "cell_line_name": self._first_value(model_group[name_col]) if name_col else None,
                            "role": role,
                            "roles": roles,
                            "dependency_probability": pd.to_numeric(model_group[dep_col], errors="coerce").max() if dep_col else None,
                            "gene_effect": pd.to_numeric(model_group[ge_col], errors="coerce").min() if ge_col else None,
                            "known_reference_compounds": known_compounds[:8],
                            "known_reference_compounds_n": len(known_compounds),
                            "available_in_laboratory": bool(lab_rows),
                            "laboratory_lines": [row.get("lab_name") for row in lab_rows if row.get("lab_name")],
                            "laboratory_roles": [row.get("laboratory_role") for row in lab_rows if row.get("laboratory_role")],
                        }
                        model_items.append(self._clean(item))
                model_items.sort(key=lambda x: (0 if x.get("available_in_laboratory") else 1, ROLE_ORDER.get(str(x.get("role")), 9), str(x.get("cell_line_name") or x.get("model_id"))))

            positive_n = sum(1 for item in model_items if item.get("role") == "positive")
            negative_n = sum(1 for item in model_items if str(item.get("role", "")).startswith("negative_"))
            lab_n = sum(1 for item in model_items if item.get("available_in_laboratory"))
            contexts.append(self._clean({
                **meta,
                "target_gene": gene,
                "priority_hypotheses_n": int(len(group)),
                "positive_models_n": positive_n,
                "negative_models_n": negative_n,
                "laboratory_models_n": lab_n,
                "models_n": len(model_items),
                "models": model_items[:30],
            }))

        contexts.sort(key=lambda x: (-int(x.get("laboratory_models_n") or 0), -int(x.get("positive_models_n") or 0), -int(x.get("priority_hypotheses_n") or 0), str(x.get("mcl_organ_ru") or ""), str(x.get("mcl_cancer_name") or "")))
        return contexts

    def detail(self, compound_id: str, *, focus_target: str | None = None) -> dict[str, Any]:
        compound_id = compound_id.strip().upper()
        registry = self.registry()
        if registry.empty:
            raise MCLDataError("PYZ registry is not available")
        hit = registry[registry["own_compound_id"].astype(str).str.upper().eq(compound_id)]
        if hit.empty:
            raise MCLDataError(f"Unknown PYZ compound: {compound_id}")
        compound = self._records(hit.head(1))[0]

        matrix = self.evidence_matrix()
        evidence = matrix[matrix["own_compound_id"].astype(str).str.upper().eq(compound_id)].copy() if not matrix.empty else pd.DataFrame()
        if not evidence.empty:
            evidence = evidence.sort_values("best_tanimoto_morgan_r2_2048", ascending=False, na_position="last")
        top_hits = self.top_hits()
        evidence_records: list[dict[str, Any]] = []
        priority_genes: list[str] = []
        if not evidence.empty:
            for rank, row in enumerate(evidence.to_dict("records"), start=1):
                sim = row.get("best_tanimoto_morgan_r2_2048")
                try:
                    similarity = float(sim) if sim is not None and not pd.isna(sim) else None
                except (TypeError, ValueError):
                    similarity = None
                gene = str(row.get("target_gene") or "").upper()
                if similarity is not None and similarity >= 0.35:
                    priority_genes.append(gene)
                include_neighbors = rank <= 3 or (similarity is not None and similarity >= 0.35) or (focus_target and gene == focus_target.upper())
                neighbors: list[dict[str, Any]] = []
                if include_neighbors and not top_hits.empty:
                    subset = top_hits[
                        top_hits["own_compound_id"].astype(str).str.upper().eq(compound_id)
                        & top_hits["target_gene"].astype(str).str.upper().eq(gene)
                    ].sort_values("rank_in_target").head(5)
                    neighbors = self._records(subset)
                cleaned = self._clean(row)
                cleaned["top_neighbors"] = neighbors
                evidence_records.append(cleaned)

        if focus_target:
            focus = focus_target.strip().upper()
            context_genes = [focus]
        elif priority_genes:
            context_genes = list(dict.fromkeys(priority_genes))
        elif evidence_records:
            context_genes = [str(evidence_records[0].get("target_gene") or "").upper()]
        else:
            context_genes = []

        contexts: list[dict[str, Any]] = []
        for gene in context_genes:
            gene_rows = [r for r in evidence_records if str(r.get("target_gene") or "").upper() == gene]
            similarity = gene_rows[0].get("best_tanimoto_morgan_r2_2048") if gene_rows else None
            target_contexts = self._contexts_for_target(gene)
            for context in target_contexts:
                context["pyz_target_similarity"] = similarity
                context["molecule_relevance_status"] = (
                    "priority_target_based_test_context" if similarity is not None and float(similarity) >= 0.35
                    else "weak_structural_link_context_only"
                )
                context["interpretation_ru"] = (
                    "Эта культура релевантна для проверки предполагаемой мишени в данном опухолевом контексте. "
                    "Это не означает, что PYZ уже показала активность на этой линии."
                )
                contexts.append(context)

        return self._clean({
            "available": True,
            "compound": compound,
            "admet_summary": self._admet_summary(),
            "target_evidence": evidence_records,
            "priority_target_genes": list(dict.fromkeys(priority_genes)),
            "relevant_contexts": contexts,
            "context_guardrail_ru": (
                "Органы, опухоли и клеточные линии здесь наследуются от биологической значимости предполагаемой мишени по Candidate v2/DepMap. "
                "До собственного cell assay нельзя говорить, что PYZ эффективна в этих культурах."
            ),
            "direct_pyz_target_evidence": "absent_not_measured",
            "docking_status": "pending_external_results",
        })

    def smiles(self, compound_id: str) -> str:
        registry = self.registry()
        hit = registry[registry["own_compound_id"].astype(str).str.upper().eq(compound_id.strip().upper())]
        if hit.empty:
            raise MCLDataError(f"Unknown PYZ compound: {compound_id}")
        return str(hit.iloc[0]["standardized_smiles"])
