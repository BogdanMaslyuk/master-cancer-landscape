from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

from .own_compounds import OwnCompoundStore
from .store import MCLDataError


class UploadedMoleculeStore:
    """Persistent local analysis for arbitrary uploaded small molecules."""

    def __init__(self, root: Path, own_compounds: OwnCompoundStore):
        self.root = Path(root).resolve()
        self.own_compounds = own_compounds
        self.runtime_dir = self.root / "data" / "runtime" / "uploaded_molecules"
        self.ligand_catalog_path = (
            self.root / "data" / "runtime" / "target_ligand_space" / "ligand_catalog.parquet"
        )
        self.ligand_manifest_path = (
            self.root / "data" / "runtime" / "target_ligand_space" / "target_ligand_space_manifest.json"
        )

    @staticmethod
    def _chemistry():
        try:
            from rdkit import Chem, DataStructs
            from rdkit.Chem import Crippen, Descriptors, Lipinski, QED
            from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors
            from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
            from rdkit.Chem.MolStandardize import rdMolStandardize
        except ImportError as exc:
            raise ValueError(
                "RDKit is required for molecule analysis. Install project chemistry extras."
            ) from exc
        return {
            "Chem": Chem,
            "DataStructs": DataStructs,
            "Crippen": Crippen,
            "Descriptors": Descriptors,
            "Lipinski": Lipinski,
            "QED": QED,
            "rdFingerprintGenerator": rdFingerprintGenerator,
            "rdMolDescriptors": rdMolDescriptors,
            "FilterCatalog": FilterCatalog,
            "FilterCatalogParams": FilterCatalogParams,
            "rdMolStandardize": rdMolStandardize,
        }

    @classmethod
    def _standardize(cls, smiles: str) -> tuple[Any, str, str]:
        chem = cls._chemistry()
        Chem = chem["Chem"]
        rdMolStandardize = chem["rdMolStandardize"]
        text = str(smiles or "").strip()
        if not text:
            raise ValueError("SMILES is required")
        mol = Chem.MolFromSmiles(text, sanitize=True)
        if mol is None:
            raise ValueError("RDKit could not parse the supplied SMILES")
        parent = rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(mol))
        canonical = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
        repaired = Chem.MolFromSmiles(canonical, sanitize=True)
        if repaired is None:
            raise ValueError("Structure could not be reparsed after standardization")
        repaired.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(repaired)
        canonical = Chem.MolToSmiles(repaired, canonical=True, isomericSmiles=True)
        inchikey = Chem.MolToInchiKey(repaired)
        if not inchikey:
            raise ValueError("InChIKey could not be generated")
        return repaired, canonical, inchikey

    @staticmethod
    def _clean(value: Any) -> Any:
        if isinstance(value, dict):
            return {str(k): UploadedMoleculeStore._clean(v) for k, v in value.items()}
        if isinstance(value, (list, tuple, set)):
            return [UploadedMoleculeStore._clean(v) for v in value]
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
    def _alerts(cls, mol: Any, name: str) -> list[str]:
        chem = cls._chemistry()
        Params = chem["FilterCatalogParams"]
        params = Params()
        params.AddCatalog(getattr(Params.FilterCatalogs, name))
        catalog = chem["FilterCatalog"](params)
        return sorted({match.GetDescription() for match in catalog.GetMatches(mol)})

    @classmethod
    def _descriptors(cls, mol: Any) -> dict[str, Any]:
        chem = cls._chemistry()
        Chem = chem["Chem"]
        Crippen = chem["Crippen"]
        Descriptors = chem["Descriptors"]
        Lipinski = chem["Lipinski"]
        QED = chem["QED"]
        rdMolDescriptors = chem["rdMolDescriptors"]
        mw = float(Descriptors.MolWt(mol))
        logp = float(Crippen.MolLogP(mol))
        hbd = int(Lipinski.NumHDonors(mol))
        hba = int(Lipinski.NumHAcceptors(mol))
        rot = int(Lipinski.NumRotatableBonds(mol))
        tpsa = float(rdMolDescriptors.CalcTPSA(mol))
        pains = cls._alerts(mol, "PAINS")
        brenk = cls._alerts(mol, "BRENK")
        try:
            from rdkit.Contrib.SA_Score import sascorer
            sa = float(sascorer.calculateScore(mol))
        except Exception:
            sa = None
        return {
            "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
            "molecular_weight": mw,
            "exact_molecular_weight": float(rdMolDescriptors.CalcExactMolWt(mol)),
            "clogp": logp,
            "tpsa": tpsa,
            "hbd": hbd,
            "hba": hba,
            "rotatable_bonds": rot,
            "rings": int(rdMolDescriptors.CalcNumRings(mol)),
            "formal_charge": int(Chem.GetFormalCharge(mol)),
            "qed": float(QED.qed(mol)),
            "lipinski_violations": int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10),
            "veber_pass": bool(tpsa <= 140 and rot <= 10),
            "synthetic_accessibility": sa,
            "pains_alerts_n": len(pains),
            "pains_alerts": pains,
            "brenk_alerts_n": len(brenk),
            "brenk_alerts": brenk,
        }

    @staticmethod
    def _support(similarity: float) -> tuple[str, str]:
        if similarity >= 0.70:
            return "high_structural_transferability", "Высокая структурная поддержка"
        if similarity >= 0.50:
            return "moderate_structural_transferability", "Умеренная структурная поддержка"
        if similarity >= 0.35:
            return "exploratory_low_structural_similarity", "Исследовательская структурная поддержка"
        return "mechanism_not_transferable_by_2d_similarity", "Слабая структурная связь"

    @lru_cache(maxsize=1)
    def _ligand_space(self) -> dict[str, Any]:
        if not self.ligand_catalog_path.exists():
            return {"available": False, "targets": {}}
        chem = self._chemistry()
        Chem = chem["Chem"]
        generator = chem["rdFingerprintGenerator"].GetMorganGenerator(radius=2, fpSize=2048)
        frame = pd.read_parquet(self.ligand_catalog_path)
        needed = {"target_gene", "ligand_id", "canonical_smiles"}
        if frame.empty or not needed.issubset(frame.columns):
            return {"available": False, "targets": {}}
        frame = frame.copy()
        frame["target_gene"] = frame["target_gene"].astype(str).str.upper()
        targets: dict[str, dict[str, Any]] = {}
        for gene, group in frame.groupby("target_gene", sort=True):
            rows, fps = [], []
            for row in group.to_dict("records"):
                mol = Chem.MolFromSmiles(str(row.get("canonical_smiles") or ""), sanitize=True)
                if mol is None:
                    continue
                rows.append(row)
                fps.append(generator.GetFingerprint(mol))
            if rows:
                targets[str(gene)] = {"rows": rows, "fps": fps}
        return {"available": bool(targets), "targets": targets, "generator": generator}

    def _target_search(self, mol: Any) -> list[dict[str, Any]]:
        space = self._ligand_space()
        if not space.get("available"):
            return []
        DataStructs = self._chemistry()["DataStructs"]
        query_fp = space["generator"].GetFingerprint(mol)
        output: list[dict[str, Any]] = []
        for gene, target in space["targets"].items():
            sims = list(DataStructs.BulkTanimotoSimilarity(query_fp, target["fps"]))
            if not sims:
                continue
            order = sorted(range(len(sims)), key=lambda i: sims[i], reverse=True)
            best_i = order[0]
            best = target["rows"][best_i]
            best_sim = float(sims[best_i])
            support_class, support_label = self._support(best_sim)
            neighbors = []
            for rank, idx in enumerate(order[:5], 1):
                row = target["rows"][idx]
                neighbors.append(self._clean({
                    "rank": rank,
                    "ligand_id": row.get("ligand_id"),
                    "tanimoto": float(sims[idx]),
                    "lowest_reported_value_nm": row.get("lowest_reported_value_nm_across_endpoints"),
                    "endpoint": row.get("lowest_value_endpoint_type"),
                    "source": row.get("lowest_value_source"),
                }))
            output.append(self._clean({
                "target_gene": gene,
                "ligands_compared_n": len(target["rows"]),
                "best_tanimoto": best_sim,
                "support_class": support_class,
                "support_label_ru": support_label,
                "best_ligand_id": best.get("ligand_id"),
                "best_lowest_reported_value_nm": best.get("lowest_reported_value_nm_across_endpoints"),
                "best_endpoint": best.get("lowest_value_endpoint_type"),
                "best_source": best.get("lowest_value_source"),
                "hits_ge_0_35_n": sum(v >= 0.35 for v in sims),
                "hits_ge_0_50_n": sum(v >= 0.50 for v in sims),
                "hits_ge_0_70_n": sum(v >= 0.70 for v in sims),
                "top_neighbors": neighbors,
                "direct_uploaded_molecule_target_evidence": "absent_not_measured",
            }))
        output.sort(key=lambda row: (-float(row.get("best_tanimoto") or -1), str(row.get("target_gene"))))
        for rank, row in enumerate(output, 1):
            row["rank"] = rank
        return output

    @staticmethod
    def _dependency_value(value: Any) -> float | None:
        """Only measured, finite CRISPR probabilities are eligible for roles."""
        try:
            number = float(value)
        except (TypeError, ValueError):
            return None
        if pd.isna(number) or not 0 <= number <= 1:
            return None
        return number

    def _full_genetic_cell_lines(
        self, supported_targets: list[dict[str, Any]]
    ) -> list[dict[str, Any]] | None:
        """Use complete independent CRISPR evidence when processed data exist.

        Returns None only when the full source is unavailable; empty list
        means the full source was queried but no relevant model was found.
        """
        full_targets: list[tuple[dict[str, Any], dict[str, Any]]] = []
        for target in supported_targets:
            gene = str(target["target_gene"])
            context = self.own_compounds.full_target_models(gene)
            if not context.get("available"):
                return None
            full_targets.append((target, context))

        out: list[dict[str, Any]] = []
        for target, context in full_targets:
            gene = str(target["target_gene"])
            similarity = float(target.get("best_tanimoto") or 0)
            laboratory_triage = context.get("laboratory_genetic_triage") or {}
            triage_by_model = {
                str(item["model_id"]): item
                for item in laboratory_triage.get("items", [])
            }
            for model in context.get("items", []):
                genetic_lab_triage = triage_by_model.get(str(model.get("model_id")), {})
                status = str(model["genetic_dependency_status"])
                if status == "dependent":
                    role = "genetically_dependent_candidate"
                    reason = "Измеренная генетическая зависимость по CRISPR; ответ на новую молекулу неизвестен."
                elif status == "nondependent_measured":
                    role = "genetically_nondependent_candidate"
                    reason = "Генетическая независимость по измеренной P(dep); это не отрицательный лекарственный контроль."
                else:
                    role = "insufficient_data"
                    reason = "Отсутствует измерение вероятности генетической зависимости."
                lab = bool(model.get("registered_in_laboratory"))
                out.append(self._clean({
                    "target_gene": gene,
                    "target_similarity": similarity,
                    "support_class": target.get("support_class"),
                    "organ": model.get("organ"),
                    "cancer_id": model.get("cancer_id"),
                    "cancer_name": model.get("cancer_name"),
                    "model_id": model.get("model_id"),
                    "cell_line_name": model.get("cell_line_name"),
                    "role": role,
                    "original_candidate_role": None,
                    "genetic_dependency_status": status,
                    "dependency_probability": model.get("dependency_probability"),
                    "gene_effect": model.get("gene_effect"),
                    "depmap_release": model.get("depmap_release"),
                    "reference_drug_evidence_status": model.get("reference_drug_evidence_status"),
                    "reference_drug_source_scope": model.get("reference_drug_source_scope"),
                    "reference_drug_observations_n": model.get("reference_drug_observations_n", 0),
                    "reference_drug_observations": model.get("reference_drug_observations") or [],
                    "direct_uploaded_molecule_response": "absent_not_measured",
                    "known_reference_compounds": model.get("known_reference_compounds") or [],
                    "available_in_laboratory": lab,
                    "laboratory_lines": model.get("laboratory_lines") or [],
                    "laboratory_genetic_role": genetic_lab_triage.get("genetic_role"),
                    "laboratory_genetic_triage_group": genetic_lab_triage.get("genetic_triage_group"),
                    "same_cancer_genetic_comparison_possible": bool(
                        genetic_lab_triage.get("same_cancer_comparison_possible", False)
                    ),
                    "physical_stock_verified": False,
                    "source_layer": "full_processed_depmap_target_context",
                    "priority_level": "exploratory" if status != "not_measured" else "insufficient_data",
                    "priority_level_ru": "Исследовательский" if status != "not_measured" else "Недостаточно данных",
                    "reasons": [reason, "Физическое наличие культуры не проверено."]
                    if lab else [reason],
                    "warnings": [
                        "Генетическая зависимость не доказывает чувствительность к PYZ.",
                        "Структурное сходство не доказывает связывание с мишенью.",
                        "Фармакологические ответы известных препаратов приведены отдельно; они не подтверждают чувствительность к PYZ.",
                    ],
                }))
        out.sort(key=lambda c: (
            0 if c["available_in_laboratory"] else 1,
            0 if c["genetic_dependency_status"] == "dependent" else 1,
            -(c.get("dependency_probability") if c.get("dependency_probability") is not None else -1),
            str(c.get("target_gene") or ""), str(c.get("cell_line_name") or ""),
        ))
        return out

    def _cell_lines(self, targets: list[dict[str, Any]]) -> list[dict[str, Any]]:
        # Weak structural neighbors remain in the target table, not in
        # evidence-backed cell-line recommendations.
        selected = [t for t in targets if float(t.get("best_tanimoto") or 0) >= 0.35]
        if selected and hasattr(self.own_compounds, "full_target_models"):
            full = self._full_genetic_cell_lines(selected)
            if full is not None:
                return full
        # Compatibility fallback only: legacy Candidate v2 is truncated and
        # its role refers to a known reference drug, NOT to uploaded PYZ.
        output: list[dict[str, Any]] = []
        seen: set[tuple[str, str]] = set()
        for target in selected[:5]:
            gene = str(target["target_gene"])
            similarity = float(target.get("best_tanimoto") or 0)
            for context in self.own_compounds._contexts_for_target(gene):
                for model in context.get("models", []):
                    model_id = str(model.get("model_id") or "")
                    if not model_id or (gene, model_id) in seen:
                        continue
                    seen.add((gene, model_id))
                    original_role = str(model.get("role") or "unclassified")
                    reference_scope = str(model.get("role_scope") or "")
                    dependency = self._dependency_value(model.get("dependency_probability"))
                    lab = bool(model.get("available_in_laboratory"))
                    known = model.get("known_reference_compounds") or []
                    reasons: list[str] = []
                    warnings: list[str] = [
                        "Структурное сходство не подтверждает связывание новой молекулы с мишенью."
                    ]
                    if reference_scope == "reference_compound_specific_candidate_v2_truncated":
                        role = "known_drug_context_only"
                        level, level_ru = "exploratory", "Исследовательский"
                        reasons.append("Роли Candidate v2 относятся только к конкретным известным препаратам.")
                        warnings.append("Неполная выборка Candidate v2: нельзя автоматически считать модель положительной или отрицательной для PYZ.")
                    elif dependency is None:
                        role = "insufficient_data"
                        level, level_ru = "insufficient_data", "Недостаточно данных"
                        warnings.append("Нет достоверного измерения CRISPR-зависимости для этой пары модель–мишень.")
                    elif original_role == "positive" and dependency > 0.5:
                        role = "positive"
                        level, level_ru = ("high", "Высокий") if lab else ("medium", "Средний")
                        reasons.append("Подтверждена CRISPR-зависимость; модель отмечена в Candidate v2")
                    elif original_role.startswith("negative_") and dependency <= 0.5:
                        role = original_role
                        level, level_ru = "exploratory", "Исследовательский"
                        reasons.append("Измеренная низкая CRISPR-зависимость; кандидат в опухолевый контроль")
                        warnings.append("Отсутствие CRISPR-зависимости не доказывает нечувствительность к веществу.")
                    else:
                        role = "discordant_or_unclassified"
                        level, level_ru = "exploratory", "Исследовательский"
                        warnings.append("Оценка CRISPR не согласуется с назначенной ролью Candidate v2.")
                    if dependency is not None:
                        reasons.append(f"CRISPR P(dep)={dependency:.2f}")
                    if model.get("gene_effect") is not None:
                        try:
                            reasons.append(f"Chronos Gene Effect={float(model['gene_effect']):.2f}")
                        except (TypeError, ValueError):
                            pass
                    if known:
                        reasons.append("Известны вещества в фармакологическом контексте мишени; эффект новой молекулы не измерен")
                    if lab:
                        reasons.append("Клеточная линия есть в лабораторном реестре")
                    if role == "positive" and not known:
                        warnings.append("Независимая фармакологическая поддержка в этих данных не представлена.")
                    output.append(self._clean({
                        "target_gene": gene,
                        "target_similarity": similarity,
                        "support_class": target.get("support_class"),
                        "organ": context.get("mcl_organ_ru") or context.get("mcl_system_ru"),
                        "cancer_id": context.get("mcl_cancer_id"),
                        "cancer_name": context.get("mcl_cancer_name"),
                        "model_id": model_id,
                        "cell_line_name": model.get("cell_line_name") or model_id,
                        "role": role,
                        "original_candidate_role": original_role,
                        "reference_drug_observations": model.get("reference_drug_observations") or [],
                        "genetic_dependency_status": (
                            "unknown" if dependency is None else
                            "dependent" if dependency > 0.5 else "nondependent_measured"
                        ),
                        "source_layer": "legacy_candidate_v2_truncated",
                        "direct_uploaded_molecule_response": "absent_not_measured",
                        "reference_drug_context_only": True,
                        "dependency_probability": dependency,
                        "gene_effect": model.get("gene_effect"),
                        "known_reference_compounds": known,
                        "available_in_laboratory": lab,
                        "laboratory_lines": model.get("laboratory_lines") or [],
                        "priority_level": level,
                        "priority_level_ru": level_ru,
                        "reasons": reasons,
                        "warnings": warnings,
                    }))
        order = {"high": 0, "medium": 1, "exploratory": 2, "insufficient_data": 3}
        output.sort(key=lambda row: (
            order.get(str(row.get("priority_level")), 9),
            0 if row.get("available_in_laboratory") else 1,
            0 if row.get("role") == "positive" else 1,
            -(row.get("dependency_probability") if row.get("dependency_probability") is not None else -1),
            str(row.get("cell_line_name") or ""),
        ))
        return output

    def _non_tumor_control(self) -> dict[str, Any] | None:
        panel = self.own_compounds.lab_panel()
        if panel.empty or "control_role" not in panel.columns:
            return None
        hit = panel[panel["control_role"].astype(str).eq("general_non_tumor_control")]
        if hit.empty:
            return None
        row = self._clean(hit.iloc[0].to_dict())
        return {
            "lab_name": row.get("lab_name") or "BJ5ta",
            "model_id": row.get("model_id"),
            "role": "general_non_tumor_control",
            "guardrail_ru": "Общий человеческий неопухолевый контроль; не является нормальной тканью конкретного органа.",
        }

    def _experiment_panel(self, targets: list[dict[str, Any]], cells: list[dict[str, Any]]) -> dict[str, Any]:
        supported = [t for t in targets if float(t.get("best_tanimoto") or 0) >= 0.35]
        if any(c.get("source_layer") == "full_processed_depmap_target_context" for c in cells):
            focus = next(
                (t for t in supported if any(
                    c["target_gene"] == t["target_gene"]
                    and c.get("genetic_dependency_status") == "dependent"
                    for c in cells
                )),
                supported[0] if supported else {},
            )
            gene = str(focus.get("target_gene") or "")
            subset = [c for c in cells if c.get("target_gene") == gene]
            dependent = [c for c in subset if c.get("genetic_dependency_status") == "dependent"]
            nondependent = [c for c in subset if c.get("genetic_dependency_status") == "nondependent_measured"]
            dependent.sort(key=lambda c: (
                0 if c.get("available_in_laboratory") else 1,
                -(c.get("dependency_probability") or 0),
            ))
            nondependent.sort(key=lambda c: (
                0 if c.get("available_in_laboratory") else 1,
                c.get("dependency_probability") if c.get("dependency_probability") is not None else 2,
            ))
            # These are *genetic comparison candidates*, never negative
            # pharmacological controls or evidence of PYZ activity.
            paired: dict[str, dict[str, Any]] = {}
            for cell in subset:
                if not cell.get("same_cancer_genetic_comparison_possible"):
                    continue
                cancer = str(cell.get("cancer_id") or "")
                if not cancer:
                    continue
                group = paired.setdefault(cancer, {
                    "cancer_id": cancer,
                    "cancer_name": cell.get("cancer_name"),
                    "genetically_dependent_models": [],
                    "genetically_nondependent_models": [],
                    "comparison_status": "provisional_same_cancer_genetic_only",
                })
                if cell.get("laboratory_genetic_role") == "genetically_dependent_hypothesis_model":
                    group["genetically_dependent_models"].append(cell)
                elif cell.get("laboratory_genetic_role") == "same_cancer_genetic_comparator_candidate":
                    group["genetically_nondependent_models"].append(cell)
            # Only retain genuine same-cancer laboratory registry pairs.
            paired_groups = [
                g for g in paired.values()
                if g["genetically_dependent_models"] and g["genetically_nondependent_models"]
            ]
            paired_groups.sort(key=lambda g: str(g["cancer_id"]))
            return {
                "available": bool(subset),
                "panel_type": "exploratory_genetic_context",
                "same_cancer_genetic_comparison_groups": paired_groups,
                "same_cancer_genetic_comparison_groups_n": len(paired_groups),
                "status": "not_ready",
                "complete": False,
                "evidence_complete": False,
                "laboratory_ready": False,
                "focus_target_gene": gene or None,
                "focus_target_tanimoto": focus.get("best_tanimoto"),
                "positive_models": [],
                "negative_tumor_control": None,
                "genetically_dependent_candidates": dependent[:10],
                "genetically_nondependent_candidates": nondependent[:10],
                "genetic_candidate_models_n": len(dependent),
                "genetic_comparator_models_n": len(nondependent),
                "non_tumor_control": self._non_tumor_control(),
                "missing_components": [
                    "experimentally_supported_uploaded_molecule_target",
                    "compound_specific_pharmacology",
                    "validated_assay_protocol",
                    "verified_laboratory_availability",
                ],
                "reason": "crispr_context_without_uploaded_molecule_pharmacology",
                "interpretation_ru": (
                    "Найдены модели с измеренной генетической зависимостью от предполагаемой мишени. "
                    "Результат PRISM известных препаратов не переносится на PYZ. "
                    "Это предварительный контекст для планирования исследования, а не готовая "
                    "мишень-ориентированная экспериментальная панель."
                ),
            }
        if not supported:
            return {
                "available": False,
                "complete": False,
                "status": "not_ready",
                "missing_components": [
                    "supported_target_hypothesis",
                    "positive_tumour_model",
                    "independent_confirmatory_model",
                    "measured_negative_tumour_control",
                ],
                "reason": "no_target_reaches_exploratory_similarity_threshold",
                "interpretation_ru": (
                    "Ни одна target-гипотеза не достигла Tanimoto 0,35. "
                    "Не сформирована мишень-ориентированная экспериментальная панель; "
                    "рассмотрите независимый поиск мишени или первичный фенотипический скрининг."
                ),
                "positive_models": [],
                "negative_tumor_control": None,
                "non_tumor_control": self._non_tumor_control(),
            }
        # A top structural neighbor without a dependent tumour model is not
        # an experimentally actionable hypothesis: consider subsequent targets.
        focus = next(
            (t for t in supported if any(
                c.get("target_gene") == t["target_gene"] and c.get("role") == "positive"
                for c in cells
            )),
            supported[0],
        )
        gene = str(focus["target_gene"])
        subset = [c for c in cells if c.get("target_gene") == gene]
        positives = [c for c in subset if c.get("role") == "positive"]
        positives.sort(key=lambda c: (
            0 if c.get("available_in_laboratory") else 1,
            -(c.get("dependency_probability") if c.get("dependency_probability") is not None else -1),
        ))
        selected = positives[:2]
        # Prioritize a comparator from the same cancer context; exclude
        # missing, discordant, and unmeasured CRISPR results.
        negatives = [
            c for c in subset
            if str(c.get("role") or "").startswith("negative_")
            and c.get("dependency_probability") is not None
            and c["dependency_probability"] <= 0.5
            and c.get("model_id") not in {p.get("model_id") for p in selected}
        ]
        cancer_ids = {p.get("cancer_id") for p in selected if p.get("cancer_id")}
        negatives.sort(key=lambda c: (
            0 if c.get("cancer_id") in cancer_ids else 1,
            0 if c.get("available_in_laboratory") else 1,
            c["dependency_probability"],
        ))
        negative = negatives[0] if negatives else None
        normal = self._non_tumor_control()
        missing = []
        if not selected:
            missing.append("positive_tumour_model")
        if len(selected) < 2:
            missing.append("independent_confirmatory_model")
        if negative is None:
            missing.append("measured_negative_tumour_control")
        if normal is None:
            missing.append("non_tumour_control")
        # Distinguish an in-silico model set from a panel that can be run
        # with the locally registered laboratory cultures.
        if selected and not any(p.get("available_in_laboratory") for p in selected):
            missing.append("laboratory_positive_tumour_model")
        if len(selected) >= 2 and sum(bool(p.get("available_in_laboratory")) for p in selected) < 2:
            missing.append("laboratory_confirmatory_model")
        if negative is not None and not negative.get("available_in_laboratory"):
            missing.append("laboratory_negative_tumour_control")
        if negative is not None and selected:
            same_cancer = any(
                p.get("cancer_id") and p.get("cancer_id") == negative.get("cancer_id")
                for p in selected
            )
            if not same_cancer:
                missing.append("same_cancer_negative_control")
        # The legacy Candidate v2 list is Top-N truncated and retains known
        # reference-drug roles. Never certify it as a complete PYZ panel.
        if cells and not any(
            c.get("source_layer") == "full_processed_depmap_target_context"
            for c in cells
        ):
            missing.append("full_depmapped_evidence_layer")
        complete = not missing
        evidence_complete = bool(
            len(selected) >= 2 and negative is not None and normal is not None
        )
        laboratory_ready = bool(
            evidence_complete
            and all(p.get("available_in_laboratory") for p in selected)
            and negative.get("available_in_laboratory")
        )
        return {
            "available": bool(selected),
            "complete": complete,
            "evidence_complete": evidence_complete,
            "laboratory_ready": laboratory_ready,
            "status": "complete" if complete else ("partial" if selected else "not_ready"),
            "missing_components": missing,
            "interpretation_ru": (
                "Предложена экспериментальная панель для проверки гипотезы; "
                "прямое действие новой молекулы на мишень не установлено."
                if complete else
                "Панель неполная. Отсутствующие модели или контроли нельзя "
                "считать проверенными и готовыми к эксперименту."
            ),
            "focus_target_gene": gene,
            "focus_target_tanimoto": focus.get("best_tanimoto"),
            "positive_models": selected,
            "negative_tumor_control": negative,
            "non_tumor_control": normal,
            "suggested_first_assay": {
                "assay": "MTT/viability pilot",
                "exposure_hours": 72,
                "concentration_points_n": 8,
                "technical_replicates": 3,
                "independent_experiments_min_n": 3,
                "controls": [
                    "vehicle",
                    "positive cytotoxic control",
                    "cell-free blank",
                    "compound + assay reagent without cells",
                ],
                "note_ru": (
                    "Диапазон концентраций задаётся после проверки растворимости и пилота; "
                    "план не является доказательством эффективности."
                ),
            },
        }

    def _source_versions(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "mcl_molecule_analysis_contract": "mcl-uploaded-molecule-analysis-v1",
            "target_ligand_space": "current_local_runtime",
            "candidate_prioritization": "v2",
            "fingerprint": "Morgan radius=2, 2048 bits",
        }
        if self.ligand_manifest_path.exists():
            try:
                manifest = json.loads(self.ligand_manifest_path.read_text(encoding="utf-8"))
                payload["target_ligand_space_manifest"] = {
                    key: manifest.get(key)
                    for key in ("contract", "built_at", "targets_n", "ligands_n", "measurements_n")
                    if manifest.get(key) is not None
                }
            except (OSError, json.JSONDecodeError):
                pass
        return payload

    def _known_pyz(self, canonical: str) -> str | None:
        registry = self.own_compounds.registry()
        if registry.empty:
            return None
        for row in registry[["own_compound_id", "standardized_smiles"]].dropna().to_dict("records"):
            try:
                _, known, _ = self._standardize(str(row["standardized_smiles"]))
            except ValueError:
                continue
            if known == canonical:
                return str(row["own_compound_id"])
        return None

    def _analysis_id(self, canonical: str) -> str:
        digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:12].upper()
        return f"MCL-UPL-{digest}"

    def _record_path(self, analysis_id: str) -> Path:
        safe = analysis_id.strip().upper()
        if not re.fullmatch(r"MCL-UPL-[A-F0-9]{12}", safe):
            raise MCLDataError(f"Unknown uploaded molecule: {analysis_id}")
        return self.runtime_dir / f"{safe}.json"

    def analyze(
        self,
        smiles: str,
        name: str | None = None,
        force_recompute: bool = False,
    ) -> dict[str, Any]:
        mol, canonical, inchikey = self._standardize(smiles)
        analysis_id = self._analysis_id(canonical)
        path = self._record_path(analysis_id)
        if path.exists() and not force_recompute:
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["reused_existing_analysis"] = True
            return payload
        targets = self._target_search(mol)
        cells = self._cell_lines(targets)
        payload = self._clean({
            "available": True,
            "analysis_id": analysis_id,
            "analysis_version": "1.0.0-mvp",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "name": (name or "").strip() or analysis_id,
            "input_smiles": smiles,
            "standardized_smiles": canonical,
            "inchikey": inchikey,
            "matched_pyz_id": self._known_pyz(canonical),
            "descriptors": self._descriptors(mol),
            "admet": {
                "status": "not_available_for_arbitrary_uploaded_molecule_in_mvp",
                "interpretation_ru": "Универсальная локальная ADMET-модель для произвольных структур ещё не подключена; данные PYZ не переносятся на новую молекулу.",
            },
            "targets": targets,
            "cell_line_recommendations": cells,
            "experiment_panel": self._experiment_panel(targets, cells),
            "normal_tissue_liability": {
                "status": "planned_next_layer",
                "interpretation_ru": "GTEx/Human Protein Atlas и расширенный normal-tissue liability слой ещё не подключены. BJ5ta используется только как общий неопухолевый контроль.",
            },
            "method": {
                "fingerprint": "Morgan radius=2, 2048 bits",
                "similarity": "Tanimoto",
                "target_space": "current local Target Ligand Space",
                "target_space_scope": "current curated MCL core targets; not a proteome-wide target search",
            },
            "source_versions": self._source_versions(),
            "guardrail_ru": "MCL выбирает наиболее информативные модели для проверки гипотезы; он не утверждает, что загруженная молекула уже действует на выбранную мишень или клеточную линию.",
            "reused_existing_analysis": False,
        })
        self.runtime_dir.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".json.part")
        tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        tmp.replace(path)
        return payload

    def detail(self, analysis_id: str) -> dict[str, Any]:
        path = self._record_path(analysis_id)
        if not path.exists():
            raise MCLDataError(f"Unknown uploaded molecule: {analysis_id}")
        return json.loads(path.read_text(encoding="utf-8"))

    def targets(self, analysis_id: str) -> dict[str, Any]:
        payload = self.detail(analysis_id)
        return {"available": True, "analysis_id": analysis_id, "targets": payload.get("targets", []), "guardrail_ru": payload.get("guardrail_ru")}

    def cell_lines(self, analysis_id: str, target_gene: str | None = None, laboratory_only: bool = False) -> dict[str, Any]:
        payload = self.detail(analysis_id)
        rows = list(payload.get("cell_line_recommendations", []))
        if target_gene:
            gene = target_gene.strip().upper()
            rows = [row for row in rows if str(row.get("target_gene") or "").upper() == gene]
        if laboratory_only:
            rows = [row for row in rows if row.get("available_in_laboratory")]
        return {"available": True, "analysis_id": analysis_id, "items": rows, "total": len(rows), "guardrail_ru": payload.get("guardrail_ru")}

    def experiment(self, analysis_id: str) -> dict[str, Any]:
        payload = self.detail(analysis_id)
        return {"available": True, "analysis_id": analysis_id, "experiment_panel": payload.get("experiment_panel", {}), "guardrail_ru": payload.get("guardrail_ru")}

    def smiles(self, analysis_id: str) -> str:
        return str(self.detail(analysis_id)["standardized_smiles"])
