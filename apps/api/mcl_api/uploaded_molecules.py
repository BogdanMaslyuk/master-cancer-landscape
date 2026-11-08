from __future__ import annotations

import hashlib
import json
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

    def _cell_lines(self, targets: list[dict[str, Any]]) -> list[dict[str, Any]]:
        selected = [t for t in targets if float(t.get("best_tanimoto") or 0) >= 0.35] or targets[:3]
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
                    role = str(model.get("role") or "unclassified")
                    lab = bool(model.get("available_in_laboratory"))
                    if similarity >= 0.35 and role == "positive" and lab:
                        level, level_ru = "high", "Высокий"
                    elif similarity >= 0.35 and role == "positive":
                        level, level_ru = "medium", "Средний"
                    else:
                        level, level_ru = "exploratory", "Исследовательский"
                    reasons = []
                    if model.get("dependency_probability") is not None:
                        reasons.append(f"CRISPR P(dep)={float(model['dependency_probability']):.2f}")
                    if model.get("gene_effect") is not None:
                        reasons.append(f"Chronos Gene Effect={float(model['gene_effect']):.2f}")
                    if role == "positive":
                        reasons.append("положительная модель Candidate v2")
                    if role.startswith("negative_"):
                        reasons.append("кандидат в отрицательный опухолевый контроль")
                    if model.get("known_reference_compounds_n"):
                        reasons.append("есть фармакологический контекст известных веществ")
                    if lab:
                        reasons.append("линия есть в лабораторном реестре")
                    output.append(self._clean({
                        "target_gene": gene,
                        "target_similarity": similarity,
                        "support_class": target.get("support_class"),
                        "organ": context.get("mcl_organ_ru") or context.get("mcl_system_ru"),
                        "cancer_id": context.get("mcl_cancer_id"),
                        "cancer_name": context.get("mcl_cancer_name"),
                        "model_id": model.get("model_id"),
                        "cell_line_name": model.get("cell_line_name") or model.get("model_id"),
                        "role": role,
                        "dependency_probability": model.get("dependency_probability"),
                        "gene_effect": model.get("gene_effect"),
                        "known_reference_compounds": model.get("known_reference_compounds") or [],
                        "available_in_laboratory": lab,
                        "laboratory_lines": model.get("laboratory_lines") or [],
                        "priority_level": level,
                        "priority_level_ru": level_ru,
                        "reasons": reasons,
                    }))
        order = {"high": 0, "medium": 1, "exploratory": 2}
        output.sort(key=lambda row: (
            order.get(str(row.get("priority_level")), 9),
            0 if row.get("available_in_laboratory") else 1,
            0 if row.get("role") == "positive" else 1,
            -float(row.get("dependency_probability") or -1),
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
        focus = next((t for t in targets if float(t.get("best_tanimoto") or 0) >= 0.35), targets[0] if targets else None)
        if focus is None:
            return {"available": False, "reason": "no_target_context"}
        gene = str(focus["target_gene"])
        subset = [c for c in cells if c.get("target_gene") == gene]
        positives = [c for c in subset if c.get("role") == "positive"]
        negatives = [c for c in subset if str(c.get("role") or "").startswith("negative_")]
        positives.sort(key=lambda c: (0 if c.get("available_in_laboratory") else 1, -float(c.get("dependency_probability") or -1)))
        negatives.sort(key=lambda c: (0 if c.get("available_in_laboratory") else 1, float(c.get("dependency_probability") or 9)))
        return {
            "available": True,
            "focus_target_gene": gene,
            "focus_target_tanimoto": focus.get("best_tanimoto"),
            "positive_models": positives[:2],
            "negative_tumor_control": negatives[0] if negatives else None,
            "non_tumor_control": self._non_tumor_control(),
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
                "note_ru": "Диапазон концентраций задаётся после проверки растворимости и пилота; это план первого эксперимента, а не доказательство эффективности.",
            },
        }

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
        if not safe.startswith("MCL-UPL-"):
            raise MCLDataError(f"Unknown uploaded molecule: {analysis_id}")
        return self.runtime_dir / f"{safe}.json"

    def analyze(self, smiles: str, name: str | None = None) -> dict[str, Any]:
        mol, canonical, inchikey = self._standardize(smiles)
        analysis_id = self._analysis_id(canonical)
        path = self._record_path(analysis_id)
        if path.exists():
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["reused_existing_analysis"] = True
            return payload
        targets = self._target_search(mol)
        cells = self._cell_lines(targets)
        payload = self._clean({
            "available": True,
            "analysis_id": analysis_id,
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
            },
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
