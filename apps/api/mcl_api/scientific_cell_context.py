from __future__ import annotations

"""Full measured DepMap gene × model context, independent of reference-drug roles.

This is an explicitly *genetic* layer. It never infers pharmacological
sensitivity of an uploaded molecule from CRISPR gene dependence.
"""

import os
from pathlib import Path
from typing import Any

import pandas as pd


def _clean(value: Any) -> Any:
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


class FullGeneModelContext:
    def __init__(self, project_root: Path):
        self.project_root = Path(project_root).resolve()
        # A read-only external source root is useful for Git worktrees where
        # large scientific Parquet inputs are kept out of version control.
        base = Path(os.environ.get("MCL_SCIENTIFIC_DATA_ROOT") or project_root).resolve()
        self.data_root = base
        self.context_path = base / "data" / "processed" / "depmap_model_target_context.parquet"
        self.atlas_path = base / "data" / "processed" / "depmap_crispr_model_atlas.parquet"
        self.laboratory_path = base / "data" / "processed" / "laboratory_panel.parquet"

    def load(self, target_gene: str) -> dict[str, Any]:
        """Return all source-measured model records (no Top-N truncation)."""
        gene = target_gene.strip().upper()
        if not self.context_path.is_file():
            return {
                "available": False,
                "status": "source_not_available",
                "target_gene": gene,
                "items": [],
                "total_measured_models": 0,
                "source": "processed_depmap_target_context",
            }
        required = {"model_id", "target_gene", "dependency_probability", "gene_effect"}
        # pandas/PyArrow predicate pushdown avoids loading the full matrix.
        frame = pd.read_parquet(
            self.context_path,
            filters=[("target_gene", "==", gene)],
        )
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"Invalid DepMap target context schema: {sorted(missing)}")
        frame = frame.copy()
        frame["model_id"] = frame["model_id"].astype(str)
        frame["target_gene"] = frame["target_gene"].astype(str).str.upper()
        if frame["model_id"].duplicated().any():
            raise ValueError(f"Duplicate ModelID for {gene} in full DepMap context")
        if self.atlas_path.is_file():
            atlas = pd.read_parquet(self.atlas_path)
            if not {"model_id", "cell_line_name"}.issubset(atlas.columns):
                raise ValueError("Missing model_id or cell_line_name in DepMap atlas")
            if atlas["model_id"].duplicated().any():
                raise ValueError("Duplicate ModelID in DepMap atlas")
            frame = frame.merge(
                atlas[["model_id", "cell_line_name"]],
                on="model_id", how="left", validate="one_to_one",
            )
        elif "cell_line_name" not in frame.columns:
            frame["cell_line_name"] = None
        laboratory: dict[str, list[str]] = {}
        if self.laboratory_path.is_file():
            lab = pd.read_parquet(self.laboratory_path)
            if {"model_id", "lab_name"}.issubset(lab.columns):
                for model_id, group in lab.loc[lab["model_id"].notna()].groupby("model_id"):
                    laboratory[str(model_id)] = sorted(
                        set(str(name) for name in group["lab_name"].dropna())
                    )
        # PRISM is an independent evidence layer. Read only known-drug
        # observations for the locally registered subset (not all DepMap models).
        from .scientific_reference_context import KnownDrugReferenceContext

        reference = KnownDrugReferenceContext(self.project_root).load(
            gene, list(laboratory)
        )
        out = []
        for row in frame.to_dict("records"):
            probability = pd.to_numeric(
                pd.Series([row.get("dependency_probability")]), errors="coerce"
            ).iloc[0]
            if pd.notna(probability) and not (0 <= float(probability) <= 1):
                raise ValueError(f"Invalid P(dep) for {gene}/{row['model_id']}")
            p = float(probability) if pd.notna(probability) else None
            gene_effect = pd.to_numeric(
                pd.Series([row.get("gene_effect")]), errors="coerce"
            ).iloc[0]
            ge = float(gene_effect) if pd.notna(gene_effect) else None
            if p is None:
                status = "not_measured"
            elif p > 0.5:
                status = "dependent"
            else:
                status = "nondependent_measured"
            model_id = row["model_id"]
            lab_names = laboratory.get(model_id, [])
            reference_observations = (
                reference["observations_by_model"].get(model_id, [])
                if reference["available"] and lab_names else []
            )
            reference_status = (
                "not_queried_non_laboratory_model" if not lab_names
                else "normalized_source_unavailable" if not reference["available"]
                else "normalized_prism_reference_observations" if reference_observations
                else "no_linked_prism_lfc_observations"
            )
            out.append({
                "target_gene": gene,
                "model_id": model_id,
                "cell_line_name": _clean(row.get("cell_line_name")) or model_id,
                "cancer_id": _clean(row.get("mcl_cancer_id")),
                "cancer_name": _clean(row.get("mcl_cancer_name")),
                "organ": _clean(row.get("mcl_organ_ru")),
                "dependency_probability": p,
                "gene_effect": ge,
                "genetic_dependency_status": status,
                "crispr_measured": p is not None,
                "depmap_release": _clean(row.get("depmap_release")),
                "laboratory_lines": lab_names,
                "registered_in_laboratory": bool(lab_names),
                "physical_stock_verified": False,
                "source_layer": "full_processed_depmap_target_context",
                "direct_uploaded_molecule_response": "absent_not_measured",
                "known_reference_compounds": sorted(set(
                    str(obs.get("preferred_name") or obs["compound_id"])
                    for obs in reference_observations
                )),
                "reference_drug_observations": reference_observations,
                "reference_drug_observations_n": len(reference_observations),
                "reference_drug_evidence_status": reference_status,
                "reference_drug_source_scope": "laboratory_registered_models_only",
            })
        out.sort(key=lambda item: (
            0 if item["registered_in_laboratory"] else 1,
            0 if item["genetic_dependency_status"] == "dependent" else
            1 if item["genetic_dependency_status"] == "nondependent_measured" else 2,
            -(item["dependency_probability"] if item["dependency_probability"] is not None else -1),
            str(item["cell_line_name"]),
        ))
        return {
            "available": True,
            "status": "full_genetic_context",
            "target_gene": gene,
            "source": "processed_depmap_target_context",
            "total_measured_models": sum(item["crispr_measured"] for item in out),
            "total_model_rows": len(out),
            "laboratory_registered_models_n": sum(item["registered_in_laboratory"] for item in out),
            "reference_drug_context": {
                "available": reference["available"],
                "status": reference["status"],
                "queried_models_n": reference["queried_models_n"],
                "observations_n": reference["observations_n"],
                "missing_source_files": reference["missing_source_files"],
                "note_ru": (
                    "Измерения PRISM приведены только для моделей, отмеченных в лабораторном "
                    "реестре, и относятся к известным веществам, а не к новой молекуле."
                ),
            },
            "items": out,
            "interpretation_ru": (
                "Показана генетическая зависимость от гена по CRISPR; "
                "чувствительность к новой молекуле не измерялась. "
                "Наличие в лабораторном реестре не подтверждает физическую готовность культуры."
            ),
        }
