from __future__ import annotations

"""Normalized PRISM observations for known compounds, kept separate from PYZ.

Each observation remains linked to a particular model, compound, assay and
source release. Target annotations are metadata, not proof that observed
cell response is mediated by that target.
"""

import os
from pathlib import Path
from typing import Any

import pandas as pd


class KnownDrugReferenceContext:
    """Read-only source reader; currently queries laboratory-registered models."""

    def __init__(self, root: Path):
        base = Path(os.environ.get("MCL_SCIENTIFIC_DATA_ROOT") or root).resolve()
        pharm = base / "data" / "runtime" / "pharmacology"
        self.links_path = pharm / "model_compound_target_links.parquet"
        self.responses_path = pharm / "responses.parquet"
        self.compounds_path = pharm / "compound_catalog.parquet"
        self.evidence_path = pharm / "target_evidence.parquet"

    @staticmethod
    def _one(value: Any) -> Any:
        if value is None:
            return None
        try:
            if pd.isna(value):
                return None
        except (TypeError, ValueError):
            pass
        if hasattr(value, "item"):
            return value.item()
        return value

    def load(self, gene: str, model_ids: list[str]) -> dict[str, Any]:
        gene = gene.strip().upper()
        ids = sorted({str(mid) for mid in model_ids if mid})
        missing = [
            str(p) for p in (
                self.links_path, self.responses_path,
                self.compounds_path, self.evidence_path
            ) if not p.is_file()
        ]
        if missing:
            return {
                "available": False, "status": "normalized_source_unavailable",
                "target_gene": gene, "queried_models_n": len(ids),
                "observations_n": None, "observations_by_model": {},
                "missing_source_files": missing,
            }
        if not ids:
            return {
                "available": True, "status": "no_models_to_query",
                "target_gene": gene, "queried_models_n": 0,
                "observations_n": 0, "observations_by_model": {},
                "missing_source_files": [],
            }

        links = pd.read_parquet(
            self.links_path, filters=[("target_gene", "==", gene),
                                      ("model_id", "in", ids)],
        )
        expected_links = {"target_gene", "model_id", "compound_id"}
        if not expected_links.issubset(links.columns):
            raise ValueError("Source target-links missing: " +
                             ", ".join(sorted(expected_links - set(links.columns))))
        links = links.loc[links["target_gene"].astype(str).str.upper().eq(gene)].copy()
        links["model_id"] = links["model_id"].astype(str)
        links["compound_id"] = links["compound_id"].astype(str)
        keys = links[["model_id", "compound_id"]].drop_duplicates()

        if keys.empty:
            return {
                "available": True, "status": "no_linked_known_drug_pairs",
                "target_gene": gene, "queried_models_n": len(ids),
                "observations_n": 0, "observations_by_model": {},
                "missing_source_files": [],
            }

        responses = pd.read_parquet(
            self.responses_path, filters=[("model_id", "in", ids)],
        )
        required = {
            "observation_id", "model_id", "compound_id", "source",
            "source_release", "source_assay_id", "endpoint", "value"
        }
        if not required.issubset(responses.columns):
            raise ValueError("Normalized PRISM responses missing columns: " +
                             ", ".join(sorted(required - set(responses.columns))))
        responses = responses.loc[
            responses["source"].astype(str).str.upper().eq("PRISM")
            & responses["endpoint"].astype(str).str.upper().eq("LFC")
        ].copy()
        responses["model_id"] = responses["model_id"].astype(str)
        responses["compound_id"] = responses["compound_id"].astype(str)
        observations = responses.merge(
            keys, on=["model_id", "compound_id"],
            how="inner", validate="many_to_one",
        )
        if observations.empty:
            return {
                "available": True, "status": "no_matching_prism_lfc",
                "target_gene": gene, "queried_models_n": len(ids),
                "observations_n": 0, "observations_by_model": {},
                "missing_source_files": [],
            }

        if observations["observation_id"].isna().any() or observations["observation_id"].duplicated().any():
            raise ValueError("Missing or duplicate PRISM observation_id in model-compound link")

        compounds = pd.read_parquet(
            self.compounds_path,
            filters=[("compound_id", "in", sorted(set(observations["compound_id"])))],
        )
        if {"compound_id", "preferred_name"}.issubset(compounds.columns):
            if compounds["compound_id"].duplicated().any():
                raise ValueError("Duplicate compound_id in known compound catalog")
            compounds = compounds[["compound_id", "preferred_name"]].copy()
            compounds["compound_id"] = compounds["compound_id"].astype(str)
            observations = observations.merge(
                compounds, on="compound_id", how="left", validate="many_to_one",
            )

        annotations = pd.read_parquet(
            self.evidence_path,
            filters=[("target_gene", "==", gene),
                     ("compound_id", "in", sorted(set(observations["compound_id"])))],
        )
        if not {"compound_id", "target_gene"}.issubset(annotations.columns):
            raise ValueError("Missing compound/target identifier in target evidence")
        annotation_by_compound = {}
        for compound_id, group in annotations.groupby("compound_id"):
            evidence_items = []
            for row in group.to_dict("records"):
                evidence_items.append({
                    "evidence_id": self._one(row.get("evidence_id")),
                    "target_gene": self._one(row.get("target_gene")),
                    "evidence_type": self._one(row.get("evidence_type")),
                    "directness": self._one(row.get("directness")),
                    "action": self._one(row.get("action")),
                    "source": self._one(row.get("source")),
                })
            annotation_by_compound[str(compound_id)] = evidence_items

        by_model: dict[str, list[dict[str, Any]]] = {}
        observations = observations.sort_values(
            ["model_id", "compound_id", "observation_id"],
        )
        for row in observations.to_dict("records"):
            compound_id = str(row["compound_id"])
            item = {
                "model_id": str(row["model_id"]),
                "compound_id": compound_id,
                "preferred_name": self._one(row.get("preferred_name")),
                "observation_id": str(row["observation_id"]),
                "source_assay_id": self._one(row.get("source_assay_id")),
                "source_release": self._one(row.get("source_release")),
                "source": "PRISM",
                "endpoint": "LFC",
                "prism_lfc": self._one(row.get("value")),
                "dose": self._one(row.get("dose")),
                "dose_unit": self._one(row.get("dose_unit")),
                "exposure_time_h": self._one(row.get("exposure_time_h")),
                "reference_model_role": None,
                "candidate_priority_status": None,
                "target_annotations": annotation_by_compound.get(compound_id, []),
                "response_is_for_uploaded_molecule": False,
                "interpretation_ru": (
                    "Ответ на конкретное известное вещество PRISM; "
                    "аннотация мишени не доказывает, что эффект опосредован этим белком."
                ),
            }
            by_model.setdefault(str(row["model_id"]), []).append(item)
        return {
            "available": True, "status": "normalized_prism_reference_observations",
            "target_gene": gene, "queried_models_n": len(ids),
            "observations_n": len(observations),
            "observations_by_model": by_model, "missing_source_files": [],
            "note": (
                "Only normalized PRISM records, not original vendor raw data. "
                "Coverage is limited to queried laboratory-registered models."
            ),
        }
