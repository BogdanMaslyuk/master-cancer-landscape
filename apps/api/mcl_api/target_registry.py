from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd


class ProteinTargetRegistry:
    """Small read-only registry resolving gene-mapped targets to protein metadata.

    Pharmacology sources frequently annotate a target only by gene symbol. MCL keeps
    that source resolution explicit and enriches the display with UniProtKB metadata
    when a unique reviewed protein mapping is available. The registry never upgrades
    a gene-level source annotation into isoform-level evidence.
    """

    def __init__(self, root: Path):
        self.root = Path(root).resolve()
        self.path = self.root / "data" / "processed" / "target_registry" / "protein_targets.parquet"

    @lru_cache(maxsize=1)
    def frame(self) -> pd.DataFrame:
        if not self.path.exists():
            return pd.DataFrame()
        frame = pd.read_parquet(self.path)
        if "target_gene" in frame.columns:
            frame["target_gene"] = frame["target_gene"].astype(str).str.upper()
        return frame

    def enrich(self, frame: pd.DataFrame, gene_col: str = "target_gene") -> pd.DataFrame:
        if frame.empty or gene_col not in frame.columns:
            return frame
        registry = self.frame()
        if registry.empty or "target_gene" not in registry.columns:
            return frame
        out = frame.copy()
        out[gene_col] = out[gene_col].astype(str).str.upper()
        keep = [c for c in (
            "target_gene", "protein_preferred_name", "protein_name_raw",
            "uniprot_primary_accession", "uniprot_accessions_json", "uniprot_entry_name",
            "protein_families", "protein_length", "protein_mapping_status",
            "source_resolution", "mapping_source", "mapping_retrieved_at", "uniprot_release",
        ) if c in registry.columns]
        return out.merge(
            registry[keep].drop_duplicates("target_gene"),
            left_on=gene_col,
            right_on="target_gene",
            how="left",
            suffixes=("", "_registry"),
        )

    def lookup(self, gene_symbol: str) -> dict[str, Any]:
        gene = gene_symbol.strip().upper()
        frame = self.frame()
        if frame.empty or "target_gene" not in frame.columns:
            return {
                "target_gene": gene,
                "protein_mapping_status": "registry_not_built",
                "source_resolution": "gene_mapped",
            }
        hit = frame[frame["target_gene"] == gene]
        if hit.empty:
            return {
                "target_gene": gene,
                "protein_mapping_status": "not_in_registry",
                "source_resolution": "gene_mapped",
            }
        row = hit.iloc[0].to_dict()
        cleaned: dict[str, Any] = {}
        for key, value in row.items():
            try:
                if pd.isna(value):
                    cleaned[str(key)] = None
                    continue
            except (TypeError, ValueError):
                pass
            if hasattr(value, "item"):
                try:
                    value = value.item()
                except (TypeError, ValueError):
                    pass
            cleaned[str(key)] = value
        return cleaned
