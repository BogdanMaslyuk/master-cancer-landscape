from __future__ import annotations

"""Compatibility entry point for Own Compound Structural Similarity v1.1.

The original v1 builder remains the scientific calculation implementation. This
entry point overlays the separately traced compound_structure_registry onto the
runtime pharmacology compound catalog before the Morgan/Tanimoto calculation.
It does not mutate the source pharmacology catalog.
"""

from pathlib import Path

import pandas as pd

import build_own_compound_similarity as impl


ROOT = Path(__file__).resolve().parents[1]
STRUCTURE_REGISTRY = ROOT / "data" / "runtime" / "pharmacology" / "compound_structure_registry.parquet"
ENRICHED_CATALOG = ROOT / "data" / "runtime" / "pharmacology" / "compound_catalog_structures_enriched.parquet"


def _text_series(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _as_text_object(series: pd.Series) -> pd.Series:
    """Return a text/object Series before string overlay assignments.

    Older parquet materializations can infer all-empty metadata columns as
    float64/NaN. Assigning SMILES/InChIKey/CID strings into those columns emits
    pandas FutureWarning and will become an error in a future pandas release.
    """
    return series.astype("object").where(series.notna(), "")


def _build_overlay() -> Path:
    if not impl.COMPOUNDS.exists():
        raise SystemExit(f"Missing pharmacology compound catalog: {impl.COMPOUNDS}")
    if not STRUCTURE_REGISTRY.exists():
        raise SystemExit(
            "Missing compound structure registry. Run:\n"
            ".\\.venv\\Scripts\\python.exe .\\scripts\\enrich_priority_compound_structures_pubchem.py"
        )

    compounds = pd.read_parquet(impl.COMPOUNDS).copy()
    registry = pd.read_parquet(STRUCTURE_REGISTRY).copy()
    if registry.empty:
        raise SystemExit("Compound structure registry is empty.")

    registry = registry[
        registry.get("resolution_status", pd.Series("", index=registry.index)).astype(str).eq("resolved")
    ].copy()
    if registry.empty:
        raise SystemExit(
            "Compound structure registry contains no resolved structures; inspect "
            "outputs\\qc\\compound_structure_registry_qc.tsv"
        )

    keep = [
        c for c in (
            "compound_id", "canonical_smiles", "isomeric_smiles", "inchikey",
            "pubchem_cid_resolved", "resolution_source", "resolution_method",
        ) if c in registry.columns
    ]
    registry = registry[keep].drop_duplicates("compound_id", keep="last")
    merged = compounds.merge(registry, on="compound_id", how="left", suffixes=("", "_resolved"))

    if "canonical_smiles" not in merged.columns:
        merged["canonical_smiles"] = pd.Series("", index=merged.index, dtype="object")
    else:
        merged["canonical_smiles"] = _as_text_object(merged["canonical_smiles"])
    resolved_smiles = _as_text_object(
        merged.get("canonical_smiles_resolved", pd.Series("", index=merged.index, dtype="object"))
    )
    missing_smiles = _text_series(merged["canonical_smiles"]).eq("")
    merged.loc[missing_smiles, "canonical_smiles"] = resolved_smiles.loc[missing_smiles]

    if "inchikey" not in merged.columns:
        merged["inchikey"] = pd.Series("", index=merged.index, dtype="object")
    else:
        merged["inchikey"] = _as_text_object(merged["inchikey"])
    resolved_inchikey = _as_text_object(
        merged.get("inchikey_resolved", pd.Series("", index=merged.index, dtype="object"))
    )
    missing_inchikey = _text_series(merged["inchikey"]).eq("")
    merged.loc[missing_inchikey, "inchikey"] = resolved_inchikey.loc[missing_inchikey]

    if "pubchem_cid" not in merged.columns:
        merged["pubchem_cid"] = pd.Series("", index=merged.index, dtype="object")
    else:
        merged["pubchem_cid"] = _as_text_object(merged["pubchem_cid"])
    resolved_cid = _as_text_object(
        merged.get("pubchem_cid_resolved", pd.Series("", index=merged.index, dtype="object"))
    )
    missing_cid = _text_series(merged["pubchem_cid"]).eq("")
    merged.loc[missing_cid, "pubchem_cid"] = resolved_cid.loc[missing_cid]

    merged["structure_overlay_source"] = _as_text_object(
        merged.get("resolution_source", pd.Series("", index=merged.index, dtype="object"))
    )
    merged["structure_overlay_method"] = _as_text_object(
        merged.get("resolution_method", pd.Series("", index=merged.index, dtype="object"))
    )

    merged = merged.drop(
        columns=[
            "canonical_smiles_resolved", "inchikey_resolved", "isomeric_smiles",
            "pubchem_cid_resolved", "resolution_source", "resolution_method",
        ],
        errors="ignore",
    )
    merged.to_parquet(ENRICHED_CATALOG, index=False, compression="zstd")

    resolved_ids = set(registry["compound_id"].astype(str))
    enriched_n = int(compounds["compound_id"].astype(str).isin(resolved_ids).sum())
    structure_n = int(_text_series(merged["canonical_smiles"]).ne("").sum())
    print("MCL Compound Structure Overlay v1.1")
    print(f"Resolved registry compounds overlapping catalog: {enriched_n}")
    print(f"Compound catalog rows with structure after overlay: {structure_n}/{len(merged)}")
    print(f"Wrote {ENRICHED_CATALOG.relative_to(ROOT)}")
    return ENRICHED_CATALOG


if __name__ == "__main__":
    impl.COMPOUNDS = _build_overlay()
    impl.main()
