from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

import build_cellular_evidence_layer as base


ROOT = Path(__file__).resolve().parents[1]
PRISM_REGISTRY = ROOT / "data" / "runtime" / "pharmacology" / "prism_structure_registry.parquet"
MANIFEST = ROOT / "data" / "runtime" / "cellular_evidence" / "cellular_evidence_manifest.json"


def _text(value: Any) -> str:
    return base._text(value)


def _build_structure_map_v11(
    ligand_catalog: pd.DataFrame,
    compounds: pd.DataFrame,
    qc: list[dict[str, Any]],
) -> pd.DataFrame:
    if not PRISM_REGISTRY.exists():
        raise SystemExit(
            "Missing data/runtime/pharmacology/prism_structure_registry.parquet. "
            "Run scripts/resolve_prism_structures_pubchem.py first."
        )

    registry = pd.read_parquet(PRISM_REGISTRY)
    if registry.empty or "resolution_status" not in registry.columns:
        raise SystemExit("PRISM structure registry is empty or malformed.")

    resolved = registry[registry["resolution_status"].astype(str).eq("resolved")].copy()
    if resolved.empty:
        raise SystemExit("PRISM structure registry contains no resolved structures.")

    prism_rows: list[dict[str, Any]] = []
    for row in resolved.to_dict("records"):
        compound_id = _text(row.get("compound_id"))
        canonical, inchikey, error = base._standardize_smiles(row.get("canonical_smiles"))
        if error:
            qc.append(
                {
                    "stage": "prism_registry_standardization",
                    "entity_id": compound_id,
                    "status": "excluded",
                    "detail": error,
                }
            )
            continue
        prism_rows.append(
            {
                "prism_compound_id": compound_id,
                "prism_preferred_name": _text(row.get("preferred_name")),
                "prism_broad_id": _text(row.get("broad_id")),
                "prism_standardized_smiles": canonical,
                "inchikey": inchikey,
                "prism_structure_resolution_source": _text(row.get("resolution_source")),
                "prism_structure_resolution_method": _text(row.get("resolution_method")),
                "prism_structure_resolution_query": _text(row.get("resolution_query")),
                "prism_pubchem_cid": _text(row.get("pubchem_cid_resolved")),
            }
        )

    prism = pd.DataFrame(prism_rows)
    if prism.empty:
        return pd.DataFrame()

    ligands = ligand_catalog.copy()
    ligands["inchikey"] = ligands["inchikey"].map(_text)
    ligands = ligands[ligands["inchikey"] != ""].copy()
    mapped = ligands.merge(prism, on="inchikey", how="inner")
    if mapped.empty:
        return mapped

    mapped["match_basis"] = "standardized_parent_inchikey_exact_via_prism_structure_registry"
    keep = [
        "target_gene",
        "ligand_id",
        "inchikey",
        "canonical_smiles",
        "sources_json",
        "source_ligand_ids_json",
        "endpoint_types_json",
        "measurements_n",
        "lowest_reported_value_nm_across_endpoints",
        "lowest_value_endpoint_type",
        "lowest_value_relation",
        "lowest_value_source",
        "prism_compound_id",
        "prism_preferred_name",
        "prism_broad_id",
        "prism_standardized_smiles",
        "prism_structure_resolution_source",
        "prism_structure_resolution_method",
        "prism_structure_resolution_query",
        "prism_pubchem_cid",
        "match_basis",
    ]
    return mapped[[c for c in keep if c in mapped.columns]].drop_duplicates().reset_index(drop=True)


def _patch_manifest() -> None:
    if not MANIFEST.exists():
        return
    payload = json.loads(MANIFEST.read_text(encoding="utf-8"))
    payload["contract"] = "mcl-cellular-evidence-v1.1"
    payload["prism_structure_registry"] = str(PRISM_REGISTRY.relative_to(ROOT))
    payload["structure_mapping"] = (
        "PRISM chemical identities are resolved from source metadata, an existing MCL structure registry, "
        "or conservative PubChem PUG REST lookups; both PRISM and Target Ligand Space structures are then "
        "standardized with RDKit Cleanup -> FragmentParent -> canonical isomeric SMILES -> InChIKey and linked "
        "only by exact standardized parent InChIKey. PubChem resolution establishes chemical identity only, "
        "not target mechanism."
    )
    MANIFEST.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def main() -> None:
    base._build_structure_map = _build_structure_map_v11
    base.main()
    _patch_manifest()


if __name__ == "__main__":
    main()
