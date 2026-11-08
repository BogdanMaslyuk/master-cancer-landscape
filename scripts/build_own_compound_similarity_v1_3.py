from __future__ import annotations

"""Own-compound similarity v1.3: compare PYZ only with strict oncology references.

The chemical calculation is unchanged (Morgan radius 2, 2048 bits, Tanimoto).
The biological comparator universe is changed deliberately: only manually curated
level A/B small-molecule oncology references are allowed. Biologics, research tools
(level C), and unsuitable/off-target relations (level D) are excluded.
"""

import json
from pathlib import Path
from typing import Any

import pandas as pd

import build_own_compound_similarity as impl
from build_own_compound_similarity_v1_2 import _standardize_smiles_v1_2

ROOT = Path(__file__).resolve().parents[1]
PHARM = ROOT / "data" / "runtime" / "pharmacology"
RUNTIME = ROOT / "data" / "runtime" / "own_compounds"
QC_DIR = ROOT / "outputs" / "qc"
REFERENCES = PHARM / "oncology_reference_structures.parquet"
OWN_CONFIG = ROOT / "config" / "own_compounds.tsv"


def _text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _load_strict_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if not OWN_CONFIG.exists():
        raise SystemExit("Missing config/own_compounds.tsv")
    if not REFERENCES.exists():
        raise SystemExit("Missing oncology_reference_structures.parquet; run resolve_oncology_reference_structures.py first.")

    own = pd.read_csv(OWN_CONFIG, sep="\t", dtype=str, keep_default_na=False)
    refs = pd.read_parquet(REFERENCES)
    refs = refs[
        refs["reference_level"].astype(str).isin({"A", "B"})
        & refs["modality"].astype(str).eq("small_molecule")
        & refs["similarity_eligible"].astype(bool)
        & refs["resolution_status"].astype(str).eq("resolved")
        & refs["canonical_smiles"].fillna("").astype(str).str.strip().ne("")
    ].copy()
    if refs.empty:
        raise SystemExit("No resolved A/B small-molecule oncology reference structures are available.")

    refs["compound_id"] = "ONCOLOGY_REFERENCE:" + refs["reference_id"].astype(str)
    refs["preferred_name"] = refs["compound_name"].astype(str)
    refs["actions_json"] = refs["action"].map(lambda x: json.dumps([_text(x)], ensure_ascii=False))
    refs["evidence_types_json"] = refs.apply(
        lambda r: json.dumps(
            [
                f"oncology_reference_level_{_text(r['reference_level'])}",
                _text(r.get("source_kind")) or "manual_curation",
                _text(r.get("target_relation")) or "target_relation_unspecified",
            ],
            ensure_ascii=False,
        ),
        axis=1,
    )
    for column in ["models_n", "dependent_models_n", "strong_dependency_models_n"]:
        refs[column] = pd.NA

    comparator_pairs = refs[[
        "compound_id", "target_gene", "preferred_name", "canonical_smiles", "inchikey",
        "actions_json", "evidence_types_json", "models_n", "dependent_models_n",
        "strong_dependency_models_n", "reference_id", "reference_level", "target_relation",
        "oncology_rationale_ru", "source_kind", "source_url",
    ]].copy()

    compound_catalog = comparator_pairs[[
        "compound_id", "preferred_name", "canonical_smiles", "inchikey"
    ]].drop_duplicates("compound_id").copy()
    for column in ["pubchem_cid", "chembl_id", "broad_id", "gdsc_id"]:
        compound_catalog[column] = ""

    return own, compound_catalog, comparator_pairs


def main() -> None:
    # Reuse the battle-tested calculation and RDKit repair path, while replacing
    # only the comparator universe and output names.
    impl._load_inputs = _load_strict_inputs
    impl._standardize_smiles = _standardize_smiles_v1_2

    impl.HITS_OUTPUT = RUNTIME / "oncology_reference_similarity_hits.parquet"
    impl.TOP_OUTPUT = RUNTIME / "oncology_reference_top_similarity_hits.parquet"
    impl.TARGET_OUTPUT = RUNTIME / "oncology_reference_target_similarity_summary.parquet"
    impl.MANIFEST = RUNTIME / "oncology_reference_similarity_manifest.json"
    impl.QC_OUTPUT = QC_DIR / "oncology_reference_similarity_qc.tsv"

    print("MCL Own Compound Structural Similarity v1.3")
    print("Comparator universe: strict oncology reference compounds, levels A+B, small molecules only")
    impl.main()

    manifest_path = impl.MANIFEST
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest.update({
            "contract": "mcl-own-compound-oncology-reference-similarity-v1.3",
            "comparator_universe": "oncology_reference_compounds levels A+B; small_molecule; resolved structure",
            "excluded_reference_levels": ["C", "D"],
            "excluded_modalities": ["biologic", "antibody_drug_conjugate"],
            "scientific_guardrail": "Structural similarity does not establish common target, binding mode, potency, selectivity, or antitumor efficacy.",
        })
        manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    refs = pd.read_parquet(REFERENCES)
    used = refs[
        refs["reference_level"].astype(str).isin({"A", "B"})
        & refs["modality"].astype(str).eq("small_molecule")
        & refs["similarity_eligible"].astype(bool)
        & refs["resolution_status"].astype(str).eq("resolved")
    ]
    print(f"Strict oncology reference structures used: {len(used)}")
    print(f"Targets represented in strict chemical universe: {used['target_gene'].nunique()}")
    missing_targets = sorted(set(refs["target_gene"].astype(str)) - set(used["target_gene"].astype(str)))
    if missing_targets:
        print("Targets without an A/B small-molecule structural comparator: " + ", ".join(missing_targets))
    print(f"Strict target summary: {impl.TARGET_OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
