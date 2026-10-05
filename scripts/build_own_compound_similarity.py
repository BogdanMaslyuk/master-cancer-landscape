from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem.MolStandardize import rdMolStandardize
except ImportError as exc:
    raise SystemExit(
        'RDKit is required for structural similarity. Install the chemistry extra:\n'
        '.\\.venv\\Scripts\\python.exe -m pip install -e ".[chemistry]"'
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "own_compounds.tsv"
PHARM = ROOT / "data" / "runtime" / "pharmacology"
RUNTIME = ROOT / "data" / "runtime" / "own_compounds"
QC_DIR = ROOT / "outputs" / "qc"

COMPOUNDS = PHARM / "compound_catalog.parquet"
TARGET_COMPOUNDS = PHARM / "target_compound_catalog.parquet"
CANDIDATES = PHARM / "candidate_hypotheses_v2.parquet"

OWN_OUTPUT = RUNTIME / "own_compounds.parquet"
HITS_OUTPUT = RUNTIME / "similarity_hits.parquet"
TOP_OUTPUT = RUNTIME / "top_similarity_hits.parquet"
TARGET_OUTPUT = RUNTIME / "target_similarity_summary.parquet"
MANIFEST = RUNTIME / "own_compound_similarity_manifest.json"
QC_OUTPUT = QC_DIR / "own_compound_similarity_qc.tsv"

FP_RADIUS = 2
FP_SIZE = 2048
PRIORITY_STATUS = "priority_for_in_vitro"
MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_SIZE)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _json_unique(values: pd.Series | list[Any]) -> str:
    result: list[str] = []
    for value in values:
        text = _text(value)
        if text and text not in result:
            result.append(text)
    return json.dumps(result, ensure_ascii=False)


def _standardize_smiles(smiles: Any) -> tuple[Any | None, str | None, str | None, str | None]:
    text = _text(smiles)
    if not text:
        return None, None, None, "missing_smiles"
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return None, None, None, "rdkit_parse_failed"
    try:
        cleaned = rdMolStandardize.Cleanup(mol)
        parent = rdMolStandardize.FragmentParent(cleaned)
        canonical = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
        inchikey = Chem.MolToInchiKey(parent)
    except Exception as exc:  # RDKit failures vary by structure; preserve the reason for QC.
        return None, None, None, f"standardization_failed:{type(exc).__name__}"
    return parent, canonical, inchikey, None


def _similarity_band(value: float) -> str:
    if value >= 0.70:
        return "high_2d_similarity"
    if value >= 0.50:
        return "moderate_2d_similarity"
    if value >= 0.35:
        return "low_notable_2d_similarity"
    return "weak_2d_similarity"


def _priority_context(priority: pd.DataFrame) -> pd.DataFrame:
    if priority.empty:
        return pd.DataFrame(
            columns=[
                "target_gene", "priority_hypotheses_n", "priority_cancers_n",
                "priority_cancers_json", "priority_compounds_json",
            ]
        )
    rows: list[dict[str, Any]] = []
    for target_gene, group in priority.groupby("target_gene", sort=True):
        rows.append(
            {
                "target_gene": str(target_gene),
                "priority_hypotheses_n": int(len(group)),
                "priority_cancers_n": int(group["mcl_cancer_id"].nunique()) if "mcl_cancer_id" in group.columns else 0,
                "priority_cancers_json": _json_unique(
                    group["mcl_cancer_name"] if "mcl_cancer_name" in group.columns else []
                ),
                "priority_compounds_json": _json_unique(
                    group["preferred_name"] if "preferred_name" in group.columns else group["compound_id"]
                ),
            }
        )
    return pd.DataFrame(rows)


def _load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    required = [CONFIG, COMPOUNDS, TARGET_COMPOUNDS, CANDIDATES]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    if missing:
        raise SystemExit("Missing own-compound similarity inputs:\n" + "\n".join(missing))

    own = pd.read_csv(CONFIG, sep="\t", dtype=str, keep_default_na=False)
    compounds = pd.read_parquet(COMPOUNDS)
    target_compounds = pd.read_parquet(TARGET_COMPOUNDS)
    candidates = pd.read_parquet(CANDIDATES)

    if "priority_status" not in candidates.columns:
        raise SystemExit("Candidate v2 input has no priority_status column.")
    if "target_gene" not in target_compounds.columns or "compound_id" not in target_compounds.columns:
        raise SystemExit("target_compound_catalog.parquet must contain target_gene and compound_id.")

    candidates["target_gene"] = candidates["target_gene"].astype(str).str.upper()
    target_compounds["target_gene"] = target_compounds["target_gene"].astype(str).str.upper()
    return own, compounds, target_compounds.merge(
        candidates[candidates["priority_status"].astype(str).eq(PRIORITY_STATUS)][
            ["target_gene"]
        ].drop_duplicates(),
        on="target_gene",
        how="inner",
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compare proprietary MCL compounds with known compounds linked to targets that already "
            "have Candidate v2 priority_for_in_vitro hypotheses."
        )
    )
    parser.add_argument("--top-per-target", type=int, default=5)
    parser.add_argument("--global-top", type=int, default=50)
    args = parser.parse_args()

    own, compound_catalog, comparator_pairs = _load_inputs()
    candidates = pd.read_parquet(CANDIDATES)
    candidates["target_gene"] = candidates["target_gene"].astype(str).str.upper()
    priority = candidates[candidates["priority_status"].astype(str).eq(PRIORITY_STATUS)].copy()
    priority_targets = set(priority["target_gene"].dropna().astype(str).str.upper())
    priority_pairs = set(
        zip(
            priority["target_gene"].astype(str).str.upper(),
            priority["compound_id"].astype(str),
        )
    )
    context = _priority_context(priority)

    # target_compound_catalog already contains structure columns in current pharmacology builds,
    # but coalesce from compound_catalog to remain compatible with older materializations.
    metadata = [
        c for c in (
            "compound_id", "preferred_name", "canonical_smiles", "inchikey",
            "pubchem_cid", "chembl_id", "broad_id", "gdsc_id",
        ) if c in compound_catalog.columns
    ]
    fallback = compound_catalog[metadata].drop_duplicates("compound_id").copy()
    comparator_pairs = comparator_pairs.merge(
        fallback,
        on="compound_id",
        how="left",
        suffixes=("", "_catalog"),
    )
    for column in ("preferred_name", "canonical_smiles", "inchikey", "pubchem_cid", "chembl_id", "broad_id", "gdsc_id"):
        catalog_column = f"{column}_catalog"
        if catalog_column in comparator_pairs.columns:
            if column not in comparator_pairs.columns:
                comparator_pairs[column] = comparator_pairs[catalog_column]
            else:
                current = comparator_pairs[column].fillna("").astype(str).str.strip()
                comparator_pairs.loc[current.eq(""), column] = comparator_pairs.loc[current.eq(""), catalog_column]
            comparator_pairs = comparator_pairs.drop(columns=catalog_column)

    comparator_pairs = comparator_pairs.drop_duplicates(["target_gene", "compound_id"]).copy()
    comparator_compounds = comparator_pairs[
        [c for c in ("compound_id", "preferred_name", "canonical_smiles", "inchikey") if c in comparator_pairs.columns]
    ].drop_duplicates("compound_id").copy()

    known_cache: dict[str, dict[str, Any]] = {}
    comparator_qc: list[dict[str, Any]] = []
    for _, row in comparator_compounds.iterrows():
        compound_id = _text(row.get("compound_id"))
        mol, canonical, inchikey, error = _standardize_smiles(row.get("canonical_smiles"))
        if error:
            comparator_qc.append(
                {"entity_type": "known_comparator", "entity_id": compound_id, "qc_status": error}
            )
            continue
        known_cache[compound_id] = {
            "mol": mol,
            "fingerprint": MORGAN.GetFingerprint(mol),
            "standardized_smiles": canonical,
            "standardized_inchikey": inchikey,
        }

    own_rows: list[dict[str, Any]] = []
    own_cache: dict[str, dict[str, Any]] = {}
    own_qc: list[dict[str, Any]] = []
    for _, row in own.iterrows():
        own_id = _text(row.get("own_compound_id"))
        mol, canonical, inchikey, error = _standardize_smiles(row.get("standardized_smiles"))
        record = dict(row)
        record.update(
            {
                "canonical_smiles_rdkit": canonical,
                "inchikey_rdkit": inchikey,
                "structure_valid": error is None,
                "structure_qc_status": error or "ok",
            }
        )
        own_rows.append(record)
        if error:
            own_qc.append({"entity_type": "own_compound", "entity_id": own_id, "qc_status": error})
            continue
        own_cache[own_id] = {
            "mol": mol,
            "fingerprint": MORGAN.GetFingerprint(mol),
            "standardized_smiles": canonical,
            "standardized_inchikey": inchikey,
        }

    own_out = pd.DataFrame(own_rows)
    if not own_cache:
        raise SystemExit("No valid proprietary structures were available after RDKit parsing/standardization.")
    if not known_cache:
        raise SystemExit(
            "No known comparator compounds with valid structures were found for Candidate v2 priority targets."
        )

    hit_rows: list[dict[str, Any]] = []
    valid_pairs = comparator_pairs[comparator_pairs["compound_id"].astype(str).isin(known_cache)].copy()
    for own_id, own_data in own_cache.items():
        for _, comparator in valid_pairs.iterrows():
            compound_id = _text(comparator.get("compound_id"))
            known = known_cache[compound_id]
            similarity = float(
                DataStructs.TanimotoSimilarity(own_data["fingerprint"], known["fingerprint"])
            )
            target_gene = _text(comparator.get("target_gene")).upper()
            hit_rows.append(
                {
                    "own_compound_id": own_id,
                    "source_molecule_id": _text(
                        own_out.loc[own_out["own_compound_id"].eq(own_id), "source_molecule_id"].iloc[0]
                    ) if "source_molecule_id" in own_out.columns else "",
                    "target_gene": target_gene,
                    "comparator_compound_id": compound_id,
                    "comparator_name": _text(comparator.get("preferred_name")) or compound_id,
                    "tanimoto_morgan_r2_2048": similarity,
                    "similarity_band": _similarity_band(similarity),
                    "exact_standardized_structure": bool(
                        own_data["standardized_inchikey"]
                        and own_data["standardized_inchikey"] == known["standardized_inchikey"]
                    ),
                    "comparator_is_priority_compound": (target_gene, compound_id) in priority_pairs,
                    "own_canonical_smiles": own_data["standardized_smiles"],
                    "comparator_canonical_smiles": known["standardized_smiles"],
                    "comparator_inchikey": known["standardized_inchikey"],
                    "actions_json": _text(comparator.get("actions_json")),
                    "evidence_types_json": _text(comparator.get("evidence_types_json")),
                    "known_models_n": comparator.get("models_n"),
                    "known_dependent_models_n": comparator.get("dependent_models_n"),
                    "known_strong_dependency_models_n": comparator.get("strong_dependency_models_n"),
                    "built_at": _now(),
                }
            )

    hits = pd.DataFrame(hit_rows)
    hits = hits.merge(context, on="target_gene", how="left")
    hits = hits.sort_values(
        ["own_compound_id", "tanimoto_morgan_r2_2048", "target_gene", "comparator_name"],
        ascending=[True, False, True, True],
    ).reset_index(drop=True)

    top_parts: list[pd.DataFrame] = []
    for own_id, group in hits.groupby("own_compound_id", sort=True):
        global_top = group.nlargest(max(1, args.global_top), "tanimoto_morgan_r2_2048").copy()
        global_top["top_scope"] = "global_priority_target_universe"
        global_top["rank_in_scope"] = range(1, len(global_top) + 1)
        top_parts.append(global_top)
        for target_gene, target_group in group.groupby("target_gene", sort=True):
            target_top = target_group.nlargest(max(1, args.top_per_target), "tanimoto_morgan_r2_2048").copy()
            target_top["top_scope"] = f"target:{target_gene}"
            target_top["rank_in_scope"] = range(1, len(target_top) + 1)
            top_parts.append(target_top)
    top_hits = pd.concat(top_parts, ignore_index=True) if top_parts else pd.DataFrame()

    summary_rows: list[dict[str, Any]] = []
    for (own_id, target_gene), group in hits.groupby(["own_compound_id", "target_gene"], sort=True):
        best = group.sort_values("tanimoto_morgan_r2_2048", ascending=False).iloc[0]
        summary_rows.append(
            {
                "own_compound_id": own_id,
                "target_gene": target_gene,
                "known_comparators_with_structure_n": int(group["comparator_compound_id"].nunique()),
                "best_comparator_compound_id": best["comparator_compound_id"],
                "best_comparator_name": best["comparator_name"],
                "best_tanimoto_morgan_r2_2048": float(best["tanimoto_morgan_r2_2048"]),
                "best_similarity_band": best["similarity_band"],
                "best_comparator_is_priority_compound": bool(best["comparator_is_priority_compound"]),
                "priority_hypotheses_n": int(best["priority_hypotheses_n"] or 0),
                "priority_cancers_n": int(best["priority_cancers_n"] or 0),
                "priority_cancers_json": best["priority_cancers_json"],
                "priority_compounds_json": best["priority_compounds_json"],
                "built_at": _now(),
            }
        )
    target_summary = pd.DataFrame(summary_rows).sort_values(
        ["own_compound_id", "best_tanimoto_morgan_r2_2048"], ascending=[True, False]
    ).reset_index(drop=True)

    RUNTIME.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    own_out.to_parquet(OWN_OUTPUT, index=False, compression="zstd")
    hits.to_parquet(HITS_OUTPUT, index=False, compression="zstd")
    top_hits.to_parquet(TOP_OUTPUT, index=False, compression="zstd")
    target_summary.to_parquet(TARGET_OUTPUT, index=False, compression="zstd")
    target_summary.to_csv(RUNTIME / "target_similarity_summary.tsv", sep="\t", index=False)

    qc = pd.DataFrame([*own_qc, *comparator_qc])
    if qc.empty:
        qc = pd.DataFrame(columns=["entity_type", "entity_id", "qc_status"])
    qc.to_csv(QC_OUTPUT, sep="\t", index=False)

    comparator_total = int(comparator_compounds["compound_id"].nunique())
    comparator_valid = int(len(known_cache))
    structure_coverage = comparator_valid / comparator_total if comparator_total else 0.0
    manifest = {
        "contract": "mcl-own-compound-structural-similarity-v1",
        "built_at": _now(),
        "own_compounds_n": int(len(own_out)),
        "own_compounds_valid_n": int(len(own_cache)),
        "priority_targets_n": int(len(priority_targets)),
        "priority_hypotheses_n": int(len(priority)),
        "known_target_linked_comparator_compounds_n": comparator_total,
        "known_comparator_compounds_with_valid_structure_n": comparator_valid,
        "known_comparator_structure_coverage_fraction": structure_coverage,
        "similarity_pairs_n": int(len(hits)),
        "fingerprint": {
            "algorithm": "RDKit Morgan bit fingerprint",
            "radius": FP_RADIUS,
            "fp_size": FP_SIZE,
            "chirality": False,
            "metric": "Tanimoto",
        },
        "structure_standardization": (
            "RDKit MolFromSmiles -> rdMolStandardize.Cleanup -> FragmentParent -> canonical isomeric SMILES. "
            "The same procedure is applied to proprietary and known comparator structures."
        ),
        "comparator_universe": (
            "All known compounds in target_compound_catalog linked to any target having at least one "
            "Candidate v2 priority_for_in_vitro hypothesis; comparators are not restricted to the 222 priority rows themselves."
        ),
        "scientific_guardrail": (
            "2D structural similarity is hypothesis-generating evidence only. A high Tanimoto value does not establish "
            "the same target, binding mode, pharmacological action, potency, selectivity, or therapeutic effect."
        ),
        "similarity_bands_descriptive_only": {
            ">=0.70": "high_2d_similarity",
            "0.50-0.70": "moderate_2d_similarity",
            "0.35-0.50": "low_notable_2d_similarity",
            "<0.35": "weak_2d_similarity",
        },
        "outputs": [
            str(path.relative_to(ROOT))
            for path in (OWN_OUTPUT, HITS_OUTPUT, TOP_OUTPUT, TARGET_OUTPUT, QC_OUTPUT)
        ],
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print("MCL Own Compound Structural Similarity v1")
    print(f"Own compounds: {len(own_out)}; valid structures: {len(own_cache)}")
    print(f"Candidate v2 priority hypotheses: {len(priority)}")
    print(f"Priority targets: {len(priority_targets)}")
    print(
        "Known target-linked comparator compounds: "
        f"{comparator_total}; valid structures: {comparator_valid} "
        f"({structure_coverage:.1%})"
    )
    print(f"Own × target-linked comparator pairs: {len(hits)}")
    print("Best target-level matches per proprietary molecule:")
    for own_id, group in target_summary.groupby("own_compound_id", sort=True):
        best = group.head(5)
        print(f"  {own_id}")
        for _, row in best.iterrows():
            print(
                f"    {row['target_gene']}: {row['best_comparator_name']} "
                f"Tanimoto={row['best_tanimoto_morgan_r2_2048']:.3f} "
                f"[{row['best_similarity_band']}]"
            )
    print(f"Wrote {RUNTIME.relative_to(ROOT)}")
    print(f"QC: {QC_OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
