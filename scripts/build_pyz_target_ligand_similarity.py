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
        'RDKit is required. Install chemistry extras: .\\.venv\\Scripts\\python.exe -m pip install -e ".[chemistry]"'
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
OWN_CONFIG = ROOT / "config" / "own_compounds.tsv"
TLS = ROOT / "data" / "runtime" / "target_ligand_space"
LIGAND_CATALOG = TLS / "ligand_catalog.parquet"
TARGET_MAP = TLS / "target_map.tsv"
RUNTIME = ROOT / "data" / "runtime" / "own_compounds"
QC_DIR = ROOT / "outputs" / "qc"

SUMMARY_OUT = RUNTIME / "target_ligand_space_similarity_summary.parquet"
SUMMARY_TSV_OUT = RUNTIME / "target_ligand_space_similarity_summary.tsv"
TOP_HITS_OUT = RUNTIME / "target_ligand_space_similarity_top_hits.parquet"
TOP_HITS_TSV_OUT = RUNTIME / "target_ligand_space_similarity_top_hits.tsv"
MANIFEST_OUT = RUNTIME / "target_ligand_space_similarity_manifest.json"
QC_OUT = QC_DIR / "target_ligand_space_similarity_qc.tsv"

FP_RADIUS = 2
FP_SIZE = 2048
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


def _standardize(smiles: Any) -> tuple[Any | None, str | None, str | None, str | None]:
    text = _text(smiles)
    if not text:
        return None, None, None, "missing_smiles"
    try:
        mol = Chem.MolFromSmiles(text, sanitize=True)
        if mol is None:
            return None, None, None, "rdkit_parse_failed"
        cleaned = rdMolStandardize.Cleanup(mol)
        parent = rdMolStandardize.FragmentParent(cleaned)
        canonical = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
        repaired = Chem.MolFromSmiles(canonical, sanitize=True)
        if repaired is None:
            return None, None, None, "post_standardization_reparse_failed"
        repaired.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(repaired)
        Chem.GetSymmSSSR(repaired)
        canonical = Chem.MolToSmiles(repaired, canonical=True, isomericSmiles=True)
        inchikey = Chem.MolToInchiKey(repaired)
        fp = MORGAN.GetFingerprint(repaired)
        return fp, canonical, inchikey, None
    except Exception as exc:
        return None, None, None, f"fingerprint_preparation_failed:{type(exc).__name__}"


def _band(value: float | None) -> str:
    if value is None:
        return "no_comparator"
    if value >= 0.70:
        return "high_2d_similarity"
    if value >= 0.50:
        return "moderate_2d_similarity"
    if value >= 0.35:
        return "low_notable_2d_similarity"
    return "weak_2d_similarity"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare PYZ-001...PYZ-065 with the experimental Target Ligand Space v1."
    )
    parser.add_argument("--top-per-target", type=int, default=20)
    args = parser.parse_args()
    if args.top_per_target < 1:
        raise SystemExit("--top-per-target must be >= 1")

    missing = [p for p in (OWN_CONFIG, LIGAND_CATALOG, TARGET_MAP) if not p.exists()]
    if missing:
        raise SystemExit("Missing similarity inputs:\n" + "\n".join(str(p.relative_to(ROOT)) for p in missing))

    own = pd.read_csv(OWN_CONFIG, sep="\t", dtype=str, keep_default_na=False)
    catalog = pd.read_parquet(LIGAND_CATALOG)
    target_map = pd.read_csv(TARGET_MAP, sep="\t", dtype=str, keep_default_na=False)
    target_map["target_gene"] = target_map["target_gene"].astype(str).str.upper()
    catalog["target_gene"] = catalog["target_gene"].astype(str).str.upper()

    RUNTIME.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    qc: list[dict[str, Any]] = []

    own_cache: dict[str, dict[str, Any]] = {}
    for row in own.to_dict("records"):
        own_id = _text(row.get("own_compound_id"))
        fp, canonical, inchikey, error = _standardize(row.get("standardized_smiles"))
        if error:
            qc.append({"entity_type": "own_compound", "entity_id": own_id, "target_gene": "", "status": error})
            continue
        own_cache[own_id] = {"fp": fp, "canonical_smiles": canonical, "inchikey": inchikey}
    if not own_cache:
        raise SystemExit("No valid PYZ structures were available.")

    ligand_rows: list[dict[str, Any]] = []
    ligand_fps: dict[str, Any] = {}
    for row in catalog.to_dict("records"):
        ligand_id = _text(row.get("ligand_id"))
        fp, canonical, inchikey, error = _standardize(row.get("canonical_smiles"))
        if error:
            qc.append(
                {
                    "entity_type": "target_ligand",
                    "entity_id": ligand_id,
                    "target_gene": _text(row.get("target_gene")),
                    "status": error,
                }
            )
            continue
        record = dict(row)
        record["canonical_smiles_similarity"] = canonical
        record["inchikey_similarity"] = inchikey
        ligand_rows.append(record)
        ligand_fps[ligand_id] = fp

    valid_catalog = pd.DataFrame(ligand_rows)
    if valid_catalog.empty:
        raise SystemExit("No valid Target Ligand Space structures were available.")

    target_groups: dict[str, pd.DataFrame] = {
        gene: group.reset_index(drop=True)
        for gene, group in valid_catalog.groupby("target_gene", sort=True)
    }

    summary_rows: list[dict[str, Any]] = []
    top_rows: list[dict[str, Any]] = []
    all_targets = target_map["target_gene"].dropna().astype(str).str.upper().drop_duplicates().tolist()

    for own_id, own_data in sorted(own_cache.items()):
        for gene in all_targets:
            group = target_groups.get(gene)
            if group is None or group.empty:
                summary_rows.append(
                    {
                        "own_compound_id": own_id,
                        "target_gene": gene,
                        "ligands_compared_n": 0,
                        "hits_ge_0_35_n": 0,
                        "hits_ge_0_50_n": 0,
                        "hits_ge_0_70_n": 0,
                        "best_ligand_id": "",
                        "best_tanimoto_morgan_r2_2048": None,
                        "best_similarity_band": "no_comparator",
                        "best_ligand_sources_json": "[]",
                        "best_endpoint_types_json": "[]",
                        "best_lowest_reported_value_nm_across_endpoints": None,
                        "best_lowest_value_endpoint_type": "",
                        "best_lowest_value_source": "",
                        "built_at": _now(),
                    }
                )
                continue

            ids = group["ligand_id"].astype(str).tolist()
            fps = [ligand_fps[x] for x in ids]
            similarities = list(DataStructs.BulkTanimotoSimilarity(own_data["fp"], fps))
            order = sorted(range(len(similarities)), key=lambda i: similarities[i], reverse=True)
            best_idx = order[0]
            best = group.iloc[best_idx]
            best_similarity = float(similarities[best_idx])

            summary_rows.append(
                {
                    "own_compound_id": own_id,
                    "target_gene": gene,
                    "ligands_compared_n": int(len(group)),
                    "hits_ge_0_35_n": int(sum(v >= 0.35 for v in similarities)),
                    "hits_ge_0_50_n": int(sum(v >= 0.50 for v in similarities)),
                    "hits_ge_0_70_n": int(sum(v >= 0.70 for v in similarities)),
                    "best_ligand_id": _text(best.get("ligand_id")),
                    "best_tanimoto_morgan_r2_2048": best_similarity,
                    "best_similarity_band": _band(best_similarity),
                    "best_ligand_sources_json": _text(best.get("sources_json")),
                    "best_endpoint_types_json": _text(best.get("endpoint_types_json")),
                    "best_lowest_reported_value_nm_across_endpoints": best.get("lowest_reported_value_nm_across_endpoints"),
                    "best_lowest_value_endpoint_type": _text(best.get("lowest_value_endpoint_type")),
                    "best_lowest_value_source": _text(best.get("lowest_value_source")),
                    "built_at": _now(),
                }
            )

            for rank, idx in enumerate(order[: args.top_per_target], start=1):
                ligand = group.iloc[idx]
                similarity = float(similarities[idx])
                top_rows.append(
                    {
                        "own_compound_id": own_id,
                        "target_gene": gene,
                        "rank_in_target": rank,
                        "ligand_id": _text(ligand.get("ligand_id")),
                        "tanimoto_morgan_r2_2048": similarity,
                        "similarity_band": _band(similarity),
                        "ligand_canonical_smiles": _text(ligand.get("canonical_smiles")),
                        "ligand_inchikey": _text(ligand.get("inchikey")),
                        "sources_json": _text(ligand.get("sources_json")),
                        "source_ligand_ids_json": _text(ligand.get("source_ligand_ids_json")),
                        "endpoint_types_json": _text(ligand.get("endpoint_types_json")),
                        "measurements_n": int(ligand.get("measurements_n") or 0),
                        "lowest_reported_value_nm_across_endpoints": ligand.get("lowest_reported_value_nm_across_endpoints"),
                        "lowest_value_endpoint_type": _text(ligand.get("lowest_value_endpoint_type")),
                        "lowest_value_relation": _text(ligand.get("lowest_value_relation")),
                        "lowest_value_source": _text(ligand.get("lowest_value_source")),
                        "built_at": _now(),
                    }
                )

    summary = pd.DataFrame(summary_rows).sort_values(
        ["own_compound_id", "best_tanimoto_morgan_r2_2048", "target_gene"],
        ascending=[True, False, True],
        na_position="last",
    ).reset_index(drop=True)
    top_hits = pd.DataFrame(top_rows).sort_values(
        ["own_compound_id", "target_gene", "rank_in_target"]
    ).reset_index(drop=True)

    summary.to_parquet(SUMMARY_OUT, index=False, compression="zstd")
    summary.to_csv(SUMMARY_TSV_OUT, sep="\t", index=False)
    top_hits.to_parquet(TOP_HITS_OUT, index=False, compression="zstd")
    top_hits.to_csv(TOP_HITS_TSV_OUT, sep="\t", index=False)
    pd.DataFrame(qc, columns=["entity_type", "entity_id", "target_gene", "status"]).to_csv(
        QC_OUT, sep="\t", index=False
    )

    manifest = {
        "contract": "mcl-pyz-target-ligand-similarity-v1",
        "built_at": _now(),
        "own_compounds_n": int(len(own_cache)),
        "targets_n": int(len(all_targets)),
        "valid_target_ligands_n": int(len(valid_catalog)),
        "fingerprint": {"type": "Morgan bit fingerprint", "radius": FP_RADIUS, "bits": FP_SIZE, "chirality": False},
        "metric": "Tanimoto",
        "top_hits_saved_per_pyz_target": args.top_per_target,
        "threshold_counts": [0.35, 0.50, 0.70],
        "similarity_bands": {
            "high_2d_similarity": ">=0.70",
            "moderate_2d_similarity": "0.50-0.70",
            "low_notable_2d_similarity": "0.35-0.50",
            "weak_2d_similarity": "<0.35",
        },
        "scientific_guardrail": (
            "Structural similarity to an experimentally measured ligand generates a target hypothesis only. "
            "It does not establish binding, mechanism, affinity, selectivity or antitumor activity for a PYZ compound."
        ),
        "potency_guardrail": (
            "lowest_reported_value_nm_across_endpoints is descriptive only; Ki, Kd and IC50 are retained separately and are not treated as interchangeable."
        ),
    }
    MANIFEST_OUT.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    comparable = summary[summary["best_tanimoto_morgan_r2_2048"].notna()].copy()
    global_best = comparable.sort_values("best_tanimoto_morgan_r2_2048", ascending=False).head(30)

    print("MCL PYZ × Target Ligand Space Similarity v1")
    print(f"PYZ compounds: {len(own_cache)}")
    print(f"Targets: {len(all_targets)}")
    print(f"Unique experimental target-ligand structures: {len(valid_catalog)}")
    print(f"PYZ × target summaries: {len(summary)}")
    print(f"Saved top ligand neighbors: {len(top_hits)}")
    print("Global best PYZ × target ligand neighbors:")
    for row in global_best.to_dict("records"):
        value = float(row["best_tanimoto_morgan_r2_2048"])
        potency = row.get("best_lowest_reported_value_nm_across_endpoints")
        endpoint = _text(row.get("best_lowest_value_endpoint_type"))
        potency_text = f"; lowest reported {endpoint}={float(potency):g} nM" if pd.notna(potency) and endpoint else ""
        print(
            f"  {row['own_compound_id']} -> {row['target_gene']}: "
            f"Tanimoto={value:.3f} [{row['best_similarity_band']}] "
            f"ligand={row['best_ligand_id']}{potency_text}"
        )
    print(f"Summary: {SUMMARY_TSV_OUT.relative_to(ROOT)}")
    print(f"Top hits: {TOP_HITS_TSV_OUT.relative_to(ROOT)}")
    print(f"QC: {QC_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
