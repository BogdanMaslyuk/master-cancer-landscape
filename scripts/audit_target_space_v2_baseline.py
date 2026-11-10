#!/usr/bin/env python3
"""Read-only scientific baseline audit for Target Ligand Space v2.

Only the dedicated audit report is written. No legacy scientific file is changed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INPUTS = {
    "candidate_v2": "data/runtime/pharmacology/candidate_hypotheses_v2.parquet",
    "ligand_v1": "data/runtime/target_ligand_space/ligand_catalog.parquet",
    "ligand_manifest_v1": "data/runtime/target_ligand_space/target_ligand_space_manifest.json",
    "own_compounds": "config/own_compounds.tsv",
    "full_depmap_context": "data/processed/depmap_model_target_context.parquet",
}
DEFAULT_OUTPUT = Path("outputs/reports/target_ligand_space_v2/baseline_audit.json")
SMILES_REGRESSION = ("CCO", "c1ccccc1", "CC(=O)O")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def audit(source: Path) -> dict:
    source = source.resolve()
    files = {}
    for key, relative in INPUTS.items():
        path = source / relative
        files[key] = {
            "relative_path": relative,
            "present": path.is_file(),
            "bytes": path.stat().st_size if path.is_file() else None,
            "sha256": digest(path) if path.is_file() else None,
        }

    core = {"status": "not_available", "target_genes": [], "target_genes_n": None}
    candidate = source / INPUTS["candidate_v2"]
    if candidate.is_file():
        frame = pd.read_parquet(candidate, columns=["target_gene", "priority_status"])
        priority = frame.loc[
            frame["priority_status"].astype(str).eq("priority_for_in_vitro"),
            "target_gene",
        ]
        genes = sorted({str(g).strip().upper() for g in priority.dropna() if str(g).strip()})
        core = {
            "status": "measured_from_candidate_v2",
            "target_genes": genes,
            "target_genes_n": len(genes),
            "expected_reference_count": 13,
            "matches_expected_reference_count": len(genes) == 13,
        }

    own = {"status": "not_available", "compound_ids_n": None}
    registry = source / INPUTS["own_compounds"]
    if registry.is_file():
        frame = pd.read_csv(registry, sep="\t", dtype=str)
        if "own_compound_id" not in frame.columns:
            raise ValueError("own_compounds.tsv missing own_compound_id")
        ids = frame["own_compound_id"].dropna().str.strip()
        own = {
            "status": "measured_from_own_compounds",
            "compound_ids_n": int(ids.nunique()),
            "duplicate_ids_n": int(ids.duplicated().sum()),
            "expected_reference_count": 65,
            "matches_expected_reference_count": ids.nunique() == 65 and not ids.duplicated().any(),
        }

    readiness = {
        "core_count_check": core.get("matches_expected_reference_count") is True,
        "pyz_count_check": own.get("matches_expected_reference_count") is True,
        "ligand_v1_present": files["ligand_v1"]["present"],
        "full_depmap_context_present": files["full_depmap_context"]["present"],
    }
    return {
        "contract": "mcl-expanded-target-space-baseline-audit-v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "science_root": str(source),
        "files": files,
        "core_targets": core,
        "own_compounds": own,
        "regression_smiles": list(SMILES_REGRESSION),
        "checks": readiness,
        "ready_for_extended_build": all(readiness.values()),
        "limitations": [
            "File hashes and counts are not a scientific validation of source measurements.",
            "SMILES entries are fixed regression inputs; this audit does not execute molecule analysis.",
            "Missing local files are reported as missing, never interpreted as empty scientific datasets.",
        ],
    }


def write_report(report: dict, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    temporary.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    temporary.replace(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, default=ROOT / DEFAULT_OUTPUT)
    parser.add_argument("--strict", action="store_true", help="Exit nonzero if baseline is incomplete.")
    args = parser.parse_args()
    result = audit(args.data_root)
    write_report(result, args.output)
    print(f"Audit: {args.output}; ready={result['ready_for_extended_build']}")
    for label, ok in result["checks"].items():
        print(f"  {'OK' if ok else 'MISSING/FAIL'} {label}")
    if args.strict and not result["ready_for_extended_build"]:
        raise SystemExit("Baseline incomplete; no scientific source files were changed.")


if __name__ == "__main__":
    main()
