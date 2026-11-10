#!/usr/bin/env python3
"""Build an independent, auditable human protein target universe (P0.2).

This is a *registry*, not a target prediction engine. It deliberately does not
filter by Candidate v2 / priority_for_in_vitro. All outputs live in v2 paths.
The optional UniProt TSV must be an official, frozen reviewed-human export.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_V1 = ROOT / "data/processed/target_registry/protein_targets.parquet"
DEFAULT_OUTPUT = ROOT / "data/processed/target_ligand_space_v2"
ACCESSION = re.compile(r"^[A-Z0-9]{6,10}$")
GENE = re.compile(r"^[A-Z][A-Z0-9.-]*$")


def _value(item: object) -> str:
    if item is None or pd.isna(item):
        return ""
    return str(item).strip()


def _first(row: dict, *fields: str) -> str:
    for name in fields:
        value = _value(row.get(name))
        if value:
            return value
    return ""


def _digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _read_frame(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    if path.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    raise ValueError(f"Unsupported input: {path}; expected .tsv or .parquet")


def collect(source_registry: pd.DataFrame, uniprot: pd.DataFrame | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return eligible targets and QC cases, preserving source-level provenance."""
    if not {"target_gene", "protein_mapping_status", "uniprot_primary_accession"}.issubset(source_registry):
        raise ValueError("v1 registry lacks target_gene/protein_mapping_status/uniprot_primary_accession")
    evidence: list[dict] = []
    rejected: list[dict] = []

    for row in source_registry.to_dict("records"):
        gene = _value(row.get("target_gene")).upper()
        acc = _value(row.get("uniprot_primary_accession")).upper()
        status = _value(row.get("protein_mapping_status"))
        if status != "unique_swissprot":
            rejected.append({"source": "existing_reviewed_registry", "target_gene": gene,
                             "uniprot_accession": acc, "reason_code": "mapping_not_unique"})
            continue
        evidence.append({
            "target_gene": gene, "uniprot_accession": acc,
            "protein_name": _value(row.get("protein_preferred_name")),
            "source": "existing_reviewed_registry",
            "uniprot_release": _value(row.get("uniprot_release")),
            "reviewed": True, "organism_id": "9606",
        })

    if uniprot is not None:
        required = {"Entry", "Reviewed", "Gene Names (primary)"}
        if not required.issubset(uniprot.columns):
            raise ValueError(f"UniProt TSV missing required columns: {sorted(required - set(uniprot))}")
        if "Organism (ID)" not in uniprot.columns and "Organism ID" not in uniprot.columns:
            raise ValueError("UniProt TSV must explicitly contain Organism (ID)")
        for row in uniprot.to_dict("records"):
            acc = _first(row, "Entry").upper()
            gene = _first(row, "Gene Names (primary)").upper()
            organism = _first(row, "Organism (ID)", "Organism ID")
            reviewed = _first(row, "Reviewed").lower()
            if organism != "9606":
                reason = "non_human"
            elif reviewed not in {"reviewed", "yes", "true"}:
                reason = "not_reviewed"
            elif len(gene.split()) != 1:
                reason = "ambiguous_primary_gene"
            else:
                reason = ""
            if reason:
                rejected.append({"source": "official_uniprot_snapshot", "target_gene": gene,
                                 "uniprot_accession": acc, "reason_code": reason})
                continue
            evidence.append({
                "target_gene": gene, "uniprot_accession": acc,
                "protein_name": _first(row, "Protein names", "Protein Names"),
                "source": "official_uniprot_snapshot", "uniprot_release": "",
                "reviewed": True, "organism_id": organism,
            })

    valid: list[dict] = []
    for row in evidence:
        gene, acc = row["target_gene"], row["uniprot_accession"]
        if not GENE.fullmatch(gene) or not ACCESSION.fullmatch(acc) or "-" in acc:
            rejected.append({"source": row["source"], "target_gene": gene,
                             "uniprot_accession": acc, "reason_code": "invalid_identifier"})
        else:
            valid.append(row)

    # Conflicting mappings are never resolved by choosing the first row.
    by_gene: dict[str, set[str]] = {}
    by_acc: dict[str, set[str]] = {}
    for row in valid:
        by_gene.setdefault(row["target_gene"], set()).add(row["uniprot_accession"])
        by_acc.setdefault(row["uniprot_accession"], set()).add(row["target_gene"])
    accepted: list[dict] = []
    for row in valid:
        gene, acc = row["target_gene"], row["uniprot_accession"]
        if len(by_gene[gene]) != 1 or len(by_acc[acc]) != 1:
            rejected.append({"source": row["source"], "target_gene": gene,
                             "uniprot_accession": acc, "reason_code": "conflicting_gene_protein_mapping"})
        else:
            accepted.append(row)

    groups: dict[tuple[str, str], list[dict]] = {}
    for row in accepted:
        groups.setdefault((row["target_gene"], row["uniprot_accession"]), []).append(row)
    rows = []
    for (gene, acc), group in sorted(groups.items()):
        named = [r["protein_name"] for r in group if r["protein_name"]]
        releases = sorted({r["uniprot_release"] for r in group if r["uniprot_release"]})
        sources = sorted({r["source"] for r in group})
        rows.append({
            "target_gene": gene,
            "uniprot_accession": acc,
            "protein_name": named[0] if named else "",
            "organism_id": "9606",
            "mapping_status": "unique_reviewed_human",
            "reviewed": True,
            "source_names_json": json.dumps(sources),
            "uniprot_releases_json": json.dumps(releases),
            "ligand_coverage_status": "not_evaluated",
            "depmap_coverage_status": "not_evaluated",
            "oncology_context_status": "not_evaluated",
            "direct_target_evidence_status": "not_evaluated",
            "source_resolution": "gene_to_reviewed_protein_mapping",
        })
    cols = ["target_gene", "uniprot_accession", "protein_name", "organism_id",
            "mapping_status", "reviewed", "source_names_json", "uniprot_releases_json",
            "ligand_coverage_status", "depmap_coverage_status", "oncology_context_status",
            "direct_target_evidence_status", "source_resolution"]
    qc_cols = ["source", "target_gene", "uniprot_accession", "reason_code"]
    return pd.DataFrame(rows, columns=cols), pd.DataFrame(rejected, columns=qc_cols)


def build(v1_path: Path, uniprot_path: Path | None, output_dir: Path,
          uniprot_release: str | None = None) -> dict:
    if not v1_path.is_file():
        raise FileNotFoundError(f"Missing locally verified v1 registry: {v1_path}")
    if uniprot_path is not None and not uniprot_path.is_file():
        raise FileNotFoundError(f"Missing pinned UniProt TSV: {uniprot_path}")
    if uniprot_path is not None and not uniprot_release:
        raise ValueError("--uniprot-release is required with --uniprot-tsv")
    v1 = _read_frame(v1_path)
    extra = _read_frame(uniprot_path) if uniprot_path is not None else None
    frame, qc = collect(v1, extra)
    if frame.empty:
        raise ValueError("No uniquely mapped reviewed human proteins; nothing published")
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / "protein_target_universe.tsv"
    qc_path = output_dir / "protein_target_universe_qc.tsv"
    manifest_path = output_dir / "protein_target_universe_manifest.json"
    # Stage temporary files: only publish a successful, internally validated build.
    frame.to_csv(output.with_name(output.name + ".tmp"), sep="\t", index=False)
    qc.to_csv(qc_path.with_name(qc_path.name + ".tmp"), sep="\t", index=False)
    manifest = {
        "contract": "mcl-protein-target-universe-v2-p0",
        "built_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "reviewed_human_single_protein_candidates",
        "universe_is_exhaustive": False,
        "uniprot_snapshot_included": uniprot_path is not None,
        "uniprot_release": uniprot_release,
        "source_files": [{"path": str(p.resolve()), "sha256": _digest(p)}
                         for p in [v1_path] + ([uniprot_path] if uniprot_path else [])],
        "target_genes_n": len(frame),
        "rejected_source_rows_n": len(qc),
        "candidate_v2_filter_applied": False,
        "ligand_or_depmap_coverage_confirmed": False,
        "limitations": [
            "This stage contains mapped proteins, not verified ChEMBL/BindingDB direct assays.",
            "No oncology or DepMap dependency claims can be made from registry membership.",
            "An optional frozen UniProt snapshot is needed to expand beyond previously annotated targets.",
        ],
    }
    manifest_path.with_name(manifest_path.name + ".tmp").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    for path in [output, qc_path, manifest_path]:
        path.with_name(path.name + ".tmp").replace(path)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v1-registry", type=Path, default=DEFAULT_V1)
    parser.add_argument("--uniprot-tsv", type=Path)
    parser.add_argument("--uniprot-release", type=str)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    manifest = build(args.v1_registry, args.uniprot_tsv, args.output_dir, args.uniprot_release)
    print(f"Unique reviewed human mapped proteins: {manifest['target_genes_n']}")
    print("Ligand, DepMap and oncology coverage: NOT YET VERIFIED")
    print(f"Manifest: {args.output_dir / 'protein_target_universe_manifest.json'}")


if __name__ == "__main__":
    main()
