#!/usr/bin/env python3
"""P0.3.6: deterministic 13-core + 5–10 expanded reviewed-human pilot panel.

Selecting a protein for the pilot does NOT mean it has a verified ligand.
Observed rows from an uncurated source are marked separately and cannot be
used as proof of physical binding, efficacy, or final catalogue coverage.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as src:
        for part in iter(lambda: src.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def preliminary_counts(path: Path) -> dict[str, int]:
    if not path.is_file():
        raise FileNotFoundError(path)
    frame = pq.ParquetFile(path)
    if not {"target_gene", "source_record_id"}.issubset(frame.schema_arrow.names):
        raise ValueError("Preliminary measurements lack target_gene/source_record_id")
    counts: dict[str, int] = {}
    for batch in frame.iter_batches(batch_size=8192, columns=["target_gene", "source_record_id"]):
        for row in batch.to_pylist():
            gene = str(row["target_gene"] or "").strip().upper()
            if gene and str(row["source_record_id"] or "").strip():
                counts[gene] = counts.get(gene, 0) + 1
    return counts


def create_panel(
    baseline_path: Path, registry_path: Path, *, additional_n: int = 7,
    extra_genes: tuple[str, ...] = (), measured_rows: Path | None = None,
) -> tuple[list[dict[str, str]], dict]:
    if not (5 <= additional_n <= 10):
        raise ValueError("Additional target count must be 5..10")
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    if baseline.get("contract") != "mcl-expanded-target-space-baseline-audit-v1":
        raise ValueError("Unrecognized P0.1 baseline contract")
    core = baseline.get("core_targets", {}).get("target_genes") or []
    core = sorted({str(g).strip().upper() for g in core})
    if len(core) != 13 or baseline.get("own_compounds", {}).get("compound_ids_n") != 65:
        raise ValueError("Real P0.1 audit must confirm 13 core targets and 65 PYZ")
    if not baseline.get("core_targets", {}).get("matches_expected_reference_count") or not baseline.get("own_compounds", {}).get("matches_expected_reference_count"):
        raise ValueError("P0.1 count checks have not passed")
    reg = pd.read_csv(registry_path, sep="\t", dtype=str, keep_default_na=False)
    required = {"target_gene", "uniprot_accession", "mapping_status", "organism_id"}
    if not required.issubset(reg.columns):
        raise ValueError("P0.2 registry schema incomplete")
    good = reg.loc[
        reg["mapping_status"].eq("unique_reviewed_human")
        & reg["organism_id"].eq("9606")
    ].copy()
    if good["target_gene"].duplicated().any() or good["uniprot_accession"].duplicated().any():
        raise ValueError("Ambiguous reviewed target gene/UniProt mapping")
    mapping = dict(zip(good["target_gene"], good["uniprot_accession"]))
    missing = sorted(set(core) - set(mapping))
    if missing:
        raise ValueError("Some core targets are not uniquely mapped in P0.2: " + ",".join(missing))
    outside = sorted(set(mapping) - set(core))
    counts = preliminary_counts(measured_rows) if measured_rows is not None else None

    explicit = [str(g).strip().upper() for g in extra_genes]
    if explicit:
        if len(explicit) != additional_n or len(set(explicit)) != additional_n:
            raise ValueError("Explicit extra genes must be distinct and match additional_n")
        if any(g not in outside for g in explicit):
            raise ValueError("Explicit extra gene is core or absent from pinned reviewed universe")
        expanded = sorted(explicit)
        basis = "researcher_selected_reviewed_protein"
    elif counts is not None:
        expanded = sorted(outside, key=lambda gene: (-counts.get(gene, 0), gene))[:additional_n]
        basis = "observed_source_rows_only_not_curated"
    else:
        expanded = outside[:additional_n]
        basis = "deterministic_gene_order_no_ligand_evidence"
    if len(expanded) != additional_n:
        raise ValueError(
            f"Only {len(expanded)} eligible non-core proteins; requested {additional_n}"
        )
    rows: list[dict[str, str]] = []
    for gene in core + expanded:
        status = "source_rows_not_evaluated" if counts is not None else "not_evaluated"
        rows.append({
            "target_gene": gene,
            "uniprot_accession": mapping[gene],
            "scope": "core" if gene in core else "expanded",
            "selected_by": "p01_measured_core" if gene in core else basis,
            "source_observed_measurement_rows": str(counts.get(gene, 0)) if counts is not None else "",
            "direct_assay_evidence_status": status,
            "verified_ligand_count": "",
            "scientific_note": (
                "Control from P0.1; no response to any PYZ inferred" if gene in core
                else "Pilot candidate only; actual ligand availability remains unverified"
            ),
        })
    report = {
        "contract": "mcl-tls-v2-p036-pilot-panel-v1",
        "core_genes_n": len(core),
        "expanded_genes_n": len(expanded),
        "total_proteins_n": len(rows),
        "expanded_selection_basis": basis,
        "baseline_sha256": sha256(baseline_path),
        "registry_sha256": sha256(registry_path),
        "preliminary_measurements_sha256": sha256(measured_rows) if measured_rows else None,
        "expanded_gene_symbols": expanded,
        "evidence_not_verified": True,
        "note": "Counts of source rows cannot establish independently confirmed A/B ligands.",
    }
    return rows, report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--additional-n", type=int, default=7)
    parser.add_argument("--extra-genes", nargs="*", default=[])
    parser.add_argument("--preliminary-measurements", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.resolve().is_relative_to(args.baseline.resolve().parent):
        raise SystemExit("Pilot outputs must not overwrite the source audit directory.")
    rows, report = create_panel(
        args.baseline, args.registry, additional_n=args.additional_n,
        extra_genes=tuple(args.extra_genes), measured_rows=args.preliminary_measurements
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out = args.output_dir / "pilot_panel.tsv"
    if out.exists() or (args.output_dir / "pilot_panel_manifest.json").exists():
        raise FileExistsError("The pilot panel is immutable; choose a new run directory")
    with out.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    report["pilot_panel_sha256"] = sha256(out)
    (args.output_dir / "pilot_panel_manifest.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"Pilot panel: {len(rows)} reviewed proteins (13 core, {args.additional_n} expanded)")
    print("Experimental ligand evidence: NOT YET VERIFIED")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
