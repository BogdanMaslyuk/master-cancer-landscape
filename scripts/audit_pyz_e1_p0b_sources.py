"""E1-P0-B read-only inventory of existing PYZ reference-ligand assay provenance.

This does NOT approve assays or infer experimental PYZ activity. It links
existing v1 rows only through exact target gene + standardized ligand InChIKey,
then exposes primary-source gaps for manual review. No online requests.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

FOCUS = (("PYZ-035", "CDK4"), ("PYZ-055", "BCL2L1"))
TOP_FILES = (
    "data/runtime/own_compounds/target_ligand_space_similarity_top_hits.parquet",
    "data/runtime/own_compounds/target_ligand_space_similarity_top_hits.tsv",
)
CATALOG_FILES = (
    "data/runtime/target_ligand_space/ligand_catalog.parquet",
    "data/runtime/target_ligand_space/ligand_catalog.tsv",
)
MEASUREMENT_FILES = (
    "data/runtime/target_ligand_space/ligand_measurements.parquet",
    "data/runtime/target_ligand_space/ligand_measurements.tsv",
)
SKIP_DIRS = {".git", ".venv", "venv", "node_modules", ".next", "__pycache__",
             "browser_profiles", "EdgeOptimizationGuideModelsManifest", "ui"}
SOURCE_FILE_RX = re.compile(
    r"(bindingdb|chembl|assay|rsid_eaids|source_manifest|ligand_measurements|"
    r"ligand_catalog|target_ligand_space_similarity_top_hits|manual_scientific_review)",
    flags=re.I,
)
NEIGHBOR_FIELDS = [
    "own_compound_id", "target_gene", "rank_in_target", "ligand_id",
    "tanimoto", "ligand_inchikey", "catalog_inchikey", "measurement_rows_n",
    "assay_ids_present_n", "assay_ids_missing_n", "link_status",
]
MEASUREMENT_FIELDS = [
    "own_compound_id", "target_gene", "rank_in_target", "ligand_id",
    "ligand_inchikey", "source", "source_ligand_id", "source_record_id",
    "uniprot_accession", "endpoint_type", "relation", "value_nm",
    "assay_id", "document_id", "source_url", "source_evidence_class",
    "provenance_status", "measurement_owner", "scientific_review_status",
]
FILE_FIELDS = ["path", "bytes", "sha256", "sha256_status"]


def _str(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (ValueError, TypeError):
        pass
    return str(value).strip()


def _find(root: Path, relative_paths: tuple[str, ...]) -> Path | None:
    return next((root / rel for rel in relative_paths if (root / rel).is_file()), None)


def _load(path: Path, *, genes: set[str] | None = None) -> pd.DataFrame:
    if path.suffix.lower() == ".parquet":
        if genes and path.name.startswith("ligand_measurements"):
            try:
                return pd.read_parquet(path, filters=[("target_gene", "in", sorted(genes))])
            except (ValueError, TypeError, KeyError):
                # Fall back for nonstandard Parquet metadata or older Arrow.
                pass
        frame = pd.read_parquet(path)
    else:
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if genes and "target_gene" in frame.columns:
        return frame.loc[frame["target_gene"].astype(str).str.upper().isin(genes)].copy()
    return frame


def _require_columns(frame: pd.DataFrame, fields: set[str], source: str) -> None:
    absent = sorted(fields - set(frame.columns))
    if absent:
        raise ValueError(f"{source}: missing required columns: {absent}")


def _sha_small(path: Path, limit: int = 32 * 1024 * 1024) -> tuple[str, str]:
    if path.stat().st_size > limit:
        return "", "not_computed_large_source"
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(2 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest(), "computed"


def discover_sources(roots: list[Path], *, max_entries: int = 500) -> list[dict[str, Any]]:
    """List likely original archives/staging reports without reading their content."""
    found: dict[str, dict[str, Any]] = {}
    for root in roots:
        if not root.exists():
            continue
        locations = [root] if root.is_file() else [
            root / name for name in ("data", "outputs", "documents", "sources")
            if (root / name).is_dir()
        ]
        if not locations and root.is_dir():
            locations = [root]
        for location in locations:
            if location.is_file():
                candidates = [location]
            else:
                candidates = []
                for directory, dirs, filenames in os.walk(location):
                    depth = len(Path(directory).relative_to(location).parts)
                    dirs[:] = [d for d in dirs if d not in SKIP_DIRS] if depth < 9 else []
                    for name in filenames:
                        if SOURCE_FILE_RX.search(name):
                            candidates.append(Path(directory) / name)
                            if len(candidates) >= max_entries:
                                break
                    if len(candidates) >= max_entries:
                        break
            for file in candidates:
                if not SOURCE_FILE_RX.search(file.name):
                    continue
                key = str(file.resolve()).casefold()
                if key in found:
                    continue
                size = file.stat().st_size
                # Preserve full path for operator; large files are not hashed on discovery.
                found[key] = {"path": str(file.resolve()), "bytes": size,
                              "sha256": "", "sha256_status": "not_computed_discovery"}
    return sorted(found.values(), key=lambda x: x["path"].casefold())


def prepare(
    hits: pd.DataFrame,
    catalog: pd.DataFrame | None,
    measurements: pd.DataFrame | None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    _require_columns(hits, {"own_compound_id", "target_gene", "ligand_id", "rank_in_target"}, "top hits")
    if catalog is not None:
        _require_columns(catalog, {"ligand_id", "target_gene", "inchikey"}, "ligand catalog")
    if measurements is not None:
        _require_columns(measurements, {"target_gene", "inchikey", "source", "endpoint_type"}, "ligand measurements")

    by_catalog: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    if catalog is not None:
        for raw in catalog.to_dict("records"):
            by_catalog[(_str(raw.get("target_gene")).upper(), _str(raw.get("ligand_id")))].append(raw)

    by_structure: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    if measurements is not None:
        for raw in measurements.to_dict("records"):
            gene, ikey = _str(raw.get("target_gene")).upper(), _str(raw.get("inchikey"))
            if gene and ikey:
                by_structure[(gene, ikey)].append(raw)

    neighbors: list[dict[str, Any]] = []
    evidence: list[dict[str, Any]] = []
    selected_n = {}
    for compound, gene in FOCUS:
        subset = hits.loc[
            hits["own_compound_id"].astype(str).str.upper().eq(compound)
            & hits["target_gene"].astype(str).str.upper().eq(gene)
        ].copy()
        subset["_rank"] = pd.to_numeric(subset["rank_in_target"], errors="coerce")
        subset = subset.loc[subset["_rank"].notna()].sort_values(["_rank", "ligand_id"]).head(5)
        selected_n[f"{compound}/{gene}"] = int(len(subset))
        for hit in subset.to_dict("records"):
            ligand_id = _str(hit.get("ligand_id"))
            hit_ikey = _str(hit.get("ligand_inchikey"))
            catalog_rows = by_catalog.get((gene, ligand_id), [])
            catalog_ikey = _str(catalog_rows[0].get("inchikey")) if len(catalog_rows) == 1 else ""
            base = {
                "own_compound_id": compound,
                "target_gene": gene,
                "rank_in_target": int(hit["_rank"]),
                "ligand_id": ligand_id,
                "tanimoto": _str(hit.get("tanimoto_morgan_r2_2048")),
                "ligand_inchikey": hit_ikey,
                "catalog_inchikey": catalog_ikey,
                "measurement_rows_n": 0,
                "assay_ids_present_n": 0,
                "assay_ids_missing_n": 0,
                "link_status": "",
            }
            if catalog is None:
                base["link_status"] = "CATALOG_UNAVAILABLE"
            elif not catalog_rows:
                base["link_status"] = "CATALOG_LIGAND_NOT_FOUND"
            elif len(catalog_rows) > 1 or not catalog_ikey:
                base["link_status"] = "AMBIGUOUS_CATALOG_KEY"
            elif not hit_ikey:
                base["link_status"] = "TOP_HIT_STRUCTURE_ID_MISSING"
            elif hit_ikey != catalog_ikey:
                base["link_status"] = "INCHIKEY_MISMATCH_NO_LINK"
            elif measurements is None:
                base["link_status"] = "MEASUREMENTS_UNAVAILABLE"
            else:
                records = by_structure.get((gene, catalog_ikey), [])
                if not records:
                    base["link_status"] = "NO_EXACT_MEASUREMENT"
                else:
                    base["measurement_rows_n"] = len(records)
                    for record in records:
                        assay_id = _str(record.get("assay_id"))
                        if assay_id:
                            base["assay_ids_present_n"] += 1
                        else:
                            base["assay_ids_missing_n"] += 1
                        evidence.append({
                            "own_compound_id": compound,
                            "target_gene": gene,
                            "rank_in_target": int(hit["_rank"]),
                            "ligand_id": ligand_id,
                            "ligand_inchikey": catalog_ikey,
                            "source": _str(record.get("source")),
                            "source_ligand_id": _str(record.get("source_ligand_id")),
                            "source_record_id": _str(record.get("source_record_id")),
                            "uniprot_accession": _str(record.get("uniprot_accession")),
                            "endpoint_type": _str(record.get("endpoint_type")),
                            "relation": _str(record.get("relation")),
                            "value_nm": _str(record.get("value_nm")),
                            "assay_id": assay_id,
                            "document_id": _str(record.get("document_id")),
                            "source_url": _str(record.get("source_url")),
                            "source_evidence_class": _str(record.get("evidence_class")),
                            "provenance_status": "SOURCE_ASSAY_ID_UNREVIEWED" if assay_id else "ASSAY_ID_MISSING",
                            "measurement_owner": "known_reference_ligand_not_pyz",
                            "scientific_review_status": "PENDING_PRIMARY_SOURCE_VERIFICATION",
                        })
                    base["link_status"] = "EXACT_STRUCTURE_MATCH_PENDING_SOURCE_REVIEW"
            neighbors.append(base)
    stats = {
        "selected_neighbors_by_focus": selected_n,
        "neighbors_total": len(neighbors),
        "neighbor_status_counts": dict(sorted(Counter(x["link_status"] for x in neighbors).items())),
        "measurement_candidates_n": len(evidence),
        "assay_id_present_rows_n": sum(bool(x["assay_id"]) for x in evidence),
        "assay_id_missing_rows_n": sum(not x["assay_id"] for x in evidence),
        "measurement_endpoint_counts": dict(sorted(Counter(x["endpoint_type"] for x in evidence).items())),
        "measurement_source_counts": dict(sorted(Counter(x["source"] for x in evidence).items())),
        "primary_source_approved_n": 0,
    }
    return neighbors, evidence, stats


def _write_tsv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--science-root", type=Path, required=True, help="PYZ v2 scientific worktree")
    parser.add_argument("--ligand-root", type=Path, action="append", default=[],
                        help="Candidate worktree for same-root ligand_catalog + ligand_measurements; repeatable")
    parser.add_argument("--scan-root", type=Path, action="append", default=[],
                        help="Existing local BindingDB/ChEMBL source/archive directory; repeatable")
    parser.add_argument("--out-dir", type=Path, required=True, help="New isolated output folder")
    args = parser.parse_args()

    science = args.science_root.resolve()
    if not science.is_dir():
        raise SystemExit("Science root does not exist: " + str(science))
    roots = list(dict.fromkeys([science] + [p.resolve() for p in args.ligand_root]))
    output = args.out_dir.resolve()
    if output.exists():
        raise SystemExit("Output already exists; refusing to overwrite: " + str(output))
    for root in roots:
        if output.is_relative_to(root):
            raise SystemExit("Output must be outside source worktrees")
    top_path = _find(science, TOP_FILES)
    pair_root = next((r for r in roots if _find(r, CATALOG_FILES) and _find(r, MEASUREMENT_FILES)), None)
    cat_path = _find(pair_root, CATALOG_FILES) if pair_root else None
    m_path = _find(pair_root, MEASUREMENT_FILES) if pair_root else None

    source_files = {}
    for label, path in [("top_hits", top_path), ("ligand_catalog", cat_path), ("ligand_measurements", m_path)]:
        if path is None:
            source_files[label] = {"available": False, "path": None}
        else:
            digest, status = _sha_small(path)
            source_files[label] = {
                "available": True, "path": str(path), "size_bytes": path.stat().st_size,
                "sha256": digest, "sha256_status": status,
            }

    hits = _load(top_path) if top_path else pd.DataFrame(
        columns=["own_compound_id", "target_gene", "ligand_id", "rank_in_target"])
    catalog = _load(cat_path, genes={"CDK4", "BCL2L1"}) if cat_path else None
    measurements = _load(m_path, genes={"CDK4", "BCL2L1"}) if m_path else None
    neighbors, evidence, stats = prepare(hits, catalog, measurements)

    scan_roots = [*roots, *(p.resolve() for p in args.scan_root)]
    files = discover_sources(scan_roots)
    # Do not claim scientific approval even if an assay identifier is present.
    state = "READY_FOR_MANUAL_REVIEW_OF_SOURCE_ROWS" if (
        len(neighbors) == 10
        and all(row["link_status"] == "EXACT_STRUCTURE_MATCH_PENDING_SOURCE_REVIEW" for row in neighbors)
    ) else "DATA_GAPS_OR_UNRESOLVED_SOURCE_LINKS"

    report = {
        "contract": "mcl-pyz-e1-p0b-source-inventory-v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_status": "NO_PRIMARY_ASSAYS_APPROVED",
        "inventory_status": state,
        "source_root": str(science),
        "ligand_catalog_measurements_root": str(pair_root) if pair_root else None,
        "cross_worktree_ligand_source": bool(pair_root and pair_root != science),
        "cross_worktree_note": (
            "Exact InChIKey joins do not prove source release parity; review manifests and source hashes."
            if pair_root and pair_root != science else ""
        ),
        "source_files": source_files,
        "stats": stats,
        "discovered_files_n": len(files),
        "limitations": [
            "No primary publication or assay has been approved by this audit.",
            "ChEMBL assay IDs from v1 are source metadata, not manual verification.",
            "BindingDB v1 REST rows commonly have blank assay_id; do not invent it.",
            "Top hits are structural analogs of PYZ; listed Ki/Kd/IC50 belong to known ligands.",
            "The v1 ligand catalog merges measurements by target+InChIKey; a source record may not have unique assay metadata.",
            "The P0.3.6 Articles/Assays/rsid_eaids join must be evaluated separately, by original reactant_set ID.",
        ],
    }
    output.mkdir(parents=True, exist_ok=False)
    _write_tsv(output / "neighbors_audit.tsv", NEIGHBOR_FIELDS, neighbors)
    _write_tsv(output / "measurement_source_candidates.tsv", MEASUREMENT_FIELDS, evidence)
    _write_tsv(output / "local_source_inventory.tsv", FILE_FIELDS, files)
    (output / "inventory_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("E1-P0-B source inventory: " + state)
    print("Science root: " + str(science))
    print("Ligand measurements root: " + str(pair_root))
    print(f"Focus neighbors: {len(neighbors)} (expected 10)")
    print(f"Exact-source candidate rows: {len(evidence)}")
    print(f"Rows with assay ID (UNREVIEWED): {stats['assay_id_present_rows_n']}")
    print(f"Rows without assay ID: {stats['assay_id_missing_rows_n']}")
    print(f"Potential archive/staging paths: {len(files)}")
    print("Scientific approval: NONE; manual primary-source review required")
    print("Report: " + str(output / "inventory_report.json"))


if __name__ == "__main__":
    main()
