"""MCL E1-P0-B3: read-only evidence triage of the E1-P0-B2 BindingDB result ZIP.

Does NOT approve assays, infer PYZ activity or interpret ADMET data.
Only Python stdlib required. Data remain in the existing MCL_Project.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

REQUIRED_MEMBERS = (
    "bindingdb_full_report.json", "bindingdb_full_matches.tsv",
    "bindingdb_full_neighbor_coverage.tsv", "run_manifest.json",
)
REGISTER_COLUMNS = [
    "record_id", "own_compound_id", "target_gene", "ligand_id",
    "rank_in_target", "tanimoto", "review_priority",
    "bindingdb_reactant_set_id", "bindingdb_monomer_id", "source_row_number",
    "endpoint", "raw_nm", "relation", "value_nm", "measurement_owner",
    "target_uniprot", "all_target_chain_uniprot", "target_name", "target_complexity",
    "expected_inchikey", "bindingdb_inchikey", "computed_inchikey",
    "structure_match", "chemical_identity_flag", "source_key_agrees_with_computed",
    "article_doi", "pmid", "curation_source", "potential_duplicate_group",
    "possible_duplicate_rows_in_group", "biological_interpretation_guard",
    "source_record_status", "primary_method_review", "structure_review",
    "independent_assay_review", "scientific_approval",
]
PAPER_COLUMNS = [
    "article_doi", "pmid", "record_count", "unique_ligands",
    "targets", "endpoint_types", "priority_1_records", "review_status",
    "method_check", "compound_identity_check", "publication_reference_check",
    "independence_check", "notes",
]
CHECK_COLUMNS = [
    "record_id", "target_gene", "ligand_id", "bindingdb_reactant_set_id",
    "article_doi", "priority", "source_structure_check", "original_paper_structure_check",
    "stereo_isomer_check", "recombinant_complex_check", "assay_protocol_check",
    "concentration_and_units_check", "replicates_error_check", "source_record_alignment_check",
    "independent_assay_check", "mechanism_vs_target_claim_check",
    "reviewer", "review_date_utc", "evidence_location", "decision", "rationale",
]
REQUIRED_ROWS = {
    "own_compound_id", "target_gene", "ligand_id", "rank_in_target", "tanimoto",
    "expected_inchikey", "bindingdb_inchikey", "computed_inchikey",
    "structure_match", "target_uniprot", "all_target_chain_uniprot",
    "target_complexity", "target_name", "bindingdb_reactant_set_id",
    "bindingdb_monomer_id", "endpoint", "raw_nm", "relation", "value_nm",
    "article_doi", "pmid", "curation_source", "source_row_number",
    "measurement_owner", "scientific_approval",
}
INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
FOCUS = {("PYZ-035", "CDK4"), ("PYZ-055", "BCL2L1")}
TARGETS = {"CDK4": "P11802", "BCL2L1": "Q07817"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_zip_json(z: zipfile.ZipFile, name: str) -> dict:
    value = json.loads(z.read(name).decode("utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {name}")
    return value


def _read_zip_tsv(z: zipfile.ZipFile, name: str, required: set[str]) -> list[dict]:
    reader = csv.DictReader(io.StringIO(z.read(name).decode("utf-8-sig")), delimiter="\t")
    if not required.issubset(reader.fieldnames or []):
        raise ValueError(f"Missing columns in {name}: {sorted(required - set(reader.fieldnames or []))}")
    rows = list(reader)
    if any(None in item for item in rows):
        raise ValueError(f"Bad tabular fields in {name}")
    return rows


def load(source_zip: Path):
    with zipfile.ZipFile(source_zip) as z:
        names = z.namelist()
        for required in REQUIRED_MEMBERS:
            if names.count(required) != 1:
                raise ValueError(f"Expected exactly one file {required} in ZIP")
        report = _read_zip_json(z, "bindingdb_full_report.json")
        manifest = _read_zip_json(z, "run_manifest.json")
        rows = _read_zip_tsv(z, "bindingdb_full_matches.tsv", REQUIRED_ROWS)
        coverage = _read_zip_tsv(z, "bindingdb_full_neighbor_coverage.tsv", {
            "own_compound_id", "target_gene", "ligand_id", "rank_in_target", "exact_endpoint_rows"})
    if report.get("contract") != "mcl-pyz-e1-p0b-full-bindingdb-readonly-v2":
        raise ValueError("Unsupported full BindingDB report contract")
    if report.get("scientific_status") != "NO_PRIMARY_ASSAYS_APPROVED" or report.get("primary_assays_approved") != 0:
        raise ValueError("Input scientific status is not the expected unapproved baseline")
    if manifest.get("scientific_status") != "NO_PRIMARY_ASSAYS_APPROVED":
        raise ValueError("Unexpected run_manifest scientific status")
    if manifest.get("stage") != "E1_P0B2":
        raise ValueError("Unexpected source run stage")
    input_hashes = report.get("input_sha256", {})
    source_hash = str(manifest.get("bindingdb_sha256", "")).lower()
    if len(source_hash) != 64 or source_hash not in {str(v).lower() for v in input_hashes.values()}:
        raise ValueError("Original full BindingDB SHA-256 disagrees with report")
    neighbors_hash = str(manifest.get("neighbors_sha256", "")).lower()
    if len(neighbors_hash) != 64 or neighbors_hash not in {str(v).lower() for v in input_hashes.values()}:
        raise ValueError("Original neighbors SHA-256 disagrees with report")
    if int(report.get("neighbors", -1)) != len(coverage) or len(coverage) != 10:
        raise ValueError("Expected 10 source neighbors")
    if int(report.get("stats", {}).get("matched_endpoint_rows", -1)) != len(rows):
        raise ValueError("Endpoint row count does not match source report")
    if len(rows) != 11:
        raise ValueError("This E1-P0-B3 baseline requires 11 original endpoint rows")
    pairs = Counter((r["own_compound_id"], r["target_gene"]) for r in coverage)
    if pairs != Counter({("PYZ-035", "CDK4"): 5, ("PYZ-055", "BCL2L1"): 5}):
        raise ValueError("Incorrect source neighbor pairs")
    neighbor_keys = [(r["own_compound_id"], r["target_gene"], r["ligand_id"]) for r in coverage]
    if len(neighbor_keys) != len(set(neighbor_keys)):
        raise ValueError("Duplicated neighbor IDs in coverage")
    if any((r["target_gene"], r["own_compound_id"]) not in
           {("CDK4", "PYZ-035"), ("BCL2L1", "PYZ-055")} for r in rows):
        raise ValueError("Unexpected compound-target mapping")
    measured = Counter((r["own_compound_id"], r["target_gene"], r["ligand_id"]) for r in rows)
    for row in coverage:
        key = (row["own_compound_id"], row["target_gene"], row["ligand_id"])
        if int(row["exact_endpoint_rows"]) != measured[key]:
            raise ValueError(f"Coverage endpoint count mismatch for {key}")
    if sum(bool(measured[key]) for key in neighbor_keys) != int(report.get("neighbors_exact_endpoint_coverage", -1)):
        raise ValueError("Exact neighbor coverage differs from source report")
    return report, manifest, rows, coverage


def analyze(source_zip: Path) -> tuple[list[dict], list[dict], list[dict], dict]:
    report, manifest, rows, coverage = load(source_zip)
    source_hash = manifest["bindingdb_sha256"].lower()
    monomer_groups = Counter((r["target_gene"], r["bindingdb_monomer_id"], r["endpoint"]) for r in rows)
    papers = defaultdict(list)
    registry = []
    seen = set()
    neighbor_keys = {(r["own_compound_id"], r["target_gene"], r["ligand_id"]) for r in coverage}
    stats = Counter()
    for row in rows:
        pair = (row["own_compound_id"], row["target_gene"])
        if pair not in FOCUS or row["target_uniprot"].upper() != TARGETS[row["target_gene"]]:
            raise ValueError(f"Unexpected PYZ/target/UniProt relation: {pair}")
        if (*pair, row["ligand_id"]) not in neighbor_keys:
            raise ValueError("Measurement not linked to original neighbor coverage")
        if row["endpoint"] not in ("Ki", "Kd", "IC50") or row["relation"] not in ("=", "<", ">", "<=", ">="):
            raise ValueError("Unexpected endpoint or relation")
        if row["target_uniprot"].upper() not in row["all_target_chain_uniprot"].upper().split(";"):
            raise ValueError("Target accession not found among assay chain identifiers")
        if row["measurement_owner"] != "known_reference_ligand_not_pyz":
            raise ValueError("PYZ experimental ownership is forbidden")
        if row["scientific_approval"] != "NOT_APPROVED_MANUAL_REVIEW_REQUIRED":
            raise ValueError("Unexpected source scientific approval")
        if not all(INCHIKEY.fullmatch(row[k]) for k in ("expected_inchikey", "computed_inchikey")):
            raise ValueError("Missing or invalid standardized InChIKey")
        if row["computed_inchikey"] != row["expected_inchikey"] or row["structure_match"] != "EXACT_STANDARDIZED_INCHIKEY":
            raise ValueError("Source row is not an exact computed structural match")
        if row["target_gene"] == "CDK4" and row["endpoint"] != "IC50":
            raise ValueError("Unexpected CDK4 endpoint in frozen baseline")
        if row["target_gene"] == "BCL2L1" and row["endpoint"] != "Ki":
            raise ValueError("Unexpected BCL2L1 endpoint in frozen baseline")
        if not row["article_doi"].strip() or not row["bindingdb_reactant_set_id"].isdigit() or not row["source_row_number"].isdigit():
            raise ValueError("Missing source provenance identifiers")
        key = (row["bindingdb_reactant_set_id"], row["target_gene"], row["endpoint"], row["source_row_number"])
        if key in seen:
            raise ValueError(f"Duplicate source record: {key}")
        seen.add(key)
        rid = "E1B3-" + hashlib.sha256("|".join([source_hash, *key]).encode()).hexdigest()[:16].upper()
        source_key_match = row["bindingdb_inchikey"] == row["computed_inchikey"]
        chemical = "SOURCE_INCHIKEY_MATCHES_COMPUTED" if source_key_match else "SOURCE_INCHIKEY_DIFF_REVIEW_STEREO_OR_STANDARDIZATION"
        dup_key = (row["target_gene"], row["bindingdb_monomer_id"], row["endpoint"])
        dup_count = monomer_groups[dup_key]
        priority = "P1_TARGETED_SOURCE_REVIEW" if int(row["rank_in_target"]) == 1 else (
            "P2_REFERENCE_CROSSCHECK" if dup_count > 1 else "P3_REMAINING")
        interpretation = ("KINASE_INHIBITION_NOT_DIRECT_BINDING_OR_PYZ_ACTIVITY" if row["target_gene"] == "CDK4"
                          else "BCLXL_KI_ASSAY_MECHANISM_UNREVIEWED_NOT_PYZ_ACTIVITY")
        record = {
            **{col: row.get(col, "") for col in REGISTER_COLUMNS},
            "record_id": rid, "review_priority": priority,
            "chemical_identity_flag": chemical,
            "source_key_agrees_with_computed": "YES" if source_key_match else "NO",
            "potential_duplicate_group": ("MONOMER_TARGET_ENDPOINT_SHARED" if dup_count > 1 else "NONE_DETECTED_IN_11_ROWS"),
            "possible_duplicate_rows_in_group": dup_count,
            "biological_interpretation_guard": interpretation,
            "source_record_status": "SOURCE_ROW_TRACED_UNREVIEWED",
            "primary_method_review": "NOT_REVIEWED",
            "structure_review": "NOT_REVIEWED",
            "independent_assay_review": "NOT_REVIEWED",
            "scientific_approval": "NO_PRIMARY_ASSAYS_APPROVED",
        }
        registry.append(record)
        papers[row["article_doi"].lower()].append(record)
        stats["records"] += 1
        stats["source_inchikey_discrepancy"] += not source_key_match
        stats["multichain_flagged"] += row["target_complexity"] == "MULTI_CHAIN_REVIEW_REQUIRED"
        stats["potential_duplicate_rows"] += dup_count > 1
    registry.sort(key=lambda r: (r["own_compound_id"], int(r["rank_in_target"]), r["bindingdb_reactant_set_id"]))
    paper_queue = []
    for doi, grouped in sorted(papers.items()):
        paper_queue.append({
            "article_doi": doi,
            "pmid": ";".join(sorted({r["pmid"] for r in grouped if r["pmid"]})),
            "record_count": len(grouped),
            "unique_ligands": len({r["bindingdb_monomer_id"] for r in grouped}),
            "targets": ";".join(sorted({r["target_gene"] for r in grouped})),
            "endpoint_types": ";".join(sorted({r["endpoint"] for r in grouped})),
            "priority_1_records": sum(r["review_priority"] == "P1_TARGETED_SOURCE_REVIEW" for r in grouped),
            "review_status": "NOT_REVIEWED",
            "method_check": "PENDING_FULL_TEXT_OR_SUPPORTING_INFO",
            "compound_identity_check": "PENDING_PUBLISHED_STRUCTURE",
            "publication_reference_check": "PENDING_LINKED_METADATA_RECONCILIATION",
            "independence_check": "PENDING_EXPERIMENT_LEVEL_RECONCILIATION",
            "notes": "One DOI is not proof of independent assay replicates",
        })
    checklist = []
    for record in registry:
        item = {
            "record_id": record["record_id"], "target_gene": record["target_gene"],
            "ligand_id": record["ligand_id"],
            "bindingdb_reactant_set_id": record["bindingdb_reactant_set_id"],
            "article_doi": record["article_doi"], "priority": record["review_priority"],
            "reviewer": "", "review_date_utc": "", "evidence_location": "",
            "decision": "NOT_REVIEWED", "rationale": "",
        }
        for k in CHECK_COLUMNS:
            if k.endswith("_check"):
                item[k] = "NOT_REVIEWED"
        checklist.append(item)
    stats["papers"] = len(papers)
    stats["distinct_reference_ligands"] = len({(r["target_gene"], r["bindingdb_monomer_id"]) for r in registry})
    stats["reference_ligands_above_035"] = len({(r["target_gene"], r["ligand_id"]) for r in registry if float(r["tanimoto"]) >= 0.35})
    stats["primary_assays_approved"] = 0
    result = {
        "contract": "mcl-pyz-e1-p0b3-assay-evidence-review-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_artifact_sha256": sha256(source_zip),
        "full_bindingdb_source_sha256": source_hash,
        "source_neighbors_sha256": manifest["neighbors_sha256"].lower(),
        "source_release_file": report["source_release_file"],
        "source_measurement_owner": "known_reference_ligand_not_pyz",
        "scientific_status": "NO_PRIMARY_ASSAYS_APPROVED",
        "statistics": dict(sorted(stats.items())),
        "limitations": [
            "Review flags are metadata triage, not a primary publication or assay acceptance.",
            "Computed standardized InChIKey may disagree with source InChIKey and original paper stereochemistry.",
            "A shared monomer and endpoint is only a potential duplicate, not proof of identical experiments.",
            "Paper abstracts do not verify concentration, constructs, protocol details or replication.",
            "Only 11 records in this source artifact were reviewed; no conclusion about other BindingDB entries.",
            "No ADMET scores or claims are used.",
            "No known-ligand endpoint is evidence of an experimental PYZ result.",
        ],
    }
    return registry, paper_queue, checklist, result


def write_tsv(path: Path, fields: list[str], values: list[dict]):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(values)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-zip", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    if not args.source_zip.is_file():
        raise FileNotFoundError(args.source_zip)
    destination = args.out_dir.resolve()
    if destination.exists():
        raise FileExistsError(f"Output folder exists: {destination}")
    if args.source_zip.resolve().is_relative_to(destination):
        raise ValueError("Output must not contain the input source ZIP")
    reg, papers, checklist, report = analyze(args.source_zip)
    destination.mkdir(parents=True)
    write_tsv(destination / "assay_evidence_register.tsv", REGISTER_COLUMNS, reg)
    write_tsv(destination / "publication_review_queue.tsv", PAPER_COLUMNS, papers)
    write_tsv(destination / "manual_source_review_checklist.tsv", CHECK_COLUMNS, checklist)
    (destination / "evidence_review_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"records": len(reg), "papers": len(papers),
                      "source_key_discrepancies": report["statistics"]["source_inchikey_discrepancy"],
                      "approved_primary_assays": 0}, ensure_ascii=False))
    return report


if __name__ == "__main__":
    main()