"""Read-only full BindingDB screening for PYZ reference ligands (E1-P0-B2).

Measured values belong to known reference ligands, NEVER to the PYZ compounds.
No primary assay or publication is approved by this script.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import sys
import zipfile
from collections import Counter, defaultdict
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.error")
TARGETS = {"CDK4": "P11802", "BCL2L1": "Q07817"}
FOCUS = {("PYZ-035", "CDK4"), ("PYZ-055", "BCL2L1")}
ENDPOINTS = {"Ki": "Ki (nM)", "Kd": "Kd (nM)", "IC50": "IC50 (nM)"}
INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
NUMERIC = re.compile(r"^\s*(<=|>=|<|>|=|≤|≥)?\s*(\d+(?:\.\d*)?|\.\d+)(?:[eE]([+-]?\d+))?\s*$")
MATCH_FIELDS = [
    "own_compound_id", "target_gene", "ligand_id", "rank_in_target", "tanimoto",
    "expected_inchikey", "bindingdb_inchikey", "computed_inchikey", "structure_match",
    "target_uniprot", "all_target_chain_uniprot", "target_complexity",
    "target_name", "organism", "bindingdb_reactant_set_id",
    "bindingdb_monomer_id", "endpoint", "raw_nm", "relation", "value_nm",
    "endpoint_quality", "article_doi", "pmid", "curation_source",
    "source_row_number", "reference_compound_smiles", "source_inchikey_vs_computed",
    "measurement_owner", "scientific_approval",
]
COVER_FIELDS = [
    "own_compound_id", "target_gene", "ligand_id", "rank_in_target", "tanimoto",
    "ligand_inchikey", "exact_endpoint_rows", "connectivity_only_endpoint_rows", "status",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def standardize_inchikey(smiles: str) -> str:
    if not smiles:
        return ""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        parent = rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(mol))
        canonical = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
        prepared = Chem.MolFromSmiles(canonical)
        return Chem.MolToInchiKey(prepared) if prepared is not None else ""
    except Exception:
        return ""


def parse_nm(raw: str) -> tuple[str, str, str]:
    match = NUMERIC.fullmatch(raw.strip())
    if match is None:
        return "", "", "INVALID_RAW_ENDPOINT"
    relation = {"≤": "<=", "≥": ">="}.get(match.group(1), match.group(1) or "=")
    number = float(match.group(2) + ("e" + match.group(3) if match.group(3) else ""))
    if not math.isfinite(number) or number <= 0:
        return relation, "", "INVALID_RAW_ENDPOINT"
    return relation, format(number, ".12g"), "PARSED_RAW_ENDPOINT"


def read_neighbors(path: Path, require_ten: bool = False):
    required = {"own_compound_id", "target_gene", "ligand_id", "ligand_inchikey", "rank_in_target"}
    result = {}
    with path.open("r", encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError("neighbors_audit.tsv is missing required columns")
        for row in reader:
            key = (row["own_compound_id"], row["target_gene"], row["ligand_id"])
            if key[:2] not in FOCUS:
                continue
            if key in result:
                raise ValueError(f"Duplicate neighbor: {key}")
            inchikey = row["ligand_inchikey"].strip().upper()
            if not INCHIKEY.fullmatch(inchikey):
                raise ValueError(f"Invalid InChIKey for {key}")
            row["ligand_inchikey"] = inchikey
            result[key] = row
    if not result or len(result) > 10:
        raise ValueError(f"Expected between 1 and 10 neighbors, found {len(result)}")
    if require_ten:
        pair_counts = Counter((k[0], k[1]) for k in result)
        if len(result) != 10 or any(pair_counts[p] != 5 for p in FOCUS):
            raise ValueError("Expected exactly 5 neighbors for each of PYZ-035/CDK4 and PYZ-055/BCL2L1")
    return result


@contextmanager
def tsv_stream(path: Path):
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            candidates = [item for item in archive.infolist()
                          if not item.is_dir() and item.filename.lower().endswith(".tsv")]
            if len(candidates) != 1:
                raise ValueError("BindingDB ZIP must contain exactly one primary TSV")
            with archive.open(candidates[0]) as source:
                with io.TextIOWrapper(source, encoding="utf-8-sig", newline="") as text:
                    yield csv.reader(text, delimiter="\t")
    else:
        with path.open("r", encoding="utf-8-sig", newline="") as source:
            yield csv.reader(source, delimiter="\t")


def unique_column(header, title, required=True):
    matches = [i for i, name in enumerate(header) if name.strip() == title]
    if len(matches) > 1 or (required and not matches):
        raise ValueError(f"Ambiguous or missing BindingDB column {title!r}: {len(matches)}")
    return matches[0] if matches else None


def scan(neighbors, db_path: Path):
    try:
        csv.field_size_limit(sys.maxsize)
    except OverflowError:
        csv.field_size_limit(2**31 - 1)
    by_target = defaultdict(list)
    for neighbor in neighbors.values():
        by_target[neighbor["target_gene"]].append(neighbor)
    matches = []
    counts, coverage_counts = Counter(), Counter()
    with tsv_stream(db_path) as records:
        header = next(records, None)
        if not header:
            raise ValueError("Empty BindingDB TSV")
        # BindingDB repeats the same chain column for multichain complexes.
        chain_cols = [i for i, h in enumerate(header)
                      if h.strip().startswith("UniProt (SwissProt) Primary ID of Target Chain")]
        if not chain_cols:
            raise ValueError("Missing UniProt target-chain columns")
        required = ("BindingDB Reactant_set_id", "Ligand SMILES", "Ligand InChI Key",
                    "BindingDB MonomerID", "Target Name", "Article DOI")
        columns = {k: unique_column(header, k) for k in required}
        columns.update({k: unique_column(header, v) for k, v in ENDPOINTS.items()})
        optional = {
            "organism": "Target Source Organism According to Curator or DataSource",
            "pmid": "PMID",
            "curation": "Curation/DataSource",
            "n_chains": "Number of Protein Chains in Target (>1 implies a multichain complex)",
        }
        columns.update({k: unique_column(header, v, required=False) for k, v in optional.items()})
        for number, fields in enumerate(records, start=2):
            counts["rows_scanned"] += 1
            if len(fields) > len(header):
                counts["rows_schema_mismatch"] += 1
                continue
            if len(fields) < len(header):
                fields.extend([""] * (len(header) - len(fields)))
            accessions = sorted({fields[i].strip().upper() for i in chain_cols if fields[i].strip()})
            genes = [gene for gene, uniprot in TARGETS.items() if uniprot in accessions]
            if not genes:
                continue
            counts["rows_relevant_uniprot"] += 1
            smile = fields[columns["Ligand SMILES"]].strip()
            computed = standardize_inchikey(smile)
            if not computed:
                counts["rows_relevant_invalid_smiles"] += 1
                continue
            original = fields[columns["Ligand InChI Key"]].strip().upper()
            if original and original != computed:
                counts["rows_source_vs_standardized_key_differ"] += 1
            raw_nchains = fields[columns["n_chains"]] if columns["n_chains"] is not None else ""
            n_chains = max(len(accessions), int(raw_nchains) if raw_nchains.isdigit() else 0)
            for gene in genes:
                for ref in by_target[gene]:
                    expected = ref["ligand_inchikey"]
                    if computed == expected:
                        relation = "EXACT_STANDARDIZED_INCHIKEY"
                    elif computed.split("-")[0] == expected.split("-")[0]:
                        relation = "CONNECTIVITY_ONLY_STEREO_UNRESOLVED"
                    else:
                        continue
                    counts["rows_with_candidate_structure"] += 1
                    for endpoint in ENDPOINTS:
                        raw = fields[columns[endpoint]].strip()
                        if not raw:
                            continue
                        qualifier, numeric, quality = parse_nm(raw)
                        if quality == "INVALID_RAW_ENDPOINT":
                            counts["unparseable_endpoint_candidates"] += 1
                        else:
                            coverage_counts[(ref["own_compound_id"], gene, ref["ligand_id"], relation)] += 1
                        matches.append({
                            "own_compound_id": ref["own_compound_id"],
                            "target_gene": gene, "ligand_id": ref["ligand_id"],
                            "rank_in_target": ref["rank_in_target"], "tanimoto": ref.get("tanimoto", ""),
                            "expected_inchikey": expected, "bindingdb_inchikey": original,
                            "computed_inchikey": computed, "structure_match": relation,
                            "target_uniprot": TARGETS[gene],
                            "all_target_chain_uniprot": ";".join(accessions),
                            "target_complexity": "MULTI_CHAIN_REVIEW_REQUIRED" if n_chains > 1 else "SINGLE_CHAIN_OR_UNRESOLVED",
                            "target_name": fields[columns["Target Name"]],
                            "organism": fields[columns["organism"]] if columns["organism"] is not None else "",
                            "bindingdb_reactant_set_id": fields[columns["BindingDB Reactant_set_id"]],
                            "bindingdb_monomer_id": fields[columns["BindingDB MonomerID"]],
                            "endpoint": endpoint, "raw_nm": raw, "relation": qualifier,
                            "value_nm": numeric, "endpoint_quality": quality,
                            "article_doi": fields[columns["Article DOI"]],
                            "pmid": fields[columns["pmid"]] if columns["pmid"] is not None else "",
                            "curation_source": fields[columns["curation"]] if columns["curation"] is not None else "",
                            "source_row_number": number, "reference_compound_smiles": smile,
                            "source_inchikey_vs_computed": "SAME" if original == computed else "DIFFERENT_OR_MISSING_REVIEW_REQUIRED",
                            "measurement_owner": "known_reference_ligand_not_pyz",
                            "scientific_approval": "NOT_APPROVED_MANUAL_REVIEW_REQUIRED",
                        })
                        counts["matched_endpoint_rows"] += 1
    coverage = []
    for key, ref in neighbors.items():
        exact = coverage_counts[(*key, "EXACT_STANDARDIZED_INCHIKEY")]
        stereo = coverage_counts[(*key, "CONNECTIVITY_ONLY_STEREO_UNRESOLVED")]
        coverage.append({
            **{k: ref.get(k, "") for k in ("own_compound_id", "target_gene", "ligand_id",
                                           "rank_in_target", "tanimoto", "ligand_inchikey")},
            "exact_endpoint_rows": exact, "connectivity_only_endpoint_rows": stereo,
            "status": "UNREVIEWED_EXACT_MATCH" if exact else
                      ("STEREO_REVIEW_REQUIRED" if stereo else "NO_VALID_MATCH_IN_THIS_RELEASE"),
        })
    return coverage, matches, dict(sorted(counts.items()))


def write_tsv(path: Path, columns, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--neighbors", required=True, type=Path)
    parser.add_argument("--bindingdb", required=True, type=Path)
    parser.add_argument("--expected-sha256", default="")
    parser.add_argument("--require-ten", action="store_true")
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    for path in (args.neighbors, args.bindingdb):
        if not path.is_file():
            raise FileNotFoundError(path)
    out = args.out_dir.resolve()
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    for path in (args.neighbors, args.bindingdb):
        if out.is_relative_to(path.resolve().parent):
            raise ValueError("Output cannot be inside scientific source directory")
    db_sha = sha256(args.bindingdb)
    if args.expected_sha256 and db_sha.lower() != args.expected_sha256.strip().lower():
        raise ValueError("BindingDB SHA256 mismatch")
    neighbors = read_neighbors(args.neighbors, require_ten=args.require_ten)
    coverage, matched, stats = scan(neighbors, args.bindingdb)
    report = {
        "contract": "mcl-pyz-e1-p0b-full-bindingdb-readonly-v2",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256": {str(args.neighbors.resolve()): sha256(args.neighbors),
                         str(args.bindingdb.resolve()): db_sha},
        "source_release_file": args.bindingdb.name, "targets": TARGETS,
        "neighbors": len(neighbors), "stats": stats,
        "neighbors_exact_endpoint_coverage": sum(int(r["exact_endpoint_rows"]) > 0 for r in coverage),
        "neighbors_connectivity_only_coverage": sum(int(r["connectivity_only_endpoint_rows"]) > 0 for r in coverage),
        "primary_assays_approved": 0,
        "scientific_status": "NO_PRIMARY_ASSAYS_APPROVED",
        "limitations": [
            "Target-chain UniProt match does not establish assay specificity or a single-protein target.",
            "BindingDB and imported ChEMBL records may duplicate underlying experiments.",
            "No experimental activity of PYZ is established by reference-ligand measurements.",
            "A partial InChIKey match is not proof of stereochemical identity.",
            "A missing result is limited to this release, not all scientific sources.",
            "A parsed endpoint is still not a manually approved primary assay.",
        ],
    }
    out.mkdir(parents=True)
    write_tsv(out / "bindingdb_full_matches.tsv", MATCH_FIELDS, matched)
    write_tsv(out / "bindingdb_full_neighbor_coverage.tsv", COVER_FIELDS, coverage)
    (out / "bindingdb_full_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"rows_scanned": stats.get("rows_scanned", 0),
                      "endpoint_rows": len(matched),
                      "exact_neighbors": report["neighbors_exact_endpoint_coverage"],
                      "primary_assays_approved": 0}, ensure_ascii=False))
    return report


if __name__ == "__main__":
    main()
