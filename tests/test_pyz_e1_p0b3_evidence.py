"""Self-contained synthetic tests for the E1-P0-B3 evidence audit.

No real science dataset, network, services or scientific-source writes needed.
"""
from __future__ import annotations

import csv
import io
import json
import tempfile
import unittest
import zipfile
from pathlib import Path

from scripts.audit_pyz_e1_p0b3_evidence import main, analyze, sha256


def tsv(rows, fields):
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def synthetic_zip(path, transform=None):
    fullhash, nhash = "a" * 64, "b" * 64
    coverage, matches = [], []
    fields_coverage = ["own_compound_id", "target_gene", "ligand_id", "rank_in_target", "exact_endpoint_rows"]
    fields_matches = [
        "own_compound_id", "target_gene", "ligand_id", "rank_in_target", "tanimoto",
        "expected_inchikey", "bindingdb_inchikey", "computed_inchikey", "structure_match",
        "target_uniprot", "all_target_chain_uniprot", "target_complexity", "target_name",
        "bindingdb_reactant_set_id", "bindingdb_monomer_id", "endpoint", "raw_nm", "relation",
        "value_nm", "article_doi", "pmid", "curation_source", "source_row_number",
        "measurement_owner", "scientific_approval",
    ]
    for gene, pyz, acc, letter in [("CDK4", "PYZ-035", "P11802", "A"),
                                   ("BCL2L1", "PYZ-055", "Q07817", "B")]:
        for rank in range(1, 6):
            ligand = f"TLS-{letter}{rank}"
            n = 2 if (gene == "CDK4" and rank == 2) else 1
            coverage.append({"own_compound_id": pyz, "target_gene": gene,
                             "ligand_id": ligand, "rank_in_target": rank, "exact_endpoint_rows": n})
            for duplicate_index in range(n):
                number = len(matches) + 1
                key = "ABCDEFGHIJKLMN-UHFFFAOYSA-N"
                matches.append({
                    "own_compound_id": pyz, "target_gene": gene, "ligand_id": ligand,
                    "rank_in_target": rank, "tanimoto": "0.4" if rank == 1 else "0.2",
                    "expected_inchikey": key, "bindingdb_inchikey": key if number <= 3 else "ABCDEFGHIJKLMN-AAAAAAAAAA-N",
                    "computed_inchikey": key, "structure_match": "EXACT_STANDARDIZED_INCHIKEY",
                    "target_uniprot": acc, "all_target_chain_uniprot": acc + (";P24385" if gene == "CDK4" else ""),
                    "target_complexity": "MULTI_CHAIN_REVIEW_REQUIRED" if gene == "CDK4" else "SINGLE_CHAIN_OR_UNRESOLVED",
                    "target_name": gene, "bindingdb_reactant_set_id": str(1000 + number),
                    "bindingdb_monomer_id": f"M-{gene}-{rank}",
                    "endpoint": "IC50" if gene == "CDK4" else "Ki",
                    "raw_nm": str(number * 100), "relation": "=", "value_nm": str(number * 100),
                    "article_doi": "10.1111/source" if gene == "CDK4" else "10.2222/source",
                    "pmid": "123456" if gene == "CDK4" else "987654",
                    "curation_source": "ChEMBL", "source_row_number": str(2000 + number),
                    "measurement_owner": "known_reference_ligand_not_pyz",
                    "scientific_approval": "NOT_APPROVED_MANUAL_REVIEW_REQUIRED",
                })
    report = {"contract": "mcl-pyz-e1-p0b-full-bindingdb-readonly-v2",
              "scientific_status": "NO_PRIMARY_ASSAYS_APPROVED", "primary_assays_approved": 0,
              "neighbors": 10, "stats": {"matched_endpoint_rows": 11},
              "neighbors_exact_endpoint_coverage": 10,
              "input_sha256": {"bindingdb.zip": fullhash, "neighbors.tsv": nhash},
              "source_release_file": "BindingDB_All_test_tsv.zip"}
    manifest = {"stage": "E1_P0B2", "scientific_status": "NO_PRIMARY_ASSAYS_APPROVED",
                "bindingdb_sha256": fullhash, "neighbors_sha256": nhash}
    if transform:
        transform(report, manifest, matches, coverage)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("bindingdb_full_report.json", json.dumps(report))
        z.writestr("run_manifest.json", json.dumps(manifest))
        z.writestr("bindingdb_full_matches.tsv", tsv(matches, fields_matches))
        z.writestr("bindingdb_full_neighbor_coverage.tsv", tsv(coverage, fields_coverage))


class EvidenceAuditTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / "result.zip"
        self.out = self.root / "new_review"

    def test_realistic_synthetic_register_and_source_immutable(self):
        synthetic_zip(self.source)
        digest = sha256(self.source)
        result = main(["--source-zip", str(self.source), "--out-dir", str(self.out)])
        self.assertEqual(sha256(self.source), digest)
        self.assertEqual(result["statistics"]["records"], 11)
        self.assertEqual(result["statistics"]["source_inchikey_discrepancy"], 8)
        self.assertEqual(result["statistics"]["potential_duplicate_rows"], 2)
        self.assertEqual(result["statistics"]["reference_ligands_above_035"], 2)
        self.assertEqual(result["statistics"]["primary_assays_approved"], 0)
        self.assertEqual(result["scientific_status"], "NO_PRIMARY_ASSAYS_APPROVED")
        with (self.out / "assay_evidence_register.tsv").open(encoding="utf-8-sig") as f:
            register = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(len(register), 11)
        self.assertEqual(len({x["record_id"] for x in register}), 11)
        self.assertTrue(all(x["scientific_approval"] == "NO_PRIMARY_ASSAYS_APPROVED" for x in register))
        with (self.out / "manual_source_review_checklist.tsv").open(encoding="utf-8-sig") as f:
            checklist = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(len(checklist), 11)
        self.assertTrue(all(x["decision"] == "NOT_REVIEWED" and x["reviewer"] == "" for x in checklist))
        with self.assertRaises(FileExistsError):
            main(["--source-zip", str(self.source), "--out-dir", str(self.out)])

    def test_reject_changed_source_hash(self):
        synthetic_zip(self.source, lambda report, manifest, rows, coverage: manifest.update(bindingdb_sha256="c" * 64))
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            analyze(self.source)
        self.assertFalse(self.out.exists())

    def test_reject_pyz_owned_measurement(self):
        synthetic_zip(self.source, lambda report, manifest, rows, coverage: rows[0].update(measurement_owner="measured_pyz"))
        with self.assertRaisesRegex(ValueError, "PYZ experimental ownership"):
            analyze(self.source)

    def test_reject_wrong_crosswalk_coverage(self):
        synthetic_zip(self.source, lambda report, manifest, rows, coverage: coverage[0].update(exact_endpoint_rows=3))
        with self.assertRaisesRegex(ValueError, "Coverage endpoint count mismatch"):
            analyze(self.source)

    def test_reject_unapproved_to_approved_status(self):
        synthetic_zip(self.source, lambda report, manifest, rows, coverage: rows[0].update(scientific_approval="APPROVED"))
        with self.assertRaisesRegex(ValueError, "source scientific approval"):
            analyze(self.source)

    def test_reject_target_misassignment(self):
        synthetic_zip(self.source, lambda report, manifest, rows, coverage: rows[0].update(target_uniprot="P50750"))
        with self.assertRaisesRegex(ValueError, "PYZ/target/UniProt"):
            analyze(self.source)


if __name__ == "__main__":
    unittest.main()