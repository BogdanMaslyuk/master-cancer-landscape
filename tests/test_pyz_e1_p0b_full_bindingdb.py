"""Synthetic regression tests for the full BindingDB reader.

No large DB download, source mutation or assay approval is required.
"""
from __future__ import annotations

import csv
import zipfile

import pytest
from scripts.audit_pyz_e1_p0b_full_bindingdb import (
    main,
    parse_nm,
    standardize_inchikey,
)


def fixture_files(tmp_path):
    neighbors = tmp_path / "neighbors_audit.tsv"
    base = standardize_inchikey("CCO")
    stereo = standardize_inchikey("C[C@H](F)Cl")
    reverse = standardize_inchikey("C[C@@H](F)Cl")
    assert stereo[:14] == reverse[:14] and stereo != reverse
    with neighbors.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["own_compound_id", "target_gene", "ligand_id",
                         "ligand_inchikey", "rank_in_target", "tanimoto"])
        writer.writerows([
            ["PYZ-035", "CDK4", "TLS-E", base, "1", "0.37"],
            ["PYZ-055", "BCL2L1", "TLS-S", stereo, "1", "0.45"],
        ])
    header = [
        "BindingDB Reactant_set_id", "Ligand SMILES", "Ligand InChI Key",
        "BindingDB MonomerID", "Target Name", "Article DOI",
        "Ki (nM)", "Kd (nM)", "IC50 (nM)",
        "Number of Protein Chains in Target (>1 implies a multichain complex)",
        "UniProt (SwissProt) Primary ID of Target Chain",
        "UniProt (SwissProt) Primary ID of Target Chain",
    ]
    rows = [
        ["100", "CCO", base, "10", "CDK4/cyclin D1", "10.1111/example",
         "", "", "<450", "2", "P11802", "P24385"],
        ["101", "C[C@@H](F)Cl", reverse, "11", "BCL2L1", "10.2222/example",
         "80", "", "", "1", "Q07817", ""],
        ["102", "CCO", base, "12", "CDK9", "10.3333/example",
         "", "", "120", "1", "P50750", ""],
        ["103", "bad SMILES", "", "13", "CDK4", "10.4444/example",
         "", "", "100", "1", "P11802", ""],
    ]
    database = tmp_path / "BindingDB_All_202611_tsv.zip"
    with zipfile.ZipFile(database, "w") as archive:
        archive.writestr(
            "BindingDB_All.tsv",
            "\n".join("\t".join(row) for row in [header, *rows]) + "\n",
        )
    return neighbors, database


def test_exact_stereo_and_false_target(tmp_path):
    neighbors, database = fixture_files(tmp_path)
    output = tmp_path / "out"
    report = main(["--neighbors", str(neighbors), "--bindingdb", str(database),
                   "--out-dir", str(output)])
    assert report["neighbors"] == 2
    assert report["stats"]["rows_scanned"] == 4
    assert report["stats"]["rows_relevant_invalid_smiles"] == 1
    assert report["neighbors_exact_endpoint_coverage"] == 1
    assert report["neighbors_connectivity_only_coverage"] == 1
    assert report["primary_assays_approved"] == 0
    with (output / "bindingdb_full_matches.tsv").open(encoding="utf-8-sig", newline="") as source:
        matches = list(csv.DictReader(source, delimiter="\t"))
    assert len(matches) == 2
    assert matches[0]["raw_nm"] == "<450"
    assert matches[0]["target_complexity"] == "MULTI_CHAIN_REVIEW_REQUIRED"
    assert matches[0]["measurement_owner"] == "known_reference_ligand_not_pyz"
    assert matches[1]["structure_match"] == "CONNECTIVITY_ONLY_STEREO_UNRESOLVED"


def test_ten_neighbor_guard_and_no_overwrite(tmp_path):
    neighbors, database = fixture_files(tmp_path)
    output = tmp_path / "out"
    with pytest.raises(ValueError, match="exactly 5"):
        main(["--neighbors", str(neighbors), "--bindingdb", str(database),
              "--out-dir", str(output), "--require-ten"])
    assert not output.exists()
    main(["--neighbors", str(neighbors), "--bindingdb", str(database),
          "--out-dir", str(output)])
    with pytest.raises(FileExistsError):
        main(["--neighbors", str(neighbors), "--bindingdb", str(database),
              "--out-dir", str(output)])


def test_sha_and_schema_guards(tmp_path):
    neighbors, database = fixture_files(tmp_path)
    output = tmp_path / "out"
    with pytest.raises(ValueError, match="SHA256"):
        main(["--neighbors", str(neighbors), "--bindingdb", str(database),
              "--expected-sha256", "0" * 64, "--out-dir", str(output)])
    assert not output.exists()
    with zipfile.ZipFile(database, "w") as archive:
        archive.writestr("wrong.tsv", "BindingDB Reactant_set_id\tLigand SMILES\n1\tCCO\n")
    with pytest.raises(ValueError, match="UniProt"):
        main(["--neighbors", str(neighbors), "--bindingdb", str(database),
              "--out-dir", str(output)])
    assert not output.exists()


def test_invalid_endpoint_never_counts_as_coverage():
    assert parse_nm("<0")[2] == "INVALID_RAW_ENDPOINT"
    assert parse_nm("n/a")[2] == "INVALID_RAW_ENDPOINT"
    assert parse_nm(">100") == (">", "100", "PARSED_RAW_ENDPOINT")
