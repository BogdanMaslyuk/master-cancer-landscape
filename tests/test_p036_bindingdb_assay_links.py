"""P0.3.6 BindingDB real-schema linkage tests with wholly synthetic records."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import io
import json
import zipfile
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "link_bindingdb_assays_p036.py"
SPEC = importlib.util.spec_from_file_location("link_bindingdb_assays_p036", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def zip_tsv(path: Path, inner: str, headers: list[str], rows: list[list[str]]) -> None:
    text = io.StringIO(newline="")
    writer = csv.writer(text, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE)
    writer.writerow(headers)
    writer.writerows(rows)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(inner, text.getvalue())


def make_input(tmp_path):
    source = tmp_path / "sources"
    source.mkdir()
    panel_path = tmp_path / "panel" / "pilot_panel.tsv"
    panel_path.parent.mkdir()
    panel = pd.DataFrame([{
        "scope": "core" if i < 13 else "expanded",
        "target_gene": "G" + str(i),
        "uniprot_accession": "Q" + str(i).zfill(5),
    } for i in range(20)])
    panel.to_csv(panel_path, sep="\t", index=False)
    articles = source / "articles.zip"
    mapping = source / "mapping.zip"
    assays = source / "assays.zip"
    zip_tsv(articles, "BindingDB_BindingDB_Articles.tsv", [
        "BindingDB Reactant_set_id", "BindingDB MonomerID", "Ligand SMILES",
        "UniProt (SwissProt) Primary ID of Target Chain 1", "Ki (nM)",
        "Kd (nM)", "IC50 (nM)", "Article DOI",
        "Number of Protein Chains in Target (>1 implies a multichain complex)",
        "UniProt (SwissProt) Primary ID of Target Chain 2",
    ], [
        ["100", "L1", "CCO", "Q00000", "10", "", "", "10.1000/a", "1", ""],
        ["101", "L2", "CCN", "Q00000", "", "", "20", "10.1000/a", "1", ""],
        ["102", "L3", "CCC", "Q00013", "", "20", "", "10.1000/b", "2", "Q00000"],
        ["103", "L4", "CCCl", "Q00013", "", "", "30", "10.1000/b", "1", ""],
        ["104", "L5", "CCBr", "Q00013", "50", "", "", "10.1000/b", "1", ""],
        ["105", "L6", "CCF", "Q00000", "", "", "<100", "10.1000/a", "1", ""],
        ["106", "L7", "CCO", "Q00000", "", "5", "", "10.1000/a", "1", ""],
        ["107", "L8", "CCO", "Q00013", "100", "", "", "10.1000/b", "1", ""],
        ["108", "L9", "CCO", "Q00000", "20", "", "", "10.1000/a", "1", ""],
        ["108", "L10", "CCO", "Q00000", "25", "", "", "10.1000/a", "1", ""],
        ["109", "L11", "CCO", "Q00013", "200", "", "", "10.1000/b", "1", ""],
        ["110", "L12", "CCO", "Q00000", "", "", "", "10.1000/a", "1", ""],
        ["111", "L13", "CCO", "P99999", "100", "", "", "10.1000/c", "1", ""],
    ])
    zip_tsv(mapping, "BindingDB_rsid_eaids.tsv", [
        "REACTANT_SET_ID", "ENTRYID_ASSAYID",
    ], [
        ["100", "285_1"], ["101", "285_2"], ["104", "285_1"],
        ["104", "285_2"], ["105", "900_1"], ["106", "777_1"],
        ["107", "BAD"], ["108", "285_1"],
    ])
    zip_tsv(assays, "BindingDB_Assays.tsv", [
        "ENTRYID", "ASSAYID", "ASSAY_NAME", "DESCRIPTION",
    ], [
        ["285", "1", "Binding assay", "Single protein competition assay"],
        ["285", "2", "Functional assay", ""],
        ["777", "1", "Assay variant", "Method A"],
        ["777", "1", "Assay variant", "Method B"],
    ])
    return {
        "articles": articles, "mapping": mapping, "assays": assays,
        "panel": panel_path,
        "checksums": {
            k: MODULE.checksum(p) for k,p in
            (("articles",articles),("mapping",mapping),("assays",assays))
        },
    }


def test_actual_header_join_and_fail_closed_review_status(tmp_path):
    inp = make_input(tmp_path)
    report, rows, coverage = MODULE.link(
        inp["articles"], inp["mapping"], inp["assays"],
        inp["panel"], inp["checksums"])
    indexed = {r["reactant_set_id"]: r for r in rows}
    assert len(rows) == 10
    assert report["article_counts"]["excluded_multichain_rows"] == 1
    assert indexed["100"]["entryid_assayid"] == "285_1"
    assert indexed["100"]["assay_description"] == "Single protein competition assay"
    assert indexed["100"]["link_status"] == "LINKED_DESCRIPTION_PRESENT_UNREVIEWED"
    assert indexed["101"]["link_status"] == "ASSAY_DESCRIPTION_MISSING"
    assert indexed["103"]["link_status"] == "NO_ASSAY_MAPPING"
    assert indexed["104"]["link_status"] == "MULTIPLE_ASSAY_MAPPINGS"
    assert indexed["105"]["link_status"] == "NO_ASSAY_METADATA"
    assert indexed["106"]["link_status"] == "CONFLICTING_ASSAY_METADATA"
    assert indexed["107"]["link_status"] == "INVALID_ASSAY_MAPPING"
    assert indexed["108"]["link_status"] == "DUPLICATE_REACTANT_SET_ID_IN_ARTICLES"
    assert all(r["review_decision"] == "" for r in rows)
    assert all(r["reviewer"] == "" for r in rows)
    assert all(r["approved_measurements"] == 0 for r in coverage)
    assert report["scientific_status"] == "NO_ASSAYS_APPROVED_PENDING_MANUAL_REVIEW"


def test_archive_checksum_is_strict(tmp_path):
    inp = make_input(tmp_path)
    inp["checksums"]["mapping"] = "0"*64
    with pytest.raises(ValueError, match="SHA-256"):
        MODULE.link(inp["articles"], inp["mapping"], inp["assays"],
                    inp["panel"], inp["checksums"])


def test_cli_writes_immutable_outputs_outside_source_dirs(tmp_path, monkeypatch):
    inp = make_input(tmp_path)
    out = tmp_path / "analysis"
    argv = [
        "link_bindingdb_assays_p036.py", "--articles-zip", str(inp["articles"]),
        "--mapping-zip", str(inp["mapping"]), "--assays-zip", str(inp["assays"]),
        "--pilot-panel", str(inp["panel"]), "--out-dir", str(out),
    ]
    for kind in ("articles", "mapping", "assays"):
        argv += ["--" + kind + "-sha256", inp["checksums"][kind]]
    monkeypatch.setattr("sys.argv", argv)
    MODULE.main()
    assert (out / "assay_review_queue.tsv").is_file()
    assert (out / "assay_link_coverage.tsv").is_file()
    assert json.loads((out / "assay_link_report.json").read_text())["scientific_status"] == "NO_ASSAYS_APPROVED_PENDING_MANUAL_REVIEW"
    with pytest.raises(FileExistsError, match="Immutable"):
        MODULE.main()
    assert inp["articles"].is_file()
