"""P0.3.6 BindingDB ZIP subset audit: synthetic data only, NOT experimental evidence."""
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

SCRIPT=Path(__file__).resolve().parents[1]/"scripts"/"screen_bindingdb_articles_p036.py"
SPEC=importlib.util.spec_from_file_location("screen_bindingdb_articles_p036",SCRIPT)
MODULE=importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def make_zip(path: Path, inner: str, columns: list[str], rows: list[list[str]]) -> None:
    stream=io.StringIO(newline="")
    writer=csv.writer(stream,delimiter="\t",lineterminator="\n",quoting=csv.QUOTE_NONE)
    writer.writerow(columns)
    writer.writerows(rows)
    with zipfile.ZipFile(path,"w",compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(inner,stream.getvalue())


def inputs(tmp_path: Path) -> dict:
    source=tmp_path/"source"
    source.mkdir()
    panel=tmp_path/"pilot"/"pilot_panel.tsv"
    panel.parent.mkdir()
    core=["GENE"+str(i).zfill(2) for i in range(13)]
    extra=["AKT1","AURKA","BRAF","CDK6","MET","MTOR","PARP1"]
    pd.DataFrame([
        {"target_gene":gene,"uniprot_accession":"Q"+str(i+1).zfill(5),
         "scope":"core" if i<13 else "expanded"}
        for i,gene in enumerate(core+extra)
    ]).to_csv(panel,sep="\t",index=False)
    cols=[
        "BindingDB Reactant_set_id","BindingDB MonomerID","Ligand SMILES",
        "UniProt (SwissProt) Primary ID of Target Chain",
        "Ki (nM)","Kd (nM)","IC50 (nM)","Article DOI",
        "Number of Protein Chains in Target (>1 implies a multichain complex)",
    ]
    records=[
        ["100","M1","CCO","Q00001","10","","","10.1234/a","1"],
        ["101","M1","CCO","Q00001","","","<200","10.1234/a","1"],
        ["102","M2","CCN","Q00014","","20","","10.1234/b","2"],
        ["103","M3","CCC","P99999","500","","","10.1234/c","1"],
        ["104","M4","CCCl","Q00014","","","","","1"],
    ]
    articles=source/"BindingDB_BindingDB_Articles_202611_tsv.zip"
    assays=source/"BindingDB_Assays_202611_tsv.zip"
    mapping=source/"BindingDB_rsid_eaids_202611_tsv.zip"
    make_zip(articles,"BindingDB_BindingDB_Articles.tsv",cols,records)
    make_zip(assays,"BindingDB_Assays.tsv",
             ["Entry ID","Assay ID","Assay Description"],[["E1","A1","Test only"]])
    make_zip(mapping,"BindingDB_rsid_eaids.tsv",
             ["Reactant_set_id","EntryID_AssayID"],[["100","E1_A1"]])
    return {"articles":articles,"assays":assays,"mapping":mapping,"panel":panel}


def test_articles_only_20_target_observation_not_a_b_approval(tmp_path):
    inp=inputs(tmp_path)
    report,rows=MODULE.audit(inp["articles"],inp["assays"],inp["mapping"],inp["panel"])
    assert len(rows)==20
    assert report["counts"]["total_article_rows_scanned"]==5
    assert report["counts"]["pilot_articles_rows"]==4
    assert report["proteins_with_any_articles_endpoint_rows"]==2
    assert report["scientific_status"]=="EXPLORATORY_NOT_A_B_APPROVED"
    first=rows[0]
    assert first["target_gene"]=="GENE00"
    assert first["raw_Ki_rows"]==1
    assert first["raw_IC50_rows"]==1
    assert first["distinct_BindingDB_monomer_ids"]==1
    second=next(row for row in rows if row["target_gene"]=="AKT1")
    assert second["raw_Kd_rows"]==1
    assert second["multichain_rows"]==1
    assert second["rows_with_endpoint"]==1
    assert not any(row["experimental_quality_status"] != "UNREVIEWED_SUBSET" for row in rows)


def test_checksum_mismatch_rejected_before_scientific_scan(tmp_path):
    inp=inputs(tmp_path)
    with pytest.raises(ValueError,match="checksum mismatch"):
        MODULE.audit(inp["articles"],inp["assays"],inp["mapping"],inp["panel"],
                     expected_hashes={"articles":"0"*64})


def test_report_writes_outside_sources_and_cannot_overwrite(tmp_path,monkeypatch):
    inp=inputs(tmp_path)
    out=tmp_path/"report"
    cmd=[
        "screen_bindingdb_articles_p036.py",
        "--articles-zip",str(inp["articles"]),
        "--assays-zip",str(inp["assays"]),
        "--mapping-zip",str(inp["mapping"]),
        "--pilot-panel",str(inp["panel"]),
        "--out-dir",str(out),
    ]
    monkeypatch.setattr("sys.argv",cmd)
    MODULE.main()
    assert (out/"pilot_bindingdb_articles_coverage.tsv").is_file()
    assert (out/"bindingdb_articles_schema_report.json").is_file()
    assert json.loads((out/"bindingdb_articles_schema_report.json").read_text())["counts"]["pilot_articles_rows"]==4
    with pytest.raises(FileExistsError,match="Immutable exploratory"):
        MODULE.main()
    assert inp["articles"].is_file()
