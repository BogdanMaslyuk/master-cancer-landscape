"""Pinned ChEMBL SQLite extraction tests with a miniature relational database."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import sqlite3
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1]/"scripts"/"export_tls_v2_pinned_sources.py"
SPEC = importlib.util.spec_from_file_location("export_tls_v2_pinned_sources",SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def source(tmp_path):
    db = tmp_path/"chembl_37_test.db"
    con = sqlite3.connect(db)
    con.executescript("""
    CREATE TABLE version (release INTEGER);
    INSERT INTO version VALUES (37);
    CREATE TABLE activities (
      activity_id INTEGER,assay_id INTEGER,molregno INTEGER,standard_type TEXT,
      standard_value TEXT,standard_units TEXT,standard_relation TEXT,
      data_validity_comment TEXT,potential_duplicate INTEGER
    );
    CREATE TABLE assays (
      assay_id INTEGER,chembl_id TEXT,tid INTEGER,doc_id INTEGER,
      confidence_score INTEGER,relationship_type TEXT,assay_type TEXT,
      bao_format TEXT,variant_id INTEGER,description TEXT
    );
    CREATE TABLE target_dictionary(tid INTEGER,chembl_id TEXT,tax_id INTEGER,target_type TEXT);
    CREATE TABLE target_components(tid INTEGER,component_id INTEGER);
    CREATE TABLE component_sequences(component_id INTEGER,accession TEXT);
    CREATE TABLE molecule_dictionary(molregno INTEGER,chembl_id TEXT);
    CREATE TABLE compound_structures(molregno INTEGER,canonical_smiles TEXT);
    CREATE TABLE bioassay_ontology(bao_id TEXT,label TEXT);
    CREATE TABLE docs(doc_id INTEGER,chembl_id TEXT,doi TEXT,pubmed_id TEXT);
    INSERT INTO target_dictionary VALUES(1,'CHEMBL-T1',9606,'SINGLE PROTEIN');
    INSERT INTO target_dictionary VALUES(2,'CHEMBL-T2',10090,'SINGLE PROTEIN');
    INSERT INTO target_components VALUES(1,10);
    INSERT INTO target_components VALUES(2,11);
    INSERT INTO component_sequences VALUES(10,'P11802');
    INSERT INTO component_sequences VALUES(11,'P99999');
    INSERT INTO molecule_dictionary VALUES(5,'CHEMBL-M1');
    INSERT INTO compound_structures VALUES(5,'CCO');
    INSERT INTO bioassay_ontology VALUES('BAO1','single protein format');
    INSERT INTO docs VALUES(9,'CHEMBL-D1','10.1234/verified','12345');
    INSERT INTO assays VALUES(10,'CHEMBL-A1',1,9,9,'D','B','BAO1',NULL,'Ki binding');
    INSERT INTO assays VALUES(11,'CHEMBL-A2',2,9,9,'D','B','BAO1',NULL,'Mouse assay');
    INSERT INTO assays VALUES(12,'CHEMBL-A3',1,9,8,'D','B','BAO1',NULL,'Lower-confidence');
    INSERT INTO assays VALUES(13,'CHEMBL-A4',1,9,9,'D','B','BAO1',77,'Variant');
    INSERT INTO activities VALUES(100,10,5,'Ki','100','nM','=','',0);
    INSERT INTO activities VALUES(101,11,5,'Ki','110','nM','=','',0);
    INSERT INTO activities VALUES(102,12,5,'Ki','120','nM','=','',0);
    INSERT INTO activities VALUES(103,13,5,'Ki','130','nM','=','',0);
    """)
    con.commit()
    con.close()
    registry = tmp_path/"protein_target_universe.tsv"
    registry.write_text(
        "target_gene\tuniprot_accession\tmapping_status\torganism_id\n"
        "CDK4\tP11802\tunique_reviewed_human\t9606\n"
    )
    return db,registry


def test_extract_filters_unmatched_species_homologues_and_variants(tmp_path):
    db,registry = source(tmp_path)
    h = hashlib.sha256(db.read_bytes()).hexdigest()
    out=tmp_path/"chembl.tsv"
    result=MODULE.export_sqlite(db,registry,out,release="ChEMBL_37",expected_db_sha256=h)
    assert result["rows"]==1
    with out.open(newline="") as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    assert len(rows)==1
    assert rows[0]["activity_id"]=="100"
    assert rows[0]["uniprot_accession"]=="P11802"
    assert rows[0]["doi"]=="10.1234/verified"
    assert rows[0]["bao_label"]=="single protein format"
    assert result["sql_sha256"]
    with pytest.raises(FileExistsError):
        MODULE.export_sqlite(db,registry,out,release="ChEMBL_37",expected_db_sha256=h)


def test_release_and_checksum_mismatch_fail_closed(tmp_path):
    db,registry=source(tmp_path)
    out=tmp_path/"chembl.tsv"
    with pytest.raises(ValueError,match="SHA-256 mismatch"):
        MODULE.export_sqlite(db,registry,out,release="ChEMBL_37",expected_db_sha256="0"*64)
    h=hashlib.sha256(db.read_bytes()).hexdigest()
    with pytest.raises(ValueError,match="not corroborated"):
        MODULE.export_sqlite(db,registry,out,release="ChEMBL_38",expected_db_sha256=h)
    assert not out.exists()
