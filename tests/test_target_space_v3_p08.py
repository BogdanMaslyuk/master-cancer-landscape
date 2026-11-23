"""Synthetic offline tests; no actual user science data or internet required."""
from __future__ import annotations
import csv
from io import StringIO
import json
from pathlib import Path
import sqlite3
import sys
import zipfile
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_target_space_v3_p08 import (CONTRACT, P0_CONTRACT, P06_CONTRACT,FOCUS6,
    ReadOnlyZip, sha_bytes,shafile,parse_nm,rows_from_pilot,chembl_local,bindingdb_local,
    open_binding_tsv,run,single_header)


def zip_from(folder:Path,out:Path):
    with zipfile.ZipFile(out,'w') as z:
        for p in folder.iterdir():z.write(p,arcname=p.name)
    return out


@pytest.fixture
def input_archive(tmp_path):
    p0dir=tmp_path/'p0';p0dir.mkdir()
    (p0dir/'target_universe.tsv').write_text('gene_symbol\n'+ '\n'.join(FOCUS6)+'\n')
    (p0dir/'coverage_summary.json').write_text(json.dumps({'contract':P0_CONTRACT,'technical_status':'REVIEW_REQUIRED'}))
    (p0dir/'source_manifest.json').write_text('{}')
    p0zip=zip_from(p0dir,tmp_path/'P0.zip')
    sha={p.name:shafile(p) for p in p0dir.iterdir()}
    p06dir=tmp_path/'p06';p06dir.mkdir()
    (p06dir/'p06_summary.json').write_text(json.dumps({'contract':P06_CONTRACT,'sample_size':6,
        'uniprot_snapshot':{'status':'RETRIEVED_REVIEWED_HUMAN'},
        'source_p0_sha256_before_after':sha}))
    rows=[]
    for i,g in enumerate(FOCUS6):
        acc=['P78527','Q07817','P08069','P21860','P49356','P20711'][i]
        rows.append({'gene_symbol':g,'uniprot_accession':acc,
                     'hgnc_id':'HGNC:'+str(i+1), 'mapping_status':'EXACT_REVIEWED_WITH_HGNC',
                     'sample_stratum':'MAPPED_OTHER'})
    for file in ('p06_selection.tsv','p06_external_coverage.tsv'):
        with (p06dir/file).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
    p06zip=zip_from(p06dir,tmp_path/'P06.zip')
    p07dir=tmp_path/'p07';p07dir.mkdir()
    (p07dir/'p07_summary.json').write_text(json.dumps({'scientific_acceptance':'NOT_APPROVED',
        'source_archives_sha256':{p0zip.name:shafile(p0zip),p06zip.name:shafile(p06zip)}}))
    p07zip=zip_from(p07dir,tmp_path/'P07.zip')
    return p0zip,p06zip,p07zip


@pytest.fixture
def mock_sqlite(tmp_path):
    p=tmp_path/'chembl_36.db'
    db=sqlite3.connect(p)
    db.executescript('''
    CREATE TABLE target_dictionary(tid INTEGER PRIMARY KEY,chembl_id TEXT,target_type TEXT,organism TEXT);
    CREATE TABLE target_components(tid INT,component_id INT);
    CREATE TABLE component_sequences(component_id INT,accession TEXT);
    CREATE TABLE assays(assay_id INTEGER PRIMARY KEY,tid INT,assay_type TEXT,confidence_score INT,chembl_id TEXT);
    CREATE TABLE activities(activity_id INTEGER PRIMARY KEY,assay_id INT,molregno INT, standard_type TEXT,
      standard_units TEXT,standard_relation TEXT,standard_value FLOAT,potential_duplicate INT,
      data_validity_comment TEXT,doc_id INT);
    CREATE TABLE molecule_dictionary(molregno INT PRIMARY KEY,chembl_id TEXT);
    CREATE TABLE compound_structures(molregno INT PRIMARY KEY,canonical_smiles TEXT);
    CREATE TABLE docs(doc_id INTEGER PRIMARY KEY,chembl_id TEXT,doi TEXT);
    INSERT INTO target_dictionary VALUES (1,'CHEMBL100','SINGLE PROTEIN','Homo sapiens');
    INSERT INTO target_components VALUES (1,100);
    INSERT INTO component_sequences VALUES (100,'P78527');
    INSERT INTO assays VALUES (10,1,'B',9,'CHEMBLASSAY10');
    INSERT INTO molecule_dictionary VALUES (1,'CHEMBLMOL1');
    INSERT INTO compound_structures VALUES (1,'CCO');
    INSERT INTO docs VALUES (1,'CHEMBLDOC1','10.1000/sample');
    INSERT INTO activities VALUES (900,10,1,'Ki','nM','=',12,0,NULL,1);
    INSERT INTO activities VALUES (901,10,1,'Kd','nM','>',5000,0,NULL,1);
    INSERT INTO activities VALUES (902,10,1,'IC50','nM','=',50,1,NULL,1);
    INSERT INTO activities VALUES (903,10,1,'IC50','nM','=',51,0,'Potential author error',1);
    INSERT INTO activities VALUES (904,10,1,'IC50','uM','=',50,0,NULL,1);
    ''')
    db.commit();db.close()
    return p


@pytest.fixture
def mock_bindingdb(tmp_path):
    p=tmp_path/'BindingDB_All_202611.tsv'
    cols=['BindingDB Reactant_set_id','Ligand SMILES','Ligand InChI Key','BindingDB MonomerID','Target Name',
          'Ki (nM)','Kd (nM)','IC50 (nM)',
          'UniProt (SwissProt) Primary ID of Target Chain',
          'UniProt (SwissProt) Primary ID of Target Chain',
          'Number of Protein Chains in Target (>1 implies a multichain complex)',
          'Article DOI','PMID','Curation/DataSource','Target Source Organism According to Curator or DataSource']
    rows=[['10','CCO','KEY','MON1','Kinase','12','','','P78527','','1','10.1000/sample','111','article','Homo sapiens'],
          ['11','CCC','KEY2','MON2','Complex','','','20','P78527','P11802','2','','222','BindingDB','Homo sapiens'],
          ['12','CNC','KEY3','MON3','Kinase','<30','','','P78527','','1','','333','BindingDB','Homo sapiens'],
          ['13','CCN','KEY4','MON4','Kinase','','5','','P99999','','1','','444','BindingDB','Homo sapiens'],
          ['14','CO','KEY5','MON5','Other','50','','','P78527','','1','','555','BindingDB','Mus musculus']]
    with p.open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(cols);w.writerows(rows)
    return p


def test_nm_normal_and_censored():
    assert parse_nm('3e2')==('=','300','EXACT_NUMERIC')
    assert parse_nm('< 1.5')==('<','1.5','CENSORED')
    assert parse_nm('≥10')==('>=','10','CENSORED')


@pytest.mark.parametrize('v',['','0','nan','-1','abc','4 nm'])
def test_bad_binding_nm(v):
    assert parse_nm(v)[2]=='INVALID_OR_MISSING'


def test_focus6_p0_p07_provenance(input_archive):
    p0,p06,p07=input_archive
    panel,meta=rows_from_pilot(p06,p0,p07,'focus6')
    assert [x['gene_symbol'] for x in panel]==list(FOCUS6)
    assert meta['p0_sha256']==shafile(p0)


def test_sample40_expands_to_available_sample(input_archive):
    p0,p06,p07=input_archive
    panel,_=rows_from_pilot(p06,p0,p07,'sample40')
    assert len(panel)==6


def test_provenance_mismatch_rejected(input_archive,tmp_path):
    p0,p06,p07=input_archive
    incorrect=tmp_path/'broken_p0.zip'
    with zipfile.ZipFile(p0) as src,zipfile.ZipFile(incorrect,'w') as dst:
        for n in src.namelist(): dst.writestr(n,b'changed' if n=='target_universe.tsv' else src.read(n))
    with pytest.raises(ValueError,match='provenance mismatch'):
        rows_from_pilot(p06,incorrect,None,'focus6')


def test_readonly_chembl_matches_exact_protein(input_archive,mock_sqlite):
    _,p06,_=input_archive
    panel,_=rows_from_pilot(p06,None,None,'focus6')
    obs,cov,meta=chembl_local(mock_sqlite,panel,100)
    assert len(obs)==4
    assert next(x for x in obs if x['activity_id']==900)['evidence_status']=='ASSAY_METADATA_MATCHED_PRIMARY_REVIEW_REQUIRED'
    assert next(x for x in obs if x['activity_id']==900)['document_chembl_id']=='CHEMBLDOC1'
    assert next(x for x in obs if x['activity_id']==900)['doi']=='10.1000/sample'
    assert next(x for x in obs if x['activity_id']==901)['evidence_status']=='REVIEW_CENSORED_OR_RELATION'
    assert next(x for x in obs if x['activity_id']==902)['evidence_status']=='REVIEW_CURATOR_FLAGGED_OR_DUPLICATE'
    assert next(x for x in obs if x['activity_id']==903)['evidence_status']=='REVIEW_CURATOR_FLAGGED_OR_DUPLICATE'
    assert cov[0]['status']=='LOCAL_ASSAYS_FOUND'
    assert meta['confirmed_pyZ_binding'] is False


def test_chembl_limit_reports_truncated(input_archive,mock_sqlite):
    panel,_=rows_from_pilot(input_archive[1],None,None,'focus6')
    obs,cov,_=chembl_local(mock_sqlite,panel,1)
    assert len(obs)==1 and cov[0]['truncated'] and cov[0]['candidate_rows']==4


def test_chembl_reject_missing_schemas(input_archive,tmp_path):
    p=tmp_path/'invalid.db';sqlite3.connect(p).close()
    panel,_=rows_from_pilot(input_archive[1],None,None,'focus6')
    with pytest.raises(ValueError,match='schema invalid'):
        chembl_local(p,panel,100)


def test_sqlite_source_no_modify(input_archive,mock_sqlite):
    before=shafile(mock_sqlite)
    panel,_=rows_from_pilot(input_archive[1],None,None,'focus6')
    chembl_local(mock_sqlite,panel,100)
    assert shafile(mock_sqlite)==before


def test_bindingdb_duplicate_chain_columns_and_exact_uniprot(input_archive,mock_bindingdb):
    panel,_=rows_from_pilot(input_archive[1],None,None,'focus6')
    rows,cov,meta=bindingdb_local(mock_bindingdb,panel,100)
    assert meta['rows_scanned']==5
    assert len(rows)==4
    assert rows[0]['gene_symbol']=='PRKDC'
    assert rows[0]['identity_status']=='SWISSPROT_PRIMARY_SINGLE_CHAIN_MATCH_ASSAY_REVIEW'
    assert rows[1]['identity_status']=='REVIEW_MULTI_CHAIN'
    assert rows[2]['measurement_status']=='CENSORED_NOT_EXACT'
    assert rows[3]['identity_status']=='REVIEW_SPECIES_MISMATCH'
    assert cov[0]['chain_match_rows']==4


def test_bindingdb_zip_stream(input_archive,mock_bindingdb,tmp_path):
    z=tmp_path/'db.zip'
    with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as f:f.write(mock_bindingdb,arcname=mock_bindingdb.name)
    panel,_=rows_from_pilot(input_archive[1],None,None,'focus6')
    with open_binding_tsv(z) as r: assert next(r)[0]=='BindingDB Reactant_set_id'
    rows,cov,_=bindingdb_local(z,panel,100)
    assert len(rows)==4


def test_bindingdb_reject_no_primary_chain(input_archive,tmp_path):
    p=tmp_path/'bad.tsv';p.write_text('Target Name\tKi (nM)\nKinase\t10\n')
    panel,_=rows_from_pilot(input_archive[1],None,None,'focus6')
    with pytest.raises(ValueError,match='SwissProt'):
        bindingdb_local(p,panel,100)


def test_bindingdb_limit_reports_overlimit_not_invalid(input_archive,mock_bindingdb):
    panel,_=rows_from_pilot(input_archive[1],None,None,'focus6')
    rows,cov,_=bindingdb_local(mock_bindingdb,panel,1)
    assert len(rows)==1 and cov[0]['truncated']
    assert cov[0]['rows_over_limit']==3
    assert cov[0]['rows_invalid']==0


def test_readonly_full_run_reproducible(input_archive,mock_sqlite,mock_bindingdb,tmp_path):
    p0,p06,p07=input_archive
    prior={p:shafile(p) for p in [p0,p06,p07,mock_sqlite,mock_bindingdb]}
    one=run(p06_zip=p06,p0_zip=p0,p07_zip=p07,chembl_sqlite=mock_sqlite,
            bindingdb_tsv=mock_bindingdb,out=tmp_path/'audit1',release_chembl='chembl_36',release_bindingdb='202611')
    two=run(p06_zip=p06,p0_zip=p0,p07_zip=p07,chembl_sqlite=mock_sqlite,
            bindingdb_tsv=mock_bindingdb,out=tmp_path/'audit2',release_chembl='chembl_36',release_bindingdb='202611')
    assert one==two
    assert one['technical_status']=='REVIEW_REQUIRED'
    assert one['chembl_assay_metadata_eligible_n']==1
    assert one['bindingdb_exact_single_chain_rows_n']==2 # exact and censored, not PYZ
    assert one['scientific_verified_binding_to_pyz_n']==0
    assert one['external_requests_n']==0
    assert prior=={p:shafile(p) for p in prior}
    assert (tmp_path/'audit1'/'p08_chembl_observations.tsv').read_bytes()==(tmp_path/'audit2'/'p08_chembl_observations.tsv').read_bytes()


def test_no_local_sources_are_not_assessed(input_archive,tmp_path):
    p0,p06,p07=input_archive
    r=run(p06_zip=p06,p0_zip=p0,p07_zip=p07,out=tmp_path/'preflight')
    assert r['technical_status']=='INPUT_REQUIRED'
    assert r['chembl_rows_reviewed_n']==0
    assert r['source_status']['chembl']=='SOURCE_NOT_PROVIDED_NOT_ASSESSED'
    assert r['scientific_acceptance']=='NOT_APPROVED'


def test_explicit_release_required(input_archive,mock_sqlite,tmp_path):
    with pytest.raises(ValueError,match='release'):
        run(p06_zip=input_archive[1],chembl_sqlite=mock_sqlite,out=tmp_path/'out')


def test_output_cannot_overwrite_previous_report(input_archive,tmp_path):
    out=tmp_path/'report';out.mkdir()
    with pytest.raises(ValueError,match='Output must be new'):
        run(p06_zip=input_archive[1],out=out)


def test_bad_zip_duplicate_member_refused(tmp_path):
    p=tmp_path/'bad.zip'
    with zipfile.ZipFile(p,'w') as z:
        z.writestr('test','x');z.writestr('test','y')
    with pytest.raises(ValueError,match='duplicate'):
        ReadOnlyZip(p)


def test_source_artifact_cannot_be_selected(input_archive,tmp_path):
    # Existing real P0.6 selection never includes the artifact UNNAMED:0.
    p0,p06,p07=input_archive
    panel,_=rows_from_pilot(p06,p0,p07,'focus6')
    assert all(x['gene_symbol']!='UNNAMED: 0' for x in panel)