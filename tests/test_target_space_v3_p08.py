"""Synthetic offline tests; no actual user science data or internet required."""
from __future__ import annotations

import csv
import json
import sqlite3
import sys
import zipfile
from pathlib import Path

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_target_space_v3_p08 import (
    FOCUS6,
    P0_CONTRACT,
    P06_CONTRACT,
    P07_CONTRACT,
    ReadOnlyZip,
    bindingdb_local,
    chembl_local,
    open_binding_tsv,
    parse_nm,
    rows_from_pilot,
    run,
    shafile,
    source_metadata,
)


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
    (p07dir/'p07_summary.json').write_text(json.dumps({'contract':P07_CONTRACT,
        'scientific_acceptance':'NOT_APPROVED',
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
    CREATE TABLE docs(doc_id INTEGER PRIMARY KEY,chembl_id TEXT,doi TEXT,pubmed_id TEXT);
    INSERT INTO target_dictionary VALUES (1,'CHEMBL100','SINGLE PROTEIN','Homo sapiens');
    INSERT INTO target_components VALUES (1,100);
    INSERT INTO component_sequences VALUES (100,'P78527');
    INSERT INTO assays VALUES (10,1,'B',9,'CHEMBLASSAY10');
    INSERT INTO molecule_dictionary VALUES (1,'CHEMBLMOL1');
    INSERT INTO compound_structures VALUES (1,'CCO');
    INSERT INTO docs VALUES (1,'CHEMBLDOC1','10.1000/sample','111');
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
          ['14','CO','KEY5','MON5','Other','50','','','P78527','','1','','555','BindingDB','Mus musculus'],
          ['15','CCCl','KEY6','MON6','Missing species','70','','','P78527','','1','','666','BindingDB',''],
          ['16','CCBr','KEY7','MON7','Mouse complex','80','','','P78527','P11802','2','','777','BindingDB','Mus musculus'],
          ['17','CCI','KEY8','MON8','Unknown chain count','90','','','P78527','','','','888','BindingDB','Homo sapiens'],
          ['18','CCF','KEY9','MON9','Missing-species complex','100','','','P78527','P11802','2','','999','BindingDB','']]
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


def test_p07_contract_is_strictly_required(input_archive,tmp_path):
    p0,p06,p07=input_archive
    broken=tmp_path/'P07-wrong-contract.zip'
    with zipfile.ZipFile(p07) as src,zipfile.ZipFile(broken,'w') as dst:
        meta=json.loads(src.read('p07_summary.json'));meta['contract']='wrong-contract'
        dst.writestr('p07_summary.json',json.dumps(meta))
    with pytest.raises(ValueError,match='P0.7 contract mismatch'):
        rows_from_pilot(p06,p0,broken,'focus6')


def test_p07_recorded_source_sha_must_be_canonical(input_archive,tmp_path):
    p0,p06,p07=input_archive
    broken=tmp_path/'P07-bad-sha.zip'
    with zipfile.ZipFile(p07) as src,zipfile.ZipFile(broken,'w') as dst:
        meta=json.loads(src.read('p07_summary.json'))
        meta['source_archives_sha256'][p06.name]='not-a-sha256'
        dst.writestr('p07_summary.json',json.dumps(meta))
    with pytest.raises(ValueError,match='lowercase SHA-256'):
        rows_from_pilot(p06,p0,broken,'focus6')


def test_sample40_expands_to_available_sample(input_archive):
    p0,p06,p07=input_archive
    panel,_=rows_from_pilot(p06,p0,p07,'sample40')
    assert len(panel)==6


def test_provenance_mismatch_rejected(input_archive,tmp_path):
    p0,p06,_p07=input_archive
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
    assert next(x for x in obs if x['activity_id']==900)['pmid']=='111'
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
    assert meta['rows_scanned']==9
    assert len(rows)==8
    assert rows[0]['gene_symbol']=='PRKDC'
    assert rows[0]['chain_status']=='SWISSPROT_PRIMARY_SINGLE_CHAIN_MATCH'
    assert rows[0]['species_status']=='HUMAN_SPECIES_MATCH'
    assert rows[1]['chain_status']=='REVIEW_MULTI_CHAIN'
    assert rows[1]['species_status']=='HUMAN_SPECIES_MATCH'
    assert rows[2]['measurement_status']=='CENSORED_NOT_EXACT'
    assert rows[3]['chain_status']=='SWISSPROT_PRIMARY_SINGLE_CHAIN_MATCH'
    assert rows[3]['species_status']=='REVIEW_SPECIES_MISMATCH'
    assert cov[0]['chain_match_rows']==8


def test_bindingdb_species_and_chain_statuses_are_independent(input_archive,mock_bindingdb):
    panel,_=rows_from_pilot(input_archive[1],None,None,'focus6')
    rows,cov,_=bindingdb_local(mock_bindingdb,panel,100)
    by_id={x['bindingdb_reactant_set_id']:x for x in rows}
    assert by_id['15']['chain_status']=='SWISSPROT_PRIMARY_SINGLE_CHAIN_MATCH'
    assert by_id['15']['species_status']=='REVIEW_SPECIES_MISSING'
    assert by_id['16']['chain_status']=='REVIEW_MULTI_CHAIN'
    assert by_id['16']['species_status']=='REVIEW_SPECIES_MISMATCH'
    assert by_id['17']['chain_status']=='REVIEW_CHAIN_COUNT_UNKNOWN'
    assert by_id['17']['species_status']=='HUMAN_SPECIES_MATCH'
    assert by_id['18']['chain_status']=='REVIEW_MULTI_CHAIN'
    assert by_id['18']['species_status']=='REVIEW_SPECIES_MISSING'
    assert by_id['15']['raw_ki_nm']=='70' and by_id['15']['raw_kd_nm']=='' and by_id['15']['raw_ic50_nm']==''
    prkdc=next(x for x in cov if x['gene_symbol']=='PRKDC')
    assert prkdc['single_chain_objects_n']==4
    assert prkdc['multi_chain_objects_n']==3
    assert prkdc['unknown_chain_count_objects_n']==1
    assert prkdc['species_missing_objects_n']==2
    assert prkdc['species_mismatch_objects_n']==2


def test_bindingdb_preserves_all_raw_affinity_fields_including_invalid(input_archive,mock_bindingdb):
    with mock_bindingdb.open('a',newline='') as f:
        csv.writer(f,delimiter='\t').writerow(
          ['19','CCX','KEY10','MON10','Raw values',' >= 4 ','5e1','~7','P78527','','1',
           '10.1000/raw','1010','BindingDB','Homo sapiens'])
    panel,_=rows_from_pilot(input_archive[1],None,None,'focus6')
    rows,cov,_=bindingdb_local(mock_bindingdb,panel,100)
    raw=[x for x in rows if x['bindingdb_reactant_set_id']=='19']
    assert len(raw)==3
    assert {(x['raw_ki_nm'],x['raw_kd_nm'],x['raw_ic50_nm']) for x in raw}=={(' >= 4 ','5e1','~7')}
    assert {x['measurement_value_status'] for x in raw}=={'CENSORED','EXACT_NUMERIC','INVALID_OR_MISSING'}
    assert next(x for x in cov if x['gene_symbol']=='PRKDC')['rows_invalid']==1


def test_bindingdb_zip_stream(input_archive,mock_bindingdb,tmp_path):
    z=tmp_path/'db.zip'
    with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as f:f.write(mock_bindingdb,arcname=mock_bindingdb.name)
    panel,_=rows_from_pilot(input_archive[1],None,None,'focus6')
    with open_binding_tsv(z) as r: assert next(r)[0]=='BindingDB Reactant_set_id'
    rows,_cov,_=bindingdb_local(z,panel,100)
    assert len(rows)==8


def test_bindingdb_reject_no_primary_chain(input_archive,tmp_path):
    p=tmp_path/'bad.tsv';p.write_text('Target Name\tKi (nM)\nKinase\t10\n')
    panel,_=rows_from_pilot(input_archive[1],None,None,'focus6')
    with pytest.raises(ValueError,match='SwissProt'):
        bindingdb_local(p,panel,100)


def test_bindingdb_limit_reports_overlimit_not_invalid(input_archive,mock_bindingdb):
    panel,_=rows_from_pilot(input_archive[1],None,None,'focus6')
    rows,cov,_=bindingdb_local(mock_bindingdb,panel,1)
    assert len(rows)==1 and cov[0]['truncated']
    assert cov[0]['rows_over_limit']==7
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
    assert one['bindingdb_single_chain_human_observations_n']==2 # exact and censored, not PYZ
    assert one['cross_source_potential_duplicate_pairs_n']==1
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
    with (tmp_path/'preflight'/'p08_target_evidence_counts.tsv').open(newline='') as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==6
    assert all(x['chembl_status']=='SOURCE_NOT_PROVIDED_NOT_ASSESSED' for x in rows)
    assert all(x['chembl_exact_numeric_n']=='' for x in rows)
    assert all(x['cross_source_duplicate_status']=='NOT_ASSESSED_SOURCE_MISSING' for x in rows)


def test_explicit_release_required(input_archive,mock_sqlite,tmp_path):
    with pytest.raises(ValueError,match='release'):
        run(p06_zip=input_archive[1],chembl_sqlite=mock_sqlite,out=tmp_path/'out')


def test_publisher_sha_release_and_license_metadata(mock_sqlite):
    digest=shafile(mock_sqlite)
    meta=source_metadata(mock_sqlite,'ChEMBL','chembl_36',digest,publisher_sha256=digest,
                         license_id='ChEMBL license',license_url='https://example.test/license')
    assert meta['publisher_sha256_status']=='MATCH'
    assert meta['filename_release_status']=='MATCH'
    assert meta['license_status']=='PROVIDED_UNVERIFIED'


def test_publisher_sha_mismatch_is_rejected(mock_sqlite):
    with pytest.raises(ValueError,match='publisher SHA-256 mismatch'):
        source_metadata(mock_sqlite,'ChEMBL','chembl_36',shafile(mock_sqlite),publisher_sha256='0'*64)


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
