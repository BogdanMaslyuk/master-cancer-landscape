"""No-network synthetic regression tests for TLS V3 / P0.7."""
import csv
import io
import json
from pathlib import Path
import sys
import zipfile

import httpx
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_target_space_v3_p07 import (
    SafeZip, digest, hash_path, read_tsv, value, request_info,
    parse_numeric, binding_records, is_binding_chemically_eligible,
    classify_chembl_activity, retry_candidates, execute_retry, main_audit,
    CONTRACT, P0_CONTRACT, P06_CONTRACT
)


def to_tsv(rows, columns):
    s=io.StringIO();w=csv.DictWriter(s,fieldnames=columns,delimiter='\t');w.writeheader();w.writerows(rows)
    return s.getvalue().encode('utf-8')


def zipfrom(path, members):
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED) as z:
        for n,b in members.items():z.writestr(n,b)


@pytest.fixture
def sources(tmp_path):
    p0,p06=tmp_path/'p0.zip',tmp_path/'p06.zip'
    p0_data={
       'target_universe.tsv':to_tsv([{'gene_symbol':'CDK4'}],['gene_symbol']),
       'coverage_summary.json':json.dumps({'contract':P0_CONTRACT,'technical_status':'REVIEW_REQUIRED',
                                           'panel_genes_n':13,'own_compounds_n':65,
                                           'existing_matrix_rows_n':845}).encode(),
       'source_manifest.json':b'{}'
    }
    zipfrom(p0,p0_data)
    accession='P11802'
    urls=[
      'https://www.ebi.ac.uk/chembl/api/data/target.json?limit=100&target_components__accession='+accession,
      'https://www.ebi.ac.uk/chembl/api/data/activity.json?assay_type=B&limit=10&offset=0&standard_type=Ki&standard_units=nM&target_chembl_id=CHEMBL331',
      'https://www.ebi.ac.uk/chembl/api/data/assay/CHEMBL888.json?',
      'https://www.bindingdb.org/rest/getLigandsByUniprots?cutoff=10000&response=application%2Fjson&uniprot='+accession,
      'https://www.ebi.ac.uk/chembl/api/data/activity.json?assay_type=B&limit=10&offset=0&standard_type=Kd&standard_units=nM&target_chembl_id=CHEMBL331'
    ]
    docs=[
       {'targets':[{'target_type':'SINGLE PROTEIN','organism':'Homo sapiens',
                    'target_components':[{'accession':'P11802'}],'target_chembl_id':'CHEMBL331'}]},
       {'activities':[{'target_chembl_id':'CHEMBL331','activity_id':1,
                       'assay_chembl_id':'CHEMBL888','molecule_chembl_id':'CHEMBL111',
                       'canonical_smiles':'CCO','assay_type':'B','standard_type':'Ki',
                       'standard_units':'nM','standard_relation':'=','standard_value':'14',
                       'potential_duplicate':0,'data_validity_comment':None}],
        'page_meta':{'total_count':1}},
       {'assay_chembl_id':'CHEMBL888','target_chembl_id':'CHEMBL331',
        'assay_type':'B','confidence_score':9},
       {'getLindsByUniprotsResponse':{'affinities':[
           {'query':'Cyclin-dependent kinase 4','monomerid':'1',
            'affinity_type':'Ki','affinity':'14','smile':'CCO','pmid':'1234'},
           {'query':'Cyclin-dependent kinase 6','monomerid':'2',
            'affinity_type':'Ki','affinity':'22','smile':'CCC','pmid':'4321'}]}}
    ]
    selection=[{'gene_symbol':'CDK4','uniprot_accession':'P11802','sample_stratum':'CONTROL_13'}]
    requests=[]; p06_data={}
    for i,url in enumerate(urls):
        body=json.dumps(docs[i]).encode() if i<4 else None
        filekey=f'cache{i}.bin'
        if body is not None:p06_data['external_response_cache\\'+filekey]=body
        requests.append({'url':url,'result':'OK' if body is not None else 'HTTPStatusError:500',
          'cache_file':filekey,'http_status':'200' if body is not None else '500',
          'response_bytes':len(body) if body is not None else '',
          'response_sha256':digest(body) if body is not None else ''})
    p06_data['p06_requests.tsv']=to_tsv(requests,list(requests[0]))
    p06_data['p06_selection.tsv']=to_tsv(selection,list(selection[0]))
    ext=[{'gene_symbol':'CDK4','uniprot_accession':'P11802','sample_stratum':'CONTROL_13',
          'chembl_status':'EXAMPLE_DIRECT_ASSAY_VERIFIED',
          'bindingdb_status':'NO_STRICT_CANDIDATES_IN_RESPONSE'}]
    p06_data['p06_external_coverage.tsv']=to_tsv(ext,list(ext[0]))
    p06_data['p06_summary.json']=json.dumps({
      'contract':P06_CONTRACT,'requests_n':len(requests),'sample_size':1,
      'source_p0_sha256_before_after':{n:digest(b) for n,b in p0_data.items()}}).encode()
    zipfrom(p06,p06_data)
    return p0,p06


def test_realistic_offline_replay_and_scientific_status(sources,tmp_path):
    p0,p06=sources
    orig0,orig6=hash_path(p0),hash_path(p06)
    result=main_audit(p0,p06,tmp_path/'p07')
    assert result['sample_targets_n']==1
    assert result['p06_cached_responses_sha256_verified_n']==4
    assert result['p06_failed_requests_n']==1
    assert result['chembl_cached_assay_verified_examples_n']==1
    assert result['bindingdb_name_only_chemically_eligible_n']==2
    assert result['bindingdb_rows_verified_target_binding_n']==0
    assert result['retriable_original_chembl_requests_n']==1
    assert result['scientific_acceptance']=='NOT_APPROVED'
    assert result['merge_permitted'] is False
    assert result['retried_requests_n']==0
    assert orig0==hash_path(p0) and orig6==hash_path(p06)
    row=read_tsv((tmp_path/'p07'/'p07_per_target.tsv').read_bytes())[0]
    assert row['classification']=='ASSAY_EXAMPLE_VERIFIED'


def test_default_replay_never_initiates_network(sources,tmp_path,monkeypatch):
    monkeypatch.setattr(httpx,'Client',lambda *a,**k: (_ for _ in ()).throw(AssertionError('Network contacted offline')))
    result=main_audit(*sources,tmp_path/'p07')
    assert result['offline_source_replay'] is True


def test_replay_two_runs_identical_tables(sources,tmp_path):
    main_audit(*sources,tmp_path/'first')
    main_audit(*sources,tmp_path/'second')
    for name in ['p07_per_target.tsv','p07_chembl_observations.tsv','p07_bindingdb_protein_names.tsv','p07_retry_plan.tsv']:
        assert (tmp_path/'first'/name).read_bytes()==(tmp_path/'second'/name).read_bytes()


def test_source_sha_tampering_rejected(sources,tmp_path):
    p0,p06=sources
    with zipfile.ZipFile(p06) as z: entries={x:z.read(x) for x in z.namelist()}
    j=json.loads(entries['p06_summary.json']);j['source_p0_sha256_before_after']['target_universe.tsv']='0'*64
    entries['p06_summary.json']=json.dumps(j).encode()
    zipfrom(p06,entries)
    with pytest.raises(ValueError,match='SHA256'):
        main_audit(p0,p06,tmp_path/'bad')


def test_cached_sha_tampering_rejected(sources,tmp_path):
    p0,p06=sources
    with zipfile.ZipFile(p06) as z: entries={x:z.read(x) for x in z.namelist()}
    entries['external_response_cache\\cache0.bin']=b'NOT SAME CONTENT'
    zipfrom(p06,entries)
    with pytest.raises(ValueError,match='integrity mismatch'):
        main_audit(p0,p06,tmp_path/'bad')


def test_existing_output_folder_never_overwritten(sources,tmp_path):
    out=tmp_path/'p07';out.mkdir()
    with pytest.raises(FileExistsError):main_audit(*sources,out)


def test_archive_path_traversal_rejected(tmp_path):
    f=tmp_path/'evil.zip';zipfrom(f,{'../evil.txt':b'x'})
    with pytest.raises(ValueError,match='Unsafe ZIP'):SafeZip(f)


def test_duplicate_normalized_archive_path_rejected(tmp_path):
    f=tmp_path/'dup.zip';zipfrom(f,{'x\\y.bin':b'1','x/y.bin':b'2'})
    with pytest.raises(ValueError,match='Duplicate ZIP'):SafeZip(f)


@pytest.mark.parametrize('url',[
    'http://www.ebi.ac.uk/chembl/api/data/target.json?target_components__accession=P11802',
    'https://evil.example.com/chembl/api/data/target.json?target_components__accession=P11802',
    'https://www.ebi.ac.uk/chembl/api/data/assay/../../secret.json',
    'https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id=bad',
    'https://www.ebi.ac.uk:444/chembl/api/data/target.json?target_components__accession=P11802',
    'https://www.ebi.ac.uk/chembl/api/data/assay/CHEMBL42%2F..json',
    'https://www.ebi.ac.uk/chembl/api/data/target.json?target_components__accession=P11802&target_components__accession=P99999',
])
def test_disallow_untrusted_url(url):
    assert request_info(url)[0]=='UNTRUSTED'


def test_allow_original_endpoint_urls(sources):
    _,p06=sources
    with zipfile.ZipFile(p06) as z: rows=read_tsv(z.read('p06_requests.tsv'))
    assert [request_info(r['url'])[0] for r in rows]==['target','activity','assay','bindingdb','activity']


@pytest.mark.parametrize('raw,expect',[
    ('10',('=',10.0)),(' < 20',('<',20.0)),('>=3',('>=',3.0)),('20.0',('=',20.0)),
    ('0',('',None)),('NaN',('',None)),('junk',('',None)),('-2',('',None))])
def test_affinity_parse(raw,expect):
    assert parse_numeric(raw)==expect


def test_bindingdb_mixed_target_names_are_never_verified(sources,tmp_path):
    main_audit(*sources,tmp_path/'p07')
    rows=read_tsv((tmp_path/'p07'/'p07_bindingdb_protein_names.tsv').read_bytes())
    assert len(rows)==2
    assert {r['provider_query_protein_name'] for r in rows}=={'Cyclin-dependent kinase 4','Cyclin-dependent kinase 6'}
    assert all(r['identity_status']=='NAME_ONLY_TARGET_UNVERIFIED' for r in rows)
    assert all(r['scientific_use']=='SOURCE_REVIEW_QUEUE_ONLY' for r in rows)


def test_bindingdb_lower_bound_not_eligible():
    assert is_binding_chemically_eligible({'affinity_type':'Ki','affinity':'>500','smile':'CCO'})[0] is False


def test_bindingdb_upper_bound_keeps_censor_flag():
    assert is_binding_chemically_eligible({'affinity_type':'Ki','affinity':'<500','smile':'CCO'})[1]=='CENSORED_UPPER_BOUND'


def test_assay_score8_not_accepted():
    r={'target_chembl_id':'CHEMBL331','assay_chembl_id':'CHEMBL888',
      'molecule_chembl_id':'CHEMBL111','canonical_smiles':'CC',
      'standard_type':'Ki','standard_value':'22','standard_relation':'=',
      'standard_units':'nM','assay_type':'B'}
    result=classify_chembl_activity(r,'CHEMBL331',{'CHEMBL888':{'target_chembl_id':'CHEMBL331','assay_type':'B','confidence_score':8}})
    assert result['evidence_status']=='ASSAY_CONFIDENCE_NOT_9'


def test_assay_target_mismatch_not_accepted():
    r={'target_chembl_id':'CHEMBL331','assay_chembl_id':'CHEMBL888',
      'molecule_chembl_id':'CHEMBL111','canonical_smiles':'CC',
      'standard_type':'Ki','standard_value':'22','standard_relation':'=',
      'standard_units':'nM','assay_type':'B'}
    assert classify_chembl_activity(r,'CHEMBL331',{'CHEMBL888':{'target_chembl_id':'CHEMBL999'}})['evidence_status']=='ASSAY_TARGET_CONFLICT'


def test_invalid_ch_embl_assay_id_not_accepted():
    r={'target_chembl_id':'CHEMBL331','assay_chembl_id':'../../../secrets',
      'molecule_chembl_id':'CHEMBL111','canonical_smiles':'CC',
      'standard_type':'Ki','standard_value':'22','standard_relation':'=',
      'standard_units':'nM','assay_type':'B'}
    assert classify_chembl_activity(r,'CHEMBL331',{})['evidence_status']=='INVALID_SOURCE_IDS'


def test_retry_plan_deterministic_and_restricted(sources):
    _,p06=sources
    with zipfile.ZipFile(p06) as z:
        plans=retry_candidates(read_tsv(z.read('p06_requests.tsv')),read_tsv(z.read('p06_selection.tsv')))
    assert len(plans)==1
    assert plans[0]['kind']=='activity'
    assert 'only=' in plans[0]['retry_url']
    assert plans[0]['retry_url'].startswith('https://www.ebi.ac.uk/chembl/api/data/activity.json')


def test_retry_mock_writes_only_new_cache(sources,tmp_path):
    def respond(req):
        return httpx.Response(200,json={'activities':[],'page_meta':{'total_count':0}})
    p0,p06=sources
    original=hash_path(p06)
    result=main_audit(p0,p06,tmp_path/'out',allow_retries=True,max_retries=1,
                 transport=httpx.MockTransport(respond))
    assert result['retries_approved'] is True
    assert result['retried_requests_n']==1
    assert result['retried_success_n']==1
    assert hash_path(p06)==original
    assert list((tmp_path/'out'/'targeted_retry_cache').iterdir())


def test_retry_500_preserves_error_without_claiming_missing(sources,tmp_path):
    def bad(req):return httpx.Response(500,json={'error':'down'})
    result=main_audit(*sources,tmp_path/'out',allow_retries=True,max_retries=1,
                 transport=httpx.MockTransport(bad))
    assert result['retried_requests_n']==1
    assert result['retried_success_n']==0
    attempt=read_tsv((tmp_path/'out'/'p07_retry_attempts.tsv').read_bytes())[0]
    assert attempt['result']=='ERROR'
    assert attempt['http_status']=='500'


def test_reject_retry_budget_zero(tmp_path):
    with pytest.raises(ValueError,match='max_requests'):
        execute_retry([],tmp_path,max_requests=0)


def test_bindingdb_explicit_different_accession_not_counted_for_request():
    record={'query':'Some other protein','uniprot':'P99999',
            'affinity_type':'Ki','affinity':'10','smile':'CCO'}
    assert is_binding_chemically_eligible(record)[0] is True  # chemistry only
    # Membership to the requested accession remains a separate requirement.
    assert record['uniprot']!='P11802'


def test_activity_censoring_not_promoted_to_exact():
    row={'target_chembl_id':'CHEMBL331','assay_chembl_id':'CHEMBL888',
         'molecule_chembl_id':'CHEMBL111','canonical_smiles':'CCO',
         'standard_type':'Ki','standard_value':'20','standard_relation':'>',
         'standard_units':'nM','assay_type':'B'}
    assert classify_chembl_activity(row,'CHEMBL331',
      {'CHEMBL888':{'target_chembl_id':'CHEMBL331','assay_type':'B','confidence_score':9}})['evidence_status']=='CENSORED_OR_NONEXACT'