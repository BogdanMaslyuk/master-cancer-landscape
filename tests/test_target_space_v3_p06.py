"""TLS V3/P0.6 unit and fixture integration tests. Network is never used."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys

import httpx
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from pilot_target_space_v3_p06 import (
    ACCESSION_RE, Limits, ControlledClient, binding_quality, deterministic_order,
    exact_chembl_target, istrue, parse_uniprot_tsv, probe_one, read_tsv, resolve_gene,
    run, sample_panel, sha_file, stratum, valid_chembl_activity, write_tsv,
)


def _row(g, kind='CRISPR_UNMAPPED'):
    return {
        'gene_symbol':g,'source_artifact':'False',
        'path_C':'True' if kind=='CONTROL_13' else 'False',
        'candidate_v2_any_status':'True' if kind=='NONPRIORITY_MAPPED' else 'False',
        'candidate_v2_priority':'False',
        'protein_unique_reviewed_human':'True' if kind in {'NONPRIORITY_MAPPED','MAPPED_OTHER','CONTROL_13'} else 'False',
        'has_crispr':'True',
        'uniprot_accession':'P11802' if g=='CDK4' else '',
    }


def test_bool_coercion():
    assert istrue('true') and istrue('Yes') and istrue(1)
    assert not istrue('False')


def test_groups_are_exclusive():
    for k in ['CONTROL_13','NONPRIORITY_MAPPED','MAPPED_OTHER','CRISPR_UNMAPPED']:
        assert stratum(_row('ABCD',k))==k


def test_placeholder_is_excluded_even_if_crispr():
    r=_row('UNNAMED: 0')
    assert stratum(r)=='EXCLUDED_SOURCE_ARTIFACT'
    assert not sample_panel([r])


def test_sample_deterministic_across_order():
    rows=[]
    for kind in ['CONTROL_13','NONPRIORITY_MAPPED','MAPPED_OTHER','CRISPR_UNMAPPED']:
        rows += [_row(f'{kind}{i:03d}', kind) for i in range(20)]
    a=sample_panel(rows)
    b=sample_panel(rows[::-1])
    assert [(x['gene_symbol'],x['sample_stratum']) for x in a] == [(x['gene_symbol'],x['sample_stratum']) for x in b]
    assert len(a)==40
    assert sum(x['sample_stratum']=='CONTROL_13' for x in a)==13
    assert sum(x['sample_stratum']=='NONPRIORITY_MAPPED' for x in a)==12
    assert sum(x['sample_stratum']=='MAPPED_OTHER' for x in a)==7
    assert sum(x['sample_stratum']=='CRISPR_UNMAPPED' for x in a)==8


def test_sample_small_panel_no_duplicates():
    panel=sample_panel([_row('CDK4','CONTROL_13'),_row('EGFR','MAPPED_OTHER')],max_targets=40)
    assert {x['gene_symbol'] for x in panel}=={'CDK4','EGFR'}


@pytest.mark.parametrize('size',[0,51,-2])
def test_reject_invalid_sample_size(size):
    with pytest.raises(ValueError):sample_panel([_row('CDK4')],max_targets=size)


def test_reject_duplicate_gene_key():
    with pytest.raises(ValueError,match='duplicate'):
        sample_panel([_row('CDK4'),_row('CDK4')])


def test_tsv_roundtrip(tmp_path):
    p=tmp_path/'test.tsv';write_tsv(p,[{'gene_symbol':'CDK4','hgnc_id':'HGNC:1771'}],['gene_symbol','hgnc_id'])
    assert read_tsv(p)[0]['gene_symbol']=='CDK4'


def _snapshot(n=15000, extra=''):
    rows=['Entry\tGene Names (primary)\tHGNC\tProtein names']
    rows.append('P11802\tCDK4\tHGNC:1771\tCyclin-dependent kinase 4')
    for i in range(n-1):
        # Synthetic accession format is structurally valid but not an actual UniProt ID.
        # Use Q7*****-style base36 keys to prevent accidental duplicate keys.
        digits='0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        k=i//10
        chars=[]
        for _ in range(3):
            chars.append(digits[k%36]);k//=36
        acc='Q7'+''.join(reversed(chars))+str(i%10)
        rows.append(f'{acc}\tSYN{i}\tHGNC:{100000+i}\tSynthetic protein')
    if extra:rows.append(extra)
    return ('\n'.join(rows)+'\n').encode()


def test_uniprot_full_reviewed_snapshot_parses():
    rows, mapping, issues=parse_uniprot_tsv(_snapshot())
    assert len(rows)==15000 and mapping['CDK4'][0]['hgnc_ids']=='HGNC:1771'
    assert not issues


def test_uniprot_missing_schema_fails():
    with pytest.raises(ValueError,match='schema changed'):
        parse_uniprot_tsv(b'Entry\tGene Names (primary)\nP11802\tCDK4\n')


def test_uniprot_small_response_rejected():
    with pytest.raises(ValueError,match='snapshot size'):
        parse_uniprot_tsv(b'Entry\tGene Names (primary)\tHGNC\nP11802\tCDK4\tHGNC:1771\n')


def test_uniprot_duplicate_accession_rejected():
    with pytest.raises(ValueError,match='duplicate/invalid accession'):
        parse_uniprot_tsv(_snapshot(extra='P11802\tCDK4\tHGNC:1771\tDuplicate'))


def test_resolve_uniprot_exact_and_conflict():
    _, m, _=parse_uniprot_tsv(_snapshot())
    assert resolve_gene('CDK4',m,'P11802')['mapping_status']=='EXACT_REVIEWED_WITH_HGNC'
    assert resolve_gene('CDK4',m,'P12345')['mapping_status']=='LOCAL_UNIPROT_CONFLICT'
    assert resolve_gene('FAKE',m,'')['mapping_status']=='NO_EXACT_PRIMARY_GENE'


def test_resolve_duplicate_primary_gene_is_ambiguous():
    _, m, _=parse_uniprot_tsv(_snapshot())
    m['CDK4'].append({'uniprot_accession':'Q99999','hgnc_ids':'HGNC:1771'})
    assert resolve_gene('CDK4',m,'')['mapping_status']=='AMBIGUOUS_PRIMARY_GENE'


def test_uniprot_missing_hgnc_does_not_claim_full_mapping():
    _, m, _=parse_uniprot_tsv(_snapshot())
    m['CDK4'][0]['hgnc_ids']=''
    assert resolve_gene('CDK4',m,'P11802')['mapping_status']=='EXACT_REVIEWED_HGNC_UNRESOLVED'


def _target_doc(acc='P11802', kind='SINGLE PROTEIN', organism='Homo sapiens', target='CHEMBL331'):
    return {'targets':[{'target_chembl_id':target,'organism':organism,'target_type':kind,
                        'target_components':[{'accession':acc}]}]}


def test_chembl_strict_exact_match():
    assert exact_chembl_target(_target_doc(), 'P11802') == ('CHEMBL331','EXACT_SINGLE_PROTEIN_HUMAN')
    assert not exact_chembl_target(_target_doc(kind='PROTEIN COMPLEX'), 'P11802')[0]
    assert not exact_chembl_target(_target_doc(organism='Mus musculus'), 'P11802')[0]
    assert not exact_chembl_target(_target_doc(), 'P00000')[0]


def test_chembl_duplicate_exact_targets_flagged():
    d=_target_doc()
    d['targets'].append({**d['targets'][0],'target_chembl_id':'CHEMBL999'})
    assert exact_chembl_target(d,'P11802')[1]=='AMBIGUOUS_CHEMBL_TARGET'


def _activity():
    return {'target_chembl_id':'CHEMBL331','standard_type':'Ki','standard_units':'nM',
            'standard_relation':'=','standard_value':'23','assay_type':'B',
            'assay_chembl_id':'CHEMBL123','molecule_chembl_id':'CHEMBL_M1','activity_id':444}


def test_activity_direct_candidate_strict_rules():
    r=_activity();assert valid_chembl_activity(r,'CHEMBL331','Ki')
    for key,value in [('standard_units','uM'),('standard_relation','<'),('assay_type','F'),
                      ('standard_value','0'),('target_chembl_id','CHEMBL_OTHER')]:
        test={**r,key:value}
        assert not valid_chembl_activity(test,'CHEMBL331','Ki')


def test_bindingdb_observation_not_treated_as_confirmed():
    doc={'result':{'records':[{'affinity_type':'Ki','affinity':'10','SMILE':'CCC',
                              'query':'P11802','monomerid':'A'}]}}
    result=binding_quality(doc,'P11802')
    assert result['bindingdb_candidate_records']==1
    assert result['bindingdb_experiment_verified'] is False


def test_bindingdb_unknown_response_not_assumed_absent():
    result=binding_quality({'unknown_wrapper':[]},'P11802')
    assert result['bindingdb_status']=='NO_PARSEABLE_RECORDS_OR_EMPTY'
    assert not result['bindingdb_experiment_verified']


def test_bindingdb_wrong_protein_rejected():
    doc={'result':{'affinitytype':'Kd','affinity':'50','smile':'CO','uniprot':'P99999'}}
    assert binding_quality(doc,'P11802')['bindingdb_candidate_records']==0


def test_block_outside_host_and_http(tmp_path):
    cl=ControlledClient(tmp_path)
    try:
        with pytest.raises(ValueError,match='allowlisted'):
            cl.get('https://127.0.0.1/private',{})
        with pytest.raises(ValueError,match='allowlisted'):
            cl.get('http://rest.uniprot.org/uniprotkb/stream',{})
    finally:cl.close()


def test_response_limit(tmp_path):
    transport=httpx.MockTransport(lambda request: httpx.Response(200,content=b'x'*110))
    cl=ControlledClient(tmp_path,transport=transport)
    try:
        with pytest.raises(RuntimeError,match='(SINGLE_RESPONSE_BUDGET_EXCEEDED|RESPONSE_BYTE_BUDGET_EXCEEDED)'):
            cl.get('https://rest.uniprot.org/uniprotkb/stream',{},max_bytes=50)
        assert cl.records[-1]['result']!='OK'
    finally: cl.close()


def test_request_budget_limits_count(tmp_path):
    transport=httpx.MockTransport(lambda request:httpx.Response(200,content=b'{}'))
    cl=ControlledClient(tmp_path,Limits(max_requests=1,min_delay_s=0),transport=transport)
    try:
        cl.get('https://rest.uniprot.org/uniprotkb/stream',{'q':'1'})
        with pytest.raises(RuntimeError,match='REQUEST_BUDGET_EXHAUSTED'):
            cl.get('https://rest.uniprot.org/uniprotkb/stream',{'q':'2'})
    finally: cl.close()


def test_cache_reuse_reduces_requests(tmp_path):
    transport=httpx.MockTransport(lambda request:httpx.Response(200,content=b'{}'))
    cl=ControlledClient(tmp_path,Limits(min_delay_s=0),transport=transport)
    try:
        a=cl.get('https://rest.uniprot.org/uniprotkb/stream',{'q':'same'})
        b=cl.get('https://rest.uniprot.org/uniprotkb/stream',{'q':'same'})
        assert a[0]==b[0] and cl.count==1 and b[1]['cache']=='REUSED'
    finally:cl.close()


@pytest.fixture
def p0(tmp_path):
    src=tmp_path/'P0';src.mkdir()
    rows=[_row('CDK4','CONTROL_13'),_row('EGFR','NONPRIORITY_MAPPED'),
          _row('UNNAMED: 0')]
    write_tsv(src/'target_universe.tsv',rows,list(rows[0]))
    (src/'coverage_summary.json').write_text(json.dumps({
        'contract':'mcl-tls-v3-p0-inventory-v1','technical_status':'REVIEW_REQUIRED',
        'valid_gene_candidates_n':2,'suspect_source_placeholders_n':1}))
    (src/'source_manifest.json').write_text('{}')
    return src


def test_realistic_offline_pilot_uses_no_network(p0,tmp_path):
    class FailTransport(httpx.BaseTransport):
        def handle_request(self,request):
            raise AssertionError('Offline request must never happen')
    before={p.name:sha_file(p) for p in p0.iterdir()}
    summary=run(p0,tmp_path/'out',max_targets=40,transport=FailTransport())
    assert summary['sample_size']==2 and summary['requests_n']==0
    assert summary['external_mode']=='OFFLINE_ONLY'
    assert summary['scientific_acceptance']=='NOT_APPROVED'
    assert before=={p.name:sha_file(p) for p in p0.iterdir()}
    selection=read_tsv(tmp_path/'out'/'p06_selection.tsv')
    assert all(row['gene_symbol']!='UNNAMED: 0' for row in selection)


def test_old_p0_is_rejected(p0,tmp_path):
    p=p0/'coverage_summary.json'
    meta=json.loads(p.read_text());meta.pop('valid_gene_candidates_n');p.write_text(json.dumps(meta))
    with pytest.raises(ValueError,match='rerun current'):
        run(p0,tmp_path/'out')


def test_p0_input_cannot_be_overwritten(p0):
    with pytest.raises(ValueError):run(p0,p0/'out')


def test_technically_blocked_p0_cannot_be_piloted(p0,tmp_path):
    p=p0/'coverage_summary.json';meta=json.loads(p.read_text());meta['technical_status']='BLOCKED';p.write_text(json.dumps(meta))
    with pytest.raises(ValueError,match='BLOCKED'):
        run(p0,tmp_path/'out')


def test_full_mocked_external_pilot_and_provenance(p0,tmp_path):
    def respond(request):
        path=request.url.path
        if path.endswith('/stream'):
            return httpx.Response(200,content=_snapshot(),headers={'x-uniprot-release':'2026_04','x-total-results':'15000'})
        if path.endswith('/target.json'):
            return httpx.Response(200,json=_target_doc())
        if path.endswith('/activity.json'):
            typ=request.url.params.get('standard_type','')
            return httpx.Response(200,json={'activities':[_activity()] if typ=='Ki' else [],
                                            'page_meta':{'total_count':1 if typ=='Ki' else 0}})
        if path.endswith('/assay/CHEMBL123.json'):
            return httpx.Response(200,json={'assay_type':'B','confidence_score':9,
                                            'target_chembl_id':'CHEMBL331'})
        if path.endswith('/getLigandsByUniprots'):
            return httpx.Response(200,json={'items':[{'affinitytype':'Ki','affinity':'23',
                                                     'smile':'CC','query':'P11802'}]})
        return httpx.Response(404)
    transport=httpx.MockTransport(respond)
    summary=run(p0,tmp_path/'out',enable_external_pilot=True,max_targets=1,
                transport=transport,limits=Limits(min_delay_s=0))
    assert summary['external_mode']=='EXPLICIT_BOUNDED_PILOT'
    assert summary['sample_size']==1
    assert summary['mcl_exact_reviewed_with_hgnc_n']==1
    assert summary['sample_chembl_exact_assay_examples_n']==1
    assert summary['requests_n']==7  # one UniProt, one ChEMBL target, three activities, one assay, one BindingDB
    row=read_tsv(tmp_path/'out'/'p06_external_coverage.tsv')[0]
    assert row['chembl_status']=='EXAMPLE_DIRECT_ASSAY_VERIFIED'
    assert row['bindingdb_experiment_verified']=='False'
    assert (tmp_path/'out'/'uniprot_reviewed_human.tsv').is_file()
    assert (tmp_path/'out'/'mcl_reviewed_uniprot_mapping.tsv').is_file()


def test_remote_404_is_source_error_not_negative(tmp_path):
    def bad(req):return httpx.Response(404,json={'error':'down'})
    cl=ControlledClient(tmp_path,Limits(min_delay_s=0),transport=httpx.MockTransport(bad))
    try:
        result=probe_one(cl,'CDK4','P11802')
        assert result['chembl_status'].startswith('ERROR_')
        assert result['bindingdb_status'].startswith('ERROR_')
    finally:cl.close()


def test_uniprot_total_header_consistency_checked(p0,tmp_path):
    def wrong(req):
        if req.url.path.endswith('/stream'):
            return httpx.Response(200,content=_snapshot(),headers={'x-total-results':'15001'})
        return httpx.Response(404)
    summary=run(p0,tmp_path/'out',enable_external_pilot=True,max_targets=1,
                transport=httpx.MockTransport(wrong),limits=Limits(min_delay_s=0))
    assert summary['uniprot_snapshot']['status']=='ERROR_ValueError'
    assert summary['requests_n']==1
    assert summary['sample_chembl_exact_assay_examples_n']==0

def test_10_character_uniprot_accession_pattern():
    assert ACCESSION_RE.fullmatch('A0A0B4J2D5')
    assert ACCESSION_RE.fullmatch('P11802')
    assert not ACCESSION_RE.fullmatch('invalid123')

def test_bindingdb_smiles_alias_allowed_but_not_verified():
    doc={'result': {'affinitytype':'Kd','affinity':'12','smiles':'CCCO','query':'P11802'}}
    row=binding_quality(doc,'P11802')
    assert row['bindingdb_candidate_records']==1
    assert not row['bindingdb_experiment_verified']


def test_chembl_injected_assay_url_is_not_requested(tmp_path):
    calls=[]
    def mock(req):
        calls.append(req.url.path)
        if req.url.path.endswith('/target.json'):
            return httpx.Response(200,json=_target_doc())
        if req.url.path.endswith('/activity.json'):
            row={**_activity(),'assay_chembl_id':'CHEMBL123/../../other'}
            typ=req.url.params.get('standard_type','')
            return httpx.Response(200,json={'activities':[row] if typ=='Ki' else [],
                                            'page_meta':{'total_count':1 if typ=='Ki' else 0}})
        if req.url.path.endswith('/getLigandsByUniprots'):
            return httpx.Response(200,json={'items':[]})
        return httpx.Response(404)
    client=ControlledClient(tmp_path,Limits(min_delay_s=0),transport=httpx.MockTransport(mock))
    try:
        row=probe_one(client,'CDK4','P11802')
        assert row['chembl_status']=='CANDIDATES_EXIST_NOT_INDEPENDENTLY_VERIFIED'
        assert not any('/assay/' in x for x in calls)
    finally:client.close()

def test_bindingdb_real_provider_name_only_fields_are_not_direct_accession_evidence():
    """BindingDB can return different proteins under one requested UniProt."""
    doc={'getLindsByUniprotsResponse':{'affinities':[
        {'query':'Protein farnesyltransferase subunit beta',
         'affinity_type':'Kd','affinity':'475','smile':'CCO','monomerid':'1'},
        {'query':'Protein farnesyltransferase/geranylgeranyltransferase type-1 subunit alpha',
         'affinity_type':'IC50','affinity':'8.9','smile':'CCN','monomerid':'2'},
        {'query':'Protein farnesyltransferase subunit beta',
         'affinity_type':'IC50','affinity':'>30000','smile':'CCC','monomerid':'3'}
    ]}}
    r=binding_quality(doc,'P49356')
    assert r['bindingdb_records_parsed']==3
    assert r['bindingdb_name_only_candidate_rows']==2
    assert r['bindingdb_distinct_query_names_n']==2
    assert r['bindingdb_mixed_query_names'] is True
    assert r['bindingdb_candidate_records']==0
    assert r['bindingdb_status']=='NAME_ONLY_EVIDENCE_TARGET_UNVERIFIED'
    assert r['bindingdb_experiment_verified'] is False


def test_bindingdb_explicit_wrong_accession_is_not_name_only_evidence():
    doc={'x':[{'query':'Protein X','uniprot':'P99999',
               'affinity_type':'Ki','affinity':'15',
               'smile':'CCO'}]}
    r=binding_quality(doc,'P11802')
    assert r['bindingdb_candidate_records']==0
    assert r['bindingdb_name_only_candidate_rows']==0


def test_bindingdb_name_only_does_not_upgrade_direct_evidence():
    doc={'x':[{'query':'Receptor tyrosine-protein kinase erbB-3',
               'affinity_type':'IC50','affinity':'10',
               'smile':'CCO','doi':'10.1000/example'}]}
    r=binding_quality(doc,'P21860')
    assert r['bindingdb_name_only_candidate_rows']==1
    assert r['bindingdb_candidate_records']==0
    assert r['bindingdb_experiment_verified'] is False
