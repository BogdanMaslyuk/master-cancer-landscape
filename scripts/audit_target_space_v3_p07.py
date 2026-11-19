#!/usr/bin/env python3
"""MCL TLS V3/P0.7: independent source-integrity and assay-identity audit.

Default: strictly offline replay of two immutable P0/P0.6 ZIPs. Optional,
explicitly approved, rate-limited, bounded ChEMBL retries; never rewrite
input ZIPs or the 65 x 13 scientific baseline.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path, PurePosixPath
import re
import time
from urllib.parse import parse_qs, urlencode, urlparse
import zipfile

import httpx

CONTRACT = 'mcl-tls-v3-p07-source-reliability-v1'
P0_CONTRACT = 'mcl-tls-v3-p0-inventory-v1'
P06_CONTRACT = 'mcl-tls-v3-p06-target-coverage-pilot-v1'
CHEMBL_ROOT = 'https://www.ebi.ac.uk/chembl/api/data/'
ACCESSION = re.compile(r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z][A-Z0-9]{2}[0-9])?)$')
CHEMBL_ID = re.compile(r'^CHEMBL[0-9]{1,12}$')
ASSAY_RE = re.compile(r'^/chembl/api/data/assay/(CHEMBL[0-9]{1,12})\.json$')
ENDPOINTS = ('Ki', 'Kd', 'IC50')
RETRIABLE_CODES = {'500','502','503','504','408','429'}
CSV_FIELDS = { 'p06_requests.tsv','p06_selection.tsv','p06_external_coverage.tsv' }
P0_INTEGRITY = ('target_universe.tsv','coverage_summary.json','source_manifest.json')
MAX_ZIP_MEMBERS = 400
MAX_MEMBER = 8_000_000
MAX_ZIP_UNCOMPRESSED = 90_000_000


def digest(body: bytes) -> str:
    return hashlib.sha256(body).hexdigest()


def hash_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def value(v) -> str:
    return str(v if v is not None else '').strip()


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        w.writeheader();w.writerows(rows)


def read_tsv(body: bytes) -> list[dict[str,str]]:
    reader=csv.DictReader(io.StringIO(body.decode('utf-8-sig')), delimiter='\t')
    if not reader.fieldnames or len(reader.fieldnames)!=len(set(reader.fieldnames)):
        raise ValueError('TSV missing or duplicate headers')
    return list(reader)


class SafeZip:
    def __init__(self,path:Path):
        self.path=path.resolve()
        self.z=zipfile.ZipFile(path)
        records={}
        total=0
        for x in self.z.infolist():
            name=x.filename.replace('\\','/')
            parts=PurePosixPath(name).parts
            if x.is_dir():continue
            if not parts or name.startswith('/') or '..' in parts or ':' in parts[0] or '\x00' in name:
                raise ValueError('Unsafe ZIP entry: '+name)
            if name in records:raise ValueError('Duplicate ZIP member: '+name)
            if x.file_size>MAX_MEMBER:raise ValueError('Too large ZIP member: '+name)
            total+=x.file_size
            if total>MAX_ZIP_UNCOMPRESSED or len(records)>=MAX_ZIP_MEMBERS:
                raise ValueError('ZIP exceeds audit size/member bounds')
            records[name]=x
        self.members=records

    def read(self,name:str) -> bytes:
        key=name.replace('\\','/')
        if key not in self.members:raise KeyError('Missing archive member '+key)
        info=self.members[key]
        with self.z.open(info) as stream:
            body=stream.read(MAX_MEMBER+1)
        if len(body)!=info.file_size or len(body)>MAX_MEMBER:
            raise ValueError('Unexpected ZIP member size '+key)
        return body

    def close(self):self.z.close()


def request_info(url: str) -> tuple[str,str,dict]:
    u=urlparse(url)
    if u.scheme!='https' or u.username or u.password or u.port not in (None,443) or u.fragment:
        return 'UNTRUSTED','',{}
    q=parse_qs(u.query,keep_blank_values=True)
    if any(len(v)!=1 for v in q.values()):return 'UNTRUSTED','',{}
    p=u.path
    if u.hostname=='rest.uniprot.org' and p=='/uniprotkb/stream':return 'uniprot','',q
    if u.hostname=='www.bindingdb.org' and p=='/rest/getLigandsByUniprots':
        ac=q.get('uniprot',[''])[0]
        return ('bindingdb',ac,q) if ACCESSION.fullmatch(ac) else ('UNTRUSTED','',q)
    if u.hostname!='www.ebi.ac.uk':return 'UNTRUSTED','',q
    if p=='/chembl/api/data/target.json':
        ac=q.get('target_components__accession',[''])[0]
        return ('target',ac,q) if ACCESSION.fullmatch(ac) else ('UNTRUSTED','',q)
    if p=='/chembl/api/data/activity.json':
        tar=q.get('target_chembl_id',[''])[0]
        return ('activity',tar,q) if CHEMBL_ID.fullmatch(tar) else ('UNTRUSTED','',q)
    match=ASSAY_RE.fullmatch(p)
    if match:return 'assay',match.group(1),q
    return 'UNTRUSTED','',q


def parse_numeric(raw: str) -> tuple[str, float|None]:
    m=re.fullmatch(r'\s*(<=|>=|<|>|=)?\s*(\d+(?:\.\d+)?)\s*',value(raw).replace(',',''))
    if not m:return '',None
    x=float(m.group(2))
    return (m.group(1) or '=',x) if 0<x<math.inf else ('',None)


def binding_records(payload:dict) -> list[dict]:
    records=[]
    def visit(node,depth=0):
        if depth>12:return
        if isinstance(node,dict):
            keys={re.sub('[^a-z0-9]','',str(k).lower()) for k in node}
            if any(k.endswith('affinitytype') for k in keys) and any(k.endswith('affinity') for k in keys) and ('smile' in keys or 'smiles' in keys):
                records.append(node)
            for v in node.values():visit(v,depth+1)
        elif isinstance(node,list):
            for v in node:visit(v,depth+1)
    visit(payload)
    return records


def field(record:dict,*aliases):
    fields={re.sub('[^a-z0-9]','',str(k).lower()):v for k,v in record.items()}
    for a in aliases:
        for k,v in fields.items():
            if k==a or k.endswith(a):return value(v)
    return ''


def is_binding_chemically_eligible(rec:dict) -> tuple[bool,str,str]:
    metric=field(rec,'affinitytype').upper()
    if metric not in {'KI','KD','IC50'}:return False,'OTHER_ENDPOINT',metric
    relation,x=parse_numeric(field(rec,'affinity'))
    if x is None or x>10000:return False,'VALUE_MISSING_OR_OUT_OF_SCOPE',metric
    if relation in ('>','>='):return False,'LOWER_BOUND',metric
    if not field(rec,'smile','smiles'):return False,'NO_SMILES',metric
    return True, 'CENSORED_UPPER_BOUND' if relation in ('<','<=') else 'EXACT_VALUE',metric


def classify_chembl_activity(row:dict,target_id:str,assay_by_id:dict) -> dict:
    assay_id=value(row.get('assay_chembl_id'))
    status='ACTIVITY_CANDIDATE'
    if value(row.get('target_chembl_id'))!=target_id:
        status='TARGET_MISMATCH'
    elif value(row.get('standard_type')) not in ENDPOINTS:
        status='ENDPOINT_NOT_ALLOWED'
    elif value(row.get('standard_units'))!='nM':
        status='UNITS_UNSUPPORTED'
    elif value(row.get('standard_relation'))!='=':
        status='CENSORED_OR_NONEXACT'
    elif value(row.get('assay_type')).upper()!='B':
        status='NON_BINDING_ASSAY'
    else:
        try:
            measured=float(row.get('standard_value'))
            if not math.isfinite(measured) or measured<=0:status='INVALID_MEASUREMENT'
        except (ValueError,TypeError):status='INVALID_MEASUREMENT'
        if not CHEMBL_ID.fullmatch(assay_id) or not CHEMBL_ID.fullmatch(value(row.get('molecule_chembl_id'))):
            status='INVALID_SOURCE_IDS'
    if status=='ACTIVITY_CANDIDATE':
        assay=assay_by_id.get(assay_id)
        if assay is None:status='ASSAY_NOT_CACHED'
        elif value(assay.get('target_chembl_id'))!=target_id:status='ASSAY_TARGET_CONFLICT'
        elif value(assay.get('assay_type')).upper()!='B':status='ASSAY_NOT_BINDING'
        elif value(assay.get('confidence_score'))!='9':status='ASSAY_CONFIDENCE_NOT_9'
        elif not value(row.get('canonical_smiles')):status='SMILES_NOT_IN_ACTIVITY'
        elif value(row.get('data_validity_comment')):status='ACTIVITY_VALIDITY_FLAG'
        elif value(row.get('potential_duplicate')) in {'1','True'}:status='POTENTIAL_DUPLICATE'
        else:status='CACHED_ASSAY_EXAMPLE_VERIFIED'
    return {'activity_id':value(row.get('activity_id')),
            'assay_chembl_id':assay_id, 'molecule_chembl_id':value(row.get('molecule_chembl_id')),
            'standard_type':value(row.get('standard_type')),
            'standard_value':value(row.get('standard_value')),
            'standard_relation':value(row.get('standard_relation')),
            'standard_units':value(row.get('standard_units')),
            'document_chembl_id':value(row.get('document_chembl_id')),
            'potential_duplicate':value(row.get('potential_duplicate')),
            'data_validity_comment':value(row.get('data_validity_comment')),
            'evidence_status':status}


def retry_projection(kind:str,q:dict) -> dict:
    """Reduce payload size without losing identity fields or relaxing assay QC."""
    v={k:xs[0] for k,xs in q.items()}
    if kind=='target':
        return {'target_components__accession':v['target_components__accession'], 'limit':'10',
                'only':'target_chembl_id,target_type,organism,tax_id,target_components'}
    if kind=='activity':
        return {'target_chembl_id':v['target_chembl_id'],'assay_type':'B',
                'standard_type':v['standard_type'],'standard_units':'nM',
                'limit':'10','offset':'0',
                'only':'activity_id,assay_chembl_id,molecule_chembl_id,standard_type,standard_units,standard_relation,standard_value,assay_type,target_chembl_id,canonical_smiles,document_chembl_id,data_validity_comment,potential_duplicate'}
    return {'only':'assay_chembl_id,target_chembl_id,assay_type,confidence_score'}


def retry_candidates(requests:list[dict], selection:list[dict]) -> list[dict]:
    priorities={x['uniprot_accession']:i for i,x in enumerate(selection) if ACCESSION.fullmatch(x['uniprot_accession'])}
    raw=[]
    for r in requests:
        if r.get('result')=='OK':continue
        kind,key,params=request_info(r['url'])
        if kind not in ('target','activity','assay'):continue
        status=value(r.get('http_status'))
        message=value(r.get('result'))
        if status not in RETRIABLE_CODES and not any(k in message for k in ('ReadTimeout','ConnectError','ConnectTimeout','BYTE_BUDGET','SINGLE_RESPONSE_BUDGET')):
            continue
        if kind=='target' and key not in priorities:continue
        if kind=='activity' and (not CHEMBL_ID.fullmatch(key) or params.get('standard_type',[''])[0] not in ENDPOINTS):continue
        if kind=='assay' and not CHEMBL_ID.fullmatch(key):continue
        transformed=retry_projection(kind,params)
        original=urlparse(r['url'])
        u=original.scheme+'://'+original.hostname+original.path+'?'+urlencode(sorted(transformed.items()))
        raw.append({'original_url':r['url'],'retry_url':u,'kind':kind,'original_error':message,
                    'original_http':status,'source_priority':priorities.get(key,999),
                    'safe_to_retry':True,'retry_requires_explicit_opt_in':True})
    # Only retry requests from the original logged failed set; deterministic.
    order={'assay':0,'activity':1,'target':2}
    return sorted(raw,key=lambda r:(order[r['kind']],r['source_priority'],r['retry_url']))


def execute_retry(plans:list[dict], out:Path, *, max_requests:int=12,
                  transport=None, sleep_between:float=0.75) -> tuple[list[dict],dict[str,bytes]]:
    if not 1<=max_requests<=24:raise ValueError('max_requests must be 1..24')
    results=[];bodies={};bytes_total=0
    cache=out/'targeted_retry_cache';cache.mkdir()
    client=httpx.Client(follow_redirects=False,timeout=16.0,transport=transport,
                        headers={'User-Agent':'MCL-TLS-V3-P07-source-audit/1.0','Accept':'application/json'})
    try:
        for i,plan in enumerate(plans[:max_requests]):
            kind,key,q=request_info(plan['retry_url'])
            if kind!=plan['kind'] or kind not in ('target','activity','assay'):
                raise ValueError('Unsafe targeted retry endpoint')
            if i and sleep_between:time.sleep(sleep_between)
            status='';detail='';body=b'';cachefile='';http_code=''
            try:
                with client.stream('GET',plan['retry_url']) as response:
                    http_code=str(response.status_code)
                    response.raise_for_status()
                    declared=response.headers.get('Content-Length')
                    if declared and int(declared)>900_000:raise ValueError('RETRY_RESPONSE_TOO_LARGE')
                    pieces=[];n=0
                    for piece in response.iter_bytes():
                        n+=len(piece)
                        if n>900_000 or bytes_total+n>8_000_000:raise ValueError('RETRY_BYTE_CAP')
                        pieces.append(piece)
                    body=b''.join(pieces)
                parsed=json.loads(body)
                if not isinstance(parsed,dict):raise ValueError('RETRY_NOT_JSON_OBJECT')
                # Fail closed if the provider projection omits the identity field.
                if kind=='target' and not isinstance(parsed.get('targets'),list):raise ValueError('RETRY_TARGET_SCHEMA')
                if kind=='activity' and not isinstance(parsed.get('activities'),list):raise ValueError('RETRY_ACTIVITY_SCHEMA')
                if kind=='assay' and not parsed.get('target_chembl_id'):raise ValueError('RETRY_ASSAY_SCHEMA')
                bytes_total+=len(body)
                cachefile=digest(plan['retry_url'].encode('utf-8'))+'.json'
                (cache/cachefile).write_bytes(body)
                bodies[plan['retry_url']]=body
                status='OK'
            except (httpx.HTTPError,ValueError,json.JSONDecodeError) as ex:
                status='ERROR';detail=f'{type(ex).__name__}:{str(ex)[:160]}'
            results.append({**plan,'result':status,'http_status':http_code,
                            'response_sha256':digest(body) if status=='OK' else '',
                            'response_bytes':len(body) if status=='OK' else 0,
                            'cache_file':cachefile,'error_detail':detail})
    finally:
        client.close()
    return results,bodies


def main_audit(p0zip:Path,p06zip:Path,out:Path,*,allow_retries:bool=False,
               max_retries:int=12,transport=None) -> dict:
    p0zip=p0zip.resolve();p06zip=p06zip.resolve();out=out.resolve()
    if p0zip==p06zip:raise ValueError('P0 and P0.6 inputs must differ')
    if not p0zip.is_file() or not p06zip.is_file():raise FileNotFoundError('Both original P0/P0.6 ZIPs are required')
    for p in (p0zip,p06zip):
        if out==p.parent or out in p.parents or p in out.parents:raise ValueError('Output overlaps input path')
    if out.exists():raise FileExistsError('Never overwrite an existing audit folder: '+str(out))
    za=SafeZip(p0zip);zb=SafeZip(p06zip)
    try:
        p0_summary=json.loads(za.read('coverage_summary.json'))
        p06_summary=json.loads(zb.read('p06_summary.json'))
        if p0_summary.get('contract')!=P0_CONTRACT or p06_summary.get('contract')!=P06_CONTRACT:
            raise ValueError('P0/P0.6 science source contract mismatch')
        if p0_summary.get('technical_status')=='BLOCKED':raise ValueError('P0 was technically blocked')
        for name in P0_INTEGRITY:
            expected=p06_summary.get('source_p0_sha256_before_after',{}).get(name)
            if not expected or digest(za.read(name))!=expected:
                raise ValueError('P0/P0.6 SHA256 source disagreement: '+name)
        requests=read_tsv(zb.read('p06_requests.tsv'))
        selection=read_tsv(zb.read('p06_selection.tsv'))
        external=read_tsv(zb.read('p06_external_coverage.tsv'))
        if len(selection)!=len(external) or len(selection)!=p06_summary['sample_size']:
            raise ValueError('Selection vs coverage size mismatch')
        if len({r['gene_symbol'] for r in selection})!=len(selection):raise ValueError('Duplicate selected gene')
        if {r['gene_symbol'] for r in selection}!={r['gene_symbol'] for r in external}:
            raise ValueError('Selection and external gene sets mismatch')
        if len(requests)!=p06_summary.get('requests_n'):raise ValueError('Request count mismatch')
        # Verify every original successful cached response by SHA; no silent loss.
        good={};integrity=[]
        failures=[]
        for r in requests:
            kind,key,params=request_info(r['url'])
            if kind=='UNTRUSTED':raise ValueError('Untrusted request in P0.6 log: '+r['url'][:80])
            row={'source':'P06','source_url':r['url'],'endpoint':kind,'key':key,
                 'http_status':r['http_status'],'result':r['result'],
                 'cache_verified':False,'response_sha256':r['response_sha256'],
                 'source_cache_member':r['cache_file']}
            if r['result']=='OK':
                name='external_response_cache/'+r['cache_file']
                body=zb.read(name)
                if digest(body)!=r['response_sha256'] or len(body)!=int(r['response_bytes']):
                    raise ValueError('P0.6 response integrity mismatch for '+name)
                row['cache_verified']=True
                good[r['url']]=body
            else:failures.append(row)
            integrity.append(row)
        plans=retry_candidates(requests,selection)
        out.mkdir(parents=True,exist_ok=False)
        retry_rows=[];retried={}
        if allow_retries:
            retry_rows,retried=execute_retry(plans,out,max_requests=max_retries,transport=transport,
                                             sleep_between=0 if transport is not None else 0.75)
            for r in retry_rows:
                integrity.append({'source':'P07_RETRY','source_url':r['retry_url'],
                                  'endpoint':r['kind'],'key':request_info(r['retry_url'])[1],
                                  'http_status':r['http_status'],'result':r['result'],
                                  'cache_verified':r['result']=='OK',
                                  'response_sha256':r['response_sha256'],
                                  'source_cache_member':r['cache_file']})
        # Target accession mappings from known external per-target rows, but do NOT
        # mistake that mapping for validated experimental binding.
        selected={r['uniprot_accession']:r['gene_symbol'] for r in external
                  if ACCESSION.fullmatch(r['uniprot_accession'])}
        bygene={r['gene_symbol']:r for r in external}
        # First pass verifies ChEMBL SINGLE PROTEIN Homo sapiens exact accession.
        target_to_gene={}
        target_conflicts=[]
        for url,raw in {**good,**retried}.items():
            kind,key,_=request_info(url)
            if kind!='target' or key not in selected:continue
            data=json.loads(raw)
            for target in data.get('targets',[]):
                if not isinstance(target,dict):continue
                if value(target.get('target_type'))!='SINGLE PROTEIN':continue
                if value(target.get('organism')).lower()!='homo sapiens' and value(target.get('tax_id'))!='9606':continue
                components=target.get('target_components') or []
                entries={value(x.get('accession')) for x in components if isinstance(x,dict)}
                if entries!={key}:continue
                tid=value(target.get('target_chembl_id'))
                if CHEMBL_ID.fullmatch(tid):
                    if tid in target_to_gene and target_to_gene[tid]!=selected[key]:
                        target_conflicts.append(tid)
                    target_to_gene[tid]=selected[key]
        for tid in target_conflicts:target_to_gene.pop(tid,None)
        # Cache assay records separately to check exactly matching assay IDs.
        assays={}
        for url,raw in {**good,**retried}.items():
            kind,key,_=request_info(url)
            if kind=='assay':
                doc=json.loads(raw)
                if value(doc.get('assay_chembl_id'))==key:assays[key]=doc
        activity_evidence=[]
        for url,raw in {**good,**retried}.items():
            kind,target,params=request_info(url)
            if kind!='activity':continue
            if target not in target_to_gene:continue
            doc=json.loads(raw)
            for entry in doc.get('activities',[]):
                if not isinstance(entry,dict):continue
                quality=classify_chembl_activity(entry,target,assays)
                activity_evidence.append({'gene_symbol':target_to_gene[target],
                    'uniprot_accession':bygene[target_to_gene[target]]['uniprot_accession'],
                    'target_chembl_id':target,'request_url':url,
                    'raw_source_sha256':digest(raw),**quality})
        # Retried projected endpoints can overlap original first pages; avoid
        # double counting identical experiment IDs when source fields agree.
        unique_rows={}
        for row in activity_evidence:
            key=(row['gene_symbol'],row['activity_id'])
            if key in unique_rows and row['evidence_status']!='CACHED_ASSAY_EXAMPLE_VERIFIED':continue
            unique_rows[key]=row
        activity_evidence=sorted(unique_rows.values(),key=lambda r:(r['gene_symbol'],r['activity_id']))
        # BindingDB observation grouping: NO name-only row is a verified target.
        bdb_groups=defaultdict(lambda:{'raw':0,'eligible':0,'exact_value':0,'upper_bound':0,'identity_declared':0,'examples':[]})
        bdb_details=[]
        empty_bindingdb_responses=[]
        for url,raw in good.items():
            kind,acc,_=request_info(url)
            if kind!='bindingdb' or acc not in selected:continue
            if not raw:
                empty_bindingdb_responses.append(selected[acc])
                bdb_details.append({'gene_symbol':selected[acc],'accession':acc,'url':url,
                                    'response_sha256':digest(raw),'parsed_rows':0,
                                    'source_status':'EMPTY_HTTP_200_NOT_NEGATIVE'})
                continue
            try:
                doc=json.loads(raw)
            except json.JSONDecodeError as ex:
                raise ValueError(f'Invalid successful BindingDB cache payload for {acc}') from ex
            records=binding_records(doc)
            for rec in records:
                pname=field(rec,'query') or '[NAME_MISSING]'
                group=bdb_groups[(selected[acc],acc,pname)]
                group['raw']+=1
                eligible,status,metric=is_binding_chemically_eligible(rec)
                explicitly_declared=field(rec,'uniprot')
                if explicitly_declared and explicitly_declared!=acc:
                    eligible=False   # Provider explicitly names another accession.
                if eligible:
                    group['eligible']+=1
                    if status=='EXACT_VALUE':group['exact_value']+=1
                    else:group['upper_bound']+=1
                    # The provider's `query` is a protein name, not an accession.
                    # If an explicit, separate UniProt field exists, record it
                    # for manual review but do not assume the assay is confirmed.
                    uniprot=field(rec,'uniprot')
                    if uniprot and uniprot==acc:group['identity_declared']+=1
                    if len(group['examples'])<2:
                        group['examples'].append({'monomerid':field(rec,'monomerid','monomer'),
                            'pmid':field(rec,'pmid'),'doi':field(rec,'doi'),
                            'metric':metric,'value':field(rec,'affinity')})
            bdb_details.append({'gene_symbol':selected[acc],'accession':acc,'url':url,
                               'response_sha256':digest(raw),'parsed_rows':len(records),
                               'source_status':'PARSED_PROVIDER_JSON'})
        bdb_rows=[]
        for (g,acc,name),v in sorted(bdb_groups.items()):
            bdb_rows.append({'gene_symbol':g,'requested_uniprot':acc,'provider_query_protein_name':name,
                'raw_rows':v['raw'],'chemically_eligible_rows':v['eligible'],
                'exact_value_rows':v['exact_value'],'upper_bound_rows':v['upper_bound'],
                'accession_declared_in_record_rows':v['identity_declared'],
                'source_record_examples_json':json.dumps(v['examples'],ensure_ascii=False),
                'identity_status':'NAME_ONLY_TARGET_UNVERIFIED',
                'scientific_use':'SOURCE_REVIEW_QUEUE_ONLY'})
        # Per-target interpreted status; never classify request failure as no ligands.
        byact=defaultdict(list)
        for row in activity_evidence:byact[row['gene_symbol']].append(row)
        bybdb=defaultdict(list)
        for row in bdb_rows:bybdb[row['gene_symbol']].append(row)
        eligible_activity_states={'ASSAY_NOT_CACHED','ASSAY_CONFIDENCE_NOT_9',
                                  'CACHED_ASSAY_EXAMPLE_VERIFIED',
                                  'ASSAY_TARGET_CONFLICT','ASSAY_NOT_BINDING',
                                  'SMILES_NOT_IN_ACTIVITY','ACTIVITY_VALIDITY_FLAG',
                                  'POTENTIAL_DUPLICATE'}
        per_target=[]
        for r in external:
            g=r['gene_symbol'];a=byact[g];b=bybdb[g]
            n=sum(x['evidence_status']=='CACHED_ASSAY_EXAMPLE_VERIFIED' for x in a)
            verified_target=bool(n)
            cache_bdb=sum(x['chemically_eligible_rows'] for x in b)
            per_target.append({'gene_symbol':g,'uniprot_accession':r['uniprot_accession'],
                'sample_stratum':r['sample_stratum'],
                'p06_chembl_status':r['chembl_status'],
                'p06_bindingdb_status':r['bindingdb_status'],
                'chembl_activity_rows_cached':len(a),
                'chembl_exact_numeric_rows_cached':sum(x['evidence_status'] in eligible_activity_states for x in a),
                'chembl_assay_verified_examples_n':n,
                'chembl_verified_example_not_full_coverage':verified_target,
                'bindingdb_raw_rows_cached':sum(x['raw_rows'] for x in b),
                'bindingdb_name_only_chemical_rows':cache_bdb,
                'bindingdb_identity_verified':False,
                'unassessed_source_is_negative':False,
                'classification':'ASSAY_EXAMPLE_VERIFIED' if verified_target else
                 'CHEMBL_ACTIVITY_CANDIDATES' if any(x['evidence_status'] in eligible_activity_states for x in a) else
                 'CHEMBL_ONLY_CENSORED_OR_UNUSABLE' if a else
                 'BINDINGDB_NAME_ONLY_UNVERIFIED' if cache_bdb else
                 'NOT_ESTABLISHED_FROM_CACHED_DATA'})
        write_tsv(out/'p07_per_target.tsv',per_target,list(per_target[0]))
        write_tsv(out/'p07_chembl_observations.tsv',activity_evidence,
                  ['gene_symbol','uniprot_accession','target_chembl_id','activity_id',
                   'assay_chembl_id','molecule_chembl_id','standard_type','standard_value',
                   'standard_relation','standard_units','document_chembl_id',
                   'potential_duplicate','data_validity_comment','evidence_status',
                   'request_url','raw_source_sha256'])
        write_tsv(out/'p07_bindingdb_protein_names.tsv',bdb_rows,
                  ['gene_symbol','requested_uniprot','provider_query_protein_name',
                   'raw_rows','chemically_eligible_rows','exact_value_rows','upper_bound_rows',
                   'accession_declared_in_record_rows','source_record_examples_json',
                   'identity_status','scientific_use'])
        write_tsv(out/'p07_bindingdb_cache_index.tsv',bdb_details,
                  ['gene_symbol','accession','url','response_sha256','parsed_rows','source_status'])
        write_tsv(out/'p07_request_integrity.tsv',integrity,
                  ['source','source_url','endpoint','key','http_status','result',
                   'cache_verified','response_sha256','source_cache_member'])
        write_tsv(out/'p07_retry_plan.tsv',plans,
                  ['kind','source_priority','original_http','original_error','original_url','retry_url',
                   'safe_to_retry','retry_requires_explicit_opt_in'])
        write_tsv(out/'p07_retry_attempts.tsv',retry_rows,
                  ['kind','original_url','retry_url','original_http','original_error','result',
                   'http_status','response_sha256','response_bytes','cache_file','error_detail'])
        qc=[{'severity':'REVIEW','code':'SAMPLE_NOT_PROTEOME','detail':'40 preselected targets; no population inference'},
            {'severity':'REVIEW','code':'BINDINGDB_NAME_ONLY','detail':'Target identity not proven by protein names'},
            {'severity':'REVIEW','code':'CHEMBL_FIRST_PAGE_ONLY','detail':'Only cached first-page activity rows and assays, not complete activity universe'},
            {'severity':'REVIEW','code':'RETRY_NOT_NEEDED_FOR_BASELINE','detail':'Original 65x13 and scientific sources not modified'}]
        if empty_bindingdb_responses:
            qc.append({'severity':'REVIEW','code':'BINDINGDB_HTTP200_EMPTY_BODIES',
                       'detail':f'{len(empty_bindingdb_responses)} empty 200 OK bodies are not evidence of absence'})
        if failures:qc.append({'severity':'REVIEW','code':'P06_SOURCE_ERRORS','detail':str(len(failures))+' failed original requests; not negative evidence'})
        if target_conflicts:qc.append({'severity':'BLOCKER','code':'CHEMBL_TARGET_ID_CONFLICT','detail':','.join(sorted(set(target_conflicts)))})
        write_tsv(out/'p07_QC.tsv',qc,['severity','code','detail'])
        status='BLOCKED' if any(x['severity']=='BLOCKER' for x in qc) else 'REVIEW_REQUIRED'
        summary={'contract':CONTRACT,'technical_status':status,
                 'scientific_acceptance':'NOT_APPROVED','merge_permitted':False,
                 'offline_source_replay':not allow_retries,
                 'sample_targets_n':len(selection),
                 'p06_requests_n':len(requests),
                 'p06_cached_responses_sha256_verified_n':sum(r['cache_verified'] for r in integrity if r['source']=='P06'),
                 'p06_failed_requests_n':len(failures),
                 'chembl_target_mappings_exact_n':len(target_to_gene),
                 'chembl_cached_activity_rows_n':len(activity_evidence),
                 'chembl_cached_exact_numeric_rows_n':sum(x['evidence_status'] in eligible_activity_states for x in activity_evidence),
                 'chembl_cached_assay_verified_examples_n':sum(x['evidence_status']=='CACHED_ASSAY_EXAMPLE_VERIFIED' for x in activity_evidence),
                 'chembl_targets_with_verified_examples_n':sum(x['chembl_assay_verified_examples_n']>0 for x in per_target),
                 'bindingdb_parsed_cached_rows_n':sum(x['raw_rows'] for x in bdb_rows),
                 'bindingdb_name_only_chemically_eligible_n':sum(x['chemically_eligible_rows'] for x in bdb_rows),
                 'bindingdb_name_groups_n':len(bdb_rows),
                 'bindingdb_success_empty_response_targets_n':len(empty_bindingdb_responses),
                 'bindingdb_rows_verified_target_binding_n':0,
                 'retriable_original_chembl_requests_n':len(plans),
                 'retries_approved':allow_retries,'retried_requests_n':len(retry_rows),
                 'retried_success_n':sum(x['result']=='OK' for x in retry_rows),
                 'p0_scientific_baseline_panel_genes':p0_summary.get('panel_genes_n'),
                 'p0_scientific_baseline_pyz':p0_summary.get('own_compounds_n'),
                 'p0_scientific_baseline_matrix_rows':p0_summary.get('existing_matrix_rows_n'),
                 'target_pyz_binding_established':False,
                 'source_archives_sha256':{p0zip.name:hash_path(p0zip),p06zip.name:hash_path(p06zip)},
                 'scientific_source_unchanged':True}
        (out/'p07_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        (out/'p07_source_manifest.json').write_text(json.dumps({
            'contract':CONTRACT,'input_p0_file':str(p0zip),'input_p06_file':str(p06zip),
            'input_sha256':summary['source_archives_sha256'],
            'p0_internal_sha256_verified':{n:digest(za.read(n)) for n in P0_INTEGRITY},
            'source_paths_are_readonly':True,'network_opt_in':allow_retries,
            'source_urls':['https://www.ebi.ac.uk/chembl/api/data/','https://www.bindingdb.org/rest/getLigandsByUniprots'],
            'cache_entry_limit_bytes':MAX_MEMBER,'total_zip_uncompressed_limit_bytes':MAX_ZIP_UNCOMPRESSED},
            ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        report=f'''# MCL TLS V3 / P0.7 — научная проверка экспериментальных источников\n\n**Статус:** {status}. **Научная приёмка:** NOT_APPROVED. **Слияние:** запрещено.\n\n- Исходная выборка: **{len(selection)}** белков. Исходных запросов P0.6: **{len(requests)}**; ответов, проверенных по SHA256: **{summary['p06_cached_responses_sha256_verified_n']}**; ошибок: **{len(failures)}**.\n- Проверенные соответствия ChEMBL SINGLE PROTEIN по UniProt: **{len(target_to_gene)}** идентификаторов (локальный кэш, не полный реестр).\n- Строк активности ChEMBL из кэшированных ответов: **{len(activity_evidence)}** (из них **{summary['chembl_cached_exact_numeric_rows_n']}** с точными положительными числами и первично допустимым типом анализа); примеров со сверенной записью assay и confidence_score=9: **{summary['chembl_cached_assay_verified_examples_n']}**. Не полное покрытие ChEMBL.\n- BindingDB: **{summary['bindingdb_parsed_cached_rows_n']}** прочитанных строк; **{summary['bindingdb_name_only_chemically_eligible_n']}** предварительно подходящих химических наблюдений, но **0 подтверждённых прямых измерений для конкретного UniProt**. Смешанные названия белков в одном запросе требуют ручной проверки.\n- Запросов ChEMBL, пригодных для контролируемого повтора: **{len(plans)}**. Сетевой режим: **{'явно разрешён' if allow_retries else 'выключен'}**; попыток: **{len(retry_rows)}**, успешных: **{summary['retried_success_n']}**.\n\n## Научные ограничения\n\nКэш — лишь часть экспериментов из ChEMBL и BindingDB, к тому же с неполными ответами. Исходная выборка не репрезентативна для всего протеома. Наличие лигандов известного вещества не доказывает связывание PYZ. Для BindingDB совпадение названия белка с UniProt не равняется экспериментально подтверждённой идентичности белкового объекта.\n\n## Решение\n\nНе переносить P0.7 в рабочую матрицу PYZ и не утверждать панель P1 без разрешения конфликтов идентификаторов и научного QC экспериментальных записей. Исходные 65×13 не изменяются.\n'''
        (out/'P07_REPORT_RU.md').write_text(report,encoding='utf-8')
        # Integrity check again after all calculations.
        if digest(zb.read('p06_summary.json'))!=digest(json.dumps(p06_summary,ensure_ascii=False,indent=2).encode('utf-8')+b'\n'):
            # Compare original file rather than guessed serialization if producer differed.
            if hash_path(p06zip)!=summary['source_archives_sha256'][p06zip.name]:
                raise RuntimeError('Source P0.6 changed during audit')
        for path in (p0zip,p06zip):
            if hash_path(path)!=summary['source_archives_sha256'][path.name]:
                raise RuntimeError('Source changed during audit: '+str(path))
        print(json.dumps({k:summary[k] for k in ('technical_status','sample_targets_n','p06_failed_requests_n',
            'chembl_cached_assay_verified_examples_n','bindingdb_name_only_chemically_eligible_n',
            'retries_approved')},ensure_ascii=False),flush=True)
        return summary
    finally:
        za.close();zb.close()


def main():
    ap=argparse.ArgumentParser(description='MCL TLS V3/P0.7 read-only replay and bounded ChEMBL source verification')
    ap.add_argument('--p0-zip',type=Path,required=True)
    ap.add_argument('--p06-zip',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--enable-targeted-retry',action='store_true')
    ap.add_argument('--max-requests',type=int,default=12)
    opts=ap.parse_args()
    main_audit(opts.p0_zip,opts.p06_zip,opts.out,allow_retries=opts.enable_targeted_retry,
               max_retries=opts.max_requests)

if __name__=='__main__':main()