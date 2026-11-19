#!/usr/bin/env python3
"""TLS V3/P0.6: bounded, reproducible, evidence-conscious target-space pilot.

OFFLINE is default. Only --enable-external-pilot permits HTTP requests to exact
allowlisted public UniProt, ChEMBL, and BindingDB REST hosts. No science writes,
no PYZ scoring, no inference of target binding, no automatic merge/acceptance.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
import io
import json
from pathlib import Path
import re
import time
from typing import Any
from urllib.parse import urlencode, urlparse

import httpx

CONTRACT = 'mcl-tls-v3-p06-target-coverage-pilot-v1'
UNIPROT_URL = 'https://rest.uniprot.org/uniprotkb/stream'
CHEMBL_TARGET_URL = 'https://www.ebi.ac.uk/chembl/api/data/target.json'
CHEMBL_ACTIVITY_URL = 'https://www.ebi.ac.uk/chembl/api/data/activity.json'
CHEMBL_ASSAY = 'https://www.ebi.ac.uk/chembl/api/data/assay'
BINDINGDB_URL = 'https://www.bindingdb.org/rest/getLigandsByUniprots'
HOSTS = {'rest.uniprot.org', 'www.ebi.ac.uk', 'www.bindingdb.org'}
ENDPOINTS = ('Ki', 'Kd', 'IC50')
UNIPROT_FIELDS = 'accession,gene_primary,xref_hgnc,protein_name'
UNIPROT_QUERY = 'organism_id:9606 AND reviewed:true'
SOURCE_PLACEHOLDER_RE = re.compile(r'^UNNAMED:\s*\d+$', re.I)
ACCESSION_RE = re.compile(r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z][A-Z0-9]{2}[0-9])?)(?:-[0-9]+)?$')
HGNC_RE = re.compile(r'HGNC:\d+')
QUOTAS = [('CONTROL_13', 13), ('NONPRIORITY_MAPPED', 12),
          ('MAPPED_OTHER', 7), ('CRISPR_UNMAPPED', 8)]
MAX_TARGETS_ALLOWED = 50


def strval(value: Any) -> str:
    return str(value or '').strip()


def istrue(value: Any) -> bool:
    return str(value).strip().lower() in {'true', '1', 'yes'}


def sha_file(path: Path) -> str:
    h = sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open('r', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames:
            raise ValueError(f'Empty TSV header: {path}')
        return list(reader)


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def deterministic_order(rows: list[dict], seed: str) -> list[dict]:
    return sorted(rows, key=lambda r: (sha256(f'{seed}:{r["gene_symbol"]}'.encode()).hexdigest(), r['gene_symbol']))


def stratum(row: dict[str, str]) -> str:
    if istrue(row.get('source_artifact')) or SOURCE_PLACEHOLDER_RE.fullmatch(strval(row.get('gene_symbol'))):
        return 'EXCLUDED_SOURCE_ARTIFACT'
    if istrue(row.get('path_C')):
        return 'CONTROL_13'
    if istrue(row.get('candidate_v2_any_status')) and not istrue(row.get('candidate_v2_priority')) and istrue(row.get('protein_unique_reviewed_human')):
        return 'NONPRIORITY_MAPPED'
    if istrue(row.get('protein_unique_reviewed_human')):
        return 'MAPPED_OTHER'
    if istrue(row.get('has_crispr')):
        return 'CRISPR_UNMAPPED'
    return 'OTHER'


def sample_panel(rows: list[dict[str, str]], *, max_targets: int = 40, seed: str = 'tls-v3-p06-v1') -> list[dict]:
    if not 1 <= max_targets <= MAX_TARGETS_ALLOWED:
        raise ValueError(f'max_targets must be 1..{MAX_TARGETS_ALLOWED}')
    symbols = [strval(r.get('gene_symbol')).upper() for r in rows]
    if len(set(symbols)) != len(symbols) or any(not v for v in symbols):
        raise ValueError('P0 target_universe has duplicate or blank gene symbols')
    candidates = []
    for row in rows:
        g = strval(row['gene_symbol']).upper()
        if SOURCE_PLACEHOLDER_RE.fullmatch(g) or istrue(row.get('source_artifact')):
            continue
        r = dict(row)
        r['gene_symbol'] = g
        r['sample_stratum'] = stratum(r)
        candidates.append(r)
    buckets = {name: deterministic_order([r for r in candidates if r['sample_stratum'] == name], seed)
               for name in [x[0] for x in QUOTAS] + ['OTHER']}
    chosen = []
    seen = set()
    for name, quota in QUOTAS:
        for r in buckets[name][:min(quota, max_targets - len(chosen))]:
            chosen.append(r)
            seen.add(r['gene_symbol'])
    remaining = deterministic_order([r for r in candidates if r['gene_symbol'] not in seen], seed + ':fill')
    chosen += remaining[:max(0, max_targets - len(chosen))]
    return chosen


def parse_uniprot_tsv(data: bytes) -> tuple[list[dict], dict[str, list[dict]], list[dict]]:
    """UniProt reviewed human snapshot, strict schema, exact gene symbol join.

    Full reviewed human UniProtKB entry count != unique genes or drug targets.
    """
    text = data.decode('utf-8-sig')
    reader = csv.DictReader(io.StringIO(text), delimiter='\t')
    columns = reader.fieldnames or []
    def col(match: str) -> str:
        options = [c for c in columns if c.lower().strip() == match.lower()]
        if len(options) != 1:
            raise ValueError(f'UniProt TSV schema changed: required {match}; got {columns}')
        return options[0]
    entry, symbol, hgnc = col('Entry'), col('Gene Names (primary)'), col('HGNC')
    all_rows: list[dict] = []
    by_gene: dict[str, list[dict]] = {}
    issues: list[dict] = []
    accessions: set[str] = set()
    for index, row in enumerate(reader, start=2):
        accession = strval(row.get(entry)).upper()
        raw_gene = strval(row.get(symbol)).upper()
        ids = sorted(set(HGNC_RE.findall(strval(row.get(hgnc)).upper())))
        if not ACCESSION_RE.fullmatch(accession) or accession in accessions:
            raise ValueError(f'UniProt duplicate/invalid accession at row {index}: {accession}')
        accessions.add(accession)
        item = {'uniprot_accession': accession, 'gene_symbol': raw_gene,
                'hgnc_ids': ';'.join(ids), 'reviewed_human_from_query': True}
        all_rows.append(item)
        if not raw_gene or ' ' in raw_gene:
            issues.append({'gene_symbol': raw_gene, 'source': 'UniProt',
                           'issue': 'NO_UNAMBIGUOUS_PRIMARY_GENE', 'detail': accession})
        else:
            by_gene.setdefault(raw_gene, []).append(item)
        if len(ids) > 1:
            issues.append({'gene_symbol': raw_gene, 'source': 'UniProt',
                           'issue': 'MULTIPLE_HGNC_IDS', 'detail': accession})
    if not 15000 <= len(all_rows) <= 30000:
        raise ValueError(f'Unexpected reviewed-human snapshot size {len(all_rows)}; reject partial API response')
    return all_rows, by_gene, issues


def resolve_gene(gene: str, by_gene: dict[str, list[dict]], local_accession: str) -> dict:
    matches = by_gene.get(gene, [])
    if len(matches) != 1:
        return {'uniprot_accession': '', 'hgnc_id': '',
                'mapping_status': 'NO_EXACT_PRIMARY_GENE' if not matches else 'AMBIGUOUS_PRIMARY_GENE',
                'local_registry_match': 'NOT_ASSESSED'}
    item = matches[0]
    ids = item['hgnc_ids'].split(';') if item['hgnc_ids'] else []
    hgnc = ids[0] if len(ids) == 1 else ''
    acc = item['uniprot_accession']
    local_status = ('MATCH' if acc == local_accession else 'CONFLICT'
                    if local_accession else 'NOT_IN_LOCAL_REGISTRY')
    status = 'EXACT_REVIEWED_WITH_HGNC' if hgnc else 'EXACT_REVIEWED_HGNC_UNRESOLVED'
    if local_status == 'CONFLICT':
        status = 'LOCAL_UNIPROT_CONFLICT'
    return {'uniprot_accession':acc, 'hgnc_id':hgnc,'mapping_status':status,
            'local_registry_match':local_status}


@dataclass
class Limits:
    max_requests: int = 320
    max_total_bytes: int = 50_000_000
    max_single_bytes: int = 5_000_000
    min_delay_s: float = 0.30
    timeout_s: float = 25.0


class ControlledClient:
    def __init__(self, out: Path, limits: Limits | None = None, transport: httpx.BaseTransport | None = None):
        self.limits = limits or Limits()
        self.out = out
        self.out.mkdir(parents=True, exist_ok=True)
        self.count = 0
        self.total_bytes = 0
        self.last_request: float | None = None
        self.records: list[dict] = []
        self.client = httpx.Client(follow_redirects=False, timeout=self.limits.timeout_s,
                                   transport=transport, headers={'User-Agent':'MCL-TLS-V3-P06-research-inventory/1.0',
                                                                 'Accept':'application/json, text/tab-separated-values'})

    def close(self):
        self.client.close()

    def get(self, url: str, params: dict[str, str | int], *, max_bytes: int | None = None) -> tuple[bytes, dict]:
        u = urlparse(url)
        if u.scheme != 'https' or u.hostname not in HOSTS or u.username or u.password or u.port not in (None, 443):
            raise ValueError('Blocked non-allowlisted external request')
        key = url + '?' + urlencode(sorted((str(k), str(v)) for k, v in params.items()))
        keyhash = sha256(key.encode()).hexdigest()
        budget = min(self.limits.max_single_bytes, max_bytes or self.limits.max_single_bytes)
        cache = self.out / f'{keyhash}.bin'
        if cache.exists():
            body = cache.read_bytes()
            if len(body) > budget: raise ValueError('Cached response exceeds size cap')
            return body, {'cache':'REUSED','sha256':sha256(body).hexdigest(), 'headers':{}}
        if self.count >= self.limits.max_requests:
            raise RuntimeError('REQUEST_BUDGET_EXHAUSTED')
        if self.last_request is not None:
            delay = self.limits.min_delay_s - (time.monotonic() - self.last_request)
            if delay > 0: time.sleep(delay)
        self.last_request = time.monotonic()
        self.count += 1
        ts = datetime.now(timezone.utc).isoformat()
        record = {'url':key,'requested_at_utc':ts, 'http_status':'', 'response_bytes':'',
                  'response_sha256':'', 'result':'ERROR', 'cache_file':cache.name,
                  'headers_date':'','headers_etag':'','headers_uniprot_release':'','headers_total_results':''}
        try:
            # Read only a bounded response, never pass through redirects from source data.
            with self.client.stream('GET', url, params=params) as response:
                record['http_status'] = str(response.status_code)
                for name, field in [('date','headers_date'),('etag','headers_etag'),
                                    ('x-uniprot-release','headers_uniprot_release'),('x-total-results','headers_total_results')]:
                    record[field] = response.headers.get(name,'')
                response.raise_for_status()
                length_header = response.headers.get('content-length')
                if length_header and int(length_header) > budget:
                    raise RuntimeError('SINGLE_RESPONSE_BUDGET_EXCEEDED')
                chunks=[]; n=0
                for chunk in response.iter_bytes():
                    n += len(chunk)
                    if n > budget or self.total_bytes + n > self.limits.max_total_bytes:
                        raise RuntimeError('RESPONSE_BYTE_BUDGET_EXCEEDED')
                    chunks.append(chunk)
                body = b''.join(chunks)
            self.total_bytes += len(body)
            cache.write_bytes(body)
            record.update(result='OK',response_bytes=str(len(body)),
                          response_sha256=sha256(body).hexdigest())
            self.records.append(record)
            return body, {'cache':'NEW','sha256':record['response_sha256'],
                          'headers':{k:record[k] for k in ('headers_date','headers_etag','headers_uniprot_release','headers_total_results')}}
        except (httpx.HTTPError, RuntimeError, ValueError) as e:
            record['result'] = f'{type(e).__name__}:{str(e)[:180]}'
            self.records.append(record)
            raise

    def json(self, url: str, params: dict[str,str|int], *, max_bytes: int | None = None) -> dict:
        body, _ = self.get(url, params, max_bytes=max_bytes)
        parsed = json.loads(body)
        if not isinstance(parsed, dict): raise ValueError('Non-object API JSON response')
        return parsed


def exact_chembl_target(doc: dict, accession: str) -> tuple[str, str]:
    rows = doc.get('targets')
    if not isinstance(rows, list): return '', 'MALFORMED_TARGET_RESPONSE'
    good = []
    for r in rows:
        if not isinstance(r,dict) or strval(r.get('target_type')).upper() != 'SINGLE PROTEIN':
            continue
        if strval(r.get('organism')).lower() != 'homo sapiens' and strval(r.get('tax_id')) != '9606':
            continue
        components = r.get('target_components') or []
        accessions = {strval(c.get('accession')).upper() for c in components if isinstance(c,dict)}
        if accessions == {accession} and strval(r.get('target_chembl_id')):
            good.append(strval(r['target_chembl_id']))
    good = sorted(set(good))
    if len(good) == 1: return good[0], 'EXACT_SINGLE_PROTEIN_HUMAN'
    return '', 'NO_EXACT_HUMAN_SINGLE_PROTEIN' if not good else 'AMBIGUOUS_CHEMBL_TARGET'


def valid_chembl_activity(row: dict, target_id: str, endpoint: str) -> bool:
    if strval(row.get('target_chembl_id')) != target_id:
        return False
    if strval(row.get('standard_type')) != endpoint or strval(row.get('standard_units')) != 'nM':
        return False
    if strval(row.get('standard_relation')) != '=':
        return False
    if strval(row.get('assay_type')).upper() != 'B':
        return False
    try:
        v = float(row.get('standard_value'))
        return 0 < v < float('inf') and bool(strval(row.get('assay_chembl_id'))) and bool(strval(row.get('molecule_chembl_id')))
    except (ValueError, TypeError): return False


def read_binding_records(payload: dict) -> list[dict]:
    """Collect candidate objects; unknown/empty structure => UNVERIFIED, not zero."""
    rows: list[dict] = []
    def field(r: dict, suffix: str):
        for key, value in r.items():
            name = re.sub(r'[^a-z0-9]', '', key.lower())
            if name.endswith(suffix): return value
        return None
    def walk(x: Any, depth: int = 0):
        if depth > 15: return
        if isinstance(x,dict):
            if field(x,'affinitytype') is not None and field(x,'affinity') is not None and (
                field(x,'smile') is not None or field(x,'smiles') is not None):
                rows.append(x)
            for v in x.values(): walk(v,depth+1)
        elif isinstance(x,list):
            for v in x:walk(v,depth+1)
    walk(payload)
    return rows


def binding_quality(payload: dict, accession: str) -> dict:
    records=read_binding_records(payload)
    if not records:
        return {'bindingdb_status':'NO_PARSEABLE_RECORDS_OR_EMPTY','bindingdb_records_parsed':0,
                'bindingdb_candidate_records':0, 'bindingdb_experiment_verified':False}
    candidates=0
    for r in records:
        fields={re.sub('[^a-z0-9]','',k.lower()):strval(v) for k,v in r.items()}
        def val(suffix):
            return next((v for k,v in fields.items() if k.endswith(suffix)), '')
        endpoint=val('affinitytype').upper()
        if endpoint not in {'KI','KD','IC50'}: continue
        raw=val('affinity').replace(',','')
        m=re.fullmatch(r'\s*(<=|>=|<|>|=)?\s*(\d+(?:\.\d+)?)\s*',raw)
        if not m or not (0 < float(m[2]) <= 10000): continue
        if m[1] not in (None,'=','<','<='):continue
        declared=val('uniprot') or val('query')
        seen=set(re.findall(r'\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])\b',declared.upper()))
        if not seen or seen!={accession}: continue # fail closed without exact reported target
        if val('smile') or val('smiles'): candidates+=1
    return {'bindingdb_status':'CANDIDATE_ROWS_NEED_SOURCE_ASSAY_REVIEW' if candidates else 'NO_STRICT_CANDIDATES_IN_RESPONSE',
            'bindingdb_records_parsed':len(records), 'bindingdb_candidate_records':candidates,
            'bindingdb_experiment_verified':False}


def probe_one(client: ControlledClient, gene: str, accession: str) -> dict:
    out={'chembl_status':'NOT_QUERIED','chembl_target_id':'','chembl_activity_candidates_total':'',
         'chembl_exact_assay_example_ids':'','bindingdb_status':'NOT_QUERIED',
         'bindingdb_records_parsed':'','bindingdb_candidate_records':'',
         'bindingdb_experiment_verified':False}
    if not accession:
        out['chembl_status']=out['bindingdb_status']='MAPPING_NOT_RESOLVED'
        return out
    try:
        doc=client.json(CHEMBL_TARGET_URL,{'target_components__accession':accession,'limit':100},max_bytes=400_000)
        target,status=exact_chembl_target(doc, accession)
        out['chembl_status']=status;out['chembl_target_id']=target
        if target:
            total=0;examples=[];queried=0
            for typ in ENDPOINTS:
                activities=client.json(CHEMBL_ACTIVITY_URL,
                    {'target_chembl_id':target,'assay_type':'B','standard_type':typ,
                     'standard_units':'nM','limit':10,'offset':0},max_bytes=400_000)
                rows=activities.get('activities');meta=activities.get('page_meta')
                if not isinstance(rows,list) or not isinstance(meta,dict):
                    raise ValueError('Unexpected ChEMBL activity page schema')
                n=meta.get('total_count')
                if not isinstance(n,int) or n<0:
                    raise ValueError('No ChEMBL total_count')
                total+=n;queried+=1
                examples.extend([x for x in rows if isinstance(x,dict) and valid_chembl_activity(x,target,typ)])
            out['chembl_activity_candidates_total']=total
            # Verify up to two distinct assay IDs, never treat a mere result count
            # as experimental binding confirmed.
            seen_assays=set()
            verified=[]
            for row in examples:
                aid=strval(row.get('assay_chembl_id'))
                if not re.fullmatch(r'CHEMBL[0-9]{1,12}', aid):
                    continue  # never interpolate arbitrary provider text in URL paths
                if aid in seen_assays or len(seen_assays)>=2:continue
                seen_assays.add(aid)
                assay=client.json(f'{CHEMBL_ASSAY}/{aid}.json',{},max_bytes=100_000)
                if strval(assay.get('assay_type')).upper() == 'B' and \
                    strval(assay.get('target_chembl_id')) == target and \
                    strval(assay.get('confidence_score')) == '9':
                    verified.append(f'{strval(row.get("activity_id"))}:{aid}:{strval(row.get("molecule_chembl_id"))}')
            out['chembl_exact_assay_example_ids']=';'.join(verified)
            out['chembl_status']='EXAMPLE_DIRECT_ASSAY_VERIFIED' if verified else \
                'CANDIDATES_EXIST_NOT_INDEPENDENTLY_VERIFIED' if total else 'ZERO_IN_FILTERED_QUERY'
    except (httpx.HTTPError, ValueError, RuntimeError, json.JSONDecodeError) as e:
        out['chembl_status']='ERROR_'+type(e).__name__
    try:
        binding=client.json(BINDINGDB_URL,
            {'uniprot':accession,'cutoff':10000,'response':'application/json'},max_bytes=2_000_000)
        out.update(binding_quality(binding,accession))
    except (httpx.HTTPError, ValueError, RuntimeError, json.JSONDecodeError) as e:
        out['bindingdb_status']='ERROR_'+type(e).__name__
    return out


def run(p0: Path, out: Path, *, enable_external_pilot: bool=False,
        max_targets: int=40, seed: str='tls-v3-p06-v1', transport: httpx.BaseTransport | None=None,
        limits: Limits | None=None) -> dict:
    p0=p0.resolve();out=out.resolve()
    if not p0.is_dir() or not (p0/'target_universe.tsv').is_file() or not (p0/'coverage_summary.json').is_file():
        raise ValueError('Existing P0 report folder required; expected target_universe.tsv and coverage_summary.json')
    if out==p0 or p0 in out.parents or out in p0.parents or out.exists():
        raise ValueError('P0.6 output must be a new independent folder, not the P0 report')
    meta=json.loads((p0/'coverage_summary.json').read_text(encoding='utf-8-sig'))
    if meta.get('contract') != 'mcl-tls-v3-p0-inventory-v1':
        raise ValueError('P0 report schema unsupported')
    if meta.get('technical_status')=='BLOCKED':
        raise ValueError('Cannot pilot a technically BLOCKED P0 snapshot')
    # Require the corrected P0 audit: old T4 is not evidence-linked and old
    # placeholder counts cannot be accepted as a scientific universe.
    if 'valid_gene_candidates_n' not in meta or 'suspect_source_placeholders_n' not in meta:
        raise ValueError('P0 audit is older than source-artifact/T4 correction: rerun current TLS V3 P0 first')
    input_files=[p0/'target_universe.tsv',p0/'coverage_summary.json',p0/'source_manifest.json']
    if any(not x.is_file() for x in input_files):raise ValueError('Required P0 provenance file missing')
    before={p.name:sha_file(p) for p in input_files}
    entries=read_tsv(input_files[0]);selection=sample_panel(entries,max_targets=max_targets,seed=seed)
    out.mkdir(parents=True,exist_ok=False)
    sample_fields=['gene_symbol','sample_stratum','uniprot_accession','protein_unique_reviewed_human',
                   'has_crispr','candidate_v2_any_status','candidate_v2_priority','path_C']
    write_tsv(out/'p06_selection.tsv',selection,sample_fields)
    problems=[]
    external_rows=[]
    snapshot_rows=[]
    request_rows=[]
    snapshot_meta={'status':'NOT_REQUESTED_OFFLINE','entries_n':None,'sha256':None,
                   'release_header':'','source_url':UNIPROT_URL}
    by_gene:dict[str,list[dict]]={}
    client=None
    if enable_external_pilot:
        client=ControlledClient(out/'external_response_cache',limits,transport)
        try:
            print('[P0.6] Explicit external pilot enabled: fetching one bounded reviewed-human UniProt TSV snapshot',flush=True)
            body, info=client.get(UNIPROT_URL,
                {'query':UNIPROT_QUERY,'fields':UNIPROT_FIELDS,'format':'tsv'},max_bytes=5_000_000)
            snapshot_rows,by_gene,issues=parse_uniprot_tsv(body)
            problems.extend(issues)
            count_header=info['headers'].get('headers_total_results','')
            if count_header and (not count_header.isdecimal() or int(count_header)!=len(snapshot_rows)):
                raise ValueError(f'UniProt response count mismatch: header={count_header}; rows={len(snapshot_rows)}')
            snapshot_meta.update(status='RETRIEVED_REVIEWED_HUMAN', entries_n=len(snapshot_rows),
                                 sha256=info['sha256'], release_header=info['headers'].get('headers_uniprot_release',''),
                                 total_count_header=count_header,
                                 completeness_status='COUNT_HEADER_MATCH' if count_header else 'COUNT_HEADER_UNAVAILABLE',
                                 retrieved_at_utc=datetime.now(timezone.utc).isoformat())
            write_tsv(out/'uniprot_reviewed_human.tsv',snapshot_rows,
                      ['uniprot_accession','gene_symbol','hgnc_ids','reviewed_human_from_query'])
        except Exception as e:
            snapshot_meta['status']='ERROR_'+type(e).__name__
            snapshot_meta['error']=str(e)[:200]
            problems.append({'gene_symbol':'','source':'UniProt', 'issue':'PROTEOME_SOURCE_UNAVAILABLE','detail':str(e)[:250]})
    # Build a whole-MCL candidate mapping from the single reviewed-human snapshot.
    # These are gene-to-reviewed-protein crosswalks, not binding-target claims.
    mapped_all=[]
    if snapshot_meta['status']=='RETRIEVED_REVIEWED_HUMAN':
        for row in entries:
            symbol=strval(row['gene_symbol']).upper()
            if SOURCE_PLACEHOLDER_RE.fullmatch(symbol) or istrue(row.get('source_artifact')):
                mapped_all.append({'gene_symbol':symbol,'mapping_status':'EXCLUDED_SOURCE_PLACEHOLDER',
                                   'uniprot_accession':'','hgnc_id':'','local_registry_match':'NOT_ASSESSED'})
                continue
            mapped_all.append({'gene_symbol':symbol, **resolve_gene(symbol,by_gene,strval(row.get('uniprot_accession')).upper())})
        write_tsv(out/'mcl_reviewed_uniprot_mapping.tsv',mapped_all,
                  ['gene_symbol','uniprot_accession','hgnc_id','mapping_status','local_registry_match'])
    print(f'[P0.6] Deterministic panel: {len(selection)} entries across 5 strata; internet={enable_external_pilot}',flush=True)
    try:
        for i,row in enumerate(selection,1):
            g=row['gene_symbol']
            local=strval(row.get('uniprot_accession')).upper()
            ident=resolve_gene(g,by_gene,local) if snapshot_meta['status']=='RETRIEVED_REVIEWED_HUMAN' else {
                'uniprot_accession':'','hgnc_id':'','mapping_status':'NOT_ASSESSED_UNIPROT',
                'local_registry_match':'NOT_ASSESSED'}
            if ident['mapping_status'] in {'AMBIGUOUS_PRIMARY_GENE','LOCAL_UNIPROT_CONFLICT','EXACT_REVIEWED_HGNC_UNRESOLVED'}:
                problems.append({'gene_symbol':g,'source':'UniProt','issue':ident['mapping_status'],
                                 'detail':local+' -> '+ident['uniprot_accession']})
            if client and snapshot_meta['status']=='RETRIEVED_REVIEWED_HUMAN':
                # Conflicting or HGNC-less mapping stays a manual review,
                # never sent to target-ligand endpoints as validated identity.
                acc=ident['uniprot_accession'] if ident['mapping_status']=='EXACT_REVIEWED_WITH_HGNC' else ''
                print(f'[P0.6] {i}/{len(selection)} {g} UniProt={acc or "UNRESOLVED"}',flush=True)
                probe=probe_one(client,g,acc)
            else:
                probe={'chembl_status':'OFFLINE_NOT_QUERIED','bindingdb_status':'OFFLINE_NOT_QUERIED'}
            external_rows.append({'gene_symbol':g,'sample_stratum':row['sample_stratum'],
                                  'local_uniprot_accession':local,**ident,**probe})
    finally:
        if client:
            request_rows=client.records
            client.close()
    write_tsv(out/'p06_external_coverage.tsv',external_rows,
              ['gene_symbol','sample_stratum','local_uniprot_accession',
               'uniprot_accession','hgnc_id','mapping_status','local_registry_match',
               'chembl_status','chembl_target_id','chembl_activity_candidates_total',
               'chembl_exact_assay_example_ids','bindingdb_status','bindingdb_records_parsed',
               'bindingdb_candidate_records','bindingdb_experiment_verified'])
    write_tsv(out/'p06_mapping_issues.tsv',problems,['gene_symbol','source','issue','detail'])
    write_tsv(out/'p06_requests.tsv',request_rows,
              ['url','requested_at_utc','http_status','response_bytes','response_sha256','result',
               'cache_file','headers_date','headers_etag','headers_uniprot_release','headers_total_results'])
    after={p.name:sha_file(p) for p in input_files}
    if before!=after:
        raise RuntimeError('SOURCE_CHANGED_DURING_PILOT')
    strata={name:sum(r['sample_stratum']==name for r in selection)
            for name in [x[0] for x in QUOTAS]+['OTHER']}
    verified=[r for r in external_rows if r.get('chembl_status')=='EXAMPLE_DIRECT_ASSAY_VERIFIED']
    source_errors=[r for r in external_rows if strval(r.get('chembl_status')).startswith('ERROR_')
                   or strval(r.get('bindingdb_status')).startswith('ERROR_')]
    summary={
        'contract':CONTRACT,'technical_status':'REVIEW_REQUIRED',
        'scientific_acceptance':'NOT_APPROVED','merge_permitted':False,
        'sample_size':len(selection),'sample_seed':seed,'strata_counts':strata,
        'external_mode':'EXPLICIT_BOUNDED_PILOT' if enable_external_pilot else 'OFFLINE_ONLY',
        'uniprot_snapshot':snapshot_meta,
        'mcl_candidate_gene_rows_reviewed_n':len(mapped_all) if mapped_all else None,
        'mcl_exact_reviewed_with_hgnc_n':sum(x['mapping_status']=='EXACT_REVIEWED_WITH_HGNC' for x in mapped_all) if mapped_all else None,
        'mcl_uniprot_local_mapping_conflicts_n':sum(x['mapping_status']=='LOCAL_UNIPROT_CONFLICT' for x in mapped_all) if mapped_all else None,
        'sample_chembl_exact_assay_examples_n':len(verified),
        'sample_source_error_targets_n':len(source_errors),
        'requests_n':len(request_rows),'sample_not_proteome_estimate':True,
        'old_panel_preserved':True,'source_p0_sha256_before_after':before,
        'target_pyz_binding_established':False,
        'sampling_limitations':['Stratified small convenience panel; no unbiased estimate of entire proteome',
         'ChEMBL queried only first 10 candidate rows per endpoint and up to two assays per target',
         'BindingDB response capped to 2 MB per target, and source assay directness remains unverified',
         'UniProt primary symbol lookup does not substitute full HGNC ID mapping for all MCL genes',
         'Raw response cache is stored in output and never used to mutate scientific data']}
    (out/'p06_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    (out/'p06_source_manifest.json').write_text(json.dumps({'input_dir':str(p0), 'input_sha256':before,
      'runtime_snapshot':snapshot_meta,'source_urls':[UNIPROT_URL,CHEMBL_TARGET_URL,CHEMBL_ACTIVITY_URL,BINDINGDB_URL],
      'network_allowlist':sorted(HOSTS),'request_count':len(request_rows),'read_only_scientific_data':True},
      ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    qc=[]
    qc.append({'severity':'REVIEW','code':'SAMPLE_NOT_POPULATION','detail':'No extrapolation to full protein universe'})
    if snapshot_meta['status']!='RETRIEVED_REVIEWED_HUMAN':
        qc.append({'severity':'REVIEW','code':'UNIPROT_FULL_MAPPING_NOT_ASSESSED','detail':snapshot_meta['status']})
    if enable_external_pilot and snapshot_meta.get('completeness_status') == 'COUNT_HEADER_UNAVAILABLE':
        qc.append({'severity':'REVIEW','code':'UNIPROT_SNAPSHOT_TOTAL_NOT_INDEPENDENTLY_VERIFIED',
                   'detail':'UniProt response does not expose X-Total-Results; reviewed-human snapshot may need separate count validation'})
    if source_errors:
        qc.append({'severity':'REVIEW','code':'REMOTE_SOURCE_ERRORS','detail':str(len(source_errors))})
    qc.append({'severity':'REVIEW','code':'BINDINGDB_DIRECTNESS_NOT_VERIFIED','detail':'Candidate rows only; no full assay verification'})
    write_tsv(out/'p06_QC.tsv',qc,['severity','code','detail'])
    text=f'''# TLS V3 / P0.6 — проверка расширяемости пространства мишеней\n\n**Режим:** {summary['external_mode']}. **Научная приёмка:** NOT_APPROVED. **Слияние:** запрещено.\n\n- Проверено в стратифицированной выборке: **{len(selection)}** генов; структура: {strata}.\n- Снимок reviewed human UniProt: **{snapshot_meta['status']}**; записей: **{snapshot_meta.get('entries_n')}**.\n- Точные reviewed human UniProt + HGNC сопоставления со всем генным каталогом MCL: **{summary.get('mcl_exact_reviewed_with_hgnc_n')}**; конфликтов с локальным UniProt: **{summary.get('mcl_uniprot_local_mapping_conflicts_n')}**. Это идентификация, не доказательство лекарственной пригодности.
- Число примеров химических измерений ChEMBL с отдельно проверенным assay (не полное покрытие): **{len(verified)}** мишеней.\n- Целевых запросов с ошибками удалённых источников: **{len(source_errors)}**.\n- Всего запросов: **{len(request_rows)}**.\n\n## Ограничения\n\nТекущий пилот не создаёт матрицу PYZ, не измеряет связывание и не подтверждает эффективность. 40 отобранных белков НЕ представляют случайную выборку полного протеома. Даже при нулевых результатах источник мог быть неполным или доступным только вне проверенного набора. По BindingDB сохраняются только кандидаты без повышения до подтверждённых прямых измерений. Никакого научно принятого числа новых мишеней для P1 нет.\n\nРезультаты первичного P0 сохранены без изменений и проверены по SHA-256.\n'''
    (out/'P06_REPORT_RU.md').write_text(text,encoding='utf-8')
    print(json.dumps({k:summary[k] for k in ('external_mode','sample_size','sample_chembl_exact_assay_examples_n','requests_n','scientific_acceptance')},ensure_ascii=False),flush=True)
    return summary


def main():
    ap=argparse.ArgumentParser(description='TLS V3 / P0.6 explicit opt-in scientific target-space external pilot')
    ap.add_argument('--p0-report',required=True,type=Path)
    ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--max-targets',type=int,default=40)
    ap.add_argument('--seed',default='tls-v3-p06-v1')
    flags=ap.add_mutually_exclusive_group()
    flags.add_argument('--offline',action='store_true')
    flags.add_argument('--enable-external-pilot',action='store_true')
    args=ap.parse_args()
    run(args.p0_report,args.out,enable_external_pilot=args.enable_external_pilot,
        max_targets=args.max_targets,seed=args.seed)

if __name__=='__main__':
    main()