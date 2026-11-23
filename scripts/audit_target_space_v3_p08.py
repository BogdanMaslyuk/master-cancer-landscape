#!/usr/bin/env python3
"""MCL TLS V3/P0.8: read-only local ChEMBL SQLite + BindingDB TSV audit.

Inputs are explicitly supplied version-pinned files. No HTTP, no modification of
scientific datasets, no PYZ activity or causal target assertions.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from contextlib import contextmanager
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sqlite3
import sys
from typing import Any
import zipfile

CONTRACT = 'mcl-tls-v3-p08-local-source-validation-v1'
P0_CONTRACT = 'mcl-tls-v3-p0-inventory-v1'
P06_CONTRACT = 'mcl-tls-v3-p06-target-coverage-pilot-v1'
FOCUS6 = ('PRKDC', 'BCL2L1', 'IGF1R', 'ERBB3', 'FNTB', 'DDC')
ENDPOINTS = ('Ki', 'Kd', 'IC50')
VALUE = re.compile(r'^\s*(<=|>=|<|>|=|≤|≥)?\s*(\d+(?:\.\d*)?|\.\d+)(?:[eE]([+-]?\d+))?\s*$')
ACCESSION = re.compile(r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z][A-Z0-9]{2}[0-9])?)(?:-[0-9]+)?$')
REQUIRED_CHEMBL = {
    'target_dictionary': {'tid', 'chembl_id', 'target_type', 'organism'},
    'target_components': {'tid', 'component_id'},
    'component_sequences': {'component_id', 'accession'},
    'assays': {'assay_id', 'tid', 'assay_type', 'confidence_score', 'chembl_id'},
    'activities': {'activity_id', 'assay_id', 'molregno', 'standard_type', 'standard_units', 'standard_relation', 'standard_value'},
    'molecule_dictionary': {'molregno', 'chembl_id'},
    'compound_structures': {'molregno', 'canonical_smiles'},
}
P0_INTEGRITY = ('target_universe.tsv', 'coverage_summary.json', 'source_manifest.json')
OBS_CHEMBL = ['gene_symbol','uniprot_accession','chembl_target_id','activity_id','assay_chembl_id',
              'molecule_chembl_id','canonical_smiles','standard_type','standard_relation','standard_value',
              'standard_units','assay_type','confidence_score','potential_duplicate','data_validity_comment',
              'document_chembl_id','doi','evidence_status','scientific_acceptance']
OBS_BINDING = ['gene_symbol','uniprot_accession','source_row','bindingdb_reactant_set_id','bindingdb_monomer_id',
               'all_chain_accessions','declared_chain_count','target_name','organism','endpoint',
               'raw_value_nm','relation','numeric_value_nm','ligand_smiles','ligand_inchikey','article_doi',
               'pmid','curation_source','identity_status','measurement_status','scientific_acceptance']


def shafile(path: Path) -> str:
    sha=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b''):sha.update(chunk)
    return sha.hexdigest()


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def t(v: Any) -> str:
    return '' if v is None else str(v).strip()


def as_true(v: Any) -> bool:
    return t(v).lower() in ('1','yes','true')


def read_tsv_bytes(data: bytes) -> list[dict[str,str]]:
    with io.StringIO(data.decode('utf-8-sig'),newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        if not reader.fieldnames or len(reader.fieldnames)!=len(set(reader.fieldnames)):
            raise ValueError('TSV header missing or duplicated')
        return list(reader)


def output_tsv(path: Path,rows: list[dict],fields: list[str]) -> None:
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=fields,extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


class ReadOnlyZip:
    def __init__(self,path: Path):
        self.zip=zipfile.ZipFile(path)
        names=[i.filename for i in self.zip.infolist() if not i.is_dir()]
        if len(names)!=len(set(names)) or len(names)>500:
            self.zip.close();raise ValueError('ZIP duplicate members or too many members')
        for info in self.zip.infolist():
            if info.file_size>100_000_000 or info.filename.startswith('/') or '..' in Path(info.filename).parts:
                self.zip.close();raise ValueError('Unsafe source ZIP member')
        self.names=set(names)
    def read(self,name:str) -> bytes:
        if name not in self.names:
            raise ValueError('Source ZIP missing required member '+name)
        return self.zip.read(name)
    def close(self):self.zip.close()


def parse_nm(raw:Any) -> tuple[str,str,str]:
    m=VALUE.fullmatch(t(raw))
    if not m:return '','','INVALID_OR_MISSING'
    relation={'≤':'<=','≥':'>='}.get(m.group(1),m.group(1) or '=')
    try:num=float(m.group(2)+('e'+m.group(3) if m.group(3) else ''))
    except (ValueError,OverflowError):return '','','INVALID_OR_MISSING'
    if not (math.isfinite(num) and num>0):return '','','INVALID_OR_MISSING'
    return relation,format(num,'.12g'),'EXACT_NUMERIC' if relation=='=' else 'CENSORED'


def rows_from_pilot(p06zip:Path,p0zip:Path|None,p07zip:Path|None, panel:str) -> tuple[list[dict],dict]:
    z=ReadOnlyZip(p06zip)
    try:
        meta=json.loads(z.read('p06_summary.json'))
        if meta.get('contract')!=P06_CONTRACT:raise ValueError('P0.6 contract mismatch')
        selection=read_tsv_bytes(z.read('p06_selection.tsv'))
        coverage=read_tsv_bytes(z.read('p06_external_coverage.tsv'))
        if len(selection)!=meta.get('sample_size') or len(coverage)!=len(selection):
            raise ValueError('P0.6 sample count mismatch')
        if len({x['gene_symbol'] for x in coverage})!=len(coverage):raise ValueError('P0.6 duplicate genes')
        if {x['gene_symbol'] for x in coverage}!={x['gene_symbol'] for x in selection}:
            raise ValueError('P0.6 sample vs coverage gene mismatch')
        if meta.get('uniprot_snapshot',{}).get('status')!='RETRIEVED_REVIEWED_HUMAN':
            raise ValueError('P0.6 UniProt snapshot not validated')
        selected={x['gene_symbol']:x for x in coverage}
        wanted=list(FOCUS6) if panel=='focus6' else [x['gene_symbol'] for x in selection]
        if any(x not in selected for x in wanted):raise ValueError('Requested focus genes missing in P0.6 sample')
        chosen=[]
        for g in wanted:
            row=selected[g]
            acc=t(row.get('uniprot_accession')).upper()
            status=t(row.get('mapping_status'))
            usable=(status=='EXACT_REVIEWED_WITH_HGNC' and bool(ACCESSION.fullmatch(acc)))
            chosen.append({'gene_symbol':g,'uniprot_accession':acc if usable else '',
                           'mapping_status':status,'uniprot_hgnc':t(row.get('hgnc_id')),
                           'sample_stratum':t(row.get('sample_stratum')),
                           'source_identity_status':'VERIFIED_UNIPROT_HGNC' if usable else 'NOT_ASSESSED_MAPPING_UNRESOLVED'})
        linkage={'p06_sha256':shafile(p06zip),'p0_present':bool(p0zip),'p07_present':bool(p07zip)}
        if p0zip:
            p0=ReadOnlyZip(p0zip)
            try:
                p0summary=json.loads(p0.read('coverage_summary.json'))
                if p0summary.get('contract')!=P0_CONTRACT:raise ValueError('P0 contract mismatch')
                for n in P0_INTEGRITY:
                    if sha_bytes(p0.read(n))!=meta.get('source_p0_sha256_before_after',{}).get(n):
                        raise ValueError('P0 and P0.6 provenance mismatch: '+n)
            finally:p0.close()
            linkage['p0_sha256']=shafile(p0zip)
        if p07zip:
            p07=ReadOnlyZip(p07zip)
            try:
                p07meta=json.loads(p07.read('p07_summary.json'))
                if p07meta.get('scientific_acceptance')!='NOT_APPROVED':
                    raise ValueError('Unexpected P0.7 scientific acceptance state')
                recorded=p07meta.get('source_archives_sha256',{})
                if recorded.get(p06zip.name)!=linkage['p06_sha256']:
                    raise ValueError('P0.7 and P0.6 provenance mismatch')
                if p0zip and recorded.get(p0zip.name)!=linkage['p0_sha256']:
                    raise ValueError('P0.7 and P0 provenance mismatch')
            finally:p07.close()
            linkage['p07_sha256']=shafile(p07zip)
        return chosen,linkage
    finally:z.close()


def check_sqlite_schema(db:sqlite3.Connection) -> dict[str,set[str]]:
    columns={}
    for table,required in REQUIRED_CHEMBL.items():
        # Table names originate in fixed internal constants, not user input.
        info=db.execute(f'PRAGMA table_info({table})').fetchall()
        got={r[1] for r in info}
        if not required.issubset(got):
            raise ValueError(f'ChEMBL schema invalid for {table}; missing {sorted(required-got)}')
        columns[table]=got
    for table in ('docs',):
        columns[table]={r[1] for r in db.execute(f'PRAGMA table_info({table})').fetchall()}
    return columns


def chembl_local(path:Path,panel:list[dict],max_rows:int) -> tuple[list[dict],list[dict],dict]:
    if max_rows<1 or max_rows>10000:raise ValueError('max_rows_per_target must be 1..10000')
    # file URI with mode=ro prevents creating or writing any SQLite database.
    db=sqlite3.connect(path.resolve().as_uri()+'?mode=ro',uri=True,timeout=5)
    db.row_factory=sqlite3.Row
    try:
        db.execute('PRAGMA query_only=ON')
        schema=check_sqlite_schema(db)
        extra=schema['activities']
        act_opt=lambda col: f'a.{col}' if col in extra else f'NULL AS {col}'
        doi_col='d.doi' if 'doi' in schema['docs'] else 'NULL AS doi'
        doc_chembl='d.chembl_id AS document_chembl_id' if {'chembl_id','doc_id'}.issubset(schema['docs']) else 'NULL AS document_chembl_id'
        use_docs={'chembl_id','doc_id'}.issubset(schema['docs']) and 'doc_id' in extra
        doc_join='LEFT JOIN docs d ON a.doc_id=d.doc_id' if use_docs else ''
        if not use_docs:doi_col='NULL AS doi';doc_chembl='NULL AS document_chembl_id'
        results=[];coverage=[]
        for entry in panel:
            gene,accession=entry['gene_symbol'],entry['uniprot_accession']
            if not accession:
                coverage.append({'gene_symbol':gene,'status':'NOT_ASSESSED_UNRESOLVED_MAPPING',
                                 'target_matches':0,'candidate_rows':0,'rows_retained':0,'truncated':False})
                continue
            matches=db.execute('''SELECT t.tid,t.chembl_id,t.target_type,t.organism
                FROM target_dictionary t JOIN target_components tc ON tc.tid=t.tid
                JOIN component_sequences cs ON cs.component_id=tc.component_id
                WHERE cs.accession=?''',(accession,)).fetchall()
            exact=[]
            for row in matches:
                if t(row['target_type'])!='SINGLE PROTEIN' or t(row['organism']).lower()!='homo sapiens':continue
                components=db.execute('''SELECT c.accession FROM target_components x
                    JOIN component_sequences c ON x.component_id=c.component_id WHERE x.tid=?''',(row['tid'],)).fetchall()
                if len(components)==1 and t(components[0][0])==accession:exact.append(row)
            if len(exact)!=1:
                coverage.append({'gene_symbol':gene,'status':'NO_EXACT_SINGLE_TARGET' if not exact else 'AMBIGUOUS_EXACT_TARGETS',
                                 'target_matches':len(exact),'candidate_rows':None,'rows_retained':0,'truncated':False})
                continue
            target=exact[0]
            where=""" WHERE s.tid=? AND s.assay_type='B' AND s.confidence_score=9
                AND a.standard_type IN ('Ki','Kd','IC50') AND a.standard_units='nM' """
            total=db.execute('SELECT COUNT(*) FROM activities a JOIN assays s ON a.assay_id=s.assay_id'+where,
                             (target['tid'],)).fetchone()[0]
            query=f'''SELECT a.activity_id, a.standard_type, a.standard_relation,a.standard_value,a.standard_units,
                    {act_opt('potential_duplicate')}, {act_opt('data_validity_comment')},
                    s.chembl_id AS assay_chembl_id, s.assay_type,s.confidence_score,
                    m.chembl_id AS molecule_chembl_id,c.canonical_smiles,
                    {doc_chembl}, {doi_col}
                    FROM activities a JOIN assays s ON a.assay_id=s.assay_id
                    LEFT JOIN molecule_dictionary m ON a.molregno=m.molregno
                    LEFT JOIN compound_structures c ON a.molregno=c.molregno
                    {doc_join} {where}
                    ORDER BY a.activity_id LIMIT ?'''
            chosen=db.execute(query,(target['tid'],max_rows)).fetchall()
            for row in chosen:
                raw=t(row['standard_value']);rel=t(row['standard_relation'])
                try: val=float(raw)
                except (TypeError,ValueError):val=float('nan')
                isexact=rel=='=' and math.isfinite(val) and val>0
                is_valid=not t(row['data_validity_comment']) and t(row['potential_duplicate']) not in ('1','True','true')
                has_structure=bool(t(row['canonical_smiles']))
                if not has_structure:status='REVIEW_MISSING_STRUCTURE'
                elif not math.isfinite(val) or val<=0:status='REVIEW_INVALID_NUMBER'
                elif rel!='=':status='REVIEW_CENSORED_OR_RELATION'
                elif not is_valid:status='REVIEW_CURATOR_FLAGGED_OR_DUPLICATE'
                elif isexact:status='ASSAY_METADATA_MATCHED_PRIMARY_REVIEW_REQUIRED'
                else:status='REVIEW_OTHER'
                results.append({'gene_symbol':gene,'uniprot_accession':accession,
                    'chembl_target_id':t(target['chembl_id']),
                    **{k:row[k] for k in ('activity_id','assay_chembl_id','molecule_chembl_id',
                       'canonical_smiles','standard_type','standard_relation','standard_value','standard_units',
                       'assay_type','confidence_score','potential_duplicate','data_validity_comment',
                       'document_chembl_id','doi')},'evidence_status':status,'scientific_acceptance':'NOT_APPROVED'})
            coverage.append({'gene_symbol':gene,'status':'LOCAL_ASSAYS_FOUND' if total else 'NO_ROWS_IN_STRICT_LOCAL_QUERY',
                             'target_matches':1,'candidate_rows':total,'rows_retained':len(chosen),'truncated':total>len(chosen)})
        return results,coverage,{'schema_checked':True,'rows_retrieved':len(results),
                                 'confirmed_pyZ_binding':False}
    finally:db.close()


@contextmanager
def open_binding_tsv(path:Path):
    try:csv.field_size_limit(sys.maxsize)
    except OverflowError:csv.field_size_limit(2**31-1)
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as z:
            members=[i for i in z.infolist() if not i.is_dir() and i.filename.lower().endswith('.tsv')]
            if len(members)!=1:raise ValueError('BindingDB ZIP must contain exactly one TSV')
            info=members[0]
            if info.file_size>12_000_000_000 or info.filename.startswith('/') or '..' in Path(info.filename).parts:
                raise ValueError('Unsafe BindingDB TSV in ZIP')
            if info.compress_size>0 and info.file_size/info.compress_size>150:
                raise ValueError('BindingDB ZIP expansion ratio exceeds limit')
            with z.open(info) as raw,io.TextIOWrapper(raw,encoding='utf-8-sig',newline='') as text:
                yield csv.reader(text,delimiter='\t')
    else:
        with path.open(encoding='utf-8-sig',newline='') as text:
            yield csv.reader(text,delimiter='\t')


def single_header(header:list[str],name:str,required:bool=True) -> int|None:
    ids=[i for i,x in enumerate(header) if x.strip()==name]
    if len(ids)!=1 and required:raise ValueError(f'BindingDB header must include exactly one {name!r} column')
    return ids[0] if len(ids)==1 else None


def bindingdb_local(path:Path,panel:list[dict],max_rows:int) -> tuple[list[dict],list[dict],dict]:
    accessions={x['uniprot_accession']:x['gene_symbol'] for x in panel if x['uniprot_accession']}
    matched=Counter();saved=Counter();multichain=Counter();excluded=Counter();overlimit=Counter();observations=[]
    count=0; schema_mismatch=0
    with open_binding_tsv(path) as reader:
        header=next(reader,None)
        if not header:raise ValueError('Empty BindingDB TSV')
        chain_cols=[i for i,c in enumerate(header) if c.strip().startswith('UniProt (SwissProt) Primary ID of Target Chain')]
        if not chain_cols:raise ValueError('BindingDB missing SwissProt primary target chain columns')
        fields=['BindingDB Reactant_set_id','Ligand SMILES','Ligand InChI Key','BindingDB MonomerID',
                'Target Name','Ki (nM)','Kd (nM)','IC50 (nM)']
        columns={x:single_header(header,x) for x in fields}
        optional=['Target Source Organism According to Curator or DataSource','PMID','Article DOI','Curation/DataSource',
                  'Number of Protein Chains in Target (>1 implies a multichain complex)']
        columns.update({x:single_header(header,x,False) for x in optional})
        def at(row,name):
            col=columns[name]
            return t(row[col]) if col is not None else ''
        for n,row in enumerate(reader,start=2):
            count+=1
            if count%500000==0: print(f'[P0.8] BindingDB scanned {count:,} TSV rows',flush=True)
            if len(row)>len(header):schema_mismatch+=1;continue
            if len(row)<len(header):row.extend(['']*(len(header)-len(row)))
            chains=[t(row[i]).upper() for i in chain_cols if t(row[i])]
            related=sorted(set(chains)&set(accessions))
            if not related:continue
            chain_count=at(row,optional[4]);chain_n=int(chain_count) if chain_count.isdigit() else len(chains)
            unique_chains=sorted(set(chains))
            organ=at(row,optional[0]);human=not organ or organ.lower() in ('homo sapiens','human')
            for acc in related:
                gene=accessions[acc];matched[gene]+=1
                id_status=('REVIEW_MULTI_CHAIN' if chain_n>1 or len(chains)>1 else
                           'REVIEW_SPECIES_MISMATCH' if not human else
                           'SWISSPROT_PRIMARY_SINGLE_CHAIN_MATCH_ASSAY_REVIEW')
                if id_status!='SWISSPROT_PRIMARY_SINGLE_CHAIN_MATCH_ASSAY_REVIEW':multichain[gene]+=1
                for endpoint in ENDPOINTS:
                    raw=at(row,endpoint+' (nM)')
                    if not raw:continue
                    relation,numeric,qual=parse_nm(raw)
                    if qual=='INVALID_OR_MISSING':excluded[gene]+=1;continue
                    if saved[gene]>=max_rows:overlimit[gene]+=1;continue
                    status=('EXACT_NUMERIC_STRUCTURE_PRESENT' if qual=='EXACT_NUMERIC' and at(row,'Ligand SMILES') else
                            'CENSORED_NOT_EXACT' if qual=='CENSORED' else 'REVIEW_MISSING_STRUCTURE')
                    observations.append({'gene_symbol':gene,'uniprot_accession':acc,'source_row':n,
                        'bindingdb_reactant_set_id':at(row,'BindingDB Reactant_set_id'),
                        'bindingdb_monomer_id':at(row,'BindingDB MonomerID'),
                        'all_chain_accessions':';'.join(unique_chains),'declared_chain_count':chain_count,
                        'target_name':at(row,'Target Name'),'organism':organ,'endpoint':endpoint,'raw_value_nm':raw,
                        'relation':relation,'numeric_value_nm':numeric,'ligand_smiles':at(row,'Ligand SMILES'),
                        'ligand_inchikey':at(row,'Ligand InChI Key'),'article_doi':at(row,'Article DOI'),
                        'pmid':at(row,'PMID'),'curation_source':at(row,'Curation/DataSource'),
                        'identity_status':id_status,'measurement_status':status,
                        'scientific_acceptance':'NOT_APPROVED'})
                    saved[gene]+=1
    coverage=[{'gene_symbol':x['gene_symbol'],'status':'ROWS_MATCHED' if matched[x['gene_symbol']] else 'NO_UNIPROT_CHAIN_MATCH_IN_LOCAL_FILE',
               'chain_match_rows':matched[x['gene_symbol']],'multi_or_species_review_rows':multichain[x['gene_symbol']],
               'rows_retained':saved[x['gene_symbol']],'rows_over_limit':overlimit[x['gene_symbol']],'rows_invalid':excluded[x['gene_symbol']],
               'truncated':overlimit[x['gene_symbol']]>0} for x in panel]
    return observations,coverage,{'rows_scanned':count,'rows_schema_mismatch':schema_mismatch,
                                'strict_primary_chain_name_matching':True,'source_assay_verified':False}


def run(*,p06_zip:Path,p0_zip:Path|None=None,p07_zip:Path|None=None,
        chembl_sqlite:Path|None=None,bindingdb_tsv:Path|None=None,out:Path,
        panel:str='focus6',max_rows_per_target:int=1000,release_chembl:str='',release_bindingdb:str='') -> dict:
    if panel not in ('focus6','sample40'):raise ValueError('Invalid panel')
    paths=[x.resolve() for x in (p06_zip,p0_zip,p07_zip,chembl_sqlite,bindingdb_tsv) if x is not None]
    if not p06_zip.is_file():raise FileNotFoundError('P0.6 ZIP missing')
    if chembl_sqlite and not chembl_sqlite.is_file():raise FileNotFoundError('ChEMBL SQLite missing')
    if bindingdb_tsv and not bindingdb_tsv.is_file():raise FileNotFoundError('BindingDB TSV/ZIP missing')
    out=out.resolve()
    if out.exists() or any(out==p or out in p.parents for p in paths):
        raise ValueError('Output must be new and outside input folders')
    if chembl_sqlite and not re.fullmatch(r'chembl_\d+',release_chembl):
        raise ValueError('Explicit ChEMBL release is required e.g. chembl_36')
    if bindingdb_tsv and not re.fullmatch(r'\d{6}',release_bindingdb):
        raise ValueError('Explicit BindingDB YYYYMM release required')
    panel_rows,linkage=rows_from_pilot(p06_zip,p0_zip,p07_zip,panel)
    if not 1<=max_rows_per_target<=10000:raise ValueError('Max rows per target must be 1..10000')
    inputs={str(x):{'size_bytes':x.stat().st_size,'sha256_before':shafile(x)} for x in paths}
    qc=[];status={'chembl':'SOURCE_NOT_PROVIDED_NOT_ASSESSED','bindingdb':'SOURCE_NOT_PROVIDED_NOT_ASSESSED'}
    chem_rows=[];chem_cov=[];bdb_rows=[];bdb_cov=[];stats={}
    if chembl_sqlite:
        try:
            chem_rows,chem_cov,stats['chembl']=chembl_local(chembl_sqlite,panel_rows,max_rows_per_target)
            status['chembl']='LOCAL_SOURCE_CHECKED'
        except (ValueError,sqlite3.Error) as e:
            status['chembl']='SOURCE_SCHEMA_OR_READ_ERROR'
            qc.append({'severity':'BLOCKER','source':'ChEMBL','code':type(e).__name__,'detail':t(e)[:500]})
    else:qc.append({'severity':'REVIEW','source':'ChEMBL','code':'MISSING_NOT_ASSESSED','detail':'Provide pinned ChEMBL SQLite explicitly'})
    if bindingdb_tsv:
        try:
            bdb_rows,bdb_cov,stats['bindingdb']=bindingdb_local(bindingdb_tsv,panel_rows,max_rows_per_target)
            status['bindingdb']='LOCAL_SOURCE_CHECKED'
        except (ValueError,zipfile.BadZipFile,csv.Error,UnicodeError) as e:
            status['bindingdb']='SOURCE_SCHEMA_OR_READ_ERROR'
            qc.append({'severity':'BLOCKER','source':'BindingDB','code':type(e).__name__,'detail':t(e)[:500]})
    else:qc.append({'severity':'REVIEW','source':'BindingDB','code':'MISSING_NOT_ASSESSED','detail':'Provide pinned full BindingDB TSV/ZIP explicitly'})
    for x in panel_rows:
        if not x['uniprot_accession']:
            qc.append({'severity':'REVIEW','source':'UniProt','code':'UNRESOLVED_MAPPING',
                       'detail':x['gene_symbol']+' : '+x['mapping_status']})
    for path in paths:
        digest=shafile(path)
        inputs[str(path)]['sha256_after']=digest
        if digest!=inputs[str(path)]['sha256_before']:
            raise RuntimeError('Source changed during run: '+str(path))
    # Both source types report experimental ligand activity only; PYZ not assessed.
    technical=('BLOCKED' if any(x['severity']=='BLOCKER' for x in qc) else
               'INPUT_REQUIRED' if not chembl_sqlite or not bindingdb_tsv else 'REVIEW_REQUIRED')
    summary={'contract':CONTRACT,'technical_status':technical,'scientific_acceptance':'NOT_APPROVED',
             'merge_permitted':False,'external_requests_n':0,'panel_mode':panel,'panel_genes_n':len(panel_rows),
             'source_status':status,'releases':{'chembl':release_chembl or None,'bindingdb':release_bindingdb or None},
             'chembl_rows_reviewed_n':len(chem_rows),
             'chembl_assay_metadata_eligible_n':sum(x['evidence_status']=='ASSAY_METADATA_MATCHED_PRIMARY_REVIEW_REQUIRED' for x in chem_rows),
             'bindingdb_rows_reviewed_n':len(bdb_rows),
             'bindingdb_exact_single_chain_rows_n':sum(x['identity_status']=='SWISSPROT_PRIMARY_SINGLE_CHAIN_MATCH_ASSAY_REVIEW' for x in bdb_rows),
             'scientific_verified_binding_to_pyz_n':0,'source_input_hashes_unchanged':True,
             'baseline_65x13_untouched':True,'linkage':linkage,'source_stats':stats,
             'limitations':['No extrapolation beyond chosen panel or scanned release',
                            'ChEMBL assay confidence is curation metadata, not independent primary literature review',
                            'BindingDB SwissProt target-chain identity is not a verified biochemical assay mechanism',
                            'Censored values and duplicate/reference records cannot be aggregated into a single potency',
                            'No PYZ activity, binding, cell selectivity, or experimental validation inferred']}
    out.mkdir(parents=True,exist_ok=False)
    output_tsv(out/'p08_target_panel.tsv',panel_rows,list(panel_rows[0]))
    output_tsv(out/'p08_chembl_observations.tsv',chem_rows,OBS_CHEMBL)
    output_tsv(out/'p08_bindingdb_observations.tsv',bdb_rows,OBS_BINDING)
    output_tsv(out/'p08_chembl_coverage.tsv',chem_cov,
               ['gene_symbol','status','target_matches','candidate_rows','rows_retained','truncated'])
    output_tsv(out/'p08_bindingdb_coverage.tsv',bdb_cov,
               ['gene_symbol','status','chain_match_rows','multi_or_species_review_rows','rows_retained','rows_over_limit','rows_invalid','truncated'])
    output_tsv(out/'p08_QC.tsv',qc,['severity','source','code','detail'])
    with (out/'p08_summary.json').open('w',encoding='utf-8') as f:json.dump(summary,f,indent=2,ensure_ascii=False);f.write('\n')
    manifest={'contract':CONTRACT,'input_paths_sha256':inputs,'no_network':True,
              'read_only_scientific_data':True,'exact_target_uniprot_required':True,
              'no_automatic_scientific_acceptance':True}
    with (out/'p08_source_manifest.json').open('w',encoding='utf-8') as f:json.dump(manifest,f,indent=2,ensure_ascii=False);f.write('\n')
    report=f'''# TLS V3 / P0.8 — локальная проверка экспериментальных источников

**Статус:** {technical}. **Научная приёмка:** NOT_APPROVED. **Слияние:** запрещено.

- Панель: {panel}; проверено генов {len(panel_rows)}. Сеть: 0 запросов.
- ChEMBL: {status['chembl']}, release {release_chembl or 'НЕ ПРЕДОСТАВЛЕН'}; строк {len(chem_rows)}; предварительных точных измерений с проверенной метаинформацией assay: {summary['chembl_assay_metadata_eligible_n']}.
- BindingDB: {status['bindingdb']}, release {release_bindingdb or 'НЕ ПРЕДОСТАВЛЕН'}; строк {len(bdb_rows)}; с точным single-chain SwissProt: {summary['bindingdb_exact_single_chain_rows_n']}.
- Связь с 65 PYZ, докинг, причинность и селективность **не устанавливались**.

Статус `SOURCE_NOT_PROVIDED_NOT_ASSESSED` не означает отрицательный результат. Подтверждение по UniProt, ChEMBL или BindingDB не заменяет первичные протоколы и оценку качества измерений. Нельзя переносить эти строки в действующую матрицу PYZ до научной приёмки.
'''
    (out/'P08_REPORT_RU.md').write_text(report,encoding='utf-8')
    print(json.dumps({k:summary[k] for k in ('technical_status','panel_genes_n','source_status',
          'chembl_rows_reviewed_n','bindingdb_rows_reviewed_n','scientific_acceptance')},ensure_ascii=False),flush=True)
    return summary


def main()->None:
    p=argparse.ArgumentParser(description='MCL TLS V3 P0.8 read-only local ChEMBL+BindingDB validation')
    p.add_argument('--p06-zip',type=Path,required=True)
    p.add_argument('--p0-zip',type=Path)
    p.add_argument('--p07-zip',type=Path)
    p.add_argument('--chembl-sqlite',type=Path)
    p.add_argument('--bindingdb-tsv',type=Path)
    p.add_argument('--chembl-release',default='')
    p.add_argument('--bindingdb-release',default='')
    p.add_argument('--panel',choices=('focus6','sample40'),default='focus6')
    p.add_argument('--max-rows-per-target',type=int,default=1000)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    run(p06_zip=a.p06_zip,p0_zip=a.p0_zip,p07_zip=a.p07_zip,
        chembl_sqlite=a.chembl_sqlite,bindingdb_tsv=a.bindingdb_tsv,
        release_chembl=a.chembl_release,release_bindingdb=a.bindingdb_release,
        panel=a.panel,max_rows_per_target=a.max_rows_per_target,out=a.out)

if __name__=='__main__':main()