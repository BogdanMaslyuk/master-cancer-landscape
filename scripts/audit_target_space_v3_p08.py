#!/usr/bin/env python3
"""MCL TLS V3/P0.8: read-only local ChEMBL SQLite + BindingDB TSV audit.

Inputs are explicitly supplied version-pinned files. No HTTP, no modification of
scientific datasets, no PYZ activity or causal target assertions.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import sqlite3
import sys
import zipfile
from collections import Counter, defaultdict
from contextlib import contextmanager
from pathlib import Path
from typing import Any

CONTRACT = 'mcl-tls-v3-p08-local-source-validation-v1'
P0_CONTRACT = 'mcl-tls-v3-p0-inventory-v1'
P06_CONTRACT = 'mcl-tls-v3-p06-target-coverage-pilot-v1'
P07_CONTRACT = 'mcl-tls-v3-p07-source-reliability-v1'
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
              'document_chembl_id','doi','pmid','measurement_value_status','evidence_status',
              'primary_experiment_review_status','cross_source_potential_duplicate','scientific_acceptance']
OBS_BINDING = ['gene_symbol','uniprot_accession','source_row','bindingdb_reactant_set_id','bindingdb_monomer_id',
               'raw_chain_accessions','all_chain_accessions','populated_chain_count','declared_chain_count',
               'target_name','organism',
               'raw_ki_nm','raw_kd_nm','raw_ic50_nm','endpoint','raw_value_nm','relation','numeric_value_nm',
               'ligand_smiles','ligand_inchikey','article_doi','pmid','curation_source','chain_status',
               'species_status','measurement_value_status','measurement_status',
               'primary_experiment_review_status','cross_source_potential_duplicate','scientific_acceptance']

CHEMBL_COVERAGE = ['gene_symbol','status','target_matches','candidate_rows','rows_retained','exact_numeric_n',
                   'censored_n','single_component_protein_objects_n','multi_component_protein_objects_n',
                   'species_missing_objects_n','species_mismatch_objects_n','doi_records_n','pmid_records_n',
                   'primary_experiment_review_records_n','curator_potential_duplicate_records_n',
                   'rows_truncated_n','truncated']
BINDING_COVERAGE = ['gene_symbol','status','chain_match_rows','single_chain_objects_n','multi_chain_objects_n',
                    'unknown_chain_count_objects_n','species_missing_objects_n','species_mismatch_objects_n',
                    'exact_numeric_n','censored_n','doi_records_n','pmid_records_n',
                    'primary_experiment_review_records_n','rows_retained','rows_over_limit','rows_invalid',
                    'rows_truncated_n','truncated']
TARGET_EVIDENCE = ['gene_symbol','uniprot_accession','chembl_status','chembl_exact_numeric_n',
                   'chembl_censored_n','chembl_single_component_protein_objects_n',
                   'chembl_multi_component_protein_objects_n','chembl_species_missing_objects_n',
                   'chembl_species_mismatch_objects_n','chembl_doi_records_n','chembl_pmid_records_n',
                   'chembl_primary_experiment_review_records_n','chembl_curator_potential_duplicate_records_n',
                   'chembl_rows_truncated_n','bindingdb_status','bindingdb_exact_numeric_n',
                   'bindingdb_censored_n','bindingdb_single_chain_objects_n','bindingdb_multi_chain_objects_n',
                   'bindingdb_unknown_chain_count_objects_n','bindingdb_species_missing_objects_n',
                   'bindingdb_species_mismatch_objects_n','bindingdb_doi_records_n','bindingdb_pmid_records_n',
                   'bindingdb_primary_experiment_review_records_n','bindingdb_rows_truncated_n',
                   'cross_source_duplicate_status','cross_source_potential_duplicate_pairs_n']


def shafile(path: Path) -> str:
    sha=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b''):sha.update(chunk)
    return sha.hexdigest()


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def require_sha256(value: Any, label: str) -> str:
    digest=t(value)
    if not re.fullmatch(r'[0-9a-f]{64}',digest):
        raise ValueError(label+' must be a lowercase SHA-256 hexadecimal digest')
    return digest


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
                if p07meta.get('contract')!=P07_CONTRACT:
                    raise ValueError('P0.7 contract mismatch')
                if p07meta.get('scientific_acceptance')!='NOT_APPROVED':
                    raise ValueError('Unexpected P0.7 scientific acceptance state')
                recorded=p07meta.get('source_archives_sha256',{})
                if not isinstance(recorded,dict):
                    raise TypeError('P0.7 source archive SHA-256 map missing')
                if require_sha256(recorded.get(p06zip.name),'P0.7 recorded P0.6 SHA-256')!=linkage['p06_sha256']:
                    raise ValueError('P0.7 and P0.6 provenance mismatch')
                if p0zip and require_sha256(recorded.get(p0zip.name),'P0.7 recorded P0 SHA-256')!=linkage['p0_sha256']:
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
        pmid_col='d.pubmed_id AS pmid' if 'pubmed_id' in schema['docs'] else 'NULL AS pmid'
        doc_chembl='d.chembl_id AS document_chembl_id' if {'chembl_id','doc_id'}.issubset(schema['docs']) else 'NULL AS document_chembl_id'
        use_docs={'chembl_id','doc_id'}.issubset(schema['docs']) and 'doc_id' in extra
        doc_join='LEFT JOIN docs d ON a.doc_id=d.doc_id' if use_docs else ''
        if not use_docs:
            doi_col='NULL AS doi';pmid_col='NULL AS pmid';doc_chembl='NULL AS document_chembl_id'
        results=[];coverage=[]
        for entry in panel:
            gene,accession=entry['gene_symbol'],entry['uniprot_accession']
            if not accession:
                coverage.append({'gene_symbol':gene,'status':'NOT_ASSESSED_UNRESOLVED_MAPPING',
                                 'target_matches':0,'candidate_rows':'','rows_retained':'','truncated':''})
                continue
            matches=db.execute('''SELECT DISTINCT t.tid,t.chembl_id,t.target_type,t.organism
                FROM target_dictionary t JOIN target_components tc ON tc.tid=t.tid
                JOIN component_sequences cs ON cs.component_id=tc.component_id
                WHERE cs.accession=?''',(accession,)).fetchall()
            exact=[];single_objects=0;multi_objects=0;species_missing=0;species_mismatch=0
            for row in matches:
                components=db.execute('''SELECT c.accession FROM target_components x
                    JOIN component_sequences c ON x.component_id=c.component_id WHERE x.tid=?''',(row['tid'],)).fetchall()
                component_accessions=[t(x[0]) for x in components]
                is_single=(t(row['target_type'])=='SINGLE PROTEIN' and len(components)==1 and
                           component_accessions[0]==accession)
                if is_single:single_objects+=1
                elif len(components)>1:multi_objects+=1
                organism=t(row['organism'])
                if not organism:species_missing+=1
                elif organism.casefold() not in ('homo sapiens','human'):species_mismatch+=1
                if is_single and organism.casefold() in ('homo sapiens','human'):exact.append(row)
            if len(exact)!=1:
                coverage.append({'gene_symbol':gene,'status':'NO_EXACT_SINGLE_TARGET' if not exact else 'AMBIGUOUS_EXACT_TARGETS',
                                 'target_matches':len(exact),'candidate_rows':'','rows_retained':'',
                                 'single_component_protein_objects_n':single_objects,
                                 'multi_component_protein_objects_n':multi_objects,
                                 'species_missing_objects_n':species_missing,
                                 'species_mismatch_objects_n':species_mismatch,'truncated':''})
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
                    {doc_chembl}, {doi_col}, {pmid_col}
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
                value_status=('EXACT_NUMERIC' if isexact else 'CENSORED' if rel and rel!='=' and
                              math.isfinite(val) and val>0 else 'INVALID_OR_MISSING')
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
                       'document_chembl_id','doi','pmid')},'measurement_value_status':value_status,
                    'evidence_status':status,'primary_experiment_review_status':'REVIEW_REQUIRED',
                    'cross_source_potential_duplicate':'NOT_ASSESSED','scientific_acceptance':'NOT_APPROVED'})
            gene_rows=results[-len(chosen):] if chosen else []
            coverage.append({'gene_symbol':gene,'status':'LOCAL_ASSAYS_FOUND' if total else 'NO_ROWS_IN_STRICT_LOCAL_QUERY',
                             'target_matches':1,'candidate_rows':total,'rows_retained':len(chosen),
                             'exact_numeric_n':sum(x['measurement_value_status']=='EXACT_NUMERIC' for x in gene_rows),
                             'censored_n':sum(x['measurement_value_status']=='CENSORED' for x in gene_rows),
                             'single_component_protein_objects_n':single_objects,
                             'multi_component_protein_objects_n':multi_objects,
                             'species_missing_objects_n':species_missing,
                             'species_mismatch_objects_n':species_mismatch,
                             'doi_records_n':sum(bool(t(x['doi'])) for x in gene_rows),
                             'pmid_records_n':sum(bool(t(x['pmid'])) for x in gene_rows),
                             'primary_experiment_review_records_n':len(gene_rows),
                             'curator_potential_duplicate_records_n':sum(as_true(x['potential_duplicate']) for x in gene_rows),
                             'rows_truncated_n':total-len(chosen),'truncated':total>len(chosen)})
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
    matched=Counter();saved=Counter();singlechain=Counter();multichain=Counter();unknownchain=Counter()
    species_missing=Counter();species_mismatch=Counter();doi_rows=Counter();pmid_rows=Counter()
    exact=Counter();censored=Counter();invalid=Counter();overlimit=Counter();observations=[]
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
        def raw_at(row,name):
            col=columns[name]
            return '' if col is None or row[col] is None else str(row[col])
        for n,row in enumerate(reader,start=2):
            count+=1
            if count%500000==0: print(f'[P0.8] BindingDB scanned {count:,} TSV rows',flush=True)
            if len(row)>len(header):schema_mismatch+=1;continue
            if len(row)<len(header):row.extend(['']*(len(header)-len(row)))
            raw_chains=[str(row[i]) for i in chain_cols if row[i] is not None and str(row[i])!='']
            chains=[t(value).upper() for value in raw_chains if t(value)]
            related=sorted(set(chains)&set(accessions))
            if not related:continue
            raw_chain_count=raw_at(row,optional[4]);chain_count=t(raw_chain_count)
            unique_chains=sorted(set(chains))
            if (chain_count.isdigit() and int(chain_count)>1) or len(chains)>1:
                chain_status='REVIEW_MULTI_CHAIN'
            elif chain_count=='1' and len(chains)==1:
                chain_status='SWISSPROT_PRIMARY_SINGLE_CHAIN_MATCH'
            else:
                chain_status='REVIEW_CHAIN_COUNT_UNKNOWN'
            raw_organ=raw_at(row,optional[0]);organ=t(raw_organ)
            if not organ:species_status='REVIEW_SPECIES_MISSING'
            elif organ.casefold() in ('homo sapiens','human'):species_status='HUMAN_SPECIES_MATCH'
            else:species_status='REVIEW_SPECIES_MISMATCH'
            for acc in related:
                gene=accessions[acc];matched[gene]+=1
                if chain_status=='SWISSPROT_PRIMARY_SINGLE_CHAIN_MATCH':singlechain[gene]+=1
                elif chain_status=='REVIEW_MULTI_CHAIN':multichain[gene]+=1
                else:unknownchain[gene]+=1
                if species_status=='REVIEW_SPECIES_MISSING':species_missing[gene]+=1
                elif species_status=='REVIEW_SPECIES_MISMATCH':species_mismatch[gene]+=1
                if at(row,'Article DOI'):doi_rows[gene]+=1
                if at(row,'PMID'):pmid_rows[gene]+=1
                raw_values={endpoint:raw_at(row,endpoint+' (nM)') for endpoint in ENDPOINTS}
                for endpoint in ENDPOINTS:
                    raw=raw_values[endpoint]
                    if not raw:continue
                    relation,numeric,qual=parse_nm(raw)
                    if saved[gene]>=max_rows:overlimit[gene]+=1;continue
                    if qual=='EXACT_NUMERIC':exact[gene]+=1
                    elif qual=='CENSORED':censored[gene]+=1
                    else:invalid[gene]+=1
                    status=('EXACT_NUMERIC_STRUCTURE_PRESENT' if qual=='EXACT_NUMERIC' and at(row,'Ligand SMILES') else
                            'CENSORED_NOT_EXACT' if qual=='CENSORED' else
                            'REVIEW_INVALID_OR_MISSING_VALUE' if qual=='INVALID_OR_MISSING' else
                            'REVIEW_MISSING_STRUCTURE')
                    observations.append({'gene_symbol':gene,'uniprot_accession':acc,'source_row':n,
                        'bindingdb_reactant_set_id':at(row,'BindingDB Reactant_set_id'),
                        'bindingdb_monomer_id':at(row,'BindingDB MonomerID'),
                        'raw_chain_accessions':';'.join(raw_chains),
                        'all_chain_accessions':';'.join(unique_chains),'populated_chain_count':len(chains),
                        'declared_chain_count':raw_chain_count,'target_name':at(row,'Target Name'),'organism':raw_organ,
                        'raw_ki_nm':raw_values['Ki'],'raw_kd_nm':raw_values['Kd'],
                        'raw_ic50_nm':raw_values['IC50'],'endpoint':endpoint,'raw_value_nm':raw,
                        'relation':relation,'numeric_value_nm':numeric,'ligand_smiles':at(row,'Ligand SMILES'),
                        'ligand_inchikey':at(row,'Ligand InChI Key'),'article_doi':at(row,'Article DOI'),
                        'pmid':at(row,'PMID'),'curation_source':at(row,'Curation/DataSource'),
                        'chain_status':chain_status,'species_status':species_status,
                        'measurement_value_status':qual,'measurement_status':status,
                        'primary_experiment_review_status':'REVIEW_REQUIRED',
                        'cross_source_potential_duplicate':'NOT_ASSESSED',
                        'scientific_acceptance':'NOT_APPROVED'})
                    saved[gene]+=1
    coverage=[{'gene_symbol':x['gene_symbol'],'status':'ROWS_MATCHED' if matched[x['gene_symbol']] else 'NO_UNIPROT_CHAIN_MATCH_IN_LOCAL_FILE',
               'chain_match_rows':matched[x['gene_symbol']],'single_chain_objects_n':singlechain[x['gene_symbol']],
               'multi_chain_objects_n':multichain[x['gene_symbol']],
               'unknown_chain_count_objects_n':unknownchain[x['gene_symbol']],
               'species_missing_objects_n':species_missing[x['gene_symbol']],
               'species_mismatch_objects_n':species_mismatch[x['gene_symbol']],
               'exact_numeric_n':exact[x['gene_symbol']],'censored_n':censored[x['gene_symbol']],
               'doi_records_n':doi_rows[x['gene_symbol']],'pmid_records_n':pmid_rows[x['gene_symbol']],
               'primary_experiment_review_records_n':matched[x['gene_symbol']],
               'rows_retained':saved[x['gene_symbol']],'rows_over_limit':overlimit[x['gene_symbol']],
               'rows_invalid':invalid[x['gene_symbol']],'rows_truncated_n':overlimit[x['gene_symbol']],
               'truncated':overlimit[x['gene_symbol']]>0} for x in panel]
    return observations,coverage,{'rows_scanned':count,'rows_schema_mismatch':schema_mismatch,
                                'strict_primary_chain_name_matching':True,'source_assay_verified':False}


def source_metadata(path:Path,source:str,release:str,actual_sha256:str,*,publisher_sha256:str='',
                    license_id:str='',license_url:str='') -> dict:
    expected=t(publisher_sha256)
    if expected:
        expected=require_sha256(expected,f'{source} publisher SHA-256')
        if expected!=actual_sha256:
            raise ValueError(f'{source} publisher SHA-256 mismatch')
        checksum_status='MATCH'
    else:checksum_status='NOT_PROVIDED_NOT_ASSESSED'
    if source=='ChEMBL':
        found=re.search(r'chembl_(\d+)',path.name,re.IGNORECASE)
        filename_release='chembl_'+found.group(1) if found else ''
    else:
        found=re.search(r'BindingDB_All_(\d{6})',path.name,re.IGNORECASE)
        filename_release=found.group(1) if found else ''
    if filename_release and filename_release.casefold()!=release.casefold():
        raise ValueError(f'{source} release label disagrees with source filename')
    return {'release':release,'filename_release':filename_release or None,
            'filename_release_status':'MATCH' if filename_release else 'FILENAME_NOT_INFORMATIVE_NOT_ASSESSED',
            'actual_sha256':actual_sha256,'publisher_sha256':expected or None,
            'publisher_sha256_status':checksum_status,'license_id':t(license_id) or None,
            'license_url':t(license_url) or None,
            'license_status':'PROVIDED_UNVERIFIED' if t(license_id) or t(license_url) else 'NOT_PROVIDED_NOT_ASSESSED'}


def _normal_doi(value:Any) -> str:
    doi=t(value).casefold()
    for prefix in ('https://doi.org/','http://doi.org/','doi:'):
        if doi.startswith(prefix):doi=doi[len(prefix):].strip()
    return doi


def _numeric_text(value:Any) -> str:
    try:number=float(t(value))
    except (TypeError,ValueError):return ''
    return format(number,'.12g') if math.isfinite(number) and number>0 else ''


def cross_source_duplicates(chem_rows:list[dict],bdb_rows:list[dict],assessed:bool) -> list[dict]:
    if not assessed:return []
    for row in chem_rows+bdb_rows:row['cross_source_potential_duplicate']='FALSE'
    chem_by_key=defaultdict(list)
    for row in chem_rows:
        key=(row['gene_symbol'],_normal_doi(row.get('doi')),t(row.get('standard_type')),
             t(row.get('standard_relation')) or '=',_numeric_text(row.get('standard_value')),
             t(row.get('canonical_smiles')))
        if all(key[1:]):chem_by_key[key].append(row)
    pairs=[]
    for bdb in bdb_rows:
        key=(bdb['gene_symbol'],_normal_doi(bdb.get('article_doi')),t(bdb.get('endpoint')),
             t(bdb.get('relation')),_numeric_text(bdb.get('numeric_value_nm')),t(bdb.get('ligand_smiles')))
        if not all(key[1:]):continue
        for chem in chem_by_key.get(key,[]):
            chem['cross_source_potential_duplicate']='TRUE'
            bdb['cross_source_potential_duplicate']='TRUE'
            pairs.append({'gene_symbol':key[0],'doi':key[1],'endpoint':key[2],'relation':key[3],
                          'numeric_value_nm':key[4],'ligand_smiles':key[5],
                          'chembl_activity_id':chem['activity_id'],
                          'bindingdb_source_row':bdb['source_row'],
                          'assessment':'POTENTIAL_DUPLICATE_EXACT_FIELDS_PRIMARY_REVIEW_REQUIRED'})
    return pairs


def unassessed_coverage(panel:list[dict],fields:list[str],status:str) -> list[dict]:
    return [{**{field:'' for field in fields},'gene_symbol':x['gene_symbol'],'status':status} for x in panel]


def target_evidence_rows(panel:list[dict],chem_cov:list[dict],bdb_cov:list[dict],
                         duplicate_pairs:list[dict],duplicates_assessed:bool) -> list[dict]:
    chem={x['gene_symbol']:x for x in chem_cov};bdb={x['gene_symbol']:x for x in bdb_cov}
    duplicate_counts=Counter(x['gene_symbol'] for x in duplicate_pairs)
    rows=[]
    for target in panel:
        gene=target['gene_symbol'];c=chem[gene];b=bdb[gene]
        rows.append({'gene_symbol':gene,'uniprot_accession':target['uniprot_accession'],
          'chembl_status':c['status'],'chembl_exact_numeric_n':c.get('exact_numeric_n',''),
          'chembl_censored_n':c.get('censored_n',''),
          'chembl_single_component_protein_objects_n':c.get('single_component_protein_objects_n',''),
          'chembl_multi_component_protein_objects_n':c.get('multi_component_protein_objects_n',''),
          'chembl_species_missing_objects_n':c.get('species_missing_objects_n',''),
          'chembl_species_mismatch_objects_n':c.get('species_mismatch_objects_n',''),
          'chembl_doi_records_n':c.get('doi_records_n',''),'chembl_pmid_records_n':c.get('pmid_records_n',''),
          'chembl_primary_experiment_review_records_n':c.get('primary_experiment_review_records_n',''),
          'chembl_curator_potential_duplicate_records_n':c.get('curator_potential_duplicate_records_n',''),
          'chembl_rows_truncated_n':c.get('rows_truncated_n',''),'bindingdb_status':b['status'],
          'bindingdb_exact_numeric_n':b.get('exact_numeric_n',''),'bindingdb_censored_n':b.get('censored_n',''),
          'bindingdb_single_chain_objects_n':b.get('single_chain_objects_n',''),
          'bindingdb_multi_chain_objects_n':b.get('multi_chain_objects_n',''),
          'bindingdb_unknown_chain_count_objects_n':b.get('unknown_chain_count_objects_n',''),
          'bindingdb_species_missing_objects_n':b.get('species_missing_objects_n',''),
          'bindingdb_species_mismatch_objects_n':b.get('species_mismatch_objects_n',''),
          'bindingdb_doi_records_n':b.get('doi_records_n',''),'bindingdb_pmid_records_n':b.get('pmid_records_n',''),
          'bindingdb_primary_experiment_review_records_n':b.get('primary_experiment_review_records_n',''),
          'bindingdb_rows_truncated_n':b.get('rows_truncated_n',''),
          'cross_source_duplicate_status':'ASSESSED_EXACT_FIELDS_ONLY' if duplicates_assessed else 'NOT_ASSESSED_SOURCE_MISSING',
          'cross_source_potential_duplicate_pairs_n':duplicate_counts[gene] if duplicates_assessed else ''})
    return rows


def run(*,p06_zip:Path,p0_zip:Path|None=None,p07_zip:Path|None=None,
        chembl_sqlite:Path|None=None,bindingdb_tsv:Path|None=None,out:Path,
        panel:str='focus6',max_rows_per_target:int=1000,release_chembl:str='',release_bindingdb:str='',
        chembl_publisher_sha256:str='',bindingdb_publisher_sha256:str='',
        chembl_license_id:str='',bindingdb_license_id:str='',
        chembl_license_url:str='',bindingdb_license_url:str='') -> dict:
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
    source_meta={}
    if chembl_sqlite:
        source_meta['chembl']=source_metadata(chembl_sqlite,'ChEMBL',release_chembl,
            inputs[str(chembl_sqlite.resolve())]['sha256_before'],publisher_sha256=chembl_publisher_sha256,
            license_id=chembl_license_id,license_url=chembl_license_url)
    if bindingdb_tsv:
        source_meta['bindingdb']=source_metadata(bindingdb_tsv,'BindingDB',release_bindingdb,
            inputs[str(bindingdb_tsv.resolve())]['sha256_before'],publisher_sha256=bindingdb_publisher_sha256,
            license_id=bindingdb_license_id,license_url=bindingdb_license_url)
    qc=[];status={'chembl':'SOURCE_NOT_PROVIDED_NOT_ASSESSED','bindingdb':'SOURCE_NOT_PROVIDED_NOT_ASSESSED'}
    chem_rows=[];chem_cov=[];bdb_rows=[];bdb_cov=[];stats={}
    if chembl_sqlite:
        try:
            chem_rows,chem_cov,stats['chembl']=chembl_local(chembl_sqlite,panel_rows,max_rows_per_target)
            status['chembl']='LOCAL_SOURCE_CHECKED'
        except (ValueError,sqlite3.Error) as e:
            status['chembl']='SOURCE_SCHEMA_OR_READ_ERROR'
            chem_cov=unassessed_coverage(panel_rows,CHEMBL_COVERAGE,'SOURCE_SCHEMA_OR_READ_ERROR_NOT_ASSESSED')
            qc.append({'severity':'BLOCKER','source':'ChEMBL','code':type(e).__name__,'detail':t(e)[:500]})
    else:
        chem_cov=unassessed_coverage(panel_rows,CHEMBL_COVERAGE,'SOURCE_NOT_PROVIDED_NOT_ASSESSED')
        qc.append({'severity':'REVIEW','source':'ChEMBL','code':'MISSING_NOT_ASSESSED','detail':'Provide pinned ChEMBL SQLite explicitly'})
    if bindingdb_tsv:
        try:
            bdb_rows,bdb_cov,stats['bindingdb']=bindingdb_local(bindingdb_tsv,panel_rows,max_rows_per_target)
            status['bindingdb']='LOCAL_SOURCE_CHECKED'
        except (ValueError,zipfile.BadZipFile,csv.Error,UnicodeError) as e:
            status['bindingdb']='SOURCE_SCHEMA_OR_READ_ERROR'
            bdb_cov=unassessed_coverage(panel_rows,BINDING_COVERAGE,'SOURCE_SCHEMA_OR_READ_ERROR_NOT_ASSESSED')
            qc.append({'severity':'BLOCKER','source':'BindingDB','code':type(e).__name__,'detail':t(e)[:500]})
    else:
        bdb_cov=unassessed_coverage(panel_rows,BINDING_COVERAGE,'SOURCE_NOT_PROVIDED_NOT_ASSESSED')
        qc.append({'severity':'REVIEW','source':'BindingDB','code':'MISSING_NOT_ASSESSED','detail':'Provide pinned full BindingDB TSV/ZIP explicitly'})
    for x in panel_rows:
        if not x['uniprot_accession']:
            qc.append({'severity':'REVIEW','source':'UniProt','code':'UNRESOLVED_MAPPING',
                       'detail':x['gene_symbol']+' : '+x['mapping_status']})
    for path in paths:
        digest=shafile(path)
        inputs[str(path)]['sha256_after']=digest
        if digest!=inputs[str(path)]['sha256_before']:
            raise RuntimeError('Source changed during run: '+str(path))
    duplicates_assessed=(status['chembl']=='LOCAL_SOURCE_CHECKED' and status['bindingdb']=='LOCAL_SOURCE_CHECKED')
    duplicate_pairs=cross_source_duplicates(chem_rows,bdb_rows,duplicates_assessed)
    evidence_rows=target_evidence_rows(panel_rows,chem_cov,bdb_cov,duplicate_pairs,duplicates_assessed)
    # Both source types report experimental ligand activity only; PYZ not assessed.
    technical=('BLOCKED' if any(x['severity']=='BLOCKER' for x in qc) else
               'INPUT_REQUIRED' if not chembl_sqlite or not bindingdb_tsv else 'REVIEW_REQUIRED')
    summary={'contract':CONTRACT,'technical_status':technical,'scientific_acceptance':'NOT_APPROVED',
             'merge_permitted':False,'external_requests_n':0,'panel_mode':panel,'panel_genes_n':len(panel_rows),
             'source_status':status,'releases':{'chembl':release_chembl or None,'bindingdb':release_bindingdb or None},
             'chembl_rows_reviewed_n':len(chem_rows),
             'chembl_assay_metadata_eligible_n':sum(x['evidence_status']=='ASSAY_METADATA_MATCHED_PRIMARY_REVIEW_REQUIRED' for x in chem_rows),
             'bindingdb_rows_reviewed_n':len(bdb_rows),
             'bindingdb_single_chain_human_observations_n':sum(
                 x['chain_status']=='SWISSPROT_PRIMARY_SINGLE_CHAIN_MATCH' and
                 x['species_status']=='HUMAN_SPECIES_MATCH' for x in bdb_rows),
             'cross_source_potential_duplicate_pairs_n':len(duplicate_pairs),
             'scientific_verified_binding_to_pyz_n':0,'source_input_hashes_unchanged':True,
             'baseline_65x13_untouched':True,'linkage':linkage,'source_stats':stats,
             'source_metadata':source_meta,
             'limitations':['No extrapolation beyond chosen panel or scanned release',
                            'ChEMBL assay confidence is curation metadata, not independent primary literature review',
                            'BindingDB SwissProt target-chain identity is not a verified biochemical assay mechanism',
                            'Censored values and duplicate/reference records cannot be aggregated into a single potency',
                            'Cross-source duplicate screening uses exact DOI, endpoint, relation, numeric value, and SMILES only',
                            'Provided license metadata is recorded but not independently validated offline',
                            'No PYZ activity, binding, cell selectivity, or experimental validation inferred']}
    out.mkdir(parents=True,exist_ok=False)
    output_tsv(out/'p08_target_panel.tsv',panel_rows,list(panel_rows[0]))
    output_tsv(out/'p08_chembl_observations.tsv',chem_rows,OBS_CHEMBL)
    output_tsv(out/'p08_bindingdb_observations.tsv',bdb_rows,OBS_BINDING)
    output_tsv(out/'p08_chembl_coverage.tsv',chem_cov,CHEMBL_COVERAGE)
    output_tsv(out/'p08_bindingdb_coverage.tsv',bdb_cov,BINDING_COVERAGE)
    output_tsv(out/'p08_target_evidence_counts.tsv',evidence_rows,TARGET_EVIDENCE)
    output_tsv(out/'p08_cross_source_potential_duplicates.tsv',duplicate_pairs,
               ['gene_symbol','doi','endpoint','relation','numeric_value_nm','ligand_smiles',
                'chembl_activity_id','bindingdb_source_row','assessment'])
    output_tsv(out/'p08_QC.tsv',qc,['severity','source','code','detail'])
    with (out/'p08_summary.json').open('w',encoding='utf-8') as f:json.dump(summary,f,indent=2,ensure_ascii=False);f.write('\n')
    manifest={'contract':CONTRACT,'input_paths_sha256':inputs,'source_metadata':source_meta,'no_network':True,
              'read_only_scientific_data':True,'exact_target_uniprot_required':True,
              'no_automatic_scientific_acceptance':True}
    with (out/'p08_source_manifest.json').open('w',encoding='utf-8') as f:json.dump(manifest,f,indent=2,ensure_ascii=False);f.write('\n')
    def shown(value:Any)->str:return 'NOT_ASSESSED' if value=='' or value is None else str(value)
    target_lines=[]
    for row in evidence_rows:
        target_lines.append(
          f"- **{row['gene_symbol']} ({row['uniprot_accession'] or 'mapping unresolved'})** — "
          f"ChEMBL `{row['chembl_status']}`: exact/censored "
          f"{shown(row['chembl_exact_numeric_n'])}/{shown(row['chembl_censored_n'])}, "
          f"single/multi protein objects {shown(row['chembl_single_component_protein_objects_n'])}/"
          f"{shown(row['chembl_multi_component_protein_objects_n'])}, species missing/mismatch "
          f"{shown(row['chembl_species_missing_objects_n'])}/{shown(row['chembl_species_mismatch_objects_n'])}, "
          f"DOI/PMID {shown(row['chembl_doi_records_n'])}/{shown(row['chembl_pmid_records_n'])}, "
          f"primary review {shown(row['chembl_primary_experiment_review_records_n'])}, truncated "
          f"{shown(row['chembl_rows_truncated_n'])}; BindingDB `{row['bindingdb_status']}`: exact/censored "
          f"{shown(row['bindingdb_exact_numeric_n'])}/{shown(row['bindingdb_censored_n'])}, "
          f"single/multi/unknown-chain {shown(row['bindingdb_single_chain_objects_n'])}/"
          f"{shown(row['bindingdb_multi_chain_objects_n'])}/{shown(row['bindingdb_unknown_chain_count_objects_n'])}, "
          f"species missing/mismatch {shown(row['bindingdb_species_missing_objects_n'])}/"
          f"{shown(row['bindingdb_species_mismatch_objects_n'])}, DOI/PMID "
          f"{shown(row['bindingdb_doi_records_n'])}/{shown(row['bindingdb_pmid_records_n'])}, primary review "
          f"{shown(row['bindingdb_primary_experiment_review_records_n'])}, truncated "
          f"{shown(row['bindingdb_rows_truncated_n'])}; cross-source duplicate pairs "
          f"{shown(row['cross_source_potential_duplicate_pairs_n'])} "
          f"(`{row['cross_source_duplicate_status']}`).")
    target_report='\n'.join(target_lines)
    report=f'''# TLS V3 / P0.8 — локальная проверка экспериментальных источников

**Статус:** {technical}. **Научная приёмка:** NOT_APPROVED. **Слияние:** запрещено.

- Панель: {panel}; проверено генов {len(panel_rows)}. Сеть: 0 запросов.
- ChEMBL: {status['chembl']}, release {release_chembl or 'НЕ ПРЕДОСТАВЛЕН'}; строк {len(chem_rows)}; предварительных точных измерений с проверенной метаинформацией assay: {summary['chembl_assay_metadata_eligible_n']}.
- BindingDB: {status['bindingdb']}, release {release_bindingdb or 'НЕ ПРЕДОСТАВЛЕН'}; строк {len(bdb_rows)}; наблюдений одновременно с однозначной цепью SwissProt и подтверждённым видом Homo sapiens: {summary['bindingdb_single_chain_human_observations_n']}.
- Связь с 65 PYZ, докинг, причинность и селективность **не устанавливались**.

## Показатели по мишеням

{target_report}

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
    p.add_argument('--chembl-publisher-sha256',default='')
    p.add_argument('--bindingdb-publisher-sha256',default='')
    p.add_argument('--chembl-license-id',default='')
    p.add_argument('--bindingdb-license-id',default='')
    p.add_argument('--chembl-license-url',default='')
    p.add_argument('--bindingdb-license-url',default='')
    p.add_argument('--panel',choices=('focus6','sample40'),default='focus6')
    p.add_argument('--max-rows-per-target',type=int,default=1000)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    run(p06_zip=a.p06_zip,p0_zip=a.p0_zip,p07_zip=a.p07_zip,
        chembl_sqlite=a.chembl_sqlite,bindingdb_tsv=a.bindingdb_tsv,
        release_chembl=a.chembl_release,release_bindingdb=a.bindingdb_release,
        chembl_publisher_sha256=a.chembl_publisher_sha256,
        bindingdb_publisher_sha256=a.bindingdb_publisher_sha256,
        chembl_license_id=a.chembl_license_id,bindingdb_license_id=a.bindingdb_license_id,
        chembl_license_url=a.chembl_license_url,bindingdb_license_url=a.bindingdb_license_url,
        panel=a.panel,max_rows_per_target=a.max_rows_per_target,out=a.out)

if __name__=='__main__':main()
