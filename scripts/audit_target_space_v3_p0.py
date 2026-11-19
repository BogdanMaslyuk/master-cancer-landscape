#!/usr/bin/env python3
"""MCL TLS V3/P0: isolated, offline, READ-ONLY target coverage inventory.

This module counts only what can be established from locally existing snapshots.
It never treats a missing ligand outside the old 13-target catalog as no activity.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from typing import Any

import pandas as pd

VERSION = 'mcl-tls-v3-p0-inventory-v1'
FILES = {
    'gene_catalog': 'data/runtime/explorer/gene_catalog.parquet',
    'gene_reference': 'data/processed/gene_explorer/gene_reference.parquet',
    'protein_registry': 'data/processed/target_registry/protein_targets.parquet',
    'protein_registry_manifest': 'data/processed/target_registry/protein_target_registry_manifest.json',
    'old_target_map': 'data/runtime/target_ligand_space/target_map.tsv',
    'old_ligand_catalog': 'data/runtime/target_ligand_space/ligand_catalog.parquet',
    'old_ligand_measurements': 'data/runtime/target_ligand_space/ligand_measurements.parquet',
    'old_tls_manifest': 'data/runtime/target_ligand_space/target_ligand_space_manifest.json',
    'old_pyz_matrix': 'data/runtime/pyz_target_evidence/pyz_target_evidence_matrix.parquet',
    'own_compounds': 'config/own_compounds.tsv',
    'candidate_v2': 'data/runtime/pharmacology/candidate_hypotheses_v2.parquet',
    'gene_dependency_summary': 'data/runtime/explorer/gene_dependency_summary.parquet',
    'gene_context_metrics': 'data/runtime/explorer/gene_context_metrics.parquet',
    'expression_axes': 'data/runtime/explorer/model_layers/expression.json',
    'copy_number_axes': 'data/runtime/explorer/model_layers/copy_number.json',
    'gene_effect_axes': 'data/runtime/explorer/model_layers/gene_effect.json',
    'depmap_input_inventory': 'data/processed/depmap_input_inventory.tsv',
}
GENE_RE = re.compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]*$')
ENDPOINTS = {'Ki', 'Kd', 'IC50'}
OUTPUTS = ('target_universe.parquet', 'target_universe.tsv', 'ligand_coverage.tsv',
           'biology_coverage.tsv', 'mapping_issues.tsv', 'filter_loss_audit.tsv',
           'source_manifest.json', 'coverage_summary.json', 'TLS_V3_P0_REPORT_RU.md', 'QC.tsv')


def hash_file(path: Path, chunk: int = 4 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for buf in iter(lambda: stream.read(chunk), b''):
            h.update(buf)
    return h.hexdigest()


def text(value: Any) -> str:
    if value is None:
        return ''
    try:
        if pd.isna(value):
            return ''
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def gene(value: Any) -> str:
    v = text(value).upper()
    return v if GENE_RE.fullmatch(v) else ''


def num(value: Any) -> float | None:
    try:
        x = float(value)
        return x if pd.notna(x) and x != float('inf') and x != float('-inf') else None
    except (TypeError, ValueError):
        return None


def source_path(key: str, roots: list[Path]) -> Path | None:
    relative = Path(FILES[key])
    if relative.is_absolute() or '..' in relative.parts:
        raise ValueError('Unsafe source path')
    for root in roots:
        p = root / relative
        if p.is_file():
            return p
    return None


def load_frame(path: Path | None, required: tuple[str, ...] = ()) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame()
    if path.suffix == '.parquet':
        frame = pd.read_parquet(path)
    elif path.suffix == '.tsv':
        frame = pd.read_csv(path, sep='\t', dtype=str, low_memory=False)
    else:
        raise ValueError('Unsupported tabular input')
    missing = set(required) - set(frame.columns)
    if missing:
        raise ValueError(f'{path.name}: missing {sorted(missing)}')
    return frame


def axes_genes(path: Path | None) -> set[str]:
    if path is None:
        return set()
    val = json.loads(path.read_text(encoding='utf-8-sig'))
    seq = val.get('genes')
    if not isinstance(seq, list):
        raise ValueError(f'{path.name}: genes axis missing')
    cleaned = [gene(x) for x in seq]
    if any(not x for x in cleaned) or len(set(cleaned)) != len(cleaned):
        raise ValueError(f'{path.name}: invalid gene axis')
    return set(cleaned)


def validate_output(out: Path, roots: list[Path]) -> None:
    out = out.resolve()
    for root in roots:
        if out == root or root in out.parents or out in root.parents:
            raise ValueError('Audit output must be isolated from scientific roots')
    if out.exists():
        raise FileExistsError(f'Output folder exists: {out}')


def evaluate_coverage(rows: pd.DataFrame, measurements: pd.DataFrame,
                      known_drugs: pd.DataFrame) -> pd.DataFrame:
    """Keep absent/unknown distinct and never synthesize target potency."""
    out = rows.copy()
    out['T0'] = True
    out['T1'] = out['protein_unique_reviewed_human'].astype(bool)
    out['T2'] = out['T1'] & out['has_direct_ligand'].astype(bool)
    out['T3'] = out['T2'] & out['has_crispr'].astype(bool)
    out['T4'] = out['T3'] & out['has_measured_known_drug_context'].astype(bool)
    out['coverage_tier'] = out.apply(lambda r: 'T4' if r['T4'] else 'T3' if r['T3']
                                     else 'T2' if r['T2'] else 'T1' if r['T1'] else 'T0', axis=1)
    return out


def audit(science_root: Path, reference_root: Path, out: Path, *,
          release_snapshot: str | None = None) -> dict[str, Any]:
    science_root, reference_root, out = (x.resolve() for x in (science_root, reference_root, out))
    if not science_root.is_dir() or not reference_root.is_dir():
        raise FileNotFoundError('Scientific and reference roots must already exist')
    roots = list(dict.fromkeys((science_root, reference_root)))
    validate_output(out, roots)
    paths = {key: source_path(key, roots) for key in FILES}
    qc = []
    source_manifest = []
    for key, p in paths.items():
        if p is None:
            qc.append({'severity': 'REVIEW', 'component': key, 'code': 'SOURCE_MISSING',
                       'detail': 'No existing source in either root'})
        source_manifest.append({'name': key, 'path': str(p) if p else None,
                                'present': p is not None, 'bytes': p.stat().st_size if p else None,
                                'sha256': hash_file(p) if p else None})
    print(f'[TLS V3/P0] Sources: {sum(p is not None for p in paths.values())}/{len(paths)} found; reading local tables', flush=True)
    # No inference of unknown protein mapping: local registry is restricted to
    # previous pharmacology targets and has no complete proteome coverage.
    catalog = load_frame(paths['gene_catalog'], ('gene_symbol',))
    if catalog.empty:
        raise ValueError('gene_catalog is required for P0 target universe')
    genes = {gene(x) for x in catalog['gene_symbol']}
    genes.discard('')
    if not genes:
        raise ValueError('No valid genes in catalog')
    issues = []
    for i, raw in enumerate(catalog['gene_symbol']):
        if not gene(raw):
            issues.append({'source':'gene_catalog', 'row':i+1, 'gene_symbol':text(raw),
                           'issue':'INVALID_GENE_IDENTIFIER'})
    reg = load_frame(paths['protein_registry'])
    old = load_frame(paths['old_target_map'])
    lig = load_frame(paths['old_ligand_catalog'])
    ms = load_frame(paths['old_ligand_measurements'])
    deps = load_frame(paths['gene_dependency_summary'])
    candidates = load_frame(paths['candidate_v2'])
    for req, frame in [('protein_registry',reg),('old_target_map',old),('old_ligand_catalog',lig),
                       ('old_ligand_measurements',ms),('gene_dependency_summary',deps),('candidate_v2',candidates)]:
        if not frame.empty and 'target_gene' not in frame.columns and req in {'protein_registry','old_target_map','old_ligand_catalog','old_ligand_measurements','candidate_v2'}:
            qc.append({'severity':'BLOCKER','component':req,'code':'MISSING_TARGET_GENE','detail':'Cannot join target keys'})
    old_genes = set(gene(x) for x in old.get('target_gene', pd.Series(dtype=str)));old_genes.discard('')
    # Existing panel may be incomplete in local snapshot. Never assume 13 were
    # checked unless they are actually present.
    pyz = load_frame(paths['own_compounds'])
    old_matrix = load_frame(paths['old_pyz_matrix'])
    if len(old_genes) != 13:
        qc.append({'severity':'REVIEW', 'component':'old_target_map', 'code':'PANEL_13_NOT_CONFIRMED',
                   'detail': f'Found {len(old_genes)} distinct targets (expected 13)'})
    if len(pyz) != 65 or ('own_compound_id' in pyz and pyz['own_compound_id'].nunique() != 65):
        qc.append({'severity':'REVIEW','component':'own_compounds','code':'PYZ_65_NOT_CONFIRMED',
                   'detail':f'Found {len(pyz)} rows (expected 65)'})
    if len(old_matrix) != 845:
        qc.append({'severity':'REVIEW','component':'old_pyz_matrix','code':'MATRIX_845_NOT_CONFIRMED',
                   'detail':f'Found {len(old_matrix)} rows (expected 845)'})
    # Candidate paths: A = genes with measured CRISPR; B = locally measured
    # direct ligand (only the previously scoped reference library); C = old panel.
    crispr_genes = axes_genes(paths['gene_effect_axes']) if paths['gene_effect_axes'] else set()
    expr_genes = axes_genes(paths['expression_axes']) if paths['expression_axes'] else set()
    cn_genes = axes_genes(paths['copy_number_axes']) if paths['copy_number_axes'] else set()
    summary_genes = {gene(x) for x in deps.get('gene_symbol', pd.Series(dtype=str))}
    summary_genes.discard('')
    biological = crispr_genes | summary_genes
    ligand_counts: dict[str,dict[str,int]] = {}
    for x in lig.to_dict('records'):
        g = gene(x.get('target_gene'))
        if not g:
            issues.append({'source':'ligand_catalog','row':None,'gene_symbol':text(x.get('target_gene')),
                           'issue':'INVALID_TARGET_GENE'})
            continue
        ligand_counts.setdefault(g, {'catalog_rows':0,'valid_structures':0})['catalog_rows']+=1
        if text(x.get('inchikey')) and text(x.get('canonical_smiles')):
            ligand_counts[g]['valid_structures']+=1
    # Binding measurements belong to known ligands, not to PYZ. Preserve endpoint
    # identities and censoring; count direct eligible records with exact target.
    meas_stats: dict[str,dict[str,int]] = {}
    for x in ms.to_dict('records'):
        g=gene(x.get('target_gene'))
        if not g:continue
        z=meas_stats.setdefault(g, {'records':0,'Ki':0,'Kd':0,'IC50':0,'eligible_direct':0,'censored':0})
        z['records']+=1
        e=text(x.get('endpoint_type') or x.get('endpoint') or x.get('standard_type'))
        relation=text(x.get('relation') or x.get('standard_relation') or '=')
        if e in ENDPOINTS:z[e]+=1
        if relation not in {'=',''}:z['censored']+=1
        # No positive classification when reported measurement is incomplete.
        value=num(x.get('value_nm') if 'value_nm' in x else x.get('standard_value'))
        if e in ENDPOINTS and value is not None and value>0 and relation in {'=','<','<=','>','>='}:
            z['eligible_direct']+=1
    reg_by:dict[str,list[dict]]= {}
    for i,x in enumerate(reg.to_dict('records')):
        g=gene(x.get('target_gene'))
        if g:reg_by.setdefault(g,[]).append(x)
        else:issues.append({'source':'protein_registry','row':i+1,'gene_symbol':text(x.get('target_gene')),'issue':'INVALID_TARGET_GENE'})
    candidate_by: dict[str, list[dict]]={}
    for x in candidates.to_dict('records'):
        g=gene(x.get('target_gene'))
        if g:candidate_by.setdefault(g,[]).append(x)
    # The intended T4 describes measured *known-drug* context, not PYZ.
    def known_active(records:list[dict]) -> bool:
        return any(num(x.get('active_models_n')) is not None and num(x.get('active_models_n'))>0
                   for x in records)
    universe = genes | biological | set(ligand_counts) | old_genes
    out_rows=[]; ligand_rows=[]; biology_rows=[]
    for g in sorted(universe):
        r=reg_by.get(g,[])
        good=[x for x in r if text(x.get('protein_mapping_status'))=='unique_swissprot'
              and text(x.get('uniprot_primary_accession'))]
        unique = len(good)==1 and len(r)==1
        if len(r)>1 or (r and not unique):
            issues.append({'source':'protein_registry','row':None,'gene_symbol':g,
                           'issue':'AMBIGUOUS_OR_UNRESOLVED_PROTEIN_MAPPING'})
        accession=text(good[0]['uniprot_primary_accession']) if unique else ''
        stats=meas_stats.get(g, {'records':0,'Ki':0,'Kd':0,'IC50':0,'eligible_direct':0,'censored':0})
        liginfo=ligand_counts.get(g, {'catalog_rows':0,'valid_structures':0})
        present_ligand = stats['eligible_direct']>0 and liginfo['valid_structures']>0
        known=candidate_by.get(g,[])
        pathA=g in biological;pathB=g in meas_stats and stats['eligible_direct']>0;pathC=g in old_genes
        knownctx=known_active(known)
        panel_coverage= 'CHECKED_OLD_PANEL' if g in old_genes else 'NOT_ASSESSED_OUTSIDE_OLD_PANEL'
        base={'gene_symbol':g,'hgnc_id':text(good[0].get('hgnc_id')) if unique else '',
              'uniprot_accession':accession,'protein_mapping_status':
               'unique_swissprot' if unique else 'unresolved_or_not_in_local_registry',
              'protein_unique_reviewed_human':unique,'path_A':pathA,'path_B':pathB,'path_C':pathC,
              'has_crispr':g in crispr_genes,'has_expression':g in expr_genes,
              'has_copy_number':g in cn_genes,'has_direct_ligand':present_ligand,
              'has_measured_known_drug_context':knownctx,
              'ligand_evidence_scope':panel_coverage,
              'missing_ligand_is_negative':False,
              'candidate_v2_any_status':bool(known),
              'candidate_v2_priority':any(text(x.get('priority_status'))=='priority_for_in_vitro' for x in known)}
        out_rows.append(base)
        ligand_rows.append({'gene_symbol':g,'source_scope':panel_coverage,'measurement_rows':stats['records'],
                            'eligible_direct_measurements':stats['eligible_direct'],
                            'ligand_catalog_rows':liginfo['catalog_rows'],
                            'ligand_structures_usable':liginfo['valid_structures'],
                            'Ki_records':stats['Ki'],'Kd_records':stats['Kd'],
                            'IC50_records':stats['IC50'],'censored_records':stats['censored'],
                            'not_assessed_outside_local_scope':g not in old_genes})
        biology_rows.append({'gene_symbol':g,'crispr_available':g in crispr_genes,
                             'dependency_summary_available':g in summary_genes,
                             'expression_available':g in expr_genes,'copy_number_available':g in cn_genes,
                             'candidate_v2_context_available':bool(known),
                             'known_drug_active_context':knownctx})
    print(f'[TLS V3/P0] Universe rows: {len(out_rows)}; building coverage and filter audit', flush=True)
    universe_df=evaluate_coverage(pd.DataFrame(out_rows),ms,candidates)
    # Audit the non-causal selection funnel separately: observed exclusion vs.
    # inferred selection bias, never assuming unseen ligands are absent.
    losses=[]
    # Index once: filtering a 20,000-row dataframe separately for each gene
    # is quadratic and can make a real local run unnecessarily expensive.
    indexed_universe=universe_df.set_index('gene_symbol',verify_integrity=True)
    for g in sorted(universe):
        r=indexed_universe.loc[g]
        losses.append({'gene_symbol':g,'in_old_13':bool(r['path_C']),
                       'in_candidate_v2':bool(r['candidate_v2_any_status']),
                       'priority_for_in_vitro':bool(r['candidate_v2_priority']),
                       'known_direct_ligand_locally':bool(r['has_direct_ligand']),
                       'gene_not_eligible_due_to_priority_filter': bool(r['candidate_v2_any_status']) and not bool(r['candidate_v2_priority']),
                       'outside_panel_reason':'PRIORITY_FILTER' if bool(r['candidate_v2_any_status']) and not bool(r['candidate_v2_priority'])
                          else 'COVERAGE_NOT_ASSESSED' if not bool(r['path_C']) else 'INCLUDED'})
    # Explicit cross-root provenance: snapshot source hashes and verify again.
    changes=[]
    for x in source_manifest:
        if x['present']:
            changed=hash_file(Path(x['path'])) != x['sha256']
            if changed:changes.append(x['name'])
    if changes:
        qc.append({'severity':'BLOCKER','component':'source_manifest','code':'SOURCE_CHANGED_DURING_AUDIT',
                   'detail':','.join(changes)})
    # Record release claims, do not certify raw source from local labels.
    release_claims=set()
    inv=paths['depmap_input_inventory']
    if inv:
        for v in load_frame(inv).get('depmap_release',pd.Series(dtype=str)).dropna():
            if text(v):release_claims.add(text(v))
    if release_snapshot and release_claims and release_snapshot not in release_claims:
        qc.append({'severity':'BLOCKER','component':'depmap_release','code':'DEPMAP_RELEASE_CONFLICT',
                   'detail':f'Requested {release_snapshot} but inventory says {sorted(release_claims)}'})
    if len(release_claims)>1:
        qc.append({'severity':'BLOCKER','component':'depmap_release','code':'MULTIPLE_DEPMAP_RELEASES',
                   'detail':','.join(sorted(release_claims))})
    if not release_claims:
        qc.append({'severity':'REVIEW','component':'depmap_release','code':'UNVERIFIED_RELEASE',
                   'detail':'No local DepMap input inventory release'})
    status='BLOCKED' if any(x['severity']=='BLOCKER' for x in qc) else 'REVIEW_REQUIRED' if any(x['severity']=='REVIEW' for x in qc) else 'PASS_TECHNICAL'
    summary={'contract':VERSION,'technical_status':status,'scientific_acceptance':'NOT_APPROVED',
             'merge_permitted':False,'external_coverage':'EXTERNAL_COVERAGE_NOT_ASSESSED',
             'full_chembl_bindingdb_coverage':'NOT_ASSESSED',
             'depmap_release_declared':sorted(release_claims),
             'universe_genes_n':len(genes),'universe_entities_n':len(universe),
             'panel_genes_n':len(old_genes),'own_compounds_n':len(pyz),
             'existing_matrix_rows_n':len(old_matrix),
             'counts':{f'T{i}':int(universe_df[f'T{i}'].sum()) for i in range(5)},
             'with_crispr_n':int(universe_df['has_crispr'].sum()),
             'with_local_direct_ligand_n':int(universe_df['has_direct_ligand'].sum()),
             'outside_panel_not_assessed_n':int((~universe_df['path_C']).sum()),
             'filter_excluded_candidate_genes_n':sum(x['gene_not_eligible_due_to_priority_filter'] for x in losses),
             'recommendation':'P1_NO_AUTOMATIC_GO: review identity gaps and acquire complete licensed ChEMBL/BindingDB source snapshots before setting new target count',
             'read_only':True,'evidence_limitations':[
                 'Local ligand collection is restricted to previous target panel; not entire ChEMBL/BindingDB.',
                 'Missing evidence outside panel is NOT an inactive target.',
                 'CRISPR knockout does not establish target engagement or efficacy of PYZ.',
                 'The local registry covers pharmacology-linked targets, not the entire human proteome.']}
    print(f'[TLS V3/P0] Writing report tables to: {out}', flush=True)
    out.mkdir(parents=True,exist_ok=False)
    universe_df.to_parquet(out/'target_universe.parquet',index=False)
    universe_df.to_csv(out/'target_universe.tsv',sep='\t',index=False)
    pd.DataFrame(ligand_rows).to_csv(out/'ligand_coverage.tsv',sep='\t',index=False)
    pd.DataFrame(biology_rows).to_csv(out/'biology_coverage.tsv',sep='\t',index=False)
    pd.DataFrame(issues,columns=['source','row','gene_symbol','issue']).to_csv(out/'mapping_issues.tsv',sep='\t',index=False)
    pd.DataFrame(losses).to_csv(out/'filter_loss_audit.tsv',sep='\t',index=False)
    pd.DataFrame(qc,columns=['severity','component','code','detail']).to_csv(out/'QC.tsv',sep='\t',index=False)
    (out/'source_manifest.json').write_text(json.dumps({'contract':VERSION,'science_root':str(science_root),
        'reference_root':str(reference_root),'sources':source_manifest,'read_only':True},
        ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    (out/'coverage_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    report=f'''# MCL / TLS V3 P0 — локальная инвентаризация мишеней

**Статус:** {status}. **Научная приёмка:** NOT_APPROVED. Слияние веток запрещено.

## Наблюдения

- Генов в исходном каталоге: **{len(genes)}**.
- Уникальных кандидатов (гены ∪ локальные мишени): **{len(universe)}**.
- Старых мишеней из карты: **{len(old_genes)}**.
- PYZ в реестре: **{len(pyz)}**; строк действующей матрицы: **{len(old_matrix)}**.
- Локальные классы покрытия: {', '.join(f'T{i}: {summary['counts'][f'T{i}']}' for i in range(5))}.
- CRISPR Gene Effect: **{summary['with_crispr_n']}** генов.
- Локальные прямые лигандные данные: **{summary['with_local_direct_ligand_n']}** мишеней.
- Вне старой панели: **{summary['outside_panel_not_assessed_n']}** генов/кандидатов; лигандные данные НЕ оценены.
- Отобраны фильтром Candidate v2 без статуса priority_for_in_vitro: **{summary['filter_excluded_candidate_genes_n']}** кандидатов (если данные доступны).
- Заявленный релиз DepMap в локальной инвентаризации: {', '.join(sorted(release_claims)) or 'НЕ ИЗВЕСТЕН'}.

## Научные ограничения

Локальная библиотека лигандов сформирована прежним узким отбором. Это НЕ полный ChEMBL/BindingDB: нельзя из локальных нулей утверждать отсутствие экспериментальных лигандов для остальных белков. T0–T4 характеризуют **покрытие**, а не вероятность связывания PYZ. Внешний поиск выключен; полное покрытие источников не оценено.

## Следующее решение

Нельзя обоснованно назначить число мишеней для P1 без ревизии первичных источников и разрешения HGNC/UniProt для полного протеома. Не изменять 65×13, PYZ, API, сайт, исходные научные данные. Перед научной приёмкой исследователь изучает QC, source_manifest и проверяет P1-панель вручную.
'''
    (out/'TLS_V3_P0_REPORT_RU.md').write_text(report,encoding='utf-8')
    return summary


def main() -> None:
    parser=argparse.ArgumentParser(description='MCL TLS V3/P0 read-only offline target inventory')
    parser.add_argument('--science-root',type=Path,required=True)
    parser.add_argument('--reference-root',type=Path,default=None)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--release-snapshot',default=None)
    parser.add_argument('--offline',action='store_true',default=True)
    parser.add_argument('--enable-external-pilot',action='store_true')
    parser.add_argument('--max-external-targets',type=int,default=50)
    args=parser.parse_args()
    if args.enable_external_pilot:
        parser.error('External pilot is not implemented in P0 offline runner: no requests were sent')
    summary=audit(args.science_root,args.reference_root or args.science_root,
                  args.out,release_snapshot=args.release_snapshot)
    print(json.dumps({k:summary[k] for k in ('technical_status','universe_entities_n','panel_genes_n',
                                         'scientific_acceptance','merge_permitted')},ensure_ascii=False))


if __name__=='__main__':
    main()