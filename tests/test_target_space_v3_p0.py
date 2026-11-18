"""Offline, synthetic coverage tests: no external queries or scientific writes."""
import json
from pathlib import Path
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_target_space_v3_p0 import (
    FILES, OUTPUTS, audit, axes_genes, evaluate_coverage, gene, hash_file,
    load_frame, num, source_path, text, validate_output
)


def test_gene_normalization():
    assert gene(' cdK4 ') == 'CDK4'


@pytest.mark.parametrize('raw', ['', 'bad symbol', '.', '🙂'])
def test_bad_gene_rejected(raw):
    assert gene(raw) == ''


def test_gene_with_numeric_and_dots():
    assert gene('HLA-DRA') == 'HLA-DRA'
    assert gene('MIR-1.1') == 'MIR-1.1'


def test_null_text():
    assert text(None) == ''
    assert text(float('nan')) == ''


@pytest.mark.parametrize('x', ['', None, 'abc', float('inf'), float('nan')])
def test_unusable_numeric(x):
    assert num(x) is None


def test_numeric_zero_is_not_missing():
    assert num(0) == 0


def test_hash_tracks_bytes(tmp_path):
    p=tmp_path/'source'
    p.write_bytes(b'a')
    one=hash_file(p)
    p.write_bytes(b'b')
    assert one!=hash_file(p)


def test_source_path_preserves_root_order(tmp_path):
    a,b=tmp_path/'a',tmp_path/'b'
    for root in (a,b):
        p=root/FILES['own_compounds'];p.parent.mkdir(parents=True);p.write_text('sample')
    assert source_path('own_compounds',[a,b])==a/FILES['own_compounds']


def test_source_path_missing(tmp_path):
    assert source_path('old_target_map',[tmp_path]) is None


def test_reject_source_output_overlap(tmp_path):
    a=tmp_path/'science';a.mkdir()
    with pytest.raises(ValueError):validate_output(a/'reports',[a])
    with pytest.raises(ValueError):validate_output(tmp_path,[a])


def test_output_not_overwritten(tmp_path):
    s=tmp_path/'s';s.mkdir();o=tmp_path/'out';o.mkdir()
    with pytest.raises(FileExistsError):validate_output(o,[s])


def test_axes_reject_duplicate_genes(tmp_path):
    p=tmp_path/'axis.json';p.write_text(json.dumps({'genes':['CDK4','CDK4']}))
    with pytest.raises(ValueError):axes_genes(p)


def test_axes_reject_bad_shape(tmp_path):
    p=tmp_path/'axis.json';p.write_text(json.dumps({'genes':'CDK4'}))
    with pytest.raises(ValueError):axes_genes(p)


def test_axes_valid(tmp_path):
    p=tmp_path/'axis.json';p.write_text(json.dumps({'genes':['CDK4','EGFR']}))
    assert axes_genes(p)=={'CDK4','EGFR'}


def test_load_frame_missing_returns_empty():
    assert load_frame(None).empty


def test_load_frame_rejects_unknown(tmp_path):
    p=tmp_path/'a.txt';p.write_text('a')
    with pytest.raises(ValueError):load_frame(p)


def test_tiers_are_nested():
    data=pd.DataFrame([{'protein_unique_reviewed_human':True,'has_direct_ligand':True,
                        'has_crispr':True,'has_measured_known_drug_context':True},
                       {'protein_unique_reviewed_human':False,'has_direct_ligand':True,
                        'has_crispr':True,'has_measured_known_drug_context':True}])
    rows=evaluate_coverage(data,pd.DataFrame(),pd.DataFrame())
    assert list(rows.coverage_tier)==['T4','T0']
    assert rows.T4.sum()<=rows.T3.sum()<=rows.T2.sum()<=rows.T1.sum()


@pytest.fixture
def mock_science(tmp_path,monkeypatch):
    science,central=tmp_path/'science',tmp_path/'central'
    science.mkdir();central.mkdir()
    frames={
        'gene_catalog':pd.DataFrame({'gene_symbol':['CDK4','EGFR','TP53']}),
        'protein_registry':pd.DataFrame([
            {'target_gene':'CDK4','protein_mapping_status':'unique_swissprot','uniprot_primary_accession':'P11802'},
            {'target_gene':'EGFR','protein_mapping_status':'multiple_swissprot','uniprot_primary_accession':'P00533'}]),
        'old_target_map':pd.DataFrame({'target_gene':['CDK4']}),
        'old_ligand_catalog':pd.DataFrame({'target_gene':['CDK4'],
             'inchikey':['ABCDEFGHIJ-KLMNOPQRST-U'],'canonical_smiles':['CCO']}),
        'old_ligand_measurements':pd.DataFrame({'target_gene':['CDK4'],'endpoint_type':['Ki'],
             'value_nm':[10],'relation':['=']}),
        'candidate_v2':pd.DataFrame([{'target_gene':'CDK4','priority_status':'priority_for_in_vitro',
                                    'active_models_n':2},
                                   {'target_gene':'EGFR','priority_status':'exploratory_hypothesis',
                                    'active_models_n':0}]),
        'gene_dependency_summary':pd.DataFrame({'gene_symbol':['CDK4','EGFR']}),
        'own_compounds':pd.DataFrame({'own_compound_id':[f'PYZ-{i:03d}' for i in range(1,66)]}),
        'old_pyz_matrix':pd.DataFrame({'own_compound_id':[f'PYZ-{i:03d}' for i in range(1,66)]*13,
                                       'target_gene':['CDK4']*845})
    }
    for key in frames:
        loc=(science if key in ('gene_catalog','old_target_map','own_compounds','old_pyz_matrix') else central)/FILES[key]
        loc.parent.mkdir(parents=True,exist_ok=True)
        if loc.suffix=='.tsv':frames[key].to_csv(loc,sep='\t',index=False)
        else:loc.write_bytes(b'synthetic dummy parquet')
    for key,genes in [('gene_effect_axes',['CDK4','EGFR']),('expression_axes',['CDK4']),
                       ('copy_number_axes',['CDK4'])]:
        p=science/FILES[key];p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(json.dumps({'genes':genes}))
    inv=central/FILES['depmap_input_inventory'];inv.parent.mkdir(parents=True,exist_ok=True)
    inv.write_text('depmap_release\tlogical_name\n26Q1\tgene_effect\n')
    raw=pd.read_parquet
    def fake_read(path,*args,**kwargs):
        for key,frame in frames.items():
            if str(path).replace('\\','/').endswith(FILES[key].replace('\\','/')):
                return frame.copy()
        return raw(path,*args,**kwargs)
    monkeypatch.setattr(pd,'read_parquet',fake_read)
    # No pyarrow in sandbox, use pickle writer only for testing output creation.
    monkeypatch.setattr(pd.DataFrame,'to_parquet',lambda self,path,**kw:Path(path).write_bytes(b'test-only-parquet'))
    return science,central


def test_audit_never_claims_complete_ligand_space(mock_science,tmp_path):
    s,c=mock_science
    r=audit(s,c,tmp_path/'report')
    assert r['full_chembl_bindingdb_coverage']=='NOT_ASSESSED'
    assert r['external_coverage']=='EXTERNAL_COVERAGE_NOT_ASSESSED'


def test_audit_never_auto_approves(mock_science,tmp_path):
    s,c=mock_science;r=audit(s,c,tmp_path/'report')
    assert r['scientific_acceptance']=='NOT_APPROVED' and not r['merge_permitted']


def test_missing_source_exposed_not_fabricated(mock_science,tmp_path):
    s,c=mock_science;r=audit(s,c,tmp_path/'report')
    qc=pd.read_csv(tmp_path/'report'/'QC.tsv',sep='\t')
    assert 'SOURCE_MISSING' in set(qc.code)


def test_t0_t1_t2_t3_t4_counts(mock_science,tmp_path):
    s,c=mock_science;r=audit(s,c,tmp_path/'report')
    assert r['counts']['T4']==1
    assert r['counts']['T0']>=r['counts']['T1']>=r['counts']['T2']>=r['counts']['T3']>=r['counts']['T4']


def test_nonpanel_not_assessed(mock_science,tmp_path):
    s,c=mock_science;audit(s,c,tmp_path/'report')
    frame=pd.read_csv(tmp_path/'report'/'target_universe.tsv',sep='\t')
    assert set(frame.loc[frame.gene_symbol=='EGFR','ligand_evidence_scope'])=={'NOT_ASSESSED_OUTSIDE_OLD_PANEL'}


def test_candidate_nonpriority_filter_shown(mock_science,tmp_path):
    s,c=mock_science;audit(s,c,tmp_path/'report')
    frame=pd.read_csv(tmp_path/'report'/'filter_loss_audit.tsv',sep='\t')
    assert frame.loc[frame.gene_symbol.eq('EGFR'),'gene_not_eligible_due_to_priority_filter'].iloc[0]


def test_ligand_endpoints_not_combined(mock_science,tmp_path):
    s,c=mock_science;audit(s,c,tmp_path/'report')
    frame=pd.read_csv(tmp_path/'report'/'ligand_coverage.tsv',sep='\t')
    row=frame.loc[frame.gene_symbol.eq('CDK4')].iloc[0]
    assert row.Ki_records==1 and row.Kd_records==0 and row.IC50_records==0


def test_sources_unchanged(mock_science,tmp_path):
    s,c=mock_science
    files=[p for root in (s,c) for p in root.rglob('*') if p.is_file()]
    before={str(p):hash_file(p) for p in files}
    audit(s,c,tmp_path/'report')
    assert before=={str(p):hash_file(p) for p in files}


def test_required_ten_outputs(mock_science,tmp_path):
    s,c=mock_science;audit(s,c,tmp_path/'report')
    assert set(OUTPUTS).issubset({p.name for p in (tmp_path/'report').iterdir()})


def test_release_contradiction_blocks(mock_science,tmp_path):
    s,c=mock_science
    r=audit(s,c,tmp_path/'report',release_snapshot='25Q4')
    assert r['technical_status']=='BLOCKED'


def test_fresh_independent_run_reproducible(mock_science,tmp_path):
    s,c=mock_science
    one=audit(s,c,tmp_path/'r1')
    two=audit(s,c,tmp_path/'r2')
    assert one==two
    assert (tmp_path/'r1'/'target_universe.tsv').read_bytes()==(tmp_path/'r2'/'target_universe.tsv').read_bytes()