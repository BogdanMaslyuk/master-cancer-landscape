"""Synthetic AC-02 regression tests: no user's DepMap data needed."""
from __future__ import annotations

import hashlib
import json
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_pyz_e0_ac02 import (  # noqa: E402
    CONTRACT, LAYERS, axis_quality, canon, inspect_layer, pick_indices,
    recorded_hashes, safe_path, sample_gene_indices, sha256_file, verify_values,
)


def test_contract_and_scope():
    assert CONTRACT.startswith("mcl-e0-ac02")
    assert LAYERS == ("gene_effect", "expression", "copy_number")


def test_traversal_and_symlink_refused(tmp_path):
    with pytest.raises(ValueError):
        safe_path(tmp_path, "../outside.npy")
    with pytest.raises(ValueError):
        safe_path(tmp_path, str(tmp_path / "absolute.npy"))
    file = tmp_path / "file.txt"
    file.write_text("scientific")
    link = tmp_path / "alias"
    try:
        link.symlink_to(file)
    except OSError:
        pytest.skip("Symlink creation unavailable")
    with pytest.raises(ValueError):
        safe_path(tmp_path, "alias")


def test_axis_duplicate_and_empty():
    assert axis_quality(["model1", "model1"], ["CDK4"]) == ["DUPLICATED_MODEL_ID"]
    assert axis_quality(["m1"], ["cdk4", "CDK4"]) == ["DUPLICATED_GENE_SYMBOL"]
    assert "BLANK_AXIS_KEY" in axis_quality([""], ["EGFR"])
    assert "EMPTY_AXIS" in axis_quality([], [])


def test_deterministic_sampling():
    assert pick_indices(1) == [0]
    assert pick_indices(0) == []
    genes = [f"G{i}" for i in range(38)] + ["CDK4", "TP53"]
    a = sample_gene_indices(genes)
    assert a == sample_gene_indices(genes)
    assert 38 in a and 39 in a
    assert max(a) < len(genes)


def test_values_float32_and_nan_equal():
    models = ["ACH-001", "ACH-002", "ACH-003"]
    genes = ["CDK4", "EGFR"]
    source = pd.DataFrame({"model_id": models,
                           "CDK4": [-0.1, float("nan"), -0.8],
                           "EGFR": [0.33, -0.44, 0.12]})
    matrix = np.asarray(source[["CDK4", "EGFR"]].to_numpy(dtype=np.float32).T)
    match, count = verify_values(matrix, source, genes, models, [0, 1])
    assert (match, count) == (6, 6)
    matrix[1, 2] = 5
    match, count = verify_values(matrix, source, genes, models, [0, 1])
    assert (match, count) == (5, 6)


def test_model_order_is_critical():
    source = pd.DataFrame({"model_id": ["ACH-002", "ACH-001"], "CDK4": [1, 2]})
    matrix = np.asarray([[1, 2]], dtype=np.float32)
    assert verify_values(matrix, source, ["CDK4"], ["ACH-001", "ACH-002"], [0]) == (0, -1)


def test_historical_sha_file_parsing_and_tampering(tmp_path):
    contents = b"scientific input v1"
    p = tmp_path / "model.npy"
    p.write_bytes(contents)
    digest = hashlib.sha256(contents).hexdigest()
    record = tmp_path / "file_checksums.sha256"
    record.write_text(f"{digest}  data/runtime/explorer/model_layers/gene_effect.npy\n")
    expected = recorded_hashes(record)
    assert expected["data/runtime/explorer/model_layers/gene_effect.npy"] == sha256_file(p)
    p.write_bytes(b"tampered")
    assert expected["data/runtime/explorer/model_layers/gene_effect.npy"] != sha256_file(p)


def test_historical_sha_rejects_traversal(tmp_path):
    record = tmp_path / "sha"
    record.write_text("a" * 64 + "  ../secret.dat\n")
    with pytest.raises(ValueError):
        recorded_hashes(record)


def make_layer(tmp_path, *, bad_shape=False, bad_orientation=False, source=True):
    folder = tmp_path / "data/runtime/explorer/model_layers"
    folder.mkdir(parents=True)
    processed = tmp_path / "data/processed"
    processed.mkdir(parents=True)
    matrix = np.asarray([[1.0, np.nan], [2.5, -3.0]], dtype=np.float32)
    if bad_shape:
        matrix = matrix[:1, :]
    np.save(folder / "gene_effect.npy", matrix)
    meta = {"schema_version": "1.0", "layer": "gene_effect",
            "orientation": "model_by_gene" if bad_orientation else "gene_by_model",
            "dtype": "float32", "models_n": 2, "genes_n": 2,
            "models": ["ACH-001", "ACH-002"], "genes": ["CDK4", "EGFR"],
            "source_file": "depmap_model_gene_effect.parquet", "matrix_file": "gene_effect.npy"}
    (folder / "gene_effect.json").write_text(json.dumps(meta))
    if source:
        (processed / "depmap_model_gene_effect.parquet").write_bytes(b"synthetic-mocked-parquet")
    return matrix


def test_inspect_layer_missing_source_is_not_claimed_verified(tmp_path):
    make_layer(tmp_path, source=False)
    issues = []
    result = inspect_layer(tmp_path, "gene_effect", lambda *a: issues.append(a), [], baseline={})
    assert result["axis_status"] == "AXES_AND_SHAPE_CHECKED"
    assert result["source_parquet_status"] == "MISSING"
    assert any(i[2] == "MISSING_PROCESSED_MATRIX" and i[0] == "REVIEW" for i in issues)


def test_inspect_layer_wrong_orientation_blocks(tmp_path):
    make_layer(tmp_path, source=False, bad_orientation=True)
    issues = []
    inspect_layer(tmp_path, "gene_effect", lambda *a: issues.append(a), [], baseline={})
    assert any(i[2] == "WRONG_LAYER_CONTRACT" and i[0] == "BLOCKER" for i in issues)


def test_inspect_layer_bad_shape_blocks(tmp_path):
    make_layer(tmp_path, source=False, bad_shape=True)
    issues = []
    inspect_layer(tmp_path, "gene_effect", lambda *a: issues.append(a), [], baseline={})
    assert any(i[2] == "ARRAY_SHAPE_OR_DTYPE" and i[0] == "BLOCKER" for i in issues)


def test_inspect_layer_parquet_values_and_axis_crosscheck(tmp_path, monkeypatch):
    source = pd.DataFrame({"model_id": ["ACH-001", "ACH-002"], "CDK4": [1, np.nan],
                           "EGFR": [2.5, -3.0]})
    make_layer(tmp_path)
    fakepyarrow = types.ModuleType("pyarrow")
    fakepyarrow.__path__ = []
    fakepq = types.ModuleType("pyarrow.parquet")
    fakepq.ParquetFile = lambda unused: types.SimpleNamespace(
        schema_arrow=types.SimpleNamespace(names=list(source.columns)))
    fakepyarrow.parquet = fakepq
    monkeypatch.setitem(sys.modules, "pyarrow", fakepyarrow)
    monkeypatch.setitem(sys.modules, "pyarrow.parquet", fakepq)
    monkeypatch.setattr(pd, "read_parquet", lambda _, columns: source[columns])
    issues = []
    result = inspect_layer(tmp_path, "gene_effect", lambda *a: issues.append(a), [], baseline={})
    assert result["source_parquet_status"] == "ORDER_AND_SAMPLE_MATCH"
    assert result["sample_compared"] == result["sample_matched"] == 4
    assert not [i for i in issues if i[0] == "BLOCKER"]
    source.loc[1, "EGFR"] = 99
    issues = []
    result = inspect_layer(tmp_path, "gene_effect", lambda *a: issues.append(a), [], baseline={})
    assert any(i[2] == "SAMPLED_VALUES_MISMATCH" for i in issues)


def test_powershell_utf8_bom_and_no_science_rebuild():
    script = Path(__file__).resolve().parents[1] / "scripts/run-pyz-e0-ac02.ps1"
    raw = script.read_bytes()
    assert raw.startswith(b"\xef\xbb\xbf")
    text = raw.decode("utf-8-sig")
    assert "build-explorer.ps1" not in text
    assert "run-pyz-e0-ac02.ps1" not in text  # no recursive launch
    assert "Compress-Archive" in text


def synthetic_science(tmp_path, monkeypatch):
    """Small in-memory Parquet reader mock + real NPY/sidecars, no pyarrow dependency."""
    from audit_pyz_e0_ac02 import TABLES, BASE
    index = tmp_path / BASE
    layers = index / "model_layers"
    processed = tmp_path / "data/processed"
    layers.mkdir(parents=True)
    processed.mkdir(parents=True)
    frames = {}
    for name, columns in TABLES.items():
        if name == "gene_catalog":
            data = {column: (['CDK4','EGFR'] if column == 'gene_symbol' else ['[]','[]']) for column in columns}
        elif name == "gene_context_metrics":
            data = {'gene_symbol':['CDK4','EGFR'], 'comparison_id':['C-1','C-1'], 'delta_gene_effect':[0.1,0.3]}
        elif name == "gene_annotations":
            data = {'gene_symbol':['CDK4','EGFR'], 'annotation_type':['pathway','pathway'],
                    'annotation_id':['p1','p2'], 'source':['GO','Reactome'],
                    'source_id':['GO:0001','R-HSA-001']}
        else:
            data = {col: ['CDK4','EGFR'] if col=='gene_symbol' else
                    [1,1] if col=='dependency_models_n' else
                    [2,2] if col=='dependency_models_total_n' else
                    [0.5,0.5] if col=='dependency_fraction' else ['checked','checked'] for col in columns}
        frame = pd.DataFrame(data)
        p = index / f"{name}.parquet"
        p.write_bytes(b'parquet-stub')
        frames[str(p)] = frame
    m = {'index_contract':'mcl-gene-explorer-runtime-v1','schema_version':'1.0',
         'runtime_root':'data/runtime/explorer','genes_n':2,'context_rows_n':2,
         'built_at':'2026-11-16T00:00:00Z','source':'synthetic','builder':'test'}
    (index/'manifest.json').write_text(json.dumps(m))
    (index/'gene_dependency_summary.json').write_text(json.dumps({'dependency_call_method':'probability_threshold'}))
    (processed/'depmap_model_multiomics_manifest.json').write_text(json.dumps({'depmap_release':'TEST_ONLY'}))
    for name in ('gene_effect','expression','copy_number'):
        matrix = np.array([[1.0, np.nan],[2.5, -3.0]], dtype=np.float32)
        np.save(layers/f'{name}.npy',matrix)
        (layers/f'{name}.json').write_text(json.dumps({
           'schema_version':'1.0','layer':name,'orientation':'gene_by_model','dtype':'float32',
           'models_n':2,'genes_n':2,'models':['ACH-001','ACH-002'],
           'genes':['CDK4','EGFR'],'source_file':f'depmap_model_{name}.parquet',
           'matrix_file':f'{name}.npy'}))
        p=processed/f'depmap_model_{name}.parquet'
        p.write_bytes(b'parquet-stub')
        frames[str(p)] = pd.DataFrame({'model_id':['ACH-001','ACH-002'],
                                        'CDK4':[1.0,np.nan], 'EGFR':[2.5,-3.0]})
    fakepyarrow = types.ModuleType('pyarrow')
    fakepyarrow.__path__=[]
    fakepq=types.ModuleType('pyarrow.parquet')
    def pf(path):
        frame = frames[str(path)]
        return types.SimpleNamespace(schema_arrow=types.SimpleNamespace(names=list(frame.columns)),
                                     metadata=types.SimpleNamespace(num_rows=len(frame)))
    fakepq.ParquetFile=pf
    fakepyarrow.parquet=fakepq
    monkeypatch.setitem(sys.modules,'pyarrow',fakepyarrow)
    monkeypatch.setitem(sys.modules,'pyarrow.parquet',fakepq)
    monkeypatch.setattr(pd,'read_parquet',lambda p,columns:frames[str(p)][columns].copy())
    sha = tmp_path/'baseline.sha256'
    from audit_pyz_e0_ac02 import TABLES, LAYERS
    names = ([f'data/runtime/explorer/{name}.parquet' for name in TABLES] +
             [f'data/runtime/explorer/model_layers/{name}.npy' for name in LAYERS])
    sha.write_text(''.join(sha256_file(tmp_path/name)+'  '+name+'\n' for name in names))
    return sha, frames


def test_full_synthetic_ac02_report_has_evidence_but_not_approval(tmp_path,monkeypatch):
    from audit_pyz_e0_ac02 import main
    baseline, _ = synthetic_science(tmp_path,monkeypatch)
    outputs = tmp_path.parent / (tmp_path.name + '_reports')
    assert main(['--root',str(tmp_path),'--out',str(outputs),'--baseline-hashes',str(baseline)]) == 0
    result = json.loads((outputs/'AC02_REPORT.json').read_text(encoding='utf-8'))
    assert result['seven_sources_found']==7
    assert result['source_hashes_unchanged']==7
    assert result['baseline_matches']==7
    assert result['sidecars_unchanged']==6
    assert result['sidecar_present']==6
    assert result['sidecar_missing']==0
    assert result['technical_status']=='REVIEW_REQUIRED'
    assert result['ac02_scientific_acceptance']=='NOT_APPROVED'
    assert result['merge_permitted'] is False
    assert all(x['source_parquet_status']=='ORDER_AND_SAMPLE_MATCH' for x in result['layers'])
    assert sorted(p.name for p in outputs.iterdir())==sorted([
        'AC02_REPORT.json','AC02_SOURCES.tsv','AC02_SIDECARS.tsv','AC02_LAYERS.tsv',
        'AC02_FINDINGS.tsv','AC02_REVIEW_RU.md'])
    with pytest.raises(FileExistsError):
        main(['--root',str(tmp_path),'--out',str(outputs),'--baseline-hashes',str(baseline)])


def test_full_synthetic_detects_model_axis_swap(tmp_path,monkeypatch):
    from audit_pyz_e0_ac02 import audit
    baseline, frames = synthetic_science(tmp_path,monkeypatch)
    key=str(tmp_path/'data/processed/depmap_model_gene_effect.parquet')
    frames[key] = frames[key].iloc[::-1].reset_index(drop=True)
    report, _, issues, _, _ = audit(tmp_path,baseline)
    assert report['technical_status']=='BLOCKED'
    assert any(i['code']=='MODEL_AXIS_OR_GENE_COLUMN_MISMATCH' for i in issues)


def test_full_synthetic_detects_mutated_scientific_cell(tmp_path,monkeypatch):
    from audit_pyz_e0_ac02 import audit
    baseline, frames = synthetic_science(tmp_path,monkeypatch)
    key=str(tmp_path/'data/processed/depmap_model_copy_number.parquet')
    frames[key].loc[1,'EGFR']=42.0
    report, _, issues, _, _ = audit(tmp_path,baseline)
    assert report['technical_status']=='BLOCKED'
    assert any(i['code']=='SAMPLED_VALUES_MISMATCH' for i in issues)


def test_annotations_same_coarse_key_distinct_sources_are_not_duplicates(tmp_path, monkeypatch):
    """One functional label may be justified by different reference terms."""
    from audit_pyz_e0_ac02 import audit
    baseline, frames = synthetic_science(tmp_path, monkeypatch)
    table = str(tmp_path / 'data/runtime/explorer/gene_annotations.parquet')
    base = frames[table].copy()
    # The same gene/type/class is supported by two distinct source identifiers.
    base.loc[0, 'annotation_id'] = 'functional:kinase'
    repeated = base.iloc[[0]].copy()
    repeated.loc[:, 'source'] = 'Reactome'
    repeated.loc[:, 'source_id'] = 'R-HSA-003'
    frames[table] = pd.concat([base, repeated], ignore_index=True)
    result, _, issues, _, _ = audit(tmp_path, baseline)
    anno = next(x for x in result['tables'] if x['table'] == 'gene_annotations')
    assert anno['duplicate_keys'] == 0
    assert anno['coarse_annotation_repeats'] == 1
    assert not any(x['code'] == 'DUPLICATE_ANNOTATION_EVIDENCE' for x in issues)


def test_annotations_same_full_provenance_key_is_review_not_false_pass(tmp_path, monkeypatch):
    from audit_pyz_e0_ac02 import audit
    baseline, frames = synthetic_science(tmp_path, monkeypatch)
    table = str(tmp_path / 'data/runtime/explorer/gene_annotations.parquet')
    frames[table] = pd.concat([frames[table], frames[table].iloc[[0]]], ignore_index=True)
    result, _, issues, _, _ = audit(tmp_path, baseline)
    anno = next(x for x in result['tables'] if x['table'] == 'gene_annotations')
    assert anno['duplicate_keys'] == 1
    assert any(x['code'] == 'DUPLICATE_ANNOTATION_EVIDENCE' and x['severity'] == 'REVIEW' for x in issues)


def test_missing_multiomics_manifest_is_not_reported_as_modified(tmp_path, monkeypatch):
    from audit_pyz_e0_ac02 import audit
    baseline, _ = synthetic_science(tmp_path, monkeypatch)
    (tmp_path / 'data/processed/depmap_model_multiomics_manifest.json').unlink()
    result, _, issues, _, sidecars = audit(tmp_path, baseline)
    assert result['sidecar_present'] == 5
    assert result['sidecar_missing'] == 1
    assert result['sidecars_unchanged'] == 5
    assert any(x['code'] == 'MANIFEST_MISSING' for x in issues)
    assert not any(x['code'] == 'MANIFEST_UNAVAILABLE' and x['component'] == 'multiomics_manifest' for x in issues)
    assert any(not x['exists'] for x in sidecars)
