"""E0 closeout fixtures. No scientific source or actual local server required."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "close_pyz_e0.py"
spec = importlib.util.spec_from_file_location("pyz_e0_closeout", SCRIPT)
closeout = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(closeout)


def test_baseline_hashes_and_no_mutation(tmp_path):
    root = tmp_path / "mcl"
    source = root / "config" / "own_compounds.tsv"
    source.parent.mkdir(parents=True)
    source.write_text("own_compound_id\tstandardized_smiles\nPYZ-035\tCCO\n", encoding="utf-8")
    report = tmp_path / "e0"
    report.mkdir()
    original = closeout.sha256(source)
    (report / "file_checksums.sha256").write_text(
        f"{original}  config/own_compounds.tsv\n", encoding="utf-8")
    check = closeout.verify_source_hashes(root, closeout.read_baseline_hashes(report))
    assert check[0]["status"] == "UNCHANGED"
    assert source.read_bytes().startswith(b"own_compound_id")
    source.write_text(source.read_text() + "PYZ-051\tCCN\n", encoding="utf-8")
    assert closeout.verify_source_hashes(root, closeout.read_baseline_hashes(report))[0]["status"] == "CHANGED"


def test_refuses_path_traversal_in_old_report(tmp_path):
    r = tmp_path / "baseline"
    r.mkdir()
    (r / "file_checksums.sha256").write_text("a" * 64 + "  ../outside\n")
    import pytest
    with pytest.raises(ValueError, match="Unsafe"):
        closeout.read_baseline_hashes(r)


def test_numeric_comparison_correct_missing_and_nearly_equal():
    assert closeout.numeric_equal(None, None)
    assert closeout.numeric_equal(float('nan'), None) is False  # raw NaN must not pretend to be zero.
    assert closeout.numeric_equal(0.395833333, 0.3958333331)
    assert not closeout.numeric_equal(0.3958, 0.4958)
    assert not closeout.numeric_equal(None, 0)


def test_scientific_reconciliation_source_vs_api_and_no_pyz_activity(tmp_path, monkeypatch):
    root = tmp_path / "mcl"
    config = root / "config"
    config.mkdir(parents=True)
    reg = pd.DataFrame({"own_compound_id": ["PYZ-001", "PYZ-035", "PYZ-051", "PYZ-055", "PYZ-065"],
                        "standardized_smiles": ["CCO", "CCN", "CCC", "COC", "CNC"]})
    reg.to_csv(config / "own_compounds.tsv", sep="\t", index=False)
    admet = pd.DataFrame({"own_compound_id": reg.own_compound_id, "qed": [0.1, 0.2, 0.3, 0.4, 0.5],
                          "integrated_risk": ["check"] * 5})
    admet.to_csv(config / "pyz_admet_snapshot.tsv", sep="\t", index=False)
    matrix_dir = root / "data" / "runtime" / "pyz_target_evidence"
    matrix_dir.mkdir(parents=True)
    m = pd.DataFrame([{"own_compound_id": cid, "target_gene": "CDK4",
                       "best_ligand_id": "CHEMBL-REF", "best_tanimoto_morgan_r2_2048": 0.39,
                       "best_lowest_reported_value_nm_across_endpoints": 18.0,
                       "best_lowest_value_endpoint_type": "IC50"} for cid in reg.own_compound_id])
    m.to_csv(matrix_dir / "pyz_target_evidence_matrix.tsv", sep="\t", index=False)
    neighbors_dir = root / "data" / "runtime" / "own_compounds"
    neighbors_dir.mkdir(parents=True)
    pd.DataFrame([{"own_compound_id": "PYZ-035", "target_gene": "CDK4",
                   "rank_in_target": 1, "ligand_id": "CHEMBL-REF",
                   "lowest_value_endpoint_type": "IC50", "lowest_value_relation": ">",
                   "lowest_reported_value_nm_across_endpoints": 18.0,
                   "measurements_n": 1}]).to_csv(neighbors_dir / "target_ligand_space_similarity_top_hits.tsv", sep="\t", index=False)
    def api(url, timeout=20.0):
        cid = url.rsplit('/', 1)[-1]
        ix = list(reg.own_compound_id).index(cid)
        return 200, {"profile_v2": {"compound_id": cid,
            "identity": {"standardized_smiles": reg.iloc[ix].standardized_smiles},
            "properties": {"qed": admet.iloc[ix].qed},
            "admet": {"integrated_risk": "check"},
            "target_hypotheses": {"items": [{"target_gene": "CDK4", "tanimoto": 0.39,
                  "best_ligand_id": "CHEMBL-REF", "reference_endpoint": "IC50",
                  "reference_value_nm": 18.0, "reference_measurement_owner": "known_ligand_not_pyz"}]},
            "scientific_summary": {"own_target_binding": "not_measured"}}}, {}, 1.0
    monkeypatch.setattr(closeout, "request", api)
    report = tmp_path / "report"
    report.mkdir()
    hashes_before = closeout.sha256(matrix_dir / "pyz_target_evidence_matrix.tsv")
    rows, ref, cell, failures = closeout.reconcile_scientific(root, 'http://127.0.0.1:8100', report)
    assert not failures
    assert len(rows) == 5 * 7  # smiles, QED, risk, similarity, ligand, endpoint, measurement.
    assert all(r["status"] == "MATCH" for r in rows)
    assert ref[0]["relation"] == ">"
    assert ref[0]["primary_assay_id"] is None
    assert ref[0]["primary_assay_status"] == "REQUIRES_PRIMARY_SOURCE"
    assert not cell
    assert closeout.sha256(matrix_dir / "pyz_target_evidence_matrix.tsv") == hashes_before
    assert (report / 'source_api_reconciliation.tsv').is_file()


def test_status_never_auto_approves_expert_checks():
    hashes = [{"path": "a", "status": "UNCHANGED"}]
    rows = [{"status": "MATCH"}]
    a = closeout.status_table(hashes, rows, [],
        {"status": "AUTO_PASS_MANUAL_REVIEW_PENDING", "n_screenshots": 9},
        {"status": "PASS"}, [{"primary_assay_id": None}], {})
    assert a['AC-08']['status'] == 'AUTO_PASS'
    assert a['AC-05']['status'] == 'AUTO_CAPTURED_REVIEW_PENDING'
    assert a['AC-07']['status'] == 'NOT_REVIEWED'
    assert a['AC-09']['status'] == 'NOT_APPROVED'
    assert all(row['status'] != 'APPROVED' for row in a.values())


def test_missing_browser_falls_back_to_manual(tmp_path):
    result = closeout.run_ui(tmp_path, "http://127.0.0.1:3100", tmp_path, tmp_path)
    assert result["status"] == "NOT_RUN"
    assert 'unavailable' in result['reason']


def test_complete_closeout_does_not_write_to_science_root(tmp_path, monkeypatch):
    from argparse import Namespace
    root = tmp_path / 'science'
    reg = root / 'config' / 'own_compounds.tsv'
    reg.parent.mkdir(parents=True)
    reg.write_text('own_compound_id\tstandardized_smiles\nPYZ-035\tCCO\n', encoding='utf-8')
    auditroot = tmp_path / 'audit'
    (auditroot / 'scripts').mkdir(parents=True)
    monkeypatch.setattr(closeout, '__file__', str(auditroot / 'scripts' / 'close_pyz_e0.py'))
    baseline = auditroot / 'outputs' / 'audits' / 'pyz_e0' / '20261114T085303Z-abcdef'
    baseline.mkdir(parents=True)
    (baseline / 'E0_ACCEPTANCE.md').write_text('# not accepted\n', encoding='utf-8')
    (baseline / 'environment.json').write_text(json.dumps({'root': str(root)}), encoding='utf-8')
    (baseline / 'file_checksums.sha256').write_text(
        closeout.sha256(reg) + '  config/own_compounds.tsv\n', encoding='utf-8')
    monkeypatch.setattr(closeout, 'request', lambda *a, **kw: (200, {'mcl_root': str(root)}, {}, 1.0))
    monkeypatch.setattr(closeout, 'reconcile_scientific', lambda *a, **kw: (
        [{'status': 'MATCH'}], [], [], []))
    before = closeout.sha256(reg)
    out, code = closeout.execute(Namespace(root=str(root), baseline=str(baseline),
        api_base='http://127.0.0.1:8100', ui_base='http://127.0.0.1:3100',
        ui_tools='', skip_ui=True, skip_cold_api=True, overwrite=False))
    assert code == 0
    assert out.is_relative_to(auditroot)
    assert (out / 'source_sha_entry.tsv').is_file()
    assert (out / 'source_sha_after.tsv').is_file()
    assert (out / 'E0_ACCEPTANCE_CLOSEOUT.md').is_file()
    assert 'NOT_APPROVED' in (out / 'E0_ACCEPTANCE_CLOSEOUT.md').read_text(encoding='utf-8')
    assert closeout.sha256(reg) == before
    assert not (root / 'outputs').exists()