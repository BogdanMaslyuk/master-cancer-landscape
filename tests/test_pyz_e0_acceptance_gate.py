"""Offline regression tests: synthetic E0 evidence only, no science data edits."""
from __future__ import annotations

import csv
import hashlib
import json
import sys
import zipfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_pyz_e0_acceptance_gate import Evidence, analyze, main


def tsv(columns, rows):
    from io import StringIO
    v = StringIO()
    w = csv.writer(v, delimiter='\t', lineterminator='\n')
    w.writerow(columns)
    w.writerows(rows)
    return v.getvalue().encode()


def fixture(tmp_path, closeout=True, mismatch=False, bad_matrix=False, corrupt_png=False):
    root = tmp_path / 'e0'
    root.mkdir()
    source_sha = hashlib.sha256(b'raw scientific data').hexdigest()
    (root / 'E0_ACCEPTANCE.md').write_text('НЕ ПРИНЯТО / AC-09\n', encoding='utf-8')
    (root / 'file_checksums.sha256').write_text(f'{source_sha}  config/own_compounds.tsv\n')
    (root / 'input_inventory.tsv').write_bytes(tsv(['status'], [['NOT_REVIEWED'], ['present']]))
    (root / 'api_cases.tsv').write_bytes(tsv(['case','pass'], [['PYZ-035', 'true']]))
    (root / 'performance.tsv').write_bytes(tsv(['route','p95_ms'], [['/api/pyz/summary','120']]))
    (root / 'schema_audit.json').write_text(json.dumps({'issues': []}))
    coverage = [[f'PYZ-{i:03d}',f'GENE-{j:02d}','present'] for i in range(1,66) for j in range(1,14)]
    if bad_matrix:
        coverage.pop()
    (root / 'coverage_matrix.tsv').write_bytes(tsv(['compound_id','target_gene','state'], coverage))
    (root / 'ui_checklist.tsv').write_bytes(tsv(['status'], [['NOT_REVIEWED']]*9))
    if closeout:
        out = root / 'closeout'
        (out/'ui/screenshots').mkdir(parents=True)
        for tag in range(9):
            (out / f'ui/screenshots/screen-{tag}.png').write_bytes(b'x' if corrupt_png else b'\x89PNG\r\n\x1a\nimage')
        for phase in ['entry','after']:
            (out / f'source_sha_{phase}.tsv').write_bytes(tsv(
                ['path','baseline_sha256','current_sha256','status'],
                [['config/own_compounds.tsv',source_sha,'a'*64 if mismatch else source_sha,'CHANGED' if mismatch else 'UNCHANGED']]
            ))
        (out / 'source_api_reconciliation.tsv').write_bytes(tsv(['status'],[['MATCH']]))
        (out / 'cold_api.json').write_text(json.dumps({'status':'PASS'}))
        (out / 'browser_status.json').write_text(json.dumps({'n_screenshots':9,'status':'PASS'}))
        (out / 'ui/browser_results.json').write_text(json.dumps({'screenshots':[{'status':'AUTO_PASS'}]*9}))
        (out / 'E0_ACCEPTANCE_CLOSEOUT.md').write_text('AC-09 NOT APPROVED\n')
    return root


def as_zip(root, path):
    with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for f in root.rglob('*'):
            if f.is_file():
                z.write(f, f.relative_to(root).as_posix())
    return path


def test_realistic_synthetic_closeout_pending_owner(tmp_path):
    root = fixture(tmp_path)
    z = as_zip(root,tmp_path/'E0.zip')
    source = Evidence(z)
    try:
        r, rows, findings, hashes = analyze(source)
    finally:
        source.close()
    assert r['criteria_count'] == 9
    assert r['final_status'] == 'PENDING_OWNER_ACCEPTANCE'
    assert r['merge_permitted'] is False
    assert r['scientific_status'] == 'NOT_ACCEPTED_MANUAL_REVIEW_REQUIRED'
    assert r['screenshot_files'] == 9
    assert len(hashes) == 2
    assert next(x for x in rows if x['criterion']=='AC-09')['manual_decision']=='PENDING'


def test_baseline_without_closeout_is_incomplete(tmp_path):
    root = fixture(tmp_path,closeout=False)
    r,*_ = analyze(Evidence(root))
    assert r['final_status'] == 'INCOMPLETE_EVIDENCE'
    assert r['merge_permitted'] is False
    assert r['closeout_missing_files']


def test_hash_mismatch_blocks(tmp_path):
    r,*_ = analyze(Evidence(fixture(tmp_path,mismatch=True)))
    assert r['final_status']=='BLOCKED' and r['blocking_findings'] > 0


def test_matrix_coverage_mismatch_blocks(tmp_path):
    r,*_ = analyze(Evidence(fixture(tmp_path,bad_matrix=True)))
    assert r['final_status']=='BLOCKED'


def test_api_mismatch_blocks(tmp_path):
    root = fixture(tmp_path)
    (root/'closeout/source_api_reconciliation.tsv').write_bytes(tsv(['status'],[['MISMATCH']]))
    r,*_ = analyze(Evidence(root))
    assert r['final_status']=='BLOCKED'


def test_invalid_png_blocks(tmp_path):
    r,*_ = analyze(Evidence(fixture(tmp_path,corrupt_png=True)))
    assert r['final_status']=='BLOCKED'


def test_zip_traversal_rejected(tmp_path):
    z=tmp_path/'evil.zip'
    with zipfile.ZipFile(z,'w') as output:
        output.writestr('E0_ACCEPTANCE.md','x')
        output.writestr('../danger.txt','bad')
    with pytest.raises(ValueError,match='Unsafe'):
        Evidence(z)


def test_zip_duplicate_entry_rejected(tmp_path):
    z=tmp_path/'duplicate.zip'
    with zipfile.ZipFile(z,'w') as output:
        output.writestr('E0_ACCEPTANCE.md','x')
        output.writestr('E0_ACCEPTANCE.md','y')
    with pytest.raises(ValueError,match='Duplicate'):
        Evidence(z)


def test_never_overwrite_or_edit_input(tmp_path):
    root=fixture(tmp_path)
    z=as_zip(root,tmp_path/'E0.zip')
    original=hashlib.sha256(z.read_bytes()).hexdigest()
    out=tmp_path/'report'
    result=main(['--e0',str(z),'--out-dir',str(out)])
    assert result['merge_permitted'] is False
    assert hashlib.sha256(z.read_bytes()).hexdigest()==original
    assert (out/'acceptance_matrix.tsv').is_file()
    with pytest.raises(FileExistsError):
        main(['--e0',str(z),'--out-dir',str(out)])


def test_missing_baseline_or_forged_approval_never_accepted(tmp_path):
    root=fixture(tmp_path)
    (root/'E0_ACCEPTANCE.md').write_text('ПРИНЯТО, ВСЁ ОТЛИЧНО')
    r,*_=analyze(Evidence(root))
    assert r['merge_permitted'] is False
    (root/'E0_ACCEPTANCE.md').unlink()
    with pytest.raises(ValueError,match='exactly one'):
        Evidence(root)


def test_live_original_sources_verified_and_mismatch_blocks(tmp_path):
    root=fixture(tmp_path)
    science=tmp_path/'science'
    (science/'config').mkdir(parents=True)
    raw=science/'config/own_compounds.tsv'
    raw.write_bytes(b'raw scientific data')
    r,*_=analyze(Evidence(root),science)
    assert r['fresh_live_sha_checked'] is True
    assert r['fresh_live_sha_rows']==1
    assert r['final_status']=='PENDING_OWNER_ACCEPTANCE'
    raw.write_bytes(b'tampered')
    r,*_=analyze(Evidence(root),science)
    assert r['final_status']=='BLOCKED'
    assert r['fresh_live_sha_checked'] is False


def test_e0_review_bundle_preserves_original_png_and_exposes_warning_reason(tmp_path):
    root = fixture(tmp_path)
    original_png = root / 'closeout/ui/screenshots/screen-0.png'
    original_hash = hashlib.sha256(original_png.read_bytes()).hexdigest()
    browser_results = {
        'screenshots': [
            {'route': '/pyz', 'width': 1440, 'height': 900, 'saved': True,
             'http_status': 200, 'heading_matches': True,
             'horizontal_overflow': False, 'keyboard_focus': False,
             'console_errors': [], 'status': 'AUTO_WARN', 'screenshot': 'screenshots/screen-0.png'}
        ]
    }
    (root / 'closeout/ui/browser_results.json').write_text(json.dumps(browser_results), encoding='utf-8')
    out = tmp_path / 'final_report'
    report = main(['--e0', str(root), '--out-dir', str(out)])
    assert report['merge_permitted'] is False
    copy = out / 'review_evidence/closeout/ui/screenshots/screen-0.png'
    assert copy.read_bytes() == original_png.read_bytes()
    assert hashlib.sha256(copy.read_bytes()).hexdigest() == original_hash
    manifest = json.loads((out / 'review_evidence_manifest.json').read_text(encoding='utf-8'))
    assert manifest['screenshots'] == 9
    entry = next(x for x in manifest['files'] if x['name'].endswith('screen-0.png'))
    assert entry['sha256'] == original_hash
    assert (out / 'review_evidence/closeout/ui/browser_results.json').is_file()
    with (out / 'ui_diagnostics.tsv').open(encoding='utf-8-sig') as handle:
        items = list(csv.DictReader(handle, delimiter='\t'))
    assert len(items) == 1
    assert items[0]['warning_causes'] == 'KEYBOARD_FOCUS_MISSING'
    assert items[0]['route'] == '/pyz'
    assert hashlib.sha256(original_png.read_bytes()).hexdigest() == original_hash


def test_e0_review_bundle_reads_zip_without_unpacking_source(tmp_path):
    root = fixture(tmp_path)
    archive = as_zip(root, tmp_path / 'source_e0.zip')
    before = hashlib.sha256(archive.read_bytes()).hexdigest()
    out = tmp_path / 'evidence_review'
    main(['--e0', str(archive), '--out-dir', str(out)])
    assert (out / 'review_evidence/closeout/ui/screenshots/screen-8.png').is_file()
    assert (out / 'review_evidence/closeout/performance.tsv').exists() is False
    assert (out / 'review_evidence/performance.tsv').is_file()
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == before



def test_viewport_only_screenshots_are_not_fake_browser_failures(tmp_path):
    root = fixture(tmp_path)
    source = root / 'closeout/ui/browser_results.json'
    source.write_text(json.dumps({
        'status': 'PARTIAL_SCREENSHOTS_MANUAL_REVIEW_PENDING',
        'method': 'headless Edge/Chrome CLI, viewport-only',
        'screenshots': [
            {'route': '/pyz', 'width': 390, 'height': 844,
             'status': 'CAPTURED', 'saved': True,
             'path': 'ui\\screenshots\\screen-0.png',
             'note': 'Only screenshot: DOM checks not executed'}
            for _ in range(9)
        ],
        'manual_review_required': True,
    }), encoding='utf-8')
    source_zip = as_zip(root, tmp_path / 'source.zip')
    result_dir = tmp_path / 'review'
    report = main(['--e0', str(source_zip), '--out-dir', str(result_dir)])
    assert report['browser_capture_only_count'] == 9
    assert report['browser_executed_check_count'] == 0
    assert report['merge_permitted'] is False
    assert report['final_status'] == 'PENDING_OWNER_ACCEPTANCE'
    findings = list(csv.DictReader(
        (result_dir / 'findings.tsv').open(encoding='utf-8-sig'), delimiter='\t'))
    ac5 = [row['detail'] for row in findings if row['criterion'] == 'AC-05']
    assert any('without HTTP/DOM/keyboard checks: 9' in s for s in ac5)
    assert all('checks with warnings: 9' not in s for s in ac5)
    with (result_dir / 'ui_diagnostics.tsv').open(encoding='utf-8-sig') as file:
        screenshots = list(csv.DictReader(file, delimiter='\t'))
    assert len(screenshots) == 9
    assert all(row['warning_causes'] == 'NOT_TESTED_VIEWPORT_ONLY' for row in screenshots)
    assert all(row['http_status'] == '' for row in screenshots)
    assert all(row['screenshot'].endswith('screen-0.png') for row in screenshots)


def test_executed_browser_failure_is_not_downgraded_to_capture_only(tmp_path):
    root = fixture(tmp_path)
    source = root / 'closeout/ui/browser_results.json'
    source.write_text(json.dumps({
        'method': 'Playwright',
        'screenshots': [{
            'route': '/pyz', 'width': 390, 'height': 844,
            'status': 'AUTO_WARN', 'saved': True, 'http_status': 200,
            'heading_matches': True, 'horizontal_overflow': True,
            'keyboard_focus': True, 'console_errors': [],
            'screenshot': 'screenshots/screen-0.png'}]
    }), encoding='utf-8')
    out = tmp_path / 'review'
    report = main(['--e0', str(root), '--out-dir', str(out)])
    assert report['browser_capture_only_count'] == 0
    assert report['browser_executed_check_count'] == 1
    with (out / 'ui_diagnostics.tsv').open(encoding='utf-8-sig') as file:
        rows = list(csv.DictReader(file, delimiter='\t'))
    assert rows[0]['warning_causes'] == 'HORIZONTAL_OVERFLOW'
