"""Synthetic tests for immutable B3 -> B3A evidence overlay (no real files required)."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import zipfile
from pathlib import Path

import pytest

from scripts.validate_pyz_e1_p0b3a_literature import validate, main

BASE = {
    'record_id': 'E1B3-ONE', 'own_compound_id': 'PYZ-035',
    'target_gene': 'CDK4', 'ligand_id': 'TLS-A',
    'bindingdb_reactant_set_id': '101', 'bindingdb_monomer_id': '500',
    'endpoint': 'IC50', 'raw_nm': '420', 'article_doi': '10.0000/test',
}


def make_fixture(tmp_path):
    zipped = tmp_path / 'b3.zip'
    source_io = io.StringIO()
    w = csv.DictWriter(source_io, list(BASE), delimiter='\t')
    w.writeheader()
    w.writerow(BASE)
    with zipfile.ZipFile(zipped, 'w') as z:
        z.writestr('assay_evidence_register.tsv', source_io.getvalue())
        z.writestr('evidence_review_report.json', json.dumps({'scientific_status': 'NO_PRIMARY_ASSAYS_APPROVED'}))
    row = {
        **BASE, 'publication_check': 'PUBLICATION_AND_BINDINGDB_MATCHED',
        'method_evidence': 'BINDINGDB_ASSAY_METADATA_ONLY',
        'chemical_identity_evidence': 'STEREO_AND_PAPER_SCHEME_PENDING',
        'independence_review': 'NOT_REVIEWED', 'review_status': 'PARTIAL_SOURCE_CHECK',
        'primary_assay_approved': 'NO', 'pyz_activity_measured': 'NO',
        'source_urls': 'https://bindingdb.org/test', 'remaining_check': 'review primary article',
    }
    ledger = tmp_path / 'evidence.tsv'
    def save(values):
        with ledger.open('w', newline='', encoding='utf-8-sig') as stream:
            writer = csv.DictWriter(stream, row, delimiter='\t')
            writer.writeheader()
            writer.writerow(values)
    save(row)
    manifest = tmp_path / 'manifest.json'
    metadata = {
        'contract': 'mcl-pyz-e1-p0b3a-source-literature-check-v1',
        'parent_zip_sha256': hashlib.sha256(zipped.read_bytes()).hexdigest(),
        'total_records': 1, 'source_checked_records': 1,
        'approved_primary_assays': 0, 'external_admet_not_integrated': True,
    }
    manifest.write_text(json.dumps(metadata), encoding='utf-8')
    return zipped, ledger, manifest, row, save, metadata


def test_source_overlay_validation(tmp_path):
    a, b, c, *_ = make_fixture(tmp_path)
    data = validate(a, b, c, expected_count=1)
    assert data['records'] == 1
    assert data['review_status_counts'] == {'PARTIAL_SOURCE_CHECK': 1}
    assert data['scientific_status'] == 'NO_PRIMARY_ASSAYS_APPROVED'
    out = tmp_path / 'out'
    main(['--b3-zip', str(a), '--ledger', str(b), '--manifest', str(c), '--expected-count', '1', '--out-dir', str(out)])
    assert (out / 'P0B3A_VALIDATION.json').is_file()


def test_reject_false_approval(tmp_path):
    a, b, c, row, save, _ = make_fixture(tmp_path)
    save({**row, 'primary_assay_approved': 'YES'})
    with pytest.raises(ValueError, match='Forbidden assay'):
        validate(a, b, c, expected_count=1)


def test_reject_mutated_endpoints(tmp_path):
    a, b, c, row, save, _ = make_fixture(tmp_path)
    save({**row, 'raw_nm': '100'})
    with pytest.raises(ValueError, match='immutable field changed'):
        validate(a, b, c, expected_count=1)


def test_reject_changed_source_hash(tmp_path):
    a, b, c, row, save, meta = make_fixture(tmp_path)
    c.write_text(json.dumps({**meta, 'parent_zip_sha256': '0' * 64}), encoding='utf-8')
    with pytest.raises(ValueError, match='SHA-256'):
        validate(a, b, c, expected_count=1)


def test_reject_no_public_source(tmp_path):
    a, b, c, row, save, meta = make_fixture(tmp_path)
    save({**row, 'source_urls': ''})
    with pytest.raises(ValueError, match='missing verifiable'):
        validate(a, b, c, expected_count=1)


def test_reject_existing_output(tmp_path):
    a, b, c, *_ = make_fixture(tmp_path)
    out = tmp_path / 'previous_result'
    out.mkdir()
    with pytest.raises(FileExistsError):
        main(['--b3-zip', str(a), '--ledger', str(b), '--manifest', str(c), '--expected-count', '1', '--out-dir', str(out)])