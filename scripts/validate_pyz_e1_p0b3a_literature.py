"""MCL E1-P0-B3A: verify a literature evidence overlay against original B3 ZIP.

Read-only scientific validation. NO experimental PYZ values or assay approvals.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

COMPARE_FIELDS = (
    'record_id', 'own_compound_id', 'target_gene', 'ligand_id',
    'bindingdb_reactant_set_id', 'bindingdb_monomer_id', 'endpoint',
    'raw_nm', 'article_doi',
)
REQUIRED_LEDGER = set(COMPARE_FIELDS) | {
    'publication_check', 'method_evidence', 'chemical_identity_evidence',
    'independence_review', 'review_status', 'primary_assay_approved',
    'pyz_activity_measured', 'source_urls', 'remaining_check',
}
VALID_REVIEW_STATUSES = {'PARTIAL_SOURCE_CHECK', 'PENDING_SOURCE_REVIEW'}


def checksum(p: Path) -> str:
    h = hashlib.sha256()
    with p.open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def records(stream):
    reader = csv.DictReader(stream, delimiter='\t')
    return reader, list(reader)


def validate(b3_zip: Path, ledger_path: Path, manifest_path: Path, expected_count: int = 11):
    for p in (b3_zip, ledger_path, manifest_path):
        if not p.is_file():
            raise FileNotFoundError(p)
    expected_hash = checksum(b3_zip)
    meta = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    if meta.get('contract') != 'mcl-pyz-e1-p0b3a-source-literature-check-v1':
        raise ValueError('Unknown evidence overlay contract')
    if meta.get('parent_zip_sha256', '').lower() != expected_hash:
        raise ValueError('Parent B3 source SHA-256 mismatch')
    with zipfile.ZipFile(b3_zip) as z:
        if z.testzip() is not None:
            raise ValueError('Damaged B3 source ZIP')
        if not {'assay_evidence_register.tsv', 'evidence_review_report.json'}.issubset(z.namelist()):
            raise ValueError('Required B3 source members missing')
        source_report = json.loads(z.read('evidence_review_report.json').decode('utf-8-sig'))
        if source_report.get('scientific_status') != 'NO_PRIMARY_ASSAYS_APPROVED':
            raise ValueError('B3 contains an unexpected scientific approval')
        source_reader, orig = records(io.StringIO(z.read('assay_evidence_register.tsv').decode('utf-8-sig')))
    with ledger_path.open(encoding='utf-8-sig', newline='') as stream:
        dest_reader, overlay = records(stream)
    if not REQUIRED_LEDGER.issubset(dest_reader.fieldnames or []):
        raise ValueError('Missing ledger columns')
    if len(orig) != expected_count or len(overlay) != expected_count:
        raise ValueError('Expected record count mismatch')
    by_orig = {r['record_id']: r for r in orig}
    by_out = {r['record_id']: r for r in overlay}
    if len(by_orig) != len(orig) or len(by_out) != len(overlay) or by_orig.keys() != by_out.keys():
        raise ValueError('Record_id cardinality or identity mismatch')
    counter = Counter()
    for rec_id, entry in by_out.items():
        src = by_orig[rec_id]
        for field in COMPARE_FIELDS:
            if entry[field] != src.get(field, ''):
                raise ValueError(f'{rec_id} immutable field changed: {field}')
        if entry['primary_assay_approved'] != 'NO' or entry['pyz_activity_measured'] != 'NO':
            raise ValueError('Forbidden assay/PYZ approval in literature ledger')
        if entry['review_status'] not in VALID_REVIEW_STATUSES:
            raise ValueError(f'Unknown review status: {entry["review_status"]}')
        if entry['review_status'] == 'PARTIAL_SOURCE_CHECK':
            links = [s for s in entry['source_urls'].split(';') if s]
            if not links or not all(x.startswith('https://') for x in links):
                raise ValueError('Partially checked record missing verifiable evidence links')
            if entry['method_evidence'] == 'NOT_REVIEWED':
                raise ValueError('Partially checked record lacks method evidence class')
        counter[entry['review_status']] += 1
    if meta.get('total_records') != len(overlay) or meta.get('source_checked_records') != counter['PARTIAL_SOURCE_CHECK']:
        raise ValueError('Manifest coverage mismatch')
    if meta.get('approved_primary_assays') != 0 or meta.get('external_admet_not_integrated') is not True:
        raise ValueError('Manifest scientific constraints violated')
    return {
        'contract': 'mcl-pyz-e1-p0b3a-verified-literature-overlay-v1',
        'verified_utc': datetime.now(timezone.utc).isoformat(),
        'input_sha256': {str(p): checksum(p) for p in (b3_zip, ledger_path, manifest_path)},
        'records': len(overlay), 'review_status_counts': dict(sorted(counter.items())),
        'scientific_status': 'NO_PRIMARY_ASSAYS_APPROVED', 'pyz_activity_measured': False,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--b3-zip', type=Path, required=True)
    parser.add_argument('--ledger', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--expected-count', type=int, default=11)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    out = args.out_dir.resolve()
    if out.exists():
        raise FileExistsError(out)
    if any(p.resolve().is_relative_to(out) for p in (args.b3_zip, args.ledger, args.manifest)):
        raise ValueError('Output directory must not contain the source data')
    report = validate(args.b3_zip, args.ledger, args.manifest, args.expected_count)
    out.mkdir(parents=True, exist_ok=False)
    (out / 'P0B3A_VALIDATION.json').write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
    print(f'Verified {report["records"]} records; zero approved primary assays.')
    return report

if __name__ == '__main__':
    main()