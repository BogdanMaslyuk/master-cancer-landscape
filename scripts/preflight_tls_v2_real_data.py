#!/usr/bin/env python3
"""P0.3.6 read-only preflight: verify *real* scientific sources before a pilot.

Writes a diagnostic JSON report in an isolated output directory. Does NOT
download scientific data, change source files, approve an assay or publish v2.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for part in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(part)
    return digest.hexdigest()


def _check_file(path: Path | None, label: str, checks: dict) -> bool:
    present = bool(path and path.is_file())
    checks[label + "_available"] = present
    return present


def _report_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def inspect(science_root: Path, baseline_file: Path, registry: Path,
            registry_manifest: Path, *, source_manifest: Path | None = None,
            review_file: Path | None = None) -> dict:
    """The report is diagnostic; only all checks True means ready for strict run."""
    science_root = science_root.resolve()
    checks: dict[str, bool] = {}
    issues: list[str] = []
    evidence: dict = {}
    baseline = None
    if _check_file(baseline_file, "p01_baseline", checks):
        try:
            baseline = _report_json(baseline_file)
            if baseline.get("contract") != "mcl-expanded-target-space-baseline-audit-v1":
                issues.append("p01_unrecognized_baseline_contract")
                checks["p01_contract"] = False
            else:
                checks["p01_contract"] = True
            checks["p01_core13"] = (
                baseline.get("core_targets", {}).get("target_genes_n") == 13
                and baseline.get("core_targets", {}).get("matches_expected_reference_count") is True
            )
            checks["p01_pyz65"] = (
                baseline.get("own_compounds", {}).get("compound_ids_n") == 65
                and baseline.get("own_compounds", {}).get("matches_expected_reference_count") is True
            )
            checks["p01_source_root"] = (
                Path(baseline.get("science_root", "")).resolve() == science_root
            )
            checks["p01_input_files_unchanged"] = True
            for logical_name, item in baseline.get("files", {}).items():
                path = science_root / str(item.get("relative_path", ""))
                expected = item.get("sha256")
                if not item.get("present") or not expected or not path.is_file():
                    checks["p01_input_files_unchanged"] = False
                    issues.append("p01_unavailable_source:" + logical_name)
                    continue
                if sha256(path) != expected:
                    checks["p01_input_files_unchanged"] = False
                    issues.append("p01_source_changed:" + logical_name)
            evidence["p01_core_genes"] = baseline.get("core_targets", {}).get("target_genes", [])
            evidence["p01_original_source_hashes"] = {
                k: v.get("sha256") for k,v in baseline.get("files", {}).items()
            }
        except (OSError, ValueError, KeyError, TypeError) as exc:
            issues.append("p01_invalid_baseline:" + type(exc).__name__)
            checks["p01_readable"] = False
    else:
        issues.append("p01_baseline_absent")
    snapshot = None
    if _check_file(registry, "p02_registry", checks) and _check_file(
        registry_manifest, "p02_manifest", checks
    ):
        try:
            snapshot = _report_json(registry_manifest)
            checks["p02_contract"] = snapshot.get("contract") == "mcl-protein-target-universe-v2-p0"
            frame = pd.read_csv(registry, sep="\t", dtype=str, keep_default_na=False)
            needed = {"target_gene", "uniprot_accession", "mapping_status", "organism_id"}
            checks["p02_schema"] = needed.issubset(frame.columns)
            checks["p02_unique_mapping"] = bool(
                checks["p02_schema"] and len(frame) > 0
                and frame["target_gene"].is_unique
                and frame["uniprot_accession"].is_unique
                and (frame["mapping_status"] == "unique_reviewed_human").all()
                and (frame["organism_id"] == "9606").all()
                and int(snapshot.get("target_genes_n", -1)) == len(frame)
            )
            checks["p02_pinned_uniprot_snapshot"] = (
                snapshot.get("uniprot_snapshot_included") is True
                and bool(snapshot.get("uniprot_release"))
                and "REPLACE" not in str(snapshot.get("uniprot_release")).upper()
                and "ACTUAL_" not in str(snapshot.get("uniprot_release")).upper()
            )
            uniprot_url = urlparse(str(snapshot.get("uniprot_source_url") or ""))
            checks["p02_official_source_declared"] = (
                uniprot_url.scheme == "https"
                and (uniprot_url.hostname or "").endswith("uniprot.org")
            )
            checks["p02_source_hashes"] = True
            for source in snapshot.get("source_files", []):
                path = Path(source.get("path", ""))
                if not path.is_file() or sha256(path) != source.get("sha256"):
                    checks["p02_source_hashes"] = False
                    issues.append("p02_source_changed_or_missing:" + str(path))
            if not snapshot.get("source_files"):
                checks["p02_source_hashes"] = False
            if checks["p02_schema"]:
                evidence["p02_proteins_n"] = len(frame)
                evidence["p02_core_genes_mapped_n"] = len(
                    set(evidence.get("p01_core_genes", [])) & set(frame["target_gene"])
                )
                checks["p02_core13_mapped"] = evidence["p02_core_genes_mapped_n"] == 13
            evidence["p02_registry_sha256"] = sha256(registry)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            checks["p02_readable"] = False
            issues.append("p02_invalid_registry:" + type(exc).__name__)
    else:
        issues.append("p02_registry_or_manifest_absent")
    if source_manifest is not None:
        if _check_file(source_manifest, "p03_source_manifest", checks):
            try:
                record = _report_json(source_manifest)
                checks["p03_manifest_contract"] = (
                    record.get("contract") == "mcl-tls-v2-pinned-extracts-v1"
                )
                rows = record.get("sources", [])
                checks["p03_exact_providers"] = (
                    len(rows) == 2 and {x.get("source") for x in rows} == {"ChEMBL","BindingDB"}
                )
                checks["p03_source_files_hash"] = bool(checks["p03_exact_providers"])
                checks["p03_source_provenance"] = True
                for source in rows:
                    path = Path(source.get("path", ""))
                    if not path.is_absolute():
                        path = source_manifest.parent / path
                    if not path.is_file() or sha256(path) != source.get("sha256"):
                        checks["p03_source_files_hash"] = False
                        issues.append("p03_source_changed_or_missing:" + str(path))
                    if not all(str(source.get(field, "")).strip()
                               for field in ("release","url","license")):
                        checks["p03_source_provenance"] = False
                    parsed = urlparse(str(source.get("url") or ""))
                    hostname = parsed.hostname or ""
                    official = (
                        (source.get("source") == "ChEMBL"
                         and (hostname == "ebi.ac.uk" or hostname.endswith(".ebi.ac.uk")))
                        or (source.get("source") == "BindingDB"
                            and (hostname == "bindingdb.org" or hostname.endswith(".bindingdb.org")))
                    )
                    if parsed.scheme != "https" or not official:
                        checks["p03_source_provenance"] = False
                    if any(
                        placeholder in str(source.get("license", "")).upper()
                        for placeholder in ("REPLACE","UNKNOWN","TODO","VERIFY")
                    ):
                        checks["p03_source_provenance"] = False
                evidence["p03_source_releases"] = {
                    s.get("source"):s.get("release") for s in rows
                }
            except (OSError, ValueError, KeyError, TypeError) as exc:
                checks["p03_manifest_readable"] = False
                issues.append("p03_invalid_manifest:" + type(exc).__name__)
    else:
        checks["p03_source_manifest_available"] = False
        issues.append("p03_sources_not_yet_pinned")
    if review_file is not None and _check_file(review_file, "assay_review", checks):
        try:
            table = pd.read_csv(review_file, sep="\t", dtype=str, keep_default_na=False)
            required = {"source","source_release","assay_id","method_class",
                        "reviewer","evidence_url","verification_status"}
            checks["assay_review_schema"] = required.issubset(table.columns)
            approved = table.loc[
                table["verification_status"].eq("approved")
            ] if checks["assay_review_schema"] else pd.DataFrame()
            checks["assay_review_has_approved_records"] = not approved.empty
            checks["assay_review_metadata_complete"] = (
                not approved.empty
                and approved[["source","source_release","assay_id","method_class",
                              "reviewer","evidence_url"]].ne("").all().all()
            )
            evidence["assay_review_approved_rows_n"] = len(approved)
            issues.append("assay_review_is_human_assertion_not_independent_source_verification")
        except (OSError, ValueError, KeyError, TypeError) as exc:
            checks["assay_review_readable"] = False
            issues.append("assay_review_invalid:" + type(exc).__name__)
    else:
        checks["assay_review_available"] = False
        issues.append("assay_review_not_confirmed")
    for check, passed in checks.items():
        if passed is not True:
            issues.append("failed_check:" + check)
    return {
        "contract": "mcl-tls-v2-p036-real-data-preflight-v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "data_root": str(science_root),
        "checks": checks,
        "ready_for_full_real_data_pilot": bool(checks) and all(checks.values()),
        "evidence": evidence,
        "blocking_or_informational_issues": sorted(set(issues)),
        "limitations": [
            "A successful preflight does not confirm binding methods or biological activity.",
            "No external archives are downloaded or fabricated.",
            "Source URLs are declarations; this check does not establish archive authenticity without official checksums.",
            "Scientific A/B signoff and P0.3.6 final acceptance require independent human review.",
        ],
    }


def _version(package: str) -> str | None:
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return None


def run_context(output_dir: Path) -> dict:
    parent = output_dir.resolve()
    while not parent.exists():
        parent = parent.parent
    usage = shutil.disk_usage(parent)
    git = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT,
        capture_output=True, text=True, check=False, timeout=20
    )
    return {
        "contract": "mcl-tls-v2-p036-run-context-v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_sha": git.stdout.strip() if git.returncode == 0 else None,
        "python": sys.version.split()[0], "os": platform.platform(),
        "versions": {name: _version(name) for name in
                     ("rdkit","pandas","pyarrow","numpy","master-cancer-landscape")},
        "output_volume_free_bytes_at_preflight": usage.free,
        "output_volume_total_bytes": usage.total,
    }


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--science-root", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--registry-manifest", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path)
    parser.add_argument("--assay-review", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()
    if args.out_dir.resolve().is_relative_to(args.science_root.resolve()):
        raise SystemExit("Preflight reports must be outside the original read-only scientific root.")
    result = inspect(
        args.science_root, args.baseline, args.registry, args.registry_manifest,
        source_manifest=args.source_manifest, review_file=args.assay_review
    )
    write_json(args.out_dir/"preflight_report.json", result)
    write_json(args.out_dir/"run_context.json", run_context(args.out_dir))
    print(f"Preflight: {'READY' if result['ready_for_full_real_data_pilot'] else 'BLOCKED'}")
    for label, ok in result["checks"].items():
        print(f"  {'OK' if ok else 'BLOCKED'}: {label}")
    if args.strict and not result["ready_for_full_real_data_pilot"]:
        raise SystemExit("Real-data P0.3.6 preflight blocked; no scientific data changed.")


if __name__ == "__main__":
    main()
