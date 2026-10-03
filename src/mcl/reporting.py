from __future__ import annotations

import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from mcl.utils.hash import sha256_file, sha256_text


def _git_commit(root: Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"], stderr=subprocess.DEVNULL, text=True
        ).strip()
    except Exception:
        return None


def _configuration_hash(root: Path) -> str:
    files = [
        root / "config/project.yaml",
        root / "config/source_versions.yaml",
        root / "config/cancer_contexts.yaml",
        root / "config/field_ownership.yaml",
        root / "config/thresholds.yaml",
    ]
    payload = "\n".join(f"{p.name}:{sha256_file(p)}" for p in files if p.exists())
    return sha256_text(payload)


def create_manifest(
    root: Path,
    seed_path: Path,
    hgnc_path: Path,
    withdrawn_path: Path | None,
    identities_count: int,
    qc_errors: int,
    outputs: list[str],
) -> dict:
    return {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "code_version": "0.1.0",
        "git_commit": _git_commit(root),
        "python_version": platform.python_version(),
        "dependency_spec_sha256": sha256_file(root / "pyproject.toml"),
        "configuration_sha256": _configuration_hash(root),
        "input_target_seed": str(seed_path),
        "input_target_seed_sha256": sha256_file(seed_path),
        "hgnc_snapshot": str(hgnc_path),
        "hgnc_snapshot_sha256": sha256_file(hgnc_path),
        "withdrawn_snapshot": str(withdrawn_path) if withdrawn_path else None,
        "withdrawn_snapshot_sha256": sha256_file(withdrawn_path) if withdrawn_path else None,
        "targets_n": identities_count,
        "qc_errors_n": qc_errors,
        "qc_result": "PASS" if qc_errors == 0 else "FAIL",
        "outputs": outputs,
    }


def create_opentargets_manifest(
    root: Path,
    pairs_path: Path,
    target_identifiers_path: Path,
    source_release: str,
    endpoint: str,
    associations_count: int,
    qc_errors: int,
    qc_warnings: int,
    outputs: list[str],
) -> dict:
    return {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "code_version": "0.2.0",
        "milestone": "2A+2B Open Targets Wave 1",
        "git_commit": _git_commit(root),
        "python_version": platform.python_version(),
        "dependency_spec_sha256": sha256_file(root / "pyproject.toml"),
        "configuration_sha256": _configuration_hash(root),
        "input_pairs": str(pairs_path),
        "input_pairs_sha256": sha256_file(pairs_path),
        "target_identifiers": str(target_identifiers_path),
        "target_identifiers_sha256": sha256_file(target_identifiers_path),
        "opentargets_release": source_release,
        "opentargets_endpoint": endpoint,
        "cancer_target_pairs_n": associations_count,
        "qc_errors_n": qc_errors,
        "qc_warnings_n": qc_warnings,
        "qc_result": "PASS" if qc_errors == 0 else "FAIL",
        "outputs": outputs,
    }


def create_depmap_manifest(
    root: Path,
    release: str,
    files: dict,
    pairs_path: Path,
    evidence_n: int,
    context_audit_n: int,
    qc_errors: int,
    qc_warnings: int,
    outputs: list[str],
) -> dict:
    return {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "code_version": "0.3.5",
        "milestone": "3 DepMap functional dependency Wave 1",
        "git_commit": _git_commit(root),
        "python_version": platform.python_version(),
        "dependency_spec_sha256": sha256_file(root / "pyproject.toml"),
        "configuration_sha256": _configuration_hash(root),
        "depmap_release": release,
        "source_url": "https://depmap.org/portal/data_page/",
        "input_pairs": str(pairs_path),
        "input_pairs_sha256": sha256_file(pairs_path),
        "depmap_files": files,
        "evidence_rows_n": evidence_n,
        "context_audit_rows_n": context_audit_n,
        "qc_errors_n": qc_errors,
        "qc_warnings_n": qc_warnings,
        "qc_result": "PASS" if qc_errors == 0 else "FAIL",
        "outputs": outputs,
    }
