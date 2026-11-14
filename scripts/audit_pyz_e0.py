"""E0: read-only, evidence-producing audit of the local MCL/PYZ review worktree.

The only writes are to outputs/audits/pyz_e0/<UTC>-<git-sha>[-N].
No raw/processed/runtime data, API, or frontend code is modified.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import re
import statistics
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

AUDITOR_VERSION = "e0-audit-1.0"
PROFILE_VERSION = "mcl-pyz-profile-v2-p0"
EXPECTED_IDS = [f"PYZ-{i:03d}" for i in range(1, 66)]
SOURCE_FILES = {
    "registry": ("config/own_compounds.tsv", True),
    "admet_individual": ("config/pyz_admet_snapshot.tsv", False),
    "admet_series": ("config/pyz_admet_summary.json", False),
    "matrix": ("data/runtime/pyz_target_evidence/pyz_target_evidence_matrix.parquet", True),
    "neighbors": ("data/runtime/own_compounds/target_ligand_space_similarity_top_hits.parquet", True),
    "candidate": ("data/runtime/pharmacology/candidate_hypotheses_v2.parquet", False),
    "candidate_models": ("data/runtime/pharmacology/candidate_hypothesis_models_v2.parquet", False),
    "laboratory": ("data/processed/laboratory_panel.parquet", False),
    "genes_manifest": ("data/runtime/explorer/manifest.json", False),
    "genes_catalog": ("data/runtime/explorer/gene_catalog.parquet", False),
    "genes_context": ("data/runtime/explorer/gene_context_metrics.parquet", False),
    "genes_annotations": ("data/runtime/explorer/gene_annotations.parquet", False),
    "genes_dependency": ("data/runtime/explorer/gene_dependency_summary.parquet", False),
    "genes_dependency_json": ("data/runtime/explorer/gene_dependency_summary.json", False),
    "gene_effect_npy": ("data/runtime/explorer/model_layers/gene_effect.npy", False),
    "gene_effect_meta": ("data/runtime/explorer/model_layers/gene_effect.json", False),
    "expression_npy": ("data/runtime/explorer/model_layers/expression.npy", False),
    "expression_meta": ("data/runtime/explorer/model_layers/expression.json", False),
    "copy_number_npy": ("data/runtime/explorer/model_layers/copy_number.npy", False),
    "copy_number_meta": ("data/runtime/explorer/model_layers/copy_number.json", False),
}
METADATA_ONLY_PARQUET = {"genes_catalog", "genes_context", "genes_annotations", "genes_dependency"}
REQUIRED_COLUMNS = {
    "registry": ("own_compound_id", "standardized_smiles"),
    "matrix": ("own_compound_id", "target_gene", "best_tanimoto_morgan_r2_2048"),
    "neighbors": ("own_compound_id", "target_gene", "rank_in_target"),
    "candidate": ("hypothesis_id", "target_gene"),
    "candidate_models": ("hypothesis_id",),
}
# Exit states: 0 = audit outputs generated; they are NOT an acceptance decision;
# 2 = unable to audit critical inputs or wrong MCL_ROOT; 3 = operational failure.


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def run_cmd(args: list[str], cwd: Path, timeout: int = 90) -> tuple[int, str]:
    try:
        p = subprocess.run(args, cwd=cwd, text=True, encoding="utf-8", errors="replace",
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout, check=False)
        return p.returncode, p.stdout.strip()
    except (OSError, subprocess.TimeoutExpired) as exc:
        return 127, f"{type(exc).__name__}: {exc}"


def git(root: Path, *args: str) -> str:
    code, value = run_cmd(["git", "-C", str(root), *args], root)
    return value if code == 0 else f"GIT_UNAVAILABLE: {value}"


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as stream:
        stream.write(value)
    os.replace(tmp, path)


def save_json(path: Path, obj: Any) -> None:
    atomic_text(path, json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False, default=str) + "\n")


def save_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    from io import StringIO
    buffer = StringIO(newline="")
    writer = csv.DictWriter(buffer, delimiter="\t", fieldnames=fields, extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: json.dumps(row[key], ensure_ascii=False, default=str) if isinstance(row.get(key), (list, dict)) else ("" if row.get(key) is None else row[key]) for key in fields})
    atomic_text(path, buffer.getvalue())


def defect(issues: list[dict], category: str, severity: str, code: str, message: str,
           source: str = "", evidence: str = "", repair_stage: str = "E1") -> None:
    issues.append({"id": f"E0-{len(issues) + 1:03d}", "category": category,
                   "severity": severity, "code": code, "message": message,
                   "source": source, "evidence": evidence, "repair_stage": repair_stage})


def read_frame(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    if path.suffix == ".tsv":
        return pd.read_csv(path, sep="\t", dtype={"own_compound_id": "string"})
    raise ValueError("Not a tabular source")


def source_path(root: Path, path: str) -> Path:
    primary = root / path
    if primary.exists():
        return primary
    if primary.suffix == ".parquet":
        alt = primary.with_suffix(".tsv")
        if alt.exists():
            return alt
    return primary


def inventory_sources(root: Path, issues: list[dict]) -> tuple[list[dict], dict[str, pd.DataFrame], dict[str, str]]:
    inventory: list[dict] = []
    frames: dict[str, pd.DataFrame] = {}
    protected_hashes: dict[str, str] = {}
    for label, (relative, required) in SOURCE_FILES.items():
        path = source_path(root, relative)
        actual = path.relative_to(root).as_posix()
        row: dict[str, Any] = {"label": label, "path": actual, "required": required,
              "exists": path.is_file(), "format": path.suffix.lstrip("."), "bytes": "",
              "modified_utc": "", "sha256": "", "rows": "", "columns_n": "",
              "column_types": "", "nulls_key": "", "duplicate_keys": "", "status": "NOT_AVAILABLE", "error": ""}
        if not path.is_file():
            if required:
                defect(issues, "DATA", "blocker", "MISSING_REQUIRED", f"Отсутствует обязательный файл {actual}", actual, repair_stage="ENV")
            inventory.append(row)
            continue
        try:
            stat = path.stat()
            row.update(bytes=stat.st_size,
                       modified_utc=datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
                       sha256=sha256(path))
            protected_hashes[actual] = row["sha256"]
            if path.suffix in (".tsv", ".parquet"):
                if label in METADATA_ONLY_PARQUET and path.suffix == ".parquet":
                    import pyarrow.parquet as pq
                    metadata = pq.ParquetFile(path)
                    row["rows"] = metadata.metadata.num_rows
                    row["columns_n"] = len(metadata.schema_arrow)
                    row["column_types"] = {field.name: str(field.type) for field in metadata.schema_arrow}
                    row["status"] = "NOT_REVIEWED"  # Full gene-index validation is delegated below.
                else:
                    frame = read_frame(path)
                    frames[label] = frame
                    row["rows"] = len(frame)
                    row["columns_n"] = len(frame.columns)
                    row["column_types"] = {name: str(dtype) for name, dtype in frame.dtypes.items()}
                    keys = [x for x in REQUIRED_COLUMNS.get(label, ()) if x in frame.columns]
                    row["nulls_key"] = {x: int(frame[x].isna().sum()) for x in keys}
                    key_grains = {
                        "registry": ["own_compound_id"],
                        "admet_individual": ["own_compound_id"],
                        "matrix": ["own_compound_id", "target_gene"],
                        "neighbors": ["own_compound_id", "target_gene", "rank_in_target"],
                    }
                    grain = key_grains.get(label, [])
                    if grain and set(grain).issubset(frame.columns):
                        row["duplicate_keys"] = int(frame.duplicated(subset=grain).sum())
                    row["status"] = "PASS"
            elif path.suffix == ".json":
                json.loads(path.read_text(encoding="utf-8-sig"))
                row["status"] = "PASS"
            elif path.suffix == ".npy":
                # Do not read potentially multi-GB arrays; hash + companion schema suffice at E0.
                row["status"] = "NOT_REVIEWED"
            else:
                row["status"] = "NOT_REVIEWED"
        except (Exception) as exc:
            row["status"] = "FAIL"
            row["error"] = f"{type(exc).__name__}: {exc}"
            defect(issues, "DATA", "blocker" if required else "major", "UNREADABLE_INPUT",
                   f"Невозможно прочитать {actual}: {row['error']}", actual, repair_stage="ENV")
        inventory.append(row)
    return inventory, frames, protected_hashes


def require_columns(frame: pd.DataFrame, required: tuple[str, ...], label: str,
                    issues: list[dict]) -> bool:
    absent = sorted(set(required) - set(frame.columns))
    if absent:
        defect(issues, "DATA", "blocker" if label in ("registry", "matrix") else "major",
               "MISSING_COLUMN", f"{label}: нет обязательных колонок: {', '.join(absent)}", label, repair_stage="ENV")
        return False
    return True


def valid_similarity(series: pd.Series) -> pd.Series:
    numbers = pd.to_numeric(series, errors="coerce")
    return numbers.notna() & numbers.map(math.isfinite) & numbers.between(0, 1)


def inspect_data(frames: dict[str, pd.DataFrame], issues: list[dict]) -> tuple[dict, list[dict]]:
    result: dict[str, Any] = {"expected_registry_n": 65, "expected_target_n": 13,
                              "expected_pairs_if_complete": 845, "auditor_version": AUDITOR_VERSION}
    coverage: list[dict] = []
    for label, required in REQUIRED_COLUMNS.items():
        if label in frames:
            require_columns(frames[label], required, label, issues)
    registry = frames.get("registry")
    matrix = frames.get("matrix")
    if registry is None or not {"own_compound_id", "standardized_smiles"}.issubset(registry.columns):
        result["registry_status"] = "FAIL"
        return result, coverage

    ids = registry["own_compound_id"].fillna("").astype(str).str.strip().tolist()
    duplicates = registry["own_compound_id"].duplicated(keep=False).sum()
    unknown = sorted(set(ids) - set(EXPECTED_IDS) - {""})
    missing = sorted(set(EXPECTED_IDS) - set(ids))
    empty_smiles = int(registry["standardized_smiles"].fillna("").astype(str).str.strip().eq("").sum())
    result["registry"] = {"rows": len(registry), "unique_ids": len(set(ids) - {""}),
                          "duplicate_ids_rows": int(duplicates), "empty_smiles": empty_smiles,
                          "missing_ids": missing, "unexpected_ids": unknown}
    if duplicates or unknown or missing or empty_smiles or "" in ids:
        defect(issues, "DATA", "blocker", "REGISTRY_QC",
               "Несоответствие регистра 65 PYZ: дубли/пропуски ID или SMILES", "config/own_compounds.tsv",
               evidence=json.dumps(result["registry"], ensure_ascii=False), repair_stage="ENV")
    try:
        from rdkit import Chem
        invalid = [str(row.own_compound_id) for row in registry.itertuples() if str(getattr(row, "standardized_smiles", "")).strip()
                   and Chem.MolFromSmiles(str(row.standardized_smiles)) is None]
        result["registry"]["rdkit_invalid_smiles"] = invalid
        if invalid:
            defect(issues, "DATA", "major", "INVALID_SMILES", "RDKit не смог разобрать SMILES", "registry", str(invalid), "E1")
    except ImportError:
        result["registry"]["rdkit_smiles_validation"] = "NOT_REVIEWED: RDKit unavailable"

    if matrix is None or not {"own_compound_id", "target_gene", "best_tanimoto_morgan_r2_2048"}.issubset(matrix.columns):
        result["matrix_status"] = "FAIL"
        return result, coverage
    m = matrix.copy()
    m["own_compound_id"] = m["own_compound_id"].fillna("").astype(str).str.strip()
    m["target_gene"] = m["target_gene"].fillna("").astype(str).str.strip().str.upper()
    m["_valid"] = valid_similarity(m["best_tanimoto_morgan_r2_2048"])
    m["_value"] = pd.to_numeric(m["best_tanimoto_morgan_r2_2048"], errors="coerce")
    targets = sorted(set(m["target_gene"]) - {""})
    pairs = m.groupby(["own_compound_id", "target_gene"], dropna=False, sort=True)
    counts = pairs.size().to_dict()
    m_by_pair = {(cid, gene): sub for (cid, gene), sub in pairs}
    invalid_n = int((~m["_valid"]).sum())
    duplicate_pairs = int(sum(n - 1 for n in counts.values() if n > 1))
    unknown_id = sorted(set(m["own_compound_id"]) - set(ids) - {""})
    result["matrix"] = {"rows": len(m), "unique_pairs": len(counts), "targets_n": len(targets),
                       "targets": targets, "duplicates_above_unique": duplicate_pairs,
                       "invalid_similarity_rows": invalid_n, "unknown_pyz_ids": unknown_id,
                       "empty_target_rows": int(m["target_gene"].eq("").sum())}
    if len(targets) != 13:
        defect(issues, "DATA", "major", "TARGET_COUNT", f"Фактически {len(targets)} мишеней вместо ожидаемых 13", "matrix", repair_stage="ENV")
    if invalid_n:
        defect(issues, "DATA", "major", "INVALID_TANIMOTO", f"В матрице {invalid_n} недопустимых значений", "matrix", repair_stage="E1")
    if duplicate_pairs:
        defect(issues, "DATA", "major", "DUPLICATE_PAIRS", f"В матрице {duplicate_pairs} дополнительных строк для повторных пар", "matrix", repair_stage="ENV")
    if unknown_id:
        defect(issues, "DATA", "major", "FOREIGN_PYZ_ID", f"Мишени с неизвестными PYZ: {unknown_id}", "matrix", repair_stage="ENV")
    for cid in EXPECTED_IDS:
        for gene in targets:
            subset = m_by_pair.get((cid, gene))
            n = 0 if subset is None else len(subset)
            if n == 0:
                state, value = "missing", None
            elif n > 1:
                state, value = "duplicate", None
            elif not bool(subset.iloc[0]["_valid"]):
                state, value = "invalid", None
            else:
                state, value = "present", float(subset.iloc[0]["_value"])
            coverage.append({"compound_id": cid, "target_gene": gene, "rows": n, "state": state,
                             "tanimoto": value, "source": "matrix"})
    result["coverage"] = dict(Counter(row["state"] for row in coverage))
    result["coverage"]["possible_pairs_for_actual_target_set"] = len(EXPECTED_IDS) * len(targets)
    result["coverage"]["complete_65_by_13"] = (
        len(targets) == 13 and result["coverage"].get("present", 0) == 845 and
        len(m) == 845 and not invalid_n and not duplicate_pairs and not unknown_id
    )
    if result["coverage"].get("missing", 0):
        defect(issues, "DATA", "major", "MISSING_PAIRS",
               f"Отсутствует {result['coverage']['missing']} сочетаний для фактического набора мишеней", "matrix", repair_stage="ENV")

    neighbors = frames.get("neighbors")
    if neighbors is not None and {"own_compound_id", "target_gene", "rank_in_target"}.issubset(neighbors.columns):
        known_ids = set(ids)
        known_pairs = set(counts)
        pairs_n = list(zip(neighbors["own_compound_id"].fillna("").astype(str),
                           neighbors["target_gene"].fillna("").astype(str).str.upper()))
        foreign_ids = int(sum(cid not in known_ids for cid, _ in pairs_n))
        foreign_pairs = int(sum(pair not in known_pairs for pair in pairs_n))
        invalid_rank = int((~pd.to_numeric(neighbors["rank_in_target"], errors="coerce").ge(1)).sum())
        duplicates_nb = int(neighbors.duplicated(subset=["own_compound_id", "target_gene", "rank_in_target"]).sum())
        result["neighbors"] = {"rows": len(neighbors), "foreign_pyz_rows": foreign_ids,
                               "pairs_absent_from_matrix": foreign_pairs, "invalid_rank_rows": invalid_rank,
                               "duplicate_rank_rows": duplicates_nb,
                               "ligand_id_column": "ligand_id" if "ligand_id" in neighbors.columns else "NOT_AVAILABLE"}
        if any((foreign_ids, foreign_pairs, invalid_rank, duplicates_nb)):
            defect(issues, "DATA", "major", "NEIGHBORS_QC",
                   "В соседях есть незнакомые PYZ/пары, неверный ранг или дубли", "neighbors",
                   json.dumps(result["neighbors"], ensure_ascii=False), "ENV")
    if "admet_individual" in frames:
        admet = frames["admet_individual"]
        if "own_compound_id" in admet.columns:
            ids_admet = admet["own_compound_id"].fillna("").astype(str).tolist()
            result["admet"] = {"rows": len(admet), "unique_ids": len(set(ids_admet)),
                               "duplicates": int(admet["own_compound_id"].duplicated().sum()),
                               "missing_registry_ids": sorted(set(ids) - set(ids_admet)),
                               "unknown_ids": sorted(set(ids_admet) - set(ids))}
            if "qed" in admet.columns:
                result["admet"]["qed_out_of_range"] = int((~valid_similarity(admet["qed"]) & admet["qed"].notna()).sum())
        else:
            defect(issues, "DATA", "major", "ADMET_ID_MISSING", "ADMET без колонки own_compound_id", "admet_individual", repair_stage="E1")
    hypothesis = frames.get("candidate")
    models = frames.get("candidate_models")
    if hypothesis is not None and models is not None and "hypothesis_id" in hypothesis and "hypothesis_id" in models:
        known_hyp = set(hypothesis["hypothesis_id"].dropna().astype(str))
        foreign_hyp = sorted(set(models["hypothesis_id"].dropna().astype(str)) - known_hyp)
        result["candidate"] = {"hypotheses": len(hypothesis), "model_rows": len(models),
                               "unknown_hypothesis_ids_n": len(foreign_hyp), "sample_unknown": foreign_hyp[:10]}
        if foreign_hyp:
            defect(issues, "DATA", "major", "UNLINKED_MODEL", "Клеточные наблюдения без гипотезы", "candidate_models", str(foreign_hyp[:10]), "E1")
        if "model_id" in models and "preferred_name" in models:
            conflicts = models.groupby("model_id", dropna=True)["preferred_name"].nunique()
            result["candidate"]["models_with_multiple_reference_drugs_n"] = int((conflicts > 1).sum())
    return result, coverage


def request_json(url: str, timeout: float) -> tuple[int, Any, float, int, str, str]:
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "MCL-E0-Audit/1.0"})
    started = time.perf_counter()
    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
            raw = response.read()
            return response.status, json.loads(raw), (time.perf_counter() - started) * 1000, len(raw), response.headers.get("X-MCL-Process-Time-Ms", ""), response.headers.get("Cache-Control", "")
    except urllib.error.HTTPError as err:
        raw = err.read()
        try:
            payload = json.loads(raw)
        except ValueError:
            payload = {"detail": raw[:300].decode(errors="replace")}
        return err.code, payload, (time.perf_counter() - started) * 1000, len(raw), "", ""
    except (urllib.error.URLError, OSError, TimeoutError, ValueError) as exc:
        return 0, {"error": f"{type(exc).__name__}: {exc}"}, (time.perf_counter() - started) * 1000, 0, "", ""


def safe_api_sample(key: str, data: Any) -> dict[str, Any]:
    if not isinstance(data, dict):
        return {"response_type": type(data).__name__}
    if key == "health":
        return {k: data.get(k) for k in ("status", "mcl_root")}
    if key == "summary":
        return {k: data.get(k) for k in ("compounds_n", "targets_n", "matrix_rows_n", "shortlist_rows_n", "target_genes")}
    if key == "catalog":
        return {"total": data.get("total"), "items_n": len(data.get("items") or []),
                "first_ids": [i.get("own_compound_id") for i in (data.get("items") or [])[:6]]}
    if key == "matrix":
        return {"items_n": len(data.get("items") or []), "targets_n": len(data.get("targets") or []),
                "compounds_n": len(data.get("compounds") or [])}
    if key == "unknown":
        return {"detail": str(data.get("detail", ""))[:300]}
    profile = data.get("profile_v2") or {}
    return {"compound_id": profile.get("compound_id"), "schema_version": profile.get("schema_version"),
            "identity_id": (profile.get("identity") or {}).get("own_compound_id"),
            "target_items_n": len((profile.get("target_hypotheses") or {}).get("items") or []),
            "context_items_n": len((profile.get("cell_context") or {}).get("items") or []),
            "has_smiles": bool((profile.get("identity") or {}).get("standardized_smiles")),
            "own_target_binding": (profile.get("scientific_summary") or {}).get("own_target_binding"),
            "own_cellular_activity": (profile.get("scientific_summary") or {}).get("own_cellular_activity")}


def audit_api(root: Path, base: str, inventory: dict[str, pd.DataFrame],
              report: Path, issues: list[dict], measure: bool,
              timeout: float) -> tuple[list[dict], list[dict], list[dict]]:
    paths = {"health": "/health", "summary": "/api/pyz/summary",
             "catalog": "/api/pyz?limit=100&offset=0", "matrix": "/api/pyz/matrix",
             "unknown": "/api/pyz/PYZ-999"}
    for cid in EXPECTED_IDS:
        paths[cid] = f"/api/pyz/{cid}"
    cases: list[dict] = []
    observations: list[dict] = []
    reconciliation: list[dict] = []
    expected_registry = {}
    reg = inventory.get("registry")
    if reg is not None and {"own_compound_id", "standardized_smiles"}.issubset(reg.columns):
        expected_registry = dict(zip(reg["own_compound_id"].astype(str), reg["standardized_smiles"].fillna("").astype(str)))
    for key, path in paths.items():
        code, data, duration, bytes_n, header, cache_header = request_json(base.rstrip("/") + path, timeout)
        ok = (code == 404 if key == "unknown" else code == 200)
        entry = {"case": key, "path": path, "http_status": code, "pass": ok,
                 "duration_ms": round(duration, 3), "response_bytes": bytes_n,
                 "process_ms": header, "cache_control": cache_header, "error": ""}
        if code == 0:
            entry["error"] = str(data.get("error", "unavailable"))
        if key == "health" and code == 200:
            actual_root = Path(str(data.get("mcl_root", ""))).resolve()
            if actual_root != root.resolve() or data.get("status") != "ok":
                ok = False
                entry["pass"] = False
                entry["error"] = "MCL_ROOT отличается от аудируемой папки или status != ok"
                defect(issues, "ENV", "blocker", "WRONG_ROOT", entry["error"], path, str(actual_root), "ENV")
                cases.append(entry)
                break  # Do not audit a foreign instance.
        if code == 200 and key == "summary" and isinstance(data, dict):
            registry_df = inventory.get("registry")
            matrix_df = inventory.get("matrix")
            if registry_df is not None and "compounds_n" in data and int(data["compounds_n"]) != len(registry_df):
                entry["pass"] = False
                entry["error"] = "summary compounds_n disagrees with registered source rows"
            if matrix_df is not None and "matrix_rows_n" in data and int(data["matrix_rows_n"]) != len(matrix_df):
                entry["pass"] = False
                entry["error"] = "summary matrix_rows_n disagrees with source rows"
        if code == 200 and key == "catalog" and isinstance(data, dict):
            registry_df = inventory.get("registry")
            if registry_df is not None and data.get("total") != len(registry_df):
                entry["pass"] = False
                entry["error"] = "catalog total disagrees with registry source rows"
        if code == 200 and key == "matrix" and isinstance(data, dict):
            matrix_df = inventory.get("matrix")
            if matrix_df is not None and len(data.get("items") or []) != len(matrix_df):
                entry["pass"] = False
                entry["error"] = "matrix API items_n disagrees with source rows"
        if code != 0 and code == 200 and key.startswith("PYZ-"):
            p = (data or {}).get("profile_v2") or {}
            actual_id = p.get("compound_id")
            actual_smiles = (p.get("identity") or {}).get("standardized_smiles")
            summary_data = p.get("scientific_summary") or {}
            entry["pass"] = bool(p.get("schema_version") == PROFILE_VERSION and actual_id == key and
                                 (p.get("identity") or {}).get("own_compound_id") == key and
                                 (not expected_registry or expected_registry.get(key) == actual_smiles) and
                                 summary_data.get("own_target_binding") == "not_measured" and
                                 summary_data.get("own_cellular_activity") == "not_measured")
            if key in ("PYZ-001", "PYZ-035", "PYZ-051", "PYZ-065"):
                reconciliation.append({"compound_id": key, "field": "identity.standardized_smiles",
                                       "source": "registry", "api_matches": expected_registry.get(key) == actual_smiles,
                                       "ui_status": "NOT_REVIEWED", "notes": "Нужно проверить изображение структуры и SMILES на экране"})
        if not entry["pass"]:
            defect(issues, "API", "major" if code else "blocker", "API_CASE_FAIL",
                   f"API {key} вернул HTTP {code} или некорректный контракт", path, entry["error"], "ENV" if code == 0 else "E1")
        cases.append(entry)
        if key in ("health", "summary", "catalog", "matrix", "unknown", "PYZ-001", "PYZ-035", "PYZ-051", "PYZ-065"):
            save_json(report / "api_samples" / f"{key.lower()}.json", safe_api_sample(key, data))
        if code == 0:
            break  # Avoid 65 identical connection timeouts on a stopped server.
    if measure and cases and cases[0]["pass"]:
        endpoints = {"summary": paths["summary"], "catalog": paths["catalog"], "card": paths["PYZ-035"], "matrix": paths["matrix"]}
        for label, path in endpoints.items():
            endpoint = base.rstrip("/") + path
            # The first isolated request (cold-ish) is retained separately; true process cold-start is manual.
            cold = request_json(endpoint, timeout)
            observations.append({"route": label, "phase": "first_sample_not_process_cold",
                                 "n": 0, "http_status": cold[0], "client_ms": round(cold[2], 3),
                                 "response_bytes": cold[3], "process_ms": cold[4]})
            for _ in range(3):
                request_json(endpoint, timeout)
            for n in range(1, 21):
                status, _, dur, size, process_ms, _ = request_json(endpoint, timeout)
                observations.append({"route": label, "phase": "warm", "n": n, "http_status": status,
                                     "client_ms": round(dur, 3), "response_bytes": size, "process_ms": process_ms})
    return cases, observations, reconciliation


def percentile(values: list[float], q: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    pos = (len(ordered) - 1) * q
    lo, hi = math.floor(pos), math.ceil(pos)
    return ordered[lo] + (ordered[hi] - ordered[lo]) * (pos - lo)


def write_reports(root: Path, out: Path, git_sha: str, git_status_before: str,
                  inventory: list[dict], frames: dict[str, pd.DataFrame], hashes: dict[str, str],
                  data_summary: dict, coverage: list[dict], issues: list[dict],
                  api_cases: list[dict], perf: list[dict], reconciliation: list[dict],
                  offline: bool, ui_base: str, api_base: str) -> None:
    save_tsv(out / "input_inventory.tsv", ["label", "path", "required", "exists", "format", "bytes",
             "modified_utc", "sha256", "rows", "columns_n", "column_types", "nulls_key",
             "duplicate_keys", "status", "error"], inventory)
    atomic_text(out / "file_checksums.sha256", "".join(f"{h}  {path}\n" for path, h in hashes.items()))
    save_json(out / "schema_audit.json", {"version": AUDITOR_VERSION, "summary": data_summary,
                                         "issues": issues, "schemas": {k: {name: str(dtype) for name, dtype in df.dtypes.items()} for k, df in frames.items()}})
    save_tsv(out / "coverage_matrix.tsv", ["compound_id", "target_gene", "rows", "state", "tanimoto", "source"], coverage)
    save_tsv(out / "api_cases.tsv", ["case", "path", "http_status", "pass", "duration_ms", "response_bytes", "process_ms", "cache_control", "error"], api_cases)
    save_tsv(out / "performance.tsv", ["route", "phase", "n", "http_status", "client_ms", "response_bytes", "process_ms"], perf)
    save_tsv(out / "reconciliation.tsv", ["compound_id", "field", "source", "api_matches", "ui_status", "notes"], reconciliation)
    for file, content in {
        "data_qc.md": ["# E0 — проверка данных", "", "Аудитор проверил фактические локальные файлы только на чтение.",
                       "Сводка: ```json", json.dumps(data_summary, ensure_ascii=False, indent=2), "```"],
        "performance_method.md": ["# E0 — методика скорости", "", "Режим: клиентское время HTTP-запроса, в том числе сериализация и передача.",
                                  "Для каждого маршрута: 1 первый образец, 3 прогрева (не входят в статистику), 20 последовательных тёплых измерений.",
                                  "Первый образец не считается холодным запуском процесса; полноценный cold start требует ручного измерения PS1/PS2.",
                                  "Кэш клиента/Next.js и браузера необходимо отдельно проверить вручную. API URL: " + api_base,
                                  "", "| Маршрут | Успешных замеров | p50, мс | p95, мс |", "|---|---:|---:|---:|"],
    }.items():
        if file == "performance_method.md":
            for name in ("summary", "catalog", "card", "matrix"):
                arr = [float(i["client_ms"]) for i in perf if i["route"] == name and i["phase"] == "warm" and i["http_status"] == 200]
                content.append(f"| {name} | {len(arr)} | {percentile(arr, 0.5)} | {percentile(arr, 0.95)} |")
        atomic_text(out / file, "\n".join(content) + "\n")
    defects = ["# E0 — реестр наблюдений и дефектов", "", "Не означает автоматического принятия Э0.", ""]
    for row in issues:
        defects.extend([f"## {row['id']} · {row['category']} · {row['severity']} · {row['code']}",
                        f"- Факт: {row['message']}", f"- Источник: `{row['source']}`",
                        f"- Доказательство: {row['evidence']}", f"- Этап устранения: {row['repair_stage']}",
                        "- Статус: OPEN", ""])
    if not issues:
        defects.append("Автоматически выявленных дефектов нет. Ручная и научная приёмка ещё требуются.\n")
    atomic_text(out / "defects.md", "\n".join(defects))
    critical = [x for x in issues if x["severity"] == "blocker"]
    manual_routes = ["/pyz", "/pyz/PYZ-035", "/pyz/matrix"]
    screen_sizes = [(1440, 900), (768, 1024), (390, 844)]
    checklist = [{"route": route, "width": w, "height": h, "status": "NOT_REVIEWED",
                  "screenshot": "", "notes": "Проверить прокрутку, доступность с клавиатуры, данные, клики"}
                 for route in manual_routes for w, h in screen_sizes]
    save_tsv(out / "ui_checklist.tsv", ["route", "width", "height", "status", "screenshot", "notes"], checklist)
    (out / "screenshots").mkdir(exist_ok=True)
    acceptance = ["# E0 — акт приёмки (ЗАПОЛНИТЬ ПОСЛЕ РУЧНОЙ ПРОВЕРКИ)", "",
                  "**Статус: НЕ ПРИНЯТО / ожидает ручной и научной проверки.**", "",
                  f"Базовый Git SHA: `{git_sha}`", f"Рабочий каталог: `{root}`",
                  f"API режим: {'OFFLINE' if offline else 'LIVE'}. URL: `{api_base}`",
                  "", "| Критерий | Текущее состояние | Что требуется |", "|---|---|---|",
                  "| AC-01 | PARTIAL | Проверить среды, SHA и исходные хеши |",
                  "| AC-02 | PARTIAL | Проверить полноту используемых таблиц и схем |",
                  "| AC-03 | PARTIAL | Сверить число уникальных пар и покрытие |",
                  "| AC-04 | " + ("NOT_RUN" if offline else "PARTIAL") + " | API все 65 + PYZ-999; контрольные ответы |",
                  "| AC-05 | NOT_REVIEWED | Добавить 9 исходных скриншотов и ручные проверки |",
                  "| AC-06 | " + ("NOT_RUN" if offline or not perf else "PARTIAL") + " | Подтвердить p50/p95, cold start и UI |",
                  "| AC-07 | NOT_REVIEWED | Ручная научная сверка источников/методов |",
                  "| AC-08 | PARTIAL | Проверить входные SHA до/после |",
                  "| AC-09 | NOT_APPROVED | Утверждение владельцем проекта |",
                  "", f"Автоматически выявлено блокеров: {len(critical)}.",
                  "", "**Начало Э1 не разрешено автоматически.** Требуется устранить блокеры среды и утвердить дефекты Э1."]
    atomic_text(out / "E0_ACCEPTANCE.md", "\n".join(acceptance) + "\n")
    e1 = ["# E1 — входные наблюдения", "", "Требуют проверки перед оформлением отдельного научного ТЗ.", "",
          "По коду нынешней карточки известны потенциальные риски: смешение Ki/Kd/IC50 в агрегате известных лигандов,",
          "упрощение ролей клеточных линий, отсутствие индивидуального Ames/DILI в profile_v2, жестко заданные числа в UI.",
          "Это **не** подтверждает конкретные ошибки в локальных данных до ручной сверки.", ""]
    for row in issues:
        if row["repair_stage"] == "E1":
            e1.append(f"- {row['id']}: {row['message']}")
    atomic_text(out / "E1_INPUTS.md", "\n".join(e1) + "\n")
    before = set(git_status_before.splitlines())
    after_status = git(root, "status", "--porcelain", "--untracked-files=no")
    after = set(after_status.splitlines())
    atomic_text(out / "git_baseline.txt", "\n".join([f"SHA={git_sha}",
                   "STATUS_BEFORE:", git_status_before, "STATUS_AFTER_TRACKED:", after_status,
                   "NEW_TRACKED_CHANGES:", *sorted(after - before)]) + "\n")
    save_json(out / "environment.json", {"version": AUDITOR_VERSION, "when_utc": utc_now(),
         "root": str(root), "git_sha": git_sha, "git_branch": git(root, "branch", "--show-current"),
         "python": sys.version, "python_executable": sys.executable,
         "platform": platform.platform(), "api_base": api_base, "ui_base": ui_base,
         "offline": offline, "git_status_before": git_status_before,
         "source_sha_before": hashes, "issues_n": len(issues)})


def audit(args: argparse.Namespace) -> tuple[Path, int]:
    root = Path(args.root).expanduser().resolve()
    if not (root / "config" / "own_compounds.tsv").exists() or not (root / "apps" / "api").is_dir():
        raise ValueError(f"Not an MCL root: {root}")
    sha = git(root, "rev-parse", "HEAD")
    short = sha[:8] if re.fullmatch(r"[0-9a-f]{40}", sha) else "gitunknown"
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    report_root = Path(getattr(args, "report_root", None) or Path(__file__).resolve().parents[1]).expanduser().resolve()
    if not (report_root / "scripts" / "audit_pyz_e0.py").is_file() and not getattr(args, "report_root", None):
        raise ValueError(f"Cannot identify audit code root: {report_root}")
    parent = report_root / "outputs" / "audits" / "pyz_e0"
    parent.mkdir(parents=True, exist_ok=True)
    out = parent / f"{stamp}-{short}"
    count = 1
    while out.exists():
        out = parent / f"{stamp}-{short}-{count}"
        count += 1
    out.mkdir()
    status = git(root, "status", "--porcelain", "--untracked-files=no")
    issues: list[dict] = []
    inventory, frames, hashes = inventory_sources(root, issues)
    summary, coverage = inspect_data(frames, issues)
    cases: list[dict] = []
    perf: list[dict] = []
    reconciliation: list[dict] = []
    if not args.offline:
        cases, perf, reconciliation = audit_api(root, args.api_base, frames, out, issues, not args.skip_performance, args.timeout)
    else:
        reconciliation = [{"compound_id": x, "field": "identity.standardized_smiles",
                           "source": "registry", "api_matches": "NOT_RUN", "ui_status": "NOT_REVIEWED",
                           "notes": "Offline audit"} for x in ("PYZ-001", "PYZ-035", "PYZ-051", "PYZ-065")]
    # Hash check: read-only evidence. Avoid claims when source could not be hashed.
    changed = []
    for rel, expected in hashes.items():
        path = root / rel
        if not path.is_file() or sha256(path) != expected:
            changed.append(rel)
    if changed:
        defect(issues, "DATA", "blocker", "INPUT_MUTATED", "Входные файлы изменились во время аудита", ", ".join(changed), repair_stage="ENV")
    # Existing Gene Explorer validator: read-only, capture real status.
    verifier = root / "scripts" / "verify_runtime_indexes.py"
    if verifier.exists():
        code, output = run_cmd([sys.executable, str(verifier)], root)
        atomic_text(out / "gene_index_verify.log", f"exit={code}\n{output}\n")
        if code:
            defect(issues, "DATA", "major", "GENE_INDEX_VERIFY", "Не пройдена проверка индексов Gene Explorer", "scripts/verify_runtime_indexes.py", output[-1500:], "ENV")
    write_reports(root, out, sha, status, inventory, frames, hashes, summary, coverage,
                  issues, cases, perf, reconciliation, args.offline, args.ui_base, args.api_base)
    return out, (2 if any(x["severity"] == "blocker" for x in issues) else 0)


def main() -> int:
    parser = argparse.ArgumentParser(description="MCL PYZ E0 baseline audit (strictly read-only inputs)")
    parser.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--report-root", default=None,
                        help="Separate audit-code worktree where reports are created (default: script root)")
    parser.add_argument("--api-base", default="http://127.0.0.1:8100")
    parser.add_argument("--ui-base", default="http://127.0.0.1:3100")
    parser.add_argument("--offline", action="store_true", help="Skip all live HTTP requests")
    parser.add_argument("--skip-performance", action="store_true")
    parser.add_argument("--timeout", type=float, default=12.0)
    args = parser.parse_args()
    try:
        out, code = audit(args)
    except (OSError, ValueError) as exc:
        print(f"E0 cannot start safely: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 3
    print(f"E0 outputs: {out}")
    print(f"E0 audit script exit code: {code}. E0 is NOT automatically accepted; see E0_ACCEPTANCE.md")
    return code


if __name__ == "__main__":
    raise SystemExit(main())