"""MCL E0 AC-02: immutable genetic source provenance and axis audit.

No source modifications, no API traffic, no rebuilding scientific data.
The report distinguishes technical data alignment from biological validity.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

CONTRACT = "mcl-e0-ac02-genetic-source-audit-v1"
BASE = Path("data/runtime/explorer")
TABLES = {
    "gene_catalog": ("gene_symbol", "mcl_domains_json", "protein_classes_json", "compartments_json", "hallmarks_json"),
    "gene_context_metrics": ("gene_symbol", "comparison_id", "delta_gene_effect"),
    "gene_annotations": ("gene_symbol", "annotation_type", "annotation_id"),
    "gene_dependency_summary": ("gene_symbol", "dependency_models_n", "dependency_models_total_n",
                                "dependency_fraction", "dependency_type", "best_cancer_name",
                                "best_cancer_dependency_fraction", "specificity_score", "specificity_label_ru"),
}
LAYERS = ("gene_effect", "expression", "copy_number")
CONTROL_GENES = ("CDK4", "BCL2L1", "EGFR", "TP53", "KRAS")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def safe_path(root: Path, relative: str | Path) -> Path:
    rel = Path(relative)
    if rel.is_absolute() or ".." in rel.parts:
        raise ValueError(f"Unsafe relative path: {relative}")
    path = root / rel
    if path.is_symlink() or not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Unsafe source path: {relative}")
    return path


def read_json(path: Path, *, max_size: int = 16_000_000) -> dict[str, Any]:
    if not path.is_file() or path.stat().st_size > max_size:
        raise ValueError(f"Missing or oversized metadata: {path.name}")
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(data, dict):
        raise ValueError(f"Invalid metadata object: {path.name}")
    return data


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        w = csv.DictWriter(handle, delimiter="\t", fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({col: "" if row.get(col) is None else str(row.get(col)) for col in columns})


INVENTORY_COMPARE_FIELDS = (
    "depmap_release", "logical_name", "filename", "size_bytes", "sha256"
)


def read_inventory_rows(path: Path) -> list[tuple[str, ...]]:
    """Parse the five scientific-identity fields; ignore machine-specific paths."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not set(INVENTORY_COMPARE_FIELDS).issubset(reader.fieldnames):
            raise ValueError("DepMap inventory lacks required scientific-identity columns")
        rows: list[tuple[str, ...]] = []
        for entry in reader:
            values = tuple(str(entry.get(k) or "").strip() for k in INVENTORY_COMPARE_FIELDS)
            if not all(values):
                raise ValueError("DepMap inventory has blank scientific-identity fields")
            if not values[3].isdigit() or int(values[3]) < 1:
                raise ValueError("DepMap inventory has invalid source file size")
            if not re.fullmatch(r"[a-fA-F0-9]{64}", values[4]):
                raise ValueError("DepMap inventory has invalid source SHA256")
            rows.append((*values[:3], str(int(values[3])), values[4].lower()))
    if not rows:
        raise ValueError("DepMap inventory is empty")
    return rows


def compare_inventories(scientific: Path, reference: Path) -> dict[str, Any]:
    """Compare source identities irrespective of TSV ordering or absolute paths."""
    left = Counter(read_inventory_rows(scientific))
    right = Counter(read_inventory_rows(reference))
    return {
        "scientific_rows": sum(left.values()),
        "reference_rows": sum(right.values()),
        "matched_rows": sum((left & right).values()),
        "scientific_only_rows": sum((left - right).values()),
        "reference_only_rows": sum((right - left).values()),
        "scientific_duplicate_rows": sum(n - 1 for n in left.values()),
        "reference_duplicate_rows": sum(n - 1 for n in right.values()),
        "semantic_match": left == right,
        "full_tsv_sha_match": sha256_file(scientific) == sha256_file(reference),
        "compared_fields": list(INVENTORY_COMPARE_FIELDS),
    }


def recorded_hashes(file: Path | None) -> dict[str, str]:
    if file is None:
        return {}
    if not file.is_file():
        raise FileNotFoundError(f"E0 reference hashes not found: {file}")
    ret: dict[str, str] = {}
    for line in file.read_text(encoding="utf-8-sig").splitlines():
        if not line.strip():
            continue
        m = re.fullmatch(r"([0-9a-fA-F]{64})  (.+)", line)
        if not m:
            raise ValueError("Malformed baseline SHA-256 listing")
        rel = m.group(2).replace("\\", "/")
        safe_path(Path.cwd(), rel)  # reject path traversal, without needing source to exist
        if rel in ret:
            raise ValueError("Duplicate baseline SHA path")
        ret[rel] = m.group(1).lower()
    return ret


def pick_indices(length: int, count: int = 11) -> list[int]:
    if length <= 0:
        return []
    n = min(length, count)
    return sorted({(i * (length - 1)) // max(1, n - 1) for i in range(n)})


def sample_gene_indices(genes: list[str]) -> list[int]:
    keep = set(pick_indices(len(genes), 11))
    gene_map = {x: i for i, x in enumerate(genes)}
    keep.update(gene_map[g] for g in CONTROL_GENES if g in gene_map)
    return sorted(keep)


def canon(values: list[Any]) -> list[str]:
    return [str(x).strip().upper() for x in values]


def axis_quality(models: list[Any], genes: list[Any]) -> list[str]:
    errors = []
    model_ids = [str(x).strip() for x in models]
    gene_ids = canon(genes)
    if not model_ids or not gene_ids:
        errors.append("EMPTY_AXIS")
    if any(not x for x in model_ids) or any(not x for x in gene_ids):
        errors.append("BLANK_AXIS_KEY")
    if len(model_ids) != len(set(model_ids)):
        errors.append("DUPLICATED_MODEL_ID")
    if len(gene_ids) != len(set(gene_ids)):
        errors.append("DUPLICATED_GENE_SYMBOL")
    return errors


def verify_values(matrix: np.ndarray, source: pd.DataFrame, genes: list[str],
                  models: list[str], indices: list[int]) -> tuple[int, int]:
    """Independent numerical crosscheck vs builder source (float32 conversion).

    Only selected Parquet columns and selected model rows are inspected.
    Pure NaN/NaN equals; infinity is never accepted as trustworthy data.
    """
    if source["model_id"].astype(str).tolist() != models:
        return 0, -1
    model_points = pick_indices(len(models), 9)
    good = total = 0
    for pos in indices:
        gene = genes[pos]
        # Column names may differ in case, but must be unique after normalization.
        matching = [c for c in source.columns if c != "model_id" and str(c).strip().upper() == gene]
        if len(matching) != 1:
            return good, -2
        observed = source[matching[0]].to_numpy(dtype=np.float32, copy=False)
        expected = matrix[pos, model_points]
        for a, b in zip(observed[model_points], expected):
            total += 1
            equal = (np.isnan(a) and np.isnan(b)) or (
                np.isfinite(a) and np.isfinite(b) and a == b)
            good += int(equal)
    return good, total


def inspect_layer(root: Path, name: str, add_issue, inventory, *,
                  baseline: dict[str, str], comparison_root: Path | None = None) -> dict[str, Any]:
    prefix = BASE / "model_layers" / name
    rel_npy = f"{prefix}.npy".replace("\\", "/")
    rel_json = f"{prefix}.json".replace("\\", "/")
    matrix_path = safe_path(root, rel_npy)
    json_path = safe_path(root, rel_json)
    result: dict[str, Any] = {"layer": name, "axis_status": "NOT_REVIEWED", "models_n": 0,
                              "genes_n": 0, "dtype": "", "orientation": "", "sample_compared": 0,
                              "sample_matched": 0, "source_parquet_status": "NOT_REVIEWED"}
    if not matrix_path.is_file() or not json_path.is_file():
        add_issue("BLOCKER", name, "MISSING_LAYER_OR_AXES", "Layer array or its ordered axes JSON is missing")
        return result
    try:
        meta = read_json(json_path)
        models = [str(x) for x in meta.get("models", [])]
        genes = canon(meta.get("genes", []))
        for issue in axis_quality(models, genes):
            add_issue("BLOCKER", name, issue, "Duplicated/empty ordered gene or model identifiers")
        if meta.get("orientation") != "gene_by_model" or str(meta.get("schema_version")) != "1.0":
            add_issue("BLOCKER", name, "WRONG_LAYER_CONTRACT", "Expected schema 1.0, gene_by_model")
        if str(meta.get("layer")) != name:
            add_issue("BLOCKER", name, "LAYER_NAME_MISMATCH", "Metadata layer does not match array filename")
        if meta.get("dtype") != "float32":
            add_issue("BLOCKER", name, "METADATA_DTYPE_MISMATCH", "Expected float32 in builder contract")
        if int(meta.get("models_n", -1)) != len(models) or int(meta.get("genes_n", -1)) != len(genes):
            add_issue("BLOCKER", name, "METADATA_AXIS_COUNTS", "Axis count differs from ordered IDs")
        matrix = np.load(matrix_path, mmap_mode="r", allow_pickle=False)
        result.update(models_n=len(models), genes_n=len(genes), dtype=str(matrix.dtype),
                      orientation=str(meta.get("orientation")), shape=list(matrix.shape))
        if matrix.ndim != 2 or matrix.shape != (len(genes), len(models)) or matrix.dtype != np.float32:
            add_issue("BLOCKER", name, "ARRAY_SHAPE_OR_DTYPE", "NumPy array is not float32 genes × models")
            return result
        result["axis_status"] = "AXES_AND_SHAPE_CHECKED"
        expected_source_name = f"depmap_model_{name}.parquet"
        if meta.get("source_file") != expected_source_name or meta.get("matrix_file") != matrix_path.name:
            add_issue("BLOCKER", name, "SOURCE_FILE_METADATA", "Layer source_file/matrix_file metadata mismatches builder")
        source_path = safe_path(comparison_root or root, Path("data/processed") / expected_source_name)
        if not source_path.is_file():
            add_issue("REVIEW", name, "MISSING_PROCESSED_MATRIX", "Cannot independently verify ordered axes/values without builder source")
            result["source_parquet_status"] = "MISSING"
            return result
        # The source is intentionally not exported and never written.
        import pyarrow.parquet as pq
        pf = pq.ParquetFile(source_path)
        cols = list(pf.schema_arrow.names)
        if "model_id" not in cols:
            add_issue("BLOCKER", name, "MISSING_MODEL_ID", "Builder source has no model_id")
            return result
        source_genes = [str(c).strip().upper() for c in cols if c != "model_id"]
        if genes != source_genes:
            add_issue("BLOCKER", name, "GENE_COLUMN_ORDER_MISMATCH", "Array gene axis differs from source Parquet order")
            result["source_parquet_status"] = "GENE_ORDER_MISMATCH"
            return result
        indices = sample_gene_indices(genes)
        selected = [cols[1 + i] for i in indices]  # requires model_id to be first in source schema
        if cols[0] != "model_id":
            selected = [c for c in cols if c != "model_id" and c.upper() in {genes[i] for i in indices}]
        sample_source = pd.read_parquet(source_path, columns=["model_id", *selected])
        # Full model axis is checked without reading all gene columns.
        good, count = verify_values(matrix, sample_source, genes, models, indices)
        result["sample_compared"] = count if count >= 0 else 0
        result["sample_matched"] = good
        if count < 0:
            add_issue("BLOCKER", name, "MODEL_AXIS_OR_GENE_COLUMN_MISMATCH",
                      "Builder source model order or sample gene column is inconsistent")
        elif good != count:
            add_issue("BLOCKER", name, "SAMPLED_VALUES_MISMATCH", f"Float32 scientific values differ: {good}/{count} matched")
        else:
            result["source_parquet_status"] = (
                "ORDER_AND_SAMPLE_MATCH_CROSS_ROOT"
                if comparison_root is not None and comparison_root.resolve() != root.resolve()
                else "ORDER_AND_SAMPLE_MATCH"
            )
        return result
    except Exception as ex:
        add_issue("BLOCKER", name, "UNREADABLE_LAYER", f"Unable to validate layer or its axes: {type(ex).__name__}: {str(ex)[:130]}")
        return result


def audit(root: Path, baseline_file: Path | None = None, *,
          comparison_root: Path | None = None) -> tuple[dict[str, Any], list[dict], list[dict], list[dict]]:
    root = root.resolve()
    comparison_root = comparison_root.resolve() if comparison_root is not None else root
    if not comparison_root.is_dir():
        raise FileNotFoundError("Comparison source root does not exist")
    if not safe_path(root, BASE).is_dir():
        raise FileNotFoundError("Explorer runtime data directory is missing")
    baseline = recorded_hashes(baseline_file)
    inventory: list[dict] = []
    issues: list[dict] = []
    def issue(severity: str, component: str, code: str, detail: str):
        issues.append({"severity": severity, "component": component, "code": code, "detail": detail})
    expected = [*(f"{BASE / k}.parquet" for k in TABLES), *(f"{BASE / 'model_layers' / k}.npy" for k in LAYERS)]
    for entry in expected:
        rel = str(entry).replace("\\", "/")
        path = safe_path(root, rel)
        if not path.is_file():
            issue("BLOCKER", rel, "MISSING_FILE", "One of seven E0 source files is missing")
            inventory.append({"path": rel, "exists": False, "bytes": 0, "sha256": "", "baseline_match": "NOT_AVAILABLE"})
            continue
        h = sha256_file(path)
        baseline_status = "MATCH" if baseline.get(rel) == h else "MISMATCH" if rel in baseline else "NOT_AVAILABLE"
        if baseline_status == "MISMATCH":
            issue("BLOCKER", rel, "E0_HASH_CHANGED", "Source SHA differs from previously verified E0 baseline")
        inventory.append({"path": rel, "exists": True, "bytes": path.stat().st_size, "sha256": h,
                          "baseline_match": baseline_status})

    # Axis sidecars are part of the scientific interpretation, despite not being
    # among the seven large datasets; detect concurrent changes explicitly.
    sidecar_paths = [
        BASE / "manifest.json", BASE / "gene_dependency_summary.json",
        *(BASE / "model_layers" / f"{name}.json" for name in LAYERS),
        Path("data/processed/depmap_model_multiomics_manifest.json"),
    ]
    sidecar_rows: list[dict[str, Any]] = []
    for relative in sidecar_paths:
        rel = str(relative).replace("\\", "/")
        p = safe_path(root, rel)
        sidecar_rows.append({"path": rel, "exists": p.is_file(),
                             "sha256": sha256_file(p) if p.is_file() else "", "unchanged_after": False})
    index_manifest_path = safe_path(root, BASE / "manifest.json")
    dependency_meta_path = safe_path(root, BASE / "gene_dependency_summary.json")
    processed_multiomics_manifest = safe_path(root, "data/processed/depmap_model_multiomics_manifest.json")
    manifest: dict[str, Any] = {}
    depmeta: dict[str, Any] = {}
    multiomics: dict[str, Any] = {}
    for name, file in (("explorer_manifest", index_manifest_path),
                       ("dependency_manifest", dependency_meta_path),
                       ("multiomics_manifest", processed_multiomics_manifest)):
        if not file.is_file():
            issue("REVIEW" if name == "multiomics_manifest" else "BLOCKER", name,
                  "MANIFEST_MISSING", f"Manifest file is absent: {file.name}")
            continue
        try:
            data = read_json(file)
            if name == "explorer_manifest": manifest = data
            elif name == "dependency_manifest": depmeta = data
            else: multiomics = data
        except (ValueError, OSError, json.JSONDecodeError) as exc:
            issue("REVIEW" if name == "multiomics_manifest" else "BLOCKER", name, "MANIFEST_UNAVAILABLE",
                  f"Cannot read manifest: {type(exc).__name__}")

    if manifest:
        if (manifest.get("index_contract") != "mcl-gene-explorer-runtime-v1" or
                str(manifest.get("schema_version")) != "1.0" or
                manifest.get("runtime_root") != "data/runtime/explorer"):
            issue("BLOCKER", "explorer_manifest", "INDEX_CONTRACT_MISMATCH", "Unexpected index contract/version/root")
        if not manifest.get("built_at"):
            issue("REVIEW", "explorer_manifest", "BUILD_TIMESTAMP_MISSING", "No UTC build provenance")
    if depmeta:
        method = depmeta.get("dependency_call_method")
        if not method:
            issue("REVIEW", "dependency_manifest", "DEPENDENCY_METHOD_MISSING", "No probability/fallback method label")
        elif method == "gene_effect_descriptive_fallback":
            issue("REVIEW", "dependency_manifest", "DESCRIPTIVE_FALLBACK",
                  "Gene Effect cutoff used instead of probability-of-dependency calls; must be labelled in science UI")
    tables: list[dict] = []
    catalog_genes: set[str] = set()
    counts: dict[str, int] = {}
    for name, columns in TABLES.items():
        path = safe_path(root, BASE / f"{name}.parquet")
        if not path.is_file():
            tables.append({"table": name, "rows": 0, "genes": 0, "null_genes": 0, "duplicate_keys": 0})
            continue
        try:
            import pyarrow.parquet as pq
            meta = pq.ParquetFile(path)
            available = set(meta.schema_arrow.names)
            absent = set(columns) - available
            if absent:
                issue("BLOCKER", name, "MISSING_COLUMNS", ", ".join(sorted(absent)))
                continue
            key_columns = ["gene_symbol"]
            if name == "gene_context_metrics": key_columns.append("comparison_id")
            if name == "gene_annotations":
                key_columns.extend(["annotation_type", "annotation_id"])
                if {"source", "source_id"}.issubset(available):
                    key_columns.extend(["source", "source_id"])
                else:
                    issue("REVIEW", name, "ANNOTATION_PROVENANCE_FIELDS_MISSING",
                          "Cannot distinguish repeated mappings by original source and source_id")
            if name == "gene_dependency_summary": key_columns.extend(["dependency_models_n", "dependency_models_total_n", "dependency_fraction"])
            frame = pd.read_parquet(path, columns=key_columns)
            raw_genes = frame["gene_symbol"].astype("string")
            keys = raw_genes.str.strip().str.upper()
            null_n = int(keys.isna().sum() + keys.eq("").sum())
            genes = set(keys.dropna().astype(str)) - {""}
            dup_key_columns = [keys.rename("gene_symbol")]
            if name == "gene_context_metrics": dup_key_columns.append(frame["comparison_id"])
            coarse_annotation_repeats = 0
            if name == "gene_annotations":
                dup_key_columns.extend([frame["annotation_type"], frame["annotation_id"]])
                coarse_annotation_repeats = int(pd.concat(dup_key_columns, axis=1).duplicated().sum())
                if {"source", "source_id"}.issubset(frame.columns):
                    # This is the *actual materialization key* used by
                    # AnnotatedGeneExplorerStore.materialize_indexes().
                    # Multiple GO/Reactome/other terms may legitimately
                    # map to the same (gene, annotation_type, annotation_id).
                    dup_key_columns.extend([frame["source"], frame["source_id"]])
            duplicate = int(pd.concat(dup_key_columns, axis=1).duplicated().sum())
            if null_n:
                issue("BLOCKER", name, "EMPTY_GENE_SYMBOL", f"Null/blank gene symbols: {null_n}")
            if duplicate:
                issue("REVIEW" if name == "gene_annotations" else "BLOCKER", name,
                      "DUPLICATE_ANNOTATION_EVIDENCE" if name == "gene_annotations" else "DUPLICATE_KEYS",
                      f"Repeated full materialization identifier keys: {duplicate}")
            if name == "gene_catalog": catalog_genes = genes
            if name == "gene_dependency_summary":
                d = pd.to_numeric(frame["dependency_models_n"], errors="coerce")
                n = pd.to_numeric(frame["dependency_models_total_n"], errors="coerce")
                f = pd.to_numeric(frame["dependency_fraction"], errors="coerce")
                invalid = ((d < 0) | (n < 0) | (d > n) | ((f < 0) | (f > 1)) & f.notna()).fillna(False)
                if int(invalid.sum()):
                    issue("BLOCKER", name, "INVALID_DEPENDENCY_FRACTION", f"Out-of-range dependency records: {int(invalid.sum())}")
            tables.append({"table": name, "rows": int(meta.metadata.num_rows), "genes": len(genes),
                           "null_genes": null_n, "duplicate_keys": duplicate,
                           "coarse_annotation_repeats": coarse_annotation_repeats})
            counts[name] = len(genes)
            if name == "gene_catalog" and manifest:
                if manifest.get("genes_n") != len(genes):
                    issue("BLOCKER", name, "GENE_COUNT_MANIFEST", "Gene catalog size differs from manifest genes_n")
            if name == "gene_context_metrics" and manifest and manifest.get("context_rows_n") != len(frame):
                issue("BLOCKER", name, "CONTEXT_ROWS_MANIFEST", "Context row count differs from manifest")
        except Exception as exc:
            issue("BLOCKER", name, "UNREADABLE_PARQUET", f"Cannot check Parquet: {type(exc).__name__}: {str(exc)[:120]}")
    if catalog_genes:
        for name in ("gene_context_metrics", "gene_annotations", "gene_dependency_summary"):
            path = safe_path(root, BASE / f"{name}.parquet")
            if not path.is_file(): continue
            try:
                symbols = pd.read_parquet(path, columns=["gene_symbol"])["gene_symbol"].astype(str).str.upper()
                orphan = len(set(symbols) - catalog_genes)
                if orphan:
                    issue("REVIEW", name, "GENES_OUTSIDE_CATALOG", f"Genes not in catalog: {orphan} distinct")
            except Exception: pass
    # Supplemental provenance for the processed Parquet files used as comparison reference.
    # This is NOT a claim that the independent working copy is a certified DepMap snapshot.
    comparison_files: list[dict[str, Any]] = []
    for rel_path in [
        *(f"data/processed/depmap_model_{name}.parquet" for name in LAYERS),
        "data/processed/depmap_model_multiomics_manifest.json",
        "data/processed/depmap_input_inventory.tsv",
    ]:
        p = safe_path(comparison_root, rel_path)
        present = p.is_file()
        comparison_files.append({
            "path": rel_path, "exists": present,
            "bytes": p.stat().st_size if present else 0,
            "sha256": sha256_file(p) if present else "",
            "unchanged_after": False,
        })
    layer_checks = [inspect_layer(root, name, issue, inventory, baseline=baseline,
                                  comparison_root=comparison_root) for name in LAYERS]
    comparison_manifest_path = safe_path(
        comparison_root, "data/processed/depmap_model_multiomics_manifest.json")
    comparison_manifest: dict[str, Any] = {}
    if comparison_manifest_path.is_file():
        try:
            comparison_manifest = read_json(comparison_manifest_path)
        except (ValueError, OSError, json.JSONDecodeError):
            issue("REVIEW", "comparison_provenance", "COMPARISON_MANIFEST_INVALID",
                  "Comparison-root DepMap multiomics manifest could not be read")
    elif comparison_root != root:
        issue("REVIEW", "comparison_provenance", "COMPARISON_MANIFEST_MISSING",
              "Comparison-root DepMap multiomics manifest does not exist")

    inventory_comparison: dict[str, Any] = {"status": "NOT_CHECKED"}
    if comparison_root != root:
        # Even matching data across two worktrees is not official release verification.
        issue("REVIEW", "comparison_provenance", "CROSS_ROOT_LINEAGE_REVIEW_REQUIRED",
              "Reference matrices reside in a different worktree; their release lineage needs manual confirmation")
        layer_manifest = comparison_manifest.get("layers", {})
        for layer_info in layer_checks:
            name = layer_info["layer"]
            declared = layer_manifest.get(name, {}) if isinstance(layer_manifest, dict) else {}
            if declared and layer_info["source_parquet_status"].startswith("ORDER_AND_SAMPLE_MATCH"):
                expected_models = layer_info["models_n"]
                expected_genes = layer_info["genes_n"]
                if (declared.get("available") is not True or
                    declared.get("models_n") != expected_models or
                    declared.get("genes_n") != expected_genes):
                    issue("BLOCKER", name, "COMPARISON_MANIFEST_AXES_MISMATCH",
                          "Comparison manifest model/gene counts differ from verified runtime array")
            elif not declared:
                issue("REVIEW", name, "COMPARISON_LAYER_PROVENANCE_MISSING",
                      "No per-layer evidence in comparison-root manifest")
        scientific_inventory = safe_path(root, "data/processed/depmap_input_inventory.tsv")
        reference_inventory = safe_path(comparison_root, "data/processed/depmap_input_inventory.tsv")
        if scientific_inventory.is_file() and reference_inventory.is_file():
            try:
                inventory_comparison = compare_inventories(
                    scientific_inventory, reference_inventory)
                if not inventory_comparison["semantic_match"]:
                    issue("REVIEW", "comparison_provenance", "INPUT_INVENTORY_CONTENT_DIFFERENT",
                          "DepMap inventory source identities (release/file/size/SHA) differ")
                if inventory_comparison["scientific_duplicate_rows"] or inventory_comparison["reference_duplicate_rows"]:
                    issue("REVIEW", "comparison_provenance", "INPUT_INVENTORY_DUPLICATE_ENTRIES",
                          "Repeated identical source-identity rows in a DepMap inventory")
            except (OSError, ValueError, csv.Error) as ex:
                issue("REVIEW", "comparison_provenance", "INPUT_INVENTORY_UNREADABLE",
                      f"Unable to verify source identities: {type(ex).__name__}")
        else:
            issue("REVIEW", "comparison_provenance", "INPUT_INVENTORY_MISSING",
                  "Cannot cross-check both DepMap input inventories")

    for row in comparison_files:
        p = safe_path(comparison_root, row["path"])
        current = sha256_file(p) if p.is_file() else ""
        row["unchanged_after"] = row["exists"] and current == row["sha256"]
        if row["exists"] and not row["unchanged_after"]:
            issue("BLOCKER", row["path"], "REFERENCE_CHANGED_DURING_AUDIT",
                  "The comparison Parquet or manifest was modified during audit")
    # Explicit release provenance, but no unaudited inference about downloaded raw source.
    release = multiomics.get("depmap_release") if multiomics else None
    provenance = {"index_contract": manifest.get("index_contract"), "index_schema_version": manifest.get("schema_version"),
                  "index_builder": manifest.get("builder"), "index_built_at": manifest.get("built_at"),
                  "reference_source_description": manifest.get("source"),
                  "depmap_release_declared": release,
                  "depmap_release_comparison_declared": comparison_manifest.get("depmap_release"),
                  "depmap_release_declared_source": (
                      "scientific_root_manifest" if release else
                      "comparison_root_manifest" if comparison_manifest.get("depmap_release") else "UNAVAILABLE"),
                  "dependency_call_method": depmeta.get("dependency_call_method"),
                  "external_source_authenticity": "NOT_INDEPENDENTLY_VERIFIED",
                  "notes": "A declared release or local checksum does not establish independent official source authenticity"}
    if not release and not comparison_manifest.get("depmap_release"):
        issue("REVIEW", "provenance", "DEPMAP_RELEASE_UNDOCUMENTED",
              "DepMap release missing in scientific and comparison manifests")
    if release and comparison_manifest.get("depmap_release") and (
        release != comparison_manifest["depmap_release"]
    ):
        issue("BLOCKER", "provenance", "DEPMAP_RELEASE_CONTRADICTION",
              "Scientific and reference-root manifests declare different DepMap releases")
    if not manifest.get("source"):
        issue("REVIEW", "provenance", "INDEX_SOURCE_DESCRIPTION_MISSING", "Index source description missing")
    # Recheck the exact seven source files after source/matrix comparisons.
    for row in inventory:
        if not row["exists"]: continue
        rel = row["path"]
        current = sha256_file(safe_path(root, rel))
        if current != row["sha256"]:
            issue("BLOCKER", rel, "SOURCE_CHANGED_DURING_AUDIT", "Input changed between checks")
        row["unchanged_after"] = current == row["sha256"]
    for row in sidecar_rows:
        p = safe_path(root, row["path"])
        current = sha256_file(p) if p.is_file() else ""
        row["unchanged_after"] = current == row["sha256"] and bool(current)
        if row["exists"] and not row["unchanged_after"]:
            issue("BLOCKER", row["path"], "AXIS_METADATA_CHANGED_DURING_AUDIT",
                  "Manifest or ordered axis JSON changed during validation")
    status = "BLOCKED" if any(x["severity"] == "BLOCKER" for x in issues) else "REVIEW_REQUIRED" if any(x["severity"] == "REVIEW" for x in issues) or provenance["external_source_authenticity"] != "VERIFIED" else "PASS_TECHNICAL"
    report = {"contract": CONTRACT, "utc": datetime.now(timezone.utc).isoformat(),
              "source_root": str(root), "seven_sources_found": sum(bool(x["exists"]) for x in inventory),
              "source_files_total": 7, "source_hashes_unchanged": sum(bool(x.get("unchanged_after")) for x in inventory),
              "baseline_entries_compared": sum(x.get("baseline_match") != "NOT_AVAILABLE" for x in inventory),
              "baseline_matches": sum(x.get("baseline_match") == "MATCH" for x in inventory),
              "tables": tables, "layers": layer_checks,
              "comparison_root": str(comparison_root),
              "comparison_cross_worktree": comparison_root != root,
              "comparison_sources": comparison_files,
              "inventory_comparison": inventory_comparison,
              "comparison_source_files_found": sum(int(x["exists"]) for x in comparison_files),
              "comparison_source_files_unchanged": sum(int(x["unchanged_after"]) for x in comparison_files),
              "sidecar_files": len(sidecar_rows),
              "sidecar_present": sum(bool(x["exists"]) for x in sidecar_rows),
              "sidecar_missing": sum(not x["exists"] for x in sidecar_rows),
              "sidecars_unchanged": sum(bool(x["unchanged_after"]) for x in sidecar_rows),
              "provenance": provenance,
              "counts_by_severity": dict(Counter(x["severity"] for x in issues)),
              "technical_status": status, "ac02_scientific_acceptance": "NOT_APPROVED",
              "merge_permitted": False, "no_pyz_activity_inferred": True,
              "no_admet_integration": True}
    return report, inventory, issues, layer_checks, sidecar_rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--baseline-hashes", type=Path, default=None)
    parser.add_argument("--comparison-root", type=Path, default=None,
                        help="Independent existing repository containing processed DepMap Parquet; read-only")
    args = parser.parse_args(argv)
    root, out = args.root.resolve(), args.out.resolve()
    if out.exists():
        raise FileExistsError("Refusing to overwrite an earlier AC-02 report")
    if out.is_relative_to(root) or root.is_relative_to(out):
        raise ValueError("Output must be isolated outside source root")
    comparison_root = args.comparison_root.resolve() if args.comparison_root is not None else root
    if out.is_relative_to(comparison_root) or comparison_root.is_relative_to(out):
        raise ValueError("Output must be isolated outside comparison root")
    report, inventory, issues, layer_checks, sidecar_rows = audit(
        root, args.baseline_hashes, comparison_root=comparison_root)
    out.mkdir(parents=True)
    (out / "AC02_REPORT.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_tsv(out / "AC02_SOURCES.tsv", ["path", "exists", "bytes", "sha256", "baseline_match", "unchanged_after"], inventory)
    write_tsv(out / "AC02_FINDINGS.tsv", ["severity", "component", "code", "detail"], issues)
    write_tsv(out / "AC02_SIDECARS.tsv", ["path", "exists", "sha256", "unchanged_after"], sidecar_rows)
    (out / "AC02_INVENTORY_COMPARISON.json").write_text(
        json.dumps(report["inventory_comparison"], indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8")
    write_tsv(out / "AC02_COMPARISON_SOURCES.tsv",
              ["path", "exists", "bytes", "sha256", "unchanged_after"],
              report["comparison_sources"])
    write_tsv(out / "AC02_LAYERS.tsv", ["layer", "axis_status", "models_n", "genes_n", "shape", "dtype", "orientation", "sample_compared", "sample_matched", "source_parquet_status"], layer_checks)
    lines = ["# MCL Э0 / AC-02 — аудит генетических источников", "",
             f"**Статус:** {report['technical_status']} (научная приёмка НЕ УТВЕРЖДЕНА)", "",
             f"- Входных источников на месте: {report['seven_sources_found']}/7",
             f"- SHA-256 повторно совпали: {report['source_hashes_unchanged']}/7",
             f"- Сопроводительные файлы обнаружены: {report['sidecar_present']}/{report['sidecar_files']}",
             f"- Из найденных сопроводительных файлов неизменны: {report['sidecars_unchanged']}/{report['sidecar_present']}",
             f"- Совпадений с исторической базой Э0: {report['baseline_matches']}/{report['baseline_entries_compared']}",
             f"- Использована отдельная копия для сравнения: {'ДА' if report['comparison_cross_worktree'] else 'НЕТ'}",
             f"- Источники сравнения целы: {report['comparison_source_files_unchanged']}/{report['comparison_source_files_found']}",
             f"- DepMap release в научной копии: {report['provenance']['depmap_release_declared'] or 'НЕ УКАЗАН'}",
             f"- DepMap release в копии сравнения: {report['provenance']['depmap_release_comparison_declared'] or 'НЕ УКАЗАН'}",
             f"- Подтверждение официального первоисточника: {report['provenance']['external_source_authenticity']}",
             "", "## Массивы gene × model", "",
             "| Слой | Генов | Моделей | Оси | Источник/значения |", "|---|---:|---:|---|---|"]
    for layer in layer_checks:
        lines.append(f"| {layer['layer']} | {layer['genes_n']} | {layer['models_n']} | {layer['axis_status']} | {layer['source_parquet_status']} ({layer['sample_matched']}/{layer['sample_compared']}) |")
    lines += ["", "## Зафиксированные вопросы", ""]
    lines += [f"- **{x['severity']} / {x['component']} / {x['code']}**: {x['detail']}" for x in issues] or ["- Автоматических замечаний нет, но научная приёмка всё равно ручная."]
    lines += ["", "## Ограничение выводов", "",
              "Верификация идентификаторов и массивов не подтверждает биологическую валидность DepMap/UniProt и не означает, что PYZ связывается с белком.",
              "Источник и версия должны быть документально проверены по официальному релизу до научной приёмки.",
              "Э0 не принимается автоматически; ни один PR не сливается; ADMET не обновляется.", ""]
    (out / "AC02_REVIEW_RU.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"status": report["technical_status"], "sources": report["seven_sources_found"],
                      "blockers": report["counts_by_severity"].get("BLOCKER", 0),
                      "reviews": report["counts_by_severity"].get("REVIEW", 0), "merge_permitted": False}, ensure_ascii=False))
    return 2 if report["technical_status"] == "BLOCKED" else 0


if __name__ == "__main__":
    raise SystemExit(main())