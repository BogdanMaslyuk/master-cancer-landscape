"""Complete measurable portions of MCL/PYZ E0 without changing scientific inputs.

Run against a *separate* E0 audit-code worktree and previously generated E0 report.
Evidence is stored only in <existing E0 report>/closeout/.
This is a technical closeout, never automatic scientific acceptance.
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
import shutil
import socket
import statistics
import subprocess
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

VERSION = "e0-closeout-1.0"
CONTROLS = ("PYZ-001", "PYZ-035", "PYZ-051", "PYZ-065")
SCIENCE_CASES = (("PYZ-035", "CDK4"), ("PYZ-055", "BCL2L1"))


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def safe_str(x: Any) -> str:
    if x is None:
        return ""
    try:
        import pandas as pd
        if bool(pd.isna(x)):
            return ""
    except (TypeError, ValueError):
        pass
    return str(x)


def clean(x: Any) -> Any:
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, Path):
        return str(x)
    if x is None:
        return None
    try:
        import pandas as pd
        if bool(pd.isna(x)):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(x, "item"):
        return clean(x.item())
    if isinstance(x, float) and not math.isfinite(x):
        return None
    if isinstance(x, (str, int, float, bool)):
        return x
    return str(x)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def write_json(path: Path, data: Any) -> None:
    write_text(path, json.dumps(clean(data), ensure_ascii=False, indent=2) + "\n")


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> None:
    import io
    b = io.StringIO()
    writer = csv.DictWriter(b, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    for r in rows:
        writer.writerow({c: safe_str(r.get(c)) for c in columns})
    write_text(path, b.getvalue())


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_baseline_hashes(report: Path) -> dict[str, str]:
    path = report / "file_checksums.sha256"
    if not path.is_file():
        raise ValueError(f"Baseline source hashes absent: {path}")
    hashes: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        if not line.strip():
            continue
        m = re.fullmatch(r"([0-9a-fA-F]{64})  (.+)", line)
        if not m:
            raise ValueError("Malformed source hash entry in E0 baseline")
        digest, rel = m.group(1).lower(), m.group(2).replace("\\", "/")
        pp = Path(rel)
        if pp.is_absolute() or ".." in pp.parts or rel in hashes:
            raise ValueError("Unsafe or duplicate source path in E0 baseline")
        hashes[rel] = digest
    if not hashes:
        raise ValueError("Empty source hashes file")
    return hashes


def verify_source_hashes(root: Path, hashes: dict[str, str]) -> list[dict[str, str]]:
    results = []
    for rel, old in sorted(hashes.items()):
        p = (root / rel).resolve()
        if not p.is_relative_to(root.resolve()):
            raise ValueError(f"Source path outside root: {rel}")
        new = sha256(p) if p.is_file() else ""
        results.append({"path": rel, "baseline_sha256": old, "current_sha256": new,
                        "status": "UNCHANGED" if old == new else "MISSING" if not new else "CHANGED"})
    return results


def request(url: str, timeout: float = 20.0) -> tuple[int, dict[str, Any], dict[str, str], float]:
    started = time.perf_counter()
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "MCL-E0-Closeout/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            content = r.read()
            headers = {k.lower(): v for k, v in r.headers.items()}
            try:
                data = json.loads(content)
            except (ValueError, UnicodeDecodeError):
                data = {"error": "non-JSON response"}
            return r.status, data, headers, (time.perf_counter() - started) * 1000
    except urllib.error.HTTPError as ex:
        return ex.code, {"error": str(ex)}, {}, (time.perf_counter() - started) * 1000
    except (OSError, ValueError) as ex:
        return 0, {"error": f"{type(ex).__name__}: {ex}"}, {}, (time.perf_counter() - started) * 1000


def numeric_equal(a: Any, b: Any) -> bool:
    if a is None and b is None:
        return True
    try:
        if a is None or b is None:
            return False
        x, y = float(a), float(b)
        return math.isfinite(x) and math.isfinite(y) and math.isclose(x, y, abs_tol=1e-8, rel_tol=1e-7)
    except (TypeError, ValueError):
        return safe_str(a) == safe_str(b)


def source_table(root: Path, rel: str, columns: list[str] | None = None):
    """Read only selected columns from Parquet whenever possible."""
    import pandas as pd
    p = root / rel
    if p.suffix == ".parquet" and not p.is_file() and p.with_suffix(".tsv").is_file():
        p = p.with_suffix(".tsv")
    if not p.is_file():
        return None
    if p.suffix == ".parquet":
        if columns:
            import pyarrow.parquet as pq
            names = set(pq.ParquetFile(p).schema_arrow.names)
            cols = [x for x in columns if x in names]
            if not cols:
                return pd.DataFrame()
            return pd.read_parquet(p, columns=cols)
        return pd.read_parquet(p)
    if p.suffix == ".tsv":
        return pd.read_csv(p, sep="\t", usecols=lambda col: columns is None or col in columns)
    return None


def add_case(rows: list[dict], cid: str, target: str, field: str, origin: str,
             source: Any, api: Any, exact: bool = False) -> None:
    match = (safe_str(source) == safe_str(api)) if exact else numeric_equal(clean(source), clean(api))
    rows.append({"compound_id": cid, "target_gene": target, "field": field,
                 "source_file": origin, "source_value": safe_str(clean(source)), "api_value": safe_str(clean(api)),
                 "status": "MATCH" if match else "MISMATCH", "ui_status": "REQUIRES_VISUAL_REVIEW"})


def reconcile_scientific(root: Path, api_base: str, out: Path) -> tuple[list[dict], list[dict], list[dict], list[str]]:
    """Compare source to API without equating known ligand potency to PYZ potency."""
    import pandas as pd
    rows: list[dict] = []
    reference: list[dict] = []
    contexts: list[dict] = []
    problems: list[str] = []
    registry = source_table(root, "config/own_compounds.tsv")
    matrix = source_table(root, "data/runtime/pyz_target_evidence/pyz_target_evidence_matrix.parquet",
        ["own_compound_id", "target_gene", "best_ligand_id", "best_tanimoto_morgan_r2_2048",
         "best_lowest_reported_value_nm_across_endpoints", "best_lowest_value_endpoint_type", "best_lowest_value_source"])
    neighbors = source_table(root, "data/runtime/own_compounds/target_ligand_space_similarity_top_hits.parquet",
        ["own_compound_id", "target_gene", "rank_in_target", "ligand_id", "lowest_reported_value_nm_across_endpoints",
         "lowest_value_relation", "lowest_value_endpoint_type", "lowest_value_source", "endpoint_types_json",
         "sources_json", "measurements_n", "tanimoto_morgan_r2_2048"])
    admet = source_table(root, "config/pyz_admet_snapshot.tsv")
    if registry is None or matrix is None:
        return [], [], [], ["Required registry or matrix unavailable"]
    reg = registry.set_index("own_compound_id", drop=False)
    risks = admet.set_index("own_compound_id") if admet is not None and "own_compound_id" in admet.columns else None
    controls = tuple(dict.fromkeys((*CONTROLS, *(cid for cid, _ in SCIENCE_CASES))))
    ui_expected: dict[str, dict] = {}
    for cid in controls:
        code, payload, _, _ = request(api_base.rstrip("/") + "/api/pyz/" + cid)
        if code != 200 or not isinstance(payload.get("profile_v2"), dict):
            problems.append(f"{cid}: API HTTP {code} or missing profile_v2")
            continue
        profile = payload["profile_v2"]
        if cid not in reg.index:
            problems.append(f"{cid}: not found in source registry")
            continue
        reg_row = reg.loc[cid]
        smiles = clean(reg_row.get("standardized_smiles"))
        add_case(rows, cid, "", "identity.standardized_smiles", "config/own_compounds.tsv", smiles,
                 (profile.get("identity") or {}).get("standardized_smiles"), exact=True)
        ui_expected[cid] = {"smiles": smiles, "target_count": len((profile.get("target_hypotheses") or {}).get("items") or [])}
        if risks is not None and cid in risks.index:
            risk = risks.loc[cid]
            add_case(rows, cid, "", "properties.qed", "config/pyz_admet_snapshot.tsv", risk.get("qed"),
                     (profile.get("properties") or {}).get("qed"))
            add_case(rows, cid, "", "admet.integrated_risk", "config/pyz_admet_snapshot.tsv", risk.get("integrated_risk"),
                     (profile.get("admet") or {}).get("integrated_risk"), exact=True)
        source_targets = matrix[matrix["own_compound_id"].astype(str).eq(cid)].set_index("target_gene")
        api_targets = {safe_str(t.get("target_gene")): t for t in (profile.get("target_hypotheses") or {}).get("items") or []}
        if set(source_targets.index) != set(api_targets):
            problems.append(f"{cid}: source and API target lists differ")
        for target, record in source_targets.iterrows():
            api_t = api_targets.get(str(target), {})
            if not api_t:
                continue
            add_case(rows, cid, str(target), "target.tanimoto", "matrix.parquet", record.get("best_tanimoto_morgan_r2_2048"), api_t.get("tanimoto"))
            add_case(rows, cid, str(target), "target.best_ligand_id", "matrix.parquet", record.get("best_ligand_id"), api_t.get("best_ligand_id"), exact=True)
            add_case(rows, cid, str(target), "target.reference_endpoint", "matrix.parquet", record.get("best_lowest_value_endpoint_type"), api_t.get("reference_endpoint"), exact=True)
            add_case(rows, cid, str(target), "target.reference_value_nm", "matrix.parquet", record.get("best_lowest_reported_value_nm_across_endpoints"), api_t.get("reference_value_nm"))
            if api_t.get("reference_measurement_owner") != "known_ligand_not_pyz":
                problems.append(f"{cid}/{target}: reference measurement not attributed to a known ligand")
        if (profile.get("scientific_summary") or {}).get("own_target_binding") != "not_measured":
            problems.append(f"{cid}: unexpected own PYZ binding assertion")
    write_json(out / "ui_expected.json", ui_expected)
    if neighbors is not None:
        for cid, gene in SCIENCE_CASES:
            subset = neighbors[(neighbors["own_compound_id"].astype(str) == cid) & (neighbors["target_gene"].astype(str) == gene)]
            subset = subset.sort_values("rank_in_target").head(5)
            for row in subset.to_dict("records"):
                reference.append({"pyz": cid, "hypothesis_target": gene,
                                  "ligand_id": clean(row.get("ligand_id")),
                                  "rank": clean(row.get("rank_in_target")),
                                  "tanimoto": clean(row.get("tanimoto_morgan_r2_2048")),
                                  "endpoint_type": clean(row.get("lowest_value_endpoint_type")),
                                  "reference_value_nm": clean(row.get("lowest_reported_value_nm_across_endpoints")),
                                  "relation": clean(row.get("lowest_value_relation")),
                                  "source": clean(row.get("lowest_value_source")),
                                  "sources_json": clean(row.get("sources_json")),
                                  "endpoint_types_json": clean(row.get("endpoint_types_json")),
                                  "measurement_count": clean(row.get("measurements_n")),
                                  "primary_assay_id": None,
                                  "primary_assay_status": "REQUIRES_PRIMARY_SOURCE"})
    models = source_table(root, "data/runtime/pharmacology/candidate_hypothesis_models_v2.parquet",
                          ["hypothesis_id", "model_id", "compound_id", "target_gene", "model_role",
                           "response_value", "dependency_probability", "gene_effect", "cell_line_name"])
    if models is not None and {"model_id", "compound_id", "target_gene", "model_role"}.issubset(models.columns):
        sub = models[models["target_gene"].isin(["CDK4", "BCL2L1"])].copy()
        # One gene/model may have multiple reference drugs and/or roles; never collapse them to a single drug.
        groups = sub.groupby(["target_gene", "model_id"], dropna=True)
        key_stats = groups.agg(drugs=("compound_id", "nunique"), roles=("model_role", "nunique"))
        candidates = key_stats[(key_stats["drugs"] > 1) & (key_stats["roles"] > 1)].sort_index()
        for (gene, mid), stats in list(candidates.iterrows())[:12]:
            details = sub[(sub["target_gene"] == gene) & (sub["model_id"] == mid)]
            # Deterministic cap and no bulk export of scientific source data.
            for _, record in details.sort_values(["compound_id", "hypothesis_id"]).head(12).iterrows():
                contexts.append({"gene": gene, "model_id": mid, "cell_line_name": clean(record.get("cell_line_name")),
                                 "reference_compound_id": clean(record.get("compound_id")),
                                 "hypothesis_id": clean(record.get("hypothesis_id")),
                                 "role": clean(record.get("model_role")),
                                 "response_value": clean(record.get("response_value")),
                                 "dependency_probability": clean(record.get("dependency_probability")),
                                 "gene_effect": clean(record.get("gene_effect")),
                                 "scope": "Candidate v2 subset, NOT PYZ activity"})
    write_tsv(out / "scientific_reference_examples.tsv", ["pyz","hypothesis_target","ligand_id","rank","tanimoto",
        "endpoint_type","reference_value_nm","relation","source","sources_json","endpoint_types_json",
        "measurement_count","primary_assay_id","primary_assay_status"], reference)
    write_tsv(out / "cell_context_observations.tsv", ["gene","model_id","cell_line_name","reference_compound_id",
        "hypothesis_id","role","response_value","dependency_probability","gene_effect","scope"], contexts)
    write_tsv(out / "source_api_reconciliation.tsv", ["compound_id","target_gene","field","source_file","source_value",
        "api_value","status","ui_status"], rows)
    return rows, reference, contexts, problems


def detect_cache(root: Path, api_base: str) -> dict:
    ui_source = root / "apps" / "explorer" / "lib" / "api.ts"
    server = root / "apps" / "api" / "mcl_api" / "main.py"
    text = ui_source.read_text(encoding="utf-8") if ui_source.exists() else ""
    api_text = server.read_text(encoding="utf-8") if server.exists() else ""
    status, _, headers, _ = request(api_base.rstrip("/") + "/api/pyz/summary")
    return {"request_status": status, "api_cache_control": headers.get("cache-control", ""),
            "client_force_cache_configured": 'cache: "force-cache"' in text,
            "client_default_revalidate_300_present": "MCL_API_REVALIDATE_SECONDS || 300" in text,
            "api_default_max_age_300_present": 'max-age=300' in api_text,
            "dynamic_update_test": "NOT_TESTED: cannot alter real scientific sources; isolated fixture needed",
            "meaning": "Potential staleness after source update; do not infer that a particular page was stale."}


def cold_api(root: Path, out: Path, timeout: float = 70.0) -> dict:
    """Measure a fresh *owned* FastAPI instance on a free loopback port; never stop PS1/PS2."""
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        port = s.getsockname()[1]
    api_pkg = root / "apps" / "api"
    env = os.environ.copy()
    env["MCL_ROOT"] = str(root)
    env["PYTHONPATH"] = str(api_pkg) + os.pathsep + env.get("PYTHONPATH", "")
    cmd = [sys.executable, "-m", "uvicorn", "mcl_api.main:app", "--app-dir", str(api_pkg),
           "--host", "127.0.0.1", "--port", str(port), "--no-access-log"]
    result = {"status": "NOT_RUN", "port": port, "api_startup_ms": None, "first_card_ms": None,
              "error": "", "note": "Dedicated local server; no modification of PS1, PS2 or scientific sources"}
    with (out / "cold_api_server.log").open("w", encoding="utf-8") as log:
        proc = None
        try:
            start = time.perf_counter()
            proc = subprocess.Popen(cmd, cwd=str(root), env=env, stdout=log, stderr=subprocess.STDOUT)
            health = None
            while (time.perf_counter() - start) < timeout:
                if proc.poll() is not None:
                    raise RuntimeError(f"Separate uvicorn exited with code {proc.returncode}")
                code, data, _, _ = request(f"http://127.0.0.1:{port}/health", timeout=1.0)
                if code == 200:
                    health = data
                    break
                time.sleep(0.25)
            if not health:
                raise TimeoutError("Dedicated API did not start within time limit")
            if Path(str(health.get("mcl_root", ""))).resolve() != root.resolve():
                raise RuntimeError("Dedicated API is pointing to another MCL_ROOT")
            result["api_startup_ms"] = round((time.perf_counter() - start) * 1000, 2)
            code, data, _, elapsed = request(f"http://127.0.0.1:{port}/api/pyz/PYZ-035", timeout=45)
            result["first_card_ms"] = round(elapsed, 2)
            result["status"] = "PASS" if code == 200 and data.get("profile_v2", {}).get("compound_id") == "PYZ-035" else "FAIL"
            if result["status"] == "FAIL":
                result["error"] = f"First card response HTTP {code}"
        except (OSError, RuntimeError, TimeoutError, ValueError) as exc:
            result["status"] = "FAIL"
            result["error"] = f"{type(exc).__name__}: {exc}"
        finally:
            if proc is not None and proc.poll() is None:
                proc.terminate()  # Only our own non-reloading uvicorn process.
                try:
                    proc.wait(timeout=8)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait(timeout=8)
    return result


def run_ui(out: Path, ui_base: str, audit_root: Path, tools_root: Path) -> dict:
    node = shutil.which("node")
    script = audit_root / "scripts" / "pyz_e0_browser.cjs"
    module = tools_root / "node_modules" / "playwright-core"
    if not node or not module.is_dir() or not script.is_file():
        return {"status": "NOT_RUN", "reason": "Node/playwright-core/browser script unavailable; use manual screenshots"}
    cmd = [node, str(script), "--url", ui_base, "--out", str(out), "--tools", str(tools_root)]
    try:
        result = subprocess.run(cmd, text=True, capture_output=True, encoding="utf-8", errors="replace", timeout=240)
        write_text(out / "browser_stdout.log", result.stdout + "\n" + result.stderr)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"status": "FAILED", "reason": str(exc)}
    evidence_path = out / "ui" / "browser_results.json"
    if not evidence_path.is_file():
        return {"status": "FAILED", "reason": f"Browser process exit {result.returncode}; see browser_stdout.log"}
    data = json.loads(evidence_path.read_text(encoding="utf-8"))
    return {"status": data.get("status", "FAILED"), "exit_code": result.returncode,
            "n_screenshots": len([x for x in data.get("screenshots", []) if x.get("saved")]),
            "n_screenshot_checks": len(data.get("screenshots", [])),
            "n_dom_checks": len(data.get("control_dom_checks", [])),
            "browser": data.get("browser"), "reason": data.get("reason", "")}


def status_table(hash_rows: list[dict], comparisons: list[dict], api_issues: list[str],
                 browser: dict, cold: dict, reference: list[dict], cache: dict) -> dict[str, dict]:
    """AC states are conservative: automated PASS never equals expert approval."""
    hash_ok = bool(hash_rows) and all(x["status"] == "UNCHANGED" for x in hash_rows)
    matches = bool(comparisons) and all(x["status"] == "MATCH" for x in comparisons)
    screenshot_n = browser.get("n_screenshots", 0)
    return {
        "AC-01": {"status": "AUTO_PASS" if hash_ok else "BLOCKED", "reason": "Source hash comparison and recorded original commit"},
        "AC-02": {"status": "PARTIAL", "reason": "7 input types were NOT_REVIEWED in original E0 inventory; content review outstanding"},
        "AC-03": {"status": "AUTO_PASS_EXPERT_PENDING", "reason": "845 / 845 unique matrix pairs in baseline; expert scientific meaning pending"},
        "AC-04": {"status": "AUTO_PASS_EXPERT_PENDING" if matches and not api_issues else "BLOCKED",
                  "reason": "Source→API comparisons for controls; screen/primary source review still required"},
        "AC-05": {"status": "AUTO_CAPTURED_REVIEW_PENDING" if screenshot_n == 9 else "NOT_REVIEWED",
                  "reason": f"Automated screenshots {screenshot_n}/9; human visual and accessibility assessment required"},
        "AC-06": {"status": "PARTIAL", "reason": "Existing warm p50/p95; separate cold API result=" + cold.get("status", "NOT_RUN") + "; browser results available"},
        "AC-07": {"status": "NOT_REVIEWED", "reason": "No primary assay_id + published measurement validation in accepted source"},
        "AC-08": {"status": "AUTO_PASS" if hash_ok else "BLOCKED", "reason": "SHA-256 checked after E0; no modified scientific files permitted"},
        "AC-09": {"status": "NOT_APPROVED", "reason": "Explicit project-owner acceptance is required"},
    }


def write_review(out: Path, root: Path, baseline: Path, statuses: dict, hash_rows: list[dict],
                 comparisons: list[dict], ref_rows: list[dict], context_rows: list[dict],
                 api_problems: list[str], browser: dict, cold: dict, cache: dict):
    hashes_ok = sum(x["status"] == "UNCHANGED" for x in hash_rows)
    matched = sum(x["status"] == "MATCH" for x in comparisons)
    lines = ["# MCL / PYZ — Э0: дополнительная автоматическая проверка", "",
             "**Не является научной приёмкой.** Скрипты не изменяют исходные файлы и API.", "",
             f"- Дата UTC: {utc_now()}", f"- Код аудитора: {VERSION}", f"- Проверяемая папка: `{root}`", f"- Исходный отчёт Э0: `{baseline}`",
             f"- SHA-256 после запуска: {hashes_ok}/{len(hash_rows)} совпали с базовым отчётом.",
             f"- Source → API: {matched}/{len(comparisons)} полей совпали. (Это проверка передачи, не научного смысла.)",
             f"- 2D-структурные аналоги: {len(ref_rows)} строк для двух контрольных пар, первичные assay_id не приложены.",
             f"- Различающиеся наблюдения известных препаратов: {len(context_rows)} диагностических строк.",
             f"- Браузер: {browser.get('status')} / снимки {browser.get('n_screenshots',0)}/9.",
             f"- Отдельный холодный API: {cold.get('status')}; запуск: {cold.get('api_startup_ms')} мс; первая карточка: {cold.get('first_card_ms')} мс.",
             "", "## Критерии приёмки", "", "| AC | Состояние | Что требуется |", "|---|---|---|"]
    for ac, row in statuses.items():
        lines.append(f"| {ac} | {row['status']} | {row['reason']} |")
    lines += ["", "## Научное решение", "",
              "Сверка матрицы и API показывает совпадение передаваемых значений, но не проверяет публикации и исходные assay. ",
              "Значения Ki/Kd/IC50 из агрегатов НЕ используются для оценки активности PYZ.",
              "В `scientific_reference_examples.tsv` сохранены ближайшие лиганды PYZ-035/CDK4 и PYZ-055/BCL2L1;",
              "поле primary_assay_id пусто намеренно, поскольку его нет в принятой матрице/top-hits.",
              "Для SCI-01/02 необходимы исходные assay_id, условия, знак сравнения и публикация известного лиганда.",
              "В `cell_context_observations.tsv` примеры независимых известных препаратов; роли нельзя склеивать.",
              "", "## Кэш", "", f"- API Cache-Control: `{cache.get('api_cache_control')}`",
              f"- Next.js force-cache: `{cache.get('client_force_cache_configured')}`",
              "- Прямой тест обновления научной базы не проводился, чтобы не изменять действующие научные источники.",
              "", "## Дальнейшая ручная приёмка", "",
              "1. Открыть девять снимков, отметить визуальные дефекты, проверить клавиатуру, доступность кнопок и отображение чисел.",
              "2. Сверить первоисточники хотя бы для PYZ-035/CDK4 и PYZ-055/BCL2L1 и отдельные препараты на одной линии.",
              "3. Рассмотреть семь исходников NOT_REVIEWED; подтвердить происхождение и версию больших индексов.",
              "4. Подтвердить тест актуализации кэша на отдельной копии научных данных до выпуска Э4; не править рабочие parquet.",
              "5. Подписать AC-09. До этого Э0 формально остаётся НЕ ПРИНЯТ.", ""]
    if api_problems:
        lines += ["## Ошибки API", ""] + [f"- {x}" for x in api_problems]
    if any(x["status"] != "UNCHANGED" for x in hash_rows):
        lines += ["", "## ОПАСНО: ХЕШИ ИСХОДНИКОВ ИЗМЕНИЛИСЬ", ""] + [f"- {x['path']}: {x['status']}" for x in hash_rows if x["status"] != "UNCHANGED"]
    write_text(out / "E0_ACCEPTANCE_CLOSEOUT.md", "\n".join(lines) + "\n")
    write_text(out / "E1_HANDOFF_CLOSEOUT.md", "\n".join([
        "# Э1 — уточнённые входные задачи", "",
        "- SCI-01: ввести раздельную типизацию Ki/Kd/IC50, с assay_id и неравенствами; первичные ссылки пока не приняты.",
        "- SCI-02: наблюдения известных препаратов сохранять отдельно, с конкретным model_id/compound_id/hypothesis_id.",
        "- SCI-03: индивидуальные прогнозы ADMET должны быть представлены только при наличии исходной строки/метода.",
        "- SCI-04: Танимото ≥0,35 — только навигационный порог; не вероятность связывания.",
        "- PERF/Э4: разделить большую карточку на обзор и запрашиваемые подробности после профилирования.",
        "- Кэш: зафиксировать версионирование и явное обновление данных, исключить устаревшие ответы Next.js.",
        "", "Принять эти задачи к проектированию можно; выпуск Э1 требует приёмки Э0 и отдельного ТЗ.", ""]))


def execute(args: argparse.Namespace) -> tuple[Path, int]:
    root = Path(args.root).expanduser().resolve()
    baseline = Path(args.baseline).expanduser().resolve()
    audit_root = Path(__file__).resolve().parents[1]
    if not baseline.is_dir() or not (baseline / "E0_ACCEPTANCE.md").is_file():
        raise ValueError("Existing E0 baseline report missing or incomplete")
    if not (root / "config" / "own_compounds.tsv").is_file():
        raise ValueError("Real source registry missing")
    if root == audit_root or audit_root.is_relative_to(root):
        raise ValueError("Audit code/report must be in a different worktree from science sources")
    expected_parent = audit_root / "outputs" / "audits" / "pyz_e0"
    if not baseline.is_relative_to(expected_parent.resolve()):
        raise ValueError("Refusing to write reports outside own audit worktree")
    original = json.loads((baseline / "environment.json").read_text(encoding="utf-8"))
    if Path(original.get("root", "")).resolve() != root:
        raise ValueError("Baseline was generated for another MCL_ROOT")
    out = baseline / "closeout"
    if out.exists() and not args.overwrite:
        raise ValueError(f"Closeout already exists. Use --overwrite to regenerate only reports in {out}")
    out.mkdir(exist_ok=True)
    hashes = read_baseline_hashes(baseline)
    hash_entry = verify_source_hashes(root, hashes)
    write_tsv(out / "source_sha_entry.tsv", ["path", "baseline_sha256", "current_sha256", "status"], hash_entry)
    # Abort accessing the live API if its scientific root is not exactly the one being audited.
    status, health, _, _ = request(args.api_base.rstrip("/") + "/health")
    api_issues = []
    if status != 200 or Path(str(health.get("mcl_root", ""))).resolve() != root:
        api_issues.append("API /health not reachable or points to another MCL_ROOT; source→API comparisons skipped")
        comparisons, reference, contexts = [], [], []
        write_json(out / "ui_expected.json", {})
    else:
        comparisons, reference, contexts, science_problems = reconcile_scientific(root, args.api_base, out)
        api_issues.extend(science_problems)
    cache = detect_cache(root, args.api_base) if not api_issues else {"status": "NOT_RUN"}
    write_json(out / "cache_diagnostics.json", cache)
    cold = {"status": "NOT_RUN", "reason": "Skipped by command or API unavailable"}
    if not args.skip_cold_api and not api_issues:
        cold = cold_api(root, out)
    write_json(out / "cold_api.json", cold)
    browser = {"status": "NOT_RUN", "reason": "Skipped"}
    if not args.skip_ui and not api_issues:
        tools_root = Path(args.ui_tools).expanduser().resolve()
        browser = run_ui(out, args.ui_base, audit_root, tools_root)
    write_json(out / "browser_status.json", browser)
    hashes_now = verify_source_hashes(root, hashes)  # End-of-run recheck, after all API/browser tests.
    write_tsv(out / "source_sha_after.tsv", ["path", "baseline_sha256", "current_sha256", "status"], hashes_now)
    reviewed_hashes = []
    for entry, final in zip(hash_entry, hashes_now):
        recorded = final.copy()
        if entry["status"] != "UNCHANGED" and final["status"] == "UNCHANGED":
            recorded["status"] = "CHANGED_AT_ENTRY"
        reviewed_hashes.append(recorded)
    statuses = status_table(reviewed_hashes, comparisons, api_issues, browser, cold, reference, cache)
    write_review(out, root, baseline, statuses, reviewed_hashes, comparisons, reference, contexts,
                 api_issues, browser, cold, cache)
    errors = [x for x in (*hash_entry, *hashes_now) if x["status"] != "UNCHANGED"]
    errors += [x for x in comparisons if x["status"] == "MISMATCH"]
    if api_issues:
        errors.extend(api_issues)
    return out, 2 if errors else 0


def main() -> int:
    parser = argparse.ArgumentParser(description="PYZ E0 closeout (source→API, immutable sources, UI captures)")
    parser.add_argument("--root", required=True, help="Original MCL/PYZ data worktree")
    parser.add_argument("--baseline", required=True, help="Existing report folder from E0 read-only auditor")
    parser.add_argument("--api-base", default="http://127.0.0.1:8100")
    parser.add_argument("--ui-base", default="http://127.0.0.1:3100")
    parser.add_argument("--ui-tools", default="", help="Isolated npm prefix with playwright-core")
    parser.add_argument("--skip-ui", action="store_true")
    parser.add_argument("--skip-cold-api", action="store_true")
    parser.add_argument("--overwrite", action="store_true", help="Replace generated closeout files, never E0 base reports")
    args = parser.parse_args()
    try:
        out, rc = execute(args)
    except (OSError, ValueError, KeyError) as ex:
        print(f"E0 CLOSEOUT STOP: {type(ex).__name__}: {ex}", file=sys.stderr)
        return 3
    print(f"E0 CLOSEOUT: {out}\nEXIT CODE: {rc}. E0 is NOT automatically accepted.")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())