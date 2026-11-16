"""Consolidate MCL PYZ E0 acceptance evidence, without accepting science or changing inputs.

Inputs are an existing E0 baseline/closeout ZIP or directory. Python stdlib only.
Outputs are an independent report, with all nine criteria explicitly requiring
human scientific review. Does not access API or start/stop servers.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import stat
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

VERSION = "mcl-e0-acceptance-gate-v1"
AC_IDS = tuple(f"AC-{i:02d}" for i in range(1, 10))
MAX_ENTRY_BYTES = 50 * 1024 * 1024
MAX_ARCHIVE_BYTES = 256 * 1024 * 1024
SHA = re.compile(r"^[0-9a-fA-F]{64}$")
SHA_LINE = re.compile(r"^([0-9a-fA-F]{64})  (.+)$")
REQUIRED_BASE = (
    "E0_ACCEPTANCE.md", "file_checksums.sha256", "input_inventory.tsv",
    "coverage_matrix.tsv", "schema_audit.json", "ui_checklist.tsv",
    "api_cases.tsv", "performance.tsv",
)
REQUIRED_CLOSEOUT = (
    "closeout/source_sha_entry.tsv", "closeout/source_sha_after.tsv",
    "closeout/source_api_reconciliation.tsv", "closeout/cold_api.json",
    "closeout/browser_status.json", "closeout/E0_ACCEPTANCE_CLOSEOUT.md",
)
STATUSES = ("AUTO_EVIDENCE", "PENDING_MANUAL", "MISSING_EVIDENCE", "BLOCKED")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def norm(path: str) -> str:
    """Reject zip traversal, symlinks, Windows drive prefixes and path ambiguity."""
    path = path.replace("\\", "/")
    pp = PurePosixPath(path)
    if (not path or path.startswith("/") or re.match(r"^[A-Za-z]:", path)
            or any(x in {"", ".", ".."} for x in path.split("/"))):
        raise ValueError("Unsafe archive path")
    return str(pp)


class Evidence:
    def __init__(self, path: Path):
        self.path = path.resolve()
        self.zip: zipfile.ZipFile | None = None
        self.paths: dict[str, object] = {}
        if path.is_file() and zipfile.is_zipfile(path):
            self.zip = zipfile.ZipFile(path)
            infos = [i for i in self.zip.infolist() if not i.is_dir()]
            if not infos:
                raise ValueError("Empty E0 evidence ZIP")
            total = 0
            for info in infos:
                p = norm(info.filename)
                mode = (info.external_attr >> 16) & 0o170000
                if mode == stat.S_IFLNK:
                    raise ValueError("Symlink is forbidden in ZIP")
                if p in self.paths:
                    raise ValueError("Duplicate archive member")
                if info.file_size > MAX_ENTRY_BYTES:
                    raise ValueError("Archive member exceeds size limit")
                total += info.file_size
                if total > MAX_ARCHIVE_BYTES:
                    raise ValueError("Archive uncompressed size exceeds limit")
                self.paths[p] = info
        elif path.is_dir():
            for p in path.rglob("*"):
                if p.is_symlink():
                    raise ValueError("Symlinks forbidden in evidence directory")
                if p.is_file():
                    relative = norm(p.relative_to(path).as_posix())
                    if relative in self.paths:
                        raise ValueError("Duplicate evidence member")
                    self.paths[relative] = p
        else:
            raise FileNotFoundError(f"E0 report ZIP/directory not found: {path}")
        prefix_candidates = [p[:-len("E0_ACCEPTANCE.md")] for p in self.paths if p.endswith("E0_ACCEPTANCE.md")]
        if len(prefix_candidates) != 1:
            raise ValueError("Expected exactly one E0_ACCEPTANCE.md at report root")
        prefix = prefix_candidates[0]
        if any(not name.startswith(prefix) for name in self.paths):
            raise ValueError("E0 archive contains mixed unrelated roots")
        self.paths = {name[len(prefix):]: value for name, value in self.paths.items()}

    def has(self, name: str) -> bool:
        return name in self.paths

    def names(self):
        return set(self.paths)

    def read(self, name: str) -> bytes:
        if name not in self.paths:
            raise FileNotFoundError(name)
        info = self.paths[name]
        if self.zip:
            with self.zip.open(info) as f:
                data = f.read(MAX_ENTRY_BYTES + 1)
        else:
            with Path(info).open("rb") as f:
                data = f.read(MAX_ENTRY_BYTES + 1)
        if len(data) > MAX_ENTRY_BYTES:
            raise ValueError("Oversized report member")
        return data

    def text(self, name: str) -> str:
        return self.read(name).decode("utf-8-sig")

    def tsv(self, name: str, required: set[str] | None = None) -> list[dict[str, str]]:
        reader = csv.DictReader(io.StringIO(self.text(name)), delimiter="\t")
        if required and not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{name}: missing columns {sorted(required - set(reader.fieldnames or []))}")
        rows = list(reader)
        if any(None in row for row in rows):
            raise ValueError(f"{name}: malformed TSV rows")
        return rows

    def json(self, name: str):
        return json.loads(self.text(name))

    def close(self):
        if self.zip:
            self.zip.close()


def baseline_hashes(source: Evidence) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for raw in source.text("file_checksums.sha256").splitlines():
        if not raw.strip():
            continue
        m = SHA_LINE.fullmatch(raw)
        if not m:
            raise ValueError("Malformed baseline source SHA-256 entry")
        rel = norm(m.group(2))
        if rel in hashes:
            raise ValueError("Duplicate baseline source path")
        hashes[rel] = m.group(1).lower()
    if not hashes:
        raise ValueError("No scientific source SHA-256 hashes in baseline")
    return hashes


def compare_recorded_hashes(source: Evidence, expected: dict[str, str]):
    details = []
    failures = []
    for phase in ("entry", "after"):
        key = f"closeout/source_sha_{phase}.tsv"
        if not source.has(key):
            continue
        rows = source.tsv(key, {"path", "baseline_sha256", "current_sha256", "status"})
        mapped = {}
        for r in rows:
            path = norm(r["path"])
            if path in mapped:
                failures.append(f"Duplicated source path in {key}: {path}")
            mapped[path] = r
        if set(mapped) != set(expected):
            failures.append(f"Source path set differs from baseline in {key}")
        for path, sha in expected.items():
            row = mapped.get(path)
            if not row:
                failures.append(f"Missing source checksum {phase}: {path}")
                continue
            before, current = row["baseline_sha256"].lower(), row["current_sha256"].lower()
            status = row["status"]
            good = before == sha and current == sha and status == "UNCHANGED"
            details.append({"phase": phase, "path": path, "status": status,
                            "sha_matches": str(good).lower()})
            if not good:
                failures.append(f"Scientific SHA mismatch {phase}: {path} [{status}]")
    return details, failures


def verify_live_hashes(root: Path, expected: dict[str, str]):
    """Independently hash current scientific files only when a trusted root is supplied."""
    root = root.resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Scientific source root not found: {root}")
    rows, failures = [], []
    for rel, known in sorted(expected.items()):
        path = root / rel
        if path.is_symlink() or not path.resolve().is_relative_to(root):
            raise ValueError(f"Unsafe scientific path: {rel}")
        observed = digest(path) if path.is_file() else ""
        same = observed == known
        rows.append({"phase": "live", "path": rel,
                     "status": "UNCHANGED" if same else "MISSING" if not observed else "CHANGED",
                     "sha_matches": str(same).lower()})
        if not same:
            failures.append(f"Live scientific file SHA mismatch: {rel}")
    return rows, failures


def add(rows, code, status, evidence, note):
    assert code in AC_IDS and status in STATUSES
    rows.append({"criterion": code, "technical_status": status,
                 "evidence": evidence, "required_manual_review": note,
                 "manual_decision": "PENDING", "reviewer": "", "review_date": ""})


def analyze(source: Evidence, science_root: Path | None = None) -> tuple[dict, list[dict], list[dict], list[dict]]:
    missing_base = [x for x in REQUIRED_BASE if not source.has(x)]
    findings: list[dict] = []
    def issue(severity, code, msg):
        findings.append({"severity": severity, "criterion": code, "detail": msg})
    for path in missing_base:
        issue("MISSING", "AC-01", f"Missing baseline evidence: {path}")
    if not source.has("file_checksums.sha256"):
        hashes = {}
    else:
        hashes = baseline_hashes(source)
    hash_details, errors = compare_recorded_hashes(source, hashes) if hashes else ([], [])
    live_count = 0
    if science_root and hashes:
        live_rows, live_errors = verify_live_hashes(science_root, hashes)
        hash_details.extend(live_rows)
        errors.extend(live_errors)
        live_count = len(live_rows)
    for x in errors:
        issue("BLOCKER", "AC-08", x)
    missing_closeout = [x for x in REQUIRED_CLOSEOUT if not source.has(x)]
    for x in missing_closeout:
        issue("MISSING", "AC-09", f"Missing closeout evidence: {x}")
    ac = []
    if errors:
        s = "BLOCKED"
    elif hashes and not any(x.startswith("closeout/source_sha_") for x in missing_closeout):
        s = "AUTO_EVIDENCE"
    else:
        s = "MISSING_EVIDENCE"
    add(ac, "AC-01", s, "file_checksums.sha256 + closeout/source_sha_entry.tsv", "Сверить версию исходников и подпись владельца")

    if source.has("schema_audit.json"):
        schema = source.json("schema_audit.json")
        for item in schema.get("issues", []):
            if isinstance(item, dict) and str(item.get("severity", "")).lower() == "blocker":
                issue("BLOCKER", "AC-02", f"Baseline schema blocker: {item.get('id', item.get('message', '?'))}")
    if source.has("input_inventory.tsv"):
        inv = source.tsv("input_inventory.tsv", {"status"})
        unresolved = sum(r["status"].strip().upper() in {"NOT_REVIEWED", "MISSING", "INVALID", "ERROR"} for r in inv)
        if unresolved:
            issue("MANUAL", "AC-02", f"Input inventory unresolved rows: {unresolved}")
        status = "BLOCKED" if any(x["criterion"] == "AC-02" and x["severity"] == "BLOCKER" for x in findings) else "PENDING_MANUAL" if unresolved else "AUTO_EVIDENCE"
    else:
        status = "MISSING_EVIDENCE"
    add(ac, "AC-02", status, "input_inventory.tsv", "Проверить происхождение и версии семи ранее отмеченных источников")

    if source.has("coverage_matrix.tsv"):
        rows = source.tsv("coverage_matrix.tsv", {"compound_id", "target_gene", "state"})
        keys = [(r["compound_id"], r["target_gene"]) for r in rows]
        invalid = [r for r in rows if r["state"] not in {"present", "missing", "duplicate", "invalid"}]
        problems = len(keys) != len(set(keys)) or bool(invalid)
        n_comp = len({r["compound_id"] for r in rows})
        n_gene = len({r["target_gene"] for r in rows})
        if len(rows) != 845 or n_comp != 65 or n_gene != 13:
            problems = True
            issue("BLOCKER", "AC-03", f"Unexpected PYZ/target matrix shape: {len(rows)} rows, {n_comp} PYZ, {n_gene} genes")
        if any(r["state"] in {"missing", "duplicate", "invalid"} for r in rows):
            issue("MANUAL", "AC-03", "Some E0 coverage entries are missing/duplicate/invalid; inspect context")
        if len(keys) != len(set(keys)) or invalid:
            issue("BLOCKER", "AC-03", "Duplicate pair or unrecognized coverage status")
        status = "BLOCKED" if problems else "AUTO_EVIDENCE"
    else:
        status = "MISSING_EVIDENCE"
    add(ac, "AC-03", status, "coverage_matrix.tsv + schema_audit.json", "Проверить биологический смысл покрытия 65 × 13")

    if source.has("closeout/source_api_reconciliation.tsv"):
        records = source.tsv("closeout/source_api_reconciliation.tsv", {"status"})
        errors_api = [r for r in records if r["status"] != "MATCH"]
        if errors_api:
            issue("BLOCKER", "AC-04", f"Source→API mismatch/non-match rows: {len(errors_api)}")
        status = "BLOCKED" if errors_api else "AUTO_EVIDENCE" if records else "MISSING_EVIDENCE"
    else:
        status = "MISSING_EVIDENCE"
    add(ac, "AC-04", status, "closeout/source_api_reconciliation.tsv", "Вручную проверить PYZ-001, -035, -051, -065, -055 на экране")

    shots = [x for x in source.names() if x.startswith("closeout/ui/screenshots/") and x.lower().endswith(".png")]
    invalid_png = [p for p in shots if not source.read(p).startswith(b"\x89PNG\r\n\x1a\n")]
    if invalid_png:
        issue("BLOCKER", "AC-05", f"Screenshot files without PNG signatures: {len(invalid_png)}")
    if source.has("closeout/ui/browser_results.json"):
        browser_results = source.json("closeout/ui/browser_results.json")
        ui_warnings = [row for row in browser_results.get("screenshots", [])
                       if row.get("status") != "AUTO_PASS"]
        if ui_warnings:
            issue("MANUAL", "AC-05", f"Browser technical checks with warnings: {len(ui_warnings)}")
    if source.has("closeout/browser_status.json"):
        browser = source.json("closeout/browser_status.json")
        reported_n = browser.get("n_screenshots", 0)
    else:
        browser = {}
        reported_n = 0
    if len(shots) < 9 or reported_n != 9:
        issue("MANUAL", "AC-05", f"Screenshots: {len(shots)} files; browser reported {reported_n}/9")
    add(ac, "AC-05", "BLOCKED" if invalid_png else "PENDING_MANUAL" if len(shots) >= 9 and reported_n == 9 else "MISSING_EVIDENCE",
        "closeout/ui/screenshots/*.png + closeout/browser_status.json", "Оценить 9 PNG, клавиатуру, прокрутку, данные и доступность")

    if source.has("closeout/cold_api.json"):
        cold = source.json("closeout/cold_api.json")
        cold_status = str(cold.get("status", "NOT_RUN"))
        if cold_status == "NOT_RUN":
            issue("MANUAL", "AC-06", "Cold API was not run")
        status = "PENDING_MANUAL"  # no asserted latency budget in E0 requirements
    else:
        cold_status = "MISSING"
        status = "MISSING_EVIDENCE"
    add(ac, "AC-06", status, "performance.tsv + closeout/cold_api.json", "Проверить p50/p95, холодный запуск, кэш и критерии задержки")

    add(ac, "AC-07", "PENDING_MANUAL", "scientific_reference_examples.tsv + Э1-P0-B3A", "Проверить структуры, методы Ki/Kd/IC50, независимость опытов; не переносить активность на PYZ")
    add(ac, "AC-08", "BLOCKED" if errors else "AUTO_EVIDENCE" if len(hash_details) >= 2 * len(hashes) and hashes else "MISSING_EVIDENCE",
        "closeout/source_sha_entry.tsv + source_sha_after.tsv", "Решить расхождения и подписать контроль неизменности научных источников")
    add(ac, "AC-09", "BLOCKED" if errors else "PENDING_MANUAL", "E0_ACCEPTANCE.md + closeout/E0_ACCEPTANCE_CLOSEOUT.md", "Только владелец проекта утверждает акт и разрешает объединение веток")

    counts = Counter(x["technical_status"] for x in ac)
    blockers = [x for x in findings if x["severity"] == "BLOCKER"]
    report = {
        "contract": VERSION,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_path": str(source.path),
        "input_sha256": digest(source.path) if source.path.is_file() else None,
        "criteria_count": len(ac),
        "scientific_source_paths": len(hashes),
        "scientific_hash_comparison_rows": len(hash_details),
        "fresh_live_sha_rows": live_count,
        "fresh_live_sha_checked": bool(science_root and live_count == len(hashes) and not errors),
        "screenshot_files": len(shots),
        "browser_reported_screenshots": reported_n,
        "cold_api_status": cold_status,
        "counts_by_technical_status": dict(counts),
        "blocking_findings": len(blockers),
        "closeout_missing_files": missing_closeout,
        "scientific_status": "NOT_ACCEPTED_MANUAL_REVIEW_REQUIRED",
        "merge_permitted": False,
        "final_status": "BLOCKED" if blockers else "INCOMPLETE_EVIDENCE" if any(x["technical_status"] == "MISSING_EVIDENCE" for x in ac) else "PENDING_OWNER_ACCEPTANCE",
        "notes": [
            "AUTO_EVIDENCE is not a signed scientific acceptance.",
            "Archived SHA rows confirm internal consistency; fresh local rehash is optional and marked explicitly.",
            "Missing E0 evidence is not inferred to have passed from E1 results or GitHub CI.",
            "ADMET values are not used for ranking or safety conclusions.",
        ],
    }
    return report, ac, findings, hash_details


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict]):
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)



# Whitelisted E0 *evidence*, not original scientific source tables/Parquet.
REVIEW_MATERIALS = {
    "E0_ACCEPTANCE.md", "input_inventory.tsv", "schema_audit.json",
    "coverage_matrix.tsv", "api_cases.tsv", "performance.tsv",
    "ui_checklist.tsv", "defects.md", "environment.json",
    "closeout/E0_ACCEPTANCE_CLOSEOUT.md", "closeout/browser_status.json",
    "closeout/cold_api.json", "closeout/cache_diagnostics.json",
    "closeout/ui/browser_results.json",
    "closeout/ui/ui_checklist_automated.tsv",
    "closeout/source_api_reconciliation.tsv",
    "closeout/ui_expected.json",
    "closeout/scientific_reference_examples.tsv",
    "closeout/cell_context_observations.tsv",
}


def export_review_evidence(source: Evidence, out: Path) -> dict:
    """Copy only whitelisted E0 evidence into new report, never source files.

    Preserve original bytes and store their SHA-256 for manual AC-02/05/06 review.
    Evidence() has already rejected unsafe archive names, duplicate entries and
    symlinks. We apply an explicit allowlist again before writing output.
    """
    src_names = sorted(name for name in source.names()
                       if name in REVIEW_MATERIALS or
                       (name.startswith("closeout/ui/screenshots/") and
                        name.lower().endswith(".png")))
    records = []
    root = out.resolve()
    for name in src_names:
        rel = norm(name)
        dest = (out / "review_evidence" / rel).resolve()
        if not dest.is_relative_to(root):
            raise ValueError("Evidence path escaped the report output directory")
        raw = source.read(name)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(raw)
        records.append({"name": rel, "size_bytes": len(raw),
                        "sha256": hashlib.sha256(raw).hexdigest()})
    manifest = {"scope": "whitelisted E0 source report files, no scientific input edits",
                "files": records, "screenshots": sum(
                    item["name"].startswith("closeout/ui/screenshots/") for item in records)}
    (out / "review_evidence_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    diagnostic_rows = []
    if source.has("closeout/ui/browser_results.json"):
        browser = source.json("closeout/ui/browser_results.json")
        for item in browser.get("screenshots", []):
            causes = []
            if not item.get("saved"): causes.append("SCREENSHOT_NOT_SAVED")
            if item.get("http_status") != 200: causes.append("HTTP_NOT_200")
            if item.get("heading_matches") is not True: causes.append("HEADING_MISMATCH")
            if item.get("horizontal_overflow") is True: causes.append("HORIZONTAL_OVERFLOW")
            if item.get("keyboard_focus") is not True: causes.append("KEYBOARD_FOCUS_MISSING")
            if item.get("console_errors"): causes.append("BROWSER_ERRORS")
            if item.get("status") != "AUTO_PASS" and not causes:
                causes.append("UNCLASSIFIED_WARNING")
            diagnostic_rows.append({
                "route": item.get("route", ""), "width": item.get("width", ""),
                "height": item.get("height", ""), "status": item.get("status", ""),
                "http_status": item.get("http_status", ""),
                "heading_matches": item.get("heading_matches", ""),
                "horizontal_overflow": item.get("horizontal_overflow", ""),
                "keyboard_focus": item.get("keyboard_focus", ""),
                "console_errors": json.dumps(item.get("console_errors", []), ensure_ascii=False),
                "screenshot": item.get("screenshot", ""),
                "warning_causes": "|".join(causes),
            })
    write_tsv(out / "ui_diagnostics.tsv", (
        "route", "width", "height", "status", "http_status",
        "heading_matches", "horizontal_overflow", "keyboard_focus",
        "console_errors", "screenshot", "warning_causes"), diagnostic_rows)
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--e0", type=Path, required=True, help="Existing E0 report ZIP or folder")
    parser.add_argument("--out-dir", type=Path, required=True, help="New directory under MCL_Project/reports")
    parser.add_argument("--science-root", type=Path, default=None, help="Optional independent rehash of original source files")
    args = parser.parse_args(argv)
    src = args.e0.resolve()
    out = args.out_dir.resolve()
    if out.exists():
        raise FileExistsError(f"Will not overwrite report: {out}")
    if out.is_relative_to(src) or src.is_relative_to(out):
        raise ValueError("Output and input paths must be separate")
    evidence = Evidence(src)
    try:
        report, criteria, findings, hashes = analyze(evidence, args.science_root)
    finally:
        evidence.close()
    out.mkdir(parents=True)
    review_source = Evidence(src)
    try:
        export_review_evidence(review_source, out)
    finally:
        review_source.close()
    write_tsv(out / "acceptance_matrix.tsv", ("criterion", "technical_status", "evidence", "required_manual_review", "manual_decision", "reviewer", "review_date"), criteria)
    write_tsv(out / "findings.tsv", ("severity", "criterion", "detail"), findings)
    write_tsv(out / "source_hash_reconciliation.tsv", ("phase", "path", "status", "sha_matches"), hashes)
    (out / "gate_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    manual = [
        "# Э0: консолидированный протокол оценки готовности", "",
        f"Источник: {src}", f"Статус: {report['final_status']}",
        "", "**Э0 НЕ УТВЕРЖДЁН автоматически. Слияние PR запрещено до ручного решения.**", "",
        "## Критерии AC-01…AC-09", "",
        "| Критерий | Технические данные | Ручная проверка |", "|---|---|---|",
    ]
    for r in criteria:
        manual.append(f"| {r['criterion']} | {r['technical_status']} | {r['required_manual_review']} |")
    manual.extend(["", "## Выявленные вопросы", ""])
    if not findings:
        manual.append("- Автоматических проблем не выявлено; научная приёмка всё равно требуется.")
    else:
        manual.extend(f"- **{x['severity']} / {x['criterion']}** — {x['detail']}" for x in findings)
    manual.extend(["", "## Следующая процедура", "",
                   "1. Сверить исходный локальный отчёт Э0 с содержимым этого архива.",
                   "2. Проверить вручную девять скриншотов и первичные научные источники.",
                   "3. Зафиксировать решения по каждому AC и подпись владельца.",
                   "4. Только после принятия Э0 рассмотреть PR №5; Э1 PR №6–9 оценивать по очереди и отдельно.",
                   "5. ADMET не включать в скоринг; будущий внешний модуль рассматривать отдельным этапом.", ""])
    (out / "E0_ACCEPTANCE_GATE_RU.md").write_text("\n".join(manual), encoding="utf-8")
    print(json.dumps({"status": report["final_status"], "criteria": len(criteria), "blockers": len([x for x in findings if x["severity"] == "BLOCKER"]), "merge_permitted": False}, ensure_ascii=False))
    return report


if __name__ == "__main__":
    main()