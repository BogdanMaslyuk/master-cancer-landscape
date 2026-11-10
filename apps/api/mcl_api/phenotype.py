from __future__ import annotations

"""Phenotype-first experiment intake, P0/P1.

Only own measurements belong here. PRISM, CRISPR and chemical-neighbor
hypotheses must never populate compound-specific viability endpoints.
No fit or selectivity result is fabricated from insufficient raw data.
"""

import csv
import hashlib
import io
import json
import math
import os
import re
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ALLOWED_CONTROLS = frozenset({
    "treated", "vehicle", "untreated", "positive_control", "blank",
    "compound_optical_blank", "compound_mtt_blank", "baseline_n0",
})
MODEL_ID_RE = re.compile(r"^[A-Za-z0-9_-]{1,64}$")
PLATE_RE = re.compile(r"^[A-Za-z0-9_-]{1,64}$")
WELL_RE = re.compile(r"^[A-H](?:0[1-9]|1[0-2])$")
EXPERIMENT_ID_RE = re.compile(r"^PHENO-[a-f0-9]{12}$")
REQUIRED_COLUMNS = (
    "plate_id", "well_id", "biological_run_id", "cell_line_id",
    "compound_id", "batch_id", "concentration_value", "concentration_unit",
    "exposure_h", "assay_type", "control_type", "raw_signal",
)
UNIT_TO_UM = {"nm": 0.001, "um": 1.0, "µm": 1.0, "μm": 1.0, "mm": 1000.0}


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _atomic_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    try:
        temp.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def _number(value: str, *, field: str, row: int, required: bool = True) -> float | None:
    if not value.strip():
        if required:
            raise ValueError(f"CSV row {row}: missing {field}")
        return None
    try:
        v = float(value)
    except ValueError as exc:
        raise ValueError(f"CSV row {row}: {field} must be numeric") from exc
    if not math.isfinite(v):
        raise ValueError(f"CSV row {row}: {field} must be finite")
    return v


class PhenotypeStore:
    def __init__(self, root: Path, runtime_root: Path | None = None):
        self.root = Path(root).resolve()
        # The safe default is the current checkout. Never point this at
        # MCL_SCIENTIFIC_DATA_ROOT or the user's original scientific archive.
        self.runtime = Path(runtime_root).resolve() if runtime_root else self.root / "data" / "runtime" / "phenotype"
        self.cell_registry = self.root / "config" / "phenotype_cell_models.tsv"
        self.compound_registry = self.root / "config" / "own_compounds.tsv"

    @staticmethod
    def _tsv(path: Path) -> list[dict[str, str]]:
        if not path.is_file():
            return []
        with path.open(encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def models(self) -> dict[str, Any]:
        rows = self._tsv(self.cell_registry)
        return {
            "items": [
                {**r, "assay_ready": r.get("assay_ready", "").lower() == "true"}
                for r in rows
            ],
            "total": len(rows),
            "note_ru": (
                "Запись в реестре не подтверждает наличие флакона, идентичность, "
                "отсутствие микоплазмы или готовность к опыту."
            ),
        }

    def compounds(self) -> dict[str, Any]:
        rows = self._tsv(self.compound_registry)
        return {
            "items": [
                {
                    "compound_id": r["own_compound_id"],
                    "preferred_name": r.get("preferred_name") or r["own_compound_id"],
                    "standardized_smiles": r.get("standardized_smiles"),
                    "physical_batch_status": "unverified",
                    "direct_hek_mrc_activity": "not_measured",
                }
                for r in rows if r.get("own_compound_id")
            ],
            "total": sum(bool(r.get("own_compound_id")) for r in rows),
            "note_ru": "Структура в каталоге не означает физическое наличие или подтверждённую чистоту вещества.",
        }

    def _directory(self, experiment_id: str) -> Path:
        if not EXPERIMENT_ID_RE.fullmatch(experiment_id):
            raise ValueError("Invalid phenotype experiment ID")
        return self.runtime / "experiments" / experiment_id

    def create(self, payload: dict[str, Any]) -> dict[str, Any]:
        models = {r["cell_line_id"] for r in self.models()["items"]}
        compounds = {r["compound_id"] for r in self.compounds()["items"]}
        model_ids = payload.get("cell_line_ids")
        compound_ids = payload.get("compound_ids")
        if not isinstance(model_ids, list) or not model_ids or len(model_ids) > 20:
            raise ValueError("Select 1–20 exact cell line IDs")
        if not isinstance(compound_ids, list) or not compound_ids or len(compound_ids) > 65:
            raise ValueError("Select 1–65 compound IDs")
        if not all(isinstance(s, str) for s in model_ids + compound_ids):
            raise ValueError("Model and compound IDs must be strings")
        if len(set(model_ids)) != len(model_ids) or len(set(compound_ids)) != len(compound_ids):
            raise ValueError("Duplicate cell or compound ID")
        unknown_models = sorted(set(model_ids) - models)
        unknown_compounds = sorted(set(compound_ids) - compounds)
        if unknown_models or unknown_compounds:
            raise ValueError(f"Unknown models: {unknown_models}; unknown PYZ: {unknown_compounds}")
        title = str(payload.get("title") or "").strip()
        protocol_id = str(payload.get("protocol_id") or "").strip()
        if not title or len(title) > 180:
            raise ValueError("A short experiment title is required")
        if not protocol_id or len(protocol_id) > 100:
            raise ValueError("Approved protocol identifier is required; use DRAFT for an unapproved plan")
        assay_type = str(payload.get("assay_type") or "MTT").strip().upper()
        if assay_type not in {"MTT", "CELL_COUNT"}:
            raise ValueError("Assay must be MTT or CELL_COUNT in this milestone")
        exposure = float(payload.get("exposure_h") or 0)
        if not math.isfinite(exposure) or not 0 < exposure <= 720:
            raise ValueError("Exposure must be within 0–720 hours")
        experiment_id = "PHENO-" + uuid.uuid4().hex[:12]
        directory = self._directory(experiment_id)
        directory.mkdir(parents=True, exist_ok=False)
        item = {
            "experiment_id": experiment_id,
            "title": title,
            "created_at": _utc(),
            "cell_line_ids": model_ids,
            "compound_ids": compound_ids,
            "protocol_id": protocol_id,
            "assay_type": assay_type,
            "exposure_h": exposure,
            "experiment_status": "draft",
            "readiness": "unverified",
            "interpretation_ru": (
                "План эксперимента, не измеренная цитотоксичность. "
                "Данные HEK293T/FT не подтверждают эффективность против рака почки."
            ),
        }
        _atomic_json(directory / "experiment.json", item)
        return item

    def detail(self, experiment_id: str) -> dict[str, Any]:
        directory = self._directory(experiment_id)
        path = directory / "experiment.json"
        if not path.is_file():
            raise ValueError("Phenotype experiment not found")
        experiment = json.loads(path.read_text(encoding="utf-8"))
        imports = sorted(
            [json.loads(x.read_text(encoding="utf-8")) for x in (directory / "imports").glob("*.json")],
            key=lambda row: (row["created_at"], row["sha256"]),
        ) if (directory / "imports").is_dir() else []
        statuses = [x["qc"]["status"] for x in imports]
        status = (
            "qc_fail" if "qc_fail" in statuses else
            "qc_pending" if imports else "draft"
        )
        return {
            **experiment,
            "experiment_status": status,
            "imports": imports,
            "imported_files_n": len(imports),
            "raw_wells_n": sum(x["wells_n"] for x in imports),
            "ic50_status": "not_estimable_no_validated_concentration_fit",
            "gr50_status": "not_available_no_validated_cell_counts",
            "selectivity_status": "not_measured",
            "verified_result": False,
        }

    def _parse_csv(self, csv_text: str, experiment: dict[str, Any]) -> list[dict[str, Any]]:
        if not csv_text.strip():
            raise ValueError("Empty CSV")
        reader = csv.DictReader(io.StringIO(csv_text.lstrip("\ufeff")), restkey="__extra__")
        headers = reader.fieldnames or []
        if len(headers) != len(set(headers)) or "__extra__" in headers:
            raise ValueError("Duplicate/reserved CSV column")
        missing = sorted(set(REQUIRED_COLUMNS) - set(headers))
        if missing:
            raise ValueError(f"Missing columns: {', '.join(missing)}")
        out: list[dict[str, Any]] = []
        seen: set[tuple[str, str]] = set()
        for index, row in enumerate(reader, 2):
            if index > 50002:
                raise ValueError("CSV exceeds 50000 rows")
            if not row or not any(str(v or "").strip() for v in row.values()):
                continue
            if row.get("__extra__") is not None or any(v is None for v in row.values()):
                raise ValueError(f"CSV row {index}: malformed column count")
            r = {str(k): str(v).strip() for k, v in row.items()}
            plate = r["plate_id"]
            well = r["well_id"].upper()
            if not PLATE_RE.fullmatch(plate) or not WELL_RE.fullmatch(well):
                raise ValueError(f"CSV row {index}: invalid plate or 96-well address")
            unique = plate, well
            if unique in seen:
                raise ValueError(f"CSV row {index}: duplicate well {plate}/{well}")
            seen.add(unique)
            if not MODEL_ID_RE.fullmatch(r["biological_run_id"]):
                raise ValueError(f"CSV row {index}: biological_run_id is required")
            cell = r["cell_line_id"]
            if cell not in experiment["cell_line_ids"]:
                raise ValueError(f"CSV row {index}: unexpected cell line {cell!r}")
            if r["assay_type"].upper() != experiment["assay_type"]:
                raise ValueError(f"CSV row {index}: assay differs from experiment")
            ctype = r["control_type"].lower()
            if ctype not in ALLOWED_CONTROLS:
                raise ValueError(f"CSV row {index}: unknown control_type {ctype!r}")
            compound = r["compound_id"]
            if compound and compound not in experiment["compound_ids"]:
                raise ValueError(f"CSV row {index}: unknown compound {compound!r}")
            if ctype in {"treated", "compound_mtt_blank", "compound_optical_blank"}:
                if not compound or not r["batch_id"]:
                    raise ValueError(f"CSV row {index}: compound and batch are required")
            signal = _number(r["raw_signal"], field="raw_signal", row=index)
            dose_value = _number(r["concentration_value"], field="concentration_value", row=index, required=False)
            unit = r["concentration_unit"].lower()
            if unit not in UNIT_TO_UM:
                raise ValueError(f"CSV row {index}: unknown concentration unit {r['concentration_unit']!r}")
            if dose_value is not None and dose_value < 0:
                raise ValueError(f"CSV row {index}: negative dose")
            if ctype in {"treated", "compound_mtt_blank", "compound_optical_blank"} and (dose_value is None or dose_value <= 0):
                raise ValueError(f"CSV row {index}: compound exposure needs a positive dose")
            exposure = _number(r["exposure_h"], field="exposure_h", row=index)
            if exposure != experiment["exposure_h"]:
                raise ValueError(f"CSV row {index}: exposure differs from experiment")
            if "reader_wavelength_nm" in r and r["reader_wavelength_nm"]:
                wavelength = _number(r["reader_wavelength_nm"], field="reader_wavelength_nm", row=index)
                if wavelength <= 0:
                    raise ValueError(f"CSV row {index}: invalid wavelength")
            else:
                wavelength = None
            out.append({
                **r,
                "well_id": well,
                "control_type": ctype,
                "assay_type": experiment["assay_type"],
                "raw_signal": signal,
                "concentration_um": dose_value * UNIT_TO_UM[unit] if dose_value is not None else None,
                "reader_wavelength_nm": wavelength,
                "source_csv_row": index,
            })
        if not out:
            raise ValueError("CSV has no measurements")
        return out

    @staticmethod
    def _qc(rows: list[dict[str, Any]], assay_type: str) -> dict[str, Any]:
        flags: list[dict[str, str]] = []
        groups: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
        for row in rows:
            key = row["plate_id"], row["biological_run_id"], row["cell_line_id"]
            groups.setdefault(key, []).append(row)
            if row["raw_signal"] < 0:
                flags.append({"severity": "fail", "code": "negative_raw_signal", "well": row["well_id"]})
        for (plate, run, cell), group in groups.items():
            controls = {r["control_type"] for r in group}
            scope = f"{plate}/{run}/{cell}"
            if "vehicle" not in controls:
                flags.append({"severity": "fail", "code": "missing_vehicle_control", "scope": scope})
            if "blank" not in controls:
                flags.append({"severity": "pending", "code": "missing_media_blank", "scope": scope})
            if "positive_control" not in controls:
                flags.append({"severity": "pending", "code": "missing_positive_control", "scope": scope})
            if assay_type == "MTT" and any(r["control_type"] == "treated" for r in group):
                treated_compounds = {r["compound_id"] for r in group if r["control_type"] == "treated"}
                for compound in sorted(treated_compounds):
                    present = {r["control_type"] for r in group if r["compound_id"] == compound}
                    for required in ("compound_optical_blank", "compound_mtt_blank"):
                        if required not in present:
                            flags.append({
                                "severity": "pending", "code": "interference_control_missing",
                                "scope": scope, "compound_id": compound, "control": required,
                            })
        # Without approved SOP thresholds and a scientifically reviewed blank
        # correction, numerical MTT interference cannot be declared absent.
        return {
            "status": "qc_fail" if any(x["severity"] == "fail" for x in flags) else "qc_pending",
            "flags": flags,
            "group_count": len(groups),
            "wells_n": len(rows),
            "sop_thresholds_verified": False,
            "assay_interference": "not_excluded" if assay_type == "MTT" else "not_assessed",
            "verified_result": False,
            "interpretation_ru": "Сохранён первичный сигнал; научная готовность и клеточная селективность пока не подтверждены.",
        }

    def import_csv(
        self, experiment_id: str, csv_text: str | None, filename: str,
        *, raw_bytes: bytes | None = None,
    ) -> dict[str, Any]:
        # Prefer the exact original browser bytes over a decoded/re-encoded
        # string. This preserves the source SHA-256, CRLF and UTF-8 BOM.
        if raw_bytes is None:
            if not isinstance(csv_text, str):
                raise ValueError("Missing raw CSV content")
            raw = csv_text.encode("utf-8")
        else:
            if not isinstance(raw_bytes, bytes):
                raise ValueError("Raw CSV must be bytes")
            raw = raw_bytes
        if not raw or len(raw) > 8_000_000:
            raise ValueError("CSV must be UTF-8 and smaller than 8 MB")
        try:
            parsed_text = raw.decode("utf-8-sig")
        except UnicodeDecodeError as exc:
            raise ValueError("CSV must be UTF-8 encoded") from exc
        if "\x00" in parsed_text:
            raise ValueError("CSV contains NUL bytes")
        if not str(filename).lower().endswith(".csv"):
            raise ValueError("P1 supports CSV only; XLSX ingestion is not implemented yet")
        experiment = self.detail(experiment_id)
        rows = self._parse_csv(parsed_text, experiment)
        old_plate_ids = {
            plate for entry in experiment["imports"] for plate in entry["plate_ids"]
        }
        new_plate_ids = {row["plate_id"] for row in rows}
        sha = hashlib.sha256(raw).hexdigest()
        directory = self._directory(experiment_id) / "imports"
        same = next((x for x in experiment["imports"] if x["sha256"] == sha), None)
        if same is not None:
            return {**same, "reused_existing_import": True}
        overlap = old_plate_ids & new_plate_ids
        if overlap:
            raise ValueError("Plate already imported; previous immutable measurements remain intact: " + ", ".join(sorted(overlap)))
        qc = self._qc(rows, experiment["assay_type"])
        directory.mkdir(parents=True, exist_ok=True)
        raw_path = directory / f"{sha}.csv"
        normalized_path = directory / f"{sha}.rows.jsonl"
        try:
            with raw_path.open("xb") as stream:
                stream.write(raw)
        except FileExistsError:
            if hashlib.sha256(raw_path.read_bytes()).hexdigest() != sha:
                raise ValueError("Stored raw file failed SHA-256 integrity check")
        # The derived row ledger is separately versioned and never modifies
        # the uploaded original. Internal line-dose calculations use this
        # explicit normalized representation, not a second implicit parser.
        if normalized_path.exists():
            raise ValueError("A normalized ledger already exists without matching import manifest")
        with normalized_path.open("x", encoding="utf-8", newline="\n") as output:
            for row in rows:
                output.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + "\n")
        normalized_sha = hashlib.sha256(normalized_path.read_bytes()).hexdigest()
        record = {
            "normalization_version": "mcl-phenotype-csv-v1",
            "normalized_rows_sha256": normalized_sha,
            "normalized_rows_file": normalized_path.name,
            "sha256": sha,
            "filename": str(filename).replace("\\", "/").split("/")[-1],
            "created_at": _utc(),
            "plate_ids": sorted(new_plate_ids),
            "wells_n": len(rows),
            "biological_runs_n": len({r["biological_run_id"] for r in rows}),
            "cell_line_ids": sorted({r["cell_line_id"] for r in rows}),
            "compound_ids": sorted({r["compound_id"] for r in rows if r["compound_id"]}),
            "qc": qc,
            "fit_status": "not_performed",
            "selectivity_status": "not_measured",
            "reused_existing_import": False,
        }
        _atomic_json(directory / f"{sha}.json", record)
        return record
