from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from mcl_api.phenotype import PhenotypeStore


def sample_store(tmp_path: Path) -> PhenotypeStore:
    config = tmp_path / "config"
    config.mkdir(exist_ok=True)
    (config / "phenotype_cell_models.tsv").write_text(
        "cell_line_id\tdisplay_name\tmodel_role\tinventory_status\tassay_ready\n"
        "HEK293\tHEK293\thuman_technical\tunverified\tfalse\n"
        "HEK293T\tHEK293T\thuman_technical\tregistered_only\tfalse\n"
        "HEK293FT\tHEK293FT\thuman_technical\tregistered_only\tfalse\n"
        "MRC-5\tMRC-5\thuman_normal_diploid\tplanned_unverified\tfalse\n",
        encoding="utf-8",
    )
    (config / "own_compounds.tsv").write_text(
        "own_compound_id\tpreferred_name\tstandardized_smiles\n"
        "PYZ-035\tPYZ-035\tCCO\n"
        "PYZ-055\tPYZ-055\tCCN\n",
        encoding="utf-8",
    )
    return PhenotypeStore(tmp_path)


def create(store: PhenotypeStore, **overrides):
    req = {
        "title": "Pilot HEK293T vs MRC-5",
        "protocol_id": "DRAFT",
        "cell_line_ids": ["HEK293T", "MRC-5"],
        "compound_ids": ["PYZ-035"],
        "exposure_h": 72.0,
        "assay_type": "MTT",
    }
    req.update(overrides)
    return store.create(req)


def record(*, plate="P1", well="A01", run="RUN1", cell="HEK293T",
           compound="PYZ-035", batch="B1", dose="10", unit="uM",
           ctrl="treated", signal="0.65", method="MTT", hours="72") -> str:
    return ",".join([plate, well, run, cell, compound, batch, dose, unit,
                     hours, method, ctrl, signal])


HEADER = ("plate_id,well_id,biological_run_id,cell_line_id,compound_id,"
          "batch_id,concentration_value,concentration_unit,exposure_h,"
          "assay_type,control_type,raw_signal")


def csv_file(*rows: str) -> str:
    return HEADER + "\n" + "\n".join(rows) + "\n"


def complete_controls(cell="HEK293T", plate="P1") -> list[str]:
    return [
        record(plate=plate, well="A01", cell=cell, ctrl="vehicle",
               compound="", batch="", dose="0", signal="1.2"),
        record(plate=plate, well="A02", cell=cell, ctrl="blank",
               compound="", batch="", dose="0", signal="0.04"),
        record(plate=plate, well="A03", cell=cell, ctrl="positive_control",
               compound="", batch="", dose="0", signal="0.15"),
        record(plate=plate, well="A04", cell=cell, ctrl="compound_optical_blank", signal="0.08"),
        record(plate=plate, well="A05", cell=cell, ctrl="compound_mtt_blank", signal="0.12"),
        record(plate=plate, well="A06", cell=cell, ctrl="treated", signal="0.65"),
    ]


def test_cell_model_variants_are_distinct_and_mrc_is_not_verified(tmp_path: Path):
    store = sample_store(tmp_path)
    models = {m["cell_line_id"]: m for m in store.models()["items"]}
    assert set(models) == {"HEK293", "HEK293T", "HEK293FT", "MRC-5"}
    assert models["MRC-5"]["inventory_status"] == "planned_unverified"
    assert models["MRC-5"]["assay_ready"] is False
    assert models["HEK293T"]["model_role"] == "human_technical"
    assert models["HEK293FT"]["model_role"] == "human_technical"
    assert models["HEK293T"]["cell_line_id"] != models["HEK293FT"]["cell_line_id"]


def test_experiment_creation_is_phenotype_first_without_tanimoto_or_depmap(tmp_path: Path):
    store = sample_store(tmp_path)
    created = create(store)
    assert created["experiment_id"].startswith("PHENO-")
    assert created["readiness"] == "unverified"
    assert store.detail(created["experiment_id"])["selectivity_status"] == "not_measured"
    assert store.detail(created["experiment_id"])["ic50_status"].startswith("not_estimable")
    assert store.detail(created["experiment_id"])["imported_files_n"] == 0
    assert "target_gene" not in created


def test_reject_unknown_cell_compound_duplicate_variants_and_invalid_assay(tmp_path: Path):
    store = sample_store(tmp_path)
    for change in (
        {"cell_line_ids": ["HEK293_unknown"]},
        {"compound_ids": ["PYZ-999"]},
        {"cell_line_ids": ["HEK293T", "HEK293T"]},
        {"assay_type": "UNTYPED"},
        {"exposure_h": 0},
    ):
        with pytest.raises(ValueError):
            create(store, **change)


def test_csv_import_keeps_immutable_source_and_does_not_invent_ic50(tmp_path: Path):
    store = sample_store(tmp_path)
    experiment_id = create(store)["experiment_id"]
    raw = csv_file(*complete_controls())
    first = store.import_csv(experiment_id, raw, "plate.csv")
    second = store.import_csv(experiment_id, raw, "plate.csv")
    assert first["reused_existing_import"] is False
    assert second["reused_existing_import"] is True
    assert first["sha256"] == hashlib.sha256(raw.encode()).hexdigest()
    path = tmp_path / "data" / "runtime" / "phenotype" / "experiments" / experiment_id / "imports"
    assert (path / f"{first['sha256']}.csv").read_text(encoding="utf-8") == raw
    assert first["wells_n"] == 6
    assert first["qc"]["status"] == "qc_pending"
    assert first["qc"]["sop_thresholds_verified"] is False
    detail = store.detail(experiment_id)
    assert detail["raw_wells_n"] == 6
    assert detail["imported_files_n"] == 1
    assert detail["ic50_status"] == "not_estimable_no_validated_concentration_fit"
    assert detail["selectivity_status"] == "not_measured"
    assert detail["verified_result"] is False


def test_import_rejects_duplicate_well_and_keeps_directory_clean(tmp_path: Path):
    store = sample_store(tmp_path)
    experiment_id = create(store)["experiment_id"]
    bad = csv_file(record(), record(signal="0.1"))
    with pytest.raises(ValueError, match="duplicate well"):
        store.import_csv(experiment_id, bad, "plate.csv")
    assert store.detail(experiment_id)["imported_files_n"] == 0


@pytest.mark.parametrize("bad, expected", [
    (csv_file(record(unit="huh")), "unknown concentration unit"),
    (csv_file(record(well="Z99")), "invalid plate"),
    (csv_file(record(cell="HEK293FT")), "unexpected cell line"),
    (csv_file(record(compound="PYZ-999")), "unknown compound"),
    (csv_file(record(dose="-1")), "negative dose"),
    (csv_file(record(signal="nan")), "must be finite"),
    (csv_file(record(hours="120")), "exposure differs"),
])
def test_csv_bad_data_are_explicitly_rejected(tmp_path: Path, bad: str, expected: str):
    store = sample_store(tmp_path)
    experiment_id = create(store)["experiment_id"]
    with pytest.raises(ValueError, match=expected):
        store.import_csv(experiment_id, bad, "plate.csv")
    assert store.detail(experiment_id)["raw_wells_n"] == 0


def test_missing_vehicle_and_negative_signal_are_qc_fail_not_silent(tmp_path: Path):
    store = sample_store(tmp_path)
    experiment_id = create(store)["experiment_id"]
    raw = csv_file(record(signal="-0.5"))
    result = store.import_csv(experiment_id, raw, "missing_controls.csv")
    assert result["qc"]["status"] == "qc_fail"
    assert {"negative_raw_signal", "missing_vehicle_control"}.issubset(
        {f["code"] for f in result["qc"]["flags"]}
    )
    assert store.detail(experiment_id)["experiment_status"] == "qc_fail"
    assert store.detail(experiment_id)["verified_result"] is False


def test_mtt_without_optical_and_chemical_blanks_remains_pending(tmp_path: Path):
    store = sample_store(tmp_path)
    experiment_id = create(store)["experiment_id"]
    raw = csv_file(
        record(well="A01", ctrl="vehicle", compound="", batch="", dose="0"),
        record(well="A02", ctrl="treated"),
    )
    result = store.import_csv(experiment_id, raw, "plate.csv")
    assert result["qc"]["status"] == "qc_pending"
    assert result["qc"]["assay_interference"] == "not_excluded"
    assert sum(f["code"] == "interference_control_missing" for f in result["qc"]["flags"]) == 2


def test_same_plate_changed_file_is_rejected_instead_of_overwriting(tmp_path: Path):
    store = sample_store(tmp_path)
    experiment_id = create(store)["experiment_id"]
    store.import_csv(experiment_id, csv_file(*complete_controls()), "plate.csv")
    with pytest.raises(ValueError, match="already imported"):
        store.import_csv(
            experiment_id, csv_file(*complete_controls(cell="MRC-5")), "changed.csv"
        )
    assert store.detail(experiment_id)["imported_files_n"] == 1


def test_separate_biological_plate_does_not_merge_hekt_with_heft(tmp_path: Path):
    store = sample_store(tmp_path)
    experiment_id = create(store, cell_line_ids=["HEK293T", "HEK293FT", "MRC-5"])["experiment_id"]
    store.import_csv(experiment_id, csv_file(*complete_controls()), "t.csv")
    store.import_csv(experiment_id, csv_file(*complete_controls(cell="HEK293FT", plate="P2")), "ft.csv")
    detail = store.detail(experiment_id)
    assert detail["imported_files_n"] == 2
    assert set(detail["cell_line_ids"]) == {"HEK293T", "HEK293FT", "MRC-5"}
    assert detail["selectivity_status"] == "not_measured"


def test_real_registry_is_not_mutated_and_contains_65_pyz():
    root = Path(__file__).resolve().parents[3]
    store = PhenotypeStore(root)
    assert store.compounds()["total"] == 65
    assert store.models()["total"] >= 5
    assert next(m for m in store.models()["items"] if m["cell_line_id"] == "MRC-5")["assay_ready"] is False


def test_api_has_explicit_response_schema_and_no_claimed_ic50():
    from mcl_api.main import app
    schema = app.openapi()
    for path in ("/api/phenotype/models", "/api/phenotype/compounds",
                 "/api/phenotype/experiments/{experiment_id}"):
        response = schema["paths"][path]["get"]["responses"]["200"]["content"]["application/json"]["schema"]
        assert response.get("$ref", "").startswith("#/components/schemas/Phenotype")
    assert schema["paths"]["/api/phenotype/experiments"]["post"]["requestBody"]["content"]["application/json"]["schema"]["$ref"].endswith("PhenotypeExperimentCreateRequest")
