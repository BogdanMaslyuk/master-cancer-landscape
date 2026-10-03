from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import typer
import yaml

from mcl.config import project_root
from mcl.models import MappingStatus
from mcl.normalize.targets import load_target_seeds, normalize_targets
from mcl.qc.targets import target_normalization_qc
from mcl.reporting import create_manifest, create_opentargets_manifest, create_depmap_manifest
from mcl.sources.hgnc import fetch_current_snapshot
from mcl.sources.depmap_sync import sync_depmap
from mcl.analysis.opentargets import collect_opentargets_wave1, load_wave1_pairs, validate_opentargets_diseases
from mcl.analysis.depmap import (
    analyze_depmap_wave1,
    analyze_depmap_kras_sensitivity,
    validate_depmap_inputs,
)
from mcl.analysis.depmap_genomewide import analyze_depmap_genomewide
from mcl.qc.opentargets import opentargets_qc, opentargets_disease_mapping_qc
from mcl.qc.depmap import depmap_qc, depmap_kras_sensitivity_qc
from mcl.qc.depmap_genomewide import depmap_genomewide_qc
from mcl.normalize.disease_mapping import load_opentargets_disease_mappings
from mcl.utils.io import write_json, write_tsv
from mcl.utils.hash import sha256_file
from mcl.utils.logging import configure_logging
from mcl.visualization.depmap_explorer import build_depmap_explorer_html

app = typer.Typer(no_args_is_help=True, help="Master Cancer Landscape pipeline")


def _root(path: Path | None) -> Path:
    return (path or project_root()).resolve()


@app.command("status")
def status(root: Path | None = typer.Option(None, help="Project root")) -> None:
    root = _root(root)
    configure_logging(root)
    seed = root / "data/input/targets_seed.tsv"
    n = len(load_target_seeds(seed)) if seed.exists() else 0
    typer.echo("Master Cancer Landscape Computational Pipeline v0.4.0")
    typer.echo("Milestone 0: IMPLEMENTED")
    m1_done = (root / "data/processed/target_identifiers.tsv").exists()
    m2_done = (root / "data/processed/opentargets_associations.tsv").exists()
    m3_done = (root / "data/processed/depmap_evidence.tsv").exists()
    m31_done = (root / "data/processed/depmap_kras_sensitivity.tsv").exists()
    m32_dir = root / "data/processed/depmap_genomewide"
    m32_done = m32_dir.exists() and any(m32_dir.glob("*_genes.parquet"))
    typer.echo(f"Milestone 1: IMPLEMENTED / {'executed' if m1_done else 'awaiting execution'}")
    typer.echo(f"Milestone 2A+2B: IMPLEMENTED / {'executed' if m2_done else 'awaiting execution'}")
    typer.echo(f"Milestone 3 DepMap: IMPLEMENTED / {'executed' if m3_done else 'awaiting pinned local release files'}")
    typer.echo(f"Milestone 3.1 KRAS sensitivity: IMPLEMENTED / {'executed' if m31_done else 'awaiting execution'}")
    typer.echo(f"Milestone 3.2 Genome-wide Dependency Explorer: IMPLEMENTED / {'executed' if m32_done else 'awaiting execution'}")
    typer.echo(f"Target seeds: {n}")
    pairs = root / "data/input/cancer_target_pairs_wave1.tsv"
    pair_n = len(load_wave1_pairs(pairs)) if pairs.exists() else 0
    typer.echo(f"Wave 1 Cancer×Target pairs: {pair_n}")
    typer.echo("Next execution gate: analyze-depmap-genome-wide --cancer-id <ID> --comparison <mode>")


@app.command("fetch-hgnc")
def fetch_hgnc(
    root: Path | None = typer.Option(None, help="Project root"),
    snapshot_date: str | None = typer.Option(None, help="YYYY-MM-DD; defaults to current UTC date"),
) -> None:
    root = _root(root)
    configure_logging(root)
    meta = fetch_current_snapshot(root / "data/raw", snapshot_date=snapshot_date)
    typer.echo(json.dumps(meta, ensure_ascii=False, indent=2))


@app.command("normalize-targets")
def normalize_targets_cmd(
    root: Path | None = typer.Option(None, help="Project root"),
    seed: Path | None = typer.Option(None, help="Target seed TSV"),
    hgnc_tsv: Path | None = typer.Option(None, help="Pinned HGNC complete-set TSV"),
    withdrawn_tsv: Path | None = typer.Option(None, help="Pinned HGNC withdrawn TSV"),
    source_release: str | None = typer.Option(None, help="Human-readable HGNC snapshot/release label"),
) -> None:
    root = _root(root)
    configure_logging(root)
    seed = seed or (root / "data/input/targets_seed.tsv")

    if hgnc_tsv is None:
        candidates = sorted((root / "data/raw/hgnc").glob("*/hgnc_complete_set.txt"))
        if not candidates:
            raise typer.BadParameter("No HGNC snapshot found. Run `mcl fetch-hgnc` or pass --hgnc-tsv.")
        hgnc_tsv = candidates[-1]
        if withdrawn_tsv is None:
            maybe = hgnc_tsv.parent / "withdrawn.txt"
            withdrawn_tsv = maybe if maybe.exists() else None

    identities, provenance = normalize_targets(seed, hgnc_tsv, withdrawn_tsv, source_release)
    seeds = load_target_seeds(seed)
    qc = target_normalization_qc(seeds, identities)

    processed = root / "data/processed"
    qc_dir = root / "outputs/qc"
    reports = root / "outputs/reports"
    processed.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)

    identity_rows = [x.model_dump(mode="json") for x in identities]
    provenance_rows = [x.model_dump(mode="json") for x in provenance]
    qc_rows = [x.model_dump(mode="json") for x in qc]

    id_tsv = processed / "target_identifiers.tsv"
    prov_tsv = processed / "provenance.tsv"
    write_tsv(id_tsv, identity_rows, list(identity_rows[0].keys()))
    write_tsv(prov_tsv, provenance_rows, list(provenance_rows[0].keys()))

    # Parquet is the canonical processed-data format; TSV is retained for human inspection.
    pd.DataFrame(identity_rows).to_parquet(processed / "target_identifiers.parquet", index=False)
    pd.DataFrame(provenance_rows).to_parquet(processed / "provenance.parquet", index=False)

    write_json(qc_dir / "target_normalization_qc.json", qc_rows)
    write_tsv(qc_dir / "target_normalization_qc.tsv", qc_rows, list(qc_rows[0].keys()))

    errors = sum(1 for x in qc if x.severity == "ERROR")
    manifest = create_manifest(
        root=root,
        seed_path=Path(seed),
        hgnc_path=Path(hgnc_tsv),
        withdrawn_path=Path(withdrawn_tsv) if withdrawn_tsv else None,
        identities_count=len(identities),
        qc_errors=errors,
        outputs=[str(id_tsv), str(prov_tsv), str(qc_dir / "target_normalization_qc.json")],
    )
    write_json(reports / "run_manifest.json", manifest)

    typer.echo(f"Targets normalized: {len(identities)}")
    typer.echo(f"QC errors: {errors}")
    typer.echo(f"Mapping statuses: {dict(pd.Series([x.mapping_status.value for x in identities]).value_counts())}")
    if errors:
        raise typer.Exit(code=2)


@app.command("qc-targets")
def qc_targets(
    root: Path | None = typer.Option(None, help="Project root"),
) -> None:
    root = _root(root)
    configure_logging(root)
    qc_path = root / "outputs/qc/target_normalization_qc.json"
    if not qc_path.exists():
        raise typer.BadParameter("No QC output. Run normalize-targets first.")
    rows = json.loads(qc_path.read_text(encoding="utf-8"))
    errors = [x for x in rows if x["severity"] == "ERROR"]
    warnings = [x for x in rows if x["severity"] == "WARNING"]
    typer.echo(f"ERROR: {len(errors)} | WARNING: {len(warnings)} | TOTAL: {len(rows)}")
    if errors:
        raise typer.Exit(code=2)


def _opentargets_config(root: Path) -> dict:
    cfg = yaml.safe_load((root / "config/source_versions.yaml").read_text(encoding="utf-8")) or {}
    return cfg.get("opentargets", {})


@app.command("validate-ot-diseases")
def validate_ot_diseases(
    root: Path | None = typer.Option(None, help="Project root"),
    source_release: str | None = typer.Option(None, help="Open Targets release label"),
    endpoint: str | None = typer.Option(None, help="GraphQL endpoint"),
) -> None:
    root = _root(root)
    configure_logging(root)
    cfg = _opentargets_config(root)
    source_release = source_release or str(cfg.get("release") or "26.09")
    endpoint = endpoint or str(cfg.get("graphql_endpoint") or "https://api.platform.opentargets.org/api/v4/graphql")
    rows, provenance, meta = validate_opentargets_diseases(
        root=root, source_release=source_release, endpoint=endpoint
    )
    qc_dir = root / "outputs/qc"
    processed = root / "data/processed"
    qc_dir.mkdir(parents=True, exist_ok=True)
    processed.mkdir(parents=True, exist_ok=True)
    write_json(qc_dir / "opentargets_disease_mapping.json", rows)
    write_tsv(qc_dir / "opentargets_disease_mapping.tsv", rows, list(rows[0].keys()))
    write_json(qc_dir / "opentargets_platform_meta.json", meta)
    write_tsv(
        processed / "opentargets_disease_mappings_resolved.tsv",
        rows,
        list(rows[0].keys()),
    )
    prov_rows = [x.model_dump(mode="json") for x in provenance]
    write_tsv(processed / "opentargets_disease_mapping_provenance.tsv", prov_rows, list(prov_rows[0].keys()))
    errors = [r for r in rows if r["status"] == "ERROR"]
    warnings = [r for r in rows if r["status"] == "WARNING"]
    platform = meta.get("platform", {})
    typer.echo(
        f"Open Targets API: {platform.get('api_version') or 'unknown'} | "
        f"data release: {platform.get('data_version') or 'unknown'}"
    )
    typer.echo(f"Disease mappings checked: {len(rows)}")
    typer.echo(f"ERROR: {len(errors)} | WARNING: {len(warnings)} | PASS: {len(rows)-len(errors)-len(warnings)}")
    for r in rows:
        typer.echo(
            f"{r['cancer_id']}: {r['configured_disease_id']} => "
            f"{r.get('resolved_disease_id')} | {r.get('returned_disease_name')} "
            f"({r.get('resolution_method')}) [{r['status']}]"
        )
    if errors:
        raise typer.Exit(code=2)


@app.command("fetch-opentargets")
def fetch_opentargets(
    root: Path | None = typer.Option(None, help="Project root"),
    source_release: str | None = typer.Option(None, help="Open Targets data release label, e.g. 26.09"),
    endpoint: str | None = typer.Option(None, help="GraphQL endpoint"),
    page_size: int | None = typer.Option(None, help="Disease association page size"),
) -> None:
    root = _root(root)
    configure_logging(root)
    cfg = _opentargets_config(root)
    source_release = source_release or str(cfg.get("release") or "26.09")
    endpoint = endpoint or str(cfg.get("graphql_endpoint") or "https://api.platform.opentargets.org/api/v4/graphql")
    page_size = page_size or int(cfg.get("page_size") or 500)

    target_ids = root / "data/processed/target_identifiers.tsv"
    if not target_ids.exists():
        raise typer.BadParameter("Missing target_identifiers.tsv. Complete Milestone 1 first.")
    pairs_path = root / "data/input/cancer_target_pairs_wave1.tsv"
    pairs = load_wave1_pairs(pairs_path)

    rows, provenance, meta = collect_opentargets_wave1(
        root=root, source_release=source_release, endpoint=endpoint, page_size=page_size
    )
    qc = opentargets_qc(rows, expected_n=len(pairs))
    mappings = load_opentargets_disease_mappings(root / "config/cancer_contexts.yaml")
    qc.extend(opentargets_disease_mapping_qc(meta.get("diseases", {}), mappings))

    processed = root / "data/processed"
    qc_dir = root / "outputs/qc"
    reports = root / "outputs/reports"
    processed.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)

    row_dicts = [x.model_dump(mode="json") for x in rows]
    prov_dicts = [x.model_dump(mode="json") for x in provenance]
    qc_dicts = [x.model_dump(mode="json") for x in qc]

    assoc_tsv = processed / "opentargets_associations.tsv"
    assoc_parquet = processed / "opentargets_associations.parquet"
    ot_prov_tsv = processed / "opentargets_provenance.tsv"
    ot_prov_parquet = processed / "opentargets_provenance.parquet"
    write_tsv(assoc_tsv, row_dicts, list(row_dicts[0].keys()))
    write_tsv(ot_prov_tsv, prov_dicts, list(prov_dicts[0].keys()))
    pd.DataFrame(row_dicts).to_parquet(assoc_parquet, index=False)
    pd.DataFrame(prov_dicts).to_parquet(ot_prov_parquet, index=False)

    # Append to canonical provenance while preserving Milestone 1 records.
    canonical_prov_tsv = processed / "provenance.tsv"
    canonical_prov_parquet = processed / "provenance.parquet"
    existing = []
    if canonical_prov_tsv.exists():
        existing = pd.read_csv(canonical_prov_tsv, sep="\t", dtype=str).fillna("").to_dict("records")
    merged = existing + prov_dicts
    if merged:
        write_tsv(canonical_prov_tsv, merged, list(merged[0].keys()))
        pd.DataFrame(merged).to_parquet(canonical_prov_parquet, index=False)

    write_json(qc_dir / "opentargets_qc.json", qc_dicts)
    write_tsv(qc_dir / "opentargets_qc.tsv", qc_dicts, list(qc_dicts[0].keys()))
    write_json(qc_dir / "opentargets_disease_validation.json", meta)

    errors = sum(1 for x in qc if x.severity == "ERROR")
    warnings = sum(1 for x in qc if x.severity == "WARNING")
    manifest = create_opentargets_manifest(
        root=root, pairs_path=pairs_path, target_identifiers_path=target_ids,
        source_release=source_release, endpoint=endpoint, associations_count=len(rows),
        qc_errors=errors, qc_warnings=warnings,
        outputs=[str(assoc_tsv), str(assoc_parquet), str(ot_prov_tsv), str(qc_dir / "opentargets_qc.json")],
    )
    write_json(reports / "run_manifest_opentargets.json", manifest)

    found = sum(1 for x in rows if x.association_found)
    typer.echo(f"Open Targets Wave 1 rows: {len(rows)}")
    typer.echo(f"Direct associations found: {found}/{len(rows)}")
    typer.echo(f"QC ERROR: {errors} | WARNING: {warnings} | TOTAL: {len(qc)}")
    if errors:
        raise typer.Exit(code=2)


@app.command("qc-opentargets")
def qc_opentargets_cmd(root: Path | None = typer.Option(None, help="Project root")) -> None:
    root = _root(root)
    configure_logging(root)
    qc_path = root / "outputs/qc/opentargets_qc.json"
    if not qc_path.exists():
        raise typer.BadParameter("No Open Targets QC output. Run fetch-opentargets first.")
    rows = json.loads(qc_path.read_text(encoding="utf-8"))
    errors = [x for x in rows if x["severity"] == "ERROR"]
    warnings = [x for x in rows if x["severity"] == "WARNING"]
    info = [x for x in rows if x["severity"] == "INFO"]
    typer.echo(f"ERROR: {len(errors)} | WARNING: {len(warnings)} | INFO: {len(info)} | TOTAL: {len(rows)}")
    if errors:
        raise typer.Exit(code=2)


@app.command("sync-depmap")
def sync_depmap_cmd(
    root: Path | None = typer.Option(None, help="Project root"),
    release: str | None = typer.Option(None, help="Pinned DepMap release, e.g. 26Q1"),
    depmap_dir: Path | None = typer.Option(None, help="Override local DepMap cache directory"),
    offline: bool = typer.Option(False, "--offline", help="Do not contact the DepMap metadata endpoint"),
) -> None:
    """Prepare/check a pinned DepMap release and record SHA256 provenance."""
    root = _root(root)
    configure_logging(root)
    cfg = _depmap_config(root)
    release = release or str(cfg.get("release") or "26Q1")
    result = sync_depmap(
        root=root,
        release=release,
        depmap_dir=depmap_dir,
        check_catalog=not offline,
    )

    processed = root / "data/processed"
    reports = root / "outputs/reports"
    processed.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    rows = result["files"]
    write_json(reports / "depmap_sync_manifest.json", result)
    write_tsv(processed / "depmap_sync_inventory.tsv", rows, list(rows[0].keys()))

    typer.echo(f"DepMap pinned release: {release}")
    if result.get("latest_catalog_release"):
        suffix = " (newer release exists; pin is unchanged)" if result.get("newer_release_available") else ""
        typer.echo(f"Latest release visible in catalog: {result['latest_catalog_release']}{suffix}")
    elif result.get("catalog_error"):
        typer.echo(f"Catalog: unavailable ({result['catalog_error']})")
    else:
        typer.echo("Catalog: not checked")

    for row in rows:
        mark = "✓" if row["status"] == "ready" else "!"
        local = row.get("local_filename") or row["expected_names"]
        details = f" | {row['header_check']}" if row.get("header_check") else ""
        typer.echo(f"{mark} {row['logical_name']}: {local} -> {row['status']}{details}")
        if row.get("note") and row["status"] != "ready":
            typer.echo(f"  {row['note']}")

    ready_n = sum(1 for x in rows if x["status"] == "ready")
    typer.echo(f"Ready files: {ready_n}/{len(rows)}")
    typer.echo(f"Manifest: {reports / 'depmap_sync_manifest.json'}")
    if result["ready"]:
        typer.echo("DepMap cache is ready. Next: mcl validate-depmap-inputs")
    else:
        typer.echo("DepMap cache is incomplete. Download only the missing files through the official DepMap portal, then rerun sync-depmap.")
        raise typer.Exit(code=2)


def _depmap_config(root: Path) -> dict:
    cfg = yaml.safe_load((root / "config/source_versions.yaml").read_text(encoding="utf-8")) or {}
    return cfg.get("depmap", {})


@app.command("validate-depmap-inputs")
def validate_depmap_inputs_cmd(
    root: Path | None = typer.Option(None, help="Project root"),
    release: str | None = typer.Option(None, help="Pinned DepMap release, e.g. 26Q1"),
    depmap_dir: Path | None = typer.Option(None, help="Directory containing official DepMap release files"),
) -> None:
    root = _root(root)
    configure_logging(root)
    cfg = _depmap_config(root)
    release = release or str(cfg.get("release") or "26Q1")
    try:
        inventory, qc, meta = validate_depmap_inputs(root, release, depmap_dir)
    except FileNotFoundError as exc:
        typer.echo(str(exc))
        raise typer.Exit(code=2)
    qc_dir = root / "outputs/qc"
    reports = root / "outputs/reports"
    processed = root / "data/processed"
    qc_dir.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    processed.mkdir(parents=True, exist_ok=True)
    write_tsv(processed / "depmap_input_inventory.tsv", inventory, list(inventory[0].keys()))
    write_json(qc_dir / "depmap_input_qc.json", [x.model_dump(mode="json") for x in qc])
    write_tsv(qc_dir / "depmap_input_qc.tsv", [x.model_dump(mode="json") for x in qc], list(qc[0].model_dump().keys()))
    write_json(qc_dir / "depmap_input_meta.json", meta)
    errors = [x for x in qc if x.severity == "ERROR"]
    warnings = [x for x in qc if x.severity == "WARNING"]
    typer.echo(f"DepMap release: {release}")
    typer.echo(f"Files checked: {len(inventory)}")
    typer.echo(f"ERROR: {len(errors)} | WARNING: {len(warnings)} | INFO: {len(qc)-len(errors)-len(warnings)}")
    for row in qc:
        if row.severity != "INFO" or row.check in {"context_rule_feasibility", "context_rule_diagnostics"}:
            entity = f" {row.entity_id}" if row.entity_id else ""
            details = f" | observed: {row.observed}" if row.observed else ""
            typer.echo(f"[{row.severity}] {row.check}{entity}: {row.message}{details}")
    if errors:
        raise typer.Exit(code=2)


@app.command("analyze-depmap")
def analyze_depmap_cmd(
    root: Path | None = typer.Option(None, help="Project root"),
    release: str | None = typer.Option(None, help="Pinned DepMap release, e.g. 26Q1"),
    depmap_dir: Path | None = typer.Option(None, help="Directory containing official DepMap release files"),
) -> None:
    root = _root(root)
    configure_logging(root)
    cfg = _depmap_config(root)
    release = release or str(cfg.get("release") or "26Q1")

    # Fail early on schema/source problems before scientific calculations.
    inventory, input_qc, _ = validate_depmap_inputs(root, release, depmap_dir)
    input_errors = [x for x in input_qc if x.severity == "ERROR"]
    if input_errors:
        for row in input_errors:
            typer.echo(f"[ERROR] {row.check}: {row.message}")
        raise typer.Exit(code=2)

    rows, audit, provenance, inventory, meta = analyze_depmap_wave1(root, release, depmap_dir)
    expected_n = len(load_wave1_pairs(root / "data/input/cancer_target_pairs_wave1.tsv"))
    qc = depmap_qc(rows, audit, expected_pairs_n=expected_n)

    processed = root / "data/processed"
    qc_dir = root / "outputs/qc"
    reports = root / "outputs/reports"
    processed.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)

    row_dicts = [x.model_dump(mode="json") for x in rows]
    audit_dicts = [x.model_dump(mode="json") for x in audit]
    prov_dicts = [x.model_dump(mode="json") for x in provenance]
    qc_dicts = [x.model_dump(mode="json") for x in qc]

    evidence_tsv = processed / "depmap_evidence.tsv"
    evidence_parquet = processed / "depmap_evidence.parquet"
    audit_tsv = processed / "depmap_context_audit.tsv"
    audit_parquet = processed / "depmap_context_audit.parquet"
    prov_tsv = processed / "depmap_provenance.tsv"
    prov_parquet = processed / "depmap_provenance.parquet"
    inventory_tsv = processed / "depmap_input_inventory.tsv"

    write_tsv(evidence_tsv, row_dicts, list(row_dicts[0].keys()))
    write_tsv(audit_tsv, audit_dicts, list(audit_dicts[0].keys()))
    write_tsv(prov_tsv, prov_dicts, list(prov_dicts[0].keys()))
    write_tsv(inventory_tsv, inventory, list(inventory[0].keys()))
    pd.DataFrame(row_dicts).to_parquet(evidence_parquet, index=False)
    pd.DataFrame(audit_dicts).to_parquet(audit_parquet, index=False)
    pd.DataFrame(prov_dicts).to_parquet(prov_parquet, index=False)

    # Preserve provenance from earlier milestones.
    canonical_prov = processed / "provenance.tsv"
    existing = []
    if canonical_prov.exists():
        existing = pd.read_csv(canonical_prov, sep="\t", dtype=str).fillna("").to_dict("records")
    merged = existing + prov_dicts
    if merged:
        write_tsv(canonical_prov, merged, list(merged[0].keys()))
        pd.DataFrame(merged).to_parquet(processed / "provenance.parquet", index=False)

    write_json(qc_dir / "depmap_qc.json", qc_dicts)
    write_tsv(qc_dir / "depmap_qc.tsv", qc_dicts, list(qc_dicts[0].keys()))
    write_json(qc_dir / "depmap_analysis_meta.json", meta)

    errors = sum(1 for x in qc if x.severity == "ERROR")
    warnings = sum(1 for x in qc if x.severity == "WARNING")
    manifest = create_depmap_manifest(
        root=root,
        release=release,
        files=meta["files"],
        pairs_path=root / "data/input/cancer_target_pairs_wave1.tsv",
        evidence_n=len(rows),
        context_audit_n=len(audit),
        qc_errors=errors,
        qc_warnings=warnings,
        outputs=[str(evidence_tsv), str(evidence_parquet), str(audit_tsv), str(qc_dir / "depmap_qc.json")],
    )
    write_json(reports / "run_manifest_depmap.json", manifest)

    groups = {}
    for x in audit:
        groups.setdefault(x.cancer_id, {"context": 0, "comparator": 0, "excluded": 0})[x.assigned_group] += 1
    typer.echo(f"DepMap release: {release}")
    typer.echo(f"Wave 1 evidence rows: {len(rows)}")
    for cancer_id, counts in groups.items():
        typer.echo(f"{cancer_id}: context={counts['context']} | comparator={counts['comparator']} | excluded={counts['excluded']}")
    typer.echo(f"QC ERROR: {errors} | WARNING: {warnings} | TOTAL: {len(qc)}")
    if errors:
        raise typer.Exit(code=2)


@app.command("qc-depmap")
def qc_depmap_cmd(root: Path | None = typer.Option(None, help="Project root")) -> None:
    root = _root(root)
    configure_logging(root)
    qc_path = root / "outputs/qc/depmap_qc.json"
    if not qc_path.exists():
        raise typer.BadParameter("No DepMap QC output. Run analyze-depmap first.")
    rows = json.loads(qc_path.read_text(encoding="utf-8"))
    errors = [x for x in rows if x["severity"] == "ERROR"]
    warnings = [x for x in rows if x["severity"] == "WARNING"]
    info = [x for x in rows if x["severity"] == "INFO"]
    typer.echo(f"ERROR: {len(errors)} | WARNING: {len(warnings)} | INFO: {len(info)} | TOTAL: {len(rows)}")
    for row in errors + warnings:
        typer.echo(f"[{row['severity']}] {row['check']} {row.get('entity_id') or ''}: {row['message']}")
    if errors:
        raise typer.Exit(code=2)


@app.command("analyze-depmap-kras-sensitivity")
def analyze_depmap_kras_sensitivity_cmd(
    root: Path | None = typer.Option(None, help="Project root"),
    release: str | None = typer.Option(None, help="Pinned DepMap release, e.g. 26Q1"),
    depmap_dir: Path | None = typer.Option(None, help="Directory containing official DepMap release files"),
) -> None:
    """M3.1: split KRAS exact-variant contexts into WT-proxy and other-KRAS comparators."""
    root = _root(root)
    configure_logging(root)
    cfg = _depmap_config(root)
    release = release or str(cfg.get("release") or "26Q1")

    # M3.1 uses exactly the same pinned files and schema gate as primary M3.
    _, input_qc, _ = validate_depmap_inputs(root, release, depmap_dir)
    input_errors = [x for x in input_qc if x.severity == "ERROR"]
    if input_errors:
        for row in input_errors:
            typer.echo(f"[ERROR] {row.check}: {row.message}")
        raise typer.Exit(code=2)

    rows, audit, meta = analyze_depmap_kras_sensitivity(root, release, depmap_dir)
    pair_table = pd.read_csv(root / "data/input/cancer_target_pairs_wave1.tsv", sep="\t", dtype=str).fillna("")
    eligible_contexts = set(meta["eligible_contexts"])
    eligible_pairs_n = int(pair_table["Cancer_ID"].isin(eligible_contexts).sum())
    comparison_types_n = len(meta["comparison_types"])
    expected_rows_n = eligible_pairs_n * comparison_types_n
    qc = depmap_kras_sensitivity_qc(rows, audit, expected_rows_n=expected_rows_n)

    processed = root / "data/processed"
    qc_dir = root / "outputs/qc"
    reports = root / "outputs/reports"
    processed.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)

    row_dicts = [x.model_dump(mode="json") for x in rows]
    audit_dicts = [x.model_dump(mode="json") for x in audit]
    qc_dicts = [x.model_dump(mode="json") for x in qc]

    evidence_tsv = processed / "depmap_kras_sensitivity.tsv"
    evidence_parquet = processed / "depmap_kras_sensitivity.parquet"
    audit_tsv = processed / "depmap_kras_sensitivity_audit.tsv"
    audit_parquet = processed / "depmap_kras_sensitivity_audit.parquet"
    write_tsv(evidence_tsv, row_dicts, list(row_dicts[0].keys()))
    write_tsv(audit_tsv, audit_dicts, list(audit_dicts[0].keys()))
    parquet_written = True
    try:
        pd.DataFrame(row_dicts).to_parquet(evidence_parquet, index=False)
        pd.DataFrame(audit_dicts).to_parquet(audit_parquet, index=False)
    except ImportError:
        parquet_written = False
    meta["parquet_written"] = parquet_written
    write_json(qc_dir / "depmap_kras_sensitivity_qc.json", qc_dicts)
    write_tsv(qc_dir / "depmap_kras_sensitivity_qc.tsv", qc_dicts, list(qc_dicts[0].keys()))
    write_json(qc_dir / "depmap_kras_sensitivity_meta.json", meta)

    errors = [x for x in qc if x.severity == "ERROR"]
    warnings = [x for x in qc if x.severity == "WARNING"]
    manifest = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "code_version": "0.3.5",
        "milestone": "3.1 DepMap KRAS comparator sensitivity",
        "depmap_release": release,
        "input_pairs": str(root / "data/input/cancer_target_pairs_wave1.tsv"),
        "input_pairs_sha256": sha256_file(root / "data/input/cancer_target_pairs_wave1.tsv"),
        "eligible_cancer_target_pairs_n": eligible_pairs_n,
        "comparison_types_n": comparison_types_n,
        "evidence_rows_n": len(rows),
        "audit_rows_n": len(audit),
        "qc_errors_n": len(errors),
        "qc_warnings_n": len(warnings),
        "qc_result": "PASS" if not errors else "FAIL",
        "outputs": [str(evidence_tsv), str(audit_tsv), str(qc_dir / "depmap_kras_sensitivity_qc.json")],
        "analysis_meta": meta,
    }
    write_json(reports / "run_manifest_depmap_kras_sensitivity.json", manifest)

    typer.echo(f"DepMap release: {release}")
    typer.echo(f"M3.1 KRAS sensitivity rows: {len(rows)}")
    for cancer_id in sorted({x.cancer_id for x in audit}):
        subset = [x for x in audit if x.cancer_id == cancer_id]
        counts: dict[str, int] = {}
        for x in subset:
            counts[x.assigned_group] = counts.get(x.assigned_group, 0) + 1
        typer.echo(
            f"{cancer_id}: target={counts.get('target_variant', 0)} | "
            f"KRAS-WT-proxy={counts.get('kras_wildtype_proxy', 0)} | "
            f"other-KRAS-driver/hotspot={counts.get('other_kras_driver_hotspot', 0)} | "
            f"excluded={counts.get('excluded', 0)}"
        )
    typer.echo(f"QC ERROR: {len(errors)} | WARNING: {len(warnings)} | TOTAL: {len(qc)}")
    if errors:
        raise typer.Exit(code=2)


@app.command("qc-depmap-kras-sensitivity")
def qc_depmap_kras_sensitivity_cmd(root: Path | None = typer.Option(None, help="Project root")) -> None:
    root = _root(root)
    configure_logging(root)
    qc_path = root / "outputs/qc/depmap_kras_sensitivity_qc.json"
    if not qc_path.exists():
        raise typer.BadParameter("No M3.1 KRAS sensitivity QC output. Run analyze-depmap-kras-sensitivity first.")
    rows = json.loads(qc_path.read_text(encoding="utf-8"))
    errors = [x for x in rows if x["severity"] == "ERROR"]
    warnings = [x for x in rows if x["severity"] == "WARNING"]
    info = [x for x in rows if x["severity"] == "INFO"]
    typer.echo(f"ERROR: {len(errors)} | WARNING: {len(warnings)} | INFO: {len(info)} | TOTAL: {len(rows)}")
    for row in errors + warnings:
        typer.echo(f"[{row['severity']}] {row['check']} {row.get('entity_id') or ''}: {row['message']}")
    if errors:
        raise typer.Exit(code=2)



if __name__ == "__main__":
    app()


@app.command("analyze-depmap-genome-wide")
def analyze_depmap_genome_wide_cmd(
    cancer_id: str = typer.Option(..., "--cancer-id", help="Cancer context ID, e.g. CANCER-001"),
    comparison: str = typer.Option(
        "primary",
        "--comparison",
        help="primary | kras-wt | other-kras",
    ),
    root: Path | None = typer.Option(None, help="Project root"),
    release: str | None = typer.Option(None, help="Pinned DepMap release, e.g. 26Q1"),
    depmap_dir: Path | None = typer.Option(None, help="Directory containing official DepMap release files"),
    top_n_heatmap: int = typer.Option(30, min=5, max=100, help="Genes shown in the explorer heatmap"),
) -> None:
    """Run genome-wide CRISPR dependency analysis and build an offline HTML explorer."""
    root = _root(root)
    configure_logging(root)
    cfg = _depmap_config(root)
    release = release or str(cfg.get("release") or "26Q1")
    comparison = comparison.strip().lower()

    # Reuse the same strict source/schema gate as M3 and M3.1.
    _, input_qc, _ = validate_depmap_inputs(root, release, depmap_dir)
    input_errors = [x for x in input_qc if x.severity == "ERROR"]
    if input_errors:
        for row in input_errors:
            typer.echo(f"[ERROR] {row.check}: {row.message}")
        raise typer.Exit(code=2)

    try:
        results, cohort, gene_effect, dependency, meta = analyze_depmap_genomewide(
            root=root,
            release=release,
            cancer_id=cancer_id,
            comparison=comparison,
            depmap_dir=depmap_dir,
        )
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc

    qc = depmap_genomewide_qc(results, meta)
    key = f"{cancer_id}__{comparison}".replace("/", "_").replace("\\", "_")
    processed = root / "data/processed/depmap_genomewide"
    qc_dir = root / "outputs/qc"
    reports = root / "outputs/reports"
    processed.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)

    genes_tsv = processed / f"{key}_genes.tsv"
    genes_parquet = processed / f"{key}_genes.parquet"
    cohort_tsv = processed / f"{key}_cohort.tsv"
    cohort_parquet = processed / f"{key}_cohort.parquet"
    ge_parquet = processed / f"{key}_gene_effect_matrix.parquet"
    dp_parquet = processed / f"{key}_dependency_matrix.parquet"
    qc_json = qc_dir / f"depmap_genomewide_{key}_qc.json"
    qc_tsv = qc_dir / f"depmap_genomewide_{key}_qc.tsv"
    meta_json = reports / f"depmap_genomewide_{key}_meta.json"
    explorer_html = reports / f"depmap_explorer_{key}.html"

    result_records = results.astype(object).where(pd.notna(results), None).to_dict("records")
    cohort_records = cohort.astype(object).where(pd.notna(cohort), None).to_dict("records")
    write_tsv(genes_tsv, result_records, list(results.columns))
    results.to_parquet(genes_parquet, index=False)
    if cohort_records:
        write_tsv(cohort_tsv, cohort_records, list(cohort.columns))
    cohort.to_parquet(cohort_parquet, index=False)
    gene_effect.reset_index().to_parquet(ge_parquet, index=False)
    dependency.reset_index().to_parquet(dp_parquet, index=False)

    qc_dicts = [x.model_dump(mode="json") for x in qc]
    write_json(qc_json, qc_dicts)
    write_tsv(qc_tsv, qc_dicts, list(qc_dicts[0].keys()))
    meta = {
        **meta,
        "outputs": {
            "genes_tsv": str(genes_tsv),
            "genes_parquet": str(genes_parquet),
            "cohort_tsv": str(cohort_tsv),
            "gene_effect_matrix_parquet": str(ge_parquet),
            "dependency_matrix_parquet": str(dp_parquet),
            "explorer_html": str(explorer_html),
        },
    }
    write_json(meta_json, meta)
    build_depmap_explorer_html(
        results=results,
        cohort=cohort,
        gene_effect=gene_effect,
        meta=meta,
        output_path=explorer_html,
        top_n_heatmap=top_n_heatmap,
    )

    errors = [x for x in qc if x.severity == "ERROR"]
    warnings = [x for x in qc if x.severity == "WARNING"]
    fdr_n = int(results["fdr_0_05"].fillna(False).astype(bool).sum())
    broad_n = int(results["broad_dependency_warning"].fillna(False).astype(bool).sum())
    typer.echo(f"DepMap release: {release}")
    typer.echo(f"Genome-wide context: {cancer_id} | comparison={comparison}")
    typer.echo(
        f"Models: context={meta['context_models_n']} | comparator={meta['comparator_models_n']}"
    )
    typer.echo(f"Genes analyzed: {len(results)}")
    typer.echo(f"FDR q<0.05: {fdr_n} | broad-dependency flags: {broad_n}")
    typer.echo(f"Explorer: {explorer_html}")
    typer.echo(f"QC ERROR: {len(errors)} | WARNING: {len(warnings)} | TOTAL: {len(qc)}")
    if errors:
        raise typer.Exit(code=2)


@app.command("qc-depmap-genome-wide")
def qc_depmap_genome_wide_cmd(
    cancer_id: str = typer.Option(..., "--cancer-id", help="Cancer context ID"),
    comparison: str = typer.Option("primary", "--comparison", help="primary | kras-wt | other-kras"),
    root: Path | None = typer.Option(None, help="Project root"),
) -> None:
    root = _root(root)
    configure_logging(root)
    key = f"{cancer_id}__{comparison.strip().lower()}".replace("/", "_").replace("\\", "_")
    qc_path = root / "outputs/qc" / f"depmap_genomewide_{key}_qc.json"
    if not qc_path.exists():
        raise typer.BadParameter(
            "No genome-wide QC output for this context/comparison. Run analyze-depmap-genome-wide first."
        )
    rows = json.loads(qc_path.read_text(encoding="utf-8"))
    errors = [x for x in rows if x["severity"] == "ERROR"]
    warnings = [x for x in rows if x["severity"] == "WARNING"]
    info = [x for x in rows if x["severity"] == "INFO"]
    typer.echo(
        f"ERROR: {len(errors)} | WARNING: {len(warnings)} | INFO: {len(info)} | TOTAL: {len(rows)}"
    )
    for row in errors + warnings:
        typer.echo(
            f"[{row['severity']}] {row['check']} {row.get('entity_id') or ''}: {row['message']}"
        )
