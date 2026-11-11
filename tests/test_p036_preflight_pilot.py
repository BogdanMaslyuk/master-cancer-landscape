"""P0.3.6: no fabricated readiness, stable panels, exact-archive integrity."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str):
    path = ROOT/"scripts"/name
    spec = importlib.util.spec_from_file_location(path.stem,path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


BASE = load_script("audit_target_space_v2_baseline.py")
PREFLIGHT = load_script("preflight_tls_v2_real_data.py")
PANEL = load_script("build_tls_v2_pilot_panel.py")


def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture(tmp_path: Path):
    root = tmp_path/"science"
    root.mkdir(parents=True)
    for relative in BASE.INPUTS.values():
        path = root/relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"fixture data")
    genes = ["GENE" + str(i).zfill(2) for i in range(13)]
    pd.DataFrame([
        {"target_gene": gene, "priority_status": "priority_for_in_vitro"}
        for gene in genes
    ]).to_parquet(root/BASE.INPUTS["candidate_v2"])
    compounds = pd.DataFrame([
        {"own_compound_id": "PYZ-" + str(i).zfill(3)} for i in range(1,66)
    ])
    compounds.to_csv(root/BASE.INPUTS["own_compounds"], sep="\t", index=False)
    baseline = BASE.audit(root)
    baseline_path = tmp_path/"audit"/"baseline_audit.json"
    write_json(baseline_path, baseline)
    registry = tmp_path/"p02"/"protein_target_universe.tsv"
    registry.parent.mkdir()
    all_genes = genes + ["EXTRA" + str(i).zfill(2) for i in range(10)]
    pd.DataFrame([
        {"target_gene":g,"uniprot_accession":"P" + str(i+1).zfill(5),
         "mapping_status":"unique_reviewed_human","organism_id":"9606"}
        for i,g in enumerate(all_genes)
    ]).to_csv(registry, sep="\t", index=False)
    original = tmp_path/"p02"/"official_reviewed_human.tsv"
    original.write_text("This is a synthetic source fixture, not official UniProt.")
    manifest = tmp_path/"p02"/"protein_target_universe_manifest.json"
    write_json(manifest, {
        "contract":"mcl-protein-target-universe-v2-p0",
        "target_genes_n":len(all_genes),
        "uniprot_snapshot_included":True,
        "uniprot_release":"synthetic-fixture-do-not-publish",
        "uniprot_source_url":"https://rest.uniprot.org/uniprotkb/stream",
        "source_files":[{"path":str(original),"sha256":sha(original)}],
    })
    return root,baseline_path,registry,manifest,genes


def test_p01_p02_preflight_requires_sources_and_assay_review(tmp_path):
    root,baseline,registry,manifest,genes=fixture(tmp_path)
    result=PREFLIGHT.inspect(root,baseline,registry,manifest)
    assert not result["ready_for_full_real_data_pilot"]
    assert result["checks"]["p01_core13"]
    assert result["checks"]["p01_pyz65"]
    assert result["checks"]["p01_input_files_unchanged"]
    assert result["checks"]["p02_core13_mapped"]
    assert not result["checks"]["p03_source_manifest_available"]
    assert "p03_sources_not_yet_pinned" in result["blocking_or_informational_issues"]


def test_source_tampering_blocks_preflight(tmp_path):
    root,baseline,registry,manifest,_=fixture(tmp_path)
    (root/BASE.INPUTS["ligand_v1"]).write_bytes(b"not the fixture hash")
    result=PREFLIGHT.inspect(root,baseline,registry,manifest)
    assert result["checks"]["p01_input_files_unchanged"] is False
    assert not result["ready_for_full_real_data_pilot"]


def test_strict_source_manifest_and_approval_checks(tmp_path):
    root,baseline,registry,manifest,_=fixture(tmp_path)
    chembl=tmp_path/"chembl_joined.tsv"
    binding=tmp_path/"bindingdb.tsv"
    chembl.write_text("synthetic fixture only")
    binding.write_text("synthetic fixture only")
    sources=tmp_path/"source_manifest.json"
    write_json(sources,{
        "contract":"mcl-tls-v2-pinned-extracts-v1",
        "sources":[
            {"source":"ChEMBL","release":"test","path":str(chembl),"sha256":sha(chembl),
             "url":"https://ftp.ebi.ac.uk/pub/databases/chembl/fixture","license":"fixture terms"},
            {"source":"BindingDB","release":"test","path":str(binding),"sha256":sha(binding),
             "url":"https://bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp","license":"fixture terms"},
        ],
    })
    review=tmp_path/"assays.tsv"
    review.write_text(
        "source\tsource_release\tassay_id\tmethod_class\treviewer\tevidence_url\tverification_status\n"
        "ChEMBL\ttest\tA1\tbinding_direct\ttest-reviewer\thttps://example.org/fixture\tapproved\n"
    )
    result=PREFLIGHT.inspect(root,baseline,registry,manifest,
                             source_manifest=sources,review_file=review)
    assert result["ready_for_full_real_data_pilot"]
    # This is only a file-integrity gate, not an independent validation of
    # whether the scientist really checked the primary paper.
    assert "assay_review_is_human_assertion_not_independent_source_verification" in result["blocking_or_informational_issues"]
    binding.write_text("changed")
    result=PREFLIGHT.inspect(root,baseline,registry,manifest,
                             source_manifest=sources,review_file=review)
    assert not result["ready_for_full_real_data_pilot"]
    assert not result["checks"]["p03_source_files_hash"]


def test_pilot_selects_13_core_plus_7_expanded_without_claiming_activity(tmp_path):
    _,baseline,registry,_,core=fixture(tmp_path)
    rows,report=PANEL.create_panel(baseline,registry)
    assert len(rows)==20
    assert sum(r["scope"]=="core" for r in rows)==13
    assert sum(r["scope"]=="expanded" for r in rows)==7
    assert set(core)=={r["target_gene"] for r in rows if r["scope"]=="core"}
    assert all(not r["verified_ligand_count"] for r in rows)
    assert report["evidence_not_verified"]
    assert report["expanded_selection_basis"]=="deterministic_gene_order_no_ligand_evidence"


def test_pilot_rejects_unmapped_core_and_invalid_extras(tmp_path):
    _,baseline,registry,_,core=fixture(tmp_path)
    frame=pd.read_csv(registry,sep="\t")
    frame=frame[frame.target_gene.ne(core[0])]
    frame.to_csv(registry,sep="\t",index=False)
    with pytest.raises(ValueError,match="not uniquely mapped"):
        PANEL.create_panel(baseline,registry)
    _,baseline,registry,_,_=fixture(tmp_path/"another")
    with pytest.raises(ValueError,match="Explicit extra genes"):
        PANEL.create_panel(baseline,registry,extra_genes=("EXTRA01",),additional_n=7)


def test_pilot_ranks_only_observed_source_rows_not_verified_ligands(tmp_path):
    _,baseline,registry,_,_=fixture(tmp_path)
    raw=tmp_path/"unreviewed.parquet"
    pd.DataFrame([
        {"target_gene":"EXTRA09","source_record_id":str(i)} for i in range(6)
    ]+[{"target_gene":"EXTRA02","source_record_id":"e1"}]).to_parquet(raw)
    rows,report=PANEL.create_panel(baseline,registry,measured_rows=raw,additional_n=5)
    extra=[r for r in rows if r["scope"]=="expanded"]
    assert extra[0]["target_gene"]=="EXTRA09"
    assert extra[0]["source_observed_measurement_rows"]=="6"
    assert all(x["verified_ligand_count"]=="" for x in extra)
    assert report["expanded_selection_basis"]=="observed_source_rows_only_not_curated"


def test_unofficial_source_urls_fail_preflight(tmp_path):
    root,baseline,registry,manifest,_=fixture(tmp_path)
    document=json.loads(manifest.read_text(encoding="utf-8"))
    document["uniprot_source_url"]="https://example.com/not-official"
    manifest.write_text(json.dumps(document),encoding="utf-8")
    result=PREFLIGHT.inspect(root,baseline,registry,manifest)
    assert not result["checks"]["p02_official_source_declared"]
    assert not result["ready_for_full_real_data_pilot"]


def test_65_compound_count_cannot_hide_replaced_pyz_identifier(tmp_path):
    root,baseline_file,registry,manifest,_=fixture(tmp_path)
    compounds=root/BASE.INPUTS["own_compounds"]
    frame=pd.read_csv(compounds,sep="\t",dtype=str)
    frame.loc[0,"own_compound_id"]="PYZ-999"
    frame.to_csv(compounds,sep="\t",index=False)
    # Recalculate the baseline to show count=65 by itself does not prove
    # preservation of the 65 canonical PYZ identities.
    baseline_file.write_text(json.dumps(BASE.audit(root)))
    result=PREFLIGHT.inspect(root,baseline_file,registry,manifest)
    assert result["checks"]["p01_pyz65"]
    assert not result["checks"]["p01_exact_pyz001_065_ids"]
    assert not result["ready_for_full_real_data_pilot"]


def test_ligand_baseline_does_not_require_full_depmap_context(tmp_path):
    """A missing P1.2 cell-context artifact must not stop P0.3 ligand assembly."""
    root, baseline_path, registry, manifest, _ = fixture(tmp_path)
    (root / BASE.INPUTS["full_depmap_context"]).unlink()
    (root / BASE.INPUTS["ligand_manifest_v1"]).unlink()
    updated = BASE.audit(root)
    assert updated["ready_for_extended_build"] is True
    assert updated["ready_for_depmap_cell_integration"] is False
    assert updated["deferred_capabilities"]["full_depmap_context_present"] is False
    baseline_path.write_text(json.dumps(updated), encoding="utf-8")
    result = PREFLIGHT.inspect(root, baseline_path, registry, manifest)
    assert result["checks"]["p01_input_files_unchanged"] is True
    assert any("p01_deferred_optional_source:full_depmap_context" == entry
               for entry in result["blocking_or_informational_issues"])
    assert result["checks"]["p01_core13"] is True
    assert result["checks"]["p01_pyz65"] is True
    # Real-data preflight is still blocked by missing ChEMBL/BindingDB
    # and missing documented assay approvals.
    assert result["ready_for_full_real_data_pilot"] is False


def test_pilot_output_under_audit_report_parent_is_allowed(tmp_path, monkeypatch):
    """Regression: MCL_TLS_V2_CHECK/pilot beside baseline_audit.json is safe."""
    _, baseline, registry, _, _ = fixture(tmp_path)
    output = baseline.parent / "pilot"
    monkeypatch.setattr("sys.argv", [
        "build_tls_v2_pilot_panel.py",
        "--baseline", str(baseline),
        "--registry", str(registry),
        "--additional-n", "7",
        "--output-dir", str(output),
    ])
    PANEL.main()
    panel = pd.read_csv(output / "pilot_panel.tsv", sep="\t")
    assert len(panel) == 20
    assert baseline.is_file()
    assert (output / "pilot_panel_manifest.json").is_file()


def test_pilot_output_inside_original_science_root_is_forbidden(tmp_path, monkeypatch):
    root, baseline, registry, _, _ = fixture(tmp_path)
    output = root / "unsafe_pilot_results"
    monkeypatch.setattr("sys.argv", [
        "build_tls_v2_pilot_panel.py",
        "--baseline", str(baseline),
        "--registry", str(registry),
        "--output-dir", str(output),
    ])
    with pytest.raises(SystemExit, match="outside the original scientific data"):
        PANEL.main()
    assert not output.exists()
