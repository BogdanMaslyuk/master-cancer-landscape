"""P0.3.6 gate tests on synthetic fixtures. NEVER real scientific evidence."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def script(name):
    path=ROOT/"scripts"/name
    spec=importlib.util.spec_from_file_location(path.stem,path)
    instance=importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(instance)
    return instance


AUDIT=script("audit_target_space_v2_baseline.py")
BUILDER=script("build_expanded_target_index.py")
ACCEPT=script("validate_tls_v2_scientific_pilot.py")


def hash_of(p:Path):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def synthetic_pilot(tmp_path:Path):
    root=tmp_path/"legacy"
    root.mkdir()
    for relative in AUDIT.INPUTS.values():
        p=root/relative
        p.parent.mkdir(parents=True,exist_ok=True)
        p.write_bytes(b"synthetic legacy fixture")
    genes=["GENE"+str(i).zfill(2) for i in range(13)]
    pd.DataFrame([
        {"target_gene":g,"priority_status":"priority_for_in_vitro"} for g in genes
    ]).to_parquet(root/AUDIT.INPUTS["candidate_v2"])
    pd.DataFrame([
        {"own_compound_id":"PYZ-"+str(i).zfill(3)} for i in range(1,66)
    ]).to_csv(root/AUDIT.INPUTS["own_compounds"],sep="\t",index=False)
    base_file=tmp_path/"baseline_audit.json"
    base_file.write_text(json.dumps(AUDIT.audit(root)))

    panel=tmp_path/"pilot_panel.tsv"
    entries = [
        {"target_gene":g,"uniprot_accession":"P"+str(i).zfill(5),
         "scope":"core","direct_assay_evidence_status":"not_evaluated"}
        for i,g in enumerate(genes,1)
    ] + [
        {"target_gene":g,"uniprot_accession":"P"+str(i).zfill(5),
         "scope":"expanded","direct_assay_evidence_status":"not_evaluated"}
        for i,g in enumerate(["CDK4","EGFR","EXTRA03","EXTRA04","EXTRA05"],14)
    ]
    pd.DataFrame(entries).to_csv(panel,sep="\t",index=False)
    source=tmp_path/"p03"/"staging"/"fixture"
    source.mkdir(parents=True)
    rows=[
        {
            "source":"ChEMBL","source_release":"fixture",
            "source_record_id":"101","target_gene":"CDK4",
            "uniprot_accession":"P11802","original_smiles":"CCO",
            "endpoint":"Ki","assay_id":"A1",
            "quality_tier":"candidate_A","relation":"=","numeric_value_nm":"350",
            "measurement_status":"exact_eligible",
            "original_chembl_activity_id":"","doi":"",
        },
        {
            "source":"ChEMBL","source_release":"fixture",
            "source_record_id":"102","target_gene":"EGFR",
            "uniprot_accession":"P00533","original_smiles":"CCN",
            "endpoint":"IC50","assay_id":"A2",
            "quality_tier":"candidate_B","relation":"=","numeric_value_nm":"500",
            "measurement_status":"exact_eligible",
            "original_chembl_activity_id":"","doi":"",
        },
    ]
    normalized=source/"measurements_normalized.parquet"
    pd.DataFrame(rows).to_parquet(normalized,index=False)
    registry=tmp_path/"protein_target_universe.tsv"
    pd.DataFrame([
        {"target_gene":"CDK4","uniprot_accession":"P11802"},
        {"target_gene":"EGFR","uniprot_accession":"P00533"},
    ]).to_csv(registry,sep="\t",index=False)
    ingest={
        "contract":"mcl-tls-v2-normalized-staging-v1",
        "status":"staging_not_published",
        "measurement_rows":len(rows),
        "normalized_sha256":hash_of(normalized),
        "registry_sha256":hash_of(registry),
        "selection_mode":"reviewed_pilot_panel",
        "pilot_panel_sha256":hash_of(panel),
    }
    (source/"ingest_manifest.json").write_text(json.dumps(ingest))
    review=tmp_path/"assay_review.tsv"
    review.write_text(
        "source\tsource_release\tassay_id\tmethod_class\treviewer\tevidence_url\tverification_status\n"
        "ChEMBL\tfixture\tA1\tbinding_direct\tfixture-reviewer\thttps://example.org/1\tapproved\n"
        "ChEMBL\tfixture\tA2\tbiochemical_function\tfixture-reviewer\thttps://example.org/2\tapproved\n"
    )
    a=BUILDER.build(source,registry,review,tmp_path/"runs_one")
    b=BUILDER.build(source,registry,review,tmp_path/"runs_two")
    stage1=tmp_path/"runs_one"/"staging"/a["build_id"]
    stage2=tmp_path/"runs_two"/"staging"/b["build_id"]

    preflight=tmp_path/"preflight_report.json"
    preflight.write_text(json.dumps({
        "contract":"mcl-tls-v2-p036-real-data-preflight-v1",
        "ready_for_full_real_data_pilot":True,
        "checks":{"fixture_only":True},
    }))

    manual=tmp_path/"manual_scientific_review.json"
    manual.write_text(json.dumps({
        "build_id":a["build_id"],
        "reviewer":"SYNTHETIC-TEST-NOT-A-REAL-SCIENTIST",
        "review_date":"2026-11-11",
        "decision":"approved",
        "sampling_method":"synthetic test case, not genuine review",
        "issues_unresolved_n":0,
        "reviewed_cases":[
            {"source":"ChEMBL","source_record_id":"101","endpoint":"Ki",
             "evidence_url":"https://example.org/fake-paper-1","decision":"approved"},
            {"source":"ChEMBL","source_record_id":"102","endpoint":"IC50",
             "evidence_url":"https://example.org/fake-paper-2","decision":"approved"},
        ],
    }))
    return preflight,panel,stage1,stage2,root,base_file,source,manual


def test_missing_real_review_and_repeat_build_block_pilot(tmp_path):
    preflight,panel,stage1,stage2,root,baseline,ingest,manual=synthetic_pilot(tmp_path)
    result=ACCEPT.inspect(preflight,panel,stage1,root,baseline,
                          ingest_stage=ingest,second_stage=stage2)
    assert result["status"]=="BLOCKED"
    assert result["checks"]["independent_repeat_build_identical"]
    assert result["checks"]["same_pinned_pilot_ingest"]
    assert not result["checks"]["manual_scientific_signoff"]
    assert result["evidence"]["verified_primary_measurements_n"]==2
    assert result["checks"]["legacy_archive_unchanged"]


def test_fully_matching_synthetic_gate_only_proves_software_behavior(tmp_path):
    preflight,panel,stage1,stage2,root,baseline,ingest,manual=synthetic_pilot(tmp_path)
    result=ACCEPT.inspect(preflight,panel,stage1,root,baseline,
                          ingest_stage=ingest,second_stage=stage2,
                          manual_review=manual)
    assert result["all_pilot_gates_met"]
    assert result["status"]=="PASS_FOR_SCIENTIFIC_REVIEW"
    assert result["checks"]["scratch_pointer_rollback"]
    assert result["evidence"]["verified_by_target"]=={"CDK4":1,"EGFR":1}
    assert len(result["evidence"]["pilot_coverage"])==18
    assert not (root/"current.json").exists()


def test_mismatched_stage_or_modified_legacy_blocks_pilot(tmp_path):
    preflight,panel,stage1,stage2,root,baseline,ingest,manual=synthetic_pilot(tmp_path)
    raw=ingest/"measurements_normalized.parquet"
    raw.write_bytes(b"tampered")
    (root/AUDIT.INPUTS["ligand_v1"]).write_bytes(b"changed")
    result=ACCEPT.inspect(preflight,panel,stage1,root,baseline,
                          ingest_stage=ingest,second_stage=stage2,
                          manual_review=manual)
    assert result["status"]=="BLOCKED"
    assert not result["checks"]["same_pinned_pilot_ingest"]
    assert not result["checks"]["legacy_archive_unchanged"]


def test_missing_repeat_stage_always_blocks(tmp_path):
    preflight,panel,stage1,stage2,root,baseline,ingest,manual=synthetic_pilot(tmp_path)
    result=ACCEPT.inspect(preflight,panel,stage1,root,baseline,
                          ingest_stage=ingest,manual_review=manual)
    assert not result["checks"]["independent_repeat_build_identical"]
    assert not result["all_pilot_gates_met"]
