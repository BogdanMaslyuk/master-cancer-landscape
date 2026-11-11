#!/usr/bin/env python3
"""P0.3.6 scientific pilot gate. Fail closed: CI success != real-data PASS.

Verifies the archived index, pinned pilot, unchanged legacy files, two
independent content-identical builds and a source-linked manual review record.
Produces a reviewable report; NEVER switches a production API or index.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from validate_expanded_target_index import validate as validate_index  # noqa: E402


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for part in iter(lambda:stream.read(1024*1024),b""):
            h.update(part)
    return h.hexdigest()


def read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def dry_run_pointer_rollback() -> bool:
    """An isolated file-pointer dry-run, NOT a switch of MCL production state."""
    with tempfile.TemporaryDirectory(prefix="mcl-tls-p036-") as folder:
        root=Path(folder)
        pointer=root/"current.json"
        pointer.write_text(json.dumps({"version":"legacy-v1"}),encoding="utf-8")
        before=sha256(pointer)
        temporary=root/"candidate.tmp"
        temporary.write_text(json.dumps({"version":"candidate-v2"}),encoding="utf-8")
        os.replace(temporary,pointer)
        if read(pointer).get("version")!="candidate-v2":
            return False
        restored=root/"restore.tmp"
        restored.write_text(json.dumps({"version":"legacy-v1"}),encoding="utf-8")
        os.replace(restored,pointer)
        return sha256(pointer)==before


def legacy_unchanged(baseline_file: Path, science_root: Path) -> tuple[bool,list[str]]:
    baseline=read(baseline_file)
    problems=[]
    if baseline.get("contract")!="mcl-expanded-target-space-baseline-audit-v1":
        return False,["invalid_baseline_contract"]
    if Path(baseline.get("science_root","")).resolve()!=science_root.resolve():
        problems.append("science_root_mismatch")
    for name,entry in baseline.get("files",{}).items():
        path=science_root/str(entry.get("relative_path",""))
        if not path.is_file() or sha256(path)!=entry.get("sha256"):
            problems.append("legacy_changed_or_missing:"+name)
    return not problems,problems


def review_gate(file: Path | None, build_id: str, actual_ids: set[tuple[str,str,str]]) -> tuple[bool,dict]:
    if file is None or not file.is_file():
        return False,{"reason":"manual_scientific_review_not_supplied"}
    try:
        record=read(file)
        required={"build_id","reviewer","review_date","decision",
                  "sampling_method","reviewed_cases","issues_unresolved_n"}
        if not required.issubset(record):
            return False,{"reason":"manual_review_missing_fields"}
        if record["build_id"]!=build_id or record["decision"]!="approved":
            return False,{"reason":"manual_review_wrong_build_or_not_approved"}
        if not str(record["reviewer"]).strip() or not str(record["sampling_method"]).strip():
            return False,{"reason":"manual_review_missing_attribution"}
        if int(record["issues_unresolved_n"])!=0:
            return False,{"reason":"manual_review_unresolved_issues"}
        cases=record["reviewed_cases"]
        if not isinstance(cases,list) or len(cases)<min(5,len(actual_ids)):
            return False,{"reason":"too_few_reviewed_primary_measurements"}
        reviewed=set()
        for sample in cases:
            key=(str(sample.get("source","")),str(sample.get("source_record_id","")),
                 str(sample.get("endpoint","")))
            if key not in actual_ids or key in reviewed:
                return False,{"reason":"manual_review_contains_unmatched_or_repeated_case"}
            if (sample.get("decision")!="approved" or
                    not str(sample.get("evidence_url","")).startswith("https://")):
                return False,{"reason":"manual_review_missing_source_or_decision"}
            reviewed.add(key)
        return True,{
            "reviewer":record["reviewer"],
            "review_date":record["review_date"],
            "sampled_verified_cases_n":len(reviewed),
            "sample_strategy":record["sampling_method"],
            "review_file_sha256":sha256(file),
            "note":"This verifies reviewer assertions and matching IDs; it cannot independently certify the primary experimental methods.",
        }
    except (ValueError,TypeError,KeyError,OSError) as exc:
        return False,{"reason":"invalid_manual_review:"+type(exc).__name__}


def compare_stages(first: Path, second: Path | None) -> tuple[bool,dict]:
    if second is None or not second.is_dir() or first.resolve()==second.resolve():
        return False,{"reason":"independent_second_build_missing"}
    other=read(second/"build_manifest.json")
    original=read(first/"build_manifest.json")
    keys=("input_normalized_sha256","input_registry_sha256","assay_review_sha256",
          "rdkit_version","fingerprint","unique_target_structure_pairs_n")
    differences=[key for key in keys if other.get(key)!=original.get(key)]
    a={entry["path"]:entry["sha256"] for entry in original.get("artifacts",[])}
    b={entry["path"]:entry["sha256"] for entry in other.get("artifacts",[])}
    if a!=b:
        differences.append("artifact_content_sha256")
    return not differences,{
        "matching_content":not differences,
        "differences":differences,
        "first_build_id":original.get("build_id"),
        "second_build_id":other.get("build_id"),
        "note":"Service timestamps are excluded from scientific content comparison.",
    }


def inspect(
    preflight_path: Path, pilot_panel: Path, stage: Path, science_root: Path,
    baseline_path: Path, *, second_stage: Path | None=None,
    manual_review: Path | None=None,
) -> dict:
    checks:dict[str,bool]={}
    notes=[]
    evidence={}
    preflight=read(preflight_path)
    checks["preflight_real_sources_and_review"] = (
        preflight.get("contract")=="mcl-tls-v2-p036-real-data-preflight-v1"
        and preflight.get("ready_for_full_real_data_pilot") is True
        and all(x is True for x in preflight.get("checks",{}).values())
    )
    panel=pd.read_csv(pilot_panel,sep="\t",dtype=str,keep_default_na=False)
    mandatory={"target_gene","uniprot_accession","scope","direct_assay_evidence_status"}
    if not mandatory.issubset(panel.columns):
        raise ValueError("Pilot panel has no required UniProt/scope evidence columns")
    checks["pilot_13_core"] = len(panel[panel["scope"]=="core"])==13
    checks["pilot_5_to_10_expanded"] = 5<=len(panel[panel["scope"]=="expanded"])<=10
    checks["pilot_unique_targets"] = (
        panel["target_gene"].is_unique and panel["uniprot_accession"].is_unique
    )
    core=set(panel.loc[panel.scope=="core","target_gene"])
    extended=set(panel.loc[panel.scope=="expanded","target_gene"])
    stage_manifest=read(stage/"build_manifest.json")
    validation=None
    try:
        validation=validate_index(stage)
        checks["index_internal_consistency"]=True
    except (ValueError,OSError,AssertionError,KeyError) as exc:
        checks["index_internal_consistency"]=False
        notes.append("index_validation_failed:"+type(exc).__name__+":"+str(exc))
    checks["same_pinned_pilot_ingest"] = False
    raw_ingest=stage.parent.parent.parent if False else None
    # The index stage records the ingest *content SHA*, but the original
    # ingest_manifest and normalized table are outside the index stage.
    # Compare a complete scientific chain only if an explicit ingest stage
    # path is passed in a future signed release record.
    checks["assay_source_linkage_present"]=stage_manifest.get("approved_assays_n",0)>0
    curated=stage/"measurements_curated.parquet"
    counts=Counter()
    by_source=Counter()
    by_endpoint=Counter()
    independent=set()
    actual_ids=set()
    target_counts=Counter()
    examples=[]
    if curated.is_file():
        for batch in pq.ParquetFile(curated).iter_batches(batch_size=8192):
            for row in batch.to_pylist():
                tier=row.get("curated_quality_tier") or ""
                counts[tier]+=1
                if tier in {"A","B"} and row.get("independence_status")=="primary":
                    key=(row.get("source") or "",row.get("source_record_id") or "",
                         row.get("endpoint") or "")
                    if key in actual_ids:
                        notes.append("duplicate_primary_measurement_id:"+str(key))
                    actual_ids.add(key)
                    independent.add(key)
                    by_source[row["source"]]+=1
                    by_endpoint[row["endpoint"]]+=1
                    target_counts[row["target_gene"]]+=1
                    if len(examples)<10:
                        examples.append({
                            "target_gene":row["target_gene"],
                            "source":row["source"],
                            "source_record_id":row["source_record_id"],
                            "endpoint":row["endpoint"],
                            "relation":row.get("relation"),
                            "value_nm":row.get("numeric_value_nm"),
                            "assay_id":row.get("assay_id"),
                            "method":row.get("review_method"),
                            "original_smiles":row.get("original_smiles"),
                        })
    checks["curated_verified_A_B_nonempty"]=len(independent)>0
    checks["curated_primary_ids_unique"]=not any(x.startswith("duplicate_primary_measurement_id") for x in notes)
    checks["scientific_count_consistency"]=(
        int(stage_manifest.get("verified_measurements_n",-1))==len(independent)
        and int(stage_manifest.get("measurement_rows_n",-1))==sum(counts.values())
    )
    checks["pilot_only_target_scope"]=set(target_counts).issubset(core|extended)
    if set(target_counts)-(core|extended):
        notes.append("validated_records_outside_pilot:"+
                     ",".join(sorted(set(target_counts)-(core|extended))))
    checks["legacy_archive_unchanged"],legacy_issues=legacy_unchanged(
        baseline_path,science_root
    )
    notes.extend(legacy_issues)
    checks["independent_repeat_build_identical"], repeat=compare_stages(
        stage,second_stage
    )
    checks["scratch_pointer_rollback"]=dry_run_pointer_rollback()
    checks["manual_scientific_signoff"], manual=review_gate(
        manual_review,stage_manifest.get("build_id",""),actual_ids
    )
    if not checks["manual_scientific_signoff"]:
        notes.append(str(manual.get("reason")))
    evidence.update({
        "index_integrity":validation,
        "build_id":stage_manifest.get("build_id"),
        "curated_tier_counts":dict(counts),
        "verified_primary_measurements_n":len(independent),
        "verified_by_source":dict(by_source),
        "verified_by_endpoint":dict(by_endpoint),
        "verified_by_target":dict(sorted(target_counts.items())),
        "pilot_coverage": [{
            "target_gene":row["target_gene"],"scope":row["scope"],
            "uniprot_accession":row["uniprot_accession"],
            "verified_primary_measurements_in_pilot":target_counts.get(row["target_gene"],0),
            "interpretation":"Only in the pinned and reviewed processed source; not universal absence"
        } for row in panel.to_dict("records")],
        "reviewer":manual,
        "repeatability":repeat,
        "examples":examples,
        "raw_input_index_manifest_sha256":sha256(stage/"build_manifest.json"),
        "baseline_sha256":sha256(baseline_path),
        "pilot_panel_sha256":sha256(pilot_panel),
    })
    # The old site is never modified by this validator, regardless of status.
    accepted=bool(checks) and all(checks.values())
    return {
        "contract":"mcl-tls-v2-p036-scientific-pilot-acceptance-v1",
        "generated_at_utc":datetime.now(timezone.utc).isoformat(),
        "status":"PASS_FOR_SCIENTIFIC_REVIEW" if accepted else "BLOCKED",
        "all_pilot_gates_met":accepted,
        "checks":checks,"evidence":evidence,
        "issues":sorted(set(notes+[k for k,v in checks.items() if not v])),
        "scope_note":"A successful pilot does not authorize P1 API switching or merge PR #3.",
        "limitations":[
            "Full-dataset build and stratified source review are separate from this pilot.",
            "Human assay approval must be genuine; matching review IDs is not independent evidence verification.",
            "Rollback check uses a scratch version pointer; no active site version is modified.",
        ],
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight",type=Path,required=True)
    parser.add_argument("--pilot-panel",type=Path,required=True)
    parser.add_argument("--stage",type=Path,required=True)
    parser.add_argument("--science-root",type=Path,required=True)
    parser.add_argument("--baseline",type=Path,required=True)
    parser.add_argument("--second-stage",type=Path)
    parser.add_argument("--manual-review",type=Path)
    parser.add_argument("--report",type=Path,required=True)
    parser.add_argument("--strict",action="store_true")
    args=parser.parse_args()
    if args.report.resolve().is_relative_to(args.science_root.resolve()):
        raise SystemExit("Report output must be outside the archived scientific source")
    report=inspect(args.preflight,args.pilot_panel,args.stage,args.science_root,
                   args.baseline,second_stage=args.second_stage,
                   manual_review=args.manual_review)
    args.report.parent.mkdir(parents=True,exist_ok=True)
    temp=args.report.with_suffix(args.report.suffix+".tmp")
    temp.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    temp.replace(args.report)
    print(f"P0.3.6 pilot gate: {report['status']}")
    print(f"Verified independent measurement rows: {report['evidence']['verified_primary_measurements_n']}")
    if args.strict and not report["all_pilot_gates_met"]:
        raise SystemExit("Scientific pilot is blocked; no active MCL files modified.")


if __name__=="__main__":
    main()
