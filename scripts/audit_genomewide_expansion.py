from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
OUTPUT_DIR = ROOT / "outputs" / "qc"


def _read_audit() -> pd.DataFrame:
    parquet = PROCESSED / "depmap_context_audit.parquet"
    tsv = PROCESSED / "depmap_context_audit.tsv"
    if parquet.exists():
        return pd.read_parquet(parquet)
    if tsv.exists():
        return pd.read_csv(tsv, sep="\t", low_memory=False)
    raise SystemExit("Missing data/processed/depmap_context_audit.tsv/parquet")


def main() -> None:
    config = yaml.safe_load((ROOT / "config" / "cancer_contexts.yaml").read_text(encoding="utf-8")) or {}
    contexts = config.get("contexts") or {}
    pathway_config = yaml.safe_load((ROOT / "config" / "pathways.yaml").read_text(encoding="utf-8")) or {}
    active = {(str(row.get("cancer_id")), str(row.get("comparison"))) for row in pathway_config.get("comparisons") or []}
    audit = _read_audit()
    required = {"cancer_id", "model_id", "assigned_group"}
    missing = required - set(audit.columns)
    if missing:
        raise SystemExit(f"Context audit is missing columns: {sorted(missing)}")

    rows: list[dict[str, object]] = []
    for cancer_id, cfg in contexts.items():
        sub = audit[audit["cancer_id"].astype(str) == str(cancer_id)].copy()
        groups = sub["assigned_group"].astype(str) if not sub.empty else pd.Series(dtype=str)
        context_n = int((groups == "context").sum())
        comparator_n = int((groups == "comparator").sum())
        excluded_n = int((groups == "excluded").sum())
        threshold = int((cfg or {}).get("minimum_context_n_warning") or 5)

        if context_n == 0:
            status = "NO_CONTEXT"
            recommendation = "Do not run: no target/context models are available."
        elif comparator_n == 0:
            status = "NO_COMPARATOR"
            recommendation = "Do not claim differential dependency. Redesign/expand the comparator before genome-wide testing."
        elif context_n < threshold:
            status = "EXPLORATORY_LOW_N"
            recommendation = (
                f"Genome-wide calculation may be run only as exploratory evidence; context n={context_n} is below the prespecified threshold {threshold}. "
                "Keep low_sample_size=true and exclude from primary pathway discovery."
            )
        elif comparator_n < threshold:
            status = "EXPLORATORY_LOW_COMPARATOR_N"
            recommendation = (
                f"Comparator n={comparator_n} is below the prespecified threshold {threshold}; keep inference exploratory."
            )
        else:
            status = "READY"
            recommendation = "Eligible for a primary context-vs-comparator genome-wide run, subject to normal QC."

        rows.append(
            {
                "cancer_id": cancer_id,
                "name": (cfg or {}).get("name") or cancer_id,
                "context_models_n": context_n,
                "comparator_models_n": comparator_n,
                "excluded_models_n": excluded_n,
                "minimum_n_threshold": threshold,
                "readiness_status": status,
                "already_in_pathway_config": any(cid == str(cancer_id) for cid, _ in active),
                "recommendation": recommendation,
            }
        )

    frame = pd.DataFrame(rows)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    tsv_path = OUTPUT_DIR / "gene_explorer_genomewide_expansion_readiness.tsv"
    json_path = OUTPUT_DIR / "gene_explorer_genomewide_expansion_readiness.json"
    frame.to_csv(tsv_path, sep="\t", index=False)
    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": "data/processed/depmap_context_audit + config/cancer_contexts.yaml",
        "principle": "No molecular context is promoted into primary genome-wide discovery without an explicit comparator and prespecified minimum sample-size check.",
        "contexts": frame.to_dict("records"),
    }
    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    print(frame[["cancer_id", "context_models_n", "comparator_models_n", "minimum_n_threshold", "readiness_status"]].to_string(index=False))
    print(f"Wrote {tsv_path.relative_to(ROOT)}")
    print(f"Wrote {json_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
