from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "data" / "raw" / "pharmacology" / "prism_24q2"
DEFAULT_OUTPUT = ROOT / "data" / "raw" / "pharmacology" / "normalized" / "prism_24q2"
SOURCE = "PRISM"
RELEASE = "PRISM Repurposing Public 24Q2"


def _find(directory: Path, tokens: tuple[str, ...]) -> Path | None:
    for path in sorted(directory.glob("*")):
        key = path.name.lower().replace("-", "_")
        if all(token.lower().replace("-", "_") in key for token in tokens):
            return path
    return None


def _col(frame: pd.DataFrame, *names: str) -> str | None:
    lookup = {
        str(c).strip().lower().replace(".", "_").replace("-", "_"): str(c)
        for c in frame.columns
    }
    for name in names:
        key = name.strip().lower().replace(".", "_").replace("-", "_")
        if key in lookup:
            return lookup[key]
    return None


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _targets(value: object) -> list[str]:
    raw = _text(value)
    if not raw:
        return []
    parts = re.split(r"[,;|/]", raw)
    output: list[str] = []
    for part in parts:
        gene = part.strip().upper()
        if not gene or gene in {"NA", "NAN", "NONE", "UNKNOWN"}:
            continue
        # Repurposing Hub target annotations are normally compact gene symbols.
        # Verbose mechanism text remains in the separate MOA/action field.
        if " " in gene or len(gene) > 30:
            continue
        if gene not in output:
            output.append(gene)
    return output


def _load_csv(path: Path) -> pd.DataFrame:
    print(f"Reading {path.name} ...")
    return pd.read_csv(path, low_memory=False)


def _resolve_files(directory: Path) -> tuple[Path, Path]:
    lfc = _find(directory, ("lfc",))
    treatment = _find(directory, ("treatment", "info"))
    if lfc and treatment:
        return lfc, treatment
    lfc = lfc or _find(directory, ("data", "matrix"))
    treatment = treatment or _find(directory, ("compound", "list"))
    if not lfc or not treatment:
        raise SystemExit(
            "Could not resolve PRISM response and treatment files. Expected LFC + Treatment_Info "
            "or Extended_Primary_Data_Matrix + Extended_Primary_Compound_List."
        )
    return lfc, treatment


def _long_lfc(frame: pd.DataFrame) -> pd.DataFrame:
    row_col = _col(frame, "row_id", "depmap_id", "model_id")
    profile_col = _col(frame, "profile_id", "broad_id", "compound_id")
    lfc_col = _col(frame, "LFC", "logfold_change", "log_fold_change", "value")
    if not row_col or not profile_col or not lfc_col:
        # Fallback for an extended wide matrix: ACH-* columns are models and the
        # remaining identifier column identifies the treatment/compound profile.
        ach_cols = [c for c in frame.columns if str(c).upper().startswith("ACH-")]
        if not ach_cols:
            raise SystemExit(
                f"Unsupported PRISM response schema. Columns: {list(frame.columns)[:20]}"
            )
        id_candidates = [c for c in frame.columns if c not in ach_cols]
        id_col = id_candidates[0]
        melted = frame[[id_col, *ach_cols]].melt(
            id_vars=[id_col], var_name="model_id", value_name="LFC"
        )
        melted = melted.rename(columns={id_col: "profile_id"})
        melted["profile_id"] = melted["profile_id"].map(_text)
        return melted.dropna(subset=["LFC"])

    out = pd.DataFrame()
    row_values = frame[row_col].astype(str)
    if row_values.str.contains("::", regex=False).any():
        out["model_id"] = row_values.str.split("::").str[0]
    else:
        out["model_id"] = row_values
    out["profile_id"] = frame[profile_col].map(_text)
    out["LFC"] = pd.to_numeric(frame[lfc_col], errors="coerce")
    pass_col = _col(frame, "PASS", "passed_qc", "pass")
    if pass_col:
        truth = frame[pass_col].astype(str).str.lower().isin(
            {"true", "1", "yes", "pass", "passed"}
        )
        out = out[truth].copy()
    return out.dropna(subset=["LFC"])


def _compound_id(profile_id: object) -> str:
    profile = _text(profile_id)
    return f"PRISM24Q2:{profile}" if profile else ""


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Normalize PRISM Repurposing Public 24Q2 for MCL Pharmacology Layer v1."
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    response_path, treatment_path = _resolve_files(input_dir)
    raw = _load_csv(response_path)
    treatment = _load_csv(treatment_path)
    response = _long_lfc(raw)

    profile_t = _col(treatment, "profile_id", "compound_id", "broad_id")
    broad_t = _col(treatment, "broad_id")
    name_t = _col(treatment, "name", "compound_name", "drug_name")
    moa_t = _col(treatment, "moa", "mechanism_of_action", "mechanism")
    target_t = _col(treatment, "target", "targets", "gene_targets")
    smiles_t = _col(treatment, "smiles", "canonical_smiles")
    if not profile_t:
        raise SystemExit(
            "PRISM treatment metadata has no profile_id/compound_id/broad_id column"
        )

    meta = treatment.copy()
    meta["_profile"] = meta[profile_t].map(_text)
    meta = meta[meta["_profile"] != ""].copy()
    meta["_compound_id"] = meta["_profile"].map(_compound_id)
    meta["_broad_id"] = meta[broad_t].map(_text) if broad_t else ""
    keep_cols = ["_profile", "_compound_id", "_broad_id"] + [
        c for c in (name_t, moa_t, target_t, smiles_t) if c
    ]
    meta = meta[keep_cols].drop_duplicates("_profile", keep="first")

    response = response[response["profile_id"].map(_text) != ""].copy()
    response = response.merge(meta, left_on="profile_id", right_on="_profile", how="left")
    response["compound_id"] = response["profile_id"].map(_compound_id)

    compounds = pd.DataFrame(
        {
            "compound_id": meta["_compound_id"],
            "preferred_name": meta[name_t].map(_text) if name_t else meta["_profile"],
            "canonical_smiles": meta[smiles_t].map(_text) if smiles_t else "",
            "inchikey": "",
            "pubchem_cid": "",
            "chembl_id": "",
            "broad_id": meta["_broad_id"],
            "gdsc_id": "",
            "source_ids_json": [
                json.dumps(
                    {
                        "prism_profile_id": profile,
                        "broad_id": broad or None,
                    },
                    ensure_ascii=False,
                )
                for profile, broad in zip(meta["_profile"], meta["_broad_id"])
            ],
        }
    ).drop_duplicates("compound_id", keep="first")

    aggregated = response.groupby(["model_id", "compound_id"], as_index=False)["LFC"].mean()
    aggregated["observation_id"] = [
        f"PRISM24Q2-{i:08d}" for i in range(1, len(aggregated) + 1)
    ]
    responses = pd.DataFrame(
        {
            "observation_id": aggregated["observation_id"],
            "model_id": aggregated["model_id"].astype(str),
            "compound_id": aggregated["compound_id"].astype(str),
            "source": SOURCE,
            "source_release": RELEASE,
            "source_assay_id": "Repurposing Public 24Q2",
            "assay_type": "PRISM pooled cell-line viability",
            "endpoint": "LFC",
            "value": aggregated["LFC"],
            "unit": "log2 fold-change vs DMSO",
            "auc": None,
            "ic50": None,
            "ec50": None,
            "gi50": None,
            "viability": None,
            "dose": 2.5,
            "dose_unit": "uM",
            "exposure_time_h": 120.0,
            "replicate_n": None,
            "quality_flag": "PASS rows aggregated where PASS was available",
            "publication": "PRISM Repurposing Public 24Q2",
        }
    )

    target_rows: list[dict[str, object]] = []
    if target_t:
        unique_meta = meta.drop_duplicates("_compound_id", keep="first")
        for _, row in unique_meta.iterrows():
            compound_id = str(row["_compound_id"])
            moa = _text(row.get(moa_t)) if moa_t else ""
            for gene in _targets(row.get(target_t)):
                target_rows.append(
                    {
                        "evidence_id": f"PRISM24Q2-TGT-{len(target_rows)+1:08d}",
                        "compound_id": compound_id,
                        "target_gene": gene,
                        "action": moa,
                        "evidence_type": "curated_repurposing_hub_annotation",
                        "activity_type": "",
                        "activity_value": None,
                        "activity_unit": "",
                        "source": "Broad Repurposing Hub / PRISM",
                        "source_assay_id": "",
                        "publication": "PRISM Repurposing Public 24Q2",
                        "confidence": "curated_annotation",
                        "directness": "compound_target_annotation_not_cell_line_mechanism",
                    }
                )
    targets = pd.DataFrame(
        target_rows,
        columns=[
            "evidence_id", "compound_id", "target_gene", "action", "evidence_type",
            "activity_type", "activity_value", "activity_unit", "source", "source_assay_id",
            "publication", "confidence", "directness",
        ],
    )

    compounds.to_csv(output_dir / "compounds.tsv", sep="\t", index=False)
    responses.to_csv(output_dir / "responses.tsv", sep="\t", index=False)
    targets.to_csv(output_dir / "target_evidence.tsv", sep="\t", index=False)

    print("PRISM 24Q2 normalized for MCL Pharmacology Layer v1")
    print(f"Models: {responses['model_id'].nunique()}")
    print(f"Treatment profiles / compounds: {compounds['compound_id'].nunique()}")
    print(f"Response observations: {len(responses)}")
    print(f"Compound-target annotations: {len(targets)}")
    duplicated_broad = compounds.loc[
        compounds["broad_id"].astype(str).ne("")
        & compounds["broad_id"].astype(str).duplicated(keep=False),
        "broad_id",
    ].nunique()
    print(f"Duplicated non-empty broad_id values kept separate by profile_id: {duplicated_broad}")
    print(f"Wrote {output_dir.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
