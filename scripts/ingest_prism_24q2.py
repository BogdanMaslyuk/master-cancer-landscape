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


def _norm_col_name(value: object) -> str:
    return str(value).strip().lower().replace(".", "_").replace("-", "_").replace(" ", "_")


def _col(frame: pd.DataFrame, *names: str) -> str | None:
    lookup = {_norm_col_name(c): str(c) for c in frame.columns}
    for name in names:
        key = _norm_col_name(name)
        if key in lookup:
            return lookup[key]
    return None


def _col_contains(frame: pd.DataFrame, *tokens: str) -> str | None:
    wanted = tuple(_norm_col_name(x) for x in tokens)
    for column in frame.columns:
        key = _norm_col_name(column)
        if all(token in key for token in wanted):
            return str(column)
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
        ach_cols = [c for c in frame.columns if str(c).upper().startswith("ACH-")]
        if not ach_cols:
            raise SystemExit(
                f"Unsupported PRISM response schema. Columns: {list(frame.columns)[:20]}"
            )
        id_candidates = [c for c in frame.columns if c not in ach_cols]
        if not id_candidates:
            raise SystemExit("PRISM wide matrix has ACH-* columns but no compound identifier column")
        id_col = id_candidates[0]
        print(f"Detected wide PRISM matrix: compound ID column = {id_col!r}; model columns = {len(ach_cols)}")
        melted = frame[[id_col, *ach_cols]].melt(
            id_vars=[id_col], var_name="model_id", value_name="LFC"
        )
        melted = melted.rename(columns={id_col: "profile_id"})
        melted["profile_id"] = melted["profile_id"].map(_text)
        melted["LFC"] = pd.to_numeric(melted["LFC"], errors="coerce")
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


def _infer_profile_column(treatment: pd.DataFrame, response_profiles: set[str]) -> tuple[str | None, int]:
    """Resolve the compound ID column by actual overlap with the PRISM response matrix.

    Extended PRISM exports have changed header names across releases. Matching values
    is safer than hard-coding one header and also surfaces malformed metadata clearly.
    """
    explicit = _col(
        treatment,
        "profile_id", "compound_id", "broad_id", "id", "ids", "pert_id", "perturbation_id",
    )
    if explicit:
        values = set(treatment[explicit].map(_text))
        overlap = len(response_profiles & values)
        if overlap:
            return explicit, overlap

    best_column: str | None = None
    best_overlap = 0
    for column in treatment.columns:
        values = set(treatment[column].map(_text))
        overlap = len(response_profiles & values)
        if overlap > best_overlap:
            best_column = str(column)
            best_overlap = overlap
    return best_column, best_overlap


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

    response_profiles = set(response["profile_id"].map(_text)) - {""}
    profile_t, profile_overlap = _infer_profile_column(treatment, response_profiles)
    if not profile_t or profile_overlap == 0:
        raise SystemExit(
            "Could not map PRISM compound metadata to the response matrix. "
            f"Response compound IDs: {len(response_profiles)}; metadata columns: {list(treatment.columns)}"
        )

    print(
        f"Detected compound metadata ID column = {profile_t!r}; "
        f"matched {profile_overlap}/{len(response_profiles)} matrix compound IDs"
    )

    broad_t = _col(treatment, "broad_id", "broad id") or _col_contains(treatment, "broad", "id")
    name_t = _col(
        treatment,
        "name", "compound_name", "drug_name", "pert_iname", "compound", "drug",
    ) or _col_contains(treatment, "compound", "name")
    moa_t = _col(
        treatment, "moa", "mechanism_of_action", "mechanism", "mechanism of action"
    ) or _col_contains(treatment, "mechanism")
    target_t = _col(
        treatment, "target", "targets", "gene_targets", "gene_target", "target_gene"
    ) or _col_contains(treatment, "target")
    smiles_t = _col(
        treatment, "smiles", "canonical_smiles", "canonical smiles"
    ) or _col_contains(treatment, "smiles")

    print(
        "Detected metadata fields: "
        f"name={name_t!r}, broad_id={broad_t!r}, moa={moa_t!r}, target={target_t!r}, smiles={smiles_t!r}"
    )

    meta = treatment.copy()
    meta["_profile"] = meta[profile_t].map(_text)
    meta = meta[meta["_profile"] != ""].copy()
    meta["_compound_id"] = meta["_profile"].map(_compound_id)
    if broad_t:
        meta["_broad_id"] = meta[broad_t].map(_text)
    else:
        # The extended primary matrix itself commonly uses BRD identifiers as the
        # treatment profile ID. Preserve that external ID when it is recognizable.
        meta["_broad_id"] = meta["_profile"].map(
            lambda x: x if str(x).upper().startswith(("BRD:", "BRD-")) else ""
        )

    keep_cols = ["_profile", "_compound_id", "_broad_id"] + [
        c for c in (name_t, moa_t, target_t, smiles_t) if c
    ]
    keep_cols = list(dict.fromkeys(keep_cols))
    meta = meta[keep_cols].drop_duplicates("_profile", keep="first")

    response = response[response["profile_id"].map(_text) != ""].copy()
    response = response.merge(meta, left_on="profile_id", right_on="_profile", how="left")
    response["compound_id"] = response["profile_id"].map(_compound_id)

    mapped_profiles = set(meta["_profile"])
    unmatched_profiles = sorted(response_profiles - mapped_profiles)
    print(f"Unmatched matrix compound IDs after metadata join: {len(unmatched_profiles)}")
    if unmatched_profiles:
        print(f"  examples: {unmatched_profiles[:5]}")

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
                    {"prism_profile_id": profile, "broad_id": broad or None},
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
