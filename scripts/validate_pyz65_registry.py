from __future__ import annotations

from pathlib import Path

import pandas as pd
from rdkit import Chem


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "config" / "own_compounds.tsv"
EXPECTED_IDS = [f"PYZ-{i:03d}" for i in range(1, 66)]
REQUIRED_COLUMNS = [
    "own_compound_id",
    "source_molecule_id",
    "preferred_name",
    "standardized_smiles",
]


def main() -> None:
    if not REGISTRY.exists():
        raise SystemExit(f"Missing approved PYZ registry: {REGISTRY}")

    frame = pd.read_csv(REGISTRY, sep="\t", dtype=str, keep_default_na=False)
    missing_columns = [column for column in REQUIRED_COLUMNS if column not in frame.columns]
    if missing_columns:
        raise SystemExit(f"PYZ registry is missing required columns: {missing_columns}")

    ids = frame["own_compound_id"].tolist()
    if ids != EXPECTED_IDS:
        raise SystemExit(
            "Approved own-compound registry must contain exactly PYZ-001...PYZ-065 in order. "
            f"Observed rows={len(ids)}; first={ids[:3]}; last={ids[-3:] if ids else []}"
        )

    for column in ("source_molecule_id", "preferred_name"):
        values = frame[column].tolist()
        if values != EXPECTED_IDS:
            raise SystemExit(f"Column {column} must match PYZ-001...PYZ-065 exactly.")

    smiles = frame["standardized_smiles"].astype(str).str.strip()
    if smiles.eq("").any():
        bad = frame.loc[smiles.eq(""), "own_compound_id"].tolist()
        raise SystemExit(f"Empty SMILES in approved PYZ registry: {bad}")
    if smiles.duplicated().any():
        duplicates = frame.loc[smiles.duplicated(keep=False), ["own_compound_id", "standardized_smiles"]]
        raise SystemExit(
            "Duplicate literal SMILES in approved PYZ registry:\n" + duplicates.to_string(index=False)
        )

    invalid: list[str] = []
    canonical: list[str] = []
    for own_id, value in zip(frame["own_compound_id"], smiles):
        mol = Chem.MolFromSmiles(value)
        if mol is None:
            invalid.append(own_id)
            canonical.append("")
        else:
            canonical.append(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True))
    if invalid:
        raise SystemExit(f"RDKit could not parse approved PYZ structures: {invalid}")

    canonical_series = pd.Series(canonical)
    if canonical_series.duplicated().any():
        dup_idx = canonical_series[canonical_series.duplicated(keep=False)].index
        duplicates = frame.loc[dup_idx, ["own_compound_id", "standardized_smiles"]]
        raise SystemExit(
            "Chemically duplicate structures detected after RDKit canonicalization:\n"
            + duplicates.to_string(index=False)
        )

    print("MCL approved own-compound registry")
    print("Registry: PYZ-001...PYZ-065")
    print(f"Rows: {len(frame)}")
    print("RDKit parse: 65/65 PASS")
    print("Canonical-structure duplicates: 0")


if __name__ == "__main__":
    main()
