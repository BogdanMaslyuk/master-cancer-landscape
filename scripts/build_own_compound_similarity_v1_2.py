from __future__ import annotations

"""Robust entry point for Own Compound Structural Similarity v1.2.

Scientific calculation remains the v1 implementation. This wrapper adds two
compatibility/QC protections discovered on the first real structure-enriched run:

1. metadata columns are coerced to object/text-compatible dtype before fallback
   assignment, avoiding pandas incompatible-dtype warnings;
2. every standardized RDKit molecule is reparsed/sanitized, its ring information
   is initialized explicitly, and Morgan fingerprint generation is pre-validated.

A comparator that still cannot produce a valid Morgan fingerprint is rejected by
structure QC rather than aborting the entire run. Similarity parameters and the
biological comparator universe are unchanged.
"""

from typing import Any

import pandas as pd
from rdkit import Chem

import build_own_compound_similarity as impl
from build_own_compound_similarity_v1_1 import _build_overlay


_ORIGINAL_LOAD_INPUTS = impl._load_inputs
_ORIGINAL_STANDARDIZE = impl._standardize_smiles
_ORIGINAL_MORGAN = impl.MORGAN
_TEXT_METADATA_COLUMNS = (
    "preferred_name",
    "canonical_smiles",
    "inchikey",
    "pubchem_cid",
    "chembl_id",
    "broad_id",
    "gdsc_id",
)


def _text_object(series: pd.Series) -> pd.Series:
    return series.astype("object").where(series.notna(), "")


def _load_inputs_v1_2() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    own, compounds, comparator_pairs = _ORIGINAL_LOAD_INPUTS()
    compounds = compounds.copy()
    comparator_pairs = comparator_pairs.copy()
    for frame in (compounds, comparator_pairs):
        for column in _TEXT_METADATA_COLUMNS:
            if column in frame.columns:
                frame[column] = _text_object(frame[column])
    return own, compounds, comparator_pairs


def _standardize_smiles_v1_2(
    smiles: Any,
) -> tuple[Any | None, str | None, str | None, str | None]:
    mol, canonical, inchikey, error = _ORIGINAL_STANDARDIZE(smiles)
    if error or mol is None or not canonical:
        return mol, canonical, inchikey, error

    try:
        # Reparse the standardized representation with full sanitization. Some
        # molecules returned by MolStandardize/FragmentParent can otherwise have
        # an uninitialized RingInfo cache in specific RDKit builds on Windows.
        repaired = Chem.MolFromSmiles(canonical, sanitize=True)
        if repaired is None:
            return None, None, None, "post_standardization_reparse_failed"
        repaired.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(repaired)
        Chem.GetSymmSSSR(repaired)
        Chem.AssignStereochemistry(repaired, cleanIt=True, force=True)

        canonical_repaired = Chem.MolToSmiles(
            repaired, canonical=True, isomericSmiles=True
        )
        inchikey_repaired = Chem.MolToInchiKey(repaired)

        # Pre-validate the exact fingerprint operation used by the scientific
        # implementation. Failure is converted into structure QC, not a global
        # pipeline crash.
        _ORIGINAL_MORGAN.GetFingerprint(repaired)
    except Exception as exc:
        return None, None, None, f"fingerprint_preparation_failed:{type(exc).__name__}"

    return repaired, canonical_repaired, inchikey_repaired, None


def main() -> None:
    impl.COMPOUNDS = _build_overlay()
    impl._load_inputs = _load_inputs_v1_2
    impl._standardize_smiles = _standardize_smiles_v1_2
    impl.main()


if __name__ == "__main__":
    main()
