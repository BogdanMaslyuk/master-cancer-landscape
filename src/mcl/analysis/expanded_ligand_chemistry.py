"""P0.3.5: strict chemical identity and compatible Morgan binary fingerprints.

Two stereoisomers can share a non-chiral fingerprint, but they must NEVER be
merged as one experimentally measured structure. Full structure identity is
computed before fragment-parent extraction.
"""
from __future__ import annotations

import hashlib
from functools import lru_cache

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize

FP_RADIUS = 2
FP_BITS = 2048
FP_CHIRALITY = False
FP_BYTES = FP_BITS // 8

# A conservative pilot definition. Other modalities are retained in QC,
# not silently treated as identical organic small molecules.
SMALL_MOLECULE_ELEMENTS = {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I", "B", "Si"}
GENERATOR = rdFingerprintGenerator.GetMorganGenerator(
    radius=FP_RADIUS, fpSize=FP_BITS, includeChirality=FP_CHIRALITY
)


def _identity(prefix: str, smiles: str) -> str:
    return hashlib.sha256((prefix + "|" + smiles).encode("utf-8")).hexdigest()[:24]


@lru_cache(maxsize=50000)
def characterize_smiles(smiles: str) -> tuple[dict[str, str], bytes | None]:
    raw = str(smiles or "").strip()
    base: dict[str, str] = {
        "exact_structure_id": "",
        "parent_connectivity_id": "",
        "canonical_smiles": "",
        "parent_smiles": "",
        "full_inchikey": "",
        "parent_inchikey": "",
        "chemistry_status": "excluded",
        "chemistry_reasons": "",
        "rdkit_version": rdBase.rdkitVersion,
    }
    if not raw:
        base["chemistry_reasons"] = "missing_structure"
        return base, None
    try:
        mol = Chem.MolFromSmiles(raw, sanitize=True)
        if mol is None:
            base["chemistry_reasons"] = "rdkit_parse_failed"
            return base, None
        # This structure, including salts and stereochemistry, is the exact
        # experimental substance annotation. Parent extraction is separate.
        full_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        full_inchikey = Chem.MolToInchiKey(mol)
        if not full_inchikey:
            base["chemistry_reasons"] = "no_full_inchikey"
            return base, None
        fragments = Chem.GetMolFrags(mol)
        cleaned = rdMolStandardize.Cleanup(mol)
        parent = rdMolStandardize.FragmentParent(cleaned)
        parent_smiles = Chem.MolToSmiles(parent, isomericSmiles=True, canonical=True)
        parent_inchikey = Chem.MolToInchiKey(parent)
        if not parent_smiles or not parent_inchikey:
            base["chemistry_reasons"] = "no_valid_parent"
            return base, None
        no_chirality = Chem.Mol(parent)
        Chem.RemoveStereochemistry(no_chirality)
        connectivity_smiles = Chem.MolToSmiles(
            no_chirality, isomericSmiles=False, canonical=True
        )
        reasons = []
        if len(fragments) != 1:
            reasons.append("salt_or_multifragment_requires_review")
        if any(atom.GetSymbol() not in SMALL_MOLECULE_ELEMENTS for atom in mol.GetAtoms()):
            reasons.append("unsupported_element_or_metal")
        if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
            reasons.append("nonorganic_or_no_carbon")
        if mol.GetNumHeavyAtoms() > 120 or Descriptors.MolWt(mol) > 1500:
            reasons.append("outside_pilot_small_molecule_scope")
        if any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()):
            reasons.append("radical_requires_review")
        if any(flag == "?" for _, flag in Chem.FindMolChiralCenters(
            mol, includeUnassigned=True
        )):
            reasons.append("unassigned_stereocenter")
        if Chem.MolFromSmiles(parent_smiles) is None:
            reasons.append("parent_not_parseable")
        base.update({
            "exact_structure_id": "EX-" + _identity("full-v1", full_smiles),
            "parent_connectivity_id": "PC-" + _identity(
                "connectivity-v1", connectivity_smiles
            ),
            "canonical_smiles": full_smiles,
            "parent_smiles": parent_smiles,
            "full_inchikey": full_inchikey,
            "parent_inchikey": parent_inchikey,
            "chemistry_status": "verified_pilot_scope" if not reasons else "manual_review",
            "chemistry_reasons": "|".join(sorted(reasons)),
        })
        if reasons:
            return base, None
        fp = GENERATOR.GetFingerprint(parent)
        binary = DataStructs.BitVectToBinaryText(fp)
        if len(binary) != FP_BYTES:
            raise RuntimeError(f"Unexpected binary Morgan length: {len(binary)}")
        return base, binary
    except Exception as exc:
        base["chemistry_status"] = "excluded"
        base["chemistry_reasons"] = "rdkit_failure:" + type(exc).__name__
        return base, None
