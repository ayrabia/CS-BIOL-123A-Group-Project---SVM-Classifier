"""
features.py — RDKit molecular descriptor and fingerprint computation.

Used primarily for featurizing novel PubChem compounds during screening (Task 6).
The instructor-provided training CSV already contains precomputed descriptors.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import numpy as np
import pandas as pd


def smiles_to_morgan(smiles: str, radius: int = 2, n_bits: int = 2048) -> np.ndarray | None:
    """Compute Morgan fingerprint (ECFP) for a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    return np.array(fp)


def smiles_to_descriptors(smiles: str) -> dict | None:
    """Compute a set of Lipinski/physicochemical descriptors for a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "MW":   Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "HBD":  rdMolDescriptors.CalcNumHBD(mol),
        "HBA":  rdMolDescriptors.CalcNumHBA(mol),
        "TPSA": Descriptors.TPSA(mol),
        "RotBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
    }


def featurize_smiles_list(smiles_list: list[str]) -> pd.DataFrame:
    """Return a DataFrame of Morgan bits + descriptors for a list of SMILES."""
    rows = []
    for smi in smiles_list:
        fp = smiles_to_morgan(smi)
        desc = smiles_to_descriptors(smi)
        if fp is None or desc is None:
            continue
        row = {f"fp_{i}": fp[i] for i in range(len(fp))}
        row.update(desc)
        rows.append(row)
    return pd.DataFrame(rows)
