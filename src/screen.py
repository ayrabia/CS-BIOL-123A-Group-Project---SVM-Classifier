"""
screen.py — PubChem compound retrieval and SVC-based screening (Task 6).

Retrieves candidate compounds from PubChem relevant to BRAF V600E inhibition,
computes RDKit descriptors, and applies the trained SVC model.
"""

import requests
import pandas as pd
from src.features import featurize_smiles_list

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def get_cids_by_name(query: str, max_results: int = 100) -> list[int]:
    """Search PubChem for compound CIDs matching a query string."""
    url = f"{PUBCHEM_BASE}/compound/name/{query}/cids/JSON"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    cids = resp.json().get("IdentifierList", {}).get("CID", [])
    return cids[:max_results]


def get_smiles_for_cids(cids: list[int]) -> dict[int, str]:
    """Fetch isomeric SMILES for a list of PubChem CIDs."""
    cid_str = ",".join(map(str, cids))
    url = f"{PUBCHEM_BASE}/compound/cid/{cid_str}/property/IsomericSMILES/JSON"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    props = resp.json().get("PropertyTable", {}).get("Properties", [])
    return {p["CID"]: p["IsomericSMILES"] for p in props if "IsomericSMILES" in p}


def screen_compounds(model, query: str = "BRAF inhibitor", max_results: int = 100) -> pd.DataFrame:
    """
    Full screening pipeline:
    1. Retrieve CIDs from PubChem
    2. Fetch SMILES
    3. Compute RDKit features
    4. Apply trained SVC
    Returns a DataFrame with CID, SMILES, and predicted class/probability.
    """
    print(f"Fetching CIDs for '{query}'...")
    cids = get_cids_by_name(query, max_results)
    print(f"  {len(cids)} CIDs retrieved.")

    print("Fetching SMILES...")
    cid_smiles = get_smiles_for_cids(cids)

    smiles_list = list(cid_smiles.values())
    cid_list = list(cid_smiles.keys())

    print("Computing RDKit features...")
    X = featurize_smiles_list(smiles_list)

    preds = model.predict(X)
    probs = model.predict_proba(X)[:, 1]

    results = pd.DataFrame({
        "CID":         cid_list[: len(X)],
        "SMILES":      smiles_list[: len(X)],
        "predicted":   preds,
        "prob_active": probs,
    })
    return results.sort_values("prob_active", ascending=False)
