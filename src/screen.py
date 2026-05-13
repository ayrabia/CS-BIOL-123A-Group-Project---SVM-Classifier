"""
screen.py — PubChem compound retrieval and SVC-based screening (Task 6).

Retrieves candidate compounds from PubChem via 2D similarity search around
approved BRAF inhibitor seed compounds, computes RDKit descriptors, and
applies the trained screening SVC.
"""

import time
import requests
import numpy as np
import pandas as pd
from src.features import compute_rdkit_features

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

SEED_NAMES = ["vemurafenib", "dabrafenib", "encorafenib", "sorafenib"]


def get_cid_by_name(name: str) -> int | None:
    """Return the top PubChem CID for a compound name."""
    url = f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/cids/JSON"
    resp = requests.get(url, timeout=30)
    if resp.status_code != 200:
        return None
    cids = resp.json().get("IdentifierList", {}).get("CID", [])
    return cids[0] if cids else None


def get_similar_cids(seed_cid: int, threshold: int = 85, max_results: int = 150) -> list[int]:
    """Return CIDs with 2D Tanimoto similarity >= threshold to seed_cid."""
    url = (
        f"{PUBCHEM_BASE}/compound/fastsimilarity_2d/cid/{seed_cid}/cids/JSON"
        f"?Threshold={threshold}&MaxRecords={max_results}"
    )
    resp = requests.get(url, timeout=60)
    if resp.status_code != 200:
        return []
    return resp.json().get("IdentifierList", {}).get("CID", [])


def get_smiles_batch(cids: list[int], batch_size: int = 100) -> dict[int, str]:
    """Fetch canonical SMILES for a list of CIDs in batches."""
    result = {}
    for i in range(0, len(cids), batch_size):
        batch = cids[i:i + batch_size]
        cid_str = ",".join(map(str, batch))
        url = f"{PUBCHEM_BASE}/compound/cid/{cid_str}/property/SMILES,MolecularWeight/JSON"
        resp = requests.get(url, timeout=30)
        if resp.status_code != 200:
            continue
        for p in resp.json().get("PropertyTable", {}).get("Properties", []):
            if "SMILES" in p:
                result[p["CID"]] = p["SMILES"]
        time.sleep(0.2)
    return result


def build_candidate_library(training_cids: list[int] | None = None) -> list[int]:
    """
    Retrieve candidate CIDs via similarity search around the four approved
    BRAF inhibitors (vemurafenib, dabrafenib, encorafenib, sorafenib).
    Excludes any CIDs already in training_cids.
    """
    seed_cids = {}
    for name in SEED_NAMES:
        cid = get_cid_by_name(name)
        if cid:
            seed_cids[name] = cid
            print(f"  {name} → CID {cid}")
        time.sleep(0.3)

    candidate_set = set(seed_cids.values())
    for name, seed_cid in seed_cids.items():
        similar = get_similar_cids(seed_cid)
        print(f"  Similarity to {name}: {len(similar)} compounds")
        candidate_set.update(similar)
        time.sleep(0.5)

    if training_cids:
        candidate_set -= set(training_cids)

    return list(candidate_set)


def screen_compounds(model, training_cids: list[int] | None = None) -> pd.DataFrame:
    """
    Full screening pipeline:
    1. Build candidate library via similarity search around seed BRAF inhibitors
    2. Fetch SMILES for all candidates
    3. Compute RDKit features (2048-bit Morgan ECFP4 + 6 physicochemical)
    4. Apply trained screening SVC
    Returns a DataFrame with CID, SMILES, predicted class, and decision score.
    """
    print("Building candidate library...")
    candidate_cids = build_candidate_library(training_cids)
    print(f"Total unique candidates (excluding training set): {len(candidate_cids)}")

    print("Fetching SMILES...")
    cid_smiles = get_smiles_batch(candidate_cids)

    rows, cids_out, smiles_out = [], [], []
    for cid, smi in cid_smiles.items():
        feat = compute_rdkit_features(smi)
        if feat is not None:
            rows.append(feat)
            cids_out.append(cid)
            smiles_out.append(smi)

    X = np.array(rows)
    preds = model.predict(X)
    scores = model.decision_function(X)

    results = pd.DataFrame({
        "CID":            cids_out,
        "SMILES":         smiles_out,
        "predicted":      preds,
        "decision_score": scores,
    })
    return results.sort_values("decision_score", ascending=False)
