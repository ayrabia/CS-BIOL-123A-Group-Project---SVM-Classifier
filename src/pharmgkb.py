"""
pharmgkb.py — PharmGKB annotation lookup for predicted active compounds (Task 7).

Cross-references each predicted inhibitor with PharmGKB to assess
pharmacogenomic relevance, known drug interactions, and existing annotations.
"""

import requests

PHARMGKB_BASE = "https://api.pharmgkb.org/v1/data"


def lookup_by_name(drug_name: str) -> dict | None:
    """
    Query PharmGKB for a drug by name.
    Returns the first matching entry or None if not found.
    """
    url = f"{PHARMGKB_BASE}/chemical"
    params = {"name": drug_name, "view": "base"}
    resp = requests.get(url, params=params, timeout=15)
    if resp.status_code != 200:
        return None
    data = resp.json().get("data", [])
    return data[0] if data else None


def lookup_by_pubchem_cid(cid: int) -> dict | None:
    """
    Query PharmGKB for a compound by PubChem CID via cross-reference.
    Returns the first matching entry or None if not found.
    """
    url = f"{PHARMGKB_BASE}/chemical"
    params = {"crossReferences.resource": "PubChem Compound", "crossReferences.resourceId": str(cid), "view": "base"}
    resp = requests.get(url, params=params, timeout=15)
    if resp.status_code != 200:
        return None
    data = resp.json().get("data", [])
    return data[0] if data else None


def annotate_predictions(predictions_df):
    """
    Add PharmGKB annotations to a predictions DataFrame (output of screen.py).
    Adds columns: pharmgkb_id, pharmgkb_name, pharmgkb_found.
    """
    ids, names, found = [], [], []
    for cid in predictions_df["CID"]:
        result = lookup_by_pubchem_cid(cid)
        if result:
            ids.append(result.get("id", ""))
            names.append(result.get("name", ""))
            found.append(True)
        else:
            ids.append(None)
            names.append(None)
            found.append(False)

    predictions_df = predictions_df.copy()
    predictions_df["pharmgkb_id"]   = ids
    predictions_df["pharmgkb_name"] = names
    predictions_df["pharmgkb_found"] = found
    return predictions_df
