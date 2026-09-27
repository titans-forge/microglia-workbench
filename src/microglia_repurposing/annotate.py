"""Annotation receipts. Flags are not automatic drops and not treatment advice."""
from __future__ import annotations

import pandas as pd

from microglia_repurposing.ingest import load_pubchem_gold_standard, load_pert_info

TOXIC_MOA_FLAG = {"hdac inhibitor", "proteasome inhibitor", "tubulin inhibitor", "topoisomerase inhibitor"}
KNOWN_MOA = {
    "vorinostat": "HDAC inhibitor",
    "sirolimus": "mTOR inhibitor",
    "minocycline": "tetracycline antibiotic",
    "pioglitazone": "PPAR receptor agonist",
    "donepezil": "acetylcholinesterase inhibitor",
    "memantine": "NMDA receptor antagonist",
    "pexidartinib": "CSF1R inhibitor",
}


def annotate_hits(ranked: pd.DataFrame) -> pd.DataFrame:
    out = ranked.copy()
    out["moa"] = [KNOWN_MOA.get(str(name).lower(), "unknown") for name in out.index.astype(str)]
    out["tox_flag"] = out["moa"].str.lower().isin(TOXIC_MOA_FLAG)
    out["bbb_class"] = "unknown"
    out["cell_context"] = "lincs_l1000_chem_pert_consensus"
    pubchem = {str(item.get("Title", "")).lower(): item for item in load_pubchem_gold_standard().get("compounds", [])}
    inchikeys = []
    cids = []
    for name in out.index.astype(str):
        hit = pubchem.get(name.lower(), {})
        inchikeys.append(hit.get("InChIKey"))
        cids.append(hit.get("CID"))
    out["pubchem_inchikey"] = inchikeys
    out["pubchem_cid"] = cids
    try:
        pert = load_pert_info()
        pert_map = pert.groupby(pert["pert_iname"].str.lower())["pert_id"].apply(lambda s: sorted(set(s.astype(str)))).to_dict()
        out["lincs_pert_ids"] = [",".join(pert_map.get(str(name).lower(), [])) or None for name in out.index.astype(str)]
    except Exception:
        out["lincs_pert_ids"] = None
    return out


def annotate_context_hits(ranked: pd.DataFrame) -> pd.DataFrame:
    """Attach mechanism flags using names while preserving LINCS identity."""
    out = ranked.copy()
    out["moa"] = [KNOWN_MOA.get(str(name).lower(), "unknown") for name in out["pert_iname"]]
    out["tox_flag"] = out["moa"].str.lower().isin(TOXIC_MOA_FLAG)
    out["tox_annotation_status"] = np_where_known(out["moa"])
    out["bbb_class"] = "unknown"
    out["cell_context"] = "lincs_l1000_chem_pert_instance"
    return out


def np_where_known(moa: pd.Series) -> list[str]:
    return ["mechanism_flag_only_not_viability_evidence" if v != "unknown" else "unknown" for v in moa]
