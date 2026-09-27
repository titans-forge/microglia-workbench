"""Descriptive control-panel coverage, not validated discovery performance."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from microglia_repurposing.contracts import require

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_GOLD = PACKAGE_ROOT / "data" / "control_panel_ad.csv"


def load_gold_standard(path: str | Path | None = None) -> pd.DataFrame:
    target = Path(path) if path is not None else DEFAULT_GOLD
    require(target.is_file(), "gold_standard_missing")
    table = pd.read_csv(target, dtype=str)
    required = {"disease_efo", "drug_inchikey", "cmap_name", "label", "label_class", "source", "pubchem_cid"}
    require(required.issubset(table.columns), "gold_standard_schema_invalid")
    table["label"] = table["label"].astype(int)
    require(set(table["label"]).issubset({0, 1}), "gold_standard_label_invalid")
    require(table["drug_inchikey"].str.fullmatch(r"[A-Z0-9]{14}-[A-Z0-9]{10}-[A-Z0-9]").all(), "gold_standard_inchikey_invalid")
    return table


def _auroc(scores: np.ndarray, labels: np.ndarray) -> float | None:
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if pos.size == 0 or neg.size == 0:
        return None
    greater = np.sum(pos[:, None] > neg[None, :])
    equal = np.sum(pos[:, None] == neg[None, :])
    return float((greater + 0.5 * equal) / (pos.size * neg.size))


def evaluate_ranked(ranked: pd.DataFrame, gold: pd.DataFrame, *, score_col: str = "cosine_reversal") -> dict[str, Any]:
    usable = ranked
    if "overlap" in ranked.columns:
        usable = ranked[ranked["overlap"] >= 3]
    by_name = {str(name).lower(): i for i, name in enumerate(usable.index.astype(str))}
    rows = []
    missing = []
    for item in gold.itertuples(index=False):
        key = str(item.cmap_name).lower()
        if str(item.cmap_name).lower() not in {str(n).lower() for n in ranked.index}:
            missing.append({"cmap_name": item.cmap_name, "reason": "missing_lincs_consensus_signature"})
            continue
        if key not in by_name:
            missing.append({"cmap_name": item.cmap_name, "reason": "insufficient_gene_overlap"})
            continue
        score = float(usable.iloc[by_name[key]][score_col])
        rank = int(usable[score_col].rank(ascending=False, method="min").iloc[by_name[key]])
        rows.append({
            "cmap_name": item.cmap_name,
            "label": int(item.label),
            "label_class": item.label_class,
            "source": item.source,
            "drug_inchikey": item.drug_inchikey,
            "score": score,
            "rank": rank,
        })
    matched = pd.DataFrame(rows)
    result = {
        "matched_count": int(len(matched)),
        "missing": missing,
        "auroc": None,
        "positive_median_rank": None,
        "negative_median_rank": None,
        "rows": rows,
    }
    if not matched.empty and matched["label"].nunique() == 2:
        result["auroc"] = _auroc(matched["score"].to_numpy(dtype="float64"), matched["label"].to_numpy(dtype="int64"))
        result["positive_median_rank"] = float(matched.loc[matched["label"] == 1, "rank"].median())
        result["negative_median_rank"] = float(matched.loc[matched["label"] == 0, "rank"].median())
    return qualify_evaluation(result, gold)


def qualify_evaluation(result: dict[str, Any], controls: pd.DataFrame) -> dict[str, Any]:
    """Retain legacy diagnostics without promoting weak labels to validation."""
    result["panel_type"] = "curated_exploratory_controls"
    result["validation_status"] = "not_validated"
    result["total_controls"] = len(controls)
    result["matched_fraction"] = result["matched_count"] / len(controls) if len(controls) else 0
    result["descriptive_legacy_label_auroc"] = result.get("auroc")
    result["auroc"] = None
    if "auroc_max" in result:
        result["auroc_max"] = None
    result["limitations"] = [
        "Legacy labels combine mechanism, modulation, trial participation and symptomatic use; they are not rescue truth.",
        "Assay-level label evidence and an independent held-out benchmark have not been adjudicated.",
        "Descriptive legacy-label AUROC is not discovery or clinical validation.",
    ]
    evidence = controls.set_index("cmap_name").to_dict(orient="index")
    for row in result["rows"]:
        row["label_evidence_status"] = evidence[row["cmap_name"]].get("evidence_status", "unadjudicated")
    return result
