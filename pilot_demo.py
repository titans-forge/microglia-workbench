"""Deprecated synthetic demo. Use run_engine.py on cached public data."""
import sys
from pathlib import Path

print("The synthetic planted-reverser demo has been replaced by the real engine.")
print("Run: PYTHONPATH=src python run_engine.py")
sys.exit(0)

import numpy as np
import pandas as pd

from microglia_repurposing.signatures import build_synthetic_dam_signature, normalize_signature, consensus_signature
from microglia_repurposing.scoring import cosine_reversal, xsum_score, rank_aggregate, driver_genes
from microglia_repurposing.annotate import annotate_hits

# 1. Disease signature (synthetic DAM-like)
dis = build_synthetic_dam_signature()
d = normalize_signature(dis)
print(f"Disease vector: {len(d)} genes | top up: {list(d.nlargest(5).index)} | top down: {list(d.nsmallest(5).index)}")

# 2. Synthetic drug panel: 8 drugs, 2 true reversers, 6 noise
rng = np.random.default_rng(42)
genes = d.index.tolist()
def make_drug(name, pid, reversal_strength=0.0, moa="unknown"):
    base = rng.normal(0, 1, size=len(genes))
    if reversal_strength:
        base = base * 0.4 - d.to_numpy() * reversal_strength  # anti-correlated
    rows = [{"pert_id": pid, "cmap_name": name, "moa": moa,
             "cell_id": c, "dose": "10uM", "time": "24h",
             "gene_symbol": g, "zscore": z}
            for c in ["THP-1", "NEU"] for g, z in zip(genes, base + rng.normal(0, 0.2, len(genes)))]
    return rows

panel = []
panel += make_drug("sirolimus-like_mTOR", "BRD-001", 1.2, "mTOR inhibitor")
panel += make_drug("pioglitazone-like_PPARg", "BRD-002", 1.0, "PPAR receptor agonist")
panel += make_drug("vorinostat-like_HDAC", "BRD-003", 0.7, "HDAC inhibitor")
for i in range(4, 9):
    panel += make_drug(f"noise-{i}", f"BRD-00{i}", 0.0, "unknown")

drug_long = pd.DataFrame(panel)
R = consensus_signature(drug_long)

# 3. Score
scores = pd.DataFrame({
    "cosine_reversal": cosine_reversal(d, R),
    "xsum_reversal": xsum_score(d, R, top_n=50),
})
scores["mean_rank"] = rank_aggregate(scores[["cosine_reversal", "xsum_reversal"]])
ranked = scores.sort_values("mean_rank")

chembl_meta = pd.DataFrame([
    {"pert_id": "BRD-001", "cmap_name": "sirolimus-like_mTOR", "target": "MTOR", "phase": 4},
    {"pert_id": "BRD-002", "cmap_name": "pioglitazone-like_PPARg", "target": "PPARG", "phase": 4},
    {"pert_id": "BRD-003", "cmap_name": "vorinostat-like_HDAC", "target": "HDAC1", "phase": 4},
])
ranked = annotate_hits(ranked, chembl_meta)

print("\nRanked candidates (higher cosine/xsum = stronger reversal):")
print(ranked[["cosine_reversal", "xsum_reversal", "mean_rank"]].round(3))

# 4. Explain top hit
top = ranked.index[0]
print(f"\nTop hit: {top} | driver genes (disease_z x drug_z reversed):")
print(driver_genes(d, R.loc[top]).round(3))
print("\nOK — synthetic pilot recovered reversers at top. Next: replace synthetic with Mathys Mic1 + L1000 NEU/THP-1.")
