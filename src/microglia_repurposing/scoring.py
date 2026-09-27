"""Reversal scoring with a permutation null. Not a clinical ranking."""
from __future__ import annotations

import numpy as np
import pandas as pd


def cosine_reversal(d: pd.Series, R: pd.DataFrame) -> pd.Series:
    """Higher = stronger anti-correlation on the shared gene universe."""
    common = [gene for gene in R.columns if gene in d.index]
    if not common:
        return pd.Series(0.0, index=R.index, name="cosine_reversal")
    dvec = d.reindex(common).to_numpy(dtype="float64")
    rmat = R.loc[:, common].to_numpy(dtype="float64")
    if not np.isfinite(dvec).all() or not np.isfinite(rmat).all():
        raise ValueError("nonfinite_signature")
    # Scale before squaring so even large finite input cannot overflow.
    dvec = dvec / max(float(np.max(np.abs(dvec))), 1e-300)
    scale = np.maximum(np.max(np.abs(rmat), axis=1), 1e-300)
    rmat = rmat / scale[:, None]
    denom = np.sqrt(np.einsum("ij,ij->i", rmat, rmat)) * np.sqrt(np.sum(dvec ** 2))
    numerator = np.einsum("ij,j->i", rmat, dvec, optimize=False)
    cos = np.divide(numerator, denom, out=np.zeros_like(numerator), where=denom > 0)
    return pd.Series(np.clip(-cos, -1.0, 1.0), index=R.index, name="cosine_reversal")


def xsum_score(d: pd.Series, R: pd.DataFrame) -> pd.Series:
    """Sum of drug signs on disease-up (should be negative) vs disease-down (should be positive)."""
    up = d[d > 0].index
    down = d[d < 0].index
    up_sum = R.reindex(columns=up, fill_value=0).sum(axis=1) if len(up) else pd.Series(0.0, index=R.index)
    down_sum = R.reindex(columns=down, fill_value=0).sum(axis=1) if len(down) else pd.Series(0.0, index=R.index)
    return (down_sum - up_sum).rename("xsum_reversal")


def _ks(hit_ranks: np.ndarray, n_total: int) -> float:
    hits = np.sort(hit_ranks)
    n = hits.size
    if n == 0 or n == n_total:
        return 0.0
    hit_step = 1.0 / n
    miss_step = 1.0 / (n_total - n)
    running = 0.0
    best_pos = 0.0
    best_neg = 0.0
    last = -1
    for rank in hits:
        gap = int(rank) - last - 1
        if gap:
            running -= miss_step * gap
            best_neg = min(best_neg, running)
        running += hit_step
        best_pos = max(best_pos, running)
        last = int(rank)
    tail = (n_total - 1) - last
    if tail > 0:
        running -= miss_step * tail
        best_neg = min(best_neg, running)
    return best_pos if best_pos >= -best_neg else best_neg


def cmap_connectivity(d: pd.Series, ranked_genes: list[str]) -> float:
    """Classic CMap connectivity of a disease up/down query vs one ranked drug list."""
    index = {gene: i for i, gene in enumerate(ranked_genes)}
    up = np.array([index[g] for g in d[d > 0].index if g in index], dtype=int)
    down = np.array([index[g] for g in d[d < 0].index if g in index], dtype=int)
    ks_up = _ks(up, len(ranked_genes))
    ks_down = _ks(down, len(ranked_genes))
    if ks_up * ks_down < 0:
        return (ks_up - ks_down) / 2.0
    return 0.0


def cmap_reversal(d: pd.Series, drug_bundle: dict, universe: list[str]) -> pd.Series:
    """Compatibility entrypoint for tie-aware binary set reversal, not CLUE.

    No within-sign drug ranking is inferred from alphabetic/query order.
    """
    from microglia_repurposing.contexts import ContextKS
    from microglia_repurposing.signatures import drug_matrix
    matrix = drug_matrix(drug_bundle, universe)
    return pd.Series(ContextKS(d, matrix.T).scores(d), index=matrix.index, name="set_reversal")


def bh_fdr(pvalues: pd.Series) -> pd.Series:
    p = pvalues.to_numpy(dtype="float64")
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n, dtype="float64")
    out[order] = np.clip(q, 0.0, 1.0)
    return pd.Series(out, index=pvalues.index, name="fdr")


def permutation_null(d: pd.Series, R: pd.DataFrame, observed: pd.Series, *, n_perm: int = 200, seed: int = 0) -> pd.Series:
    """Gene-label permutation of the disease vector. Returns one-sided p for high reversal."""
    if n_perm <= 0:
        raise ValueError("n_perm_must_be_positive")
    rng = np.random.default_rng(seed)
    values = d.to_numpy(dtype="float64")
    greater = np.zeros(len(observed), dtype="int64")
    R = R.reindex(observed.index)
    names = list(observed.index)
    aligned = observed
    for _ in range(n_perm):
        shuffled = pd.Series(rng.permutation(values), index=d.index)
        null = cosine_reversal(shuffled, R).reindex(aligned.index).fillna(0)
        greater += (null.to_numpy() >= aligned.to_numpy() - 1e-12).astype("int64")
    p = (greater + 1) / (n_perm + 1)
    return pd.Series(p, index=names, name="cosine_p").reindex(observed.index)


def rank_aggregate(scores: pd.DataFrame) -> pd.Series:
    ranks = scores.rank(ascending=False, method="min")
    return ranks.mean(axis=1).rename("mean_rank")


def driver_genes(d: pd.Series, r: pd.Series, top_n: int = 20) -> pd.DataFrame:
    common = d.index.intersection(r.index)
    frame = pd.DataFrame({"disease_sign": d.loc[common], "drug_sign": r.loc[common]})
    frame["reversal_contrib"] = -(frame["disease_sign"] * frame["drug_sign"])
    return frame.sort_values(["reversal_contrib", "disease_sign"], ascending=[False, False]).head(top_n)
