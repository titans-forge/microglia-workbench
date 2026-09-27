"""Cell-context-resolved reversal scoring from the Harmonizome LINCS CMap matrix.

Research-planning only. Every LINCS CMap context is a cell-line assay
(MCF7, PC3, THP1, NPC, ...) — none of them is a microglia assay. Cell-class
labels below are assay annotations for triage, not microglia claims.
"""
from __future__ import annotations

import gzip
import json
import os
import tempfile
from pathlib import Path
from dataclasses import asdict
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from microglia_repurposing.contracts import require, sha256_identity
from microglia_repurposing.identity import pin_symbol
from microglia_repurposing.ingest import cache_file
from microglia_repurposing.scoring import bh_fdr

HARMONIZOME_MATRIX = "harmonizome_lincs_cmap_chemical_gene_attribute_matrix.txt.gz"
CELL_INFO = "GSE92742_Broad_LINCS_cell_info.txt.gz"

# Assay-context classes for triage. THP1/U937 are monocytic leukemia lines,
# the closest LINCS analog to microglia; NPC/NEU/ASC are CNS-lineage cells.
# These are annotations of the LINCS assay, not statements about microglia.
MONOCYTIC_CELLS = frozenset({"THP1", "U937"})
CNS_CELLS = frozenset({"NPC", "NEU", "ASC"})


@dataclass(frozen=True)
class ContextMeta:
    pert_id: str
    pert_iname: str
    cell_id: str
    time: str
    dose: str


def cell_class(cell_id: str) -> str:
    if cell_id in MONOCYTIC_CELLS:
        return "monocytic"
    if cell_id in CNS_CELLS:
        return "cns"
    return "other"


def parse_context_name(name: str) -> ContextMeta:
    """Parse 'BRD-K81418486_Vorinostat_A375_24.0_h_10.0_um' style headers.

    Split from both ends: pert_id is the first token, the trailing five
    tokens are dose, dose_unit, time, time_unit, cell. Everything between
    is the perturbation name (it may itself contain underscores).
    """
    parts = str(name).split("_")
    require(len(parts) >= 7, f"context_name_unparseable:{name}")
    pert_id = parts[0]
    cell_id = parts[-5]
    time = f"{parts[-4]} {parts[-3]}"
    dose = f"{parts[-2]} {parts[-1]}"
    pert_iname = "_".join(parts[1:-5])
    require(pert_id.startswith("BRD-"), f"context_pert_id_invalid:{name}")
    return ContextMeta(pert_id=pert_id, pert_iname=pert_iname, cell_id=cell_id, time=time, dose=dose)


def load_cell_ids() -> set[str]:
    path = cache_file(CELL_INFO)
    table = pd.read_csv(path, sep="\t", dtype=str, low_memory=False, keep_default_na=False)
    require("cell_id" in table.columns, "cell_info_schema_invalid")
    return set(table["cell_id"])


def _load_harmonizome_contexts(
    universe_hgnc: dict[str, str],
) -> tuple[pd.DataFrame, list[ContextMeta], dict[str, Any]]:
    """Load the Harmonizome matrix restricted to an HGNC-pinned disease universe.

    Returns (context frame [genes x contexts of +/-1/0], context metadata,
    load report). Fails closed if the cache is missing. Duplicate HGNC rows
    are dropped (first wins) and reported, never merged or invented.
    """
    path = cache_file(HARMONIZOME_MATRIX)
    wanted: dict[str, str] = {}
    for symbol, hgnc_id in universe_hgnc.items():
        wanted[hgnc_id] = symbol
    rows: dict[str, np.ndarray] = {}
    duplicate_hgnc: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        handle.readline()  # description row
        context_names = header[3:]
        n_ctx = len(context_names)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            symbol = parts[0]
            hit = pin_symbol(symbol)
            if hit is None or hit["hgnc_id"] not in wanted:
                continue
            official = wanted[hit["hgnc_id"]]
            if official in rows:
                duplicate_hgnc.append(hit["hgnc_id"])
                continue
            values = np.fromstring("\t".join(parts[3:]), sep="\t", dtype=np.float64)
            require(values.size == n_ctx, f"harmonizome_row_length_invalid:{symbol}")
            require(np.isfinite(values).all() and np.isin(values, [-1, 0, 1]).all(), f"harmonizome_values_invalid:{symbol}")
            rows[official] = values

    require(rows, "harmonizome_universe_empty")
    frame = pd.DataFrame(rows).T  # genes x contexts
    frame.columns = context_names
    known_cells = load_cell_ids()
    meta: list[ContextMeta] = []
    keep: list[int] = []
    unrecognized: list[str] = []
    for i, name in enumerate(context_names):
        parsed = parse_context_name(name)
        if parsed.cell_id not in known_cells:
            unrecognized.append(name)
            continue
        meta.append(parsed)
        keep.append(i)
    require(meta, "harmonizome_no_valid_contexts")
    frame = frame.iloc[:, keep]
    report = {
        "file": HARMONIZOME_MATRIX,
        "genes_matched": int(frame.shape[0]),
        "genes_requested": len(universe_hgnc),
        "unmatched_symbols": sorted(set(universe_hgnc.keys()) - set(frame.index)),
        "contexts_total": len(context_names),
        "contexts_valid": len(meta),
        "contexts_unrecognized_cell": len(unrecognized),
        "duplicate_hgnc_dropped": sorted(set(duplicate_hgnc)),
    }
    return frame, meta, report


def context_source_block() -> dict[str, Any]:
    from microglia_repurposing.contracts import file_sha256

    return {
        "library": "Harmonizome LINCS L1000 CMap chemical gene attribute matrix",
        "url": "https://maayanlab.cloud/static/hdfs/harmonizome/data/lincscmapchemical/gene_attribute_matrix.txt.gz",
        "sha256": file_sha256(cache_file(HARMONIZOME_MATRIX)),
        "note": (
            "Per-context (perturbagen x cell line x time x dose) binarized +/-1 "
            "signatures. THP1/U937 are monocytic leukemia lines and the closest "
            "LINCS analog to microglia; NPC/NEU/ASC are CNS-lineage cells. "
            "This is not a microglia assay."
        ),
    }


class ContextKS:
    """Vectorized signed set enrichment against a disease-weight ranking.

    Local directional_set_ks/v2; not native CLUE WTCS, NCS or tau.
    Both extrema are evaluated at tied-weight group boundaries.

    The ranking is the disease genes sorted by continuous weight (log2FC),
    not a binary block ordering. Drug-down genes should sit at the
    disease-up end (high weight) and drug-up genes at the disease-down end.
    Index matrices are precomputed once so permutations only re-rank.
    """

    def __init__(self, d: pd.Series, ctx: pd.DataFrame) -> None:
        require(d.index.is_unique and ctx.index.is_unique, "duplicate_gene_identity")
        require(np.isfinite(d.to_numpy()).all() and np.isfinite(ctx.to_numpy()).all(), "nonfinite_signature")
        order = list(d.sort_values(ascending=False).index)
        self.n = len(order)
        self.gene_pos = {gene: i for i, gene in enumerate(order)}
        genes = [g for g in ctx.index if g in self.gene_pos]
        mat = ctx.loc[genes].to_numpy(dtype=np.float64)  # genes x contexts
        self.gene_rank_idx = np.array([self.gene_pos[g] for g in genes], dtype=np.int64)
        self.up_idx, self.up_sizes = self._pad(mat < 0)  # drug down-regulates
        self.down_idx, self.down_sizes = self._pad(mat > 0)  # drug up-regulates

    def _pad(self, sets: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        sizes = sets.sum(axis=0).astype(np.int64)
        width = int(sizes.max()) if sizes.size and sizes.max() > 0 else 0
        n_ctx = sets.shape[1]
        if width == 0:
            return np.zeros((n_ctx, 0), dtype=np.int64), sizes
        padded = np.zeros((n_ctx, width), dtype=np.int64)
        for c in range(n_ctx):
            idx = np.flatnonzero(sets[:, c])
            # Store disease-universe positions, not matrix row indices:
            # rank_of in scores() is indexed by universe position.
            padded[c, : idx.size] = self.gene_rank_idx[idx]
        return padded, sizes

    def _positions(self, idx_matrix: np.ndarray, sizes: np.ndarray, rank_of: np.ndarray) -> np.ndarray:
        n_ctx, width = idx_matrix.shape
        if width == 0:
            return np.full((n_ctx, 0), self.n, dtype=np.float64)
        positions = rank_of[idx_matrix]
        j = np.arange(width)[None, :]
        positions = np.where(j < sizes[:, None], positions, self.n)
        return np.sort(positions, axis=1)

    @staticmethod
    def _extremes(positions: np.ndarray, sizes: np.ndarray, n: int, starts: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
        """Max/min KS deviation for padded sorted position matrices."""
        if positions.shape[1] == 0:
            zeros = np.zeros(positions.shape[0])
            return zeros, zeros
        n_ctx, width = positions.shape
        j = np.arange(width)[None, :]
        hits = j < sizes[:, None]
        safe_sizes = np.maximum(sizes, 1)[:, None]
        safe_miss = np.maximum(n - sizes, 1)[:, None]
        # positions are the final ranks of weight-tie groups. Only group
        # boundaries may contribute; gene order within ties has no meaning.
        first = np.ones_like(hits)
        last = np.ones_like(hits)
        first[:, 1:] = positions[:, 1:] != positions[:, :-1]
        last[:, :-1] = positions[:, :-1] != positions[:, 1:]
        before_positions = positions if starts is None else starts[np.minimum(positions.astype(int), n - 1)]
        after = (j + 1) / safe_sizes - (positions - j) / safe_miss
        before = j / safe_sizes - (before_positions - j) / safe_miss
        max_dev = np.maximum(0, np.max(np.where(hits & last, after, -np.inf), axis=1))
        min_dev = np.minimum(0, np.min(np.where(hits & first, before, np.inf), axis=1))
        empty = (sizes == 0) | (sizes == n)
        max_dev = np.where(empty, 0.0, max_dev)
        min_dev = np.where(empty, 0.0, min_dev)
        return max_dev, min_dev

    def scores(self, d: pd.Series) -> np.ndarray:
        """Connectivity per context for one disease vector."""
        require(set(d.index) == set(self.gene_pos), "permutation_universe_mismatch")
        ordered = d.sort_values(ascending=False, kind="stable")
        values = ordered.to_numpy()
        require(np.isfinite(values).all(), "nonfinite_signature")
        group_starts = np.r_[0, np.flatnonzero(values[1:] != values[:-1]) + 1]
        group_ends = np.r_[group_starts[1:] - 1, self.n - 1]
        lengths = group_ends - group_starts + 1
        starts = np.repeat(group_starts, lengths)
        ends = np.repeat(group_ends, lengths)
        rank_of = np.empty(self.n, dtype=np.int64)
        for i, gene in enumerate(ordered.index):
            rank_of[self.gene_pos[gene]] = ends[i]
        def enrichment(idx, sizes):
            hi, lo = self._extremes(self._positions(idx, sizes, rank_of), sizes, self.n, starts)
            return np.where(np.isclose(hi, -lo, atol=1e-12, rtol=0), 0, np.where(hi > -lo, hi, lo))
        up = enrichment(self.up_idx, self.up_sizes)
        down = enrichment(self.down_idx, self.down_sizes)
        return np.where(up * down < 0, (up - down) / 2.0, 0.0)


def weighted_connectivity(d: pd.Series, ctx: pd.DataFrame) -> pd.Series:
    """Independent dense-walk reference for directional_set_ks/v2.

    Historical function name retained for callers; enrichment is unweighted
    on sets, with continuous disease weights defining tie groups.
    """
    ordered = d.sort_values(ascending=False, kind="stable")
    values = ordered.to_numpy()
    endpoints = np.r_[np.flatnonzero(values[1:] != values[:-1]), len(values) - 1]
    def enrichment(mask):
        k = int(mask.sum())
        if k == 0 or k == len(mask):
            return 0.0
        walk = np.cumsum(np.where(mask, 1.0 / k, -1.0 / (len(mask) - k)))
        walk = np.r_[0.0, walk[endpoints]]
        hi, lo = float(walk.max()), float(walk.min())
        if np.isclose(hi, -lo, atol=1e-12, rtol=0):
            return 0.0
        return hi if hi > -lo else lo
    scores = []
    for name in ctx.columns:
        signs = ctx[name].reindex(ordered.index).fillna(0).to_numpy()
        up, down = enrichment(signs < 0), enrichment(signs > 0)
        scores.append((up - down) / 2.0 if up * down < 0 else 0.0)
    return pd.Series(scores, index=list(ctx.columns), name="weighted_connectivity")


def aggregate_context_scores(
    per_context: pd.DataFrame,
    meta: list[ContextMeta],
    *,
    min_hits: int = 3,
) -> pd.DataFrame:
    """Per-perturbation aggregation with cell-class breakdown. Deterministic.

    Contexts overlapping fewer than `min_hits` disease genes are excluded
    from the statistics (reported in counts, never scored). Perturbations
    with no scoring context get status=insufficient_coverage.
    """
    require(len(per_context) == len(meta), "context_meta_mismatch")
    frame = per_context.copy()
    frame["pert_id"] = [m.pert_id for m in meta]
    frame["pert_iname"] = [m.pert_iname for m in meta]
    frame["cell_id"] = [m.cell_id for m in meta]
    frame["time"] = [m.time for m in meta]
    frame["dose"] = [m.dose for m in meta]
    frame["cell_class"] = [cell_class(m.cell_id) for m in meta]
    if "hits" not in frame.columns:
        frame["hits"] = 0
    scoring = frame[frame["hits"] >= min_hits].copy()
    grouped = frame.groupby("pert_id", sort=True)
    scored_grouped = scoring.groupby("pert_id", sort=True)
    out = pd.DataFrame({
        "pert_iname": grouped["pert_iname"].first(),
        "n_contexts": grouped.size(),
        "n_contexts_scored": scored_grouped.size().reindex(grouped.size().index).fillna(0).astype(int),
    })
    out["connectivity_mean"] = scored_grouped["weighted_connectivity"].mean()
    out["connectivity_max"] = scored_grouped["weighted_connectivity"].max()
    out["overlap_min"] = scored_grouped["hits"].min()
    out["overlap_median"] = scored_grouped["hits"].median()
    out["positive_context_fraction"] = scored_grouped["weighted_connectivity"].agg(lambda x: float((x > 0).mean()))
    cell_means = scoring.groupby(["pert_id", "cell_id"])["weighted_connectivity"].mean()
    cells = cell_means.groupby(level=0)
    out["n_cells_scored"] = cells.size()
    out["cell_balanced_mean"] = cells.mean()
    out["cell_mean_min"] = cells.min()
    out["cell_mean_max"] = cells.max()
    def leave_cell_min(values):
        if len(values) < 2:
            return np.nan
        return float(((values.sum() - values) / (len(values) - 1)).min())
    out["leave_one_cell_out_min"] = cells.agg(leave_cell_min)
    out["context_robustness_status"] = [
        "unavailable" if pd.isna(n) else "single_cell_only" if n < 2 else "positive_across_cells" if low > 0 else "cell_dependent_or_nonpositive"
        for n, low in zip(out["n_cells_scored"], out["cell_mean_min"])
    ]
    if "cosine_reversal" in scoring.columns:
        out["cosine_mean"] = scored_grouped["cosine_reversal"].mean()
    best = scoring.loc[scoring.groupby("pert_id")["weighted_connectivity"].idxmax()]
    best = best.set_index("pert_id").reindex(out.index)
    out["best_context"] = best.apply(
        lambda row: None if pd.isna(row.get("cell_id")) else f"{row.cell_id} {row.time} {row.dose}", axis=1
    )
    out["best_context_cell"] = best["cell_id"].where(best["cell_id"].notna()).to_numpy()
    out["best_context_time"] = best["time"].where(best["time"].notna()).to_numpy()
    out["best_context_dose"] = best["dose"].where(best["dose"].notna()).to_numpy()
    out["best_context_score"] = best["weighted_connectivity"].to_numpy()
    for label in ("monocytic", "cns", "other"):
        part = scoring[scoring["cell_class"] == label].groupby("pert_id")["weighted_connectivity"].mean()
        out[f"connectivity_{label}_mean"] = part.reindex(out.index)
    out["status"] = ["scored" if n > 0 else "insufficient_coverage" for n in out["n_contexts_scored"]]
    return out.sort_values(["connectivity_mean"], ascending=False, na_position="last")


def connectivity_null(
    ks: ContextKS,
    d: pd.Series,
    pert_ids: list[str],
    observed_mean: pd.Series,
    *,
    n_perm: int = 200,
    seed: int = 0,
    positions: np.ndarray | None = None,
) -> pd.Series:
    """Gene-label permutation null for per-perturbation mean connectivity.

    pert_ids are per scoring context; positions are the matching indices
    into the full context axis (all contexts when None). Permutation
    shuffles the disease weights across the gene universe; the drug sign
    sets stay fixed. One-sided p for high mean reversal.
    """
    require(n_perm > 0, "n_perm_must_be_positive")
    rng = np.random.default_rng(seed)
    codes_arr = np.array(pert_ids)
    drug_index = {pid: i for i, pid in enumerate(pd.unique(codes_arr))}
    codes = np.array([drug_index[p] for p in pert_ids])
    counts = np.bincount(codes, minlength=len(drug_index))
    require((counts > 0).all(), "null_counts_invalid")
    if positions is None:
        positions = np.arange(len(pert_ids))
    positions = np.asarray(positions, dtype=np.int64)
    obs = observed_mean.reindex(list(drug_index.keys())).fillna(0).to_numpy(dtype=np.float64)
    greater = np.zeros(len(drug_index), dtype=np.int64)
    values = d.to_numpy(dtype=np.float64)
    for _ in range(n_perm):
        shuffled = pd.Series(rng.permutation(values), index=d.index)
        null_scores = ks.scores(shuffled)[positions]
        null_means = np.bincount(codes, weights=null_scores, minlength=len(drug_index)) / counts
        greater += (null_means >= obs - 1e-12).astype(np.int64)
    p = (greater + 1) / (n_perm + 1)
    return pd.Series(p, index=list(drug_index.keys()), name="connectivity_p")


def fdr_from_pvalues(pvalues: pd.Series) -> pd.Series:
    return bh_fdr(pvalues)

def load_harmonizome_contexts(universe_hgnc: dict[str, str]):
    """Verified, content-keyed normalized cache; no pickle deserialization."""
    from microglia_repurposing.provenance import fingerprint
    from microglia_repurposing.identity import DEFAULT_HGNC
    source = cache_file(HARMONIZOME_MATRIX)
    key = sha256_identity({"schema": 2, "universe": universe_hgnc,
                           "matrix": fingerprint(source), "hgnc": fingerprint(DEFAULT_HGNC),
                           "cells": fingerprint(cache_file(CELL_INFO)),
                           "loader": fingerprint(Path(__file__)),
                           "identity": fingerprint(Path(__file__).with_name("identity.py"))})
    folder = source.parent / "normalized"
    folder.mkdir(exist_ok=True)
    target = folder / (key + ".npz")
    receipt = target.with_suffix(".sha256")
    if target.exists() and receipt.exists():
        require(fingerprint(target) == receipt.read_text().strip(), "normalized_cache_hash_mismatch")
        with np.load(target, allow_pickle=False) as data:
            frame = pd.DataFrame(data["values"], index=data["genes"].tolist(), columns=data["columns"].tolist())
            meta = [ContextMeta(**m) for m in json.loads(str(data["meta"]))]
            report = json.loads(str(data["report"]))
        report["normalized_cache"] = "hit"
        return frame, meta, report
    frame, meta, report = _load_harmonizome_contexts(universe_hgnc)
    report["normalized_cache_key"] = key
    with tempfile.NamedTemporaryFile(dir=folder, suffix=".npz", delete=False) as handle:
        temporary = Path(handle.name)
        np.savez_compressed(handle, values=frame.to_numpy(dtype=np.int8),
                            genes=np.array(frame.index, dtype=str), columns=np.array(frame.columns, dtype=str),
                            meta=json.dumps([asdict(m) for m in meta]), report=json.dumps(report))
    os.replace(temporary, target)
    receipt.write_text(fingerprint(target) + "\n")
    report["normalized_cache"] = "miss"
    return frame, meta, report
