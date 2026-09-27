import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from microglia_repurposing.contracts import CLAIM_BOUNDARY
from microglia_repurposing.contexts import (
    ContextKS,
    aggregate_context_scores,
    cell_class,
    connectivity_null,
    context_source_block,
    load_harmonizome_contexts,
    parse_context_name,
    weighted_connectivity,
)

CACHE = Path(__file__).resolve().parents[1] / "data" / "cache"
HARMONIZOME = CACHE / "harmonizome_lincs_cmap_chemical_gene_attribute_matrix.txt.gz"
MATHYS = CACHE / "mathys2019_mic_pathology_de.tsv"


def test_parse_context_name_handles_real_headers():
    meta = parse_context_name("BRD-K81418486_Vorinostat_A375_24.0_h_10.0_um")
    assert meta.pert_id == "BRD-K81418486"
    assert meta.pert_iname == "Vorinostat"
    assert meta.cell_id == "A375"
    assert meta.time == "24.0 h"
    assert meta.dose == "10.0 um"
    spaced = parse_context_name("BRD-U94846492_Quinine_hemisulfate_salt_monohydrate_SW948_6.0_h_10.0_um")
    assert spaced.pert_iname == "Quinine_hemisulfate_salt_monohydrate"
    assert spaced.cell_id == "SW948"


def test_parse_context_name_fails_closed_on_junk():
    with pytest.raises(ValueError):
        parse_context_name("not_a_context_header")
    with pytest.raises(ValueError):
        parse_context_name("XX123_A375_24.0_h_10.0_um")


def test_cell_class_labels_are_assay_annotations():
    assert cell_class("THP1") == "monocytic"
    assert cell_class("U937") == "monocytic"
    assert cell_class("NPC") == "cns"
    assert cell_class("NEU") == "cns"
    assert cell_class("ASC") == "cns"
    assert cell_class("MCF7") == "other"


def _toy_disease_and_contexts():
    genes = [f"G{i}" for i in range(20)]
    weights = pd.Series(np.linspace(2.0, -2.0, 20), index=genes, name="disease")
    rng = np.random.default_rng(7)
    ctx = pd.DataFrame(0.0, index=genes, columns=["reverser", "mimic", "noise"])
    # Reverser: -1 on disease-up genes, +1 on disease-down genes.
    ctx.loc[genes[:8], "reverser"] = -1.0
    ctx.loc[genes[12:], "reverser"] = 1.0
    # Mimic: same direction as the disease.
    ctx.loc[genes[:8], "mimic"] = 1.0
    ctx.loc[genes[12:], "mimic"] = -1.0
    # Noise: random sparse signs.
    noise = rng.choice([-1.0, 0.0, 1.0], size=(20,), p=[0.3, 0.4, 0.3])
    ctx["noise"] = noise
    return weights, ctx


def test_weighted_connectivity_separates_reverser_from_mimic():
    d, ctx = _toy_disease_and_contexts()
    scores = weighted_connectivity(d, ctx)
    assert scores["reverser"] > 0.5
    # Signed local set score distinguishes mimicry from no signal.
    assert scores["mimic"] == pytest.approx(-1.0)
    assert scores["reverser"] > scores["noise"]


def test_context_ks_vectorized_matches_scalar_reference():
    d, ctx = _toy_disease_and_contexts()
    ks = ContextKS(d, ctx)
    vec = ks.scores(d)
    ref = weighted_connectivity(d, ctx)
    assert np.allclose(vec, ref.to_numpy())
    shuffled = pd.Series(np.random.default_rng(3).permutation(d.to_numpy()), index=d.index)
    assert np.allclose(ks.scores(shuffled), weighted_connectivity(shuffled, ctx).to_numpy())


def test_connectivity_null_is_deterministic():
    d, ctx = _toy_disease_and_contexts()
    ks = ContextKS(d, ctx)
    pert_ids = ["BRD-1", "BRD-1", "BRD-2"]
    observed = pd.Series({"BRD-1": 0.9, "BRD-2": 0.0})
    left = connectivity_null(ks, d, pert_ids, observed, n_perm=20, seed=1, positions=[0, 1, 2])
    right = connectivity_null(ks, d, pert_ids, observed, n_perm=20, seed=1, positions=[0, 1, 2])
    assert left.equals(right)
    assert ((left >= 1.0 / 21) & (left <= 1.0)).all()


def test_aggregate_reports_cell_class_and_coverage():
    from microglia_repurposing.contexts import ContextMeta

    d, ctx = _toy_disease_and_contexts()
    per_context = pd.DataFrame({
        "weighted_connectivity": [0.8, 0.4, 0.0, 0.0],
        "cosine_reversal": [0.5, 0.2, 0.0, 0.0],
        "hits": [10, 10, 1, 1],
    }, index=["reverser", "mimic", "noise", "extra"])
    meta = [
        ContextMeta("BRD-A", "drugone", "THP1", "6 h", "10 um"),
        ContextMeta("BRD-A", "drugone", "MCF7", "24 h", "10 um"),
        ContextMeta("BRD-B", "drugtwo", "NPC", "24 h", "10 um"),
        ContextMeta("BRD-C", "drugthree", "A549", "24 h", "10 um"),
    ]
    agg = aggregate_context_scores(per_context, meta, min_hits=3)
    row_a = agg.loc["BRD-A"]
    assert row_a["status"] == "scored"
    assert row_a["n_contexts_scored"] == 2
    assert row_a["connectivity_mean"] == pytest.approx(0.6)
    assert row_a["connectivity_monocytic_mean"] == pytest.approx(0.8)
    assert row_a["best_context"].startswith("THP1")
    assert agg.loc["BRD-C", "status"] == "insufficient_coverage"
    assert np.isnan(agg.loc["BRD-C", "connectivity_mean"])


def test_claim_boundary_is_preserved_on_context_outputs():
    block = context_source_block()
    assert "not a microglia assay" in block["note"]
    assert "monocytic" in block["note"]
    assert "not diagnosis" in CLAIM_BOUNDARY


def test_harmonizome_loader_fails_closed_when_cache_missing(monkeypatch):
    def missing(name):
        raise ValueError(f"cache_missing:{name}")

    import microglia_repurposing.contexts as contexts

    monkeypatch.setattr(contexts, "cache_file", missing)
    with pytest.raises(ValueError, match="cache_missing"):
        contexts.load_harmonizome_contexts({"APOE": "HGNC:613"})


@pytest.mark.skipif(not HARMONIZOME.is_file() or not MATHYS.is_file(), reason="public cache not present")
def test_harmonizome_load_resolves_disease_universe():
    from microglia_repurposing.signatures import build_disease_signature

    disease = build_disease_signature(source="mathys")
    universe = {g["gene_symbol"]: g["hgnc_id"] for g in disease["genes"]}
    frame, meta, report = load_harmonizome_contexts(universe)
    assert report["genes_matched"] == 24  # mixed-model default; unavailable genes remain explicit
    assert report["contexts_valid"] >= 30000
    assert set(report["unmatched_symbols"]) <= set(universe)
    assert len(meta) == frame.shape[1]
    assert all(m.cell_id for m in meta)
    assert not report["duplicate_hgnc_dropped"]


@pytest.mark.skipif(not HARMONIZOME.is_file() or not MATHYS.is_file(), reason="public cache not present")
def test_context_engine_runs_on_real_cached_data():
    from microglia_repurposing.context_engine import run_context_engine
    from microglia_repurposing.evaluate import load_gold_standard

    names = list(load_gold_standard()["cmap_name"])
    result = run_context_engine(names=names, n_perm=50, output=None)
    assert result["schema"] == "titans-med-context-reversal-result/v2"
    assert result["research_planning_only"] is True
    assert result["claim_eligibility_mutation"] is False
    assert result["disease"]["cell_context"] == "human_brain_microglia"
    unresolved = {i["cmap_name"] for i in result["context_panel"]["counts"]["unresolved_names"]}
    assert "pexidartinib" in unresolved  # missing in LINCS, reported not filled
    ranked = result["_ranked"]
    assert "scored" in set(ranked["status"])
    scored = ranked[ranked["status"] == "scored"]
    assert (scored["n_contexts_scored"] > 0).all()
    assert scored["connectivity_p"].between(0, 1).all()
    assert scored["fdr"].between(0, 1).all()
    # Known HDAC probe positive control is resolved by name and identity.
    if "BRD-K81418486" in scored.index:
        assert scored.loc["BRD-K81418486", "pert_iname"].lower() == "vorinostat"
    packet = result["packet"]
    assert packet["claim_ids"] == []
    blockers = packet.get("blocking_reasons") or (packet.get("next_experiment") or {}).get("blocking_reasons") or []
    assert "drug_signature_is_not_microglia" in blockers


@pytest.mark.skipif(not HARMONIZOME.is_file() or not MATHYS.is_file(), reason="public cache not present")
def test_context_engine_output_is_valid_json():
    from microglia_repurposing.context_engine import run_context_engine

    result = run_context_engine(names=["vorinostat", "donepezil"], n_perm=20, output=None)
    clean = {k: v for k, v in result.items() if not k.startswith("_")}
    text = json.dumps(clean, ensure_ascii=True, allow_nan=False)
    assert "research_planning_only" in text
    assert "not a microglia assay" in text