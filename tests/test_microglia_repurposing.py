import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from microglia_repurposing.contracts import CLAIM_BOUNDARY
from microglia_repurposing.evaluate import evaluate_ranked, load_gold_standard
from microglia_repurposing.scoring import cmap_connectivity, cosine_reversal, permutation_null, xsum_score


def test_gold_standard_uses_real_inchikeys_and_does_not_invent_groups():
    gold = load_gold_standard()
    assert gold["drug_inchikey"].str.startswith("group_").sum() == 0
    assert set(gold["cmap_name"]) >= {"vorinostat", "donepezil", "memantine", "pexidartinib"}
    assert {"hdac_probe_positive_control", "symptomatic_non_reversal_negative_control"} <= set(gold["label_class"])


def test_cosine_reversal_recovers_an_anti_correlated_vector():
    d = pd.Series({"APOE": 1.0, "TREM2": 1.0, "P2RY12": -1.0, "CX3CR1": -1.0})
    R = pd.DataFrame(
        {
            "APOE": [-1.0, 1.0],
            "TREM2": [-1.0, 1.0],
            "P2RY12": [1.0, -1.0],
            "CX3CR1": [1.0, -1.0],
        },
        index=["reverser", "mimic"],
    )
    scores = cosine_reversal(d, R)
    assert scores["reverser"] > scores["mimic"]
    assert scores["reverser"] > 0.9


def test_cmap_connectivity_is_zero_when_up_and_down_agree():
    ranked = ["A", "B", "C", "D"]
    d = pd.Series({"A": 1.0, "B": 1.0})
    assert cmap_connectivity(d, ranked) == 0.0
    d_same = pd.Series({"A": 1.0, "B": -1.0})
    cs = cmap_connectivity(d_same, ranked)
    assert -1 <= cs <= 1


def test_xsum_prefers_opposite_signs():
    d = pd.Series({"UP": 1.0, "DOWN": -1.0})
    R = pd.DataFrame({"UP": [-1.0, 1.0], "DOWN": [1.0, -1.0]}, index=["reverser", "mimic"])
    scores = xsum_score(d, R)
    assert scores["reverser"] > scores["mimic"]


def test_permutation_null_is_deterministic():
    d = pd.Series({"A": 1.0, "B": -1.0, "C": 1.0, "D": -1.0})
    R = pd.DataFrame(np.eye(4), index=["w", "x", "y", "z"], columns=["A", "B", "C", "D"])
    observed = cosine_reversal(d, R)
    left = permutation_null(d, R, observed, n_perm=20, seed=1)
    right = permutation_null(d, R, observed, n_perm=20, seed=1)
    assert left.equals(right)


def test_claim_boundary_forbids_treatment_language():
    assert "not diagnosis" in CLAIM_BOUNDARY
    assert "dosing" in CLAIM_BOUNDARY


CACHE = Path(__file__).resolve().parents[1] / "data" / "cache"


@pytest.mark.skipif(not (CACHE / "mathys2019_mic_pathology_de.tsv").is_file(), reason="Mathys DE cache missing")
def test_mathys_signature_is_directional_and_hgnc_pinned(monkeypatch):
    from microglia_repurposing.signatures import build_disease_signature

    monkeypatch.setattr("microglia_repurposing.signatures.spine_overlay", lambda symbol: {"id": "gene:HGNC:613"} if symbol == "APOE" else None)
    disease = build_disease_signature(source="mathys")
    symbols = {gene["gene_symbol"] for gene in disease["genes"]}
    assert disease["source"]["doi"] == "10.1038/s41586-019-1195-2"
    assert disease["cell_state"] == "pathology_vs_no_pathology"
    assert "APOE" in symbols and "SPP1" in symbols
    assert disease["counts"]["up"] >= 10 and disease["counts"]["down"] >= 5
    assert all(gene["hgnc_id"].startswith("HGNC:") for gene in disease["genes"])
    assert any(gene.get("spine_entity_id") == "gene:HGNC:613" for gene in disease["genes"])


@pytest.mark.skipif(
    not (CACHE / "LINCS_L1000_Chem_Pert_Consensus_Sigs.gmt").is_file()
    or not (CACHE / "mathys2019_mic_pathology_de.tsv").is_file(),
    reason="public cache not present",
)
def test_engine_runs_on_real_cached_data():
    from microglia_repurposing.engine import run_engine

    result = run_engine(names=["vorinostat", "donepezil", "memantine", "pioglitazone"], n_perm=20, output=None)
    assert result["research_planning_only"] is True
    assert result["claim_eligibility_mutation"] is False
    assert result["disease"]["cell_context"] == "human_brain_microglia"
    assert result["disease"]["cell_state"] == "pathology_vs_no_pathology"
    assert result["drug_panel"]["cell_context"] == "lincs_l1000_chem_pert_consensus"
    assert result["disease"]["counts"]["up"] >= 10
    assert all(gene["hgnc_id"].startswith("HGNC:") for gene in result["disease"]["genes"])
    ranked = result["_ranked"]
    assert "scored" in set(ranked["status"])
    assert int(ranked.loc[ranked["status"] == "scored", "overlap"].min()) >= 3
    packet = result["packet"]
    assert packet["claim_ids"] == []
    blockers = packet.get("blocking_reasons") or (packet.get("next_experiment") or {}).get("blocking_reasons") or []
    assert "drug_signature_is_not_microglia" in blockers
    if packet.get("next_experiment") is None:
        assert packet.get("empty_reason") in {
            "no_eligible_positive_reversal",
            "no_significant_reversal",
        }
