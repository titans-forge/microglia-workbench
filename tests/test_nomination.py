import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from microglia_repurposing.contracts import CLAIM_BOUNDARY
from microglia_repurposing.ingest import (
    MATHYS_CELL_TYPES,
    MATHYS_SUBPOP_CONTRASTS,
    load_mathys_celltype_contrast,
    load_mathys_mic_selective_de,
    load_mathys_subpopulation_markers,
)
from microglia_repurposing.signatures import build_disease_signature, disease_vector
from microglia_repurposing.ensemble import (
    ENSEMBLE,
    ENSEMBLE_SCHEMA,
    PROTOCOL_SCHEMA,
    default_protocol,
    ensemble_test,
    member_ids_in_order,
    nominate,
    score_member,
    verify_protocol,
)

PACKAGE_ROOT = Path(__file__).resolve().parents[1]


# --- microglia selectivity rule -------------------------------------------------

@pytest.fixture(scope="module")
def selective_mix():
    return load_mathys_mic_selective_de(model="mix")


def test_selective_rule_is_mic_degs_minus_other_cell_type_degs(selective_mix):
    sel = selective_mix["selectivity"]
    assert sel["rule"].startswith("microglia DEGs for this contrast minus the union")
    assert sel["flag_column"] == "deg_ind_mix"
    assert sel["contrast"] == "no_pathology_vs_pathology"

    mic = load_mathys_celltype_contrast("Mic", model="mix")
    mic_degs = set(mic.loc[mic["deg_ind_mix"] == 1, "gene_symbol"])
    shared = set()
    for cell in MATHYS_CELL_TYPES:
        if cell == "Mic":
            continue
        other = load_mathys_celltype_contrast(cell, model="mix")
        shared |= (mic_degs & set(other.loc[other["deg_ind_mix"] == 1, "gene_symbol"]))

    assert sel["microglia_degs"] == len(mic_degs)
    assert set(sel["shared_with_any_other_cell_type"]) == shared
    assert set(sel["selective_genes"]) == (mic_degs - shared)
    # Every reported selective gene really is absent from all five other types.
    for cell, count in sel["shared_counts"].items():
        assert count == len(mic_degs & set(
            load_mathys_celltype_contrast(cell, model="mix")
            .query("deg_ind_mix == 1")["gene_symbol"]))
    assert not (set(sel["selective_genes"]) & shared)


def test_selective_rule_honours_the_model_flag():
    mix = load_mathys_mic_selective_de(model="mix")
    ind = load_mathys_mic_selective_de(model="ind")
    assert mix["flag_column"] == "deg_ind_mix" and ind["flag_column"] == "deg_ind"
    # The cell-level flag calls more DEGs, so it shares more with other cell types
    # and leaves a different selective set. Neither is a subset of the other.
    assert mix["selectivity"]["microglia_degs"] != ind["selectivity"]["microglia_degs"]
    assert set(mix["selectivity"]["selective_genes"]) != set(ind["selectivity"]["selective_genes"])


def test_selectivity_carries_an_explicit_non_claim(selective_mix):
    sel = selective_mix["selectivity"]
    assert "not evidence of" in sel["not_selectivity"]
    assert "not_selectivity" in sel


def test_selectivity_is_not_a_coverage_improvement(selective_mix):
    """The honest trade-off: selectivity costs genes and therefore panel coverage."""
    unrestricted = build_disease_signature(source="mathys", model="mix")
    selective = build_disease_signature(source="mathys_selective", model="mix")
    assert len(selective["genes"]) < len(unrestricted["genes"])
    assert selective["selectivity"]["microglia_degs"] > len(selective["genes"])


# --- subpopulation member -------------------------------------------------------

def test_subpopulation_contrast_is_declared_and_signed():
    assert MATHYS_SUBPOP_CONTRASTS["mic1_vs_mic0"] == ("Mic1", "Mic0")
    raw = load_mathys_subpopulation_markers(contrast="mic1_vs_mic0")
    assert raw["cell_state"] == "Mic1_vs_Mic0_subpopulation_markers"
    assert raw["selectivity"]["up_subpopulation"] == "Mic1"
    assert raw["selectivity"]["down_subpopulation"] == "Mic0"
    assert raw["selectivity"]["up_markers"] == len(raw["up_symbols"])
    assert raw["selectivity"]["down_markers"] == len(raw["down_symbols"])
    # Up genes carry a positive weight, down genes a negative one.
    assert all(raw["weights"][g] > 0 for g in raw["up_symbols"])
    assert all(raw["weights"][g] < 0 for g in raw["down_symbols"])
    # A gene cannot be in both directions.
    assert not (set(raw["up_symbols"]) & set(raw["down_symbols"]))
    assert raw["selectivity"]["ambiguous_direction_genes"] == []


def test_subpopulation_contrast_rejects_unknown_contrasts():
    with pytest.raises(ValueError):
        load_mathys_subpopulation_markers(contrast="mic9_vs_mic0")


def test_subpopulation_member_is_large_enough_to_matter():
    """The subpopulation member exists to fix a power problem; assert it is bigger."""
    case_control = build_disease_signature(source="mathys", model="mix")
    subpop = build_disease_signature(source="mathys_subpop", contrast="mic1_vs_mic0")
    assert len(subpop["genes"]) > 3 * len(case_control["genes"])


# --- ensemble integrity ---------------------------------------------------------

def test_ensemble_members_are_unique_and_ordered():
    ids = member_ids_in_order()
    assert len(ids) == len(set(ids)) == len(ENSEMBLE)
    assert ids[0] == "mathys_mix"
    assert "mathys_subpop_mic1" in ids
    roles = [m.role for m in ENSEMBLE]
    assert len(roles) == len(set(roles))


def test_ensemble_contains_the_declared_primary_and_an_independent_source():
    by_id = {m.member_id: m for m in ENSEMBLE}
    assert by_id["mathys_mix"].role == "primary_default"
    assert by_id["cellmarker_region"].source == "cellmarker"
    assert by_id["cellmarker_region"].weight_basis == "sign_only"


def test_every_member_builds_a_signed_valid_signature():
    from microglia_repurposing.replication import validate_signature

    for member in ENSEMBLE:
        bundle = member.bundle()
        validate_signature(bundle)
        vector = disease_vector(bundle)
        assert (vector > 0).sum() >= 3 and (vector < 0).sum() >= 3
        assert vector.index.is_unique


# --- protocol freeze ------------------------------------------------------------

def test_default_protocol_is_self_consistent_and_verifies():
    protocol = default_protocol(n_perm=100, seed=0, min_hits=3)
    assert protocol["schema"] == PROTOCOL_SCHEMA
    assert protocol["revision"] == 2 and protocol["supersedes_revision"] == 1
    assert protocol["declared_before_results"] is True
    assert protocol["research_planning_only"] is True
    assert verify_protocol(protocol) is protocol
    assert [m["member_id"] for m in protocol["members"]] == member_ids_in_order()


def test_tier_b_thresholds_are_fractions_of_the_ensemble_not_tuned_counts():
    protocol = default_protocol(n_perm=100)
    tier_b = protocol["tier_b"]
    n = tier_b["n_members"]
    assert n == len(ENSEMBLE)
    assert tier_b["min_members_scored"] == int(np.ceil(0.75 * n))
    assert tier_b["min_members_positive"] == int(np.ceil(0.8 * n))
    assert tier_b["min_members_positive"] <= n
    # Tier B must state plainly that it is not Tier A.
    non_claim = tier_b["not_a_substitute_for_tier_a"].lower()
    assert "does not establish significance under the primary analysis" in non_claim
    assert "cannot be reported as if it did" in non_claim


def test_protocol_hash_covers_content():
    protocol = default_protocol(n_perm=100)
    tampered = dict(protocol)
    tampered["tier_a"] = dict(protocol["tier_a"], fdr_max=0.5)
    with pytest.raises(ValueError):
        verify_protocol(tampered)


def test_protocol_refuses_a_changed_member_list():
    protocol = default_protocol(n_perm=100)
    tampered = dict(protocol)
    tampered["members"] = protocol["members"][:-1]
    tampered.pop("content_hash")
    from microglia_repurposing.contracts import sha256_identity
    tampered["content_hash"] = sha256_identity(tampered)
    with pytest.raises(ValueError):
        verify_protocol(tampered)


def test_protocol_rejects_a_stale_revision():
    protocol = default_protocol(n_perm=100)
    tampered = {k: v for k, v in protocol.items() if k != "content_hash"}
    tampered["revision"] = 1
    from microglia_repurposing.contracts import sha256_identity
    tampered["content_hash"] = sha256_identity(tampered)
    with pytest.raises(ValueError):
        verify_protocol(tampered)


def test_default_protocol_rejects_invalid_configuration():
    with pytest.raises(ValueError):
        default_protocol(n_perm=0)
    with pytest.raises(ValueError):
        default_protocol(min_hits=0)


# --- the agreement test ---------------------------------------------------------

@pytest.fixture(scope="module")
def two_members():
    return [score_member(ENSEMBLE[0], min_hits=3), score_member(ENSEMBLE[5], min_hits=3)]


def test_agreement_is_counted_only_where_a_member_is_coverage_eligible(two_members):
    statistics = ensemble_test(two_members, n_perm=25, seed=0)
    scored = statistics["members_scored"]
    agreement = statistics["observed_agreement"]
    assert (agreement <= scored).all()
    assert (scored >= 1).all()
    assert scored.max() <= len(two_members)
    # Every candidate on the union axis is covered by at least one member.
    assert statistics["members_scored"].min() >= 1
    assert len(statistics["candidate_ids"]) == len(set(statistics["candidate_ids"]))
    # The union axis is at least as large as either member's own axis.
    for result in two_members:
        assert len(statistics["candidate_ids"]) >= len(result.member_ids)


def test_members_are_aligned_on_the_union_axis_not_their_own(two_members):
    statistics = ensemble_test(two_members, n_perm=25, seed=0)
    n_global = len(statistics["candidate_ids"])
    for index in range(len(two_members)):
        assert statistics["observed_mean"][index].shape == (n_global,)
        assert statistics["member_p"][index].shape == (n_global,)
    # A candidate scored by only one member has no score in the other.
    alone = [i for i in range(n_global) if statistics["members_scored"][i] == 1]
    assert alone, "expected some candidate covered by a single member"
    i = alone[0]
    missing = [m for m in range(len(two_members)) if np.isnan(statistics["observed_mean"][m][i])]
    assert len(missing) == len(two_members) - 1


def test_agreement_p_is_a_probability_and_respects_the_permutation_floor(two_members):
    n_perm = 25
    statistics = ensemble_test(two_members, n_perm=n_perm, seed=0)
    p = statistics["agreement_p"]
    assert ((p > 0) & (p <= 1)).all()
    assert p.min() >= 1.0 / (n_perm + 1)
    assert statistics["p_resolution"] == pytest.approx(1.0 / (n_perm + 1))


def test_ensemble_test_is_deterministic_for_a_seed(two_members):
    a = ensemble_test(two_members, n_perm=25, seed=7)
    b = ensemble_test(two_members, n_perm=25, seed=7)
    assert np.array_equal(a["agreement_p"], b["agreement_p"])
    assert np.array_equal(a["observed_agreement"], b["observed_agreement"])
    for i in range(len(two_members)):
        assert np.array_equal(a["member_p"][i], b["member_p"][i], equal_nan=True)


def test_ensemble_test_rejects_invalid_permutation_counts(two_members):
    with pytest.raises(ValueError):
        ensemble_test(two_members, n_perm=0, seed=0)
    with pytest.raises(ValueError):
        ensemble_test([], n_perm=10, seed=0)


# --- nomination document --------------------------------------------------------

@pytest.fixture(scope="module")
def nomination(tmp_path_factory):
    protocol = default_protocol(n_perm=40, seed=0, min_hits=3)
    output = tmp_path_factory.mktemp("nomination") / "report.json"
    document = nominate(protocol=protocol, output=output)
    return document, output


def test_nomination_document_is_serializable_without_nan(nomination):
    document, output = nomination
    text = output.read_text(encoding="utf-8")
    json.loads(text)
    assert "NaN" not in text and "Infinity" not in text
    assert document["schema"] == ENSEMBLE_SCHEMA


def test_nomination_keeps_the_claim_boundary(nomination):
    document, _ = nomination
    assert document["research_planning_only"] is True
    assert document["claim_eligibility_mutation"] is False
    assert document["claim_boundary"] == CLAIM_BOUNDARY


def test_nomination_scores_every_declared_member(nomination):
    document, _ = nomination
    assert list(document["members"]) == member_ids_in_order()
    for member_id, report in document["members"].items():
        assert report["genes_in_signature"] >= 6
        assert report["eligible_contexts"] > 0
        assert report["perturbations_eligible"] > 0
        assert report["fdr_family_size"] == report["perturbations_eligible"]
        assert report["signature_hash"]


def test_nomination_rows_carry_per_member_detail(nomination):
    document, _ = nomination
    rows = document["ranking"]
    assert rows
    for row in rows[:50]:
        assert row["members_scored"] <= len(ENSEMBLE)
        assert row["members_positive"] <= row["members_scored"]
        assert 0 < row["agreement_p"] <= 1
        assert isinstance(row["tier_a"], bool) and isinstance(row["tier_b"], bool)
        for member_id in member_ids_in_order():
            entry = row["members"][member_id]
            assert "scored" in entry
            if entry["scored"]:
                assert entry["n_contexts_scored"] > 0
                assert entry["fdr"] is None or entry["fdr"] >= entry["p"]


def test_tier_a_is_the_primary_member_at_the_declared_bar(nomination):
    document, _ = nomination
    fdr_max = document["protocol"]["tier_a"]["fdr_max"]
    primary = document["protocol"]["tier_a"]["primary_member"]
    for row in document["ranking"]:
        entry = row["members"][primary]
        expected = bool(row["named"] and entry.get("positive")
                        and entry.get("fdr") is not None and entry["fdr"] <= fdr_max)
        assert row["tier_a"] is expected
        if row["tier_a"]:
            assert row["primary_fdr"] <= fdr_max


def test_tier_b_requires_agreement_and_membership(nomination):
    document, _ = nomination
    tier_b = document["protocol"]["tier_b"]
    for row in document["ranking"]:
        expected = bool(
            row["named"]
            and row["members_scored"] >= tier_b["min_members_scored"]
            and row["members_positive"] >= tier_b["min_members_positive"]
            and row["agreement_p"] <= tier_b["ensemble_p_max"]
        )
        assert row["tier_b"] is expected


def test_tier_b_is_never_reported_as_tier_a(nomination):
    document, _ = nomination
    for row in document["review_queue"]:
        if row["tier_b"]:
            assert row["tier_a"] is False
    assert "not significance under the" in document["tiers"]["interpretation"]


def test_exclusion_reason_is_absent_for_nominated_rows_and_ordered(nomination):
    document, _ = nomination
    for row in document["ranking"]:
        if row["tier_a"] or row["tier_b"]:
            assert row["exclusion_reason"] is None
        else:
            assert row["exclusion_reason"] in document["protocol"]["exclusions"]


def test_review_queue_respects_the_declared_limit_and_ordering(nomination):
    document, _ = nomination
    limit = document["protocol"]["queue_limit"]
    queue = document["review_queue"]
    assert len(queue) <= limit
    # Tier A rows come before Tier B rows, and the queue is the tier order.
    tiers = [row["tier_a"] for row in queue]
    assert tiers == sorted(tiers, reverse=True)
    assert len(queue) == min(limit, document["tiers"]["tier_a_count"] + document["tiers"]["tier_b_count"])


def test_empty_result_is_a_valid_document(nomination):
    """The primary analysis may legitimately nominate nothing."""
    document, _ = nomination
    if document["tiers"]["tier_a_count"] == 0:
        assert document["tiers"]["tier_a_empty"] is True
        assert all(not row["tier_a"] for row in document["ranking"])
        # Every row that reached neither tier carries a reason; tier B rows do not.
        assert all(row["exclusion_reason"] for row in document["ranking"]
                   if not (row["tier_a"] or row["tier_b"]))
    assert "Neither tier is experimental rescue." in document["tiers"]["interpretation"]


def test_nomination_is_reproducible(nomination, tmp_path):
    document, _ = nomination
    protocol = default_protocol(n_perm=40, seed=0, min_hits=3)
    again = nominate(protocol=protocol, output=tmp_path / "again.json")
    assert again["content_hash"] == document["content_hash"]


def test_nomination_refuses_a_tampered_protocol():
    protocol = default_protocol(n_perm=40)
    protocol["tier_a"] = dict(protocol["tier_a"], fdr_max=0.99)
    with pytest.raises(ValueError):
        nominate(protocol=protocol)


def test_nomination_records_that_the_agreement_null_is_not_independence(nomination):
    document, _ = nomination
    statistics = document["statistics"]
    assert statistics["validation_status"] == "exploratory"
    assert statistics["per_member_family"] == "separate BH family per member; never pooled"
    assert statistics["subgroup_and_best_context_tests"] == "not_performed"
    assert "not independent replication" in " ".join(document["protocol"]["notes"])
