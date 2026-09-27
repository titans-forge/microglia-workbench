"""Declared signature ensemble and cross-signature robustness nomination.

Research-planning only. Every LINCS CMap context is a cell-line assay, not a
microglia assay. A candidate that is positive under several signature
constructions is better supported than one that is positive under a single
construction, but this is still not experimental rescue.

The tested statistic is unchanged: each member is scored with the same
directional_set_ks/v2 connectivity, its own HGNC-pinned gene universe, its own
coverage floor and its own BH family. What is new is that agreement across
members is itself tested against a permutation null instead of being asserted.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from microglia_repurposing.contexts import (
    ContextKS,
    aggregate_context_scores,
    context_source_block,
    fdr_from_pvalues,
    load_harmonizome_contexts,
)
from microglia_repurposing.contracts import require, sha256_identity
from microglia_repurposing.signatures import build_disease_signature, disease_vector

LINCS_MISSING = "-666"
ENSEMBLE_SCHEMA = "titans-med-signature-ensemble-result/v1"
PROTOCOL_SCHEMA = "titans-med-nomination-protocol/v1"


@dataclass(frozen=True)
class EnsembleMember:
    """One declared signature construction. Declared before any result is read."""

    member_id: str
    source: str
    model: str | None
    role: str
    weight_basis: str
    contrast: str | None = None

    def bundle(self) -> dict[str, Any]:
        if self.model is None:
            return build_disease_signature(source=self.source, contrast=self.contrast)
        return build_disease_signature(source=self.source, model=self.model,
                                       contrast=self.contrast)


# The ensemble is fixed here, in code, so it cannot be edited after a run
# without changing the code fingerprint recorded in the report.
#
# Revision 2 adds the MOESM8 subpopulation member. It is the only member with
# substantially more genes than the case-control signature, so it is the first
# construction that reaches a useful share of the LINCS panel. The other MOESM8
# contrasts (Mic3 vs Mic0, Mic2 vs Mic0) are declared in the loader but
# deliberately not members: they are additional activated-state constructions
# from the same table, so they would inflate the agreement count without adding
# an independent line of evidence.
ENSEMBLE: tuple[EnsembleMember, ...] = (
    EnsembleMember("mathys_mix", "mathys", "mix", "primary_default", "log2fc"),
    EnsembleMember("mathys_ind", "mathys", "ind", "sensitivity_cell_level_model", "log2fc"),
    EnsembleMember("mathys_selective_mix", "mathys_selective", "mix", "microglia_selectivity", "log2fc"),
    EnsembleMember("mathys_selective_ind", "mathys_selective", "ind", "microglia_selectivity_sensitivity", "log2fc"),
    EnsembleMember("cellmarker_region", "cellmarker", None, "independent_marker_source", "sign_only"),
    EnsembleMember("mathys_subpop_mic1", "mathys_subpop", None, "independent_cell_state_contrast",
                   "log2fc", "mic1_vs_mic0"),
)


def default_protocol(*, n_perm: int = 2000, seed: int = 0, min_hits: int = 3) -> dict[str, Any]:
    """The nomination rule, declared before it is applied.

    Tier A is the project's existing bar: a named compound with a positive
    aggregate reversal in the primary analysis at FDR <= 0.25. Tier B is a
    declared secondary tier for candidates that are positive under almost every
    signature construction and whose agreement count is itself unusual under a
    permutation null. Tier B is not a substitute for Tier A and is labelled as
    such everywhere it appears. An empty result is a valid published outcome.

    Tier B thresholds are declared as fractions of the ensemble and resolved to
    counts here, so the rule keeps its meaning when the ensemble changes instead
    of being re-tuned to a member count.
    """
    import math

    require(n_perm > 0 and min_hits > 0, "invalid_run_configuration")
    n_members = len(ENSEMBLE)
    require(n_members >= 3, "ensemble_too_small")
    scored_fraction = 0.75
    positive_fraction = 0.8
    min_members_scored = math.ceil(scored_fraction * n_members)
    min_members_positive = math.ceil(positive_fraction * n_members)
    protocol = {
        "schema": PROTOCOL_SCHEMA,
        "revision": 2,
        "supersedes_revision": 1,
        "revision_note": (
            "Revision 1 declared a five-member ensemble and was frozen and run before the MOESM8 "
            "subpopulation member was added. Revision 2 adds that member and restates the Tier B "
            "thresholds as fractions of the ensemble. Revision 1 is kept unmodified and is not "
            "rewritten to look as if it had always contained six members."
        ),
        "research_planning_only": True,
        "declared_before_results": True,
        "members": [asdict(member) for member in ENSEMBLE],
        "configuration": {"n_perm": n_perm, "seed": seed, "min_hits": min_hits},
        "tier_a": {
            "label": "primary_analysis_significant",
            "rule": "named compound, positive aggregate reversal in the primary analysis, primary FDR <= fdr_max",
            "primary_member": "mathys_mix",
            "fdr_max": 0.25,
        },
        "tier_b": {
            "label": "robust_across_signature_constructions",
            "rule": (
                "named compound, scored in at least min_members_scored members, positive in at "
                "least min_members_positive members, and an agreement count at least this high in "
                "no more than ensemble_p_max of independent per-member weight permutations"
            ),
            "n_members": n_members,
            "min_members_scored_fraction": scored_fraction,
            "min_members_positive_fraction": positive_fraction,
            "min_members_scored": min_members_scored,
            "min_members_positive": min_members_positive,
            "ensemble_p_max": 0.05,
            "not_a_substitute_for_tier_a": (
                "Tier B is a cross-construction agreement test. It does not establish significance "
                "under the primary analysis and cannot be reported as if it did."
            ),
        },
        "queue_limit": 25,
        "exclusions": [
            "unnamed_perturbation",
            "insufficient_member_coverage",
            "not_positive_in_primary",
            "insufficient_agreement",
            "agreement_not_unusual_under_null",
        ],
        "notes": [
            "Members share one LINCS panel and one source dataset family; agreement is not "
            "independent replication.",
            "Per-member FDR families are separate and are never pooled.",
            "The replication column is descriptive and read from the frozen replication report; "
            "it is not recomputed here.",
            "Positive means the aggregate mean connectivity is greater than zero. It is a sign "
            "consistency test, not a magnitude or significance test.",
        ],
    }
    protocol["content_hash"] = sha256_identity(protocol)
    return protocol


def verify_protocol(protocol: dict[str, Any]) -> dict[str, Any]:
    require(protocol.get("schema") == PROTOCOL_SCHEMA, "nomination_protocol_schema_invalid")
    require(protocol.get("content_hash") == sha256_identity(
        {k: v for k, v in protocol.items() if k != "content_hash"}), "nomination_protocol_hash_mismatch")
    require([m["member_id"] for m in protocol["members"]] == [m.member_id for m in ENSEMBLE],
            "nomination_protocol_members_changed")
    require(protocol.get("revision") == 2, "nomination_protocol_revision_stale")
    return protocol


def load_or_freeze_protocol(path: Any, **configuration: Any) -> dict[str, Any]:
    """Read the frozen protocol, or write it on first run. Never rewritten silently."""
    import json

    target = Path(path)
    if target.is_file():
        return verify_protocol(json.loads(target.read_text(encoding="utf-8")))
    protocol = default_protocol(**configuration)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(protocol, indent=2, ensure_ascii=True, allow_nan=False), encoding="utf-8")
    return protocol


# Identity and name resolution are imported from the existing engine rather than
# reimplemented, so a compound can never be named one way here and another way
# in the primary analysis. The engine module is read, never modified.
from microglia_repurposing.context_engine import _pert_identity_map, _pert_name_map  # noqa: E402


@dataclass
class MemberResult:
    member: EnsembleMember
    bundle: dict[str, Any]
    frame: pd.DataFrame
    meta: list[Any]
    ks: ContextKS
    vector: pd.Series
    ranked: pd.DataFrame
    eligible_positions: np.ndarray
    eligible_codes: np.ndarray
    eligible_counts: np.ndarray
    member_ids: list[str]
    load_report: dict[str, Any]

    @property
    def eligible_contexts(self) -> int:
        return int(self.eligible_positions.size)


def score_member(member: EnsembleMember, *, min_hits: int) -> MemberResult:
    """Score one declared member on the cached panel with the unchanged statistic."""
    bundle = member.bundle()
    d = disease_vector(bundle)
    universe = {item["gene_symbol"]: item["hgnc_id"] for item in bundle["genes"]}
    frame, meta, load_report = load_harmonizome_contexts(universe)
    d = d.reindex(frame.index)
    require(d.index.is_unique and len(d) > 0, f"member_universe_empty:{member.member_id}")
    ks = ContextKS(d, frame)
    connectivity = ks.scores(d)
    hits = (frame.to_numpy() != 0).sum(axis=0)
    per_context = pd.DataFrame({"weighted_connectivity": connectivity, "hits": hits},
                               index=list(frame.columns))
    ranked = aggregate_context_scores(per_context, meta, min_hits=min_hits)
    eligible = hits >= min_hits
    positions = np.flatnonzero(eligible)
    require(positions.size > 0, f"member_no_eligible_context:{member.member_id}")
    pert_ids = np.array([meta[i].pert_id for i in positions])
    unique_ids = list(pd.unique(pert_ids))
    index_of = {pid: i for i, pid in enumerate(unique_ids)}
    codes = np.array([index_of[pid] for pid in pert_ids], dtype=np.int64)
    counts = np.bincount(codes, minlength=len(unique_ids))
    observed_mean = ranked["connectivity_mean"].reindex(unique_ids).to_numpy(dtype=np.float64)
    require(np.isfinite(observed_mean).all(), f"member_observed_mean_invalid:{member.member_id}")
    return MemberResult(
        member=member, bundle=bundle, frame=frame, meta=meta, ks=ks, vector=d, ranked=ranked,
        eligible_positions=positions, eligible_codes=codes, eligible_counts=counts,
        member_ids=unique_ids, load_report=load_report,
    )


def _member_means(result: MemberResult, scores: np.ndarray) -> np.ndarray:
    return np.bincount(
        result.eligible_codes, weights=scores[result.eligible_positions],
        minlength=len(result.member_ids),
    ) / result.eligible_counts


def _global_axis(results: list[MemberResult]) -> tuple[list[str], list[np.ndarray]]:
    """Union candidate axis plus, per member, its local-to-global scatter index."""
    ids: list[str] = []
    seen: set[str] = set()
    for result in results:
        for pid in result.member_ids:
            if pid not in seen:
                seen.add(pid)
                ids.append(pid)
    position = {pid: i for i, pid in enumerate(ids)}
    scatter = [np.array([position[pid] for pid in result.member_ids], dtype=np.int64)
               for result in results]
    return ids, scatter


def ensemble_test(results: list[MemberResult], *, n_perm: int, seed: int) -> dict[str, Any]:
    """Per-member p-values and a tested cross-member agreement count.

    One permutation loop serves both tests. Every member's weight vector is
    permuted independently, so the agreement null describes how often a
    compound with this coverage pattern would look positive in this many
    members by chance. Members are aligned on the union candidate axis so a
    member that does not cover a compound contributes nothing rather than a
    zero.
    """
    require(results, "ensemble_empty")
    require(n_perm > 0, "n_perm_must_be_positive")
    global_ids, scatter = _global_axis(results)
    n_members = len(results)
    n_global = len(global_ids)
    observed_mean = np.full((n_members, n_global), np.nan, dtype=np.float64)
    observed_positive = np.zeros((n_members, n_global), dtype=bool)
    for i, result in enumerate(results):
        means = _member_means(result, result.ks.scores(result.vector))
        observed_mean[i, scatter[i]] = means
        observed_positive[i, scatter[i]] = means > 0
    scored = ~np.isnan(observed_mean)
    members_scored = scored.sum(axis=0)
    observed_agreement = observed_positive.sum(axis=0)

    rng = np.random.default_rng(seed)
    greater = np.zeros((n_members, n_global), dtype=np.int64)
    agreement_greater = np.zeros(n_global, dtype=np.int64)
    for _ in range(n_perm):
        null_positive = np.zeros((n_members, n_global), dtype=bool)
        for i, result in enumerate(results):
            values = result.vector.to_numpy(dtype=np.float64)
            shuffled = pd.Series(rng.permutation(values), index=result.vector.index)
            means = _member_means(result, result.ks.scores(shuffled))
            index = scatter[i]
            greater[i, index] += (means >= observed_mean[i, index] - 1e-12).astype(np.int64)
            null_positive[i, index] = means > 0
        agreement_greater += (null_positive.sum(axis=0) >= observed_agreement).astype(np.int64)

    member_p = np.full((n_members, n_global), np.nan, dtype=np.float64)
    member_fdr = np.full((n_members, n_global), np.nan, dtype=np.float64)
    for i, result in enumerate(results):
        index = scatter[i]
        p = (greater[i, index] + 1) / (n_perm + 1)
        member_p[i, index] = p
        member_fdr[i, index] = fdr_from_pvalues(
            pd.Series(p, index=result.member_ids)).to_numpy()
    return {
        "candidate_ids": global_ids,
        "n_members": n_members,
        "n_perm": n_perm,
        "p_resolution": 1.0 / (n_perm + 1),
        "observed_mean": observed_mean,
        "observed_positive": observed_positive,
        "members_scored": members_scored,
        "observed_agreement": observed_agreement,
        "member_p": member_p,
        "member_fdr": member_fdr,
        "agreement_p": (agreement_greater + 1) / (n_perm + 1),
    }


def load_replication_concordance(path: Any) -> dict[str, dict[str, Any]]:
    """Descriptive cross-cohort column read from the frozen replication report.

    Recomputed nowhere. States share a cohort and the LINCS panel, so agreement
    counts are descriptive and never a joint significance claim.
    """
    import json

    target = Path(path)
    if not target.is_file():
        return {}
    report = json.loads(target.read_text(encoding="utf-8"))
    concordance: dict[str, dict[str, Any]] = {}
    for state in report.get("states", []):
        comparison = state.get("comparison")
        if not comparison:
            continue
        for record in comparison.get("records", []):
            entry = concordance.setdefault(record["pert_id"], {"states": {}, "counts": {}})
            entry["states"][state["state"]] = record["comparison"]
            entry["counts"][record["comparison"]] = entry["counts"].get(record["comparison"], 0) + 1
    return concordance


def nominate(*, protocol: dict[str, Any], replication_path: Any = None, output: Any = None) -> dict[str, Any]:
    """Run every declared member, test agreement, and build the review queue."""
    import json

    verify_protocol(protocol)
    configuration = protocol["configuration"]
    min_hits = int(configuration["min_hits"])
    results = [score_member(member, min_hits=min_hits) for member in ENSEMBLE]
    statistics = ensemble_test(results, n_perm=int(configuration["n_perm"]), seed=int(configuration["seed"]))
    member_order = [result.member.member_id for result in results]
    global_ids = statistics["candidate_ids"]
    position = {pid: i for i, pid in enumerate(global_ids)}
    local = [{pid: j for j, pid in enumerate(result.member_ids)} for result in results]

    name_map = _pert_name_map()
    identity_map = _pert_identity_map()
    tier_a_rule = protocol["tier_a"]
    tier_b_rule = protocol["tier_b"]
    primary = member_order.index(tier_a_rule["primary_member"])
    replication = load_replication_concordance(replication_path) if replication_path is not None else {}

    rows: list[dict[str, Any]] = []
    for pid in global_ids:
        g = position[pid]
        name = name_map.get(pid, LINCS_MISSING)
        row: dict[str, Any] = {
            "pert_id": pid,
            "pert_iname": name,
            "named": name != LINCS_MISSING,
            "pubchem_inchikey": identity_map.get(pid, {}).get("pubchem_inchikey"),
            "pubchem_cid": identity_map.get(pid, {}).get("pubchem_cid"),
            "members_scored": int(statistics["members_scored"][g]),
            "members_positive": int(statistics["observed_agreement"][g]),
            "agreement_p": _finite(statistics["agreement_p"][g]),
            "agreement_scope": "members where this compound is coverage-eligible",
            "members": {},
        }
        for i, result in enumerate(results):
            member_id = result.member.member_id
            j = local[i].get(pid)
            if j is None:
                row["members"][member_id] = {"scored": False}
                continue
            ranked_row = result.ranked.loc[pid]
            row["members"][member_id] = {
                "scored": True,
                "connectivity_mean": _finite(statistics["observed_mean"][i, g]),
                "positive": bool(statistics["observed_positive"][i, g]),
                "p": _finite(statistics["member_p"][i, g]),
                "fdr": _finite(statistics["member_fdr"][i, g]),
                "n_contexts_scored": int(ranked_row["n_contexts_scored"]),
                "best_context": _text(ranked_row.get("best_context")),
                "best_context_cell": _text(ranked_row.get("best_context_cell")),
                "connectivity_monocytic_mean": _finite(ranked_row.get("connectivity_monocytic_mean")),
                "connectivity_cns_mean": _finite(ranked_row.get("connectivity_cns_mean")),
                "context_robustness_status": _text(ranked_row.get("context_robustness_status")),
            }
        primary_entry = row["members"][tier_a_rule["primary_member"]]
        row["primary_positive"] = bool(primary_entry.get("positive", False))
        row["primary_fdr"] = primary_entry.get("fdr")
        row["primary_score"] = primary_entry.get("connectivity_mean")
        row["tier_a"] = bool(
            row["named"] and row["primary_positive"]
            and row["primary_fdr"] is not None and row["primary_fdr"] <= tier_a_rule["fdr_max"]
        )
        row["tier_b"] = bool(
            row["named"]
            and row["members_scored"] >= tier_b_rule["min_members_scored"]
            and row["members_positive"] >= tier_b_rule["min_members_positive"]
            and row["agreement_p"] is not None and row["agreement_p"] <= tier_b_rule["ensemble_p_max"]
        )
        row["exclusion_reason"] = _exclusion(row, tier_b_rule)
        if replication:
            entry = replication.get(pid)
            row["replication"] = {
                "states": entry["states"] if entry else {},
                "counts": entry["counts"] if entry else {},
                "note": (
                    "descriptive read of the frozen replication report; states share a cohort "
                    "and the LINCS panel, so agreement is not independent replication"
                ),
            }
        rows.append(row)

    tier_a = sorted([r for r in rows if r["tier_a"]],
                    key=lambda r: (r["primary_fdr"], -r["members_positive"], r["pert_id"]))
    tier_b = sorted([r for r in rows if r["tier_b"] and not r["tier_a"]],
                    key=lambda r: (r["agreement_p"], -r["members_positive"],
                                   r["primary_fdr"] if r["primary_fdr"] is not None else 1.0, r["pert_id"]))
    queue = (tier_a + tier_b)[: int(protocol["queue_limit"])]

    member_report = {}
    for i, result in enumerate(results):
        member_report[result.member.member_id] = {
            "role": result.member.role,
            "weight_basis": result.member.weight_basis,
            "signature_hash": result.bundle["content_hash"],
            "cell_state": result.bundle["cell_state"],
            "genes_in_signature": len(result.bundle["genes"]),
            "genes_resolved_in_panel": int(result.frame.shape[0]),
            "genes_unresolved": result.load_report["unmatched_symbols"],
            "contexts": int(result.frame.shape[1]),
            "eligible_contexts": result.eligible_contexts,
            "perturbations_eligible": len(result.member_ids),
            "fdr_family_size": int((result.ranked["status"] == "scored").sum()),
            "selectivity": result.bundle.get("selectivity"),
        }

    document: dict[str, Any] = {
        "schema": ENSEMBLE_SCHEMA,
        "research_planning_only": True,
        "claim_eligibility_mutation": False,
        "claim_boundary": (
            "research-planning-only signature reversal; not diagnosis, dosing, "
            "pathogenicity, treatment, or a clinical ranking"
        ),
        "protocol": protocol,
        "protocol_hash": protocol["content_hash"],
        "panel_source": context_source_block(),
        "members": member_report,
        "statistics": {
            "tested_statistic": "connectivity_mean",
            "method": "directional_set_ks/v2",
            "per_member_family": "separate BH family per member; never pooled",
            "agreement_statistic": "number of declared members with a positive aggregate mean",
            "agreement_null": "independent per-member gene-label permutation of the disease weights",
            "agreement_coverage_scope": "each member contributes only where it is coverage-eligible",
            "p_resolution": statistics["p_resolution"],
            "n_perm": int(configuration["n_perm"]),
            "subgroup_and_best_context_tests": "not_performed",
            "validation_status": "exploratory",
        },
        "tiers": {
            "tier_a_count": len(tier_a),
            "tier_b_count": len(tier_b),
            "tier_a_empty": len(tier_a) == 0,
            "interpretation": (
                "Tier A applies the project's existing FDR bar to the primary analysis alone. "
                "Tier B is a cross-construction agreement tier and is not significance under the "
                "primary analysis. Neither tier is experimental rescue."
            ),
        },
        "ranking": rows,
        "review_queue": queue,
        "excluded_counts": _exclusion_counts(rows),
        "replication_concordance": {
            "available": bool(replication),
            "source": str(replication_path) if replication_path is not None else None,
            "note": (
                "Descriptive external-cohort direction agreement read from the frozen replication "
                "report. States share a cohort and the LINCS panel; no joint any-state FDR."
            ),
        },
    }
    document["content_hash"] = sha256_identity(
        {k: v for k, v in document.items() if k != "content_hash"})
    if output is not None:
        target = Path(output)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(document, indent=2, ensure_ascii=True, allow_nan=False), encoding="utf-8")
    document["_results"] = results
    return document


def _finite(value: Any) -> float | None:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        return None
    return float(value)


def _text(value: Any) -> str | None:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        return None
    return str(value)


def _exclusion(row: dict[str, Any], tier_b_rule: dict[str, Any]) -> str | None:
    if row["tier_a"] or row["tier_b"]:
        return None
    if not row["named"]:
        return "unnamed_perturbation"
    if row["members_scored"] < tier_b_rule["min_members_scored"]:
        return "insufficient_member_coverage"
    if not row["primary_positive"]:
        return "not_positive_in_primary"
    if row["members_positive"] < tier_b_rule["min_members_positive"]:
        return "insufficient_agreement"
    return "agreement_not_unusual_under_null"


def _exclusion_counts(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        reason = row["exclusion_reason"]
        if reason:
            counts[reason] = counts.get(reason, 0) + 1
    return counts


def queue_summary(document: dict[str, Any]) -> dict[str, Any]:
    """Small, printable summary of a nomination document."""
    return {
        "tier_a": document["tiers"]["tier_a_count"],
        "tier_b": document["tiers"]["tier_b_count"],
        "queue": len(document["review_queue"]),
        "excluded": document["excluded_counts"],
        "members": {
            key: {
                "role": value["role"],
                "signature_genes": value["genes_in_signature"],
                "resolved_in_panel": value["genes_resolved_in_panel"],
                "eligible_contexts": value["eligible_contexts"],
                "perturbations": value["perturbations_eligible"],
            }
            for key, value in document["members"].items()
        },
    }


def member_ids_in_order() -> list[str]:
    return [member.member_id for member in ENSEMBLE]


def iter_queue(document: dict[str, Any]) -> Iterable[dict[str, Any]]:
    for row in document["review_queue"]:
        yield row
