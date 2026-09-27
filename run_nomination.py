#!/usr/bin/env python3
"""Nominate candidates from a declared signature ensemble.

Research-planning only. Every LINCS CMap context is a cell-line assay, not a
microglia assay. Tier A applies the project's existing FDR bar to the primary
analysis. Tier B is a cross-construction agreement tier and is not significance
under the primary analysis. An empty result is a valid published outcome.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from microglia_repurposing.ensemble import (  # noqa: E402
    ENSEMBLE,
    load_or_freeze_protocol,
    member_ids_in_order,
    nominate,
    queue_summary,
)

PACKAGE_ROOT = Path(__file__).resolve().parent
# Revision 1 lives at data/nomination/protocol.json and is kept unmodified.
# The active revision is 2, which adds the MOESM8 subpopulation member.
DEFAULT_PROTOCOL = PACKAGE_ROOT / "data" / "nomination" / "v2" / "protocol.json"
DEFAULT_OUTPUT = PACKAGE_ROOT / "data" / "nomination" / "v2" / "report.json"
DEFAULT_REPLICATION = PACKAGE_ROOT / "data" / "replication" / "report.json"


def _print_members(summary: dict) -> None:
    print("Ensemble members (declared before scoring):")
    for member in ENSEMBLE:
        entry = summary["members"][member.member_id]
        print(
            f"  {member.member_id:<22} {entry['role']:<34} "
            f"genes={entry['signature_genes']:>3} resolved={entry['resolved_in_panel']:>3} "
            f"eligible_contexts={entry['eligible_contexts']:>5} perturbations={entry['perturbations']:>5}"
        )
        if member.weight_basis != "log2fc":
            print(f"    weight basis: {member.weight_basis} (not a continuous effect size)")


def _print_queue(document: dict) -> None:
    queue = document["review_queue"]
    if not queue:
        print("Review queue: empty. No candidate met either declared tier.")
        return
    print(f"Review queue ({len(queue)}):")
    header = f"{'#':>3} {'compound':<26} {'tier':<6} {'pos/scored':>10} {'agree_p':>8} {'primary_fdr':>11} {'primary_score':>13}"
    print(header)
    print("-" * len(header))
    for i, row in enumerate(queue, start=1):
        tier = "A" if row["tier_a"] else "B"
        fdr = row["primary_fdr"]
        score = row["primary_score"]
        print(
            f"{i:>3} {str(row['pert_iname'])[:26]:<26} {tier:<6} "
            f"{row['members_positive']:>4}/{row['members_scored']:<5} "
            f"{row['agreement_p']:>8.4f} "
            f"{('n/a' if fdr is None else format(fdr, '.4f')):>11} "
            f"{('n/a' if score is None else format(score, '+.4f')):>13}"
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-perm", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--min-hits", type=int, default=3)
    parser.add_argument("--protocol", default=str(DEFAULT_PROTOCOL))
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument(
        "--replication",
        default=str(DEFAULT_REPLICATION),
        help="Frozen replication report read for a descriptive concordance column. Never recomputed here.",
    )
    parser.add_argument("--no-replication", action="store_true")
    args = parser.parse_args(argv)

    requested = {"n_perm": args.n_perm, "seed": args.seed, "min_hits": args.min_hits}
    protocol_path = Path(args.protocol)
    existed = protocol_path.is_file()
    protocol = load_or_freeze_protocol(protocol_path, **requested)
    if existed:
        frozen = protocol["configuration"]
        if any(frozen[key] != requested[key] for key in requested):
            print(
                "Frozen protocol wins. Requested "
                + ", ".join(f"{k}={v}" for k, v in requested.items())
                + " but "
                + ", ".join(f"{k}={frozen[k]}" for k in requested)
                + f" is frozen in {protocol_path}."
            )
        print(f"Protocol {protocol['content_hash'][:16]} verified (frozen, not rewritten).")
    else:
        print(f"Protocol {protocol['content_hash'][:16]} written to {protocol_path} before scoring.")

    replication_path = None if args.no_replication else Path(args.replication)
    if replication_path is not None and not replication_path.is_file():
        print(f"Replication report not found at {replication_path}; concordance column will be absent.")
        replication_path = None

    document = nominate(protocol=protocol, replication_path=replication_path, output=args.output)
    summary = queue_summary(document)

    print(document["claim_boundary"])
    print(
        f"Tested statistic: {document['statistics']['method']} "
        f"({document['statistics']['tested_statistic']}); "
        f"n_perm={document['statistics']['n_perm']} "
        f"p_resolution={document['statistics']['p_resolution']:.5f}"
    )
    print(f"Agreement null: {document['statistics']['agreement_null']}")
    _print_members(summary)

    tiers = document["tiers"]
    tier_b = protocol["tier_b"]
    print(f"Tier A (primary analysis, FDR <= {protocol['tier_a']['fdr_max']}): {tiers['tier_a_count']}")
    print(
        f"Tier B (scored in >= {tier_b['min_members_scored']} and positive in >= "
        f"{tier_b['min_members_positive']} of {tier_b['n_members']} members, "
        f"agreement p <= {tier_b['ensemble_p_max']}): {tiers['tier_b_count']}"
    )
    if tiers["tier_a_empty"]:
        print("Tier A is empty. The primary analysis nominates nothing at the declared bar.")
    print(f"Excluded: {summary['excluded']}")
    print(f"Replication concordance attached: {document['replication_concordance']['available']}")
    _print_queue(document)
    print(f"Wrote {args.output}")
    print(f"Member order: {', '.join(member_ids_in_order())}")

    # Fail loudly if a non-finite value reached the serialized document.
    json.dumps({k: v for k, v in document.items() if not k.startswith("_")}, allow_nan=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
