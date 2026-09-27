#!/usr/bin/env python3
"""Run the cell-context-resolved reversal engine on cached public data.

Research-planning only. LINCS CMap contexts are cell-line assays, not
microglia assays.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from microglia_repurposing.context_engine import run_context_engine
from microglia_repurposing.evaluate import load_gold_standard


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--all-harmonizome",
        action="store_true",
        help="Score every Harmonizome LINCS CMap perturbation context, not just the exploratory control panel.",
    )
    parser.add_argument("--n-perm", type=int, default=2000)
    parser.add_argument(
        "--output",
        default=str(Path(__file__).resolve().parent / "data" / "cache" / "context_reversal_latest.json"),
    )
    parser.add_argument("--model", choices=["ind", "mix"], default="mix")
    parser.add_argument("--min-hits", type=int, default=3)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)
    names = None if args.all_harmonizome else list(load_gold_standard()["cmap_name"])
    result = run_context_engine(names=names, n_perm=args.n_perm, output=args.output, model=args.model, min_hits=args.min_hits, seed=args.seed)
    ranked = result["_ranked"]
    print(result["claim_boundary"])
    print(
        f"Disease {result['disease']['cell_context']} / {result['disease']['cell_state']} "
        f"up={result['disease']['counts']['up']} down={result['disease']['counts']['down']}"
    )
    panel = result["context_panel"]
    print(
        f"Context panel {panel['counts']['perturbations']} perturbations, "
        f"{panel['counts']['contexts']} contexts, scored={panel['counts']['perturbations_scored']}"
    )
    for item in panel["counts"]["unresolved_names"]:
        print(f"  unresolved: {item['cmap_name']}: {item['reason']}")
    ev = result["evaluation"]
    print(f"Controls matched={ev['matched_count']} missing={len(ev['missing'])} auroc={ev['auroc']} auroc_max={ev['auroc_max']}")
    for item in ev["missing"]:
        print(f"  {item['cmap_name']}: {item['reason']}")
    cols = [
        c
        for c in (
            "pert_iname", "n_contexts", "n_contexts_scored", "connectivity_mean",
            "connectivity_max", "connectivity_monocytic_mean", "connectivity_cns_mean",
            "best_context", "connectivity_p", "fdr", "status",
        )
        if c in ranked.columns
    ]
    print(ranked[cols].head(15).round(4).to_string())
    nxt = result["packet"].get("next_experiment")
    if nxt is None:
        print(f"No next experiment: {result['packet'].get('empty_reason')} blockers={result['packet'].get('blocking_reasons')}")
    else:
        print(
            f"Next experiment candidate: {nxt.get('top_pert_iname')} "
            f"[{nxt.get('top_best_context')}] fdr={nxt.get('fdr')} blockers={nxt.get('blocking_reasons')}"
        )
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())