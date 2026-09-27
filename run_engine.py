#!/usr/bin/env python3
"""Run the real signature-reversal engine on cached public data."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from microglia_repurposing.engine import run_engine
from microglia_repurposing.evaluate import load_gold_standard


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all-lincs-consensus", action="store_true", help="Score every Enrichr LINCS consensus compound, not just the exploratory control panel.")
    parser.add_argument("--n-perm", type=int, default=2000)
    parser.add_argument("--output", default=str(Path(__file__).resolve().parent / "data" / "cache" / "signature_reversal_latest.json"))
    parser.add_argument("--model", choices=["ind", "mix"], default="mix")
    parser.add_argument("--min-hits", type=int, default=3)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)
    names = None if args.all_lincs_consensus else list(load_gold_standard()["cmap_name"])
    result = run_engine(names=names, n_perm=args.n_perm, output=args.output, model=args.model, min_hits=args.min_hits, seed=args.seed)
    ranked = result["_ranked"]
    print(result["claim_boundary"])
    print(
        f"Disease {result['disease']['cell_context']} / {result['disease']['cell_state']} "
        f"up={result['disease']['counts']['up']} down={result['disease']['counts']['down']}"
    )
    print(f"Drug panel {result['drug_panel']['cell_context']} n={result['drug_panel']['counts']['drugs']}")
    ev = result["evaluation"]
    print(f"Exploratory controls matched={ev['matched_count']} missing={len(ev['missing'])} auroc={ev['auroc']}")
    if ev["missing"]:
        print("Control-panel gaps:")
        for item in ev["missing"]:
            print(f"  {item['cmap_name']}: {item['reason']}")
    cols = [c for c in ("cosine_reversal", "xsum_reversal", "set_reversal", "mean_rank", "fdr", "tox_flag") if c in ranked.columns]
    print(ranked[cols].head(12).round(4).to_string())
    nxt = result["packet"].get("next_experiment")
    if nxt is None:
        print(f"No next experiment: {result['packet'].get('empty_reason')} blockers={result['packet'].get('blocking_reasons')}")
    else:
        print(f"Next experiment candidate: {nxt['top_cmap_name']} blockers={nxt['blocking_reasons']}")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
