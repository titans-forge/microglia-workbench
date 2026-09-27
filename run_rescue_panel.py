#!/usr/bin/env python3
"""Report the frozen DAM rescue panel. Does not invent observations or scores."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from microglia_repurposing.rescue import workspace_status, RESCUE_DIR


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=RESCUE_DIR / "report.json")
    args = parser.parse_args()
    status = workspace_status()
    report = status["report"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": status["status"],
        "protocol_hash": status["protocol_hash"],
        "eligible_rescue_controls": report["eligible_rescue_controls"],
        "rescue_auroc": report["rescue_auroc"],
        "n_observed": report["n_observed"],
        "n_entries": report["n_entries"],
        "report": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
