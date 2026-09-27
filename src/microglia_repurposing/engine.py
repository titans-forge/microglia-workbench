"""Run identity-pinned signature reversal on cached public data."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from microglia_repurposing.annotate import annotate_hits
from microglia_repurposing.contracts import CLAIM_BOUNDARY, RESULT_SCHEMA, sha256_identity
from microglia_repurposing.evaluate import evaluate_ranked, load_gold_standard
from microglia_repurposing.packet import build_packet
from microglia_repurposing.provenance import run_provenance
from microglia_repurposing.scoring import bh_fdr, cmap_reversal, cosine_reversal, permutation_null, rank_aggregate, xsum_score
from microglia_repurposing.signatures import build_disease_signature, build_drug_signatures, disease_vector, drug_matrix

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = PACKAGE_ROOT / "data" / "cache" / "signature_reversal_latest.json"


def _jsonable(frame: pd.DataFrame) -> list[dict[str, Any]]:
    clean = frame.head(25).reset_index().rename(columns={"index": "cmap_name"})
    records = json.loads(clean.to_json(orient="records"))
    return records


def run_engine(*, names: Iterable[str] | None = None, n_perm: int = 2000, output: str | Path | None = DEFAULT_OUTPUT, disease_source: str = "mathys", model: str = "mix", min_hits: int = 3, seed: int = 0) -> dict[str, Any]:
    if n_perm <= 0 or min_hits < 1:
        raise ValueError("invalid_run_configuration")
    disease = build_disease_signature(source=disease_source, model=model)
    d = disease_vector(disease)
    drugs = build_drug_signatures(names=names)
    universe = list(d.index)
    matrix = drug_matrix(drugs, universe)
    scores = pd.DataFrame({
        "cosine_reversal": cosine_reversal(d, matrix),
        "xsum_reversal": xsum_score(d, matrix),
        "set_reversal": cmap_reversal(d, drugs, universe),
        "overlap": (matrix != 0).sum(axis=1),
    })
    min_overlap = min_hits
    scores["status"] = ["scored" if int(v) >= min_overlap else "insufficient_overlap" for v in scores["overlap"]]
    scored = scores[scores["status"] == "scored"]
    if scored.empty:
        scores["mean_rank"] = float("nan")
        scores["cosine_p"] = float("nan")
        scores["fdr"] = float("nan")
    else:
        scores.loc[scored.index, "mean_rank"] = rank_aggregate(scored[["cosine_reversal", "xsum_reversal", "set_reversal"]])
        scores.loc[scored.index, "cosine_p"] = permutation_null(d, matrix.loc[scored.index], scored["cosine_reversal"], n_perm=n_perm, seed=seed)
        scores.loc[scored.index, "fdr"] = bh_fdr(scores.loc[scored.index, "cosine_p"])
    scores["status_order"] = scores["status"].map({"scored": 0, "insufficient_overlap": 1}).fillna(2)
    ranked = scores.sort_values(["status_order", "cosine_reversal"], ascending=[True, False], na_position="last")
    ranked = ranked.drop(columns=["status_order"])
    ranked = annotate_hits(ranked)
    gold = load_gold_standard()
    evaluation = evaluate_ranked(ranked, gold)
    packet = build_packet(disease, ranked, matrix, d)
    result = {
        "schema": RESULT_SCHEMA,
        "research_planning_only": True,
        "claim_eligibility_mutation": False,
        "claim_boundary": CLAIM_BOUNDARY,
        "disease_signature_hash": disease["content_hash"],
        "drug_signature_hash": drugs["content_hash"],
        "disease": {
            "cell_context": disease["cell_context"],
            "cell_state": disease["cell_state"],
            "counts": disease["counts"],
            "source": disease["source"],
            "genes": disease["genes"],
            "l1000_note": disease.get("l1000_note"),
        },
        "drug_panel": {
            "cell_context": drugs["cell_context"],
            "counts": drugs["counts"],
            "source": drugs["source"],
        },
        "configuration": {"model": model, "min_hits": min_hits, "seed": seed, "n_perm": n_perm},
        "statistics": {"tested_statistic": "cosine_reversal", "ranking_statistic": "cosine_reversal",
                       "secondary_scores": "set_reversal and mean_rank are descriptive only",
                       "p_resolution": 1 / (n_perm + 1), "fdr_family_size": len(scored),
                       "fdr_scope": "all coverage-eligible consensus compounds", "triage_fdr_threshold": 0.25,
                       "validation_status": "exploratory"},
        "provenance": run_provenance(),
        "evaluation": evaluation,
        "ranking": json.loads(ranked.rename_axis("cmap_name").reset_index().to_json(orient="records")),
        "packet": packet,
        "top": _jsonable(ranked),
    }
    result["content_hash"] = sha256_identity({k: v for k, v in result.items() if k != "content_hash"})
    if output is not None:
        target = Path(output)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(result, indent=2, ensure_ascii=True, allow_nan=False), encoding="utf-8")
    result["_ranked"] = ranked
    result["_disease"] = disease
    result["_drugs"] = drugs
    return result
