"""Cell-context-resolved reversal run on cached public data.

Research-planning only. LINCS CMap contexts are cell-line assays, not
microglia assays. Hits are hypotheses for a declared microglia experiment.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from microglia_repurposing.contexts import (
    ContextKS,
    aggregate_context_scores,
    connectivity_null,
    context_source_block,
    fdr_from_pvalues,
    load_harmonizome_contexts,
)
from microglia_repurposing.contracts import (
    CLAIM_BOUNDARY,
    RESULT_SCHEMA_CONTEXT,
    require,
    sha256_identity,
)
from microglia_repurposing.evaluate import _auroc, load_gold_standard, qualify_evaluation
from microglia_repurposing.annotate import annotate_context_hits
from microglia_repurposing.provenance import run_provenance
from microglia_repurposing.ingest import load_pert_info
from microglia_repurposing.packet import build_packet
from microglia_repurposing.scoring import cosine_reversal
from microglia_repurposing.signatures import build_disease_signature, disease_vector

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = PACKAGE_ROOT / "data" / "cache" / "context_reversal_latest.json"
LINCS_MISSING = "-666"


def _pert_name_map() -> dict[str, str]:
    """pert_id -> canonical LINCS pert_iname. Real names win over '-666'."""
    pert = load_pert_info()
    mapping: dict[str, str] = {}
    for row in pert[["pert_id", "pert_iname"]].itertuples(index=False):
        pid, iname = str(row.pert_id), str(row.pert_iname)
        current = mapping.get(pid)
        if current is None or (current == LINCS_MISSING and iname != LINCS_MISSING):
            mapping[pid] = iname
    return mapping


def _pert_identity_map() -> dict[str, dict[str, str]]:
    """pert_id -> LINCS-pinned InChIKey / PubChem CID, never invented."""
    pert = load_pert_info()
    mapping: dict[str, dict[str, str]] = {}
    for row in pert[["pert_id", "inchi_key", "pubchem_cid"]].itertuples(index=False):
        pid = str(row.pert_id)
        inchi = str(row.inchi_key) if str(row.inchi_key) not in {"", LINCS_MISSING} else None
        cid = str(row.pubchem_cid) if str(row.pubchem_cid) not in {"", LINCS_MISSING} else None
        current = mapping.get(pid)
        if current is None or (current.get("pubchem_inchikey") is None and inchi):
            mapping[pid] = {"pubchem_inchikey": inchi, "pubchem_cid": cid}
    return mapping


def _resolve_names(names: Iterable[str], pert_name_map: dict[str, str]) -> tuple[set[str], list[dict[str, str]]]:
    """Map cmap-style names to LINCS pert_ids. Unresolved names are reported."""
    lowered = {str(name).lower(): str(name) for name in names}
    wanted: set[str] = set()
    missing: list[dict[str, str]] = []
    by_pert = {}
    for pert_id, pert_iname in pert_name_map.items():
        key = pert_iname.lower()
        if key in lowered:
            wanted.add(pert_id)
            by_pert.setdefault(key, []).append(pert_id)
    for key, original in lowered.items():
        if key not in by_pert:
            missing.append({"cmap_name": original, "reason": "missing_lincs_perturbation"})
    return wanted, missing


def _jsonable(ranked: pd.DataFrame, limit: int = 25) -> list[dict[str, Any]]:
    columns = [
        "pert_id", "pert_iname", "tox_flag", "tox_annotation_status", "moa", "n_contexts", "n_contexts_scored", "connectivity_mean",
        "connectivity_max", "cosine_mean", "best_context", "best_context_cell",
        "best_context_score", "connectivity_monocytic_mean", "connectivity_cns_mean",
        "connectivity_other_mean", "pubchem_inchikey", "pubchem_cid",
        "identity_source", "status", "connectivity_p", "fdr", "overlap_min", "overlap_median",
        "positive_context_fraction", "cell_balanced_mean", "n_cells_scored", "leave_one_cell_out_min", "context_robustness_status",
    ]
    frame = ranked.rename_axis("pert_id").reset_index()
    columns = [c for c in columns if c in frame.columns]
    clean = json.loads(frame[columns].head(limit).to_json(orient="records"))
    return clean


def evaluate_against_gold(
    ranked: pd.DataFrame,
    gold: pd.DataFrame,
    pert_name_map: dict[str, str],
    *,
    score_col: str = "connectivity_mean",
) -> dict[str, Any]:
    """Gold-standard evaluation on the per-perturbation aggregate.

    A gold cmap_name may map to several pert_ids (salts/batches); the name
    score is the mean of its perturbations' aggregate scores. Ties and
    missing perturbations are reported, never filled in.
    """
    rows = []
    missing = []
    for item in gold.itertuples(index=False):
        name = str(item.cmap_name).lower()
        pert_ids = [pid for pid, iname in pert_name_map.items() if iname.lower() == name]
        if not pert_ids:
            missing.append({"cmap_name": item.cmap_name, "reason": "missing_lincs_perturbation"})
            continue
        present = [pid for pid in pert_ids if pid in ranked.index]
        if not present:
            missing.append({"cmap_name": item.cmap_name, "reason": "no_context_in_harmonizome"})
            continue
        scores = ranked.loc[present, score_col].astype(float)
        if scores.isna().all():
            missing.append({"cmap_name": item.cmap_name, "reason": "insufficient_context_coverage"})
            continue
        rows.append({
            "cmap_name": item.cmap_name,
            "label": int(item.label),
            "label_class": item.label_class,
            "source": item.source,
            "drug_inchikey": item.drug_inchikey,
            "score": float(scores.mean()),
            "max_score": float(scores.max()),
        })
    matched = pd.DataFrame(rows)
    result = {
        "score_column": score_col,
        "matched_count": int(len(matched)),
        "missing": missing,
        "auroc": None,
        "auroc_max": None,
        "positive_median_score": None,
        "negative_median_score": None,
        "rows": rows,
    }
    if not matched.empty and matched["label"].nunique() == 2:
        labels = matched["label"].to_numpy(dtype=np.int64)
        result["auroc"] = _auroc(matched["score"].to_numpy(dtype=np.float64), labels)
        result["auroc_max"] = _auroc(matched["max_score"].to_numpy(dtype=np.float64), labels)
        pos = matched.loc[matched["label"] == 1, "score"]
        neg = matched.loc[matched["label"] == 0, "score"]
        result["positive_median_score"] = float(pos.median())
        result["negative_median_score"] = float(neg.median())
    return qualify_evaluation(result, gold)


def run_context_engine(
    *,
    names: Iterable[str] | None = None,
    n_perm: int = 2000,
    output: str | Path | None = DEFAULT_OUTPUT,
    disease_source: str = "mathys",
    min_hits: int = 3,
    model: str = "mix",
    seed: int = 0,
    disease_bundle: dict[str, Any] | None = None,
) -> dict[str, Any]:
    require(disease_source == "mathys", "context_engine_requires_mathys_source")
    require(n_perm > 0 and min_hits > 0, "invalid_run_configuration")
    disease = disease_bundle if disease_bundle is not None else build_disease_signature(source=disease_source, model=model)
    if disease_bundle is not None:
        from microglia_repurposing.replication import validate_signature
        validate_signature(disease)
    d = disease_vector(disease)
    universe = {item["gene_symbol"]: item["hgnc_id"] for item in disease["genes"]}
    frame, meta, load_report = load_harmonizome_contexts(universe)
    # Permute only genes represented in the assay, never absent genes as zeros.
    d = d.reindex(frame.index)
    pert_name_map = _pert_name_map()
    missing_names: list[dict[str, str]] = []
    if names is not None:
        wanted, missing_names = _resolve_names(names, pert_name_map)
        keep = [i for i, m in enumerate(meta) if m.pert_id in wanted]
        require(keep, "context_panel_empty_after_name_filter")
        frame = frame.iloc[:, keep]
        meta = [meta[i] for i in keep]
    ks = ContextKS(d, frame)
    connectivity = ks.scores(d)
    cosine = cosine_reversal(d, frame.T).reindex(frame.columns).to_numpy()
    hits = (frame.to_numpy() != 0).sum(axis=0)
    per_context = pd.DataFrame({
        "weighted_connectivity": connectivity,
        "cosine_reversal": cosine,
        "hits": hits,
    }, index=list(frame.columns))
    ranked = aggregate_context_scores(per_context, meta, min_hits=min_hits)
    scored = ranked[ranked["status"] == "scored"]
    if scored.empty:
        ranked["connectivity_p"] = np.nan
        ranked["fdr"] = np.nan
    else:
        scoring_mask = hits >= min_hits
        scoring_pert_ids = [m.pert_id for m, keep in zip(meta, scoring_mask) if keep]
        scoring_positions = np.flatnonzero(scoring_mask)
        p_values = connectivity_null(
            ks, d, scoring_pert_ids, scored["connectivity_mean"],
            n_perm=n_perm, seed=seed, positions=scoring_positions,
        )
        ranked["connectivity_p"] = p_values.reindex(ranked.index)
        ranked["fdr"] = fdr_from_pvalues(p_values).reindex(ranked.index)
    ranked["status_order"] = ranked["status"].map({"scored": 0, "insufficient_coverage": 1}).fillna(2)
    ranked = ranked.sort_values(
        ["status_order", "fdr", "connectivity_mean"],
        ascending=[True, True, False],
        na_position="last",
    )
    ranked = ranked.drop(columns=["status_order"])
    ranked["pert_iname"] = [
        pert_name_map.get(pid, LINCS_MISSING) for pid in ranked.index
    ]
    # Real compound identity from LINCS pert metadata where present.
    identity_map = _pert_identity_map()
    ranked["pubchem_inchikey"] = [identity_map.get(pid, {}).get("pubchem_inchikey") for pid in ranked.index]
    ranked["pubchem_cid"] = [identity_map.get(pid, {}).get("pubchem_cid") for pid in ranked.index]
    ranked["identity_source"] = ["lincs_pert_info" if identity_map.get(pid, {}).get("pubchem_inchikey") else "unknown" for pid in ranked.index]
    ranked = annotate_context_hits(ranked)
    # Canonical LINCS pert names where available; never invented.
    named = ranked[ranked["pert_iname"] != LINCS_MISSING].copy()
    dupes = named.groupby(named["pert_iname"].str.lower()).size()
    dupes = dupes[dupes > 1]
    gold = load_gold_standard()
    evaluation = evaluate_against_gold(ranked, gold, pert_name_map)
    # Best-context vectors for driver-gene explanation of the top hit.
    best_vectors = pd.DataFrame(0.0, index=list(ranked.index), columns=list(d.index))
    context_lookup = {}
    for i, m in enumerate(meta):
        context_lookup.setdefault((m.pert_id, m.cell_id, m.time, m.dose), i)
    for pid, row in ranked.iterrows():
        key = (pid, row.get("best_context_cell"), row.get("best_context_time"), row.get("best_context_dose"))
        position = context_lookup.get(key)
        if position is not None:
            best_vectors.loc[pid] = frame.iloc[:, position].to_numpy()
    # Packet candidacy: named perturbations only. A next-experiment packet for
    # an unnamed LINCS perturbation (-666) is not actionable; unnamed rows stay
    # in the ranked list and report, never in the packet.
    named = ranked[ranked["pert_iname"] != LINCS_MISSING].copy()
    packet_input = named
    packet = build_packet(disease, packet_input, best_vectors, d, mode="context")
    top_name = packet.get("next_experiment", {}) or {}
    if top_name.get("top_cmap_name"):
        pert_id = top_name["top_cmap_name"]
        packet["next_experiment"]["top_pert_iname"] = pert_name_map.get(pert_id, LINCS_MISSING)
        packet["next_experiment"]["top_best_context"] = ranked.loc[pert_id, "best_context"]
        packet["next_experiment"]["top_context_cell_class"] = None
        cell = ranked.loc[pert_id, "best_context_cell"]
        if isinstance(cell, str) and cell:
            from microglia_repurposing.contexts import cell_class
            packet["next_experiment"]["top_context_cell_class"] = cell_class(cell)
        packet["next_experiment"]["assay"] = "signature_reversal_in_declared_microglia_state"
        packet["next_experiment"]["falsification_links"] = [
            "no_reversal_of_declared_microglia_state_signature",
            "reversal_only_in_non_microglia_lincs_line",
        ]
    # Secondary triage: named hits with positive monocytic or CNS-lineage
    # context evidence (closer assay analogs to microglia). Hypotheses only.
    context_hits = []
    for pid, row in named.iterrows():
        if row.get("status") != "scored":
            continue
        fdr = row.get("fdr")
        if fdr is None or pd.isna(fdr) or float(fdr) > 0.25:
            continue
        monocytic = row.get("connectivity_monocytic_mean")
        cns = row.get("connectivity_cns_mean")
        positive = [c for c in (monocytic, cns) if c is not None and not pd.isna(c) and float(c) > 0]
        if not positive:
            continue
        context_hits.append({
            "pert_id": pid,
            "pert_iname": row["pert_iname"],
            "connectivity_mean": float(row["connectivity_mean"]),
            "connectivity_monocytic_mean": None if pd.isna(monocytic) else float(monocytic),
            "connectivity_cns_mean": None if pd.isna(cns) else float(cns),
            "best_context": row["best_context"],
            "fdr": float(fdr),
            "fdr_scope": "global mean; selected subgroup is descriptive and not separately tested",
        })
    context_hits.sort(key=lambda item: (-item["connectivity_mean"], item["pert_id"]))
    packet["cns_or_monocytic_evidence_hits"] = context_hits[:10]
    packet["context_hits_note"] = (
        "Named perturbations with FDR<=0.25 and positive reversal in monocytic "
        "(THP1/U937) or CNS-lineage (NPC/NEU/ASC) LINCS contexts. Hypotheses for "
        "a declared microglia experiment, not treatment candidates."
    )
    packet["content_hash"] = sha256_identity(
        {k: v for k, v in packet.items() if k != "content_hash"}
    )
    result = {
        "schema": RESULT_SCHEMA_CONTEXT,
        "research_planning_only": True,
        "claim_eligibility_mutation": False,
        "claim_boundary": CLAIM_BOUNDARY,
        "disease_signature_hash": disease["content_hash"],
        "disease": {
            "cell_context": disease["cell_context"],
            "cell_state": disease["cell_state"],
            "counts": disease["counts"],
            "source": disease["source"],
            "genes": disease["genes"],
            "l1000_note": disease.get("l1000_note"),
        },
        "context_panel": {
            "source": context_source_block(),
            "load_report": load_report,
            "panel": "named_gold_standard" if names is not None else "all_harmonizome_perturbations",
            "counts": {
                "contexts": len(meta),
                "perturbations": int(ranked.shape[0]),
                "perturbations_scored": int((ranked["status"] == "scored").sum()),
                "unresolved_names": missing_names,
            },
            "min_hits": min_hits,
            "n_perm": n_perm,
            "duplicate_pert_inames": sorted(dupes.index.tolist())[:20],
            "cell_class_note": (
                "monocytic=THP1,U937 (monocytic leukemia lines, closest LINCS analog "
                "to microglia); cns=NPC,NEU,ASC (CNS-lineage cells); other=rest. "
                "No LINCS context is a microglia assay."
            ),
        },
        "configuration": {"model": model, "min_hits": min_hits, "seed": seed, "n_perm": n_perm},
        "statistics": {"tested_statistic": "connectivity_mean", "method": "directional_set_ks/v2",
                       "p_resolution": 1 / (n_perm + 1), "fdr_family_size": len(scored),
                       "fdr_scope": "all coverage-eligible perturbation means", "triage_fdr_threshold": 0.25,
                       "gene_universe": list(d.index), "subgroup_and_best_context_tests": "not_performed",
                       "cell_balanced_and_leave_one_cell_out": "descriptive sensitivity only; not separately tested",
                       "validation_status": "exploratory"},
        "provenance": run_provenance(),
        "evaluation": evaluation,
        "ranking": json.loads(ranked.rename_axis("pert_id").reset_index().to_json(orient="records")),
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
    result["_meta"] = meta
    result["_frame"] = frame
    result["_per_context"] = per_context
    return result
