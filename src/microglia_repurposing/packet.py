"""Research-planning packet for one next microglia signature experiment."""
from __future__ import annotations

from typing import Any

from microglia_repurposing.contracts import CLAIM_BOUNDARY, PACKET_SCHEMA, sha256_identity
from microglia_repurposing.scoring import driver_genes


def build_packet(disease: dict[str, Any], ranked: Any, matrix: Any, disease_vec: Any, *, mode: str = "consensus") -> dict[str, Any]:
    context = mode == "context"
    drug_context = "lincs_l1000_chem_pert_instance" if context else "lincs_l1000_chem_pert_consensus"
    primary = "connectivity_mean" if context else "cosine_reversal"
    scored = ranked[ranked.get("status", "scored") == "scored"] if "status" in ranked.columns else ranked
    scored = scored[scored[primary] > 0] if primary in scored else scored
    significant = scored
    if "fdr" in scored.columns and not scored.empty:
        significant = scored[scored["fdr"].apply(lambda v: v is not None and not pd_na(v) and float(v) <= 0.25)]
    identity_links = []
    for gene in disease["genes"]:
        if gene.get("spine_entity_id"):
            identity_links.append({
                "entity_id": gene["spine_entity_id"],
                "kind": "gene",
                "label": gene["gene_symbol"],
                "hgnc_id": gene["hgnc_id"],
                "relation": "disease_signature_member",
                "epistemic_status": "hypothesis",
                "evidence_required": True,
            })
    if scored.empty or significant.empty:
        packet = {
            "schema": PACKET_SCHEMA,
            "research_planning_only": True,
            "claim_eligibility_mutation": False,
            "claim_boundary": CLAIM_BOUNDARY,
            "program": "neuro",
            "empty_reason": "no_eligible_positive_reversal" if scored.empty else "no_significant_reversal",
            "disease_cell_context": disease["cell_context"],
            "drug_cell_context": drug_context,
            "next_experiment": None,
            "claim_ids": [],
            "identity_links": identity_links,
            "blocking_reasons": (
                ["insufficient_overlap", "drug_signature_is_not_microglia", "core_dam_markers_mostly_outside_l1000_landmarks"]
                if scored.empty
                else ["no_significant_reversal", "drug_signature_is_not_microglia"]
            ),
        }
        packet["content_hash"] = sha256_identity({k: v for k, v in packet.items() if k != "content_hash"})
        return packet
    top_name = str(significant.index[0])
    top = significant.iloc[0]
    drivers = driver_genes(disease_vec, matrix.loc[top_name])
    packet = {
        "schema": PACKET_SCHEMA,
        "research_planning_only": True,
        "claim_eligibility_mutation": False,
        "claim_boundary": CLAIM_BOUNDARY,
        "program": "neuro",
        "disease_cell_context": disease["cell_context"],
        "disease_cell_state": disease["cell_state"],
        "drug_cell_context": drug_context,
        "next_experiment": {
            "assay": "signature_reversal_in_declared_microglia_state",
            "matrix": "human_iPSC_microglia_or_snRNAseq_DAM_state",
            "population": "research_planning_only",
            "query_signature_hash": disease["content_hash"],
            "top_cmap_name": top_name,
            "pert_id": top_name if context else None,
            "score": {"name": primary, "value": float(top[primary]),
                      "aggregation": "mean_over_coverage_eligible_contexts" if context else "consensus_signature",
                      "method": "directional_set_ks/v2" if context else "cosine_reversal/v2"},
            "cosine_reversal": nullable_float(top.get("cosine_mean" if context else "cosine_reversal")),
            "set_reversal": None if context else nullable_float(top.get("set_reversal")),
            "cmap_reversal": None,
            "tox_flag": bool(top.get("tox_flag", False)),
            "tox_annotation_status": top.get("tox_annotation_status", "mechanism_annotation_only"),
            "moa": top.get("moa", "unknown"),
            "pubchem_inchikey": top.get("pubchem_inchikey"),
            "pubchem_cid": top.get("pubchem_cid"),
            "evidence_status": "exploratory_hypothesis_not_validated",
            "coverage": {"overlap": nullable_float(top.get("overlap")),
                         "contexts_scored": nullable_float(top.get("n_contexts_scored"))},
            "driver_explanation": "best_context_only; significance tests the aggregate, not this selected context" if context else "consensus_signature",
            "fdr": None if pd_na(top.get("fdr")) else float(top.get("fdr")),
            "falsification_links": [
                "no_reversal_of_declared_microglia_state_signature",
                "signal_only_in_non_microglia_lincs_line",
            ],
            "blocking_reasons": _blockers(top, disease),
        },
        "claim_ids": [],
        "identity_links": identity_links,
        "driver_genes": drivers.reset_index().rename(columns={"index": "gene_symbol"}).head(12).to_dict(orient="records"),
    }
    packet["content_hash"] = sha256_identity({k: v for k, v in packet.items() if k != "content_hash"})
    return packet


def pd_na(value: Any) -> bool:
    try:
        import pandas as pd
        return value is None or bool(pd.isna(value))
    except Exception:
        return value is None


def _blockers(top: Any, disease: dict[str, Any]) -> list[str]:
    reasons = ["drug_signature_is_not_microglia", "no_independent_rescue_validation", "viability_not_tested"]
    if bool(top.get("tox_flag")):
        reasons.append("toxic_moa_flag")
    if pd_na(top.get("pubchem_inchikey")):
        reasons.append("pubchem_identity_missing")
    if disease["cell_context"] != "human_brain_microglia":
        reasons.append("disease_context_not_human_microglia")
    return reasons


def nullable_float(value: Any) -> float | None:
    return None if pd_na(value) else float(value)
