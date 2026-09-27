"""Declared DAM-state rescue assay. Unobserved until an assay receipt exists."""
from __future__ import annotations

import json
import math
import copy
from pathlib import Path
from typing import Any, Mapping


from microglia_repurposing.contracts import (
    CLAIM_BOUNDARY,
    file_sha256,
    require,
    sha256_identity,
    validate_gene_row,
)
from microglia_repurposing.identity import pin_symbol

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
RESCUE_DIR = PACKAGE_ROOT / "data" / "rescue" / "v2"
PROTOCOL_SCHEMA = "titans-med-dam-rescue-protocol/v1"
PANEL_SCHEMA = "titans-med-dam-rescue-panel/v1"
OBSERVATION_SCHEMA = "titans-med-dam-rescue-observation-log/v1"
REPORT_SCHEMA = "titans-med-dam-rescue-report/v1"
INCHIKEY = r"[A-Z0-9]{14}-[A-Z0-9]{10}-[A-Z0-9]"
ROLES = frozenset({
    "assay_vehicle",
    "assay_induction",
    "assay_viability_gate",
    "proposed_negative",
    "proposed_depletion_confounder",
    "proposed_cytotoxic_confounder",
    "legacy_not_rescue",
    "proposed_positive_hypothesis",
})
OUTCOMES = frozenset({"rescued", "not_rescued", "depleted", "cytotoxic", "inconclusive"})
ELIGIBLE_OUTCOMES = frozenset({"rescued", "not_rescued"})
RANKING_ROLES = frozenset({"proposed_negative", "proposed_positive_hypothesis"})


def verify_document(document: Mapping[str, Any]) -> dict[str, Any]:
    require(document.get("content_hash") == sha256_identity(
        {k: v for k, v in document.items() if k != "content_hash"}), "document_hash_mismatch")
    return dict(document)


def seal(document: dict[str, Any]) -> dict[str, Any]:
    document["content_hash"] = sha256_identity({k: v for k, v in document.items() if k != "content_hash"})
    return document


def _load_sealed(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"rescue_file_missing:{path.name}")
    return verify_document(json.loads(path.read_text(encoding="utf-8")))


def load_protocol(path: str | Path | None = None) -> dict[str, Any]:
    protocol = _load_sealed(Path(path) if path is not None else RESCUE_DIR / "protocol.json")
    return validate_protocol(protocol)


def validate_protocol(protocol):
    verify_document(protocol)
    require(protocol.get("schema") == PROTOCOL_SCHEMA, "rescue_protocol_schema_invalid")
    require(protocol.get("research_planning_only") is True, "rescue_protocol_claim_invalid")
    require(protocol.get("claim_boundary") == CLAIM_BOUNDARY, "rescue_protocol_boundary_invalid")
    require(protocol.get("assay", {}).get("matrix") == "human_iPSC_microglia", "rescue_matrix_invalid")
    require(protocol.get("assay", {}).get("cell_state") == "DAM_vs_homeostatic", "rescue_state_invalid")
    require(protocol.get("assay", {}).get("species") == "human", "rescue_species_invalid")
    genes = protocol["query_signature"]["genes"]
    require(len(genes) >= 6, "rescue_signature_too_small")
    seen = set()
    up = down = 0
    for gene in genes:
        validate_gene_row(gene)
        require(gene["hgnc_id"] not in seen, "rescue_gene_duplicate")
        seen.add(gene["hgnc_id"])
        require(gene["sign"] in (1, -1), "rescue_gene_sign_invalid")
        identity = pin_symbol(gene["gene_symbol"])
        require(identity is not None and identity["hgnc_id"] == gene["hgnc_id"], "rescue_gene_identity_mismatch")
        require(identity["gene_symbol"] == gene["gene_symbol"], "rescue_gene_symbol_mismatch")
        up += gene["sign"] == 1
        down += gene["sign"] == -1
    require(up >= 3 and down >= 3, "rescue_signature_unbalanced")
    require(protocol["query_signature"]["counts"] == {"up": up, "down": down}, "rescue_signature_counts_mismatch")
    primary = protocol["endpoints"]["primary"]
    require(primary["name"] == "DAM_signature_reversal" and primary["observed"] is False, "rescue_primary_not_unobserved")
    require(protocol["endpoints"]["viability"]["required"] is True, "rescue_viability_not_required")
    require(protocol["holdout"]["legacy_controls_eligible"] is False, "rescue_holdout_legacy_allowed")
    return protocol


def load_panel(path: str | Path | None = None, *, protocol: Mapping[str, Any] | None = None) -> dict[str, Any]:
    protocol = protocol or load_protocol()
    panel = _load_sealed(Path(path) if path is not None else RESCUE_DIR / "panel.json")
    return validate_panel(panel, protocol)


def validate_panel(panel, protocol):
    verify_document(panel)
    require(panel.get("schema") == PANEL_SCHEMA, "rescue_panel_schema_invalid")
    require(panel.get("protocol_hash") == protocol["content_hash"], "rescue_panel_protocol_mismatch")
    require(panel.get("research_planning_only") is True, "rescue_panel_claim_invalid")
    entries = panel["entries"]
    require(entries, "rescue_panel_empty")
    seen = set()
    for entry in entries:
        require(entry["compound_id"] not in seen, "rescue_panel_duplicate")
        seen.add(entry["compound_id"])
        require(entry["role"] in ROLES, "rescue_panel_role_invalid")
        require(entry["observed"] is False, "rescue_panel_observed_without_log")
        require(entry["benchmark_label"] is None, "rescue_panel_label_not_null")
        require(entry["eligible_for_rescue_benchmark"] is False, "rescue_panel_premature_eligibility")
        require(entry["validation_split"] in {"development", "holdout", "assay_control"}, "rescue_split_invalid")
        if entry["role"] in RANKING_ROLES and entry["validation_split"] == "development":
            require(entry["already_examined"] is True, "rescue_development_unexamined")
        if entry["validation_split"] == "holdout":
            require(entry["already_examined"] is False, "rescue_holdout_already_examined")
        if entry.get("pubchem_inchikey"):
            import re
            require(re.fullmatch(INCHIKEY, entry["pubchem_inchikey"]), "rescue_inchikey_invalid")
            require(entry["identity_status"] == "pinned", "rescue_identity_status_mismatch")
        else:
            require(entry["identity_status"] in {"unresolved", "assay_reagent", "biological_control"},
                    "rescue_identity_status_invalid")
            require(entry["eligible_for_rescue_benchmark"] is False, "unresolved_identity_not_eligible")
    require(panel["eligible_rescue_controls"] == 0, "rescue_panel_counts_nonzero")
    require(panel["rescue_auroc"] is None, "rescue_panel_auroc_not_null")
    return panel


def load_observations(path: str | Path | None = None, *, protocol: Mapping[str, Any] | None = None) -> dict[str, Any]:
    protocol = protocol or load_protocol()
    log = _load_sealed(Path(path) if path is not None else RESCUE_DIR / "observations.json")
    return validate_log(log, protocol["content_hash"])


def validate_log(log, protocol_hash):
    verify_document(log)
    require(log.get("schema") == OBSERVATION_SCHEMA, "rescue_observation_schema_invalid")
    require(log.get("protocol_hash") == protocol_hash, "rescue_observation_protocol_mismatch")
    require(isinstance(log.get("records"), list), "rescue_observation_records_invalid")
    ids = set()
    for record in log["records"]:
        validate_record(record, protocol_hash)
        if record.get("observation_id"):
            require(record["observation_id"] not in ids, "observation_duplicate")
            ids.add(record["observation_id"])
    return log


def validate_record(record, protocol_hash):
    require(isinstance(record, dict), "observation_record_invalid")
    require(record.get("protocol_hash") == protocol_hash, "observation_protocol_mismatch")
    require(record.get("cell_model") == "human_iPSC_microglia" and record.get("cell_state") == "DAM_vs_homeostatic", "observation_assay_mismatch")
    require(record.get("outcome") in OUTCOMES, "observation_outcome_invalid")
    require(type(record.get("viability_pass")) is bool, "observation_viability_invalid")
    for key in ("compound_id", "assay_receipt"):
        require(isinstance(record.get(key), str) and record[key].strip(), "observation_" + key + "_invalid")
    for key in ("observation_id", "donor_id", "replicate_id", "batch_id", "induction_arm", "dose"):
        if record.get(key) is not None:
            require(isinstance(record[key], str) and record[key].strip(), "observation_" + key + "_invalid")
    for key in ("time_h", "dam_reversal_score", "viability_loss_fraction", "induction_shift", "homeostatic_restoration", "functional_effect"):
        if record.get(key) is not None:
            require(_finite(record[key]), "observation_" + key + "_invalid")
    if record.get("time_h") is not None:
        require(record["time_h"] > 0, "observation_time_h_invalid")
    if record.get("viability_loss_fraction") is not None:
        require(record["viability_loss_fraction"] <= 1, "observation_viability_fraction_invalid")
    require(record.get("benchmark_label") is None and "eligible_for_rescue_benchmark" not in record, "observation_cannot_set_eligibility")


def ingest_observation(log: dict[str, Any], panel: Mapping[str, Any], observation: Mapping[str, Any]) -> dict[str, Any]:
    """Append one wet-lab record. Does not edit the frozen panel or confer a ranking label."""
    validate_log(log, panel["protocol_hash"])
    require(observation.get("protocol_hash") == log["protocol_hash"], "observation_protocol_mismatch")
    require(observation.get("protocol_hash") == panel["protocol_hash"], "observation_panel_mismatch")
    verify_document(panel)
    compound_id = observation.get("compound_id")
    entry = next((item for item in panel["entries"] if item["compound_id"] == compound_id), None)
    require(entry is not None, "observation_compound_unknown")
    require(isinstance(observation.get("assay_receipt"), str) and observation["assay_receipt"].strip(),
            "observation_receipt_missing")
    require(observation.get("outcome") in OUTCOMES, "observation_outcome_invalid")
    require(isinstance(observation.get("viability_pass"), bool), "observation_viability_missing")
    require(observation.get("cell_model") == "human_iPSC_microglia", "observation_matrix_mismatch")
    require(observation.get("cell_state") == "DAM_vs_homeostatic", "observation_state_mismatch")
    require(observation.get("benchmark_label") is None, "observation_cannot_set_benchmark_label")
    require("eligible_for_rescue_benchmark" not in observation, "observation_cannot_set_eligibility")
    record = {
        "compound_id": compound_id,
        "protocol_hash": observation["protocol_hash"],
        "assay_receipt": observation["assay_receipt"].strip(),
        "outcome": observation["outcome"],
        "viability_pass": observation["viability_pass"],
        "cell_model": observation["cell_model"],
        "cell_state": observation["cell_state"],
        "dose": observation.get("dose"),
        "time_h": observation.get("time_h"),
        "dam_reversal_score": observation.get("dam_reversal_score"),
        "note": observation.get("note") or "",
    }
    score = record["dam_reversal_score"]
    if score is not None:
        require(isinstance(score, (int, float)) and not isinstance(score, bool) and math.isfinite(score),
                "observation_score_invalid")
        record["dam_reversal_score"] = float(score)
    # Keep structured measurements and identifiers instead of silently dropping them.
    for key in ("observation_id", "donor_id", "replicate_id", "batch_id", "induction_arm",
                "viability_loss_fraction", "induction_shift", "homeostatic_restoration",
                "functional_effect", "evidence", "adjudication"):
        record[key] = copy.deepcopy(observation.get(key))
    if record["observation_id"] is not None:
        require(isinstance(record["observation_id"], str) and record["observation_id"].strip(), "observation_id_invalid")
        require(not any(r.get("observation_id") == record["observation_id"] for r in log["records"]), "observation_duplicate")
    for key in ("time_h", "viability_loss_fraction", "induction_shift", "homeostatic_restoration", "functional_effect"):
        if record[key] is not None:
            require(_finite(record[key]), f"observation_{key}_invalid")
    validate_record(record, log["protocol_hash"])
    updated = dict(log)
    updated["records"] = list(log["records"]) + [record]
    return seal(updated)


def _finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def protocol_blockers(protocol):
    reasons = []
    design = protocol.get("design", {})
    if design.get("status") != "qualified":
        reasons.append("assay_design_not_qualified")
    for key in ("induction_minimum_shift", "homeostatic_minimum_restoration",
                "maximum_viability_loss_fraction"):
        if not _finite(design.get(key)):
            reasons.append("design_missing:" + key)
    for key in ("minimum_independent_donors", "minimum_replicates_per_donor"):
        value = design.get(key)
        if type(value) is not int or value < 2:
            reasons.append("design_missing:" + key)
    for key in ("dose_time_plan", "adjudication_plan"):
        if not design.get(key): reasons.append("design_missing:" + key)
    for arm in ("abeta_oligomer", "lps_ifng"):
        spec = design.get("arms", {}).get(arm, {})
        if not all(spec.get(k) for k in ("preparation", "concentration", "duration_h")):
            reasons.append("induction_arm_unqualified:" + arm)
    for section, text_key, number_key in (("functional_endpoint", "method", "minimum_effect"),
                                           ("signature_score", "formula", "minimum_reversal")):
        spec = design.get(section, {})
        if not spec.get(text_key) or not _finite(spec.get(number_key)):
            reasons.append("design_missing:" + section)
    cap = design.get("maximum_viability_loss_fraction")
    if _finite(cap) and not 0 <= cap <= 0.2: reasons.append("viability_limit_invalid")
    functional = design.get("functional_endpoint", {})
    if functional.get("required") is not True or functional.get("beneficial_direction") != "increase":
        reasons.append("functional_endpoint_invalid")
    for gene in protocol["query_signature"]["genes"]:
        source = gene.get("provenance", {})
        if source.get("status") != "verified" or not all(source.get(k) for k in
                ("source_doi", "source_locator", "source_species", "human_transfer_justification")):
            reasons.append("gene_provenance_unverified:" + gene["gene_symbol"])
    return reasons


def _evidence_valid(evidence):
    # Receipts bind measurements to a real, local, immutable-by-hash artifact.
    if not isinstance(evidence, dict): return False
    try:
        path = (PACKAGE_ROOT / evidence["path"]).resolve()
        if not path.is_relative_to(PACKAGE_ROOT / "data" / "rescue"): return False
        return path.is_file() and file_sha256(path) == evidence["sha256"]
    except (KeyError, TypeError, OSError, ValueError):
        return False


def observation_identity(observation):
    return sha256_identity({k: v for k, v in observation.items() if k != "adjudication"})


def eligibility(entry, observation, protocol=None):
    reasons = []
    if observation is None: reasons.append("unobserved")
    if entry["identity_status"] != "pinned": reasons.append("identity_not_pinned")
    if entry["role"] not in RANKING_ROLES: reasons.append("role_not_ranking_label")
    if entry["validation_split"] == "assay_control": reasons.append("assay_control_not_ranking_label")
    if protocol is None: reasons.append("protocol_required")
    else: reasons.extend(protocol_blockers(protocol))
    if observation is not None:
        if observation.get("outcome") not in ELIGIBLE_OUTCOMES: reasons.append("outcome_not_rescued_or_not_rescued")
        if observation.get("viability_pass") is not True: reasons.append("viability_gate_failed")
        for key in ("observation_id", "donor_id", "replicate_id", "batch_id", "dose", "assay_receipt"):
            if not isinstance(observation.get(key), str) or not observation[key].strip():
                reasons.append("missing:" + key)
        for key in ("time_h", "dam_reversal_score", "viability_loss_fraction", "induction_shift",
                    "homeostatic_restoration", "functional_effect"):
            if not _finite(observation.get(key)): reasons.append("missing_measurement:" + key)
        if not _evidence_valid(observation.get("evidence")): reasons.append("measurement_evidence_unbound")
        review = observation.get("adjudication") or {}
        if not isinstance(review, dict) or not review.get("reviewer") or not _evidence_valid(review.get("evidence")):
            reasons.append("adjudication_evidence_unbound")
        if not isinstance(review, dict) or review.get("observation_hash") != observation_identity(observation) or review.get("decision") != observation.get("outcome"):
            reasons.append("adjudication_observation_mismatch")
        design = (protocol or {}).get("design", {})
        if observation.get("induction_arm") not in design.get("arms", {}): reasons.append("induction_arm_invalid")
        checks = [("viability_loss_fraction", design.get("maximum_viability_loss_fraction"), "maximum", "viability_gate_failed"),
                  ("induction_shift", design.get("induction_minimum_shift"), "minimum", "induction_gate_failed")]
        if observation.get("outcome") == "not_rescued":
            tests = [("dam_reversal_score", design.get("signature_score", {}).get("minimum_reversal")),
                     ("homeostatic_restoration", design.get("homeostatic_minimum_restoration")),
                     ("functional_effect", design.get("functional_endpoint", {}).get("minimum_effect"))]
            if all(_finite(observation.get(k)) and _finite(t) and observation[k] >= t for k,t in tests):
                reasons.append("negative_label_conflicts_with_measurements")
        if observation.get("outcome") == "rescued":
            checks += [("dam_reversal_score", design.get("signature_score", {}).get("minimum_reversal"), "minimum", "signature_gate_failed"),
                       ("homeostatic_restoration", design.get("homeostatic_minimum_restoration"), "minimum", "homeostasis_gate_failed"),
                       ("functional_effect", design.get("functional_endpoint", {}).get("minimum_effect"), "minimum", "function_gate_failed")]
        for key, threshold, direction, reason in checks:
            value = observation.get(key)
            if _finite(value) and _finite(threshold):
                if (direction == "maximum" and value > threshold) or (direction == "minimum" and value < threshold): reasons.append(reason)
    return {"eligible_for_rescue_benchmark": not reasons,
            "benchmark_label": (int(observation["outcome"] == "rescued") if not reasons else None),
            "blocking_reasons": sorted(set(reasons)), "observed": observation is not None}


def evaluate_rescue(protocol=None, panel=None, observations=None):
    protocol = load_protocol() if protocol is None else validate_protocol(protocol)
    panel = load_panel(protocol=protocol) if panel is None else validate_panel(panel, protocol)
    observations = load_observations(protocol=protocol) if observations is None else validate_log(observations, protocol["content_hash"])
    require(panel["protocol_hash"] == protocol["content_hash"], "rescue_panel_protocol_mismatch")
    require(observations["protocol_hash"] == protocol["content_hash"], "rescue_observation_protocol_mismatch")
    known = {e["compound_id"] for e in panel["entries"]}
    ids = set()
    for record in observations["records"]:
        require(record.get("compound_id") in known, "observation_compound_unknown")
        require(record.get("protocol_hash") == protocol["content_hash"], "observation_protocol_mismatch")
        require(record.get("cell_model") == "human_iPSC_microglia" and record.get("cell_state") == "DAM_vs_homeostatic", "observation_assay_mismatch")
        if record.get("observation_id"):
            require(record["observation_id"] not in ids, "observation_duplicate")
            ids.add(record["observation_id"])
    rows = []
    for entry in panel["entries"]:
        records = [r for r in observations["records"] if r["compound_id"] == entry["compound_id"]]
        statuses = [eligibility(entry, r, protocol) for r in records] or [eligibility(entry, None, protocol)]
        reasons = sorted({reason for status in statuses for reason in status["blocking_reasons"]})
        # Never let a later well hide an earlier failure or a different condition.
        conditions = {(r.get("induction_arm"), str(r.get("dose")), r.get("time_h")) for r in records}
        outcomes = {r.get("outcome") for r in records}
        if len(conditions) > 1: reasons.append("multiple_conditions_require_separate_analysis")
        if len(outcomes) > 1: reasons.append("discordant_observations")
        donors = {r.get("donor_id") for r in records if r.get("donor_id")}
        design = protocol.get("design", {})
        minimum = design.get("minimum_independent_donors")
        replicates = design.get("minimum_replicates_per_donor")
        if type(minimum) is not int or len(donors) < minimum: reasons.append("insufficient_independent_donors")
        if type(replicates) is not int or any(len({r.get("replicate_id") for r in records if r.get("donor_id") == donor and r.get("replicate_id")}) < replicates for donor in donors):
            reasons.append("insufficient_replicates")
        eligible = bool(records) and not reasons
        rows.append({**{k: entry.get(k) for k in ("compound_id", "name", "role", "validation_split", "identity_status", "pubchem_inchikey")},
            "observed": bool(records), "n_observations": len(records), "n_donors": len(donors),
            "observations": records, "outcome": next(iter(outcomes)) if len(outcomes) == 1 else None,
            "viability_pass": all(r.get("viability_pass") is True for r in records) if records else None,
            "dam_reversal_score": None, "benchmark_label": int(outcomes == {"rescued"}) if eligible else None,
            "eligible_for_rescue_benchmark": eligible, "blocking_reasons": sorted(set(reasons))})
    return seal({"schema": REPORT_SCHEMA, "research_planning_only": True, "claim_boundary": CLAIM_BOUNDARY,
        "protocol_hash": protocol["content_hash"], "panel_hash": panel["content_hash"], "observation_hash": observations["content_hash"],
        "n_entries": len(rows), "n_observed": sum(r["observed"] for r in rows), "n_observations": len(observations["records"]),
        "eligible_rescue_controls": sum(r["eligible_for_rescue_benchmark"] for r in rows),
        "rescue_auroc": None, "benchmark_status": "unavailable_no_registered_independent_predictions",
        "holdout_status": "Prospective assay hypotheses only; development and holdout are not pooled into performance metrics.",
        "design_blockers": protocol_blockers(protocol), "rows": rows,
        "limitations": ["Signature reversal alone is not functional rescue.",
            "Assay qualification and gene-level human transfer remain explicit requirements.",
            "All observations are retained; mixed conditions or outcomes require separate adjudication.",
            "Prediction AUROC is disabled until a separately frozen prediction benchmark exists; assay scores are never substituted."]})


def workspace_status() -> dict[str, Any]:
    protocol = load_protocol()
    panel = load_panel(protocol=protocol)
    observations = load_observations(protocol=protocol)
    report = evaluate_rescue(protocol, panel, observations)
    return {
        "status": "frozen_unobserved" if report["n_observed"] == 0 else "observations_present",
        "protocol_hash": protocol["content_hash"],
        "panel_hash": panel["content_hash"],
        "report": report,
        "protocol": {k: protocol[k] for k in (
            "schema", "revision", "purpose", "assay", "design", "endpoints", "falsifiers", "holdout", "freeze_note") if k in protocol},
        "query_signature": protocol["query_signature"],
        "panel": {k: panel[k] for k in (
            "schema", "eligible_rescue_controls", "rescue_auroc", "holdout_status", "entries") if k in panel},
        "source_hashes": {
            "protocol.json": file_sha256(RESCUE_DIR / "protocol.json"),
            "panel.json": file_sha256(RESCUE_DIR / "panel.json"),
            "observations.json": file_sha256(RESCUE_DIR / "observations.json"),
        },
    }
