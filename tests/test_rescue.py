import copy
import json
from pathlib import Path

import pytest

from microglia_repurposing.rescue import (
    evaluate_rescue,
    ingest_observation,
    load_observations,
    load_panel,
    load_protocol,
    verify_document,
    workspace_status,
)

ROOT = Path(__file__).resolve().parents[1]


def test_frozen_panel_has_no_rescue_labels_and_null_auroc():
    status = workspace_status()
    report = status["report"]
    assert status["status"] == "frozen_unobserved"
    assert report["eligible_rescue_controls"] == 0
    assert report["rescue_auroc"] is None
    assert report["n_observed"] == 0
    protocol = load_protocol()
    assert protocol["assay"]["matrix"] == "human_iPSC_microglia"
    assert protocol["query_signature"]["counts"] == {"up": 12, "down": 10}
    panel = load_panel(protocol=protocol)
    assert all(e["benchmark_label"] is None and not e["eligible_for_rescue_benchmark"] for e in panel["entries"])
    assert any(e["compound_id"] == "sirolimus" and e["validation_split"] == "holdout" for e in panel["entries"])
    assert all(e["already_examined"] is True for e in panel["entries"] if e["validation_split"] == "development")


def test_panel_and_observation_hashes_fail_closed(tmp_path):
    protocol = json.loads((ROOT / "data/rescue/protocol.json").read_text())
    protocol["freeze_note"] = "tampered"
    with pytest.raises(ValueError, match="document_hash"):
        verify_document(protocol)
    observations = json.loads((ROOT / "data/rescue/observations.json").read_text())
    observations["records"] = [{"compound_id": "donepezil"}]
    path = tmp_path / "observations.json"
    path.write_text(json.dumps(observations))
    with pytest.raises(ValueError, match="document_hash"):
        load_observations(path)


def test_observation_requires_receipt_and_does_not_mutate_panel():
    protocol = load_protocol()
    panel = load_panel(protocol=protocol)
    log = load_observations(protocol=protocol)
    before = json.dumps(panel)
    with pytest.raises(ValueError, match="receipt"):
        ingest_observation(log, panel, {
            "compound_id": "donepezil", "protocol_hash": protocol["content_hash"],
            "assay_receipt": "  ", "outcome": "not_rescued", "viability_pass": True,
            "cell_model": "human_iPSC_microglia", "cell_state": "DAM_vs_homeostatic",
        })
    with pytest.raises(ValueError, match="cannot_set_eligibility"):
        ingest_observation(log, panel, {
            "compound_id": "donepezil", "protocol_hash": protocol["content_hash"],
            "assay_receipt": "assay-1", "outcome": "not_rescued", "viability_pass": True,
            "cell_model": "human_iPSC_microglia", "cell_state": "DAM_vs_homeostatic",
            "eligible_for_rescue_benchmark": True,
        })
    assert json.dumps(panel) == before


def test_unresolved_identity_and_viability_loss_are_not_eligible():
    protocol = load_protocol()
    panel = load_panel(protocol=protocol)
    log = load_observations(protocol=protocol)
    log = ingest_observation(log, panel, {
        "compound_id": "GW3965", "protocol_hash": protocol["content_hash"],
        "assay_receipt": "assay-gw", "outcome": "rescued", "viability_pass": True,
        "cell_model": "human_iPSC_microglia", "cell_state": "DAM_vs_homeostatic",
        "dam_reversal_score": 0.9,
    })
    log = ingest_observation(log, panel, {
        "compound_id": "vorinostat", "protocol_hash": protocol["content_hash"],
        "assay_receipt": "assay-vor", "outcome": "rescued", "viability_pass": False,
        "cell_model": "human_iPSC_microglia", "cell_state": "DAM_vs_homeostatic",
        "dam_reversal_score": 0.8,
    })
    log = ingest_observation(log, panel, {
        "compound_id": "vehicle_dmso", "protocol_hash": protocol["content_hash"],
        "assay_receipt": "assay-veh", "outcome": "not_rescued", "viability_pass": True,
        "cell_model": "human_iPSC_microglia", "cell_state": "DAM_vs_homeostatic",
        "dam_reversal_score": 0.0,
    })
    report = evaluate_rescue(protocol, panel, log)
    by_id = {row["compound_id"]: row for row in report["rows"]}
    assert by_id["GW3965"]["eligible_for_rescue_benchmark"] is False
    assert "identity_not_pinned" in by_id["GW3965"]["blocking_reasons"]
    assert by_id["vorinostat"]["eligible_for_rescue_benchmark"] is False
    assert "viability_gate_failed" in by_id["vorinostat"]["blocking_reasons"]
    assert by_id["vehicle_dmso"]["eligible_for_rescue_benchmark"] is False
    assert report["rescue_auroc"] is None


def test_self_reported_outcomes_do_not_become_labels_or_auroc():
    protocol = load_protocol()
    panel = load_panel(protocol=protocol)
    frozen = copy.deepcopy(panel)
    log = load_observations(protocol=protocol)
    log = ingest_observation(log, panel, {
        "compound_id": "sirolimus", "protocol_hash": protocol["content_hash"],
        "assay_receipt": "assay-sir", "outcome": "rescued", "viability_pass": True,
        "cell_model": "human_iPSC_microglia", "cell_state": "DAM_vs_homeostatic",
        "dam_reversal_score": 0.7,
    })
    log = ingest_observation(log, panel, {
        "compound_id": "donepezil", "protocol_hash": protocol["content_hash"],
        "assay_receipt": "assay-don", "outcome": "not_rescued", "viability_pass": True,
        "cell_model": "human_iPSC_microglia", "cell_state": "DAM_vs_homeostatic",
        "dam_reversal_score": -0.1,
    })
    report = evaluate_rescue(protocol, panel, log)
    assert report["eligible_rescue_controls"] == 0
    assert report["rescue_auroc"] is None
    assert panel == frozen
    workspace = workspace_status()
    assert workspace["report"]["rescue_auroc"] is None


def qualified_fixture(tmp_path, monkeypatch):
    """Synthetic engineering fixture, never written to the real observation log."""
    import microglia_repurposing.rescue as rescue
    p = copy.deepcopy(load_protocol())
    p['design'].update(status='qualified', induction_minimum_shift=0.3,
        homeostatic_minimum_restoration=0.2, minimum_independent_donors=2,
        minimum_replicates_per_donor=2, dose_time_plan='synthetic fixed condition',
        adjudication_plan='synthetic independent review')
    for arm in p['design']['arms'].values():
        arm.update(preparation='synthetic', concentration='synthetic', duration_h=24)
    p['design']['functional_endpoint'].update(method='synthetic clearance', minimum_effect=0.2)
    p['design']['signature_score'].update(formula='synthetic fixed formula', minimum_reversal=0.3)
    for g in p['query_signature']['genes']:
        g['provenance'].update(status='verified', source_doi='synthetic', source_locator='synthetic',
            source_species='human', human_transfer_justification='synthetic test only')
    rescue.seal(p)
    panel = copy.deepcopy(load_panel())
    panel['protocol_hash'] = p['content_hash']; rescue.seal(panel)
    log = rescue.seal(dict(schema=rescue.OBSERVATION_SCHEMA, protocol_hash=p['content_hash'], records=[]))
    monkeypatch.setattr(rescue, 'PACKAGE_ROOT', tmp_path)
    artifact = tmp_path/'data/rescue/fixture.json'; artifact.parent.mkdir(parents=True)
    artifact.write_text('{"synthetic": true}')
    evidence = dict(path='data/rescue/fixture.json', sha256=rescue.file_sha256(artifact))
    for donor in ('d1', 'd2'):
        for rep in ('r1', 'r2'):
            log = ingest_observation(log, panel, dict(compound_id='sirolimus', protocol_hash=p['content_hash'],
                assay_receipt='synthetic', outcome='rescued', viability_pass=True,
                cell_model='human_iPSC_microglia', cell_state='DAM_vs_homeostatic',
                observation_id=donor+rep, donor_id=donor, replicate_id=rep, batch_id='b1',
                dose='synthetic fixed dose', time_h=24, induction_arm='abeta_oligomer',
                dam_reversal_score=0.7, viability_loss_fraction=0.1, induction_shift=0.5,
                homeostatic_restoration=0.5, functional_effect=0.5, evidence=evidence))
            record = log['records'][-1]
            record['adjudication'] = dict(reviewer='synthetic reviewer', evidence=evidence,
                observation_hash=rescue.observation_identity(record), decision='rescued')
            rescue.seal(log)
    return p, panel, log, rescue


def test_complete_replicated_evidence_can_qualify_but_never_self_score_auroc(tmp_path, monkeypatch):
    p, panel, log, rescue = qualified_fixture(tmp_path, monkeypatch)
    report = evaluate_rescue(p, panel, log)
    assert report['eligible_rescue_controls'] == 1
    assert report['rescue_auroc'] is None
    assert report['n_observations'] == 4
    assert next(r for r in report['rows'] if r['compound_id']=='sirolimus')['n_donors'] == 2


@pytest.mark.parametrize('field,value,reason', [
    ('viability_loss_fraction',0.21,'viability_gate_failed'),
    ('functional_effect',0.1,'function_gate_failed'),
    ('induction_shift',0.1,'induction_gate_failed'),
    ('homeostatic_restoration',0.1,'homeostasis_gate_failed'),
    ('dam_reversal_score',0.1,'signature_gate_failed'),
    ('functional_effect',None,'missing_measurement:functional_effect'),
    ('dose','another dose','multiple_conditions_require_separate_analysis'),
    ('outcome','not_rescued','discordant_observations'),
])
def test_bad_earlier_well_cannot_be_hidden_by_latest_success(tmp_path, monkeypatch, field, value, reason):
    p, panel, log, rescue = qualified_fixture(tmp_path, monkeypatch)
    log['records'][0][field] = value
    rescue.seal(log)
    report = evaluate_rescue(p, panel, log)
    row = next(r for r in report['rows'] if r['compound_id']=='sirolimus')
    assert not row['eligible_for_rescue_benchmark']
    assert reason in row['blocking_reasons']
    assert len(row['observations']) == 4


def test_evidence_tamper_and_duplicate_and_wrong_matrix_fail_closed(tmp_path, monkeypatch):
    p, panel, log, rescue = qualified_fixture(tmp_path, monkeypatch)
    (tmp_path/'data/rescue/fixture.json').write_text('changed')
    assert evaluate_rescue(p,panel,log)['eligible_rescue_controls'] == 0
    duplicate = copy.deepcopy(log); duplicate['records'].append(copy.deepcopy(duplicate['records'][0]))
    rescue.seal(duplicate)
    with pytest.raises(ValueError, match='duplicate'): evaluate_rescue(p,panel,duplicate)
    log['records'][0]['cell_model'] = 'human_chondrocyte'; rescue.seal(log)
    with pytest.raises(ValueError, match='assay_mismatch'): evaluate_rescue(p,panel,log)


def test_single_donor_and_reused_replicate_do_not_establish_replication(tmp_path, monkeypatch):
    p, panel, log, rescue = qualified_fixture(tmp_path, monkeypatch)
    for record in log['records']:
        record['donor_id'] = 'd1'; record['replicate_id'] = 'r1'
    rescue.seal(log)
    row = next(r for r in evaluate_rescue(p,panel,log)['rows'] if r['compound_id']=='sirolimus')
    assert 'insufficient_independent_donors' in row['blocking_reasons']
    assert 'insufficient_replicates' in row['blocking_reasons']


def test_resealed_disk_log_still_requires_valid_record(tmp_path):
    from microglia_repurposing.rescue import seal
    log = load_observations(); log['records'] = [{'compound_id':'sirolimus'}]; seal(log)
    path=tmp_path/'bad.json'; path.write_text(json.dumps(log))
    with pytest.raises(ValueError, match='protocol_mismatch'): load_observations(path)
