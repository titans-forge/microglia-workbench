#!/usr/bin/env python3
"""Emit the frozen DAM rescue artifacts. Run once; later edits need a new protocol."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from microglia_repurposing.contracts import CLAIM_BOUNDARY
from microglia_repurposing.rescue import seal

OUT = ROOT / "data" / "rescue"
if any((OUT / name).exists() for name in ("protocol.json", "panel.json", "observations.json")):
    raise SystemExit("Frozen v1 artifacts already exist. Create a new revision; never reset observations.")
OUT.mkdir(parents=True, exist_ok=True)

DAM_UP = [
    ("APOE", "HGNC:613"), ("TREM2", "HGNC:17761"), ("TYROBP", "HGNC:12449"),
    ("SPP1", "HGNC:11255"), ("ITGAX", "HGNC:6152"), ("CLEC7A", "HGNC:14558"),
    ("CST7", "HGNC:2479"), ("CD9", "HGNC:1709"), ("LPL", "HGNC:6677"),
    ("GPNMB", "HGNC:4462"), ("CD68", "HGNC:1693"), ("CTSD", "HGNC:2529"),
]
DAM_DOWN = [
    ("P2RY12", "HGNC:18124"), ("TMEM119", "HGNC:27884"), ("CX3CR1", "HGNC:2558"),
    ("P2RY13", "HGNC:4537"), ("SELPLG", "HGNC:10722"), ("GPR34", "HGNC:4490"),
    ("HEXB", "HGNC:4879"), ("OLFML3", "HGNC:24956"), ("TGFBR1", "HGNC:11772"),
    ("SALL1", "HGNC:10524"),
]
genes = (
    [{"gene_symbol": s, "hgnc_id": h, "sign": 1, "status": "Approved"} for s, h in DAM_UP]
    + [{"gene_symbol": s, "hgnc_id": h, "sign": -1, "status": "Approved"} for s, h in DAM_DOWN]
)

protocol = seal({
    "schema": "titans-med-dam-rescue-protocol/v1",
    "created_at": "2026-09-10T14:00:00Z",
    "research_planning_only": True,
    "claim_boundary": CLAIM_BOUNDARY,
    "purpose": "Prespecified human iPSC-microglia DAM-state rescue assay; not a LINCS ranking and not a clinical trial.",
    "assay": {
        "name": "DAM_state_rescue_in_human_iPSC_microglia",
        "species": "human",
        "matrix": "human_iPSC_microglia",
        "cell_state": "DAM_vs_homeostatic",
        "induction": "Declared DAM-like challenge (Aβ oligomers or LPS/IFN-γ) versus matched vehicle in the same donor line.",
        "note": "LINCS THP1/U937/NPC/NEU lines are not this assay. Mathys pathology-vs-no-pathology DE is a related but distinct contrast.",
    },
    "query_signature": {
        "name": "core_DAM_vs_homeostatic",
        "source": {
            "library": "Keren-Shaul 2017 DAM core intersected with human microglia markers",
            "doi": "10.1038/nature22999",
            "supporting_doi": ["10.1016/j.immuni.2017.08.008", "10.1038/s41586-019-1195-2"],
            "note": "Frozen gene identities, not a newly fit mixed model. Human official symbols only.",
        },
        "genes": genes,
        "counts": {"up": len(DAM_UP), "down": len(DAM_DOWN)},
    },
    "endpoints": {
        "primary": {
            "name": "DAM_signature_reversal",
            "description": "Signed reversal of the frozen DAM vs homeostatic query on iPSC-microglia RNA-seq after induction.",
            "observed": False,
        },
        "viability": {
            "required": True,
            "description": "Cell count or ATP assay. Viability loss >20% vs vehicle blocks a rescued label.",
            "observed": False,
        },
        "secondary": [
            {"name": "homeostatic_marker_restoration", "examples": ["P2RY12", "TMEM119"], "observed": False},
            {"name": "phagocytosis", "observed": False, "required": False},
        ],
    },
    "falsifiers": [
        "reversal_only_with_viability_loss",
        "induction_control_fails_to_shift_DAM_signature",
        "homeostatic_genes_not_restored",
        "signal_only_in_non_microglia_lincs_line",
        "label_from_clinical_trial_or_depletion_paper",
    ],
    "holdout": {
        "legacy_controls_eligible": False,
        "rule": "Compounds already examined in control_panel_ad.csv remain development references and cannot be relabeled as a blinded holdout. Holdout hypotheses were reserved before any iPSC-microglia observations.",
    },
    "freeze_note": "Frozen before any wet-lab DAM rescue observation. Literature roles are hypotheses. A changed assay or gene set requires a new protocol hash.",
})

def entry(**kwargs):
    row = {
        "observed": False,
        "benchmark_label": None,
        "eligible_for_rescue_benchmark": False,
        "already_examined": False,
        "pubchem_cid": None,
        "pubchem_inchikey": None,
        "identity_status": "unresolved",
        "reference_doi": None,
        "finding": None,
        "boundary": "No DAM-state rescue observation in human iPSC-microglia has been recorded for this protocol.",
    }
    row.update(kwargs)
    return row

panel = seal({
    "schema": "titans-med-dam-rescue-panel/v1",
    "protocol_hash": protocol["content_hash"],
    "research_planning_only": True,
    "claim_boundary": CLAIM_BOUNDARY,
    "eligible_rescue_controls": 0,
    "rescue_auroc": None,
    "holdout_status": "No iPSC-microglia observations yet. Legacy exploratory compounds cannot become a blinded holdout. Holdout hypotheses remain unobserved.",
    "entries": [
        entry(compound_id="vehicle_dmso", name="Matched DMSO vehicle", role="assay_vehicle",
              validation_split="assay_control", identity_status="assay_reagent",
              finding="Assay baseline. Not a ranking label.",
              boundary="Vehicle defines the comparison, it does not adjudicate drug rescue."),
        entry(compound_id="lps_ifng_induction", name="LPS plus IFN-γ induction", role="assay_induction",
              validation_split="assay_control", identity_status="biological_control",
              finding="Expected to shift the culture toward a DAM-like signature if the assay is working.",
              boundary="Induction success is a quality-control outcome, not a drug-rescue positive."),
        entry(compound_id="abeta_oligomer_induction", name="Aβ oligomer induction", role="assay_induction",
              validation_split="assay_control", identity_status="biological_control",
              finding="AD-relevant DAM-like challenge. Exact oligomer prep must be recorded on the assay receipt.",
              boundary="Induction success is a quality-control outcome, not a drug-rescue positive."),
        entry(compound_id="viability_gate", name="Viability and cell-count readout", role="assay_viability_gate",
              validation_split="assay_control", identity_status="assay_reagent",
              finding="Required gate. Cell loss can masquerade as signature reversal.",
              boundary="A rescued label is forbidden when the viability gate fails."),
        entry(compound_id="donepezil", name="donepezil", role="proposed_negative",
              validation_split="development", already_examined=True, identity_status="pinned",
              pubchem_cid="3152", pubchem_inchikey="ADEBPBSSDYVVLD-UHFFFAOYSA-N",
              reference_doi=None, cmap_name="donepezil",
              finding="Symptomatic AChE inhibitor. Proposed not to reverse a DAM-state signature in iPSC-microglia.",
              boundary="Clinical use is not a negative rescue observation until this assay is run."),
        entry(compound_id="memantine", name="memantine", role="proposed_negative",
              validation_split="development", already_examined=True, identity_status="pinned",
              pubchem_cid="4054", pubchem_inchikey="BUGYDGFZZOZRHP-UHFFFAOYSA-N",
              reference_doi=None, cmap_name="memantine",
              finding="Symptomatic NMDAR antagonist. Proposed not to reverse a DAM-state signature in iPSC-microglia.",
              boundary="Clinical use is not a negative rescue observation until this assay is run."),
        entry(compound_id="pexidartinib", name="pexidartinib", role="proposed_depletion_confounder",
              validation_split="development", already_examined=True, identity_status="pinned",
              pubchem_cid="25151352", pubchem_inchikey="JGWRKYUXBBNENE-UHFFFAOYSA-N",
              reference_doi="10.1016/j.neuron.2014.02.040", cmap_name="pexidartinib",
              finding="CSF1R inhibition depletes microglia in mice. Expected viability loss or cell dropout, not functional rescue.",
              boundary="Depletion is a competing explanation and cannot become a rescued label."),
        entry(compound_id="vorinostat", name="vorinostat", role="proposed_cytotoxic_confounder",
              validation_split="development", already_examined=True, identity_status="pinned",
              pubchem_cid="5311", pubchem_inchikey="WAEXFXRVDQXREF-UHFFFAOYSA-N",
              reference_doi=None, cmap_name="vorinostat",
              finding="HDAC inhibitor often ranks in CMap via stress/cytotoxicity. Not a DAM rescue positive.",
              boundary="A cytotoxic signature is not rescue. Viability gate required."),
        entry(compound_id="minocycline", name="minocycline", role="legacy_not_rescue",
              validation_split="development", already_examined=True, identity_status="pinned",
              pubchem_cid="54675783", pubchem_inchikey="FFTVPQUHLQBXQZ-KVUCHLLUSA-N",
              reference_doi=None, cmap_name="minocycline",
              finding="MADE did not delay impairment in mild Alzheimer disease. Clinical outcome is not this assay.",
              boundary="Cannot be promoted to a rescue positive or blinded holdout from the legacy label."),
        entry(compound_id="pioglitazone", name="pioglitazone", role="legacy_not_rescue",
              validation_split="development", already_examined=True, identity_status="pinned",
              pubchem_cid="4829", pubchem_inchikey="HYAFETHFCAUJAY-UHFFFAOYSA-N",
              reference_doi="10.1016/S1474-4422(21)00043-0", cmap_name="pioglitazone",
              finding="TOMMORROW tested delay of cognitive impairment, not microglia DAM rescue.",
              boundary="Trial participation does not establish a rescue benchmark label."),
        entry(compound_id="sirolimus", name="sirolimus", role="proposed_positive_hypothesis",
              validation_split="holdout", already_examined=False, identity_status="pinned",
              pubchem_cid="5284616", pubchem_inchikey="QFJCIRLUMZQUOT-HPLJOQBZSA-N",
              cmap_name="sirolimus",
              finding="mTOR inhibition is a literature hypothesis for immune/microglia state change. Unobserved in this assay.",
              boundary="Holdout hypothesis reserved before iPSC-microglia data. LINCS reversal is not adjudication."),
        entry(compound_id="GW3965", name="GW3965", role="proposed_positive_hypothesis",
              validation_split="holdout", already_examined=False, identity_status="unresolved",
              finding="LXR agonism is a literature hypothesis for DAM/lipid-microglia biology. Compound identity is not pinned in the local PubChem cache.",
              boundary="Unresolved identity cannot become an eligible rescue control."),
        entry(compound_id="MCC950", name="MCC950", role="proposed_positive_hypothesis",
              validation_split="holdout", already_examined=False, identity_status="unresolved",
              finding="NLRP3 inhibition is a literature inflammasome hypothesis. Compound identity is not pinned in the local PubChem cache.",
              boundary="Unresolved identity cannot become an eligible rescue control."),
    ],
})

observations = seal({
    "schema": "titans-med-dam-rescue-observation-log/v1",
    "protocol_hash": protocol["content_hash"],
    "research_planning_only": True,
    "records": [],
})

for name, document in {"protocol.json": protocol, "panel.json": panel, "observations.json": observations}.items():
    path = OUT / name
    path.write_text(json.dumps(document, indent=2, ensure_ascii=True, allow_nan=False) + "\n", encoding="utf-8")
    print(name, document["content_hash"])
print("panel_protocol", panel["protocol_hash"] == protocol["content_hash"])
