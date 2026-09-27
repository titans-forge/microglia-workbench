"""Hash-pinned, research-planning-only contracts for signature reversal."""
from __future__ import annotations

import hashlib
import json
import re
from typing import Any, Mapping

DISEASE_SCHEMA = "titans-med-microglia-disease-signature/v1"
DRUG_SCHEMA = "titans-med-lincs-drug-signature/v1"
RESULT_SCHEMA = "titans-med-signature-reversal-result/v2"
RESULT_SCHEMA_CONTEXT = "titans-med-context-reversal-result/v2"
PACKET_SCHEMA = "titans-med-microglia-repurposing-packet/v2"
CLAIM_BOUNDARY = (
    "research-planning-only signature reversal; not diagnosis, dosing, "
    "pathogenicity, treatment, or a clinical ranking"
)
CELL_CONTEXTS = frozenset({
    "human_brain_microglia",
    "lincs_l1000_chem_pert_consensus",
    "lincs_l1000_chem_pert_instance",
    "synthetic_fixture",
})
UTC = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")
SHA = re.compile(r"^[0-9a-f]{64}$")
HGNC = re.compile(r"^HGNC:[0-9]+$")
PMC = re.compile(r"^PMC[0-9]+$")


def canonical_json_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False) + "\n").encode()


def sha256_identity(value: Any) -> str:
    raw = value if isinstance(value, (bytes, bytearray)) else canonical_json_bytes(value)
    return hashlib.sha256(raw).hexdigest()


def file_sha256(path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, code: str) -> None:
    if not condition:
        raise ValueError(code)


def validate_gene_row(row: Mapping[str, Any]) -> dict[str, Any]:
    require(isinstance(row, Mapping), "gene_row_invalid")
    require(isinstance(row.get("gene_symbol"), str) and row["gene_symbol"].strip(), "gene_symbol_invalid")
    require(isinstance(row.get("hgnc_id"), str) and HGNC.fullmatch(row["hgnc_id"]), "hgnc_id_invalid")
    require(row.get("status") == "Approved", "hgnc_status_not_approved")
    return dict(row)
