"""Identity-pinned disease and drug signatures from public caches."""
from __future__ import annotations

from typing import Any, Iterable

import numpy as np
import pandas as pd

from microglia_repurposing.contracts import (
    CLAIM_BOUNDARY,
    DISEASE_SCHEMA,
    DRUG_SCHEMA,
    require,
    sha256_identity,
)
from microglia_repurposing.identity import pin_symbol, pin_symbols, spine_overlay
from microglia_repurposing.ingest import (
    load_cellmarker_microglia,
    load_landmark_genes,
    load_lincs_consensus,
    load_mathys_mic_pathology_de,
    load_mathys_mic_selective_de,
    load_mathys_subpopulation_markers,
    load_pert_info,
)


def _pin_signed(symbols: Iterable[str], sign: int, *, overlay: bool = False, weights: dict[str, float] | None = None) -> tuple[list[dict[str, Any]], list[str]]:
    pinned, missing = pin_symbols(symbols)
    official_weight: dict[str, float] = {}
    for original, value in (weights or {}).items():
        hit = pin_symbol(original)
        if hit is not None:
            official_weight[hit["gene_symbol"]] = float(value)
    rows = []
    for item in pinned:
        weight = official_weight.get(item["gene_symbol"], float(sign))
        row = {**item, "sign": 1 if weight > 0 else -1, "weight": float(weight), "spine_entity_id": None}
        if overlay:
            hit = spine_overlay(item["gene_symbol"])
            if hit is not None:
                row["spine_entity_id"] = hit.get("id")
        rows.append(row)
    return rows, missing


def build_disease_signature(*, source: str = "mathys", model: str = "mix",
                            contrast: str | None = None) -> dict[str, Any]:
    if source == "cellmarker":
        raw = load_cellmarker_microglia()
        weights = None
    elif source == "mathys":
        raw = load_mathys_mic_pathology_de(model=model)
        weights = raw["weights"]
    elif source == "mathys_selective":
        raw = load_mathys_mic_selective_de(model=model)
        weights = raw["weights"]
    elif source == "mathys_subpop":
        raw = load_mathys_subpopulation_markers(
            contrast=contrast or "mic1_vs_mic0")
        weights = raw["weights"]
    else:
        raise ValueError("disease_source_invalid")
    up_rows, up_missing = _pin_signed(raw["up_symbols"], 1, overlay=True, weights=weights)
    down_rows, down_missing = _pin_signed(raw["down_symbols"], -1, overlay=True, weights=weights)
    require(len(up_rows) >= 3 and len(down_rows) >= 3, "disease_signature_too_small")
    genes = up_rows + down_rows
    genes.sort(key=lambda item: item["gene_symbol"])
    source_block = {
        "library": raw["library"],
        "url": raw["source_url"],
        "sha256": raw["source_sha256"],
        "note": raw["note"],
    }
    if "source_doi" in raw:
        source_block["doi"] = raw["source_doi"]
        source_block["model"] = raw.get("model")
    if "terms" in raw:
        source_block["terms"] = raw["terms"]
    bundle = {
        "schema": DISEASE_SCHEMA,
        "research_planning_only": True,
        "claim_boundary": CLAIM_BOUNDARY,
        "cell_context": raw["cell_context"],
        "cell_state": raw["cell_state"],
        "species": raw["species"],
        "source": source_block,
        "genes": genes,
        "unmapped_symbols": sorted(set(up_missing + down_missing)),
        "ambiguous_symbols": raw.get("ambiguous_symbols") or [],
        "identity_overlay": {"attached": sum(bool(g["spine_entity_id"]) for g in genes),
                             "missing_or_unavailable": sum(not bool(g["spine_entity_id"]) for g in genes)},
        "counts": {"up": len(up_rows), "down": len(down_rows), "unmapped": len(set(up_missing + down_missing))},
        "l1000_note": "Query genes are HGNC-pinned Mathys microglia DE, not restricted to the 978 L1000 landmarks. LINCS consensus gene sets can include inferred genes.",
    }
    if "selectivity" in raw:
        bundle["selectivity"] = raw["selectivity"]
    bundle["content_hash"] = sha256_identity({k: v for k, v in bundle.items() if k != "content_hash"})
    return bundle


def disease_vector(bundle: dict[str, Any]) -> pd.Series:
    return pd.Series(
        {item["gene_symbol"]: float(item.get("weight", item["sign"])) for item in bundle["genes"]},
        dtype="float64",
        name="disease_sign",
    )


def build_drug_signatures(*, names: Iterable[str] | None = None) -> dict[str, Any]:
    consensus = load_lincs_consensus()
    pert = load_pert_info()
    landmarks = set(load_landmark_genes()["gene_symbol"].str.upper())
    wanted = None if names is None else {str(name).strip().lower() for name in names if str(name).strip()}
    pert_lookup = pert.groupby(pert["pert_iname"].str.lower())["pert_id"].agg(lambda v: sorted(set(v.astype(str)))).to_dict()
    records = []
    skipped = []
    for cmap_name, directions in consensus["drugs"].items():
        if wanted is not None and cmap_name.lower() not in wanted:
            continue
        up_rows, up_missing = _pin_signed(directions.get("up") or [], 1)
        down_rows, down_missing = _pin_signed(directions.get("down") or [], -1)
        if not up_rows and not down_rows:
            skipped.append({"cmap_name": cmap_name, "reason": "no_pinned_genes"})
            continue
        conflicts = {g["hgnc_id"] for g in up_rows} & {g["hgnc_id"] for g in down_rows}
        up_rows = [g for g in up_rows if g["hgnc_id"] not in conflicts]
        down_rows = [g for g in down_rows if g["hgnc_id"] not in conflicts]
        pert_ids = pert_lookup.get(cmap_name.lower(), [])
        records.append({
            "cmap_name": cmap_name,
            "pert_ids": pert_ids,
            "cell_context": consensus["cell_context"],
            "up": up_rows,
            "down": down_rows,
            "ambiguous_direction_hgnc": sorted(conflicts),
            "unmapped_symbols": sorted(set(up_missing + down_missing)),
            "landmark_overlap": sum(item["gene_symbol"].upper() in landmarks for item in up_rows + down_rows),
        })
    require(records, "drug_signatures_empty")
    records.sort(key=lambda item: item["cmap_name"].lower())
    bundle = {
        "schema": DRUG_SCHEMA,
        "research_planning_only": True,
        "claim_boundary": CLAIM_BOUNDARY,
        "cell_context": consensus["cell_context"],
        "source": {
            "library": consensus["library"],
            "url": consensus["source_url"],
            "sha256": consensus["source_sha256"],
            "note": "LINCS L1000 chemical-perturbation consensus gene sets. This is not a microglia assay.",
        },
        "records": records,
        "counts": {"drugs": len(records), "skipped": len(skipped)},
        "skipped": skipped[:50],
    }
    bundle["content_hash"] = sha256_identity({k: v for k, v in bundle.items() if k != "content_hash"})
    return bundle


def drug_matrix(bundle: dict[str, Any], universe: Iterable[str]) -> pd.DataFrame:
    columns = list(universe)
    index = [item["cmap_name"] for item in bundle["records"]]
    mat = np.zeros((len(index), len(columns)), dtype="float64")
    col_index = {gene: i for i, gene in enumerate(columns)}
    for row_i, item in enumerate(bundle["records"]):
        for gene in item["up"]:
            j = col_index.get(gene["gene_symbol"])
            if j is not None:
                mat[row_i, j] = 1.0
        for gene in item["down"]:
            j = col_index.get(gene["gene_symbol"])
            if j is not None:
                mat[row_i, j] = -1.0
    return pd.DataFrame(mat, index=index, columns=columns)
