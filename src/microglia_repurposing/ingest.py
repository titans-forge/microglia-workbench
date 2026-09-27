"""Bounded ingest of public, already-cached signature sources."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from microglia_repurposing.contracts import file_sha256, require

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
CACHE = PACKAGE_ROOT / "data" / "cache"

CELLMARKER_HOMEOSTATIC = "Homeostatic Microglial Cell Brain Human"
CELLMARKER_ENTORHINAL = "Microglial Cell Entorhinal Cortex Human"
CELLMARKER_SFG = "Microglial Cell Superior Frontal Gyrus Human"
LINCS_CONSENSUS = "LINCS_L1000_Chem_Pert_Consensus_Sigs.gmt"
LINCS_UP = "LINCS_L1000_Chem_Pert_up.gmt"
LINCS_DOWN = "LINCS_L1000_Chem_Pert_down.gmt"
CELLMARKER_GMT = "CellMarker_2024.gmt"
LANDMARK = "GSE92742_Broad_LINCS_gene_info_delta_landmark.txt.gz"
PERT_INFO = "GSE92742_Broad_LINCS_pert_info.txt.gz"


def cache_file(name: str) -> Path:
    path = CACHE / name
    require(path.is_file() and path.stat().st_size > 0, f"cache_missing:{name}")
    from microglia_repurposing.provenance import verify_source
    verify_source(path)
    return path


def parse_gmt(path: Path) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        require(len(parts) >= 3 and bool(parts[0]), "gmt_row_invalid")
        name = parts[0]
        require(name not in sets, f"gmt_duplicate_term:{name}")
        genes = []
        seen = set()
        for item in parts[2:]:
            symbol = item.split(",")[0].strip()
            if symbol and symbol not in seen:
                seen.add(symbol)
                genes.append(symbol)
        sets[name] = genes
    require(sets, "gmt_empty")
    return sets


MATHYS_MIC_DE = "mathys2019_mic_pathology_de.tsv"
MATHYS_MIC_XLSX = "mathys2019_MOESM4.xlsx"
MATHYS_URL = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1195-2/MediaObjects/41586_2019_1195_MOESM4_ESM.xlsx"
MATHYS_DOI = "10.1038/s41586-019-1195-2"

# MOESM4 holds the same three contrasts side by side for six cell types. Block
# offsets are positional; the expected header is asserted so a changed sheet
# layout fails closed instead of silently reading the wrong contrast.
MATHYS_CELL_TYPES = ("Ex", "In", "Ast", "Oli", "Opc", "Mic")
MATHYS_MICROGLIA_SHEET = "Mic"
MATHYS_CONTRASTS = {
    "no_pathology_vs_pathology": 0,
    "no_pathology_vs_early_pathology": 11,
    "early_vs_late_pathology": 22,
}
MATHYS_BLOCK_COLUMNS = (
    "gene_symbol", "ind_adj_pval", "no_pathology_mean", "pathology_mean",
    "log2fc", "mixed_z", "mixed_p", "deg_ind", "deg_ind_mix",
)
MATHYS_BLOCK_HEADER = (
    "IndModel.adj.pvals", "no.pathology.mean", "pathology.mean", "IndModel.FC",
    "MixedModel.z", "MixedModel.p", "DEGs.Ind.Model", "DEGs.Ind.Mix.models",
)
MATHYS_MODEL_FLAG = {"ind": "deg_ind", "mix": "deg_ind_mix"}


def load_mathys_celltype_contrast(cell_type: str, *, model: str = "mix",
                                  contrast: str = "no_pathology_vs_pathology") -> pd.DataFrame:
    """One cell type's block for one declared MOESM4 contrast.

    Columns are read positionally and the header is asserted, so a publisher
    layout change raises rather than returning a different contrast.
    """
    require(cell_type in MATHYS_CELL_TYPES, f"mathys_cell_type_unknown:{cell_type}")
    require(model in MATHYS_MODEL_FLAG, "mathys_model_invalid")
    require(contrast in MATHYS_CONTRASTS, f"mathys_contrast_unknown:{contrast}")
    path = cache_file(MATHYS_MIC_XLSX)
    start = MATHYS_CONTRASTS[contrast]
    raw = pd.read_excel(path, sheet_name=cell_type, header=None)
    require(raw.shape[1] >= start + len(MATHYS_BLOCK_COLUMNS), "mathys_sheet_too_narrow")
    header = [str(value).strip() for value in raw.iloc[1, start + 1:start + len(MATHYS_BLOCK_COLUMNS)]]
    require(tuple(header) == MATHYS_BLOCK_HEADER, f"mathys_block_header_invalid:{cell_type}:{contrast}")
    block = raw.iloc[2:, start:start + len(MATHYS_BLOCK_COLUMNS)].copy()
    block.columns = list(MATHYS_BLOCK_COLUMNS)
    block = block.dropna(subset=["gene_symbol"])
    block["gene_symbol"] = block["gene_symbol"].astype(str).str.strip()
    block = block[block["gene_symbol"].ne("") & block["gene_symbol"].ne("nan")]
    block = block.drop_duplicates("gene_symbol", keep="first")
    for column in ("ind_adj_pval", "no_pathology_mean", "pathology_mean", "log2fc", "mixed_z", "mixed_p", "deg_ind", "deg_ind_mix"):
        block[column] = pd.to_numeric(block[column], errors="coerce")
    require(not block["log2fc"].isna().all(), f"mathys_effect_missing:{cell_type}")
    return block.reset_index(drop=True)


def load_mathys_mic_selective_de(*, model: str = "mix",
                                 contrast: str = "no_pathology_vs_pathology") -> dict[str, Any]:
    """Microglia DEGs for one contrast that are not DEGs in any other MOESM4 cell type.

    The rule is declared before scoring: a gene enters the selective signature
    only when the same contrast and the same model flag call it differentially
    expressed in microglia and in none of Ex, In, Ast, Oli, Opc. Excluded genes
    and per-cell-type overlaps are reported, never silently dropped.
    """
    require(model in MATHYS_MODEL_FLAG, "mathys_model_invalid")
    flag = MATHYS_MODEL_FLAG[model]
    path = cache_file(MATHYS_MIC_XLSX)
    contrasts = {cell: load_mathys_celltype_contrast(cell, model=model, contrast=contrast) for cell in MATHYS_CELL_TYPES}
    mic = contrasts[MATHYS_MICROGLIA_SHEET]
    mic_degs = set(mic.loc[mic[flag] == 1, "gene_symbol"])
    require(mic_degs, "mathys_selective_microglia_degs_empty")
    overlaps: dict[str, int] = {}
    shared_any: set[str] = set()
    for cell in MATHYS_CELL_TYPES:
        if cell == MATHYS_MICROGLIA_SHEET:
            continue
        other = contrasts[cell]
        other_degs = set(other.loc[other[flag] == 1, "gene_symbol"])
        overlaps[cell] = len(mic_degs & other_degs)
        shared_any |= (mic_degs & other_degs)
    selective = mic_degs - shared_any
    require(selective, "mathys_selective_signature_empty")
    rows = mic[mic["gene_symbol"].isin(selective)]
    up = rows[rows["log2fc"] > 0].sort_values("log2fc", ascending=False)
    down = rows[rows["log2fc"] < 0].sort_values("log2fc")
    return {
        "library": f"Mathys2019_Nature_MOESM4_{MATHYS_MICROGLIA_SHEET}_selective",
        "source_url": MATHYS_URL,
        "source_doi": MATHYS_DOI,
        "source_sha256": file_sha256(path),
        "model": model,
        "flag_column": flag,
        "contrast": contrast,
        "up_symbols": list(up["gene_symbol"].astype(str)),
        "down_symbols": list(down["gene_symbol"].astype(str)),
        "weights": {str(row.gene_symbol): float(row.log2fc) for row in rows.itertuples(index=False)},
        "ambiguous_symbols": [],
        "cell_context": "human_brain_microglia",
        "cell_state": f"{contrast}_microglia_selective",
        "species": "human",
        "brain_region": "prefrontal_cortex",
        "selectivity": {
            "rule": (
                "microglia DEGs for this contrast minus the union of DEGs called by the same "
                "model flag in Ex, In, Ast, Oli and Opc; identity resolution happens after the "
                "rule, so the selective set is fixed before any LINCS coverage is inspected"
            ),
            "contrast": contrast,
            "model": model,
            "flag_column": flag,
            "microglia_degs": len(mic_degs),
            "shared_with_any_other_cell_type": sorted(shared_any),
            "shared_counts": overlaps,
            "selective_genes": sorted(selective),
            "not_selectivity": (
                "Selectivity here is a DEG-calling filter on one dataset, not evidence of "
                "microglia-specific function, causality, or disease driving. A gene called "
                "selective can still act in other cell types."
            ),
        },
        "note": (
            "Mathys et al. 2019 Nature 570:332-337, MOESM4 Mic sheet, no-pathology vs "
            "pathology, restricted to genes not called as DEGs in the other five MOESM4 "
            "cell types for the same contrast. ROSMAP snRNA-seq microglia; not iPSC-microglia "
            "and not a LINCS cell line. Selectivity costs coverage: fewer genes reach the "
            "LINCS panel than in the unrestricted signature."
        ),
    }


MATHYS_SUBPOP_XLSX = "mathys2019_MOESM8.xlsx"
MATHYS_SUBPOP_SHEET = "Subpopulation_markers"
MATHYS_SUBPOP_COLUMNS = ("gene_symbol", "adj_pval", "lfc", "subpopulation")
MATHYS_SUBPOP_HEADER = ("gene.name", "adj.pvals", "lFC", "subpopulation")
# Declared subpopulation contrasts. Mic1 is the disease-associated microglial
# subpopulation reported by the authors; Mic0 is the homeostatic one. Markers
# are over-expression relative to the other subpopulations of the same type,
# so this is a cell-state identity contrast, not a case-control disease contrast.
MATHYS_SUBPOP_CONTRASTS = {
    "mic1_vs_mic0": ("Mic1", "Mic0"),
    "mic3_vs_mic0": ("Mic3", "Mic0"),
    "mic2_vs_mic0": ("Mic2", "Mic0"),
}


def load_mathys_subpopulation_markers(*, contrast: str = "mic1_vs_mic0") -> dict[str, Any]:
    """Microglial subpopulation markers for one declared MOESM8 contrast.

    The up direction is the first subpopulation of the declared pair and the
    down direction is the second. A gene marked in both is dropped as ambiguous
    and reported, never forced into a direction. Non-finite log fold changes
    (the publisher records some as inf) are dropped and counted.
    """
    require(contrast in MATHYS_SUBPOP_CONTRASTS, f"mathys_subpop_contrast_unknown:{contrast}")
    up_state, down_state = MATHYS_SUBPOP_CONTRASTS[contrast]
    path = cache_file(MATHYS_SUBPOP_XLSX)
    raw = pd.read_excel(path, sheet_name=MATHYS_SUBPOP_SHEET, header=None)
    require(raw.shape[1] >= len(MATHYS_SUBPOP_COLUMNS), "mathys_subpop_sheet_too_narrow")
    header = [str(value).strip() for value in raw.iloc[1, :len(MATHYS_SUBPOP_COLUMNS)]]
    require(tuple(header) == MATHYS_SUBPOP_HEADER, "mathys_subpop_header_invalid")
    table = raw.iloc[2:, :len(MATHYS_SUBPOP_COLUMNS)].copy()
    table.columns = list(MATHYS_SUBPOP_COLUMNS)
    table = table.dropna(subset=["gene_symbol", "subpopulation"])
    table["gene_symbol"] = table["gene_symbol"].astype(str).str.strip()
    table["subpopulation"] = table["subpopulation"].astype(str).str.strip()
    table = table[table["gene_symbol"].ne("") & table["gene_symbol"].ne("nan")]
    table["lfc"] = pd.to_numeric(table["lfc"], errors="coerce")
    require(up_state in set(table["subpopulation"]), f"mathys_subpop_state_missing:{up_state}")
    require(down_state in set(table["subpopulation"]), f"mathys_subpop_state_missing:{down_state}")
    nonfinite = int((~np.isfinite(table["lfc"].to_numpy(dtype=np.float64))).sum())
    table = table[np.isfinite(table["lfc"])]
    up_rows = table[table["subpopulation"] == up_state].drop_duplicates("gene_symbol", keep="first")
    down_rows = table[table["subpopulation"] == down_state].drop_duplicates("gene_symbol", keep="first")
    require(len(up_rows) >= 3 and len(down_rows) >= 3, "mathys_subpop_signature_too_small")
    up_symbols = set(up_rows["gene_symbol"])
    down_symbols = set(down_rows["gene_symbol"])
    ambiguous = sorted(up_symbols & down_symbols)
    up_rows = up_rows[~up_rows["gene_symbol"].isin(ambiguous)]
    down_rows = down_rows[~down_rows["gene_symbol"].isin(ambiguous)]
    require(len(up_rows) >= 3 and len(down_rows) >= 3, "mathys_subpop_signature_too_small")
    weights = {str(row.gene_symbol): float(row.lfc) for row in up_rows.itertuples(index=False)}
    weights.update({str(row.gene_symbol): -float(row.lfc) for row in down_rows.itertuples(index=False)})
    return {
        "library": f"Mathys2019_Nature_MOESM8_{MATHYS_SUBPOP_SHEET}",
        "source_url": MATHYS_URL,
        "source_doi": MATHYS_DOI,
        "source_sha256": file_sha256(path),
        "contrast": contrast,
        "up_symbols": list(up_rows["gene_symbol"]),
        "down_symbols": list(down_rows["gene_symbol"]),
        "weights": weights,
        "ambiguous_symbols": ambiguous,
        "cell_context": "human_brain_microglia",
        "cell_state": f"{up_state}_vs_{down_state}_subpopulation_markers",
        "species": "human",
        "brain_region": "prefrontal_cortex",
        "selectivity": {
            "rule": (
                "MOESM8 subpopulation marker table; the up direction is the markers of the first "
                "declared subpopulation and the down direction is the markers of the second, each "
                "being over-expression relative to the other microglial subpopulations. Genes "
                "marked in both directions are dropped as ambiguous. Non-finite log fold changes "
                "are dropped before any LINCS coverage is inspected."
            ),
            "contrast": contrast,
            "up_subpopulation": up_state,
            "down_subpopulation": down_state,
            "up_markers": len(up_rows),
            "down_markers": len(down_rows),
            "dropped_nonfinite_lfc": nonfinite,
            "ambiguous_direction_genes": ambiguous,
            "not_selectivity": (
                "This is a subpopulation marker contrast within microglia. It is not a case-control "
                "disease contrast and it does not establish that either subpopulation drives disease."
            ),
        },
        "note": (
            "Mathys et al. 2019 Nature 570:332-337, MOESM8 subpopulation markers, "
            f"{up_state} versus {down_state}. Markers are defined against the other subpopulations "
            "of the same cell type, so this member tests a cell-state identity contrast rather than "
            "a pathology contrast. ROSMAP snRNA-seq microglia; not iPSC-microglia and not a LINCS "
            "cell line. This member has substantially more genes than the case-control signature, "
            "so it reaches more of the LINCS panel."
        ),
    }


def load_cellmarker_microglia() -> dict[str, Any]:
    path = cache_file(CELLMARKER_GMT)
    sets = parse_gmt(path)
    for required in (CELLMARKER_HOMEOSTATIC, CELLMARKER_ENTORHINAL, CELLMARKER_SFG):
        require(required in sets, f"cellmarker_term_missing:{required}")
    homeostatic = set(sets[CELLMARKER_HOMEOSTATIC])
    activated = set(sets[CELLMARKER_ENTORHINAL]) | set(sets[CELLMARKER_SFG])
    up = sorted(activated - homeostatic)
    down = sorted(homeostatic - activated)
    ambiguous = sorted(activated & homeostatic)
    return {
        "library": "CellMarker_2024",
        "source_url": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=CellMarker_2024",
        "source_sha256": file_sha256(path),
        "terms": {
            "homeostatic": CELLMARKER_HOMEOSTATIC,
            "entorhinal": CELLMARKER_ENTORHINAL,
            "superior_frontal_gyrus": CELLMARKER_SFG,
        },
        "up_symbols": up,
        "down_symbols": down,
        "ambiguous_symbols": ambiguous,
        "cell_context": "human_brain_microglia",
        "cell_state": "region_associated_vs_homeostatic",
        "species": "human",
        "note": "Not Mathys Mic1 DE and not a mouse DAM table. Ambiguous genes are excluded, not forced into a direction.",
    }


def load_mathys_mic_pathology_de(*, model: str = "ind") -> dict[str, Any]:
    """Mathys 2019 Nature microglia no-pathology vs pathology DE (MOESM4 Mic sheet).

    model='ind' uses DEGs.Ind.Model (cell-level FDR).
    model='mix' uses DEGs.Ind.Mix.models (individual-level mixed-model consensus).
    """
    require(model in {"ind", "mix"}, "mathys_model_invalid")
    path = cache_file(MATHYS_MIC_DE)
    table = pd.read_csv(path, sep="\t")
    required = {"gene_symbol", "log2fc", "ind_adj_pval", "deg_ind", "deg_ind_mix"}
    require(required.issubset(table.columns), "mathys_de_schema_invalid")
    flag = "deg_ind" if model == "ind" else "deg_ind_mix"
    de = table[pd.to_numeric(table[flag], errors="coerce") == 1].copy()
    de["log2fc"] = pd.to_numeric(de["log2fc"], errors="coerce")
    de = de.dropna(subset=["log2fc"])
    require(len(de) >= 10, "mathys_de_too_small")
    up = de[de["log2fc"] > 0].sort_values("log2fc", ascending=False)
    down = de[de["log2fc"] < 0].sort_values("log2fc")
    return {
        "library": "Mathys2019_Nature_MOESM4_Mic",
        "source_url": MATHYS_URL,
        "source_doi": MATHYS_DOI,
        "source_sha256": file_sha256(path),
        "model": model,
        "flag_column": flag,
        "up_symbols": list(up["gene_symbol"].astype(str)),
        "down_symbols": list(down["gene_symbol"].astype(str)),
        "weights": {str(row.gene_symbol): float(row.log2fc) for row in de.itertuples(index=False)},
        "ambiguous_symbols": [],
        "cell_context": "human_brain_microglia",
        "cell_state": "pathology_vs_no_pathology",
        "species": "human",
        "brain_region": "prefrontal_cortex",
        "note": (
            "Mathys et al. 2019 Nature 570:332-337, MOESM4 Mic sheet, "
            "no-pathology vs pathology. ROSMAP snRNA-seq microglia, not iPSC-microglia "
            "and not a LINCS cell line."
        ),
    }


def load_landmark_genes() -> pd.DataFrame:
    path = cache_file(LANDMARK)
    table = pd.read_csv(path, sep="\t", dtype=str)
    require({"pr_gene_id", "pr_gene_symbol"}.issubset(table.columns), "landmark_schema_invalid")
    table = table.rename(columns={"pr_gene_id": "entrez_id", "pr_gene_symbol": "gene_symbol"})
    return table


def load_pert_info() -> pd.DataFrame:
    path = cache_file(PERT_INFO)
    table = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    require({"pert_id", "pert_iname"}.issubset(table.columns), "pert_info_schema_invalid")
    return table


def load_lincs_consensus() -> dict[str, dict[str, list[str]]]:
    path = cache_file(LINCS_CONSENSUS)
    sets = parse_gmt(path)
    drugs: dict[str, dict[str, list[str]]] = {}
    for term, genes in sets.items():
        if term.endswith(" Up"):
            name = term[:-3]
            drugs.setdefault(name, {"up": [], "down": []})["up"] = genes
        elif term.endswith(" Down"):
            name = term[:-5]
            drugs.setdefault(name, {"up": [], "down": []})["down"] = genes
    require(drugs, "lincs_consensus_empty")
    return {
        "library": "LINCS_L1000_Chem_Pert_Consensus_Sigs",
        "source_url": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=LINCS_L1000_Chem_Pert_Consensus_Sigs",
        "source_sha256": file_sha256(path),
        "cell_context": "lincs_l1000_chem_pert_consensus",
        "drugs": drugs,
    }


def load_lincs_instances_for_names(names: list[str]) -> list[dict[str, Any]]:
    """Optional cell-line-resolved instances for named compounds only."""
    wanted = {name.lower() for name in names}
    rows = []
    for direction, filename in (("up", LINCS_UP), ("down", LINCS_DOWN)):
        path = cache_file(filename)
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            term = line.split("\t", 1)[0]
            lowered = term.lower()
            if not any(name in lowered for name in wanted):
                continue
            genes = [item.split(",")[0].strip() for item in line.split("\t")[2:] if item.strip()]
            rows.append({"term": term, "direction": direction, "genes": genes})
    return rows


def load_pubchem_gold_standard() -> dict[str, Any]:
    path = cache_file("pubchem_gold_standard.json")
    return json.loads(path.read_text(encoding="utf-8"))
