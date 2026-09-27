"""Pin gene symbols to official HGNC identity. Do not invent identifiers."""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping

import pandas as pd

from microglia_repurposing.contracts import require

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_HGNC = PACKAGE_ROOT / "data" / "cache" / "hgnc_complete_set.txt"


@lru_cache(maxsize=1)
def load_hgnc_table(path: str | Path | None = None) -> pd.DataFrame:
    target = Path(path) if path is not None else DEFAULT_HGNC
    require(target.is_file(), "hgnc_table_missing")
    from microglia_repurposing.provenance import verify_source
    verify_source(target)
    table = pd.read_csv(target, sep="\t", dtype=str, low_memory=False, keep_default_na=False)
    require("symbol" in table.columns and "hgnc_id" in table.columns, "hgnc_table_schema_invalid")
    if "status" in table.columns:
        table = table[table["status"] == "Approved"].copy()
    return table.reset_index(drop=True)


def _split_aliases(value: str) -> list[str]:
    if not value:
        return []
    return [part.strip() for part in str(value).replace(",", "|").split("|") if part.strip()]


@lru_cache(maxsize=1)
def _symbol_index(path: str | None = None) -> dict[str, dict[str, str]]:
    table = load_hgnc_table(path)
    index: dict[str, dict[str, str]] = {}
    ambiguous = set()
    official = set(table["symbol"].str.upper())
    for row in table.itertuples(index=False):
        record = {
            "gene_symbol": str(getattr(row, "symbol")),
            "hgnc_id": str(getattr(row, "hgnc_id")),
            "name": str(getattr(row, "name", "")),
            "status": str(getattr(row, "status", "Approved")),
            "entrez_id": str(getattr(row, "entrez_id", "")),
            "ensembl_gene_id": str(getattr(row, "ensembl_gene_id", "")),
            "uniprot_ids": str(getattr(row, "uniprot_ids", "")),
        }
        keys = {record["gene_symbol"].upper()}
        for field in ("alias_symbol", "prev_symbol"):
            if hasattr(row, field):
                keys.update(alias.upper() for alias in _split_aliases(getattr(row, field)))
        for key in keys:
            # Official symbol wins if an alias collides.
            current = index.get(key)
            if key in ambiguous and key not in official:
                continue
            if current and current["hgnc_id"] != record["hgnc_id"] and key not in official:
                ambiguous.add(key)
                index.pop(key, None)
            elif current is None or current["gene_symbol"].upper() != key:
                index[key] = record
    return index


def pin_symbol(symbol: str, *, hgnc_path: str | Path | None = None) -> dict[str, str] | None:
    text = str(symbol or "").strip()
    if not text:
        return None
    index = _symbol_index(None if hgnc_path is None else str(hgnc_path))
    hit = index.get(text.upper())
    if hit is None:
        return None
    return dict(hit)


def pin_symbols(symbols: Iterable[str], *, hgnc_path: str | Path | None = None) -> tuple[list[dict[str, str]], list[str]]:
    pinned = []
    missing = []
    seen = set()
    for symbol in symbols:
        row = pin_symbol(symbol, hgnc_path=hgnc_path)
        if row is None:
            missing.append(str(symbol))
            continue
        if row["hgnc_id"] in seen:
            continue
        seen.add(row["hgnc_id"])
        pinned.append(row)
    pinned.sort(key=lambda item: item["gene_symbol"])
    return pinned, sorted(set(missing))


def spine_overlay(symbol: str) -> dict[str, Any] | None:
    """Attach Titans Med spine identity when the gene is already catalogued."""
    try:
        from pathlib import Path as _Path
        import sys
        vault = _Path(__file__).resolve().parents[3]
        if str(vault) not in sys.path:
            sys.path.insert(0, str(vault))
        from forge_genomics_protein_spine import resolve_identity_query
        hits = [hit for hit in resolve_identity_query(symbol, limit=8) if hit.get("kind") == "gene" and hit.get("score", 0) >= 80]
    except Exception:
        return None
    return hits[0] if hits else None
