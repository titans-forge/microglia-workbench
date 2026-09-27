"""Source and runtime receipts. Hash caching invalidates on file metadata change."""
from __future__ import annotations

import json
import platform
from functools import lru_cache
from pathlib import Path

from microglia_repurposing.contracts import file_sha256, require

ROOT = Path(__file__).resolve().parents[2]


@lru_cache(maxsize=128)
def _fingerprint(path: str, size: int, mtime_ns: int, ctime_ns: int) -> str:
    return file_sha256(path)


def fingerprint(path: Path) -> str:
    stat = path.stat()
    return _fingerprint(str(path.resolve()), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


def verify_source(path: Path) -> str:
    actual = fingerprint(path)
    sidecar = path.with_name(path.name + ".provenance.json")
    if sidecar.is_file():
        expected = json.loads(sidecar.read_text()).get("sha256")
        if expected:
            require(actual == expected, f"source_hash_mismatch:{path.name}")
    return actual


def run_provenance() -> dict:
    import numpy, pandas
    sources = {}
    for path in sorted((ROOT / "data/cache").iterdir()):
        if not path.is_file() or path.suffix == ".json":
            continue
        digest = verify_source(path)
        sidecar = path.with_name(path.name + ".provenance.json")
        metadata = json.loads(sidecar.read_text()) if sidecar.is_file() else {}
        sources[path.name] = {"sha256": digest, "upstream_hash_verified": bool(metadata.get("sha256")),
                              "source": metadata}
    return {"python": platform.python_version(), "numpy": numpy.__version__, "pandas": pandas.__version__,
            "code": {str(p.relative_to(ROOT)): fingerprint(p) for p in sorted((ROOT / "src").rglob("*.py"))},
            "sources": sources, "control_panel_sha256": fingerprint(ROOT / "data/control_panel_ad.csv")}
