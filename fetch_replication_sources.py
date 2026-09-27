#!/usr/bin/env python3
"""Fetch only the six declared public aggregate tables; never individual-cell data."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
BASE = 'https://sea-ad-single-cell-profiling.s3.amazonaws.com/'

def main():
    folder = ROOT / 'data/replication'
    protocol = json.loads((folder / 'protocol.json').read_text())
    raw = folder / 'sources'
    raw.mkdir(exist_ok=True)
    for state in protocol['states']:
        key = ('MTG/RNAseq/Supplementary Information/Nebula Results/All Donors/'
               'Continuous_Pseudo-progression_Score/Micro-PVM_' + state +
               '_across_Continuous_Pseudo-progression_Score_DE.csv')
        url = BASE + quote(key)
        path = raw / (state + '.csv')
        receipt = path.with_suffix('.provenance.json')
        if path.exists():
            expected = json.loads(receipt.read_text())['sha256']
            if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise ValueError('Cached source checksum mismatch: ' + state)
            print('Verified cached source:', state, flush=True)
            continue
        with urlopen(url, timeout=60) as response:
            data = response.read(12_000_001)
            if len(data) > 12_000_000:
                raise ValueError('Aggregate source exceeds declared 12 MB limit')
            metadata = dict(url=url, retrieved_at=datetime.now(timezone.utc).isoformat(),
                            sha256=hashlib.sha256(data).hexdigest(), bytes=len(data),
                            etag=response.headers.get('ETag'),
                            last_modified=response.headers.get('Last-Modified'),
                            upstream_checksum_verified=False,
                            data_kind='public aggregate differential expression',
                            protocol_hash=protocol['content_hash'])
        temp = path.with_suffix('.tmp')
        temp.write_bytes(data)
        receipt.write_text(json.dumps(metadata, indent=2))
        temp.replace(path)
        print('Saved aggregate source:', state, len(data), 'bytes', flush=True)

if __name__ == '__main__':
    main()
