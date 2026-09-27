#!/usr/bin/env python3
"""Evaluate all frozen SEA-AD states and atomically publish a complete report."""
import json
import sys
from pathlib import Path
from datetime import datetime, timezone
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'src'))
from microglia_repurposing.replication import verify_document, seal, seaad_signature, compare_rankings
from microglia_repurposing.contracts import file_sha256, require
from microglia_repurposing.context_engine import run_context_engine


def write_immutable(folder, document):
    verify_document(document)
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / (document['content_hash'] + '.json')
    if path.exists():
        verify_document(json.loads(path.read_text()))
        return path
    with path.open('x') as f:
        json.dump(document, f, indent=2, allow_nan=False)
    return path


def main():
    folder = ROOT / 'data/replication'
    protocol = verify_document(json.loads((folder / 'protocol.json').read_text()))
    baseline = verify_document(json.loads((ROOT / 'data/review/context_mix.json').read_text()))
    require(baseline['content_hash'] == protocol['baseline_run_hash'], 'baseline_changed_after_freeze')
    for path, expected in {**protocol['source_code_hashes'], **protocol['baseline_inputs']}.items():
        require(file_sha256(ROOT / path) == expected, 'frozen_input_changed:' + path)
    # Pin the actual perturbation cache to the already frozen baseline run as well.
    for name, receipt in baseline['provenance']['sources'].items():
        require(file_sha256(ROOT / 'data/cache' / name) == receipt['sha256'], 'panel_input_changed:' + name)
    implementation = seal({'schema':'replication-implementation/v1','protocol_hash':protocol['content_hash'],
        'files':{str(p.relative_to(ROOT)):file_sha256(p) for p in [Path(__file__), ROOT/'src/microglia_repurposing/replication.py']}})
    write_immutable(folder/'implementations', implementation)
    baseline_path = write_immutable(folder / 'runs', baseline)
    entries = []
    config = protocol['scoring']
    for state in protocol['states']:
        print('Evaluating', state, flush=True)
        signature = seaad_signature(folder/'sources'/(state+'.csv'), state, protocol)
        write_immutable(folder/'signatures', signature)
        entry = {'state':state, 'signature_hash':signature['content_hash'], 'counts':signature['counts'],
                 'import_quality':signature['import_quality'], 'source':signature['source']}
        if min(signature['counts']['up'], signature['counts']['down']) < protocol['selection']['min_genes_per_direction']:
            entry.update(status='insufficient_signature', run_hash=None, comparison=None,
                         reason='Frozen selection yielded fewer than three genes in at least one direction; settings unchanged.')
        else:
            result = run_context_engine(disease_bundle=signature, model='external', output=None,
                                        n_perm=config['n_perm'], min_hits=config['min_hits'], seed=config['seed'])
            run = {k:v for k,v in result.items() if not k.startswith('_')}
            write_immutable(folder/'runs', run)
            entry.update(status='evaluated', run_hash=run['content_hash'],
                         comparison=compare_rankings(baseline, run), statistics=run['statistics'],
                         eligible=sum(r['status']=='scored' and (r.get('connectivity_mean') or 0)>0 and
                                      r.get('fdr') is not None and r['fdr']<=config['fdr_threshold'] for r in run['ranking']))
        entries.append(entry)
    report = seal({'schema':'titans-replication-report/v1', 'created_at':datetime.now(timezone.utc).isoformat(),
        'protocol_hash':protocol['content_hash'], 'implementation_hash':implementation['content_hash'],
        'baseline_hash':baseline['content_hash'], 'states':entries,
        'independent_replication_established':False,
        'limitations':['Donor-level crosswalk is not available locally; distinct recruitment sources are documented, but donor non-overlap is not verified.',
                       'SEA-AD MTG state-specific continuous progression differs from Mathys prefrontal binary pathology contrast.',
                       'All six SEA-AD states share a cohort and all runs share the LINCS panel; neither is six independent validations.',
                       'Per-state FDR families differ with coverage. No joint any-state FDR or biological rescue claim.',
                       'Original published aggregate coefficients reused; no new donor-level model fitted.'],
        'claim_boundary':'External cohort transfer evaluation only; no experimental rescue or clinical validation.'})
    write_immutable(folder/'reports', report)
    # Publish last: the dashboard cannot observe half-written state results.
    pointer = folder/'report.tmp'
    pointer.write_text(json.dumps(report, indent=2, allow_nan=False))
    pointer.replace(folder/'report.json')
    print(json.dumps({'report':report['content_hash'],'states':[{'state':e['state'],'status':e['status'],
            'counts':e['counts'],'comparison':e['comparison']['counts'] if e['comparison'] else None} for e in entries]},indent=2))

if __name__ == '__main__':
    main()
