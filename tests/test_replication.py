import copy
import json
from pathlib import Path
import pandas as pd
import pytest
from microglia_repurposing.replication import (bh_adjust, compare_rankings,
    seaad_signature, seal, validate_signature, verify_document)
from microglia_repurposing.contracts import file_sha256

ROOT = Path(__file__).resolve().parents[1]

def test_bh_accounts_for_the_whole_family_and_ties():
    assert bh_adjust([.01,.04,.03,.5]).tolist() == pytest.approx([.04,.053333333,.053333333,.5])
    assert bh_adjust([.01,.01,1]).tolist() == pytest.approx([.015,.015,1])
    with pytest.raises(ValueError):bh_adjust([float('nan')])
    with pytest.raises(ValueError):bh_adjust([-0.1])

def test_comparison_preserves_opposition_neutral_and_missing_in_full_union():
    def row(key, score, status='scored'):
        return dict(pert_id=key,pert_iname=key,connectivity_mean=score,status=status,fdr=.8)
    a={'ranking':[row('opposed',.4),row('positive',.2),row('negative',-.3),row('neutral',0),row('missing',.5)]}
    b={'ranking':[row('opposed',-.4),row('positive',.1),row('negative',-.1),row('neutral',.2),row('missing',None,'insufficient_coverage'),row('new',.9)]}
    comparison=compare_rankings(a,b)
    assert comparison['counts']==dict(opposite_direction=1,both_positive=1,both_negative=1,neutral=1,missing_coverage=2)
    missing=next(r for r in comparison['records'] if r['pert_id']=='missing')
    assert missing['external_score'] is None and missing['external_fdr'] is None
    new=next(r for r in comparison['records'] if r['pert_id']=='new')
    assert new['baseline_status']=='not_present' and new['baseline_score'] is None

def test_published_csv_uses_symbol_index_not_numeric_gene_id_and_keeps_sparse_result():
    protocol=verify_document(json.loads((ROOT/'data/replication/protocol.json').read_text()))
    signature=seaad_signature(ROOT/'data/replication/sources/Micro-PVM_1.csv','Micro-PVM_1',protocol)
    assert signature['counts']==dict(up=3,down=1,unmapped=0)
    assert all(not g['gene_symbol'].isdigit() for g in signature['genes'])
    with pytest.raises(ValueError,match='too_small'):validate_signature(signature)

def test_source_tampering_and_bad_probability_fail_closed(tmp_path):
    protocol=json.loads((ROOT/'data/replication/protocol.json').read_text())
    source=ROOT/'data/replication/sources/Micro-PVM_1.csv'
    p=tmp_path/'state.csv';p.write_bytes(source.read_bytes())
    receipt=json.loads(source.with_suffix('.provenance.json').read_text())
    p.with_suffix('.provenance.json').write_text(json.dumps(receipt))
    with p.open('a') as f:f.write('corrupt')
    with pytest.raises(ValueError,match='source_hash'):seaad_signature(p,'Micro-PVM_1',protocol)
    table=pd.read_csv(source);table.loc[0,'p_Continuous_Pseudo-progression_Score']=-1
    table.to_csv(p,index=False);receipt['sha256']=file_sha256(p)
    p.with_suffix('.provenance.json').write_text(json.dumps(receipt))
    with pytest.raises(ValueError,match='invalid_p_value'):seaad_signature(p,'Micro-PVM_1',protocol)

def test_complete_report_keeps_every_frozen_state_and_no_replication_claim():
    folder=ROOT/'data/replication'
    protocol=verify_document(json.loads((folder/'protocol.json').read_text()))
    report=verify_document(json.loads((folder/'report.json').read_text()))
    assert report['protocol_hash']==protocol['content_hash']
    assert [s['state'] for s in report['states']]==protocol['states']
    assert sum(s['status']=='evaluated' for s in report['states'])==4
    assert report['independent_replication_established'] is False
    for entry in report['states']:
        if entry['status']=='evaluated':
            run=verify_document(json.loads((folder/'runs'/(entry['run_hash']+'.json')).read_text()))
            assert run['configuration']==dict(model='external',n_perm=2000,min_hits=3,seed=0)
            assert sum(entry['comparison']['counts'].values())==len(entry['comparison']['records'])
        else:
            assert entry['run_hash'] is None and entry['comparison'] is None
    benchmark=verify_document(json.loads((folder/'benchmarks.json').read_text()))
    assert benchmark['rescue_auroc'] is None and benchmark['eligible_rescue_controls']==0
    assert all(e['benchmark_label'] is None and not e['eligible_for_rescue_benchmark'] for e in benchmark['entries'])
