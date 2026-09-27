import importlib.util
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('dashboard_server',ROOT/'dashboard_server.py')
server=importlib.util.module_from_spec(spec);spec.loader.exec_module(server)

@pytest.fixture
def workspace(tmp_path):
    return server.Workspace(ROOT/'data/review',tmp_path/'reviews.sqlite3')

def test_reviews_survive_reopen_without_changing_scores(workspace):
    run=workspace.runs['context_mix'];key='BRD-K81418486';before=dict(run['lookup'][key])
    body=dict(run='context_mix',run_hash=run['content_hash'],candidate_id=key,status='Needs evidence',note='Synthetic test',reference='',outcome='')
    workspace.save(body)
    reopened=server.ReviewStore(workspace.store.path)
    assert reopened.get(run['content_hash'],key)['note']=='Synthetic test'
    assert run['lookup'][key]==before
    assert reopened.all()[0]['biological_validation_changed'] is False

def test_wrong_hash_unknown_candidate_and_unreferenced_outcome_rejected(workspace):
    body=dict(run='context_mix',run_hash=workspace.runs['context_mix']['content_hash'],candidate_id='BRD-K81418486',status='Observed result',note='Synthetic',outcome='Negative',reference='')
    with pytest.raises(ValueError,match='require'):workspace.save(body)
    body['reference']='synthetic-assay';body['run_hash']='stale'
    with pytest.raises(ValueError,match='changed'):workspace.save(body)
    body['run_hash']=workspace.runs['context_mix']['content_hash'];body['candidate_id']='not-present'
    with pytest.raises(ValueError,match='changed'):workspace.save(body)

def test_exact_id_model_comparison_keeps_missing_coverage_explicit(workspace):
    result=workspace.candidate('context_mix','BRD-K07691486')
    assert next(r for r in result['comparisons'] if r['run']=='context_mix')['rank'] is None
    assert next(r for r in result['comparisons'] if r['run']=='context_ind')['rank']==2
    assert workspace.manifest()['runs'][0]['eligible']==0

def test_corrupt_analysis_is_rejected(tmp_path):
    (tmp_path/'summary.json').write_text(json.dumps({'runs':{'context_mix':{}}}))
    (tmp_path/'context_mix.json').write_text(json.dumps({'content_hash':'wrong'}))
    with pytest.raises(ValueError,match='hash mismatch'):server.Workspace(tmp_path,tmp_path/'review.db')

@pytest.mark.parametrize('origin,token',[('https://unrelated.example','valid'),('http://127.0.0.1:8879','invalid')])
def test_mutating_api_rejects_cross_origin_and_invalid_session(workspace,origin,token):
    from types import SimpleNamespace
    instance=object.__new__(server.handler(workspace));instance.server=SimpleNamespace(server_port=8879)
    instance.headers={'Host':'127.0.0.1:8879','Origin':origin,'X-Review-Token':workspace.token if token=='valid' else 'wrong'}
    responses=[];instance.send_json=lambda status,body:responses.append(status)
    instance.do_POST()
    assert responses==[403]
    assert workspace.store.all()==[]

def test_archive_retains_original_scores_and_allows_original_review(workspace):
    import copy
    name='context_mix';old=workspace.runs[name];digest=old['content_hash'];key='BRD-K81418486'
    expected=dict(old['lookup'][key])
    newer=copy.deepcopy(old);newer['content_hash']='f'*64;newer['lookup'][key]['connectivity_mean']=-.99
    workspace.runs[name]=newer
    detail=workspace.candidate(name,key,digest)
    assert detail['historical'] is True
    assert detail['row']==expected
    body=dict(run=name,run_hash=digest,candidate_id=key,status='Needs evidence',note='Original run review',reference='',outcome='')
    workspace.save(body)
    assert workspace.store.get(digest,key)['note']=='Original run review'
    assert workspace.store.get(newer['content_hash'],key) is None
    assert workspace.candidate(name,key)['row']['connectivity_mean']==-.99

def test_refresh_loads_validated_snapshot_and_preserves_token_and_reviews(workspace):
    name='context_mix';run=workspace.runs[name];key='BRD-K81418486'
    workspace.save(dict(run=name,run_hash=run['content_hash'],candidate_id=key,status='Needs evidence',note='Refresh canary',reference='',outcome=''))
    replacement=workspace.refreshed()
    assert replacement.token==workspace.token
    assert replacement.candidate(name,key)['review']['note']=='Refresh canary'
    assert replacement.loaded_fingerprint==workspace.loaded_fingerprint
    assert replacement.replication_summary()['report']['independent_replication_established'] is False
    assert len(replacement.replication_summary()['report']['states'])==6
    assert replacement.rescue_summary()['report']['rescue_auroc'] is None
    assert replacement.rescue_summary()['report']['eligible_rescue_controls']==0

def test_failed_refresh_keeps_previous_snapshot_and_reviews(workspace,monkeypatch):
    from types import SimpleNamespace
    handler_type=server.handler(workspace)
    instance=object.__new__(handler_type);instance.server=SimpleNamespace(server_port=8879)
    instance.headers={'Host':'127.0.0.1:8879','X-Review-Token':workspace.token};instance.path='/api/refresh'
    old=workspace.loaded_fingerprint
    def fail():raise ValueError('Corrupt staged source')
    monkeypatch.setattr(workspace,'refreshed',fail)
    responses=[];instance.send_json=lambda status,body:responses.append((status,body))
    instance.do_POST()
    assert responses[0][0]==400
    assert workspace.loaded_fingerprint==old
    instance.path='/api/manifest';instance.do_GET()
    assert responses[1][0]==200 and responses[1][1]['snapshot']==old
