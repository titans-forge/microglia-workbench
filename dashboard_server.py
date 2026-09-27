#!/usr/bin/env python3
"""Loopback-only evidence review dashboard. Reviews never alter scientific scores."""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import secrets
import sqlite3
import sys
import threading
from urllib.parse import urlparse, parse_qs

ROOT = Path(__file__).resolve().parent
if str(ROOT / 'src') not in sys.path:
    sys.path.insert(0, str(ROOT / 'src'))
APP = 'titans-repurposing-dashboard'
VERSION = '0.6.0'
LABELS = {'context_mix':'Mixed model · primary', 'context_ind':'Original model · sensitivity',
          'context_mix_overlap5':'Higher overlap · sensitivity', 'consensus_controls':'Consensus controls',
          'consensus_all':'Full consensus'}
STATUSES = {'Needs evidence','Plan experiment','Deprioritize','Observed result'}


def canonical_hash(data):
    raw=(json.dumps(data,sort_keys=True,separators=(',', ':'),ensure_ascii=True,allow_nan=False)+'\n').encode()
    return hashlib.sha256(raw).hexdigest()


def verified(document):
    if canonical_hash({k:v for k,v in document.items() if k!='content_hash'}) != document.get('content_hash'):
        raise ValueError('Document hash mismatch')
    return document


def build_fingerprint():
    paths=[ROOT/'dashboard_server.py', ROOT/'src/microglia_repurposing/rescue.py', *sorted((ROOT/'dashboard').glob('*'))]
    return canonical_hash({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()})


class ReviewStore:
    def __init__(self, path):
        self.path=path
        path.parent.mkdir(parents=True,exist_ok=True)
        with self.connect() as db:
            db.execute('CREATE TABLE IF NOT EXISTS reviews (run_hash TEXT, candidate_id TEXT, body TEXT NOT NULL, PRIMARY KEY(run_hash,candidate_id))')
            db.execute('CREATE TABLE IF NOT EXISTS history (id INTEGER PRIMARY KEY, body TEXT NOT NULL)')
            db.execute('CREATE TABLE IF NOT EXISTS analyses (hash TEXT PRIMARY KEY, name TEXT NOT NULL, body TEXT NOT NULL, contexts TEXT NOT NULL)')

    def archive(self, name, run, contexts):
        raw={k:v for k,v in run.items() if k!='lookup'}
        verified(raw)
        with self.connect() as db:
            db.execute('INSERT OR IGNORE INTO analyses VALUES(?,?,?,?)',
                       (raw['content_hash'], name, json.dumps(raw,allow_nan=False), json.dumps(contexts,allow_nan=False)))

    def analysis(self, digest):
        with self.connect() as db:
            row=db.execute('SELECT name,body,contexts FROM analyses WHERE hash=?',(digest,)).fetchone()
        if row is None: raise ValueError('Original analysis snapshot is unavailable')
        return row[0], verified(json.loads(row[1])), json.loads(row[2])

    def connect(self):
        return sqlite3.connect(self.path,timeout=10)

    def get(self, run_hash, candidate_id):
        with self.connect() as db:
            row=db.execute('SELECT body FROM reviews WHERE run_hash=? AND candidate_id=?',(run_hash,candidate_id)).fetchone()
        return json.loads(row[0]) if row else None

    def all(self):
        with self.connect() as db:
            return [json.loads(r[0]) for r in db.execute('SELECT body FROM reviews ORDER BY rowid DESC')]

    def save(self, body):
        text=json.dumps(body,allow_nan=False)
        with self.connect() as db:
            db.execute('INSERT INTO reviews VALUES(?,?,?) ON CONFLICT(run_hash,candidate_id) DO UPDATE SET body=excluded.body',
                       (body['run_hash'],body['candidate_id'],text))
            db.execute('INSERT INTO history(body) VALUES(?)',(text,))
        return body


class Workspace:
    def __init__(self, folder, database):
        self.folder=folder.resolve();self.store=ReviewStore(database);self.runs={};self.token=secrets.token_urlsafe(32)
        self.build=build_fingerprint()
        self.labels=dict(LABELS)
        self.replication=None;self.datasets=None;self.benchmarks=None;self.rescue=None
        before=self.disk_fingerprint()
        summary=json.loads((folder/'summary.json').read_text())
        for name in summary['runs']:
            if name not in LABELS: continue
            run=json.loads((folder/(name+'.json')).read_text())
            digest=canonical_hash({k:v for k,v in run.items() if k!='content_hash'})
            if digest!=run['content_hash']: raise ValueError(f'Analysis hash mismatch: {name}')
            if not run['schema'].endswith('/v2') or not run['research_planning_only']: raise ValueError('Unsupported analysis')
            run['lookup']={}
            for i,row in enumerate(run['ranking']):
                row=dict(row);row['rank']=i+1;key=row.get('pert_id',row.get('cmap_name'))
                run['lookup'][key]=row
            self.runs[name]=run
        if 'context_mix' not in self.runs: raise ValueError('Primary analysis is missing')
        self.comparisons=summary.get('comparisons',{})
        self.contexts={}
        path=folder/'context_mix_contexts.csv'
        if path.exists():
            with path.open() as stream:
                for row in csv.DictReader(stream): self.contexts.setdefault(row['pert_id'],[]).append(row)
        replication_folder=self.folder.parent/'replication'
        if (replication_folder/'report.json').exists():
            self.replication=verified(json.loads((replication_folder/'report.json').read_text()))
            protocol=verified(json.loads((replication_folder/'protocol.json').read_text()))
            if self.replication['protocol_hash']!=protocol['content_hash']:raise ValueError('Replication protocol mismatch')
            baseline=verified(json.loads((replication_folder/'runs'/(self.replication['baseline_hash']+'.json')).read_text()))
            if baseline['content_hash']!=self.replication['baseline_hash']:raise ValueError('Replication baseline mismatch')
            if baseline['content_hash']!=self.runs['context_mix']['content_hash']:
                self.store.archive('context_mix',baseline,{})
            self.datasets=verified(json.loads((replication_folder/'datasets.json').read_text()))
            self.benchmarks=verified(json.loads((replication_folder/'benchmarks.json').read_text()))
            for entry in self.replication['states']:
                if entry['status']!='evaluated':continue
                run=verified(json.loads((replication_folder/'runs'/(entry['run_hash']+'.json')).read_text()))
                if run['content_hash']!=entry['run_hash']:raise ValueError('Replication run mismatch')
                name='seaad_'+entry['state'];self.labels[name]='SEA-AD · '+entry['state']
                run['lookup']={r['pert_id']:{**r,'rank':i+1} for i,r in enumerate(run['ranking'])}
                self.runs[name]=run
        rescue_folder=self.folder.parent/'rescue'/'v2'
        if (rescue_folder/'protocol.json').exists():
            from microglia_repurposing.rescue import workspace_status
            self.rescue=workspace_status()
        self.loaded_fingerprint=self.disk_fingerprint()
        if self.loaded_fingerprint!=before:raise ValueError('Analysis files changed while loading; retry refresh')
        for name,run in self.runs.items():self.store.archive(name,run,self.contexts if name=='context_mix' else {})

    def disk_fingerprint(self):
        paths=sorted(self.folder.glob('*.json'))+[self.folder/'context_mix_contexts.csv']
        paths += [self.folder.parent/'replication'/name for name in ['report.json','protocol.json','datasets.json','benchmarks.json']]
        paths += [self.folder.parent/'rescue'/'v2'/name for name in ['protocol.json','panel.json','observations.json','report.json']]
        return canonical_hash({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()})

    def refreshed(self):
        result=Workspace(self.folder,self.store.path)
        result.token=self.token
        # A data refresh cannot claim the loaded Python process changed its code.
        result.build=self.build
        return result

    def replication_summary(self):
        if self.replication is None:return {'status':'not_run'}
        report={k:v for k,v in self.replication.items() if k!='states'}
        report['states']=[{**{k:v for k,v in e.items() if k!='comparison'},
                           'comparison':{'counts':e['comparison']['counts']} if e.get('comparison') else None} for e in self.replication['states']]
        return {'status':'available','report':report,'datasets':self.datasets,'benchmarks':self.benchmarks,
                'baseline_is_current':self.replication['baseline_hash']==self.runs['context_mix']['content_hash']}

    def rescue_summary(self):
        if self.rescue is None:return {'status':'not_run'}
        return self.rescue

    def get_run(self, name, digest=None):
        if digest and (name not in self.runs or self.runs[name]['content_hash']!=digest):
            original,run,contexts=self.store.analysis(digest)
            if original!=name:raise ValueError('Snapshot does not belong to this analysis')
            run['lookup']={r.get('pert_id',r.get('cmap_name')):{**r,'rank':i+1} for i,r in enumerate(run['ranking'])}
            return run,contexts
        return self.runs[name], self.contexts if name=='context_mix' else {}

    def manifest(self):
        runs=[]
        for key,run in self.runs.items():
            scored=[r for r in run['ranking'] if r['status']=='scored']
            positives=[r for r in scored if (r.get('connectivity_mean',r.get('cosine_reversal')) or 0)>0 and r.get('fdr') is not None and r['fdr']<=.25]
            runs.append({'id':key,'label':self.labels[key],'total':len(run['ranking']),'scored':len(scored),'eligible':len(positives),
                         'hash':run['content_hash'],'configuration':run['configuration'],
                         'gene_count':sum(run['disease']['counts'].get(k,0) for k in ['up','down'])})
        return {'app':APP,'version':VERSION,'runs':runs,'comparisons':self.comparisons,'review_token':self.token,
                'snapshot':self.loaded_fingerprint,'refresh_available':self.disk_fingerprint()!=self.loaded_fingerprint}

    def candidate(self, name, key, digest=None):
        run,contexts=self.get_run(name,digest);row=run['lookup'][key];comparisons=[]
        historical=name not in self.runs or run['content_hash']!=self.runs[name]['content_hash']
        if name.startswith('context_') and not historical:
            for other,r in self.runs.items():
                if not other.startswith('context_'):continue
                found=r['lookup'].get(key)
                comparisons.append({'analysis':LABELS[other],'run':other,'status':found['status'] if found else 'missing',
                                    'rank':found['rank'] if found and found['status']=='scored' else None,
                                    'score':found.get('connectivity_mean') if found else None,'fdr':found.get('fdr') if found else None})
        gaps=['Independent rescue validation not established','Viability validation not established']
        if row.get('status')!='scored':gaps.insert(0,'Insufficient gene overlap')
        if row.get('n_cells_scored',0)==1:gaps.append('Evidence comes from one cell line')
        if row.get('tox_flag'):gaps.append('Mechanism toxicity flag requires follow-up')
        if row.get('moa')=='unknown':gaps.append('Mechanism annotation unavailable')
        if any(c['status']!='scored' for c in comparisons):gaps.append('Coverage changes across model or overlap choices')
        external=[]
        if self.replication and name=='context_mix' and run['content_hash']==self.replication['baseline_hash']:
            for state in self.replication['states']:
                match=next((r for r in (state.get('comparison') or {}).get('records',[]) if r['pert_id']==key),None)
                external.append({'state':state['state'],'status':state['status'],'comparison':match})
        return {'row':row,'run_hash':run['content_hash'],'comparisons':comparisons,'gaps':gaps,
                'contexts':contexts.get(key,[]), 'historical':historical,'external':external,
                'review':self.store.get(run['content_hash'],key),'statistics':run['statistics']}

    def save(self, body):
        try:run,_=self.get_run(body['run'],body.get('run_hash'))
        except ValueError:raise ValueError('Candidate or analysis changed; original snapshot unavailable')
        key=body['candidate_id']
        if key not in run['lookup'] or body.get('run_hash')!=run['content_hash']:raise ValueError('Candidate or analysis changed; reload before saving')
        if body.get('status') not in STATUSES:raise ValueError('Invalid review status')
        for field in ['note','reference','outcome']:
            if not isinstance(body.get(field,''),str) or len(body.get(field,''))>8000:raise ValueError('Invalid review field')
        if body['status']=='Observed result' and (not body.get('reference','').strip() or body.get('outcome') not in ['Positive','Negative','Inconclusive']):
            raise ValueError('Observed results require an assay receipt or reference and an outcome')
        record={k:body.get(k,'') for k in ['run','run_hash','candidate_id','status','note','reference','outcome']}
        record.update(schema='repurposing-review/v2',saved_at=datetime.now(timezone.utc).isoformat(),biological_validation_changed=False)
        return self.store.save(record)


def handler(workspace):
    refresh_lock=threading.Lock()
    class Handler(BaseHTTPRequestHandler):
        def send_json(self, status, value):
            data=json.dumps(value,allow_nan=False).encode();self.send_response(status)
            self.send_header('Content-Type','application/json');self.send_header('Content-Length',str(len(data)))
            self.send_header('Cache-Control','no-store');self.end_headers();self.wfile.write(data)
        def valid_host(self):
            return self.headers.get('Host')==f'127.0.0.1:{self.server.server_port}'
        def do_GET(self):
            if not self.valid_host():return self.send_json(403,{'error':'Loopback host required'})
            parsed=urlparse(self.path);q=parse_qs(parsed.query)
            try:
                if parsed.path=='/health':return self.send_json(200,{'app':APP,'version':VERSION,'workspace':str(workspace.folder),
                    'build':workspace.build,'snapshot':workspace.loaded_fingerprint})
                if parsed.path=='/api/replication':return self.send_json(200,workspace.replication_summary())
                if parsed.path=='/api/replication/comparisons':
                    if workspace.replication is None:raise ValueError('Replication has not run')
                    state=next((s for s in workspace.replication['states'] if s['state']==q.get('state',[''])[0]),None)
                    if state is None:raise ValueError('Unknown external state')
                    rows=(state.get('comparison') or {}).get('records',[])
                    term=q.get('q',[''])[0].casefold();category=q.get('category',['all'])[0]
                    rows=[r for r in rows if term in (r['name']+' '+r['pert_id']).casefold() and (category=='all' or r['comparison']==category)]
                    offset=max(0,int(q.get('offset',['0'])[0]))
                    return self.send_json(200,{'total':len(rows),'rows':rows[offset:offset+30]})
                if parsed.path=='/api/rescue':return self.send_json(200,workspace.rescue_summary())
                if parsed.path=='/api/manifest':return self.send_json(200,workspace.manifest())
                if parsed.path=='/api/reviews':return self.send_json(200,workspace.store.all())
                if parsed.path=='/api/candidate':return self.send_json(200,workspace.candidate(q['run'][0],q['id'][0],q.get('hash',[None])[0]))
                if parsed.path=='/api/candidates':
                    run=workspace.runs[q.get('run',['context_mix'])[0]];rows=list(run['lookup'].values())
                    query=q.get('q',[''])[0].casefold();coverage=q.get('coverage',['scored'])[0]
                    rows=[r for r in rows if (coverage=='all' or r['status']=='scored') and query in ' '.join(str(r.get(k,'')) for k in ['pert_iname','cmap_name','pert_id']).casefold()]
                    if q.get('filter',['all'])[0]=='flagged':rows=[r for r in rows if r.get('tox_flag')]
                    if q.get('filter',['all'])[0]=='multicell':rows=[r for r in rows if (r.get('n_cells_scored') or 0)>=2]
                    offset=max(0,int(q.get('offset',['0'])[0]));return self.send_json(200,{'total':len(rows),'rows':rows[offset:offset+30]})
                if parsed.path.startswith('/analysis/'):
                    name=parsed.path.removeprefix('/analysis/').removesuffix('.json')
                    run,_=workspace.get_run(name,q.get('hash',[None])[0])
                    return self.send_json(200,{k:v for k,v in run.items() if k!='lookup'})
                files={'/':'index.html','/app.js':'app.js','/style.css':'style.css'}
                if parsed.path not in files:return self.send_json(404,{'error':'Not found'})
                file=ROOT/'dashboard'/files[parsed.path];data=file.read_bytes();self.send_response(200)
                self.send_header('Content-Type',{'html':'text/html; charset=utf-8','js':'text/javascript','css':'text/css'}[file.suffix[1:]])
                self.send_header('Content-Length',str(len(data)));self.send_header('Cache-Control','no-store')
                self.send_header('Content-Security-Policy',"default-src 'self'; script-src 'self'; style-src 'self'; connect-src 'self'; object-src 'none'; frame-ancestors 'none'")
                self.end_headers();self.wfile.write(data)
            except (KeyError,ValueError,OSError) as error:self.send_json(400,{'error':str(error)})
        def do_POST(self):
            nonlocal workspace
            origin=self.headers.get('Origin')
            if not self.valid_host() or (origin and origin!=f'http://127.0.0.1:{self.server.server_port}') or self.headers.get('X-Review-Token')!=workspace.token:
                return self.send_json(403,{'error':'Review session invalid'})
            try:
                if self.path=='/api/refresh':
                    with refresh_lock:
                        replacement=workspace.refreshed()
                        workspace=replacement
                    return self.send_json(200,workspace.manifest())
                if self.path!='/api/reviews':return self.send_json(404,{'error':'Not found'})
                length=int(self.headers.get('Content-Length','0'))
                if not 0<length<=32768:raise ValueError('Invalid request size')
                body=json.loads(self.rfile.read(length));return self.send_json(200,workspace.save(body))
            except (KeyError,ValueError,TypeError,OSError) as error:self.send_json(400,{'error':str(error)})
    return Handler


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--port',type=int,default=8879)
    parser.add_argument('--data',type=Path,default=ROOT/'data/review');parser.add_argument('--database',type=Path,default=ROOT/'data/reviews.sqlite3')
    args=parser.parse_args();workspace=Workspace(args.data,args.database)
    server=ThreadingHTTPServer(('127.0.0.1',args.port),handler(workspace))
    print(f'{APP} {VERSION}: http://127.0.0.1:{server.server_port}/',flush=True);server.serve_forever()


if __name__=='__main__':main()
