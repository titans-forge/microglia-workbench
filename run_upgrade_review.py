#!/usr/bin/env python3
"""Regenerate rankings and sensitivity comparisons from cached public inputs."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))
from microglia_repurposing.context_engine import run_context_engine
from microglia_repurposing.engine import run_engine
from microglia_repurposing.evaluate import load_gold_standard


def compare(left, right):
    a = left['_ranked']; b = right['_ranked']
    common = a.index.intersection(b.index)
    common = [k for k in common if a.loc[k,'status'] == 'scored' and b.loc[k,'status'] == 'scored']
    top_a = list(a[a.status == 'scored'].head(25).index)
    top_b = list(b[b.status == 'scored'].head(25).index)
    changes = []
    ra = {key:i+1 for i,key in enumerate(a.index)}; rb = {key:i+1 for i,key in enumerate(b.index)}
    for key in dict.fromkeys(top_a+top_b):
        changes.append({'pert_id':key, 'name':str(a.loc[key,'pert_iname']) if key in a.index else str(b.loc[key,'pert_iname']),
                        'left_rank':ra.get(key),'right_rank':rb.get(key),
                        'left_status':str(a.loc[key,'status']) if key in a.index else 'missing',
                        'right_status':str(b.loc[key,'status']) if key in b.index else 'missing'})
    # Pearson of average ranks = Spearman, without optional scipy dependency.
    correlation = a.loc[common,'connectivity_mean'].rank().corr(b.loc[common,'connectivity_mean'].rank()) if len(common)>1 else None
    return {'common_scored':len(common),'score_rank_correlation':correlation,'top25_overlap':len(set(top_a)&set(top_b)),
            'top_changes':changes, 'interpretation':'Sensitivity within one source dataset; not independent validation.'}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent/'data/review')
    ap.add_argument('--n-perm',type=int,default=2000)
    args=ap.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    summary={'n_perm':args.n_perm,'runs':{},'comparisons':{}}
    results={}
    scenarios=[('context_mix',{'model':'mix','min_hits':3}),('context_ind',{'model':'ind','min_hits':3}),('context_mix_overlap5',{'model':'mix','min_hits':5})]
    for name,config in scenarios:
        start=time.perf_counter()
        result=run_context_engine(n_perm=args.n_perm,output=args.output_dir/(name+'.json'),**config)
        results[name]=result
        summary['runs'][name]={'seconds':round(time.perf_counter()-start,3),'counts':result['context_panel']['counts'],
                              'top':result['top'][:5], 'statistics':result['statistics']}
        (args.output_dir/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
        print(name, summary['runs'][name]['seconds'],summary['runs'][name]['counts'],flush=True)
    summary['comparisons']['model_ind_vs_mix']=compare(results['context_ind'],results['context_mix'])
    summary['comparisons']['overlap3_vs_5']=compare(results['context_mix'],results['context_mix_overlap5'])
    # Read-only reuse of per-context scores permits filter/inspection without rerunning nulls.
    base=results['context_mix']; per=base['_per_context'].copy()
    per['pert_id']=[m.pert_id for m in base['_meta']];per['cell_id']=[m.cell_id for m in base['_meta']]
    per['time']=[m.time for m in base['_meta']];per['dose']=[m.dose for m in base['_meta']]
    per.rename_axis('context_id').reset_index().to_csv(args.output_dir/'context_mix_contexts.csv',index=False)
    # Consensus inputs have no cell axis; report their rankings separately.
    names=list(load_gold_standard()['cmap_name'])
    for name,panel in [('consensus_controls',names),('consensus_all',None)]:
        start=time.perf_counter()
        result=run_engine(names=panel,n_perm=args.n_perm,output=args.output_dir/(name+'.json'))
        summary['runs'][name]={'seconds':round(time.perf_counter()-start,3),'counts':result['drug_panel']['counts'],'top':result['top'][:5]}
        print(name,summary['runs'][name]['seconds'],flush=True)
    (args.output_dir/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
    from microglia_repurposing.review import write_review
    write_review(args.output_dir)
    print('Review:',args.output_dir/'index.html',flush=True)


if __name__ == '__main__':
    main()
