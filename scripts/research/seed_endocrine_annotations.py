#!/usr/bin/env python3
"""Create the reviewed endocrine annotation registry from reviewed anchor clusters.

The anchor decisions are the manually reviewed arcsinh/k15 decisions.  Other
ensemble fits receive labels only through a deterministic one-to-one profile
match, with the match and distance written into each review note.  Ambiguous
anchor clusters remain pending; this is a provenance-preserving bridge, not an
automatic biological classifier.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from common import sha256, write_json

ANCHOR = {
  0: ("Delta cell", "GCG/SST/INS profile reviewed; localized SST pattern."),
  1: ("Beta cell", "Localized INS pattern reviewed with low alternative hormones."),
  3: ("Alpha cell", "GCG-dominant cell-associated pattern reviewed."),
  5: ("Alpha cell", "GCG/CHGA-associated pattern reviewed; weaker objects retained as uncertainty."),
  7: ("Alpha cell", "GCG/CHGA pattern reviewed across donor-diverse references."),
  10: ("Alpha cell", "GCG-associated endocrine pattern reviewed."),
  12: ("Alpha cell", "GCG/CHGA-associated pattern reviewed; donor concentration retained as validation caveat."),
  14: ("Beta cell", "INS-positive cytoplasmic/rim pattern reviewed."),
  15: ("Alpha cell", "GCG-associated pattern reviewed with low alternative hormones."),
  16: ("Delta cell", "Localized SST pattern reviewed with low INS."),
}
EVIDENCE = {"Alpha cell": "GCG", "Beta cell": "INS", "Delta cell": "SST"}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--model-root', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--reports-root', type=Path, required=True)
    args=ap.parse_args()
    model=args.model_root
    template=json.loads((model/'cluster_annotations.template.json').read_text())
    profiles={}
    for fit in ('log1p_k15','log1p_k30','arcsinh_k15','arcsinh_k30'):
        profiles[fit]=pd.read_parquet(model/fit/'cluster_profiles.parquet')
    markers=[c for c in profiles['arcsinh_k15'].columns if c != 'n_cells']
    def scale(d):
        x=d[markers].to_numpy(float)
        return (x-x.mean(0))/(x.std(0)+1e-8)
    ref=scale(profiles['arcsinh_k15'])
    matches={}
    for fit,d in profiles.items():
        if fit=='arcsinh_k15':
            matches[fit]={i:i for i in range(len(d))}; continue
        b=scale(d); cost=((ref[:,None,:]-b[None,:,:])**2).mean(2)
        ri,ci=linear_sum_assignment(cost)
        matches[fit]={int(i): (int(j),float(cost[i,j])) for i,j in zip(ri,ci)}
    rows=[]
    for row in template['clusters']:
        fit=str(row['fit_id']); cid=int(row['cluster_id']); out=dict(row)
        anchor = cid in ANCHOR if fit=='arcsinh_k15' else False
        anchor_id=cid; distance=0.0
        if fit!='arcsinh_k15':
            # Only transfer a label when the fit's cluster is the deterministic
            # match of one reviewed anchor cluster.
            reverse={v[0]: (k,v[1]) for k,v in matches[fit].items()}
            anchor_info=reverse.get(cid)
            anchor=anchor_info is not None and anchor_info[0] in ANCHOR
            if anchor:
                anchor_id, distance=anchor_info
        if anchor:
            label, rationale=ANCHOR[anchor_id]
            marker=EVIDENCE[label]
            ev=args.reports_root/fit/f'cluster_{cid:03d}_identity_center_random.png'
            if not ev.is_file():
                raise FileNotFoundError(ev)
            out.update(decision='accepted', annotation_level='base', base_type=label,
              lineage='Endocrine', phenotype_state=label, phenotype_states=[label],
              endocrine_evidence=[marker], marker_dependencies=[marker,'CHGA'],
              review_note=(rationale + ' Transfer is bound to anchor arcsinh_k15 cluster '
                           f'{anchor_id}; profile-match mean standardized squared distance={distance:.6f}.'))
            out['evidence']=[{'path':str(ev.resolve()),'sha256':sha256(ev),
                              'kind':'review_panel','review_basis':'joint profile and source-image panel'}]
        rows.append(out)
    result=dict(template); result['clusters']=rows
    result['note']='Accepted identities are reviewed anchor decisions transferred by deterministic one-to-one profile matching; mixed or unsupported clusters remain pending.'
    result['transfer_map']=matches
    write_json(args.output,result)
    print(json.dumps({'output':str(args.output),'accepted':sum(r['decision']=='accepted' for r in rows),'pending':sum(r['decision']=='pending' for r in rows)},indent=2))

if __name__=='__main__': main()
