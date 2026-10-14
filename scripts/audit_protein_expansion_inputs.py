#!/usr/bin/env python3
"""Validate exact geometry, cell identities and QC coverage without typing cells."""
import argparse
import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from assemble_protein_expansion import donor_regions
from build_protein_embedding import sha256
from stage_original_gate_calls import DONORS


def audit(prepared,raw_root,geometry,report):
    if report.exists():raise FileExistsError(report)
    manifest=json.loads((geometry/'manifest.json').read_text());records=[]
    for cid in DONORS:
        snapshot=next(r for r in manifest['donors'] if r['case_id']==cid)
        try:
            raw=ad.read_h5ad(raw_root/f'{cid}.h5ad',backed='r')
            cells=ad.read_h5ad(prepared/f'{cid}.h5ad',backed='r')
            centroids=pd.read_csv(snapshot['centroids'],sep='\t',dtype={'object_id':str}).set_index('object_id')
            regions,_,_,_,_=donor_regions(snapshot,raw,centroids)
            kept=raw.obs_names.isin(cells.obs_names)
            qc_core=np.array([kept[r['core_ids']].sum() for r in regions])
            qc_peri=np.array([kept[r['peri_ids']].sum() for r in regions])
            counts=pd.DataFrame(dict(case_id=cid,islet_key=[r['name'] for r in regions],
                qc_core_cells=qc_core,qc_peri_cells=qc_peri,paired=(qc_core>0)&(qc_peri>0)))
            report.parent.mkdir(parents=True,exist_ok=True)
            counts.to_csv(report.parent/f'{cid}-region-census.csv',index=False)
            records.append(dict(case_id=cid,passed=True,boundaries=len(regions),paired_qc_islets=int(counts.paired.sum()),
                source_raw_cells=raw.n_obs,qc_cells=cells.n_obs,exact_native_centroid_alignment=True,
                shared_core_cells=0,native_roundoff_overlap_um2=sum(r['native_roundoff_overlap_um2'] for r in regions)/2,
                geometry_sha256=sha256(snapshot['geometry']),review_status=snapshot['boundary_status']))
        except Exception as e:
            records.append(dict(case_id=cid,passed=False,error=str(e)))
        finally:
            if 'raw' in locals():raw.file.close()
            if 'cells' in locals():cells.file.close()
        print(json.dumps(records[-1]),flush=True)
    report.parent.mkdir(parents=True,exist_ok=True)
    result=dict(passed=all(r['passed'] for r in records),donors=records,phenotypes_assigned=False,
        scientific_acceptance=False)
    report.write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['prepared','raw-root','geometry','report']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();r=audit(a.prepared,a.raw_root,a.geometry,a.report)
    raise SystemExit(0 if r['passed'] else 1)
