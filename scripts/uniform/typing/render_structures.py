#!/usr/bin/env python3
"""Render each established structure over its stain and DAPI, against its previous geometry.

Reads only. Output is a whole-image view, individual native-resolution zooms and a quantitative
before/after table. A render is not an accepted review; reviewers record their findings separately.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape

sys.path.insert(0,str(Path(__file__).resolve().parent))
import collagen_structures as C
import structure_thresholds as ST
import codex_typing as CT

def one(run, baseline, out, tag):
    dapi=C.load(run,tag,'DAPI');tissue=ST.tissue_mask(run,tag,dapi.shape)
    tissue_path=run/'inputs'/'tissue'/f'{tag}.geojson'
    tissue_geom=shapely.union_all([shape(f['geometry']) for f in json.loads(tissue_path.read_text())['features']])
    records=[]
    for cls,marker in ST.CLASSES.items():
        current=run/'structures'/cls/f'{cls}__{tag}.geojson'
        old=baseline/'structures'/cls/f'{cls}__{tag}.geojson'
        binary=ST.object_mask(current,dapi.shape);previous=ST.object_mask(old,dapi.shape)
        geoms=[shape(f['geometry']) for f in json.loads(current.read_text())['features']]
        union=shapely.union_all(geoms)
        missing=union.difference(tissue_geom).area*CT.PIXEL_UM**2
        if missing > 1e-5:raise ValueError(f'{tag}/{cls}: {missing} um2 outside final Tissue')
        C.render(tag,C.load(run,tag,marker),dapi,tissue,binary,np.zeros_like(tissue),previous,out/cls)
        records.append({'image':tag,'class':cls,'marker':marker,'objects':len(geoms),
            'union_area_um2':union.area*CT.PIXEL_UM**2,'tissue_area_um2':tissue_geom.area*CT.PIXEL_UM**2,
            'area_density':union.area/tissue_geom.area,'outside_tissue_um2':missing,
            'pixel_dice_vs_baseline':2*int((binary&previous).sum())/max(int(binary.sum())+int(previous.sum()),1),
            'structure_sha256':hashlib.sha256(current.read_bytes()).hexdigest(),
            'tissue_sha256':hashlib.sha256(tissue_path.read_bytes()).hexdigest()})
    print(tag,'three structure classes rendered',flush=True)
    return records

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--baseline',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--images',nargs='*');p.add_argument('--workers',type=int,default=2)
    a=p.parse_args();tags=a.images or [r['image_name'] for r in CT.cohort_rows()]
    a.out.mkdir(parents=True,exist_ok=True)
    with ProcessPoolExecutor(a.workers) as pool:
        records=list(pool.map(one,[a.run]*len(tags),[a.baseline]*len(tags),[a.out]*len(tags),tags))
    pd.DataFrame([r for rs in records for r in rs]).to_csv(a.out/'checks.csv',index=False)
if __name__=='__main__':main()
