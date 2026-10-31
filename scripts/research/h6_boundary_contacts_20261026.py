#!/usr/bin/env python3
"""Use existing cell polygons for the original H6 contact-definition sensitivity."""
from pathlib import Path
import json,mmap,re
import numpy as np
import pandas as pd
from shapely.geometry import shape
from shapely.strtree import STRtree
from h6_permutation_measurements_20261026 import evaluate_graph

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/h5_h10_original_plan_20261026/h6'
SOURCE=ROOT/'data/development/uniform-20donor-20261015/cell_geojson'

def selected_geometries(path,wanted):
    geometries={};pattern=re.compile(rb'\{"type":"Feature","id":"([^"]+)"')
    with path.open('rb') as fh,mmap.mmap(fh.fileno(),0,access=mmap.ACCESS_READ) as mm:
        previous=None
        for match in pattern.finditer(mm):
            if previous is not None:
                start,uid=previous
                if uid in wanted:
                    feature=json.JSONDecoder().raw_decode(mm[start:match.start()].decode())[0]
                    geometries[uid]=shape(feature['geometry'])
            previous=(match.start(),match.group(1).decode())
        if previous is not None and previous[1] in wanted:
            feature=json.JSONDecoder().raw_decode(mm[previous[0]:].decode())[0]
            geometries[previous[1]]=shape(feature['geometry'])
    assert set(geometries)==wanted
    return geometries

def main():
    graphout=OUT/'boundary_graphs_by_donor';graphout.mkdir(exist_ok=True)
    rows=[];nulls=[];ids=[]
    for donor in sorted((OUT/'graphs_by_donor').iterdir()):
        loaded=[(p,np.load(p)) for p in sorted(donor.glob('*.npz'))]
        wanted={v for _,g in loaded for v in g['object_ids']}
        path=list(SOURCE.glob(f'*_{int(donor.name):04d}.geojson'))[0]
        geoms=selected_geometries(path,wanted)
        for path,g in loaded:
            polys=[geoms[v] for v in g['object_ids']]
            pairs=STRtree(polys).query(polys,predicate='intersects').T
            edges=pairs[pairs[:,0]<pairs[:,1]].astype(np.int32)
            record={'islet_id':path.stem,'case_id':int(donor.name),'boundary_edges':len(edges),
                    'invalid_source_polygons':sum(not p.is_valid for p in polys)}
            if len(edges):
                measurement,null=evaluate_graph((path.stem,int(donor.name),g['labels'],edges))
                record.update(measurement);nulls.append(null);ids.append(path.stem)
                target=graphout/donor.name;target.mkdir(exist_ok=True)
                np.savez_compressed(target/path.name,edges=edges,labels=g['labels'],object_ids=g['object_ids'])
            rows.append(record)
        print(f'Boundary graphs: donor {donor.name}, {len(loaded)} objects',flush=True)
    d=pd.DataFrame(rows);d.to_csv(OUT/'boundary_contact_measurements.csv',index=False)
    if nulls:np.savez_compressed(OUT/'boundary_permutation_null_draws.npz',islet_ids=np.array(ids),draws=np.stack(nulls))
    print('Boundary-edge objects',int(d.boundary_edges.gt(0).sum()),'of',len(d),flush=True)

if __name__=='__main__':main()
