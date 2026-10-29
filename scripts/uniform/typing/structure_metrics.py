#!/usr/bin/env python3
"""Structural union area in Tissue/endocrine objects and independent Collagen IV cell distances.

These endpoints never modify cell context, positivity or identity. Component counts are diagnostics;
collagen amount is union area divided by the region's area. Distances are centroid-to-area in um,
zero inside an object, missing when there are no objects.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape

sys.path.insert(0, str(Path(__file__).resolve().parent))
import codex_typing as CT
import collagen_structures as C

CLASSES = ('Capillary', 'Nerve', 'Lymphatic', 'Collagen IV')
MANIFEST = 'structure_metrics_manifest.json'

def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()

def input_paths(tag):
    return [f'inputs/tissue/{tag}.geojson', f'islets/{tag}.geojson', f'cells/{tag}.parquet',
            *[f'structures/{cls}/{cls}__{tag}.geojson' for cls in CLASSES]]

def output_paths(tag):
    return [f'cell_structure_distances/{tag}.parquet', f'cell_structure_distances/{tag}.json']

def file_hashes(run, paths):
    return {path: sha256(run / path) for path in paths}

def validate_manifest(run, expected_images=None):
    """Fail before assembly writes anything if geometry, cells, code or endpoint outputs are stale."""
    run = Path(run)
    path = run / MANIFEST
    if not path.is_file():
        raise ValueError('Missing structure metrics manifest; run typing/structure_metrics.py first')
    manifest = json.loads(path.read_text())
    if manifest.get('schema') != 'uniform.structure_metrics/1':
        raise ValueError('Unsupported structure metrics manifest')
    if manifest.get('metric_code_sha256') != sha256(__file__):
        raise ValueError('Structure metrics code changed; rebuild endpoints')
    if manifest.get('pixel_um') != CT.PIXEL_UM or manifest.get('classes') != list(CLASSES):
        raise ValueError('Structure metrics definitions changed; rebuild endpoints')
    tags = manifest.get('images', [])
    if len(tags) != len(set(tags)) or (expected_images is not None and set(tags) != set(expected_images)):
        raise ValueError('Structure metrics image set does not match the requested cohort')
    expected_inputs = {'structures/collagen_iv.json'} | {p for tag in tags for p in input_paths(tag)}
    expected_outputs = {'structure_region_metrics.csv'} | {p for tag in tags for p in output_paths(tag)}
    for kind, expected in [('inputs', expected_inputs), ('outputs', expected_outputs)]:
        recorded = manifest.get(kind, {})
        if set(recorded) != expected:
            raise ValueError(f'Structure metrics {kind} inventory is incomplete')
        for relative, digest in recorded.items():
            target = run / relative
            if not target.is_file() or sha256(target) != digest:
                raise ValueError(f'Stale structure metrics {kind}: {relative}')
    C.validate_current(run,tags)
    return manifest

def write_manifest(run, results, global_inputs):
    """Publish only after every result and every input has been checked again for concurrent changes."""
    tags = [item['image'] for item in results]
    inputs = dict(global_inputs)
    for item in results:
        inputs.update(item['inputs'])
    for relative, digest in inputs.items():
        if sha256(run / relative) != digest:
            raise ValueError(f'Input changed while measuring structures: {relative}')
    C.validate_current(run,tags)
    outputs = ['structure_region_metrics.csv'] + [p for tag in tags for p in output_paths(tag)]
    manifest = {'schema': 'uniform.structure_metrics/1', 'source_run': str(run.resolve()),
                'metric_code_sha256': sha256(__file__), 'pixel_um': CT.PIXEL_UM,
                'classes': list(CLASSES), 'images': tags, 'inputs': inputs,
                'outputs': file_hashes(run, outputs), 'rows': sum(len(item['records']) for item in results)}
    part = run / (MANIFEST + '.part')
    part.write_text(json.dumps(manifest, indent=2) + '\n')
    os.replace(part, run / MANIFEST)
    return manifest

def features(path):
    raw=json.loads(Path(path).read_text())
    return raw['features'] if isinstance(raw,dict) else raw

def valid(g):
    return g if g.is_valid else shapely.make_valid(g)

def area_metrics(region, geoms, tree=None):
    if region.is_empty or region.area == 0:
        return {'area_um2':0.0,'region_um2':0.0,'area_density':np.nan,'count':0}
    geoms=np.asarray(geoms,dtype=object)
    near=(tree or shapely.STRtree(geoms)).query(region,predicate='intersects') if len(geoms) else []
    intersections=shapely.intersection(geoms[near],region) if len(near) else []
    areas=shapely.area(intersections) if len(intersections) else np.array([])
    union=shapely.union_all(intersections)
    area=float(union.area)
    if area > region.area + 1e-6: raise ValueError('union area exceeds region area')
    return {'area_um2':area*CT.PIXEL_UM**2,'region_um2':region.area*CT.PIXEL_UM**2,
            'area_density':min(1.0,area/region.area),'count':int((areas>0).sum())}

def distances(xy_um, geoms):
    if not len(geoms): return np.full(len(xy_um),np.nan)
    points=shapely.points(np.asarray(xy_um)/CT.PIXEL_UM)
    pairs,values=shapely.STRtree(geoms).query_nearest(points,return_distance=True,all_matches=False)
    result=np.full(len(points),np.nan)
    result[pairs[0]]=values*CT.PIXEL_UM
    return result

def one(run, row):
    tag=row['image_name']; tissue_path=run/'inputs'/'tissue'/f'{tag}.geojson'
    inputs=file_hashes(run,input_paths(tag))
    tissue=shapely.union_all([valid(shape(f['geometry'])) for f in features(tissue_path)])
    regions=[('Tissue','tissue','Tissue',tissue)]
    for f in features(run/'islets'/f'{tag}.geojson'):
        p=f.get('properties',{}); cls=p.get('classification',{}).get('name')
        if cls in CT.ENDOCRINE_OBJECTS and f.get('geometry'):
            regions.append((p.get('name',f.get('id')),'endocrine_object',cls,valid(shape(f['geometry'])).intersection(tissue)))
    records=[]; hashes={}
    for cls in CLASSES:
        path=run/'structures'/cls/f'{cls}__{tag}.geojson'
        hashes[cls]=inputs[str(path.relative_to(run))]
        geoms=np.asarray([valid(shape(f['geometry'])).intersection(tissue) for f in features(path) if f.get('geometry')],dtype=object)
        geoms=geoms[~shapely.is_empty(geoms)]
        tree=shapely.STRtree(geoms)
        for region,rtype,object_class,g in regions:
            records.append({'image':tag,'case_id':row['case_id'],'region':region,'region_type':rtype,
                            'object_class':object_class,'class':cls,**area_metrics(g,geoms,tree)})
        if cls=='Collagen IV':
            cells=pd.read_parquet(run/'cells'/f'{tag}.parquet',columns=['Object ID','Centroid X µm','Centroid Y µm'])
            if cells['Object ID'].duplicated().any():raise ValueError(tag+': duplicate cell IDs')
            d=pd.DataFrame({'object_id':cells['Object ID'],
                            'distance_to_collagen_iv_um':distances(cells.iloc[:,1:].to_numpy(),geoms)})
            dest=run/'cell_structure_distances'/f'{tag}.parquet';dest.parent.mkdir(parents=True,exist_ok=True)
            d.to_parquet(dest,index=False)
            dest.with_suffix('.json').write_text(json.dumps({'image':tag,'objects':len(geoms),'cells':len(cells),
                'definition':'centroid to nearest Collagen IV area, um; zero inside; NA if no objects',
                'structure_sha256':hashes[cls],'tissue_sha256':inputs[str(tissue_path.relative_to(run))],
                'cells_sha256':inputs[f'cells/{tag}.parquet']},indent=2)+'\n')
    print(tag,len(records),'structure regions',flush=True)
    return {'image':tag,'records':records,'inputs':inputs}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--workers',type=int,default=3)
    a=p.parse_args();a.run=a.run.resolve();rows=CT.cohort_rows()
    C.validate_current(a.run,[r['image_name'] for r in rows])
    global_inputs=file_hashes(a.run,['structures/collagen_iv.json'])
    with ProcessPoolExecutor(a.workers) as pool:
        tables=list(pool.map(one,[a.run]*len(rows),rows))
    csv_path=a.run/'structure_region_metrics.csv'
    part=csv_path.with_suffix('.csv.part')
    pd.DataFrame([r for result in tables for r in result['records']]).to_csv(part,index=False)
    os.replace(part,csv_path)
    write_manifest(a.run,tables,global_inputs)
if __name__=='__main__':main()
