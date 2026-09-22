#!/usr/bin/env python3
"""Rebuild canonical UUID cell inputs and run manifest-bound REDSEA for 21 donors.

The earlier H5AD row-index identifiers are archival, not QuPath object IDs.
The original native measurement table and detection UUID define the new cell
universe. Every older row table is preserved alongside a reconciliation report.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import subprocess
import sys
import time

for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ.setdefault(name,'1')
REPO=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'external/Phenocycler_Analysis'))
sys.path.insert(0,str(REPO/'scripts/codex_phenocycler'))
from phenocycler.config import PipelineConfig
import numpy as np
import pandas as pd
import tifffile
from codex_panel import ome_channel_names, expected_channels, DONORS
from common import MARKERS, sha256, write_json


def prepare(args):
    root=args.input_root.resolve(strict=True)
    cohort=pd.read_csv(root/'cohort.csv',dtype=str)
    records=[]
    for row in cohort.to_dict('records'):
        cid=row['case_id'];path=root/'raw'/f'{cid}-cells.parquet'
        archive=root/'provenance/archival_raw'/path.name
        archive.parent.mkdir(parents=True,exist_ok=True)
        if archive.exists():
            raise FileExistsError('Canonical preparation is immutable; use its recorded outputs: '+str(archive))
        os.link(path,archive)
        prior=pd.read_parquet(archive)
        if cid in DONORS:
            source=args.native_root/'cells'/f'donor_id={cid}'/'data_0.parquet'
            native=pd.read_parquet(source)
            raw=native.rename(columns={'centroid_x_um':'X_centroid','centroid_y_um':'Y_centroid',
                'cell_area_um2':'Cell Area','nucleus_area_um2':'Nucleus Area',
                'classification':'archival_phenotype','parent':'archival_parent'}).copy()
            raw['case_id']=cid
            geometry=args.native_root/'redsea_scratch/geojson'/f'cells__{int(cid):04d}.geojson'
            native_metadata=dict(source=str(source),source_sha256=sha256(source),
                geometry=str(geometry.resolve(strict=True)),geometry_sha256=sha256(geometry))
        else:
            raw=prior.copy();raw['object_id']=raw.object_id.str.removeprefix(cid+':')
            geometry=root/'native_cells'/cid/'cells.geojson'
            proof=root/'native_cells'/cid/'manifest.json'
            native_metadata=json.loads(proof.read_text())
            if sha256(geometry)!=native_metadata['geometry_sha256']:
                raise ValueError(f'{cid}: native polygon export changed')
            centroids=pd.read_csv(root/'native_cells'/cid/'native_centroids.tsv',sep='\t',dtype={'object_id':str})
            if set(raw.object_id)!=set(centroids.object_id):
                raise ValueError(f'{cid}: native and raw UUID universes differ')
            raw=raw.merge(centroids,on='object_id',validate='one_to_one')
            raw['X_centroid']=raw.x_source_px*.3774;raw['Y_centroid']=raw.y_source_px*.3774
        if raw.object_id.duplicated().any() or raw.object_id.isna().any():
            raise ValueError(f'{cid}: missing/duplicate native UUID')
        raw['marker_measurements_complete']=np.isfinite(raw[MARKERS].to_numpy(float)).all(axis=1)
        raw['measurement_exclusion_reason']=np.where(raw.marker_measurements_complete,'','incomplete_raw_marker_measurements')
        raw['image']=Path(row['image']).name
        intern=root/'preprocessing/cells'/f'donor_id={cid}'/'data_0.parquet'
        intern.parent.mkdir(parents=True,exist_ok=True);raw.to_parquet(intern,index=False)
        raw['object_id']=cid+':'+raw.object_id.astype(str)
        # Atomic replacement breaks the staging hardlink; the previous candidate
        # and archived H5AD-derived raw inputs remain byte-for-byte unchanged.
        temp=path.with_suffix('.new.parquet');raw.to_parquet(temp,index=False);temp.replace(path)
        records.append(dict(case_id=cid,image=row['image'],geometry=str(geometry.resolve(strict=True)),
            source_raw=str(intern),raw_sha256=sha256(intern),native= native_metadata,
            archived_rows=len(prior),canonical_rows=len(raw),
            incomplete_raw_cells=int((~raw.marker_measurements_complete).sum()),
            raw_universe_reason='native UUID detections retain objects omitted from earlier H5AD' if cid in DONORS else 'unchanged UUID universe'))
    write_json(root/'preprocessing/inputs.json',dict(donors=records,markers=MARKERS,
        archived_inputs_preserved=True,raw_ids='case_id:QuPath_UUID',
        source_implementation_sha256=sha256(REPO/'external/Phenocycler_Analysis/phenocycler/redsea.py')))
    pd.DataFrame([{k:v for k,v in r.items() if k not in ['native']} for r in records]).to_csv(root/'preprocessing/cell_universe.csv',index=False)
    print('Prepared',sum(r['canonical_rows'] for r in records),'native cells',flush=True)


def correct_one(args):
    from phenocycler import redsea
    root=args.input_root.resolve(strict=True)
    input_manifest=json.loads((root/'preprocessing/inputs.json').read_text())
    row=next(r for r in input_manifest['donors'] if r['case_id']==args.donor)
    cid=args.donor;image=Path(row['image']);image_stat=image.stat()
    if sha256(row['source_raw'])!=row['raw_sha256']:
        raise ValueError('Raw measurements changed')
    geometry_hash=sha256(row['geometry'])
    expected=row['native']['geometry_sha256']
    if geometry_hash!=expected:raise ValueError('Native cell polygons changed')
    def resolve(cfg,donor):
        if donor!=cid:raise ValueError('Cross-donor REDSEA request')
        return image.name,Path(row['geometry']),image
    channel_record={}
    def channels(path,level):
        if Path(path)!=image or level!=0:raise ValueError('Only declared full-resolution image is supported')
        tf=tifffile.TiffFile(path);series=tf.series[0]
        sidecar=path.with_suffix('').with_suffix('.ome.xml')
        xml=sidecar.read_text() if sidecar.exists() else tf.ome_metadata
        if not xml:raise ValueError('Named OME channel metadata required')
        names=ome_channel_names(xml)
        if len(names)!=len(set(names)) or not set(MARKERS).issubset(names):
            raise ValueError('Duplicate or missing panel markers')
        if series.axes!='CYX' or len(names)!=series.shape[0]:raise ValueError('Image channel/dimension mismatch')
        spec=json.loads((root/'boundaries/sources'/cid/'source.json').read_text())
        # Conversion may reorder actual channel pixels; raw acquisition order
        # alone cannot validate converted OME labels. The independent UUID/raw
        # pixel alignment audit is the final channel-identity acceptance gate.
        if tuple(series.shape[-2:])!=tuple(spec['image_shape_cyx'][-2:]):raise ValueError('Native geometry and image dimensions differ')
        channel_record.update(channel_names=names,image_shape=list(series.shape),pixel_um=.3774)
        return tf,series.pages,names,series.shape[-2:],.3774
    redsea.resolve_paths=resolve;redsea.qptiff_channels=channels
    cfg=PipelineConfig(data_dir=root/'preprocessing',images_dir=root/'images')
    params=redsea.RedseaParams(downsample=1,edge_radius=0,comp_mode=0,alpha=1.,gap_bridge=1,
                              save_intermediates=True,keep_mask=False)
    started=time.monotonic()
    redsea.process_donor(cid,cfg,params)
    corrected=pd.read_parquet(cfg.cells_redsea_dir/f'donor_id={cid}'/'data_0.parquet')
    raw=pd.read_parquet(root/'raw'/f'{cid}-cells.parquet')
    corrected['object_id']=cid+':'+corrected.object_id.astype(str)
    if set(raw.object_id)!=set(corrected.object_id):raise ValueError('Corrected cell UUID universe differs')
    measured=raw.marker_measurements_complete.to_numpy(bool)
    merged=raw.drop(columns=MARKERS).merge(corrected.drop(columns=['donor_id']),on='object_id',validate='one_to_one',sort=False)
    # Preserve the documented missing measurements, even when geometry pixels
    # could supply a new estimate. That is an explicit measurement revision.
    merged.loc[~merged.marker_measurements_complete,MARKERS]=np.nan
    finite=np.isfinite(merged[MARKERS].to_numpy(float)).all(axis=1)
    merged['corrected_measurements_complete']=finite & merged.marker_measurements_complete
    merged.loc[merged.marker_measurements_complete & ~finite,'measurement_exclusion_reason']='no_rasterized_cell_pixels'
    merged['marker_measurements_complete']=merged.corrected_measurements_complete
    dest=root/'corrected'/f'{cid}-cells.parquet'
    temp=dest.with_suffix('.new.parquet');merged.to_parquet(temp,index=False);temp.replace(dest)
    if image.stat().st_size!=image_stat.st_size or image.stat().st_mtime_ns!=image_stat.st_mtime_ns:
        raise ValueError('Image changed during correction')
    write_json(root/'preprocessing/checkpoints'/f'{cid}.json',dict(case_id=cid,status='complete',
        image=str(image),image_bytes=image_stat.st_size,image_mtime_ns=image_stat.st_mtime_ns,
        raw_sha256=row['raw_sha256'],native_geometry_sha256=geometry_hash,
        corrected_sha256=sha256(dest),cells=len(merged),measured_cells=int(finite.sum()),
        excluded=int((~merged.marker_measurements_complete).sum()),channel_record=channel_record,
        params=vars(params),algorithm_sha256=sha256(Path(redsea.__file__)),
        runtime_seconds=time.monotonic()-started))
    print(f'{cid}: correction complete',flush=True)


def sync_centroids(args):
    """Bind exact native centroids before starting expensive correction jobs."""
    root=args.input_root.resolve(strict=True)
    path=root/'preprocessing/inputs.json';manifest=json.loads(path.read_text())
    if list((root/'preprocessing/checkpoints').glob('*.json')):
        raise ValueError('Centroid dependencies must be fixed before correction starts')
    for row in manifest['donors']:
        cid=row['case_id'];folder=root/'native_cells'/cid
        proof=json.loads((folder/'manifest.json').read_text())
        centroids=folder/'native_centroids.tsv'
        if sha256(centroids)!=proof['centroid_sha256']:raise ValueError('Native centroids changed')
        native=pd.read_csv(centroids,sep='\t',dtype={'object_id':str})
        raw=pd.read_parquet(row['source_raw'])
        if not set(raw.object_id).issubset(set(native.object_id)):
            raise ValueError(f'{cid}: canonical cells absent from native hierarchy')
        # Native extra fragments remain in REDSEA rasterization; a missing
        # measurement record is never fabricated from a neighboring cell.
        native=native.set_index('object_id').reindex(raw.object_id)
        raw['x_source_px']=native.x_source_px.to_numpy();raw['y_source_px']=native.y_source_px.to_numpy()
        raw['X_centroid']=raw.x_source_px*.3774;raw['Y_centroid']=raw.y_source_px*.3774
        raw.to_parquet(row['source_raw'],index=False)
        row['raw_sha256']=sha256(row['source_raw'])
        row['native_centroids_sha256']=proof['centroid_sha256']
        row['native_source_sha256']=proof['source_sha256']
        row['native_only_objects']=len(native.index.difference(pd.Index(raw.object_id)))
        raw['object_id']=cid+':'+raw.object_id.astype(str)
        output=root/'raw'/f'{cid}-cells.parquet';temp=output.with_suffix('.new.parquet')
        raw.to_parquet(temp,index=False);temp.replace(output)
    manifest['exact_native_centroids']=True;write_json(path,manifest)
    print('Exact native centroids reconciled for',len(manifest['donors']),'donors',flush=True)


def run(args):
    root=args.input_root.resolve(strict=True)
    rows=json.loads((root/'preprocessing/inputs.json').read_text())['donors']
    if args.donor:rows=[r for r in rows if r['case_id'] in args.donor]
    def work(row):
        cid=row['case_id'];checkpoint=root/'preprocessing/checkpoints'/f'{cid}.json'
        if checkpoint.exists():
            record=json.loads(checkpoint.read_text())
            if record['raw_sha256']!=row['raw_sha256'] or sha256(root/'corrected'/f'{cid}-cells.parquet')!=record['corrected_sha256']:
                raise ValueError(f'{cid}: stale checkpoint')
            return cid
        command=[sys.executable,str(Path(__file__).resolve()),'one','--input-root',str(root),'--donor',cid]
        with (root/'logs'/f'redsea_{cid}.log').open('a') as log:
            subprocess.run(command,cwd=REPO,stdout=log,stderr=subprocess.STDOUT,check=True)
        return cid
    failed=[]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures={pool.submit(work,r):r['case_id'] for r in rows}
        for f in as_completed(futures):
            try:print('Corrected',f.result(),flush=True)
            except Exception as exc:
                failed.append(dict(case_id=futures[f],error=str(exc)));print('FAILED',futures[f],str(exc),flush=True)
    write_json(root/'preprocessing/run_status.json',dict(requested=len(rows),completed=len(rows)-len(failed),failures=failed))
    if failed:raise SystemExit(1)


def alignment_audit(args):
    """Compare every declared channel to UUID-matched native raw measurements."""
    root=args.input_root.resolve(strict=True)
    donors=json.loads((root/'preprocessing/inputs.json').read_text())['donors']
    def one(row):
        cid=row['case_id'];raw=pd.read_parquet(row['source_raw']).set_index('object_id')
        with np.load(root/'preprocessing/redsea_scratch/intermediates'/f'{cid}.npz',allow_pickle=True) as z:
            selected=pd.Index(z['oids'].astype(str)).get_indexer(raw.index)
            if np.any(selected<0):raise ValueError('UUID absent from pixel measurement')
            cols=list(z['cols']);sizes=z['sizes'][selected]
            check=json.loads((root/'preprocessing/checkpoints'/f'{cid}.json').read_text())
            if check.get('channel_identity_revision')=='HLADR_PDPN_native_acquisition_reconciliation_v1':
                a,b=cols.index('HLADR'),cols.index('PDPN');cols[a],cols[b]=cols[b],cols[a]
            values=z['data'][selected][:,[cols.index(m) for m in MARKERS]]/np.where(sizes>0,sizes,np.nan)[:,None]
        original=raw[MARKERS].to_numpy(float);records=[]
        for i,marker in enumerate(MARKERS):
            ok=np.isfinite(original[:,i])&np.isfinite(values[:,i]);a=original[ok,i];b=values[ok,i]
            corr=float(np.corrcoef(a,b)[0,1]) if len(a)>2 and a.std()>0 and b.std()>0 else None
            log_corr=float(np.corrcoef(np.log1p(a),np.log1p(b))[0,1]) if corr is not None else None
            ratio=float(np.median(b[a>0]/a[a>0])) if (a>0).any() else None
            records.append(dict(case_id=cid,marker=marker,matched_measured_cells=int(ok.sum()),
                pearson_r=corr,log1p_pearson_r=log_corr,median_pixel_to_export_ratio=ratio,
                source_export='native_cell_mean',pixel_measurement='uncorrected_full_resolution_UUID_mask_mean',
                review_required=corr is None or corr<.98 or ratio is None or not .95<=ratio<=1.05))
        return records
    rows=[]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for f in as_completed([pool.submit(one,row) for row in donors]):rows.extend(f.result())
    table=pd.DataFrame(rows).sort_values(['case_id','marker']);table.to_csv(root/'preprocessing/channel_alignment.csv',index=False)
    write_json(root/'preprocessing/channel_alignment.json',dict(comparisons=len(table),
        minimum_pearson=float(table.pearson_r.min()),review_required=int(table.review_required.sum()),
        interpretation='Correlation and scale diagnostics; not biological cell-typing validation'))
    print(table.loc[table.review_required,['case_id','marker','pearson_r','median_pixel_to_export_ratio']].to_string(index=False),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare');a.add_argument('--native-root',type=Path,required=True);a.set_defaults(func=prepare)
    b=sub.add_parser('one');b.add_argument('--donor',required=True);b.set_defaults(func=correct_one)
    c=sub.add_parser('run');c.add_argument('--donor',action='append');c.add_argument('--jobs',type=int,default=12);c.set_defaults(func=run)
    d=sub.add_parser('sync-centroids');d.set_defaults(func=sync_centroids)
    e=sub.add_parser('alignment-audit');e.add_argument('--jobs',type=int,default=4);e.set_defaults(func=alignment_audit)
    for parser in [a,b,c,d,e]:parser.add_argument('--input-root',type=Path,required=True)
    args=p.parse_args();args.func(args)


if __name__=='__main__':main()
