#!/usr/bin/env python3
"""Reconcile candidate OME labels to UUID-matched acquisition measurements.

No pixel intensities are changed. The two channels may be relabeled only when
the source acquisition and matched raw-cell means independently support the
same exact two-channel permutation. Source/frozen images remain unchanged.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
import os
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET
for _key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ.setdefault(_key,'1')
import numpy as np
import pandas as pd
import tifffile
from common import MARKERS,sha256,write_json


def match_correlations(root,cid):
    raw=pd.read_parquet(root/'preprocessing/cells'/f'donor_id={cid}'/'data_0.parquet').set_index('object_id')
    with np.load(root/'preprocessing/redsea_scratch/intermediates'/f'{cid}.npz',allow_pickle=True) as z:
        idx=pd.Index(z['oids']).get_indexer(raw.index)
        if (idx<0).any():raise ValueError('Missing native UUID')
        cols=list(z['cols']);sizes=z['sizes'][idx]
        values=z['data'][idx][:,[cols.index(m) for m in ['HLADR','PDPN']]]/np.where(sizes>0,sizes,np.nan)[:,None]
    matrix=[]
    for marker in ['HLADR','PDPN']:
        a=raw[marker].to_numpy(float);record=[]
        for b in values.T:
            good=np.isfinite(a)&np.isfinite(b)
            record.append(float(np.corrcoef(a[good],b[good])[0,1]))
        matrix.append(record)
    return np.asarray(matrix)


def reconcile(root,row):
    cid=row['case_id'];path=Path(row['image']);spec=json.loads(Path(row['boundary_source']).read_text())
    check=root/'preprocessing/channel_reconciliation'/f'{cid}.json'
    if check.exists():
        record=json.loads(check.read_text())
        if sha256(root/'corrected'/f'{cid}-cells.parquet')!=record['corrected_sha256']:
            raise ValueError('Reconciled corrected table changed')
        return record
    with tifffile.TiffFile(path) as tf:
        xml=tf.ome_metadata
        tree=ET.fromstring(xml)
        channels=tree.findall('.//{*}Channel');names=[c.attrib['Name'].replace('PGP9_5','PGP9.5') for c in channels]
    expected=spec['source_channel_names']
    archive=root/'provenance/channel_identity_before'/cid
    archived_xml=archive/'original_ome.xml'
    prior_xml=archived_xml.read_text() if archived_xml.exists() else xml
    prior_names=[c.attrib['Name'] for c in ET.fromstring(prior_xml).findall('.//{*}Channel')]
    current_indices=[prior_names.index(m) for m in ['HLADR','PDPN']]
    image_already_reconciled=[names.index(m) for m in ['HLADR','PDPN']]==[expected.index(m) for m in ['HLADR','PDPN']]
    expected_indices=[expected.index(m) for m in ['HLADR','PDPN']]
    matrix=match_correlations(root,cid)
    cross_supported=matrix[0,1]>.98 and matrix[1,0]>.98 and matrix[0,1]>matrix[0,0]+.2 and matrix[1,0]>matrix[1,1]+.2
    direct_supported=matrix[0,0]>.9 and matrix[1,1]>.9 and matrix[0,0]>matrix[0,1]+.2 and matrix[1,1]>matrix[1,0]+.2
    swap=cross_supported
    if swap and current_indices!=expected_indices[::-1]:raise ValueError(f'{cid}: unexpected acquisition permutation')
    if not swap and not direct_supported:
        raise ValueError(f'{cid}: independent UUID measurements do not confirm the acquisition swap: {matrix}')
    if not swap:
        return dict(case_id=cid,permutation_required=False,raw_pixel_correlations=matrix.tolist(),
                    image_channel_names=names,
                    conversion_reordered_channels=current_indices!=expected_indices,
                    interpretation='Current OME channel labels match independent native UUID measurements; a source acquisition order difference alone is not an error.')
    oldstat=path.stat()
    archive.mkdir(parents=True,exist_ok=True)
    # Independent new inode before changing embedded TIFF metadata. Replacing
    # the candidate hardlink must never modify any frozen/source image.
    temporary=path.with_name(path.name+'.reconciled')
    original_xml=prior_xml
    token='__CODEX_HLADR_PDPN_SWAP__'
    xml=prior_xml.replace('Name="HLADR"','Name="'+token+'"').replace('Name="PDPN"','Name="HLADR"').replace('Name="'+token+'"','Name="PDPN"')
    if len(xml.encode())!=len(original_xml.encode()):raise ValueError('Unexpected OME description length change')
    if not image_already_reconciled:
        print(f'{cid}: copying an independent image for metadata reconciliation',flush=True)
        shutil.copyfile(path,temporary)
        if os.path.samefile(path,temporary):raise ValueError('Candidate image copy is not independent')
        tifffile.tiffcomment(temporary,xml)
    with tifffile.TiffFile(path if image_already_reconciled else temporary) as tf:
        checknames=[c.attrib['Name'] for c in ET.fromstring(tf.ome_metadata).findall('.//{*}Channel')]
        if [checknames.index(m) for m in ['HLADR','PDPN']]!=expected_indices:raise ValueError('TIFF metadata roundtrip failed')
    if path.stat().st_mtime_ns!=oldstat.st_mtime_ns or path.stat().st_size!=oldstat.st_size:raise ValueError('Source image changed during copy')
    if not archived_xml.exists():archived_xml.write_text(original_xml)
    if not image_already_reconciled:temporary.replace(path)
    corrected=root/'corrected'/f'{cid}-cells.parquet'
    if not (archive/corrected.name).exists():shutil.copy2(corrected,archive/corrected.name)
    frame=pd.read_parquet(archive/corrected.name);frame[['HLADR','PDPN']]=frame[['PDPN','HLADR']].to_numpy()
    temp=corrected.with_suffix('.new.parquet');frame.to_parquet(temp,index=False);temp.replace(corrected)
    checkpoint=root/'preprocessing/checkpoints'/f'{cid}.json';record=json.loads(checkpoint.read_text())
    if not (archive/'redsea_checkpoint.json').exists():shutil.copy2(checkpoint,archive/'redsea_checkpoint.json')
    record['corrected_sha256']=sha256(corrected)
    record['image_mtime_ns']=path.stat().st_mtime_ns
    record['channel_record']['channel_names']=checknames
    record['channel_identity_revision']='HLADR_PDPN_native_acquisition_reconciliation_v1'
    record['redsea_pixel_rerun_required']=False
    record['relabeling_equivalence']='Scalar REDSEA correction operates independently on each channel; swapping output columns exactly matches rerunning with corrected channel names.'
    write_json(checkpoint,record)
    result=dict(case_id=cid,permutation_required=True,source_indices=current_indices,
        corrected_indices=expected_indices,raw_pixel_correlations=matrix.tolist(),
        acquisition_metadata=str(Path(row['boundary_source'])),acquisition_metadata_sha256=sha256(row['boundary_source']),
        image=str(path),image_sha256=sha256(path),image_bytes=path.stat().st_size,
        image_mtime_ns=path.stat().st_mtime_ns,corrected_sha256=record['corrected_sha256'],
        source_image_inode_preserved=True,pixels_modified=False,metadata_only_change=True,
        original_inputs_preserved=True,image_channel_names=checknames,
        redsea_channelwise_equivalence=record['relabeling_equivalence'])
    write_json(check,result)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-root',type=Path,required=True);p.add_argument('--workers',type=int,default=4)
    a=p.parse_args();root=a.input_root.resolve(strict=True)
    rows=pd.read_csv(root/'cohort.csv',dtype=str).to_dict('records');records=[]
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(reconcile,root,row) for row in rows]):
            r=f.result();records.append(r);print(r['case_id'],r['permutation_required'],flush=True)
    records.sort(key=lambda r:int(r['case_id']))
    write_json(root/'preprocessing/channel_reconciliation.json',dict(status='channel_identity_reconciled',
        donors=records,corrected_donors=sum(r['permutation_required'] for r in records),
        interpretation='Metadata labels reconciled to independent native UUID measurements and acquisition metadata. No pixel intensities or original image files modified.'))


if __name__=='__main__':main()
