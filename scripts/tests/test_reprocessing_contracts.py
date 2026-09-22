"""Regression coverage for failures found during the real 21-donor rebuild."""
import importlib.util
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest
import shapely
import tifffile

RESEARCH=Path(__file__).resolve().parents[1]/'research'
sys.path.insert(0,str(RESEARCH))
from common import sha256,write_json
from build_candidate import positive_area
from reconcile_channel_identity import reconcile
from rebuild_spatial import overlap_flags,duct_distance


def channel_fixture(tmp_path,converted_reorder=False):
    root=tmp_path/'candidate';root.mkdir()
    rng=np.random.default_rng(321);native=rng.uniform(10,10000,(120,2))
    pixels=np.arange(2*16*16,dtype=np.uint16).reshape(2,16,16)
    source=tmp_path/'frozen.ome.tiff'
    tifffile.imwrite(source,pixels,ome=True,metadata={'axes':'CYX','Channel':{'Name':['HLADR','PDPN']}})
    image=root/'test.ome.tiff';os.link(source,image)
    raw=root/'preprocessing/cells/donor_id=112/data_0.parquet';raw.parent.mkdir(parents=True)
    ids=np.array([f'uuid{i}' for i in range(len(native))])
    pd.DataFrame({'object_id':ids,'HLADR':native[:,0],'PDPN':native[:,1]}).to_parquet(raw,index=False)
    intermediate=root/'preprocessing/redsea_scratch/intermediates';intermediate.mkdir(parents=True)
    values=native if converted_reorder else native[:,::-1]
    np.savez(intermediate/'112.npz',oids=ids,cols=['HLADR','PDPN'],sizes=np.ones(len(ids)),data=values)
    corrected=root/'corrected/112-cells.parquet';corrected.parent.mkdir()
    pd.DataFrame({'object_id':ids,'case_id':'112','HLADR':values[:,0]*.8,'PDPN':values[:,1]*.8}).to_parquet(corrected,index=False)
    spec=root/'source.json';write_json(spec,{'source_channel_names':['PDPN','HLADR']})
    write_json(root/'preprocessing/checkpoints/112.json',{'channel_record':{'channel_names':['HLADR','PDPN']}})
    return root,source,image,native,pixels,{'case_id':'112','image':str(image),'boundary_source':str(spec)}


def test_relabeling_preserves_frozen_pixels_and_is_idempotent(tmp_path):
    root,source,image,native,pixels,row=channel_fixture(tmp_path)
    original_hash=sha256(source)
    result=reconcile(root,row)
    assert result['permutation_required']
    assert sha256(source)==original_hash
    assert not os.path.samefile(source,image)
    np.testing.assert_array_equal(tifffile.imread(image),pixels)
    corrected=pd.read_parquet(root/'corrected/112-cells.parquet')
    np.testing.assert_allclose(corrected[['HLADR','PDPN']],native*.8)
    assert reconcile(root,row)==result


def test_legitimate_pixel_reordering_is_not_mislabeled(tmp_path):
    root,source,image,native,pixels,row=channel_fixture(tmp_path,converted_reorder=True)
    original=sha256(image)
    result=reconcile(root,row)
    assert not result['permutation_required']
    assert result['conversion_reordered_channels']
    assert sha256(image)==original and os.path.samefile(source,image)


def test_recovery_after_image_metadata_write_before_table_relabel(tmp_path):
    root,source,image,native,pixels,row=channel_fixture(tmp_path)
    with tifffile.TiffFile(image) as tf:old_xml=tf.ome_metadata
    archive=root/'provenance/channel_identity_before/112';archive.mkdir(parents=True)
    (archive/'original_ome.xml').write_text(old_xml)
    fresh=root/'recovered.ome.tiff'
    tifffile.imwrite(fresh,pixels,ome=True,metadata={'axes':'CYX','Channel':{'Name':['PDPN','HLADR']}})
    fresh.replace(image)
    result=reconcile(root,row)
    assert result['permutation_required']
    np.testing.assert_allclose(pd.read_parquet(root/'corrected/112-cells.parquet')[['HLADR','PDPN']],native*.8)


def test_absent_structure_is_missing_but_assessed_empty_is_zero():
    region=shapely.box(0,0,10,10)
    assert np.isnan(positive_area(region,None))
    empty=(np.array([],dtype=object),shapely.STRtree([]))
    assert positive_area(region,empty)==0


def test_spatial_overlap_depends_on_exact_region_not_larger_sensitivity_band():
    cores=np.array([shapely.box(40,40,50,50),shapely.box(95,40,105,50)],dtype=object)
    tissue=shapely.box(0,0,150,150)
    assert not overlap_flags(cores,1).any()
    union=shapely.union_all(cores)
    primary=np.array([g.buffer(20).intersection(tissue).difference(union) for g in cores],dtype=object)
    sensitivity=np.array([g.buffer(30).intersection(tissue).difference(union) for g in cores],dtype=object)
    assert not overlap_flags(primary,1).any()
    assert overlap_flags(sensitivity,1).all()


def test_duct_distance_censors_unobserved_field_without_substituting_zero():
    core=shapely.box(20,40,30,50);tissue=shapely.box(0,0,100,100)
    def index(geoms):return np.array(geoms,dtype=object),shapely.STRtree(geoms)
    result=duct_distance(core,tissue,index([shapely.box(70,40,75,50)]),1)
    assert result['duct_distance_censored'] and np.isnan(result['duct_distance_um'])
    assert result['duct_distance_observed_um']==40 and result['duct_distance_lower_bound_um']==20
    assert duct_distance(core,tissue,index([shapely.box(35,40,38,50)]),1)['duct_distance_um']==5
    assert duct_distance(core,tissue,index([]),1)['duct_distance_lower_bound_um']==20
    assert np.isnan(duct_distance(core,tissue,None,1)['duct_distance_lower_bound_um'])
