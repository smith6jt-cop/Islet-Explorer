"""Regression checks for source separation, signal labels, and boundary geometry."""
import sys
from pathlib import Path

import numpy as np
import pytest
import tifffile
from shapely.geometry import Polygon, box

sys.path.insert(0,str(Path(__file__).parents[1]/'cohort_expansion'))
from refine_endocrine_pixel_classifier import (
    ORIGINAL_DONORS,balanced_sample,clean_regions,normalization,normalize_features,
    rasterize,read_endocrine,refuse_invalid_samples,training_masks,validate_training_donors,
)


def test_target_donor_cannot_enter_training():
    refs=[{'case_id':d} for d in ORIGINAL_DONORS]
    validate_training_donors(refs)
    with pytest.raises(ValueError):validate_training_donors(refs+[{'case_id':'78'}])
    with pytest.raises(ValueError):validate_training_donors(refs[:-1]+[{'case_id':'78'}])


@pytest.mark.parametrize('byteorder',['>','<'])
def test_original_and_new_tiff_byte_orders_give_the_same_pixels(tmp_path,byteorder):
    raw=np.arange(3*32*40,dtype=np.uint16).reshape(3,32,40)
    path=tmp_path/'source.tif'
    tifffile.imwrite(path,raw,byteorder=byteorder,metadata={'axes':'CYX'},photometric='minisblack')
    spec=dict(source_image=str(path),source_channel_names=['INS','GCG','SST'],image_shape_cyx=list(raw.shape))
    actual,_=read_endocrine(spec)
    expected=raw.astype(np.float32).reshape(3,16,2,20,2).mean(axis=(2,4)).transpose(1,2,0)
    np.testing.assert_array_equal(actual,expected)


def test_invalid_draft_samples_cannot_be_reused(tmp_path):
    (tmp_path/'quality_status.json').write_text('{"status":"invalid_image_reader"}')
    with pytest.raises(ValueError):refuse_invalid_samples(tmp_path)


def test_dark_sam_interiors_do_not_become_foreground_labels():
    mask=np.zeros((320,400),bool);mask[20:220,20:220]=True
    signal=np.ones((320,400,3),np.float32);signal[30:140,30:140,0]=100
    positive,negative,interior,background,_=training_masks(mask,np.ones_like(mask),signal)
    assert not positive[180,180] and interior[180,180]
    assert not negative[180,180]
    assert positive[80,80] and negative[270,300] and background[270,300]


def test_raster_centers_holes_and_overlapping_objects():
    ring=Polygon([(0,0),(20,0),(20,20),(0,20)],holes=[[(4,4),(16,4),(16,16),(4,16)]])
    raster=rasterize([ring,box(8,8,12,12)],(10,10),downsample=2)
    assert raster[0,0] and raster[4,4]
    assert not raster[2,2]
    assert raster.sum()==68


def test_normalization_handles_acquisition_scale_without_labels():
    rng=np.random.default_rng(1);raw=rng.uniform(20,200,(30,40,3)).astype('float32')
    tissue=np.ones((30,40),bool)
    f=np.concatenate([raw,raw*1.1],axis=-1)
    a=normalize_features(f,normalization(raw,tissue))
    scale=np.array([2,4,8]);offset=np.array([100,20,40])
    b=normalize_features(f*np.tile(scale,2)+np.tile(offset,2),normalization(raw*scale+offset,tissue))
    np.testing.assert_allclose(a,b,atol=2e-6)


def test_bounded_morphology_preserves_large_gaps_and_keeps_small_foci():
    mask=np.zeros((140,200),bool);mask[20:110,20:110]=True
    mask[35:39,35:39]=False  # 16 um2 enclosed hole may be filled
    mask[55:80,55:80]=False  # 625 um2 void must remain
    mask[30:75,130:170]=True;mask[120:123,180:183]=True
    tissue=np.ones_like(mask);tissue[:,:22]=False
    regions,foci=clean_regions(mask,tissue,pixel_um=1.)
    assert regions[36,36] and not regions[65,65]
    assert not regions[50,120] and not regions[:,:22].any()
    assert foci[121,181] and not regions[121,181]


def test_group_sampling_is_deterministic_and_covers_small_regions():
    mask=np.ones((10,100),bool);groups=np.zeros(mask.shape,int);groups[:,:3]=1
    a=balanced_sample(mask,groups,100,np.random.default_rng(10))
    b=balanced_sample(mask,groups,100,np.random.default_rng(10))
    np.testing.assert_array_equal(a,b)
    assert len(np.unique(a))==100 and (groups.ravel()[a]==1).sum()==30
    assert balanced_sample(np.zeros_like(mask),groups,10,np.random.default_rng(1)).size==0
