import gzip
import json
import sys
from pathlib import Path
import numpy as np
import pytest
from shapely.geometry import shape
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'research'))
from structural_masks import exact_pixel_geojson,duct_masks,threshold_proposal,accept
from common import sha256


def test_positive_pixel_polygons_preserve_holes_and_exact_area(tmp_path):
    mask=np.zeros((260,12),bool);mask[2:250,2:10]=True;mask[5:245,4:8]=False
    record=exact_pixel_geojson(mask,tmp_path/'positive.geojson.gz','ColIV')
    with gzip.open(record['path'],'rt') as stream:doc=json.load(stream)
    geoms=[shape(f['geometry']) for f in doc['features']]
    assert all(g.is_valid for g in geoms)
    assert sum(g.area for g in geoms)==int(mask.sum())
    assert all(f['properties']['classification']['name']=='ColIV' for f in doc['features'])


def test_empty_mask_has_zero_area_without_claiming_assessability(tmp_path):
    record=exact_pixel_geojson(np.zeros((5,5),bool),tmp_path/'empty.geojson.gz','Duct')
    assert record['positive_pixels']==0
    assert record['features']==0
    assert 'validated' not in record


def test_duct_envelope_never_replaces_positive_pixel_area():
    ck=np.zeros((60,60));ck[10:50,10:50]=100;ck[20:40,20:40]=0
    pan=ck.copy();threshold={'threshold':10,'strong_threshold':200}
    positive,envelope=duct_masks(ck,pan,threshold,threshold,1.)
    assert positive.sum()==1200
    assert envelope.sum()==1600
    assert not positive[30,30] and envelope[30,30]


def test_degenerate_marker_cannot_become_a_positive_classifier():
    result=threshold_proposal(np.zeros(3000))
    assert not result['supported']


def image_review_fixture(tmp_path):
    folder=tmp_path/'structural'/'78';folder.mkdir(parents=True)
    files={}
    for name in ['mask','review','companion']:
        path=folder/name;path.write_text(name);files[name]=path
    mask={'path':str(files['mask']),'sha256':sha256(files['mask']),'validated':False,'assessed_empty':False}
    evidence={'path':str(files['review']),'sha256':sha256(files['review']),
              'companion_path':str(files['companion']),'companion_sha256':sha256(files['companion'])}
    (folder/'proposal.json').write_text(json.dumps({'masks':{'Duct':mask},'reviews':{'Duct':evidence}}))
    registry=tmp_path/'structural/78.json';registry.write_text(json.dumps({'masks':{'Duct':mask}}))
    decision={'case_id':'78','structure':'Duct','mask_sha256':mask['sha256'],
              'review_sha256':evidence['sha256'],'companion_sha256':evidence['companion_sha256'],
              'decision':'accepted','rationale':'Known synthetic source evidence'}
    review=tmp_path/'review.json';review.write_text(json.dumps({'reviewer':'test','method':'source images','decisions':[decision]}))
    return files,registry,review


@pytest.mark.parametrize('changed',['mask','review','companion'])
def test_current_image_review_rejects_changed_evidence(tmp_path,changed):
    files,registry,review=image_review_fixture(tmp_path)
    files[changed].write_text('Changed after visual review')
    with pytest.raises(ValueError,match='changed'):accept(tmp_path,review)
    assert not json.loads(registry.read_text())['masks']['Duct']['validated']


def test_current_image_review_needs_companion_marker_and_preserves_empty_uncertainty(tmp_path):
    _,registry,review=image_review_fixture(tmp_path)
    document=json.loads(review.read_text());document['decisions'][0].pop('companion_sha256')
    review.write_text(json.dumps(document))
    with pytest.raises(ValueError,match='Companion-marker'):accept(tmp_path,review)
    assert not json.loads(registry.read_text())['masks']['Duct']['validated']
    _,registry,review=image_review_fixture(tmp_path/'complete')
    accept(tmp_path/'complete',review)
    result=json.loads(registry.read_text())['masks']['Duct']
    assert result['validated'] and not result['assessed_empty']
