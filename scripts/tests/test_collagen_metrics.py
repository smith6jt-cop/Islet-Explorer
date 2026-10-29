"""Union area and independent structural distance semantics."""
from pathlib import Path
import json
import subprocess
import sys
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import mapping
from shapely.geometry import box
from shapely.geometry import Polygon, MultiPolygon

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'uniform'/'typing'))
import structure_metrics as M

def test_overlaps_and_clipping_count_area_once():
    region=box(0,0,10,10)
    m=M.area_metrics(region,[box(-5,0,7,10),box(3,0,15,10)])
    assert m['area_density']==1
    assert m['count']==2
    assert np.isclose(m['area_um2'],100*M.CT.PIXEL_UM**2)

def test_empty_mask_and_region():
    m=M.area_metrics(box(0,0,10,10),[])
    assert m['count']==0 and m['area_density']==0 and m['area_um2']==0
    assert np.isnan(M.distances([[1,2]],[])[0])
    assert np.isnan(M.area_metrics(box(0,0,0,0),[])['area_density'])

def test_centroid_to_area_distance_units():
    xy=np.array([[5,5],[12,5]])*M.CT.PIXEL_UM
    d=M.distances(xy,[box(0,0,10,10)])
    assert np.allclose(d,[0,2*M.CT.PIXEL_UM])

def test_clipping_discards_line_contacts_and_keeps_polygon_parts():
    import collagen_structures as C
    tissue=box(0,0,10,10).union(box(20,0,30,10))
    result=C.clip_polygons([box(-1,0,20,5)],tissue,1.0)
    assert len(result)==1 and result[0].geom_type=='Polygon'
    assert result[0].area==50


@pytest.fixture
def measured_run(tmp_path,monkeypatch):
    # These tests exercise endpoint provenance. Full mask-input validation has separate focused tests.
    monkeypatch.setattr(M.C,'validate_current',lambda run,tags: {})
    def write(path, data):
        target=tmp_path/path;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_text(json.dumps(data))
    def feature(geometry, cls, name=None):
        return {'id':name or cls,'geometry':mapping(geometry),
                'properties':{'classification':{'name':cls},'name':name}}
    write('inputs/tissue/image.geojson',{'features':[feature(box(0,0,10,10),'Tissue')]})
    write('islets/image.geojson',{'features':[feature(box(2,2,4,4),'EndocrineSingle','Single_1')]})
    for cls in M.CLASSES:
        write(f'structures/{cls}/{cls}__image.geojson',{'features':[feature(box(0,0,5,10),cls)]})
    write('structures/collagen_iv.json',{'selected_method':'all-tissue'})
    (tmp_path/'cells').mkdir()
    pd.DataFrame({'Object ID':['a','b'],'Centroid X µm':np.array([3.,7.])*M.CT.PIXEL_UM,
                  'Centroid Y µm':np.array([5.,5.])*M.CT.PIXEL_UM}).to_parquet(tmp_path/'cells/image.parquet',index=False)
    # Independent cell identities are a sentinel: computing endpoints must not alter them.
    write('typing/identity/identity_sentinel.json',{'a':'Beta cell','b':'T cell'})
    before=(tmp_path/'typing/identity/identity_sentinel.json').read_bytes()
    result=M.one(tmp_path,{'image_name':'image','case_id':'1'})
    pd.DataFrame(result['records']).to_csv(tmp_path/'structure_region_metrics.csv',index=False)
    M.write_manifest(tmp_path,[result],M.file_hashes(tmp_path,['structures/collagen_iv.json']))
    assert (tmp_path/'typing/identity/identity_sentinel.json').read_bytes()==before
    return tmp_path,result


def test_metrics_manifest_covers_singles_cells_and_outputs(measured_run):
    run,result=measured_run
    manifest=M.validate_manifest(run,['image'])
    assert manifest['rows']==8
    assert {r['object_class'] for r in result['records']}=={'Tissue','EndocrineSingle'}
    distances=pd.read_parquet(run/'cell_structure_distances/image.parquet')
    assert distances.object_id.tolist()==['a','b']
    assert np.allclose(distances.distance_to_collagen_iv_um,[0,2*M.CT.PIXEL_UM])
    assert 'cells/image.parquet' in manifest['inputs']
    assert 'structure_region_metrics.csv' in manifest['outputs']


@pytest.mark.parametrize('relative',[
    'inputs/tissue/image.geojson','islets/image.geojson','cells/image.parquet',
    'structures/Collagen IV/Collagen IV__image.geojson','structures/Capillary/Capillary__image.geojson',
    'structures/collagen_iv.json','structure_region_metrics.csv',
    'cell_structure_distances/image.parquet','cell_structure_distances/image.json'])
def test_manifest_rejects_changed_input_or_output(measured_run,relative):
    run,_=measured_run
    with (run/relative).open('ab') as stream:stream.write(b'\n')
    with pytest.raises(ValueError,match='Stale structure metrics'):
        M.validate_manifest(run,['image'])


def test_manifest_rejects_code_change_missing_outputs_and_wrong_cohort(measured_run):
    run,_=measured_run
    with pytest.raises(ValueError,match='image set'):
        M.validate_manifest(run,['image','missing'])
    path=run/M.MANIFEST;manifest=json.loads(path.read_text())
    saved=manifest['metric_code_sha256'];manifest['metric_code_sha256']='old-code'
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError,match='code changed'):
        M.validate_manifest(run,['image'])
    manifest['metric_code_sha256']=saved
    manifest['outputs'].pop('cell_structure_distances/image.json')
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError,match='inventory is incomplete'):
        M.validate_manifest(run,['image'])


def test_manifest_not_published_when_input_changes_during_measurement(measured_run):
    run,result=measured_run
    (run/M.MANIFEST).unlink()
    with (run/'cells/image.parquet').open('ab') as stream:stream.write(b'\n')
    with pytest.raises(ValueError,match='Input changed while measuring'):
        M.write_manifest(run,[result],M.file_hashes(run,['structures/collagen_iv.json']))
    assert not (run/M.MANIFEST).exists()


def test_addition_structure_fingerprint_uses_actual_threshold_values(tmp_path):
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'uniform'))
    import redraw_structures as R
    (tmp_path/'structures').mkdir()
    path=tmp_path/'structures/thresholded.json'
    path.write_text(json.dumps({'reference_raw_threshold':{'Nerve':12.,'Capillary':24.},'sigma_px':1.}))
    before=R.threshold_inputs(tmp_path,{'cohort':'addition'})
    path.write_text(json.dumps({'reference_raw_threshold':{'Nerve':13.,'Capillary':24.},'sigma_px':1.}))
    assert R.threshold_inputs(tmp_path,{'cohort':'addition'})!=before
    assert R.threshold_inputs(tmp_path,{'cohort':'original'}) is None


def test_transferred_structure_contours_are_clipped_to_exact_tissue():
    import structure_thresholds as S
    tissue = box(0.25, 0.25, 9.75, 9.75)
    # A raster contour extends beyond a subpixel Tissue boundary; a touching line is not an object.
    result = S.clip_to_tissue([box(0, 0, 5, 5), box(9.75, 2, 12, 4)], tissue, 1.0)
    assert len(result) == 1
    assert result[0].equals(box(0.25, 0.25, 5, 5))
    assert tissue.covers(result[0])


def test_polygon_hole_does_not_erase_a_separate_tissue_island():
    import structure_thresholds as S
    outer=Polygon([(1,1),(19,1),(19,19),(1,19)],holes=[[(4,4),(16,4),(16,16),(4,16)]])
    island=box(7,7,13,13)
    first=np.zeros((22,22),np.uint8);second=np.zeros_like(first);separate=np.zeros_like(first)
    S.fill(first,mapping(MultiPolygon([island,outer])))
    S.fill(second,mapping(MultiPolygon([outer,island])))
    S.fill(separate,mapping(island));S.fill(separate,mapping(outer))
    assert np.array_equal(first,second) and np.array_equal(first,separate)
    assert first[10,10]==1 and first[5,5]==0 and first[2,2]==1
    assert first.sum()==(19*19-13*13)+7*7


def test_polygon_union_fill_keeps_subpixel_coordinates_and_clips_negative_bounds():
    import structure_thresholds as S
    import cv2
    geometry=Polygon([(-2.3,-1.2),(9.8,-1.2),(9.8,9.7),(-2.3,9.7)],
                     holes=[[(2.25,2.25),(6.75,2.25),(6.75,6.75),(2.25,6.75)]])
    expected=np.zeros((12,12),np.uint8)
    cv2.fillPoly(expected,[np.round(np.asarray(geometry.exterior.coords)*16).astype(np.int32)],1,shift=4)
    cv2.fillPoly(expected,[np.round(np.asarray(geometry.interiors[0].coords)*16).astype(np.int32)],0,shift=4)
    actual=np.zeros_like(expected);S.fill(actual,mapping(geometry))
    assert np.array_equal(actual,expected)


def test_classifier_refuses_wrong_embedded_tissue_before_inference(tmp_path):
    qupath=Path('/home/smith6jt/QuPath/bin/QuPath')
    if not qupath.is_file():pytest.skip('QuPath is unavailable')
    expected=tmp_path/'expected.geojson'
    expected.write_text(json.dumps({'type':'FeatureCollection','features':[{
        'type':'Feature','geometry':mapping(box(0,0,12,10)),
        'properties':{'objectType':'annotation','classification':{'name':'Tissue','color':[102,153,153]}}}]}))
    fixture=tmp_path/'check.groovy'
    fixture.write_text(r'''
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClass
import qupath.lib.roi.ROIs
import qupath.lib.regions.ImagePlane
def tissue = PathObjects.createAnnotationObject(ROIs.createRectangleROI(0, 0, 10, 10, ImagePlane.getDefaultPlane()),
    PathClass.fromString('Tissue'))
def local = new Binding([args: [args[1], 'Nerve', args[2]] as String[], getAnnotationObjects: { -> [tissue] }])
boolean rejected = false
try { new GroovyShell(local).evaluate(new File(args[0])) }
catch (IllegalStateException e) {
    assert e.getMessage().contains('Embedded QuPath Tissue differs')
    rejected = true
}
assert rejected
println('embedded Tissue mismatch rejected before inference')
''')
    script=Path(__file__).resolve().parents[1]/'uniform/qupath/structures_in_tissue.groovy'
    command=[str(qupath),'script','-a',str(script),'-a',str(tmp_path/'output'),'-a',str(expected),str(fixture)]
    result=subprocess.run(command,capture_output=True,text=True,timeout=60)
    assert result.returncode==0,result.stdout+result.stderr
    assert 'embedded Tissue mismatch rejected before inference' in result.stdout
