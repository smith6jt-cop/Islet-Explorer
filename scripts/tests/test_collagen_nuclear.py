"""Local nuclear-reference calibration, identity independence and mask freshness."""
import json
from pathlib import Path
import sys

import cv2
import numpy as np
import pytest
from shapely.geometry import Polygon, box, mapping

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'uniform'/'typing'))
import collagen_structures as C
import structure_metrics as M


def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value))


@pytest.fixture
def small_blocks(monkeypatch):
    monkeypatch.setattr(C,'BLOCK_UM',20*C.CT.PIXEL_UM)
    monkeypatch.setattr(C,'NUCLEAR_MIN_PIXELS',10)


def test_nucleus_union_preserves_holes_overlaps_skips_zero_area_and_ignores_identity(tmp_path):
    path=tmp_path/'cells.geojson'
    donut=Polygon([(2,2),(22,2),(22,22),(2,22)],holes=[[(8,8),(16,8),(16,16),(8,16)]])
    other=box(10,10,20,20)
    feats=[{'nucleusGeometry':mapping(g),'properties':{'classification':{'name':'arbitrary'}}}
           for g in [donut,other]]
    feats += [{'nucleusGeometry':{'type':'Polygon','coordinates':[[[0,0]]*5]}},{'geometry':mapping(box(0,0,25,25))}]
    write(path,{'features':feats});before=path.read_bytes()
    tissue=np.ones((28,28),bool);tissue[:4]=False
    ref,stats=C.nucleus_reference(path,tissue.shape,tissue)
    assert path.read_bytes()==before
    assert stats['nuclei_used']==2 and stats['degenerate_nucleus_skipped']==1 and stats['missing_nucleus']==1
    assert not ref[0,0] and not ref[9,9] and ref[14,14] and not ref[:4].any()
    for f in feats:f['properties']={'classification':{'name':'different label'},'measurements':{'irrelevant':999}}
    write(path,{'features':list(reversed(feats))})
    other_ref,_=C.nucleus_reference(path,tissue.shape,tissue)
    assert np.array_equal(ref,other_ref)


def test_lower_side_scale_is_unchanged_by_upper_tail_magnitude(small_blocks):
    values=np.tile(np.array([0,1,2,3,4,5,6,7,10,12],float),40).reshape(20,20)
    sm=np.expm1(values).astype(np.float32);ref=np.ones(sm.shape,bool)
    t1,s1=C.nuclear_threshold_map(sm,ref,ref)
    values[values>=10]=20
    t2,s2=C.nuclear_threshold_map(np.expm1(values).astype(np.float32),ref,ref)
    expected=np.expm1(4.5+2.5/C.LOWER_MAD_NORMALIZER)
    assert np.allclose(t1,expected,rtol=1e-6) and np.array_equal(t1,t2)
    assert s1['specificity_claim'] is False and s2['lower_scale_multiplier']==1


def test_missing_and_zero_scale_blocks_borrow_with_explicit_diagnostics(small_blocks):
    sm=np.tile(np.expm1(np.linspace(.1,2,20)),(20,3)).astype(np.float32)
    tissue=np.ones(sm.shape,bool);ref=np.zeros(sm.shape,bool);ref[:,:20]=True;ref[:,40:]=True
    sm[:,40:]=0
    threshold,stats=C.nuclear_threshold_map(sm,ref,tissue)
    assert threshold.shape==sm.shape and np.isfinite(threshold).all()
    assert stats['unsupported_tissue_blocks']==2 and stats['degenerate_scale_blocks']==1
    assert stats['nearest_reference_distance_blocks_p90_max_in_tissue'][1]==2
    assert 'reference_borrowing_beyond_one_block' in stats['qc_flags']
    assert [r['status'] for r in stats['block_diagnostics']]==[
        'own_reference','insufficient_reference','degenerate_lower_scale']
    assert np.allclose(threshold,threshold[0,0])


def test_empty_reference_or_flat_background_fails_instead_of_floor(small_blocks):
    tissue=np.ones((20,20),bool)
    for ref,sm in [(np.zeros_like(tissue),np.ones(tissue.shape)),(tissue,np.zeros(tissue.shape)),
                   (tissue,np.full(tissue.shape,np.nan))]:
        with pytest.raises(ValueError,match='no block'):
            C.nuclear_threshold_map(sm,ref,tissue)


def test_partial_edge_block_keeps_native_grid_position(small_blocks):
    # Last block is only half width. Its measured center remains x=50, not x=41.67.
    signal=np.tile(np.linspace(.1,1,20)[:,None],(1,50))
    signal[:,20:40]+=1;signal[:,40:]+=2
    sm=np.expm1(signal).astype(np.float32);ref=np.ones(sm.shape,bool)
    tmap,stats=C.nuclear_threshold_map(sm,ref,ref)
    last=np.expm1(stats['block_diagnostics'][-1]['threshold_log'])
    assert tmap[10,42]<.9*last
    assert tmap[10,49]>.95*last


@pytest.fixture
def collagen_run(tmp_path):
    image=tmp_path/'raw.tif';image.write_bytes(b'raw image fixture')
    write(tmp_path/'project_manifest.json',{'entries':[{'image_name':'image','entry_id':1,'image_path':str(image)}]})
    write(tmp_path/'qupath/data/1/server.json',{'metadata':{'channels':[{'name':'ColIV'},{'name':'DAPI'}]}})
    write(tmp_path/'inputs/tissue/image.geojson',{'features':[{'geometry':mapping(box(0,0,10,10)),
          'properties':{'classification':{'name':'Tissue'}}}]})
    write(tmp_path/'cell_geojson/image.geojson',{'features':[{'nucleusGeometry':mapping(box(2,2,8,8))}]})
    write(tmp_path/'structures/Capillary/Capillary__image.geojson',{'features':[]})
    output=tmp_path/'structures/Collagen IV/Collagen IV__image.geojson'
    write(output,{'features':[]})
    record={'provenance':C.provenance(tmp_path,'image'),'output_sha256':C.sha256(output)}
    write(tmp_path/'structures/collagen_iv.json',{'selected_method':'all-tissue','method_version':C.METHOD_VERSION,
          'images':{'image':record}})
    return tmp_path,record


@pytest.mark.parametrize('relative',[
    'cell_geojson/image.geojson','inputs/tissue/image.geojson','qupath/data/1/server.json',
    'structures/Collagen IV/Collagen IV__image.geojson','raw.tif'])
def test_current_mask_rejects_changed_geometry_channel_source_or_output(collagen_run,relative):
    run,_=collagen_run
    C.validate_current(run,['image'])
    with (run/relative).open('ab') as stream:stream.write(b'\n')
    with pytest.raises(ValueError,match='Stale Collagen IV'):
        C.validate_current(run,['image'])


def test_current_mask_rejects_old_method_alternate_nuclei_and_changed_code(collagen_run,tmp_path):
    run,record=collagen_run
    alternate=tmp_path/'alternate';alternate.mkdir()
    (alternate/'image.geojson').write_bytes((run/'cell_geojson/image.geojson').read_bytes())
    record['provenance']=C.provenance(run,'image',cell_dir=alternate)
    with pytest.raises(ValueError,match='Stale Collagen IV inputs'):
        C.validate_record(run,'image',record)
    record['provenance']=C.provenance(run,'image')
    record['provenance']['code'][str(Path(C.__file__).resolve())]='old code'
    with pytest.raises(ValueError,match='Stale Collagen IV inputs'):
        C.validate_record(run,'image',record)
    manifest=run/'structures/collagen_iv.json'
    write(manifest,{'selected_method':'all-tissue','images':{'image':record}})
    with pytest.raises(ValueError,match='method'):
        C.validate_current(run,['image'])


def test_changed_loaded_dependency_cannot_be_stamped_as_current(collagen_run,monkeypatch):
    run,_=collagen_run
    snapshot=dict(C.LOADED_CODE_SHA256)
    snapshot[str(Path(C.ST.__file__).resolve())]='different code loaded earlier'
    monkeypatch.setattr(C,'LOADED_CODE_SHA256',snapshot)
    with pytest.raises(ValueError,match='changed after module import'):
        C.provenance(run,'image')


def test_metrics_cli_rejects_stale_collagen_before_writes(collagen_run,monkeypatch):
    run,_=collagen_run
    with (run/'cell_geojson/image.geojson').open('ab') as stream:stream.write(b'\n')
    monkeypatch.setattr(M.CT,'cohort_rows',lambda:[{'image_name':'image'}])
    monkeypatch.setattr(sys,'argv',['structure_metrics.py','--run',str(run),'--workers','1'])
    with pytest.raises(ValueError,match='Stale Collagen IV'):
        M.main()
    assert not (run/'structure_region_metrics.csv').exists()
    assert not (run/'cell_structure_distances').exists()


@pytest.mark.parametrize('change_during_render',[False,True])
def test_one_publishes_fresh_checkpoint_without_changing_source_cells(collagen_run,monkeypatch,change_during_render):
    run,_=collagen_run
    monkeypatch.setattr(C,'NUCLEAR_MIN_PIXELS',5)
    arr=np.tile(np.arange(20,dtype=np.float32),(20,1))
    monkeypatch.setattr(C,'load',lambda run,tag,marker: arr.copy())
    source=run/'cell_geojson/image.geojson';before=source.read_bytes()
    sentinel=run/'typing/identity/sentinel.json';write(sentinel,{'cell':'Beta cell'});identity=sentinel.read_bytes()
    renders=run/'renders'
    def render(tag,col,dapi,tissue,binary,capillary,previous,dest):
        dest.mkdir()
        cv2.imwrite(str(dest/f'{tag}_whole.png'),np.zeros((20,20,3),np.uint8))
        write(dest/f'{tag}_tiles.json',[])
        if change_during_render:
            source.write_bytes(before+b'\n')
    monkeypatch.setattr(C,'render',render)
    output=run/'structures/Collagen IV/Collagen IV__image.geojson';old_output=output.read_bytes()
    if change_during_render:
        with pytest.raises(ValueError,match='inputs changed during processing'):
            C.one(run,'image',renders)
        assert output.read_bytes()==old_output and not C.checkpoint_path(run,'image').exists()
    else:
        result=C.one(run,'image',renders)
        C.validate_record(run,'image',result,check_renders=True)
        assert source.read_bytes()==before and result['specificity_claim'] is False
        assert C.checkpoint_path(run,'image').is_file()
    assert sentinel.read_bytes()==identity
