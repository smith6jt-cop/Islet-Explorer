"""Explicit decision partition, lossless artifacts, atomic readiness and native round-trip tests."""
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

UNIFORM=Path(__file__).resolve().parents[1]/'uniform'
sys.path.insert(0,str(UNIFORM))
import new_detection_filter as F


def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value))


@pytest.fixture
def raw(tmp_path):
    source=tmp_path/'raw';tag='image';area=source/'augmentation'/tag
    ids=['old','keep','reject','defer'];new=ids[1:]
    image=tmp_path/'raw_image.tif';image.write_bytes(b'immutable pixels')
    write(source/'project_manifest.json',{'entries':[{'image_name':tag,'entry_id':1,'image_path':str(image),'cohort':'original'}]})
    write(source/'qupath/project.qpproj',{'images':[{'imageName':tag,'entryID':1}]})
    write(source/'qupath/data/1/server.json',{'metadata':{'channels':[{'name':'DAPI'}]}})
    (source/'qupath/data/1/data.qpdata').write_bytes(b'unfiltered qpdata')
    frame=pd.DataFrame({'Object ID':ids,'Cell: DAPI: Mean':[np.nan,-0.,2.,3.],
                        'Cell: Area µm^2':[10.,11.,12.,13.],'Parent':['Tissue']*4})
    for folder in ('cells','final_cells','context'):
        (source/folder).mkdir();frame.to_parquet(source/folder/f'{tag}.parquet',index=False)
    features=[{'type':'Feature','id':oid,'geometry':{'type':'Polygon','coordinates':[[[i,0],[i+1,0],[i+1,1],[i,0]]]},
               'nucleusGeometry':{'type':'Point','coordinates':[i+.5,.5]},'properties':{'name':oid}} for i,oid in enumerate(ids)]
    write(source/'cell_geojson/image.geojson',{'type':'FeatureCollection','features':features})
    write(source/'islets/image.geojson',{'features':[{'id':'endo','properties':{'classification':{'name':'Islet'},'measurements':{'Whole cells':1}}}]})
    write(source/'inputs/tissue/image.geojson',{'features':[]})
    (source/'annotations').mkdir();pd.DataFrame({'Object ID':['tissue'],'Num Detections':[4],'Whole cells':[1]}).to_csv(source/'annotations/image.tsv',sep='\t',index=False)
    write(area/'original_inventory.json',{'cells':{'old':'original hash'},'annotations':{}})
    write(area/'new_object_ids.json',new)
    write(area/'verification.json',{'original_cells':1,'new_cells':3,'total_cells':4})
    mapping={'augment_detect':['inputs/tissue/image.geojson','augmentation/image/original_inventory.json','augmentation/image/new_object_ids.json'],
             'augment_C_cells':['cells/image.parquet'],'augment_C_geojson':['cell_geojson/image.geojson'],
             'augment_J_cells':['final_cells/image.parquet'],'augment_J_context':['context/image.parquet'],
             'augment_J_annotations':['annotations/image.tsv'],'augment_J_islets':['islets/image.geojson'],
             'augment_verify':['augmentation/image/verification.json']}
    state={'schema':'uniform.tissue_augmentation/1','inputs':{'outline_sha256':F.sha(source/'inputs/tissue/image.geojson')},
           'last_qpdata_sha256':F.sha(source/'qupath/data/1/data.qpdata'),
           'steps':{step:{'fingerprint':'fixture','outputs':{str(source/name):F.sha(source/name) for name in mapping.get(step,[])}} for step in F.REQUIRED_STEPS}}
    write(area/'state.json',state)
    evidence=tmp_path/'review.txt';evidence.write_text('synthetic reviewed decisions')
    decision={'schema':F.DECISIONS_SCHEMA,'rule':{'id':'fixture','version':'1','parameters':{},'decision_basis':'explicit fixture'},
              'images':{tag:{'source_seal':F.source_seal(source,tag),'review':{'status':'reviewed','reviewer':'test',
                         'evidence':[{'path':str(evidence),'sha256':F.sha(evidence)}]},
                         'decisions':[{'object_id':i,'decision':d,'reason':'fixture'} for i,d in zip(new,['accept','reject','defer'])]}}}
    decision_path=tmp_path/'decisions.json';write(decision_path,decision)
    return source,tmp_path/'qc',decision_path


def fake_native(stage,tag,original,raw_new,excluded):
    area=stage/'filter'/tag
    count=len(original|raw_new)-len(excluded)
    result={'ok':True,'retained_full_signatures_match':True,'annotation_full_signatures_match':True,
            'excluded_cells':len(excluded),'retained_cells':count,'annotation_detection_counts':{'tissue':count},
            'annotation_export_detection_counts':{'tissue':count}}
    write(area/'native_before.json',{'fixture':'full native behavior independently exercised below'})
    write(area/'native_apply.json',result);write(area/'native_verify.json',result)
    if excluded:
        qp=F.AU.qpdata(stage,tag);qp.write_bytes(qp.read_bytes()+b'filtered')
        p=stage/'annotations'/f'{tag}.tsv';df=pd.read_csv(p,sep='\t');df['Num Detections']=count;df.to_csv(p,sep='\t',index=False)
    return result


def test_complete_partition_rejects_original_unknown_duplicate_and_pending(raw):
    source,out,path=raw;d=F.read_json(path);seal=F.source_seal(source,'image')
    _,a,r,u=F.validate_decisions(d,'image',seal)
    assert (a,r,u)==(['keep'],['reject'],['defer'])
    for replacement in ('old','unknown','keep'):
        bad=json.loads(json.dumps(d));bad['images']['image']['decisions'][1]['object_id']=replacement
        with pytest.raises(ValueError):F.validate_decisions(bad,'image',seal)
    bad=json.loads(json.dumps(d));bad['images']['image']['review']['status']='pending'
    with pytest.raises(ValueError,match='reviewed scientific'):F.validate_decisions(bad,'image',seal)


def test_filter_preserves_source_and_retained_float_bits_and_defer(raw,monkeypatch):
    source,out,path=raw;before=F.source_seal(source,'image');monkeypatch.setattr(F,'run_native',fake_native)
    result=F.process(source,out,path,'image');assert F.source_seal(source,'image')==before
    assert result['accepted_new_ids']==['keep'] and result['rejected_new_ids']==['reject'] and result['deferred_new_ids']==['defer']
    assert result['counts']['total_cells']==2 and result['counts']['excluded_new_cells']==2
    assert list(pd.read_parquet(out/'cells/image.parquet')['Object ID'])==['old','keep']
    assert np.signbit(pd.read_parquet(out/'cells/image.parquet')['Cell: DAPI: Mean'].iloc[1])
    assert not os.path.samefile(source/'qupath/data/1/data.qpdata',out/'qupath/data/1/data.qpdata')
    assert F.accepted_inputs(out,'image')['total_cells']==2
    # Unrelated image additions to the same document do not invalidate this image's reviewed decision projection.
    d=F.read_json(path);d['images']['other']={};write(path,d)
    assert F.accepted_inputs(out,'image')['total_cells']==2
    assert F.process(source,out,path,'image')==result


def test_noop_is_byte_identical_and_readiness_rejects_changed_output(raw,monkeypatch):
    source,out,path=raw;d=F.read_json(path)
    for item in d['images']['image']['decisions']:item['decision']='accept'
    write(path,d);monkeypatch.setattr(F,'run_native',fake_native);result=F.process(source,out,path,'image')
    assert result['checks']['noop_byte_identical']
    for key in ('cells','final_cells','context','geometry','annotations','islets','tissue','qpdata'):
        assert F.sha(Path(result['source_seal']['files'][key]['path']))==result['outputs'][key]['sha256']
    (out/'cell_geojson/image.geojson').write_text('{}')
    with pytest.raises(ValueError,match='filtered output changed'):F.accepted_inputs(out,'image')


def test_annotation_export_counts_follow_spatial_measurements_not_ancestry(raw,monkeypatch):
    source,out,path=raw
    def spatial_native(*args):
        report=fake_native(*args)
        report['annotation_detection_counts']['tissue']=0
        return report
    monkeypatch.setattr(F,'run_native',spatial_native)
    F.process(source,out,path,'image')
    report={'annotation_detection_counts':{'tissue':0},'annotation_export_detection_counts':{'tissue':1}}
    with pytest.raises(ValueError,match='dynamic counts differ'):
        F.verify_annotations(source,out,'image',{'reject','defer'},report)


def test_failed_native_attempt_never_releases_and_retry_uses_raw_source(raw,monkeypatch):
    source,out,path=raw;before=F.source_seal(source,'image')
    def fail(*args):raise RuntimeError('interrupted native write')
    monkeypatch.setattr(F,'run_native',fail)
    with pytest.raises(RuntimeError):F.process(source,out,path,'image')
    assert F.accepted_inputs(out,'image') is None and F.source_seal(source,'image')==before
    assert len(list((out/'filter/.staging').iterdir()))==1
    monkeypatch.setattr(F,'run_native',fake_native);F.process(source,out,path,'image')
    assert F.accepted_inputs(out,'image')['total_cells']==2
    assert len(list((out/'filter/.staging').iterdir()))==2


def test_changed_review_source_or_retained_value_is_rejected(raw,monkeypatch):
    source,out,path=raw;monkeypatch.setattr(F,'run_native',fake_native);F.process(source,out,path,'image')
    d=F.read_json(path);d['images']['image']['decisions'][0]['decision']='defer';write(path,d)
    with pytest.raises(ValueError,match='decisions changed'):F.accepted_inputs(out,'image')
    # The source seal rejects mutations even when no completion is being consumed.
    (source/'qupath/data/1/data.qpdata').write_bytes(b'changed')
    with pytest.raises(ValueError,match='qpdata changed'):F.source_seal(source,'image')


def test_source_change_during_work_withholds_completion(raw,monkeypatch):
    source,out,path=raw
    def change(stage,*args):
        result=fake_native(stage,*args)
        (source/'qupath/data/1/data.qpdata').write_bytes(b'external mutation')
        return result
    monkeypatch.setattr(F,'run_native',change)
    with pytest.raises(ValueError,match='qpdata changed'):F.process(source,out,path,'image')
    assert F.accepted_inputs(out,'image') is None


def test_refuses_unowned_output_and_changed_pending_decisions(raw,monkeypatch):
    source,out,path=raw;target=out/'cells/image.parquet';target.parent.mkdir(parents=True);target.write_bytes(b'user data')
    monkeypatch.setattr(F,'run_native',fake_native)
    with pytest.raises(ValueError,match='unowned existing'):F.process(source,out,path,'image')
    assert target.read_bytes()==b'user data'
    target.unlink()
    def fail(*args):raise RuntimeError('failed')
    monkeypatch.setattr(F,'run_native',fail)
    with pytest.raises(RuntimeError):F.process(source,out,path,'image')
    d=F.read_json(path);d['images']['image']['decisions'][0]['decision']='defer';write(path,d)
    with pytest.raises(ValueError,match='Pending filter attempt'):F.process(source,out,path,'image')


def test_structure_reuse_pins_inputs_without_forging_redraw_record(raw,monkeypatch):
    source,out,path=raw;monkeypatch.setattr(F,'run_native',fake_native);F.process(source,out,path,'image')
    inputs={'fixture':'unchanged pixel inference'}
    monkeypatch.setattr(F,'_structure_inputs',lambda *args:dict(inputs))
    classes={}
    for cls in ('Nerve','Capillary','Lymphatic'):
        target=source/'structures'/cls/f'{cls}__image.geojson';write(target,{'features':[]})
        classes[cls]={'sha256':F.sha(target),'edge_touching_retained':True}
    write(source/'structures/redraw_records/image.json',{'inputs':dict(inputs),'classes':classes})
    certificate=F.certify_structures(out,'image')
    assert F.certify_structures(out,'image')==certificate
    assert not (out/'structures/redraw_records/image.json').exists()
    assert F.structure_reuse_current(out,'image') is not None
    inputs['fixture']='changed classifier'
    with pytest.raises(ValueError,match='classifier/threshold/code'):F.structure_reuse_current(out,'image')


def test_native_deletion_and_reload_preserve_parents_measurements_and_annotations(tmp_path):
    if not F.QUPATH.is_file():pytest.skip('QuPath unavailable')
    fixture=tmp_path/'fixture.groovy';fixture.write_text(r'''
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClass
import qupath.lib.objects.hierarchy.PathObjectHierarchy
import qupath.lib.roi.ROIs
import qupath.lib.regions.ImagePlane
import qupath.lib.images.ImageData
import qupath.lib.images.servers.WrappedBufferedImageServer
import qupath.lib.io.PathIO
import java.awt.image.BufferedImage
import com.google.gson.Gson
def helper=new GroovyShell().evaluate(new File(args[0]))
new GroovyShell().parse(new File(args[1]))
def rect={x,y,w,h -> ROIs.createRectangleROI(x,y,w,h,ImagePlane.getDefaultPlane())}
def h=new PathObjectHierarchy()
def t=PathObjects.createAnnotationObject(rect(0,0,100,100),PathClass.fromString('Tissue'))
def islet=PathObjects.createAnnotationObject(rect(0,0,50,50),PathClass.fromString('Islet'))
islet.setName('Islet_1');islet.getMeasurementList().put('Whole cells',7)
def old=PathObjects.createCellObject(rect(1,1,5,5),rect(2,2,2,2),null,null)
def keep=PathObjects.createCellObject(rect(10,10,5,5),rect(11,11,2,2),null,null)
def removed=PathObjects.createCellObject(rect(20,20,5,5),rect(21,21,2,2),null,null)
def rootCell=PathObjects.createCellObject(rect(70,70,5,5),rect(71,71,2,2),null,null)
old.getMeasurementList().put('Cell: DAPI: Mean',12.25)
old.getMeasurementList().put('Signed distance to annotation with Islet µm',-3.0)
old.getMetadata().put('raw','preserve');old.setColor(0xFF00FF)
h.addObject(t);t.addChildObject(islet);islet.addChildObjects([old,keep,removed])
h.getRootObject().addChildObject(rootCell)
def ids={objects -> objects.collect{it.getID().toString()} as Set}
def originals=ids([old,rootCell]), rawNew=ids([keep,removed]), excluded=ids([removed])
boolean caught=false
try {helper.remove(h,originals,rawNew,originals)} catch(IllegalStateException e) {caught=true}
assert caught
def before=helper.remove(h,originals,rawNew,excluded)
assert h.getDetectionObjects().size()==3
assert helper.verify(h,before,excluded).annotation_detection_counts[islet.getID().toString()]==2
assert helper.verify(h,before,excluded).annotation_detection_counts[t.getID().toString()]==2
old.getMeasurementList().put('Signed distance to annotation with Islet µm',-4.0)
caught=false
try {helper.verify(h,before,excluded)} catch(IllegalStateException e){caught=true}
assert caught
old.getMeasurementList().put('Signed distance to annotation with Islet µm',-3.0)
t.addChildObject(old)
caught=false
try {helper.verify(h,before,excluded)} catch(IllegalStateException e){caught=true}
assert caught
islet.addChildObject(old)
def server=new WrappedBufferedImageServer('fixture',new BufferedImage(100,100,BufferedImage.TYPE_BYTE_GRAY))
def imageData=new ImageData(server,h,ImageData.ImageType.FLUORESCENCE)
def saved=new File(args[2]);PathIO.writeImageData(saved,imageData)
def reloaded=PathIO.readImageData(saved,server)
assert helper.verify(reloaded.getHierarchy(),before,excluded).retained_cells==3
assert reloaded.getHierarchy().getAnnotationObjects().find{it.getName()=='Islet_1'}.getMeasurementList().get('Whole cells')==7
// Exercise the actual apply and verify script using a fresh independent ImageData and a native save.
def reportDir=new File(args[3]);reportDir.mkdirs()
def cfg=[helper:args[0],original_ids:originals.toList(),raw_new_ids:[keep.getID().toString()],excluded_ids:[keep.getID().toString()],
 before_inventory:new File(reportDir,'before.json').getPath(),annotation_output:new File(reportDir,'ann.tsv').getPath(),
 mode:'apply',report:new File(reportDir,'apply.json').getPath()]
def config=new File(reportDir,'config.json');config.text=new Gson().toJson(cfg)
def destination=new File(reportDir,'after.qpdata')
def entry=new Expando(getImageName:{->'fixture'},saveImageData:{data->PathIO.writeImageData(destination,data)})
def local=new Binding([args:[config.getPath()] as String[],getCurrentHierarchy:{->reloaded.getHierarchy()},
 getCurrentImageData:{->reloaded},getProjectEntry:{->entry},saveAnnotationMeasurements:{String path->new File(path).text='export invoked'}])
new GroovyShell(local).evaluate(new File(args[1]))
assert new File(cfg.annotation_output).text=='export invoked'
def applied=new Gson().fromJson(new File(cfg.report).text,Map)
assert applied.annotation_detection_counts[t.getID().toString()]==1
assert applied.annotation_export_detection_counts[t.getID().toString()]==2
assert applied.annotation_export_detection_counts[islet.getID().toString()]==1
def finalData=PathIO.readImageData(destination,server)
cfg.mode='verify';cfg.report=new File(reportDir,'verify.json').getPath();config.text=new Gson().toJson(cfg)
local.setVariable('getCurrentHierarchy',{->finalData.getHierarchy()});local.setVariable('getCurrentImageData',{->finalData})
new GroovyShell(local).evaluate(new File(args[1]))
assert finalData.getHierarchy().getDetectionObjects().size()==2
def verified=new Gson().fromJson(new File(cfg.report).text,Map)
assert verified.ok
assert verified.annotation_detection_counts[t.getID().toString()]==1
assert verified.annotation_export_detection_counts[t.getID().toString()]==2
assert verified.annotation_export_detection_counts==applied.annotation_export_detection_counts
assert finalData.getHierarchy().getDetectionObjects().find{it.getID()==rootCell.getID()}.getParent().isRootObject()
println('native retained raw/geometry/parent/annotation signatures survive deletion and reload')
''')
    command=[str(F.QUPATH),'script','-a',str(F.HELPER),'-a',str(F.SCRIPT),'-a',str(tmp_path/'saved.qpdata'),'-a',str(tmp_path/'reports'),str(fixture)]
    result=subprocess.run(command,capture_output=True,text=True,timeout=90,
                          env={**os.environ,'JAVA_TOOL_OPTIONS':'-Xmx1g -XX:ActiveProcessorCount=2'})
    assert result.returncode==0,result.stdout+result.stderr
    assert 'signatures survive deletion and reload' in result.stdout
