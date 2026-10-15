#!/usr/bin/env python3
"""Disease-blinded structural pixel candidates and explicit image audit.

Mixture thresholds generate proposals, never certified masks. Each donor/target
requires a signed image-audit decision that binds the candidate and review
montage hashes. Positive-pixel polygons and filled duct envelopes are separate.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import gzip
import hashlib
import json
import os
from pathlib import Path
import time
import xml.etree.ElementTree as ET
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='1'
import cv2
import numpy as np
import pandas as pd
from PIL import Image,ImageDraw
from scipy.ndimage import binary_fill_holes
from shapely import box,coverage_union_all,get_parts,geometry,contains_xy
from shapely.geometry import mapping
from sklearn.mixture import GaussianMixture
import tifffile
from common import sha256,write_json

MARKERS=['ColIV','CK19','panCK']
CONTEXT=['ECAD','BCatenin']
INHERITED={'Nerve':'B3TUBB','Capillary':'CD31','Lymphatic':'PDPN'}
COMPANIONS={'ColIV':['ColIV','ECAD','panCK'], 'Duct':['CK19','panCK','ECAD','BCatenin'],
            'Nerve':['B3TUBB','PGP9_5','GAP43'], 'Capillary':['CD31','CD34','SMA'],
            'Lymphatic':['PDPN','CD31','HLADR']}
SEED=20260919


def source_info(image):
    with tifffile.TiffFile(image) as tif:
        root=ET.fromstring(tif.ome_metadata)
        pixels=root.find('.//{*}Pixels')
        names=[node.attrib.get('Name','') for node in pixels.findall('{*}Channel')]
        if len(set(names))!=len(names) or not set(MARKERS+list(INHERITED.values())).issubset(names):
            raise ValueError('Missing/duplicate structural channels')
        units={'µm':1.,'um':1.,'micrometer':1.,'mm':1000.,'nm':.001}
        try:
            pixel=float(pixels.attrib['PhysicalSizeX'])*units[pixels.attrib.get('PhysicalSizeXUnit','µm')]
            pixel_y=float(pixels.attrib['PhysicalSizeY'])*units[pixels.attrib.get('PhysicalSizeYUnit','µm')]
        except KeyError:raise ValueError('Unsupported physical unit')
        if abs(pixel-pixel_y)>1e-8:raise ValueError('Non-square pixels require explicit calibration')
        return names,pixel,tuple(tif.series[0].shape[-2:])


def read_channel(image,names,marker):
    with tifffile.TiffFile(image) as tif:
        return tif.series[0].pages[names.index(marker)].asarray(maxworkers=2)


def threshold_proposal(values,seed=SEED,posterior=.95):
    """Fit an intensity-only proposal. Overlapping/unimodal signals fail closed."""
    intensity=np.asarray(values,float)
    # Exact zero is a point mass (camera/background censoring), not a Gaussian
    # positive population. Fit the continuous nonzero signal distribution.
    v=np.log1p(intensity[np.isfinite(intensity)&(intensity>0)])
    if len(v)<1000 or np.std(v)<1e-8:return dict(supported=False,reason='insufficient_intensity_information')
    model=GaussianMixture(2,n_init=3,random_state=seed,reg_covar=1e-4).fit(v[:,None])
    order=np.argsort(model.means_[:,0]);low,high=order
    separation=float((model.means_[high,0]-model.means_[low,0])/np.sqrt(model.covariances_[low,0,0]+model.covariances_[high,0,0]))
    grid=np.linspace(float(np.min(v)),float(np.max(v)),8192)
    probability=model.predict_proba(grid[:,None])[:,high]
    positive=(grid>model.means_[low,0])&(probability>=posterior)
    if not positive.any() or separation<.75:
        # An intensity quantile supplies an image-review proposal only; it does
        # not establish a positive/background split or certify any mask.
        return dict(supported=True,statistically_separated=False,reason='unresolved_positive_background_split_quantile_review_only',separation=separation,
            threshold=float(np.expm1(np.quantile(v,.995))),strong_threshold=float(np.expm1(np.quantile(v,.999))),
            interpretation='Conservative bright-signal proposal; mandatory image audit, no statistical positivity claim')
    threshold=float(np.expm1(grid[np.flatnonzero(positive)[0]]))
    strong=(grid>model.means_[high,0])&(probability>=.999)
    return dict(supported=True,statistically_separated=True,threshold=threshold,strong_threshold=float(np.expm1(grid[np.flatnonzero(strong)[0]])) if strong.any() else float(np.expm1(grid[-1])),
        separation=separation,log_means=model.means_[:,0].tolist(),log_variances=model.covariances_[:,0,0].tolist(),weights=model.weights_.tolist(),
        positive_posterior=posterior,interpretation='Unsupervised proposal requiring image audit; not validated by mixture separation')


def exact_pixel_geojson(mask,path,label):
    """Union half-open pixel squares, retaining exact positive area and holes."""
    features=[]
    for start in range(0,len(mask),128):
        tile=mask[start:start+128]
        change=np.diff(np.pad(tile.astype(np.int8),((0,0),(1,1))),axis=1)
        ys,left=np.nonzero(change==1);ye,right=np.nonzero(change==-1)
        if not np.array_equal(ys,ye):raise ValueError('Invalid binary run ordering')
        if not len(left):continue
        polygons=coverage_union_all(box(left,ys+start,right,ys+start+1))
        for poly in get_parts(polygons):
            colors={'ColIV':[190,80,220],'Duct':[255,210,50],'DuctEnvelope':[255,160,30]}
            features.append(dict(type='Feature',geometry=mapping(poly),properties=dict(objectType='annotation',name=f'{label}_{len(features)+1}',
                classification=dict(name=label,color=colors.get(label,[100,200,200])),units='source_pixels',positive_pixels=poly.area)))
        if start%512==0:print(label,'polygonized rows',start+len(tile),'of',len(mask),flush=True)
    payload=dict(type='FeatureCollection',features=features)
    with gzip.open(path,'wt') as stream:json.dump(payload,stream,separators=(',',':'),allow_nan=False)
    area=sum(feature['properties']['positive_pixels'] for feature in features)
    if not np.isclose(area,int(mask.sum()),rtol=0,atol=1e-6):raise ValueError('Positive mask polygon area differs from pixel count')
    return dict(path=str(path),sha256=sha256(path),positive_pixels=int(mask.sum()),features=len(features))


def minimum_objects(mask,min_pixels):
    n,labels,stats,_=cv2.connectedComponentsWithStats(mask.astype(np.uint8),connectivity=8)
    keep=stats[:,cv2.CC_STAT_AREA]>=min_pixels;keep[0]=False
    return keep[labels]


def duct_masks(ck,pan,ck_threshold,pan_threshold,pixel):
    positive=(ck>ck_threshold['threshold'])&((pan>pan_threshold['threshold'])|(ck>ck_threshold['strong_threshold']))
    # This cleanup defines a classifier object envelope. It cannot replace the
    # measured positive-pixel area, which retains only originally positive pixels.
    radius=max(1,round(2/pixel));kernel=cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(2*radius+1,2*radius+1))
    envelope=cv2.morphologyEx(positive.astype(np.uint8),cv2.MORPH_CLOSE,kernel)>0
    envelope=binary_fill_holes(envelope)
    envelope=minimum_objects(envelope,max(1,round(20/pixel**2)))
    return positive&envelope,envelope


def locations(image,tissue,seed,count=4,size=256):
    rng=np.random.default_rng(seed);shape=image.shape
    ys=rng.integers(size//2,shape[0]-size//2,20000);xs=rng.integers(size//2,shape[1]-size//2,20000)
    keep=contains_xy(tissue,xs,ys);ys,xs=ys[keep],xs[keep]
    intensity=image[ys,xs];order=np.argsort(intensity)
    picked=[]
    for fraction in [.60,.90,.99,.999]:
        location=int(order[min(len(order)-1,int(fraction*(len(order)-1)))])
        picked.append((int(xs[location]),int(ys[location])))
    return picked


def display_channel(data,lo=None,hi=None):
    x=np.maximum(data.astype(float),0);lo=float(np.quantile(x,.02)) if lo is None else lo;hi=float(np.quantile(x,.995)) if hi is None else hi
    return np.uint8(255*np.clip((x-lo)/max(hi-lo,1e-6),0,1))


def render_review(path,label,marker,channel,mask,tissue,seed,context=None,source_polygons=None):
    size=256;locations_=locations(channel,tissue,seed);canvas=Image.new('RGB',(size*3,size*4+72),'white');draw=ImageDraw.Draw(canvas)
    draw.text((12,8),f'{label} | raw {marker} | proposed overlay | epithelial / companion context',fill='black')
    draw.text((12,30),'Disease-blinded donor; low, intermediate, high and bright sampled signals. Cyan = candidate.',fill='black')
    records=[]
    for row,(x,y) in enumerate(locations_):
        x0=x-size//2;y0=y-size//2;sl=np.s_[y0:y0+size,x0:x0+size]
        raw=display_channel(channel[sl]);rgb=np.stack([raw]*3,axis=-1);overlay=rgb.copy()
        if mask is not None:
            selected=mask[sl]
            edges=cv2.morphologyEx(selected.astype(np.uint8),cv2.MORPH_GRADIENT,np.ones((3,3),np.uint8))>0
            overlay[edges]=[0,255,255]
        if source_polygons is not None:
            p=Image.fromarray(overlay);pen=ImageDraw.Draw(p)
            for poly in source_polygons:
                if not poly.intersects(box(x0,y0,x0+size,y0+size)):continue
                for part in get_parts(poly):
                    if part.geom_type!='Polygon':continue
                    pen.line([(xx-x0,yy-y0) for xx,yy in part.exterior.coords],fill=(0,255,255),width=1)
            overlay=np.asarray(p)
        ctx=np.zeros_like(rgb)
        if context:
            for index,data in enumerate(context[:3]):ctx[...,index]=display_channel(data[sl])
        else:ctx=rgb
        for col,img in enumerate([rgb,overlay,ctx]):canvas.paste(Image.fromarray(img),(col*size,72+row*size))
        draw.text((3,72+row*size+3),f'x={x0}, y={y0}',fill=(255,80,80))
        records.append(dict(x=x0,y=y0,width=size,height=size))
    canvas.save(path)
    return dict(path=str(path),sha256=sha256(path),crops=records,marker=marker)


def propose_one(args):
    root,cid,sample_count=args;root=Path(root);start=time.monotonic();cv2.setNumThreads(2)
    cohort=pd.read_csv(root/'cohort.csv',dtype=str);row=cohort[cohort.case_id.eq(cid)].iloc[0]
    image=Path(row.image);names,ome_pixel,shape=source_info(image)
    calibration=json.loads((root/'boundaries/sources'/cid/'source.json').read_text())
    pixel=float(calibration['pixel_size_um'])
    source=root/'boundaries/sources'/cid/'valid_tissue.geojson'
    tissue_doc=json.loads(source.read_text())
    geoms=[geometry.shape(f['geometry']) for f in tissue_doc['features']] if tissue_doc.get('type')=='FeatureCollection' else [geometry.shape(tissue_doc.get('geometry',tissue_doc))]
    from shapely import union_all
    tissue=union_all(geoms)
    exclusion_path=root/'typing/eligibility/technical_exclusions.json'
    if exclusion_path.exists():
        exclusions=json.loads(exclusion_path.read_text()).get('regions',[])
        excluded=[geometry.shape(region['geometry']) for region in exclusions if str(region['case_id'])==cid]
        if excluded:tissue=tissue.difference(union_all(excluded))
    out=root/'structural'/cid;out.mkdir(parents=True,exist_ok=True)
    if (out/'proposal.json').exists():return dict(case_id=cid,status='already_proposed')
    rng=np.random.default_rng(SEED+int(cid));ys=rng.integers(0,shape[0],sample_count*3);xs=rng.integers(0,shape[1],sample_count*3)
    keep=contains_xy(tissue,xs+.5,ys+.5);ys,xs=ys[keep][:sample_count],xs[keep][:sample_count]
    arrays={m:read_channel(image,names,m) for m in MARKERS}
    posterior={'ColIV':.95,'CK19':.90,'panCK':.50}
    thresholds={m:threshold_proposal(arrays[m][ys,xs],SEED+int(cid),posterior[m]) for m in MARKERS}
    pd.DataFrame({m:arrays[m][ys,xs] for m in MARKERS}).to_parquet(out/'intensity_samples.parquet',index=False)
    reviews={};masks={}
    if thresholds['ColIV']['supported']:
        mask=arrays['ColIV']>thresholds['ColIV']['threshold']
        reviews['ColIV']=render_review(out/'review_ColIV.png','ColIV','ColIV',arrays['ColIV'],mask,tissue,SEED+int(cid),[arrays['CK19'],arrays['panCK'],arrays['ColIV']])
        info=exact_pixel_geojson(mask,out/'ColIV_positive.geojson.gz','ColIV');info['positive_area_um2']=int(mask.sum())*pixel**2
        masks['ColIV']=dict(info,validated=False,assessed_empty=False,definition='thresholded positive pixels; no filled envelope')
        del mask
    if all(thresholds[m]['supported'] for m in ['CK19','panCK']):
        positive,envelope=duct_masks(arrays['CK19'],arrays['panCK'],thresholds['CK19'],thresholds['panCK'],pixel)
        info=exact_pixel_geojson(positive,out/'Duct_positive.geojson.gz','Duct');envelope_info=exact_pixel_geojson(envelope,out/'Duct_envelope.geojson.gz','DuctEnvelope')
        masks['Duct']=dict(info,validated=False,assessed_empty=False,envelope=envelope_info,definition='CK19-positive pixels with panCK support or stringent CK19; morphology only defines duct envelope',positive_area_um2=int(positive.sum())*pixel**2)
        masks['DuctEnvelope']=dict(envelope_info,validated=False,assessed_empty=False,definition='filled 2-um closed duct envelope for distance only; never structural positive area')
        reviews['Duct']=render_review(out/'review_Duct.png','Duct','CK19',arrays['CK19'],envelope,tissue,SEED+int(cid)+1,[arrays['CK19'],arrays['panCK'],arrays['ColIV']])
        reviews['DuctEnvelope']=reviews['Duct']
        del positive,envelope
    del arrays
    with gzip.open(row.previous_geometry,'rt') as stream:legacy=json.load(stream)
    for label,marker in INHERITED.items():
        features=[f for f in legacy['features'] if f.get('properties',{}).get('classification',{}).get('name')==label]
        polygons=[geometry.shape(f['geometry']) for f in features]
        channel=read_channel(image,names,marker)
        reviews[label]=render_review(out/f'review_{label}.png',label,marker,channel,None,tissue,SEED+int(cid)+len(label),source_polygons=polygons)
        path=out/f'{label}_inherited.geojson.gz'
        with gzip.open(path,'wt') as stream:json.dump(dict(type='FeatureCollection',features=features),stream,separators=(',',':'))
        masks[label]=dict(path=str(path),sha256=sha256(path),validated=False,assessed_empty=False,features=len(features),source=str(row.previous_geometry),source_sha256=sha256(row.previous_geometry),definition='inherited classifier outlines; positive-pixel fidelity requires image audit')
        del channel
    image_stat=image.stat()
    record=dict(case_id=cid,source_image=str(image),image_bytes=image_stat.st_size,image_mtime_ns=image_stat.st_mtime_ns,
        channel_names=names,pixel_um=pixel,analysis_pixel_size_um=pixel,acquired_ome_pixel_size_um=ome_pixel,
        calibration_source=str(root/'boundaries/sources'/cid/'source.json'),shape_yx=shape,threshold_proposals=thresholds,masks=masks,reviews=reviews,
        tissue_source=str(source),tissue_sha256=sha256(source),technical_exclusions_sha256=sha256(exclusion_path) if exclusion_path.exists() else None,
        disease_labels_used=False,seconds=time.monotonic()-start,
        validation='pending image-reference review; proposal generation is not acceptance')
    write_json(out/'proposal.json',record)
    write_json(root/'structural'/f'{cid}.json',dict(case_id=cid,masks=masks,proposal_sha256=sha256(out/'proposal.json')))
    print(cid,'structural proposals complete',flush=True)
    return dict(case_id=cid,seconds=record['seconds'],supported_proposals=list(masks))


def accept(root,review_path):
    review=json.loads(Path(review_path).read_text())
    if not review.get('reviewer') or not review.get('method'):raise ValueError('Image reviewer and method required')
    pending={}
    for decision in review['decisions']:
        cid=str(decision['case_id']);label=decision['structure'];path=root/'structural'/cid/'proposal.json'
        proposal=json.loads(path.read_text());mask=proposal['masks'][label];montage=proposal['reviews'][label]
        if decision['mask_sha256']!=mask['sha256'] or decision['review_sha256']!=montage['sha256']:raise ValueError('Stale image-review evidence')
        if decision.get('companion_sha256')!=montage.get('companion_sha256') or not montage.get('companion_sha256'):
            raise ValueError('Companion-marker image evidence is missing or stale')
        if decision.get('object_sha256')!=montage.get('object_sha256'):
            raise ValueError('Object-targeted image evidence is stale')
        for path_key,hash_key in [('path','sha256'),('companion_path','companion_sha256'),('object_path','object_sha256')]:
            if montage.get(hash_key) and sha256(montage[path_key])!=montage[hash_key]:
                raise ValueError('Rendered image-review evidence changed')
        for evidence in decision.get('supporting_evidence',[]):
            if sha256(evidence['path'])!=evidence['sha256']:
                raise ValueError('Supporting source-image evidence changed')
        if sha256(mask['path'])!=mask['sha256']:raise ValueError('Structural mask changed after image review')
        if decision['decision'] not in {'accepted','rejected','unavailable'}:raise ValueError('Invalid image decision')
        if not decision.get('rationale'):raise ValueError('Review rationale required')
        record_path=root/'structural'/f'{cid}.json'
        record=pending.setdefault(record_path,json.loads(record_path.read_text()))
        record['masks'][label].update(validated=decision['decision']=='accepted',assessed_empty=decision.get('assessed_empty',False),
            review=decision,reviewer=review['reviewer'],review_registry_sha256=sha256(review_path))
    # Reject a stale decision before applying any part of a multi-donor review.
    for record_path,record in pending.items():write_json(record_path,record)


def normalize_existing(root,cid):
    """Complete QuPath properties and bind geometry to native analysis calibration.

    The sub-0.02% OME calibration difference cannot change the rounded duct
    morphology radii/object cutoffs. Fail if it could; do not relabel geometry.
    """
    path=root/'structural'/cid/'proposal.json';record=json.loads(path.read_text())
    registry_path=root/'structural'/f'{cid}.json';registry=json.loads(registry_path.read_text())
    if any(m.get('validated') for m in registry['masks'].values()):raise ValueError('Cannot normalize accepted evidence')
    calibration_path=root/'boundaries/sources'/cid/'source.json'
    pixel=float(json.loads(calibration_path.read_text())['pixel_size_um']);old=record['pixel_um']
    if round(2/pixel)!=round(2/old) or round(20/pixel**2)!=round(20/old**2):
        raise ValueError('Calibration changes duct morphology; regenerate masks')
    colors={'ColIV':[190,80,220],'Duct':[255,210,50],'DuctEnvelope':[255,160,30]}
    for label,color in colors.items():
        if label not in record['masks']:continue
        info=record['masks'][label];target=Path(info['path'])
        with gzip.open(target,'rt') as stream:doc=json.load(stream)
        changed=False
        for index,feature in enumerate(doc['features']):
            prop=feature.setdefault('properties',{});classification=prop.setdefault('classification',{})
            for key,value in [('name',label),('color',color)]:
                if classification.get(key)!=value:classification[key]=value;changed=True
            for key,value in [('objectType','annotation'),('name',f'{label}_{index+1}')]:
                if key not in prop:prop[key]=value;changed=True
        if changed:
            info['previous_unaccepted_sha256']=info['sha256']
            with gzip.open(target,'wt') as stream:json.dump(doc,stream,separators=(',',':'),allow_nan=False)
        info['sha256']=sha256(target)
        info['positive_area_um2']=info['positive_pixels']*pixel**2
    if 'Duct' in record['masks'] and 'DuctEnvelope' in record['masks']:
        record['masks']['Duct']['envelope']=dict(record['masks']['DuctEnvelope'])
    record.update(pixel_um=pixel,analysis_pixel_size_um=pixel,
                  acquired_ome_pixel_size_um=record.get('acquired_ome_pixel_size_um',old),
                  calibration_source=str(calibration_path),calibration_sha256=sha256(calibration_path))
    write_json(path,record)
    registry.update(masks=record['masks'],proposal_sha256=sha256(path));write_json(registry_path,registry)
    return dict(case_id=cid,status='normalized_native_calibration_and_QuPath_properties')


def companion_review(args):
    root,cid=args;root=Path(root);path=root/'structural'/cid/'proposal.json';record=json.loads(path.read_text())
    names=record['channel_names'];patches={}
    for marker in sorted(set(sum(COMPANIONS.values(),[]))):
        data=read_channel(record['source_image'],names,marker)
        for label,markers in COMPANIONS.items():
            if marker not in markers:continue
            for index,crop in enumerate(record['reviews'][label]['crops']):
                x,y,w,h=[crop[k] for k in ['x','y','width','height']]
                patches[label,index,marker]=np.array(data[y:y+h,x:x+w])
        del data
    for label,markers in COMPANIONS.items():
        montage=Image.open(record['reviews'][label]['path']).convert('RGB')
        canvas=Image.new('RGB',((len(markers)+1)*256,1096),'white');pen=ImageDraw.Draw(canvas)
        pen.text((8,8),f'Donor {cid} / {label}: candidate overlay and independent native companion channels',fill='black')
        pen.text((8,30),'Fixed source locations; per-crop display ranges. Masks remain unvalidated pending review.',fill='black')
        for j,marker in enumerate(['Candidate overlay']+markers):pen.text((j*256+8,53),marker,fill='black')
        for index in range(4):
            canvas.paste(montage.crop((256,72+index*256,512,72+(index+1)*256)),(0,72+index*256))
            for j,marker in enumerate(markers):
                patch=display_channel(patches[label,index,marker]);canvas.paste(Image.fromarray(patch).convert('RGB'),((j+1)*256,72+index*256))
        destination=path.parent/f'review_{label}_companions.png';canvas.save(destination)
        record['reviews'][label]['companion_path']=str(destination)
        record['reviews'][label]['companion_sha256']=sha256(destination)
        record['reviews'][label]['companion_markers']=markers
    if 'DuctEnvelope' in record['reviews']:record['reviews']['DuctEnvelope']=dict(record['reviews']['Duct'])
    fresh=json.loads(path.read_text());fresh['reviews']=record['reviews'];write_json(path,fresh)
    registry_path=root/'structural'/f'{cid}.json';registry=json.loads(registry_path.read_text())
    registry['proposal_sha256']=sha256(path);write_json(registry_path,registry)
    print(cid,'companion review complete',flush=True)
    return dict(case_id=cid,status='companion_reviews_rendered')


def object_review(args):
    """Target rare inherited mask objects, not only random tissue intensities."""
    root,cid,labels=args;root=Path(root);path=root/'structural'/cid/'proposal.json';record=json.loads(path.read_text())
    from shapely import union_all,STRtree
    tissue_doc=json.loads(Path(record['tissue_source']).read_text())
    tissue=union_all([geometry.shape(f['geometry']) for f in tissue_doc['features']]) if tissue_doc.get('type')=='FeatureCollection' else geometry.shape(tissue_doc.get('geometry',tissue_doc))
    exclusion_path=root/'typing/eligibility/technical_exclusions.json'
    if exclusion_path.exists():
        excluded=[geometry.shape(x['geometry']) for x in json.loads(exclusion_path.read_text()).get('regions',[]) if str(x['case_id'])==cid]
        if excluded:tissue=tissue.difference(union_all(excluded))
    reviews={}
    for label in labels:
        with gzip.open(record['masks'][label]['path'],'rt') as stream:doc=json.load(stream)
        polygons=[geometry.shape(f['geometry']) for f in doc['features']]
        assessable=[(i,p.intersection(tissue)) for i,p in enumerate(polygons) if p.intersects(tissue)]
        assessable=[(i,p) for i,p in assessable if not p.is_empty and p.area>0]
        assessable.sort(key=lambda item:item[1].area)
        if not assessable:continue
        selected=[assessable[min(len(assessable)-1,round((len(assessable)-1)*q))] for q in [.5,.9,1]]
        markers=COMPANIONS[label];patches={};crops=[];size=256;h,w=record['shape_yx']
        for index,poly in selected:
            point=poly.representative_point();x=max(0,min(w-size,round(point.x)-128));y=max(0,min(h-size,round(point.y)-128))
            crops.append(dict(x=x,y=y,width=size,height=size,feature_index=index,assessable_object_area_px=poly.area))
        for marker in markers:
            plane=read_channel(record['source_image'],record['channel_names'],marker)
            for index,crop in enumerate(crops):patches[index,marker]=np.array(plane[crop['y']:crop['y']+size,crop['x']:crop['x']+size])
            del plane
        canvas=Image.new('RGB',((len(markers)+1)*size,3*size+72),'white');pen=ImageDraw.Draw(canvas)
        pen.text((8,8),f'Donor {cid} / {label}: objects at median, 90th percentile and largest assessable area',fill='black')
        pen.text((8,30),'Candidate outlines are compared with native primary and companion markers; disease labels unused.',fill='black')
        for j,marker in enumerate(['Candidate overlay']+markers):pen.text((j*size+8,53),marker,fill='black')
        tree=STRtree(polygons)
        for index,crop in enumerate(crops):
            x,y=crop['x'],crop['y'];raw=display_channel(patches[index,markers[0]])
            overlay=Image.fromarray(raw).convert('RGB');draw=ImageDraw.Draw(overlay)
            for k in tree.query(box(x,y,x+size,y+size),predicate='intersects'):
                for part in get_parts(polygons[k]):
                    if part.geom_type!='Polygon':continue
                    for ring in [part.exterior]+list(part.interiors):
                        draw.line([(xx-x,yy-y) for xx,yy in ring.coords],fill=(0,255,255),width=1)
            canvas.paste(overlay,(0,72+index*size))
            for j,marker in enumerate(markers):canvas.paste(Image.fromarray(display_channel(patches[index,marker])).convert('RGB'),((j+1)*size,72+index*size))
        destination=path.parent/f'review_{label}_objects.png';canvas.save(destination)
        reviews[label]=dict(object_path=str(destination),object_sha256=sha256(destination),object_crops=crops,object_sampling='median,90th,maximum assessable mask-object area')
    fresh=json.loads(path.read_text())
    for label,review in reviews.items():fresh['reviews'][label].update(review)
    write_json(path,fresh)
    registry_path=root/'structural'/f'{cid}.json';registry=json.loads(registry_path.read_text());registry['proposal_sha256']=sha256(path);write_json(registry_path,registry)
    print(cid,'object review complete',flush=True)
    return dict(case_id=cid,structures=labels,status='object_review_rendered')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=['propose','accept','normalize','companions','objects'])
    parser.add_argument('--input-root',type=Path,required=True);parser.add_argument('--donor',action='append');parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--sample-pixels',type=int,default=100000);parser.add_argument('--review',type=Path)
    parser.add_argument('--structure',action='append',choices=list(INHERITED))
    args=parser.parse_args();root=args.input_root.resolve(strict=True)
    if args.command=='accept':accept(root,args.review)
    else:
        if args.workers<1 or args.workers>4:raise ValueError('Maximum four workers with two TIFF threads each')
        cohort=pd.read_csv(root/'cohort.csv',dtype=str);donors=args.donor or cohort.case_id.tolist()
        if set(donors)-set(cohort.case_id):raise ValueError('Unknown donor')
        (root/'structural').mkdir(exist_ok=True)
        results=[]
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            action=propose_one if args.command=='propose' else object_review if args.command=='objects' else companion_review
            pending={pool.submit(normalize_existing,root,cid) if args.command=='normalize' else
                     pool.submit(action,(root,cid,args.sample_pixels) if args.command=='propose' else
                                 (root,cid,args.structure or ['Lymphatic']) if args.command=='objects' else (root,cid)):cid for cid in donors}
            for future in as_completed(pending):
                try:results.append(future.result())
                except Exception as error:results.append(dict(case_id=pending[future],failed=str(error)));print(results[-1],flush=True)
        write_json(root/'structural'/f'{args.command}_run.json',dict(donors=results,workers=args.workers,scientific_acceptance=False))
        if any('failed' in r for r in results):raise SystemExit(1)
