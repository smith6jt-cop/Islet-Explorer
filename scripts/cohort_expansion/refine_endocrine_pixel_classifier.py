#!/usr/bin/env python3
"""Fit an endocrine pixel classifier from signal-supported original-cohort regions.

Original islet annotations locate references, not pixel-level ground truth. SAM
geometry is never a prediction input. Evaluation reports reference agreement,
not biological accuracy. All destinations must be new directories.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import multiprocessing
from pathlib import Path
import platform
import shutil
import time

import cv2
import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage as ndi
import shapely
from shapely.geometry import shape
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, roc_auc_score
import tifffile

from prepare_endocrine_candidates import digest, write_json

ORIGINAL_DONORS = {'112','6356','6450','6479','6505','6516','6521','6533',
                   '6534','6538','6548','6549','6550','6551','6563'}
MARKERS = ['INS', 'GCG', 'SST']
SEED = 20260917
DOWNSAMPLE = 2
PIXEL_UM = .3774
CONFIG = dict(markers=MARKERS, source_pixel_size_um=PIXEL_UM, downsample=DOWNSAMPLE,
              gaussian_sigmas_pixels=[1.,4.], foreground_erode_um=2., background_gap_um=10.,
              signal_reference_background_quantile=.995, training_pixels_per_class_per_donor=6000,
              forest_trees=160, forest_max_depth=12, forest_min_sample_count=40,
              seed=SEED, probability_cutoff=.5, minimum_region_area_um2=80.,
              closing_radius_um=2., maximum_filled_hole_um2=80.)
MODEL_CONFIG = dict(CONFIG, intensity_normalization='Per-image tissue median and IQR, sampled every second 2x image pixel; scale floor 1 intensity unit')


def donor_seed(cid):
    return (int(hashlib.sha256(str(cid).encode()).hexdigest()[:8],16)+SEED) % (2**32)


def sources(args):
    """Extract verified image/tissue inputs without SAM prompts or cell gates."""
    out=args.output_root.resolve();out.mkdir(parents=True,exist_ok=False)
    proof_path=args.source_root/'source-image-verification.json'
    proof=json.loads(proof_path.read_text());records=[]
    for check in proof['images']:
        cid=check['case_id'];prior=args.source_root/cid
        spec=json.loads((prior/'prompts.json').read_text())
        path=Path(spec['source_image']).resolve(strict=True);stat=path.stat()
        if str(path)!=check['path'] or stat.st_size!=check['bytes'] or stat.st_mtime_ns!=check['mtime_ns']:
            raise ValueError('Source differs from reconciled acquisition')
        if spec['pixel_size_um']!=PIXEL_UM:raise ValueError('Source calibration differs')
        folder=out/cid;folder.mkdir()
        clean={k:spec[k] for k in ['case_id','source_image','image_shape_cyx','pixel_size_um','source_channel_names']}
        clean.update(source_bytes=stat.st_size,source_mtime_ns=stat.st_mtime_ns,
            previously_verified_image_sha256=check['previously_verified_sha256'],full_image_hash_recomputed=False,
            source_identity_verification=str(proof_path.resolve()),source_identity_verification_sha256=digest(proof_path),
            prior_image_metadata_sha256=digest(prior/'prompts.json'),SAM_prompts_included=False,cell_gates_included=False)
        shutil.copyfile(prior/'valid_tissue.geojson',folder/'valid_tissue.geojson')
        clean['valid_tissue_sha256']=digest(folder/'valid_tissue.geojson')
        write_json(folder/'source.json',clean);records.append(clean)
    write_json(out/'source_manifest.json',dict(status='image_and_tissue_inputs_verified',donors=records))


def validate_training_donors(records):
    ids=[str(r['case_id']) for r in records]
    if len(ids)!=15 or set(ids)!=ORIGINAL_DONORS:
        raise ValueError('Training requires exactly the original 15 donors; target donors are forbidden')


def read_endocrine(spec):
    """Read only three channels; use a native 2x pyramid or exact 2x area average."""
    names=spec['source_channel_names']
    indices=[names.index(m) for m in MARKERS]
    with tifffile.TiffFile(spec['source_image']) as tf:
        series=tf.series[0]
        if series.axes!='CYX' or tuple(series.shape)!=tuple(spec['image_shape_cyx']):
            raise ValueError('Source axes/dimensions differ from the manifest')
        if series.shape[1] % DOWNSAMPLE or series.shape[2] % DOWNSAMPLE:
            raise ValueError('Exact integer coordinate scaling requires even dimensions')
        if len(series.levels)>1 and series.levels[1].shape[1:]==tuple(x//2 for x in series.shape[1:]):
            raw=np.stack([series.levels[1].pages[i].asarray() for i in indices],-1).astype(np.float32)
            reader='native TIFF 2x pyramid'
        else:
            mapped=tifffile.memmap(spec['source_image'],mode='r')
            # OpenCV ignores NumPy's non-native byte-order flag. Convert before
            # calling it: old ImageJ TIFFs are big-endian, target images may not be.
            raw=np.stack([cv2.resize(np.asarray(mapped[i],dtype=np.float32),(mapped.shape[2]//2,mapped.shape[1]//2),
                                     interpolation=cv2.INTER_AREA) for i in indices],-1).astype(np.float32)
            reader='read-only TIFF memory map; native-endian float32 conversion before 2x area averaging'
    return raw,reader


def refuse_invalid_samples(root):
    if (root/'quality_status.json').exists():
        status=json.loads((root/'quality_status.json').read_text())['status']
        if status.startswith(('failed','withdrawn','invalid')):raise ValueError(f'Invalid samples: {root} ({status})')


def features(raw):
    return np.concatenate([cv2.GaussianBlur(raw,(0,0),sigmaX=s,sigmaY=s,
                                           borderType=cv2.BORDER_REFLECT_101)
                           for s in CONFIG['gaussian_sigmas_pixels']],axis=-1)


def rasterize(geometries,shape_yx,downsample=DOWNSAMPLE,center_offset=None):
    """Rasterize each object separately so one object's hole cannot erase another."""
    output=np.zeros(shape_yx,np.uint8)
    for geom in geometries:
        geom=shapely.make_valid(geom)
        x0,y0,x1,y1=geom.bounds
        left=max(0,int(np.floor(x0/downsample))-1);top=max(0,int(np.floor(y0/downsample))-1)
        right=min(shape_yx[1],int(np.ceil(x1/downsample))+1);bottom=min(shape_yx[0],int(np.ceil(y1/downsample))+1)
        if right<=left or bottom<=top:continue
        # Pixel centers correspond to the 2x area-averaged image, not raw array indices.
        yy,xx=np.mgrid[top:bottom,left:right]
        offset=downsample/2 if center_offset is None else center_offset
        hit=shapely.intersects_xy(geom,xx*downsample+offset,yy*downsample+offset)
        output[top:bottom,left:right] |= hit.astype(np.uint8)
    return output.astype(bool)


def normalization(raw_overview,tissue_overview):
    q25,median,q75=np.quantile(raw_overview[tissue_overview],[.25,.5,.75],axis=0)
    return dict(markers=MARKERS,offset=median.tolist(),scale=np.maximum(q75-q25,1.).tolist(),
                method=MODEL_CONFIG['intensity_normalization'],tissue_sample_pixels=int(tissue_overview.sum()))


def normalize_features(feat,norm):
    return ((feat-np.tile(norm['offset'],2))/np.tile(norm['scale'],2)).astype(np.float32)


def training_masks(islets,tissue,smooth,pixel_um=PIXEL_UM*DOWNSAMPLE):
    inside=ndi.distance_transform_edt(islets)*pixel_um
    outside=ndi.distance_transform_edt(~islets)*pixel_um
    interior=tissue & (inside>=CONFIG['foreground_erode_um'])
    background=tissue & (outside>=CONFIG['background_gap_um'])
    if background.sum()<10000:raise ValueError('Insufficient background reference area')
    cutoff=np.quantile(smooth[background],CONFIG['signal_reference_background_quantile'],axis=0)
    cutoff=np.maximum(cutoff,1.)
    supported=(smooth>cutoff).any(axis=-1)
    positive=interior & supported
    negative=background & ~supported
    if positive.sum()<6000 or negative.sum()<6000:raise ValueError('Insufficient clear reference pixels')
    return positive,negative,interior,background,cutoff


def balanced_sample(mask,groups,n,rng):
    """Equal group representation where possible; avoid letting one large islet dominate."""
    flat=np.flatnonzero(mask)
    if not len(flat):return np.empty(0,dtype=int)
    values=groups.ravel()[flat]
    unique=np.unique(values)
    quota=max(1,n//len(unique))
    picks=[]
    for g in unique:
        candidates=flat[values==g]
        picks.extend(rng.choice(candidates,size=min(quota,len(candidates)),replace=False).tolist())
    if len(picks)>n:picks=rng.choice(picks,size=n,replace=False).tolist()
    if len(picks)<n:
        remainder=np.setdiff1d(flat,np.asarray(picks,dtype=int))
        picks.extend(rng.choice(remainder,size=min(n-len(picks),len(remainder)),replace=False).tolist())
    return np.asarray(picks,dtype=int)


def prepare(args):
    reference=json.loads((args.reference_root/'reference_manifest.json').read_text())
    validate_training_donors(reference['donors'])
    out=args.output_root.resolve();out.mkdir(parents=True,exist_ok=False)
    records=[]
    for spec in sorted(reference['donors'],key=lambda r:int(r['case_id'])):
        cid=spec['case_id'];folder=out/cid;folder.mkdir()
        if digest(Path(spec['reference_annotations']))!=spec['reference_sha256']:
            raise ValueError('Reference annotations changed')
        source=Path(spec['source_image']);stat=source.stat()
        if stat.st_size!=spec['source_image_bytes'] or stat.st_mtime_ns//1000000!=spec['source_image_mtime_ms']:
            raise ValueError('Reference source image changed')
        print(f'{cid}: reading endocrine channels',flush=True)
        raw,reader=read_endocrine(spec);feat=features(raw)
        refs=json.loads(Path(spec['reference_annotations']).read_text())['features']
        classes={k:[shape(f['geometry']) for f in refs if f['properties']['classification']['name']==k]
                 for k in ['Islet','Tissue']}
        islets=rasterize(classes['Islet'],raw.shape[:2]);tissue=rasterize(classes['Tissue'],raw.shape[:2])
        positive,negative,interior,background,cutoff=training_masks(islets,tissue,feat[...,:3])
        rng=np.random.default_rng(donor_seed(cid));n=CONFIG['training_pixels_per_class_per_donor']
        labels,_=ndi.label(islets)
        # Stratify foreground by hormone with greatest excess over background.
        dominant=np.argmax(feat[...,:3]/cutoff,axis=-1)
        chosen=[]
        for marker in range(3):
            chosen.extend(balanced_sample(positive & (dominant==marker),labels,n//3,rng).tolist())
        remaining=np.setdiff1d(np.flatnonzero(positive),np.asarray(chosen,dtype=int))
        if len(chosen)<n:chosen.extend(rng.choice(remaining,size=n-len(chosen),replace=False).tolist())
        pos=np.asarray(chosen,dtype=int)
        yy,xx=np.indices(raw.shape[:2]);tiles=(yy//256)*(raw.shape[1]//256+1)+xx//256
        neg=balanced_sample(negative,tiles,n,rng)
        train=np.concatenate([pos,neg]);target=np.r_[np.ones(n,np.int32),np.zeros(n,np.int32)]
        # A separate random sample includes ALL reference interiors/background,
        # including ambiguous signal. It exposes disagreement hidden by clear-pixel selection.
        eval_idx=np.r_[rng.choice(np.flatnonzero(interior),min(12000,interior.sum()),replace=False),
                       rng.choice(np.flatnonzero(background),min(12000,background.sum()),replace=False)]
        eval_y=interior.ravel()[eval_idx].astype(np.int32)
        np.savez_compressed(folder/'samples.npz',x=feat.reshape(-1,6)[train],y=target,
            linear_indices=train,eval_x=feat.reshape(-1,6)[eval_idx],eval_y=eval_y,
            eval_linear_indices=eval_idx,shape_yx=np.array(raw.shape[:2]),
            dominant_marker=dominant.ravel()[pos])
        # Low-resolution audit layers retain all reference pixels, not just the sample.
        Image.fromarray((islets[::2,::2]*255).astype('uint8')).save(folder/'source_islet_regions.png')
        Image.fromarray((positive[::2,::2]*255).astype('uint8')).save(folder/'supported_training_regions.png')
        np.save(folder/'endocrine_overview.npy',raw[::2,::2])
        record=dict(case_id=cid,source_reference_manifest=str(args.reference_root/cid/'reference_manifest.json'),
            reference_sha256=spec['reference_sha256'],source_native_sha256=spec['source_native_sha256'],reader=reader,
            training_foreground=n,training_background=n,positive_reference_pixels=int(positive.sum()),
            ambiguous_interior_pixels=int((interior & ~positive).sum()),
            total_reference_interior_pixels=int(interior.sum()),negative_reference_pixels=int(negative.sum()),
            ambiguous_background_pixels=int((background & ~negative).sum()),
            foreground_stratum_counts={m:int((dominant.ravel()[pos]==i).sum()) for i,m in enumerate(MARKERS)},
            reference_signal_cutoffs=dict(zip(MARKERS,cutoff.tolist())),sample_sha256=digest(folder/'samples.npz'))
        records.append(record);write_json(folder/'sample_manifest.json',record)
        print(f'{cid}: {n}+{n} sampled pixels; {record["ambiguous_interior_pixels"]}/{int(interior.sum())} interior pixels excluded from foreground training',flush=True)
        del raw,feat,islets,tissue,positive,negative,interior,background,labels,dominant,yy,xx,tiles
    write_json(out/'sampling_manifest.json',dict(status='original15_weakly_supervised_references',config=CONFIG,donors=records,
        reference_manifest_sha256=digest(args.reference_root/'reference_manifest.json'),
        training_target='Hormone-supported pixels in reference islet interiors versus clear tissue background',
        no_manual_pixel_ground_truth=True,held_out_new_donors=['78','6090','6232','6520','6539','6558']))


def harmonize(args):
    """Preserve raw samples and normalize each donor independently of its labels."""
    refuse_invalid_samples(args.samples_root)
    out=args.output_root.resolve();out.mkdir(parents=True,exist_ok=False)
    source=json.loads((args.samples_root/'sampling_manifest.json').read_text())
    validate_training_donors(source['donors'])
    if source['config']!=CONFIG:raise ValueError('Raw sampling configuration differs')
    records=[]
    for record in source['donors']:
        cid=record['case_id'];folder=out/cid;folder.mkdir()
        rawdir=args.samples_root/cid
        if digest(rawdir/'samples.npz')!=record['sample_sha256']:raise ValueError('Raw sample changed')
        ref=json.loads(Path(record['source_reference_manifest']).read_text())
        refs=json.loads(Path(ref['reference_annotations']).read_text())['features']
        tissue=[shape(f['geometry']) for f in refs if f['properties']['classification']['name']=='Tissue']
        overview=np.load(rawdir/'endocrine_overview.npy')
        mask=rasterize(tissue,overview.shape[:2],downsample=4,center_offset=1)
        norm=normalization(overview,mask)
        data=dict(np.load(rawdir/'samples.npz'))
        data['x']=normalize_features(data['x'],norm);data['eval_x']=normalize_features(data['eval_x'],norm)
        np.savez_compressed(folder/'samples.npz',**data)
        write_json(folder/'normalization.json',norm)
        updated=dict(record,raw_samples_root=str(rawdir.resolve()),raw_sample_sha256=record['sample_sha256'],
                     sample_sha256=digest(folder/'samples.npz'),normalization=norm)
        records.append(updated);write_json(folder/'sample_manifest.json',updated)
        print(f'{cid}: tissue median {norm["offset"]}, IQR {norm["scale"]}',flush=True)
    source.update(config=MODEL_CONFIG,donors=records,raw_sampling_manifest_sha256=digest(args.samples_root/'sampling_manifest.json'))
    write_json(out/'sampling_manifest.json',source)


def fit(x,y,seed=SEED):
    cv2.setNumThreads(1);cv2.setRNGSeed(int(seed))
    model=cv2.ml.RTrees_create();model.setMaxDepth(CONFIG['forest_max_depth'])
    model.setMinSampleCount(CONFIG['forest_min_sample_count']);model.setActiveVarCount(3)
    model.setCalculateVarImportance(True)
    model.setTermCriteria((cv2.TERM_CRITERIA_MAX_ITER,CONFIG['forest_trees'],0))
    model.train(np.ascontiguousarray(x,dtype=np.float32),cv2.ml.ROW_SAMPLE,np.asarray(y,dtype=np.int32))
    return model


def probability(model,x):
    votes=model.getVotes(np.ascontiguousarray(x,dtype=np.float32),0)
    col=np.flatnonzero(votes[0]==1)
    if len(col)!=1:raise ValueError('Missing foreground class in forest votes')
    return votes[1:,col[0]]/votes[1:].sum(axis=1)


def evaluate(y,p):
    tn,fp,fn,tp=confusion_matrix(y,p>=.5,labels=[0,1]).ravel()
    return dict(roc_auc=float(roc_auc_score(y,p)),balanced_accuracy=float(balanced_accuracy_score(y,p>=.5)),
        sensitivity=float(tp/(tp+fn)),specificity=float(tn/(tn+fp)),true_positive=int(tp),
        false_positive=int(fp),true_negative=int(tn),false_negative=int(fn))


def fit_fold(cid,training_ids,x,y,held,out):
    model=fit(x,y)
    clear=evaluate(held['y'],probability(model,held['x']))
    all_refs=evaluate(held['eval_y'],probability(model,held['eval_x']))
    model.save(str(out/f'held_out_{cid}_forest.json'))
    return dict(case_id=cid,training_donors=';'.join(training_ids),
        **{'clear_reference_'+k:v for k,v in clear.items()},
        **{'all_reference_'+k:v for k,v in all_refs.items()})


def train(args):
    refuse_invalid_samples(args.samples_root)
    out=args.output_root.resolve();out.mkdir(parents=True,exist_ok=False)
    sampling=json.loads((args.samples_root/'sampling_manifest.json').read_text());validate_training_donors(sampling['donors'])
    if sampling['config']!=MODEL_CONFIG:raise ValueError('Configuration changed since harmonization')
    data={}
    for record in sampling['donors']:
        path=args.samples_root/record['case_id']/'samples.npz'
        if digest(path)!=record['sample_sha256']:raise ValueError('Samples changed')
        data[record['case_id']]=dict(np.load(path))
    reports=[]
    # Independent donor fits run in separate processes; each has its own seeded
    # OpenCV RNG. Separate output files and deterministic report order avoid races.
    with ProcessPoolExecutor(max_workers=4,mp_context=multiprocessing.get_context('spawn')) as pool:
        pending=[]
        for cid,held in data.items():
            training_ids=[donor for donor in data if donor!=cid];fitting=[data[donor] for donor in training_ids]
            pending.append(pool.submit(fit_fold,cid,training_ids,np.concatenate([d['x'] for d in fitting]),
                                       np.concatenate([d['y'] for d in fitting]),held,out))
        for future in as_completed(pending):
            report=future.result();reports.append(report)
            print(f'Leave donor {report["case_id"]} out: clear-reference balanced agreement {report["clear_reference_balanced_accuracy"]:.3f}; all-region sensitivity {report["all_reference_sensitivity"]:.3f}',flush=True)
    reports.sort(key=lambda r:int(r['case_id']))
    model=fit(np.concatenate([d['x'] for d in data.values()]),np.concatenate([d['y'] for d in data.values()]))
    model.save(str(out/'endocrine_forest.json'))
    pd.DataFrame(reports).to_csv(out/'leave_one_donor_out.csv',index=False)
    manifest=dict(status='reference_trained_pixel_classifier_requires_boundary_review',config=MODEL_CONFIG,
        training_donors=list(data),training_pixels=sum(len(d['y']) for d in data.values()),
        target_donors_excluded_from_training=sampling['held_out_new_donors'],
        feature_names=[f'{m}_gaussian_sigma_{s:g}px' for s in CONFIG['gaussian_sigmas_pixels'] for m in MARKERS],
        feature_importance=model.getVarImportance().ravel().tolist(),
        forest_sha256=digest(out/'endocrine_forest.json'),sampling_manifest_sha256=digest(args.samples_root/'sampling_manifest.json'),
        validation_interpretation='Leave-one-donor-out agreement with weak reference labels; not manual boundary accuracy',
        manual_ground_truth=False,accepted_islets=0,python_version=platform.python_version(),opencv_version=cv2.__version__,
        numpy_version=np.__version__,tifffile_version=tifffile.__version__,created_at_utc=datetime.now(timezone.utc).isoformat())
    write_json(out/'training_manifest.json',manifest)


def clean_regions(foreground,tissue,pixel_um=PIXEL_UM*DOWNSAMPLE):
    """Bounded closing and small-hole filling; preserve broad exterior and interior gaps."""
    r=CONFIG['closing_radius_um']/pixel_um
    n=int(np.ceil(r));yy,xx=np.mgrid[-n:n+1,-n:n+1];kernel=(xx*xx+yy*yy<=r*r).astype(np.uint8)
    closed=cv2.morphologyEx(foreground.astype(np.uint8),cv2.MORPH_CLOSE,kernel).astype(bool) & tissue
    # Fill enclosed holes only if small; never bridge large voids with a convex hull.
    holes=ndi.binary_fill_holes(closed) & ~closed
    lab,num=ndi.label(holes);sizes=np.bincount(lab.ravel())
    fill=(sizes*pixel_um**2<=CONFIG['maximum_filled_hole_um2']);fill[0]=False
    closed |= fill[lab] & tissue
    lab,num=ndi.label(closed);sizes=np.bincount(lab.ravel())
    keep=(sizes*pixel_um**2>=CONFIG['minimum_region_area_um2']);keep[0]=False
    return keep[lab],closed & ~keep[lab]


def predict(args):
    """Save fixed foreground and a legacy uniform comparator.

    Native tracing followed by build_size_adaptive_islets.groovy is required
    before new predictions can be imported as review boundaries.
    """
    started=time.monotonic()
    out=args.output_root.resolve();out.mkdir(parents=True,exist_ok=False)
    spec=json.loads(args.specification.read_text())
    source=Path(spec['source_image']);source_stat=source.stat()
    if source_stat.st_size!=spec['source_bytes'] or source_stat.st_mtime_ns!=spec['source_mtime_ns']:
        raise ValueError('Target source image fingerprint changed')
    if digest(args.specification.parent/'valid_tissue.geojson')!=spec['valid_tissue_sha256']:
        raise ValueError('Target tissue geometry changed')
    trained=json.loads((args.model_root/'training_manifest.json').read_text())
    trees=int(trained['config']['forest_trees'])
    minimum_votes=getattr(args,'minimum_votes',None)
    minimum_votes=trees//2 if minimum_votes is None else minimum_votes
    if not 1 <= minimum_votes <= trees:
        raise ValueError(f'Minimum votes must be between 1 and {trees}')
    cutoff=minimum_votes/trees
    training_donor = str(spec['case_id']) in trained['training_donors']
    deployment = bool(getattr(args, 'deployment', False))
    if training_donor and not deployment:
        raise ValueError('Use donor-held-out model for evaluation, or explicitly select --deployment for frozen-model production prediction')
    path=args.model_root/'endocrine_forest.json'
    if digest(path)!=trained['forest_sha256']:raise ValueError('Model changed')
    model=cv2.ml.RTrees_load(str(path));cv2.setNumThreads(getattr(args,'threads',8))
    raw,reader=read_endocrine(spec)
    tissue_geom=shape(json.loads((args.specification.parent/'valid_tissue.geojson').read_text()))
    tissue=rasterize([tissue_geom],raw.shape[:2])
    norm=normalization(raw[::2,::2],tissue[::2,::2])
    feat=normalize_features(features(raw),norm)
    p=np.empty(raw.shape[:2],np.float32)
    for start in range(0,len(p),128):
        if start%1024==0:print(f'{spec["case_id"]}: classifying row {start}/{len(p)}',flush=True)
        part=feat[start:start+128];p[start:start+128]=probability(model,part.reshape(-1,6)).reshape(part.shape[:2])
    foreground=(p>=cutoff) & tissue
    regions,foci=clean_regions(foreground,tissue)
    for name,array in [('pixel_classes',foreground),('boundary_regions',regions),('small_foci',foci)]:
        Image.fromarray(array.astype('uint8')).save(out/f'{name}.png')
    np.save(out/'foreground_probability.npy',p.astype(np.float16))
    if source.stat().st_mtime_ns!=source_stat.st_mtime_ns or source.stat().st_size!=source_stat.st_size:
        raise ValueError('Source image changed during prediction')
    # Vectorize the exact class raster with QuPath in a separate read-only step.
    manifest=dict(status='refined_pixel_prediction_unreviewed',case_id=str(spec['case_id']),source_image=spec['source_image'],
        source_specification=str(args.specification.resolve()),source_specification_sha256=digest(args.specification),
        source_pixel_size_um=PIXEL_UM,prediction_downsample=DOWNSAMPLE,prediction_width=raw.shape[1],prediction_height=raw.shape[0],
        source_width=spec['image_shape_cyx'][2],source_height=spec['image_shape_cyx'][1],reader=reader,markers=MARKERS,
        model_manifest=str((args.model_root/'training_manifest.json').resolve()),model_sha256=trained['forest_sha256'],
        source_bytes=source_stat.st_size,source_mtime_ns=source_stat.st_mtime_ns,
        previously_verified_image_sha256=spec['previously_verified_image_sha256'],full_image_hash_recomputed=False,
        foreground_area_um2=float(foreground.sum()*(PIXEL_UM*DOWNSAMPLE)**2),
        raw_connected_regions=int(ndi.label(foreground)[1]),boundary_connected_regions=int(ndi.label(regions)[1]),
        small_focus_regions=int(ndi.label(foci)[1]),boundary_area_um2=float(regions.sum()*(PIXEL_UM*DOWNSAMPLE)**2),
        config=dict(MODEL_CONFIG,probability_cutoff=cutoff),normalization=norm,SAM_used=False,SAM_used_for_prediction=False,
        operating_point=dict(minimum_votes=minimum_votes,forest_trees=trees,
            selection_status='Fixed requested cutoff; donor visual acceptance is recorded separately'),
        boundary_policy=dict(active='io60_size_adaptive_v1',stage='pending_native_refinement',
            uniform_boundary_raster_role='legacy comparison only'),
        training_region_provenance='Signal-supported interiors of original15 annotations carrying SAM Quality',
        cell_phenotype_gates_used=False,accepted_islets=0,runtime_seconds=time.monotonic()-started,
        original_project_modified=False,source_native_objects_modified=False,
        prediction_purpose='production_deployment' if deployment else 'held_out_target_prediction',
        donor_present_in_training=training_donor,
        independent_validation_result=False)
    manifest['raster_sha256']={name:digest(out/name) for name in ['pixel_classes.png','boundary_regions.png','small_foci.png']}
    write_json(out/'prediction_manifest.json',manifest)
    print(json.dumps(manifest),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('sources');p.add_argument('--source-root',type=Path,required=True);p.set_defaults(func=sources)
    p=sub.add_parser('prepare');p.add_argument('--reference-root',type=Path,required=True);p.set_defaults(func=prepare)
    p=sub.add_parser('harmonize');p.add_argument('--samples-root',type=Path,required=True);p.set_defaults(func=harmonize)
    p=sub.add_parser('train');p.add_argument('--samples-root',type=Path,required=True);p.set_defaults(func=train)
    p=sub.add_parser('predict');p.add_argument('--model-root',type=Path,required=True);p.add_argument('--specification',type=Path,required=True)
    p.add_argument('--minimum-votes',type=int,help='Fixed minimum forest votes; default: half the trees')
    p.add_argument('--deployment',action='store_true',help='Apply the frozen model to production images, including training donors; never represents held-out validation')
    p.add_argument('--threads',type=int,default=8)
    p.set_defaults(func=predict)
    for p in sub.choices.values():p.add_argument('--output-root',type=Path,required=True)
    args=parser.parse_args();args.func(args)


if __name__=='__main__':main()
