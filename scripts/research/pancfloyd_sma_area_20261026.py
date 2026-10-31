#!/usr/bin/env python3
"""SMA-positive area counterparts for reference arteriole cell-count panels.

Reuse the approved H7 per-image Otsu/2 µm² component method and current
geometry. This measures marker-positive area, without asserting arteriole identity.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/pancfloyd_mpl")
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import shapely
from skimage.filters import threshold_otsu

import h7_nerve_features_20261026 as h7
import pancfloyd_slides1_4_20261026 as panels

ROOT=panels.ROOT; OUT=panels.OUT; DATA=panels.DATA; FIG=panels.FIG
SRC=ROOT/"data/development/uniform-20donor-20261015"
APP=ROOT/"data/development/uniform-20donor-20261022-v10"
CACHE=DATA/"sma_by_donor"
PX=h7.PX


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def extract_donor(task):
    tag,entry=task;case=int(entry["case_id"])
    dest=CACHE/f"{case}.csv";prov=CACHE/f"{case}.json"
    if dest.exists() and prov.exists():
        return {"case_id":case,"status":"cached"}
    census=pd.read_parquet(panels.CENSUS_PATH)
    sub=census[census.case_id.eq(case)&census.object_class.isin(["Islet","EndocrineCluster"])].copy()
    cores,expanded,tissue,invalid,gp=h7.hf.geometry(APP,case)
    assert set(sub.object_name).issubset(cores)
    ep=SRC/"islets"/f"{tag}.geojson"
    endocrine=h7.geoms(ep);etree=shapely.STRtree(endocrine)
    names=h7.ct.channels(SRC,entry["entry_id"]);channel=names.index("SMA")
    img=h7.Channels(entry["image_path"])
    raw=img.arr[channel] if img.memmap else img(channel)
    sampled=np.asarray(raw[::8,::8]);samplemask=h7.pixel_mask(tissue,sampled.shape,step=8)
    sample=sampled[samplemask];sample=sample[np.isfinite(sample)]
    threshold=float(threshold_otsu(sample))
    rows=[]
    for r in sub.itertuples():
        core=cores[r.object_name][1];context=core.buffer(42)
        x0=max(0,int(np.floor(context.bounds[0]/PX)));y0=max(0,int(np.floor(context.bounds[1]/PX)))
        x1=min(raw.shape[1],int(np.ceil(context.bounds[2]/PX))+1);y1=min(raw.shape[0],int(np.ceil(context.bounds[3]/PX))+1)
        roi=np.asarray(raw[y0:y1,x0:x1])
        tm=h7.pixel_mask(tissue.intersection(context),roi.shape,x0,y0)
        cm=h7.pixel_mask(core,roi.shape,x0,y0)&tm
        endom=h7.pixel_mask(h7.hf.union_local(etree,endocrine,context),roi.shape,x0,y0)
        peri_geom=expanded[r.object_name].difference(core)
        pm=h7.pixel_mask(peri_geom,roi.shape,x0,y0)&tm&~endom
        positive=h7.component_filter((roi>threshold)&tm&np.isfinite(roi))
        for region,mask in [("core",cm),("peri",pm)]:
            denom=int(mask.sum())*PX**2;numer=int((positive&mask).sum())*PX**2
            rows.append(dict(islet_id=r.islet_id,case_id=case,donor_status=r.donor_status,
                             object_class=r.object_class,core_area_um2=r.core_area_um2,
                             core_complete=r.core_complete,neighborhood_complete_20=r.neighborhood_complete_20,
                             region=region,numerator=numer,denominator=denom,value=100*numer/denom if denom else np.nan,
                             sma_positive_pixels=int((positive&mask).sum()),assessed_pixels=int(mask.sum()),
                             threshold_raw=threshold,pixel_size_um=PX))
    result=pd.DataFrame(rows);result.to_csv(dest,index=False)
    source=Path(entry["image_path"]);stat=source.stat()
    provenance=dict(case_id=case,tag=tag,marker="SMA",threshold_method="Otsu on every eighth image pixel inside Tissue",
                    threshold_raw=threshold,threshold_sample_pixels=len(sample),channel_index=channel,
                    source_image=str(source),source_image_bytes=stat.st_size,source_image_mtime_ns=stat.st_mtime_ns,
                    image_shape=list(raw.shape),pixel_size_um=PX,minimum_positive_component_area_um2=2,
                    connected_component_neighborhood=8,component_filter_region="Tissue within object +42 µm context",
                    core_region="current core intersect Tissue",peri_region="current expanded minus core, intersect Tissue, exclude all endocrine cores",
                    geometry_sha256=sha(gp),all_endocrine_geometry_sha256=sha(ep),objects=len(sub),records=len(result),
                    quantitative_INS_used=False,cell_retyping=False,arteriole_identity_assumed=False)
    prov.write_text(json.dumps(provenance,indent=2)+"\n")
    del raw;img.close()
    return {"case_id":case,"status":"measured","objects":len(sub),"threshold":threshold}


def main():
    CACHE.mkdir(parents=True,exist_ok=True)
    entries=sorted(h7.ct.manifest(SRC).items())
    assert len(entries)==20 and all(int(e["case_id"])!=78 for _,e in entries)
    with ProcessPoolExecutor(max_workers=2) as pool:
        for result in pool.map(extract_donor,entries):
            print(result,flush=True)
    source=pd.concat([pd.read_csv(CACHE/f"{int(e['case_id'])}.csv") for _,e in entries],ignore_index=True)
    assert source.case_id.nunique()==20 and len(source)==9812
    assert not source.duplicated(["islet_id","region"]).any()
    assert (source.numerator>=0).all() and (source.numerator<=source.denominator).all()
    coverage=[]
    for region,panel in [("core","slide03_image07"),("peri","slide03_image08")]:
        metric=f"sma_{region}_area";d=source[source.region.eq(region)].copy();peri=region=="peri"
        d["eligible"]=d.core_complete & (d.neighborhood_complete_20 if peri else True)
        d=panels.bin_data(d);used,donor,summary=panels.summarize(metric,d)
        definition=("SMA-positive area replaces the reference's SMA+ arteriole typed-cell percentage. "
          "This amendment measures marker-positive area and does not identify arterioles or count smooth-muscle cells. "
          "Original SMA images were thresholded separately per donor with Otsu on every eighth pixel inside the current Tissue annotation; "
          "8-connected positive components smaller than 2 µm² were excluded, reusing the approved H7/local-paper image-area method. "
          "Marker classification was evaluated in a 42 µm context around each object; no INS channel was used. "
          f"Numerator is SMA-positive raster area; denominator is assessed {'20 µm peri-band' if peri else 'core'} raster area, both intersected with Tissue. "
          + ("All endocrine cores are excluded from the peri band. " if peri else "")
          + "Thresholds and source-image/geometry provenance are exported for every donor. ")
        panels.curve(panel,metric,f"SMA-positive area · {region}",f"SMA-positive area / {region} area (%)",d,summary,definition,peri)
        coverage.append(dict(metric=metric,source_objects=len(d),edge_eligible_objects=int(d.eligible.sum()),
                             plotted_objects=len(used),donors=used.case_id.nunique(),outside_display=int((d.eligible&~d.display_range).sum())))
    pd.DataFrame([json.loads((CACHE/f"{int(e['case_id'])}.json").read_text()) for _,e in entries]).to_csv(DATA/"sma_threshold_provenance.csv",index=False)
    covpath=DATA/"coverage.csv";cov=pd.read_csv(covpath);cov=cov[~cov.metric.str.startswith("sma_")]
    pd.concat([cov,pd.DataFrame(coverage)],ignore_index=True).to_csv(covpath,index=False)
    mpath=OUT/"panel_manifest.csv";manifest=pd.read_csv(mpath)
    for region,panel in [("core","slide03_image07"),("peri","slide03_image08")]:
        mask=manifest.panel.eq(panel)
        manifest.loc[mask,"status"]="explicit_area_endpoint_replacement"
        manifest.loc[mask,"metric"]=f"sma_{region}_area"
        manifest.loc[mask,"endpoint"]="SMA-positive raster area percent"
        manifest.loc[mask,"replacement"]="SMA Otsu-positive area; not arteriole identity or typed-cell percent"
        manifest.loc[mask,"reason"]="Original SMA images measured using existing approved H7 Otsu + 2 µm² component-filter method."
    manifest.to_csv(mpath,index=False)
    specpath=OUT/"analysis_specification.json";spec=json.loads(specpath.read_text())
    spec["native_figure_count"]=20;spec["figures"]=[p.stem for p in sorted(FIG.glob("*.png"))]
    spec["sma_area_replacement"]={"script":str(Path(__file__).resolve()),"donors":20,"method":"per-image Otsu +2µm² component filter, approved H7 method","arteriole_identity":False}
    spec["auc_primary_domain_um"]=[25,225];spec["auc_full_available_range_companion_retained"]=True
    specpath.write_text(json.dumps(spec,indent=2)+"\n")
    hashes={str(p.relative_to(OUT)):sha(p) for p in sorted(OUT.rglob("*")) if p.is_file() and p.name!="file_manifest.json"}
    (OUT/"file_manifest.json").write_text(json.dumps(hashes,indent=2)+"\n")
    print(pd.DataFrame(coverage).to_string(index=False),flush=True)


if __name__=="__main__":
    main()
