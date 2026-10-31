#!/usr/bin/env python3
"""Analysis-only recovery of existing structure calls with Senior union overlap.

The live app phenotype is retained. No thresholds are fitted, no served data are
changed, and no endocrine/immune chain is rerun. Only structure_calls is replayed.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
import warnings

import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts/uniform/typing"))
import classify as CL
import codex_typing as CT

SOURCE=ROOT/"data/development/uniform-20donor-20261015"
CORRECTED=ROOT/"output/h1_edge_corrected_20261024"
OUTPUT=ROOT/"output/pancfloyd_20donor_20261026/structural_cells"
AREA=ROOT/"output/pancfloyd_20donor_20261026/slides1_4/data/structure_union_area_values.csv"
CENSUS=ROOT/"output/h5_h10_original_plan_20261026/h8/full_object_census.parquet"
STEMS={"Blood vessel":"vessel","Nerve":"nerve","Lymphatic":"lymphatic"}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for data in iter(lambda:f.read(8*1024*1024),b""):h.update(data)
    return h.hexdigest()


def stream_features(path):
    """Stream a GeoJSON FeatureCollection without retaining unrelated cell outlines."""
    decoder=json.JSONDecoder();buffer="";pos=0;eof=False
    with Path(path).open() as f:
        while True:
            chunk=f.read(1024*1024)
            if not chunk:raise ValueError(f"Missing features array: {path}")
            buffer+=chunk
            match=re.search(r'"features"\s*:\s*\[',buffer)
            if match:pos=match.end();break
        while True:
            while pos<len(buffer) and buffer[pos] in " \r\n\t,":pos+=1
            if pos<len(buffer) and buffer[pos]=="]":return
            try:
                feature,end=decoder.raw_decode(buffer,pos)
            except json.JSONDecodeError:
                if eof:raise
                buffer=buffer[pos:];pos=0
                chunk=f.read(1024*1024);eof=not chunk;buffer+=chunk
                continue
            yield feature
            pos=end
            if pos>1024*1024:buffer=buffer[pos:];pos=0


def union_fraction(tree,geometries,outline):
    """Senior covered_fraction: union clipped same-class fragments / whole-cell area."""
    if outline.is_empty or outline.area<=0:return np.nan
    hits=tree.query(outline,predicate="intersects")
    if not len(hits):return 0.0
    pieces=shapely.intersection(geometries[hits],outline)
    union=pieces[0] if len(pieces)==1 else shapely.union_all(pieces)
    result=float(union.area/outline.area)
    if result < -1e-8 or result > 1+1e-7:raise ValueError(f"Union fraction out of bounds: {result}")
    return result


def process_donor(item):
    tag,entry=item;case=int(entry["case_id"])
    dest=OUTPUT/"per_cell"/f"{case}.parquet";ap=OUTPUT/"donor_audits"/f"{case}.json"
    if dest.exists() and ap.exists():return json.loads(ap.read_text())
    typ=SOURCE/"typing"
    cp=typ/"calls"/f"{tag}.parquet";jp=cp.with_suffix(".json")
    call=pd.read_parquet(cp,columns=["object_id","top_expressed_marker","cell_class"])
    saved=json.loads(jp.read_text());levels=saved["levels"]
    ctx=pd.read_parquet(typ/"context"/f"{tag}.parquet").set_index("object_id").loc[call.object_id].reset_index()
    assert ctx.qc_analysis_eligible.all() and ctx.object_id.is_unique
    wanted=set(CT.IDENTITY_MARKERS)|set(m for s in CT.STRUCTURES for m in CT.SUBPANEL[s])
    aligned=pd.read_parquet(typ/"aligned"/f"{tag}.parquet",columns=["object_id",*sorted(wanted)])
    aligned=aligned.set_index("object_id").loc[call.object_id].reset_index()
    scale=CL.load_scale(SOURCE);values=CL.scaled_markers(aligned,wanted,scale)
    top=CL.identity_argmax(values)
    assert np.array_equal(top,call.top_expressed_marker.fillna("").to_numpy()),f"{tag}: highest-marker replay differs"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        old_class,_,_,old_audit=CL.structure_calls(ctx,values,levels,top)
    for name in CT.STRUCTURES:
        for key,value in old_audit[name].items():
            assert value==saved["structures"][name][key],f"{tag} {name} {key}: {value} != cached audit"
    inside_any=ctx[[f"in_{x}" for x in CT.STRUCTURES]].any(axis=1).to_numpy()
    wanted_ids=set(ctx.loc[inside_any,"object_id"]);outlines={};geometry_path=SOURCE/"cell_geojson"/f"{tag}.geojson"
    for feature in stream_features(geometry_path):
        oid=feature.get("id")
        if oid in wanted_ids:
            outlines[oid]=shapely.make_valid(shape(feature["geometry"]))
    assert set(outlines)==wanted_ids,f"{tag}: missing {len(wanted_ids-set(outlines))} candidate cell geometries"
    corrected=ctx.copy();fraction_records={};struct_paths=[]
    for name in CT.STRUCTURES:
        object_class=CT.OBJECT_CLASS[name]
        path=SOURCE/"structures"/object_class/f"{object_class}__{tag}.geojson";struct_paths.append(path)
        geometries=np.asarray([shapely.make_valid(shape(f["geometry"])) for f in stream_features(path)],dtype=object)
        tree=shapely.STRtree(geometries)
        inside=ctx[f"in_{name}"].to_numpy(bool)
        union=np.full(len(ctx),np.nan)
        for i in np.flatnonzero(inside):union[i]=union_fraction(tree,geometries,outlines[ctx.object_id.iloc[i]])
        original=ctx[CL.CROSSED[name]].to_numpy(float)
        assert np.all(union[inside]<=original[inside]+1e-6),f"{tag}: union exceeds saved sum"
        corrected.loc[inside,CL.CROSSED[name]]=union[inside]
        fraction_records[name]=union
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        new_class,_,passed,new_audit=CL.structure_calls(corrected,values,levels,top)
    accepted_path=CORRECTED/"analysis_cells_by_donor"/f"{case}.parquet"
    accepted=pd.read_parquet(accepted_path)
    accepted=accepted[accepted.edge_analysis_eligible].copy()
    assert accepted.object_id.is_unique and accepted.object_id.isin(call.object_id).all()
    retained=ctx.object_id.isin(accepted.object_id).to_numpy()
    result=ctx[["object_id","x_um","y_um"]].copy()
    result["source_cell_class"]=call.cell_class.to_numpy()
    result["top_expressed_marker"]=top
    result["original_structure_class_before_override"]=pd.Series(old_class).fillna("").to_numpy()
    result["recovered_structure_class"]=pd.Series(new_class).fillna("").to_numpy()
    audit_rows=[]
    for name in CT.STRUCTURES:
        stem=STEMS[name];inside=ctx[f"in_{name}"].to_numpy(bool);u=fraction_records[name]
        markers=list(CT.SUBPANEL[name]);marker=markers[0]
        v=values[marker].to_numpy(float);owns=np.isin(top,CT.CLASSIFIER_MARKER[name]);clears=v>=levels[name]
        old_hit=inside&(ctx[CL.CROSSED[name]].to_numpy(float)>CT.MAJORITY_AREA)&clears&owns
        hit=inside&(u>CT.MAJORITY_AREA)&clears&owns
        result[f"centroid_inside_{stem}"]=inside
        result[f"original_summed_capped_fraction_{stem}"]=ctx[CL.CROSSED[name]].to_numpy(float)
        result[f"union_fraction_{stem}"]=u
        result[f"union_fraction_measured_{stem}"]=np.isfinite(u)
        result[f"scaled_{marker}"]=v
        result[f"cached_level_{stem}"]=levels[name]
        result[f"marker_level_clear_{stem}"]=clears
        result[f"forming_marker_top_{stem}"]=owns
        result[f"original_structure_hit_{stem}"]=old_hit
        result[f"recovered_structure_hit_{stem}"]=hit
        result[f"passed_through_{stem}"]=passed[name]
        assert np.isnan(passed[name][hit]).all()
        assert np.isnan(passed[name][~inside]).all()
        assert np.allclose(passed[name][inside&~hit],u[inside&~hit],equal_nan=True)
        audit_rows.append(dict(structure=name,cached_intermediate_hits=int(saved["structures"][name]["called"]),
            original_replay_hits_before_edge=int(old_hit.sum()),union_hits_before_edge=int(hit.sum()),
            original_replay_hits_accepted_edge=int((old_hit&retained).sum()),union_hits_accepted_edge=int((hit&retained).sum()),
            final_exclusive_hits_accepted_edge=int(((new_class==name)&retained).sum()),
            original_exclusive_hits_accepted_edge=int(((old_class==name)&retained).sum()),
            geometry_measured=int(np.isfinite(u).sum()),overlap_fraction_reduced=int((inside&(u<ctx[CL.CROSSED[name]].to_numpy(float)-1e-8)).sum()),
            cached_marker_level=float(levels[name])))
    result=result.merge(accepted.drop(columns=["x_um","y_um"]),on="object_id",how="inner",validate="one_to_one")
    result=result.rename(columns={"phenotype":"app_phenotype"})
    assert len(result)==len(accepted) and result.edge_analysis_eligible.all()
    assert result.set_index("object_id").app_phenotype.equals(accepted.set_index("object_id").phenotype.reindex(result.object_id))
    result.to_parquet(dest,index=False)
    membership_path=CORRECTED/"object_memberships_by_donor"/f"{case}.parquet"
    membership=pd.read_parquet(membership_path)
    assert not membership.duplicated(["islet_id","object_id","cell_region"]).any()
    assert membership.object_id.isin(result.object_id).all()
    member=membership.merge(result,on="object_id",validate="many_to_one",suffixes=("_membership",""))
    assert (member.phenotype==member.app_phenotype).all()
    object_rows=[];passthrough_rows=[]
    for (islet,region),g in member.groupby(["islet_id","cell_region"],observed=True):
        for name in CT.STRUCTURES:
            stem=STEMS[name]
            object_rows.append(dict(islet_id=islet,region=region,structure=name,n_recovered_cells=int(g.recovered_structure_class.eq(name).sum()),n_accepted_cells=len(g)))
            fractions=g[f"passed_through_{stem}"]
            passthrough_rows.append(dict(islet_id=islet,region=region,structure=name,n_accepted_cells=len(g),
                n_centroid_inside=int(g[f"centroid_inside_{stem}"].sum()),n_passthrough_measured=int(fractions.notna().sum()),
                passthrough_fraction_sum=float(fractions.sum(min_count=1)),
                mean_passthrough_fraction_among_measured=float(fractions.mean())))
    pd.DataFrame(object_rows).to_csv(OUTPUT/"by_donor"/f"{case}_counts.csv",index=False)
    pd.DataFrame(passthrough_rows).to_csv(OUTPUT/"by_donor"/f"{case}_passthrough.csv",index=False)
    audit=dict(case_id=case,image=tag,source_analysis_eligible_cells=len(call),accepted_edge_eligible_cells=len(accepted),
        cell_geometries_loaded=len(outlines),highest_marker_replay_matches=True,cached_audit_replay_matches=True,
        accepted_app_phenotypes_unchanged=True,exclusive_priority=list(CT.STRUCTURES),structure_rows=audit_rows,
        per_cell_output=str(dest),membership_rows=len(membership),membership_unique=True,
        source_hashes={str(p):sha(p) for p in [jp,cp,typ/"context"/f"{tag}.parquet",typ/"aligned"/f"{tag}.parquet",typ/"scale.parquet",accepted_path,membership_path,*struct_paths]},
        cell_geometry=dict(path=str(geometry_path),bytes=geometry_path.stat().st_size,mtime_ns=geometry_path.stat().st_mtime_ns))
    ap.write_text(json.dumps(audit,indent=2)+"\n")
    return audit


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--workers",type=int,default=2);args=parser.parse_args()
    for folder in ["per_cell","by_donor","donor_audits"]:(OUTPUT/folder).mkdir(parents=True,exist_ok=True)
    specification=dict(status="analysis-only recovered structural calls; served app phenotype retained",
        universe="accepted edge-eligible cell IDs from corrected analysis_cells_by_donor",structures=list(CT.STRUCTURES),
        marker_mapping=CT.CLASSIFIER_MARKER,majority_fraction_strictly_greater_than=CT.MAJORITY_AREA,
        marker_level="existing cached per-donor levels; no refit",top_marker="existing highest expressed identity marker on aligned pre-scale values, replay verified",
        overlap="Senior union of same-class object fragments clipped to each whole-cell polygon / whole-cell area",
        passthrough="union overlap fraction when centroid inside and own structure call false; NA outside or own call true",
        membership="existing corrected core/peri membership; unique cell ID per endocrine object and region",
        cell_typing_pipeline_changed=False,served_labels_changed=False,quantitative_INS_endpoint=False,
        source_script=str(Path(__file__).resolve()),source_script_sha256=sha(__file__))
    (OUTPUT/"analysis_specification.json").write_text(json.dumps(specification,indent=2)+"\n")
    entries=sorted(CT.manifest(SOURCE).items());assert len(entries)==20 and all(int(e["case_id"])!=78 for _,e in entries)
    audits=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs={pool.submit(process_donor,item):item[0] for item in entries}
        for future in as_completed(jobs):
            a=future.result();audits.append(a)
            print(a["image"],a["accepted_edge_eligible_cells"],{r["structure"]:r["union_hits_accepted_edge"] for r in a["structure_rows"]},flush=True)
    raw=pd.concat([pd.read_csv(OUTPUT/"by_donor"/f"{a['case_id']}_counts.csv") for a in audits],ignore_index=True)
    census=pd.read_parquet(CENSUS);census=census[census.object_class.isin(["Islet","EndocrineCluster"])].copy()
    areas=pd.read_csv(AREA);areas=areas[areas["class"].eq("Capillary")][["islet_id","region","region_um2"]]
    geometry=census[["islet_id","case_id","donor_status","object_class","core_area_um2","core_complete","neighborhood_complete_20"]].merge(areas,on="islet_id",validate="one_to_many")
    full=geometry.merge(pd.DataFrame({"structure":list(CT.STRUCTURES)}),how="cross")
    full=full.merge(raw,on=["islet_id","region","structure"],how="left",validate="one_to_one")
    for column in ["n_recovered_cells","n_accepted_cells"]:full[column]=full[column].fillna(0).astype(int)
    full=full.rename(columns={"region_um2":"region_area_um2"})
    full["diameter_um"]=2*np.sqrt(full.core_area_um2/np.pi)
    full["edge_eligible"]=full.core_complete&(full.region.eq("core")|full.neighborhood_complete_20)
    full["count_fraction_percent"]=100*full.n_recovered_cells/full.n_accepted_cells.replace(0,np.nan)
    full["density_cells_mm2"]=1e6*full.n_recovered_cells/full.region_area_um2.replace(0,np.nan)
    assert len(full)==4906*2*3 and not full.duplicated(["islet_id","region","structure"]).any()
    assert full.n_recovered_cells.le(full.n_accepted_cells).all()
    paired=pd.read_parquet(CORRECTED/"h1_object_features.parquet")
    for region,column in [("core","n_cells_core"),("peri","n_cells_peri")]:
        check=full[(full.region==region)&(full.structure=="Blood vessel")].merge(paired[["islet_id",column]],on="islet_id")
        assert np.array_equal(check.n_accepted_cells,check[column]),f"{region}: accepted count denominator changed"
    full.to_csv(OUTPUT/"object_structure_counts.csv",index=False);full.to_parquet(OUTPUT/"object_structure_counts.parquet",index=False)
    passthrough=pd.concat([pd.read_csv(OUTPUT/"by_donor"/f"{a['case_id']}_passthrough.csv") for a in audits],ignore_index=True)
    passthrough=full.drop(columns=["n_recovered_cells","n_accepted_cells","count_fraction_percent","density_cells_mm2"]).merge(passthrough,on=["islet_id","region","structure"],how="left",validate="one_to_one")
    for c in ["n_accepted_cells","n_centroid_inside","n_passthrough_measured"]:passthrough[c]=passthrough[c].fillna(0).astype(int)
    passthrough.to_csv(OUTPUT/"object_passthrough_summaries.csv",index=False)
    audit_table=pd.DataFrame([dict(case_id=a["case_id"],image=a["image"],**r) for a in audits for r in a["structure_rows"]])
    audit_table.to_csv(OUTPUT/"structure_call_audit.csv",index=False)
    (OUTPUT/"provenance.json").write_text(json.dumps(dict(donors=20,geometric_objects=4906,object_region_structure_rows=len(full),
        accepted_cells=sum(a["accepted_edge_eligible_cells"] for a in audits),saved_intermediate_replay_matches=True,
        accepted_cell_denominators_match_corrected_features=True,area_source=str(AREA),area_sha256=sha(AREA),
        census_sha256=sha(CENSUS),donor_audits=sorted(audits,key=lambda x:x["case_id"])),indent=2)+"\n")
    print(audit_table.groupby("structure")[["cached_intermediate_hits","union_hits_before_edge","union_hits_accepted_edge"]].sum().to_string(),flush=True)


if __name__=="__main__":main()
