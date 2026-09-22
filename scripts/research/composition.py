#!/usr/bin/env python3
"""Versioned post-typing composition in classifier-defined endocrine cores.

The endocrine pixel classifier defines the islet core and therefore the
per-endocrine-cell denominator.  Downstream typing supplies phenotype
numerators and uncertainty, but cannot revoke classifier-level islet
eligibility.  Positivity and continuous expression retain independent
availability; no missing type or marker call is converted to a negative.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from common import MARKERS, sha256, write_json

KEY = ['case_id', 'object_id']
ENDOCRINE_COMPONENTS = ['Alpha cell', 'Beta cell', 'Delta cell', 'Hormone-unresolved endocrine cell']


def check_unique(frame, keys, label):
    if not set(keys).issubset(frame) or frame.duplicated(keys).any():
        raise ValueError(f'{label} requires unique {keys}')


def canonical_cells(typing, corrected,accepted_statuses=('validated',)):
    typing=typing.copy(); corrected=corrected.copy()
    for f in [typing,corrected]:
        f['case_id']=f.case_id.astype(str);f['object_id']=f.object_id.astype(str)
    check_unique(typing,KEY,'Typing');check_unique(corrected,KEY,'Measurements')
    if set(map(tuple,typing[KEY].to_numpy())) != set(map(tuple,corrected[KEY].to_numpy())):
        raise ValueError('Typing and corrected measurements do not contain exactly the same raw objects')
    required={'base_type','classification_status','projection_supported','is_endocrine'}
    if not required.issubset(typing):raise ValueError(f'Missing typing fields: {sorted(required-set(typing))}')
    if not set(MARKERS).issubset(corrected):raise ValueError('Corrected input must contain every panel marker')
    if any(m in typing for m in MARKERS):raise ValueError('Typing contains ambiguous duplicate intensity columns')
    cells=typing.merge(corrected,on=KEY,validate='one_to_one')
    values=cells[MARKERS].to_numpy(float)
    cells['measured_raw']=np.isfinite(values).all(axis=1)&(values>=0).all(axis=1)
    quality_excluded=cells.classification_status.eq('segmentation_quality_excluded')
    for column in ['segmentation_quality_excluded','segmentation_quality_excluded_x','segmentation_quality_excluded_y']:
        if column in cells:quality_excluded=quality_excluded|cells[column].fillna(False).astype(bool)
    cells['segmentation_quality_excluded']=quality_excluded
    cells['measured']=cells.measured_raw&~quality_excluded
    cells['typed']=cells.classification_status.isin(accepted_statuses)&cells.projection_supported.fillna(False).astype(bool)&cells.measured
    # This explicit validated flag is authoritative; labels alone cannot create E.
    endocrine=cells.is_endocrine.astype('boolean')
    if ((~cells.typed)&endocrine.notna()).any():
        raise ValueError('Unvalidated identity has a non-null endocrine flag')
    cells['endocrine']=cells.typed&endocrine.fillna(False)
    if cells.loc[cells.endocrine,'base_type'].isin(ENDOCRINE_COMPONENTS).eq(False).any():
        raise ValueError('Endocrine flag outside the registered expanded union')
    return cells


def aggregate_donor(islets,memberships,cells,base_types,state_names=(),include_marker_summaries=True):
    """Return feature, count and long marker-summary tables for one donor."""
    islets=islets.copy();memberships=memberships.copy()
    for f in [islets,memberships]:f['case_id']=f.case_id.astype(str)
    memberships['object_id']=memberships.object_id.astype(str)
    check_unique(islets,['islet_id'],'Islets')
    required={'area_um2','geometry_valid','core_valid_area_um2','peri_valid_area_um2',
              'classifier_defined_core','image_edge_contact','image_edge_clearance_um'}
    if not required.issubset(islets):raise ValueError(f'Missing geometry fields {sorted(required-set(islets))}')
    if not islets.classifier_defined_core.fillna(False).astype(bool).all():
        raise ValueError('Every islet must have verified endocrine pixel-classifier provenance')
    if not (islets.geometry_valid.fillna(False).astype(bool)
            & np.isfinite(islets.core_valid_area_um2) & islets.core_valid_area_um2.gt(0)).all():
        raise ValueError('Classifier cores must have valid positive geometry before composition')
    if not memberships.compartment.isin(['core','peri']).all():raise ValueError('Unknown compartment')
    if 'band_um' in memberships:
        memberships=memberships.loc[memberships.band_um.eq(20)].copy()
    check_unique(memberships,KEY+['islet_id','compartment'],'Membership')
    donor_set=set(islets.case_id)
    if set(cells.case_id)!=donor_set or set(memberships.case_id)-donor_set:
        raise ValueError('Mixed donor membership or measurements')
    joined=memberships.merge(cells,on=KEY,how='left',validate='many_to_one',indicator=True)
    if not joined._merge.eq('both').all():raise ValueError('Membership points outside canonical cells')
    if set(joined.islet_id)-set(islets.islet_id):raise ValueError('Membership points outside geometry')
    lookup=joined.groupby(['islet_id','compartment'],sort=False)
    empty=joined.iloc[:0]
    positive_supported={m:('positive__'+m in cells and cells['positive__'+m].notna().any()) for m in MARKERS}
    features=[];counts=[];summaries=[]
    for geom in islets.to_dict('records'):
        iid=geom['islet_id'];cid=geom['case_id']
        core=lookup.get_group((iid,'core')) if (iid,'core') in lookup.groups else empty
        valid_core=bool(geom['geometry_valid']) and np.isfinite(geom['core_valid_area_um2']) and geom['core_valid_area_um2']>0
        # The pixel classifier already established this as endocrine core.  E
        # is therefore the unique canonical core-cell census, independent of a
        # downstream phenotype label or marker-gate support decision.
        denominator=float(len(core))
        f=dict(geom,endocrine_count=denominator,measured_cells=int(core.measured.sum()),
               raw_cells=len(core),typed_cells=int(core.typed.sum()),unresolved_cells=int((~core.typed&core.measured).sum()),excluded_cells=int((~core.measured).sum()),
               fraction_typed=float(core.typed.sum()/core.measured.sum()) if core.measured.any() else np.nan)
        reasons=[]
        if bool(geom['image_edge_contact']):reasons.append('image_edge_contact')
        if len(core)<1:reasons.append('no_core_cells_in_classifier_area')
        for component in ENDOCRINE_COMPONENTS:
            f['endocrine_component__'+component]=int((core.endocrine&core.base_type.eq(component)).sum()) if valid_core else np.nan
        for comp in ['core','peri']:
            d=lookup.get_group((iid,comp)) if (iid,comp) in lookup.groups else empty
            valid=bool(geom['geometry_valid']) and np.isfinite(geom[f'{comp}_valid_area_um2']) and geom[f'{comp}_valid_area_um2']>0
            count_available=valid and (d.typed.any() or len(d)==0)
            # Peri-islet counts are retained with their availability flag.
            # Missing peri typing must not erase a donor from the primary
            # core-composition trajectory; spatial analyses apply their own
            # compartment availability gates.
            nt=int(d.typed.sum()); nm=int(d.measured.sum());nu=int((~d.typed&d.measured).sum());ne=int((~d.measured).sum())
            for name in base_types:
                count=float((d.typed&d.base_type.eq(name)).sum()) if count_available else np.nan
                raw=count/denominator if np.isfinite(denominator) and denominator>0 else np.nan
                plotted=float(np.log1p(raw)) if np.isfinite(raw) else np.nan
                f[f'{comp}__count__{name}']=plotted
                counts.append(dict(case_id=cid,islet_id=iid,compartment=comp,phenotype=name,count=count,
                                   available=count_available,n_typed=nt,n_measured=nm,unresolved_cells=nu,excluded_cells=ne,
                                   fraction_typed=nt/nm if nm else np.nan,endocrine_count=denominator,
                                   raw_ratio=raw,plotted_value=plotted))
            state_column='phenotype_state' if 'phenotype_state' in d else 'state'
            for name in state_names:
                states_available='phenotype_states' in d or state_column in d
                if 'phenotype_states' in d:
                    matches=d.phenotype_states.map(lambda values: name in values if isinstance(values,(list,tuple,np.ndarray)) else False)
                else:matches=d[state_column].eq(name) if state_column in d else pd.Series(False,index=d.index)
                count=float((d.typed&matches).sum()) if count_available and states_available else np.nan
                raw=count/denominator if np.isfinite(denominator) and denominator>0 else np.nan
                counts.append(dict(case_id=cid,islet_id=iid,compartment=comp,phenotype=name,count=count,
                                   available=count_available and states_available,n_typed=nt,n_measured=nm,unresolved_cells=nu,excluded_cells=ne,
                                   fraction_typed=nt/nm if nm else np.nan,endocrine_count=denominator,
                                   raw_ratio=raw,plotted_value=float(np.log1p(raw)) if np.isfinite(raw) else np.nan))
            if not include_marker_summaries:continue
            populations={'All cells':d.loc[d.measured], 'Endocrine union':d.loc[d.endocrine]}
            populations.update({name:d.loc[d.typed&d.base_type.eq(name)] for name in base_types})
            for population,pop in populations.items():
                means=np.log1p(pop[MARKERS].to_numpy(float)).mean(axis=0) if valid and len(pop) else np.full(len(MARKERS),np.nan)
                for marker,mean in zip(MARKERS,means):
                    column='positive__'+marker
                    available=valid and (len(d)==0 or nm>0) and positive_supported[marker] and pop[column].notna().all()
                    # Empty measured regions have known zero positives only when
                    # a reviewed call column exists; missing gates stay missing.
                    positive=float(pop[column].sum()) if available else np.nan
                    summaries.append(dict(case_id=cid,islet_id=iid,compartment=comp,population=population,
                                          marker=marker,n_cells=len(pop),excluded_cells=ne,unresolved_cells=nu,mean_log1p=float(mean),
                                          positive_count=positive,positive_available=available,endocrine_count=denominator,
                                          positive_ratio=positive/denominator if np.isfinite(denominator) and denominator>0 else np.nan))
        primary=[f'core__count__{name}' for name in base_types]
        f['eligible']=not reasons
        f['trajectory_feature_supported']=bool(f['eligible'] and np.isfinite([f[name] for name in primary]).all())
        f['exclusion_reason']=';'.join(dict.fromkeys(reasons));features.append(f)
    return pd.DataFrame(features),pd.DataFrame(counts),pd.DataFrame(summaries)


def process_donor(cid,islets,memberships,root,out,base_types,state_names):
    tp=root/'typing/cells'/f'{cid}.parquet';cp=root/'corrected'/f'{cid}-cells.parquet'
    cells=canonical_cells(pd.read_parquet(tp),pd.read_parquet(cp))
    f,c,summary=aggregate_donor(islets,memberships,cells,base_types,state_names)
    piece=out/'analysis/composition_parts'/f'{cid}-marker_summaries.parquet'
    piece.parent.mkdir(parents=True,exist_ok=True);summary.to_parquet(piece,index=False)
    print(f'{cid}: {len(f)} islets, {f.eligible.sum()} eligible; {cells.typed.sum()}/{len(cells)} validated cells',flush=True)
    return f,c,str(piece),{str(tp):sha256(tp),str(cp):sha256(cp)}


def concatenate_parquet(paths,target):
    import pyarrow.parquet as pq
    writer=None
    try:
        for path in paths:
            table=pq.read_table(path)
            if writer is None:writer=pq.ParquetWriter(target,table.schema,compression='zstd')
            writer.write_table(table)
    finally:
        if writer is not None:writer.close()


def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True);p.add_argument('--output-root',required=True)
    p.add_argument('--memberships');p.add_argument('--islets');p.add_argument('--catalog');p.add_argument('--workers',type=int,default=1);a=p.parse_args()
    root=Path(a.input_root);out=Path(a.output_root)
    if (out/'trajectory/coordinates.csv').exists():raise ValueError('Existing fitted trajectory is immutable; reprocess into a fresh candidate root')
    catalog_path=Path(a.catalog) if a.catalog else root/'typing/catalog.json'
    catalog=json.loads(catalog_path.read_text());base_types=[t['name'] if isinstance(t,dict) else t for t in catalog['base_types']]
    state_names=[t['name'] if isinstance(t,dict) else t for t in catalog.get('states',[])]
    if not base_types or len(set(base_types))!=len(base_types):raise ValueError('Empty or duplicate base type catalog')
    islet_path=Path(a.islets) if a.islets else root/'analysis/islets.parquet'
    member_path=Path(a.memberships) if a.memberships else root/'analysis/memberships.parquet'
    islets=pd.read_parquet(islet_path);islets['case_id']=islets.case_id.astype(str)
    eligibility_path=root/'analysis/classifier_eligibility_manifest.json'
    if not eligibility_path.is_file():raise ValueError('Classifier eligibility provenance is required')
    eligibility=json.loads(eligibility_path.read_text())
    if eligibility.get('islets_sha256')!=sha256(islet_path):raise ValueError('Classifier eligibility islet table changed')
    if eligibility.get('definition',{}).get('downstream_typing_can_revoke_eligibility') is not False:
        raise ValueError('Downstream typing must not define classifier islet eligibility')
    cohort=pd.read_csv(root/'cohort.csv',dtype={'case_id':str})
    if cohort.case_id.duplicated().any() or set(islets.case_id)!=set(cohort.case_id):raise ValueError('Islet geometry does not account for exactly the manifest donor cohort')
    memberships=pd.read_parquet(member_path);memberships['case_id']=memberships.case_id.astype(str)
    (out/'trajectory').mkdir(parents=True,exist_ok=True);(out/'analysis').mkdir(parents=True,exist_ok=True)
    from joblib import Parallel,delayed
    sources={str(p):sha256(p) for p in [catalog_path,islet_path,member_path,eligibility_path]}
    results=Parallel(n_jobs=a.workers)(delayed(process_donor)(cid,g,memberships.loc[memberships.case_id.eq(cid)],root,out,base_types,state_names) for cid,g in islets.groupby('case_id',sort=True))
    frame=pd.concat([r[0] for r in results],ignore_index=True);count=pd.concat([r[1] for r in results],ignore_index=True)
    for result in results:sources.update(result[3])
    concatenate_parquet([r[2] for r in results],out/'analysis/marker_summaries.parquet')
    import shutil
    shutil.copyfile(out/'analysis/marker_summaries.parquet',out/'trajectory/marker_summaries.parquet')
    frame.to_parquet(out/'trajectory/composition_features.parquet',index=False);frame.to_csv(out/'trajectory/composition_features.csv',index=False)
    count.to_parquet(out/'analysis/counts.parquet',index=False)
    count.assign(validated=count.available).to_parquet(out/'trajectory/counts.parquet',index=False)
    count.assign(validated=count.available).to_csv(out/'trajectory/counts.csv',index=False)
    endocrine=frame[['islet_id','case_id','endocrine_count','raw_cells','typed_cells','unresolved_cells','excluded_cells','fraction_typed']+[c for c in frame if c.startswith('endocrine_component__')]].copy()
    endocrine['endocrine_validated']=endocrine.endocrine_count.notna()
    endocrine['validated']=endocrine.endocrine_validated
    endocrine.to_parquet(out/'trajectory/endocrine_counts.parquet',index=False)
    endocrine.to_csv(out/'trajectory/endocrine_counts.csv',index=False)
    write_json(out/'trajectory/count_manifest.json',dict(schema_version=4,
               definition='validated identity numerators; unresolved cells counted separately; classifier-core cell denominator',sources=sources,
               endocrine_components=ENDOCRINE_COMPONENTS,
               endocrine_denominator='all canonical core cells inside the endocrine pixel-classifier area, counted once',
               downstream_typing_can_revoke_islet_eligibility=False,
               counts_sha256=sha256(out/'trajectory/counts.parquet')))
    write_json(out/'trajectory/marker_panel.json',dict(schema_version=4,markers=MARKERS,marker_summary_path='trajectory/marker_summaries.parquet',
               default_population='Endocrine union',mean_definition='mean log1p corrected cell intensity within selected typed population',
               positivity_definition='reviewed marker-positive cell count per classifier-defined core cell count'))
    state_columns=[]
    if state_names:
        state_rows=count.loc[count.phenotype.isin(state_names)].copy()
        state_rows['column']=state_rows.compartment+'__count__state::'+state_rows.phenotype
        state_table=state_rows.pivot(index='islet_id',columns='column',values='plotted_value').reset_index()
        state_columns=[c for c in state_table if c!='islet_id']
        state_frame=frame.merge(state_table,on='islet_id',validate='one_to_one')
        state_frame['state_feature_supported']=state_frame.eligible & np.isfinite(state_frame[state_columns]).all(axis=1)
        state_frame.to_parquet(out/'trajectory/state_composition_features.parquet',index=False)
    # The counts and all-marker heatmap tables must stay separate from the fit.
    # The primary trajectory uses the fully observed core subset. Peri-islet
    # ratios remain declared and exported for display and secondary analyses,
    # but incomplete peri typing cannot remove a donor from the primary fit.
    columns=[f'{comp}__count__{t}' for comp in ['core','peri'] for t in base_types]
    primary_columns=[f'core__count__{t}' for t in base_types]
    exclusion_counts={}
    for value in frame.loc[~frame.eligible,'exclusion_reason'].fillna(''):
        for reason in {item for item in str(value).split(';') if item}:
            exclusion_counts[reason]=exclusion_counts.get(reason,0)+1
    eligible_by_donor={str(case_id):int(count) for case_id,count in
        frame.loc[frame.eligible].groupby('case_id').size().items()}
    write_json(out/'trajectory/feature_manifest.json',dict(schema_version=4,representation='typed_composition',
               columns=columns,primary_columns=primary_columns,state_sensitivity_columns=columns+state_columns if state_columns else [],
               state_features_sha256=sha256(out/'trajectory/state_composition_features.parquet') if state_columns else None,base_types=base_types,endocrine_components=ENDOCRINE_COMPONENTS,
               formula='log1p(validated base-type count / classifier-defined core cell count)',
               minimum_classifier_core_cells=1,
               eligibility=dict(all_of=[
                   'the endocrine pixel-classifier core has no image-edge contact',
                   'at least one canonical cell centroid is inside the classifier-defined core',
               ],classifier_core_cell_definition='unique canonical cell centroid inside the endocrine pixel-classifier core',
                  endocrine_denominator_definition='all canonical classifier-core cells, counted once',
                  downstream_typing_can_revoke_eligibility=False),
               exclusion_counts=exclusion_counts,eligible_islets_by_donor=eligible_by_donor,
               feature_supported_islets=int((frame.eligible&frame.trajectory_feature_supported).sum()),
               uncertainty='Unvalidated cells remain in the classifier-core denominator, are excluded from phenotype numerators, and are reported as unresolved; typing support affects coordinates, not islet eligibility.',
               sources=sources,features_sha256=sha256(out/'trajectory/composition_features.parquet'),
               markers=MARKERS,typing_required=True,eligible_islets=int(frame.eligible.sum()),donors=sorted(frame.case_id.unique()),
               unsupported_donors=sorted(set(frame.case_id)-set(frame.loc[frame.eligible,'case_id']))))

if __name__=='__main__':main()
