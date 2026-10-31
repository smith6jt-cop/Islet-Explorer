#!/usr/bin/env python3
"""Native Prism graph-ready data and exact CSV sources for the current app run.

Reuses the existing validated PZFX writer; no analyses or graph layouts are fitted.
"""
from pathlib import Path
import hashlib,json,re,shutil,sys
import numpy as np
import pandas as pd
from lxml import etree as ET
import export_prism_20261026 as writer

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'output/app_hypotheses_20261026'
OUT=RUN/'prism'
GROUPS=['ND','Aab+','T1D']
SOURCES={};REGISTERED=[]


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def key(s):
    clean=re.sub('[^A-Za-z0-9_]+','_',str(s)).strip('_')
    return clean if len(clean)<95 else clean[:84]+'_'+hashlib.sha256(clean.encode()).hexdigest()[:9]
def read(p):return pd.read_csv(p,dtype={'case_id':str,'islet_id':str},low_memory=False)
def source(p):
    p=Path(p);SOURCES[str(p.relative_to(RUN))]=sha(p)
    dest=OUT/'source_csv'/p.relative_to(RUN);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)


def register(k,title,frame,numeric,p,desc,x=None,ci=False,paired=False):
    if not len(frame) or not numeric:return
    k=key(k);assert k not in REGISTERED,k
    a=frame[numeric].apply(pd.to_numeric,errors='raise').to_numpy(float)
    assert not np.isinf(a).any(),k
    if ci:
        assert len(numeric)==3
        v=frame[numeric].dropna()
        assert (v.iloc[:,2]<=v.iloc[:,0]).all() and (v.iloc[:,0]<=v.iloc[:,1]).all(),k
    writer.register(k,title,frame,numeric,desc,p,paired=paired,x=x,ci=ci)
    REGISTERED.append(k);source(p)


def group_table(p,d,value,title,suffix=''):
    assert d.case_id.is_unique,(p,suffix)
    d=d.sort_values(['donor_status','case_id']);z=pd.DataFrame({'Donor_ID':d.case_id.to_numpy()})
    for g in GROUPS:z[g]=np.where(d.donor_status.eq(g),d[value],np.nan)
    register('GROUP_'+str(p.relative_to(RUN))+'_'+suffix,title,z,GROUPS,p,
        'One observation per donor in its disease column; groups are unpaired. Blank cells are not zero. Values and units are those of '+value+'. See the supplied figure legend and separate statistics; no test is recomputed.')


def primary_graphs():
    base=RUN/'h1_h4/statistics'
    for h,col in [('h1','standardized_beta'),('h2','standardized_beta'),('h3','difference_density_mm2')]:
        p=base/f'{h}_donor_values.csv';group_table(p,read(p),col,h.upper()+' donor comparison')
    p=RUN/'h5/data/donor_measurements.csv';d=read(p)
    for endpoint,z in d.groupby('endpoint'):group_table(p,z,'value','Endocrine quantity: '+endpoint,endpoint)
    p=RUN/'h9_h10/h9/data/donor_slopes.csv';group_table(p,read(p),'coverage_core_immune_slope','Coverage–immune-density donor coefficients')
    p=RUN/'h9_h10/h10/data/fixed_score_and_clinical_measurements.csv';d=read(p)
    group_table(p,d,'spatial_score_all20_display','Fixed spatial score by disease')
    for label in ['clinical','clinical_disease_conditioned']:
        p=RUN/f'h9_h10/h10/data/{label}_held_out_predictions.csv';d=read(p)
        for name in ['baseline','augmented']:
            z=d[['case_id','observed','baseline_prediction','augmented_prediction']].copy()
            z.insert(0,'Donor',d.donor_status+' | '+d.case_id)
            register(label+'_'+name+'_XY',label.replace('_',' ')+' — '+name,z[['Donor','observed',name+'_prediction']],
                [name+'_prediction'],p,'One measured donor per row; X=observed ln C-peptide, Y=held-out predicted ln C-peptide. Censored donor is absent; no imputed concentration. Scaling and fits use training donors only.',x='observed')
    p=RUN/'h9_h10/h10/data/all_six_T1D_observations_and_predictions.csv';d=read(p)
    z=pd.DataFrame({'Donor_ID':d.case_id,'Measured_C_peptide_ng_mL':d.C_peptide_ng_mL.where(~d.C_peptide_censored),
        'Censored_upper_bound_ng_mL':d.C_peptide_ng_mL.where(d.C_peptide_censored),
        'Baseline_predicted_median_ng_mL':d.baseline_predicted_median_ng_mL,'Augmented_predicted_median_ng_mL':d.augmented_predicted_median_ng_mL})
    register('T1D_six_measured_and_censored','Six T1D donors: measured values, bound and predictions',z,list(z)[1:],p,
        'Rows are all six T1D donors. The measured column is blank for6550; its0.02 limit is in the separate upper-bound column, not a measured value. Predicted medians are exponentiated log predictions. Five measured donors use LOO;6550 trains on19 measured donors. Concentration units ng/mL. No squared error is assigned to6550.',paired=True)
    # Matching rows preserve each donor's paired spatial compartments.
    for p,cat,val in [(RUN/'h1_h4/data/crossing_compartment_plot_values.csv','compartment','density_cells_mm2')]:
        d=read(p)
        for g,z in d.groupby('donor_status'):
            w=z.pivot(index='case_id',columns=cat,values=val);w.index.name='Donor_ID';w=w.reset_index()
            register('PAIRED_'+p.stem+'_'+g,'Boundary compartments — '+g,w,list(w)[1:],p,'Within-donor paired compartment densities in cells/mm². Row identity defines the pairing. No cross-disease pairing.',paired=True)
    for p in (RUN/'h6_h8').rglob('*donor*.csv'):
        if 'extra_by_donor' in str(p):continue
        d=read(p)
        if not {'case_id','donor_status'}<=set(d):continue
        val=next((c for c in ['estimate','standardized_beta','value','adjacency_excess'] if c in d),None)
        if val:
            split=[c for c in ['model','endpoint','hypothesis','variant','field','width_um'] if c in d]
            for label,z in d.groupby(split,dropna=False) if split else [('all',d)]:
                if z.case_id.is_unique:group_table(p,z,val,p.stem+' '+str(label),str(label))
    p=RUN/'h6_h8/data/h6_observed_permuted_donor.csv'
    if p.exists():
        d=read(p)
        for g,q in d.groupby('donor_status'):
            z=q[['case_id','permuted','observed']].rename(columns={'case_id':'Donor_ID'})
            register('PAIRED_h6_observed_permuted_'+g,'Observed and permuted adjacency — '+g,z,['permuted','observed'],p,
                'One donor per row; arithmetic object means of permuted and observed adjacency percentages. The same donor is paired across columns. No INS measurement and no new paired test.',paired=True)


def representation_tables(paths):
    for p in paths:
        if p.stem=='h8_census_with_states':
            d=read(p);identifier='object_id' if 'object_id' in d else 'islet_id'
            d.insert(0,'Object',d.donor_status+' | '+d.case_id+' | '+d[identifier].astype(str))
            values=['centroid_y_um','leiden_0.5']
            register('XY_h8_full_census_states','Full endocrine census: physical coordinates and saved state',
                d[['Object','centroid_x_um',*values]],values,p,
                'X and first Y are physical centroid coordinates in micrometres. Saved Leiden0.5 cluster is color metadata; missing cluster is blank, not zero or absent endocrine cells. Includes all census objects and donor IDs. This is a display table, not independent donor replicates.',x='centroid_x_um')
        if not any(t in p.stem for t in ['pseudotime','cluster_occupancy','cluster_heatmap','state_objects','app_representation_overview']):continue
        d=read(p);tag=str(p.relative_to(RUN))
        if {'pt_midpoint','mean','feature','donor_status'}<=set(d):
            for feature,z in d.groupby('feature'):
                cols=['donor_status','case_id'] if 'case_id' in z else ['donor_status']
                z=z.copy();z['Series']=z[cols].astype(str).agg(' | '.join,axis=1)
                w=z.pivot(index='pt_midpoint',columns='Series',values='mean').sort_index();w=w.reset_index();w.insert(0,'Bin',np.arange(1,len(w)+1))
                register('XY_'+tag+'_'+feature,p.stem+' — '+feature,w,list(w)[2:],p,'X is the fixed app rank-pseudotime bin midpoint. Y columns preserve the source object-bin or donor-bin arithmetic means. Missing bins are blank, not zero. This is a display profile; no test or refit is performed.',x='pt_midpoint')
        if 'donor_pseudotime_core_peri_profiles' in p.stem:
            for f in [c for c in d if c.endswith('_mean')]:
                z=d.copy();z['Series']=z.donor_status+' | '+z.case_id
                w=z.pivot(index='pt_midpoint',columns='Series',values=f).sort_index().reset_index();w.insert(0,'Bin',np.arange(1,len(w)+1))
                register('XY_'+tag+'_'+f,'Donor pseudotime profile — '+f,w,list(w)[2:],p,'One Y series per donor; arithmetic object mean within a fixed app pseudotime bin. Missing is blank. Source counts and object identities are supplied separately.',x='pt_midpoint')
        if {'cluster','occupancy_percent','case_id','donor_status'}<=set(d):
            w=d.pivot(index=['donor_status','case_id'],columns='cluster',values='occupancy_percent');w.columns=['Cluster '+str(c) for c in w]
            w=w.reset_index();w.insert(0,'Donor',w.donor_status+' | '+w.case_id)
            register('HEATMAP_'+tag,'Donor composition-cluster occupancy',w[['Donor']+list(w.columns[3:])],list(w.columns[3:]),p,'One donor per row; columns are saved Leiden0.5 clusters. Values are percentages of that donor’s complete embedded objects; absent clusters are measured zeros. Rows are not object replicates.')
        if {'feature','donor_status','leiden_0.5','equal_donor_mean'}<=set(d):
            for f,z in d.groupby('feature'):
                w=z.pivot(index='donor_status',columns='leiden_0.5',values='equal_donor_mean');w.columns=['Cluster '+str(c) for c in w];w=w.reset_index()
                register('HEATMAP_'+tag+'_'+f,'Cluster means — '+f,w,list(w)[1:],p,'Each heatmap cell is the source equal-donor mean within disease and saved cluster. Missing disease/cluster measurements remain blank; values are not independent replicate observations.')
        if {'islet_id','umap_1','umap_2'}<=set(d):
            d=d.copy();d.insert(0,'Object',d.donor_status+' | '+d.case_id+' | '+d.islet_id)
            numeric=[c for c in d.select_dtypes(include='number') if c not in ['umap_1','pt_bin','leiden_0.5']]
            register('XY_'+tag,p.stem+' — saved coordinates and coloring values',d[['Object','umap_1',*numeric]],numeric,p,'X=existing UMAP1; UMAP2 and feature-color measurements are separate numeric columns. Rows retain donor and object identity. No embedding is fitted. Feature columns are color metadata, not additional coordinates or clinical replicates.',x='umap_1')


def exact_tables(paths):
    for p in paths:
        source(p);d=read(p)
        if not len(d):continue
        assert not any(re.search(r'(^|_)INS(_|$)|ins_mean|mean_beta_INS',str(c)) for c in d),str(p)
        # Full high-dimensional analysis inputs are preserved byte-for-byte in source_csv.
        if len(d.columns)>55 or len(d)>12000:continue
        if 'extra_by_donor' in str(p):continue
        labels=[c for c in ['hypothesis','model','analysis','endpoint','feature','donor_status','case_id','islet_id','comparison','contrast','term','group','width_um','compartment','pt_bin','cluster','leiden_0.5'] if c in d]
        row=d[labels].fillna('').astype(str).agg(' | '.join,axis=1) if labels else pd.Series(np.arange(1,len(d)+1).astype(str))
        numeric=[c for c in d.select_dtypes(include='number') if c not in ['case_id','npod_case_id']]
        if not numeric:continue
        z=d[numeric].copy();z.insert(0,'Row_ID',row)
        if {'C_peptide_ng_mL','C_peptide_censored'}<=set(d):
            cens=d.C_peptide_censored.astype(str).str.lower().eq('true')
            z['Measured_C_peptide_ng_mL']=d.C_peptide_ng_mL.where(~cens)
            z['Censored_upper_bound_ng_mL']=d.C_peptide_ng_mL.where(cens)
            z=z.drop(columns=['C_peptide_ng_mL']);numeric=[c for c in numeric if c!='C_peptide_ng_mL']+['Measured_C_peptide_ng_mL','Censored_upper_bound_ng_mL']
        tag=str(p.relative_to(RUN))
        register('EXACT_'+tag,p.stem+' — exact numeric table',z,numeric,p,
            'Exact numeric columns from '+tag+'. Row titles concatenate source identifiers; all categorical fields and precise original CSV are in source_csv. This is a source/statistics table, not a new replicate design. Model estimates, fitted values, limits, p values and profile means must not be treated as independent donor measurements. C-peptide limits, when present, are kept separate from measured values.')
        triples=[('estimate','lower_95','upper_95'),('estimate','lower','upper'),('mean_beta','lower_95','upper_95'),('source_z.trend','lower.CL','upper.CL'),('exposure_z.trend','lower.CL','upper.CL'),('coverage_z.trend','lower.CL','upper.CL'),('contact_z.trend','lower.CL','upper.CL')]
        for e,l,u in triples:
            if {e,l,u}<=set(d):
                ok=d[[e,l,u]].notna().all(axis=1);q=d.loc[ok]
                if len(q):
                    frame=pd.DataFrame({'Label':row[ok].to_numpy(),'X':np.arange(1,len(q)+1),'Estimate':q[e].to_numpy(),'Upper_95CI':q[u].to_numpy(),'Lower_95CI':q[l].to_numpy()})
                    register('CI_'+tag,p.stem+' — estimates and confidence limits',frame,['Estimate','Upper_95CI','Lower_95CI'],p,
                        'Estimates with actual upper/lower95% confidence limits in Prism order: estimate, upper, lower. Limits are not error widths and are not replicate observations. See source model labels and exact statistics.',x='X',ci=True)
                break


def figure_sources():
    mappings={}
    for family,filename in [('h1_h4','manifest.json'),('h6_h8','figure_manifest.json'),('h9_h10','figure_manifest.json')]:
        p=RUN/family/filename
        if not p.exists():continue
        doc=json.loads(p.read_text());rows=doc.get('figures',[]) if isinstance(doc,dict) else doc
        for row in rows:
            vals=row.get('data',row.get('data_csv',row.get('source_csv',[])))
            src=[]
            for v in vals:
                q=Path(v);q=q if q.is_absolute() else RUN/family/q
                assert q.is_file(),q
                src.append(str(q.relative_to(RUN)))
            mappings[(family,row['stem'])]=src
        dest=OUT/'manifests'/family/filename;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
        SOURCES[str(p.relative_to(RUN))]=sha(p)
    for stem in ['cell_quantity','cluster_occupancy','endocrine_composition']:
        mappings[('h5',stem)]=['h5/data/'+stem+'_points.csv','h5/data/donor_measurements.csv','h5/omnibus_tests.csv','h5/pairwise_tests.csv']
    for stem in ['endocrine_states_states','endocrine_states_clusters']:
        mappings[('h5',stem)]=['h5/data/endocrine_states_'+s+'.csv' for s in ['state_objects','pseudotime_profiles','donor_pseudotime','donor_cluster_values']]
    mappings[('h5','beta_presence_tissue_maps')]=['h5/data/beta_presence_tissue_maps.csv']
    mappings[('shared','app_representation_overview')]=['shared/data/app_representation_overview.csv']
    return mappings


def main():
    for sub in ['csv','txt','source_csv','statistics','plots','legends','support']:(OUT/sub).mkdir(parents=True,exist_ok=True)
    writer.OUT=OUT;writer.TABLES=[];writer.FIGURES=[];writer.INS_NOTICE=OUT/'.unused_historical_notice'
    shutil.copy2(ROOT/'output/prism_export_20261026/support/Prism8XMLSchema.xsd',OUT/'support/Prism8XMLSchema.xsd')
    paths=sorted(p for p in RUN.rglob('*.csv') if 'prism' not in p.relative_to(RUN).parts and 'h7_extra_by_donor' not in p.parts)
    primary_graphs();representation_tables(paths);exact_tables(paths)
    figrows=[];mapping=figure_sources()
    for p in sorted(RUN.rglob('*.png')):
        if 'prism' in p.relative_to(RUN).parts or p.parent.name!='figures':continue
        stem=key(str(p.relative_to(RUN).with_suffix('')));legend=p.with_name(p.stem+'_legend.txt')
        assert legend.is_file(),legend
        writer.figure(stem,p,legend)
        family=p.relative_to(RUN).parts[0];datafiles=mapping.get((family,p.stem),[])
        native=[t['key'] for t in writer.TABLES if str(Path(t['source']).relative_to(RUN.relative_to(ROOT))) in datafiles]
        figrows.append(dict(figure=stem,source=str(p.relative_to(RUN)),legend=str(legend.relative_to(RUN)),source_sha256=sha(p),source_csv=datafiles,native_tables=native))
        SOURCES[str(p.relative_to(RUN))]=sha(p);SOURCES[str(legend.relative_to(RUN))]=sha(legend)
    validation=writer.prism_file()
    # Replace only the legacy generic project note; table numbers and writer are unchanged.
    legacy=OUT/'Islet_Explorer_plot_data.pzfx';p=OUT/'App_Hypotheses_plot_data.pzfx';doc=ET.parse(str(legacy));ns={'p':writer.NS}
    note=doc.find('.//{'+writer.NS+'}Info/{'+writer.NS+'}Notes')
    note.text='Current app-based H1–H10 data. Native data tables only; figures, complete legends and exact sources are separate. Disease columns are unpaired donor observations. Compartment and model-prediction rows retain donor pairing where explicitly stated. Model estimates and confidence limits are not replicate observations. Blank is missing, not zero. Censored C-peptide bounds are separate from measured values. No quantitative INS endpoint, Prism analysis or graph layout is included.'
    doc.write(str(p),xml_declaration=True,encoding='UTF-8',pretty_print=True)
    ET.XMLSchema(ET.parse(str(OUT/'support/Prism8XMLSchema_plain_cells.xsd'))).assertValid(ET.parse(str(p)))
    legacy.unlink()
    table_rows=[{k:v for k,v in x.items() if k!='frame'}|{'rows':len(x['frame'])} for x in writer.TABLES]
    pd.DataFrame(table_rows).to_csv(OUT/'table_index.csv',index=False)
    pd.DataFrame(figrows).to_csv(OUT/'figure_index.csv',index=False)
    # Statistics remain distinct files, even when also exposed as native tables.
    for p in paths:
        rel=str(p.relative_to(RUN))
        if any(t in rel for t in ['statistics/','models/']) or any(t in p.name for t in ['test','comparison','contrasts','coefficients','metrics','diagnostic']):
            dest=OUT/'statistics'/p.relative_to(RUN);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
    assert all(sha(RUN/k)==v for k,v in SOURCES.items()),'Sources changed during export; rerun after generation finishes.'
    slots=set()
    for name in ['primary_comparisons.csv','primary_tests.csv']:
        for primary in RUN.rglob(name):
            if 'prism' not in primary.relative_to(RUN).parts:slots.update(read(primary).hypothesis.astype(str))
    if '--final' in sys.argv:assert slots==set('H'+str(i) for i in range(1,11)),slots
    validation.update(source_hashes_unchanged=True,source_files=len(SOURCES),figures=len(figrows),quantitative_INS_endpoints=0,
        hypothesis_slots=sorted(slots,key=lambda v:int(v[1:])),all_ten_hypotheses_present=len(slots)==10,
        CI_order='estimate, actual upper limit, actual lower limit',censored_values_separate_from_measured=True,
        exporter_sha256=sha(__file__),prism_sha256=sha(OUT/'App_Hypotheses_plot_data.pzfx'))
    (OUT/'source_hashes.json').write_text(json.dumps(SOURCES,indent=2)+'\n')
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    (OUT/'README.md').write_text('# Current app-based figures: Prism data\n\nOpen `App_Hypotheses_plot_data.pzfx` in GraphPad Prism. It contains native donor-column, paired-compartment, XY pseudotime/held-out-prediction, heatmap-value and estimate/95%CI tables. The `table_index.csv` gives each table’s source and definitions. Graph layouts and Prism analyses are not embedded.\n\nFigures and legends are in `plots` and `legends`. Exact CSV sources preserve donor/object identifiers and categorical metadata in `source_csv`; statistical results are also copied separately to `statistics`. No test is recomputed by this exporter. Disease columns contain unpaired donors; paired tables explicitly retain one donor per row. Object/bin tables are display/source data, not independent clinical replicates.\n\nBlank values remain missing. Confidence limits use estimate, actual upper, actual lower order. Censored C-peptide is a separate upper-bound field; no measured concentration or primary squared error is assigned. The native export contains no quantitative INS endpoint.\n\nThe file passes the existing stricter plain-cell subset of GraphPad’s published schema, exact numeric and row-identity roundtrips, and source-hash checks. It has not been opened in a running Prism application in this environment.\n')
    print(json.dumps(validation),flush=True)


if __name__=='__main__':main()
