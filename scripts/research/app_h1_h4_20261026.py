#!/usr/bin/env python3
"""Link unchanged approved H1–H4 endpoints/tests to the saved app state.

No hypothesis model, embedding, pseudotime, or cell classification is fitted.
The saved model outputs are reused only after checking their input bindings.
"""
from pathlib import Path
import hashlib,json,os,re,shutil,sys
os.environ.setdefault('MPLCONFIGDIR','/tmp/app_hypotheses_mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import app_hypotheses_shared_20261026 as S

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/app_hypotheses_20261026/h1_h4'
OLD=ROOT/'output/h1_h4_disease_comparisons_20261026'
PAPER=ROOT/'output/h1_h4_literature_methods_20261025'
AUDIT=ROOT/'output/h1_h10_conditional_audit_20261026/association_audit'
GROUPS=S.GROUPS;COLORS=S.COLORS
SOURCES={};FIGURES=[]
MODEL_METHODS={
 'H1':'Exposure: the core fraction of the unique core-plus-peri immune pool, among nearby-immune-exposed objects. Response: core Beta-cell count/core area in cells/mm². Exposure and continuous covariates are standardized within donor; response is standardized across included objects. Adjustments: log core area, log peri/core area ratio, log1p nearby immune count, log1p 20–100µm field immune density, and object class. The recorded correlated and independent random-slope fits are singular; the retained primary comparison uses one separately adjusted coefficient per donor, Kruskal–Wallis across disease groups and three two-sided Wilcoxon rank-sum comparisons with BH adjustment. Within-state means reuse equal donor weights, Student-t 95% intervals and two-sided one-sample t tests, with BH across the three states.',
 'H2':'Exposure: threshold-positive ColIV coverage in the symmetric ±2µm endocrine-boundary band after excluding the union of capillary pixels. Response: core immune-cell count/core area in cells/mm². Exposure and continuous covariates are standardized within donor; response is standardized across included objects. Adjustments: log1p peri immune density, log core area, log perimeter, log1p distance to tissue edge, union-raster capillary core area fraction, and object class. The recorded primary random-slope fits are singular; the retained primary comparison uses one separately adjusted coefficient per donor, Kruskal–Wallis and three two-sided Wilcoxon rank-sum comparisons with BH adjustment. Within-state mean coefficients and tests use the already approved equal-donor Student-t approach, with BH across the three states. The ±5µm and log1p-density sensitivity outputs remain separately identified.',
 'H3':'Endpoint: donor near-crossing immune density minus reference-band immune density. Both compartments lie within 10µm inside the endocrine boundary; near-crossing geometry is the core intersection with a 10µm buffer of unioned capillary–boundary crossings, clipped to that boundary band. The reference is the remainder. Counts and assessed areas are summed separately across eligible objects for each donor, then divided to obtain cells/mm². Zero immune counts remain valid. The primary disease comparison is Kruskal–Wallis plus three two-sided Wilcoxon rank-sum tests with BH adjustment. Existing paired two-sided Wilcoxon signed-rank tests compare near and reference densities within state; their BH family includes the pooled and three state tests at this width. This is an observed-compartment comparison, not a spatially randomized null.',
 'H4':'Exposure: macrophage LGALS3 signal divided by macrophage-plus-epithelial LGALS3 signal, using nonzero measured source pools. Signal is the sum of served scaled cell means, not protein mass. Macrophages retain their accepted identity; epithelial source is accepted Exocrine with ECAD or panCK ≥0.10. Response: core Beta-cell density. Exposure and continuous adjustments are standardized within donor; response is standardized across included objects. Adjustments: combined source signal per neighborhood area, macrophage density, epithelial-source density, nearby immune density, core area, and object class. The primary LMM has a donor random intercept and source slope and is nonsingular. Disease-by-source interaction and within-state trends use lmerTest/emmeans Satterthwaite inference. Pairwise slope contrasts use BH across three disease pairs; existing within-state trend tests are unadjusted. Recovered structural qualification does not replace accepted cell identities.'}
METHOD_SOURCE_FILES=['h1_disease_random_slopes_20261026.R','h2_disease_comparisons_20261026.R','h3_disease_comparisons_20261026.py','h4_disease_comparisons_20261026.R','h2_coliv_boundary_20261025.py','h3_boundary_compartments_20261025.py','h4_source_features_20261025.py','h1_h2_h5_conditional_audit_20261026.py']

def method_metadata():
    (OUT/'model_methods.txt').write_text('\n\n'.join(h+': '+v for h,v in MODEL_METHODS.items())+'\n\nMethod sources: Barlow2025 PDF p22; Bodenmiller2025 PDF p28; Richardson2023 PDF pp3/8; Heide PDF pp30/32.\n')
    spec=json.loads((OUT/'analysis_specification.json').read_text());spec['exact_endpoint_and_model_definitions']=MODEL_METHODS
    (OUT/'analysis_specification.json').write_text(json.dumps(spec,indent=2)+'\n')
    for path in (OUT/'figures').glob('*_legend.txt'):
        stem=path.stem.removesuffix('_legend')
        selected=[stem[:2].upper()] if stem[:2] in ['h1','h2','h3','h4'] else ['H3'] if stem=='crossing_compartment_comparison' else ['H1','H2','H4'] if stem=='conditional_associations' else list(MODEL_METHODS)
        text=path.read_text().split('\n\nModel and measurement definitions\n')[0]
        path.write_text(text+'\n\nModel and measurement definitions\n'+'\n\n'.join(h+': '+MODEL_METHODS[h] for h in selected)+'\n')
    for name in METHOD_SOURCE_FILES:
        path=ROOT/'scripts/research'/name;SOURCES[str(path.relative_to(ROOT))]=S.digest(path)

def read(path):
    path=Path(path);SOURCES[str(path.relative_to(ROOT))]=S.digest(path)
    return S.read(path)

def copied(path,name):
    d=read(path);d.to_csv(OUT/'statistics'/name,index=False);return d

def figure(fig,stem,title,data,statistics,methods,note=''):
    S.savefig(fig,OUT/'figures',stem)
    FIGURES.append(dict(stem=stem,title=title,data=data,statistics=statistics,methods=methods,note=note))

def bound_join(endpoint,shared,hypothesis):
    assert endpoint.islet_id.is_unique
    e=endpoint.set_index('islet_id');s=shared.set_index('islet_id').loc[e.index]
    assert e.case_id.astype(str).equals(s.case_id.astype(str))
    assert e.donor_status.equals(s.donor_status)
    for c in ['n_beta_core','n_immune_core','n_immune_peri','core_area_um2','peri_area_um2']:
        if c in e:assert np.allclose(e[c],s[c],equal_nan=True), (hypothesis,c)
    additions=e[[c for c in e if c not in s]]
    d=s.join(additions).reset_index()
    d['hypothesis']=hypothesis
    assert d.case_id.nunique()==20 and not d.case_id.eq('78').any()
    return d

def inputs():
    shared=read(ROOT/'output/app_hypotheses_20261026/shared/object_measurements.parquet')
    h1=read(PAPER/'h1/analysis_objects.csv')
    # The exact input builder is pure until its CLI main is called.
    import h1_paper_input_20261025 as H1
    generated=H1.build_frame().to_csv(index=False).encode()
    assert hashlib.sha256(generated).hexdigest()==H1.EXPECTED_CSV_SHA256
    assert generated==(PAPER/'h1/analysis_objects.csv').read_bytes()
    h2=read(PAPER/'h2/h2_analysis_objects.csv')
    h3=read(PAPER/'h3/h3_object_compartments.csv')
    h3=h3.loc[h3.width_um.eq(10)&h3.status.eq('included')].copy()
    assert len(h3)==1865
    h3['near_density_mm2']=h3.n_near*1e6/h3.near_area_um2
    h3['reference_density_mm2']=h3.n_reference*1e6/h3.reference_area_um2
    h3['difference_density_mm2']=h3.near_density_mm2-h3.reference_density_mm2
    h4=read(OLD/'h4/primary_analysis_objects.csv')
    data={h:bound_join(d,shared,h) for h,d in [('H1',h1),('H2',h2),('H3',h3),('H4',h4)]}
    assert [len(data[h]) for h in data]==[2937,2937,1865,3681]
    # Check all saved H4 source file bindings, including original cell phenotypes
    # and LGALS3/epithelial marker measurements; none is rewritten here.
    hp=PAPER/'h4/source_hashes.json';SOURCES[str(hp.relative_to(ROOT))]=S.digest(hp)
    hashes=json.loads(hp.read_text());changed=[p for p,h in hashes.items() if S.digest(ROOT/p)!=h]
    assert not changed,changed
    h3_provenance=PAPER/'h3/source_hashes.json'
    SOURCES[str(h3_provenance.relative_to(ROOT))]=S.digest(h3_provenance)
    h3_hashes=json.loads(h3_provenance.read_text())
    assert all(S.digest(Path(p))==v for p,v in h3_hashes.items())
    h2_provenance=PAPER/'h2/source_provenance.json'
    SOURCES[str(h2_provenance.relative_to(ROOT))]=S.digest(h2_provenance)
    h2_sources=json.loads(h2_provenance.read_text())
    assert S.digest(ROOT/h2_sources['feature_input'])==h2_sources['feature_input_sha256']
    for r in h2_sources['images']:
        assert S.digest(ROOT/r['current_geometry'])==r['geometry_sha256']
        assert S.digest(ROOT/r['capillary_objects'])==r['capillary_sha256']
        st=Path(r['source_image']).stat()
        assert st.st_size==r['bytes'] and st.st_mtime_ns==r['mtime_ns']
    # Record overlap of accepted phenotype and separately recovered structural
    # qualification; this audit never overwrites accepted cell identities.
    cross=[]
    for path in sorted((ROOT/'output/pancfloyd_20donor_20261026/structural_cells/per_cell').glob('*.parquet')):
        SOURCES[str(path.relative_to(ROOT))]=S.digest(path)
        d=pd.read_parquet(path,columns=['app_phenotype','recovered_structure_class','edge_analysis_eligible','case_id'])
        d=d.loc[d.edge_analysis_eligible&d.recovered_structure_class.notna()&d.recovered_structure_class.ne('')]
        cross.append(d.groupby(['case_id','app_phenotype','recovered_structure_class']).size().rename('n_cells').reset_index())
    pd.concat(cross,ignore_index=True).to_csv(OUT/'data/accepted_phenotype_structural_qualification_crosstab.csv',index=False)
    audits=[
      dict(hypothesis='H1',input_unchanged=True,model_refit=False,objects=len(h1),donors=20,
        dependency='Accepted Beta/immune counts, core/peri geometry, regional field; no structural-cell or structural-area model term.',
        verification='Exact regenerated input CSV bytes and shared corrected counts/areas match.'),
      dict(hypothesis='H2',input_unchanged=True,model_refit=False,objects=len(h2),donors=20,
        dependency='ColIV raster coverage excludes a geometric union of vessels; the capillary covariate already uses the same union raster. No structural-cell label term.',
        verification='h2_coliv_boundary_20261025.py::boundary_features calls hf.union_local before rasterization; no summed-annotation area input. All recorded geometry hashes and source-image sizes/mtime match.'),
      dict(hypothesis='H3',input_unchanged=True,model_refit=False,objects=len(h3),donors=20,
        dependency='Crossing geometry is union_all of crossings; near/reference areas are complementary compartments. No structural-cell label term.',
        verification=f'h3_boundary_compartments_20261025.py::crossing_geometry/measure_compartments; no summed-annotation area input. All {len(h3_hashes)} recorded source hashes match.'),
      dict(hypothesis='H4',input_unchanged=True,model_refit=False,objects=len(h4),donors=20,
        dependency='Accepted Macrophage source and accepted Exocrine + ECAD/panCK epithelial operation; LGALS3 source/immune densities and core area. Structural qualification stays separate.',
        verification=f'All {len(hashes)} recorded source hashes match; shared corrected counts/areas match.')]
    pd.DataFrame(audits).to_csv(OUT/'input_correction_audit.csv',index=False)
    for h,d in data.items():
        # These are count/state/source features. No quantitative INS column is exported.
        d=d[[c for c in d if 'INS' not in c and 'ins_' not in c.lower()]]
        data[h]=d;d.to_parquet(OUT/'data'/f'{h.lower()}_analysis_objects.parquet',index=False)
        d.to_csv(OUT/'data'/f'{h.lower()}_analysis_objects.csv',index=False)
    counts=pd.concat([d.groupby(['case_id','donor_status']).agg(n_objects=('islet_id','size'),n_embedded=('in_composition_embedding','sum'),n_saved_display_eligible=('eligible_plotted','sum')).reset_index().assign(hypothesis=h) for h,d in data.items()])
    counts.to_csv(OUT/'data/cohort_counts.csv',index=False)
    return data,audits

def statistics(data):
    tables={};donors={}
    mapping={
      'H1':('h1/donor_coefficients.csv','h1/donor_disease_global_test.csv','h1/donor_pairwise_disease_tests.csv'),
      'H2':('h2/primary_donor_coefficients.csv','h2/primary_disease_overall_test.csv','h2/primary_disease_pairwise_tests.csv'),
      'H3':('h3/h3_primary_donor_values.csv','h3/h3_global_test.csv','h3/h3_pairwise_disease_tests.csv')}
    for h,paths in mapping.items():
        donors[h]=copied(OLD/paths[0],f'{h.lower()}_donor_values.csv')
        tables[h]=dict(global_test=copied(OLD/paths[1],f'{h.lower()}_global_test.csv'),pairwise=copied(OLD/paths[2],f'{h.lower()}_pairwise_tests.csv'))
    slopes=copied(OLD/'h4/disease_slopes.csv','h4_all_disease_slopes.csv');slopes=slopes.loc[slopes.model.eq('primary')].copy()
    hpairs=copied(OLD/'h4/disease_slope_comparisons.csv','h4_all_pairwise_tests.csv');hpairs=hpairs.loc[hpairs.model.eq('primary')].copy()
    hanova=copied(OLD/'h4/model_anova.csv','h4_all_model_tests.csv')
    slopes.to_csv(OUT/'statistics/h4_primary_disease_slopes.csv',index=False)
    hpairs.to_csv(OUT/'statistics/h4_primary_pairwise_tests.csv',index=False)
    tables['H4']=dict(slopes=slopes,pairwise=hpairs,global_test=hanova.loc[hanova.model.eq('primary')&hanova.term.eq('donor_status:source_z')])
    donors['H4']=copied(OLD/'h4/donor_fitted_slopes.csv','h4_conditional_donor_slopes.csv')
    old_primary=read(OLD/'primary_disease_comparison_tests.csv')
    primary=[]
    for r in old_primary.itertuples():
        primary.append(dict(hypothesis=r.hypothesis,comparison='ND vs Aab+ vs T1D: primary disease contrast',test=r.test,p_value=r.p,n_donors=20,n_objects=len(data[r.hypothesis]),status='completed_reused_unchanged',source=str((OLD/'primary_disease_comparison_tests.csv').relative_to(ROOT))))
    pd.DataFrame(primary).to_csv(OUT/'primary_comparisons.csv',index=False)
    cond=read(AUDIT/'state_conditional_means.csv');cond=cond[cond.hypothesis.isin(['H1','H2'])].copy()
    cond.to_csv(OUT/'statistics/h1_h2_within_disease_conditional_means.csv',index=False)
    h3=read(PAPER/'h3/h3_summary.csv');h3=h3[h3.width_um.eq(10)&h3.group.isin(GROUPS)].copy()
    h3.to_csv(OUT/'statistics/h3_within_disease_paired_tests.csv',index=False)
    within=[]
    for r in cond.itertuples():
        within.append(dict(hypothesis=r.hypothesis,donor_status=r.donor_status,estimate=r.mean_beta,lower_95=r.lower_95,upper_95=r.upper_95,p_value=r.p_value,p_adjusted=r.p_bh_three_states,n_donors=r.donors,n_objects=int(data[r.hypothesis].donor_status.eq(r.donor_status).sum()),test=r.test,estimate_definition='Equal-donor mean adjusted standardized regression coefficient',multiplicity='BH across three states within hypothesis',status='completed_reused_unchanged'))
    for _,r in slopes.iterrows():
        within.append(dict(hypothesis='H4',donor_status=r.donor_status,estimate=r['source_z.trend'],lower_95=r['lower.CL'],upper_95=r['upper.CL'],p_value=r['p.value'],p_adjusted=np.nan,n_donors=r.donors,n_objects=r.objects,test='lmerTest/emmeans Satterthwaite two-sided slope test',estimate_definition='Adjusted standardized regression coefficient from donor intercept/slope model',multiplicity='Existing within-state p values unadjusted; direct disease pairwise BH remains separate',status='completed_reused_unchanged'))
    for r in h3.itertuples():
        within.append(dict(hypothesis='H3',donor_status=r.group,estimate=r.median_difference_mm2,lower_95=np.nan,upper_95=np.nan,p_value=r.p_value,p_adjusted=r.p_fdr_four_groups,n_donors=r.donors,n_objects=r.objects,test='Paired two-sided Wilcoxon signed-rank test of donor near/reference densities',estimate_definition='Median donor near-minus-reference immune density (cells/mm²)',multiplicity='Existing BH family: pooled and three disease groups at width 10µm',status='completed_reused_unchanged'))
    within=pd.DataFrame(within);within.to_csv(OUT/'within_disease_estimates_tests.csv',index=False)
    # Preserve technical diagnostics and relevant sensitivity output as data, without
    # running old report generators or copying their interpretation text.
    for h,files in {'h1':['model_diagnostics.csv','group_slopes.csv','pairwise_slope_differences.csv'],
                    'h2':['model_diagnostics.txt','model_tests.csv','standardized_slopes.csv','disease_slope_comparisons.csv'],
                    'h4':['model_diagnostics.json']}.items():
        for name in files:
            path=OLD/h/name;SOURCES[str(path.relative_to(ROOT))]=S.digest(path)
            shutil.copyfile(path,OUT/'statistics'/f'{h}_{name}')
    return tables,donors,within

def stars(p):return '***' if p<.001 else '**' if p<.01 else '*' if p<.05 else ''

def bracket(ax,i,j,y,h,p):
    if not np.isfinite(p) or p>=.05:return
    ax.plot([i,i,j,j],[y,y+h,y+h,y],color='black',lw=.8);ax.text((i+j)/2,y+h,stars(p),ha='center',va='bottom')

def conditional_figures(data,tables,donors,within):
    titles={'H1':'Immune localization and β-cell density','H2':'ColIV coverage and core immune density','H3':'Immune density at capillary crossings','H4':'LGALS3 source and β-cell density'}
    fig,axes=plt.subplots(1,4,figsize=(17,4.7),layout='constrained')
    plotrows=[]
    for ax,h in zip(axes,['H1','H2','H3','H4']):
        if h!='H4':
            y='difference_density_mm2' if h=='H3' else 'standardized_beta';d=donors[h]
            arr=[d.loc[d.donor_status.eq(g),y].to_numpy() for g in GROUPS]
            bp=ax.boxplot(arr,positions=range(3),widths=.55,patch_artist=True,showfliers=False)
            for box,g in zip(bp['boxes'],GROUPS):box.set(facecolor=COLORS[g],alpha=.22)
            for i,g in enumerate(GROUPS):
                z=d.loc[d.donor_status.eq(g)].sort_values('case_id').copy();x=i+np.linspace(-.12,.12,len(z))
                ax.scatter(x,z[y],c=COLORS[g],s=30,zorder=3)
                for xx,(_,r) in zip(x,z.iterrows()):plotrows.append(dict(hypothesis=h,case_id=r.case_id,donor_status=g,x=xx,estimate=r[y],lower=np.nan,upper=np.nan))
            lo=min(d[y].min(),0);hi=max(d[y].max(),0);span=max(hi-lo,.05)
            pairs=tables[h]['pairwise']
            for k,(_,r) in enumerate(pairs.iterrows()):
                a,b=(r['group1'],r['group2']) if h=='H2' else (r['group_1'],r['group_2'])
                p=r['p_BH'] if h=='H2' else r['p_bh'];bracket(ax,GROUPS.index(a),GROUPS.index(b),hi+span*(.1+.13*k),span*.035,p)
            ax.set_ylim(lo-.1*span,hi+.5*span)
        else:
            d=tables[h]['slopes'];lo=min(d['lower.CL'].min(),0);hi=max(d['upper.CL'].max(),0);span=max(hi-lo,.05)
            for i,g in enumerate(GROUPS):
                r=d.loc[d.donor_status.eq(g)].iloc[0];y=r['source_z.trend']
                ax.errorbar(i,y,yerr=[[y-r['lower.CL']],[r['upper.CL']-y]],fmt='o',color=COLORS[g],capsize=4,ms=6)
                plotrows.append(dict(hypothesis=h,case_id='',donor_status=g,x=i,estimate=y,lower=r['lower.CL'],upper=r['upper.CL']))
            for k,(_,r) in enumerate(tables[h]['pairwise'].iterrows()):
                i,j=[(0,1),(0,2),(1,2)][k];bracket(ax,i,j,hi+span*(.1+.13*k),span*.035,r.p_bh)
            ax.set_ylim(lo-.1*span,hi+.5*span)
        ax.axhline(0,c='#888888',ls='--',lw=.8);ax.set_xticks(range(3),GROUPS);ax.set_title(titles[h],fontsize=10);ax.set_xlabel('Disease state')
        ax.set_ylabel('Near − reference immune density\n(cells/mm²)' if h=='H3' else 'Adjusted standardized regression coefficient (β)')
    pd.DataFrame(plotrows).to_csv(OUT/'data/direct_disease_plot_values.csv',index=False)
    figure(fig,'direct_disease_comparisons','Direct disease comparisons',['data/direct_disease_plot_values.csv'],['primary_comparisons.csv','statistics/h1_pairwise_tests.csv','statistics/h2_pairwise_tests.csv','statistics/h3_pairwise_tests.csv','statistics/h4_primary_pairwise_tests.csv'],
        'H1/H2/H3: individual donor values with median/IQR and 1.5×IQR whiskers. H4: model trends and pointwise 95% CIs. Brackets use the preserved three-pair BH adjustments; no bracket for p≥.05.')
    fig,axes=plt.subplots(1,3,figsize=(13,4.6),layout='constrained')
    pp=within[within.hypothesis.isin(['H1','H2','H4'])].copy()
    for ax,h in zip(axes,['H1','H2','H4']):
        z=pp[pp.hypothesis.eq(h)];lo=min(z.lower_95.min(),0);hi=max(z.upper_95.max(),0);span=max(hi-lo,.05)
        for i,g in enumerate(GROUPS):
            r=z.loc[z.donor_status.eq(g)].iloc[0]
            ax.errorbar(i,r.estimate,yerr=[[r.estimate-r.lower_95],[r.upper_95-r.estimate]],fmt='o',color=COLORS[g],capsize=4,ms=7)
            if np.isfinite(r.p_adjusted) and r.p_adjusted<.05:ax.text(i,r.upper_95+.03*span,stars(r.p_adjusted),ha='center')
        ax.set(xticks=range(3),xticklabels=GROUPS,xlabel='Disease state',title=titles[h],ylim=(lo-.1*span,hi+.2*span))
        ax.axhline(0,c='#888888',ls='--',lw=.8);ax.set_ylabel('Adjusted standardized coefficient (β)\n95% confidence interval')
    pp.to_csv(OUT/'data/conditional_association_plot_values.csv',index=False)
    figure(fig,'conditional_associations','Within-disease conditional associations',['data/conditional_association_plot_values.csv'],['within_disease_estimates_tests.csv'],
        'H1/H2: existing equal-donor mean coefficients with Student-t 95% intervals (n−1 df), and existing two-sided one-sample t tests; BH across three states separately for each hypothesis. H4: existing nonsingular donor intercept/slope LMM trends, Satterthwaite tests and pointwise 95% CIs. Stars only show preserved adjusted p values below .05.')
    fig,axes=plt.subplots(1,3,figsize=(12,4.6),layout='constrained')
    p=donors['H3'].copy();pairs=[]
    for ax,g in zip(axes,GROUPS):
        z=p[p.donor_status.eq(g)].sort_values('case_id')
        for _,r in z.iterrows():
            ax.plot([0,1],[r.reference_density_mm2,r.near_density_mm2],'-o',color=COLORS[g],alpha=.65,lw=1,ms=4)
            for x,col in enumerate(['reference_density_mm2','near_density_mm2']):pairs.append(dict(case_id=r.case_id,donor_status=g,compartment='Reference' if x==0 else 'Near crossing',plot_x=x,density_cells_mm2=r[col]))
        ax.set(xticks=[0,1],xticklabels=['Reference band','Near crossings'],title=g,ylabel='Immune density (cells/mm²)')
    pd.DataFrame(pairs).to_csv(OUT/'data/crossing_compartment_plot_values.csv',index=False)
    figure(fig,'crossing_compartment_comparison','Core boundary immune compartments',['data/crossing_compartment_plot_values.csv'],['statistics/h3_within_disease_paired_tests.csv'],
        'Each connected pair is one donor: summed immune counts divided by summed observed compartment area. Primary width 10µm. Paired two-sided Wilcoxon tests and the existing BH family of pooled plus three states are reused. Zero immune counts are retained.')

def neighborhoods(data):
    rows=[]
    fields=[('Core','immune_density_core_mm2'),('Peri 0–20 µm','immune_density_peri_mm2'),('Field 20–100 µm','immune_density_field20_100_mm2')]
    for h,d in data.items():
        for compartment,col in fields:
            q=d.groupby(['case_id','donor_status'])[col].agg(['mean','count']).reset_index()
            q['hypothesis']=h;q['compartment']=compartment;q['density_cells_mm2']=q['mean'];rows.append(q)
    q=pd.concat(rows,ignore_index=True);q['log1p_density']=np.log1p(q.density_cells_mm2)
    q.to_csv(OUT/'data/immune_neighborhood_donor_values.csv',index=False)
    fig,axes=plt.subplots(1,4,figsize=(16,4.5),layout='constrained')
    titles=['Localization–Beta cohort','ColIV–immune cohort','Crossing-compartment cohort','LGALS3-source cohort']
    for ax,h,title in zip(axes,data,titles):
        z=q[q.hypothesis.eq(h)]
        for (case,g),a in z.groupby(['case_id','donor_status']):
            a=a.set_index('compartment').loc[[v for v,c in fields]]
            ax.plot(range(3),a.log1p_density,'o-',ms=3,lw=.7,alpha=.5,c=COLORS[g])
        ax.set(xticks=range(3),xticklabels=['Core','Peri\n0–20 µm','Field\n20–100 µm'],ylabel='log1p(mean object-associated\nimmune density, cells/mm²)',title=title)
    for g in GROUPS:axes[0].plot([],[],color=COLORS[g],marker='o',label=g)
    axes[0].legend(frameon=False)
    figure(fig,'immune_neighborhood_profiles','Immune density across core, peri and field',['data/immune_neighborhood_donor_values.csv'],[],
        'Each line is one donor. For each endpoint cohort, values are the arithmetic mean of object-associated density in the named compartment, then displayed as log1p. Shared peri/field cells retain existing per-object membership; these are not deduplicated whole-tissue densities. No test is added.')

def states(data):
    features={
      'H1':[('immune_localization_fraction','Core fraction of nearby immune cells'),('beta_density_mm2','β-cell density (cells/mm²)')],
      'H2':[('coliv_coverage_2um','ColIV boundary coverage (fraction)'),('immune_density_core_mm2','Core immune density (cells/mm²)')],
      'H3':[('near_density_mm2','Crossing-near immune density (cells/mm²)'),('reference_density_mm2','Reference immune density (cells/mm²)')],
      'H4':[('myeloid_fraction','Myeloid fraction of source LGALS3'),('beta_density_mm2','β-cell density (cells/mm²)')]}
    for h,d in data.items():
        result=S.state_views(d,features[h],OUT,h.lower(),f'{len(d):,} eligible endocrine objects in the original {h} endpoint cohort; 20 donors. Endpoint definitions and conditional statistics are separate.')
        for stem in result['figures']:
            FIGURES.append(dict(stem=stem,title='App composition-state and pseudotime measurements',data=[f'data/{h.lower()}_state_objects.csv',f'data/{h.lower()}_pseudotime_profiles.csv',f'data/{h.lower()}_donor_pseudotime.csv',f'data/{h.lower()}_donor_cluster_values.csv',f'data/{h.lower()}_cluster_heatmap_values.csv'],statistics=[],methods='Saved composition UMAP/Leiden0.5 and unchanged full-cohort rank DPT. Twenty fixed-bin object means; equal-donor cluster means. No new test.',n_objects=result['objects'],n_donors=result['donors'],features=features[h]))
        cd=pd.read_csv(OUT/'data'/f'{h.lower()}_donor_cluster_values.csv')
        cd.groupby(['feature','donor_status','leiden_0.5'],observed=True)['mean'].agg(['mean','count']).reset_index().rename(columns={'mean':'equal_donor_mean','count':'n_donors'}).to_csv(OUT/'data'/f'{h.lower()}_cluster_heatmap_values.csv',index=False)
    features=[('structure_blood_vessel_core_density_cells_mm2','Structure-qualified vessel cells (cells/mm²)'),
              ('structure_blood_vessel_core_mean_passthrough_fraction_among_measured','Mean vessel passthrough area fraction'),
              ('structure_capillary_core_area_percent','Capillary union area / core area (%)')]
    result=S.state_views(data['H3'],features,OUT,'h3_structure_context',
        'Descriptive companions in the H3 endpoint cohort. Vessel cells follow centroid, union-area containment and marker/highest-marker qualification. Passthrough is the cell-area union overlap fraction among measured inside cells without their own vessel call; empty pools are missing. Annotation union area remains a separate endpoint; no localization test is replaced.')
    for stem in result['figures']:
        FIGURES.append(dict(stem=stem,title='Vessel cells, passthrough and union area',data=[f'data/h3_structure_context_{x}.csv' for x in ['state_objects','pseudotime_profiles','donor_pseudotime','donor_cluster_values','cluster_heatmap_values']],statistics=[],methods='Separate recovered structural-cell density, passthrough fraction and annotation-union area; descriptive saved-state views only.',n_objects=result['objects'],n_donors=result['donors'],features=features))
    cd=pd.read_csv(OUT/'data/h3_structure_context_donor_cluster_values.csv')
    cd.groupby(['feature','donor_status','leiden_0.5'],observed=True)['mean'].agg(['mean','count']).reset_index().rename(columns={'mean':'equal_donor_mean','count':'n_donors'}).to_csv(OUT/'data/h3_structure_context_cluster_heatmap_values.csv',index=False)

def validate(data):
    for stem,h in [('h1','H1'),('h2','H2'),('h3','H3'),('h4','H4'),('h3_structure_context','H3')]:
        q=pd.read_csv(OUT/'data'/f'{stem}_state_objects.csv').set_index('islet_id')
        d=data[h].set_index('islet_id');expected=d.loc[d.in_composition_embedding]
        assert set(q.index)==set(expected.index)
        for c in ['umap_1','umap_2','pseudotime','pt_rank','leiden_0.5']:
            assert np.allclose(q[c].astype(float),d.loc[q.index,c].astype(float),equal_nan=True),(stem,c)
        profiles=pd.read_csv(OUT/'data'/f'{stem}_pseudotime_profiles.csv')
        for field,p in profiles.groupby('feature'):
            calc=q.reset_index().groupby(['donor_status','pt_bin'])[field].agg(['mean','count'])
            observed=p.set_index(['donor_status','pt_bin'])
            assert np.allclose(observed['mean'],calc.loc[observed.index,'mean'],equal_nan=True)
            assert np.array_equal(observed['count'],calc.loc[observed.index,'count'])
        cd=pd.read_csv(OUT/'data'/f'{stem}_donor_cluster_values.csv')
        heat=pd.read_csv(OUT/'data'/f'{stem}_cluster_heatmap_values.csv').set_index(['feature','donor_status','leiden_0.5'])
        calc=cd.groupby(['feature','donor_status','leiden_0.5'])['mean'].agg(['mean','count'])
        assert np.allclose(heat.equal_donor_mean,calc.loc[heat.index,'mean'],equal_nan=True)
        assert np.array_equal(heat.n_donors,calc.loc[heat.index,'count'])
    hd=pd.read_csv(OUT/'statistics/h3_donor_values.csv',dtype={'case_id':str}).set_index('case_id')
    calc=data['H3'].groupby('case_id')[['n_near','n_reference','near_area_um2','reference_area_um2']].sum()
    for c in calc:assert np.allclose(hd.loc[calc.index,c],calc[c])
    for f in FIGURES:
        for suffix in ['png','pdf','svg']:
            assert (OUT/'figures'/f"{f['stem']}.{suffix}").stat().st_size>1000
    report=dict(saved_coordinates_and_DPT_unchanged=True,state_object_ids_exact=True,
        fixed_bin_means_and_counts_verified=True,equal_donor_cluster_means_verified=True,
        H3_donor_count_area_aggregation_verified=True,figure_sets=len(FIGURES),
        quantitative_INS_column_exported=False,model_refits=0,new_scientific_tests=0)
    for d in data.values():assert not any('INS' in c or 'ins_' in c.lower() for c in d)
    (OUT/'validation.json').write_text(json.dumps(report,indent=2)+'\n')

def main():
    for folder in ['data','statistics','figures']:(OUT/folder).mkdir(parents=True,exist_ok=True)
    specification=dict(population='Existing accepted 20 donors; donor78 excluded',representation='Shared saved app composition UMAP/Leiden0.5 and exact fitted app DPT; no refit/rerank',primary='Reuse existing direct disease comparisons',conditional='Reuse existing H1/H2 equal-donor mean tests, H3 paired donor tests and H4 mixed-model trends',new_models=False,new_scientific_tests=False,INS_intensity_endpoint=False,structural_cell_identity='Separate qualification; accepted app phenotype retained',branch_judgment='Not assigned')
    (OUT/'analysis_specification.json').write_text(json.dumps(specification,indent=2)+'\n')
    data,audits=inputs();tables,donors,within=statistics(data)
    states(data);conditional_figures(data,tables,donors,within);neighborhoods(data)
    execution=read(ROOT/'output/h1_h10_conditional_audit_20261026/plan_execution_audit.csv')
    execution=execution[execution.hypothesis.isin(data)].copy();execution['current_app_representation_outputs']='Completed: saved UMAP, fixed DPT bins, saved-cluster/donor summaries, immune neighborhood displays'
    execution.to_csv(OUT/'original_plan_execution_audit.csv',index=False)
    for f in FIGURES:
        if not f['stem'].endswith(('_states','_clusters')):
            text=f['methods']+'\n\nAll quantitative values and tests are in the linked CSVs. Donors are the biological replicates. No biological interpretation or original-condition branch judgment is assigned.\n'
            for path in f['statistics']:
                if '*' in path:
                    for p in sorted(OUT.glob(path)):text+='\n'+p.name+'\n'+p.read_text()
                else:text+='\n'+path+'\n'+(OUT/path).read_text()
            (OUT/'figures'/f"{f['stem']}_legend.txt").write_text(text)
    summary=pd.DataFrame([dict(hypothesis=h,n_objects=len(d),n_donors=d.case_id.nunique(),n_mapped_objects=int(d.in_composition_embedding.sum()),n_mapped_donors=d.loc[d.in_composition_embedding,'case_id'].nunique()) for h,d in data.items()])
    summary.to_csv(OUT/'cohort_summary.csv',index=False)
    (OUT/'README.md').write_text('# App-linked H1–H4 figures and statistics\n\nExisting approved endpoints and statistical results are reused after input checks. No model, embedding or DPT was refitted; no new scientific test or quantitative INS endpoint was added. Recovered structural qualification is retained separately from accepted phenotypes. The structural correction audit and accepted-phenotype cross-tab are included.\n\nThe primary disease comparisons are in `primary_comparisons.csv`; within-disease coefficients and tests are in `within_disease_estimates_tests.csv`. `input_correction_audit.csv` records why the later structural-area/cell corrections do not alter these model inputs. H1/H2 preserve the donor-coefficient fallback after singular random-slope fits; H4 preserves the nonsingular primary donor-intercept/slope model. H3 preserves observed compartment comparisons rather than a randomized spatial null.\n\nSaved state displays use each endpoint cohort and its mapped subset. Counts are in `cohort_summary.csv` and `data/cohort_counts.csv`; smaller nonembedded objects remain in the endpoint models. Figures contain measured values; separate legends specify methods/statistics. Original conditions are recorded as conditions in `original_plan_execution_audit.csv`, with no branch judgment assigned.\n')
    validate(data)
    method_metadata()
    for f in FIGURES:
        if 'n_objects' not in f:
            selected=['H3'] if f['stem']=='crossing_compartment_comparison' else ['H1','H2','H4'] if f['stem']=='conditional_associations' else list(data)
            f['n_objects_by_hypothesis']={h:len(data[h]) for h in selected}
            f['n_donors_by_hypothesis']={h:int(data[h].case_id.nunique()) for h in selected}
    for p in [Path(__file__),Path(S.__file__),ROOT/'docs/app_based_hypothesis_plan_20261026.md',ROOT/'docs/T1D_research_hypotheses_20261022.md']:
        SOURCES[str(p.relative_to(ROOT))]=S.digest(p)
    assert all(S.digest(ROOT/p)==v for p,v in SOURCES.items()),'Source changed during run'
    manifest=dict(figures=FIGURES,cohorts=summary.to_dict('records'),sources=SOURCES,models_refitted=False,new_scientific_tests=False,branch_judgment='Not assigned',files={str(p.relative_to(OUT)):S.digest(p) for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='manifest.json'})
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(summary.to_string(index=False));print(f'{len(FIGURES)} figure sets; all source bindings match.')

if __name__=='__main__':
    if '--metadata-only' in sys.argv:
        p=OUT/'manifest.json';m=json.loads(p.read_text());SOURCES.update(m['sources'])
        method_metadata();SOURCES[str(Path(__file__).relative_to(ROOT))]=S.digest(__file__)
        # A spelling-only documentation edit may occur while collaborators work.
        # Rebind the plan only after its quoted H1–H4 conditions still match.
        plan_path=ROOT/'docs/T1D_research_hypotheses_20261022.md';plan=plan_path.read_text()
        execution=pd.read_csv(OUT/'original_plan_execution_audit.csv')
        for r in execution.itertuples():
            section=plan.split('### '+r.hypothesis+' ',1)[1].split('\n### ',1)[0]
            conditions=re.findall(r'\*\*(If[^*]+)\*\*\s*(.*?)(?=\n\n|\Z)',section,flags=re.S)
            assert conditions[0][1].strip()==r.original_if_supported
            assert conditions[1][1].strip()==r.original_if_excluded
        SOURCES[str(plan_path.relative_to(ROOT))]=S.digest(plan_path)
        assert all(S.digest(ROOT/p)==v for p,v in SOURCES.items()),'Non-metadata source changed'
        m['sources']=SOURCES
        m['files']={str(p.relative_to(OUT)):S.digest(p) for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='manifest.json'}
        p.write_text(json.dumps(m,indent=2)+'\n');print('Methods metadata refreshed; no numerical output or figure changed.')
    else:main()
