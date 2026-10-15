import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'research'))
from analysis_v2 import model_matrix,negative_binomial_fit,region_wide,positive_counts,endocrine_profiles,endpoint
from common import MARKERS
from structural_coupling import shape_preserving_proximity
from exploratory_v2 import lymphatic_localization,nested_signature,censored_c_peptide
from ins_independent_signature import predict_fold,FEATURE_MARKERS,STRUCTURAL
from spatial_followup_endpoints import donor_endpoints
from outcome_trajectory import outcome_module,held_out_smooth


class RevisedAnalysisTest(unittest.TestCase):
    def test_all_three_disease_categories_in_within_donor_design(self):
        rng=np.random.default_rng(4);n=180
        d=pd.DataFrame({'case_id':np.repeat(np.arange(18),10).astype(str),'donor_status':np.repeat(['ND','Aab+','T1D'],60),
                        'log_area':rng.normal(5,1,n),'age':np.repeat(np.arange(18)+20,10),
                        'sex':np.repeat(['male','female']*9,10),'lymphatic':rng.uniform(0,1,n)})
        x=model_matrix(d,exposure='lymphatic')
        self.assertIn('aab_status',x);self.assertEqual(int(x.aab_status.sum()),60)
        self.assertEqual(int(x.disease.sum()),60)
        np.testing.assert_allclose(x.groupby(d.case_id).lymphatic_within.mean(),0,atol=1e-15)

    def test_negative_binomial_area_offset_recovers_density(self):
        rng=np.random.default_rng(11);n=1000;x=rng.normal(size=n);area=rng.uniform(30,200,n)
        mu=area*np.exp(-2+.5*x);y=rng.negative_binomial(5,5/(5+mu))
        beta=negative_binomial_fit(np.c_[np.ones(n),x],y,np.log(area),np.ones(n))
        np.testing.assert_allclose(beta,[-2,.5],atol=.12)
        with self.assertRaisesRegex(ValueError,'Invalid'):negative_binomial_fit(np.ones((1,1)),[-1],[0],[1])

    def test_missing_structural_mask_never_becomes_zero(self):
        r=pd.DataFrame({'islet_id':['a','a'],'band_um':[20,20],'geometry_valid':True,'valid_area_um2':100,
                        'compartment':['core','peri'],'nerve_area_fraction':[0.,0.],
                        'nerve_validated':[True,False],'overlap':False})
        w=region_wide(r).iloc[0]
        self.assertEqual(w.nerve_area_fraction_core,0)
        self.assertTrue(np.isnan(w.nerve_area_fraction_peri));self.assertTrue(np.isnan(w.capillary_area_fraction_core))

    def test_core_endocrine_endpoint_survives_peri_overlap(self):
        base=pd.DataFrame({'islet_id':[f'i{i}' for i in range(12)],'case_id':[str(i%4) for i in range(12)],
            'donor_status':['ND','T1D']*6,'age':30.,'sex':['m','f']*6,'log_area':np.arange(12)+3.,
            'success':4.,'trials':10.,'primary_assessable':True,'overlap':True})
        regions=base[['islet_id']].merge(pd.DataFrame({'compartment':['core','peri']}),how='cross')
        regions['band_um']=20;regions['geometry_valid']=True;regions['valid_area_um2']=100.
        regions['overlap']=regions.compartment.eq('peri');regions['primary_assessable']=True
        frame,_=endpoint('H1',base,regions,pd.DataFrame(),pd.DataFrame(),Path('/tmp'))
        self.assertEqual(len(frame),12);self.assertTrue(frame.valid_area_um2_peri.isna().all())

    def test_marker_calls_unavailable_remain_missing(self):
        d=pd.DataFrame({'islet_id':['a','a'],'compartment':['core','peri'],'population':'All cells',
                        'marker':'CD45','positive_count':[0.,0.],'positive_available':[True,False]})
        p=positive_counts(d,'CD45').iloc[0]
        self.assertEqual(p.CD45_core,0);self.assertTrue(np.isnan(p.CD45_peri))

    def test_spatial_null_preserves_shapes_and_reports_infeasibility(self):
        import shapely
        region=shapely.box(0,0,100,100);nerve=shapely.box(10,10,12,12);capillary=shapely.box(10,10,11,11)
        a=shape_preserving_proximity(region,nerve,capillary,1.,replicates=30,seed=7)
        b=shape_preserving_proximity(region,nerve,capillary,1.,replicates=30,seed=7)
        self.assertEqual(a,b);self.assertEqual(a['null_successes'],30);self.assertEqual(a['null_status'],'supported')
        self.assertGreater(a['excess_proximity'],0)
        failure=shape_preserving_proximity(region,region,capillary,1.,replicates=30,seed=7,max_proposals=30)
        self.assertEqual(failure['null_successes'],0);self.assertTrue(np.isnan(failure['excess_proximity']))
        self.assertEqual(failure['null_status'],'insufficient_feasible_shape_preserving_placements')

    def test_lymphatic_localization_uses_actual_positive_area(self):
        d=pd.DataFrame({'case_id':np.repeat(np.arange(6).astype(str),2),'compartment':['peri','distant']*6,
                        'positive_area_um2':[20.,5.]*6,'valid_area_um2':[100.]*12,'validated':True})
        result=lymphatic_localization(d,bootstrap=20)
        self.assertAlmostEqual(result['effect'],.15)
        self.assertEqual(result['n_donors'],6)

    def test_outcome_modules_remove_associated_phenotypes(self):
        catalog={'base_types':[{'name':'CD8 T cell','lineage':'Immune','marker_dependencies':['CD3e','CD8a']},
                               {'name':'Beta cell','lineage':'Endocrine','marker_dependencies':['INS','CHGA']}]}
        omitted,types=outcome_module('CD8 T cell',catalog)
        self.assertIn('CD45',omitted);self.assertIn('CD8 T cell',types);self.assertNotIn('Beta cell',types)
        omitted,_=outcome_module('INS',catalog)
        self.assertTrue({'INS','GCG','SST','CHGA'}<=omitted)

    def test_heldout_smooth_never_uses_heldout_outcomes_to_fit(self):
        rng=np.random.default_rng(77);n=240;t=np.linspace(-2,2,n)
        f=pd.DataFrame({'case_id':[str(i%8) for i in range(n)],'islet_id':[f'i{i}' for i in range(n)],'area_um2':np.exp(t+6),'eligible':True})
        for k in range(8):f[f'core__count__type{k}']=np.log1p(np.exp(np.sin(t+k/3)+rng.normal(0,.05,n)))
        clinical=pd.DataFrame({'case_id':[str(i) for i in range(8)],'age':[25,31,45,26,38,56,33,40],
                               'sex':['m','f']*4,'donor_status':['ND','Aab+','T1D','ND','Aab+','T1D','ND','T1D']})
        y=pd.DataFrame({'islet_id':f.islet_id,'response':t**2+rng.normal(0,.1,n)})
        _,a=held_out_smooth(f,'7',y,clinical)
        altered=y.copy();altered.loc[f.case_id.eq('7'),'response']+=100
        _,b=held_out_smooth(f,'7',altered,clinical)
        np.testing.assert_allclose(a.smooth_prediction,b.smooth_prediction)
        np.testing.assert_allclose(a.pseudotime_percentile,b.pseudotime_percentile)

    def test_nested_INS_signature_does_not_fit_heldout_burden(self):
        rng=np.random.default_rng(95);n=72;donors=[str(i) for i in range(6)]
        profile=pd.DataFrame({'case_id':np.repeat(donors,12),'islet_id':[f's{i}' for i in range(n)]})
        x=rng.uniform(size=n)
        for j,m in enumerate(FEATURE_MARKERS):profile[m]=x+rng.normal(0,.1+j*.01,n)
        endpoint=profile[['case_id','islet_id']].copy();endpoint['insulin_burden']=.1+.6*x+rng.normal(0,.05,n)
        for c in STRUCTURAL:endpoint[c]=rng.uniform(size=n)
        keys={('0',)}|{tuple(sorted(('0',d))) for d in donors if d!='0'}
        frames={k:profile.copy() for k in keys}
        a,_=predict_fold('0',frames,endpoint,donors)
        changed=endpoint.copy();changed.loc[changed.case_id.eq('0'),'insulin_burden']+=100
        b,_=predict_fold('0',frames,changed,donors)
        np.testing.assert_allclose(a.held_out_signature,b.held_out_signature)
        self.assertNotIn('INS',FEATURE_MARKERS)

    def test_lymphatic_builder_uses_exact_cells_and_positive_pixels(self):
        import tempfile,json,shapely
        with tempfile.TemporaryDirectory() as path:
            root=Path(path);(root/'analysis').mkdir();(root/'typing/cells').mkdir(parents=True);(root/'corrected').mkdir()
            features=[]
            for name,g in [('Tissue',shapely.box(0,0,100,100)),('Islet',shapely.box(40,40,60,60)),('Lymphatic',shapely.box(20,44,25,48))]:
                features.append({'type':'Feature','geometry':shapely.geometry.mapping(g),'properties':{'name':'Islet_1','classification':{'name':name}}})
            geo=root/'geometry.json';geo.write_text(json.dumps({'type':'FeatureCollection','features':features}))
            names=['CD8 T cell','CD4 T cell','B cell','Macrophage'];catalog={'base_types':[{'name':x} for x in names]}
            (root/'typing/catalog.json').write_text(json.dumps(catalog))
            cells=pd.DataFrame({'case_id':'1','object_id':[str(i) for i in range(8)],'x_source_px':[35.,65.]*4,'y_source_px':46.})
            for marker in MARKERS:cells[marker]=1.
            cells.to_parquet(root/'corrected/1-cells.parquet',index=False)
            typed=cells[['case_id','object_id']].assign(base_type=np.repeat(names,2),classification_status='validated',projection_supported=True,is_endocrine=False)
            typed.to_parquet(root/'typing/cells/1.parquet',index=False)
            membership=cells[['case_id','object_id']].assign(islet_id='1_Islet_1',compartment='peri',band_um=20)
            region=pd.DataFrame([dict(case_id='1',islet_id='1_Islet_1',band_um=20,compartment='peri',geometry_valid=True,overlap=False,lymphatic_validated=True,pixel_size_um=1.)])
            membership.to_parquet(root/'analysis/memberships.parquet',index=False);region.to_parquet(root/'analysis/regions.parquet',index=False)
            counts,areas=donor_endpoints(root,dict(case_id='1',geometry=str(geo)),catalog,region,membership)
            self.assertEqual(counts.near_count.tolist(),[1.]*4);self.assertEqual(counts.total_count.tolist(),[2.]*4)
            self.assertTrue(areas.positive_area_um2.ge(0).all());self.assertTrue(areas.valid_area_um2.gt(0).all())

    def test_profiles_select_typed_endocrine_population(self):
        records=[]
        for population,value,n in [('All cells',999,100),('Endocrine union',2,5)]:
            records.extend(dict(islet_id='a',case_id='1',compartment='core',population=population,marker=m,mean_log1p=value,n_cells=n) for m in MARKERS)
        p=endocrine_profiles(pd.DataFrame(records)).iloc[0]
        self.assertEqual(p.INS,2);self.assertEqual(p.endocrine_measured_cells,5);self.assertTrue(p.profile_eligible)

if __name__=='__main__':unittest.main()
