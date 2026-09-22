import sys
import unittest
import gzip
import json
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from shapely.geometry import box, mapping
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'research'))
from common import MARKERS, sha256
from composition import canonical_cells,aggregate_donor
from materialize_classifier_eligibility import materialize
from trajectory import fit,project,reference_percentile,graph_layout
from cell_typing import TypingConfig,fit_ensemble_member,load_model,project_clusters,consensus_assignments


def fixture(n=12):
    identity=pd.DataFrame({'case_id':['1']*n,'object_id':[str(i) for i in range(n)],
                          'base_type':['Beta cell']*5+['Hormone-unresolved endocrine cell']*2+['B cell']*(n-7),
                          'classification_status':['validated']*n,'projection_supported':True,'is_endocrine':[True]*7+[False]*(n-7)})
    corrected=identity[['case_id','object_id']].copy()
    for marker in MARKERS:corrected[marker]=np.arange(n)+1
    geometry=pd.DataFrame([dict(islet_id='1_Islet_1',case_id='1',area_um2=100,geometry_valid=True,
                              core_tissue_coverage=1.,core_valid_area_um2=100.,peri_valid_area_um2=200.,
                              classifier_defined_core=True,image_edge_contact=False,image_edge_clearance_um=5.)])
    members=identity[['case_id','object_id']].assign(islet_id='1_Islet_1',compartment='core',band_um=20)
    return identity,corrected,geometry,members


class CompositionTest(unittest.TestCase):
    def test_classifier_eligibility_materializes_edge_and_core_cell_rules(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'analysis').mkdir();(root/'corrected').mkdir();native=root/'native';native.mkdir();adaptive=root/'adaptive';adaptive.mkdir()
            source=root/'source.json';source.write_text(json.dumps({'image_shape_cyx':[31,100,100],'pixel_size_um':.5}))
            candidates=adaptive/'review_candidates.geojson';candidates.write_text('{}')
            model_sha='a'*64
            (adaptive/'prediction_manifest.json').write_text(json.dumps({'boundary_config':{'edge_margin_um':5.},'model_sha256':model_sha}))
            def feature(name,geometry):
                return {'type':'Feature','geometry':mapping(geometry),'properties':{'name':name,
                    'classification':{'name':'Endocrine pixel candidate'},
                    'metadata':{'Segmentation method':'INS/GCG/SST foreground; fixed test','Model SHA256':model_sha}}}
            features=[feature('inside',box(10,10,20,20)),feature('edge',box(0,30,10,40))]
            cores=native/'cores.geojson.gz'
            with gzip.open(cores,'wt') as stream:json.dump({'type':'FeatureCollection','features':features},stream)
            members=pd.DataFrame({'case_id':['1','1'],'object_id':['a','b'],'islet_key':['inside','edge']})
            members.to_csv(native/'native_core_memberships.csv',index=False)
            pd.DataFrame({'object_id':['a','b'],'marker_measurements_complete':[True,True]}).to_parquet(root/'corrected/1-cells.parquet',index=False)
            (native/'manifest.json').write_text(json.dumps({'cores_sha256':sha256(cores),'memberships_sha256':sha256(native/'native_core_memberships.csv')}))
            pd.DataFrame({'case_id':['1'],'boundary_source':[source],'native_core_root':[native],'boundary_candidates':[candidates]}).to_csv(root/'cohort.csv',index=False)
            pd.DataFrame({'case_id':['1','1'],'islet_id':['1_inside','1_edge'],'islet_key':['inside','edge']}).to_parquet(root/'analysis/islets.parquet',index=False)
            members.assign(islet_id=['1_inside','1_edge'],compartment='core',band_um=20).drop(columns='islet_key').to_parquet(root/'analysis/memberships.parquet',index=False)
            report=materialize(root);result=pd.read_parquet(root/'analysis/islets.parquet').set_index('islet_key')
            self.assertEqual(report['eligible_islets'],1)
            self.assertTrue(result.loc['inside','classifier_eligible']);self.assertFalse(result.loc['edge','classifier_eligible'])
            self.assertEqual(result.loc['edge','classifier_exclusion_reason'],'image_edge_contact')

    def test_endocrine_union_and_assessed_empty_peri(self):
        t,c,g,m=fixture();cells=canonical_cells(t,c)
        f,counts,summary=aggregate_donor(g,m,cells,['Beta cell','Hormone-unresolved endocrine cell','B cell'])
        self.assertEqual(f.endocrine_count.iloc[0],12)
        self.assertAlmostEqual(f['core__count__B cell'].iloc[0],np.log1p(5/12))
        self.assertEqual(f['peri__count__B cell'].iloc[0],0)
        self.assertTrue(f.eligible.iloc[0])
        self.assertEqual(set(summary.marker),set(MARKERS))
        self.assertFalse(summary.positive_available.any())
        self.assertTrue(summary.loc[summary.population.eq('Endocrine union')&summary.compartment.eq('core'),'mean_log1p'].notna().all())

    def test_unknowns_are_not_zeros_or_endocrine_cells(self):
        t,c,g,m=fixture();t.loc[:,'classification_status']='unreviewed';t['is_endocrine']=pd.NA
        f,counts,_=aggregate_donor(g,m,canonical_cells(t,c),['Beta cell','B cell'])
        self.assertEqual(f.endocrine_count.iloc[0],12);self.assertTrue(f.eligible.iloc[0])
        self.assertFalse(f.trajectory_feature_supported.iloc[0])
        self.assertTrue(counts.loc[counts.compartment.eq('core'),'count'].isna().all())
        t.loc[0,'is_endocrine']=True
        with self.assertRaisesRegex(ValueError,'Unvalidated'):canonical_cells(t,c)

    def test_quality_excluded_fragments_do_not_enter_markers_or_minimum_cells(self):
        t,c,g,m=fixture();t.loc[:3,'classification_status']='segmentation_quality_excluded';t['is_endocrine']=t.is_endocrine.astype('boolean');t.loc[:3,'is_endocrine']=pd.NA
        t['segmentation_quality_excluded']=False;t.loc[:3,'segmentation_quality_excluded']=True
        cells=canonical_cells(t,c)
        f,counts,s=aggregate_donor(g,m,cells,['Beta cell','Hormone-unresolved endocrine cell','B cell'])
        self.assertEqual(f.measured_cells.iloc[0],8);self.assertEqual(f.excluded_cells.iloc[0],4)
        self.assertTrue(f.eligible.iloc[0]);self.assertEqual(f.endocrine_count.iloc[0],12)
        self.assertTrue(s.query("population == 'All cells' and compartment == 'core'").n_cells.eq(8).all())

    def test_only_image_edge_and_empty_classifier_core_control_eligibility(self):
        t,c,g,m=fixture();g['primary_assessable']=False
        f,counts,_=aggregate_donor(g,m,canonical_cells(t,c),['Beta cell','B cell'])
        self.assertTrue(f.eligible.iloc[0]);self.assertTrue(counts.loc[counts.compartment.eq('core'),'count'].notna().all())
        g['image_edge_contact']=True;g['overlap']=True
        f,_,_=aggregate_donor(g,m,canonical_cells(t,c),['Beta cell'])
        self.assertFalse(f.eligible.iloc[0]);self.assertEqual(f.exclusion_reason.iloc[0],'image_edge_contact')
        self.assertTrue(f.overlap.iloc[0])

    def test_empty_classifier_core_is_ineligible_and_denominator_stays_zero(self):
        t,c,g,m=fixture();m=m.iloc[:0].copy()
        f,counts,_=aggregate_donor(g,m,canonical_cells(t,c),['B cell'])
        self.assertEqual(f.endocrine_count.iloc[0],0);self.assertTrue(counts.raw_ratio.isna().all())
        self.assertEqual(f.exclusion_reason.iloc[0],'no_core_cells_in_classifier_area')

    def test_simultaneous_states_preserve_base_counts(self):
        t,c,g,m=fixture();t['phenotype_states']=[['Beta cell: CD44 high','Beta cell: LGALS3 high'] if i<3 else [] for i in range(len(t))]
        names=['Beta cell: CD44 high','Beta cell: LGALS3 high']
        _,counts,_=aggregate_donor(g,m,canonical_cells(t,c),['Beta cell','B cell','Hormone-unresolved endocrine cell'],names)
        core=counts.loc[counts.compartment.eq('core')].set_index('phenotype')
        self.assertEqual(core.loc['Beta cell','count'],5)
        self.assertEqual(core.loc[names[0],'count'],3);self.assertEqual(core.loc[names[1],'count'],3)

    def test_graph_layout_retains_disconnected_components_without_pseudotime(self):
        from scipy import sparse
        ids=np.array([f'islet{i}' for i in range(6)])
        f=pd.DataFrame({'islet_id':ids,'case_id':'1','eligible':True,'x':[1,1,2,-1,-1,-2],'y':[1,2,1,-1,-2,-1]})
        coordinates=f[['islet_id','case_id']].assign(component=[0,0,0,1,1,1],pseudotime_raw=[0.,.5,1.,np.nan,np.nan,np.nan])
        before=coordinates.copy(deep=True)
        model=dict(landmark_components=np.array([0,0,0,1,1,1]),graph=sparse.block_diag([np.ones((3,3))-np.eye(3)]*2).tocsr(),
            landmark_ids=ids,landmark_vectors=f[['x','y']].to_numpy(float),columns=['x','y'],center=np.zeros(2),scale=np.ones(2),k=2)
        layout=graph_layout(f,coordinates,model)
        pd.testing.assert_frame_equal(before,coordinates)
        self.assertTrue(layout[['graph_x','graph_y']].notna().all().all());self.assertEqual(set(layout.component),{0,1})
        self.assertGreater(layout.loc[layout.component.eq(1),'graph_x'].min(),layout.loc[layout.component.eq(0),'graph_x'].max())

    def test_positivity_requires_complete_reviewed_calls(self):
        t,c,g,m=fixture();t['positive__INS']=pd.Series([True]*6+[False]*6,dtype='boolean')
        _,_,s=aggregate_donor(g,m,canonical_cells(t,c),['Beta cell','B cell','Hormone-unresolved endocrine cell'])
        ins=s.query("marker == 'INS' and population == 'All cells' and compartment == 'core'").iloc[0]
        self.assertEqual(ins.positive_count,6);self.assertAlmostEqual(ins.positive_ratio,6/12)
        t.loc[0,'positive__INS']=pd.NA
        _,_,s=aggregate_donor(g,m,canonical_cells(t,c),['Beta cell'])
        self.assertFalse(s.query("marker == 'INS' and population == 'All cells' and compartment == 'core'").positive_available.iloc[0])

    def test_membership_duplicates_and_cross_donor_rejected(self):
        t,c,g,m=fixture();cells=canonical_cells(t,c)
        with self.assertRaisesRegex(ValueError,'unique'):aggregate_donor(g,pd.concat([m,m.iloc[:1]]),cells,['Beta cell'])
        m.loc[0,'case_id']='2'
        with self.assertRaisesRegex(ValueError,'Mixed'):aggregate_donor(g,m,cells,['Beta cell'])

    def test_composition_fit_is_label_invariant(self):
        rng=np.random.default_rng(44);n=160;x=np.linspace(0,3,n)
        f=pd.DataFrame({'islet_id':[str(i) for i in range(n)],'case_id':[str(i%4) for i in range(n)],'area_um2':np.exp(x+4),'eligible':True})
        for i,name in enumerate(['Alpha','Beta','B','CD8']):
            for comp in ['core','peri']:f[f'{comp}__count__{name}']=np.log1p(np.exp(np.sin(x+i)+rng.normal(0,.03,n)))
        a,model=fit(f,k=10,cap=25)
        altered=f.copy();altered['donor_status']='T1D';altered['plotted_ratio']=np.arange(n)*19
        b,other=fit(altered,k=10,cap=25)
        np.testing.assert_allclose(a.pseudotime_raw,b.pseudotime_raw,equal_nan=True)
        self.assertEqual(model['root_id'],other['root_id'])
        projected=project(f,model)
        np.testing.assert_allclose(a.pseudotime_raw,projected.pseudotime_raw,equal_nan=True)
        np.testing.assert_allclose(a.pseudotime_percentile,projected.pseudotime_percentile,equal_nan=True)
        with self.assertRaisesRegex(ValueError,'rebuilding typing'):fit(f,excluded=['INS'])

    def test_omitted_training_marker_cannot_change_reconstructed_counts(self):
        import tempfile
        rng=np.random.default_rng(301);n=120
        original=pd.DataFrame(rng.lognormal(0,.15,(n,len(MARKERS))),columns=MARKERS)
        original.loc[:59,'INS']+=30;original.loc[60:,'GCG']+=30
        original['case_id']='1';original['object_id']=[str(i) for i in range(n)]
        altered=original.copy();altered['CD45']=np.geomspace(1e-8,1e10,n)
        config=TypingConfig(omitted_markers=('CD45',),neighbor_counts=(5,10),neighbor_backend='exact')
        projected=[{},{}];mapping={}
        with tempfile.TemporaryDirectory() as directory:
            for method in ['log1p','arcsinh']:
                for k in [5,10]:
                    fit_id=f'{method}_k{k}';models=[]
                    for version,source in enumerate([original,altered]):
                        path=Path(directory)/str(version)/fit_id
                        fit_ensemble_member(source,config,method,k,1,path,f'input-{version}')
                        model=load_model(path);models.append(model)
                        self.assertNotIn('CD45',model['preprocessing']['markers'])
                        markers=model['preprocessing']['input_markers']
                        projected[version][fit_id]=project_clusters(model['index'],source[markers].to_numpy(),model['cluster_labels'],model['preprocessing'],model['support'],model['k'],model['backend'],model['projection_agreement'])
                    self.assertEqual(models[0]['preprocessing'],models[1]['preprocessing'])
                    np.testing.assert_array_equal(models[0]['cluster_labels'],models[1]['cluster_labels'])
                    labels=models[0]['cluster_labels']
                    for cluster in np.unique(labels):
                        ids=np.flatnonzero(labels==cluster)
                        name='Beta cell' if np.mean(ids<60)>=.5 else 'Alpha cell'
                        mapping[(fit_id,int(cluster))]={'decision':'accepted','base_type':name,'lineage':'Endocrine','phenotype_state':''}
            identities=[consensus_assignments(original[['case_id','object_id']],p,mapping) for p in projected]
            pd.testing.assert_frame_equal(identities[0],identities[1])
            geometry=pd.DataFrame([dict(case_id='1',islet_id='1_Islet_1',area_um2=1000.,geometry_valid=True,core_tissue_coverage=1.,core_valid_area_um2=1000.,peri_valid_area_um2=1000.,classifier_defined_core=True,image_edge_contact=False,image_edge_clearance_um=5.)])
            members=original[['case_id','object_id']].assign(islet_id='1_Islet_1',compartment='core',band_um=20)
            a,_,_=aggregate_donor(geometry,members,canonical_cells(identities[0],original),['Alpha cell','Beta cell'],include_marker_summaries=False)
            b,_,_=aggregate_donor(geometry,members,canonical_cells(identities[1],altered),['Alpha cell','Beta cell'],include_marker_summaries=False)
            pd.testing.assert_frame_equal(a,b)
            self.assertEqual(a.endocrine_count.iloc[0],n)

    def test_frozen_percentiles_handle_ties(self):
        ref=np.array([0,1,1,2.]);np.testing.assert_allclose(reference_percentile(ref,ref),[.125,.5,.5,.875])
        np.testing.assert_allclose(reference_percentile([-1,3],ref),[0,1])

if __name__=='__main__':unittest.main()
