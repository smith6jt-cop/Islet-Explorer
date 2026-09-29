import sys
import unittest
from unittest.mock import patch
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'research'))
from common import ratio, require_cohort, MARKERS, RETIRED, write_json, offline_report_root
from trajectory import fit, project, features, landmarks, excluded_for_outcome, apply_identity_review
from phenotypes import read_rules, dependencies, rule_mask, audit
from run_hypotheses import bh_ten, beta_binomial_fit, followups
from prediction import nested_predict
from build_candidate import positive_area, membership

def sample(n=120):
    rng=np.random.default_rng(12);t=np.linspace(-2,2,n)
    d=pd.DataFrame({'islet_id':[f'{i%4}_Islet_{i}' for i in range(n)],'case_id':[str(i%4) for i in range(n)],
                    'donor_status':['ND']*n,'area_um2':np.exp(t+5),'eligible':True})
    for j,c in enumerate(features()):d[c]=3+np.sin(t+j/7)+rng.normal(0,.06,n)
    return d

class ResearchContract(unittest.TestCase):
    def test_label_outcome_and_row_order_invariance(self):
        d=sample();a,m=fit(d,k=10,cap=20,legacy_intensity=True)
        changed=d.copy();changed.donor_status='T1D';changed['ratio']=np.arange(len(d))*100
        b,n=fit(changed.sample(frac=1,random_state=7),k=10,cap=20,legacy_intensity=True)
        np.testing.assert_allclose(a.sort_values('islet_id').pseudotime_raw,b.sort_values('islet_id').pseudotime_raw,equal_nan=True,atol=1e-10)
        self.assertEqual(m['root_id'],n['root_id'])
        self.assertTrue(all(c.split('__')[1] not in ['INS','GCG','SST','CHGA'] for c in m['columns']))

    def test_equal_donor_landmarks(self):
        d=sample();ids=landmarks(d,17)
        self.assertEqual(d.iloc[ids].groupby('case_id').size().tolist(),[17]*4)

    def test_projection_alignment_and_heldout(self):
        d=sample();train=d.loc[d.case_id!='3'];a,m=fit(train,k=10,cap=20,legacy_intensity=True)
        q=d.loc[d.case_id=='3'].sample(frac=1,random_state=1);p=project(q,m)
        self.assertEqual(p.islet_id.tolist(),q.islet_id.tolist())
        self.assertTrue(p.pseudotime_raw.notna().any())
        self.assertTrue(set(q.islet_id).isdisjoint(m['landmark_ids']))

    def test_disconnected_components_not_forced(self):
        d=sample(120);rng=np.random.default_rng(4)
        # Two well-separated angular clusters with more than k landmarks each.
        for j,c in enumerate(features()):d[c]=(np.where(np.arange(len(d))<60,-1,1)*(1+j/50)+rng.normal(0,.001,len(d)))
        a,m=fit(d,k=5,cap=25,legacy_intensity=True)
        self.assertGreater(len(set(m['landmark_components'])),1)
        self.assertTrue(a.pseudotime_raw.isna().any())

    def test_constant_features_dropped_and_missing_rejected(self):
        d=sample();d[features()[0]]=1
        _,m=fit(d,k=10,cap=20,legacy_intensity=True);self.assertIn(features()[0],m['dropped'])
        with self.assertRaises(ValueError):fit(d.drop(columns=features()[1]),k=10,cap=20,legacy_intensity=True)

    def test_outcome_module_and_parent_chain_exclusion(self):
        rules=read_rules(Path(__file__).parents[2]/'data/phenotype_rules.csv')
        omit=excluded_for_outcome('CD8a Tcell',rules)
        self.assertTrue({'CD3e','CD8a','CD4','INS','SMA','CD45'}<=omit)
        self.assertFalse(any(c.split('__')[1] in omit for c in features(omit)))
        with self.assertRaises(ValueError):excluded_for_outcome('unregistered')

    def test_zero_denominator_and_cohort_validation(self):
        raw,plotted=ratio([3,0,9],[0,4,np.nan])
        self.assertTrue(np.isnan(raw[0]) and np.isnan(raw[2]));self.assertEqual(plotted[1],0)
        with self.assertRaises(ValueError):require_cohort(pd.DataFrame({'case_id':['1','2']}),['1','3'])

    def test_ten_slots_with_unavailable_tests(self):
        q=bh_ten([.001,.02]+[np.nan]*8)
        np.testing.assert_allclose(q[:2],[.01,.1]);self.assertTrue(np.isnan(q[2:]).all())

    def test_parent_negative_gate_is_not_ignored(self):
        rules=read_rules(Path(__file__).parents[2]/'data/phenotype_rules.csv')
        chain=dependencies('Beta cell',rules);v=np.zeros((2,31));v[:,MARKERS.index('INS')]=2;v[1,MARKERS.index('CD3e')]=2
        mask=rule_mask(v,MARKERS,chain,{m:1 for m in MARKERS})
        self.assertEqual(mask.tolist(),[True,False])

    def test_new_assignments_do_not_require_archival_labels(self):
        rules=read_rules(Path(__file__).parents[2]/'data/phenotype_rules.csv')
        v=np.zeros((3,31));v[0,MARKERS.index('INS')]=2;v[0,MARKERS.index('GCG')]=2
        v[1,MARKERS.index('INS')]=2;v[1,MARKERS.index('CD3e')]=2;v[2,:]=np.nan
        reports,masks=audit(np.array(['','','']),v,MARKERS,rules,{m:1 for m in MARKERS},MARKERS)
        self.assertEqual(masks['Beta cell'].tolist(),[True,False,False])
        self.assertEqual(masks['Alpha cell'].tolist(),[True,False,False])
        self.assertEqual(int(np.logical_or.reduce([masks[t] for t in ['Alpha cell','Beta cell','Delta cell']]).sum()),1)
        self.assertEqual(next(r for r in reports if r['phenotype']=='Beta cell')['existing_count'],0)
        _,blocked=audit(np.array(['Beta cell']*3),v,MARKERS,rules,{m:1 for m in MARKERS},[m for m in MARKERS if m!='CD20'])
        self.assertFalse(blocked['Beta cell'].any())

    def test_missing_thresholds_are_not_reported_as_zero_rule_matches(self):
        rules=read_rules(Path(__file__).parents[2]/'data/phenotype_rules.csv')
        reports,_=audit(np.array(['']),np.zeros((1,31)),MARKERS,rules,{},[])
        self.assertTrue(all(np.isnan(r['rule_positive_count']) for r in reports))

    def test_validated_marker_name_cannot_substitute_for_missing_threshold(self):
        rules=read_rules(Path(__file__).parents[2]/'data/phenotype_rules.csv')
        gates={m:1 for m in MARKERS};gates['CD20']=np.nan
        reports,masks=audit(np.array(['Beta cell']),np.zeros((1,31)),MARKERS,rules,gates,MARKERS)
        beta=next(r for r in reports if r['phenotype']=='Beta cell')
        self.assertFalse(beta['validated']);self.assertFalse(masks['Beta cell'].any())
        self.assertIn('CD20',beta['missing_validated_gates'])

    def test_sam_alignment_ignores_collinear_vertices_but_not_shift(self):
        import shapely
        from audit_sam_candidates import geometry_key
        box=shapely.box(0,0,12,12)
        imported=shapely.Polygon([(0,0),(6,0),(12,0),(12,12),(0,12)])
        self.assertEqual(geometry_key(box),geometry_key(imported))
        self.assertNotEqual(geometry_key(box),geometry_key(shapely.box(1,0,13,12)))

    def test_candidate_counts_match_raw_cells_without_archival_labels(self):
        import tempfile,types,shapely,anndata
        import build_candidate
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'analysis').mkdir();(root/'raw').mkdir()
            values=np.zeros((5,31),dtype=np.float32)
            values[0,MARKERS.index('INS')]=3;values[0,MARKERS.index('GCG')]=3
            values[1,MARKERS.index('INS')]=3;values[2,:]=np.nan
            values[3,MARKERS.index('INS')]=3
            obs=pd.DataFrame({'imageid':['78']*5,'X_centroid':np.array([2,4,6,20,25])*.3774,
                              'Y_centroid':np.array([2,4,6,20,25])*.3774},index=[f'78:{i}' for i in range(5)])
            cells=anndata.AnnData(values,obs=obs,var=pd.DataFrame(index=MARKERS))
            source=root/'source.h5ad';cells.write_h5ad(source)
            donor=types.SimpleNamespace(case_id='78',cells=source,cell_slot='X',geometry='mocked',
                                        donor_status='T1D',donor_source='HPAP',annotation_review='pending')
            gates=pd.DataFrame({'imageid':['78']*31,'markers':MARKERS,'gate':[1.]*31,'validated':[True]*31})
            rules=read_rules(Path(__file__).parents[2]/'data/phenotype_rules.csv')
            geometry=([('Islet_1',shapely.box(0,0,10,10)),
                       ('Islet_2',shapely.Polygon([(40,40),(50,50),(40,50),(50,40),(40,40)]))],
                      {'Tissue':shapely.box(-1,-1,100,100)}, {})
            with patch.object(build_candidate,'geometry_data',return_value=geometry):
                build_candidate.process(donor,root,rules,gates,gates)
            counts=pd.read_parquet(root/'analysis/78-counts.parquet')
            beta=counts.loc[counts.phenotype.eq('Beta cell')].set_index(['islet_id','compartment'])
            self.assertEqual(beta.loc[('78_Islet_1','core'),'count'],2)
            self.assertEqual(beta.loc[('78_Islet_1','peri'),'count'],1)
            self.assertTrue(np.isnan(beta.loc[('78_Islet_2','core'),'count']))
            regions=pd.read_parquet(root/'analysis/78.parquet')
            core=regions.query("islet_id == '78_Islet_1' and band_um == 20 and compartment == 'core'").iloc[0]
            self.assertEqual(core.endocrine_count,2)  # The INS/GCG double-positive cell counts once.
            self.assertEqual(core.raw_cells,3);self.assertEqual(core.measured_cells,2)
            archived=pd.read_parquet(root/'raw/78-cells.parquet')
            self.assertEqual(len(archived),5);self.assertEqual(archived.marker_measurements_complete.sum(),4)

    def test_metadata_replacement_preserves_linked_source_and_reports_stay_offline(self):
        import tempfile,os,json
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'source.json';copy=root/'copy.json'
            write_json(source,{'version':1});os.link(source,copy);write_json(copy,{'version':2})
            self.assertEqual(json.loads(source.read_text()),{'version':1})
            self.assertFalse(os.path.samefile(source,copy))
            with self.assertRaises(ValueError):offline_report_root(root,root/'research')

    def test_identity_review_cannot_override_classifier_eligibility(self):
        frame=sample(40);review=frame[['islet_id','case_id']].copy()
        review['decision']='confirmed_islet';review['geometry_sha256']='geometry';review['review_note']='Reviewed matching marker channels and tissue context'
        hashes={str(i):'geometry' for i in range(4)}
        review.loc[0,'decision']='non_islet'
        reviewed=apply_identity_review(frame,review,hashes)
        pd.testing.assert_series_equal(reviewed.eligible,frame.eligible)
        self.assertEqual(reviewed.identity_review_decision.iloc[0],'non_islet')
        review.loc[review.case_id=='3','decision']='unreviewed'
        reviewed=apply_identity_review(frame,review,hashes)
        pd.testing.assert_frame_equal(reviewed[frame.columns],frame)
        review['decision']='confirmed_islet';review.loc[0,'geometry_sha256']='modified'
        with self.assertRaisesRegex(ValueError,'geometry differs'):apply_identity_review(frame,review,hashes)

    def test_unknown_results_do_not_choose_followup_branches(self):
        d=pd.DataFrame({'id':[f'H{i}' for i in range(1,11)],'q_value':np.nan})
        self.assertTrue(followups(d).status.str.contains('awaiting').all())

    def test_duplicate_bootstrap_donors_never_cross_prediction_folds(self):
        class Spy:
            def fit(self,X,y):
                self.donors=set(X[:,0]);self.coef_=np.zeros((1,X.shape[1]));return self
            def predict(self,X):
                if self.donors & set(X[:,0]):raise AssertionError('Donor leakage into validation fold')
                return np.zeros(len(X),int)
            def __getitem__(self,index):return self
        groups=np.repeat(np.arange(8),2);X=np.c_[groups,groups**2];y=groups%2
        with patch('prediction.estimator',lambda c:Spy()):
            prediction,_,_=nested_predict(X,y,groups)
        self.assertEqual(len(prediction),16)

    def test_structural_area_unions_overlapping_masks(self):
        import shapely
        region=shapely.box(0,0,10,10);a=shapely.box(0,0,5,10)
        parts=np.array([a,a],dtype=object)
        self.assertEqual(positive_area(region,(parts,shapely.STRtree(parts))),50)

    def test_native_polygon_edge_is_not_shifted_by_export_rounding(self):
        import shapely
        # Actual HPAP export failure: native float vertex lies on the cell;
        # the JTS GeoJSON 0.01-pixel representation moves the right edge out.
        edge=float(np.float32(7055.59))
        x=np.array([edge,edge-.1]);y=np.array([2270.365478515625]*2)
        points=shapely.STRtree(shapely.points(x,y))
        exact=shapely.box(7050,2266,edge,2275)
        rounded=shapely.box(7050,2266,7055.59,2275)
        self.assertEqual(membership(points,exact,x,y).tolist(),[1])
        self.assertEqual(membership(points,rounded,x,y).tolist(),[0,1])
        # GeometryROI explicitly includes JTS boundary points, unlike this
        # PolygonROI right-edge winding convention.
        self.assertEqual(membership(points,exact,x,y,'qupath_geometry').tolist(),[0,1])
        with self.assertRaises(ValueError):membership(points,exact,x,y,'unknown')

    def test_beta_binomial_rejects_impossible_count(self):
        with self.assertRaises(ValueError):beta_binomial_fit(np.ones((2,1)),[3,5],[2,10])

if __name__=='__main__':unittest.main()
