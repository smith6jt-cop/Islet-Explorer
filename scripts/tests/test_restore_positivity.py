import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'research'))
from restore_positivity import sampled_cells,join_compartments,NO_PAIR,cautions,adapter,rv
from common import MARKERS


class FreshRestoreTest(unittest.TestCase):
    def fixture(self):
        frame=pd.DataFrame({'case_id':'1','object_id':[f'1:{i}' for i in range(400)],'Cell Area':10.,'x_source_px':np.arange(400),'y_source_px':1.})
        for marker in MARKERS:frame[marker]=np.arange(400)+1.
        return frame

    def test_sample_is_outcome_blind_and_excludes_quality_fragments(self):
        raw=self.fixture();corrected=raw.copy();corrected.loc[:19,'Cell Area']=1.
        quality=dict(minimum_cell_area_um2=5.,technical_exclusion_regions=[])
        a,_=sampled_cells(raw,corrected,'1',quality,2)
        corrected['donor_status']='T1D';corrected['previous_phenotype']='Beta cell'
        b,_=sampled_cells(raw,corrected,'1',quality,2)
        self.assertEqual(a.object_id.tolist(),b.object_id.tolist())
        self.assertTrue(a['Cell Area'].ge(5).all())
        raw.loc[0,'object_id']='other'
        with self.assertRaisesRegex(ValueError,'universes'):sampled_cells(raw,corrected,'1',quality,2)

    def test_native_means_are_checked_by_uuid_and_channel(self):
        frame=self.fixture();native=pd.DataFrame({'object_id':[str(i) for i in range(len(frame))]})
        for marker in MARKERS:
            native['Cell__'+marker]=frame[marker]
            native['Membrane__'+marker]=frame[marker]/2
        joined=join_compartments(native.iloc[::-1],frame,frame,'1')
        np.testing.assert_array_equal(joined['raw__INS'],joined['Cell__INS'])
        native['Cell__INS']+=500
        with self.assertRaisesRegex(ValueError,'channel identity'):join_compartments(native,frame,frame,'1')

    def test_no_constitutive_reference_is_invented(self):
        self.assertEqual(set(NO_PAIR),{'BActin','NaKATPase','BCatenin','ColIV'})
        self.assertIn('new_donor_no_historical_image_acceptance',cautions('78','CHGA'))
        self.assertIn('historical_donor_channel_problem_requires_current_review',cautions('6549','PDPN'))

    def test_fresh_adapter_preserves_strict_maximum_despite_extreme_control(self):
        adapter.patch_upstream()
        controls=np.r_[np.ones(199),100.]
        statistics=rv.negative_control_statistics(controls)
        self.assertEqual(statistics['maximum'],100.)
        self.assertNotEqual(statistics['maximum'],np.quantile(controls,.995))


if __name__=='__main__':unittest.main()
