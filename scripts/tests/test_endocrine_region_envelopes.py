"""Whole-region geometry regressions; biological identity is not inferred."""
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cohort_expansion'))
from endocrine_region_envelopes import region_envelopes


class EnvelopeTests(unittest.TestCase):
    def test_enclosed_interior_does_not_split_whole_region(self):
        yy,xx=np.mgrid[:80,:80];r=np.hypot(xx-40,yy-40)
        ring=(r>=15)&(r<=20)
        result=region_envelopes(ring,[[58,40]],.5)
        self.assertEqual(result.max(),1)
        self.assertGreater(result[40,40],0)
        self.assertEqual(result[5,5],0)

    def test_scattered_cells_do_not_get_a_shared_dark_envelope(self):
        support=np.zeros((80,100),bool)
        support[10:20,10:20]=True;support[55:65,70:80]=True
        result=region_envelopes(support,[[15,15],[75,60]],.5)
        self.assertEqual(result.max(),2)
        self.assertEqual(result[35,45],0)

    def test_open_dark_crescent_is_not_filled_as_a_hole(self):
        yy,xx=np.mgrid[:80,:80];r=np.hypot(xx-40,yy-40)
        crescent=(r>=15)&(r<=20)&(xx<48)
        result=region_envelopes(crescent,[[22,40]],.5)
        self.assertEqual(result[40,40],0)

    def test_edge_signal_is_not_eroded_away(self):
        support=np.zeros((40,40),bool);support[:10,:10]=True
        result=region_envelopes(support,[[3,3]],.5)
        self.assertTrue((result[support]>0).all())

    def test_unsupported_components_are_retained_only_in_raw_signal(self):
        support=np.zeros((40,70),bool);support[5:15,5:15]=True;support[25:35,55:65]=True
        result=region_envelopes(support,[[10,10]],.5)
        self.assertGreater(result[10,10],0);self.assertEqual(result[30,60],0)
        self.assertTrue(support[30,60])


if __name__=='__main__':unittest.main()
