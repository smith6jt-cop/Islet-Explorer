"""Signal-support failure cases independent of the SAM model."""
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cohort_expansion'))
from endocrine_signal_support import (hormone_support, mask_support_metrics,
                                      signal_footprint, signal_distances,
                                      keep_seed_near_components)


class SignalSupportTests(unittest.TestCase):
    def test_few_bright_cells_do_not_validate_dark_background(self):
        mask = np.zeros((100,100),bool); mask[10:90,10:90]=True
        raw = np.zeros((100,100,3)); raw[20:25,20:25,2]=10000
        support = hormone_support(raw,[100,100,100],1,smoothing_um=0)
        metric = mask_support_metrics(mask,support,1)
        self.assertLess(metric['area_within_signal_distance_fraction'],.1)

    def test_empty_signal_has_zero_support_even_near_image_edge(self):
        mask = np.ones((20,20),bool)
        metric = mask_support_metrics(mask,np.zeros_like(mask),.5)
        self.assertEqual(metric['area_within_signal_distance_fraction'],0)
        self.assertIsNone(metric['signal_distance_p90_um'])
        self.assertFalse(signal_footprint(np.zeros_like(mask),.5).any())

    def test_large_dark_hole_cannot_be_filled(self):
        yy,xx = np.mgrid[:80,:80]
        r = np.sqrt((xx-40)**2+(yy-40)**2)
        ring = (r>=20)&(r<=24)
        result = signal_footprint(ring,.5)
        self.assertFalse(result[40,40])
        self.assertLessEqual(signal_distances(ring,.5)[result].max(),3)
        self.assertTrue(result[ring].all())

    def test_separate_structures_stay_separate_across_dark_gap(self):
        support = np.zeros((40,70),bool)
        support[10:20,5:15]=True; support[10:20,50:60]=True
        result = signal_footprint(support,.5)
        self.assertFalse(result[:,30:40].any())
        kept = keep_seed_near_components(result,[[10,15]],.5)
        self.assertTrue(kept[15,10]); self.assertFalse(kept[15,55])

    def test_thresholds_are_applied_per_raw_channel(self):
        raw = np.zeros((10,10,3)); raw[2,3,0]=5; raw[5,6,2]=5
        support = hormone_support(raw,[10,10,4],1,smoothing_um=0)
        self.assertFalse(support[2,3]); self.assertTrue(support[5,6])

    def test_mask_metrics_require_nonempty_aligned_geometry(self):
        with self.assertRaises(ValueError):
            mask_support_metrics(np.zeros((5,5),bool),np.ones((5,5),bool),1)
        with self.assertRaises(ValueError):
            mask_support_metrics(np.ones((5,5),bool),np.ones((6,5),bool),1)


if __name__ == '__main__': unittest.main()
