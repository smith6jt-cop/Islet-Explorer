"""Check prompt placement without models; run with the analysis environment."""
import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'cohort_expansion'))
from prepare_endocrine_candidates import hormone_points


class HormonePointTests(unittest.TestCase):
    def test_dark_nuclear_centre_moves_to_hormone_ring(self):
        yy, xx = np.mgrid[:31, :31]
        radius = np.sqrt((xx-15)**2+(yy-15)**2)
        rgb = np.zeros((31, 31, 3), dtype=np.uint8)
        rgb[(radius >= 5) & (radius <= 7), 1] = 200
        centres = np.array([[15., 15.]])
        point = hormone_points(rgb, centres, 8)[0]
        self.assertGreaterEqual(np.linalg.norm(point-centres[0]), 4)
        self.assertLessEqual(np.linalg.norm(point-centres[0]), 8)
        self.assertTrue(np.array_equal(centres, [[15., 15.]]))

    def test_distant_brighter_structure_cannot_attract_prompt(self):
        rgb = np.zeros((40, 40, 3), dtype=np.uint8)
        rgb[18:21, 18:21, 0] = 120
        rgb[2:8, 2:8, 0] = 255
        centres = np.array([[22., 22.], [0., 0.]])
        points = hormone_points(rgb, centres, 5)
        self.assertTrue((np.linalg.norm(points-centres, axis=1) <= 5).all())
        self.assertTrue((points >= 0).all())
        self.assertTrue(np.array_equal(points, hormone_points(rgb, centres, 5)))


if __name__ == '__main__': unittest.main()
