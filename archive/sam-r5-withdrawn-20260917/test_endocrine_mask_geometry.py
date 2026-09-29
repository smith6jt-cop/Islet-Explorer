"""Spatial safety checks for candidate masks; run in the existing SAM environment."""
import sys
from pathlib import Path
import unittest

import numpy as np
from shapely.geometry import Point

sys.path.insert(0, str(Path(__file__).resolve().parent))
from predict_endocrine_candidates import geometry, supported_components


class CandidateMaskTests(unittest.TestCase):
    def test_disconnected_supported_islands_both_survive(self):
        mask = np.zeros((30, 30), dtype=bool)
        mask[3:8, 3:8] = True
        mask[20:26, 20:26] = True
        mask[2:4, 20:22] = True
        kept, recall, removed = supported_components(mask, np.array([[5, 5], [22, 22]]))
        result = geometry(kept, 1, np.array([0, 0]))
        self.assertEqual(len(result.geoms), 2)
        self.assertEqual(recall, 1)
        self.assertEqual(removed, 1)
        self.assertTrue(result.covers(Point(5, 5)))
        self.assertTrue(result.covers(Point(22, 22)))
        self.assertFalse(result.covers(Point(20, 2)))

    def test_holes_are_not_silently_filled(self):
        mask = np.zeros((20, 20), dtype=bool)
        mask[2:18, 2:18] = True
        mask[7:13, 7:13] = False
        result = geometry(mask, 1, np.array([0, 0]))
        self.assertEqual(len(result.interiors), 1)
        self.assertFalse(result.covers(Point(9, 9)))
        self.assertTrue(result.covers(Point(3, 3)))

    def test_scale_offset_and_xy_order_align_with_source_image(self):
        mask = np.zeros((12, 18), dtype=bool)
        mask[3:8, 6:13] = True
        result = geometry(mask, 3, np.array([120, 450]))
        self.assertEqual(result.bounds, (136.5, 457.5, 157.5, 472.5))
        self.assertTrue(result.covers(Point(120+3*8, 450+3*5)))

    def test_no_prompt_support_does_not_make_an_islet(self):
        mask = np.zeros((15, 15), dtype=bool)
        mask[8:12, 8:12] = True
        kept, recall, removed = supported_components(mask, np.array([[2, 2]]))
        self.assertFalse(kept.any())
        self.assertEqual(recall, 0)
        self.assertEqual(removed, 1)
        self.assertTrue(geometry(kept, 1, np.array([0, 0])).is_empty)

    def test_image_edge_does_not_shift_coordinates(self):
        mask = np.zeros((10, 10), dtype=bool)
        mask[0:4, 0:6] = True
        result = geometry(mask, 1, np.array([0, 0]))
        self.assertTrue(result.covers(Point(0, 0)))
        self.assertEqual(result.bounds, (-.5, -.5, 5.5, 3.5))


if __name__ == '__main__': unittest.main()
