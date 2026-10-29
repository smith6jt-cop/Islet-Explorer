"""Diagnostic cell outlines must show the entire cell used for image review."""
import sys
from pathlib import Path

import cv2
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'research'))
import typing_accuracy_tiles_20261029 as T


def ring(x0, y0, x1, y1):
    return [[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]


def multipart():
    return {'type': 'MultiPolygon', 'coordinates': [
        [ring(12, 12, 16, 16)],  # Island listed before the polygon containing its hole.
        [ring(2, 2, 26, 26), ring(8, 8, 20, 20)],
        [ring(30, 3, 38, 13)]]}


def test_disconnected_components_hole_and_island_have_order_independent_cell_masks():
    components = T.polygon_components(multipart())
    own, interior = T.cell_masks(components, (48, 48), (0, 0))
    expected = np.zeros((48, 48), np.uint8)
    expected[2:27, 2:27] = 1
    expected[8:21, 8:21] = 0
    expected[12:17, 12:17] = 1
    expected[3:14, 30:39] = 1
    np.testing.assert_array_equal(own, expected)
    np.testing.assert_array_equal(T.cell_masks(components[::-1], (48, 48), (0, 0))[0], expected)
    assert T.polygon_bounds(components) == (2, 2, 38, 26)
    assert interior[5, 5] and interior[14, 14] and interior[7, 34]
    assert not interior[7, 7] and not interior[13, 13] and not interior[7, 30]
    assert not interior[9, 9] and not interior[4, 28]


def test_all_exteriors_and_hole_boundaries_are_drawn_after_crop_translation():
    components = T.polygon_components(multipart())
    image = np.zeros((48, 48, 3), np.uint8)
    T.draw_cell_outline(image, components, (-3, -4), (255, 255, 255))
    for x, y in [(2, 5), (8, 11), (12, 14), (38, 7)]:
        assert (image[y + 4, x + 3] == 255).all()
    assert not image[9 + 4, 9 + 3].any()  # Empty hole interior is not an outline.


def test_simple_polygon_display_mask_and_two_pixel_interior_remain_unchanged():
    coords = ring(4.2, 5.2, 17.2, 20.2)
    components = T.polygon_components({'type': 'Polygon', 'coordinates': [coords]})
    own, interior = T.cell_masks(components, (30, 30), (1, 2))
    legacy = np.zeros((30, 30), np.uint8)
    cv2.fillPoly(legacy, [np.round(np.array(coords) - [1, 2]).astype(np.int32)], 1)
    np.testing.assert_array_equal(own, legacy)
    np.testing.assert_array_equal(interior, cv2.erode(legacy, np.ones((3, 3), np.uint8), iterations=2).astype(bool))
