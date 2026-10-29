"""Empty-slide calibration must never sample a Tissue island inside a hole."""
import json
from pathlib import Path
import sys

import numpy as np
import pytest
from shapely.geometry import MultiPolygon, Polygon, box, mapping

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'uniform' / 'typing'))
import redsea_cells as RS


@pytest.mark.parametrize('split_features', [False, True])
def test_empty_windows_union_is_independent_of_polygon_and_feature_order(tmp_path, monkeypatch, split_features):
    monkeypatch.setattr(RS, 'EMPTY', {'window_px': 1, 'windows': 20000, 'margin_um': 0., 'seed': 0})
    donut = Polygon(box(5, 5, 95, 95).exterior.coords, holes=[box(20, 20, 80, 80).exterior.coords])
    island = box(30, 30, 50, 50)
    path = tmp_path / 'inputs' / 'tissue' / 'image.geojson'
    path.parent.mkdir(parents=True)
    outputs = []
    for polygons in ([island, donut], [donut, island]):
        geometries = polygons if split_features else [MultiPolygon(polygons)]
        path.write_text(json.dumps({'features': [{'geometry': mapping(g)} for g in geometries]}))
        ys, xs = RS.empty_windows(tmp_path, 'image', 101, 101)
        outputs.append((ys, xs))
        found = np.zeros((101, 101), bool)
        found[ys, xs] = True
        expected = np.ones_like(found)
        expected[5:96, 5:96] = False
        expected[20:81, 20:81] = True
        expected[30:51, 30:51] = False
        expected[[0, -1], :] = False
        expected[:, [0, -1]] = False
        np.testing.assert_array_equal(found, expected)
    np.testing.assert_array_equal(outputs[0], outputs[1])


def test_empty_windows_keeps_margin_and_window_fully_outside_island(tmp_path, monkeypatch):
    monkeypatch.setattr(RS, 'EMPTY', {'window_px': 3, 'windows': 20000,
                                   'margin_um': 4 * RS.CT.PIXEL_UM, 'seed': 0})
    donut = Polygon(box(-5, -5, 105, 105).exterior.coords, holes=[box(10, 10, 90, 90).exterior.coords])
    path = tmp_path / 'inputs' / 'tissue' / 'image.geojson'
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({'features': [{'geometry': mapping(box(40, 40, 60, 60))},
                                             {'geometry': mapping(donut)}]}))
    ys, xs = RS.empty_windows(tmp_path, 'image', 101, 101)
    found = np.zeros((101, 101), bool)
    found[ys, xs] = True
    assert found[25, 25] and not found[50, 50]
    assert not found[:15].any() and not found[:, :15].any()
    assert not found[86:].any() and not found[:, 86:].any()
    assert not found[50, 35:66].any() and not found[35:66, 50].any()
