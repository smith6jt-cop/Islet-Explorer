"""Geometry regressions for the content-rule tissue continuation."""
import json
import sys
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import MultiPolygon, box, mapping

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'uniform'))
import redraw_tissue as RT


def test_original_footprints_survive_edge_crop_and_internal_holes():
    tissue = box(10, 10, 100, 100).difference(box(40, 40, 60, 60))
    cell = box(5, 25, 15, 35)
    nucleus = box(4, 28, 8, 32)
    endocrine = box(45, 45, 50, 50)
    kept, counts = RT.preserve_exact_footprints(tissue, [([cell, nucleus, endocrine], ['cell', 'nucleus', 'endocrine'])])
    assert shapely.covers(kept, [cell, nucleus, endocrine]).all()
    assert kept.equals(shapely.union_all([tissue, cell, nucleus, endocrine]))
    assert not kept.intersects(box(52, 52, 58, 58))
    assert not kept.intersects(box(0, 40, 9, 50))
    assert [counts[k] for k in ('cell', 'nucleus', 'endocrine')] == [1, 1, 1]


def test_content_geometry_preservation_is_idempotent():
    geom = box(10, 10, 30, 30)
    original = box(2.21, 3.98, 15.73, 16.42)
    first, _ = RT.preserve_exact_footprints(geom, [([original], ['cell'])])
    second, counts = RT.preserve_exact_footprints(first, [([original], ['cell'])])
    assert second.equals_exact(first, 0)
    assert counts['cell'] == 0
    assert counts['added_area_um2'] == 0


def test_original_nucleus_geometry_and_multipart_endocrine_are_loaded(tmp_path):
    for name in ('cell_geojson', 'islets'):
        (tmp_path / name).mkdir()
    cell, nucleus = box(0, 0, 5, 5), box(-1, 1, 1, 2)
    endocrine = MultiPolygon([box(20, 20, 25, 25), box(30, 30, 40, 40)])
    (tmp_path / 'cell_geojson' / 'example.geojson').write_text(json.dumps({'features': [
        {'geometry': mapping(cell), 'nucleusGeometry': mapping(nucleus)}]}))
    (tmp_path / 'islets' / 'example.geojson').write_text(json.dumps({'features': [{'geometry': mapping(endocrine)}]}))
    batches = list(RT.footprint_batches(tmp_path, 'example', batch_size=1))
    assert [k[0] for _, k in batches] == ['cell', 'nucleus', 'endocrine']
    assert all(a[0].equals_exact(b, 0) for (a, _), b in zip(batches, (cell, nucleus, endocrine)))


def test_geometry_conversion_preserves_large_empty_hole():
    mask = np.zeros((40, 40), dtype=bool)
    mask[3:37, 3:37] = True
    mask[12:28, 12:28] = False
    geom = RT.to_geometry(mask, (320, 320))
    assert geom.is_valid
    assert not geom.covers(box(120, 120, 180, 180))
    assert geom.covers(box(40, 40, 75, 75))


def test_zero_area_source_nucleus_is_recorded_without_origin_tissue(tmp_path):
    for name in ('cell_geojson', 'islets'):
        (tmp_path / name).mkdir()
    cell = box(20, 20, 25, 25)
    (tmp_path / 'cell_geojson' / 'example.geojson').write_text(json.dumps({'features': [{
        'id': 'sentinel', 'geometry': mapping(cell),
        'nucleusGeometry': {'type': 'Polygon', 'coordinates': [[[0, 0]] * 5]}}]}))
    (tmp_path / 'islets' / 'example.geojson').write_text('{"features": []}')
    diagnostics = {}
    batches = list(RT.footprint_batches(tmp_path, 'example', diagnostics=diagnostics))
    kept, _ = RT.preserve_exact_footprints(box(15, 15, 30, 30), batches)
    assert kept.covers(cell)
    assert not kept.covers(shapely.Point(0, 0))
    assert diagnostics['empty_source_footprints'] == [
        {'object_id': 'sentinel', 'kind': 'nucleus', 'geometry_type': 'Polygon'}]


def test_content_rule_rejects_smooth_lumen_and_keeps_thick_bright_wall():
    # A large bright wall encloses smooth haze. Its width exceeds the local
    # opening window, so contrast alone would incorrectly remove its interior.
    yy, xx = np.mgrid[:160, :160]
    radius = np.hypot(xx - 80, yy - 80)
    stroma = np.full((160, 160), 0.1)
    stroma[radius < 50] = 0.35
    stroma[(radius > 20) & (radius < 50)] = 0.95
    dapi = np.full((160, 160), 50.0)
    dapi[(radius > 20) & (radius < 50)] = 1000
    dapi[radius < 20] = 300
    diagnostics = {}
    mask = RT.content_mask(dapi, stroma, diagnostics)
    assert not mask[80, 80]
    assert mask[80, 115]
    assert diagnostics['contrast_opening_pixels'] == 17
    assert diagnostics['rule_version'] == RT.CONTENT_RULE_VERSION


def test_positive_otsu_handles_no_local_contrast():
    assert np.isinf(RT.positive_otsu(np.zeros((10, 10))))


def test_large_lumen_survives_footprint_padding_below_small_hole_threshold():
    content = np.zeros((40, 40), bool)
    content[3:37, 3:37] = True
    content[12:24, 12:24] = False  # 144 pixels exceeds the 109-pixel cleanup threshold.
    footprint = np.zeros_like(content)
    footprint[12:17, 12:24] = True  # Remaining 84 pixels would otherwise be filled.
    diagnostics = {}
    result = RT.add_footprints_preserving_content_holes(content, footprint, diagnostics)
    np.testing.assert_array_equal(result, content | footprint)
    assert result[footprint].all()
    assert not result[17:24, 12:24].any()
    assert diagnostics['protected_pre_padding_holes'] == 1
    assert diagnostics['restored_empty_pixels_after_second_fill'] == 84


def test_small_new_hole_between_padded_footprints_still_filled():
    content = np.zeros((30, 30), bool)
    footprint = np.zeros_like(content)
    footprint[8:22, 8:22] = True
    footprint[12:18, 12:18] = False  # Not a large hole in the pre-padding content.
    result = RT.add_footprints_preserving_content_holes(content, footprint)
    assert result[8:22, 8:22].all()
    assert not result[:8].any()


def test_protected_hole_rule_requires_completed_first_cleanup():
    import pytest
    content = np.ones((20, 20), bool)
    content[8:12, 8:12] = False
    with pytest.raises(ValueError, match='below the first cleanup threshold'):
        RT.add_footprints_preserving_content_holes(content, np.zeros_like(content))


def test_corner_touching_contours_keep_shell_hole_hierarchy(monkeypatch):
    # Default make_valid(linework) fills seven background pixels in this
    # topology, reproducing the large false septum found during image review.
    rows = ['1111110101', '1011110100', '1111110010', '1110001100',
            '0111011101', '0011101110', '0101101110', '1110010101',
            '0001010101', '0011000101']
    mask = np.array([[v == '1' for v in row] for row in rows])
    monkeypatch.setattr(RT, 'EDGE_MARGIN_UM', 0)
    diagnostics = {}
    geom = RT.to_geometry(mask, (80, 80), diagnostics)
    yy, xx = np.indices(mask.shape)
    back = shapely.covers(geom, shapely.points((xx.ravel() + .5) * 8,
                                             (yy.ravel() + .5) * 8)).reshape(mask.shape)
    np.testing.assert_array_equal(back, mask)
    assert diagnostics['repair_method'] == 'structure'
    assert diagnostics['pre_crop_added_pixel_centers'] == 0
    assert diagnostics['pre_crop_removed_pixel_centers'] == 0
    assert diagnostics['shapely_version'] == shapely.__version__
