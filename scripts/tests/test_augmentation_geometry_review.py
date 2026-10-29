"""Exact containment and native-resolution field selection for augmentation review."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import shapely
from shapely.geometry import box
import tifffile

PATH = Path(__file__).resolve().parents[1] / 'research/augmentation_geometry_review_20261029.py'
spec = importlib.util.spec_from_file_location('augmentation_review', PATH)
GR = importlib.util.module_from_spec(spec)
spec.loader.exec_module(GR)


def fixture_geometry():
    ids = ['old', 'new']
    cells = np.array([box(1, 1, 3, 3), box(11, 1, 13, 3)], object)
    nuclei = np.array([shapely.GeometryCollection(), box(11.5, 1.5, 12.5, 2.5)], object)
    frame = pd.DataFrame({'Object ID': ids, 'Centroid X µm': np.array([2, 12]) * GR.CT.PIXEL_UM,
                          'Centroid Y µm': np.array([2, 2]) * GR.CT.PIXEL_UM})
    return frame, ids, cells, nuclei, {'new'}, box(0, 0, 10, 10), box(0, 0, 20, 10), []


def test_exact_containment_and_zero_area_nuclei_are_distinct():
    report, detail = GR.check_geometry(*fixture_geometry())
    assert report['geometry_checks_passed']
    assert report['new_centroids_outside_added'] == 0
    assert report['old_empty_or_zero_area_nuclei'] == 1
    assert report['new_positive_area_nuclei'] == 1
    assert report['visual_review_status'] == 'pending'
    assert list(detail.loc[detail['is_new'], 'object_id']) == ['new']


def test_centroid_inside_does_not_hide_cell_or_nucleus_area_loss():
    args = list(fixture_geometry())
    args[6] = box(0, 0, 12.1, 10)
    report, detail = GR.check_geometry(*args)
    assert report['new_centroids_outside_final'] == 0
    assert report['new_cell_positive_area_outside_final'] == 1
    assert report['new_nucleus_positive_area_outside_final'] == 1
    assert not report['geometry_checks_passed']


def test_new_centroid_must_be_in_added_area_even_if_full_cell_is_in_final():
    args = list(fixture_geometry())
    args[2][1] = box(7, 1, 9, 3)
    args[3][1] = box(7.5, 1.5, 8.5, 2.5)
    report, _ = GR.check_geometry(*args)
    assert report['new_centroids_outside_added'] == 1
    assert report['new_cell_positive_area_outside_final'] == 0
    assert not report['geometry_checks_passed']


def test_export_rounding_cannot_create_false_roi_centroid_failure():
    args = list(fixture_geometry())
    args[0].loc[1, 'Centroid X µm'] = 9.99 * GR.CT.PIXEL_UM
    report, detail = GR.check_geometry(*args)
    assert report['new_centroids_outside_added'] == 0
    assert report['new_exported_rounded_centroids_outside_added'] == 1
    assert report['geometry_checks_passed']
    assert detail.loc[1, 'centroid_x_um'] == 12 * GR.CT.PIXEL_UM


def test_nucleus_extending_outside_own_cell_is_reported():
    args = list(fixture_geometry())
    args[3][1] = box(11, 1, 14, 3)
    report, _ = GR.check_geometry(*args)
    assert report['new_nucleus_positive_area_outside_final'] == 0
    assert report['new_nucleus_positive_area_outside_own_cell'] == 1
    assert not report['geometry_checks_passed']


def test_seeded_fields_cover_space_without_repeating_cells():
    detail = pd.DataFrame({'object_id': [f'cell{i:02d}' for i in range(20)], 'is_new': True,
                           'centroid_x_um': np.arange(20) * 100., 'centroid_y_um': 10.})
    first = GR.choose_fields(detail, seed=29)
    assert first.equals(GR.choose_fields(detail.sample(frac=1, random_state=7), seed=29))
    assert len(first) == first['object_id'].nunique() == 6
    assert first['centroid_x_um'].max() - first['centroid_x_um'].min() >= 1800
    assert GR.choose_fields(detail.assign(is_new=False), seed=29).empty


def test_tiles_use_native_pixels_not_enlarged_pyramid(tmp_path):
    path = tmp_path / 'tiny.ome.tiff'
    native = np.arange(3 * 16 * 24, dtype=np.uint16).reshape(3, 16, 24)
    low = native[:, ::2, ::2].copy()
    # Native pixels deliberately differ from a replicated low-resolution display.
    native[:, 1::2, 1::2] = 0
    with tifffile.TiffWriter(path, ome=True) as writer:
        writer.write(native, subifds=1, photometric='minisblack',
                     metadata={'axes': 'CYX', 'Channel': {'Name': list(GR.CHANNELS)}})
        writer.write(low, subfiletype=1, photometric='minisblack')
    selected = pd.DataFrame({'object_id': ['new'], 'centroid_x_um': [20 * GR.CT.PIXEL_UM],
                              'centroid_y_um': [14 * GR.CT.PIXEL_UM]})
    whole, crops, windows, meta = GR.native_views(path, selected, tile_um=3, whole_level=1)
    assert whole.shape == (8, 12, 3)
    assert crops[0].shape == (8, 8, 3)
    assert windows[0]['pyramid_level'] == 0
    x, y = windows[0]['x_px'], windows[0]['y_px']
    expected = GR.scale(native[0, y:y + 8, x:x + 8], meta['display_limits_log1p']['DAPI'])
    np.testing.assert_array_equal(crops[0][:, :, 2], expected)
    assert meta['whole_pyramid_level'] == 1
    assert meta['channels_rgb'] == {'red': 'ECAD', 'green': 'BActin', 'blue': 'DAPI'}


def test_existing_review_is_never_overwritten(tmp_path):
    folder = tmp_path / 'reviews/A'
    folder.mkdir(parents=True)
    saved = folder / 'review.json'
    saved.write_text('already reviewed')
    with pytest.raises(ValueError, match='Refusing to overwrite'):
        GR.render('A', tmp_path / 'not-even-a-run', tmp_path / 'reviews')
    assert saved.read_text() == 'already reviewed'
