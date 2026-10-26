"""Physical geometry and counting contracts for the analysis-only H1 correction."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest
import shapely
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'research'))
from h1_edge_eligibility_20261024 import border_contact, label_preclip, component_ids, neighborhood_support, restrict_to_tissue, damaged_member_mask
from h1_features_20261022 import exclusive_peri_owners
from h1_analysis_20261022 import eligible


@pytest.mark.parametrize('g', [box(0, 40, 12, 50), box(40, 0, 50, 12), box(88, 40, 100, 50),
                              box(40, 88, 50, 100), box(0, 0, 10, 10), box(90, 90, 100, 100),
                              box(0, 90, 10, 100), box(90, 0, 100, 10)])
def test_four_edges_and_corners(g):
    assert border_contact(g, 100, 100)
    assert 0 < g.centroid.x < 100 and 0 < g.centroid.y < 100


def test_nucleus_contact_and_missing_geometry():
    assert not border_contact(box(5, 5, 25, 25), 100, 100)
    assert border_contact(box(1, 10, 12, 12), 100, 100)
    assert border_contact(None, 100, 100) is None


def test_source_qc_missing_nucleus_is_not_evidence_of_truncation():
    d = pd.DataFrame(dict(cell_edge=[False, True, False, False],
                          nucleus_edge=[None, None, None, False],
                          geometry_missing=[True, True, True, False],
                          qc_analysis_eligible=[False, False, True, True]))
    assert damaged_member_mask(d).tolist() == [False, True, True, False]


def test_preclip_component_propagates_to_interior_split_and_merge():
    mask = np.zeros((100, 100), bool)
    mask[30:50, :40] = True
    mask[60:70, 60:70] = True
    labels, edge = label_preclip(mask, 1, 100, 100, 1)
    # Both children are interior after clip/erosion, but have a cut ancestor.
    for child in (box(12, 32, 20, 40), box(25, 40, 35, 48)):
        assert not border_contact(child, 100, 100)
        assert edge[component_ids(child, labels, 1)].any()
    interior = box(61, 61, 68, 68)
    assert not edge[component_ids(interior, labels, 1)].any()
    merged = shapely.union_all([interior, box(25, 40, 35, 48)])
    assert edge[component_ids(merged, labels, 1)].any()


def test_inset_clipping_is_flagged_even_without_physical_contact():
    mask = np.zeros((100, 100), bool)
    mask[20:40, 4:30] = True
    labels, edge = label_preclip(mask, 1, 100, 100, 1)
    assert edge[component_ids(box(12, 25, 20, 35), labels, 1)].any()


def test_outside_tissue_bridge_does_not_make_interior_object_edge_connected():
    mask = np.zeros((100, 100), bool)
    mask[40:50, :70] = True
    tissue = shapely.union_all([box(0, 30, 20, 60), box(50, 30, 80, 60)])
    restricted = restrict_to_tissue(mask, tissue, 1)
    labels, edge = label_preclip(restricted, 1, 100, 100, 1)
    assert edge[component_ids(box(10, 42, 15, 48), labels, 1)].any()
    assert not edge[component_ids(box(55, 42, 65, 48), labels, 1)].any()


def test_neighborhood_support_is_separate_from_core():
    core, acquisition = box(30, 30, 40, 40), box(0, 0, 100, 100)
    assert not border_contact(core, 100, 100)
    complete, requested, observed = neighborhood_support(core, acquisition, 20)
    assert complete and requested == pytest.approx(observed)
    complete, requested, observed = neighborhood_support(core, acquisition, 100)
    assert not complete and observed < requested


def test_exclusive_owner_stable_on_ties_and_uses_only_valid_memberships():
    points = shapely.points([[15, 5], [25, 5]])
    geoms = np.array([box(0, 0, 10, 10), box(20, 0, 30, 10)], dtype=object)
    names = ['B', 'A']
    ids = [np.array([0]), np.array([0, 1])]
    assert exclusive_peri_owners(points, geoms, names, ids).tolist() == [1, 1]
    assert exclusive_peri_owners(points, geoms, names, ids, [True, False]).tolist() == [0, -1]
    reverse = exclusive_peri_owners(points, geoms[::-1], names[::-1], ids[::-1])
    assert np.array(names[::-1])[reverse].tolist() == ['A', 'A']


def test_partial_neighborhood_and_incomplete_core_cannot_enter_h1():
    d = pd.DataFrame(dict(case_id=['1']*3, n_immune_nearby=[2]*3, core_area_um2=[100]*3,
        peri_area_um2=[200]*3, field20_100_area_um2=[2000]*3, immune_density_field20_100_mm2=[100]*3,
        immune_localization_fraction=[.5]*3, n_beta_core=[1]*3, n_immune_core=[1]*3,
        n_endocrine_core=[2]*3, core_complete=[True, False, True], neighborhood_complete_20=[True]*3,
        neighborhood_complete_100=[True, True, False], neighborhood_complete_250=[False]*3))
    assert eligible(d).index.tolist() == [0]
    d['core_complete'] = d.core_complete.astype('boolean')
    d.loc[0, 'core_complete'] = None
    with pytest.raises(ValueError, match='Incomplete acquisition'):
        eligible(d)
