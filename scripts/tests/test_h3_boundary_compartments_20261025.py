"""Check spatial compartment assignment against simple known geometry."""
from pathlib import Path
import sys

import numpy as np
import shapely
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "research"))
import h3_boundary_compartments_20261025 as h3


def fixture_geometry():
    core = box(0, 0, 10, 10)
    vessels = [box(-2, 4, 2, 6), box(3, 3, 4, 4), box(12, 4, 13, 6)]
    cross, n = h3.crossing_geometry(core, vessels, shapely.STRtree(vessels))
    return core, cross, n


def test_only_boundary_spanning_capillary_is_a_crossing():
    _, cross, n = fixture_geometry()
    assert n == 1
    assert cross.length == 2


def test_known_cells_and_areas_partition_peripheral_band():
    core, cross, _ = fixture_geometry()
    result = h3.measure_compartments(core, cross, 1,
                                   np.array([[.5, 5], [5, .5], [5, 5], [-1, 5]]))
    assert result["n_near"] == 1
    assert result["n_reference"] == 1
    assert result["n_band"] == 2
    assert abs(result["near_area_um2"] + result["reference_area_um2"] - 36) < 1e-8


def test_zero_immune_objects_are_included():
    core, cross, _ = fixture_geometry()
    result = h3.measure_compartments(core, cross, 1, np.empty((0, 2)))
    assert result is not None
    assert result["n_near"] == result["n_reference"] == 0


def test_no_reference_area_is_not_an_observed_zero_density():
    core, cross, _ = fixture_geometry()
    assert h3.measure_compartments(core, cross, 100, np.empty((0, 2))) is None
