"""Native coordinate precision must not excuse actual overlapping regions."""
import sys
from pathlib import Path
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).parents[1] / 'cohort_expansion'))
from verify_size_adaptive_output import overlap_exceeds_rounding


def test_rounding_strip_has_no_interior_beyond_two_edge_errors():
    overlap = box(0, 0, 10, 10).intersection(box(9.999, 0, 20, 10))
    assert not overlap_exceeds_rounding(overlap, .02)


def test_actual_overlap_is_rejected_even_when_small_relative_to_islet():
    overlap = box(0, 0, 10, 10).intersection(box(9.9, 0, 20, 10))
    assert overlap_exceeds_rounding(overlap, .02)
