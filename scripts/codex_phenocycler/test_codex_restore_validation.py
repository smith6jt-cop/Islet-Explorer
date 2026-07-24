"""Tests for the CODEX RESTORE level-2 adapters -- the patch contract, pair web, and clip split.

These exercise the module-level config and the `patch_*()` functions WITHOUT touching data or the
11 GB TSV, so they run in ~1 s. The end-to-end data path is covered by
`codex_restore_validation.py --smoke` and the screen runs.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, os.environ.get("PHENOCYCLER_ROOT",
                str(Path(__file__).resolve().parents[2] / "external" / "Phenocycler_Analysis")))

import phenocycler  # noqa: E402,F401  BEFORE numpy
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
from phenocycler import reference_selection as rs  # noqa: E402
from phenocycler import restore_validation as rv  # noqa: E402

import codex_lineage as cl  # noqa: E402
import codex_reference_selection as crs  # noqa: E402
import codex_restore_validation as crv  # noqa: E402


def test_codex_screens_thirty_pairs_cd20_is_a_target():
    assert len(crv.CODEX_PAIRS) == 30
    assert ("CD20", "ECAD") in crv.CODEX_PAIRS, "CD20 is re-included as a B target at level 2"
    assert ("CHGA", "CD31") in crv.CODEX_PAIRS, "endocrine split is carried into the level-2 screen"
    # the three hormone exocrine comparators are extra candidates
    for t in ("INS", "GCG", "SST"):
        assert (t, "panCK") in crv.CODEX_PAIRS


def test_clip_groups_partition_the_pairs():
    assert len(crv.CLIP_0995_PAIRS) == 17
    assert len(crv.CLIP_UNCLIPPED_PAIRS) == 13
    assert set(crv.CLIP_0995_PAIRS).isdisjoint(crv.CLIP_UNCLIPPED_PAIRS)
    assert set(crv.CLIP_0995_PAIRS) | set(crv.CLIP_UNCLIPPED_PAIRS) == set(crv.CODEX_PAIRS)
    # a pair is unclipped iff its TARGET is clip-sensitive
    for t, r in crv.CLIP_UNCLIPPED_PAIRS:
        assert t in cl.CLIP_UNCLIPPED_TARGETS
    for t, r in crv.CLIP_0995_PAIRS:
        assert t not in cl.CLIP_UNCLIPPED_TARGETS


def test_screened_markers_are_twenty_seven_and_all_have_a_floor():
    assert len(crv.SCREENED_MARKERS) == 27
    assert "CD20" in crv.SCREENED_MARKERS
    # the excluded non-lineage markers are absent
    for m in ("BActin", "NaKATPase", "BCatenin", "ColIV", "epCAM", "KRT8-18"):
        assert m not in crv.SCREENED_MARKERS
    for m in crv.SCREENED_MARKERS:
        assert m in cl.MARKER_INPUT_FLOOR_METHODS, f"{m} has no floor method"


def test_restore_excluded_markers_is_empty_cd20_unwound():
    assert cl.RESTORE_EXCLUDED_MARKERS == {}
    assert cl.DEFERRED_EXTRACTION_MARKERS == ()


def test_patch_upstream_sets_the_codex_contract():
    crv.patch_upstream()
    assert set(rv.QUPATH_MARKERS) == set(crv.SCREENED_MARKERS)
    assert all(rv.QUPATH_MARKERS[m] == m for m in crv.SCREENED_MARKERS)
    assert set(rv.MARKER_INPUT_FLOOR_METHODS) == set(rv.QUPATH_MARKERS)
    assert set(crv.CODEX_PAIRS) <= set(rv.CANDIDATE_PAIRS)
    assert ("CD20", "ECAD") in rv.CANDIDATE_PAIRS, "CD20 pair must be admissible after the unwind"


def test_cd20_is_now_admitted_but_self_pairs_still_raise():
    # CD20 as target is now valid (the both-role exclusion is unwound).
    cl.assert_pairs_valid([("CD20", "ECAD")])          # must NOT raise
    with pytest.raises(ValueError):
        cl.assert_pairs_valid([("ECAD", "ECAD")])      # self-pair
    with pytest.raises(ValueError):
        cl.assert_pairs_valid([])                      # empty


def test_cd4_cd8a_take_ecad_not_cd3e():
    # The co-receptor guard is now vacuous (immune -> ECAD) but must still hold as an invariant.
    d = {t: r for t, r in cl.LEVEL2_PAIRS}
    assert d["CD4"] == "ECAD"
    assert d["CD8a"] == "ECAD"
    assert d["CD3e"] == "ECAD"
    cl.no_rare_counterparts()  # must not raise on the configured web


def test_faithful_feature_space_not_cell_only():
    assert rv.PairValidationConfig().feature_compartments == ("Membrane", "Cell")


def test_reference_selection_adapter_rebinds_the_import_copies():
    crs.patch_reference_selection()
    assert set(rs.CANDIDATE_PAIRS) == set(crv.CODEX_PAIRS)            # 30 candidates
    assert set(rs.BASELINE_PAIRS) == set(crs.CODEX_BASELINE_PAIRS)   # 27 one-per-target
    assert len(rs._TARGET_ORDER) == 27, "without _TARGET_ORDER the dossier emits zero figures"
    # every override references a real screened pair
    for tgt, ref in crs.CODEX_REFERENCE_OVERRIDES.items():
        assert (tgt, ref) in set(crv.CODEX_PAIRS)


def test_build_reference_matrix_accepts_a_synthetic_v8_screen():
    crs.patch_reference_selection()
    rows = []
    for donor in ("d1", "d2", "d3"):
        for tgt, ref in (("CD8a", "ECAD"), ("INS", "GCG"), ("INS", "panCK")):
            rows.append({
                "donor": donor, "target": tgt, "reference": ref,
                "reference_n": 5000, "target_n": 100,
                "reference_to_target_ratio": 50.0, "reference_undersized": False,
                "candidate_full_maximum": 1234.0, "reference_fold": 20.0,
                "threshold_positive_n": 90, "n_assigned": 10000,
            })
    matrix = rs.build_reference_matrix(pd.DataFrame(rows), overrides=crs.CODEX_REFERENCE_OVERRIDES)
    assert {"CD8a", "INS"} <= set(matrix["target"])
    # INS has two candidate references in the synthetic screen
    assert set(matrix[matrix["target"] == "INS"]["reference"]) == {"GCG", "panCK"}


def test_review_bundle_patches_apply_per_group():
    import codex_restore_pair_review as crpr
    from phenocycler import restore_pair_review as rpr
    # clipped group
    crpr.patch_review(crpr.CODEX_SHORTLIST_CLIPPED)
    assert rpr.SHORTLIST_PAIRS == crpr.CODEX_SHORTLIST_CLIPPED
    assert ("CD20", "ECAD") in rpr.SHORTLIST_PAIRS, "CD20 is now a reviewable target"
    assert "ECAD" in rpr.MARKER_INPUT_FLOOR_METHODS
    assert rpr._qptiff_channel_names.__name__ == "_codex_qptiff_channel_names"
    # unclipped group
    crpr.patch_review(crpr.CODEX_SHORTLIST_UNCLIPPED)
    assert rpr.SHORTLIST_PAIRS == crpr.CODEX_SHORTLIST_UNCLIPPED
    assert ("HLADR", "ECAD") in rpr.SHORTLIST_PAIRS
    # every shortlisted pair has a rationale
    for p in crpr.CODEX_SHORTLIST_CLIPPED + crpr.CODEX_SHORTLIST_UNCLIPPED:
        assert crpr.CODEX_PAIR_RATIONALE.get(p), f"no rationale for {p}"
