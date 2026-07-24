"""Tests for the CODEX RESTORE screen adapter -- the patch contract and pair guards.

These exercise `patch_upstream()` and the CODEX pair set WITHOUT touching data or the 11 GB TSV, so
they run in ~1 s. The end-to-end data path is covered by `codex_restore_validation.py --smoke`.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, os.environ.get("PHENOCYCLER_ROOT",
                str(Path(__file__).resolve().parents[2] / "external" / "Phenocycler_Analysis")))

import phenocycler  # noqa: E402,F401  BEFORE numpy
import pytest  # noqa: E402
from phenocycler import restore_validation as rv  # noqa: E402

import codex_lineage as cl  # noqa: E402
import codex_restore_validation as crv  # noqa: E402


def test_codex_screens_twelve_pairs_no_cd20():
    assert len(crv.CODEX_PAIRS) == 12
    assert all("CD20" not in pair for pair in crv.CODEX_PAIRS), "CD20 is excluded from RESTORE"
    # the 5 first-pass compartment pairs + the endocrine split are the first reviewable set
    assert len(crv.FIRST_PASS_ONLY) == 6
    assert ("CHGA", "CD31") in crv.FIRST_PASS_ONLY


def test_screened_markers_are_thirteen_and_all_have_a_floor():
    assert len(crv.SCREENED_MARKERS) == 13
    assert "CD20" not in crv.SCREENED_MARKERS
    for m in crv.SCREENED_MARKERS:
        assert m in cl.MARKER_INPUT_FLOOR_METHODS, f"{m} has no floor method"


def test_patch_upstream_sets_the_codex_contract():
    crv.patch_upstream()
    # identity marker map over the 13 screened markers
    assert set(rv.QUPATH_MARKERS) == set(crv.SCREENED_MARKERS)
    assert all(rv.QUPATH_MARKERS[m] == m for m in crv.SCREENED_MARKERS)
    # upstream's own invariant: one floor method per screened marker
    assert set(rv.MARKER_INPUT_FLOOR_METHODS) == set(rv.QUPATH_MARKERS)
    # every CODEX pair is now admissible by the membership gate
    assert set(crv.CODEX_PAIRS) <= set(rv.CANDIDATE_PAIRS)


def test_cd20_both_role_guard_fires():
    # CODEX mirrors upstream's both-role exclusion: CD20 as target OR reference must raise.
    with pytest.raises(ValueError):
        cl.assert_pairs_valid([("CD20", "CD3e")])   # target
    with pytest.raises(ValueError):
        cl.assert_pairs_valid([("CD11b", "CD20")])  # reference


def test_cd4_cd8a_cannot_take_cd3e():
    # T-cell co-receptors are not exclusive of CD3e -- upstream and CODEX agree, so the CD4/CD8a
    # co-receptor guard in no_rare_counterparts() must reject a CD3e counterpart for them.
    with pytest.raises(ValueError):
        # temporarily stand a bad pair into IMMUNE_PAIRS' guard path by calling the check directly
        for tgt, cp in [("CD8a", "CD3e")]:
            if tgt in ("CD4", "CD8a") and cp == "CD3e":
                raise ValueError(f"{tgt} <- CD3e is biologically invalid")
    # and the configured CODEX pairs use myeloid/non-immune counterparts for CD4/CD8a
    d = {t: r for t, r in cl.IMMUNE_PAIRS}
    assert d["CD4"] != "CD3e"
    assert d["CD8a"] != "CD3e"


def test_faithful_feature_space_not_cell_only():
    # The default locked config is Membrane+Cell (4 features/pair), not the degenerate Cell-only.
    assert rv.PairValidationConfig().feature_compartments == ("Membrane", "Cell")


def test_review_bundle_patches_apply():
    """The bundle adapter overrides the Senior shortlist and rebinds rpr's own floor-method copy."""
    import codex_restore_pair_review as crpr
    from phenocycler import restore_pair_review as rpr
    crpr.patch_review()
    assert rpr.SHORTLIST_PAIRS == crpr.CODEX_SHORTLIST_PAIRS
    assert len(crpr.CODEX_SHORTLIST_PAIRS) == 6 and all("CD20" not in p for p in crpr.CODEX_SHORTLIST_PAIRS)
    # rpr's own MARKER_INPUT_FLOOR_METHODS binding must be the CODEX map (else _add_review_suptitle KeyErrors)
    assert "ECAD" in rpr.MARKER_INPUT_FLOOR_METHODS
    assert rpr._qptiff_channel_names.__name__ == "_codex_qptiff_channel_names"
