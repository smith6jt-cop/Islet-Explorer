"""Validate codex_panel against the ACTUAL image and QuPath artifacts, not fixtures.

The failure mode this guards against is silent: a marker resolved to the wrong channel index, or
a donor id that does not join, produces plausible numbers rather than an error. So these tests
read the real OME-XML sidecars and the real QuPath TSV header.

    conda activate scvi-env && python -m pytest scripts/codex_phenocycler/test_codex_panel.py -q

Tests that need the large artifacts skip cleanly when they are absent (e.g. in CI).
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from codex_panel import (
    CHANNELS_33,
    CHANNELS_35,
    DONORS,
    NON_PANEL_CHANNELS,
    PANEL_31,
    channel_index,
    expected_channels,
    normalize_donor_id,
    ome_channel_names,
    panel_markers,
    qupath_column,
    to_panel,
)

REPO = Path(__file__).resolve().parents[2]
IMAGES = REPO / "app" / "shiny_app" / "www" / "local_images"
CELLS_TSV = Path("/home/smith6jt/panc_CODEX/results/old_data_plots/all_data.tsv")


# --- pure logic -------------------------------------------------------------------------------

def test_35ch_is_a_strict_prefix_extension_of_33ch():
    """The only reason index 0..32 is stable across donors. If this breaks, so does everything."""
    assert CHANNELS_35[:33] == CHANNELS_33
    assert len(CHANNELS_35) == 35


def test_every_panel_marker_is_an_acquired_channel():
    assert set(PANEL_31) <= set(CHANNELS_35)


def test_non_panel_channels_are_exactly_the_leftovers():
    assert set(CHANNELS_35) - set(PANEL_31) == set(NON_PANEL_CHANNELS)


def test_donor_id_normalization_strips_qupath_zero_padding():
    assert normalize_donor_id("0112") == "112"
    assert normalize_donor_id("6533") == "6533"
    assert normalize_donor_id("ND_0112") == "112"
    with pytest.raises(ValueError):
        normalize_donor_id("no-digits")


def test_channel_index_raises_rather_than_guessing():
    ch33 = expected_channels("6533")
    assert channel_index("PGP9.5", ch33) == ch33.index("PGP9.5")
    # epCAM exists on 35-ch donors only -- must raise for a 33-ch donor, never return a default
    with pytest.raises(KeyError):
        channel_index("epCAM", ch33)


def test_ki67_is_the_only_median_export():
    assert qupath_column("Ki67") == "Cell: Ki67: Median"
    assert qupath_column("INS") == "Cell: INS: Mean"


def test_panel_markers_is_donor_aware_but_panel_is_fully_covered_either_way():
    # the two extra channels are not panel markers, so every donor carries all 31
    for d in ("112", "6356"):
        assert panel_markers(d) == PANEL_31


def test_name_normalization():
    assert to_panel("PGP9_5") == "PGP9.5"
    assert to_panel("KRT8_18") == "KRT8-18"
    assert to_panel("INS") == "INS"


# --- against the real artifacts ---------------------------------------------------------------

@pytest.mark.skipif(not IMAGES.is_dir(), reason="local_images not present")
def test_ome_channel_names_match_the_declared_split_for_every_donor():
    sidecars = sorted(IMAGES.glob("*.ome.xml"))
    assert len(sidecars) == 15, f"expected 15 donors, found {len(sidecars)}"
    seen = set()
    for xml in sidecars:
        donor = normalize_donor_id(xml.name.replace(".ome.xml", "").split("_")[-1])
        seen.add(donor)
        names = ome_channel_names(xml)
        assert names == expected_channels(donor), (
            f"{donor}: OME-XML channel order diverged from codex_panel. Got {names}"
        )
    assert seen == set(DONORS)


@pytest.mark.skipif(not CELLS_TSV.exists(), reason="QuPath per-cell TSV not present")
def test_qupath_tsv_carries_a_column_for_every_panel_marker():
    header = CELLS_TSV.open(encoding="utf-8").readline().rstrip("\n").split("\t")
    cols = set(header)
    missing = [m for m in PANEL_31 if qupath_column(m) not in cols]
    assert not missing, f"panel markers with no QuPath column: {missing}"
    # Ki67's stat really is Median in the export -- guards the special case above
    assert "Cell: Ki67: Median" in cols
    assert "Cell: Ki67: Mean" not in cols


@pytest.mark.skipif(not CELLS_TSV.exists(), reason="QuPath per-cell TSV not present")
def test_qupath_marker_column_order_matches_acquisition_order():
    """A positional read of the TSV must land on the same markers as a positional image read."""
    header = CELLS_TSV.open(encoding="utf-8").readline().rstrip("\n").split("\t")
    order = [
        to_panel(c.split(": ")[1])
        for c in header
        if c.startswith("Cell: ") and c.count(": ") == 2
        and c.rsplit(": ", 1)[1] in ("Mean", "Median")
    ]
    assert order == CHANNELS_35, f"TSV marker order diverged: {order}"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
