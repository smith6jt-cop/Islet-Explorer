"""CODEX 31-plex panel definition — the single source of truth for the Phenocycler port.

WHY THIS EXISTS
---------------
`Phenocycler_Analysis` was written for the Senior 59-plex PhenoCycler cohort. Three CODEX
realities break its assumptions, and all three are handled here so nothing downstream has to
special-case them:

1. **Channel count is not uniform.** 8 of 15 CODEX donors have 35 channels, 7 have 33. The two
   extras (`epCAM`, `KRT8_18`) are APPENDED, so the 35-channel order is a strict superset of the
   33-channel order — verified across all 15 `www/local_images/*.ome.xml`. Consequence: index
   0..32 happens to be stable today, but **never key a marker by index**. Use `channel_index()`.

2. **`redsea.qptiff_channels()` cannot read CODEX channel names.** It regexes each TIFF page
   description for Akoya `<Biomarker>` / `<Name>` *elements*. CODEX OME-TIFFs store the name as
   an XML *attribute* (`<Channel ... Name="HLADR">`), so every channel would come back `ch0..chN`
   and silently fail to join the marker panel. `ome_channel_names()` is the adapter.

3. **Two image copies disagree.** `/data/QuPath_Projects/.../Images/` has per-donor marker-name
   variants (SMA/SMActin, SST/SOM, LGALS3/MAC2, CHGA/CgA, ...), a HLADR/PDPN order swap in the 11
   plain ImageJ `.tif` entries, and labels the KRT14 slot `TRAM2.tif` on donors 6521/6548/6551.
   `app/shiny_app/www/local_images/*.ome.tiff` is harmonized and carries 6-level pyramids, so it
   is the declared source of truth. Both copies were verified dimension-identical across all 15
   donors, so QuPath-derived per-cell geometry aligns with local_images pixel space.

NAMING
------
Three vocabularies for the same markers, reconciled by `to_panel()`:
    OME-XML channel  : PGP9_5 , KRT8_18
    QuPath TSV column: PGP9.5 , KRT8-18      (`Cell: <name>: Mean`)
    h5ad var_names   : PGP9.5               (the 31-marker analysis panel)
"""

from __future__ import annotations

import re
from pathlib import Path

# --- The 31-marker analysis panel (h5ad var_names order is NOT this; this is sorted) ---------
PANEL_31 = [
    "B3TUBB", "BActin", "BCatenin", "CD163", "CD20", "CD31", "CD34", "CD35", "CD3e", "CD4",
    "CD44", "CD45", "CD68", "CD8a", "CHGA", "CK19", "ColIV", "ECAD", "GAP43", "GCG", "HLADR",
    "INS", "KRT14", "LGALS3", "NaKATPase", "PDPN", "PGP9.5", "SMA", "SST", "VIM", "panCK",
]

# --- Acquisition channel order, in panel vocabulary. -----------------------------------------
# Verified identical across all 15 donors' OME-XML and against the QuPath TSV column order.
# The first 33 are present on every donor; the last two only on the 35-channel donors.
CHANNELS_33 = [
    "DAPI", "CD31", "CD8a", "CD3e", "SMA", "Ki67", "CD4", "CD34", "HLADR", "PDPN", "panCK",
    "ECAD", "CD163", "SST", "ColIV", "VIM", "CD20", "LGALS3", "B3TUBB", "GCG", "KRT14", "GAP43",
    "CD35", "CHGA", "PGP9.5", "INS", "CD44", "CD45", "NaKATPase", "BCatenin", "CD68", "BActin",
    "CK19",
]
CHANNELS_EXTRA = ["epCAM", "KRT8-18"]
CHANNELS_35 = CHANNELS_33 + CHANNELS_EXTRA

# Donors carrying the two extra channels. Derived from the OME-XML, not assumed.
DONORS_35CH = ["6356", "6450", "6479", "6505", "6516", "6538", "6550", "6563"]
DONORS_33CH = ["112", "6521", "6533", "6534", "6548", "6549", "6551"]
DONORS = sorted(DONORS_33CH + DONORS_35CH, key=int)

# --- Channels that are NOT analysis markers ---------------------------------------------------
# DAPI: nuclear stain, never a phenotype marker (Senior's marker_taxonomy EXCLUDES it too).
# Ki67: a PROCESS/proliferation marker with no mutually-exclusive counterpart, so RESTORE pairing
#       does not apply -- and it is the ONE marker QuPath exported as `Median`, not `Mean`
#       (`Cell: Ki67: Median`). Any code reading `Cell: {m}: Mean` must special-case it.
# epCAM / KRT8-18: acquired on only 8 of 15 donors and absent from the 31-marker h5ad panel, so
#       they cannot be used for any cohort-wide call without stratifying the cohort.
NON_PANEL_CHANNELS = ["DAPI", "Ki67"] + CHANNELS_EXTRA
KI67_QUPATH_STAT = "Median"   # every other marker is "Mean"

# --- Name reconciliation ----------------------------------------------------------------------
# OME-XML / QuPath spelling -> panel spelling. Only genuine divergences belong here.
_TO_PANEL = {
    "PGP9_5": "PGP9.5",     # OME-XML
    "KRT8_18": "KRT8-18",   # OME-XML (QuPath TSV already uses KRT8-18)
}


def to_panel(name: str) -> str:
    """Normalize an OME-XML or QuPath marker name to panel vocabulary."""
    return _TO_PANEL.get(name, name)


def ome_channel_names(ome_xml: str | Path) -> list[str]:
    """Channel names from an OME-XML sidecar, in acquisition order, in panel vocabulary.

    Replaces `redsea.qptiff_channels()` for CODEX: Bio-Formats writes the channel name as an
    ATTRIBUTE, which that function's element-form regex does not match (it would return
    ch0..chN for every CODEX image).
    """
    text = Path(ome_xml).read_text() if not str(ome_xml).lstrip().startswith("<") else ome_xml
    names = re.findall(r"<[^>]*Channel[^>]*\bName=\"([^\"]*)\"", text)
    if not names:
        raise ValueError(
            f"no Channel Name attributes found in {ome_xml}; the OME-XML layout changed and "
            "this adapter needs updating -- do NOT fall back to positional indexing"
        )
    return [to_panel(n) for n in names]


def channel_index(marker: str, channels: list[str]) -> int:
    """Index of `marker` within a donor's channel list. Fail loudly rather than guess.

    The whole point of this function is that CODEX donors do not share a channel count, so an
    index valid for one donor is not guaranteed valid for another. Resolve per donor, per read.
    """
    want = to_panel(marker)
    try:
        return channels.index(want)
    except ValueError:
        raise KeyError(
            f"marker {marker!r} (normalized {want!r}) is not among this donor's {len(channels)} "
            f"channels: {channels}"
        ) from None


def expected_channels(donor_id: str) -> list[str]:
    """The channel list a donor should have, from the verified 33/35 split."""
    d = str(donor_id).lstrip("0") or "0"
    if d in DONORS_35CH:
        return list(CHANNELS_35)
    if d in DONORS_33CH:
        return list(CHANNELS_33)
    raise KeyError(f"{donor_id!r} is not a CODEX donor; known: {DONORS}")


def qupath_column(marker: str) -> str:
    """The `Cell: ...` measurement column for a marker in the QuPath per-cell TSV."""
    m = to_panel(marker)
    stat = KI67_QUPATH_STAT if m == "Ki67" else "Mean"
    return f"Cell: {m}: {stat}"


def panel_markers(donor_id: str) -> list[str]:
    """The analysis markers actually available for a donor (panel ∩ that donor's channels)."""
    have = set(expected_channels(donor_id))
    return [m for m in PANEL_31 if m in have]


# Image source of truth. NOT /data/QuPath_Projects/.../Images/ -- see the module docstring.
LOCAL_IMAGES = Path(__file__).resolve().parents[2] / "app" / "shiny_app" / "www" / "local_images"


def image_path(donor_id: str, images_dir: Path | None = None) -> Path:
    """The OME-TIFF for a donor.

    Senior's `redsea.resolve_paths` derives this as `images_dir / img.split(" - ")[0]`, which
    assumes the QuPath image name IS the filename ('...qptiff - resolution #1' -> '...qptiff').
    CODEX image names are bare zero-padded donor ids, so that yields `images_dir/0112` -- a path
    that does not exist. Files here are '<status>_<padded case>.ome.tiff' (Aab_6505.ome.tiff),
    so the mapping has to go through the donor id.
    """
    root = Path(images_dir) if images_dir else LOCAL_IMAGES
    want = normalize_donor_id(donor_id)
    hits = [
        p for p in sorted(root.glob("*.ome.tiff"))
        if normalize_donor_id(p.name.replace(".ome.tiff", "").split("_")[-1]) == want
    ]
    if len(hits) != 1:
        raise FileNotFoundError(
            f"expected exactly 1 image for donor {want!r} in {root}, found {[p.name for p in hits]}"
        )
    return hits[0]


def normalize_donor_id(image_name: str) -> str:
    """QuPath `Image` values are zero-padded ('0112'); the rest of the repo uses '112'.

    Senior derives donor_id as `regexp_extract(Image,'[0-9]+')`, which preserves the padding and
    would produce a `donor_id=0112` partition that never joins the Islet-Explorer data.
    """
    digits = re.search(r"\d+", str(image_name))
    if not digits:
        raise ValueError(f"no donor id in image name {image_name!r}")
    return digits.group(0).lstrip("0") or "0"
