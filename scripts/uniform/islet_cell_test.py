#!/usr/bin/env python3
"""islet_cell_test.py -- CODEX port of the library islet_classifier_fit.py imports.

Ported from Islet-Explorer-Senior's scripts/v3/islet_cell_test.py (see ../PORTED_FROM.json). Only the
library survives; the Senior command line was already withdrawn as a workflow step there. CODEX changes:

- CHGA is the level marker in CD99's place (owner decision 2026-10-15); this panel has no CD99.
- The identity markers come from LINEAGE_MARKERS below, the Senior lineage table mapped onto the 31
  analysis markers of this CODEX panel. Markers Senior lists that this panel lacks are simply absent.
  PGP9.5 is left out (SHARED_WITH_ENDOCRINE): on this panel it is brightest in beta cells.
- Values are QuPath's raw 'Cell: <marker>: Mean' (this pipeline has no REDSEA correction). The highest
  identity marker is read on common_scale(): each marker over its own 99.9th percentile in the image, a
  hormone's never below the cohort's median. On raw CODEX counts E-cadherin outranks insulin in beta cells
  (canary 6505, 2026-10-15).

THE REFERENCE LEVEL is Senior's: the REFERENCE_QUANTILE-th percentile (99) of the level marker over cells
that cannot belong to an islet -- with a nucleus, in Tissue, away from hormone cells (and, when structure
objects are given, outside them). Fitted on each image alone, never pooled.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import shapely
from shapely import STRtree
from shapely.geometry import shape

HORMONES = ("INS", "GCG", "SST")
LEVEL_MARKER = "CHGA"
NOT_OBJECTS = ("Tissue", "pLN")
TISSUE_PARENT = "Annotation (Tissue)"
STRUCTURES = ("Nerve", "Capillary", "Lymphatic")
REFERENCE_QUANTILE = 99.0   # Senior subpanel_io.REFERENCE_QUANTILE

# Senior's lineages.LINEAGE_MARKERS with this panel's marker names. Endocrine keeps the hormones only;
# CHGA is the level marker, not an identity marker, exactly as CD99 is in Senior.
LINEAGE_MARKERS: dict[str, tuple[str, ...]] = {
    "Endocrine":   ("INS", "GCG", "SST"),
    "Exocrine":    ("panCK", "CK19", "ECAD", "KRT14"),
    "Vascular":    ("CD31", "PDPN", "CD34"),
    "Stromal":     ("SMA", "VIM", "ColIV"),
    "T_NK":        ("CD3e", "CD8a", "CD4"),
    "B_Plasma":    ("CD20",),
    "Myeloid":     ("CD68", "CD163", "LGALS3"),
    "Granulocyte": (),
    "Neural":      ("B3TUBB", "GAP43"),
}
# Markers that cannot name a lineage on this panel because islet endocrine cells carry them too (Senior's
# lineages.SUPPORTING_ONLY reasoning: "too broad or too shared to name a cell by itself"). PGP9.5 (UCHL1),
# median raw cell means on the canaries of 2026-10-15, cells matched to the earlier typing:
#   6505: beta 4,034 / delta 1,646 / alpha 628 / nerve 856;  6232: beta 4,095 / delta 2,469 / alpha 660 / nerve 4,966.
# As a Neural marker it was the highest marker of 67 % (6505) and 77 % (6232) of beta cells, and a Neural
# highest marker bars a cell from the level-marker route, so the fit learned to leave insulin out.
# GAP43 stays Neural: beta 65 / 45 against nerve 881 / 384.
SHARED_WITH_ENDOCRINE: tuple[str, ...] = ("PGP9.5",)
IDENTITY = sorted({m for markers in LINEAGE_MARKERS.values() for m in markers})
# The quantile of each marker, over the image's cells with a nucleus in Tissue, that common_scale() divides by.
SCALE_QUANTILE = 99.9
# Hormones whose top is never below the cohort's median top (cohort_scale.py). A p99.9 lands on a hormone's cells
# only where they are more than 1 per 1,000 of the image's cells. Earlier typing, 20 images (2026-10-15): beta cells
# fall below that in four T1D images (0.02-0.24 per 1,000) and delta cells in 14 (0.04-0.91), so there the image's
# own top is background; alpha cells never do (5-35 per 1,000, T1D included). A glucagon floor only penalized dimmer
# acquisitions (additions' glucagon tops 9,249-17,425 against a median of 18,288): it halved the share of earlier
# alpha cells with glucagon highest (ND_6232 53 -> 25 %, Aab_6090 57 -> 23 %) and moved them to `other`.
FLOORED_HORMONES = ("INS", "SST")


def sanitize(name: str) -> str:
    """Channel name -> marker name, as Senior's phenocycler.redsea.SANITIZE does."""
    return re.sub(r"[ /\-]", "_", name)


def subdivide(geoms, max_vertices: int = 256):
    """Cut large polygons into pieces so point-in-polygon queries stay fast (Senior uses
    build_crossed_fraction.subdivide for the same purpose)."""
    for g in geoms:
        if g.is_empty:
            continue
        for part in shapely.get_parts(g):
            yield from _split(part, max_vertices)


def _split(g, max_vertices: int):
    if shapely.get_num_coordinates(g) <= max_vertices:
        yield g
        return
    minx, miny, maxx, maxy = g.bounds
    if maxx - minx >= maxy - miny:
        mid = (minx + maxx) / 2
        halves = (shapely.box(minx, miny, mid, maxy), shapely.box(mid, miny, maxx, maxy))
    else:
        mid = (miny + maxy) / 2
        halves = (shapely.box(minx, miny, maxx, mid), shapely.box(minx, mid, maxx, maxy))
    for h in halves:
        part = g.intersection(h)
        for piece in shapely.get_parts(part):
            if piece.geom_type in ("Polygon", "MultiPolygon") and not piece.is_empty:
                yield from _split(piece, max_vertices)


def read_objects(path: Path) -> list:
    raw = json.loads(Path(path).read_text())
    feats = raw["features"] if isinstance(raw, dict) else raw
    out = []
    for f in feats:
        if not isinstance(f, dict) or not f.get("geometry"):
            continue
        if ((f.get("properties") or {}).get("classification") or {}).get("name") in NOT_OBJECTS:
            continue
        g = shape(f["geometry"])
        if not g.is_valid:
            g = shapely.make_valid(g)
        if not g.is_empty:
            out.append(g)
    return out


def centroid_in_any(x_px: np.ndarray, y_px: np.ndarray, geoms: list) -> np.ndarray:
    """True where the point lies within one of `geoms` (large objects are cut first, so the test stays fast)."""
    inside = np.zeros(len(x_px), dtype=bool)
    if not geoms:
        return inside
    pieces = list(subdivide(geoms))
    hit, _ = STRtree(pieces).query(shapely.points(x_px, y_px), predicate="within")
    inside[np.unique(hit)] = True
    return inside


def image_tops(values: pd.DataFrame, eligible: np.ndarray, quantile: float = SCALE_QUANTILE) -> dict:
    """Each column's `quantile`-th percentile over the `eligible` rows (NaN for a column with no measurement)."""
    tops = {}
    for c in values.columns:
        col = values[c].to_numpy(dtype=float)[eligible]
        col = col[np.isfinite(col)]
        tops[c] = float(np.percentile(col, quantile)) if col.size else float("nan")
    return tops


def common_scale(values: pd.DataFrame, eligible: np.ndarray, quantile: float = SCALE_QUANTILE,
                 floor: dict | None = None) -> tuple[pd.DataFrame, dict]:
    """Each column divided by its own `quantile`-th percentile over the `eligible` rows, so that the highest
    marker compares markers on one scale. CODEX substitution: Senior's fit compares REDSEA-corrected counts
    directly, and Senior's typing lane scales every marker to its p99.9 before its argmax. Here raw counts
    differ by two orders of magnitude between channels (E-cadherin p99.9 ~37,000, insulin ~800 on 6505).

    The top is the image's own, so channel differences between acquisitions cancel. `floor` (marker -> value)
    raises a top that falls below it: the fit passes the cohort's median top for FLOORED_HORMONES, because an image
    can lack a hormone's cells almost entirely (insulin in T1D) and then its own p99.9 is background, which the
    division would turn into "hormone highest" (T1D_6534, 2026-10-15: 8,255 insulin-highest cells more than 30 um
    from any earlier endocrine cell on its own top, 23 with the floor). A column whose top is not positive keeps its
    raw values. Only the highest-marker test reads the scaled values; levels stay on raw values."""
    top = image_tops(values, eligible, quantile)
    for c, f in (floor or {}).items():
        if c in top and np.isfinite(f) and not top[c] >= f:
            top[c] = float(f)
    out = values.copy()
    for c, t in top.items():
        if np.isfinite(t) and t > 0:
            out[c] = values[c] / t
    return out, top


def top_marker(values: pd.DataFrame) -> np.ndarray:
    """The highest of the columns per row; "" where none was measured (never a tie-break winner)."""
    m = values.to_numpy(dtype=float)
    none = np.isnan(m).all(axis=1)
    top = np.asarray(values.columns, dtype=object)[np.argmax(np.where(np.isnan(m), -np.inf, m), axis=1)]
    return np.where(none, "", top)


def reference_level(values: np.ndarray, quantile: float) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise SystemExit(f"{LEVEL_MARKER}: the reference population carries no measurement -- refusing to fit a level")
    level = float(np.percentile(finite, quantile))
    if not np.isfinite(level) or level <= 0:
        # a level at the floor admits every measured cell while looking fitted (Senior subpanel_io.endocrine_level)
        raise SystemExit(f"{LEVEL_MARKER}: reference level {level!r} is at the floor -- refusing to classify")
    return level
