"""Build the CODEX RESTORE pair-review bundle (5-panel PDFs + QuPath CSVs) from a completed screen.

Thin adapter over upstream `restore_pair_review.build_review_bundle`, unmodified. It injects CODEX's
pair shortlist and fixes the two image breakages the plan flagged, then builds the same bundle a
Senior reviewer scores.

CODEX-SPECIFIC PATCHES (all verified against real CODEX OME-TIFFs):
  * `_qptiff_channel_names` -- upstream regexes each TIFF page description for Akoya `<Biomarker>`
    ELEMENTS; CODEX Bio-Formats OME-TIFFs carry channel names as namespaced `<ns0:Channel Name=...>`
    ATTRIBUTES in the embedded OME-XML, so upstream returns ch0..chN. Replaced with a reader that
    pulls `series -> parent.ome_metadata` and calls the proven `codex_panel.ome_channel_names`.
  * IMAGE RESOLUTION -- upstream derives the QPTIFF path as `images_dir / Path(image).name`, and the
    CODEX `image` column is the bare padded donor id ("0112"). Rather than patch the big
    `load_marker_zoom`, `build_images_dir()` makes a directory of no-extension symlinks named exactly
    the image value -> the real `www/local_images/*.ome.tiff` (tifffile opens by content). Verified.
  * PAIR SHORTLIST -- `SHORTLIST_PAIRS`/`PAIR_RATIONALE` and `select_review_queue`'s default are the
    Senior 59-plex set. Both are overridden; the queue function is re-bound because its default arg
    was captured at definition time.

INPUTS the bundle needs (built here):
  * a ::20 sample tree + a FULL (modulus=1) compartment tree -- the full tree gives dense spatial
    context and is loaded with `n == canonical_n`, so the raw/redsea it is checked against must
    contain EXACTLY the same finite-in-both cell set. `build_bundle_inputs()` builds all three
    consistently (the ~47 all-NaN cells excluded from every tree).
  * `results_dir` -- a COMPLETE screen output (method_spec.json etc.). `run_cohort_full_robust` writes
    the full contract for BOTH level-2 clip groups, so each group's screen dir can feed a bundle.

Two bundles, one per clip group (a bundle re-thresholds at its screen's recorded clip, so it can only
shortlist that group's pairs). The 6 level-1/endocrine-split pairs are omitted (already reviewer-
accepted at level 1). The FULL-compartment tree + clean raw/redsea on disk are the stale 13-marker
level-1 ones, so the FIRST level-2 bundle must pass --force-build.

    conda activate scvi-env
    python codex_restore_pair_review.py --group clipped   --force-build   # 11 pairs, 0.995
    python codex_restore_pair_review.py --group unclipped                 # 13 pairs, 1.0 (reuses trees)
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
_PA = Path(os.environ.get("PHENOCYCLER_ROOT",
           str(Path(__file__).resolve().parents[2] / "external" / "Phenocycler_Analysis")))
sys.path.insert(0, str(_PA))
import phenocycler  # noqa: E402,F401  BEFORE numpy
from phenocycler import restore_pair_review as rpr  # noqa: E402
from phenocycler import restore_validation as rv  # noqa: E402

import codex_lineage as cl  # noqa: E402
import codex_panel  # noqa: E402
import codex_restore_validation as crv  # noqa: E402
from codex_panel import DONORS, image_path  # noqa: E402

# Level-2 review shortlists, ONE per clip group (a bundle re-thresholds at its results dir's recorded
# clip, so a bundle can only shortlist pairs from its own group). The 6 level-1/endocrine-split pairs
# (ECAD/CD31/CD45/PGP9.5/VIM<-... + CHGA<-CD31) are OMITTED: they were reviewer-accepted at level 1 and
# reproduced their exact divisors at v8 (8582/2580/2623/3411/3817/650), so they need no re-review.
#
# CLIPPED group (0.995): the new level-2 pairs whose target is not clip-sensitive.
CODEX_SHORTLIST_CLIPPED = (
    ("CD20", "ECAD"), ("CD68", "ECAD"), ("CD163", "ECAD"),
    ("CK19", "CD31"), ("panCK", "CD31"),
    ("GCG", "INS"), ("GCG", "panCK"),
    ("CD34", "ECAD"), ("GAP43", "panCK"), ("B3TUBB", "panCK"), ("SMA", "ECAD"),
)
# UNCLIPPED group (1.0): the clip-sensitive targets (bright tail above p99.5). All new.
CODEX_SHORTLIST_UNCLIPPED = (
    ("CD3e", "ECAD"), ("CD4", "ECAD"), ("CD8a", "ECAD"), ("HLADR", "ECAD"), ("CD35", "ECAD"),
    ("KRT14", "CD31"), ("SST", "INS"), ("SST", "panCK"), ("INS", "GCG"), ("INS", "panCK"),
    ("PDPN", "ECAD"), ("CD44", "ECAD"), ("LGALS3", "CD31"),
)
CODEX_PAIR_RATIONALE = {
    # --- clipped group ---
    ("CD20", "ECAD"):   "B target, re-included after the both-role unwind; watch ECAD+ epithelial stray CD20",
    ("CD68", "ECAD"):   "macrophage; ECAD counterpart (whole-scene threshold, interpreted within CD45+)",
    ("CD163", "ECAD"):  "M2 macrophage; nested under CD68 at interpretation",
    ("CK19", "CD31"):   "ductal epithelium; CD31 negative is undersized (endothelium sparse) but exclusive",
    ("panCK", "CD31"):  "epithelial confirmatory; CD31 undersized like ECAD<-CD31 (accepted precedent)",
    ("GCG", "INS"):     "alpha; INS reference undersized in T1D beta-loss -- compare the panCK comparator",
    ("GCG", "panCK"):   "alpha exocrine comparator: is panCK a cleaner reference than the undersized INS?",
    ("CD34", "ECAD"):   "capillary EC / resting fibroblast (dual-use, parent-resolved at interpretation)",
    ("GAP43", "panCK"): "sprouting axons; exocrine-restricted panCK negative (epCAM dropped, 8/15-donor)",
    ("B3TUBB", "panCK"): "pan-neural consensus with PGP9.5; panCK negative",
    ("SMA", "ECAD"):    "smooth muscle / pericyte; positive mesenchymal call (VIM is the residual gate)",
    # --- unclipped group ---
    ("CD3e", "ECAD"):   "T anchor; ECAD counterpart; CD4/CD8 read WITHIN CD3e+ at interpretation",
    ("CD4", "ECAD"):    "T-helper; interpret only within CD3e+ (monocytes are CD4-dim)",
    ("CD8a", "ECAD"):   "T-cytotoxic; unclipped (bright tail above the 0.995 clip on low-rate donors)",
    ("HLADR", "ECAD"):  "WEAKEST pair -- ductal upregulates HLA-II in T1D, so the ECAD negative is contaminated; flag conservative",
    ("CD35", "ECAD"):   "FDC/TLS -- near-zero without lymphoid aggregates (prevalence-collapse marker)",
    ("KRT14", "CD31"):  "basal/metaplastic -- near-zero in normal tissue (prevalence-collapse marker)",
    ("SST", "INS"):     "delta; INS reference undersized in T1D -- compare the panCK comparator",
    ("SST", "panCK"):   "delta exocrine comparator vs the undersized hormone reference",
    ("INS", "GCG"):     "beta; GCG reference undersized in ND (alpha<beta) -- compare the panCK comparator",
    ("INS", "panCK"):   "beta exocrine comparator vs the undersized hormone reference",
    ("PDPN", "ECAD"):   "lymphatic EC / inflammatory-CAF (dual-use; CD31 gate disambiguates at interpretation)",
    ("CD44", "ECAD"):   "activation/stemness STATE annotation, not a lineage class call",
    ("LGALS3", "CD31"): "macrophage+ductal STATE annotation; CD31 undersized -- never a class call",
}

CELLS_RAW = crv.CELLS_RAW
CELLS_REDSEA = crv.CELLS_REDSEA
LOCAL_IMAGES = codex_panel.LOCAL_IMAGES
WORKDIR = crv.WORKDIR


def _codex_qptiff_channel_names(series):
    """CODEX channel names from the embedded OME-XML attributes (not Akoya <Biomarker> elements)."""
    tif = series.levels[0].pages[0].parent
    return codex_panel.ome_channel_names(tif.ome_metadata)


def patch_review(shortlist=CODEX_SHORTLIST_CLIPPED, rationale=None) -> None:
    crv.patch_upstream()                                   # markers / floors / candidate pairs / pin
    shortlist = tuple(tuple(p) for p in shortlist)
    rationale = dict(rationale or {p: CODEX_PAIR_RATIONALE.get(p, "") for p in shortlist})
    # rpr imported MARKER_INPUT_FLOOR_METHODS into its OWN namespace at import time (used directly in
    # _add_review_suptitle), so patching rv's copy is not enough -- rebind rpr's too, to the CODEX map.
    rpr.MARKER_INPUT_FLOOR_METHODS = rv.MARKER_INPUT_FLOOR_METHODS
    rpr._qptiff_channel_names = _codex_qptiff_channel_names
    rpr.SHORTLIST_PAIRS = shortlist
    rpr.PAIR_RATIONALE = rationale
    # select_review_queue captured SHORTLIST_PAIRS as a default arg at def time -> re-bind it so the
    # queue is built for the CODEX pairs, not the stale Senior default.
    _orig_queue = rpr.select_review_queue.__wrapped__ if hasattr(rpr.select_review_queue, "__wrapped__") \
        else rpr.select_review_queue
    def _codex_queue(evaluations, *, pairs=shortlist):
        return _orig_queue(evaluations, pairs=pairs)
    rpr.select_review_queue = _codex_queue


def _duckdb():
    import duckdb
    con = duckdb.connect()
    con.execute("PRAGMA threads=8; SET enable_progress_bar=false")
    return con


def build_images_dir(out_dir: Path, donors=DONORS, *, force: bool = False) -> Path:
    """A dir of no-extension symlinks named exactly the `image` column value -> real ome.tiff."""
    if out_dir.exists() and not force:
        return out_dir
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    import pandas as pd, glob
    for donor in donors:
        f = glob.glob(f"{WORKDIR}/sample/donor_id={donor}/*.parquet")
        img = (pd.read_parquet(f[0], columns=["image"])["image"].iloc[0] if f
               else str(donor).zfill(4))
        (out_dir / str(img)).symlink_to(image_path(donor).resolve())
    return out_dir


def _build_clean_canon(full_sample: Path, out_raw: Path, out_redsea: Path, donors) -> None:
    """Build raw_view + redsea restricted to the FULL sample's cell set, so the bundle's full-cell
    audit (n == canonical_n) holds and raw/redsea key sets match."""
    ids = f"(SELECT DISTINCT CAST(object_id AS VARCHAR) AS oid FROM read_parquet('{full_sample}/donor_id=*/*.parquet', hive_partitioning=true))"
    markers = ", ".join(f'"{m}"' for m in crv.SCREENED_MARKERS)
    din = ", ".join("'" + str(d) + "'" for d in donors)
    con = _duckdb()
    for src, out, extra in ((CELLS_RAW, out_raw, "centroid_x_um AS X_centroid, centroid_y_um AS Y_centroid, "),
                            (CELLS_REDSEA, out_redsea, "")):
        if out.exists():
            shutil.rmtree(out)
        tmp = out.with_name(out.name + ".building")
        if tmp.exists():
            shutil.rmtree(tmp)
        con.execute(f"""
        COPY (
          SELECT CAST(c.donor_id AS VARCHAR) AS donor_id, CAST(c.object_id AS VARCHAR) AS object_id,
                 {extra}{markers}
          FROM read_parquet('{src}/donor_id=*/*.parquet', hive_partitioning=true) c
          INNER JOIN {ids} f ON CAST(c.object_id AS VARCHAR) = f.oid
          WHERE CAST(c.donor_id AS VARCHAR) IN ({din})
        ) TO '{tmp}' (FORMAT PARQUET, PARTITION_BY (donor_id), COMPRESSION ZSTD)
        """)
        tmp.rename(out)
    con.close()


def build_bundle_inputs(work: Path, donors=DONORS, *, force: bool = False):
    """Build ::20 sample, full compartment tree, and matching clean raw/redsea. Returns the paths."""
    work.mkdir(parents=True, exist_ok=True)
    sample = crv.build_sample(work / "sample", donors, modulus=20, force=force)
    full = crv.build_sample(work / "full_sample", donors, modulus=1, force=force)
    raw_clean, redsea_clean = work / "raw_clean", work / "redsea_clean"
    if force or not raw_clean.exists() or not redsea_clean.exists():
        _build_clean_canon(full, raw_clean, redsea_clean, donors)
    images = build_images_dir(work / "images", donors, force=force)
    return sample, full, raw_clean, redsea_clean, images


def build_bundle(results_dir: Path, output_dir: Path, *, shortlist=CODEX_SHORTLIST_CLIPPED,
                 work: Path = WORKDIR, donors=DONORS, force_build: bool = False):
    patch_review(shortlist)
    sample, full, raw_clean, redsea_clean, images = build_bundle_inputs(work, donors, force=force_build)
    return rpr.build_review_bundle(sample, full, raw_clean, redsea_clean, images, results_dir, output_dir)


_GROUPS = {
    "clipped":   (CODEX_SHORTLIST_CLIPPED,   "screen_l2_clipped",   "review_l2_clipped"),
    "unclipped": (CODEX_SHORTLIST_UNCLIPPED, "screen_l2_unclipped", "review_l2_unclipped"),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--group", choices=sorted(_GROUPS), required=True,
                    help="which clip group's bundle to build (its shortlist + its screen dir)")
    ap.add_argument("--results", type=Path, default=None,
                    help="override the screen dir (default WORKDIR/screen_l2_<group>)")
    ap.add_argument("--output", type=Path, default=None,
                    help="override the output dir (default WORKDIR/review_l2_<group>)")
    ap.add_argument("--work", type=Path, default=WORKDIR)
    ap.add_argument("--donors", nargs="+", default=None)
    ap.add_argument("--force-build", action="store_true",
                    help="rebuild the full-compartment tree + clean raw/redsea (REQUIRED the first "
                         "level-2 build -- the on-disk trees are the stale 13-marker level-1 ones)")
    args = ap.parse_args()
    donors = args.donors or DONORS
    shortlist, screen_name, review_name = _GROUPS[args.group]
    results = args.results or (WORKDIR / screen_name)
    output = args.output or (WORKDIR / review_name)
    print(f"Building CODEX {args.group} review bundle for {len(shortlist)} pairs "
          f"(results {results.name}, pin {cl.UPSTREAM_PINNED_SHA}) -> {output}")
    outputs = build_bundle(results, output, shortlist=shortlist, work=args.work, donors=donors,
                           force_build=args.force_build)
    for name, path in outputs.items():
        print(f"[review] {name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
