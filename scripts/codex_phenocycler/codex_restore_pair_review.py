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
  * `results_dir` -- a COMPLETE `evaluate_locked_cohort` output (method_spec.json etc.). Use the
    first-pass screen dir; the immune screen used the robust per-pair runner and has no method_spec.

    conda activate scvi-env
    # first-pass bundle (6 pairs, locked 0.995 clip):
    python codex_restore_pair_review.py --results /data/codex_phenocycler/restore_validation/screen_first_pass \
        --output /data/codex_phenocycler/restore_validation/review_first_pass
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

# CODEX first-pass shortlist: the 5 compartment pairs + the endocrine split. These are the reviewable
# set with a complete evaluate_locked_cohort results dir (0.995 clip). The immune pairs are a separate
# bundle (unclipped) built after the first-pass compartments are adjudicated.
CODEX_SHORTLIST_PAIRS = tuple(tuple(p) for p in (cl.FIRST_PASS_PAIRS + [cl.ENDOCRINE_SPLIT_PAIR]))
CODEX_PAIR_RATIONALE = {
    ("ECAD", "CD31"):   "epithelial compartment anchor; must exclude endothelial ECAD bleed",
    ("CD31", "ECAD"):   "endothelial anchor; epithelial adjacency risk",
    ("CD45", "ECAD"):   "pan-leukocyte immune gate; acinar-autofluorescence negative control",
    ("PGP9.5", "panCK"): "neural anchor vs exocrine-restricted panCK (epCAM unavailable cohort-wide)",
    ("VIM", "ECAD"):    "mesenchymal candidate; must exclude genuine epithelial Vimentin",
    ("CHGA", "CD31"):   "endocrine-split anchor inside Epithelial; pan-endocrine chromogranin",
}

CELLS_RAW = crv.CELLS_RAW
CELLS_REDSEA = crv.CELLS_REDSEA
LOCAL_IMAGES = codex_panel.LOCAL_IMAGES
WORKDIR = crv.WORKDIR


def _codex_qptiff_channel_names(series):
    """CODEX channel names from the embedded OME-XML attributes (not Akoya <Biomarker> elements)."""
    tif = series.levels[0].pages[0].parent
    return codex_panel.ome_channel_names(tif.ome_metadata)


def patch_review() -> None:
    crv.patch_upstream()                                   # markers / floors / candidate pairs / pin
    # rpr imported MARKER_INPUT_FLOOR_METHODS into its OWN namespace at import time (used directly in
    # _add_review_suptitle), so patching rv's copy is not enough -- rebind rpr's too, to the CODEX map.
    rpr.MARKER_INPUT_FLOOR_METHODS = rv.MARKER_INPUT_FLOOR_METHODS
    rpr._qptiff_channel_names = _codex_qptiff_channel_names
    rpr.SHORTLIST_PAIRS = CODEX_SHORTLIST_PAIRS
    rpr.PAIR_RATIONALE = dict(CODEX_PAIR_RATIONALE)
    # select_review_queue captured SHORTLIST_PAIRS as a default arg at def time -> re-bind it so the
    # queue is built for the CODEX pairs, not the stale Senior default.
    _orig_queue = rpr.select_review_queue.__wrapped__ if hasattr(rpr.select_review_queue, "__wrapped__") \
        else rpr.select_review_queue
    def _codex_queue(evaluations, *, pairs=CODEX_SHORTLIST_PAIRS):
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


def build_bundle(results_dir: Path, output_dir: Path, *, work: Path = WORKDIR, donors=DONORS,
                 force_build: bool = False):
    patch_review()
    sample, full, raw_clean, redsea_clean, images = build_bundle_inputs(work, donors, force=force_build)
    return rpr.build_review_bundle(sample, full, raw_clean, redsea_clean, images, results_dir, output_dir)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", type=Path, required=True,
                    help="a COMPLETE evaluate_locked_cohort output dir (needs method_spec.json)")
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--work", type=Path, default=WORKDIR)
    ap.add_argument("--donors", nargs="+", default=None)
    ap.add_argument("--force-build", action="store_true")
    args = ap.parse_args()
    donors = args.donors or DONORS
    print(f"Building CODEX review bundle for {len(CODEX_SHORTLIST_PAIRS)} pairs "
          f"(results {args.results.name}, pin {cl.UPSTREAM_PINNED_SHA})")
    outputs = build_bundle(args.results, args.output, work=args.work, donors=donors,
                           force_build=args.force_build)
    for name, path in outputs.items():
        print(f"[review] {name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
