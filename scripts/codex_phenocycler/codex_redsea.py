#!/usr/bin/env python
"""Run REDSEA spillover correction on the CODEX cohort.

Two functions in `phenocycler.redsea` are Senior-specific and are PATCHED here rather than
forked, so this repo does not carry a divergent copy of an upstream package under active
development:

  * `resolve_paths` derives the image file as `images_dir / img.split(" - ")[0]`, which assumes
    the QuPath image name IS the filename ('...qptiff - resolution #1' -> '...qptiff'). CODEX
    image names are bare zero-padded donor ids, so it yields `images_dir/0112` — nonexistent.
  * `qptiff_channels` reads channel names by regexing each page description for Akoya
    `<Biomarker>` / `<Name>` ELEMENTS. CODEX OME-TIFFs store the name as an XML ATTRIBUTE, so
    every channel would come back `ch0..chN`. That does not raise — it silently produces a
    corrected table whose columns cannot join the marker panel.

Both patches fail loudly. `codex_qptiff_channels` cross-checks the names it recovers against
`codex_panel.expected_channels()` for that donor and raises on any divergence, so a changed
export or a mis-mapped donor stops the run instead of quietly producing wrong columns.

    conda activate scvi-env
    python codex_redsea.py --check                 # patch + resolve every donor, read nothing
    python codex_redsea.py --donor 112             # one donor
    python codex_redsea.py --all --jobs 2

Cost, from Senior's instrumentation (2.53M cells, 59 channels, CPU): ~31 min/donor, ~66% of it
in the per-channel bincount. CODEX is 33-35 channels, so expect materially less. MEMORY is the
real constraint, not CPU: an int32 (H,W) mask is ~0.5 GB at 12672x9504, and gap_bridge>0 adds a
`distance_transform_edt(..., return_indices=True)` allocation several times that. Keep --jobs low.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# Locate the Phenocycler_Analysis package. Override with PHENOCYCLER_ROOT. The plan is to consume
# it as a git submodule + `pip install -e .`; until then this keeps the scripts runnable standalone.
_PA = Path(os.environ.get(
    "PHENOCYCLER_ROOT",
    "/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/Phenocycler_Analysis",
))
if not (_PA / "phenocycler").is_dir():
    raise SystemExit(
        f"[err] Phenocycler_Analysis not found at {_PA}. Set PHENOCYCLER_ROOT or pip install -e it."
    )
sys.path.insert(0, str(_PA))
sys.path.insert(0, str(Path(__file__).resolve().parent))

# Import phenocycler BEFORE numpy: its __init__ pins the BLAS/OpenMP thread counts to avoid a
# documented futex_wait deadlock (0 markers in 15 h). Importing numpy first defeats that.
import phenocycler  # noqa: E402,F401  (import order is load-bearing)

import argparse  # noqa: E402

import tifffile

from codex_panel import DONORS, expected_channels, image_path, normalize_donor_id

from phenocycler import redsea as _redsea
from phenocycler.config import PipelineConfig

CODEX_DATA_DIR = Path("/data/codex_phenocycler")


def codex_resolve_paths(cfg, donor: str):
    """(image_name, geojson, image_path) for a CODEX donor.

    Same contract as `redsea.resolve_paths`, but the image is looked up by donor id and the
    GeoJSON tag uses QuPath's ZERO-PADDED image name, which is what the export groovy wrote.
    """
    img = _redsea.donor_image(cfg, donor)              # QuPath `image` column, e.g. '0112'
    tag = "".join(c if c.isalnum() or c in "._-" else "_" for c in img)
    geojson = cfg.geojson_dir / f"cells__{tag}.geojson"
    if not geojson.exists():
        raise SystemExit(
            f"[err] missing geojson: {geojson}\n"
            f"      run: bash export_cells_geojson.sh {img}"
        )
    tif = image_path(donor)
    if not tif.exists():
        raise SystemExit(f"[err] missing image: {tif}")
    return img, geojson, tif


def codex_qptiff_channels(tif: Path, level: int):
    """(tiff, level_pages, channel_names, (H,W), pixel_size_um) for a CODEX OME-TIFF.

    Channel names come from the OME-XML `Name` ATTRIBUTE (via the sidecar when present, else the
    embedded page-0 description), normalized to panel vocabulary, then CROSS-CHECKED against the
    donor's expected channel list. Divergence raises: a silently mislabelled channel is the worst
    failure mode available here, because everything downstream still runs.
    """
    tf = tifffile.TiffFile(str(tif))
    series = tf.series[0]
    lv = series.levels[level]
    pages = lv.pages

    from codex_panel import ome_channel_names
    sidecar = tif.with_suffix("").with_suffix(".ome.xml")   # ND_0112.ome.tiff -> ND_0112.ome.xml
    source = sidecar if sidecar.exists() else (series.levels[0].pages[0].description or "")
    names = ome_channel_names(source)

    donor = normalize_donor_id(tif.name.replace(".ome.tiff", "").split("_")[-1])
    want = expected_channels(donor)
    if names != want:
        raise SystemExit(
            f"[err] {tif.name}: channel names do not match codex_panel.expected_channels({donor!r}).\n"
            f"      got      {names}\n      expected {want}\n"
            "      Refusing to run: mislabelled channels would produce plausible wrong output."
        )
    if len(pages) != len(names):
        raise SystemExit(
            f"[err] {tif.name}: {len(pages)} pages at level {level} but {len(names)} channel names"
        )

    H, W = lv.shape[-2], lv.shape[-1]
    px_um = None
    try:
        p0 = series.levels[0].pages[0]
        xr = p0.tags["XResolution"].value
        unit = int(p0.tags["ResolutionUnit"].value)          # 3 == cm, 2 == inch
        ppu = xr[0] / xr[1]
        per_cm = ppu if unit == 3 else ppu / 2.54
        px_um = 1e4 / per_cm * (2 ** level)
    except Exception:
        pass
    return tf, pages, names, (H, W), px_um


def install_patches() -> None:
    _redsea.resolve_paths = codex_resolve_paths
    _redsea.qptiff_channels = codex_qptiff_channels


def _run_one(donor: str, *, cfg, params):
    """Worker entry point for the per-donor pool.

    Re-installs the patches inside the worker. On Linux/CPython 3.13 the pool forks, so the
    parent's patched module state is inherited and this is redundant -- but the default start
    method is scheduled to change, and a spawned worker would re-import `phenocycler.redsea`
    clean and silently run Senior's `resolve_paths`/`qptiff_channels` instead. That failure is
    not loud: `qptiff_channels` would return ch0..chN and the run would complete with unjoinable
    columns. Cheap insurance against a very expensive silent wrong answer.
    """
    install_patches()
    return _redsea.process_donor(donor, cfg, params)


def make_config(data_dir: Path) -> PipelineConfig:
    # images_dir is unused once resolve_paths is patched, but set it truthfully anyway.
    from codex_panel import LOCAL_IMAGES
    return PipelineConfig(data_dir=data_dir, images_dir=LOCAL_IMAGES)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=CODEX_DATA_DIR)
    ap.add_argument("--donor", action="append", help="repeatable")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--jobs", type=int, default=1, help="per-donor pool; memory-bound, keep low")
    ap.add_argument("--check", action="store_true",
                    help="resolve paths and channel names for every donor; read no pixels")
    args = ap.parse_args()

    install_patches()
    cfg = make_config(args.data_dir)
    donors = DONORS if (args.all or args.check) else (args.donor or [])
    if not donors:
        ap.error("pass --donor, --all or --check")

    if args.check:
        bad = 0
        for d in donors:
            try:
                img, gj, tif = codex_resolve_paths(cfg, d)
                _, pages, names, (H, W), px = codex_qptiff_channels(tif, 0)
                print(f"  {d:>6}  image={img:<5} {tif.name:<22} {len(names)}ch {W}x{H} "
                      f"px_um={px:.4f}  geojson={gj.stat().st_size/1e6:.0f}MB")
            except SystemExit as e:
                print(f"  {d:>6}  FAIL {e}")
                bad += 1
        print(f"\n{len(donors) - bad}/{len(donors)} donors resolve cleanly")
        return 1 if bad else 0

    from phenocycler.redsea import RedseaParams
    params = RedseaParams()

    # Skip donors already done, so an interrupted cohort run resumes instead of restarting.
    out_dir = cfg.cells_redsea_dir
    todo = [d for d in donors if not list((out_dir / f"donor_id={d}").glob("*.parquet"))]
    if len(todo) < len(donors):
        print(f"Skipping {len(donors) - len(todo)} already-corrected donors: "
              f"{sorted(set(donors) - set(todo), key=int)}")
    if not todo:
        print("Nothing to do.")
        return 0

    if args.jobs > 1:
        from functools import partial

        from phenocycler.parallel import map_donors, resolve_jobs
        n = resolve_jobs(args.jobs, len(todo))
        print(f"Running {len(todo)} donors across {n} processes")
        # functools.partial over a TOP-LEVEL function, never a lambda: the callable is pickled
        # across the process boundary, and a closure fails with
        # "Can't get local object 'main.<locals>.<lambda>'".
        results = map_donors(partial(_run_one, cfg=cfg, params=params), todo,
                             n_jobs=n, ordered=True, on_error="log")
        failed = [d for d, r in zip(todo, results) if isinstance(r, Exception) or r is None]
    else:
        failed = []
        for i, d in enumerate(todo, 1):
            print(f"=== {d}  ({i}/{len(todo)}) ===", flush=True)
            try:
                _run_one(d, cfg, params)
            except Exception as e:                      # keep going; report at the end
                print(f"  FAILED {d}: {e}", flush=True)
                failed.append(d)

    print(f"\n{len(todo) - len(failed)}/{len(todo)} donors corrected -> {out_dir}")
    if failed:
        print(f"FAILED: {sorted(failed, key=int)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
