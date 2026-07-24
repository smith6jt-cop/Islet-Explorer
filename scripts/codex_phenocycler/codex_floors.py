"""Measure the donor-local input-QC arm floors RESTORE uses, and pick each marker's floor statistic.

WHY THIS EXISTS
---------------
`codex_lineage.MARKER_INPUT_FLOOR_METHODS` was a documentation registry with zero callers -- there
was no floor-estimation code in the CODEX port at all. This module is the execution path: it calls
upstream's OWN `restore_validation.estimate_marker_input_floor` (so the numbers match what RESTORE
computes internally) over the CODEX cohort, and:

  1. floors on the **RAW** parquet -- corrected 2026-07-23. Upstream's `pair_input_qc` floors on
     `raw__{marker}`, and `load_validation_sample` np.allclose-asserts `raw__{m}` == the UNCORRECTED
     QuPath `Cell: {m}: Mean`. So `OBSERVED_FLOORS_RAW` is the operative table; REDSEA values enter
     only as the divisor/separator axis. (The old "measure floors post-REDSEA" instruction was wrong.)

  2. reports `p99_5_over_median` per donor-marker -- the predictor of the p99.5-clip failure that
     produced the CD8a/HLADR/CHGA anomalies. `raw_otsu_floor` clips at p99.5 before running Otsu; a
     marker whose true positives live above the clip loses its positive mode and 2-component Otsu
     bisects a unimodal background. Every marker whose per-donor MIN p99.5/median >= DYNRANGE_MIN is
     safe on all 15 donors; below it, at least one donor's Otsu lands in the background.

  3. picks each marker's method by that guard, NOT by inheriting Senior's map (which `floor_method()`
     fails closed to prevent): max(Otsu, Triangle) for the low-dynamic-range markers (the same
     two-statistic floor upstream keeps for CD3e/CD20), plain Otsu for the high-dynamic-range ones.

The floor is the arm-QC threshold and the Step-1 anchor input; it is NOT the positivity call. For
CD8a/HLADR specifically, max(Otsu,Triangle) fixes the ARM floor but is not a validated positivity
threshold (it makes CD8a+ ~= CD3e+, and HLADR miss 93-99% of CD68+); those markers are routed to
RESTORE + image review for the actual call. See codex_lineage.PREVALENCE_COLLAPSE_MARKERS and
codex_lineage.CLIP_UNCLIPPED_TARGETS.

    conda activate scvi-env
    python codex_floors.py                 # measure + print the recommended method map
    python codex_floors.py --check         # regression-assert against OBSERVED_FLOORS_RAW
    python codex_floors.py --parquet /data/codex_phenocycler/cells_redsea   # diagnostic only
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from pathlib import Path

# codex modules first; phenocycler BEFORE numpy (BLAS pinning).
sys.path.insert(0, str(Path(__file__).resolve().parent))
_PA = Path(os.environ.get("PHENOCYCLER_ROOT",
           str(Path(__file__).resolve().parents[2] / "external" / "Phenocycler_Analysis")))
sys.path.insert(0, str(_PA))
import phenocycler  # noqa: E402,F401  BEFORE numpy
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402
from phenocycler import restore_validation as rv  # noqa: E402

import codex_lineage as cl  # noqa: E402
from codex_panel import DONORS, to_panel  # noqa: E402

RAW_PARQUET = Path("/data/codex_phenocycler/cells")

# Per-donor minimum of p99.5/median at or above which a marker's 2-component Otsu is safe on every
# donor. Below it, the p99.5 clip CAN delete the positive mode -- but ONLY for a sparse marker. An
# ABUNDANT marker (ECAD 31%, panCK 16%) also has a low p99.5/median, because its median already sits
# in the positive population, not because Otsu fails -- so low dynamic range must be paired with low
# prevalence before it means "broken". Measured cut: the five anomalous markers (CHGA/CD3e/CD4/CD8a/
# HLADR) are the only ones that are BOTH low-dynrange (<=4.6) AND sparse (<2% at the max floor).
DYNRANGE_MIN = 7.0
SPARSE_MAX_PREV = 10.0   # a marker above this prevalence is abundant, and low dynrange is expected

# Markers whose floors the RESTORE screen needs: every marker in the level-2 screen (30 pairs).
# INS/GCG/SST ARE here now -- at level 2 the endocrine subtypes are RESTORE pairs (INS<-GCG etc.),
# so the hormones are screened targets and need floors.
SCREENED_MARKERS = sorted({m for pair in cl.ALL_SCREEN_PAIRS for m in pair})


def _marker_present(donor_dir: str, marker: str) -> bool:
    f = glob.glob(f"{donor_dir}/*.parquet")[0]
    return marker in pq.ParquetFile(f).schema.names


def measure_floors(parquet_dir: Path, markers=SCREENED_MARKERS, donors=DONORS) -> pd.DataFrame:
    """Per (donor, marker): Otsu, Triangle, p99.5/median, and the prevalence at each floor.

    Uses upstream's raw_otsu_floor / raw_triangle_floor verbatim, so the numbers equal what
    RESTORE's pair_input_qc computes.
    """
    rows = []
    donor_dirs = {d.split("donor_id=")[1]: d
                  for d in sorted(glob.glob(f"{parquet_dir}/donor_id=*"))}
    for marker in markers:
        pm = to_panel(marker)
        for donor in donors:
            ddir = donor_dirs.get(donor)
            if ddir is None or not _marker_present(ddir, pm):
                continue
            v = pd.read_parquet(ddir, columns=[pm])[pm].to_numpy(float)
            v = v[np.isfinite(v)]
            if v.size == 0:
                continue
            otsu = float(rv.raw_otsu_floor(v))
            tri = float(rv.raw_triangle_floor(v))
            med = float(np.median(v))
            p995 = float(np.quantile(v, 0.995))
            rows.append(dict(
                donor=donor, marker=pm, median=med, otsu=otsu, triangle=tri,
                p99_5_over_median=(p995 / med if med > 0 else np.inf),
                prev_otsu=float((v >= otsu).mean() * 100),
                prev_max=float((v >= max(otsu, tri)).mean() * 100),
            ))
    return pd.DataFrame(rows)


def clip_sensitive(df_marker: pd.DataFrame) -> bool:
    """True if the p99.5 clip can bisect background for this marker: BOTH low-dynamic-range AND sparse.

    Low dynamic range alone is not enough -- an abundant marker (ECAD, panCK) has a low p99.5/median
    because its median already sits in the positive population, and plain clipped Otsu is correct
    there. Requiring sparseness too leaves exactly the markers whose positive tail sits above p99.5 on
    the low-rate donors (CHGA/CD3e/CD4/CD8a/HLADR). The RESOLUTION for these is NOT max(Otsu,Triangle)
    -- it is to screen them UNCLIPPED (codex_lineage.IMMUNE_SCREEN_UPPER_QUANTILE); CHGA is the one
    kept on max()@0.995 because it rides the first-pass run. So this predicate only FLAGS the markers;
    the method/clip policy lives in codex_lineage and this module echoes it rather than re-deriving.
    """
    min_dynrange = df_marker["p99_5_over_median"].min()
    prev = df_marker["prev_max"].mean()
    return bool(min_dynrange < DYNRANGE_MIN and prev < SPARSE_MAX_PREV)


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    out = []
    for m in SCREENED_MARKERS:
        pm = to_panel(m)
        d = df[df["marker"] == pm]
        method = cl.MARKER_INPUT_FLOOR_METHODS[pm]   # echo the configured policy, do not re-derive
        is_max = method == cl.MAX_LINEAR_OTSU_TRIANGLE_FLOOR
        floor_col = d.apply(lambda r: max(r.otsu, r.triangle) if is_max else r.otsu, axis=1)
        prev_col = d["prev_max"] if is_max else d["prev_otsu"]
        out.append(dict(
            marker=pm, method=method.replace("max_linear_otsu_triangle", "max").replace("linear_otsu", "otsu"),
            floor_median=round(float(floor_col.median()), 1),
            floor_spread=round(float(floor_col.max() / floor_col.min()), 1),
            prev_pct=round(float(prev_col.mean()), 2),
            min_dynrange=round(float(d["p99_5_over_median"].min()), 1),
            clip_sensitive=int(clip_sensitive(d)),
        ))
    return pd.DataFrame(out).sort_values("min_dynrange")


def check_against_observed(df: pd.DataFrame) -> list[str]:
    """Regression-assert measured RAW floors against codex_lineage.OBSERVED_FLOORS_RAW (operative)."""
    problems = []
    for m, expected in cl.OBSERVED_FLOORS_RAW.items():
        d = df[df["marker"] == to_panel(m)]
        if d.empty:
            problems.append(f"{m}: no measurement produced")
            continue
        method = cl.MARKER_INPUT_FLOOR_METHODS.get(m, cl.LINEAR_OTSU_FLOOR)
        floors = d.apply(lambda r: max(r.otsu, r.triangle)
                         if method == cl.MAX_LINEAR_OTSU_TRIANGLE_FLOOR else r.otsu, axis=1)
        got = float(floors.median())
        if abs(got - expected) / expected > 0.02:
            problems.append(f"{m}: measured RAW floor median {got:.1f} vs recorded {expected} (>2%)")
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parquet", type=Path, default=RAW_PARQUET,
                    help="cells parquet to floor on. Default is the RAW tree (the operative input).")
    ap.add_argument("--check", action="store_true",
                    help="regression-assert against OBSERVED_FLOORS_RAW and fail on >2%% drift")
    args = ap.parse_args()

    if not args.parquet.is_dir():
        print(f"[err] parquet dir not found: {args.parquet}")
        return 2
    tag = "RAW (operative)" if args.parquet == RAW_PARQUET else str(args.parquet)
    print(f"Measuring floors on {tag} for {len(SCREENED_MARKERS)} screened markers x {len(DONORS)} donors")

    df = measure_floors(args.parquet)
    summ = summarize(df)
    print("\n" + summ.to_string(index=False))
    flagged = sorted(summ.loc[summ["clip_sensitive"] == 1, "marker"])
    unclipped = sorted(cl.CLIP_UNCLIPPED_TARGETS)
    print(f"\n`clip_sensitive=1` markers: {flagged}")
    print("Their bright positive population sits above the p99.5 clip on the low-rate donors, so a")
    print("clipped 2-component Otsu bisects background. Policy (codex_lineage): CHGA rides max()@0.995")
    print(f"(first-pass); CLIP_UNCLIPPED_TARGETS = {unclipped} screen UNCLIPPED. max() is NOT used.")

    if args.check:
        problems = check_against_observed(df)
        if problems:
            print("\nCHECK FAILED:")
            for p in problems:
                print(f"  - {p}")
            return 1
        print("\nCHECK OK: measured RAW floors match OBSERVED_FLOORS_RAW within 2%.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
