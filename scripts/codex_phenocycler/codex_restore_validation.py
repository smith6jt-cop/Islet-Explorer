"""Run upstream's locked v6 RESTORE pair screen on the CODEX cohort, unmodified.

This is a thin ADAPTER, not a fork: it builds CODEX-shaped inputs, overrides three upstream
module-level constants for CODEX's marker vocabulary and pair set, and then calls upstream
`restore_validation.evaluate_locked_cohort` verbatim. No scoring is reimplemented.

WHY EACH PIECE EXISTS (all verified against upstream 92f379a):
  * `patch_upstream()` -- upstream's `QUPATH_MARKERS`/`MARKER_INPUT_FLOOR_METHODS`/`CANDIDATE_PAIRS`
    are frozen for the Senior 59-plex. `compartment_column` raises for a marker not in QUPATH_MARKERS
    and `evaluate_locked_pair` raises for a pair not in CANDIDATE_PAIRS, so CODEX's markers and its 7
    non-Senior pairs must be injected. The CD20 both-role guard is inherited from upstream unchanged.
  * `build_sample()` -- upstream's own `extract_compartment_sample` derives `donor_id` as
    `regexp_extract(Image,'[0-9]+')`, which yields '0112' for donor 112 and drops it from the join.
    This builds the same sample schema but takes donor_id from the normalized id, reads CODEX's single
    combined `CellCompartments.tsv`, and DROPS the 47 all-NaN-Cell cells that would otherwise raise
    `nonfinite whole-cell raw value` in `load_validation_sample` and kill the whole run.
  * `build_raw_view()` -- `load_validation_sample` reads the raw parquet with columns
    `["object_id","X_centroid","Y_centroid"] + markers`; CODEX's cells parquet names the centroids
    `centroid_x_um`/`centroid_y_um`, so a renamed column view is required. `cells_redsea/` already
    satisfies `["object_id"] + markers` and is used directly.

The NNMF runs on the DEFAULT locked feature space (Membrane + Cell) -- NOT the degenerate Cell-only
variant -- because CellCompartments.tsv supplies real per-compartment means. `method_version` stays
`restore-pair-validation-v6` with no CODEX suffix, because this is the faithful config.

    conda activate scvi-env
    python codex_restore_validation.py --smoke 6516 6534      # build + load 2 donors, no scoring
    python codex_restore_validation.py --donors 6516 6534     # screen 2 donors, first-pass pairs
    python codex_restore_validation.py                        # full 15 donors, all 12 pairs
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

# codex modules + phenocycler BEFORE numpy (BLAS pinning).
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
_PA = Path(os.environ.get("PHENOCYCLER_ROOT",
           str(Path(__file__).resolve().parents[2] / "external" / "Phenocycler_Analysis")))
sys.path.insert(0, str(_PA))
import phenocycler  # noqa: E402,F401  BEFORE numpy
from phenocycler import restore_validation as rv  # noqa: E402

import codex_lineage as cl  # noqa: E402
from codex_panel import DONORS, to_panel  # noqa: E402

# --- CODEX screen configuration -----------------------------------------------------------------
# The 13 markers in the 12 screened pairs. NOT CD20 (excluded from RESTORE in both roles, 92f379a);
# CD20 rides along in CellCompartments.tsv for the deferred B-cell pass but takes no part here.
SCREENED_MARKERS = sorted({m for pair in (cl.FIRST_PASS_PAIRS + cl.IMMUNE_PAIRS
                                          + [cl.ENDOCRINE_SPLIT_PAIR]) for m in pair})
CODEX_PAIRS = [tuple(p) for p in (cl.FIRST_PASS_PAIRS + cl.IMMUNE_PAIRS + [cl.ENDOCRINE_SPLIT_PAIR])]
# Two screen classes, screened at DIFFERENT input-QC clip quantiles:
#   * first pass  (5 compartment pairs + endocrine split): the LOCKED 0.995 clip. These abundant
#     markers were never broken by the clip; unclipping only perturbs their divisors.
#   * immune      (CD3e/CD68/CD163/CD4/CD8a/HLADR): UNCLIPPED (IMMUNE_SCREEN_UPPER_QUANTILE = 1.0).
#     The 0.995 clip deletes the sparse positive tail on the 2-3 low-rate donors and drives the fit
#     MODEL_UNSTABLE; unclipped it is VALID and clean on all 15. See codex_lineage's resolution note.
FIRST_PASS_ONLY = [tuple(p) for p in (cl.FIRST_PASS_PAIRS + [cl.ENDOCRINE_SPLIT_PAIR])]
IMMUNE_ONLY = [tuple(p) for p in cl.IMMUNE_PAIRS]
LOCKED_UPPER_QUANTILE = 0.995

COMPARTMENTS = ("Nucleus", "Cytoplasm", "Membrane", "Cell")
TSV = Path("/home/smith6jt/panc_CODEX/Islet-Explorer/CellCompartments.tsv")
CELLS_RAW = Path("/data/codex_phenocycler/cells")
CELLS_REDSEA = Path("/data/codex_phenocycler/cells_redsea")
WORKDIR = Path("/data/codex_phenocycler/restore_validation")


def _upstream_sha() -> str | None:
    try:
        out = subprocess.run(["git", "-C", str(_PA), "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=10)
        return out.stdout.strip() if out.returncode == 0 else None
    except Exception:
        return None


def patch_upstream() -> None:
    """Inject CODEX's marker vocabulary, floor methods and pair set into upstream's frozen globals.

    Fails loudly if upstream has drifted off the pinned commit, so a future upstream change to the
    pair/marker contract cannot silently alter CODEX's screen.
    """
    sha = _upstream_sha()
    if sha and not sha.startswith(cl.UPSTREAM_PINNED_SHA):
        raise RuntimeError(
            f"upstream Phenocycler_Analysis is at {sha}, not the pinned {cl.UPSTREAM_PINNED_SHA}. "
            "Re-audit codex_lineage against the new upstream before screening, then bump "
            "cl.UPSTREAM_PINNED_SHA."
        )
    # Sentinel: the committed both-role exclusion must be present (guards against a stale checkout
    # even when git is unavailable).
    if not hasattr(rv, "RESTORE_EXCLUDED_MARKERS") or "CD20" not in getattr(rv, "RESTORE_EXCLUDED_MARKERS", {}):
        raise RuntimeError("upstream lacks RESTORE_EXCLUDED_MARKERS['CD20']; it is behind 92f379a")
    # METHOD_VERSION catches UNCOMMITTED method changes the SHA pin cannot see (upstream is developed
    # in the working tree). It bumps only on a real method change, so this is low-noise.
    if rv.METHOD_VERSION != cl.EXPECTED_METHOD_VERSION:
        raise RuntimeError(
            f"upstream METHOD_VERSION is {rv.METHOD_VERSION!r}, not the audited "
            f"{cl.EXPECTED_METHOD_VERSION!r}. Re-audit the method change (floors/pairs/divisor/clip?) "
            "then bump cl.EXPECTED_METHOD_VERSION. A version bump alone (e.g. a call refinement) is "
            "contract-safe; a floor/pair/divisor change is not."
        )

    rv.QUPATH_MARKERS = {m: m for m in SCREENED_MARKERS}                       # identity (CODEX names)
    rv.MARKER_INPUT_FLOOR_METHODS = {m: cl.MARKER_INPUT_FLOOR_METHODS[m] for m in SCREENED_MARKERS}
    rv.CANDIDATE_PAIRS = tuple(CODEX_PAIRS)

    # Upstream's own invariant (tests/test_restore_validation.py): one floor method per screened marker.
    if set(rv.MARKER_INPUT_FLOOR_METHODS) != set(rv.QUPATH_MARKERS):
        raise RuntimeError("patched MARKER_INPUT_FLOOR_METHODS and QUPATH_MARKERS disagree")
    cl.assert_pairs_valid(CODEX_PAIRS)   # re-run the CD20 both-role guard on the CODEX set


def _duckdb():
    import duckdb
    con = duckdb.connect()
    con.execute("PRAGMA threads=8; SET enable_progress_bar=false")
    return con


def _donor_filter(donors) -> str:
    return ", ".join("'" + str(d) + "'" for d in donors)


def build_sample(out_dir: Path, donors=DONORS, modulus: int = 20, *, force: bool = False) -> Path:
    """CODEX compartment sample: upstream's schema, CODEX donor-id, all-NaN-Cell cells dropped."""
    if out_dir.exists() and not force:
        return out_dir
    import shutil
    tmp = out_dir.with_name(out_dir.name + ".building")
    if tmp.exists():
        shutil.rmtree(tmp)

    meas = ",\n        ".join(
        f'TRY_CAST(s."{comp}: {to_panel(m)}: Mean" AS FLOAT) AS "{comp}__{to_panel(m)}"'
        for comp in COMPARTMENTS for m in SCREENED_MARKERS
    )
    # Drop cells with a NULL Cell-mean on any screened marker in the NEW export.
    not_null_tsv = " AND ".join(f's."Cell: {to_panel(m)}: Mean" IS NOT NULL' for m in SCREENED_MARKERS)
    # AND require the cell to be finite in the OLD canonical parquet (the raw axis load_validation_sample
    # allclose-checks against). 47 cells cohort-wide are all-marker-NaN in the old export (9 on 6534);
    # one in the ::20 sample raises `nonfinite whole-cell raw value` and kills the run. The raw parquet
    # stays the independent source, so the allclose remains a genuine new-vs-old cross-check.
    raw_finite = " AND ".join(f'(r."{to_panel(m)}" IS NOT NULL AND NOT isnan(r."{to_panel(m)}"))'
                              for m in SCREENED_MARKERS)
    cells_glob = f"'{CELLS_RAW}/donor_id=*/*.parquet'"

    sql = f"""
    COPY (
      WITH canon_finite AS (
        SELECT CAST(object_id AS VARCHAR) AS object_id
        FROM read_parquet({cells_glob}, hive_partitioning=true) r
        WHERE {raw_finite.replace('r."', '"')}
      )
      SELECT
        CAST(CAST(regexp_extract(s."Image", '[0-9]+') AS BIGINT) AS VARCHAR) AS donor_id,
        CAST(s."Image" AS VARCHAR)     AS image,
        CAST(s."Object ID" AS VARCHAR) AS object_id,
        TRY_CAST(s."Cell: Area µm^2" AS FLOAT) AS cell_area,
        TRY_CAST(s."Centroid X µm" AS FLOAT)   AS x_um,
        TRY_CAST(s."Centroid Y µm" AS FLOAT)   AS y_um,
        'CellCompartments.tsv' AS source_csv,
        {meas}
      FROM read_csv('{TSV}', delim='\\t', header=true, all_varchar=true, sample_size=-1) s
      INNER JOIN canon_finite c ON CAST(s."Object ID" AS VARCHAR) = c.object_id
      WHERE CAST(CAST(regexp_extract(s."Image", '[0-9]+') AS BIGINT) AS VARCHAR) IN ({_donor_filter(donors)})
        AND hash(CAST(s."Object ID" AS VARCHAR)) % {modulus} = 0
        AND {not_null_tsv}
    ) TO '{tmp}' (FORMAT PARQUET, PARTITION_BY (donor_id), COMPRESSION ZSTD)
    """
    con = _duckdb()
    con.execute(sql)
    con.close()
    if out_dir.exists():
        shutil.rmtree(out_dir)
    tmp.rename(out_dir)
    return out_dir


def build_raw_view(out_dir: Path, donors=DONORS, *, force: bool = False) -> Path:
    """Raw cells with centroids renamed to X_centroid/Y_centroid, as load_validation_sample expects."""
    if out_dir.exists() and not force:
        return out_dir
    import shutil
    tmp = out_dir.with_name(out_dir.name + ".building")
    if tmp.exists():
        shutil.rmtree(tmp)
    markers = ", ".join(f'"{to_panel(m)}"' for m in SCREENED_MARKERS)
    src = f"'{CELLS_RAW}/donor_id=*/*.parquet'"
    sql = f"""
    COPY (
      SELECT CAST(donor_id AS VARCHAR) AS donor_id,
             CAST(object_id AS VARCHAR) AS object_id,
             centroid_x_um AS X_centroid,
             centroid_y_um AS Y_centroid,
             {markers}
      FROM read_parquet({src}, hive_partitioning=true)
      WHERE CAST(donor_id AS VARCHAR) IN ({_donor_filter(donors)})
    ) TO '{tmp}' (FORMAT PARQUET, PARTITION_BY (donor_id), COMPRESSION ZSTD)
    """
    con = _duckdb()
    con.execute(sql)
    con.close()
    if out_dir.exists():
        shutil.rmtree(out_dir)
    tmp.rename(out_dir)
    return out_dir


def prepare_inputs(work: Path, donors=DONORS, modulus: int = 20, *, force: bool = False):
    work.mkdir(parents=True, exist_ok=True)
    sample = build_sample(work / "sample", donors, modulus, force=force)
    raw_view = build_raw_view(work / "raw_view", donors, force=force)
    return sample, raw_view, CELLS_REDSEA


def smoke(donors, work: Path, modulus: int = 20) -> None:
    """Build inputs and prove load_validation_sample accepts them -- the allclose gate, no scoring."""
    patch_upstream()
    sample, raw_view, redsea = prepare_inputs(work, donors, modulus, force=True)
    for d in donors:
        df, audit = rv.load_validation_sample(sample, raw_view, redsea, str(d))
        print(f"  donor {d}: loaded {len(df):,} sampled cells; "
              f"membrane_cell_complete={audit.get('membrane_cell_complete_fraction', 'n/a')}")
    print("SMOKE OK: sample/raw/redsea load and the Cell-mean allclose check passes.")


def run(donors, output_dir: Path, *, pairs=None, work: Path = WORKDIR, modulus: int = 20,
        upper_quantile: float = LOCKED_UPPER_QUANTILE,
        pair_reviews: Path | None = None, expression_reviews: Path | None = None,
        force_build: bool = False):
    patch_upstream()
    pairs = [tuple(p) for p in (pairs or CODEX_PAIRS)]
    cl.assert_pairs_valid(pairs)
    sample, raw_view, redsea = prepare_inputs(work, donors, modulus, force=force_build)
    return rv.evaluate_locked_cohort(
        sample, raw_view, redsea, output_dir,
        donors=[str(d) for d in donors], pairs=pairs,
        # DEFAULT = Membrane+Cell (faithful, not Cell-only). upper_quantile is the ONE declared CODEX
        # deviation from locked v6, and only for the immune run (1.0 vs 0.995) -- see codex_lineage.
        config=rv.PairValidationConfig(input_qc_upper_quantile=upper_quantile),
        pair_reviews=pair_reviews, expression_reviews=expression_reviews,
    )


def run_robust(donors, output_dir: Path, *, pairs, work: Path = WORKDIR, modulus: int = 20,
               upper_quantile: float = LOCKED_UPPER_QUANTILE,
               pair_reviews: Path | None = None):
    """Per-(donor,pair) screen that RECORDS a degenerate pair instead of crashing the cohort.

    `evaluate_locked_cohort` runs every pair in one try block, so a single degenerate arm (e.g.
    HLADR near-absent on a donor -> <2 cells in an exclusive arm at the unclipped floor) aborts the
    whole run. A cohort screen must instead report that pair and continue. Successful pairs use
    upstream `locked_pair_metrics` verbatim; a raise becomes an `INPUT_QC_DEGENERATE` row.
    """
    import pandas as pd
    patch_upstream()
    pairs = [tuple(p) for p in pairs]
    cl.assert_pairs_valid(pairs)
    sample, raw_view, redsea = prepare_inputs(work, donors, modulus)
    review_map, _ = rv.load_pair_reviews(pair_reviews)
    cfg = rv.PairValidationConfig(input_qc_upper_quantile=upper_quantile)
    rows = []
    for donor in donors:
        donor_df, _ = rv.load_validation_sample(sample, raw_view, redsea, str(donor))
        for target, reference in pairs:
            try:
                ev = rv.evaluate_locked_pair(
                    donor_df, str(donor), target, reference, config=cfg,
                    pair_review=review_map.get((target, reference), rv.PairReview.UNREVIEWED))
                rows.append(rv.locked_pair_metrics(ev))
            except Exception as e:
                rows.append({"method_version": rv.METHOD_VERSION, "donor": str(donor),
                             "target": target, "reference": reference,
                             "state": "INPUT_QC_DEGENERATE", "state_reason": str(e)[:120]})
    df = pd.DataFrame(rows)
    output_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_dir / "pair_evaluations.csv", index=False)
    return df


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--donors", nargs="+", default=None, help="donor ids (default all 15)")
    ap.add_argument("--first-pass", action="store_true",
                    help="screen only the 5 first-pass pairs + endocrine split, at the locked 0.995 clip")
    ap.add_argument("--immune", action="store_true",
                    help="screen only the 6 immune pairs, UNCLIPPED (IMMUNE_SCREEN_UPPER_QUANTILE)")
    ap.add_argument("--smoke", nargs="+", default=None, metavar="DONOR",
                    help="build + load-check the named donors, no scoring")
    ap.add_argument("--output", type=Path, default=WORKDIR / "screen")
    ap.add_argument("--work", type=Path, default=WORKDIR)
    ap.add_argument("--modulus", type=int, default=20)
    ap.add_argument("--pair-reviews", type=Path, default=None)
    ap.add_argument("--force-build", action="store_true")
    args = ap.parse_args()

    if args.smoke:
        smoke(args.smoke, args.work, args.modulus)
        return 0

    donors = args.donors or DONORS
    if args.immune:
        pairs, uq = IMMUNE_ONLY, cl.IMMUNE_SCREEN_UPPER_QUANTILE
    elif args.first_pass:
        pairs, uq = FIRST_PASS_ONLY, LOCKED_UPPER_QUANTILE
    else:
        pairs, uq = CODEX_PAIRS, LOCKED_UPPER_QUANTILE
    print(f"Screening {len(donors)} donors x {len(pairs)} pairs at clip={uq} "
          f"(upstream {_upstream_sha()}, pin {cl.UPSTREAM_PINNED_SHA})")
    if args.immune:
        # unclipped immune floors can empty an exclusive arm on near-absent-marker donors -> use the
        # robust per-pair runner so one degenerate pair is recorded, not fatal.
        pe = run_robust(donors, args.output, pairs=pairs, work=args.work, modulus=args.modulus,
                        upper_quantile=uq, pair_reviews=args.pair_reviews)
    else:
        results = run(donors, args.output, pairs=pairs, work=args.work, modulus=args.modulus,
                      upper_quantile=uq, pair_reviews=args.pair_reviews, force_build=args.force_build)
        pe = results.get("pair_evaluations")
    if pe is not None:
        print(f"\nwrote {len(pe)} pair-evaluation rows to {args.output}")
        if "state" in pe:
            print(pe["state"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
