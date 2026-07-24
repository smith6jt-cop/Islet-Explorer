"""Run upstream's locked v8 RESTORE pair screen on the CODEX cohort, unmodified.

This is a thin ADAPTER, not a fork: it builds CODEX-shaped inputs, overrides three upstream
module-level constants for CODEX's marker vocabulary and pair set, and then calls upstream's
per-pair scorer verbatim. No scoring is reimplemented.

LEVEL 2 (2026-07-24): the screen is now the 30-pair / 27-marker level-2 set (cl.ALL_SCREEN_PAIRS).
Because the input-QC clip is a WHOLE-CONFIG scalar, the screen runs as TWO cohort passes -- a 0.995
`--clipped` group and a 1.0 `--unclipped` group (targets whose bright tail sits above the clip; see
cl.CLIP_UNCLIPPED_TARGETS). Both use `run_cohort_full_robust` (a per-(donor,pair) try/except that
still writes the full `evaluate_locked_cohort` file set, incl. method_spec.json), because level-2 has
sparse targets that empty an arm on some donors. reference_selection + positivity consume the
CONCATENATION of the two screen dirs; each review bundle is built per-group.

WHY EACH PIECE EXISTS (verified against upstream v8 / 8db37af):
  * `patch_upstream()` -- upstream's `QUPATH_MARKERS`/`MARKER_INPUT_FLOOR_METHODS`/`CANDIDATE_PAIRS`
    are frozen for the Senior 59-plex. `compartment_column` raises for a marker not in QUPATH_MARKERS
    and `evaluate_locked_pair` raises for a pair not in CANDIDATE_PAIRS, so CODEX's markers and its
    non-Senior pairs must be injected. CODEX re-includes CD20 as a target: it rebinds rv.CANDIDATE_PAIRS
    to a set containing `CD20 <- ECAD`, so upstream's both-role guard (over upstream's own pairs) never
    blocks it at runtime. The sentinel that upstream STILL excludes CD20 is kept to guard the SHA pin.
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
    python codex_restore_validation.py --clipped              # 15 donors, 0.995 group -> screen_l2_clipped
    python codex_restore_validation.py --unclipped            # 15 donors, 1.0 group -> screen_l2_unclipped
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
# The level-2 screen: 27 markers over 30 pairs (cl.ALL_SCREEN_PAIRS = FIRST_PASS_PAIRS + LEVEL2_PAIRS
# + HORMONE_COMPARATOR_PAIRS). CD20 is now a screened TARGET (CD20 <- ECAD); the CD20 both-role guard
# is unwound on the CODEX side (cl.RESTORE_EXCLUDED_MARKERS is empty), while patch_upstream still
# asserts UPSTREAM keeps its exclusion (that assertion guards the SHA pin, not the CODEX policy).
SCREENED_MARKERS = sorted({m for pair in cl.ALL_SCREEN_PAIRS for m in pair})
CODEX_PAIRS = [tuple(p) for p in cl.ALL_SCREEN_PAIRS]

# The clip is a WHOLE-CONFIG scalar in upstream (PairValidationConfig.input_qc_upper_quantile), so a
# mixed-clip screen MUST run as TWO cohort passes. A pair goes in the UNCLIPPED (1.0) group iff its
# TARGET is in cl.CLIP_UNCLIPPED_TARGETS (its bright positive tail sits above the 0.995 clip on the
# low-rate donors); otherwise the LOCKED 0.995 group. reference_selection + positivity consume the
# CONCATENATION of the two screen dirs (both v8); each review bundle is built per-group (its
# method_spec.json records one clip, and build_review_bundle re-thresholds at that clip).
LOCKED_UPPER_QUANTILE = 0.995
CLIP_UNCLIPPED_PAIRS = [tuple(p) for p in CODEX_PAIRS if p[0] in cl.CLIP_UNCLIPPED_TARGETS]
CLIP_0995_PAIRS = [tuple(p) for p in CODEX_PAIRS if p[0] not in cl.CLIP_UNCLIPPED_TARGETS]

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
    # Sentinel: UPSTREAM's committed both-role exclusion must still be present (guards the SHA pin even
    # when git is unavailable). CODEX UNWINDS this at level 2 -- it re-includes CD20 as a target and
    # rebinds rv.CANDIDATE_PAIRS below to a set that contains CD20 <- ECAD, so upstream's import-time
    # guard (which runs over upstream's OWN CANDIDATE_PAIRS) never blocks the CODEX pair at runtime.
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


def run_cohort_full_robust(donors, output_dir: Path, *, pairs, work: Path = WORKDIR, modulus: int = 20,
                           upper_quantile: float = LOCKED_UPPER_QUANTILE,
                           pair_reviews: Path | None = None, expression_reviews: Path | None = None,
                           force_build: bool = False):
    """`evaluate_locked_cohort` with a per-(donor,pair) try/except, writing the FULL output contract.

    Two facts force this at level 2:
      1. `evaluate_locked_cohort` runs every (donor,pair) in ONE try block, so a single degenerate arm
         (a sparse target -- CD20/CD35/KRT14/SST/GAP43 -- emptying an exclusive arm on a donor) aborts
         the whole cohort. A level-2 screen has such arms by construction, so it must record that pair
         and continue.
      2. The retired `run_robust` survived (1) but wrote ONLY `pair_evaluations.csv`. The review bundle
         (`build_review_bundle`) requires a COMPLETE results dir -- `method_spec.json`,
         `expression_reviews.csv`, `pair_reviews.csv` -- and reference_selection needs the v8 columns.

    So this reproduces `evaluate_locked_cohort`'s file set VERBATIM (data_audit / pair_evaluations /
    threshold_sampling / seed_stability / expression_reviews / pair_reviews / method_spec.json /
    run_manifest.json), differing only in the per-pair try/except. A degenerate pair becomes an
    `INPUT_QC_DEGENERATE` row with NaN v8 columns -- dropped by build_reference_matrix.median() and
    skipped by the bundle (not VALID). Successful pairs use upstream `locked_pair_metrics` verbatim.
    """
    import json
    import shutil
    import pandas as pd
    from phenocycler.cohort import DONOR_EXCLUSIONS, ensure_eligible_donors

    patch_upstream()
    pairs = tuple((str(t), str(r)) for t, r in pairs)
    cl.assert_pairs_valid(pairs)
    if not pairs or len(set(pairs)) != len(pairs):
        raise ValueError("pairs must be a nonempty unique sequence")
    donors = ensure_eligible_donors([str(d) for d in donors], context="CODEX level-2 cohort screen")
    sample, raw_view, redsea = prepare_inputs(work, donors, modulus, force=force_build)
    review_map, review_table = rv.load_expression_reviews(expression_reviews)
    pair_review_map, pair_review_table = rv.load_pair_reviews(pair_reviews)

    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"output already exists: {output_dir}")
    temp_dir = output_dir.with_name(output_dir.name + ".tmp")
    if temp_dir.exists():
        raise FileExistsError(f"temporary output already exists: {temp_dir}")
    temp_dir.mkdir(parents=True)

    cfg = rv.PairValidationConfig(input_qc_upper_quantile=upper_quantile)
    audits, pair_rows, threshold_rows, stability_rows = [], [], [], []
    try:
        for donor in donors:
            donor_df, audit = rv.load_validation_sample(sample, raw_view, redsea, str(donor))
            audits.append(audit)
            degen = 0
            for target, reference in pairs:
                try:
                    ev = rv.evaluate_locked_pair(
                        donor_df, str(donor), target, reference, config=cfg,
                        expression_review=review_map.get((str(donor), target), rv.ExpressionReview.UNREVIEWED),
                        pair_review=pair_review_map.get((target, reference), rv.PairReview.UNREVIEWED))
                    pair_rows.append(rv.locked_pair_metrics(ev))
                    for axis, grid in (("target", ev["target_threshold_grid"]),
                                       ("reference", ev["reference_threshold_grid"])):
                        for row in grid.to_dict("records"):
                            threshold_rows.append({"donor": str(donor), "target": target,
                                                   "reference": reference, "axis": axis, **row})
                    for row in ev["stability_rows"]:
                        stability_rows.append({"donor": str(donor), "target": target,
                                               "reference": reference, **row})
                except Exception as e:
                    degen += 1
                    pair_rows.append({"method_version": rv.METHOD_VERSION, "donor": str(donor),
                                      "target": target, "reference": reference,
                                      "state": "INPUT_QC_DEGENERATE", "state_reason": str(e)[:200]})
            print(f"[locked-evaluate-robust] donor {donor}: {len(donor_df):,} cells, "
                  f"{degen} degenerate pair(s)", flush=True)

        tables = {
            "data_audit": pd.DataFrame(audits),
            "pair_evaluations": pd.DataFrame(pair_rows),
            "threshold_sampling": pd.DataFrame(threshold_rows),
            "seed_stability": pd.DataFrame(stability_rows),
            "expression_reviews": review_table,
            "pair_reviews": pair_review_table,
        }
        for name, table in tables.items():
            table.to_csv(temp_dir / f"{name}.csv", index=False)
        (temp_dir / "method_spec.json").write_text(json.dumps(cfg.specification(), indent=2) + "\n")
        (temp_dir / "run_manifest.json").write_text(json.dumps({
            "method_version": rv.METHOD_VERSION,
            "donors": list(donors),
            "cohort_exclusions": DONOR_EXCLUSIONS,
            "pairs": [list(p) for p in pairs],
            "input_qc_upper_quantile": upper_quantile,
            "robust_wrapper": "codex_restore_validation.run_cohort_full_robust",
        }, indent=2) + "\n")
        temp_dir.rename(output_dir)
        return tables
    except Exception:
        if temp_dir.exists():
            shutil.rmtree(temp_dir)
        raise


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--donors", nargs="+", default=None, help="donor ids (default all 15)")
    grp = ap.add_mutually_exclusive_group()
    grp.add_argument("--clipped", action="store_true",
                     help=f"screen the {len(CLIP_0995_PAIRS)} pairs whose target is NOT clip-sensitive, "
                          "at the locked 0.995 clip (incl. the level-1 gates + CHGA)")
    grp.add_argument("--unclipped", action="store_true",
                     help=f"screen the {len(CLIP_UNCLIPPED_PAIRS)} pairs whose target IS clip-sensitive, "
                          "UNCLIPPED (input_qc_upper_quantile=1.0)")
    ap.add_argument("--smoke", nargs="+", default=None, metavar="DONOR",
                    help="build + load-check the named donors, no scoring")
    ap.add_argument("--output", type=Path, default=None,
                    help="output dir (default screen_l2_clipped / screen_l2_unclipped under --work)")
    ap.add_argument("--work", type=Path, default=WORKDIR)
    ap.add_argument("--modulus", type=int, default=20)
    ap.add_argument("--pair-reviews", type=Path, default=None,
                    help="pair_reviews.csv to promote ACCEPTED pairs to VALID with divisors (post-review)")
    ap.add_argument("--force-build", action="store_true")
    args = ap.parse_args()

    if args.smoke:
        smoke(args.smoke, args.work, args.modulus)
        return 0

    if not (args.clipped or args.unclipped):
        ap.error("choose a clip group: --clipped or --unclipped (or --smoke DONOR...)")

    donors = args.donors or DONORS
    if args.unclipped:
        pairs, uq = CLIP_UNCLIPPED_PAIRS, cl.IMMUNE_SCREEN_UPPER_QUANTILE
        default_out = args.work / "screen_l2_unclipped"
    else:
        pairs, uq = CLIP_0995_PAIRS, LOCKED_UPPER_QUANTILE
        default_out = args.work / "screen_l2_clipped"
    output = args.output or default_out

    print(f"Screening {len(donors)} donors x {len(pairs)} pairs at clip={uq} "
          f"(upstream {_upstream_sha()}, pin {cl.UPSTREAM_PINNED_SHA}) -> {output}")
    # Both groups use the robust runner: level-2 has sparse targets in BOTH clip groups (CD20/GAP43 in
    # the 0.995 group; CD35/KRT14/SST in the unclipped group), so a degenerate arm must be recorded,
    # not fatal -- and the runner still writes method_spec.json for the review bundle.
    tables = run_cohort_full_robust(donors, output, pairs=pairs, work=args.work, modulus=args.modulus,
                                    upper_quantile=uq, pair_reviews=args.pair_reviews,
                                    force_build=args.force_build)
    pe = tables.get("pair_evaluations")
    if pe is not None:
        print(f"\nwrote {len(pe)} pair-evaluation rows to {output}")
        if "state" in pe:
            print(pe["state"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
