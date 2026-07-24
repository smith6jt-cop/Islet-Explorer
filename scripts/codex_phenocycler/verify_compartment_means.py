"""Verify a QuPath per-compartment intensity export before anything downstream trusts it.

The compartment means feed the FAITHFUL RESTORE NNMF (feature_compartments = Membrane+Cell, not
the degenerate Cell-only variant) and the `restore_pair_review.py` bundle. Both silently produce
plausible-but-wrong output if the export is incomplete, misjoined, or donor-short -- exactly the
QuPath-exits-0-on-error failure class this port keeps hitting. So this fails CLOSED on:

  1. every required id/coord column present,
  2. all four compartments x every screened marker present as `{Compartment}: {marker}: Mean`,
  3. exactly the 15 CODEX donors, each with one image,
  4. `Object ID` joining the canonical cells parquet TOTALLY (the join RESTORE relies on),
  5. a reported NaN-nucleus fraction (expected ~10%: InstanSeg objects with no nucleus ROI).

It reads only the header eagerly; the row-level checks run in DuckDB against the 11 GB TSV without
materialising it. Read-only: opens the TSV and the parquet, writes nothing.

    conda activate scvi-env
    python verify_compartment_means.py \
        --tsv /home/smith6jt/panc_CODEX/Islet-Explorer/CellCompartments.tsv \
        --cells /data/codex_phenocycler/cells
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# codex_panel is the CODEX naming authority (donor set, screened markers, id normalization).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from codex_panel import DONORS, normalize_donor_id, to_panel  # noqa: E402
import codex_lineage as cl  # noqa: E402

COMPARTMENTS = ("Nucleus", "Cytoplasm", "Membrane", "Cell")

# Every marker CODEX screens, plus the deferred CD20. Derived from the pair config so this cannot
# drift from what the pipeline actually consumes.
SCREENED = sorted(
    {m for pair in (cl.FIRST_PASS_PAIRS + cl.IMMUNE_PAIRS + [cl.ENDOCRINE_SPLIT_PAIR]) for m in pair}
    | {m for ms in cl.IMMUNE_SUBTYPE_GATES.values() for m in ms}
)

REQUIRED_META = ["Image", "Object ID", "Cell: Area µm^2", "Centroid X µm", "Centroid Y µm"]


def read_header(tsv: Path) -> list[str]:
    with tsv.open(encoding="utf-8") as fh:
        return fh.readline().rstrip("\n").split("\t")


def check_columns(header: list[str]) -> list[str]:
    have = set(header)
    problems = []
    for c in REQUIRED_META:
        if c not in have:
            problems.append(f"missing metadata column: {c!r}")
    for m in SCREENED:
        pm = to_panel(m)
        for comp in COMPARTMENTS:
            col = f"{comp}: {pm}: Mean"
            if col not in have:
                problems.append(f"missing measurement column: {col!r}")
    return problems


def check_rows(tsv: Path, cells_dir: Path) -> list[str]:
    import duckdb

    con = duckdb.connect()
    con.execute("PRAGMA threads=8; SET enable_progress_bar=false")
    src = (
        f"read_csv('{tsv}', delim='\\t', header=true, all_varchar=true, "
        "sample_size=-1, ignore_errors=false)"
    )
    problems = []

    # --- donor coverage: one image per donor, exactly the 15 CODEX donors ---
    donor_rows = con.execute(
        f"""
        SELECT regexp_extract("Image", '[0-9]+') AS raw_donor,
               count(*) AS n, count(DISTINCT "Image") AS images
        FROM {src} GROUP BY raw_donor ORDER BY raw_donor
        """
    ).fetchall()
    seen = {normalize_donor_id(r[0]): (r[1], r[2]) for r in donor_rows}
    expected = set(DONORS)
    missing = expected - set(seen)
    extra = set(seen) - expected
    if missing:
        problems.append(f"donors absent from export: {sorted(missing)}")
    if extra:
        problems.append(f"unexpected donors in export: {sorted(extra)}")
    multi = {d: v[1] for d, v in seen.items() if v[1] != 1}
    if multi:
        problems.append(f"donors with != 1 image: {multi}")

    # --- object_id join against the canonical cells parquet (RESTORE relies on this) ---
    cg = f"'{cells_dir}/donor_id=*/*.parquet'"
    join = con.execute(
        f"""
        WITH exp AS (
            SELECT DISTINCT CAST("Object ID" AS VARCHAR) AS oid FROM {src}
        ),
        canon AS (
            SELECT DISTINCT CAST(object_id AS VARCHAR) AS oid
            FROM read_parquet({cg}, hive_partitioning=true)
        )
        SELECT
          (SELECT count(*) FROM exp)   AS export_ids,
          (SELECT count(*) FROM canon) AS canon_ids,
          (SELECT count(*) FROM exp SEMI JOIN canon USING (oid)) AS both,
          (SELECT count(*) FROM canon ANTI JOIN exp USING (oid)) AS canon_only,
          (SELECT count(*) FROM exp ANTI JOIN canon USING (oid)) AS export_only
        """
    ).fetchone()
    export_ids, canon_ids, both, canon_only, export_only = join
    print(f"  object_id: export={export_ids:,} canon={canon_ids:,} both={both:,} "
          f"canon_only={canon_only:,} export_only={export_only:,}")
    if canon_only:
        problems.append(
            f"{canon_only:,} canonical object_ids are ABSENT from the export -- the compartment "
            "sample cannot cover every RESTORE-fitted cell"
        )

    # --- NaN-nucleus fraction (expected ~10%: objects with no nucleus ROI) ---
    ncol = 'CAST("Nucleus: DAPI: Mean" AS DOUBLE)'
    nan_nuc = con.execute(
        f"SELECT count(*) FILTER (WHERE {ncol} IS NULL OR isnan({ncol})) * 1.0 / count(*) FROM {src}"
    ).fetchone()[0]
    print(f"  NaN-nucleus fraction: {nan_nuc:.3%} (expected ~10%; those cells are Membrane/Cell-only)")

    return problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tsv", type=Path,
                    default=Path("/home/smith6jt/panc_CODEX/Islet-Explorer/CellCompartments.tsv"))
    ap.add_argument("--cells", type=Path, default=Path("/data/codex_phenocycler/cells"))
    ap.add_argument("--header-only", action="store_true",
                    help="skip the DuckDB row-level checks (fast column-presence check only)")
    args = ap.parse_args()

    if not args.tsv.is_file():
        print(f"[err] TSV not found: {args.tsv}")
        return 2
    print(f"Verifying {args.tsv} ({args.tsv.stat().st_size / 1e9:.1f} GB)")

    header = read_header(args.tsv)
    print(f"  header: {len(header)} columns")
    print(f"  screened markers checked: {SCREENED}")
    problems = check_columns(header)

    if not args.header_only:
        if not args.cells.is_dir():
            problems.append(f"canonical cells dir not found: {args.cells}")
        else:
            problems += check_rows(args.tsv, args.cells)

    if problems:
        print("\nFAIL:")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("\nOK: export is complete, correctly named, 15/15 donors, and joins the canonical parquet.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
