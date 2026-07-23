#!/usr/bin/env python
"""Stage 1 for CODEX — QuPath per-cell TSV -> hive-partitioned parquet.

Replaces `phenocycler.cells_parquet` for this cohort rather than patching it. That module opens
the input with a comma `csv.reader` and calls DuckDB `read_csv(...)` with no `delim`, so CODEX's
TAB-separated export yields **zero marker columns** — silently, because `ignore_errors=true`
swallows the mismatch. Rewriting the delimiter handling upstream would also mean carrying
Senior-specific column expectations; a small explicit adapter is cheaper and auditable.

Output matches the contract Stage 3 (RESTORE) reads:

    <data_dir>/cells/donor_id=<id>/data_0.parquet
        object_id            QuPath detection UUID -- the join key for every later stage
        donor_id             unpadded ('112', not '0112')
        <marker>             one float32 column per panel marker, RAW LINEAR MFI
        + geometry/morphology/region columns (see METADATA_COLUMNS)

Three CODEX specifics, all enforced rather than assumed:

1. **Raw linear MFI only.** RESTORE applies no transform to its threshold-fitting input, and
   `_norm`/`_log2r` are ratio-scale, so a log1p input silently changes the answer. This reads
   `Cell: <marker>: Mean` straight out of QuPath.
2. **Ki67 is `Median`, not `Mean`** — the sole exception in the export. `codex_panel.qupath_column()`
   handles it; a naive f-string over markers drops it.
3. **Donor ids are unpadded.** QuPath's `Image` is `0112`; Senior's `regexp_extract(Image,'[0-9]+')`
   preserves that and produces a `donor_id=0112` partition that never joins Islet-Explorer data.

Usage:
    conda activate scvi-env
    python build_cells_parquet.py --dry-run          # schema + counts, writes nothing
    python build_cells_parquet.py --data-dir /data/codex_phenocycler
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import duckdb

from codex_panel import (
    DONORS,
    NON_PANEL_CHANNELS,
    PANEL_31,
    qupath_column,
    to_panel,
)

# Default lives on /data: /home is ~97% full and the parquet + REDSEA scratch will not fit.
DEFAULT_DATA_DIR = Path("/data/codex_phenocycler")
CELLS_TSV = Path("/home/smith6jt/panc_CODEX/results/old_data_plots/all_data.tsv")

# Carried through from the QuPath export. `Parent` is what encodes islet membership
# (`Islet_N` = core, `Islet_N_exp20um` = peri) and drives cell_region downstream.
METADATA_COLUMNS = [
    ("Object ID", "object_id", "VARCHAR"),
    ("Parent", "parent", "VARCHAR"),
    ("Classification", "classification", "VARCHAR"),
    ("Centroid X µm", "centroid_x_um", "DOUBLE"),
    ("Centroid Y µm", "centroid_y_um", "DOUBLE"),
    ("Cell: Area µm^2", "cell_area_um2", "DOUBLE"),
    ("Nucleus: Area µm^2", "nucleus_area_um2", "DOUBLE"),
    ("Cell: Circularity", "cell_circularity", "DOUBLE"),
    ("Cell: Solidity", "cell_solidity", "DOUBLE"),
    ("Cell: Max diameter µm", "cell_max_diameter_um", "DOUBLE"),
    ("Cell: Min diameter µm", "cell_min_diameter_um", "DOUBLE"),
]

# Markers written out. The 31 analysis markers plus Ki67 (a PROCESS marker that still gets a
# positivity call, via a standalone bimodal fit rather than a RESTORE pair) and DAPI (segmentation
# QC). epCAM/KRT8-18 are deliberately EXCLUDED: acquired on only 8 of 15 donors, so including them
# would produce a ragged schema and invite a cohort-wide call that cannot be made.
MARKER_COLUMNS = PANEL_31 + ["DAPI", "Ki67"]
EXCLUDED_MARKERS = [m for m in NON_PANEL_CHANNELS if m not in MARKER_COLUMNS]


def _q(name: str) -> str:
    """Quote an identifier for DuckDB (QuPath column names contain spaces, colons and µ)."""
    return '"' + name.replace('"', '""') + '"'


def build_sql(tsv: Path) -> str:
    """One scan of the TSV, projecting to the Stage-1 schema."""
    proj = [
        # Unpadded donor id. CAST through BIGINT rather than lstrip so '0112' -> '112' is total.
        "CAST(CAST(regexp_extract(\"Image\", '[0-9]+') AS BIGINT) AS VARCHAR) AS donor_id",
        # The QuPath image name VERBATIM (zero-padded: '0112'). Required by
        # `redsea.donor_image()`, which reads a column literally named `image`, and it is what
        # the GeoJSON filename tag is derived from. Keep it distinct from donor_id -- conflating
        # the two is what produces a `donor_id=0112` partition that never joins.
        '"Image" AS image',
    ]
    for src, dst, typ in METADATA_COLUMNS:
        proj.append(f"CAST({_q(src)} AS {typ}) AS {_q(dst)}")
    for m in MARKER_COLUMNS:
        proj.append(f"CAST({_q(qupath_column(m))} AS FLOAT) AS {_q(to_panel(m))}")
    cols = ",\n       ".join(proj)
    return f"""
SELECT {cols}
FROM read_csv({_q(str(tsv)).replace('"', "'")},
              delim='\\t', header=true, sample_size=-1)
WHERE {_q('Object type')} = 'Cell'
"""


def verify_header(tsv: Path) -> list[str]:
    """Fail before a 1.4 GB scan if any expected column is absent."""
    header = tsv.open(encoding="utf-8").readline().rstrip("\n").split("\t")
    have = set(header)
    missing = [s for s, _, _ in METADATA_COLUMNS if s not in have]
    missing += [qupath_column(m) for m in MARKER_COLUMNS if qupath_column(m) not in have]
    if "Image" not in have or "Object type" not in have:
        missing += [c for c in ("Image", "Object type") if c not in have]
    if missing:
        raise SystemExit(f"ERROR: {tsv} is missing expected columns: {missing}")
    return header


def report_debris(con: duckdb.DuckDBPyConnection, relation: str, min_cell_area: float = 5.0) -> None:
    """Report QuPath segmentation debris rather than letting it propagate silently.

    Measured on this export: 47 cells across 13 donors are NaN on ALL 33 markers, with mean cell
    area 0.2 um2 vs 51.0 um2 for the rest — degenerate polygons with no pixels to average. This is
    what Senior's `restore.min_cell_area = 5.0` exists to remove.

    They are deliberately NOT dropped here. Stage 1 is the faithful raw export, and `object_id` has
    to stay a total join key back to QuPath and the canonical h5ad. `raw_otsu_floor()` already
    filters non-finite values, so the floors are safe; pass `--min-cell-area` if you want them gone
    from the parquet itself.
    """
    m0 = MARKER_COLUMNS[0]
    allnan = " AND ".join([f'isnan({_q(to_panel(m))})' for m in MARKER_COLUMNS])
    n_all, n_small = con.execute(
        f"SELECT sum(CASE WHEN {allnan} THEN 1 ELSE 0 END),"
        f"       sum(CASE WHEN cell_area_um2 < {min_cell_area} THEN 1 ELSE 0 END)"
        f" FROM {relation}"
    ).fetchone()
    if n_all or n_small:
        print(f"\nQC: {n_all or 0} cells are NaN on all {len(MARKER_COLUMNS)} markers; "
              f"{n_small or 0} have cell_area < {min_cell_area} um2 (QuPath debris).")
        print("    Retained — Stage 1 is the faithful raw export and object_id must stay a total")
        print("    join key. Floors filter non-finite values already; use --min-cell-area to drop.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells-tsv", type=Path, default=CELLS_TSV)
    ap.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--dry-run", action="store_true",
                    help="report schema and per-donor counts; write nothing")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--min-cell-area", type=float, default=0.0,
                    help="drop cells below this Cell:Area um^2 (QuPath debris). Default 0 = keep "
                         "all, so object_id stays a total join key; Senior's RESTORE uses 5.0.")
    args = ap.parse_args()

    if not args.cells_tsv.exists():
        raise SystemExit(f"ERROR: {args.cells_tsv} not found")
    verify_header(args.cells_tsv)
    print(f"Header OK: all {len(METADATA_COLUMNS)} metadata + {len(MARKER_COLUMNS)} marker columns present")
    if EXCLUDED_MARKERS:
        print(f"Excluded (not cohort-wide): {EXCLUDED_MARKERS}")

    out = args.data_dir / "cells"
    con = duckdb.connect()
    con.execute(f"PRAGMA threads={args.threads}")
    sql = build_sql(args.cells_tsv)
    if args.min_cell_area > 0:
        sql = f"SELECT * FROM ({sql}) WHERE cell_area_um2 >= {args.min_cell_area}"
        print(f"Filtering to cell_area_um2 >= {args.min_cell_area}")

    if args.dry_run:
        con.execute(f"CREATE VIEW cells AS {sql}")
        report_debris(con, "cells")
        n_cols = len(con.execute("DESCRIBE cells").fetchall())
        counts = con.execute(
            "SELECT donor_id, count(*) n FROM cells GROUP BY 1 ORDER BY CAST(donor_id AS BIGINT)"
        ).fetchall()
        total = sum(n for _, n in counts)
        print(f"\nSchema: {n_cols} columns")
        print(f"Cells:  {total:,} across {len(counts)} donors")
        for d, n in counts:
            print(f"    {d:>6}  {n:>9,}")
        unexpected = {d for d, _ in counts} - set(DONORS)
        missing = set(DONORS) - {d for d, _ in counts}
        if unexpected or missing:
            raise SystemExit(f"ERROR: donor mismatch. unexpected={sorted(unexpected)} missing={sorted(missing)}")
        print("\nDonor set matches codex_panel.DONORS exactly. (dry run: nothing written)")
        return 0

    if out.exists():
        if not args.overwrite:
            raise SystemExit(f"ERROR: {out} exists; pass --overwrite")
        shutil.rmtree(out)
    out.parent.mkdir(parents=True, exist_ok=True)

    print(f"\nWriting {out} ...")
    con.execute(f"""
        COPY ({sql}) TO '{out}'
        (FORMAT PARQUET, PARTITION_BY (donor_id), COMPRESSION ZSTD, FILENAME_PATTERN 'data_{{i}}')
    """)
    parts = sorted(p.name for p in out.glob("donor_id=*"))
    print(f"Wrote {len(parts)} partitions: {parts}")

    got = {p.split("=", 1)[1] for p in parts}
    if got != set(DONORS):
        raise SystemExit(f"ERROR: partitions {sorted(got)} != expected {DONORS}")
    n = con.execute(f"SELECT count(*) FROM read_parquet('{out}/**/*.parquet')").fetchone()[0]
    print(f"Verified: {n:,} cells, donor set matches codex_panel.DONORS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
