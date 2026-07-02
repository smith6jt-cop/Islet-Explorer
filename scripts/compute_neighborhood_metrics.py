#!/usr/bin/env python3
"""
Compute per-islet peri-islet neighborhood metrics from single-cell H5AD.

Reads the phenotyped single-cell H5AD (2.6M cells) and computes metrics
quantifying the cellular microenvironment around each islet (20µm expansion zone).

Metrics computed (4 categories):
  1. Peri-islet composition — proportion & count of each phenotype in _exp20um zone
  2. Immune infiltration — immune fractions (peri & core), peri/core ratio, ratios
  3. Enrichment z-scores — Poisson z comparing peri-islet vs tissue-wide proportion
  4. Distance metrics — min distance from islet centroid to nearest immune cells

Outputs:
  - data/neighborhood_metrics.csv  (one row per qualified islet; ~5,214 × ~100 cols)
  - data/signed_distance_kde.csv   (one row per islet-associated cell: phenotype,
      case_id, donor_status, signed_dist, cell_region, is_immune) — donor-aware
      source for the Spatial-tab Card C signed-distance KDE
  - data/cell_distance.csv         (one row per ALL cells: case_id, donor_status,
      cell_region, phenotype, signed_dist) — tissue-wide source for the Card C
      density + composition views (symlog distance axis reaches bulk tissue)

Usage:
    python scripts/compute_neighborhood_metrics.py
    python scripts/compute_neighborhood_metrics.py --input path/to/single_cell.h5ad --output data/neighborhood_metrics.csv
"""

import argparse
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Constants ──────────────────────────────────────────────────────────────

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.join(SCRIPT_DIR, "..")

DEFAULT_SC_H5AD = os.path.join(
    PROJECT_ROOT, "single_cell_analysis",
    "CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad"
)
DEFAULT_ISLET_H5AD = os.path.join(
    PROJECT_ROOT, "islet_analysis", "islets_core_fixed.h5ad"
)
DEFAULT_OUTPUT = os.path.join(PROJECT_ROOT, "data", "neighborhood_metrics.csv")
DEFAULT_KDE_OUTPUT = os.path.join(PROJECT_ROOT, "data", "signed_distance_kde.csv")
DEFAULT_CELL_DISTANCE_OUTPUT = os.path.join(PROJECT_ROOT, "data", "cell_distance.csv")

# Immune phenotypes — drives the AGGREGATE metrics where "what counts as
# immune" must be defined narrowly: immune_frac_peri/core, immune_ratio,
# tcell_density_peri, immune_count_peri/core, and the min_dist_immune_mean
# summary. Do NOT extend this list — it would make those aggregates
# meaningless (e.g. "fraction of peri cells that are 'Beta cell or B cell'").
IMMUNE_TYPES = [
    "CD8a Tcell", "CD4 Tcell", "T cell", "B cell",
    "Macrophage", "APCs", "Immune"
]

# Per-type metrics (enrich_z_*, min_dist_*) are computed for ALL phenotypes
# discovered in the data, minus this noise bucket. This replaces the earlier
# 7-phenotype ENRICH_TYPES subset so the app can expose phenotype-driven
# selectors uniformly across every cell type, not a hand-curated subset.
PER_TYPE_EXCLUDE = {"Unknown"}

# Lineage markers for the Card C "Marker" mode. Positivity = scimap-rescaled X
# > 0.5 (the same gate used to phenotype cells) — CODEX's equivalent of
# Senior's RESTORE _pos. Only markers present in the panel are emitted.
POS_MARKERS_CODEX = [
    "INS", "GCG", "SST", "CD3e", "CD8a", "CD4", "CD20", "CD68", "CD163",
    "HLADR", "CD31", "CD34", "SMA", "PDPN", "ECAD", "panCK",
]

PIXEL_SIZE_UM = 0.3774  # micrometers per pixel, matches app constant


def parse_parent(parent_str):
    """Parse Parent column to extract region type and islet ID.

    Returns (islet_name, region) where region is 'core' or 'peri'.
    Returns (None, None) for non-islet cells.
    """
    if not isinstance(parent_str, str):
        return None, None
    m = re.match(r"^(Islet_\d+)(_exp20um)?$", parent_str)
    if m:
        islet_name = m.group(1)
        region = "peri" if m.group(2) else "core"
        return islet_name, region
    return None, None


def write_kde_long_table(islet_cells, per_type_phenotypes, kde_output_path):
    """Per-cell signed-distance long table for the Spatial-tab KDE (Card C-right).

    One row per islet-associated cell (core+peri) of a non-Unknown phenotype.
    ``signed_dist`` is the QuPath 'Dist to Closest Islet' (negative = inside the
    islet core, positive = outside / peri-islet zone). The ``case_id`` column lets
    the app aggregate to donor level, which removes the pseudoreplication that
    pooling all cells introduced (e.g. one large ND donor dominating the pool).
    Regenerated in lockstep with neighborhood_metrics.csv so it can never drift
    from the current phenotyping.
    """
    keep = islet_cells[islet_cells["phenotype"].isin(per_type_phenotypes)].copy()
    safe = (keep["phenotype"].astype(str)
            .str.replace(" ", "_", regex=False)
            .str.replace("+", "plus", regex=False))
    signed = pd.to_numeric(keep["Dist to Closest Islet"], errors="coerce")
    out = pd.DataFrame({
        "phenotype":    safe.values,                          # safe-name: matches neighborhood cols + app selectors
        "case_id":      keep["imageid"].astype(str).values,   # donor id == app Case ID
        "donor_status": keep["Donor Status"].astype(str).values,
        "signed_dist":  signed.values,
        "cell_region":  keep["cell_region"].astype(str).values,
        "is_immune":    keep["phenotype"].isin(IMMUNE_TYPES).astype(int).values,
    })
    out = out[np.isfinite(out["signed_dist"])].reset_index(drop=True)
    os.makedirs(os.path.dirname(kde_output_path), exist_ok=True)
    out.to_csv(kde_output_path, index=False)
    print(f"   Signed-distance KDE long table written to: {kde_output_path}")
    print(f"   Shape: {out.shape}  "
          f"(phenotypes: {out['phenotype'].nunique()}, donors: {out['case_id'].nunique()})")
    return out


def build_marker_positivity(sc):
    """Per-cell marker positivity from scimap-rescaled X (gate 0.5 — the same
    threshold that phenotyped the cells). CODEX's equivalent of Senior's RESTORE
    _pos. Returns an int8 DataFrame of {marker}_pos columns in sc row order, or
    None if no target marker is in the panel.
    """
    avail = [m for m in POS_MARKERS_CODEX if m in list(sc.var_names)]
    if not avail:
        print("   [marker positivity] none of POS_MARKERS_CODEX in panel; skipping")
        return None
    Xm = sc[:, avail].to_memory().X
    Xm = Xm.toarray() if hasattr(Xm, "toarray") else np.asarray(Xm)
    df = pd.DataFrame({f"{m}_pos": (Xm[:, i] > 0.5).astype("int8")
                       for i, m in enumerate(avail)})
    fr = ", ".join(f"{m}={df[m + '_pos'].mean():.3f}" for m in avail)
    print(f"   [marker positivity] {len(avail)} markers @ X>0.5 | pos frac: {fr}")
    return df


def write_cell_distance_table(obs, marker_pos, cell_distance_output_path):
    """All-cells per-cell signed-distance table for Spatial Card C (density +
    composition views), including per-cell marker positivity.

    One row per cell with a finite QuPath 'Dist to Closest Islet'. Unlike
    signed_distance_kde.csv (islet-associated core+peri only), this spans ALL
    cells — islet interior (negative) out to bulk tissue (large positive) — so
    the symlog distance axis reaches the full range. ``cell_region`` is
    core/peri for islet-associated cells, 'tissue' otherwise. Phenotype labels
    are the raw app names (match PHENOTYPE_COLORS in the R app); ``{marker}_pos``
    columns (int 0/1, scimap X>0.5) drive the app's "Marker" mode.
    """
    region = obs["cell_region"].where(obs["cell_region"].notna(), "tissue")
    signed = pd.to_numeric(obs["Dist to Closest Islet"], errors="coerce")
    out = pd.DataFrame({
        "case_id":      obs["imageid"].astype(str).values,
        "donor_status": obs["Donor Status"].astype(str).values,
        "cell_region":  region.astype(str).values,
        "phenotype":    obs["phenotype"].astype(str).values,
        "signed_dist":  signed.values.astype("float32"),
    })
    if marker_pos is not None:
        for c in marker_pos.columns:
            out[c] = marker_pos[c].values
    out = out[np.isfinite(out["signed_dist"])].reset_index(drop=True)
    os.makedirs(os.path.dirname(cell_distance_output_path), exist_ok=True)
    out.to_csv(cell_distance_output_path, index=False)
    print(f"   Cell-distance table written to: {cell_distance_output_path}")
    print(f"   Shape: {out.shape}  "
          f"(phenotypes: {out['phenotype'].nunique()}, donors: {out['case_id'].nunique()})")
    return out


def compute_metrics(sc_path, islet_path, output_path, kde_output_path=None,
                    cell_distance_output_path=None):
    """Main computation: read single-cell data, compute per-islet neighborhood metrics."""
    import anndata as ad

    print("=" * 60)
    print("Computing Peri-Islet Neighborhood Metrics")
    print("=" * 60)

    # ── 1. Load qualified islet list ──────────────────────────────────────
    print(f"\n1. Loading qualified islets: {islet_path}")
    islet_adata = ad.read_h5ad(islet_path)
    qualified = islet_adata.obs[["imageid", "base_islet_id"]].copy()
    qualified["imageid"] = qualified["imageid"].astype(str)
    qualified["base_islet_id"] = qualified["base_islet_id"].astype(str)
    qualified["combined_islet_id"] = (
        qualified["imageid"] + "_" + qualified["base_islet_id"]
    )
    print(f"   Qualified islets: {len(qualified)}")

    # ── 2. Load single-cell H5AD (backed to manage memory) ───────────────
    print(f"\n2. Loading single-cell H5AD: {sc_path}")
    sc = ad.read_h5ad(sc_path, backed="r")
    print(f"   Total cells: {sc.shape[0]:,}")
    print(f"   Markers: {sc.shape[1]}")

    # ── 3. Parse Parent column → islet assignment + region ───────────────
    print("\n3. Parsing cell-islet assignments from Parent column...")
    kde_obs_cols = ["Dist to Closest Islet", "Donor Status"]
    missing_kde = [c for c in kde_obs_cols if c not in sc.obs.columns]
    if missing_kde:
        raise KeyError(
            "single-cell obs is missing columns needed for the signed-distance "
            f"KDE output: {missing_kde}"
        )
    obs = sc.obs[["imageid", "Parent", "phenotype", "X_centroid", "Y_centroid",
                  "Dist to Closest Islet", "Donor Status"]].copy()
    obs["imageid"] = obs["imageid"].astype(str)
    obs["Parent"] = obs["Parent"].astype(str)
    obs["phenotype"] = obs["phenotype"].astype(str)

    parsed = obs["Parent"].apply(parse_parent)
    obs["islet_name"] = [p[0] for p in parsed]
    obs["cell_region"] = [p[1] for p in parsed]

    # Filter to islet-associated cells only
    islet_cells = obs[obs["islet_name"].notna()].copy()
    islet_cells["combined_islet_id"] = (
        islet_cells["imageid"] + "_" + islet_cells["islet_name"]
    )
    print(f"   Islet-associated cells: {len(islet_cells):,}")
    print(f"     Core: {(islet_cells['cell_region'] == 'core').sum():,}")
    print(f"     Peri: {(islet_cells['cell_region'] == 'peri').sum():,}")

    # Filter to qualified islets only
    qualified_set = set(qualified["combined_islet_id"])
    islet_cells = islet_cells[
        islet_cells["combined_islet_id"].isin(qualified_set)
    ].copy()
    print(f"   Cells in qualified islets: {len(islet_cells):,}")

    # ── 4. Compute tissue-wide phenotype proportions (baseline) ──────────
    print("\n4. Computing tissue-wide phenotype baseline...")
    all_phenotypes = sorted(obs["phenotype"].unique())
    tissue_total = len(obs)
    tissue_counts = obs["phenotype"].value_counts()
    tissue_props = (tissue_counts / tissue_total).to_dict()
    print(f"   {len(all_phenotypes)} phenotypes, {tissue_total:,} total cells")
    # Per-type metrics span all phenotypes minus PER_TYPE_EXCLUDE.
    per_type_phenotypes = [p for p in all_phenotypes if p not in PER_TYPE_EXCLUDE]
    print(f"   Per-type metrics (enrich_z, min_dist) computed for "
          f"{len(per_type_phenotypes)} phenotypes")

    # ── 5. Compute per-islet metrics ─────────────────────────────────────
    print("\n5. Computing per-islet neighborhood metrics...")
    results = []

    peri_cells = islet_cells[islet_cells["cell_region"] == "peri"]
    core_cells = islet_cells[islet_cells["cell_region"] == "core"]

    peri_grouped = peri_cells.groupby("combined_islet_id")
    core_grouped = core_cells.groupby("combined_islet_id")

    for _, row in qualified.iterrows():
        cid = row["combined_islet_id"]
        rec = {"combined_islet_id": cid}

        # --- Peri-islet composition ---
        if cid in peri_grouped.groups:
            peri = peri_grouped.get_group(cid)
            peri_total = len(peri)
            peri_pheno_counts = peri["phenotype"].value_counts()
        else:
            peri = pd.DataFrame()
            peri_total = 0
            peri_pheno_counts = pd.Series(dtype=int)

        rec["total_cells_peri"] = peri_total

        for pheno in all_phenotypes:
            safe_name = pheno.replace(" ", "_").replace("+", "plus")
            cnt = int(peri_pheno_counts.get(pheno, 0))
            rec[f"peri_count_{safe_name}"] = cnt
            rec[f"peri_prop_{safe_name}"] = cnt / peri_total if peri_total > 0 else np.nan

        # --- Core composition (for ratio metrics) ---
        if cid in core_grouped.groups:
            core = core_grouped.get_group(cid)
            core_total = len(core)
            core_pheno_counts = core["phenotype"].value_counts()
        else:
            core = pd.DataFrame()
            core_total = 0
            core_pheno_counts = pd.Series(dtype=int)

        rec["total_cells_core"] = core_total

        # --- Immune infiltration metrics ---
        peri_immune = sum(
            int(peri_pheno_counts.get(t, 0)) for t in IMMUNE_TYPES
        )
        core_immune = sum(
            int(core_pheno_counts.get(t, 0)) for t in IMMUNE_TYPES
        )

        rec["immune_count_peri"] = peri_immune
        rec["immune_count_core"] = core_immune
        rec["immune_frac_peri"] = (
            peri_immune / peri_total if peri_total > 0 else np.nan
        )
        rec["immune_frac_core"] = (
            core_immune / core_total if core_total > 0 else np.nan
        )
        # Peri/core immune ratio
        if core_immune > 0 and peri_total > 0:
            rec["immune_ratio"] = (peri_immune / peri_total) / (
                core_immune / core_total
            )
        else:
            rec["immune_ratio"] = np.nan

        # CD8/macrophage ratio in peri-islet
        cd8_peri = int(peri_pheno_counts.get("CD8a Tcell", 0))
        macro_peri = int(peri_pheno_counts.get("Macrophage", 0))
        rec["cd8_to_macro_ratio"] = (
            cd8_peri / macro_peri if macro_peri > 0 else np.nan
        )

        # T-cell density in peri-islet (T cells per 100 peri cells)
        tcell_peri = sum(
            int(peri_pheno_counts.get(t, 0))
            for t in ["CD8a Tcell", "CD4 Tcell", "T cell"]
        )
        rec["tcell_density_peri"] = (
            100 * tcell_peri / peri_total if peri_total > 0 else np.nan
        )

        # --- Enrichment z-scores (Poisson model) — per-phenotype ---
        for etype in per_type_phenotypes:
            safe_name = etype.replace(" ", "_").replace("+", "plus")
            observed = int(peri_pheno_counts.get(etype, 0))
            expected = tissue_props.get(etype, 0) * peri_total
            if expected > 0 and peri_total > 0:
                # Poisson z-score: (observed - expected) / sqrt(expected)
                rec[f"enrich_z_{safe_name}"] = (observed - expected) / np.sqrt(
                    expected
                )
            else:
                rec[f"enrich_z_{safe_name}"] = np.nan

        # --- Distance metrics (min distance from islet core centroid to
        #     nearest cell of each phenotype in the peri-islet zone) ---
        if peri_total > 0 and len(core) > 0:
            # Use core centroid as reference point
            core_cx = core["X_centroid"].astype(float).mean()
            core_cy = core["Y_centroid"].astype(float).mean()

            # Aggregate "any immune cell" min distance — uses IMMUNE_TYPES.
            peri_immune_mask = peri["phenotype"].isin(IMMUNE_TYPES)
            peri_immune_cells = peri[peri_immune_mask]
            if len(peri_immune_cells) > 0:
                dx = peri_immune_cells["X_centroid"].astype(float).values - core_cx
                dy = peri_immune_cells["Y_centroid"].astype(float).values - core_cy
                rec["min_dist_immune_mean"] = float(np.min(np.sqrt(dx**2 + dy**2)))
            else:
                rec["min_dist_immune_mean"] = np.nan

            # Per-phenotype min distances — extends across all phenotypes.
            peri_phenos = peri["phenotype"].values
            peri_x = peri["X_centroid"].astype(float).values
            peri_y = peri["Y_centroid"].astype(float).values
            for dtype in per_type_phenotypes:
                safe_name = dtype.replace(" ", "_").replace("+", "plus")
                type_mask = peri_phenos == dtype
                if type_mask.any():
                    tdx = peri_x[type_mask] - core_cx
                    tdy = peri_y[type_mask] - core_cy
                    rec[f"min_dist_{safe_name}"] = float(
                        np.min(np.sqrt(tdx**2 + tdy**2))
                    )
                else:
                    rec[f"min_dist_{safe_name}"] = np.nan
        else:
            rec["min_dist_immune_mean"] = np.nan
            for dtype in per_type_phenotypes:
                safe_name = dtype.replace(" ", "_").replace("+", "plus")
                rec[f"min_dist_{safe_name}"] = np.nan

        results.append(rec)

    df = pd.DataFrame(results)

    # ── 6. Summary and output ────────────────────────────────────────────
    print(f"\n6. Results summary:")
    print(f"   Total islets: {len(df)}")
    has_peri = df["total_cells_peri"] > 0
    print(f"   Islets with peri-islet data: {has_peri.sum()} ({100*has_peri.mean():.1f}%)")
    print(f"   Islets without peri data: {(~has_peri).sum()}")
    print(f"   Columns: {len(df.columns)}")

    # Quick sanity check
    if "immune_frac_peri" in df.columns:
        ifp = df.loc[has_peri, "immune_frac_peri"]
        print(f"   immune_frac_peri: mean={ifp.mean():.4f}, median={ifp.median():.4f}")

    # Write output
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"\n   Written to: {output_path}")
    print(f"   Shape: {df.shape}")

    # ── 7. Signed-distance KDE long table (per-cell, donor-aware) ─────────
    if kde_output_path:
        print("\n7. Writing signed-distance KDE long table...")
        write_kde_long_table(islet_cells, per_type_phenotypes, kde_output_path)

    # ── 8. All-cells cell-distance table (tissue-wide, Card C density/comp) ─
    if cell_distance_output_path:
        print("\n8. Writing all-cells cell-distance table...")
        marker_pos = build_marker_positivity(sc)
        write_cell_distance_table(obs, marker_pos, cell_distance_output_path)

    print(f"\n{'='*60}")
    print("DONE")
    print(f"{'='*60}")

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute per-islet peri-islet neighborhood metrics"
    )
    parser.add_argument(
        "--input", default=DEFAULT_SC_H5AD,
        help="Path to single-cell H5AD (phenotyped)"
    )
    parser.add_argument(
        "--islets", default=DEFAULT_ISLET_H5AD,
        help="Path to islets_core_fixed.h5ad (qualified islet list)"
    )
    parser.add_argument(
        "--output", default=DEFAULT_OUTPUT,
        help="Output path for neighborhood_metrics.csv"
    )
    parser.add_argument(
        "--kde-output", default=DEFAULT_KDE_OUTPUT,
        help="Output path for the per-cell signed-distance KDE long table"
    )
    parser.add_argument(
        "--cell-distance-output", default=DEFAULT_CELL_DISTANCE_OUTPUT,
        help="Output path for the all-cells tissue-wide cell-distance table"
    )
    args = parser.parse_args()

    for path, label in [
        (args.input, "single-cell H5AD"),
        (args.islets, "islets_core_fixed.h5ad"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    compute_metrics(args.input, args.islets, args.output, args.kde_output,
                    args.cell_distance_output)
