#!/usr/bin/env python
"""
Reaggregate islets, rebuild trajectory, and compute Leiden clustering.

Standalone pipeline script that replaces the manual notebook workflow:
  1. Load single-cell H5AD
  2. Aggregate to islet-level (core + peri, min_cells=0, require_paired=True)
  3. Merge core + peri datasets
  4. Compute trajectory on core (neighbors → PAGA → UMAP → DPT)
  5. Validate trajectory (INS correlation, donor ordering)
  6. Compute Leiden clustering on merged (4 resolutions + UMAP)
  7. Save all outputs

Usage:
    conda activate scvi-env
    python scripts/reaggregate_islets.py
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy import stats
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.join(SCRIPT_DIR, '..')

# Add islet_analysis to path
sys.path.insert(0, os.path.join(PROJECT_ROOT, 'islet_analysis'))

# Default paths
DEFAULT_SC_H5AD = os.path.join(
    PROJECT_ROOT, 'single_cell_analysis',
    'CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad'
)
DEFAULT_ISLET_DIR = os.path.join(PROJECT_ROOT, 'islet_analysis')
DEFAULT_TRAJECTORY_OUT = os.path.join(PROJECT_ROOT, 'data', 'adata_ins_root.h5ad')


def reaggregate(sc_path, islet_dir, trajectory_out, min_cells=0):
    """Run full reaggregation pipeline."""
    import scanpy as sc
    from fixed_islet_aggregation import (
        create_islet_dataset_fixed,
        merge_islet_peri_datasets,
    )

    sc.settings.verbosity = 2

    # ---------------------------------------------------------------
    # Step 1: Load single-cell data
    # ---------------------------------------------------------------
    print("=" * 60)
    print("STEP 1: Loading single-cell data")
    print("=" * 60)
    adata_sc = sc.read_h5ad(sc_path)
    print(f"Shape: {adata_sc.shape}")
    print(f"Donors: {adata_sc.obs['imageid'].nunique()}")
    print(f"Donor status: {dict(adata_sc.obs['Donor Status'].value_counts())}")

    # ---------------------------------------------------------------
    # Step 2: Aggregate to islet-level (separate core + peri)
    # ---------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 2: Aggregating islets")
    print(f"  min_cells={min_cells}, require_paired=True")
    print("=" * 60)

    adata_islet, adata_peri = create_islet_dataset_fixed(
        adata_sc,
        region='separate',
        min_cells=min_cells,
        require_paired=True,
    )
    print(f"\nCore islets: {adata_islet.n_obs}")
    print(f"Peri islets: {adata_peri.n_obs}")

    # ---------------------------------------------------------------
    # Step 3: Merge core + peri
    # ---------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 3: Merging core + peri datasets")
    print("=" * 60)

    adata_merged = merge_islet_peri_datasets(adata_islet, adata_peri)
    print(f"Merged: {adata_merged.n_obs} islets")

    # ---------------------------------------------------------------
    # Step 4: Save aggregated H5ADs
    # ---------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 4: Saving aggregated H5ADs")
    print("=" * 60)

    core_path = os.path.join(islet_dir, 'islets_core_fixed.h5ad')
    peri_path = os.path.join(islet_dir, 'islets_peri_fixed.h5ad')
    merged_path = os.path.join(islet_dir, 'islets_merged_fixed.h5ad')

    adata_islet.write_h5ad(core_path)
    print(f"Saved: {core_path} ({adata_islet.n_obs} islets)")
    adata_peri.write_h5ad(peri_path)
    print(f"Saved: {peri_path} ({adata_peri.n_obs} islets)")
    adata_merged.write_h5ad(merged_path)
    print(f"Saved: {merged_path} ({adata_merged.n_obs} islets)")

    # ---------------------------------------------------------------
    # Step 5: Compute trajectory on core
    # ---------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 5: Computing trajectory")
    print("=" * 60)

    # Clear any existing embeddings
    for key in ['X_umap', 'X_diffmap', 'X_pca']:
        if key in adata_islet.obsm:
            del adata_islet.obsm[key]
    for key in ['neighbors', 'paga']:
        if key in adata_islet.uns:
            del adata_islet.uns[key]

    # Neighbors from scVI latent means
    if 'X_scVI_mean' in adata_islet.obsm:
        print("Computing neighbors from X_scVI_mean (batch-corrected)...")
        sc.pp.neighbors(adata_islet, n_neighbors=30, use_rep='X_scVI_mean', metric='cosine')
    else:
        print("WARNING: X_scVI_mean not available, falling back to PCA")
        sc.pp.pca(adata_islet, n_comps=10)
        sc.pp.neighbors(adata_islet, n_neighbors=15)

    # PAGA for topology
    sc.tl.paga(adata_islet, groups='donor_status')
    print("PAGA connectivity:")
    print(adata_islet.uns['paga']['connectivities'].toarray())

    # PAGA-initialized UMAP
    sc.pl.paga(adata_islet, plot=False)
    sc.tl.umap(adata_islet, init_pos='paga', min_dist=0.3, spread=1.5)
    print(f"UMAP computed: {adata_islet.obsm['X_umap'].shape}")

    # Find root: ND islet with highest INS
    ins_idx = list(adata_islet.var_names).index('INS')
    nd_mask = adata_islet.obs['donor_status'] == 'ND'
    if not nd_mask.any():
        raise ValueError("No ND islets found!")

    nd_indices = np.where(nd_mask)[0]
    ins_values = adata_islet.X[nd_mask, ins_idx]
    root_idx = nd_indices[np.argmax(ins_values)]
    adata_islet.uns['iroot'] = root_idx
    print(f"Root islet: index={root_idx}, INS={adata_islet.X[root_idx, ins_idx]:.3f}")

    # Diffusion map and pseudotime
    sc.tl.diffmap(adata_islet, n_comps=10)
    sc.tl.dpt(adata_islet)
    print(f"Pseudotime range: [{adata_islet.obs['dpt_pseudotime'].min():.3f}, "
          f"{adata_islet.obs['dpt_pseudotime'].max():.3f}]")

    adata_islet.obs['DC1'] = adata_islet.obsm['X_diffmap'][:, 0]
    adata_islet.obs['DC2'] = adata_islet.obsm['X_diffmap'][:, 1]

    # ---------------------------------------------------------------
    # Step 6: Validate trajectory
    # ---------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 6: Validating trajectory")
    print("=" * 60)

    pt = adata_islet.obs['dpt_pseudotime'].values
    valid = ~np.isnan(pt)

    # Check 1: INS negatively correlated with pseudotime
    ins_expr = adata_islet.X[:, ins_idx]
    r_ins, p_ins = stats.spearmanr(pt[valid], ins_expr[valid])
    pass_ins = r_ins < -0.3
    print(f"1. INS vs pseudotime: r={r_ins:.3f}, p={p_ins:.2e} -> {'PASS' if pass_ins else 'FAIL'}")

    # Check 2: GCG positively correlated
    pass_gcg = True
    if 'GCG' in adata_islet.var_names:
        gcg_idx = list(adata_islet.var_names).index('GCG')
        gcg_expr = adata_islet.X[:, gcg_idx]
        r_gcg, p_gcg = stats.spearmanr(pt[valid], gcg_expr[valid])
        pass_gcg = r_gcg > 0.2
        print(f"2. GCG vs pseudotime: r={r_gcg:.3f}, p={p_gcg:.2e} -> {'PASS' if pass_gcg else 'FAIL'}")

    # Check 3: Donor status ordering (ND < Aab+ < T1D)
    status_means = []
    for status in ['ND', 'Aab+', 'T1D']:
        mask = adata_islet.obs['donor_status'] == status
        if mask.any():
            mean_pt = pt[mask].mean()
            std_pt = pt[mask].std()
            status_means.append(mean_pt)
            print(f"   {status:5s}: {mean_pt:.3f} +/- {std_pt:.3f} (n={mask.sum()})")
    ordering_ok = all(status_means[i] <= status_means[i+1] for i in range(len(status_means)-1))
    print(f"3. Donor ordering ND < Aab+ < T1D: {'PASS' if ordering_ok else 'FAIL'}")

    # Check 4: 6533 included
    has_6533 = any('6533' in str(x) for x in adata_islet.obs['imageid'].unique())
    print(f"4. Donor 6533 included: {'PASS' if has_6533 else 'FAIL'}")

    all_pass = pass_ins and pass_gcg and ordering_ok and has_6533
    print(f"\nOverall: {'ALL PASSED' if all_pass else 'SOME CHECKS FAILED'}")

    if not all_pass:
        print("WARNING: Some validation checks failed. Outputs will still be saved.")

    # Visualization UMAP from raw marker expression
    # The scVI latent has donor-level variation corrected out (Age+Gender covariates
    # bijectively map to donors), which reduces disease-stage separation in UMAP.
    # Raw marker means (INS, GCG, immune markers) preserve this biological signal.
    print("\nComputing visualization UMAP from raw marker expression...")
    sc.pp.pca(adata_islet, n_comps=15)
    sc.pp.neighbors(adata_islet, n_neighbors=30, use_rep='X_pca', metric='euclidean',
                    key_added='neighbors_viz')
    sc.tl.umap(adata_islet, neighbors_key='neighbors_viz', min_dist=0.5, spread=2.0)
    print(f"Visualization UMAP computed: {adata_islet.obsm['X_umap'].shape}")

    # ---------------------------------------------------------------
    # Step 7: Save trajectory H5AD
    # ---------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 7: Saving trajectory H5AD")
    print("=" * 60)

    # Add combined_islet_id for Shiny compatibility
    if 'islet_id' in adata_islet.obs.columns and 'combined_islet_id' not in adata_islet.obs.columns:
        adata_islet.obs['combined_islet_id'] = adata_islet.obs['islet_id'].copy()

    # Also re-save the core file with trajectory data
    adata_islet.write_h5ad(core_path)
    print(f"Updated: {core_path} (with trajectory)")

    adata_islet.write_h5ad(trajectory_out)
    print(f"Saved: {trajectory_out} ({adata_islet.n_obs} islets)")

    # ---------------------------------------------------------------
    # Step 8: Compute Leiden clustering on merged dataset
    # ---------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 8: Computing Leiden clustering")
    print("=" * 60)

    # Use the core dataset for Leiden (same as original notebook)
    adata_clust = adata_islet.copy()

    # Recompute neighbors for clustering (same params as trajectory)
    if 'X_scVI_mean' in adata_clust.obsm:
        sc.pp.neighbors(adata_clust, use_rep='X_scVI_mean', n_neighbors=15, metric='cosine')
    else:
        sc.pp.pca(adata_clust, n_comps=10)
        sc.pp.neighbors(adata_clust, n_neighbors=15)

    # Compute UMAP for clustering visualization (different params from trajectory)
    sc.tl.umap(adata_clust, min_dist=0.3, spread=1.0)
    adata_clust.obs['umap_1'] = adata_clust.obsm['X_umap'][:, 0]
    adata_clust.obs['umap_2'] = adata_clust.obsm['X_umap'][:, 1]

    # Leiden at 4 resolutions
    resolutions = [0.3, 0.5, 0.8, 1.0]
    for res in resolutions:
        key = f'leiden_{res}'
        sc.tl.leiden(adata_clust, resolution=res, key_added=key)
        n_clusters = adata_clust.obs[key].nunique()
        print(f"Resolution {res}: {n_clusters} clusters")

    # Save clustered H5AD
    clust_path = os.path.join(islet_dir, 'islets_core_clustered.h5ad')
    adata_clust.write_h5ad(clust_path)
    print(f"Saved: {clust_path} ({adata_clust.n_obs} islets)")

    # ---------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)
    print(f"Islets aggregated: {adata_islet.n_obs}")
    print(f"\nOutputs:")
    print(f"  {core_path}")
    print(f"  {peri_path}")
    print(f"  {merged_path}")
    print(f"  {trajectory_out}")
    print(f"  {clust_path}")
    print(f"\nNext steps:")
    print(f"  python scripts/compute_neighborhood_metrics.py")
    print(f"  python scripts/extract_per_islet_cells.py")
    print(f"  python scripts/build_h5ad_for_app.py")


def main():
    parser = argparse.ArgumentParser(
        description="Reaggregate islets, compute trajectory, and Leiden clustering."
    )
    parser.add_argument(
        '--input', default=DEFAULT_SC_H5AD,
        help='Path to single-cell H5AD (default: single_cell_analysis/...h5ad)'
    )
    parser.add_argument(
        '--islet-dir', default=DEFAULT_ISLET_DIR,
        help='Directory for islet-level H5AD outputs (default: islet_analysis/)'
    )
    parser.add_argument(
        '--trajectory-out', default=DEFAULT_TRAJECTORY_OUT,
        help='Path for trajectory H5AD output (default: data/adata_ins_root.h5ad)'
    )
    parser.add_argument(
        '--min-cells', type=int, default=0,
        help='Minimum cells per islet (default: 0, keeps all islets with >=1 cell)'
    )
    args = parser.parse_args()

    # Validate input exists
    if not os.path.exists(args.input):
        print(f"ERROR: Single-cell H5AD not found: {args.input}")
        sys.exit(1)

    # Ensure output dirs exist
    os.makedirs(args.islet_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.trajectory_out), exist_ok=True)

    reaggregate(args.input, args.islet_dir, args.trajectory_out, args.min_cells)


if __name__ == '__main__':
    main()
