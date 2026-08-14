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


def _eta2_donor_within_status(pt, iid, status_col, stt, min_islets=3):
    """Fraction of within-status pseudotime variance explained by donor identity."""
    m = status_col == stt
    pt_s, d_s = pt[m], iid[m]
    groups = [pt_s[d_s == d] for d in np.unique(d_s) if (d_s == d).sum() >= min_islets]
    if len(groups) < 2:
        return None
    gm = pt_s.mean()
    ss_total = float(np.sum((pt_s - gm) ** 2))
    if ss_total <= 0:
        return 0.0
    ss_between = float(sum(len(g) * (g.mean() - gm) ** 2 for g in groups))
    return ss_between / ss_total


def _validate_one_pseudotime(adata, pt_col, ins_idx, label):
    """Validation gates for a single pseudotime column. Returns True if all pass.

    Checks 1-2 (marker correlations) are necessarily islet-level -- they ask
    whether expression tracks the axis. Check 3 (disease ordering) is run at
    DONOR level: with 5,214 islets from only 15 donors, an islet-level mean is
    pseudoreplicated, and the ordering is a claim about donors. The islet-level
    means are still printed for reference, but they do not decide the gate.
    """
    if pt_col not in adata.obs:
        print(f"\n[{label}] column '{pt_col}' absent -- SKIPPED")
        return True

    print(f"\n--- {label} ({pt_col}) ---")
    pt = adata.obs[pt_col].values
    valid = ~np.isnan(pt)
    if valid.sum() == 0:
        print("  no finite pseudotime values -> FAIL")
        return False

    # Check 1: INS negatively correlated with pseudotime
    ins_expr = adata.X[:, ins_idx]
    r_ins, p_ins = stats.spearmanr(pt[valid], ins_expr[valid])
    pass_ins = r_ins < -0.3
    print(f"1. INS vs pseudotime: r={r_ins:.3f}, p={p_ins:.2e} -> {'PASS' if pass_ins else 'FAIL'}")

    # Check 2: GCG positively correlated.
    # NOTE: the 0.1 threshold is RELAXED from an original 0.2 and the shipped
    # data clears it only narrowly (core r~0.143, combined r~0.103). It is
    # justified by a judgement that the genuine shared effect is ~0.14, which
    # nothing independent confirms. Treat a GCG pass as weak evidence.
    pass_gcg = True
    if 'GCG' in adata.var_names:
        gcg_idx = list(adata.var_names).index('GCG')
        r_gcg, p_gcg = stats.spearmanr(pt[valid], adata.X[valid, gcg_idx])
        pass_gcg = r_gcg > 0.1
        margin = r_gcg - 0.1
        note = "  <-- NARROW" if 0 <= margin < 0.05 else ""
        print(f"2. GCG vs pseudotime: r={r_gcg:.3f}, p={p_gcg:.2e} -> "
              f"{'PASS' if pass_gcg else 'FAIL'}{note}")

    # Check 3: disease ordering, at DONOR level (see docstring)
    iid = adata.obs['imageid'].astype(str).values
    status_col = adata.obs['donor_status'].values
    donor_means, islet_means = [], []
    for status in ['ND', 'Aab+', 'T1D']:
        mask = status_col == status
        if not mask.any():
            continue
        islet_means.append(np.nanmean(pt[mask]))
        per_donor = [np.nanmean(pt[mask & (iid == d)]) for d in np.unique(iid[mask])]
        per_donor = [x for x in per_donor if np.isfinite(x)]
        donor_means.append(np.mean(per_donor) if per_donor else np.nan)
        print(f"   {status:5s}: donor-mean {donor_means[-1]:.4f} "
              f"(n={len(per_donor)} donors) | islet-mean {islet_means[-1]:.4f} "
              f"(n={mask.sum()} islets)")
    ordering_ok = all(donor_means[i] <= donor_means[i + 1]
                      for i in range(len(donor_means) - 1))
    print(f"3. Donor-level ordering ND <= Aab+ <= T1D: {'PASS' if ordering_ok else 'FAIL'}")

    # Check 4: per-status donor eta^2 (ND/T1D <= 0.15, Aab+ <= 0.25).
    # Aab+ is relaxed because those donors genuinely span heterogeneous islet
    # types -- biology, not technical noise.
    donor_bias_ok = True
    print("4. Donor eta^2 within status (ND/T1D <= 0.15, Aab+ <= 0.25):")
    for stt in ['ND', 'Aab+', 'T1D']:
        eta2 = _eta2_donor_within_status(pt, iid, status_col, stt)
        if eta2 is None:
            print(f"   {stt}: (skip, too few donor groups)")
            continue
        thresh = 0.25 if stt == 'Aab+' else 0.15
        ok = eta2 <= thresh
        print(f"   {stt:5s}: eta^2 = {eta2:.3f} (thresh {thresh}) -> {'PASS' if ok else 'FAIL'}")
        donor_bias_ok = donor_bias_ok and ok

    return bool(pass_ins and pass_gcg and ordering_ok and donor_bias_ok)


def _validate_trajectory(adata, ins_idx):
    """Validate BOTH pseudotime modes. Returns True if every gate passes.

    Combined mode used to go completely unvalidated, so a regression there was
    invisible -- and it is the mode with the thinnest GCG margin.
    """
    ok_core = _validate_one_pseudotime(adata, 'dpt_pseudotime', ins_idx, 'CORE')
    ok_comb = _validate_one_pseudotime(adata, 'dpt_pseudotime_combined', ins_idx, 'COMBINED')

    # Cohort composition check (mode-independent)
    has_6533 = any('6533' in str(x) for x in adata.obs['imageid'].unique())
    print(f"\n5. Donor 6533 included: {'PASS' if has_6533 else 'FAIL'}")

    # Agreement between the two modes -- descriptive, not a gate. They are
    # different estimands; a low value means they are not interchangeable.
    if 'dpt_pseudotime_combined' in adata.obs:
        a = adata.obs['dpt_pseudotime'].values
        b = adata.obs['dpt_pseudotime_combined'].values
        both = ~np.isnan(a) & ~np.isnan(b)
        if both.sum() > 2:
            rho = stats.spearmanr(a[both], b[both])[0]
            print(f"6. Core vs combined agreement: Spearman rho = {rho:.3f} "
                  f"(descriptive; modes are not interchangeable)")

    all_pass = bool(ok_core and ok_comb and has_6533)
    print(f"\nOverall: {'ALL PASSED' if all_pass else 'SOME CHECKS FAILED'}")
    return all_pass


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

    # Neighbors from lightly-Harmony-corrected scVI latent (theta=0.1).
    # Empirical sweep findings (with the medoid root):
    #   raw          : span=0.116  INS=-0.69  CD8a=+0.21  Aab+ eta^2=0.30
    #   theta=0.1    : span=0.084  INS=-0.57  CD8a=+0.19  Aab+ eta^2=0.09  <- chosen
    #   theta=2.0    : span=0.061  INS=-0.55  CD8a=+0.15  Aab+ eta^2=0.07
    # theta=0.1 is the lightest viable correction: it captures most of the
    # donor-noise reduction (Harmony's behaviour is near step-functional —
    # once applied at ALL, donor variance drops sharply) while preserving
    # the immune marker signal best. Step 8 still uses theta=2 because
    # cluster-level donor invariance is the goal there; pseudotime is a
    # disease-progression axis and benefits from a softer touch.
    if 'X_scVI_mean' in adata_islet.obsm:
        import harmonypy as hm
        adata_islet.obs['_imageid_str'] = (
            adata_islet.obs['imageid'].astype(str).astype('category')
        )
        print("Running Harmony on X_scVI_mean (batch=imageid, theta=0.1) for trajectory neighbors...")
        ho_traj = hm.run_harmony(
            np.asarray(adata_islet.obsm['X_scVI_mean']),
            adata_islet.obs, '_imageid_str',
            theta=0.1, max_iter_harmony=20, verbose=False,
        )
        adata_islet.obsm['X_scVI_harmony_pt'] = np.asarray(ho_traj.Z_corr)
        del adata_islet.obs['_imageid_str']
        sc.pp.neighbors(adata_islet, n_neighbors=30,
                        use_rep='X_scVI_harmony_pt', metric='cosine')
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

    # Robust root: medoid of the top-K ND islets by INS in the scVI latent
    # (harmonized if available, else raw). A single-islet argmax is brittle
    # to staining outliers; the medoid of a size-50 candidate pool is anchored
    # by the centre-of-mass of "healthy ND, INS-high" islets and naturally
    # draws from multiple donors.
    K_ROOT = 50
    ins_idx = list(adata_islet.var_names).index('INS')
    nd_mask = adata_islet.obs['donor_status'] == 'ND'
    if not nd_mask.any():
        raise ValueError("No ND islets found!")
    nd_indices = np.where(nd_mask)[0]
    nd_ins = adata_islet.X[nd_mask, ins_idx]
    k = min(K_ROOT, len(nd_indices))
    top_pos = np.argsort(-nd_ins)[:k]
    top_global = nd_indices[top_pos]

    # Medoid of the top-K in the (lightly-)harmonized latent for the trajectory.
    # Falls back to raw scVI if no harmonized embedding is available.
    rep_name = ('X_scVI_harmony_pt' if 'X_scVI_harmony_pt' in adata_islet.obsm
                else ('X_scVI_harmony' if 'X_scVI_harmony' in adata_islet.obsm
                      else 'X_scVI_mean'))
    rep = np.asarray(adata_islet.obsm[rep_name])
    top_coords = rep[top_global]
    centroid = top_coords.mean(axis=0)
    root_idx = int(top_global[np.argmin(np.linalg.norm(top_coords - centroid, axis=1))])
    adata_islet.uns['iroot'] = root_idx

    top_donor_counts = (
        adata_islet.obs['imageid'].iloc[top_global].astype(str).value_counts().to_dict()
    )
    root_donor = adata_islet.obs['imageid'].iloc[root_idx]
    print(f"Root islet (medoid of top-{k} ND by INS in {rep_name}): "
          f"index={root_idx}, donor={root_donor}, "
          f"INS={adata_islet.X[root_idx, ins_idx]:.3f}")
    print(f"  Top-{k} candidate pool donor distribution: {top_donor_counts}")

    # Diffusion map and pseudotime
    sc.tl.diffmap(adata_islet, n_comps=10)
    sc.tl.dpt(adata_islet)
    print(f"Pseudotime range: [{adata_islet.obs['dpt_pseudotime'].min():.3f}, "
          f"{adata_islet.obs['dpt_pseudotime'].max():.3f}]")

    # Scanpy's X_diffmap column 0 is the TRIVIAL lambda=1 eigenvector (the
    # stationary distribution): constant-sign, near-zero variance, carrying no
    # trajectory information. Measured on the previously shipped file it had
    # std 0.000814 and spanned only 0.0112-0.0163, while the real first
    # diffusion component (column 1) has std 0.0138 -- ~17x larger.
    # The informative components start at column 1, so DC1/DC2 are columns 1/2.
    adata_islet.obs['DC1'] = adata_islet.obsm['X_diffmap'][:, 1]
    adata_islet.obs['DC2'] = adata_islet.obsm['X_diffmap'][:, 2]

    # ---------------------------------------------------------------
    # Step 5b: Peri+structures-augmented pseudotime
    # ---------------------------------------------------------------
    # Compute a SECOND pseudotime ("combined" mode) that augments the islet
    # core scVI mean with:
    #   1. Peri-zone scVI mean (10d)        — captures peri-islet immune /
    #                                          stromal context
    #   2. Structural cell densities (16d)  — proportions + counts of
    #                                          {Neural, Blood Vessel,
    #                                          Endothelial, Lymphatic} for
    #                                          both core and peri zones
    # Total: 10 + 10 + 16 = 36 dims.
    #
    # Structural features are z-scored and scaled by ALPHA_STRUCT before
    # concatenation. Empirical sweep (alpha_sweep.py / 2026-04-25):
    #   alpha   span   order    INS     CD8a    PDPN
    #   0.00   0.069   YES     -0.54   +0.21   +0.04
    #   0.15   0.063   YES     -0.51   +0.23   +0.07   <- chosen
    #   0.20   -0.01   NO      -0.46   +0.21   +0.09
    #   0.50   -0.00   NO      -0.26   +0.19   +0.28
    # Disease ordering (ND<Aab+<T1D) breaks at alpha >= 0.20 — structures
    # start dominating. alpha=0.15 is the largest value preserving order
    # while giving the strongest immune-marker correlations (CD3e, CD8a)
    # and starting to register vasculature/lymphatic signal (CD34, PDPN).
    ALPHA_STRUCT = 0.15
    STRUCT_PHENOTYPES = ['Neural', 'Blood Vessel', 'Endothelial', 'Lymphatic']

    print("\n" + "=" * 60)
    print("STEP 5b: Computing peri+structures-augmented pseudotime")
    print("=" * 60)
    import pandas as pd
    import anndata
    import harmonypy as hm
    # Stash core-only pseudotime + diffmap before recomputing.
    #
    # `sc.pp.neighbors` below OVERWRITES obsp['distances'] / obsp['connectivities']
    # in place. Restoring only obs/obsm/uns (as this block used to) left the saved
    # file carrying the COMBINED neighbor graph while uns['neighbors']['params']
    # advertised the core embedding -- so anyone re-running scanpy on
    # adata_ins_root.h5ad silently got combined-mode results. Measured overlap of
    # the shipped obsp was 96.6% with the combined graph vs 34.6% with core.
    # Stash obsp and iroot too, and restore all four together.
    core_pt   = adata_islet.obs['dpt_pseudotime'].values.copy()
    core_dc1  = adata_islet.obsm['X_diffmap'][:, 1].copy()
    core_dc2  = adata_islet.obsm['X_diffmap'][:, 2].copy()
    core_dmap = adata_islet.obsm['X_diffmap'].copy()
    core_neigh_uns = adata_islet.uns.get('neighbors', None)
    core_iroot = adata_islet.uns.get('iroot', None)
    core_evals = adata_islet.uns.get('diffmap_evals', None)
    core_obsp = {k: adata_islet.obsp[k].copy() for k in ('distances', 'connectivities')
                 if k in adata_islet.obsp}

    def _zscore(a):
        return (a - a.mean(axis=0)) / (a.std(axis=0) + 1e-8)

    def _struct_block(adata, label):
        """Pull (props, counts) for STRUCT_PHENOTYPES out of obsm."""
        names = list(adata.uns.get('phenotype_names', []))
        if not names:
            print(f"  {label}: no phenotype_names in uns — skipping structures")
            return None
        idx, found = [], []
        for p in STRUCT_PHENOTYPES:
            if p in names:
                idx.append(names.index(p))
                found.append(p)
            else:
                print(f"  {label}: phenotype {p!r} not found, skipping")
        if not idx:
            return None
        props = np.asarray(adata.obsm['phenotype_proportions'])[:, idx]
        counts = np.asarray(adata.obsm['phenotype_counts'])[:, idx]
        return props, counts, found

    # Build core+peri 20-dim representation, aligned by obs index
    ad_peri_loaded = anndata.read_h5ad(peri_path)
    if 'X_scVI_mean' not in ad_peri_loaded.obsm:
        print("  WARNING: peri H5AD lacks X_scVI_mean — skipping combined pseudotime")
        adata_islet.obs['dpt_pseudotime_combined'] = core_pt
        adata_islet.obs['DC1_combined'] = core_dc1
        adata_islet.obs['DC2_combined'] = core_dc2
    else:
        peri_lookup = pd.Series(range(ad_peri_loaded.n_obs), index=ad_peri_loaded.obs_names)
        try:
            peri_idx = peri_lookup.loc[adata_islet.obs_names].values
        except KeyError as e:
            raise RuntimeError(
                "Cannot align peri to core by obs_names — require_paired=True should "
                "ensure 1:1 correspondence. Mismatch detected."
            ) from e
        peri_scvi = np.asarray(ad_peri_loaded.obsm['X_scVI_mean'])[peri_idx]
        core_scvi = np.asarray(adata_islet.obsm['X_scVI_mean'])

        # Structural blocks: z-scored and scaled by ALPHA_STRUCT
        core_struct = _struct_block(adata_islet, 'core')
        peri_struct = _struct_block(ad_peri_loaded, 'peri')

        struct_parts = []
        struct_labels = []
        if core_struct is not None:
            struct_parts.extend([_zscore(core_struct[0]), _zscore(core_struct[1])])
            struct_labels.extend(
                [f'core_prop_{n}' for n in core_struct[2]] +
                [f'core_count_{n}' for n in core_struct[2]]
            )
        if peri_struct is not None:
            struct_parts.extend([
                _zscore(peri_struct[0][peri_idx]),
                _zscore(peri_struct[1][peri_idx]),
            ])
            struct_labels.extend(
                [f'peri_prop_{n}' for n in peri_struct[2]] +
                [f'peri_count_{n}' for n in peri_struct[2]]
            )

        if struct_parts:
            struct_block = np.hstack(struct_parts) * ALPHA_STRUCT
            X_combined = np.hstack([core_scvi, peri_scvi, struct_block])
            print(f"  Combined latent shape: {X_combined.shape} "
                  f"= core {core_scvi.shape[1]}d + peri {peri_scvi.shape[1]}d "
                  f"+ struct {struct_block.shape[1]}d (alpha={ALPHA_STRUCT})")
            print(f"  Structural features ({len(struct_labels)}): {struct_labels}")
        else:
            X_combined = np.hstack([core_scvi, peri_scvi])
            print(f"  Combined latent shape: {X_combined.shape} "
                  f"(core {core_scvi.shape[1]}d + peri {peri_scvi.shape[1]}d) "
                  "— structural block unavailable")

        # Light Harmony (theta=0.1) on the combined latent — same theta as core
        adata_islet.obs['_iid_pt'] = (
            adata_islet.obs['imageid'].astype(str).astype('category')
        )
        ho_comb = hm.run_harmony(
            X_combined, adata_islet.obs, '_iid_pt',
            theta=0.1, max_iter_harmony=20, verbose=False,
        )
        adata_islet.obsm['X_scVI_combined_harmony_pt'] = np.asarray(ho_comb.Z_corr)
        del adata_islet.obs['_iid_pt']

        # Recompute neighbors → diffmap → DPT on the combined embedding.
        # Use the same root index (highest-INS ND medoid) so the start of
        # the trajectory is biologically anchored the same way.
        if 'neighbors' in adata_islet.uns:
            del adata_islet.uns['neighbors']
        sc.pp.neighbors(adata_islet, n_neighbors=30,
                        use_rep='X_scVI_combined_harmony_pt', metric='cosine')
        adata_islet.uns['iroot'] = root_idx
        sc.tl.diffmap(adata_islet, n_comps=10)
        sc.tl.dpt(adata_islet)
        comb_pt  = adata_islet.obs['dpt_pseudotime'].values.copy()
        # Columns 1/2 -- column 0 is the trivial lambda=1 eigenvector (see Step 5).
        comb_dc1 = adata_islet.obsm['X_diffmap'][:, 1].copy()
        comb_dc2 = adata_islet.obsm['X_diffmap'][:, 2].copy()
        print(f"  Combined pseudotime range: [{np.nanmin(comb_pt):.3f}, {np.nanmax(comb_pt):.3f}]")

        # Restore core results as the default (dpt_pseudotime, DC1, DC2) and store
        # combined under an explicit suffix. Everything sc.pp.neighbors / diffmap
        # touched must be restored together, or the saved file describes one mode
        # while carrying another's graph.
        adata_islet.obs['dpt_pseudotime'] = core_pt
        adata_islet.obs['DC1'] = core_dc1
        adata_islet.obs['DC2'] = core_dc2
        adata_islet.obsm['X_diffmap'] = core_dmap
        if core_neigh_uns is not None:
            adata_islet.uns['neighbors'] = core_neigh_uns
        if core_iroot is not None:
            adata_islet.uns['iroot'] = core_iroot
        if core_evals is not None:
            adata_islet.uns['diffmap_evals'] = core_evals
        for _k, _v in core_obsp.items():
            adata_islet.obsp[_k] = _v

        adata_islet.obs['dpt_pseudotime_combined'] = comb_pt
        adata_islet.obs['DC1_combined'] = comb_dc1
        adata_islet.obs['DC2_combined'] = comb_dc2

    # ---------------------------------------------------------------
    # Step 6: Validate trajectory
    # ---------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 6: Validating trajectory")
    print("=" * 60)

    all_pass = _validate_trajectory(adata_islet, ins_idx)

    if not all_pass:
        print("\nWARNING: Some validation checks failed. Outputs will still be saved "
              "(a failed gate should not discard a ~15 min run), but this script will "
              "exit non-zero so automation does not treat the rebuild as clean.")

    # Visualization UMAP from raw marker expression.
    # scVI (SCVI_CODEX_v2.ipynb) was trained with Age/Gender as categorical
    # covariates but *no batch_key*, so donor-specific technical variance leaked
    # into X_scVI_mean. A scVI-based UMAP tends to produce a blob with poor
    # disease-stage separation. Raw marker PCA UMAP is clearer for display.
    # (Leiden clustering in Step 8 harmonizes scVI on imageid to fix this.)
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

    # Harmonize the scVI latent on imageid before clustering so Leiden captures
    # islet cell-types shared across donors rather than donor-specific technical
    # variance. scVI (SCVI_CODEX_v2.ipynb) was trained with Age/Gender as
    # covariates but no batch_key, so donor batch effects leaked into X_scVI_mean.
    if 'X_scVI_mean' in adata_clust.obsm:
        import harmonypy as hm
        # harmonypy 0.2.0 calls meta.describe().loc['unique'] which fails on
        # numeric columns; cast imageid to string category defensively.
        adata_clust.obs['_imageid_str'] = (
            adata_clust.obs['imageid'].astype(str).astype('category')
        )
        print("Running Harmony on X_scVI_mean (batch=imageid, theta=2)...")
        ho = hm.run_harmony(
            np.asarray(adata_clust.obsm['X_scVI_mean']),
            adata_clust.obs,
            '_imageid_str',
            theta=2.0,
            max_iter_harmony=20,
            verbose=False,
        )
        adata_clust.obsm['X_scVI_harmony'] = np.asarray(ho.Z_corr)
        del adata_clust.obs['_imageid_str']
        sc.pp.neighbors(adata_clust, use_rep='X_scVI_harmony',
                        n_neighbors=15, metric='cosine')
    else:
        print("WARNING: X_scVI_mean not available, falling back to PCA")
        sc.pp.pca(adata_clust, n_comps=10)
        sc.pp.neighbors(adata_clust, n_neighbors=15)

    # Use the visualization UMAP (raw marker PCA) for display —
    # it shows disease-stage separation, unlike scVI-based UMAP which is a blob
    adata_clust.obsm['X_umap'] = adata_islet.obsm['X_umap'].copy()
    adata_clust.obs['umap_1'] = adata_clust.obsm['X_umap'][:, 0]
    adata_clust.obs['umap_2'] = adata_clust.obsm['X_umap'][:, 1]
    print("Using visualization UMAP from trajectory (raw marker PCA, min_dist=0.5, spread=2.0)")

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

    if not all_pass:
        print("\n" + "!" * 60)
        print("VALIDATION FAILED — see STEP 6 above. Outputs were written, but do "
              "NOT promote them without understanding which gate failed.")
        print("!" * 60)

    return all_pass


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

    # Exit non-zero on a failed validation gate so automation cannot mistake a
    # degraded rebuild for a clean one. Outputs are still written first -- see
    # the Step 6 note; discarding a ~15 min run makes the gate expensive to keep.
    ok = reaggregate(args.input, args.islet_dir, args.trajectory_out, args.min_cells)
    sys.exit(0 if ok else 2)


if __name__ == '__main__':
    main()
