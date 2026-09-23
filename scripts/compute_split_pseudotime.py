#!/usr/bin/env python
"""
Recompute pseudotime SEPARATELY for small (<= CELL_SPLIT cells) and large
(> CELL_SPLIT cells) islets, and write both back into data/adata_ins_root.h5ad.

WHY THIS EXISTS
---------------
The shipped `dpt_pseudotime` is computed over all 5,214 islets at once. 60.8% of
those are aggregates of <= 10 cells, and their mean scVI latent is correspondingly
noisy — an islet can sit far out in latent space purely from sampling noise while
having ordinary marker expression. Measured on the shipped file:

    spearman(dpt_pseudotime, total_cells) = -0.274
    islets with dpt > 0.6:  median 7 cells   (vs 10-11 in the bulk)
    min_cells >= 25 collapses the tail: frac(dpt > 0.6) 3.97% -> 0.44%

Because scanpy normalises DPT by dividing by the maximum reachable distance,
those few noisy islets set the top of the axis for everyone: 96.0% of islets end
up below 0.6, and the bulk (to p95) occupies only 57.4% of the drawn axis.

Splitting by cell count and recomputing within each stratum lets the two
populations be read on their own terms instead of one being compressed by the
other's noise.

>>> THE TWO PSEUDOTIMES ARE NOT ON A COMMON SCALE. <<<
Each is an independent diffusion map with its own root, its own Harmony fit and
its own divide-by-max normalisation. `small = 0.4` and `large = 0.4` do NOT mean
the same thing. Compare SHAPES and ORDERINGS within a stratum; never compare the
numbers across strata, and never merge the two columns into one axis.

METHOD — identical to reaggregate_islets.py Step 5, applied within each stratum:
  1. Harmony on obsm['X_scVI_mean'] (batch = imageid, theta = 0.1)
  2. sc.pp.neighbors(n_neighbors=30, metric='cosine') on that embedding
  3. root = medoid of the top-50 ND islets by INS, in the harmonised embedding
  4. sc.tl.diffmap(n_comps=10) -> sc.tl.dpt()

Outputs (added to data/adata_ins_root.h5ad, NaN outside the stratum):
  obs['dpt_pseudotime_small']   islets with total_cells <= CELL_SPLIT
  obs['dpt_pseudotime_large']   islets with total_cells >  CELL_SPLIT
  obs['cell_count_stratum']     'small' / 'large'
  uns['split_pseudotime']       parameters + per-stratum validation stats

Usage:
    conda activate scvi-env
    python scripts/compute_split_pseudotime.py [--split 10] [--h5ad PATH]
"""

import os
import sys
import argparse
import numpy as np
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.join(SCRIPT_DIR, '..')
DEFAULT_H5AD = os.path.join(PROJECT_ROOT, 'data', 'adata_ins_root.h5ad')

# Matches reaggregate_islets.py Step 5. Kept as module constants so a future
# change there is a visible mismatch here rather than a silent divergence.
N_NEIGHBORS = 30
METRIC = 'cosine'
HARMONY_THETA = 0.1
N_DIFFMAP_COMPS = 10
K_ROOT = 50


def _eta2_donor_within_status(pt, iid, status_col, stt, min_islets=3):
    """Fraction of within-status pseudotime variance explained by donor identity."""
    m = status_col == stt
    p, d = pt[m], iid[m]
    if p.size < min_islets:
        return np.nan
    grand = p.mean()
    ss_tot = ((p - grand) ** 2).sum()
    if ss_tot <= 0:
        return np.nan
    ss_bet = 0.0
    for dv in np.unique(d):
        g = p[d == dv]
        if g.size:
            ss_bet += g.size * (g.mean() - grand) ** 2
    return float(ss_bet / ss_tot)


def compute_stratum(adata_full, mask, label, sc, hm):
    """Recompute DPT within one cell-count stratum. Returns (pt, stats)."""
    from scipy import stats as sps

    n = int(mask.sum())
    print(f"\n{'=' * 64}\nSTRATUM '{label}': {n} islets\n{'=' * 64}")
    if n < 50:
        print(f"  SKIP: only {n} islets, too few for a diffusion map")
        return None, {'error': f'only {n} islets'}

    ad_s = adata_full[mask].copy()

    # Drop everything carried over from the all-islet run. Leaving any of these
    # in place lets scanpy silently reuse a graph built on a different island of
    # the data -- the exact class of bug documented in docs/trajectory.md, where
    # a saved obsp belonged to a different embedding than uns claimed.
    for key in ('X_umap', 'X_diffmap', 'X_pca', 'X_scVI_harmony_pt',
                'X_scVI_combined_harmony_pt'):
        if key in ad_s.obsm:
            del ad_s.obsm[key]
    for key in ('neighbors', 'paga', 'iroot', 'diffmap_evals'):
        if key in ad_s.uns:
            del ad_s.uns[key]
    ad_s.obsp = {}

    if 'X_scVI_mean' not in ad_s.obsm:
        raise RuntimeError("X_scVI_mean missing; cannot recompute pseudotime")

    # 1. Harmony (theta=0.1), batch = donor
    ad_s.obs['_imageid_str'] = ad_s.obs['imageid'].astype(str).astype('category')
    print(f"  Harmony on X_scVI_mean (batch=imageid, theta={HARMONY_THETA})...")
    ho = hm.run_harmony(
        np.asarray(ad_s.obsm['X_scVI_mean']),
        ad_s.obs, ['_imageid_str'],
        theta=HARMONY_THETA, max_iter_harmony=20, verbose=False,
    )
    ad_s.obsm['X_scVI_harmony_pt'] = np.asarray(ho.Z_corr)

    # 2. Neighbors
    sc.pp.neighbors(ad_s, n_neighbors=N_NEIGHBORS,
                    use_rep='X_scVI_harmony_pt', metric=METRIC)

    # 3. Root: medoid of top-K ND islets by INS, inside THIS stratum. Recomputing
    #    the root per stratum is required, not optional -- the all-islet root may
    #    not even be a member of the stratum, and a root outside the graph makes
    #    DPT undefined.
    ins_idx = list(ad_s.var_names).index('INS')
    nd_mask = (ad_s.obs['donor_status'] == 'ND').to_numpy()
    if not nd_mask.any():
        raise RuntimeError(f"stratum '{label}' has no ND islets; no root available")
    nd_indices = np.where(nd_mask)[0]
    X = ad_s.X
    X = X.toarray() if hasattr(X, 'toarray') else np.asarray(X)
    nd_ins = X[nd_mask, ins_idx]
    k = min(K_ROOT, len(nd_indices))
    top_global = nd_indices[np.argsort(-nd_ins)[:k]]
    rep = np.asarray(ad_s.obsm['X_scVI_harmony_pt'])
    centroid = rep[top_global].mean(axis=0)
    root_idx = int(top_global[np.argmin(np.linalg.norm(rep[top_global] - centroid, axis=1))])
    ad_s.uns['iroot'] = root_idx
    root_donor = str(ad_s.obs['imageid'].iloc[root_idx])
    print(f"  Root: idx={root_idx} donor={root_donor} INS={X[root_idx, ins_idx]:.3f} "
          f"(medoid of top-{k} ND by INS)")

    # 4. Diffusion map -> DPT
    sc.tl.diffmap(ad_s, n_comps=N_DIFFMAP_COMPS)
    sc.tl.dpt(ad_s)
    pt = ad_s.obs['dpt_pseudotime'].to_numpy(dtype=float)
    print(f"  Pseudotime range: [{np.nanmin(pt):.3f}, {np.nanmax(pt):.3f}]")

    # 5. Validation, same gates the all-islet run reports
    st = {'n_islets': n, 'root_index': root_idx, 'root_donor': root_donor,
          'n_neighbors': N_NEIGHBORS, 'metric': METRIC,
          'harmony_theta': HARMONY_THETA, 'k_root': k}
    finite = np.isfinite(pt)
    for marker in ('INS', 'GCG', 'SST'):
        if marker in ad_s.var_names:
            mi = list(ad_s.var_names).index(marker)
            r = sps.spearmanr(pt[finite], X[finite, mi]).statistic
            st[f'spearman_{marker}'] = float(r)
            print(f"  Spearman(pt, {marker}) = {r:+.3f}")

    status = ad_s.obs['donor_status'].astype(str).to_numpy()
    iid = ad_s.obs['imageid'].astype(str).to_numpy()
    for stt in ('ND', 'Aab+', 'T1D'):
        m = status == stt
        if m.any():
            st[f'islet_mean_{stt}'] = float(np.nanmean(pt[m]))
            st[f'donor_eta2_{stt}'] = _eta2_donor_within_status(pt, iid, status, stt)
    donor_means = {}
    for stt in ('ND', 'Aab+', 'T1D'):
        m = status == stt
        if m.any():
            dm = [np.nanmean(pt[(iid == d) & m]) for d in np.unique(iid[m])]
            donor_means[stt] = float(np.nanmean(dm))
    st['donor_mean_by_status'] = donor_means
    print(f"  Donor-mean pseudotime by status: "
          f"{ {k2: round(v, 4) for k2, v in donor_means.items()} }")
    st['monotone_ND_Aab_T1D'] = bool(
        len(donor_means) == 3
        and donor_means['ND'] < donor_means['Aab+'] < donor_means['T1D'])
    print(f"  Monotone ND < Aab+ < T1D: {st['monotone_ND_Aab_T1D']}")

    # Axis usage -- the reason this script exists.
    st['frac_above_0.6'] = float((pt[finite] > 0.6).mean())
    st['p95'] = float(np.percentile(pt[finite], 95))
    st['iqr_width'] = float(np.percentile(pt[finite], 75) - np.percentile(pt[finite], 25))
    print(f"  Axis usage: frac>0.6 = {st['frac_above_0.6']:.4f}, "
          f"p95 = {st['p95']:.3f}, IQR width = {st['iqr_width']:.3f}")
    return pt, st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--h5ad', default=DEFAULT_H5AD)
    ap.add_argument('--split', type=int, default=10,
                    help='cells per islet; <= split is "small", > split is "large"')
    args = ap.parse_args()

    import scanpy as sc
    import harmonypy as hm
    import anndata as ad
    sc.settings.verbosity = 1

    print(f"Reading {args.h5ad}")
    adata = ad.read_h5ad(args.h5ad)
    print(f"  {adata.n_obs} islets x {adata.n_vars} markers")

    tc = adata.obs['total_cells'].to_numpy(dtype=float)
    small = tc <= args.split
    large = tc > args.split
    print(f"  split at {args.split}: small={int(small.sum())} large={int(large.sum())}")

    pt_small, st_small = compute_stratum(adata, small, 'small', sc, hm)
    pt_large, st_large = compute_stratum(adata, large, 'large', sc, hm)

    out_small = np.full(adata.n_obs, np.nan)
    out_large = np.full(adata.n_obs, np.nan)
    if pt_small is not None:
        out_small[small] = pt_small
    if pt_large is not None:
        out_large[large] = pt_large

    adata.obs['dpt_pseudotime_small'] = out_small
    adata.obs['dpt_pseudotime_large'] = out_large
    adata.obs['cell_count_stratum'] = np.where(small, 'small', 'large')
    adata.uns['split_pseudotime'] = {
        'split_at': args.split,
        'note': ('Independent diffusion maps per stratum: separate Harmony fit, '
                 'separate root, separate divide-by-max normalisation. The two '
                 'columns are NOT on a common scale and must never be merged '
                 'into one axis or compared as absolute numbers.'),
        'embedding': 'Harmony(X_scVI_mean, batch=imageid, theta=0.1) -> X_scVI_harmony_pt',
        'small': st_small, 'large': st_large,
    }

    adata.write_h5ad(args.h5ad)
    print(f"\nWrote {args.h5ad}")
    print("  obs['dpt_pseudotime_small'] : "
          f"{int(np.isfinite(out_small).sum())} finite")
    print("  obs['dpt_pseudotime_large'] : "
          f"{int(np.isfinite(out_large).sum())} finite")


if __name__ == '__main__':
    main()
