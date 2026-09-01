"""
PAGA pseudotime + force-directed layout for the Islet Explorer trajectory tab.

Ports Barlow et al., eLife 2024;13:RP100535, Methods -> "PAGA analysis":

    "The PAGA embedding was computed using the default parameters except for the
     following: The neighborhood search was performed using cosine distance and
     15 nearest neighbors; Leiden clustering used a resolution of 1."

and their Figure 2B, "PAGA-force directed layout of islets colored by pseudotime.
Each point represents an islet."

The force-directed layout is the point of this module. A feature-vs-pseudotime
scatter piles ~85-96% of islets onto already-inked pixels because most cell-type
proportions are exactly zero (95.6% for CD8a T cells, 98.2% for B cells) and no
axis transform moves an exact zero. The PAGA layout is a 2-D graph embedding, so
the layout algorithm spreads islets apart by construction: measured on the same
occupancy grid, 8.1% overplotted vs 84-96% for the scatter, densest pixel 3
islets vs 125-129.

Called from R via reticulate (app/shiny_app/R/paga_layout.R). Kept as a real
module rather than inline R strings so it can be run and tested standalone:

    python app/shiny_app/python/paga_pseudotime.py --min-cells 10
"""

import os
import numpy as np
import pandas as pd

# Barlow's stated parameters.
N_NEIGHBORS = 15
METRIC = "cosine"
LEIDEN_RESOLUTION = 1.0
N_DIFFMAP_COMPS = 10

_FEATURES_CSV = None


def _features_path():
    """Resolve data/barlow_islet_features.csv from either the app dir or repo root."""
    global _FEATURES_CSV
    if _FEATURES_CSV and os.path.exists(_FEATURES_CSV):
        return _FEATURES_CSV
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (
        os.path.join(os.environ.get("ISLET_DATA_ROOT", ""), "barlow_islet_features.csv"),
        os.path.join(here, "..", "..", "..", "data", "barlow_islet_features.csv"),
        os.path.join(here, "..", "..", "data", "barlow_islet_features.csv"),
        os.path.join("data", "barlow_islet_features.csv"),
    ):
        if os.path.exists(cand):
            _FEATURES_CSV = os.path.abspath(cand)
            return _FEATURES_CSV
    return None


def compute(min_cells=1, max_cells=None):
    """
    Returns a dict of equal-length lists (islet_id, case_id, donor_status, leiden,
    fa_1, fa_2, pseudotime) plus scalar diagnostics, or {'error': msg}.

    R/reticulate converts this to a data.frame; returning plain lists rather than
    a DataFrame keeps the conversion predictable across reticulate versions.
    """
    path = _features_path()
    if path is None:
        return {"error": "data/barlow_islet_features.csv not found "
                         "(run scripts/build_barlow_features.py)"}
    try:
        import scanpy as sc
        import anndata as ad
    except ImportError as e:
        return {"error": f"scanpy/anndata unavailable in this interpreter: {e}"}

    try:
        mc = max(1, int(float(min_cells)))
    except (TypeError, ValueError):
        mc = 1
    try:
        xc = None if max_cells is None else int(float(max_cells))
    except (TypeError, ValueError):
        xc = None

    b = pd.read_csv(path)
    n = b["total_cells_core"]
    sel = n >= mc
    if xc is not None:
        sel &= n <= xc
    sub = b[sel].reset_index(drop=True)
    if len(sub) < 50:
        return {"error": f"only {len(sub)} islets in cell range "
                         f"[{mc}, {xc if xc is not None else 'inf'}]; "
                         "a diffusion map needs more"}

    feat = [c for c in sub.columns if c.startswith(("core_", "peri_"))]
    X = sub[feat].to_numpy(dtype=float)
    # Constant columns carry no information and make the cosine metric unstable.
    keep = X.std(axis=0) > 0
    X = X[:, keep]
    if X.shape[1] < 3:
        return {"error": "fewer than 3 informative composition features"}

    obs = sub[["islet_id", "case_id", "donor_status"]].astype(str).copy()
    # AnnData requires string observation names. Supplying them explicitly
    # avoids its ImplicitModificationWarning without changing row identity.
    obs.index = pd.Index(sub.index.astype(str), name="obs_id")
    A = ad.AnnData(X=X, obs=obs)
    sc.pp.neighbors(A, n_neighbors=N_NEIGHBORS, metric=METRIC, use_rep="X")
    # flavor="igraph" is scanpy's supported non-deprecated path and is what the
    # future default will be; n_iterations=2 matches that default.
    sc.tl.leiden(A, resolution=LEIDEN_RESOLUTION, key_added="leiden",
                 flavor="igraph", n_iterations=2, directed=False, random_state=0)
    sc.tl.paga(A, groups="leiden")
    # PAGA connectivities must be computed before they can seed the layout.
    sc.pl.paga(A, plot=False)
    sc.tl.draw_graph(A, init_pos="paga", random_state=0)
    lay_key = next(k for k in A.obsm if k.startswith("X_draw_graph"))
    xy = np.asarray(A.obsm[lay_key])

    # ---- Root: "pseudotemporal distance from the centroid of non-T1D islets" ----
    # DPT needs an actual vertex, so use the ND islet nearest the ND centroid --
    # the medoid of the healthy cloud. Candidates are restricted to ND because
    # the globally nearest islet to that centroid can be Aab+, which anchors the
    # axis in the wrong disease stage.
    sc.tl.diffmap(A, n_comps=N_DIFFMAP_COMPS)
    dm = np.asarray(A.obsm["X_diffmap"])[:, 1:4]   # skip the trivial 0th eigenvector
    nd = np.where(A.obs["donor_status"].to_numpy() == "ND")[0]
    if nd.size == 0:
        return {"error": "no ND islets above this threshold; no root available"}
    cen = dm[nd].mean(axis=0)
    root = int(nd[np.argmin(np.linalg.norm(dm[nd] - cen, axis=1))])
    A.uns["iroot"] = root
    sc.tl.dpt(A)
    pt = A.obs["dpt_pseudotime"].to_numpy(dtype=float)

    # Orient so beta-cell content DECREASES with pseudotime. The diffusion map's
    # sign is arbitrary, so without this the axis can flip between thresholds and
    # the same islet appears at opposite ends on consecutive settings.
    beta = sub["core_Beta cell"].to_numpy(dtype=float)
    fin = np.isfinite(pt) & np.isfinite(beta)
    flipped = False
    if fin.sum() > 10:
        from scipy.stats import spearmanr
        if spearmanr(pt[fin], beta[fin]).statistic > 0:
            pt = np.nanmax(pt[np.isfinite(pt)]) - pt
            flipped = True

    return {
        "islet_id": sub["islet_id"].astype(str).tolist(),
        "case_id": sub["case_id"].astype(str).tolist(),
        "donor_status": sub["donor_status"].astype(str).tolist(),
        "leiden": A.obs["leiden"].astype(str).tolist(),
        "fa_1": xy[:, 0].astype(float).tolist(),
        "fa_2": xy[:, 1].astype(float).tolist(),
        "pseudotime": pt.astype(float).tolist(),
        "n_islets": int(A.n_obs),
        "n_clusters": int(A.obs["leiden"].nunique()),
        "layout_key": str(lay_key),
        "min_cells": mc,
        "max_cells": xc,
        "root_index": root,
        "root_donor": str(A.obs["case_id"].iloc[root]),
        "flipped": bool(flipped),
    }


if __name__ == "__main__":
    import argparse, json
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-cells", type=int, default=10)
    ap.add_argument("--max-cells", type=int, default=None)
    a = ap.parse_args()
    r = compute(a.min_cells, a.max_cells)
    if "error" in r:
        print("ERROR:", r["error"])
    else:
        print(json.dumps({k: v for k, v in r.items()
                          if not isinstance(v, list)}, indent=2))
