#!/usr/bin/env python
"""Derive per-image marker gates calibrated against independent ground truth.

WHY THIS EXISTS
---------------
The shipped phenotyping used gates drawn manually in napari on ONE image
(6450), broadcast to all 15 donors: `uns['gates']` is 31x15 but 0 of 31 markers
differ across columns, and `napariGaterProvenance['manually_adjusted']` holds
only '6450'. Image 6450 is also compositionally atypical (5.6% INS+ but 33.8%
GCG+), so gates drawn there misfire elsewhere.

Two automated alternatives were tested and REJECTED:
  * Per-image GMM (scimap's `rescale(method='by_image')` rule: midpoint between
    the top two component means). On rare markers the mixture spends all its
    components on the negative population and the gate lands inside it --
    CD20 went 0.08% -> 26.8% positive, HLADR 1.27% -> 45.8%, CD3e 0.88% ->
    35.5%. 6 of 13 markers tested failed this way. Do not re-propose it.
  * Negative-mode alignment. Helps abundant structural markers (BCatenin
    spread 4.5x -> 1.3x) but does nothing for rare markers (CD3e 14.4x ->
    14.5x) and made HLADR worse (43.6x -> 68.0x).

This script instead calibrates against two INDEPENDENT ground-truth sources
that already exist in the project:

  TIER 1a -- groovy endocrine gating (`data/annotations.tsv`).
    QuPath groovy scripts classified every cell inside each islet annotation
    into mutually exclusive Num {INS,GCG,SST,Multi}_pos / Endo_Neg counts
    (verified: these partition Num Detections exactly, 5382/5382 rows).
    Gates for INS/GCG/SST are set per image so the positive-cell count matches.

  TIER 1b -- pixel-classifier structures (`data/annotations.tsv` +
    signed `Dist to Closest {Capillary,Nerve,Lymphatic}` in the H5AD obs).
    Cells with signed distance < 0 lie inside the structure. Those labels
    separate their marker essentially perfectly (AUC: CD31/Capillary 0.990,
    B3TUBB/Nerve 0.991, PDPN/Lymphatic 0.951), so the gate is set per image at
    the Youden-J optimum.

  TIER 2 -- lineage markers with no ground truth available. Manual 6450 gate
    retained and FLAGGED in the output, so consumers know these are unvalidated.

  TIER 3 -- markers not used in broad lineage assignment. Deferred; manual
    gate retained. Revisit only if a downstream rule starts using them.

VALIDATION ALREADY PERFORMED (per-islet agreement with groovy; gates are fitted
on per-IMAGE totals, so per-ISLET agreement is partially held out):

    marker  scheme       spearman   mean|err|
    INS     manual          0.860       2.08
    INS     calibrated      0.895       0.71
    GCG     manual          0.916       1.91
    GCG     calibrated      0.950       0.67
    SST     manual          0.652       2.30
    SST     calibrated      0.958       0.57

Usage:
    python scripts/derive_gates.py [--n-sample 300000] [--dry-run]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
H5AD = REPO / "single_cell_analysis" / "CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad"
ANNOTATIONS = REPO / "data" / "annotations.tsv"
# QuPath groovy exports live OUTSIDE the repo. Note build_master_excel.py still
# defaults to a stale 'data/results' path that does not exist.
GROOVY_DIR = Path("/home/smith6jt/panc_CODEX/results/groovy_exports")
OUT_CSV = REPO / "data" / "phenotype_gates_per_image.csv"
OUT_PROV = REPO / "data" / "phenotype_gates_provenance.csv"

RANDOM_STATE = 0

# Markers calibrated against groovy endocrine counts (per image).
GROOVY_MARKERS = ["INS", "GCG", "SST"]

# Markers calibrated against pixel-classifier structures (per image).
# structure name in annotations.tsv -> marker
#
# CAVEAT on PDPN: inside-vs-far intensity separation is only ~3.3x, against 78x
# for CD31/Capillary and 30x for B3TUBB/Nerve, and PDPN carries substantial
# non-lymphatic background (mesothelium, fibroblasts). Its pooled AUC (0.951)
# flatters it. Treat the PDPN gate as provisional and confirm before relying on
# any lymphatic call.
PIXEL_CLASSIFIER = {"Capillary": "CD31", "Nerve": "B3TUBB", "Lymphatic": "PDPN"}
WEAK_PIXEL_MARKERS = {"PDPN"}

# Lineage markers with no independent ground truth. Manual gate kept, flagged.
TIER2_LINEAGE = [
    "CD45", "CD68", "CD163", "CD3e", "CD4", "CD8a", "CD20", "HLADR",
    "panCK", "ECAD", "CK19", "SMA", "ColIV", "VIM", "CD34",
]

# Minimum inside-structure cells per image before trusting a per-image ROC fit.
MIN_INSIDE = 150
# Minimum AUC before a pixel classifier is accepted as ground truth.
MIN_AUC = 0.70


def groovy_gates(a, ann, images, log_core, im_core) -> list[dict]:
    """Per-image INS/GCG/SST gates taken from QuPath's own thresholds.

    `islet_composition_<case>.tsv` records the literal per-donor QuPath
    intensity thresholds as `Ins_thr` / `Glu_thr` / `Stt_thr`, constant within a
    donor and expressed in the same raw-MFI space as `adata.raw.X`. Using them
    directly is strictly better than fitting a gate to match groovy's counts:
    a fitted gate reproduces the counts by construction and so cannot be
    validated against them, whereas these are external values that either
    reproduce the counts or do not.

    They do. Applying `Ins_thr` etc. to raw MFI lands inside groovy's own
    [X_pos, X_pos + Multi_Pos] bounds for every marker/image spot-checked,
    while the broadcast manual gate misses badly in both directions (SST on
    image 6563: manual yields 17 cells against a groovy range of 810-1278).

    Thresholds are converted to log1p space because that is where the gates are
    compared downstream.
    """
    rows = []
    var_names = list(a.var_names)
    tsv_col = {"INS": "Ins_thr", "GCG": "Glu_thr", "SST": "Stt_thr"}

    # groovy per-image counts, used only as a CHECK on the external thresholds
    isl = ann[ann["Classification"] == "Islet"]
    per_img = isl.groupby("Image")[
        ["Num GCG_pos", "Num INS_pos", "Num SST_pos", "Num Multi_Pos"]
    ].sum()

    for im in images:
        tsv = GROOVY_DIR / f"islet_composition_{int(im):04d}.tsv"
        if not tsv.exists():
            print(f"  WARNING: {tsv.name} missing; INS/GCG/SST fall back to manual")
            continue
        comp = pd.read_csv(tsv, sep="\t")
        for m in GROOVY_MARKERS:
            col = tsv_col[m]
            if col not in comp.columns or comp.empty:
                continue
            raw_thr = float(comp[col].iloc[0])
            gate = float(np.log1p(raw_thr))  # layers['log'] is log1p(raw MFI)
            j = var_names.index(m)
            sel = im_core == im
            v = log_core[sel, j]
            achieved = int((v > gate).sum()) if v.size else -1

            lo = hi = -1
            if int(im) in per_img.index:
                lo = int(per_img.loc[int(im), f"Num {m}_pos"])
                hi = lo + int(per_img.loc[int(im), "Num Multi_Pos"])
                if not (lo <= achieved <= hi):
                    print(
                        f"  NOTE: {m} image {im}: {achieved} positive, outside groovy"
                        f" bounds [{lo}, {hi}] -- threshold may not transfer cleanly"
                    )
            rows.append(
                {
                    "markers": m,
                    "imageid": str(im),
                    "gate": gate,
                    "source": "qupath_explicit_threshold",
                    "tier": 1,
                    "target_n": lo,
                    "achieved_n": achieved,
                    "n_cells": int(v.size),
                    "auc": np.nan,
                    "raw_threshold": raw_thr,
                }
            )
    return rows


def pixel_gates(a, obs, images, log_all, im_all) -> list[dict]:
    """Per-image CD31/B3TUBB/PDPN gates at the Youden-J optimum."""
    from sklearn.metrics import roc_auc_score, roc_curve

    rows = []
    var_names = list(a.var_names)
    for struct, m in PIXEL_CLASSIFIER.items():
        col = f"Dist to Closest {struct}"
        if col not in obs.columns:
            print(f"  WARNING: {col} absent; skipping {m}")
            continue
        j = var_names.index(m)
        dist = obs[col].astype(float).values
        for im in images:
            sel = im_all == im
            d, v = dist[sel], log_all[sel, j]
            ok = np.isfinite(d) & np.isfinite(v)
            d, v = d[ok], v[ok]
            # Signed distance < 0 => centroid lies inside the structure.
            y = (d < 0).astype(int)
            if y.sum() < MIN_INSIDE or (len(y) - y.sum()) < MIN_INSIDE:
                continue
            auc = float(roc_auc_score(y, v))
            if auc < MIN_AUC:
                print(
                    f"  WARNING: {m} vs {struct} image {im}: AUC {auc:.3f} < {MIN_AUC}"
                    " -- not used"
                )
                continue
            fpr, tpr, th = roc_curve(y, v)
            gate = float(th[int(np.argmax(tpr - fpr))])
            rows.append(
                {
                    "markers": m,
                    "imageid": str(im),
                    "gate": gate,
                    "source": f"pixel_classifier:{struct}",
                    "tier": 1,
                    "target_n": int(y.sum()),
                    "achieved_n": int((v > gate).sum()),
                    "n_cells": int(len(v)),
                    "auc": auc,
                }
            )
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--n-sample",
        type=int,
        default=300_000,
        help="cells sampled for the pixel-classifier ROC fits (0 = all)",
    )
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    import anndata as ad

    for p in (H5AD, ANNOTATIONS):
        if not p.exists():
            print(f"ERROR: {p} not found", file=sys.stderr)
            return 1

    print(f"Reading {H5AD.name} (backed)...")
    a = ad.read_h5ad(H5AD, backed="r")
    obs = a.obs
    if "log" not in a.layers:
        print("ERROR: layers['log'] missing; gates live in log1p space", file=sys.stderr)
        return 1

    var_names = list(a.var_names)
    img = obs["imageid"].astype(str).values
    images = pd.unique(img)
    parent = obs["Parent"].astype(str).values
    core = pd.Series(parent).str.match(r"^Islet_\d+$").values
    print(f"  {a.n_obs:,} cells, {len(images)} images, {core.sum():,} islet-core cells")

    print(f"Reading {ANNOTATIONS.name}...")
    ann = pd.read_csv(ANNOTATIONS, sep="\t", low_memory=False)

    # --- Tier 1a: groovy-calibrated endocrine gates (islet-core cells only) ---
    print("\nTier 1a: calibrating INS/GCG/SST against groovy endocrine counts...")
    ci = np.where(core)[0]
    gi = [var_names.index(m) for m in GROOVY_MARKERS]
    log_core = np.asarray(a.layers["log"][ci, :][:, gi])
    # re-index columns so log_core column order matches var_names lookups below
    full_core = np.full((log_core.shape[0], len(var_names)), np.nan)
    for k, j in enumerate(gi):
        full_core[:, j] = log_core[:, k]
    rows = groovy_gates(a, ann, images, full_core, img[ci])
    print(f"  {len(rows)} gates derived")

    # --- Tier 1b: pixel-classifier gates ---
    print("\nTier 1b: calibrating CD31/B3TUBB/PDPN against pixel classifiers...")
    rng = np.random.default_rng(RANDOM_STATE)
    if args.n_sample and a.n_obs > args.n_sample:
        si = np.sort(rng.choice(a.n_obs, args.n_sample, replace=False))
    else:
        si = np.arange(a.n_obs)
    pi = [var_names.index(m) for m in PIXEL_CLASSIFIER.values()]
    log_s = np.asarray(a.layers["log"][si, :][:, pi])
    full_s = np.full((log_s.shape[0], len(var_names)), np.nan)
    for k, j in enumerate(pi):
        full_s[:, j] = log_s[:, k]
    prows = pixel_gates(a, obs.iloc[si], images, full_s, img[si])
    print(f"  {len(prows)} gates derived")
    rows += prows

    # --- Fallback: retain the manual gate wherever no calibration exists ---
    #
    # Tier is a property of the MARKER (what evidence could exist for it), not
    # of whether this particular image happened to fit. A tier-1 marker whose
    # ROC fit failed on one image is still tier 1 -- it just falls back for that
    # image, which `source` records. Conflating the two previously filed CD31 /
    # B3TUBB / PDPN under both tier 1 and tier 3.
    tier1_markers = set(GROOVY_MARKERS) | set(PIXEL_CLASSIFIER.values())

    def marker_tier(m: str) -> int:
        if m in tier1_markers:
            return 1
        return 2 if m in TIER2_LINEAGE else 3

    manual = a.uns["gates"]["6450"]
    calibrated = {(r["markers"], r["imageid"]) for r in rows}
    for m in var_names:
        for im in images:
            if (m, str(im)) in calibrated:
                continue
            rows.append(
                {
                    "markers": m,
                    "imageid": str(im),
                    "gate": float(manual[m]),
                    "source": "manual_6450_broadcast",
                    "tier": marker_tier(m),
                    "target_n": -1,
                    "achieved_n": -1,
                    "n_cells": -1,
                    "auc": np.nan,
                }
            )

    prov = pd.DataFrame(rows)
    wide = prov.pivot(index="markers", columns="imageid", values="gate")

    print("\n--- Gate summary by tier ---")
    for tier, label in [
        (1, "calibrated against ground truth"),
        (2, "lineage marker, NO ground truth (manual 6450 retained)"),
        (3, "not used in broad lineage; deferred"),
    ]:
        mk = sorted(prov[prov.tier == tier]["markers"].unique())
        print(f"  tier {tier} ({label}): {len(mk)} markers")
        if mk:
            print(f"    {', '.join(mk)}")

    fell_back = prov[(prov.tier == 1) & (prov.source == "manual_6450_broadcast")]
    if not fell_back.empty:
        print(
            f"\n  {len(fell_back)} tier-1 (marker, image) pairs FELL BACK to the manual"
            " gate because calibration was unavailable:"
        )
        for m, grp in fell_back.groupby("markers"):
            print(f"    {m}: images {sorted(grp.imageid.tolist())}")
        print("    -> these images are uncalibrated for that marker; treat with care.")

    t1 = prov[prov.tier == 1]
    if not t1.empty:
        print("\nTier-1 per-image gate spread vs the single manual gate:")
        summ = t1.groupby("markers").agg(
            gate_min=("gate", "min"),
            gate_med=("gate", "median"),
            gate_max=("gate", "max"),
            n_images=("gate", "size"),
            mean_auc=("auc", "mean"),
        )
        summ["manual_6450"] = [manual[m] for m in summ.index]
        summ["med_delta"] = summ["gate_med"] - summ["manual_6450"]
        print(summ.to_string(float_format=lambda x: f"{x:8.2f}"))

    if args.dry_run:
        print("\n--dry-run: not writing")
        return 0

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    wide.to_csv(OUT_CSV)
    prov.to_csv(OUT_PROV, index=False)
    print(f"\nWrote {OUT_CSV}  ({wide.shape[0]} markers x {wide.shape[1]} images)")
    print(f"Wrote {OUT_PROV}  (per-gate source/tier/fit diagnostics)")
    print(
        "\nNOTE: tier-2 gates are the unvalidated manual 6450 values. They are the"
        "\n      next thing to address if a lineage call looks wrong."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
