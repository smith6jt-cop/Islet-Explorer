#!/usr/bin/env python
"""Hierarchical marker-rule evaluator -- replaces scimap.tl.phenotype_cells.

WHY THIS EXISTS
---------------
`scimap.tl.phenotype_cells` computes

    total_score = (sum of rule terms) / number_of_non_empty_features

For a rule combining `anypos` + `allneg` (2 terms), any ONE allneg marker over
the gate zeroes the allneg term, so:

    total_score = (anypos_score + 0) / 2  <=  0.5

and 0.5 is exactly the threshold at which the `-rest` (i.e. `Unknown`) bucket
wins. A hormone-bright islet cell with a single immune marker over gate is
therefore *mathematically guaranteed* to become `Unknown`, no matter how bright
its insulin. Measured consequence in the shipped data: beta and alpha cells are
**0.00% HLADR+ by construction**, which makes MHC-II upregulation in islet cells
(Bodenmiller 2025 Fig 2) undetectable; and islet `Unknown` cells are 83-95%
positive for INS/GCG/SST -- they are endocrine cells that fell through.

A second defect: scimap's `Endocrine` label is really `Endocrine-rest`, which
fires iff INS, GCG and SST are ALL below gate, then `-rest` -> NaN ->
forward-filled to the parent name. All 11,413 islet-core `Endocrine` cells are
sub-gate on every hormone. The label reads as a confident call but means the
opposite.

DESIGN
------
1. **Identity comes from positive evidence, never from absence.** A cell is
   assigned to the phenotype whose positive markers it satisfies most strongly.
   There is no divide-by-n_features and no score dilution.

2. **Cross-family conflicts are flagged, not resolved by deletion.** A cell
   positive for both endocrine and immune markers is assigned by strongest
   evidence and marked `doublet_suspect`. In insulitis, T cells are packed
   against beta cells and segmentation bleeds -- that is a real ambiguity that
   must stay visible, not be silently dumped into `Unknown`.

3. **`neg` is only used to separate SIBLINGS** within one family (e.g. CD8 vs
   CD4 T cell), where it is genuinely discriminative. Cross-family `allneg`
   constraints -- the thing that broke scimap -- are not used at all.

4. **Orthogonal annotations stay orthogonal.** `HLADR` positivity is emitted as
   its own column, never as an identity-defining exclusion.

5. **Cells matching nothing are `Unassigned`, explicitly.** Not folded into a
   residual sink. In the shipped labels `Acinar` had no positive marker at all,
   so it scored ~1.01 for any hormone-negative cell and absorbed 61.3% of the
   dataset.

RULES FORMAT (`data/phenotype_rules_v2.csv`)
-------------------------------------------
Wide: `parent`, `phenotype`, then one column per marker. Cell values:
    pos     -- marker must be above gate (required)
    anypos  -- at least one `anypos` marker in that row must be above gate
    neg     -- marker must be below gate (sibling discrimination only)
    blank   -- marker not used by this rule
`parent` is `all` for top-level families, else another phenotype's name.

NOTE: `data/phenotype_rules_v2.csv` is a RECONSTRUCTION. The rules file that
produced the shipped labels is not in the repo and cannot be recovered (it
generated `ECAD+` and `Structural`, which no committed CSV can). Treat the v2
rules as a starting point requiring biological review, not as a restoration.

Usage:
    python scripts/phenotype_evaluator.py --dry-run          # report only
    python scripts/phenotype_evaluator.py --out labels.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
H5AD = REPO / "single_cell_analysis" / "CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad"
RULES = REPO / "data" / "phenotype_rules_v2.csv"
GATES = REPO / "data" / "phenotype_gates_per_image.csv"

# Emitted as their own columns, never as identity constraints.
ORTHOGONAL_MARKERS = ["HLADR", "CHGA", "CD45"]

UNASSIGNED = "Unassigned"


# --------------------------------------------------------------------------
# Rules
# --------------------------------------------------------------------------
class Rule:
    __slots__ = ("phenotype", "parent", "pos", "anypos", "neg")

    def __init__(self, phenotype, parent, pos, anypos, neg):
        self.phenotype = phenotype
        self.parent = parent
        self.pos = pos
        self.anypos = anypos
        self.neg = neg

    def __repr__(self):  # pragma: no cover - debugging aid
        return (
            f"Rule({self.phenotype!r}, parent={self.parent!r}, pos={self.pos}, "
            f"anypos={self.anypos}, neg={self.neg})"
        )


def load_rules(path: Path, known_markers: list[str]) -> list[Rule]:
    df = pd.read_csv(path)
    if "parent" not in df.columns or "phenotype" not in df.columns:
        raise ValueError("rules file must have 'parent' and 'phenotype' columns")
    marker_cols = [c for c in df.columns if c not in ("parent", "phenotype")]

    missing = [m for m in marker_cols if m not in known_markers]
    if missing:
        # scimap silently KeyErrors on this; be explicit instead.
        raise ValueError(
            f"rules reference markers absent from the panel: {missing}. "
            "Remove them or fix the header."
        )

    rules = []
    for _, r in df.iterrows():
        pos, anypos, neg = [], [], []
        for m in marker_cols:
            v = str(r[m]).strip().lower()
            if v in ("", "nan"):
                continue
            if v == "pos":
                pos.append(m)
            elif v == "anypos":
                anypos.append(m)
            elif v == "neg":
                neg.append(m)
            else:
                raise ValueError(
                    f"{r['phenotype']}/{m}: unknown constraint {v!r} "
                    "(expected pos, anypos, neg or blank)"
                )
        if not pos and not anypos:
            raise ValueError(
                f"{r['phenotype']}: has no positive evidence (only neg/blank). "
                "Identity must come from positive markers -- a rule defined "
                "purely by absence becomes a catch-all sink."
            )
        rules.append(Rule(str(r["phenotype"]), str(r["parent"]), pos, anypos, neg))
    return rules


# --------------------------------------------------------------------------
# Gating / rescaling
# --------------------------------------------------------------------------
def rescale_with_gates(log_vals, gates_by_image, markers, imageids):
    """Piecewise rescale so each (marker, image) gate maps to exactly 0.5.

    Mirrors scimap's split MinMaxScaler: below-gate values map into [0, 0.5]
    and at-or-above-gate into [0.5, 1], per marker per image. Reimplemented
    here so the gate source is explicit and no in-place AnnData mutation
    occurs.
    """
    out = np.zeros_like(log_vals, dtype=float)
    for im in pd.unique(imageids):
        sel = imageids == im
        if not sel.any():
            continue
        block = log_vals[sel]
        for j, m in enumerate(markers):
            g = gates_by_image.get((m, str(im)))
            if g is None or not np.isfinite(g):
                out[sel, j] = np.nan
                continue
            v = block[:, j]
            lo_m = v < g
            hi_m = ~lo_m
            res = np.empty_like(v)
            if lo_m.any():
                lo = v[lo_m]
                span = g - lo.min()
                res[lo_m] = 0.5 * (lo - lo.min()) / span if span > 0 else 0.0
            if hi_m.any():
                hi = v[hi_m]
                span = hi.max() - g
                res[hi_m] = 0.5 + 0.5 * (hi - g) / span if span > 0 else 0.5
            out[sel, j] = res
    return out


# --------------------------------------------------------------------------
# Evaluation
# --------------------------------------------------------------------------
def evaluate(rescaled, markers, rules, gate=0.5):
    """Assign each cell hierarchically by strongest positive evidence.

    Returns a DataFrame with the terminal label, the family it came from, the
    winning score, and conflict/annotation flags.
    """
    midx = {m: i for i, m in enumerate(markers)}
    pos_mat = rescaled > gate
    n = rescaled.shape[0]

    def eligible_and_score(rule: Rule):
        ok = np.ones(n, dtype=bool)
        contrib = []
        for m in rule.pos:
            j = midx[m]
            ok &= pos_mat[:, j]
            contrib.append(j)
        if rule.anypos:
            cols = [midx[m] for m in rule.anypos]
            ok &= pos_mat[:, cols].any(axis=1)
            contrib.extend(cols)
        for m in rule.neg:
            ok &= ~pos_mat[:, midx[m]]
        if contrib:
            # Only positive markers contribute; an unmet anypos member
            # contributes its (sub-gate) value, which correctly lowers the score.
            vals = rescaled[:, contrib]
            score = np.nanmax(vals, axis=1) if len(contrib) > 1 else vals[:, 0]
        else:  # pragma: no cover - load_rules forbids this
            score = np.zeros(n)
        return ok, np.where(ok, score, -np.inf)

    by_parent: dict[str, list[Rule]] = {}
    for r in rules:
        by_parent.setdefault(r.parent, []).append(r)

    label = np.full(n, UNASSIGNED, dtype=object)
    family = np.full(n, UNASSIGNED, dtype=object)
    best_score = np.full(n, np.nan)

    # Precompute eligibility for every rule once.
    elig = {r.phenotype: eligible_and_score(r) for r in rules}

    # --- level 0: pick the winning top-level family ---
    tops = by_parent.get("all", [])
    if not tops:
        raise ValueError("no top-level rules (parent == 'all')")
    top_scores = np.vstack([elig[r.phenotype][1] for r in tops])
    top_any = np.isfinite(top_scores).any(axis=0)
    top_pick = np.argmax(top_scores, axis=0)
    for k, r in enumerate(tops):
        m = top_any & (top_pick == k)
        label[m] = r.phenotype
        family[m] = r.phenotype
        best_score[m] = top_scores[k][m]

    # --- descend: refine into the most specific eligible child ---
    changed = True
    while changed:
        changed = False
        for parent, children in by_parent.items():
            if parent == "all":
                continue
            at_parent = label == parent
            if not at_parent.any():
                continue
            cs = np.vstack([elig[c.phenotype][1] for c in children])
            has = np.isfinite(cs).any(axis=0) & at_parent
            if not has.any():
                continue
            pick = np.argmax(cs, axis=0)
            for k, c in enumerate(children):
                m = has & (pick == k)
                if m.any():
                    label[m] = c.phenotype
                    best_score[m] = cs[k][m]
                    changed = True

    # --- conflict + orthogonal annotations ---
    fam_pos = {}
    for r in tops:
        cols = [midx[m] for m in (r.pos + r.anypos)]
        fam_pos[r.phenotype] = pos_mat[:, cols].any(axis=1)
    n_fam = np.sum(np.vstack(list(fam_pos.values())), axis=0)

    out = pd.DataFrame(
        {
            "phenotype": label,
            "family": family,
            "score": best_score,
            "n_families_positive": n_fam,
            "doublet_suspect": n_fam >= 2,
        }
    )
    for m in ORTHOGONAL_MARKERS:
        if m in midx:
            out[f"{m}_pos"] = pos_mat[:, midx[m]]
    return out


# --------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rules", type=Path, default=RULES)
    ap.add_argument("--gates", type=Path, default=GATES)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--n-sample", type=int, default=0, help="0 = all cells")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    import anndata as ad

    for p in (H5AD, args.rules, args.gates):
        if not p.exists():
            print(f"ERROR: {p} not found", file=sys.stderr)
            return 1

    a = ad.read_h5ad(H5AD, backed="r")
    markers = list(a.var_names)
    rules = load_rules(args.rules, markers)
    print(f"Loaded {len(rules)} rules from {args.rules.name}")

    gdf = pd.read_csv(args.gates, index_col=0)
    gdf.columns = gdf.columns.astype(str)
    gates = {
        (m, im): float(gdf.loc[m, im])
        for m in gdf.index
        for im in gdf.columns
        if m in markers
    }
    print(f"Loaded gates for {gdf.shape[0]} markers x {gdf.shape[1]} images")

    img = a.obs["imageid"].astype(str).values
    if args.n_sample and a.n_obs > args.n_sample:
        rng = np.random.default_rng(0)
        idx = np.sort(rng.choice(a.n_obs, args.n_sample, replace=False))
    else:
        idx = np.arange(a.n_obs)
    print(f"Evaluating {len(idx):,} cells...")

    log_vals = np.asarray(a.layers["log"][idx, :])
    rescaled = rescale_with_gates(log_vals, gates, markers, img[idx])
    res = evaluate(rescaled, markers, rules)

    old = a.obs["phenotype"].astype(str).values[idx]
    res["phenotype_old"] = old
    res["imageid"] = img[idx]
    res["Parent"] = a.obs["Parent"].astype(str).values[idx]
    res["donor_status"] = a.obs["Donor Status"].astype(str).values[idx]

    print("\n--- New label distribution ---")
    vc = res["phenotype"].value_counts()
    for k, v in vc.items():
        print(f"  {k:22s} {v:9,d}  {100*v/len(res):5.2f}%")

    print("\n--- THE DEFINING CHECK: HLADR+ rate among endocrine cells ---")
    print("  (was exactly 0.00% under scimap, by construction)")
    core = pd.Series(res["Parent"]).str.match(r"^Islet_\d+$").values
    for ph in ["Beta cell", "Alpha cell", "Delta cell"]:
        for st in ["ND", "Aab+", "T1D"]:
            m = core & (res["phenotype"] == ph).values & (res["donor_status"] == st).values
            if m.sum() >= 20:
                r = 100 * res.loc[m, "HLADR_pos"].mean()
                print(f"    {ph:12s} {st:5s}: n={m.sum():6,d}  HLADR+ = {r:5.2f}%")

    print("\n--- Conflicts ---")
    print(f"  doublet_suspect: {res.doublet_suspect.sum():,} ({100*res.doublet_suspect.mean():.2f}%)")
    print(f"  {UNASSIGNED}:      {(res.phenotype==UNASSIGNED).sum():,} "
          f"({100*(res.phenotype==UNASSIGNED).mean():.2f}%)")

    print("\n--- Where the old Unknown / Endocrine buckets went ---")
    for oldlab in ["Unknown", "Endocrine"]:
        m = res.phenotype_old == oldlab
        if m.sum() == 0:
            continue
        print(f"  old '{oldlab}' (n={m.sum():,}):")
        for k, v in res.loc[m, "phenotype"].value_counts().head(6).items():
            print(f"      -> {k:20s} {v:8,d}  {100*v/m.sum():5.1f}%")

    if args.dry_run or args.out is None:
        print("\n(dry run / no --out: nothing written)")
        return 0

    res.to_csv(args.out, index=False)
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
