# Trajectory / Pseudotime — Accuracy Reference

Companion to [`statistics.md`](statistics.md). Linked from `app/shiny_app/CLAUDE.md`.

This file is the in-repo record of the **pseudotime accuracy audit (2026-07-20/21)**, plus
the follow-up measurements of 2026-08-14. Until then the audit's findings lived only in
session plan files and agent memory, so a reader of this repository had no way to learn that
it had happened — or which of its confident-sounding predecessors had been disproved.

Every number below was **measured**, not inferred. Where a claim was tested and failed, it is
recorded as failed rather than quietly dropped, because two of them keep getting re-proposed.

---

## 1. Ground truth for the shipped file

`data/adata_ins_root.h5ad`, re-measured 2026-08-14:

```
uns['neighbors']['params'] = {use_rep: 'X_scVI_harmony_pt', n_neighbors: 30, metric: 'cosine'}
uns['iroot']               = 703            (donor 6516, INS 0.663, 99.7th pct)
diffmap_evals[:5]          = [1.0, 0.9399, 0.9268, 0.9142, 0.9060]
```

| | core (`dpt_pseudotime`) | combined (`dpt_pseudotime_combined`) |
|---|---|---|
| INS Spearman | −0.570 | −0.511 |
| GCG Spearman | 0.143 | 0.103 |
| SST Spearman | −0.206 | −0.225 |
| donor-mean ND / Aab+ / T1D | 0.3149 / 0.3668 / 0.3870 | 0.2317 / 0.2649 / 0.2853 |
| islet-mean ND / Aab+ / T1D | 0.3040 / 0.3634 / 0.3883 | 0.2245 / 0.2619 / 0.2867 |
| donor η² ND / Aab+ / T1D | 0.034 / 0.090 / 0.024 | 0.024 / 0.078 / 0.041 |

Core vs combined agreement: **Spearman ρ = 0.846** — the two modes are *not* interchangeable,
and the Trajectory tab's "Pseudotime mode" selector is a real analytical choice, not a display
preference.

> **The committed docs were stale until 2026-08-14.** `scripts/CLAUDE.md` described "raw
> `X_scVI_mean` (no Harmony)" with `n_neighbors=15`, INS −0.69 and GCG 0.274 — figures from a
> raw-scVI variant that never shipped. Verify against the file, not the prose.

---

## 2. Confirmed real

- **Donor identity dominates disease.** Donor η² = 0.0913 vs status η² = 0.0445 — a **2.05×**
  ratio. Any per-pseudotime display needs a per-donor control; Steenbuck et al. 2025 use a
  donor rug strip under the main panel (their Fig 5E).
- **ND vs Aab+ is never separable** at donor level (p = 0.111 / 0.421 / 0.155 across
  min-cells thresholds). **T1D drives the entire signal.** The UI must not imply a clean,
  evenly-spaced three-stage ordering — and indeed the Aab+→T1D step (0.020 donor-level) is
  under half the ND→Aab+ step (0.052).
- **External validation passes.** Donor-mean pseudotime vs β-cell fraction r = −0.850 (all
  islets) to **−0.949** (min_cells ≥ 25), against Bodenmiller's −0.74. Partly circular — the
  same panel produces both sides — as is theirs.

## 3. Claims from an earlier audit that measurement REFUTED

These are recorded because they are confident, plausible, and wrong. Do not re-litigate.

| Claim | Verdict |
|---|---|
| "ND<Aab+<T1D is essentially one donor" | **False.** Donor-level means stay monotonic; Spearman(donor mean, disease code) ρ = +0.680, p = 0.0053; **leave-one-donor-out: 0 of 15 break the ordering.** The validation *is* pseudoreplicated and was fixed on principle, but the conclusion holds. |
| "min_cells ≥ 25 nearly erases the separation" | **False, and backwards.** The donor-level trend *strengthens*: ρ +0.680 → +0.643 → **+0.813 (p = 0.0002, d = 3.95)**. |
| "No trajectory structure, it's a blob" | **Overstated.** λ₁/λ₂ = 1.0141 is weak, but there is a leading gap (1.74× the later ones). A caveat, not a blocker. |

## 4. Not viable — do not re-propose

- **Lead/lag (`ccf` + Granger).** Raw ccf peaks are trend artifacts (every marker trends
  monotonically); differencing *inverts* every estimate; finer binning does not converge
  (peaks wander +2, −3, +15, −1, +0). **Autocorrelation, not bin count, is binding**: AR1 =
  0.95 at 20 bins (n_eff ≈ 0), still only n_eff = 48 at 300 bins. Cell-level pseudotime would
  make this *worse* — cleaner bins mean higher autocorrelation and lower effective n. The real
  gap vs Bodenmiller is **88 donors vs our 15**; their 7,025 ROIs is only 1.35× our 5,214
  islets. If ever revisited, gate empirically on a bin-count convergence check.
- **Cell-level β-cell pseudotime (literal Bodenmiller).** 28,221 β-cells vs their >300,000,
  and only **522 in T1D** — 413 from one donor, 3 from another. Confirmed **genuine biology,
  not a gating artifact**: T1D images are brighter on BActin and comparable on
  NaKATPase/BCatenin/panCK, with GCG higher; INS is the sole marker that collapses (~58×).
  (Groovy ground truth shows scimap still *exaggerates* that collapse ~2×.)
- **Lamian.** `compute_lamian_pseudotime()` was 157 lines with zero call sites; deleted in
  `1159ebe`. The app has no Lamian code path.

## 5. Known-imperfect, accepted

**The root is a medoid, and a medoid is the wrong shape for a root.** It is the centre of mass
of the top-50 INS-high ND islets, so it sits at the **88.2nd percentile of DC1 — not a tip**.
DPT's unsigned distance therefore folds. Impact is tail-only: the 11.8% beyond the root is the
healthy ND tail and correctly receives *low* pseudotime; contamination is 3.3% of the top
quartile and 15.7% of the top 5%.

**Accepted, not fixed.** Slingshot would eliminate the problem structurally (arc length, no
root) and `scripts/compute_pseudotime.R` is already a working destiny→slingshot template, but
switching is a method change requiring full revalidation, not a bug fix. Do not re-run for
this alone.

**Upstream contaminant, flagged not fixed.** `obsm['X_scVI_mean']` averages *all* islet-core
cells with no phenotype mask (`islet_analysis/fixed_islet_aggregation.py:341-346`). 21.5% of
islet-core cells are acinar, so the embedding feeding DPT is diluted by exocrine signal —
worst in T1D, where β-cells are scarcest. See [[phenotyping_gates_ground_truth]]: phenotyping
sits upstream of all of this.

---

## 6. Defects found and fixed

### App: the pseudotime axis was unstable (fixed)

`rank(pt)/n` ran **after** the `!is.na(value)` and `min_cells` filters, so the x-axis moved
whenever either changed — including merely switching the plotted Y-axis feature. Measured
displacement at min_cells = 25: **mean |Δrank| 0.179, max 0.357.**

Fixed by ranking **once**, over the full unfiltered islet set, in `traj_with_pseudotime()`
(`mod_trajectory_server.R`). `apply_pt_scale()` is shared by `traj_data_clean()` and
`traj_multi_feature_data()` so the scatter and the multi-feature heatmap cannot drift apart.
**Never re-rank a filtered subset.**

A Rank/Raw-DPT toggle was added, with the caveat that on the rank axis LOESS curvature and bin
density are artifacts of the transform (ranks are uniform by construction).

### App: the Trajectory CSV lied about its own filter (fixed 2026-08-14)

The provenance header read `input$traj_min_cells`, but the UI defines `ns("min_cells")`. The
input never existed, so the `%||% 1` fallback fired every time and the header **always**
reported "Min cells per islet: 1" regardless of the value applied to the rows. Same class as
the Statistics tab's silent min-cells no-op: a filter and its provenance line disagreeing.

Fixed, and both filter sites plus the header now route through `traj_min_cells_value()`, which
also closes the `NA`-from-a-cleared-`numericInput` crash (`%||%` catches NULL but not NA, so
`if (NA > 1)` aborted the reactive) that had already been fixed on the Statistics tab.

### Pipeline: two H5AD-integrity bugs (fixed 2026-08-14)

Neither affected the live app — it reads `dpt_pseudotime` from `.obs`, and nothing in
`app/shiny_app/R/` or `build_h5ad_for_app.py` consumes `DC1` or `X_diffmap` — but both
misled anyone working from the saved file.

1. **`DC1`/`DC2` stored the trivial λ=1 eigenvector.** `X_diffmap[:,0]` is scanpy's stationary
   component: measured std **0.000814**, strictly positive (0.0112→0.0163), carrying no
   trajectory information. The real first component (column 1) has std 0.0138, ~17× larger.
   Now columns 1 and 2.
2. **`obsp` held the *combined* graph while `uns` advertised core.** Step 5b's
   `sc.pp.neighbors` overwrites `obsp` in place, and the restore block put back `obs`, `obsm`
   and `uns` but never `obsp` — measured overlap 96.6% with the combined graph vs 34.6% with
   core. **Re-running scanpy on the shipped file silently produced combined-mode results.**
   `obsp`, `iroot` and `diffmap_evals` are now stashed and restored with the rest.

The superseded `notebooks/rebuild_trajectory.ipynb` still contains the original DC1 pattern.
It is not part of the pipeline; `scripts/reaggregate_islets.py` is authoritative.

### Pipeline: Step 6 validation was itself pseudoreplicated (fixed 2026-08-14)

The gate aggregated **islets**, not donors — the exact error the Statistics tab exists to
avoid — **never validated `dpt_pseudotime_combined` at all**, and printed a warning while
saving anyway, so a failure was invisible to automation.

Now: `_validate_trajectory()` runs `_validate_one_pseudotime()` over **both** modes; the
disease-ordering gate is computed on donor means (islet means still printed for reference);
core-vs-combined agreement is reported; and `main()` exits **2** on failure. Outputs are still
written first — discarding a ~15-minute run makes the gate too expensive to keep.

⚠️ **The GCG gate is the weak one.** The 0.1 threshold is relaxed from an original 0.2, and
the shipped data clears it only narrowly — core 0.143, **combined 0.103**. The relaxation
rests on a code comment asserting the genuine shared effect is ≈0.14; nothing independent
confirms that. The validator now prints `<-- NARROW` when the margin is under 0.05, and both
modes currently trip it.

---

## 7. Still open

- **`ALPHA_STRUCT = 0.15` is conditionally invalidated** by the pending CODEX re-phenotyping:
  3 of the 4 `STRUCT_PHENOTYPES` names disappear, and the structural block **shrinks silently**
  (a `print`, not an error). The sweep would need re-running — but `alpha_sweep.py`, cited at
  `reaggregate_islets.py:224` as its source, **has never existed** in the worktree or in git
  history. Re-deriving α means rewriting it from scratch.
- **The UI does not carry the audit's own conclusion.** Donor identity explains 2× the
  variance disease does and ND vs Aab+ is not separable, yet `docs/user_guide.md` still
  presents the trajectory as a clean three-stage progression. A donor rug strip under the
  scatter is the established convention.
- **`traj_validate()` is a misleading name** — it only checks that the selected feature exists
  in `var_names`. There is no in-app trajectory validation.
- **Combined mode's GCG margin (0.103)** is 0.003 above the gate. Worth watching on any rebuild.

## 8. How to re-verify

```bash
# Ground truth straight from the shipped file
python -c "
import anndata as ad, numpy as np
from scipy import stats
a = ad.read_h5ad('data/adata_ins_root.h5ad')
print(a.uns['neighbors']['params'])
pt = a.obs['dpt_pseudotime'].values; v = ~np.isnan(pt)
for m in ['INS','GCG']:
    i = list(a.var_names).index(m)
    print(m, round(stats.spearmanr(pt[v], a.X[v,i])[0], 3))
print('DC1 std', round(a.obsm['X_diffmap'][:,0].std(), 6))  # ~0.0138 after the fix
"

# Run the validation gates against any H5AD without rebuilding
python -c "
import importlib.util, anndata as ad
s = importlib.util.spec_from_file_location('r','scripts/reaggregate_islets.py')
m = importlib.util.module_from_spec(s); s.loader.exec_module(m)
a = ad.read_h5ad('data/adata_ins_root.h5ad')
print('PASS' if m._validate_trajectory(a, list(a.var_names).index('INS')) else 'FAIL')
"
```
