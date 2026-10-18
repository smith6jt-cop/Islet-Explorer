# Trajectory / Pseudotime — Accuracy Reference

> **2026-10-03 implementation:** the expanded cohort preserves the original source data and uses Harmony on protein PCA for spatial UMAP and Leiden. The live trajectory uses endocrine-normalized composition features; marker expression is the mean gate-rescaled protein signal. Earlier scVI descriptions below are historical. See the [current implementation record](../../../docs/protein_expansion_20261003.md).

Companion to [`statistics.md`](statistics.md). Linked from `app/shiny_app/CLAUDE.md`.

### Retired-feature correction (2026-10-14)

The live trajectory previously selected **all 42 core/peri phenotype ratios**,
even though nine phenotype categories had been retired elsewhere in the app.
This created a cohort artifact: Acinar, ECAD+ and Structural are absent by
construction in the added donors, while their unresolved cells accumulate in
Unknown. Original donors retain their historical labels. For islets with
25–1902 cells, original/added cohort membership explained 99.3% of raw DPT
variation, and every islet in the highest pseudotime decile was an addition.
At 25–200 cells the graph split exactly by cohort, leaving all 178 additions
unreachable from the original-donor root.

Both Scanpy and destiny now use the app's **12 retained phenotypes**, in core
and peri separately (24 candidate columns, followed by the existing constant
column removal). This is the existing phenotype policy, not a feature set
chosen to maximize insulin correlation. The stored endocrine denominator,
source data, ND-centroid root rule, beta-ratio orientation convention,
15-neighbor cosine graph and range-dependent fitting are preserved. PAGA
uses the same corrected graph.

Measured on expanded-v2, before/after this feature correction:

| Cell range | Islets in range | INS Spearman before → after | Cohort share of raw-DPT variance before → after |
|---|---:|---:|---:|
| 10–1902 | 2,696 | −0.169 → −0.288 | 44.6% → 2.7% |
| 25–1902 | 1,120 | −0.265 → −0.801 | 99.3% → 1.7% |
| 25–200 | 957 | −0.497 → −0.783 | Before: disconnected; after: 2.8% |

The old 25–200 correlation includes only 779 reachable original islets;
the corrected graph reaches all 957. These are descriptive islet-level
diagnostics, not donor-level significance tests. The original/addition
labeling and intensity-scale differences still require scientific validation;
removing retired categories does not establish full method equivalence.

The recent rendering/deferred-PAGA changes were checked separately: the old
and new calculation code produced identical DPT on identical inputs. The live
CSV matched a direct calculation within 9e-7, marker means matched within
3e-8, and the original donors' stored marker means were exactly preserved.
The saved September 23 composition-axis export already had a high-insulin
tail at the 10-cell threshold. Older INS-rooted scVI columns describe a
different method and must not be used as expected values for this fit.

DPT measures distance from its chosen root; eigenvector signs cannot reverse
that distance. Its current root and beta-ratio convention do not guarantee a
monotonic INS curve, particularly for small islets. The correction changes
computed pseudotime values; the historical fitted values below are not the
expected results of the corrected 24-feature calculation.

Reproduce the comparison without modifying data:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=2 \
  python scripts/tests/audit_live_trajectory_inputs.py \
  --data-root data/development/protein-expanded-20261003-v2 \
  --output tmp/trajectory-audit-20261014/feature-audit.json
```

### Linked islet inspection (2026-10-14)

The pseudotime scatter and both UMAPs share a donor + full islet-name key.
Hovering an islet previews it in all three plots; clicking pins it and opens
the existing segmentation/details card. Moving away restores the pinned
highlight. **Clear selection** and the details card's **Close** button clear it.
Display changes preserve the pin while the islet remains in the plotted data;
recalculating a range that excludes it clears the selection. An islet without
UMAP coordinates can still be highlighted in the scatter. PAGA and heatmaps
are not part of this interaction.

The UMAPs are interactive Plotly widgets, while their PNG downloads still use
the original ggplots. Browser-local SVG rings follow the actual plotted
coordinates, including scatter jitter, zoom and resize. Hover events do not
reach Shiny, recompute pseudotime, redraw the point cloud or change axis ranges.
Clicks use one server-validated selection route; coordinate-nearest matching
is not used. The interaction check is
[`check-trajectory-links.mjs`](../../../deploy/browser/check-trajectory-links.mjs).

This file is the in-repo record of the **pseudotime accuracy audit (2026-07-20/21)**, plus
the follow-up measurements of 2026-08-14. Until then the audit's findings lived only in
session plan files and agent memory, so a reader of this repository had no way to learn that
it had happened — or which of its confident-sounding predecessors had been disproved.

Every number below was **measured**, not inferred. Where a claim was tested and failed, it is
recorded as failed rather than quietly dropped, because two of them keep getting re-proposed.

> **Current-release routing:** Sections 1–10 audit legacy H5AD pseudotime columns and
> experiments; those axes are retained as provenance artifacts but are no longer selectable
> in the reviewer UI. Section 11 is the authority for current reviewer behavior.

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

---

## 9. What features the pseudotime is actually computed on (Aug 2026)

Asked often enough to belong at the top of any pseudotime discussion, and the answer is not what
the axis labels suggest.

**Pseudotime is computed on a 10-dimensional scVI latent embedding, not on the 31 markers.**

```
single cells (2.6M x 31 markers, scVI-trained with batch covariates)
   |  mean of the 10-d scVI latent over the cells in each islet core
   v
obsm['X_scVI_mean']            (5,214 x 10)   <- the ONLY input to the trajectory
   |  Harmony, batch = imageid (donor), theta = 0.1
   v
obsm['X_scVI_harmony_pt']      (5,214 x 10)
   |  sc.pp.neighbors(n_neighbors = 30, metric = 'cosine')
   |  sc.tl.diffmap(n_comps = 10)  ->  sc.tl.dpt()
   v
obs['dpt_pseudotime']
```

Consequences worth stating plainly:

- **The 31 markers enter only through scVI.** They are not fed to the diffusion map directly and
  they are not weighted or selected. Whatever scVI encoded is what the trajectory sees.
- **INS is used for the ROOT, not for the axis.** The root is the medoid, in the harmonised latent,
  of the top-50 ND islets ranked by mean INS (`K_ROOT = 50`). A single-islet `argmax(INS)` root was
  replaced because it is brittle to staining outliers. After that, INS has no further role.
- **The marker picked as the Y axis is display only.** It comes from `.X` (mean marker expression
  per islet) and does not participate in the pseudotime at all. `Spearman(pt, INS) = -0.570` on the
  shipped core column is therefore a *validation* statistic, not a definition.
- **Cell-type proportions, morphology and structures do not enter the core pseudotime.**
  `obsm['phenotype_proportions']`, `obsm['phenotype_counts']` and `obsm['morphology']` are carried
  in the file but unused by Step 5. The `combined` mode is the exception: it concatenates the peri
  scVI latent plus z-scored structural blocks scaled by `ALPHA_STRUCT = 0.15`.
- **Each islet is a MEAN over its cells.** An islet of 3 cells and one of 1,902 cells are both a
  single point, weighted equally in the neighbour graph. This is the origin of the noise problem
  below.

## 10. Cell-count-stratified pseudotime (Aug 2026)

### Why

The shipped `dpt_pseudotime` is computed over all 5,214 islets at once, and 60.8% are aggregates of
<= 10 cells whose mean scVI latent is correspondingly noisy. Measured:

```
spearman(dpt_pseudotime, total_cells) = -0.274      smaller islets -> higher pseudotime
islets with dpt > 0.6: median 7 cells               vs 10-11 in the bulk
min_cells >= 25 collapses the tail: frac(dpt>0.6) 3.97% -> 0.44%, max 1.00 -> 0.76
```

Because scanpy normalises DPT by dividing by the maximum reachable distance, that tail sets the top
of the axis for everyone: **96.0% of islets land below 0.6**, and the bulk (to p95) occupies only
57.4% of the drawn axis.

### What was done

`scripts/compute_split_pseudotime.py` recomputes DPT *within* each cell-count stratum, repeating
Step 5 exactly (Harmony theta=0.1 on `X_scVI_mean`, n_neighbors=30 cosine, diffmap n_comps=10),
including **re-deriving the root inside each stratum** — the all-islet root is not necessarily a
member of a stratum, and a root outside the graph makes DPT undefined.

New columns in `data/adata_ins_root.h5ad` (NaN outside their stratum, which is what subsets the
plot — the app's ordinary `is.finite()` filters do it, so the axis and the point set cannot
disagree):

| column | stratum | n |
|---|---|---|
| `dpt_pseudotime_small` | `total_cells <= 10` | 3,171 |
| `dpt_pseudotime_large` | `total_cells > 10` | 2,043 |
| `cell_count_stratum` | `'small'` / `'large'` | 5,214 |

Parameters and per-stratum validation land in `uns['split_pseudotime']`. The shipped
`dpt_pseudotime` column is left **bit-identical** (verified by md5 against a backup); only the three
columns above are added and `obsm` is untouched.

### Results

```
                    all (shipped)   small (<=10)   large (>10)
n                        5,214          3,171          2,043
Spearman(pt, INS)       -0.570         -0.507         -0.664
Spearman(pt, GCG)       +0.143         +0.386         +0.169
Spearman(pt, SST)       -0.206         -0.211         -0.188
donor-mean ND           0.3149         0.3912         0.1654
donor-mean Aab+         0.3668         0.4662         0.1874
donor-mean T1D          0.3870         0.4890         0.2396
monotone ND<Aab+<T1D      yes            yes            yes
frac(pt > 0.6)          3.97%         17.75%          1.62%
bulk-to-p95 axis use    57.4%          67.0%          36.1%
```

**The real gain is signal, not axis usage.** Large islets give `Spearman(pt, INS) = -0.664` against
the all-islet -0.570 — the cleanest INS correlation this pipeline has produced, consistent with the
small-islet aggregates having diluted it. Small islets separately show a much stronger GCG
correlation (+0.386 vs +0.143), worth a look on its own. Both strata preserve the monotone
ND < Aab+ < T1D donor ordering.

**Splitting did NOT fix the axis compression, and for large islets it made it worse** (bulk-to-p95
36.1% vs 57.4%). That is the honest result: the compression is a property of scanpy's divide-by-max
normalisation, which manufactures a long right tail in *any* stratum — the small islets were not
solely responsible. The axis-limits fix (clipping the *view* to ~p99 with `coord_cartesian()`,
leaving the data untouched) is orthogonal and still outstanding.

### >>> The strata are NOT on a common scale <<<

Each is an independent diffusion map: separate Harmony fit, separate root islet, separate
divide-by-max normalisation. `small = 0.4` and `large = 0.4` do not mean the same thing. The app
enforces this rather than trusting the reader — `cmp_scale_warning` in `mod_trajectory_server.R`
renders a red banner whenever the two pseudotime cards show different strata, an amber one when
either shows a stratum, and the comparison card's CSV header carries the same warning.
**Compare shapes and orderings within a stratum; never merge the two columns onto one axis.**

## 11. The trajectory tab as it now stands (Aug 2026)

Consolidates what were sections 11-14 during the Aug 2026 rework. Read this, not the git history of
those sections.

### One pseudotime, no modes

The tab briefly carried five selectable pseudotime "modes" (core / core+peri scVI, two cell-count
strata, and a composition one). That is gone. Every extra mode was another axis whose numbers were
not comparable with the others', which is a garden-of-forking-paths surface rather than a feature,
and it invited exactly the cross-axis comparison the warning banners existed to forbid.

The one that survived is the **islet-composition diffusion map**, on measurement: rho -0.70 to -0.73
against core beta-cell content above 25 cells/islet, versus -0.57 for the scVI variants, with the
least compressed axis of any option tried. It is *not* named after its source paper in the UI --
the method is described, not branded.

Method (Barlow et al., eLife 2024;13:RP100535, Methods -> "PAGA analysis"), all parameters theirs:

| step | value |
|---|---|
| features | cell-type counts, islet core + 20 µm peri, ÷ endocrine cells in the core, log1p |
| neighbours | cosine distance, 15 nearest |
| clustering | Leiden, resolution 1 |
| embedding | PAGA, then `sc.tl.draw_graph(init_pos = "paga")` = ForceAtlas2 |
| pseudotime axis | scanpy diffusion map + DPT; destiny DPT fallback when scanpy is unavailable |
| root | ND islet nearest the ND centroid in diffusion space |

`app/shiny_app/python/paga_pseudotime.py` (real scanpy) with `R/barlow_pseudotime.R` (destiny) as
the fallback when scanpy is unreachable from reticulate. Both are diffusion pseudotime over the same
graph on the same features.

The main scatter computes the selected cohort's neighbor graph and DPT first. Leiden clustering,
PAGA abstraction, and ForceAtlas2 run only after the reviewer selects **Build and show PAGA layout**.
The selected feature table and neighbor graph are cached and reused between those steps (up to two
recent cell-count ranges), so opening the PAGA card does not rebuild the neighbor graph or DPT.

Two deliberate deviations, stated not hidden: their hand-picked "temporarily omit 25 islets in the
middle of the map" path-isolation step is **not** ported (unauditable, and undefinable on our data),
and the root is restricted to ND islets because the globally nearest islet to the ND centroid was
Aab+ at several thresholds, which anchors the axis in the wrong stage.

### Cells-per-islet is a RANGE, and it defines the axis

> **Superseded 2026-10-17 / 2026-10-18.** The `cell_range` slider and Recalculate button were removed at the
> owner's direction (2026-10-17). For one day the endocrine-object classes chosen above the tabs (Islet /
> Cluster / Single) defined the fitted cohort; fitted on all three, ρ(pseudotime, core beta-cell ratio) was only
> −0.20 with the status medians barely separated, because one- to few-cell objects carry almost no composition.
> Since 2026-10-18 (owner: "Use 10 or more cells as the cutoff instead of a category") the trajectory fits every
> object with at least `TRAJ_MIN_CELLS` = 10 core cells, whatever its class, and the class selector is hidden on
> the Trajectory tab. On the uniform v3 root that is 1,265 of 9,496 objects (1,242 Islets, 23 Clusters):
> ρ = −0.70, median pseudotime ND 0.31 < Aab+ 0.39 < T1D 0.75; PAGA takes 10 s. The section below is the
> historical record.

The old `min_cells` numericInput is replaced by a two-ended slider, `cell_range`, defaulting to
**10 - <largest islet>** (1,902 cells; read from the data by `traj_max_islet_cells()`, with a
literal fallback). 10 is Barlow's own exclusion threshold, and below it an islet's composition
vector is mostly sampling noise -- 60.8% of our islets are <= 10 cells.

The upper bound matters as much as the lower: the pseudotime is **refitted on whatever islets the
applied range selects**, so capping it isolates a size band instead of letting the few 1,000+
cell islets dominate the embedding. The slider moves in five-cell steps without recomputing;
**Recalculate trajectory** applies the pending range. Loaded and fitted islet counts appear next
to the slider. Until the button is pressed, plots and downloads retain the previous fit.

```
range        n      median pseudotime
10-1902    2,310         0.2525
25-1902      908         0.3825
10-200     2,181         0.2386
```

This is not the §6 axis-instability defect. That was the axis moving as a side effect of changing
the *Y-axis feature*, which nobody asked for; it still cannot happen and is regression-tested
(`all.equal` across a feature switch). Here the axis moves only when the islet set that defines it
changes on recalculation; a progress bar appears during the fit.

### Y axis: markers AND cell types, never one replacing the other

`traj_y_mode` switches the measure; both remain available at all times.

- **Marker expression** -- scVI-denoised mean from `.X`, the 31 markers.
- **Cell-type ratio** -- `log1p(cells of type T / endocrine cells in the islet core)`, with a region
  control for **core (default)** / peri / core + peri. This is the default Y-axis measure.

Only the **retained** phenotypes are offered as ratio numerators: `barlow_ratio_types()` runs the
list through `phenotype_keep()`, so the nine retired labels (`Endocrine, Acinar, Endothelial,
T cell, Immune, Unknown, ECAD+, Structural, APCs`) never reappear -- 12 selectable types. Note
`Endocrine` remains the ratio **denominator**; it is retired as a plotted label, not as a
denominator.

Separate input ids per mode (`traj_feature` vs `traj_ratio_type`) so switching cannot leave a marker
name selected against ratio data, which would silently plot nothing. `traj_validate()` branches
before checking `var_names`, because cell-type names are deliberately not in the marker list.

The multi-feature heatmap initially selects INS, GCG, SST, CHGA, CD45, CD31, panCK and B3TUBB;
all 31 markers remain selectable. The scatter uses smaller points and wider vertical display
jitter to separate overlapping islets. Hover values and LOESS trends use the unjittered data;
the jitter does not change the fitted pseudotime or exported measurements.

**A ratio does not cure the zero pile-up** -- `log1p(0) = 0`, and Barlow's own Fig 2F has the same
problem. Counting core + peri or raising the range floor can reduce it:

```
                 core, all   core+peri, all   core+peri, >=25 cells
Macrophage         0.644         0.235              0.042
Beta cell          0.615         0.610              0.307
CD8a Tcell         0.956         0.864              0.718
CD4 Tcell          0.993         0.969              0.927
```

The zero-percentage callout was removed from the tab; the measured comparison remains here for reference.

### PAGA force-directed layout

Barlow's Fig 2B/2E, and the plot that actually separates islets. A feature-vs-pseudotime scatter
cannot be de-overplotted by any transform because the pile-up is at exact zero; a 2-D graph
embedding spreads points by construction. Measured on the same 133x83 occupancy grid at n = 2,310:

```
                                   coverage   densest pixel   overplotted
PAGA force-directed layout           17.1%          3            34.8%
scatter: INS vs pseudotime           14.0%          8            55.0%
scatter: prop_Beta cell (Plot tab)   12.4%        125            84.0%
scatter: prop_CD8a (Plot tab)         2.8%        129            96.4%
```

Colour by Pseudotime, Selected feature, Donor status, Leiden cluster or Donor. `coord_fixed()` and
blanked ticks: layout coordinates are unitless and isotropic, and only relative position means
anything. Rebuilds with the cell range, cached per range.

### Removed from this tab

- **Outlier exclusion** (checkbox, z slider, ghost-× layer, excluded-outlier table). With the
  pseudotime now refitted on the selected range, a second independent islet-dropping rule was a
  hidden interaction: which islets counted as outliers changed with the range, silently changing the
  embedding's input set. `trajectory_server()` still *accepts* `remove_outliers` /
  `outlier_threshold` because `app.R` passes them to every analysis module, but ignores them; the
  `traj_remove_outliers` / `traj_outlier_z` sync observers were removed from `app.R`.
- **The UMAP show/hide toggle.** The UMAP card is back to always visible, as before.
- **The second (comparison) pseudotime card.** Superseded by there being only one pseudotime.

### scanpy runtime

Production installs the exact Python environment from `deploy/requirements-python.lock` and
selects it through `RETICULATE_PYTHON`. Local development may use another compatible environment.
Without scanpy, the axis falls back to destiny and the PAGA card reports that it is unavailable.
Without `fa2-modified`, scanpy can fall back to Fruchterman–Reingold; the card reports the resulting
layout key rather than hiding the change.

## 12. Restoration audit — 2026-09-23

Sections 1–11 above are restored verbatim from `000e970`. The default app's
R/Python sources and entry point are byte-identical to that revision. The
ignored `R/00_00_candidate_profile.R` had silently selected the September 18
schema-3 candidate: its browser showed 8,568 selected islets, zero trajectory
points, and no UMAP panels. That profile and the candidate-only helpers/tests
are now outside Shiny's auto-loaded directories. Schema 4 remains separate;
the owner requires exact preservation of original islet identities and values
before expansion, and revision/evaluation of the candidate contract.

The three canonical app artifacts remain unchanged and contain 15 donors.
Independent backup copies and hashes are in
`data/_pre_restore_backup/20260923/manifest.json`. Original-source aggregation
reproduces all 5,214 islets' marker means, latent means, phenotype counts and
proportions, morphology, cell counts and donor metadata exactly. The original
enriched-H5AD builder reproduces all shared values; the April enriched file
predates the August split columns. The original Barlow builder reproduces all
5,214 × 90 CSV values. No rebuilding or `.bak` substitution was necessary.

The read-only `scripts/tests/audit_frozen_trajectory.py` verifies the documented
H5AD statistics and immutable file hashes. Re-measured donor/status eta-squared
are 0.09128898/0.04445529; core INS/GCG/SST correlations are
−0.570334/0.143115/−0.205998. The documented ND/Aab+ p-values use an exact
**one-sided** donor-level Mann–Whitney test; their two-sided values are twice
as large and also nonsignificant.

**Contradiction to the historical metadata-fix claim in section 6:** the
preserved trajectory file and its August backup both still store `DC1`/`DC2`
from `X_diffmap` columns 0/1, not 1/2. This is a pre-existing documentation/data
discrepancy, not verified evidence that the on-disk fix was applied. Neither
the restored live composition axis nor its UMAP/heatmaps consumes these
columns. They were preserved rather than silently repaired or refitted.

The September 2025 Excel fallback also differs from current original TSVs in
historical `islet_union` target rows; core/band rows and its other three sheets
match. Rebuilding would discard 13,140 saved union rows and change some union
measurements. The external donor-key workbook's primary columns now describe
a different cohort, so a blind Excel rebuild is unsafe for original metadata.

Live verification on `127.0.0.1:7777` passes all nine browser checks, including
all six rendered trajectory outputs and marker/ratio switches that leave the
axis unchanged. Raw composition-pseudotime results reproduce section 11:

| Core-cell range | Islets | Median raw pseudotime |
|---|---:|---:|
| 10–1902 | 2,310 | 0.2524542212 |
| 25–1902 | 908 | 0.3825007379 |
| 10–200 | 2,181 | 0.2385783792 |

The median tolerance is half the last documented decimal (`0.00005`). The
six export checks compare UMAP coordinates and INS/GCG values directly to
the frozen H5AD, and ratio values to original core-plus-peri counts divided
by floored endocrine counts; maximum errors are below `5e-14`. Independent
heatmap bin means/counts and clamped z-scores agree within `5e-15`.
The frozen-data audit passes 46 checks; testthat passes 192 assertions and
the relevant Python trajectory/data-loading suite passes 52 tests. Browser
and export harnesses are `deploy/browser/check-frozen-trajectory.mjs` and
`scripts/tests/audit_frozen_trajectory_exports.py`. Local reports, exports,
screenshots and test logs are under `tmp/restoration-20260923/`.

The ignored `www/local_images` symlink also pointed to the September 18
candidate. Its old target is quarantined alongside the profile. The local
static image route now uses the frozen `release-5809be4/data/images` bundle
under `tmp/ufrc-release-20260908/artifacts/`, containing exactly the original
15 donor images. These workstation routes are local configuration; no public
service deployment or candidate promotion is part of this checkpoint.
