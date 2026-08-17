# Handoff — Statistics correctness, nPOD metadata, Inflammation tab (2026-07-27/28)

Scope: `app/shiny_app/` plus one decision on the CODEX Phenocycler workstream. Seven commits on
`claude/spatial-cardc-per-cell-distance`, `5106c55..63f07a8`, 17 files, +3,454/−162. **All deployed**
(`touch app/shiny_app/restart.txt`; both `:3838` and the `:8080` nginx proxy verified 200).

Previous handoff archived to `archive/docs/HANDOFF_2026-07-25_shiny-correctness-pass.md` — its one
open app item (the three uncommitted Phenocycler files) is closed below.

---

## What prompted it

"Investigate whether the Statistics tab computes and presents its tests accurately for all inputs.
Min cells/islet doesn't change the p-value for Microenvironment=Capillary."

That reproduced, and the audit that followed found five more defects. The session then widened to a
second statistical question (does the framework account for islet number/size differences between
groups?), a methods comparison against a published cohort, and a new tab.

**Method used throughout: drive the real modules headlessly with `shiny::testServer` and measure.**
Most of these defects are invisible to code review because they are silent no-ops or size-dependent
method switches. See the memory note `shiny_testserver_audit_harness` for the two-stage harness.

---

## Commits

| SHA | What |
|---|---|
| `37ee4b7` | Six Statistics accuracy defects (below) |
| `d403563` | Group spread, islet-size companions, per-bin coverage |
| `3ed107e` | Pairwise gated on the omnibus; Section 3 layout clipping |
| `5923df9` | nPOD donor metadata, age-group companion, Inflammation tab |
| `6a385d7` | Infiltration score: whole-slide false negatives |
| `6925ec4` | Two open infiltration issues; track `islet_analysis/CLAUDE.md` |
| `63f07a8` | Reject the adaptive RESTORE divisor; preserve + record |

---

## 1. Six Statistics-tab defects (`37ee4b7`)

1. **"Min cells/islet" was a silent no-op in Microenvironment mode.** `prep_data()` built
   `targets_all` with a fixed `select()` that never carried `total_cells_core`/`total_cells_peri`,
   and both the Plot and Statistics filters are guarded on those column names — so the filter
   *skipped* rather than erroring. Capillary held at 5,380 islets and p = 0.3476 from min-cells 1
   through 200 while the CSV header still recorded the value set. Fixed by joining the counts onto
   `targets_all`; this also fixed the identical skip on the Plot tab.
2. **The headline p was always scientific**: 0.3476 rendered as `3.48e-01` — one misread order of
   magnitude from 0.0348, and almost certainly the origin of the report. `format_p_value()` is now
   the single rule everywhere.
3. **emmeans silently switched inference method at 3,000 islets.** Above `pbkrtest.limit` it drops
   to asymptotic z (df = Inf): p(ND−T1D) = 0.00007 at 5,212 islets vs 0.0027 at 2,874. Now requests
   Satterthwaite explicitly and records which method ran.
4. **The autoantibody Kruskal-Wallis was pseudoreplicated** — p = 1.6e-15 from 1,666 islets
   belonging to 5 donors split 1/3/1. Donor-level: p = 0.67.
5. **AUC ignored every Statistics filter and obeyed hidden Plot controls** (bin-width slider,
   Mean/Median radio). Now derived from this tab's own frame; `summary_df` removed from the module
   signature so the coupling cannot return.
6. **Clearing the min-cells box crashed the run** (`NA` from an emptied `numericInput`).

## 2. Group comparability (`d403563`)

Answering "does this account for islet number/size differences between groups?" — measured, and the
answer is partly reassuring:

- **Islet-count imbalance is already neutralised.** Raw counts run 121–981 per donor (8.1:1), but
  effective observations saturate at ~1/ICC: at ICC = 0.139 that is 6.87–7.17 per donor, an
  effective ratio of **1.04:1**. Donor counts are balanced 5/5/5.
- **Residual heteroscedasticity is present but harmless to the tests** (51:1 variance ratio;
  `varIdent` refit moves the omnibus 0.0039 → 0.0036). Contrast SEs are driven by between-donor
  variance. **No variance structure was added** — it buys nothing.
- **Effect sizes were the real casualty.** Cohen's d pools a 50:1 variance gap: d = 6.15 with
  CI [3.18, 9.12] while Glass's Δ ranges 4.39–31.03 for the same contrast. The forest plot now shows
  the **raw difference with a Welch CI**; d and the Glass range moved to the table with a `†` marker.
- **Islet size is a strong predictor whose distribution differs by group** (KS p = 2.9e-06), and
  `donor_status × diameter` is p < 2e-16 for β-cell and Ins_any. Two companion models now report
  this beside the primary p, with a red box when the interaction is significant.

## 3. Per-bin gating + layout (`3ed107e`)

`bin_coverage_ok()` is now the single definition of an analysable bin; ANOVA, Kendall and pairwise
all gate on it, **so no pairwise post-hoc runs without an omnibus** (user decision). The per-bin
pairwise BH family shrank 16 → 15; no conclusion flipped.

Section 3's cards were clipping text. Diagnosed by measuring the DOM through chromote: the binding
constraint is that these are `column(4)` cards, ~401 px. A 453 px title clipped at *both* ends (that
was the "cut off on top"), and rotated two-line tick labels pushed the axis title 23 px off-canvas.

## 4. nPOD metadata, age model, Inflammation tab (`5923df9`)

Three methods ported from **Steenbuck et al. 2025** (bioRxiv 2025.03.05.641526, Bodenmiller/Atkinson;
88 donors, 75 ROIs each). Their primary framework — LMM with donor random intercept, donor-mean
Wilcoxon fallback, FDR — is *already* what this tab does, so only the additive parts were taken.

- `scripts/build_donor_metadata.py` → `data/donor_metadata.csv` (committed; the 741-donor nPOD xlsx
  stays gitignored). Adds disease duration, BMI, HbA1c, C-peptide, AAb panel.
- **Age-group companion** `value ~ donor_status + age_group + (1|donor)`. Adjustment barely moves the
  effect and age never carries it.
- **Inflammation tab** (new, sixth): donor cohort table, edgeR differential abundance, infiltration
  score.

## 5. Infiltration score — whole-slide false negatives (`6a385d7`, `6925ec4`)

Prompted by "our data are whole-slide, theirs were single-islet ROIs". Three mechanisms, all measured:

- **Denominator** was the dominant defect. Their per-image denominator is islet-local by
  construction; ours was the whole section, and the all-cells:islet-proximal ratio ranges **4.1× to
  23.2×** across donors. Not confounded with stage (KW p = 0.53) — variance inflation, not bias.
  Switching to islet+peri: T cell 0.1845 → 0.0081, Immune 0.1119 → 0.0207, B cell 0.0652 → 0.0155,
  CD8a 0.0271 → 0.0108. **Four of seven phenotypes cross from null to significant on the denominator
  alone.**
- **Hard power ceiling.** At 5 donors/stage the best attainable KW p is 0.00193 and pairwise Wilcoxon
  0.00794. BH across 21 pairwise contrasts floors at q = 0.167 — **no pairwise contrast could ever be
  significant**. Those columns are now reported raw and labelled exploratory, with the ceiling printed
  beside the table so a null reads as underpowered rather than as no effect.
- **Pre-specified family** `INFLAM_FDR_FAMILY` (7 immune phenotypes, fixed in source, floor q = 0.014
  vs 0.041 across 21). Out-of-family phenotypes display uncorrected and marked `‡`.
- **Numerator/denominator mismatch** — 29% of numerator mass sits outside the islet+peri denominator.
  Exposed as a *Match numerator to denominator* control (off by default: matching confines the `1/d`
  gradient to the 20 µm peri band).
- **Nearest-islet attribution — checked, not a confounder.** Islet crowding and count don't differ by
  stage (p = 0.83, 0.53); score-vs-geometry |rho| ≤ 0.35, and negative.

## 6. Adaptive RESTORE divisor — REJECTED (`63f07a8`)

Closes the previous handoff's open item. The three uncommitted `codex_phenocycler` files proposed an
outlier-robust divisor. Reverted; preserved verbatim at
`scripts/codex_phenocycler/docs/rejected/adaptive_restore_divisor.patch`. Full reasoning is in that
workstream's `HANDOFF.md` (2026-07-28 entry). Short version: measured against the shipped screens, it
**does not fire on any of the donors it was written to fix** (HLADR 6521/6534/6551 at ratios
1.46/1.47/1.94; PDPN 6479/6538 at 2.12/1.00, all far below its k=5 trigger), it **does** change
markers nobody flagged (SMA fires on over half the donors), it contradicts two reviewed decisions,
and it bypasses the project's own divisor tripwire by patching CODEX-side so `METHOD_VERSION` never
moves.

---

## Archiving done this session

- `docs/HANDOFF.md` → `archive/docs/HANDOFF_2026-07-25_shiny-correctness-pass.md` (`git mv`, history
  preserved). Superseded: this session reworked the same Statistics code and closed its open item.
- `scripts/codex_phenocycler/docs/rejected/adaptive_restore_divisor.patch` — rejected work kept
  retrievable rather than discarded.

**Data: nothing was archived, deliberately.** Everything under `data/` is live:

| File | Status |
|---|---|
| `master_results.xlsx` (Sep 2025) | live Excel fallback for `load_master_auto()` |
| `adata_ins_root.h5ad`, `islet_explorer.h5ad` | live pipeline inputs |
| `islet_spatial_lookup.csv` | live — **but holds 10,736 records for 5,214 analysed islets** (~2 entries each); fine for lookups, do not use it for absolute geometry |
| `cell_distance.csv` (2.65M rows) | live; both Inflammation cards read it |
| `All-Download_From_nPOD_*.xlsx` (repo root, gitignored) | **keep** — needed to regenerate `donor_metadata.csv` |
| `CellCompartments.tsv` (repo root, 11.5 GB) | **live input**, hard-coded at `codex_restore_validation.py:81`. The previous handoff said "delete if disk pressure appears" — that is misleading; it would break the Phenocycler screen. `/home` is at 80% with 693 GB free, so no pressure. |

Two `.gitignore` defects fixed along the way: `data/` and `islet_analysis/` were directory-form
ignores, and **git never descends into an excluded directory, so a `!path` negation is inert** — which
is why two gate CSVs had to be force-added with no record of why. Both are now `dir/*` plus named
exceptions, which surfaced that `data/DATA_PROVENANCE.md` and `islet_analysis/CLAUDE.md` had **never
been tracked** despite `CLAUDE.md` citing both as canonical. Both are now in the repo.

---

## Still open

- **`app_redesign.R`** (repo root, untracked, created 2026-07-28). Not mine; never read or touched.
- **Phenocycler workstream** continues at its own handoff — Phases H–K, the CD8a/HLADR prevalence
  problem, and the container that has never been built.
- ~~**The Statistics tab's remaining known issues**~~ — **all closed 2026-08-14.** The
  `eta_sq`/epsilon-squared misnomer, the pooled-SD-vs-Welch mismatch (and its reversed contrast
  labels), the silent `NULL`→placeholder revert, the covariate-card wording, the raw column
  encodings, and the `outlier_summary` ordering are all fixed; resolutions are recorded inline in
  `app/shiny_app/docs/statistics.md` § "Known, still open". One item remains genuinely open: the
  **donor-level estimand is undeclared** (unweighted donor means p = 0.00374 vs islet-count-weighted
  p = 0.00855) — a choice to state, not a bug.
- **Trajectory**: an equivalent audit had been done but never written into the repo. It now lives at
  `app/shiny_app/docs/trajectory.md`, along with four defects fixed the same day — the Trajectory
  CSV's provenance header naming a non-existent input (so it always claimed "Min cells per islet: 1"),
  the same `NA`-from-a-cleared-`numericInput` crash the Statistics tab had, `DC1`/`DC2` storing the
  trivial λ=1 eigenvector, and `obsp` carrying the combined neighbor graph while `uns` advertised
  core. Step 6 validation is now donor-level, covers both pseudotime modes, and exits non-zero.
  **`scripts/CLAUDE.md` and `DATA_PROVENANCE.md` were stale** — they described raw scVI / n=15 /
  INS −0.69, where the shipped file is Harmony θ=0.1 / n=30 / INS −0.570 / GCG 0.143; both corrected.
- **Infiltration**: 16–20% of cells sit >200 µm from any islet (44% for donor 6548) — excluded from
  the default denominator but still in the numerator; that is the 29% mismatch quantified above.

## Verification available to repeat

```bash
# headless module tests (two-stage: drive plot_server, feed its frame to the module under test)
# COMMITTED as of 2026-08-14 — was a scratchpad harness, now a standing regression net.
# 26 checks pinning all six Jul 2026 defects; reproduces the audit's own numbers.
cd app/shiny_app && Rscript tests/manual/stats_audit.R

cd app/shiny_app
Rscript -e 'for (f in c("app.R", list.files("R","[.]R$",full.names=TRUE))) parse(f)'
Rscript -e 'shiny::runApp(".", port=7777)'   # dev only; prod is the symlink + restart.txt

~/miniconda3/envs/scvi-env/bin/python scripts/build_donor_metadata.py   # 15 rows + both warnings
cd scripts/codex_phenocycler && python -m pytest test_codex_restore_validation.py -q   # 11 pass
```

Layout regressions are caught by measuring, not by eye — chromote `getBoundingClientRect` deltas;
`b$screenshot(selector=)` is unreliable here, use raw CDP `Page.captureScreenshot` with a clip rect
and sleep before `b$close()`.

## Deployment

`/srv/shiny-server/islet-explorer` → `app/shiny_app/` (symlink). `touch app/shiny_app/restart.txt`
retires the workers; a fresh one appears within seconds. Verify both
`http://127.0.0.1:3838/islet-explorer/` and the nginx proxy on `:8080`. Worker logs are
`shiny:shiny` mode 640 — this account cannot read them without sudo.
