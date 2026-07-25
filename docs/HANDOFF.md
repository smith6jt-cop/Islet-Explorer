# Handoff — Shiny app correctness pass (2026-07-24/25)

Scope: `app/shiny_app/` only. The CODEX Phenocycler port is a separate, in-flight
workstream — see `scripts/codex_phenocycler/docs/HANDOFF.md`, untouched by this work.

**Commits:** `659210b` (app fixes, 19 files, +1014/−244) and the archive/handoff commit
that follows it, both on `claude/spatial-cardc-per-cell-distance`.
**Deployed:** yes — see [Deployment](#deployment).

---

## What prompted it

A review pass named six issues. Four turned out to be genuine correctness bugs, not
cosmetics; two of those silently changed statistical results.

## 1. Zero-exclusion corrupted the z-score three ways

`flag_outliers()` (`R/utils_stats.R`) takes `exclude_zeros`. Three defects:

- **Zeros were scored against a distribution they were excluded from.** The centre/scale
  came from the non-zero subset (correct) but the score was written to *every* row, so a
  zero got `0.6745·(0 − median_nz)/MAD_nz` — a large spurious negative that surfaced in
  the outlier tables. Zeros now leave the scoring population entirely (`outlier_score = NA`,
  never flagged).
- **The Statistics tab flagged outliers BEFORE dropping its zeros**, deriving the MAD from
  a zero-inflated vector and then discarding the zeros — a different flagged set from the
  Plot tab, which has always filtered first. Reordered to match.
- **A control on the Plot tab silently changed every statistical test.** Statistics consumed
  `plot_returns$raw_df`, which had already applied the *Plot sidebar's* `exclude_zero_top`
  checkbox — while the Statistics CSV header still reported `Zero-value filter: OFF`. It now
  consumes `raw_df_prezero` (= the Plot module's `raw_df_base`), so `stats_exclude_zero` is
  the only zero filter applied to the statistics.

**Also: Spatial Card A's `exclude_zeros = TRUE` was a no-op resting on a false premise.**
`enrich_z` is a Poisson z, `(observed − expected)/sqrt(expected)`, so a phenotype absent
from an islet's peri zone yields `−sqrt(expected)` — a negative continuum, never exactly 0.
Measured across all 5,214 islets: **every `enrich_z_*` column contains zero exact zeros.**
The real zero-inflation lives in `peri_prop_*`, which that card does not plot. Removed the
dead flag; corrected five places (code comments, the global outlier tooltip, the Statistics
methods text, `docs/statistics.md`) that asserted enrichment is zero-inflated.

Three test cases added for the `exclude_zeros` branch — it previously had none.

## 2. "CD4 Tcell" was missing from the Spatial phenotype checkboxes

The choices were derived from `sort(unique(donor_cells()$phenotype))` — the **currently
loaded donor's** CSV. Donor 112 is the default (`list.files()` sorts `112.csv` first) and
has **zero CD4 Tcell cells**; so do 6533 and 6534. The checkbox simply did not exist.

A second defect compounded it: `intersect(prev, phenos)` preserved selections against a
per-donor list, so any phenotype that *newly appeared* on a donor switch stayed permanently
unchecked — and its cells were filtered out of the scatter.

Choices now come from `canonical_phenotype_order()` (`R/spatial_helpers.R`), backed by
`PHENOTYPE_COLORS`. All 21 always present; ones absent for the selected donor render greyed
with a `(0)` count.

**Related:** the Plot "Population Measure" dropdown displayed three spellings of the same
phenotype — `CD4_Tcell`/`ECADplus`/`SMAplus` in the Peri-Islet optgroup vs `CD4 Tcell`/`ECAD+`/`SMA+`
in Core and Core+Peri. Peri labels now route through `extract_phenotype_from_column()`.

## 3. Fonts — the cause was plotly settings, not the size constants

Three in-plot plotly settings were unset on **every** plot, which is why text still read
small after the May 2026 sizing pass:

| Setting | Effect when unset |
|---|---|
| `hoverlabel.font` | Tooltips render at plotly's ~13px default regardless of `layout(font=)`. Several plots (Leiden UMAP, tissue scatter) are read *primarily* via hover. |
| `uniformtext` | For bar/heatmap traces plotly silently **shrinks** in-cell text to fit, so a 12pt `textfont` could render at 8px. |
| `automargin` | Angled tick labels get clipped against a fixed margin — which reads as "the font is too small" when the label is really cut off. |

All three now apply centrally in `plotly_islet_config()` (already piped last on every plotly
output, so zero call-site churn). `uniform_text_mode` defaults to `"hide"`; the Statistics
`bin_heatmap` passes `"show"` so a significance star is never hidden — a hidden star reads
as "not significant".

**Two silent `ggplotly()` conversion losses:**
- `labs(subtitle = ...)` is **dropped**. Three Statistics plots put their statistics there,
  so the age correlation (`r`, `p`), the sex BH-correction line, and the AAb Kruskal-Wallis
  p **had never rendered at all**. Rebuilt into the title via `plotly_title_with_sub()`.
- Faceting relocates axis titles *and* strip labels into `layout.annotations`, where no
  `layout()` overlay reaches them (~14.6px). `plotly_islet_strip_fonts()` restyles both.

**`annotate(size=)` is millimetres, not points** — the five placeholders at `size = 4/5`
were 11.4/14.2pt, *below* the base. Added `ISLET_ANNOTATE_SIZE`.

Scale raised to base 18 / `dense_min` 14; cramped containers grown (see
`app/shiny_app/CLAUDE.md` § "Plot Fonts").

## 4. "Sex" not "Gender"

The sidebar filter — shared by the Plot **and** Statistics tabs — was the last live UI
violation, labelled `"Gender"` with raw `F`/`M` values. Now `"Sex"` with Female/Male.
Two Statistics CSV column headers renamed on export. The `gender` **data column**
(H5AD `.obs`) and internal input ids are unchanged.

## 5. CSV exports

`download_spatial` **matched no plot on the page.** It re-derived its own subset from
`prepared()$comp` instead of the shared `nbr_comp()` every card consumes:

| Control | On-screen cards | CSV (before) |
|---|---|---|
| `nbr_min_cells` | applied | not applied — exported islets the cards dropped |
| `total_cells_peri > 0` | not applied | applied — dropped islets the cards keep |
| `pheno_filter` | applied | not applied |
| outlier flags | applied | not exported |
| provenance header | — | none at all |

Rewritten onto `nbr_comp()`. Six CSVs added (Spatial Cards A/B/C, Trajectory heatmap matrix
+ UMAP, Plot distribution), each consuming the plot's **own reactive** — four inline
`renderPlotly` bodies were extracted into reactives (`infiltration_data()`, `distance_data()`,
`dist_data()`, `traj_multi_heatmap_matrix()`) so an export cannot drift from its figure.

Shared `write_islet_csv()` / `islet_csv_header()` helpers; every CSV now records the filters
that determined its rows.

## 6. Statistics CSV omitted the post-hoc tests for its own primary result

Computed, displayed on screen, carried on the result list — and never written:

- **`lmer_pairwise`** — the emmeans Tukey contrasts, i.e. the post-hoc for the *pre-specified
  primary test*. Only the omnibus p and ICC survived, as header comments, so per-contrast
  primary results were unrecoverable from the export. (A live run exports ND vs T1D at
  Tukey p = 5.8e-04.)
- **`bin_pairwise`** — rows 2–4 of the on-screen `bin_heatmap`; only the ANOVA and Kendall
  rows reached the file.
- **`age_model_p`** — the covariate-adjusted model.

Adjacent fixes: the header no longer reports the *Plot sidebar's* bin width and diameter
range (it bins at the fixed `STATS_BIN_WIDTH_UM` and filters on `stats_diam_range`); it now
records `stats_min_cells` and `stats_diam_range` as the UI promises; downloading **before**
clicking Run writes the "no results" fallback instead of failing (`stats_run` is an
`eventReactive(ignoreInit = TRUE)` — it *raises* before first fire, it does not return NULL);
and the Tukey-adjusted p is displayed under `p (adj)` rather than `p-value`.

Dead `pairwise_wilcox()` removed (no app caller; only its own test).

---

## Deployment

`/srv/shiny-server/islet-explorer` is a symlink **directly to the git working tree**
(`app/shiny_app/`). Files go live the moment they are saved; only the running R workers hold
stale code. Deploy = `touch app/shiny_app/restart.txt`, then hit the URL to spawn a fresh
worker.

Pre-flight before respawning: confirm **no new package dependencies** — shiny-server runs as
the `shiny` user and cannot see a personal R library, so a new dep returns HTTP 500. This
change added none (only `base`/`ggplot2`/`plotly`/`utils`, all already in use).

Verified live on production: Sex label; `CD4 Tcell (0)` on donor 112; all 9 CSV buttons;
`hover=15` / `uniformtext=14` / `base=18` on every plotly plot; all Spatial cards rendering.

## Verification available to repeat

```bash
cd app/shiny_app
Rscript -e 'shiny::runApp(".", port = 7777)'          # dev
Rscript -e 'testthat::test_dir("tests/testthat")'      # 71 pass
```

Headless checks use `chromote`. **Gotcha:** `selectInput` renders as a *selectize* widget
whose native `<select>.options` are stripped — read `el.selectize.options`, not the DOM
options, or every dropdown looks empty. Likewise, plotly nests an inner `.plotly` div with
no data; scan `document.querySelectorAll('.plotly')` and resolve the host via
`closest('.shiny-plotly-output')`, or every plot looks unrendered.

## Known tradeoff

Card B's enrichment heatmap shows in-cell values at 14pt (was 12). With all 17 phenotypes
displayed the numbers sit close together. Two existing controls tune it — the "Show values"
checkbox and the phenotype filter. Revert `dense_min` or drop that one `textfont` back to 12
if the crowding matters more than the size.

## Not done / still open

- **Three `scripts/codex_phenocycler/` files remain uncommitted** (+78 lines): an adaptive,
  outlier-robust RESTORE divisor that falls back to a control quantile when the control max
  is contaminated — fixing the "zero call for a marker present in every image" failure. Left
  out deliberately: separate workstream, mid-review, and its acceptance gate is image-based
  pair review. Commit it on its own after that review.
- The repo-root `CellCompartments.tsv` is 11.5 GB and gitignored. Regenerated by QuPath;
  delete it locally if disk pressure appears.

## Archived this session

See `archive/README.md`. Moved with `git mv` (history preserved via `git log --follow`):
`app_original.R` + `app_modular.R` (dead entrypoints sitting in the *deployed* directory),
`shiny-Ex.R` (+ its `.Rout`; auto-generated Shiny package examples, committed by accident),
`FileReview.txt` (superseded Phase-2 summary).

Two things that *look* archivable but are not: `app/docs/logo.png` (the live app logo —
`www/logo.png` is a symlink to it) and the CODEX Phenocycler `HANDOFF.md` (active).
