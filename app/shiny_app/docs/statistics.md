# Statistics Tab — Detailed Reference

Detailed reference for the Statistics tab. Linked from `app/shiny_app/CLAUDE.md`.

## Statistics Tab (Phase 6+14+16, Feb-Mar 2026)

### Streamlined Redesign — De-Risk P-Hacking (Jul 2026)

The stats engine (donor-level aggregation, mixed-effects primary test, BH correction,
effect sizes) was already sound; the redesign targets the **UI's researcher-degrees-of-freedom
surface**. Changes, all in `mod_statistics_ui.R` + `mod_statistics_server.R`:

- **Fixed α = 0.05** — the `alpha` selectInput is removed. `stats_run()` sets `alpha_num <- 0.05`.
  A user-movable threshold is a goalpost-moving lever.
- **Auto-selected test type** — the Parametric/Non-parametric `test_type` radio is removed.
  `stats_run()` derives `is_parametric` from the donor-level Shapiro-Wilk result
  (`all(norm_result$is_normal, na.rm=TRUE)`; all-NA → non-parametric). This is the single source
  of truth; the Section 1 normality readout reports the same decision, and the affected
  sensitivity test is auto-selected (never flipped by hand to chase significance).
- **Fixed bin width = 50 µm** — the `bin_width` slider is removed; `STATS_BIN_WIDTH_UM <- 50L`
  (top of `mod_statistics_server.R`). Per-bin analysis is always computed. A tunable width is a
  multiple-comparison surface.
- **Demoted data-subset knobs** — outliers (`stats_remove_outliers` + `stats_outlier_z`),
  zero-exclusion (`stats_exclude_zero`), min-cells (`stats_min_cells`), and diameter range
  (`stats_diam_range`) moved into a collapsed **Advanced / sensitivity settings** `<details>`
  block in Section 1. Same input IDs, so `stats_run()` reads them unchanged. Opening the block is
  a deliberate, off-path act; every setting is recorded in the CSV header.
- **Single primary result** — Section 2 shows one headline (`primary_headline`: mixed-effects
  omnibus p + ICC + the auto-selected sensitivity-test name) and the effect-size forest plot. The
  full hypothesis table (`test_results_table`) + donor-level tests are collapsed under a
  **Sensitivity check** `<details>`. The mixed-effects p-badge was removed from the Section-1
  `overview_banner` (now descriptive only).
- **Consolidated exploratory section** — AUC moved out of Section 2 into Section 3 alongside the
  per-bin heatmap + Kendall trend, under a prominent "exploratory, not confirmatory" caption and
  a client-side `stats_show_exploratory` display toggle.
- **Static section numbering (1–5)** — `section_offset()` and the two dynamic heading renderUIs
  (`section_confounders_heading`, `section_methods_heading`) plus the server-side `section_heading`
  helper were deleted; all five headings are static in `mod_statistics_ui.R`.
- **Honest correction scope** — the Methods text no longer claims blanket BH coverage; it states
  BH is applied *within* each family (donor-level pairwise, per-bin, sex-stratified) and *not*
  across families or across features, and flags Sections 3–4 as hypothesis-generating. A green
  "Design safeguards against p-hacking" paragraph heads the Methods body.

### Pseudoreplication Fix (Phase 16)

The app has 5,214 islets from only **15 donors** (5 per group). Islets within a donor are correlated. Phase 16 fixes this textbook pseudoreplication by:

1. **Donor-level aggregation**: All tests use `aggregate_to_donor()` to compute per-donor means (N=15), not islet-level (N=5,214). Produces realistic p-values (0.01-0.05 range, not 1e-50).
2. **Mixed-effects sensitivity**: `lmer_test_donor()` runs `lmerTest::lmer(value ~ donor_status + (1 | case_id))` on islet-level data with donor as random intercept. ICC quantifies donor-level clustering.
3. **Min-cells filter**: Available on all 4 analysis tabs — Plot, Trajectory, Statistics, Spatial. All use `total_cells_core + total_cells_peri >= threshold` (Trajectory uses `total_cells` directly). Defaults to 1.
4. **Normality testing**: Shapiro-Wilk on donor-level means per group, with auto-suggestion for test type.
5. **Exploratory display toggle** (superseded Jul 2026): the former `stats_no_binning` checkbox is replaced by a client-side `stats_show_exploratory` toggle that hides/shows the Section 3 exploratory views. Per-bin analysis is now always computed.
6. **Fixed bin width** (Jul 2026): the per-bin width is fixed at 50 µm (`STATS_BIN_WIDTH_UM`); the former `bin_width` slider was removed (a tunable width is a multiple-comparison p-hacking surface). Heatmap x-axis labels angled (-45°).
7. **Per-bin donor-level**: `per_bin_donor_anova()` and `per_bin_donor_kendall()` aggregate to donor means within each bin. Bins with <2 donors per group are skipped.
8. **Demographics redesign**: Age scatter shows all islets (not donor means) with Pearson r. Sex analysis uses box plot faceted by sex. Autoantibody analysis (Aab+ only): per-donor AAb profile table + feature by AAb count box plot. Use 'Sex' not 'Gender' in user-facing labels.

### Key Utilities (`utils_stats.R`)
- `aggregate_to_donor(rdf, agg_fn)` -- Groups by `Case ID` + `donor_status`, computes mean/median of `value`, preserves age/gender
- `lmer_test_donor(rdf)` -- Mixed-effects model, returns `list(p_value, icc, fit, error_msg)`. Uses `lmerTest::anova()` (Satterthwaite F-test). Returns diagnostic `error_msg` on failure (singular fit, convergence, missing packages).
- `per_bin_donor_anova(df, ...)` -- Donor-level per-bin ANOVA. Skips bins unless ≥2 donor-status groups each have ≥2 donors (shared coverage rule).
- `per_bin_donor_kendall(df, ...)` -- Donor-level per-bin Kendall tau. Same coverage rule as `per_bin_donor_anova()` so the heatmap and trend plot include the same bins (May 2026 harmonization).
- `code_donor_status(x)` -- Maps `ND/Aab+/T1D` → `0/1/2` for Kendall trend coding. Single source of truth, used by both `per_bin_kendall()` and `per_bin_donor_kendall()`.
- `normality_tests(ddf)` -- Shapiro-Wilk per group, returns data.frame with `W`, `p_value`, `is_normal`
- `cohens_d(x, y)` -- returns `list(d, ci_lo, ci_hi)` (NOT `ci_lower`/`ci_upper`)
- `eta_squared(fit)` -- from `anova()` output
- `flag_outliers(df, value_col, group_col, method, threshold, transform, exclude_zeros, enabled)` -- Single source of truth for outlier flagging. Appends `is_outlier` (logical) and `outlier_score` (z or fence-distance) columns. Attributes carry method/threshold/enabled for downstream summaries.
- `summarize_outlier_filter(df)` -- Compact summary list for badges: `n_total`, `n_excluded`, `n_kept`, `method`, `threshold`, `enabled`, `label`.
- `outlier_title_suffix(summary)` -- Formats the `" — N outliers excluded (label)"` string for plot titles. Returns `""` when no exclusion.
- `GLOBAL_OUTLIER_TOOLTIP` -- String constant used as `title=` attribute on every outlier control across tabs.
- `strip_plotly_color_leak(gg)` -- Scrubs each trace's `text`/`hovertemplate`/`hovertext` to remove `colour:\s*black` (and `color:\s*black`) lines, plus collapses adjacent duplicate `donor_status:\s*<value>` lines. Apply after every `ggplotly()` call where any geom uses `color = "black"` as a parameter (currently: Plot scatter, Plot distribution panel, Trajectory scatter). Defensive — no-op when no leak is present.

### Unified Outlier Handling (May 2026 — Statistics decoupled May 14; defaults OFF + dynamic z May 15)
- One non-namespaced checkbox + slider pair per tab (`sidebar_remove_outliers` + `sidebar_outlier_z`, `traj_remove_outliers` + `traj_outlier_z`, `spatial_remove_outliers` + `spatial_outlier_z`) — synced bidirectionally via six `observeEvent`s in `app.R` (three for the checkboxes, three for the sliders). Reactives `remove_outliers_global()` and `outlier_z_global()` are injected into `plot_server`, `trajectory_server`, and `spatial_server` as the `remove_outliers` and `outlier_threshold` arguments.
- **Statistics tab is independent**: namespaced controls `stats_remove_outliers` and `stats_outlier_z` live inside the Section 1 *Advanced / sensitivity settings* `<details>` block (Jul 2026 redesign; formerly a top-level control row) and are read via `input$…` inside `stats_run()`. They do NOT sync with the Plot/Trajectory/Spatial toggles. The `statistics_server()` function still accepts the `remove_outliers` argument for backward compatibility but ignores it.
- Defaults: **OFF** for all four checkboxes (raw distribution renders out of the box; users opt in to outlier exclusion). z-threshold slider default 3.0, range 0.5–10, step 0.5, on every tab. Toggling the Plot sidebar propagates to Trajectory and Spatial but leaves Statistics alone (and vice versa); same for the slider.
- Plot, Trajectory, Statistics, Spatial Card C: `flag_outliers(method = "zscore", threshold = outlier_threshold())` per donor status — threshold pulled from the live slider. Card C transforms via `log1p` to handle the long-tailed distance distribution.
- **Robust two-sided z (May 15)**: the `method = "zscore"` branch uses the Iglewicz–Hoaglin modified z-score, `0.6745·(x − median)/MAD`, instead of the classical mean/SD z. Because every metric in this app is non-negative and right-skewed, mean/SD pushed the achievable lower z toward 0 and the rule only flagged the upper tail in practice. Median/MAD is symmetric and genuinely two-sided. Fallback: when `MAD = 0` (zero-inflated groups), the function reverts per-group to classical mean/SD. The threshold slider scale is unchanged; the 0.6745 rescale keeps `|z| > 3` roughly equivalent to the old classical cutoff. The active rule is exposed on `attr(df, "outlier_statistic")` (`"modified_z"`, `"classical_z"`, or `"modified_z (classical_z fallback used for some groups)"`).
- Spatial Card A: `flag_outliers(method = "iqr", threshold = 1.5)` per donor status. Uses the Tukey fence rather than the modified z because enrichment z is already a standardised, roughly symmetric quantity; intentionally NOT routed through the slider (different scale). It previously passed `exclude_zeros = TRUE` on the belief that enrichment is zero-inflated — it isn't. `enrich_z` is a Poisson z, `(observed − expected)/sqrt(expected)`, so an absent phenotype yields `−sqrt(expected)`, never exactly 0; measured across all 5,214 islets every `enrich_z_*` column has zero exact zeros, making the flag a no-op. The real zero-inflation is in `peri_prop_*`, which this card does not plot.
- Visual surfacing: every affected plot shows a **count badge** in its title or subtitle and renders excluded points as **shape=4 grey × ghost markers** (Plot, Trajectory, Spatial Card C). Stats forest plot title shows the badge; the Statistics outlier checkbox carries a small caption beneath it reporting `N excluded last Run` (or a yellow "Filter OFF" notice). Statistics Methods Reference paragraph names the local rule and cites the exact `n_excluded` for the current Run. Download CSV header records `Outlier rule: …` and the count.
- Spatial Card B's "Clip z > 5" checkbox is a **display saturation only** (caps the colour scale) — left alone, NOT an outlier control. Documented as such in `docs/user_guide.md`.
- Trajectory bug fix: previously discarded outliers silently before plotting and logged only to console; now drawn as ghosts and counted in a badge above the scatter.
- **Spatial per-card outlier tables** (May 13 follow-up): Cards A and C now expose `infiltration_show_outlier_table` and `distance_show_outlier_table` checkboxes that surface the list of excluded islets in a yellow warning panel below each plot. Reactive state lives in module-scope `infiltration_outliers` / `distance_outliers` `reactiveVal`s, populated inside each renderPlotly after the `flag_outliers()` call. Renders three states: table, "no outliers detected" blue note, or filter-OFF amber label. Matches the Plot/Trajectory `show_outlier_table` UX exactly.
- **Plot distribution panel independence** (May 13 follow-up): `output$dist` now reads `raw_df_base()` (not `raw_df()`) and applies its own `exclude_zero_dist` filter, so left and right panels' "Exclude zero values" controls no longer cascade. The dist panel also re-runs `flag_outliers()` on its own subset and renders ghost × markers + a title-suffix badge, fixing the gap where outliers only appeared on the left panel. Both panels' `position_jitter()` calls now use explicit `seed = 17` (normal points) / `seed = 31` (ghost ×) so visuals stay stable across reactive re-renders.
- **`strip_plotly_color_leak` post-process** (May 13 follow-up): added because hard-coded `color = "black"` parameters on `geom_point(shape=21, …, color="black")` (Plot summary points) and `geom_smooth(…, color="black")` (Trajectory overall trend) were surfacing as a literal `colour: black` line in plotly hover tooltips, alongside a duplicate `donor_status` entry from ggplotly emitting both `fill` and inherited `colour` aesthetics. Applied at `mod_plot_server.R:684`, `mod_plot_server.R:891`, and `mod_trajectory_server.R:662`.

### R Package Dependencies for Phase 16
- `lmerTest` (3.2-1) -- Mixed-effects model with Satterthwaite df
- `lme4` (1.1.37) -- `VarCorr()` for ICC computation

### Shared Sidebar Architecture
The Plot sidebar (mode, feature, region, donor status, AAb, age, sex filters) is visible on both the Plot and Statistics tabs. Achieved via `conditionalPanel` condition `"input.tabs == 'Plot' || input.tabs == 'Statistics'"` and matching JS `adjustLayout()` logic in `app.R`. The Statistics module consumes `plot_returns$raw_df_prezero` only. **It no longer takes `summary_df`** (Jul 2026) — see "AUC is computed on this tab's own data" below.

**The sidebar's own filters run upstream and its controls are hidden here.** `min_cells`, `stat`, `binwidth`, and the outlier controls sit inside `conditionalPanel("input.tabs != 'Statistics'")`, but a `conditionalPanel` only hides the DOM — the inputs keep their values and keep filtering `raw_df_base()`. Setting the Plot sidebar's min-cells to 50 changes the Statistics run from 5,212 to 1,900 islets (lmer p 0.0039 → 0.0027) with nothing on the Statistics page to say so, and the sex filter can silently take the cohort from 15 donors to 7. `stats_run()` therefore captures `get_selection_description()` as `plot_subset_desc`; the overview banner renders it in a **"Data subset inherited from the Plot sidebar"** `<details>`, and the CSV header prints it under an explicit `# DATA SUBSET FROM THE PLOT SIDEBAR` heading. Corollary: the Statistics diameter slider is applied *on top of* the Plot sidebar's, so it can only ever narrow — widening it past the sidebar range does nothing.

### Zero-Value Exclusion (May 2026)
- UI: `stats_exclude_zero` checkbox in Section 1 controls row (`mod_statistics_ui.R`). Default OFF.
- Server: `stats_run()` filters `value != 0` immediately after outlier flagging and before min-cells / diameter / aggregation. State and drop count (`zero_dropped`) propagate through the result list.
- Scope: applies to mixed-effects model, donor-level ANOVA / Kruskal, pairwise Cohen's d, per-bin tests, demographics regressions, and the covariate-adjusted model. Alters AUC too (recomputed on the same subset since Jul 2026). Does NOT alter the normality table (computed before filtering would change power conclusions for a comparison user's reference).
- Surfacing: overview banner shows an orange "Zero-valued excluded (N dropped)" badge; forest plot title appends `— N zero-valued dropped`; Methods Reference inserts a dedicated paragraph; CSV header records `Zero-value filter: ON/OFF; N dropped = …`.
- Use case: zero-inflated features — rare phenotypes whose core (`prop_*`) or peri-islet (`peri_prop_*`) proportion is exactly 0 in most islets. Lets users sanity-check whether group differences survive on the non-zero subset. (Enrichment z-scores are *not* zero-inflated — see the Card A note above.)
- **Ordering (Jul 2026)**: the zero filter runs BEFORE `flag_outliers()`, so the robust median/MAD are computed on the values that actually enter the tests. It previously ran after, deriving the MAD from a zero-inflated vector and then discarding the zeros — a different outlier set from the Plot tab's, which has always filtered first.
- **Decoupled from the Plot tab (Jul 2026)**: the Statistics module consumes `plot_returns$raw_df_prezero` (= the Plot module's `raw_df_base`), not `raw_df`. Previously it consumed `raw_df`, which had already applied the *Plot sidebar's* `exclude_zero_top` checkbox — so a control on another tab silently altered every test while the CSV header still reported `Zero-value filter: OFF`. `stats_exclude_zero` is now the only zero filter applied to the statistics.

### Sex-Stratified Tests — BH Correction (May 2026)
- Pairwise tests within each sex are run UNADJUSTED (`p.adjust.method = "none"`) and then BH-corrected across all 6 contrasts (3 pairs × 2 sexes) in `gender_strat$p_adj`. Replaces the prior "minimum raw p, no correction" reporting.
- `gender_plot` subtitle reports minimum BH-adjusted q per sex with a `(BH-corrected across 6 contrasts)` annotation. The `low_power` flag (n=2-3 donors per sex × group) remains as a separate badge.
- CSV export includes a `# Post-hoc: Sex-Stratified Pairwise (BH-corrected across 6 contrasts: 3 pairs x 2 sexes)` section. Its column is written as **`sex`** (the in-memory data column stays `gender`).

### 5-Section Narrative Layout

A "About this tab" banner above Section 1 describes the tab's purpose, primary test (mixed-effects with donor random intercept), and how to use it. Each section heading carries a `(?)` icon with a native-tooltip one-line description (helper: `section_heading(step, title, subtitle, tip = NULL)` in `mod_statistics_ui.R`). Prior to the Jul 2026 redesign a duplicate copy also lived in `mod_statistics_server.R` to render dynamically-numbered headings; that server-side helper (and the dynamic numbering) was removed when section numbering became a static 1–5, so the helper now exists only in the UI file.

(Section names updated in the Jul 2026 redesign; see "Streamlined Redesign" above.)

1. **Configure Analysis** -- Descriptive overview banner (feature + N/donors only), Run/CSV buttons, read-only normality readout (drives auto test selection), and a collapsed **Advanced / sensitivity settings** `<details>` block holding the outlier toggle+slider, zero-exclusion, min-cells, and diameter-range knobs. Test type is auto-selected; α is fixed at 0.05 (no selectors).
2. **Primary Result** -- Pseudoreplication info banner, then a single headline card (`primary_headline`: mixed-effects omnibus p + ICC + auto-selected sensitivity-test name, col-5) + effect-size forest plot (col-7). The full hypothesis table + donor-level tests are collapsed under a **Sensitivity check** `<details>`. AUC no longer lives here.
3. **Exploratory & Size-Dependent Patterns** -- Gated by the client-side `stats_show_exploratory` toggle under an "exploratory, not confirmatory" caption. AUC (moved from Section 2) + donor-level stratified heatmap + Kendall trend (fixed 50 µm bins). Heatmap and trend share a single coverage rule (≥2 groups × ≥2 donors).
4. **Confounders & Demographics** -- Demographics: donor summary table, age scatter (islet-level, all points), sex box plot (faceted by sex, BH-corrected across 6 contrasts), autoantibody profile (Aab+ only), covariate model (donor-level).
5. **Methods Reference** -- "Design safeguards against p-hacking" summary, pseudoreplication, donor-level aggregation, ICC, mixed-effects, normality-driven auto test selection, fixed α / bin width, zero-value exclusion state, and honest correction scope.

### Multiple Testing Correction
- Per-bin ANOVA, Kendall τ, and pairwise tests are BH-corrected across bins via `p.adjust(method = "BH")`. Heatmap shows `-log10(q)` (corrected) with star annotations using corrected values. All three reach the CSV (see below) — the pairwise family previously did not, so 3 of the heatmap's 5 rows could not be reproduced from the export.

### Download CSV contents (Jul 2026)

The export is a stacked multi-section CSV: a `#` provenance/methodology header, then one `# <Section>` block per result table. **Every post-hoc family that appears on screen is now exported**; three were previously computed, plotted, and silently dropped:

| Section | Source | Was exported before? |
|---|---|---|
| `# Shapiro-Wilk Normality Tests` | `norm_result` | yes |
| `# Donor-Level Means` | `donor_df` (`gender` → **`sex`** on export) | yes |
| `# Post-hoc: Mixed-Effects Pairwise (emmeans, Tukey-adjusted)` | `lmer_pairwise` — `contrast, estimate, SE, df, t/z.ratio, p_tukey, sig` | **NO** — this is the post-hoc for the tab's own PRIMARY test; only the omnibus p and ICC survived, as header comments, so per-contrast primary results were unrecoverable |
| `# Post-hoc: Pairwise Comparisons (donor-level, BH across 3)` | `pairwise_df` | yes |
| `# Covariate-Adjusted Model` | `age_model_p` from `lm(value ~ donor_status + age + sex)` | **NO** |
| `# Per-Bin ANOVA` / `# Per-Bin Kendall` | `bin_anova` / `bin_kendall` | yes |
| `# Post-hoc: Per-Bin Pairwise (BH across bins)` | `bin_pairwise` — rows 2-4 of `bin_heatmap` | **NO** |
| `# AUC by Donor Group`, `# Demographics Summary`, `# Age-Feature Correlation`, `# Post-hoc: Sex-Stratified Pairwise`, `# Autoantibody Profile` | — | yes |

Header notes:
- `get_selection_description()` belongs to the **Plot** module, so its lines describe the *Plot sidebar*, which filters the data before this tab sees it. Both blocks used to emit a bare `Min cells per islet (core+peri): N` line, so one file could carry two contradictory values for the same named quantity (Plot 50 / Statistics 1) with nothing saying which bound. Both are now attributed: the Plot lines sit under `# DATA SUBSET FROM THE PLOT SIDEBAR (applied upstream of this tab)`, this tab's under `# [Statistics tab] …`. Only `Bin Width` is dropped — it controls the Plot scatter alone and no longer reaches any number here.
- `stats_run()` is an `eventReactive(..., ignoreInit = TRUE)`, so before Run it **raises** rather than returning NULL. `content()` wraps it in `tryCatch` — clicking Download before Run now writes the `# No results computed` fallback instead of failing the download.
- On screen, the emmeans `p.value` is Tukey-adjusted and is displayed under **`p (adj)`**; it was previously printed under `p-value` with `p (adj)` showing `-`, which read as uncorrected.
- Donor-level pairwise contrasts are BH-corrected across the 3 pairs.
- Sex-stratified pairwise contrasts are BH-corrected across all 6 contrasts (3 pairs × 2 sexes).
- **Scope caveat (Jul 2026)**: BH is applied *within* each family only — NOT across families, and NOT across the many features a user may test. The Methods text states this honestly (the earlier "BH correction applied to pairwise p-values" wording overstated coverage), and Sections 3–4 are labelled exploratory / hypothesis-generating. Confirmatory claims should come from the Section 2 primary test on a single pre-specified feature.

### Plotly Legend Positioning
Forest plot uses `theme(legend.position = "none")` in ggplot (prevents ggplotly auto-placement), then `layout(margin = list(b = 60), legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.25))` in plotly. The `margin(b=60)` creates room below the x-axis title to prevent overlap. The trend plot uses `margin(b=100)` with legend at `y = -0.3`.


## Jul 2026 correctness pass (statistics accuracy audit)

Six defects found by driving the real modules headlessly (`shiny::testServer`) against the live H5AD across every input and Advanced/sensitivity setting.

### 1. "Min cells/islet" was a silent no-op in Microenvironment mode
`prep_data()` built `targets_all` (the Microenvironment source) with a fixed `select()` that never carried `total_cells_core` / `total_cells_peri`; both the Plot sidebar filter and `stats_run()` are guarded on those column names, so the filter *skipped* rather than erroring. `Capillary|islet_core` at min-cells 1 → 200 stayed at 5,380 islets and p = 0.3476 throughout, while the exported CSV still recorded the value the user set.

**Fix:** `data_loading.R` joins the per-islet cell counts from `comp` onto `targets_all` by `(Case ID, islet_key)` — many-to-one, row count asserted unchanged. Microenvironment now matches Cell Populations islet-for-islet at every threshold (5,193 / 5,061 / 4,244 / 1,900 / 687 / 307). This also fixes the identical silent skip on the **Plot** tab. Where counts are genuinely unavailable (Excel fallback), `min_cells_applied = FALSE` drives a red banner badge and a `*** NOT APPLIED ***` marker in the CSV — never a silent skip again.

### 2. The headline p-value was always in scientific notation
`formatC(format = "e")` unconditionally, so p = 0.3476 rendered as **`p = 3.48e-01`** — one misread order of magnitude from 0.0348 — while the hypothesis table three lines below printed `0.3476`. `format_p_value()` (`utils_stats.R`) is now the single rule everywhere: fixed to 4 dp down to 1e-4, scientific below.

### 3. Mixed-effects pairwise contrasts silently switched inference method at 3,000 islets
emmeans caps its d.f. machinery at `pbkrtest.limit` / `lmerTest.limit` = 3000 and drops to **asymptotic z-tests (df = Inf)** above it — badly anti-conservative for a 15-donor design. Same feature, same model: p(ND−T1D) = **0.00007** at 5,212 islets (z) vs **0.0027** at 2,874 islets (t, df ≈ 12). The method flipped on dataset size alone with nothing on screen to say so.

**Fix:** request `lmer.df = "satterthwaite"` with both limits raised past `nrow(rdf)` — the same approximation `lmerTest::anova()` uses for the omnibus p, so the contrasts are now consistent with the headline (Kenward-Roger, emmeans' default, is far too expensive at this N). The method used is recorded as `lmer_pairwise_df_method`, surfaced in the Methods Reference and the CSV header, so a regression to df = Inf is visible rather than invisible.

### 4. The autoantibody Kruskal-Wallis was pseudoreplicated
AAb count is a **donor** attribute, but the test ran across islets: p = **1.6e-15** from 1,666 islets belonging to just **5 donors** split 1 / 3 / 1 across the count groups. Donor-level: p = **0.67**. The one test on a tab built to eliminate pseudoreplication that reintroduced it, uncaveated.

**Fix:** `aab_kw` is computed in `stats_run()` on donor means and carries `donors_per_group` + an `underpowered` flag (any group with < 2 donors). The subtitle now states the test level, the donors behind each box, and a descriptive-only warning; the result is exported (it never was).

### 5. AUC ignored every Statistics filter and obeyed hidden Plot controls
It integrated the Plot module's `summary_df()`, so it: ignored all four Advanced filters (with 5,212 islets cut to 131, AUC was byte-identical to the unfiltered run, and "T1D shows X% decrease" described data the rest of the page had discarded); obeyed the Plot `binwidth` slider (ND 9,639 → 10,483 → 9,124 across 10/50/75 µm) while Section 3 advertises a fixed 50 µm width; and obeyed the Plot Mean/Median radio, integrating medians while the card said "mean" (T1D 465 → 186).

**Fix:** AUC derives from `rdf_binned` — this tab's filtered frame at `STATS_BIN_WIDTH_UM`. `summary_df` is removed from the module signature and from the `app.R` call site so the coupling cannot come back.

### 6. Clearing the min-cells box crashed the whole run
`%||%` only catches NULL; an emptied `numericInput` sends `NA`, so `if (NA > 1 && TRUE)` aborted the `eventReactive` with *"missing value where TRUE/FALSE needed"*, blanking every card. Now coerced to a finite value ≥ 1.

### Verified unchanged
lmer Satterthwaite p and ICC match independent `lmerTest` / `VarCorr` to 10 significant figures; donor-level ANOVA / Kruskal-Wallis match; BH correction is correct and NA-safe (R's lazy `n = length(p)` excludes NAs, so 2-group subsets are not over-corrected); Cohen's d and CIs correct; zero-exclusion still runs before `flag_outliers()`; `eventReactive` isolation holds — only the Run button recomputes.

### Known, still open
- `eta_sq` for the Kruskal-Wallis branch is **epsilon-squared** (`H/(n−1)`), exported under the name `eta_sq`.
- Sex-stratified tests use pooled-SD `pairwise.t.test` while Section 2 uses Welch `t.test` (Males ND-vs-T1D: 0.104 vs 0.088); contrast labels are also reversed (`Aab+ vs ND` vs `ND vs Aab+`).
- `per_bin_donor_pairwise()` has no coverage rule beyond n ≥ 2 per pair, so the heatmap can show a pairwise star in a bin whose ANOVA/Kendall cells are blank — contradicting the "same rule for the heatmap and the trend plot" claim in the methods text.
- An over-restrictive filter returns `NULL`, and every card reverts to the pre-run placeholder — indistinguishable from "you haven't clicked Run".
- The covariate card prints "does NOT survive adjustment for age and sex" even when donor status was never significant (Capillary: 0.348 → 0.387).
- Raw column encodings leak into the UI as the feature name (`prop_Beta cell`, `peri_prop_CD8a_Tcell`, `Capillary|islet_core`).
- `outlier_summary` is measured before the min-cells / diameter filters, so "N islets excluded last Run" is counted against a superset of what was tested.

## Jul 2026 — group-comparability pass (islet number & size distribution)

Follow-up to the correctness pass above, answering: *the tests resist pseudoreplication, but do they account for differences in islet number and distribution between disease groups?* Measured, not assumed.

### What was already sound (leave it alone)

- **Islet-count imbalance is neutralised by the random intercept.** Raw counts run 121–981 per donor (8.1:1) and 2,300 / 1,666 / 1,246 per group. But effective independent observations saturate at ~1/ICC: `n_eff = n / (1 + (n-1)·ICC)`. At the observed ICC = 0.139 that is **6.87–7.17** per donor — an effective ratio of **1.04:1**, and effective N per group of 35.3 / 35.4 / 35.1. Even at ICC = 0.03 the effective ratio is only 1.23:1. Donor 6356's 981 islets do *not* dominate ND.
- **Donor counts are balanced 5/5/5**, so the donor-level sensitivity test is balanced by construction.
- **No detectable differential islet sampling**: islets per donor does not differ by status (Kruskal-Wallis p = 0.44; low power at n=5/group, and there is no imaged-tissue-area column to normalise by).
- **Residual heteroscedasticity does not move the tests.** Residual variance is ND 378 / Aab+ 318 / T1D 7.4 (51:1; Levene-equivalent F = 721, p < 2.2e-16). Refitting with `nlme::lme(..., weights = varIdent(~1|donor_status))` moves the omnibus 0.0039 → 0.0036 and the contrasts 0.0031 → 0.0029 and 0.0493 → 0.0470. Negligible, because contrast SEs are driven by between-donor variance, not residual. Donor-level pairwise already use Welch `t.test`. **No variance-structure change was made** — it would add a dependency and a fit for no measurable benefit.

### What was fixed

**1. Effect sizes were pooled across a 50:1 variance gap.** For prop_Beta cell the donor-level SDs are ND 4.06 / Aab+ 10.75 / T1D 0.575. Cohen's d pools these into 2.90, giving ND vs T1D d = 6.15 with CI [3.18, 9.12] — while Glass's Δ is **4.39** against the ND SD and **31.03** against T1D's. d = 6.15 was an arbitrary point in that range, plotted with a CI implying precision it does not have, on the one figure Section 2 foregrounds *over* the p-value.

- `cohens_d()` (`utils_stats.R`) now also returns `raw_diff` + a **Welch CI**, per-group means/SDs, both Glass's Δ, and `var_ratio` / `var_ratio_warn` (ratio > 4).
- The **forest plot now shows the unstandardised difference** in the feature's own units with the Welch CI. d and the Glass range moved into the Sensitivity-check table; contrasts with >4× variance ratio are marked `†`.
- The **primary card gained a per-group donors / islets / mean / SD table**, plus an amber note when spreads differ more than 4-fold.
- Worth knowing what this surfaces: for Aab+ vs T1D the Wilcoxon `p_adj` = 0.032 while the Welch CI on the mean difference is [−2.08, 24.6] — it crosses zero. A rank test and a mean-difference CI genuinely disagree at n=5 with a floor effect, and the table now shows both rather than only the one that reads as significant.

**2. Islet size: strong predictor, group-different distribution, and an unreported interaction.**

- Diameter distributions differ (ND vs T1D KS D = 0.091, p = 2.9e-06; medians 29.6 / 28.8 / 31.3 µm).
- Diameter is a strong predictor (p = 1.6e-141 for prop_Beta cell).
- Adding diameter as a covariate barely moves the omnibus — 0.0039 → 0.0049 (β-cell), 0.3476 → 0.3479 (Capillary). **The size shift is not confounding the group comparison.**
- But `donor_status × diameter` is **1.2e-93** for prop_Beta cell, 1.2e-50 for Ins_any, 1.1e-08 for Nerve|islet_band — and non-significant for Capillary (0.73), Lymphatic (0.36), peri CD8a (0.30). The slope of feature-on-size genuinely differs by status for some features. Where it does, the single omnibus p is a *marginal average over a size mix that itself differs by group* — an estimand issue, not bias.

`size_companion_models()` (`utils_stats.R`) fits `value ~ donor_status + diam_c + (1|donor)` and `value ~ donor_status * diam_c + (1|donor)` on the same filtered data, with diameter **mean-centred** so the Type III `donor_status` term is the difference at the mean islet diameter rather than at diameter zero. Results appear in the primary card, as two indented rows in the hypothesis table, and in the CSV header. When the interaction is significant a red box says the group difference depends on islet size and that Section 3 is part of the primary reading. Cost: 3 `lmer` fits instead of 1 — full `stats_run()` is ~3.1 s.

**3. Per-bin coverage was group-dependent and invisible.**

```
bin   donors ND/Aab+/T1D   islets ND/Aab+/T1D
 25       5 / 5 / 5         1915 / 1346 / 1021
175       3 / 5 / 5           12 /   35 /   22
225       2 / 5 / 3            4 /   14 /    6   ← tested, previously same visual weight as bin 25
275       1 / 3 / 2            1 /    4 /    3   ← no omnibus row, but pairwise still shown
```

- Axis tick labels carry each bin's donor / islet totals.
- Cell hover shows the per-group coverage **behind that specific cell** — for a pairwise row, only the two groups in that contrast.
- New `bin_coverage_note` output under the heatmap names every bin without an omnibus row, with the reason and its actual status.

That note exposed a real inconsistency: `per_bin_donor_anova()` required ≥2 groups present **and every present group** to have ≥2 donors, while `per_bin_donor_pairwise()` only required ≥2 donors in the two groups being compared. Bin 275 (ND = 1 donor) therefore got no ANOVA/Kendall row but still drew an `Aab+ vs T1D` pairwise star — a post-hoc with no omnibus behind it. **Resolved (see below).**

### Still open

- The donor-level estimand is still undeclared: unweighted donor means give p = 0.00374, precision-weighting by islet count gives p = 0.00855. Both defensible ("average donor" vs "average islet").
- Everything in the previous section's "Known, still open" list.

### Gotcha for future edits

`isTRUE(x)` on a **vector** collapses to a single `FALSE` — it is not vectorised. Both the `†` table marker and the forest-plot title suffix were written with `isTRUE(pw$var_ratio_warn)` and silently never fired. Use `x %in% TRUE` for element-wise NA-safe truth. Also: `renderTable(align = ...)` needs exactly one character per column; adding two effect-size columns without widening `"lccccc"` throws at render time, not at parse time.

## Jul 2026 — per-bin coverage rules unified (no pairwise without an omnibus)

Resolution of the inconsistency documented above, per an explicit call: *no pairwise should take place without the ANOVA.*

`bin_coverage_ok(group_values)` in `utils_stats.R` is now the **single definition** of "this bin is analysable" — ≥ 2 donor-status groups present AND every group present in the bin having ≥ 2 donors. All three per-bin tests gate on it:

| function | before | after |
|---|---|---|
| `per_bin_donor_anova()` | inline `nrow(grp_n) < 2 \|\| any(grp_n$n < 2)` | `bin_coverage_ok()` |
| `per_bin_donor_kendall()` | equivalent rule re-derived on the status *codes* | `bin_coverage_ok()` |
| `per_bin_donor_pairwise()` | **no bin rule** — only `length(g1) >= 2 && length(g2) >= 2` per pair | `bin_coverage_ok()` gates the whole bin before any pair is tested |

A pairwise post-hoc is therefore never reported for a bin with no omnibus behind it.

**Verified** on prop_Beta cell / Ins_any / prop_Macrophage: the three families now cover an identical bin set (25, 75, 125, 175, 225), and the orphan bin-275 `Aab+ vs T1D` contrast is gone.

**Consequence for the BH family** — the per-bin pairwise family shrank 16 → 15 tests, so every `p_adj` in it moved slightly *down*:

```
bin  contrast          p        q (before, n=16)   q (after, n=15)
 25  ND vs T1D         0.0112      0.0382             0.0358
 25  Aab+ vs T1D       0.0345      0.0656             0.0615
125  Aab+ vs T1D       0.0119      0.0382             0.0358
```

No conclusion flipped at α = 0.05, but the numbers do change — expected, and the right direction (fewer tests in the family).

**Drift assertion.** The `bin_dropped$status` field still checks whether any pairwise contrast survived in a bin without an omnibus row. It should now always read `"not tested at all"`; if it ever reports `"BUG: pairwise contrasts survived without an omnibus …"`, the three gates have diverged again. Cheap self-check — keep it.

**Not implemented:** the stricter *protected* variant, where pairwise runs only if the bin's omnibus is **significant**. What is implemented gates on the omnibus being *computable*, which is what was asked for. The protected form would further shrink the BH family and is a one-line change in `per_bin_donor_pairwise()` if wanted.

## Jul 2026 — Section 3 layout: text clipped in the col-4 cards

Reported as "text in the Stratified Tests by Islet Diameter card is cut off on top and under the x-axis labels". Diagnosed by measuring the rendered DOM through chromote rather than eyeballing, which is the only way the three separate causes separated cleanly.

**The binding constraint is card width.** Section 3's three cards are `column(4, ...)`, ~**401 px** at a 1600 px viewport. Everything below follows from that.

Measured before the fix (`getBoundingClientRect`, relative to the 520 px plot):

```
title_w 453 in a 401 px plot, left -26 / right +427   → clipped at BOTH ends
xtick_h 76 (two lines, rotated -45°)
xtitle_bottom 543 vs plot height 520                  → axis title 23 px off-canvas
```

Three causes, three fixes:

1. **In-plot titles that don't fit.** `"Stratified Significance (BH-corrected, donor-level)"` measured 453 px at 18 pt. Plotly does not wrap titles — it clips them. **Removed**; the card's own `h5` already names the plot and the caption already states the BH correction. The AUC card had the same defect (`"Integrated Area by Donor Group"`, 340 px, overflowing 21 px) and its title was removed too.
2. **No explicit `margin`.** `automargin = TRUE` does not measure HTML-styled multi-line tick labels reliably and never reserves room for a title. Both plots now set `margin` explicitly (`b = 105` heatmap, `b = 95` AUC). Tune by measurement — `xtitle_overflow_y` should be 0.
3. **Donor/islet counts stacked under each tick.** This was the Jul 2026 coverage work's idea and it does not survive 401 px: adjacent bins' second lines collide into an unreadable run (`15d / 4,4285 / 216i / 168i / 24i`) and the 76 px tick block pushed the axis title off-canvas and over the note below. **Reverted to plain bin midpoints.** The counts moved to a collapsible "Donors / islets per bin" strip above the coverage note, which has the full column width, and remain in the per-cell hover.

Two follow-on traps hit while fixing it:

- **`type = "category"` is required on the heatmap x-axis.** Once the labels became plain `"25"`, `"75"`, … plotly auto-detected a *linear* axis and drew its own ticks (0, 100, 200) that no longer corresponded to the heatmap columns. The HTML in the old labels had been masking this by making them non-numeric.
- **`tickangle = -45` was never the problem** — rotating a *two-line* label was. With single-line labels the rotation is required: five 3-digit midpoints collide flat (`125175225`).

Final state, all three Section 3 plots: `xtitle_overflow_y = 0`, no in-plot titles, ticks one-per-bin and legible.

### Method note

`chromote` measurement beats screenshots for this, but check both — the geometry said "no overflow" while the picture still showed the illegible colliding tick labels, and only the picture caught it. Working recipe:

```r
b <- ChromoteSession$new(width = 1600, height = 1400)
b$Page$navigate(url); Sys.sleep(12)
# click the tab by anchor text, click Run Statistics, wait ~30 s for the lmer fits
b$Runtime$evaluate("<JS returning JSON.stringify(getBoundingClientRect deltas)>")$result$value
```

`b$screenshot(selector=)` proved unreliable here (wrote nothing, and `region=` takes a keyword not a rect in this chromote version). Use raw CDP instead, and note it races `b$close()`:

```r
png <- b$Page$captureScreenshot(format = "png", clip = list(x=, y=, width=, height=, scale=1))
writeBin(jsonlite::base64_dec(png$data), out)
```
