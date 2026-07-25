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
The Plot sidebar (mode, feature, region, donor status, AAb, age, sex filters) is visible on both the Plot and Statistics tabs. Achieved via `conditionalPanel` condition `"input.tabs == 'Plot' || input.tabs == 'Statistics'"` and matching JS `adjustLayout()` logic in `app.R`. The Statistics module consumes `plot_returns$raw_df` and `plot_returns$summary_df` directly -- no data duplication.

### Zero-Value Exclusion (May 2026)
- UI: `stats_exclude_zero` checkbox in Section 1 controls row (`mod_statistics_ui.R`). Default OFF.
- Server: `stats_run()` filters `value != 0` immediately after outlier flagging and before min-cells / diameter / aggregation. State and drop count (`zero_dropped`) propagate through the result list.
- Scope: applies to mixed-effects model, donor-level ANOVA / Kruskal, pairwise Cohen's d, per-bin tests, demographics regressions, and the covariate-adjusted model. Does NOT alter AUC (read from the Plot tab's `summary_df`) or the normality table (computed before filtering would change power conclusions for a comparison user's reference).
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
- `get_selection_description()` belongs to the **Plot** module, so its `Bin Width` and `Diameter Range` lines describe the *Plot sidebar*. The Statistics handler strips those two and restates the values that actually applied: the fixed `STATS_BIN_WIDTH_UM` and `input$stats_diam_range`. `stats_min_cells` is recorded too — the UI promises "every setting is recorded in the downloaded CSV header" and these three were the omissions.
- `stats_run()` is an `eventReactive(..., ignoreInit = TRUE)`, so before Run it **raises** rather than returning NULL. `content()` wraps it in `tryCatch` — clicking Download before Run now writes the `# No results computed` fallback instead of failing the download.
- On screen, the emmeans `p.value` is Tukey-adjusted and is displayed under **`p (adj)`**; it was previously printed under `p-value` with `p (adj)` showing `-`, which read as uncorrected.
- Donor-level pairwise contrasts are BH-corrected across the 3 pairs.
- Sex-stratified pairwise contrasts are BH-corrected across all 6 contrasts (3 pairs × 2 sexes).
- **Scope caveat (Jul 2026)**: BH is applied *within* each family only — NOT across families, and NOT across the many features a user may test. The Methods text states this honestly (the earlier "BH correction applied to pairwise p-values" wording overstated coverage), and Sections 3–4 are labelled exploratory / hypothesis-generating. Confirmatory claims should come from the Section 2 primary test on a single pre-specified feature.

### Plotly Legend Positioning
Forest plot uses `theme(legend.position = "none")` in ggplot (prevents ggplotly auto-placement), then `layout(margin = list(b = 60), legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.25))` in plotly. The `margin(b=60)` creates room below the x-axis title to prevent overlap. The trend plot uses `margin(b=100)` with legend at `y = -0.3`.

