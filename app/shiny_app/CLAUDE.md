# CLAUDE.md - Shiny App (Islet Explorer)

Deep details for the modular Shiny app under `app/shiny_app/`. The top-level `CLAUDE.md` covers the entrypoint architecture, root-level outputs, shared reactive state, and data loading. This file documents per-tab features and conventions.

## Spatial Neighborhood Analysis

See [`docs/spatial.md`](docs/spatial.md) for full details — peri-islet metrics, tab layout (Phase 9+13+15), neighborhood-card server architecture, KDE selector race fix, tissue scatter design.

## Single-Cell Drill-Down (Phase 8, Feb 2026)

### Per-Islet Cell Data

`extract_per_islet_cells.py` reads single-cell H5AD -> outputs `data/cells/{imageid}_Islet_{N}.csv`:
- 5,214 files, ~203 MB total
- 37 columns: `X_centroid`, `Y_centroid`, `phenotype`, `cell_region` (core/peri), `Cell Area`, `Nucleus Area`, + 31 protein markers
- File naming matches `combined_islet_id` from `islet_spatial_lookup.csv`

### Islet Viewer Panel (Segmentation + Drill-Down)

When an islet is clicked, the viewer panel shows with controls on the left sidebar:
- **Title**: `Islet_N (case_id, status, age, sex)` — donor info from `prepared()$comp`
- **Left sidebar (col-2)**: Layer toggles (Single Cells, Peri Boundary, Structures), Color-by dropdown, Phenotype Palette, Show peri-islet cells.
- **Center plot (col-6)**: Segmentation map with GeoJSON boundaries + optional single-cell overlay
- **Right sidebar (col-4)**: Cell composition bar chart + cell count table (when Single Cells enabled)
- `selected_islet()` carries `case_id`, `islet_key`, `centroid_x/y`, `donor_status`, `donor_age`, `donor_gender`

### Coordinate Alignment

Cell centroids (`X_centroid`, `Y_centroid`) are in micrometers. GeoJSON polygons are in pixels. Conversion: `x_px = X_centroid / PIXEL_SIZE_UM`. Matches `render_islet_segmentation_plot()` coordinate system.

### drilldown_helpers.R

- `PHENOTYPE_COLORS` -- 21 phenotypes with distinctive hex palette (endocrine=warm, immune=cool, structural=neutral)
- `drilldown_available()` -- checks if `data/cells/` exists with CSVs
- `load_islet_cells(imageid, islet_key)` -- CSV loader with `new.env()` caching
- `render_islet_drilldown_plot(info, cells, color_by, show_peri)` -- cells over `build_segmentation_base_plot()`. Uses `fill` aesthetic (shape 21) for phenotype/marker coloring, `alpha` for core/peri (0.9/0.4). `colour` aesthetic reserved for structure boundaries.
- `render_drilldown_summary(cells)` -- horizontal bar chart of phenotype composition

### segmentation_helpers.R Refactor

`build_segmentation_base_plot(info)` extracted from `render_islet_segmentation_plot()`. Returns ggplot with GeoJSON polygon layers + coord_sf (WITHOUT crosshairs or title). Reused by both `render_islet_segmentation_plot()` and `render_islet_drilldown_plot()`.

**Structure Legend**: combines all visible polygons into one sf with a `structure_type` column, uses `aes(colour = structure_type)` + `scale_colour_manual(name = "Structures")` with `key_glyph = "path"`. Display names: Islet, Peri Boundary, Nerve, Capillary, Lymphatic. The `colour` aesthetic is reserved for structures; drilldown cell scatter uses `fill` (shape 21) to avoid scale conflict.

### macOS Native-Click Bridge + Shared Resolver (Jun 2026)

**Problem**: On macOS, plotly's click/drag discrimination suppresses `plotly_click` entirely, so the server `observeEvent(event_data("plotly_click", ...))` never fires — clicking a point did nothing (no notification). Both Plot and Trajectory drill-down were affected (Trajectory's key→customdata→coord fallbacks all live *inside* the plotly_click observer, so a suppressed event bypassed them too). Works on Linux/Windows; macOS-only — no platform-conditional code, so the cause is the macOS browser event layer.

**Fix** — `islet_click_bridge(p, input_id)` in `00_globals.R`:
- `htmlwidgets::onRender` that listens for native DOM `mousedown`/`mouseup` (which always fire, immune to plotly's drag suppression). On a click-like gesture (movement < 6px) it reads the hovered islet's key and pushes `{key, nonce}` to a Shiny input via `Shiny.setInputValue(..., {priority:'event'})`.
- Key source: the ggplot `key` aes (`"case_id|islet_key"`), read as `pt.data.key[pt.pointNumber]` from the cached `plotly_hover` event AND from plotly's live `gd._hoverdata` at mousedown (belt-and-suspenders).
- Piped after each scatter's `event_register("plotly_click")`: Plot uses `ns("plot_native_click")`, Trajectory `ns("traj_native_click")`. The legacy `plotly_click` handler is kept for platforms where it already works.
- **`gd.__isletBridge` attach-once guard**: `onRender` re-runs on every plot re-render and would stack listeners (one click → N fires). Guard at the top: `if (gd.__isletBridge) return; gd.__isletBridge = true;`.

**Gotcha**: ggplotly surfaces the `key` aes to the browser ONLY on the **individual-point marker traces** (`gd.data[i].key`), NOT on summary-point/line traces — so the Plot tab needs "Show individual points" ON for a click to resolve a key. Verified headlessly (chromote): trace keys like `112|Islet_1`; a real native gesture fires exactly one `[PLOT CLICK] native bridge` / `[TRAJ CLICK] native bridge` and renders the panel.

**Shared resolver** — `select_islet_from_key(case_id, islet_key, prepared, selected_islet)` in `segmentation_helpers.R`: the common tail of both drill-down paths (SF/lookup guards, centroid lookup with zero-padded case-id fallback, donor-metadata join, `selected_islet(list(...))`). Called by both the legacy plotly_click observers and the new native-click observers in Plot + Trajectory.

## Statistics Tab

See [`docs/statistics.md`](docs/statistics.md) for full details — the Jul 2026 streamlined redesign (fixed α, auto-selected test type, fixed 50 µm bins, Advanced/sensitivity + Sensitivity-check collapsibles, single primary result, consolidated exploratory section), pseudoreplication fix (Phase 16), `utils_stats.R` reference, unified outlier handling, zero-value exclusion, sex-stratified BH correction, 5-section narrative layout, multiple-testing correction.

## Trajectory Tab Features

### Cell-Count-Weighted Trajectory (Phase 12, Mar 2026)

#### Problem
60.8% of islets have ≤10 cells — their aggregated scVI embeddings are noisy. DPT operates on 10-dim scVI latent space but the scatter shows single-marker expression from `.X`. An islet can be "far" in latent space (unusual across ALL 31 markers → high DPT) while having normal INS expression. Real biology, but shouldn't get equal visual weight.

#### Solution: Cell-Count-Aware Visualization
- **Point sizing**: `sqrt(total_cells)` — area proportional to cell count. Default in trajectory scatter.
- **Weighted LOESS**: `weight = log1p(total_cells)` — well-measured islets drive trends. Raw counts too skewed (median=9, max=1,902) → `log1p` compresses 1,902:1 ratio to ~11:1.
- **Hover tooltip**: Shows "Cells: N".
- **UI controls**: Point size selector: Cell Count (default) / Islet Diameter / Uniform.

#### Key Details
- `total_cells` is in `adata_ins_root.h5ad` `.obs` (all 5,214 islets, min=1, max=1,902)
- Added to `traj_data_clean()` result_df with NA guard (fills missing with 1)
- `log1p(total_cells)` for LOESS weights — prevents extreme islets from distorting trend lines
- All three trend line variants weighted: overall, by-donor (donor_id coloring), by-donor (donor_status coloring)
- **Critical**: Raw `total_cells` as LOESS weight causes wild curves — a 1,902-cell islet dominates the fit. Always use `log1p()`.

### Pseudotime Mode Selector (Phase 17, Apr 2026)

`mod_trajectory_ui.R` has a "Pseudotime mode" dropdown with two options: "Core only" (default, `dpt_pseudotime`) and "Core + Peri" (`dpt_pseudotime_combined`). `mod_trajectory_server.R::traj_with_pseudotime()` switches `obs$pseudotime` based on `input$pt_mode`. The selector controls X-axis (pseudotime) only — Y-axis features remain the 31-marker panel. See `scripts/CLAUDE.md` for the upstream pipeline that produces these two pseudotimes.

### Pseudotime Heatmap: Dominant Status × Feature Intensity (May 2026)

The strip heatmap below the trajectory scatter uses a 2-D encoding per pseudotime bin:

- **Hue** = dominant donor_status, computed as the argmax of per-status NORMALIZED density (`count_in_bin / total_for_status`). Normalization corrects for cohort imbalance (ND ~2300, Aab+ ~1668, T1D ~1246 islets) so dominance reflects enrichment relative to each status's pool, not raw count.
- **Intensity** = bin-mean of selected feature, normalized via 5–95th percentile clip, blended in `[INTENSITY_FLOOR, 1]` toward the donor color so even the dimmest bin is a clear tint (never white). `INTENSITY_FLOOR = 0.50`.
- **Renderer**: ggplot + `geom_tile` + `scale_fill_identity` (NOT plotly heatmap). Guarantees bars span [0, 1]. Plotly previously misaligned bin widths with the scatter axis margins.
- **Hover**: Shows pseudotime, dominant status, raw counts per status (n_ND, n_Aab, n_T1D), normalized density per status (%), and feature mean.

Replaces the previous "average donor_status code (0/1/2)" heatmap, which conflated cohort-level imbalance with per-bin enrichment and produced a single colorscale gradient that hid expression intensity entirely.

### Multi-Feature Trajectory Heatmap
Z-scored expression along pseudotime. Clamped [-2.5, 2.5], blue-white-red, dynamic height.

## Plot Tab Features (Phase 5, Feb 2026)

### Phenotype Composition Explorer
Plot tab Cell Populations mode exposes 21 cell-type proportions (`prop_*` from H5AD `.obs`) + 21 peri-islet proportions (`peri_prop_*`) + 21 core+peri combined proportions (`combined_prop_*`, computed on the fly in `raw_df_base()`) alongside 3 hormone fractions. Grouped `selectInput` with up to 4 option groups (Hormone Positivity, Core, Peri-Islet, Core+Peri). The Combined group's choices are the intersection of `prop_*` and `peri_prop_*` stems — one entry per phenotype that has both region values.

### Age & Sex Demographic Filters
Age slider and sex checkboxes in Plot sidebar (H5AD only). `renderUI` returns `NULL` when absent.

### Segmentation + Drill-Down Pattern
Click a point -> `selected_islet()` updates -> embedded panel renders inline with Boundaries/Single Cells toggle. Close button sets `selected_islet(NULL)`.

### Metric Label Vocabulary (May 2026)
Sidebar radio labels and axis titles are unit-explicit — never use the bare word "density" without naming the denominator. Canonical phrasing:
- **Cell Populations Counts** → `"Cells per islet (core)"` / `"Cells per islet (peri)"` / `"Cells per islet (core+peri)"`.
- **Cell Populations Density** → selector `"Cells per mm² of core area"` / `"… peri area"` / `"… core+peri area"`; y-axis `"Cells per mm² (core area)"` / `"(peri area)"` / `"(core+peri area)"`. Denominator is the QuPath-measured annotation area (`core_region_um2`, `peri_region_um2`, or their sum), joined into `comp` by `prep_data()` from `targets_all` rows. Peri-islet area comes from the `_exp20um` annotation in `annotations.tsv`.
- **Microenvironment** → metric radio hidden (only one option). Y-axis names both the class and the region: `"Nerve area (% of core area)"` / `"Capillary area (% of peri-islet area)"` / `"Lymphatic area (% of core+peri area)"`. **Value computed from `100 * area_um2 / region_um2`** (% of the selected region's QuPath annotation area occupied by the structure). Counts are not exposed — these targets (Nerve / Capillary / Lymphatic) are area-valued (a single nerve may segment into several pieces), so object counts depend on segmentation morphology, not biology.
- **Core+Peri combined proportions** → computed on the fly: `value = (prop_X · total_cells_core + peri_prop_X · total_cells_peri)` divided by either `(total_cells_core + total_cells_peri)` (Percentage) or `(core_region_um2 + peri_region_um2)` (Density). NA totals/areas are treated as 0 so islets with only one zone still contribute.
- **Spatial Card C distance KDE** is **donor-aware** — it averages per-donor signed-distance densities with equal donor weight (a `shape`/`abundance` radio toggle), reading `data/signed_distance_kde.csv`. Y-axis is `"Per-donor mean density (area = 1)"` (shape, default) or `"Median cells per µm per donor"` (abundance), never a per-area density. This replaced the earlier pooled `density × N`, which pseudoreplicated large donors (ND 6356 = 71% of ND's APC cells). The `shape`/`abundance` split matches the `enrich_summary()` donor-aggregation philosophy (`mod_spatial_server.R`). Card C-right also has a `dist_view` toggle exposing a second, Senior-ported **"Per-cell distance"** view over the tissue-wide `data/cell_distance.csv` (all ~2.65M cells): a **symlog**-axis density (`cd_symlog` in `spatial_helpers.R`; cells per real µm; pooled per cell; by phenotype or all) or a **composition** (stacked bars faceted by status, zero-anchored at the islet border). Two shared toggles: **`dist_mode`** (Phenotype vs **Marker+** — per-cell `{marker}_pos` = scimap-rescaled X>0.5, CODEX's equivalent of Senior's RESTORE `_pos`; ~16 lineage markers via `cell_distance_markers()` + `CELL_DIST_MARKER_COLORS`) and **`dist_norm`** (**Proportion of class** default — normalizes by each disease group's total cell count so groups compare fairly despite different totals — vs Raw count). Marker mode overlaps (a cell may be + for several markers). Binning helpers `query_cell_distance_hist` / `query_cell_distance_composition` (`spatial_helpers.R`), pickers repopulated via `observeEvent(input$dist_mode)`. CSV + the scimap gate, no DuckDB (unlike Islet-Explorer-Senior, ported from).

A `helpText` block under the Plot metric selector restates the active mode's denominator. Keep this phrasing if you add new metrics. The y-axis class name is parsed from `input$which` (encoded as `className|region`).

## Inflammation Tab (Jul 2026)

Sixth tab, `mod_inflammation_ui.R` / `mod_inflammation_server.R`, wired in `app.R` after Spatial.
Carries two cohort-level analyses adapted from Steenbuck et al. 2025 (Bodenmiller lab, bioRxiv
2025.03.05.641526) that do not fit the Statistics tab's one-feature-at-a-time flow — both test
**every cell type at once**.

- **Card A — Donor cohort.** `data/donor_metadata.csv` per donor: age, sex, age group, T1D
  duration, BMI, HbA1c, C-peptide, AAb profile. States the cohort scope (all T1D donors
  recent-onset) and surfaces the donor-112 age discrepancy rather than hiding it.
- **Card B — edgeR differential abundance.** Per-donor cell-type counts → `filterByExpr` →
  `calcNormFactors` → `estimateDisp` → `glmQLFit` → `glmQLFTest`; three contrasts, three cell
  regions, volcano + table + CSV. **Counts, not proportions:** proportions are compositionally
  constrained, so the β-cell collapse in T1D inflates every other type's proportion; the
  library-size offset separates that artefact from a real abundance change. T1D vs ND gives
  Beta cell logFC −5.55 (FDR 0.0011), CD8a Tcell +2.71 (0.0013), Delta −2.88, Immune +3.23.
- **Card C — Infiltration score.** `1` inside the islet, `1/d` outside (d = µm from the edge,
  floored at 1 µm). See "Whole-slide false negatives" below — the denominator, the FDR family and
  the pairwise columns all had to be reworked for whole-slide data.

The tab is deliberately **absent** from the sidebar-visibility condition in `app.R`'s layout JS,
so the Plot sidebar stays hidden — these analyses read `data/cell_distance.csv` through the cached
`load_cell_distance()` in `spatial_helpers.R`, not the Plot module's islet frame, so sidebar
filters genuinely do not apply.

### Whole-slide false negatives in the infiltration score (Jul 2026)

Steenbuck et al. image **75 islet-centred ROIs per donor**; we image **whole sections**. Three
consequences, all of which suppress real signal, all measured:

**1. The denominator was the dominant defect.** Their per-image denominator is islet-local by
construction. Ours was every cell on the section — and the ratio of all cells to islet-proximal
cells ranges **4.1× to 23.2×** across our 15 donors. That variation is technical and is *not*
confounded with disease stage (Kruskal-Wallis p = 0.53 on the dilution factor), so it does not bias
direction — it inflates variance and costs power. Measured effect of switching the denominator from
whole-section to islet+peri:

```
              whole-section   islet+peri
T cell           0.1845         0.0081
Immune           0.1119         0.0207
B cell           0.0652         0.0155
CD8a Tcell       0.0271         0.0108
n surviving FDR      0             4
```

Four of seven family phenotypes cross from null to significant **on the denominator alone**. The
default is now `isletperi` (the ROI analogue); `near100` and `all` remain selectable, the last
labelled as a comparison rather than a baseline.

The **numerator needed no change**: 87.8% of the CD8a score mass already comes from within 20 µm of
an islet although only 16% of CD8a cells are there, so the `1/d` weight works as intended. Capping
the numerator at 100 µm changes p by nothing.

**2. There is a hard power ceiling, computed and displayed.** At 5 donors per stage the best
attainable Kruskal-Wallis p is **0.00193** and the best attainable pairwise Wilcoxon p is
**0.00794**, both requiring perfect separation of every donor. BH multiplies that floor by family
size:

```
KW family, m=21 phenotypes  -> best possible q = 0.041   (only a perfect result clears 0.05)
KW family, m=7  immune      -> best possible q = 0.014   (real headroom)
pairwise, m=21 contrasts    -> best possible q = 0.167   NOTHING can ever clear 0.05
```

`power_ceiling()` computes this from the actual group sizes and the card prints it beside the table,
so a null reads as *underpowered* rather than *no effect*.

**3. Therefore: pre-specified family, uncorrected pairwise.** `INFLAM_FDR_FAMILY` fixes the seven
immune phenotypes **in source** so the family cannot be resized after seeing results; phenotypes
outside it display uncorrected and marked `‡` with `q = NA`. The pairwise columns are reported
**raw** — a corrected pairwise column could never contain a significant value at this cohort size,
so printing one would be a guaranteed false negative dressed as a result.

**4. Numerator/denominator population mismatch — exposed, not hidden.** Steenbuck et al. sum the
weights over the *same* cells they divide by, so their score is a per-cell mean over one population.
Changing only the denominator makes ours a ratio of two populations: **29% of the numerator mass
lies outside the islet+peri denominator**. The unmatched form is the default *deliberately* —
matching to islet+peri confines the `1/d` gradient to the 20 µm peri band and guts the weighting the
score exists for (a T cell 25 µm out would contribute nothing). The *Match numerator to denominator*
checkbox restores the paper's exact statistic:

```
denominator   matched   CD8a     T cell   Immune   n surviving FDR
isletperi     FALSE     0.0108   0.0081   0.0207        4   <- default
isletperi     TRUE      0.0092   0.0140   0.0750        3
near100       FALSE     0.0122   0.0209   0.0358        3
near100       TRUE      0.0148   0.0209   0.0493        0
all           either    0.0271   0.1845   0.1119        0
```

Neither form is wrong; they answer slightly different questions, so the choice is a control rather
than a silent assumption.

**5. Nearest-islet attribution — checked, not a confounder.** On a whole slide `signed_dist` is to
the *nearest* islet, so in principle a peri-ductal or perivascular aggregate could be credited as
islet infiltration in a way a single-islet ROI makes impossible. Tested: islet crowding (mean
centroid nearest-neighbour distance) and islet count do not differ by stage (Kruskal-Wallis p = 0.83
and 0.53), and the score correlates weakly with islet geometry (|rho| ≤ 0.35 vs islet count, ≤ 0.26
vs crowding). What correlation exists is *negative* — more islets means a bigger denominator, not
accidental credit. Caveat: `islet_spatial_lookup.csv` holds 10,736 records for the 5,214 analysed
islets (roughly two entries per islet), so the crowding figure is valid for between-donor comparison
but not as an absolute distance.

Still unaddressed: 16–20% of our cells sit >200 µm from any islet (44% for donor 6548). Under the
default denominator they are excluded from it but still contribute to the numerator — that is the
mismatch quantified in point 4.

### Conventions this tab established (worth copying)

- **`shiny::validate(shiny::need(...))`, never the bare form.** `jsonlite` masks `validate` on this
  app's search path; the bare call resolves to the wrong function and dies with
  `is.character(txt) is not TRUE`. `app.R:420` already used the namespaced form for this reason.
- **FDR runs across everything present, not the display selection.** The infiltration card's
  phenotype picker is a display filter only; correcting across the ticked subset would let
  unticking phenotypes shrink the correction and manufacture significance. (Same lesson as the
  Spatial Card B filter: filter the summary, never the test.)
- **`renderTable(align = )` needs exactly one character per column and fails at RENDER time, not
  parse time.** Bitten twice in one session; recount when columns change.
- **Plotly legend `y` is a fraction of the PLOT AREA**, which shrinks as the bottom margin grows —
  so enlarging `margin$b` alone never fixes a legend that hangs off the canvas. Move the legend
  (reduce `|y|`) as well. Measured, not guessed.

## AI Assistant (Phase 10, Feb 2026)

### Architecture

The AI chat panel (`mod_ai_assistant_ui/server.R`) uses the University of Florida Navigator AI Toolkit via an OpenAI-compatible API. Wired as `ai_assistant_server("ai_chat")` in `app.R`.

**Key files:**
- `ai_helpers.R` -- Credential loading, API calls (streaming + non-streaming), model selection, error handling
- `mod_ai_assistant_ui.R` -- Chat panel UI with model picker, textarea, send/reset buttons
- `mod_ai_assistant_server.R` -- Chat history management, streaming callback, token budget tracking, `strip_markdown()` for plain-text display

### UF Navigator API

- **Endpoint**: `https://api.ai.it.ufl.edu/v1/chat/completions`
- **Available models** (as of Feb 2026):
  - `gpt-oss-20b` -- Fast model (default)
  - `gpt-oss-120b` -- Large reasoning model (chain-of-thought in `reasoning_content`, answer in `content`)
- **Authentication**: Bearer token from `KEY` env var (loaded from `.Renviron`)
- **Base URL**: `BASE` env var in `.Renviron` (defaults to `https://api.openai.com/v1` if unset)

### Reasoning Model Behavior (gpt-oss-120b)

The 120b model is a reasoning model (similar to DeepSeek R1):

- **Non-streaming**: Response has `message.content` (final answer) + `message.reasoning_content` (chain-of-thought). Content may be `null` if `max_tokens` is too low for reasoning + answer.
- **Streaming**: Chunks first arrive with `delta.reasoning_content` only (no `delta.content`). After reasoning completes, `delta.content` chunks appear with the actual answer. The streaming parser shows "Thinking..." during the reasoning phase.

### Credential Loading

`ai_helpers.R` searches for `.Renviron` in priority order:
1. `app/shiny_app/.Renviron` (deployment)
2. Script directory `.Renviron`
3. `R_ENVIRON_USER` env var
4. `~/.Renviron`
5. `~/.Renviron.local`
6. `$HOME/.Renviron`

Required env vars: `KEY` (API key), `BASE` (API base URL, optional).

### Model Fallback Logic

`select_openai_models()` returns a candidate list. For each candidate, `call_openai_chat()` tries streaming first, then non-streaming. If a model returns HTTP 400/401/404 with "model" in the error, it falls back to the next candidate. The fallback model defaults to `gpt-oss-120b` (env var `OPENAI_DEFAULT_MODEL`).

**Critical**: UF Navigator returns HTTP **401** (not 404) for model access denied. The fallback logic must include 401 in its status code checks.

### Important AI Conventions

- **Model names must match API**: Query `GET /v1/models` to verify available models.
- **Reasoning models need sufficient tokens**: `max_output_tokens` must be large enough for reasoning + answer.
- **UF Navigator uses `/chat/completions` only**: The code detects `api.ai.it.ufl.edu` in the base URL and skips the `/responses` endpoint.
- **Streaming parser must handle `reasoning_content`**: For reasoning models, `delta.reasoning_content` chunks arrive first — use a "Thinking..." indicator.
- **Token budget**: Controlled by `OPENAI_TOKEN_BUDGET` env var. Chat panel shows cumulative usage and blocks requests when budget exhausted.
- **Debug mode**: Set `DEBUG_CREDENTIALS=1` env var for verbose credential loading logs.
- **Plain-text display**: `strip_markdown()` in `mod_ai_assistant_server.R` removes markdown formatting from AI responses before display.

## Viewer Tab (Avivator) — iframe Persistence (Jun 2026)

The Viewer embeds a local Avivator build (`www/avivator`, Vite bundle) in an `<iframe>` to render OME-TIFFs. Files: `mod_viewer_server.R`, `mod_viewer_ui.R`, `viewer_helpers.R`.

**Problem (macOS-only flashing)**: a single `renderUI` emitted the dropdown AND the iframe together, so any reactive invalidation destroyed + recreated the `<iframe>` (a WebGL canvas). Repeated canvas add/remove flickers heavily under Firefox's WebRender compositor on macOS (Mozilla bug #1555544); Retina also doubles deck.gl's drawing buffer (`useDevicePixels`).

**Architecture (the fix)**:
- **Split outputs**: `output$image_selector` (dropdown; depends on the stable available-images list) and `output$viewer_frame` (iframe; mounts once).
- **Mount-once iframe**: `viewer_frame` depends on `has_image_stable()` (flip-only boolean, below) + `current_tab()`, and reads the URL via `isolate(viewer_info()$iframe_src)` so it does NOT re-render when the image URL changes.
- **In-place src update**: `observeEvent(viewer_info()$iframe_src, ...)` → `session$sendCustomMessage("islet_update_viewer_src", list(id, src))`. A JS handler registered once in `mod_viewer_ui.R` sets `iframe.src` only when it differs from a `data-relsrc` attribute (dedup). Switching images updates the existing node — never recreates it (verified via chromote: same DOM node, src updated).
- **`probe_viewer_asset` gated** behind `VIEWER_DEBUG_ENABLED` (it's a blocking `curl` HEAD); `session$clientData` URL reads `isolate()`d (origin is fixed per session).

**CRITICAL gotcha — flip-only boolean**: a plain `reactive()` boolean invalidates its dependents whenever its INPUTS change (Shiny propagates on invalidation, NOT value-equality). So a `has_image` reactive reading `input$selected_image` made `viewer_frame` re-render (recreate the iframe) on every image switch. Fix: mirror it into a `reactiveVal` via `observe({ has_image_stable(isTRUE(has_image_now())) })` — `reactiveVal` dedups, so it only propagates when the boolean actually FLIPS (no-image → image), keeping the iframe alive across image switches.

**Secondary mitigation** (if macOS flicker persists after the above): rebuild Avivator with `useDevicePixels: false`. The bundle only reads `image_url` + `channel_config` query params, so this can't be toggled at runtime.

## App-Specific Conventions

- **CSS overflow for cards with dropdowns**: `selectInput` menus extend below their container. Add `overflow: visible;` to card styles or dropdowns get clipped behind cards.
- **Font size minimums**: Use `h5` (not `h6`) with explicit `font-size: 15px` for panel headings. Legend items minimum 14-15px.
- **Large scatter plots (>50K points)**: Use `ggplot2::renderPlot()` NOT `plotlyOutput()` -- plotly cannot handle 100K+ points interactively. Use small `size` (0.3-0.5) and low `alpha` (0.3-0.6).
- **Tissue scatter coordinate convention**: `coord_fixed() + scale_y_reverse()` matches microscopy convention (y increases downward, spatial proportions preserved)
- **Leiden cluster mapping for cells**: Islet-level Leiden assignments map to single cells via `islet_name` column lookup. Tissue background cells without an islet get `cluster = "tissue"` (grey).
- **Per-donor tissue CSVs**: `data/donors/{imageid}.csv` with 5 columns (X_centroid, Y_centroid, phenotype, cell_region, islet_name). `cell_region` = core/peri/tissue.
- **Leiden in H5AD**: 4 resolution columns (`leiden_0.3/0.5/0.8/1.0`) + 2 UMAP coords (`leiden_umap_1/2`). Extracted via `^leiden_` regex in `data_loading.R`.
- **Leiden UMAP uses visualization coords**: `islets_core_clustered.h5ad` X_umap is copied from trajectory's raw marker PCA UMAP, NOT computed separately from scVI (which produces a blob). Leiden clustering itself still uses scVI neighbors.
- **Donor palette syncing**: 3 non-namespaced `selectInput`s (`sidebar_donor_palette`, `traj_donor_palette`, `spatial_donor_palette`) synced bidirectionally via `observeEvent` + `updateSelectInput` in `app.R`. All feed `donor_palette_name()` -> `active_donor_colors()`. Default selection: **Bright** (`DONOR_COLORS_BRIGHT` = `#3366CC` / `#FF6600` / `#109618`); the legacy "Paul Tol" palette remains available as a choice. The fallback `get_donor_color_palette()` default arg in `00_globals.R` is also `"Bright"` so headless/test contexts match.
- **Trajectory legend uses reactive colors**: `donor_colors_reactive()` NOT hardcoded hex. Changes when palette selector updates.
- **Spatial background cell coloring**: `color_background` checkbox. When checked + phenotype mode, background cells get phenotype colors at alpha=0.25.
- **Reticulate bridge crossing**: `ad$uns[[key]]` costs ~221ms per call. ALWAYS cache `uns <- ad$uns` as an R list first, then index. `reconstruct_groovy_df_from_list()` uses cached list; legacy `reconstruct_groovy_df()` is a thin wrapper.
- **Deferred heavy file loading**: `annotations.tsv` (72 MB) lazy-loaded via `.seg_lazy` environment in `segmentation_helpers.R`.
- **Synthetic union rows need density**: `prep_data()` synthesizes `islet_union` rows by summing core+band. Must compute `area_um2`, `region_um2`, and `area_density = area_um2/region_um2` — setting these to `NA` causes Microenvironment Core+Peri to drop 12/15 donors (default metric is Density, NA filtered by `!is.na(value)`).
- **`updateSelectInput` after `renderUI`**: Server-side updates targeting an input that lives inside a `renderUI`-mounted block must share a reactive dependency with the renderUI, or fire only after the UI has mounted on the client. Putting the update inside its own naked `observe({})` can race the UI mount and silently leave the dropdown empty. Pattern: nest the `updateSelectInput` inside the same observer that builds the choices reactive.
- **`reactiveVal` to stabilise a `renderUI` boolean**: a `reactive()` boolean re-fires its dependents whenever its INPUTS change, even if the value is unchanged — so a `renderUI` keyed on it rebuilds (destroying child DOM, including iframes/WebGL canvases). Mirror the boolean into a `reactiveVal` via an `observe`; `reactiveVal` dedups, so the `renderUI` only re-renders when the boolean truly flips. See § "Viewer Tab (Avivator) — iframe Persistence".
- **Robust plotly clicks**: `plotly_click` can be suppressed on macOS (clicking does nothing, no notification). Use `islet_click_bridge()` (native mousedown/mouseup → `Shiny.setInputValue`) for click-to-select, and guard `onRender` with an attach-once flag so listeners don't stack. See § "macOS Native-Click Bridge".

## Phenotype Rules — Marker Definitions in the UI (May 2026)

Source of truth: `data/phenotype_rules.csv`. Parsed once at startup in
`R/00_globals.R` into `PHENOTYPE_RULES_CSV` (19 phenotypes, marker columns
`INS, GCG, SST, SMA, CD20, CD3e, CD68, CD163, CD8a, CD4, HLADR, PDPN,
CD34, CD31, B3TUBB`). The composer walks each phenotype's parent chain
(e.g. `CD8 T cell` → `T cell` → `Immune`) and unions positive/negative
constraints, then strips redundant `anypos` groups that are supersets of a
stricter constraint downstream. Four residual phenotypes (`ECAD+`, `APCs`,
`Structural`, `Unknown`) have no CSV row — they live in `PHENOTYPE_RESIDUALS`
with hand-written descriptions.

Naming aliases: PHENOTYPE_COLORS (drilldown_helpers.R) uses
`CD8a Tcell` / `CD4 Tcell` / `Blood Vessel` / `SMA+`; the CSV uses
`CD8 T cell` / `CD4 T cell` / `Blood vessel` / `Smooth Muscle`. The
`PHENOTYPE_ALIASES` map normalises PHENOTYPE_COLORS → CSV before lookup.

**Canonical display label (Jul 2026):** the PHENOTYPE_COLORS / data
`phenotype`-column form is the single canonical label shown to users
everywhere — drill-down legend, composition chart, Plot and Spatial selectors,
and (as of Jul 2026) the Phenotype-Rules reference modal. The modal reads the
CSV, so it previously displayed `Smooth Muscle` / `Blood vessel` / `CD8 T cell`
for phenotypes the rest of the app calls `SMA+` / `Blood Vessel` / `CD8a Tcell`.
`phenotype_display_label()` (backed by `PHENOTYPE_ALIASES_REV`, the inverse of
`PHENOTYPE_ALIASES`) maps any CSV/rule-form name back to the canonical label;
`phenotype_rules_table_df()` applies it to the modal's Phenotype column so all
views agree. Note: marker `SMA` (expression / gating feature) and phenotype
`SMA+` (smooth-muscle cell type) are intentionally distinct — the `+` marks the
phenotype, not a spelling variant of the marker.

Helpers:

- `get_phenotype_rule(name)` — returns `list(parent, pos, neg, anypos_groups,
  description)` for either CSV form or PHENOTYPE_COLORS form. Handles aliases
  and residuals.
- `format_phenotype_rule(name, style = "inline" | "html" | "plain")` — formats
  the rule as `"CD3e+, CD8a+  —  CD4−, CD20−"` (inline), styled HTML
  (`<strong>` positives, muted negatives), or plain text for plot subtitles.
- `extract_phenotype_from_column(col)` — resolves a sanitised column name
  (e.g. `prop_CD8a_Tcell`, `peri_prop_Beta_cell`, `min_dist_Macrophage`,
  `enrich_z_B_cell`) back to a phenotype label by stripping known prefixes
  and reversing `_`→space / `plus`→`+`.
- `phenotype_hint_ui(name, show_modal_link = "<id>")` — small blue inline
  hint block embedded under each phenotype selector. The trailing
  "see all rules" actionLink triggers the same modal as the tab-top button.
- `phenotype_rules_button(input_id)` — pill-styled actionLink placed at the
  top of every tab. Each tab uses a unique non-namespaced id
  (`show_phenotype_rules_<tab>`); a single shared loop in `app.R`
  registers an `observeEvent` per id that fires `showModal(phenotype_rules_modal_ui())`.
- `phenotype_rules_modal_ui()` — returns the modal: a sortable table of all
  23 phenotype names (19 CSV + 4 residuals) with category, positive markers,
  negative markers, and notes.

Per-tab integration:

- **Plot tab** (`mod_plot_server.R::output$which_ui`) — when in `"Cell
  Populations"` mode, parses `input$which` via `extract_phenotype_from_column()`
  and emits `phenotype_hint_ui()` directly under the metric `selectInput`.
- **Spatial Card A** — `output$infiltration_pheno_hint` and
  `output$scatter_pheno_hint` look up the selected `safe` name via
  `phenotype_options()` and emit the hint UI.
- **Spatial Card C** — `output$distance_pheno_hint` and `output$kde_pheno_hint`
  follow the same pattern. The aggregate `"Immune (all)"` selection emits a
  custom description rather than a per-phenotype rule.
- **Drill-down composition panel** — root-level `output$islet_drilldown_rules`
  in `app.R` renders an expandable `<details>` element listing colour-coded
  swatches + rules for every phenotype present in the clicked islet. Embedded
  via `uiOutput("islet_drilldown_rules")` inside both the Plot and Trajectory
  `segmentation_viewer_panel` renderUIs.

When adding a new phenotype-selecting input anywhere in the app, follow the
existing pattern: resolve the selected value to a label via
`extract_phenotype_from_column()` (or your existing `phenotype_options()`
record), then emit `phenotype_hint_ui(label, show_modal_link = ns("…"))` and
register an `observeEvent(input$…, showModal(phenotype_rules_modal_ui()))` to
wire the inline link.

## Plot Fonts & Publication-Quality Downloads (May 2026)

Single source of truth in `R/00_globals.R`:

- `ISLET_FONT_SIZES` — named list with `base=18`, `title=18`, `subtitle=15`, `axis_title=16`, `axis_text=16`, `legend_title=16`, `legend_text=16`, `strip=16`, `dense_min=14` (heatmap rows / in-cell values), `hover=15`, `annotation=15`. (Raised from the May 2026 16/14/12 scale in Jul 2026.)
- `theme_islet(base_size = 18)` — ggplot theme extending `theme_minimal()`. Use on every `renderPlot` / inline `ggplot()` that is also handed to `ggplotly()`.
- `plotly_islet_fonts()` — returns a named list of `list(size = N)` font specs (keys: `global`, `title`, `axis_title`, `axis_tick`, `legend_title`, `legend_text`, `colorbar_title`, `hover`, `annotation`). Splice into `layout(font=..., title=list(font=...), xaxis=list(title=list(font=...), tickfont=...), yaxis=..., legend=list(font=..., title=list(font=...)))`. `ggplotly()` does **not** fully respect ggplot themes — the explicit `layout(...)` overlay is required after every `ggplotly()` call.
- `plotly_islet_config(p, filename_root, width=1400, height=900, scale=2, uniform_text_mode="hide")` — pipe with `%>%` after `layout()`. Sets `toImageButtonOptions` (modebar camera exports 2800 × 1800 px PNG), disables the plotly logo, **and applies the two in-plot text settings `layout(font=)` does not cascade to** (see below). Apply to **every** plotly output except deliberate UI strips like the trajectory pseudotime heatmap (`mod_trajectory_server.R` uses `config(displayModeBar = FALSE)` and is left as-is).

**In-plot plotly text — the settings `layout(font=)` does NOT reach (Jul 2026).** Three separate mechanisms, all previously unset on every plot in the app, and collectively the reason text still read small after the May 2026 pass:

- **`hoverlabel.font`** — tooltips render at plotly's ~13px default no matter what `layout(font=)` says. Several plots (Leiden UMAP, tissue scatter, Card A bars) are read *primarily* via hover. Now applied centrally by `plotly_islet_config()`.
- **`uniformtext`** — for `bar` and `heatmap` traces plotly silently **shrinks** in-cell text to fit the box, so a `textfont = list(size = 14)` can render at 8px. `minsize` sets the floor; `mode` decides the overflow behaviour: `"hide"` (default) drops the label — right for numeric heatmap values, which overlap illegibly — and `"show"` forces it, which is required for the Statistics `bin_heatmap` significance stars (a hidden star reads as "not significant"). Pass `uniform_text_mode = "show"` there.
- **`automargin = TRUE`** — required on any axis with `tickangle` or long category names, else plotly clips the labels against a fixed margin. Clipped labels read as "the font is too small" when the label is really just cut off. Set on the enrichment bars/heatmap, the bin heatmap, and the ggplotly-converted Statistics plots.

**Two ggplotly conversion losses** (both silent, both fixed centrally):

- **`labs(subtitle = ...)` is DROPPED.** Three Statistics plots put their correlation / BH-correction / Kruskal-Wallis p-values there, so those lines never reached the page at all. Use `plotly_title_with_sub(title, subtitle)` to rebuild the title as `"<title><br><sup><subtitle></sup>"` and pass it to `layout(title = list(text = ...))`. Add `margin = list(t = 80)` so the two-line title has room.
- **Faceting relocates axis titles into `layout.annotations`.** For an unfaceted plot ggplotly keeps them in `layout.xaxis.title`, so the usual overlay works. Add a facet and *both* the strip labels and the axis titles move to `annotations`, where the overlay is ignored and they render at ~14.6px. `plotly_islet_strip_fonts(gg)` restyles both — strips get `strip`, `annotationType == "axis"` entries get `axis_title`. Pipe it directly after `ggplotly()`, before `layout()`.

**`annotate()` / `geom_text()` `size` is MILLIMETRES, not points.** `size = 4` is 11.4pt and `size = 5` is 14.2pt — both *below* the 18pt base, which is what the five placeholder annotations in this app were using. Use `ISLET_ANNOTATE_SIZE` (= `ISLET_FONT_SIZES$base / ggplot2::.pt`) instead of a literal.

**Empty states**: use `plotly_islet_empty(msg)` rather than `plotly_empty() %>% layout(title = msg)` — the latter renders its title at plotly's default size.
- `ISLET_PNG_DIMS` (10 × 7 in @ 300 dpi → 3000 × 2100 px) and `ISLET_HEATMAP_PNG_DIMS` (12 × 8 in) for `ggsave()` in `downloadHandler()` `content` functions.
- `islet_png_filename(stem)` — `<stem>_YYYYMMDD_HHMMSS.png`.

**Pattern for ggplot downloads** — when a `renderPlot` output should be downloadable, refactor it into a `reactive()` and a thin `renderPlot({ plot_reactive() })` so both the renderer and the `downloadHandler` consume the same ggplot object. Then `ggsave(file, plot = plot_reactive(), width = ISLET_PNG_DIMS$width, height = ISLET_PNG_DIMS$height, dpi = 300, units = "in")`. Examples: `mod_trajectory_server.R` UMAP donor/feature + multi-feature heatmap (`traj_umap_donor_plot()`, `traj_umap_feature_plot()`, `traj_multi_heatmap_plot()`); `mod_spatial_server.R` tissue scatter and donor UMAP (`spatial_tissue_scatter_plot()`, `spatial_umap_donor_plot()`); `app.R` segmentation + drilldown (`islet_segmentation_view_plot()`, `islet_drilldown_summary_plot()`).

**Download button placement** — non-namespaced `downloadButton("dl_islet_segmentation_view", ...)` for root-level outputs (segmentation, drilldown) lives inside each module's `segmentation_viewer_panel` renderUI alongside the plotOutput; namespaced `downloadButton(ns("dl_..."), ...)` for module-local plots lives in the module UI immediately under the plotOutput. Render handlers for the non-namespaced download buttons sit at root in `app.R`; namespaced ones live in the module server.

## CSV Downloads (Jul 2026)

Nine CSV exports: Plot `dl_summary` + `dl_distribution_csv`; Spatial `download_spatial` + `dl_infiltration_csv` + `dl_enrichment_csv` + `dl_distance_csv`; Trajectory `dl_multi_heatmap_csv` + `dl_umap_csv`; Statistics `download_stats`. Two rules, both learned from bugs:

**1. A CSV must export the plot's own reactive, never a re-derivation.** `download_spatial` used to build its own subset from `prepared()$comp` instead of the shared `nbr_comp()` every card consumes — so it ignored `nbr_min_cells` and `pheno_filter`, applied a `total_cells_peri > 0` filter the cards do *not*, and carried no outlier flags. The result corresponded to **no plot on the page**. When a card's data is computed inline in `renderPlotly`, extract it into a reactive first and have the renderer consume it (`infiltration_data()`, `distance_data()`, `dist_data()`, `traj_multi_heatmap_matrix()` were all extracted for exactly this). Then the export cannot drift.

Watch for panels with their *own* filter state: the Plot Distribution panel reads `raw_df_base()` plus its own `exclude_zero_dist`, independent of the main scatter's `exclude_zero_top` — exporting `raw_df()` there would be wrong.

**2. Every CSV carries a `#` provenance header.** Use `write_islet_csv(file, df, header)` and `islet_csv_header(title, ...)` from `00_globals.R`. Writing to an open connection also avoids R's "appending column names to file" warning that `writeLines()` + `write.table(append = TRUE)` produces. The header must name **everything that changes which rows are in the file** — donor groups, min-cells, phenotype filter, diameter/age/sex filters, and the outlier + zero-exclusion state. `get_selection_description()` (Plot module) omitted min-cells, age, sex, outlier and zero state, so two files with different row sets looked identical in provenance.

**Column naming**: user-facing headers only. The `gender` data column is renamed to `sex` at export time (in-memory it stays `gender` so `lm()` / `aggregate_to_donor()` are untouched). Per-islet exports keep `is_outlier` / `outlier_score` and do **not** drop flagged rows, so the full distribution stays recoverable.

**Don't** mix the helper with manual `font = list(size = X)` per-key calls — the helper is the entire convention. When you need a smaller font for dense rows (e.g., heatmap text annotations) reach for `ISLET_FONT_SIZES$dense_min`, not a fresh literal, and not a `max(dense_min, other - 2)` formula (two of those existed and both evaluated to exactly `dense_min`).

**Container height is part of font sizing.** 18pt in a 260px box is the same complaint as an 11pt label. The Jul 2026 pass grew the cramped ones: Statistics AAb 260→380, age scatter / sex box 300→420, AUC 280→380, forest 350→430, bin heatmap 400→520; Spatial Leiden UMAP 300→400, cluster composition 350→450, enrichment bars 420→500, enrichment heatmap 420→520; Trajectory multi-heatmap 36→44 px per marker row. When adding a plot, size the container to the text, not the other way round.
