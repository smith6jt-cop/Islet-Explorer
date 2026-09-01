library(shiny)
library(shinyjs)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(plotly)
library(broom)
library(jsonlite)

# Pixel size constant for coordinate conversion (micrometers per pixel)
# GeoJSON polygons use pixel coordinates, while islet_spatial_lookup.csv uses micrometers
PIXEL_SIZE_UM <- 0.3774  # micrometers per pixel

# Donor status color palettes
DONOR_COLORS <- c("ND" = "#4477AA", "Aab+" = "#CC6633", "T1D" = "#228833")

DONOR_COLORS_BRIGHT <- c("ND" = "#3366CC", "Aab+" = "#FF6600", "T1D" = "#109618")

DONOR_COLORS_COLORBLIND <- c("ND" = "#0072B2", "Aab+" = "#D55E00", "T1D" = "#009E73")

DONOR_COLORS_CLASSIC <- c("ND" = "#2166AC", "Aab+" = "#B2182B", "T1D" = "#1B7837")

DONOR_PALETTES <- list(
  "Paul Tol (default)" = DONOR_COLORS,
  "Bright"             = DONOR_COLORS_BRIGHT,
  "Okabe-Ito"          = DONOR_COLORS_COLORBLIND,
  "Diverging"          = DONOR_COLORS_CLASSIC
)

get_donor_color_palette <- function(name = "Bright") {
  DONOR_PALETTES[[name]] %||% DONOR_COLORS_BRIGHT
}

# Load sf package for GeoJSON/spatial operations (islet segmentation viewer)
SF_AVAILABLE <- FALSE
suppressPackageStartupMessages({
  if (requireNamespace("sf", quietly = TRUE)) {
    library(sf)
    SF_AVAILABLE <- TRUE
    message("[SF] Package sf loaded for segmentation viewer")
  } else {
    message("[SF] Package 'sf' not found. Islet segmentation viewer will be disabled.")
    message("[SF] Install with: install.packages('sf')")
  }
})

# Check for httr2 package (required for AI assistant)
httr2_available <- requireNamespace("httr2", quietly = TRUE)
if (!httr2_available) {
  message("Package 'httr2' not found. The AI assistant will be disabled until you install it (install.packages('httr2')).")
}

# Load anndata package for trajectory analysis (install with BiocManager::install('anndata'))
suppressPackageStartupMessages({
  if (!requireNamespace("anndata", quietly = TRUE)) {
    message("Package 'anndata' not found. Install with: BiocManager::install('anndata')")
  }
})

# Load vitessceR for embedded image viewer (GitHub v0.1.0 required)
VITESSCE_AVAILABLE <- FALSE
suppressPackageStartupMessages({
  if (!requireNamespace("vitessceR", quietly = TRUE)) {
    message("Package 'vitessceR' not found. Install with: remotes::install_github('vitessce/vitessceR')")
    message("Version required: >= 0.1.0")
  } else {
    # Verify we have the GitHub version with required features
    tryCatch({
      pkg_version <- packageVersion("vitessceR")
      if (pkg_version < "0.1.0") {
        warning("vitessceR version ", pkg_version, " detected. Please upgrade to >= 0.1.0 from GitHub")
      }
      test_vc <- vitessceR::VitessceConfig$new(schema_version = "1.0.16", name = "Test")
      message("[VITESSCE] Package verified working (version ", pkg_version, ")")
      VITESSCE_AVAILABLE <<- TRUE
    }, error = function(e) {
      warning("[VITESSCE] Package installed but failed basic test: ", conditionMessage(e))
    })

    # CRITICAL FIX: Ensure all child routes from MultiImageWrapper are registered
    # This fixes segmentation layers not loading properly and appearing at wrong sizes
    if (VITESSCE_AVAILABLE) {
      tryCatch({
        vitessceR::MultiImageWrapper$set("public", "get_routes", function() {
          routes <- list()
          if (length(self$image_wrappers) == 0L) return(routes)
          for (w in self$image_wrappers) {
            wr <- w$get_routes()
            if (wr > 0L) routes <- c(routes, wr)
          }
          routes
        }, overwrite = TRUE)
        message("[VITESSCE] Applied MultiImageWrapper get_routes fix for segmentation layer routing")
      }, error = function(e) {
        warning("[VITESSCE] Failed to apply get_routes fix: ", conditionMessage(e))
      })
    }
  }
})

## for base64 encoding channel_config payload to Avivator
suppressPackageStartupMessages({
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop("Package 'base64enc' is required. Run scripts/install_shiny_deps.R to install dependencies.")
  }
})

# ---- Runtime diagnostics ----
APP_VERSION <- "modular-v1"
try({
  message("[APP LOAD] Version=", APP_VERSION, " time=", Sys.time())
  message("[APP LOAD] getwd=", tryCatch(getwd(), error = function(e) NA))
  message("[APP LOAD] list.files(.) contains app.R? ", any(grepl('^app.R$', list.files('.'))))
  # Publication containers keep the application and release mounts read-only.
  # Write a load flag only when an operator explicitly requests one, and only
  # below the process temporary directory (never below the source tree).
  flag_path <- Sys.getenv("ISLET_APP_LOAD_FLAG", unset = "")
  if (nzchar(flag_path)) {
    temp_root <- normalizePath(tempdir(), mustWork = TRUE)
    flag_path <- normalizePath(flag_path, mustWork = FALSE)
    is_below_temp <- identical(flag_path, temp_root) ||
      startsWith(flag_path, paste0(temp_root, .Platform$file.sep))
    if (!is_below_temp) {
      warning("ISLET_APP_LOAD_FLAG must resolve below tempdir(); flag not written.")
    } else {
      if (!dir.exists(dirname(flag_path))) {
        dir.create(dirname(flag_path), recursive = TRUE, showWarnings = FALSE)
      }
      writeLines(paste0('loaded ', Sys.time(), ' version=', APP_VERSION), flag_path)
    }
  }
}, silent = TRUE)

# Runtime paths. Production sets ISLET_DATA_ROOT to the mounted, versioned data
# release; a repository checkout keeps the historical ../../data default.
.repo_data_root <- file.path("..", "..", "data")
data_root <- Sys.getenv("ISLET_DATA_ROOT", unset = .repo_data_root)
data_root <- tryCatch(normalizePath(data_root, mustWork = FALSE),
                      error = function(e) data_root)
data_path <- function(...) file.path(data_root, ...)

master_path <- data_path("master_results.xlsx")
h5ad_path   <- data_path("islet_explorer.h5ad")
project_root <- tryCatch(normalizePath(file.path("..", ".."), mustWork = FALSE), error = function(e) NULL)

# Ensure reticulate/anndata can discover a Python binary when RETICULATE_PYTHON
# isn't set (shiny-server environments often lack PATH).
#
# The interpreter MUST have the python `anndata` module, otherwise read_h5ad()
# throws and load_master_auto() silently falls back to the Excel source — which
# lacks the prop_*/peri_prop_* columns, gutting the Plot tab's cell-type
# selector. We therefore can't just grab the first python on PATH: the deploying
# user's venv (which ships anndata) may live under a different home than the one
# running the app. shiny-server runs as `shiny` (home /home/shiny), so the old
# path.expand("~") lookup found nothing and fell through to a bare system
# python3 with no anndata. Search several candidate locations and pick the first
# interpreter that can actually `import anndata`.
if (!nzchar(Sys.getenv("RETICULATE_PYTHON", ""))) {
  .venv_rel <- file.path(".local", "share", "islet-explorer-py", "bin", "python")
  .py_candidates <- unique(c(
    file.path(path.expand("~"), .venv_rel),               # current user's venv
    Sys.glob(file.path("/home", "*", .venv_rel)),         # any user's venv (deployment)
    Sys.glob(file.path("/srv", "*", .venv_rel)),
    Sys.getenv("PYTHON", ""),
    Sys.which("python3"),
    Sys.which("python")
  ))
  .py_candidates <- .py_candidates[nzchar(.py_candidates) & file.exists(.py_candidates)]

  .py_has_anndata <- function(py) {
    status <- suppressWarnings(tryCatch(
      system2(py, c("-c", shQuote("import anndata")),
              stdout = FALSE, stderr = FALSE),
      error = function(e) 1L))
    identical(as.integer(status), 0L)
  }

  candidate_python <- ""
  for (.py in .py_candidates) {
    if (.py_has_anndata(.py)) { candidate_python <- .py; break }
  }
  # If none advertised anndata (e.g. the capability probe couldn't run), fall
  # back to the first existing candidate so behaviour is no worse than before.
  if (!nzchar(candidate_python) && length(.py_candidates)) {
    candidate_python <- .py_candidates[[1]]
    message("[PYTHON] WARNING: no candidate interpreter could import anndata; ",
            "H5AD loading may fail and fall back to Excel. Set RETICULATE_PYTHON ",
            "to a python that has the anndata module.")
  }

  if (nzchar(candidate_python)) {
    Sys.setenv(RETICULATE_PYTHON = candidate_python)
    message("[PYTHON] RETICULATE_PYTHON=", candidate_python)
  } else {
    message("[PYTHON] No python executable found. Set RETICULATE_PYTHON to a valid interpreter.")
  }
}

# Debug flags
DEBUG_CREDS <- isTRUE(as.logical(Sys.getenv("DEBUG_CREDENTIALS", "0")))
VIEWER_DEBUG_ENABLED <- identical(Sys.getenv("VIEWER_DEBUG", "0"), "1")

# Feature flags
#
# Inflammation tab: HIDDEN by default (Aug 2026, at the user's request). The
# module and its docs are untouched — this only controls whether the tab is
# mounted, so re-enabling is a one-line env change rather than a code revert.
# Both the tabPanel and the inflammation_server() call in app.R are gated on it;
# leaving the server running without its UI would be harmless (Shiny reactives
# are lazy) but would still be dead wiring, so both are gated together.
#   SHOW_INFLAMMATION_TAB=1 in .Renviron  ->  tab reappears
SHOW_INFLAMMATION_TAB <- identical(Sys.getenv("SHOW_INFLAMMATION_TAB", "0"), "1")

# Statistics tab: HIDDEN by default (Aug 2026, at the user's request). Same
# pattern as the Inflammation tab above — the module, its docs and its tests are
# untouched, so re-enabling is a one-line env change rather than a code revert.
# Both the tabPanel and the statistics_server() call in app.R are gated on it.
# The Plot sidebar's visibility condition still names 'Statistics' in both the
# layout JS and the conditionalPanel: with the tab absent `input.tabs` can never
# equal 'Statistics', so those branches are inert, and leaving them means
# flipping this flag back on needs no second edit.
#   SHOW_STATISTICS_TAB=1 in .Renviron  ->  tab reappears
SHOW_STATISTICS_TAB <- identical(Sys.getenv("SHOW_STATISTICS_TAB", "0"), "1")

# ---------------------------------------------------------------------------
# Trajectory: the single pseudotime definition, and the islet cell-count range.
# ---------------------------------------------------------------------------
#
# There is deliberately ONE pseudotime. The tab previously offered five
# selectable "modes" (core / core+peri scVI, two cell-count strata, and a
# composition one); every additional mode was another axis whose numbers were
# not comparable with the others', which is a garden-of-forking-paths surface
# rather than a feature. The one that survived is the islet-composition
# diffusion map -- it gave the strongest beta-cell correlation measured
# (rho -0.70 to -0.73 above 25 cells/islet, against -0.57 for the scVI variants)
# and the least compressed axis. See app/shiny_app/docs/trajectory.md.
#
# The pseudotime is computed on the CURRENT cell-count range, not on a fixed
# islet set: the range is part of the definition, so moving it refits the
# diffusion map rather than filtering a frozen axis.

# Largest islet in the cohort, used as the upper bound of the range slider.
# Read from the composition features so it tracks the data; the literal is a
# fallback for when that file is absent.
TRAJ_MAX_ISLET_CELLS_FALLBACK <- 1902L

traj_max_islet_cells <- function() {
  f <- data_path("barlow_islet_features.csv")
  if (!file.exists(f) || !requireNamespace("data.table", quietly = TRUE)) {
    return(TRAJ_MAX_ISLET_CELLS_FALLBACK)
  }
  v <- tryCatch({
    d <- data.table::fread(f, select = "total_cells_core")
    as.integer(max(suppressWarnings(as.numeric(d$total_cells_core)), na.rm = TRUE))
  }, error = function(e) NA_integer_)
  if (is.finite(v) && v > 1) v else TRAJ_MAX_ISLET_CELLS_FALLBACK
}

# Default lower bound. Barlow et al. exclude islets with fewer than ten cells
# before computing their embedding; below that an islet's composition vector is
# mostly sampling noise (60.8% of our islets are <= 10 cells).
TRAJ_DEFAULT_MIN_CELLS <- 10L

#' Normalise the cell-range slider to c(lo, hi). Guards the same NA hole the old
#' numericInput had: a cleared input sends NA, and `if (NA > 1)` aborts the whole
#' reactive with "missing value where TRUE/FALSE needed".
traj_cell_range <- function(x, max_default = NULL) {
  hi_default <- suppressWarnings(as.numeric(max_default %||% TRAJ_MAX_ISLET_CELLS_FALLBACK))[1]
  if (!is.finite(hi_default) || hi_default < 1) {
    hi_default <- TRAJ_MAX_ISLET_CELLS_FALLBACK
  }
  hi_default <- floor(hi_default)
  default_lo <- min(TRAJ_DEFAULT_MIN_CELLS, hi_default)
  v <- suppressWarnings(as.numeric(x))
  v <- v[is.finite(v)]
  if (length(v) < 2) return(c(default_lo, hi_default))
  lo <- min(max(1, floor(min(v[1:2]))), hi_default)
  hi <- min(max(lo, ceiling(max(v[1:2]))), hi_default)
  c(lo, hi)
}


# AI chat card: OFF, so it is both disabled and hidden. Gating rather than
# deleting keeps the module files intact and makes this a one-line env change
# rather than a code revert, exactly as with the Inflammation tab above. With it
# off nothing reads the UF Navigator credentials (KEY/BASE), so no key is needed
# to run the app, and the Plot tab's Distribution card widens 4 -> 6 to fill the
# row the chat card used to occupy.
#   ISLET_ENABLE_AI=TRUE in .Renviron  ->  card reappears
AI_ENABLED <- isTRUE(as.logical(Sys.getenv("ISLET_ENABLE_AI", "FALSE")))

# ============================================================================
# Plot font + theme helpers (May 2026)
# ============================================================================
# Centralised font sizing so every plot in the app reads clearly on screen and
# downloads as a presentation/publication-ready PNG. "Presentation-first"
# baseline: base 16pt, 18pt titles, 14pt legends, 12pt minimum for dense rows.

ISLET_FONT_SIZES <- list(
  base         = 18,
  title        = 18,
  subtitle     = 15,
  axis_title   = 16,
  axis_text    = 16,
  legend_title = 16,
  legend_text  = 16,
  strip        = 16,
  dense_min    = 14,  # heatmap row labels, in-cell values, dense tick labels
  hover        = 15,  # plotly tooltips — NOT covered by layout(font=)
  annotation   = 15   # plotly annotations / ggplotly facet strips
)

# ggplot's `annotate("text", size = N)` / `geom_text(size = N)` take
# MILLIMETRES, not points. Passing a point-like number (14, 16) silently
# produces enormous text; passing 4 or 5 — as this app did — produces 11.4pt /
# 14.2pt, i.e. BELOW the 18pt base. Convert explicitly.
ISLET_ANNOTATE_SIZE <- ISLET_FONT_SIZES$base / ggplot2::.pt

# Standard ggplot theme. Use directly: `theme_islet()` or extend with
# additional `theme(...)` calls for plot-specific overrides.
theme_islet <- function(base_size = ISLET_FONT_SIZES$base) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = ISLET_FONT_SIZES$title, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = ISLET_FONT_SIZES$subtitle, color = "#555"),
      axis.title    = ggplot2::element_text(size = ISLET_FONT_SIZES$axis_title),
      axis.text     = ggplot2::element_text(size = ISLET_FONT_SIZES$axis_text),
      legend.title  = ggplot2::element_text(size = ISLET_FONT_SIZES$legend_title, face = "bold"),
      legend.text   = ggplot2::element_text(size = ISLET_FONT_SIZES$legend_text),
      strip.text    = ggplot2::element_text(size = ISLET_FONT_SIZES$strip, face = "bold")
    )
}

# Named list of plotly `list(size = N)` font specs ready to splice into
# layout() calls. Keys: global, title, axis_title, axis_tick, legend_title,
# legend_text, colorbar_title.
plotly_islet_fonts <- function() {
  list(
    global         = list(size = ISLET_FONT_SIZES$base),
    title          = list(size = ISLET_FONT_SIZES$title),
    axis_title     = list(size = ISLET_FONT_SIZES$axis_title),
    axis_tick      = list(size = ISLET_FONT_SIZES$axis_text),
    legend_title   = list(size = ISLET_FONT_SIZES$legend_title),
    legend_text    = list(size = ISLET_FONT_SIZES$legend_text),
    colorbar_title = list(size = ISLET_FONT_SIZES$legend_title),
    hover          = list(size = ISLET_FONT_SIZES$hover),
    annotation     = list(size = ISLET_FONT_SIZES$annotation)
  )
}

# Apply publication-quality PNG download config AND the in-plot text defaults
# that `layout(font = ...)` does NOT cascade to.
#
# Two plotly settings control text the global font never reaches, and both were
# unset on every plot in this app:
#   * `hoverlabel.font` — tooltips render at plotly's ~13px default regardless of
#     `layout(font=)`. This is the single largest "text is tiny" surface, since
#     several plots (Leiden UMAP, tissue scatter) are read primarily via hover.
#   * `uniformtext`     — for bar/heatmap traces plotly silently SHRINKS in-cell
#     text to fit the box, so a `textfont = list(size = 14)` can render at 8px.
#     `minsize` establishes a floor; `mode` decides what happens when it will not
#     fit: "hide" drops the label (right for numeric heatmap values, which
#     overlap illegibly), "show" forces it (right for single-glyph significance
#     stars, which must never disappear).
#
# This function is already piped last on every plotly output, so applying the
# layout here covers them all with no call-site changes.
# Resulting PNG is `width * scale` x `height * scale` px (default 2800x1800).
plotly_islet_config <- function(p, filename_root = "islet_explorer_plot",
                                width = 1400, height = 900, scale = 2,
                                uniform_text_mode = "hide") {
  p <- plotly::layout(p,
    hoverlabel  = list(font = list(size = ISLET_FONT_SIZES$hover)),
    uniformtext = list(minsize = ISLET_FONT_SIZES$dense_min,
                       mode    = uniform_text_mode)
  )
  plotly::config(p,
    toImageButtonOptions = list(
      format   = "png",
      filename = filename_root,
      width    = width,
      height   = height,
      scale    = scale
    ),
    displaylogo = FALSE
  )
}

# ============================================================================
# CSV export helper
# ============================================================================
# Every CSV the app emits should carry a `#`-commented provenance block naming
# the filters that produced it — otherwise two files with different row sets are
# indistinguishable after download. Writing to an open connection (rather than
# writeLines() + write.table(append = TRUE)) also avoids R's "appending column
# names to file" warning.
#
# @param file    Path or connection from a downloadHandler `content` function.
# @param df      Data frame to write. NULL/empty writes the header alone.
# @param header  Character vector of provenance lines (the "# " is added here).
write_islet_csv <- function(file, df, header = character(0)) {
  con <- base::file(file, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  if (length(header) > 0) {
    writeLines(paste("#", header), con)
    writeLines("", con)
  }
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    writeLines("# No data for the current selection", con)
  } else {
    utils::write.csv(df, con, row.names = FALSE)
  }
  invisible(TRUE)
}

# Standard first lines of a provenance block. Callers append their own.
islet_csv_header <- function(title, ...) {
  c(title, paste("Generated:", Sys.time()), ...)
}

# Empty-state placeholder with the standard title font. `plotly_empty()` on its
# own renders its title at plotly's default size, so the nine "No data" messages
# across the app read noticeably smaller than everything around them.
plotly_islet_empty <- function(msg) {
  # Supplying an explicit trace type/mode keeps plotly_build() from guessing
  # (and warning) every time an empty-state widget is rendered.
  plotly::layout(plotly::plotly_empty(type = "scatter", mode = "markers"),
                 title = list(text = msg,
                              font = list(size = ISLET_FONT_SIZES$title)))
}

# ggplotly() silently DROPS `labs(subtitle = ...)`. Several plots put their
# p-value / correction annotations there, so those lines simply never rendered.
# Build a plotly title string that carries the subtitle as a `<sup>` second line
# and pass it to `layout(title = list(text = ...))`.
plotly_title_with_sub <- function(title, subtitle = NULL) {
  if (is.null(subtitle) || !nzchar(trimws(subtitle))) return(title)
  paste0(title, "<br><sup>", subtitle, "</sup>")
}

# Restyle the layout ANNOTATIONS that ggplotly() creates and that no `layout()`
# overlay can reach.
#
# For an unfaceted plot ggplotly keeps axis titles in `layout.xaxis.title`, so
# the usual `layout(xaxis = list(title = list(font = ...)))` works. Add a facet
# and it relocates BOTH the strip labels and the axis titles into
# `layout.annotations` — where that overlay is silently ignored and they render
# at ggplot's converted default (~14.6px). Verified with plotly_build().
#
# Strip labels have no `annotationType`; axis titles carry `"axis"`.
plotly_islet_strip_fonts <- function(gg,
                                     strip_size = ISLET_FONT_SIZES$strip,
                                     axis_title_size = ISLET_FONT_SIZES$axis_title) {
  ann <- gg$x$layout$annotations
  if (is.null(ann) || length(ann) == 0) return(gg)
  for (i in seq_along(ann)) {
    gg$x$layout$annotations[[i]]$font$size <-
      if (identical(ann[[i]]$annotationType, "axis")) axis_title_size else strip_size
  }
  gg
}

# Robust click bridge for plotly scatters used by the islet drill-down.
# On macOS, plotly's own click/drag discrimination can suppress `plotly_click`
# entirely, so the server-side observer never fires (the symptom: clicking a
# point does nothing, no notification). This attaches a native mousedown/mouseup
# listener — which always fires regardless of platform — and, on a click-like
# gesture (movement < 6px), reads the hovered point's `key` ("case_id|islet_key",
# from the ggplot key aes) and pushes it to `input_id` via Shiny.setInputValue.
# The legacy plotly_click handler is kept too, so working platforms are
# unaffected. Pipe after the plotly build:
#   p %>% islet_click_bridge(session$ns("plot_native_click"))
islet_click_bridge <- function(p, input_id) {
  js <- sprintf("
function(el, x){
  // el is the plotly graph div; fall back to a child .js-plotly-plot if not.
  var gd = (el && typeof el.on === 'function') ? el
           : (el && el.querySelector ? el.querySelector('.js-plotly-plot') : el);
  if (!gd) gd = el;
  // onRender re-runs on every plot re-render; attach listeners only once per
  // graph div (plotly reuses the same node) so a click doesn't fire N times.
  if (gd.__isletBridge) return;
  gd.__isletBridge = true;
  var lastKey = null, downKey = null, downX = 0, downY = 0;
  function keyOf(pt){
    try {
      if (pt && pt.data && pt.data.key && pt.data.key[pt.pointNumber] != null)
        return pt.data.key[pt.pointNumber];
    } catch(e){}
    return null;
  }
  // plotly's live hover state — robust even if the plotly_hover event timing
  // is off, since at mousedown the cursor is over the target point.
  function hoverKey(){
    try { if (gd && gd._hoverdata && gd._hoverdata.length) return keyOf(gd._hoverdata[0]); } catch(e){}
    return null;
  }
  if (gd && typeof gd.on === 'function') {
    gd.on('plotly_hover', function(d){
      if (d && d.points && d.points.length){ var k = keyOf(d.points[0]); if (k != null) lastKey = k; }
    });
    gd.on('plotly_unhover', function(){ lastKey = null; });
    gd.on('plotly_click', function(){ console.log('[ISLET] plotly_click fired'); });
  }
  (gd || el).addEventListener('mousedown', function(e){ downX = e.clientX; downY = e.clientY; downKey = lastKey || hoverKey(); }, true);
  (gd || el).addEventListener('mouseup', function(e){
    if (Math.abs(e.clientX - downX) > 5 || Math.abs(e.clientY - downY) > 5) return;
    var k = downKey || lastKey || hoverKey();
    if (k == null) return;
    console.log('[ISLET] native click -> ' + k);
    Shiny.setInputValue('%s', {key: k, nonce: (Date.now() + Math.random())}, {priority: 'event'});
  }, true);
}", input_id)
  htmlwidgets::onRender(p, js)
}

# Standard PNG download dimensions for ggplot ggsave() in downloadHandlers.
# Default: 10x7 in at 300 dpi -> 3000x2100 px. Heatmaps use wider 12x8.
ISLET_PNG_DIMS  <- list(width = 10, height = 7, dpi = 300, units = "in")
ISLET_HEATMAP_PNG_DIMS <- list(width = 12, height = 8, dpi = 300, units = "in")

islet_png_filename <- function(stem) {
  paste0(stem, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
}

# ============================================================================
# Phenotype rules (May 2026) — marker -> celltype definitions surfaced in UI
# ============================================================================
# Source of truth: data/phenotype_rules.csv (parsed once at startup).
# Column-1 is the parent rule, column-2 the phenotype, columns 3+ are markers
# with cell values:
#   "pos"    -> required positive
#   "anypos" -> at least one of the marker(s) in the same `anypos` group
#               (the parser groups all `anypos` cells per phenotype together)
#   "allneg" -> required negative
#
# PHENOTYPE_COLORS (drilldown_helpers.R) uses 21 names; a few don't have rows
# in the CSV (residuals from the upstream single-cell pipeline). Those get
# hand-written descriptions in PHENOTYPE_RESIDUALS below. The alias map
# normalises PHENOTYPE_COLORS keys -> CSV keys for the four known mismatches.

PHENOTYPE_RULES_PATH <- data_path("phenotype_rules.csv")

# PHENOTYPE_COLORS-form -> CSV-form (and back via reverse lookup). Used to
# resolve a UI label to its rule row.
PHENOTYPE_ALIASES <- c(
  "CD8a Tcell"   = "CD8 T cell",
  "CD4 Tcell"    = "CD4 T cell",
  "Blood Vessel" = "Blood vessel",
  "SMA+"         = "Smooth Muscle"
)

# CSV/rule-form -> canonical display label (the PHENOTYPE_COLORS / data
# `phenotype`-column form used everywhere else in the app: drill-down legend,
# composition chart, Plot and Spatial selectors). The reference modal reads the
# CSV, so without this map it would show "Smooth Muscle" / "Blood vessel" /
# "CD8 T cell" while the rest of the app shows "SMA+" / "Blood Vessel" /
# "CD8a Tcell" for the same phenotype. Inverse of PHENOTYPE_ALIASES.
PHENOTYPE_ALIASES_REV <- stats::setNames(names(PHENOTYPE_ALIASES), unname(PHENOTYPE_ALIASES))

# Resolve any phenotype name (CSV-form or already-canonical) to the app's
# canonical display label. Falls back to the input unchanged.
phenotype_display_label <- function(name) {
  if (is.null(name) || is.na(name) || !nzchar(name)) return(name)
  if (name %in% names(PHENOTYPE_ALIASES_REV)) return(unname(PHENOTYPE_ALIASES_REV[[name]]))
  name
}

# Phenotypes with no CSV row — residuals from the upstream pipeline.
# Each entry: list(parent, pos, neg, anypos_groups, description).
PHENOTYPE_RESIDUALS <- list(
  `ECAD+` = list(
    parent = "Epithelial", pos = character(), neg = character(),
    anypos_groups = list(),
    description = "Epithelial cells (E-cadherin+); residual ductal/exocrine bucket from the upstream phenotyping pipeline. Not gated by the CSV rules."
  ),
  APCs = list(
    parent = "Immune", pos = character(), neg = character(),
    anypos_groups = list(),
    description = "Antigen-presenting cells aggregate (HLADR+ myeloid lineage). Used as an aggregate immune category in neighbourhood metrics; not a CSV-rule phenotype."
  ),
  Structural = list(
    parent = "all", pos = character(), neg = character(),
    anypos_groups = list(),
    description = "Residual structural cells negative for all endocrine/immune markers; bucket for fibroblast/stroma-like cells."
  ),
  Unknown = list(
    parent = "all", pos = character(), neg = character(),
    anypos_groups = list(),
    description = "Default fallback for cells that match no rule above."
  )
)

# Parse the CSV into a named list of rule records keyed by phenotype name.
# Returns NULL if the file is missing — callers should guard.
.parse_phenotype_rules_csv <- function(path) {
  if (!file.exists(path)) {
    message("[PHENOTYPE] Rules CSV not found at ", path, " — modal/hints disabled")
    return(NULL)
  }
  raw <- utils::read.csv(path, header = FALSE, stringsAsFactors = FALSE,
                         na.strings = c("", "NA"), check.names = FALSE)
  # Row 1 is the header; columns 3..N are marker names.
  markers <- as.character(raw[1, 3:ncol(raw)])
  markers <- markers[!is.na(markers) & nzchar(markers)]
  data_rows <- raw[-1, , drop = FALSE]
  rules <- list()
  for (i in seq_len(nrow(data_rows))) {
    pheno <- as.character(data_rows[i, 2])
    if (is.na(pheno) || !nzchar(pheno)) next
    parent <- as.character(data_rows[i, 1])
    parent <- if (is.na(parent) || !nzchar(parent)) "all" else parent
    vals <- as.character(data_rows[i, 3:(2 + length(markers))])
    names(vals) <- markers
    pos <- names(vals)[!is.na(vals) & vals == "pos"]
    neg <- names(vals)[!is.na(vals) & vals == "allneg"]
    anypos_markers <- names(vals)[!is.na(vals) & vals == "anypos"]
    # CSV puts all `anypos` markers on one row as a single group.
    anypos_groups <- if (length(anypos_markers) > 0) list(anypos_markers) else list()
    rules[[pheno]] <- list(
      parent        = parent,
      pos           = pos,
      neg           = neg,
      anypos_groups = anypos_groups,
      description   = NA_character_
    )
  }
  rules
}

PHENOTYPE_RULES_CSV_RAW <- .parse_phenotype_rules_csv(PHENOTYPE_RULES_PATH)

# Compose hierarchical rules: a sub-rule (e.g. "CD8 T cell") inherits the
# positive/negative requirements of its parent chain ("T cell" -> "Immune").
# We walk the parent chain until we hit "all" (root) or a missing parent.
.compose_phenotype_rules <- function(raw) {
  if (is.null(raw)) return(NULL)
  composed <- list()
  compose <- function(name, seen = character()) {
    if (!is.null(composed[[name]])) return(composed[[name]])
    if (name %in% seen) return(raw[[name]])  # cycle guard
    if (!name %in% names(raw)) return(NULL)
    r <- raw[[name]]
    if (!is.null(r$parent) && nzchar(r$parent) && r$parent != "all"
        && r$parent %in% names(raw)) {
      p <- compose(r$parent, c(seen, name))
      if (!is.null(p)) {
        r$pos           <- unique(c(p$pos, r$pos))
        r$neg           <- unique(c(p$neg, r$neg))
        r$anypos_groups <- c(p$anypos_groups, r$anypos_groups)
      }
    }
    # Cleanup: drop anypos groups already satisfied by a required positive.
    r$anypos_groups <- Filter(function(grp) length(intersect(grp, r$pos)) == 0,
                              r$anypos_groups)
    # Deduplicate identical groups.
    if (length(r$anypos_groups) > 1) {
      keys <- vapply(r$anypos_groups, function(g) paste(sort(g), collapse = "|"), character(1))
      r$anypos_groups <- r$anypos_groups[!duplicated(keys)]
    }
    # Drop any group that is a strict SUPERSET of another group in the same rule
    # (the smaller group is the stricter, sufficient constraint).
    if (length(r$anypos_groups) > 1) {
      keep <- rep(TRUE, length(r$anypos_groups))
      for (i in seq_along(r$anypos_groups)) {
        for (j in seq_along(r$anypos_groups)) {
          if (i == j) next
          a <- r$anypos_groups[[i]]; b <- r$anypos_groups[[j]]
          if (length(b) < length(a) && all(b %in% a)) { keep[i] <- FALSE; break }
        }
      }
      r$anypos_groups <- r$anypos_groups[keep]
    }
    r$neg <- setdiff(r$neg, r$pos)
    composed[[name]] <<- r
    r
  }
  for (nm in names(raw)) compose(nm)
  composed
}

PHENOTYPE_RULES_CSV <- .compose_phenotype_rules(PHENOTYPE_RULES_CSV_RAW)

# Resolve a phenotype name (in either CSV or PHENOTYPE_COLORS form) to its rule.
# Returns NULL if the name has no rule and no residual entry.
get_phenotype_rule <- function(name) {
  if (is.null(name) || is.na(name) || !nzchar(name)) return(NULL)
  # Try residuals first (small set, includes "Unknown" fallback).
  if (name %in% names(PHENOTYPE_RESIDUALS)) return(PHENOTYPE_RESIDUALS[[name]])
  # Apply alias (PHENOTYPE_COLORS -> CSV).
  csv_name <- if (name %in% names(PHENOTYPE_ALIASES)) unname(PHENOTYPE_ALIASES[[name]]) else name
  if (!is.null(PHENOTYPE_RULES_CSV) && csv_name %in% names(PHENOTYPE_RULES_CSV)) {
    return(PHENOTYPE_RULES_CSV[[csv_name]])
  }
  NULL
}

# Format the rule for display. Styles:
#   "inline" -> plain text, e.g. "CD3e+, CD8a+  —  CD4−, CD20−"
#   "html"   -> HTML with strong positives and muted negatives
#   "plain"  -> minimal text, safe for plot subtitles
format_phenotype_rule <- function(name, style = c("inline", "html", "plain")) {
  style <- match.arg(style)
  rule <- get_phenotype_rule(name)
  if (is.null(rule)) {
    msg <- paste0("No marker rule recorded for “", name, "”.")
    if (style == "html") return(htmltools::HTML(paste0("<em>", msg, "</em>")))
    return(msg)
  }
  if (!is.na(rule$description) && nzchar(rule$description)) {
    # Residual phenotype — just return the description.
    if (style == "html") return(htmltools::HTML(paste0("<em>", htmltools::htmlEscape(rule$description), "</em>")))
    return(rule$description)
  }
  # Build positive-side string
  pos_parts <- character()
  for (m in rule$pos) pos_parts <- c(pos_parts, paste0(m, "+"))
  for (grp in rule$anypos_groups) {
    if (length(grp) == 1) {
      pos_parts <- c(pos_parts, paste0(grp, "+"))
    } else {
      pos_parts <- c(pos_parts, paste0("any of ", paste(paste0(grp, "+"), collapse = " / ")))
    }
  }
  neg_parts <- if (length(rule$neg) > 0) paste0(rule$neg, "−") else character()
  if (style == "html") {
    pos_html <- if (length(pos_parts) > 0)
      paste0("<strong>", paste(htmltools::htmlEscape(pos_parts), collapse = ", "), "</strong>")
    else "<span style='color:#888;'>(no positive requirement)</span>"
    neg_html <- if (length(neg_parts) > 0)
      paste0("<span style='color:#666;'>", paste(htmltools::htmlEscape(neg_parts), collapse = ", "), "</span>")
    else ""
    sep <- if (nzchar(neg_html)) "  —  " else ""
    return(htmltools::HTML(paste0(pos_html, sep, neg_html)))
  }
  pos_txt <- if (length(pos_parts) > 0) paste(pos_parts, collapse = ", ") else "(no positive requirement)"
  neg_txt <- if (length(neg_parts) > 0) paste(neg_parts, collapse = ", ") else ""
  if (nzchar(neg_txt)) paste0(pos_txt, "  —  ", neg_txt) else pos_txt
}

# Order in which phenotypes are shown in the reference modal — endocrine first,
# then immune (T cell subtypes grouped), then structural/endothelial, residuals last.
PHENOTYPE_DISPLAY_ORDER <- c(
  "Beta cell", "Alpha cell", "Delta cell", "Endocrine",
  "Acinar", "Ductal", "Exocrine",
  "T cell", "CD8 T cell", "CD4 T cell", "B cell",
  "Macrophage", "Dendritic cell", "Immune", "APCs",
  "Lymphatic", "Blood vessel", "Endothelial",
  "Smooth Muscle", "Neural", "ECAD+", "Structural",
  "Unknown"
)

# Build a tidy data.frame for the modal table.
# Returns columns: Phenotype, Category, Positive markers, Negative markers, Notes.
phenotype_rules_table_df <- function() {
  csv_names <- if (is.null(PHENOTYPE_RULES_CSV)) character() else names(PHENOTYPE_RULES_CSV)
  all_names <- c(csv_names, names(PHENOTYPE_RESIDUALS))
  all_names <- unique(all_names)
  # Order: known order first, then anything else alphabetical.
  ordered <- c(intersect(PHENOTYPE_DISPLAY_ORDER, all_names),
               sort(setdiff(all_names, PHENOTYPE_DISPLAY_ORDER)))
  rows <- lapply(ordered, function(nm) {
    rule <- get_phenotype_rule(nm)
    # Display the app's canonical label (SMA+, Blood Vessel, CD8a Tcell, …),
    # not the CSV-form name, so the modal matches every other view.
    disp <- phenotype_display_label(nm)
    if (is.null(rule)) {
      return(data.frame(Phenotype = disp, Category = "?", Positive = "",
                        Negative = "", Notes = "",
                        stringsAsFactors = FALSE))
    }
    pos_parts <- character()
    for (m in rule$pos) pos_parts <- c(pos_parts, paste0(m, "+"))
    for (grp in rule$anypos_groups) {
      pos_parts <- c(pos_parts,
                     if (length(grp) == 1) paste0(grp, "+")
                     else paste0("any of ", paste(paste0(grp, "+"), collapse = " / ")))
    }
    neg_parts <- if (length(rule$neg) > 0) paste0(rule$neg, "−") else character()
    data.frame(
      Phenotype = disp,
      Category  = rule$parent,
      Positive  = if (length(pos_parts) > 0) paste(pos_parts, collapse = ", ") else "—",
      Negative  = if (length(neg_parts) > 0) paste(neg_parts, collapse = ", ") else "—",
      Notes     = if (!is.na(rule$description)) rule$description else "",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Modal builder. Pure shiny — no DT dependency. Pass to showModal() from any tab's
# Phenotype-Rules button handler.
phenotype_rules_modal_ui <- function() {
  df <- phenotype_rules_table_df()
  header_row <- shiny::tags$tr(
    shiny::tags$th("Phenotype",          style = "text-align:left; padding:6px 10px;"),
    shiny::tags$th("Category",           style = "text-align:left; padding:6px 10px;"),
    shiny::tags$th("Positive markers",   style = "text-align:left; padding:6px 10px;"),
    shiny::tags$th("Negative markers",   style = "text-align:left; padding:6px 10px;"),
    shiny::tags$th("Notes",              style = "text-align:left; padding:6px 10px;")
  )
  body_rows <- lapply(seq_len(nrow(df)), function(i) {
    shiny::tags$tr(
      shiny::tags$td(shiny::tags$strong(df$Phenotype[i]), style = "padding:6px 10px;"),
      shiny::tags$td(df$Category[i],                       style = "padding:6px 10px; color:#555;"),
      shiny::tags$td(df$Positive[i],                       style = "padding:6px 10px;"),
      shiny::tags$td(df$Negative[i],                       style = "padding:6px 10px; color:#666;"),
      shiny::tags$td(df$Notes[i],                          style = "padding:6px 10px; color:#444; font-size: 13px;")
    )
  })
  shiny::modalDialog(
    title = shiny::tagList(shiny::icon("dna"), "Cell Phenotype Rules — Marker Definitions"),
    size = "l", easyClose = TRUE, fade = TRUE,
    footer = shiny::modalButton("Close"),
    shiny::tags$p(
      "Every cell in the app is classified into one of these phenotypes by the upstream ",
      shiny::tags$code("data/phenotype_rules.csv"), " gating rules. A ",
      shiny::tags$code("+"),
      " means the marker is required to be positive; a ",
      shiny::tags$code("−"),
      " means it must be negative. “any of X+ / Y+” means at least one of the listed markers must be positive."
    ),
    shiny::tags$div(style = "max-height: 60vh; overflow-y: auto; border-top: 1px solid #ddd;",
      shiny::tags$table(
        class = "table table-striped",
        style = "width: 100%; font-size: 14px; margin: 0;",
        shiny::tags$thead(header_row),
        shiny::tags$tbody(body_rows)
      )
    ),
    shiny::tags$p(style = "margin-top: 10px; font-size: 12px; color: #777;",
      "Rule source: ",
      shiny::tags$code("data/phenotype_rules.csv"),
      " (15 markers: INS, GCG, SST, SMA, CD20, CD3e, CD68, CD163, CD8a, CD4, HLADR, PDPN, CD34, CD31, B3TUBB)."
    ),
    shiny::tags$p(style = "margin-top: 6px; font-size: 12px; color: #777;",
      shiny::tags$strong("Data vs. rules: "),
      "this table is the gating spec. A few rule phenotypes have no cells in the ",
      "shipped dataset because the upstream scVI labeling merged or renamed them ",
      "(e.g. Dendritic cell folds into the ", shiny::tags$code("APCs"),
      "/Immune aggregates; Exocrine splits into Acinar/Ductal). ",
      shiny::tags$code("APCs"),
      " is an HLADR+ myeloid aggregate created upstream, not itself a ",
      shiny::tags$code("phenotype_rules.csv"), " row."
    ),
    shiny::tags$p(style = "margin-top: 6px; font-size: 12px; color: #777;",
      shiny::tags$strong("Panel note: "),
      "the original Rules1 template also listed ", shiny::tags$code("CD99"),
      " (Endocrine) and ", shiny::tags$code("CD56"),
      " (Neural), but neither antibody is in this 31-plex panel, so the ",
      "phenotyping step dropped them. The effective gates were therefore ",
      shiny::tags$code("Endocrine = any of INS/GCG/SST+"), " and ",
      shiny::tags$code("Neural = B3TUBB+"), " — reflected in the table above."
    )
  )
}

# Resolve a sanitised data-column name (e.g. "prop_CD8a_Tcell", "peri_prop_Beta_cell",
# "min_dist_Macrophage", "enrich_z_B_cell") back to a phenotype label. Returns
# NULL/the input unchanged if no known prefix is recognised.
extract_phenotype_from_column <- function(col) {
  if (is.null(col) || is.na(col) || !nzchar(col)) return(NULL)
  for (pre in c("combined_prop_", "peri_prop_", "peri_count_",
                "enrich_z_", "min_dist_", "prop_")) {
    if (startsWith(col, pre)) {
      stem <- sub(paste0("^", pre), "", col)
      return(gsub("_", " ", gsub("plus", "+", stem, fixed = TRUE), fixed = TRUE))
    }
  }
  col
}

# Human-readable label for a raw data-column / selection encoding, for use in
# plot axes, card titles and CSV headers.
#
# The Statistics tab took its feature name straight from the Plot module's
# `selection_used` attribute, so raw encodings leaked onto the page as if they
# were the feature's name: "prop_Beta cell", "peri_prop_CD8a_Tcell",
# "Capillary|islet_core". This resolves the two encodings the app actually uses:
#   - `<prefix>_<phenotype>`  (prop_ / peri_prop_ / combined_prop_ / ...)
#   - `<class>|<region>`      (Microenvironment targets)
# Anything unrecognised is returned unchanged, so this is always safe to apply.
feature_display_label <- function(col) {
  if (is.null(col) || length(col) != 1 || is.na(col) || !nzchar(col)) return(col)

  # Microenvironment encoding: "Capillary|islet_core"
  if (grepl("|", col, fixed = TRUE)) {
    parts  <- strsplit(col, "|", fixed = TRUE)[[1]]
    region <- switch(parts[2],
                     islet_core  = "islet core",
                     islet_band  = "peri-islet",
                     islet_union = "core+peri",
                     gsub("_", " ", parts[2], fixed = TRUE))
    return(paste0(parts[1], " (", region, ")"))
  }

  region <- if (startsWith(col, "combined_prop_")) "core+peri"
            else if (startsWith(col, "peri_prop_") || startsWith(col, "peri_count_")) "peri-islet"
            else if (startsWith(col, "prop_") || startsWith(col, "count_")) "core"
            else NULL

  measure <- if (grepl("^(combined_prop_|peri_prop_|prop_)", col)) "% of cells"
             else if (grepl("^(peri_count_|count_)", col)) "cells"
             else if (startsWith(col, "enrich_z_")) "enrichment z"
             else if (startsWith(col, "min_dist_")) "min distance, um"
             else NULL

  pheno <- extract_phenotype_from_column(col)
  if (is.null(pheno) || identical(pheno, col)) return(col)   # no prefix matched

  qual <- paste(c(region, measure), collapse = ", ")
  if (nzchar(qual)) paste0(pheno, " (", qual, ")") else pheno
}

# UI button that opens the Phenotype Rules modal. The caller chooses the input
# id (one per tab to avoid duplicate-DOM-id collisions across simultaneously
# rendered tabPanels) and registers a matching `observeEvent` that calls
# `showModal(phenotype_rules_modal_ui())`. Style mirrors the existing "About
# this tab" info banner.
phenotype_rules_button <- function(input_id,
                                   label = "Phenotype rules — marker definitions") {
  shiny::tags$div(
    style = "margin-bottom: 12px; padding: 8px 14px; background-color: #eef4fb; border: 1px solid #c8d8e9; border-radius: 6px; display: inline-block;",
    shiny::actionLink(input_id,
      shiny::tagList(shiny::icon("dna", style = "color:#4477AA;"),
                     shiny::tags$span(style = "margin-left:6px; color:#2a4d75; font-weight:600;",
                                      label)),
      style = "text-decoration: none;")
  )
}

# UI helper: small muted hint block describing the currently selected phenotype.
# Pass a literal phenotype name. The output is a tagList suitable for embedding
# directly under a selectInput. `show_modal_link` controls the trailing
# "see all rules" actionLink; set to NULL to omit.
phenotype_hint_ui <- function(name, prefix = "Markers:",
                              show_modal_link = "show_phenotype_rules_inline") {
  if (is.null(name) || is.na(name) || !nzchar(name)) return(NULL)
  rule_html <- format_phenotype_rule(name, "html")
  link <- if (!is.null(show_modal_link)) {
    shiny::tagList(
      shiny::tags$span(style = "color:#888;", " · "),
      shiny::actionLink(show_modal_link, "see all rules →",
                        style = "font-size: 12px; color: #4477AA;")
    )
  }
  shiny::tags$div(
    style = "margin-top: -6px; margin-bottom: 8px; padding: 6px 10px; background-color: #f7f9fc; border-left: 3px solid #4477AA; border-radius: 3px; font-size: 13px; color: #333; line-height: 1.4;",
    shiny::tags$strong(style = "color:#4477AA;", paste0(prefix, " ")),
    shiny::tags$span(rule_html),
    link
  )
}

# ============================================================================
# Retired phenotypes (Jul 2026)
# ============================================================================
# Nine labels are excluded from the app everywhere a phenotype is listed,
# coloured, selected or tested. Two reasons, and they are the same reason twice:
#
#   Endocrine, Immune, T cell, Endothelial   parents of other labels in the set,
#                                            so they are gate-failure residuals
#                                            ("immune, not resolvable further")
#                                            whose abundance tracks gate quality
#   Structural, Unknown                      `-rest` catch-alls, the same thing
#                                            without children to give them away
#   Acinar, ECAD+, APCs                      retired with them
#
# Excluded rather than deleted from PHENOTYPE_COLORS: the DATA still carries
# these labels, and every lookup in this app is by exact name, so a palette
# without them would return NA and silently grey out real cells. Keeping the
# colours and filtering at the point of USE means a stray row is dropped, not
# mis-drawn.
#
# `phenotype_keep()` is the single gate. Anything enumerating phenotypes for a
# selector, a legend, a family or a test goes through it.
PHENOTYPE_RETIRED <- c("Endocrine", "Acinar", "Endothelial", "T cell",
                       "Immune", "Unknown", "ECAD+", "Structural", "APCs")

#' Drop retired phenotypes from a character vector of labels.
phenotype_keep <- function(x) setdiff(x, PHENOTYPE_RETIRED)

#' TRUE for rows whose phenotype is still in use. NA-safe.
phenotype_is_kept <- function(x) !(as.character(x) %in% PHENOTYPE_RETIRED)

#' Drop retired phenotypes from a data frame with a phenotype column.
phenotype_filter_rows <- function(df, col = "phenotype") {
  if (is.null(df) || !is.data.frame(df) || !(col %in% names(df))) return(df)
  df[phenotype_is_kept(df[[col]]), , drop = FALSE]
}
