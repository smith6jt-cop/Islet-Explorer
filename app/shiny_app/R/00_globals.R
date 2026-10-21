# Debug modes are a deployment error in production. Keep this check before
# package loading and other startup work so a misconfigured worker fails closed
# without entering either diagnostic path.
DEBUG_CREDS <- islet_flag_enabled(Sys.getenv("DEBUG_CREDENTIALS", "0"))
VIEWER_DEBUG_ENABLED <- islet_flag_enabled(Sys.getenv("VIEWER_DEBUG", "0"))
if (identical(tolower(trimws(Sys.getenv("ISLET_RUNTIME_MODE", "development"))),
              "production") &&
    (DEBUG_CREDS || VIEWER_DEBUG_ENABLED)) {
  islet_log_event("runtime_profile_rejected")
  stop("Production does not permit credential or viewer debug modes.", call. = FALSE)
}

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

# Bright is the default selection (app.R, the three palette selectors); no choice is labelled "(default)".
DONOR_PALETTES <- list(
  "Bright"             = DONOR_COLORS_BRIGHT,
  "Paul Tol"           = DONOR_COLORS,
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
    islet_log_event("sf_ready")
  } else {
    islet_log_event("sf_unavailable")
  }
})

# Check for httr2 package (required for AI assistant)
httr2_available <- requireNamespace("httr2", quietly = TRUE)
if (!httr2_available) {
  islet_log_event("httr2_unavailable")
}

# Load anndata package for trajectory analysis (install with BiocManager::install('anndata'))
suppressPackageStartupMessages({
  if (!requireNamespace("anndata", quietly = TRUE)) {
    islet_log_event("anndata_unavailable")
  }
})

# Load vitessceR for embedded image viewer (GitHub v0.1.0 required)
VITESSCE_AVAILABLE <- FALSE
suppressPackageStartupMessages({
  if (!requireNamespace("vitessceR", quietly = TRUE)) {
    islet_log_event("vitessce_unavailable")
  } else {
    # Verify we have the GitHub version with required features
    tryCatch({
      pkg_version <- packageVersion("vitessceR")
      if (pkg_version < "0.1.0") {
        islet_log_event("vitessce_version_unsupported")
      }
      test_vc <- vitessceR::VitessceConfig$new(schema_version = "1.0.16", name = "Test")
      islet_log_event("vitessce_ready")
      VITESSCE_AVAILABLE <<- TRUE
    }, error = function(e) {
      islet_log_event("vitessce_check_failed")
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
        islet_log_event("vitessce_routes_ready")
      }, error = function(e) {
        islet_log_event("vitessce_routes_failed")
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
  islet_log_event("app_loaded")
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
      islet_log_event("app_load_flag_rejected")
    } else {
      if (!dir.exists(dirname(flag_path))) {
        dir.create(dirname(flag_path), recursive = TRUE, showWarnings = FALSE)
      }
      writeLines(paste0('loaded ', Sys.time(), ' version=', APP_VERSION), flag_path)
    }
  }
}, silent = TRUE)

# Runtime paths. Production points ISLET_DATA_ROOT at its verified release;
# development uses one selected root, defaulting to the repository data dir.
.repo_data_root <- file.path("..", "..", "data")
data_root <- Sys.getenv("ISLET_DATA_ROOT", unset = .repo_data_root)
data_root <- tryCatch(normalizePath(data_root, mustWork = FALSE),
                      error = function(e) data_root)
data_path <- function(...) file.path(data_root, ...)

master_path <- data_path("master_results.xlsx")
h5ad_path   <- data_path("islet_explorer.h5ad")
project_root <- tryCatch(normalizePath(file.path("..", ".."), mustWork = FALSE), error = function(e) NULL)

# Donors removed from the app and all analysis by owner direction. Keep in step
# with scripts/excluded_donors.csv (a pipeline test checks that they match).
# Donor 78 (HDL078 / HPAP-078) was removed on 2026-10-15.
ISLET_EXCLUDED_DONORS <- c("78")

#' Bare numeric donor ID: drops an HPAP/HDL prefix, separators and leading
#' zeros, so "0078", "HPAP-078" and "HDL078" all compare as "78".
islet_canonical_donor_id <- function(x) {
  x <- trimws(as.character(x))
  sub("^(?:HPAP|HDL)?[-_ ]?0*([0-9]+)(?:\\.0+)?$", "\\1", x, ignore.case = TRUE, perl = TRUE)
}

#' The supplied donor IDs that name an excluded donor.
islet_excluded_donors_in <- function(ids) {
  ids <- unique(as.character(ids[!is.na(ids)]))
  ids[islet_canonical_donor_id(ids) %in% ISLET_EXCLUDED_DONORS]
}

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
    islet_log_event("python_anndata_unavailable")
  }

  if (nzchar(candidate_python)) {
    Sys.setenv(RETICULATE_PYTHON = candidate_python)
    islet_log_event("python_configured")
  } else {
    islet_log_event("python_unavailable")
  }
}

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
# The pseudotime is fitted on every endocrine object with at least TRAJ_MIN_CELLS core cells, whatever its
# class (Islet / Cluster); the app-wide class selection does not apply to it. Barlow et al. 2024 fitted
# islets of ten or more cells. On the uniform root that is 1,265 objects (1,242 Islets, 23 Clusters), with
# rho(pseudotime, beta-cell ratio) -0.70; fitted on every object it was -0.20, because one- to few-cell objects
# carry almost no composition to order. The user-set cells-per-islet range slider was removed on 2026-10-17;
# the fixed cutoff replaced the class selection on 2026-10-18 (owner: "Use 10 or more cells as the cutoff
# instead of a category").
TRAJ_MIN_CELLS <- 10L


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
# URL of one of the app's own www/ files, versioned by its content. shiny-server sends them with Last-Modified
# and no Cache-Control, so a browser may keep using an old copy after the file changes (2026-10-19: both header
# logos showed because a cached islet-theme.css lacked the rule that hides one).
islet_asset <- function(path) {
  file <- file.path("www", path)
  if (!file.exists(file)) return(path)
  paste0(path, "?v=", substr(unname(tools::md5sum(file)), 1, 10))
}

# ---- Colour theme (2026-10-18: dark by default, after VS Code's Dark Modern) ----
# The page theme is the <html data-theme> attribute (www/islet-theme.css / islet-theme.js), mirrored to the server
# as input$app_theme. theme_islet() and plotly_islet_config() read it, so every plot drawn inside a render
# depends on it and redraws when the reviewer switches theme. Downloads stay white (islet_ggsave(), and the
# plotly camera button re-draws the figure in light colours).
ISLET_THEME_PALETTE <- list(
  dark  = list(paper = "#1f1f1f", ink = "#cccccc", muted = "#9d9d9d", grid = "#333333", zero = "#555555",
               line = "#6b6b6b"),
  light = list(paper = "#ffffff", ink = "#000000", muted = "#555555", grid = "#ebebeb", zero = "#bdbdbd",
               line = "#444444"))

#' The theme's foreground colour, for lines and text a plot draws in a fixed colour (medians, reference lines).
islet_ink <- function() ISLET_THEME_PALETTE[[islet_theme_mode()]]$ink

#' "dark" or "light": the reviewer's theme for the current session (dark when there is no session or input).
islet_theme_mode <- function() {
  forced <- getOption("islet.theme")
  if (!is.null(forced)) return(if (identical(forced, "light")) "light" else "dark")
  session <- shiny::getDefaultReactiveDomain()
  # Inside a module the domain is the module's proxy, whose input is namespaced; the theme is an app-level input.
  if (!is.null(session) && is.function(session$rootScope)) session <- session$rootScope()
  mode <- if (is.null(session)) NULL else tryCatch(session$input$app_theme, error = function(e) NULL)
  if (identical(mode, "light")) "light" else "dark"
}

theme_islet <- function(base_size = ISLET_FONT_SIZES$base, mode = islet_theme_mode()) {
  pal <- ISLET_THEME_PALETTE[[mode]]
  # ggplot2 4: ink / paper also set the default colours of geoms, so text and lines drawn without an explicit
  # colour follow the theme.
  ggplot2::theme_minimal(base_size = base_size, ink = pal$ink, paper = pal$paper) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = pal$paper, colour = NA),
      plot.title    = ggplot2::element_text(size = ISLET_FONT_SIZES$title, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = ISLET_FONT_SIZES$subtitle, color = pal$muted),
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
  if (identical(islet_theme_mode(), "dark")) p <- plotly_islet_dark(p)
  # The camera button is replaced by one that exports on a white background whatever the screen theme
  # (www/islet-theme.js isletDownloadLightPng); same file name, size and scale as before.
  opts <- jsonlite::toJSON(list(format = "png", filename = filename_root, width = width, height = height,
                                scale = scale), auto_unbox = TRUE)
  plotly::config(p,
    toImageButtonOptions = list(
      format   = "png",
      filename = filename_root,
      width    = width,
      height   = height,
      scale    = scale
    ),
    modeBarButtonsToRemove = list("toImage"),
    modeBarButtonsToAdd = list(list(
      name  = "Download plot as a PNG (white background)",
      icon  = htmlwidgets::JS("Plotly.Icons.camera"),
      click = htmlwidgets::JS(sprintf("function(gd) { window.isletDownloadLightPng(gd, %s); }", opts)))),
    displaylogo = FALSE
  )
}

#' Recolour a built plotly figure for the dark theme: backgrounds, text, every axis, legend, annotations and
#' colour bars. Trace colours are left alone.
plotly_islet_dark <- function(p) {
  pal <- ISLET_THEME_PALETTE$dark
  p <- plotly::plotly_build(p)
  L <- p$x$layout
  L$paper_bgcolor <- pal$paper
  L$plot_bgcolor <- pal$paper
  L$font$color <- pal$ink
  if (is.list(L$title)) L$title$font$color <- pal$ink
  for (k in grep("^[xy]axis[0-9]*$", names(L), value = TRUE)) {
    ax <- L[[k]]
    ax$gridcolor <- pal$grid; ax$zerolinecolor <- pal$zero; ax$linecolor <- pal$line; ax$tickcolor <- pal$line
    ax$tickfont$color <- pal$ink
    if (is.list(ax$title)) ax$title$font$color <- pal$ink
    L[[k]] <- ax
  }
  L$legend$bgcolor <- "rgba(0,0,0,0)"
  L$legend$font$color <- pal$ink
  if (is.list(L$legend$title)) L$legend$title$font$color <- pal$ink
  if (length(L$annotations)) {
    L$annotations <- lapply(L$annotations, function(a) { a$font$color <- pal$ink; a })
  }
  L$hoverlabel$bordercolor <- pal$line
  p$x$layout <- L
  p$x$data <- lapply(p$x$data, function(tr) {
    recolour <- function(cb) { cb$tickfont$color <- pal$ink; if (is.list(cb$title)) cb$title$font$color <- pal$ink; cb }
    if (!is.null(tr$colorbar)) tr$colorbar <- recolour(tr$colorbar)
    if (!is.null(tr$marker$colorbar)) tr$marker$colorbar <- recolour(tr$marker$colorbar)
    tr
  })
  p
}

#' ggsave() for downloadHandlers: publication PNGs stay on a white background whatever the screen theme, by
#' overriding only the theme's colours (the plot's own theme tweaks are kept).
islet_ggsave <- function(file, plot, ...) {
  pal <- ISLET_THEME_PALETTE$light
  light <- ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = pal$paper, colour = NA),
    panel.background = ggplot2::element_rect(fill = pal$paper, colour = NA),
    legend.background = ggplot2::element_rect(fill = pal$paper, colour = NA),
    legend.key = ggplot2::element_rect(fill = pal$paper, colour = NA),
    text = ggplot2::element_text(colour = pal$ink),
    axis.text = ggplot2::element_text(colour = "#4d4d4d"),
    plot.subtitle = ggplot2::element_text(colour = pal$muted),
    panel.grid = ggplot2::element_line(colour = pal$grid),
    strip.text = ggplot2::element_text(colour = pal$ink),
    geom = ggplot2::element_geom(ink = pal$ink, paper = pal$paper))
  ggplot2::ggsave(file, plot = plot + light, bg = "white", ...)
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

# Islet-level UMAP markers, shared by the Spatial and Trajectory UMAPs. Colouring a plotly scatter by a group makes
# one trace per group, drawn in turn, so the last group (T1D) was painted over the others wherever they met and hid
# where each disease state sat. Here every point is one trace in a fixed interleaved order (golden-ratio stride over
# the rows, which come grouped by donor), and an empty trace per group carries the legend, whose clicks are turned
# off because they could only hide the empty trace. The thin outline in the background colour keeps touching points
# apart; the white-background PNG export swaps it to white (www/islet-theme.js).
# Size 4 px: on the v7 layout (min_dist 3, spread 3) in a Spatial UMAP, 7 % of points lie within half a marker of a
# neighbour (13 % at 5 px; 42 % at 4 px on v6's layout).
ISLET_UMAP_MARKER <- list(size = 4, opacity = 0.9, outline = 0.6)

#' Row order for drawing UMAP points: a fixed interleave (golden-ratio stride), the same on every render.
islet_umap_order <- function(n) order((seq_len(n) * 0.6180339887498949) %% 1)

#' `colour` names a discrete column (coloured by `pal`, a named vector) or a numeric one (inferno, scaled to its
#' range, with a colour bar titled `colorbar_title`; NA grey).
islet_umap_plotly <- function(df, x, y, colour, pal = NULL, text, key = NULL, source = NULL,
                              colorbar_title = NULL, size = ISLET_UMAP_MARKER$size) {
  numeric_colour <- is.numeric(df[[colour]])
  levels <- if (numeric_colour) character(0) else if (is.factor(df[[colour]]))
    levels(droplevels(df[[colour]])) else sort(unique(as.character(df[[colour]])))
  df <- df[islet_umap_order(nrow(df)), , drop = FALSE]
  paper <- ISLET_THEME_PALETTE[[islet_theme_mode()]]$paper
  marker <- list(size = size, opacity = ISLET_UMAP_MARKER$opacity,
                 line = list(color = paper, width = ISLET_UMAP_MARKER$outline))
  if (numeric_colour) {
    v <- df[[colour]]
    lim <- if (any(is.finite(v))) range(v, na.rm = TRUE) else c(0, 1)
    stops <- viridisLite::inferno(11)
    scale <- lapply(seq_along(stops), function(i) list((i - 1) / (length(stops) - 1), stops[i]))
    marker <- c(marker, list(color = ifelse(is.finite(v), v, NA), colorscale = scale, cmin = lim[1], cmax = lim[2],
                             showscale = TRUE, colorbar = list(title = list(text = colorbar_title %||% colour))))
  } else {
    col <- unname(pal[as.character(df[[colour]])])
    col[is.na(col)] <- "#9d9d9d"
    marker$color <- col
  }
  p <- plotly::plot_ly(source = source %||% "A")
  p <- plotly::add_trace(p, x = df[[x]], y = df[[y]], type = "scatter", mode = "markers",
                         text = df[[text]], hoverinfo = "text", showlegend = FALSE,
                         key = if (!is.null(key)) df[[key]] else NULL, marker = marker)
  # The legend stubs are added after the build: plot_ly() drops an all-NA trace's name along with its data.
  p <- plotly::plotly_build(p)
  for (lv in levels) {
    lc <- unname(pal[lv])
    if (is.na(lc)) lc <- "#9d9d9d"
    p$x$data[[length(p$x$data) + 1L]] <- list(type = "scatter", mode = "markers", x = list(NA), y = list(NA),
                                              name = lv, showlegend = TRUE, hoverinfo = "skip",
                                              marker = list(color = lc, size = 10))
  }
  plotly::layout(p, legend = list(itemclick = FALSE, itemdoubleclick = FALSE))
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
# UI reference: data/phenotype_rules.csv (parsed once at startup). This is not
# verified as the complete original method that produced all shipped labels.
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
    islet_log_event("phenotype_rules_unavailable")
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
    else "<span style='color:var(--islet-fg-subtle);'>(no positive requirement)</span>"
    neg_html <- if (length(neg_parts) > 0)
      paste0("<span style='color:var(--islet-fg-muted);'>", paste(htmltools::htmlEscape(neg_parts), collapse = ", "), "</span>")
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
  "Beta cell", "Alpha cell", "Delta cell", "CHGA+ hormone-negative", "Endocrine",
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
  protein_path <- data_path("protein_phenotype_rules.csv")
  if (file.exists(protein_path)) {
    rules <- utils::read.csv(protein_path, stringsAsFactors = FALSE, na.strings = "")
    rules[is.na(rules)] <- ""
    headings <- c("Parent", "Label", "All positive", "Any positive", "All negative", "Note")
    rows <- lapply(seq_len(nrow(rules)), function(i) {
      shiny::tags$tr(lapply(unname(as.list(rules[i, ])), shiny::tags$td))
    })
    return(shiny::modalDialog(
      title = "Cell Phenotype Rules", size = "l", easyClose = TRUE,
      footer = shiny::modalButton("Close"),
      if (file.exists(data_path("uniform_release.json"))) shiny::tagList(
        shiny::tags$p("All 20 donors are typed by one workflow: the Islet-Explorer-Senior cell typing, ported to ",
          "this CODEX panel. Marker values are spillover-corrected (REDSEA), aligned between images (UniFORM) and ",
          "put on one cohort scale from background to the 99.9th percentile."),
        shiny::tags$p("Endocrine cells are called only inside an endocrine object (Islet, Cluster or Single; a ",
          "Single is then counted, not analysed, and its cell treated as a tissue cell). Insulin, glucagon and ",
          "somatostatin are mutually exclusive: the strongest names the cell only when it is at least twice the ",
          "next, otherwise the cell is excluded. A CHGA+ cell with no hormone above the image's level is CHGA+ ",
          "hormone-negative. Nerve, ",
          "capillary and lymphatic structures are areas (the Plot tab's Microenvironment layer), not cell types. ",
          "Immune rules count a marker at 0.10 or above on the cohort scale. Every other cell is Exocrine, the ",
          "residual class (acinar, ductal, stromal, endothelial and neural cells together); residual cells inside an ",
          "endocrine object and cells without a nucleus are left out of the analysis."),
        shiny::tags$p("Endocrine, T cell and Immune hold the cells left without a subtype. They count in every total ",
          "but are not offered as phenotypes beside their subtypes.")
      ) else shiny::tagList(
        shiny::tags$p("The original 15 donors retain their saved phenotype labels. ",
          "The rules below were accepted for the added donors 6090, 6232, 6520, 6539 and 6558."),
        shiny::tags$p("A marker is positive at a gate-rescaled signal of 0.5 or higher. ",
          "Rules are applied from parent to child. Conflicting lineages and unsupported cells remain Unknown; ",
          "conflicting subtypes retain their parent label. These rules do not assign Acinar, ECAD+ or Structural.")
      ),
      shiny::tags$div(style = "max-height: 60vh; overflow: auto;",
        shiny::tags$table(class = "table table-striped",
          shiny::tags$thead(shiny::tags$tr(lapply(headings, shiny::tags$th))),
          shiny::tags$tbody(rows)))
    ))
  }
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
      shiny::tags$td(df$Category[i],                       style = "padding:6px 10px; color:var(--islet-fg-muted);"),
      shiny::tags$td(df$Positive[i],                       style = "padding:6px 10px;"),
      shiny::tags$td(df$Negative[i],                       style = "padding:6px 10px; color:var(--islet-fg-muted);"),
      shiny::tags$td(df$Notes[i],                          style = "padding:6px 10px; color:var(--islet-fg); font-size: 13px;")
    )
  })
  shiny::modalDialog(
    title = shiny::tagList(shiny::icon("dna"), "Cell Phenotype Rule Reference"),
    size = "l", easyClose = TRUE, fade = TRUE,
    footer = shiny::modalButton("Close"),
    shiny::tags$p(
      "This table describes a documented upstream Rules1 marker-gating stage. ",
      "It does not fully reproduce the phenotype labels saved in the app data. A ",
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
    shiny::tags$p(style = "margin-top: 10px; font-size: 12px; color: var(--islet-fg-subtle);",
      "Reference CSV: ",
      shiny::tags$code("data/phenotype_rules.csv"),
      " (15 markers: INS, GCG, SST, SMA, CD20, CD3e, CD68, CD163, CD8a, CD4, HLADR, PDPN, CD34, CD31, B3TUBB)."
    ),
    shiny::tags$p(style = "margin-top: 6px; font-size: 12px; color: var(--islet-fg-subtle);",
      shiny::tags$strong("Data vs. rules: "),
      "the full original rule matrix has not been verified. Some rule phenotypes ",
      "have no cells in the shipped dataset because upstream labeling merged or renamed them ",
      "(e.g. Dendritic cell folds into the ", shiny::tags$code("APCs"),
      "/Immune aggregates; Exocrine splits into Acinar/Ductal). ",
      shiny::tags$code("APCs"),
      " is an HLADR+ myeloid aggregate created upstream, not itself a ",
      shiny::tags$code("phenotype_rules.csv"), " row."
    ),
    shiny::tags$p(style = "margin-top: 6px; font-size: 12px; color: var(--islet-fg-subtle);",
      shiny::tags$strong("Panel note: "),
      "the original Rules1 template also listed ", shiny::tags$code("CD99"),
      " (Endocrine) and ", shiny::tags$code("CD56"),
      " (Neural), but neither antibody is in this 31-plex panel, so the ",
      "available Rules1 implementation omitted them. In that reference, the gates become ",
      shiny::tags$code("Endocrine = any of INS/GCG/SST+"), " and ",
      shiny::tags$code("Neural = B3TUBB+"), ". These are not a complete reconstruction of shipped labels."
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
    style = "margin-bottom: 12px; padding: 8px 14px; background-color: var(--islet-info-bg); border: 1px solid var(--islet-border-strong); border-radius: 6px; display: inline-block;",
    shiny::actionLink(input_id,
      shiny::tagList(shiny::icon("dna", style = "color:var(--islet-link);"),
                     shiny::tags$span(style = "margin-left:6px; color:var(--islet-link); font-weight:600;",
                                      label)),
      style = "text-decoration: none;")
  )
}

# ---- Control panels (owner, 2026-10-17: one layout and separation in every tab) ----
# Every tab puts its tab-wide controls in two collapsible panels at the top of the tab: "Data controls" (what is
# selected, filtered or computed; open) and "View controls" (how it is drawn; closed). A selector that only
# chooses what one card plots stays on that card. Downloads sit under the plot they export.
CONTROL_PANEL_CSS <- paste0(
  ".control-panel { border: 1px solid var(--islet-border-strong); border-radius: 6px; padding: 0 12px; margin-bottom: 12px; ",
  "background: var(--islet-panel-bg); } ",
  ".control-panel > summary { cursor: pointer; list-style: none; font-weight: 600; color: var(--islet-panel-fg); ",
  "padding: 10px 0; display: flex; align-items: center; gap: 8px; } ",
  ".control-panel > summary::-webkit-details-marker { display: none; } ",
  ".control-panel > summary::marker { content: ''; } ",
  ".control-panel > summary:before { content: '▸'; color: var(--islet-link); font-size: 14px; } ",
  ".control-panel[open] > summary:before { content: '▾'; } ",
  ".control-panel > summary:after { content: 'Click to show'; margin-left: auto; font-size: 12px; ",
  "font-weight: 400; color: var(--islet-panel-hint); } ",
  ".control-panel[open] > summary:after { content: 'Click to hide'; } ",
  ".control-panel[open] > summary { border-bottom: 1px solid var(--islet-border-strong); margin-bottom: 10px; } ",
  ".control-panel .form-group { margin-bottom: 6px; } ",
  ".control-panel h5 { margin: 4px 0 6px; font-size: 14px; font-weight: 600; color: var(--islet-fg); } ",
  ".control-panel .form-group:last-child { margin-bottom: 8px; }",
  ".control-panel-wrap { position: relative; } ",
  ".control-panel-wrap > .control-panel > summary { padding-right: 34px; } ",
  ".control-panel-wrap > .card-help { position: absolute; top: 7px; right: 4px; } ",
  ".card-help { display: inline-block; position: relative; z-index: 20; } ",
  ".card-help > summary { cursor: pointer; list-style: none; color: var(--islet-link); padding: 0 4px; ",
  "line-height: 1; } ",
  ".card-help > summary::-webkit-details-marker { display: none; } ",
  ".card-help > summary::marker { content: ''; } ",
  ".card-help > summary i { font-size: 18px; } ",
  ".card-help-panel { position: absolute; right: 0; top: 100%; width: min(420px, calc(100vw - 32px)); ",
  "max-height: 70vh; overflow-y: auto; padding: 12px; background: var(--islet-surface); text-align: left; ",
  "border: 1px solid var(--islet-border-strong); box-shadow: 0 4px 16px #0002; font-size: 13px; ",
  "font-weight: 400; color: var(--islet-fg); z-index: 1050; } ",
  ".card-help-panel p:last-child { margin-bottom: 0; }")

#' A "?" icon whose popover holds the help for the card or control panel it sits on (owner, 2026-10-20:
#' help is split by card, not gathered under one icon per tab). www/card-help.js keeps the popover inside the
#' window and closes it on an outside click or Escape. Put it beside the card's title.
card_help <- function(title, ...) {
  shiny::tags$details(class = "card-help",
    shiny::tags$summary(title = title, shiny::tags$span(class = "sr-only", title),
                        shiny::tags$i(class = "far fa-circle-question", `aria-hidden` = "true")),
    shiny::tags$div(class = "card-help-panel", ...))
}

#' A card title (h5) with its help icon beside it. `help` is a list of tags; NULL gives the bare title.
card_title_help <- function(title, help = NULL, style = "margin: 0;") {
  heading <- shiny::h5(title, style = style)
  if (is.null(help)) return(heading)
  shiny::tags$div(style = "display: flex; align-items: center; gap: 4px;",
                  heading, do.call(card_help, c(list(paste(title, "help")), help)))
}

#' A collapsible control group. Use control_panel("Data controls", open = TRUE, ...) and
#' control_panel("View controls", ...) at the top of a tab or section. `help`, a list of tags, puts a "?" icon at
#' the right of the panel's title bar, outside the <summary> so opening it does not open or close the panel.
control_panel <- function(title, ..., open = FALSE, help = NULL) {
  panel <- shiny::tags$details(class = "control-panel", open = if (isTRUE(open)) NA else NULL,
                               shiny::tags$summary(title), ...)
  if (is.null(help)) return(panel)
  shiny::tags$div(class = "control-panel-wrap", panel,
                  do.call(card_help, c(list(paste(title, "help")), help)))
}

#' TRUE when the served root is the uniform 20-donor build (one typing workflow for every donor).
is_uniform_root <- function() file.exists(data_path("uniform_release.json"))

#' The cell labels a uniform root assigns (uniform_release.json "labels"); NULL for other roots.
root_phenotype_labels <- function() {
  path <- data_path("uniform_release.json")
  if (!file.exists(path)) return(NULL)
  labels <- tryCatch(jsonlite::fromJSON(path)$labels, error = function(e) NULL)
  if (is.null(labels)) NULL else as.character(labels)
}

#' The "Endocrine objects" class selector. It sits in the Data controls of the Plot tab (`object_classes`) and of
#' the Spatial tab (`spatial_object_classes`), which app.R keeps in step; both filter prepared() (owner, 2026-10-19:
#' "The selector for islets and clusters needs to go in the data controls").
object_class_selector <- function(input_id) {
  shiny::checkboxGroupInput(input_id, "Endocrine objects", choices = OBJECT_CLASSES,
                            selected = unname(OBJECT_CLASSES), inline = TRUE)
}

#' Help paragraph for the app-wide "Endocrine objects" selector, shared by every tab's help popover.
object_classes_help <- function() {
  if (!is_uniform_root()) return(NULL)
  shiny::tags$p(shiny::tags$strong("Endocrine objects."),
    " The selector in the Data controls of the Plot and Spatial tabs (the two are kept in step) chooses the ",
    "classes analysed. An Islet has two or more whole cells and is at least ",
    "35 µm across (the mean of its largest and smallest caliper diameter); a Cluster has two or more cells ",
    "and is smaller. Rows and plots use only the selected classes. The Trajectory tab instead fits every object ",
    "with ten or more core cells, whatever its class. A Single (one endocrine cell outside any islet or cluster) ",
    "is not an analysis object: it has no core or peri-islet band, and its cell is an ordinary tissue cell, which ",
    "may fall in a neighbouring object's band. Singles are only counted (the Plot tab's Endocrine objects per ",
    "donor card) and drawn on the segmentation map.")
}

#' Help paragraph on how cells are typed, shared by every tab's help popover.
cell_typing_help <- function() {
  if (!is_uniform_root()) {
    return(shiny::tags$p(shiny::tags$strong("Cell labels."),
      " The original 15 donors retain their saved cell labels; the five additions use reviewed reconstructed ",
      "phenotype rules. Phenotype comparisons span both methods."))
  }
  shiny::tags$p(shiny::tags$strong("Cell labels."),
    " Every cell of all 20 donors is typed by one workflow. The phenotypes listed are Beta, Alpha and Delta cell, ",
    "CHGA+ hormone-negative (an endocrine cell with no hormone), CD8a and CD4 T cell, B cell, Macrophage and ",
    "Exocrine (the residual class). Nerve, capillary and lymphatic structures are areas, not cell types. The ",
    "hormones are mutually exclusive: a cell whose strongest hormone is under twice the next one is excluded. Cells ",
    "left without a subtype (T cells without CD4 or CD8, other leukocytes) count in every total but are not listed. ",
    "Residual cells inside an endocrine object and cells without a nucleus are left out. Endocrine cells are called ",
    "inside Islets and Clusters, and Singles' cells keep their hormone type, so a peri-islet band holds endocrine ",
    "cells only where a Single lies in it.")
}

#' One help icon for a whole tab. Since 2026-10-20 only the hidden Statistics and Inflammation tabs use it;
#' the reviewer tabs put each card's help on that card (card_help(), control_panel(help = )).
tab_help <- function(title, rules_input_id, ...) {
  shiny::tags$div(
    style = "display: flex; justify-content: flex-end; margin-bottom: 6px;",
    shiny::tags$details(style = "position: relative; z-index: 20;",
      shiny::tags$summary(title = title,
        style = "cursor: pointer; list-style: none; color: var(--islet-link); padding: 4px 8px;",
        shiny::tags$span(class = "sr-only", title),
        shiny::tags$i(class = "far fa-circle-question", `aria-hidden` = "true",
                      style = "font-size: 24px;")),
      shiny::tags$div(
        style = paste0("position: absolute; right: 0; top: 100%; width: min(420px, calc(100vw - 32px));",
                       " max-height: 70vh; overflow-y: auto; padding: 12px; background: var(--islet-surface);",
                       " border: 1px solid var(--islet-border-strong); box-shadow: 0 4px 16px #0002; font-size: 13px;"),
        phenotype_rules_button(rules_input_id), ...)
    )
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
      shiny::tags$span(style = "color:var(--islet-fg-subtle);", " · "),
      shiny::actionLink(show_modal_link, "see all rules →",
                        style = "font-size: 12px; color: var(--islet-link);")
    )
  }
  shiny::tags$div(
    style = "margin-top: -6px; margin-bottom: 8px; padding: 6px 10px; background-color: var(--islet-surface-2); border-left: 3px solid #4477AA; border-radius: 3px; font-size: 13px; color: var(--islet-fg); line-height: 1.4;",
    shiny::tags$strong(style = "color:var(--islet-link);", paste0(prefix, " ")),
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
# "Endocrine", "T cell" and "Immune" are general classes whose subtypes are listed (Beta / Alpha / Delta cell;
# CD8a / CD4 Tcell; T cell and Macrophage), so they are not offered beside them (owner, 2026-10-17: "There should
# not be a phenotype listed that is a general class when there are subtypes"). In the uniform root they hold the
# cells left without a subtype; those cells still count in every denominator. The rest are labels of the frozen
# typing that the uniform root does not produce. Keep python/paga_pseudotime.py RETIRED_PHENOTYPES in step.
PHENOTYPE_RETIRED <- c("Endocrine", "Acinar", "Endothelial", "T cell",
                       "Immune", "Unknown", "ECAD+", "Structural", "APCs")

# Endocrine object classes analysed by the app (owner A5, 2026-10-16, with a filter). Singles (one endocrine cell)
# are not analysis objects (owner, 2026-10-18: "Singles should not be included in any islet-level analysis except
# for quantifying or showing how many there are"); a uniform root counts them in endocrine_objects_by_donor.csv
# and draws them on the segmentation map.
OBJECT_CLASSES <- c("Islet" = "Islet", "Cluster" = "EndocrineCluster")

#' Cache-key suffix for an object-class selection (NULL = all classes).
object_class_key <- function(object_classes) {
  if (is.null(object_classes)) return("")
  paste0("_", paste(sort(as.character(object_classes)), collapse = "+"))
}

#' The selection as the Python bridge takes it (NULL = all classes).
object_class_arg <- function(object_classes) {
  if (is.null(object_classes)) NULL else as.list(as.character(object_classes))
}

#' Drop retired phenotypes from a character vector of labels.
phenotype_keep <- function(x) setdiff(x, PHENOTYPE_RETIRED)

#' TRUE for rows whose phenotype is still in use. NA-safe.
phenotype_is_kept <- function(x) !(as.character(x) %in% PHENOTYPE_RETIRED)

#' Marker definition text for each retained phenotype, for help popovers.
#'
#' When the data root ships `protein_phenotype_rules.csv` (the accepted rules
#' for the five added donors), each label's rule is composed with its parents:
#' required positives and negatives accumulate down the tree, and a parent's
#' "any of" group is omitted when the child's own positive evidence already
#' satisfies it. Otherwise the Rules1 reference (`format_phenotype_rule()`) is
#' used. Returns a named character vector in palette order.
phenotype_marker_definitions <- function(path = data_path("protein_phenotype_rules.csv")) {
  kept <- phenotype_keep(names(PHENOTYPE_COLORS))
  if (!file.exists(path)) {
    out <- vapply(kept, function(n) format_phenotype_rule(n, "plain"), character(1))
    return(out[!startsWith(out, "No marker rule")])
  }
  rules <- utils::read.csv(path, stringsAsFactors = FALSE, na.strings = "", check.names = FALSE)
  rules[is.na(rules)] <- ""
  split_markers <- function(x) if (nzchar(x)) strsplit(x, "|", fixed = TRUE)[[1]] else character(0)
  by_label <- split(rules, rules$label)
  compose <- function(label) {
    row <- by_label[[label]]
    pos <- split_markers(row$all_positive); any_child <- split_markers(row$any_positive)
    neg <- split_markers(row$all_negative); groups <- list()
    if (length(any_child)) groups <- list(any_child)
    parent <- row$parent
    while (nzchar(parent) && parent != "all" && !is.null(by_label[[parent]])) {
      prow <- by_label[[parent]]
      pos <- union(pos, split_markers(prow$all_positive))
      pany <- split_markers(prow$any_positive)
      if (length(pany) && !length(intersect(pany, c(pos, any_child)))) groups <- c(groups, list(pany))
      neg <- union(neg, split_markers(prow$all_negative))
      parent <- prow$parent
    }
    parts <- c(if (length(pos)) paste0(pos, "+"),
               vapply(groups, function(g) if (length(g) == 1) paste0(g, "+")
                      else paste0("any of ", paste0(g, "+", collapse = " / ")), character(1)))
    txt <- paste(parts, collapse = ", ")
    if (length(neg)) txt <- paste0(txt, "  \u2014  ", paste0(neg, "\u2212", collapse = ", "))
    txt
  }
  kept <- kept[kept %in% names(by_label)]
  stats::setNames(vapply(kept, compose, character(1)), kept)
}

#' Help-popover block listing the retained phenotypes' marker definitions.
# The listed phenotypes of the uniform root, in words (the full rule table is in the rules modal).
UNIFORM_PHENOTYPE_DEFINITIONS <- c(
  "Beta cell" = paste("inside an Islet or Cluster, a cell admitted as endocrine whose highest hormone is insulin, at",
                      "least twice the next hormone (a cell under that margin is excluded: the hormones are mutually",
                      "exclusive, and the runner-up is mostly a touching neighbour's hormone)"),
  "Alpha cell" = "as Beta cell, with glucagon the highest hormone",
  "Delta cell" = "as Beta cell, with somatostatin the highest hormone",
  "CHGA+ hormone-negative" = paste("inside an Islet or Cluster, CHGA at or above the image's tissue CHGA level",
                                   "and insulin, glucagon and somatostatin all below the image's hormone level"),
  "CD8a Tcell" = "CD3e+ and CD8a+, CD4\u2212, no epithelial marker (E-cadherin, panCK)",
  "CD4 Tcell" = "CD3e+ and CD4+, CD8a\u2212, no epithelial marker",
  "B cell" = "CD20+ and CD45+, CD3e\u2212, CD68\u2212, no epithelial marker",
  "Macrophage" = "CD68+, CD3e\u2212, no epithelial marker",
  "Exocrine" = paste("the residual class: every other cell outside the endocrine objects (acinar, ductal, stromal,",
                     "endothelial and neural cells together)"))

phenotype_marker_definitions_ui <- function() {
  if (is_uniform_root()) {
    defs <- UNIFORM_PHENOTYPE_DEFINITIONS
    return(shiny::tagList(
      shiny::tags$h6(style = "margin-top: 12px; font-weight: 700;", "Phenotype definitions"),
      shiny::tags$p(style = "color: var(--islet-fg-muted);",
        "For the immune and epithelial markers, + means 0.10 or above on the cohort scale (0 = background, ",
        "1 = the cohort's 99.9th percentile). Hormones and structure markers use each image's own levels."),
      shiny::tags$ul(style = "padding-left: 18px; margin-bottom: 6px;",
        lapply(names(defs), function(n) shiny::tags$li(shiny::tags$strong(n), ": ", defs[[n]])))))
  }
  defs <- phenotype_marker_definitions()
  if (!length(defs)) return(NULL)
  protein <- file.exists(data_path("protein_phenotype_rules.csv"))
  shiny::tagList(
    shiny::tags$h6(style = "margin-top: 12px; font-weight: 700;", "Phenotype marker definitions"),
    shiny::tags$p(style = "color: var(--islet-fg-muted);",
      if (protein) paste0("Rules accepted for the five added donors (a marker is positive at a gate-rescaled ",
                          "signal of 0.5 or higher). The original 15 donors keep their saved labels.")
      else "Rules1 marker reference (a partial reference for the saved labels)."),
    shiny::tags$ul(style = "padding-left: 18px; margin-bottom: 6px;",
      lapply(names(defs), function(n)
        shiny::tags$li(shiny::tags$strong(n), ": ", defs[[n]])))
  )
}

#' Drop retired phenotypes from a data frame with a phenotype column.
phenotype_filter_rows <- function(df, col = "phenotype") {
  if (is.null(df) || !is.data.frame(df) || !(col %in% names(df))) return(df)
  df[phenotype_is_kept(df[[col]]), , drop = FALSE]
}
