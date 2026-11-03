# ===== Single-Cell Drill-Down Helper Functions =====
# Provides cell-level visualization for the segmentation viewer.
# Loads per-islet CSV files from data/cells/ and renders cell scatter plots
# over the GeoJSON segmentation base.

# ── Phenotype color palette (21 cell types → distinctive hex) ──────────
# Endocrine = warm, Immune = cool, Structural = neutral
PHENOTYPE_COLORS <- c(
  "Beta cell"     = "#E63946",  # red
  "Alpha cell"    = "#F4A261",  # orange
  "Delta cell"    = "#E9C46A",  # gold
  "CHGA+ hormone-negative" = "#6D597A",  # muted plum
  "Endocrine"     = "#D4A373",  # tan
  "Acinar"        = "#A8DADC",  # light teal
  "Exocrine"     = "#A8DADC",  # Senior's residual class: the Acinar color
  "Ductal"        = "#457B9D",  # steel blue
  "ECAD+"         = "#6A994E",  # olive green
  "Structural"    = "#9B9B9B",  # gray
  "SMA+"          = "#C49792",  # mauve
  "Blood Vessel"  = "#BC4749",  # brick red
  "Endothelial"   = "#F28482",  # salmon
  "Lymphatic"     = "#84A98C",  # sage
  "Neural"        = "#B5838D",  # dusty rose
  "Macrophage"    = "#264653",  # dark teal
  "APCs"          = "#2A9D8F",  # teal
  "CD8a Tcell"    = "#023E8A",  # navy
  "CD4 Tcell"     = "#0077B6",  # blue
  "T cell"        = "#0096C7",  # light blue
  "B cell"        = "#48CAE4",  # sky blue
  "Immune"        = "#3A5A8C",  # medium blue
  "Unknown"       = "#CCCCCC"   # light gray
)

# ── High Contrast palette: immune spread across purple/cyan/magenta/indigo ──
PHENOTYPE_COLORS_HIGH_CONTRAST <- c(
  "Beta cell"     = "#D32F2F",  # deep red
  "Alpha cell"    = "#FF8F00",  # amber
  "Delta cell"    = "#FDD835",  # bright yellow
  "CHGA+ hormone-negative" = "#8D6E63",  # brown (was deep purple, a shade from CD8a Tcell and Macrophage; 2026-10-21)
  "Endocrine"     = "#8D6E63",  # brown
  "Acinar"        = "#00ACC1",  # cyan
  "Exocrine"     = "#00ACC1",  # Senior's residual class: the Acinar color
  "Ductal"        = "#43A047",  # green
  "ECAD+"         = "#7CB342",  # lime green
  "Structural"    = "#9E9E9E",  # gray
  "SMA+"          = "#F06292",  # pink
  "Blood Vessel"  = "#E91E63",  # magenta-red
  "Endothelial"   = "#FF7043",  # deep orange
  "Lymphatic"     = "#26A69A",  # teal
  "Neural"        = "#AB47BC",  # purple
  "Macrophage"    = "#5C6BC0",  # indigo
  "APCs"          = "#26C6DA",  # light cyan
  "CD8a Tcell"    = "#7B1FA2",  # deep purple
  "CD4 Tcell"     = "#1565C0",  # strong blue
  "T cell"        = "#0097A7",  # dark cyan
  "B cell"        = "#EC407A",  # hot pink
  "Immune"        = "#283593",  # dark indigo
  "Unknown"       = "#BDBDBD"   # light gray
)

# ── Colorblind Safe palette: adapted from Paul Tol qualitative schemes ──
PHENOTYPE_COLORS_COLORBLIND <- c(
  "Beta cell"     = "#CC6677",  # rose
  "Alpha cell"    = "#DDCC77",  # sand
  "Delta cell"    = "#999933",  # olive
  "CHGA+ hormone-negative" = "#661100",  # dark red
  "Endocrine"     = "#AA4499",  # purple
  "Acinar"        = "#44AA99",  # teal
  "Exocrine"     = "#44AA99",  # Senior's residual class: the Acinar color
  "Ductal"        = "#332288",  # indigo
  "ECAD+"         = "#117733",  # green
  "Structural"    = "#888888",  # gray
  "SMA+"          = "#882255",  # wine
  "Blood Vessel"  = "#EE8866",  # peach
  "Endothelial"   = "#EEDD88",  # light yellow
  "Lymphatic"     = "#77AADD",  # light blue
  "Neural"        = "#FFAABB",  # pink
  "Macrophage"    = "#44BB99",  # mint
  "APCs"          = "#BBCC33",  # pear
  "CD8a Tcell"    = "#AA4455",  # dark rose
  "CD4 Tcell"     = "#4477AA",  # steel blue
  "T cell"        = "#66CCEE",  # cyan
  "B cell"        = "#EE6677",  # coral
  "Immune"        = "#225566",  # dark teal
  "Unknown"       = "#BBBBBB"   # light gray
)

# ── Maximum Distinction palette: adapted from Kelly's 22 maximally distinct ──
PHENOTYPE_COLORS_MAX_DISTINCT <- c(
  "Beta cell"     = "#BE0032",  # vivid red
  "Alpha cell"    = "#F3C300",  # vivid yellow
  "Delta cell"    = "#875692",  # strong purple
  "CHGA+ hormone-negative" = "#882D17",  # reddish brown
  "Endocrine"     = "#F38400",  # vivid orange
  "Acinar"        = "#A1CAF1",  # very light blue
  "Exocrine"     = "#A1CAF1",  # Senior's residual class: the Acinar color
  "Ductal"        = "#008856",  # vivid green
  "ECAD+"         = "#C2B280",  # dark yellowish brown
  "Structural"    = "#848482",  # medium gray
  "SMA+"          = "#E68FAC",  # purplish pink
  "Blood Vessel"  = "#F99379",  # vivid orange-yellow
  "Endothelial"   = "#604E97",  # strong violet
  "Lymphatic"     = "#0067A5",  # strong blue
  "Neural"        = "#DCD300",  # vivid greenish yellow
  "Macrophage"    = "#7F180D",  # vivid reddish brown
  "APCs"          = "#B3446C",  # deep purplish pink
  "CD8a Tcell"    = "#F6A600",  # vivid yellowish brown
  "CD4 Tcell"     = "#2B3D26",  # dark olive green
  "T cell"        = "#8DB600",  # vivid yellowish green
  "B cell"        = "#654522",  # deep yellowish brown
  "Immune"        = "#E25822",  # vivid reddish orange
  "Unknown"       = "#C0C0C0"   # silver gray
)

# ── Macrophage marker-profile labels (owner, 2026-10-27; MACROPHAGE_LABELS in 00_globals.R) ──
# Twelve shades from each palette's Macrophage color: the CD68/CD163 group turns the hue (CD68+CD163+ the base,
# CD68+CD163- one way, CD68-CD163+ the other) and LGALS3/CD4 set the lightness (LGALS3+CD4+ darkest, LGALS3-CD4-
# lightest), so the family reads as macrophages and each label differs from its neighbors in lightness, which survives
# color-vision deficiency. Each palette's single "Macrophage" (the frozen root) stays as it is.
macrophage_family <- function(base) {
  hsv0 <- grDevices::rgb2hsv(grDevices::col2rgb(base))[, 1]
  mix <- function(col, t) {          # t < 0 toward black, t > 0 toward white
    rgb <- grDevices::col2rgb(col)[, 1] / 255
    rgb <- if (t < 0) rgb * (1 + t) else rgb + (1 - rgb) * t
    grDevices::rgb(rgb[1], rgb[2], rgb[3])
  }
  heads <- c(0, 0.09, -0.09)                     # hue turn, CD68+CD163+, CD68+CD163-, CD68-CD163+
  tails <- c(-0.35, -0.08, 0.25, 0.5)            # LGALS3+CD4+, LGALS3+CD4-, LGALS3-CD4+, LGALS3-CD4-
  cols <- unlist(lapply(heads, function(dh) {
    hue <- grDevices::hsv((hsv0[["h"]] + dh) %% 1, max(hsv0[["s"]], 0.35), max(min(hsv0[["v"]], 0.85), 0.45))
    vapply(tails, function(t) mix(hue, t), character(1))
  }))
  stats::setNames(cols, MACROPHAGE_LABELS)
}
# The merged categories of every marker choice (MACROPHAGE_MARKERS, owner 2026-10-27): with all four the twelve shades
# above; with fewer, the choice's categories ramp from dark to light around the Macrophage color, in key order;
# with none the one "Macrophage" keeps its color. Only one choice is shown at a time.
macrophage_key_colors <- function(base, family) {
  out <- character(0)
  for (k in seq_along(MACROPHAGE_MARKERS)) {
    for (markers in utils::combn(MACROPHAGE_MARKERS, k, simplify = FALSE)) {
      keys <- unique(macrophage_collapse(MACROPHAGE_LABELS, markers))
      keys <- keys[order(chartr("+-", "01", keys))]
      cols <- if (k == length(MACROPHAGE_MARKERS)) unname(family[match(keys, macrophage_collapse(names(family), markers))])
              else grDevices::colorRampPalette(c(family[[1]], base, family[[length(family)]]))(max(length(keys), 2))[
                     seq_along(keys)]
      out <- c(out, stats::setNames(cols, keys))
    }
  }
  out
}
with_macrophage_family <- function(pal) {
  i <- match("Macrophage", names(pal))
  family <- macrophage_family(pal[["Macrophage"]])
  c(pal[seq_len(i)], family, macrophage_key_colors(pal[["Macrophage"]], family), pal[-seq_len(i)])
}
# Top-marker labels (root v15; TOP_MARKER_LABELS in 00_globals.R): a muted family -- low chroma, hues around the
# wheel, two alternating lightnesses -- so they read apart from each other but not as the saturated endocrine,
# immune or macrophage colors. On a root that assigns them, "SMA+" takes its family color (the frozen root's "SMA+",
# smooth muscle, keeps its own).
top_marker_family <- function() {
  n <- length(TOP_MARKER_LABELS)
  stats::setNames(grDevices::hcl(h = seq(15, 345, length.out = n), c = 26, l = rep(c(50, 68), length.out = n)),
                  TOP_MARKER_LABELS)
}
with_top_marker_family <- function(pal) {
  family <- top_marker_family()
  assigned <- isTRUE(tryCatch(root_has_top_marker_labels(), error = function(e) FALSE))
  rule <- isTRUE(tryCatch(root_cd45_is_rule_label(), error = function(e) FALSE))
  # "CD45+": from v16 the CD45 rule's class, in the palette's Immune color (that root has no "Immune"); on v15 the
  # leftover CD45 cell, in the family color of "CD45+ (top marker)".
  family[[CD45_RULE_LABEL]] <- if (rule) pal[["Immune"]] else family[[CD45_TOP_MARKER_LABEL]]
  add <- if (assigned) family else family[setdiff(names(family), names(pal))]
  pal[names(add)] <- add
  pal
}
PHENOTYPE_COLORS <- with_macrophage_family(PHENOTYPE_COLORS)
PHENOTYPE_COLORS_HIGH_CONTRAST <- with_macrophage_family(PHENOTYPE_COLORS_HIGH_CONTRAST)
PHENOTYPE_COLORS_COLORBLIND <- with_macrophage_family(PHENOTYPE_COLORS_COLORBLIND)
PHENOTYPE_COLORS_MAX_DISTINCT <- with_macrophage_family(PHENOTYPE_COLORS_MAX_DISTINCT)
PHENOTYPE_COLORS <- with_top_marker_family(PHENOTYPE_COLORS)
PHENOTYPE_COLORS_HIGH_CONTRAST <- with_top_marker_family(PHENOTYPE_COLORS_HIGH_CONTRAST)
PHENOTYPE_COLORS_COLORBLIND <- with_top_marker_family(PHENOTYPE_COLORS_COLORBLIND)
PHENOTYPE_COLORS_MAX_DISTINCT <- with_top_marker_family(PHENOTYPE_COLORS_MAX_DISTINCT)

# ── Palette registry and accessor ────────────────────────────────────────
PHENOTYPE_PALETTES <- list(
  "Original"             = PHENOTYPE_COLORS,
  "High Contrast"        = PHENOTYPE_COLORS_HIGH_CONTRAST,
  "Colorblind Safe"      = PHENOTYPE_COLORS_COLORBLIND,
  "Maximum Distinction"  = PHENOTYPE_COLORS_MAX_DISTINCT
)

get_phenotype_palette <- function(name = "Original") {
  PHENOTYPE_PALETTES[[name]] %||% PHENOTYPE_COLORS
}

# ── Cell CSV cache (env-based, same pattern as geojson_cache) ──────────
drilldown_cache <- new.env()

# Path to per-islet cell CSVs
CELLS_DIR <- data_path("cells")

#' Check if single-cell drill-down data is available
drilldown_available <- function() {
  dir.exists(CELLS_DIR) && length(list.files(CELLS_DIR, pattern = "\\.csv$")) > 0
}

#' Load single-cell data for a specific islet
#' @param imageid Case/image ID (e.g., "6505")
#' @param islet_key Islet key (e.g., "Islet_284")
#' @return data.frame with columns: X_centroid, Y_centroid, phenotype, cell_region, + markers
load_islet_cells <- function(imageid, islet_key) {
  imageid <- as.character(imageid)
  islet_key <- as.character(islet_key)
  if (!islet_safe_path_segment(imageid) || !islet_safe_path_segment(islet_key)) {
    islet_log_validation_failure("islet_key_invalid")
    return(NULL)
  }
  combined_id <- paste0(imageid, "_", islet_key)
  cache_key <- combined_id

  if (exists(cache_key, envir = drilldown_cache)) {
    return(get(cache_key, envir = drilldown_cache))
  }

  csv_path <- file.path(CELLS_DIR, paste0(combined_id, ".csv"))
  if (!file.exists(csv_path)) {
    return(NULL)
  }

  cells <- tryCatch({
    df <- read.csv(csv_path, stringsAsFactors = FALSE)
    if (nrow(df) == 0) return(NULL)
    df
  }, error = function(e) {
    islet_log_event("drilldown_cells_load_failed")
    NULL
  })

  if (!is.null(cells)) {
    assign(cache_key, cells, envir = drilldown_cache)
  }

  cells
}

#' Render single-cell drill-down plot over GeoJSON base
#' @param info List with case_id, islet_key, centroid_x, centroid_y
#' @param cells data.frame from load_islet_cells()
#' @param color_by "phenotype" or a marker column name
#' @param show_peri Logical; if FALSE, only show core cells
#' @return ggplot object
render_islet_drilldown_plot <- function(info, cells, color_by = "phenotype", show_peri = TRUE,
                                        show_peri_boundary = TRUE, show_structures = TRUE,
                                        palette = PHENOTYPE_COLORS) {
  if (is.null(info) || is.null(cells) || nrow(cells) == 0) return(NULL)

  # Build base plot: GeoJSON boundaries (islet core always drawn)
  base_plot <- build_segmentation_base_plot(info, show_peri_boundary = show_peri_boundary,
                                             show_structures = show_structures)
  if (is.null(base_plot)) base_plot <- ggplot2::ggplot() + ggplot2::theme_void()

  # Filter cells by region
  if (!show_peri && "cell_region" %in% colnames(cells)) {
    cells <- cells[cells$cell_region == "core", , drop = FALSE]
  }
  # Apply the presentation policy to points, not just the palette. Keep the
  # loaded scientific table intact for original-denominator counts elsewhere.
  if ("phenotype" %in% names(cells)) {
    visible <- if (color_by == "phenotype") phenotype_is_kept(cells$phenotype) else
      !phenotype_is_conflict(cells$phenotype)
    cells <- cells[visible, , drop = FALSE]
  }
  if (nrow(cells) == 0) return(base_plot)

  # Convert cell centroids from µm to pixels (matching GeoJSON coordinate space)
  cells$x_px <- cells$X_centroid / PIXEL_SIZE_UM
  cells$y_px <- cells$Y_centroid / PIXEL_SIZE_UM

  if (color_by == "phenotype" && "phenotype" %in% colnames(cells)) {
    # Categorical coloring by phenotype (using fill to avoid color conflict with structures)
    pheno_present <- phenotype_keep(sort(unique(cells$phenotype)))
    pal <- stats::setNames(unname(palette[pheno_present]), pheno_present)
    # Fill missing phenotypes with gray
    pal[is.na(pal)] <- "#CCCCCC"

    p <- base_plot +
      ggplot2::geom_point(
        data = cells,
        ggplot2::aes(x = x_px, y = y_px, fill = phenotype),
        shape = 21, colour = "grey30", stroke = 0.3, size = 3.0, alpha = 0.9,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_manual(values = pal, name = "Phenotype", labels = phenotype_text) +
      ggplot2::labs(title = paste(info$islet_key, "- Single Cells"))
  } else if (color_by %in% colnames(cells)) {
    # Continuous coloring by marker expression (using fill)
    cells$marker_val <- suppressWarnings(as.numeric(cells[[color_by]]))

    p <- base_plot +
      ggplot2::geom_point(
        data = cells,
        ggplot2::aes(x = x_px, y = y_px, fill = marker_val),
        shape = 21, colour = "grey30", stroke = 0.3, size = 3.0, alpha = 0.9,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_viridis_c(option = "inferno", name = color_by, na.value = "gray80") +
      ggplot2::labs(title = paste(info$islet_key, "-", color_by))
  } else {
    # Fallback: just plot cells in gray
    p <- base_plot +
      ggplot2::geom_point(
        data = cells,
        ggplot2::aes(x = x_px, y = y_px),
        color = "gray40", size = 2.5, alpha = 0.7,
        inherit.aes = FALSE
      ) +
      ggplot2::labs(title = paste(info$islet_key, "- Single Cells"))
  }

  # theme_islet() FIRST, then the local overrides. Without it this plot inherited
  # whatever base_plot carried — and on the `theme_void()` fallback path (no
  # GeoJSON) that meant ggplot's 11pt default for the entire phenotype legend,
  # the densest text in the drill-down.
  p +
    theme_islet() +
    ggplot2::theme(
      legend.position = "right",
      legend.key.size = ggplot2::unit(0.9, "cm"),
      legend.text     = ggplot2::element_text(size = ISLET_FONT_SIZES$legend_text),
      legend.title    = ggplot2::element_text(size = ISLET_FONT_SIZES$legend_title,
                                              face = "bold"),
      # Spatial map: axes are pixel coordinates and carry no meaning.
      axis.title      = ggplot2::element_blank(),
      axis.text       = ggplot2::element_blank(),
      panel.grid      = ggplot2::element_blank()
    )
}

#' Render phenotype composition summary bar chart
#' @param cells data.frame from load_islet_cells()
#' @return ggplot object
render_drilldown_summary <- function(cells, palette = PHENOTYPE_COLORS) {
  if (is.null(cells) || nrow(cells) == 0 || !"phenotype" %in% colnames(cells)) {
    return(NULL)
  }

  # General classes are not listed beside their subtypes (PHENOTYPE_RETIRED).
  cells <- cells[phenotype_is_kept(cells$phenotype), , drop = FALSE]
  if (nrow(cells) == 0) return(NULL)

  # Count phenotypes, split by region
  if ("cell_region" %in% colnames(cells)) {
    counts <- as.data.frame(table(cells$phenotype, cells$cell_region),
                            stringsAsFactors = FALSE)
    colnames(counts) <- c("phenotype", "region", "count")
    counts$region <- factor(counts$region, levels = c("core", "peri"))
  } else {
    counts <- as.data.frame(table(cells$phenotype), stringsAsFactors = FALSE)
    colnames(counts) <- c("phenotype", "count")
    counts$region <- "all"
  }

  # Order by total count descending
  total_by_pheno <- tapply(counts$count, counts$phenotype, sum)
  counts$phenotype <- factor(counts$phenotype,
                              levels = names(sort(total_by_pheno, decreasing = TRUE)))

  # Colors
  pheno_present <- levels(counts$phenotype)
  pal <- stats::setNames(unname(palette[pheno_present]), pheno_present)
  pal[is.na(pal)] <- "#CCCCCC"

  if ("cell_region" %in% colnames(cells) && length(unique(counts$region)) > 1) {
    ggplot2::ggplot(counts, ggplot2::aes(x = phenotype, y = count, fill = phenotype, colour = region)) +
      ggplot2::geom_col(position = "stack", linewidth = 0.5) +
      ggplot2::scale_fill_manual(values = pal, guide = "none") +
      ggplot2::scale_colour_manual(values = c("core" = "grey20", "peri" = "white"), name = "Region") +
      ggplot2::scale_x_discrete(labels = phenotype_text) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Cell count") +
      theme_islet() +
      ggplot2::theme(legend.position = "bottom")
  } else {
    ggplot2::ggplot(counts, ggplot2::aes(x = phenotype, y = count, fill = phenotype)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = pal, guide = "none") +
      ggplot2::scale_x_discrete(labels = phenotype_text) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Cell count") +
      theme_islet()
  }
}
