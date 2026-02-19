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
  "Endocrine"     = "#D4A373",  # tan
  "Acinar"        = "#A8DADC",  # light teal
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

# ── Cell CSV cache (env-based, same pattern as geojson_cache) ──────────
drilldown_cache <- new.env()

# Path to per-islet cell CSVs
CELLS_DIR <- file.path("..", "..", "data", "cells")

#' Check if single-cell drill-down data is available
drilldown_available <- function() {
  dir.exists(CELLS_DIR) && length(list.files(CELLS_DIR, pattern = "\\.csv$")) > 0
}

#' Load single-cell data for a specific islet
#' @param imageid Case/image ID (e.g., "6505")
#' @param islet_key Islet key (e.g., "Islet_284")
#' @return data.frame with columns: X_centroid, Y_centroid, phenotype, cell_region, + markers
load_islet_cells <- function(imageid, islet_key) {
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
    message("[DRILLDOWN] Error loading cells for ", combined_id, ": ", e$message)
    NULL
  })

  if (!is.null(cells)) {
    assign(cache_key, cells, envir = drilldown_cache)
    message("[DRILLDOWN] Cached cells for ", combined_id, " (", nrow(cells), " cells)")
  }

  cells
}

#' Render single-cell drill-down plot over GeoJSON base
#' @param info List with case_id, islet_key, centroid_x, centroid_y
#' @param cells data.frame from load_islet_cells()
#' @param color_by "phenotype" or a marker column name
#' @param show_peri Logical; if FALSE, only show core cells
#' @return ggplot object
render_islet_drilldown_plot <- function(info, cells, color_by = "phenotype", show_peri = TRUE) {
  if (is.null(info) || is.null(cells) || nrow(cells) == 0) return(NULL)

  # Build GeoJSON base plot
  base_plot <- build_segmentation_base_plot(info)
  if (is.null(base_plot)) return(NULL)

  # Filter cells by region
  if (!show_peri && "cell_region" %in% colnames(cells)) {
    cells <- cells[cells$cell_region == "core", , drop = FALSE]
  }
  if (nrow(cells) == 0) return(base_plot)

  # Convert cell centroids from µm to pixels (matching GeoJSON coordinate space)
  cells$x_px <- cells$X_centroid / PIXEL_SIZE_UM
  cells$y_px <- cells$Y_centroid / PIXEL_SIZE_UM

  # Assign shape: filled for core, open for peri
  if ("cell_region" %in% colnames(cells)) {
    cells$pt_shape <- ifelse(cells$cell_region == "core", 16, 1)
  } else {
    cells$pt_shape <- 16
  }

  if (color_by == "phenotype" && "phenotype" %in% colnames(cells)) {
    # Categorical coloring by phenotype
    pheno_present <- sort(unique(cells$phenotype))
    pal <- PHENOTYPE_COLORS[pheno_present]
    # Fill missing phenotypes with gray
    pal[is.na(pal)] <- "#CCCCCC"

    p <- base_plot +
      ggplot2::geom_point(
        data = cells,
        ggplot2::aes(x = x_px, y = y_px, color = phenotype, shape = factor(pt_shape)),
        size = 1.8, alpha = 0.8,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_color_manual(values = pal, name = "Phenotype") +
      ggplot2::scale_shape_manual(values = c("16" = 16, "1" = 1),
                                  labels = c("16" = "Core", "1" = "Peri"),
                                  name = "Region") +
      ggplot2::labs(title = paste(info$islet_key, "- Single Cells"))
  } else if (color_by %in% colnames(cells)) {
    # Continuous coloring by marker expression
    cells$marker_val <- suppressWarnings(as.numeric(cells[[color_by]]))

    p <- base_plot +
      ggplot2::geom_point(
        data = cells,
        ggplot2::aes(x = x_px, y = y_px, color = marker_val, shape = factor(pt_shape)),
        size = 1.8, alpha = 0.8,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_color_viridis_c(option = "inferno", name = color_by, na.value = "gray80") +
      ggplot2::scale_shape_manual(values = c("16" = 16, "1" = 1),
                                  labels = c("16" = "Core", "1" = "Peri"),
                                  name = "Region") +
      ggplot2::labs(title = paste(info$islet_key, "-", color_by))
  } else {
    # Fallback: just plot cells in gray
    p <- base_plot +
      ggplot2::geom_point(
        data = cells,
        ggplot2::aes(x = x_px, y = y_px),
        color = "gray40", size = 1.5, alpha = 0.7,
        inherit.aes = FALSE
      ) +
      ggplot2::labs(title = paste(info$islet_key, "- Single Cells"))
  }

  p +
    ggplot2::theme(
      legend.position = "right",
      legend.key.size = ggplot2::unit(0.4, "cm"),
      legend.text = ggplot2::element_text(size = 8)
    )
}

#' Render phenotype composition summary bar chart
#' @param cells data.frame from load_islet_cells()
#' @return ggplot object
render_drilldown_summary <- function(cells) {
  if (is.null(cells) || nrow(cells) == 0 || !"phenotype" %in% colnames(cells)) {
    return(NULL)
  }

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
  pal <- PHENOTYPE_COLORS[pheno_present]
  pal[is.na(pal)] <- "#CCCCCC"

  if ("cell_region" %in% colnames(cells) && length(unique(counts$region)) > 1) {
    ggplot2::ggplot(counts, ggplot2::aes(x = phenotype, y = count, fill = phenotype, alpha = region)) +
      ggplot2::geom_col(position = "stack") +
      ggplot2::scale_fill_manual(values = pal, guide = "none") +
      ggplot2::scale_alpha_manual(values = c("core" = 1.0, "peri" = 0.5), name = "Region") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Cell count") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")
  } else {
    ggplot2::ggplot(counts, ggplot2::aes(x = phenotype, y = count, fill = phenotype)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = pal, guide = "none") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Cell count") +
      ggplot2::theme_minimal(base_size = 11)
  }
}
