# ---------- Spatial Tab module server ----------
# Exports: spatial_server(id, prepared)
#
# Dependencies:
#   prepared()$comp — must contain peri_prop_*, immune_*, enrich_z_*, leiden_* columns
#   prepared()$neighborhood — raw neighborhood data with leiden_umap_1/2
#   PHENOTYPE_COLORS — from drilldown_helpers.R
#   donor_tissue_available(), get_available_donors(), load_donor_tissue() — from spatial_helpers.R

# Colourings offered for each Spatial UMAP (owner, 2026-10-21: "multiple choices for coloring the umaps including
# donor demographics and disease state"). Values name the source: "leiden" (the resolution in the Data controls),
# "donor" and "donor_status", or a column of prepared()$comp; spatial_umap_color_choices() drops the columns a
# root lacks.
SPATIAL_UMAP_COLOR_CHOICES <- list(
  "Clusters" = c("Leiden cluster" = "leiden"),
  "Disease state" = c("Donor status" = "donor_status", "Autoantibodies (number positive)" = "n_aab"),
  "Donor" = c("Donor" = "donor", "Age (years)" = "age", "Age group" = "age_group", "Sex" = "gender",
              "Race/ethnicity" = "race_ethnicity", "Height (cm)" = "height_cm", "Weight (kg)" = "weight_kg",
              "BMI" = "BMI", "BMI percentile" = "BMI_percentile"),
  "Clinical" = c("HbA1c (%)" = "HbA1c_percent", "C-peptide (ng/mL)" = "C_peptide_ng_mL",
                 "T1D duration (years)" = "t1d_duration_years", "Admission glucose (mg/dL)" = "admission_glucose_mg_dL",
                 "Peak glucose (mg/dL)" = "peak_glucose_mg_dL"),
  # From the nPOD export (owner, 2026-10-22; v10 donor_metadata.csv).
  "Donation and tissue" = c("Cause of death" = "cause_of_death", "Recovery type" = "recovery_type",
                            "ICU stay (days)" = "ICU_days", "Organ transit time (min)" = "transit_time_min",
                            "Pancreas weight (g)" = "pancreas_weight_g", "RNA integrity (RIN)" = "RIN",
                            "Blood group (ABO)" = "ABO_group"),
  "Object" = c("Band immune proportion" = "immune_frac_peri", "Core immune proportion" = "immune_frac_core",
               "Diameter (\u00b5m)" = "islet_diam_um", "Core cells" = "total_cells_core")
)

#' SPATIAL_UMAP_COLOR_CHOICES without the column-backed choices `comp_cols` lacks.
spatial_umap_colour_choices <- function(comp_cols) {
  available <- function(v) v %in% c("leiden", "donor") | (v == "donor_status" & "Donor Status" %in% comp_cols) |
    v %in% comp_cols
  groups <- lapply(SPATIAL_UMAP_COLOR_CHOICES, function(g) g[available(unname(g))])
  groups[lengths(groups) > 0]
}

# Wrap display text only; scale breaks, cell labels and colors keep their source names.
spatial_tissue_legend_labels <- function(labels) {
  vapply(phenotype_text(labels), function(label) paste(strwrap(label, width = 28), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}

# Two-column phenotype legends need room for each row's longest wrapped label.
# Keep roughly 500 px for the map plus its axes/title instead of shrinking it
# when the user shows every macrophage subtype. Cluster legends keep the floor.
spatial_tissue_plot_height <- function(plot) {
  scale <- plot$scales$get_scales("colour")
  if (is.null(scale) || !is.character(scale$breaks) || !length(scale$breaks)) return(800L)
  labels <- scale$get_labels(scale$breaks)
  lines <- lengths(strsplit(labels, "\n", fixed = TRUE))
  row_lines <- vapply(split(lines, ceiling(seq_along(lines) / 2)), max, integer(1))
  legend_height <- 40 + sum(pmax(28, ISLET_FONT_SIZES$legend_text * 1.25 * row_lines) + 8)
  as.integer(ceiling(max(800, 600 + legend_height)))
}

spatial_server <- function(id, prepared, palette = reactive(PHENOTYPE_COLORS),
                           donor_colors_reactive = reactive(DONOR_COLORS_BRIGHT),
                           selected_islet = reactiveVal(NULL), active_tab = reactive("Spatial"),
                           macrophage_markers = reactiveVal(character(0))) {

  moduleServer(id, function(input, output, session) {
    # The macrophage markers (owner, 2026-10-27) sit under the Macrophage heading of each grouped phenotype list
    # (Data controls, composition & enrichment, distance profile); they set the app-wide setting.
    bind_macrophage_marker_inputs(input, session, c("mac_markers_filter", "mac_markers_enrich", "mac_markers_dist"),
                                  macrophage_markers)
    output$embedding_method <- renderText({
      # Since 2026-10-19 (owner) a uniform root serves the endocrine composition of the objects with ten or
      # more core cells (scripts/uniform/composition_embedding.py); older roots serve the protein embedding.
      comp_method <- prepared()$composition_embedding
      if (!is.null(comp_method)) {
        return(sprintf(paste0(
          "Built from each object's endocrine composition: the Endocrine Beta, Alpha, Delta and CHGA⁺ hormone-negative",
          " proportions of its core endocrine cells, for the %s objects with %s or more core cells. UMAP and Leiden",
          " use one neighbor graph. No donor correction: the proportions are counted within each object, so staining",
          " differences between donors do not enter them. Smaller objects have too few cells to place, so they",
          " have no UMAP point or cluster and are drawn gray, as tissue, when the tissue scatter is colored by",
          " cluster."),
          format(as.integer(comp_method$objects), big.mark = ","), as.integer(comp_method$min_core_cells)))
      }
      method <- prepared()$protein_embedding
      if (is.null(method)) return("")
      correction <- if (identical(as.character(method$donor_correction), "harmony")) {
        " Donor-adjusted with Harmony; donor and biological effects may overlap."
      } else {
        " No donor correction; clusters may reflect donor differences."
      }
      paste0("Built from each object's mean marker values: protein PCA, with UMAP and Leiden on the same marker",
             " graph.", correction)
    })
    # Rendered inside the closed Islet-Level UMAPs help popover, so keep it current while hidden.
    outputOptions(output, "embedding_method", suspendWhenHidden = FALSE)
    ns <- session$ns

    # 20-color qualitative palette for Leiden clusters
    leiden_palette <- c(
      "#1f77b4", "#ff7f0e", "#66c2a5", "#d62728", "#8da0cb",
      "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
      "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
      "#9467bd", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"
    )

    # ---- Helpers: column identification ----
    get_nbr_columns <- function(comp_cols) {
      list(
        peri_prop  = grep("^peri_prop_", comp_cols, value = TRUE),
        peri_count = grep("^peri_count_", comp_cols, value = TRUE),
        immune     = intersect(c("immune_frac_peri", "immune_frac_core", "immune_ratio",
                                 "cd8_to_macro_ratio", "tcell_density_peri"), comp_cols),
        enrich     = grep("^enrich_z_", comp_cols, value = TRUE),
        distance   = grep("^min_dist_", comp_cols, value = TRUE),
        leiden     = grep("^leiden_[0-9]", comp_cols, value = TRUE),
        leiden_umap = intersect(c("leiden_umap_1", "leiden_umap_2"), comp_cols)
      )
    }

    has_neighborhood <- reactive({
      pd <- prepared()
      if (is.null(pd$comp)) return(FALSE)
      nbr <- get_nbr_columns(colnames(pd$comp))
      length(nbr$peri_prop) > 0
    })

    has_leiden <- reactive({
      pd <- prepared()
      if (is.null(pd$comp)) return(FALSE)
      nbr <- get_nbr_columns(colnames(pd$comp))
      length(nbr$leiden) > 0
    })

    has_tissue <- reactive({
      donor_tissue_available()
    })

    # ---- Donor selector ----
    output$donor_selector <- renderUI({
      donors <- get_available_donors()
      if (length(donors) == 0) {
        return(div(style = "color: var(--islet-fg-subtle); font-style: italic; padding-top: 8px;",
                   "No per-donor tissue data available."))
      }
      pd <- prepared()
      # Build friendly labels: imageid (Donor Status)
      donor_labels <- sapply(donors, function(d) {
        if (!is.null(pd$comp) && "Donor Status" %in% colnames(pd$comp) &&
            "Case ID" %in% colnames(pd$comp)) {
          ds_match <- pd$comp$`Donor Status`[pd$comp$`Case ID` == as.integer(d)]
          ds <- if (length(ds_match) > 0) ds_match[1] else "?"
          paste0(d, " (", ds, ")")
        } else {
          d
        }
      })
      choices <- setNames(donors, donor_labels)
      selectInput(ns("donor"), "Donor", choices = choices, selected = donors[1])
    })

    # Browser input is accepted only when it is one of the donor partitions
    # discovered by the server. This check happens before the value reaches a
    # filesystem path or download filename.
    selected_donor <- reactive({
      donor_id <- input$donor
      if (is.null(donor_id)) return(NULL)
      if (!islet_valid_scalar_choice(donor_id, get_available_donors(),
                                     "spatial_donor_invalid")) return(NULL)
      donor_id
    })

    # ---- Leiden resolution selector ----
    output$leiden_res_selector <- renderUI({
      if (!has_leiden()) return(NULL)
      nbr <- get_nbr_columns(colnames(prepared()$comp))
      # Available resolutions (e.g., leiden_0.3 -> "0.3")
      res_labels <- gsub("^leiden_", "", nbr$leiden)
      choices <- setNames(nbr$leiden, res_labels)
      selectInput(ns("leiden_res"), "Leiden Resolution",
                  choices = choices, selected = nbr$leiden[2] %||% nbr$leiden[1])
    })

    # ---- Donor cells reactive (load tissue CSV for selected donor) ----
    donor_cells <- reactive({
      donor_id <- selected_donor()
      req(donor_id)
      cells <- load_donor_tissue(donor_id)
      # macrophages merged for the checked markers (owner, 2026-10-27)
      if (!is.null(cells) && "phenotype" %in% names(cells))
        cells$phenotype <- macrophage_collapse(cells$phenotype, macrophage_markers())
      cells
    })

    # ---- Donor status for selected donor ----
    donor_status <- reactive({
      donor_id <- selected_donor()
      req(donor_id)
      pd <- prepared()
      if (!is.null(pd$comp) && "Donor Status" %in% colnames(pd$comp) &&
          "Case ID" %in% colnames(pd$comp)) {
        ds <- pd$comp$`Donor Status`[pd$comp$`Case ID` == as.integer(donor_id)]
        if (length(ds) > 0) return(ds[1])
      }
      "?"
    })

    # Leiden mapping for tissue scatter: map islet_name -> cluster for core/peri cells
    islet_leiden_map <- reactive({
      req(has_leiden())
      pd <- prepared()
      leiden_col <- input$leiden_res
      req(leiden_col, leiden_col %in% colnames(pd$comp))

      # Build islet_key -> cluster mapping for the selected donor
      donor_id <- selected_donor()
      req(donor_id)

      comp <- pd$comp
      donor_mask <- comp$`Case ID` == as.integer(donor_id)
      if (!any(donor_mask)) return(NULL)

      sub <- comp[donor_mask, c("islet_key", leiden_col), drop = FALSE]
      setNames(as.character(sub[[leiden_col]]), as.character(sub$islet_key))
    })

    # ---- Phenotype filter ----
    #
    # Choices come from the CANONICAL vocabulary (`canonical_phenotype_order()`,
    # backed by PHENOTYPE_COLORS), never from `unique(donor_cells()$phenotype)`.
    # Deriving them per-donor is why "CD4 Tcell" was missing: the default donor
    # (112) has zero CD4 cells, as do 6533 and 6534, so the checkbox simply did
    # not exist on first load. Phenotypes with no cells for the current donor are
    # still listed, annotated with their count, so the list never changes shape
    # when the donor changes.
    PHENO_HIDE_DEFAULT <- c("Acinar", "Ductal", "Endocrine", "T cell", "Unknown")

    # Canonical choices + this donor's per-phenotype counts.
    pheno_choice_info <- reactive({
      cells <- tryCatch(donor_cells(), error = function(e) NULL)
      observed <- if (is.null(cells)) character(0) else cells$phenotype
      phenos <- phenotype_partition_labels(macrophage_collapse(canonical_phenotype_order(observed), macrophage_markers()))
      counts <- table(factor(observed, levels = phenos))
      list(phenos = phenos, counts = counts)
    })

    output$phenotype_filter <- renderUI({
      info <- pheno_choice_info()
      phenos <- info$phenos
      if (length(phenos) == 0) return(NULL)
      counts <- info$counts

      # First render applies the hide-by-default set; afterwards the user's
      # selection is preserved verbatim. The old code did
      # `intersect(prev, phenos)` against a per-donor choice list, which left any
      # newly-appearing phenotype permanently unchecked. With a donor-independent
      # list there is nothing to intersect away.
      prev <- isolate(input$pheno_filter)
      selected_phenos <- if (!is.null(prev)) remap_macrophage_selection(prev, phenos)
                         else setdiff(phenos, PHENO_HIDE_DEFAULT)

      # Grouped, with All | None per group and for the list (owner, 2026-10-27). Labels carry the donor's count, so
      # a zero-cell phenotype is visibly present but empty rather than silently absent.
      phenotype_group_checkboxes(ns("pheno_filter"), phenos, selected = selected_phenos, counts = counts,
                                 marker_input_id = if (root_has_marker_profile_labels()) ns("mac_markers_filter"),
                                 markers = isolate(macrophage_markers()))
    })

    # ---- Zoom state for tissue scatter ----
    scatter_zoom <- reactiveValues(xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL)

    observeEvent(input$scatter_brush, {
      brush <- input$scatter_brush
      if (!is.null(brush)) {
        scatter_zoom$xmin <- brush$xmin
        scatter_zoom$xmax <- brush$xmax
        # y is reversed in the plot, so brush coords are already in data space
        scatter_zoom$ymin <- brush$ymin
        scatter_zoom$ymax <- brush$ymax
      }
    })

    observeEvent(input$scatter_dblclick, {
      scatter_zoom$xmin <- NULL
      scatter_zoom$xmax <- NULL
      scatter_zoom$ymin <- NULL
      scatter_zoom$ymax <- NULL
    })

    observeEvent(input$scatter_reset_zoom, {
      scatter_zoom$xmin <- NULL
      scatter_zoom$xmax <- NULL
      scatter_zoom$ymin <- NULL
      scatter_zoom$ymax <- NULL
    })

    # Reset zoom when donor changes
    observeEvent(input$donor, {
      scatter_zoom$xmin <- NULL
      scatter_zoom$xmax <- NULL
      scatter_zoom$ymin <- NULL
      scatter_zoom$ymax <- NULL
    })

    # ==== Card 2: Tissue Scatter Plot (ggplot2, rasterized) ====
    spatial_tissue_scatter_plot <- reactive({
      cells <- donor_cells()
      req(cells)

      # Copy so we can filter
      plot_df <- cells

      # Region filter: highlight selected region, dim others
      region_mode <- input$region_filter %||% "all"

      # Determine coloring mode
      color_mode <- input$color_by %||% "phenotype"
      use_leiden <- (color_mode == "leiden") && has_leiden()

      if (use_leiden) {
        # Map islet_name -> leiden cluster
        lmap <- islet_leiden_map()
        if (!is.null(lmap)) {
          plot_df$cluster <- lmap[plot_df$islet_name]
          # Tissue cells with no islet get NA cluster
          plot_df$cluster[is.na(plot_df$cluster)] <- "tissue"
        } else {
          plot_df$cluster <- "tissue"
        }
      }

      # Phenotype filter (only in phenotype mode)
      selected_phenos <- phenotype_partition_labels(input$pheno_filter)
      if (!use_leiden && !is.null(selected_phenos) && length(selected_phenos) > 0) {
        plot_df <- plot_df[plot_df$phenotype %in% selected_phenos, , drop = FALSE]
      }

      # Split into foreground (highlighted) and background (dimmed) layers
      if (region_mode == "all") {
        # All cells colored; tissue background slightly dimmed
        fg <- plot_df[plot_df$cell_region %in% c("core", "peri"), , drop = FALSE]
        bg <- plot_df[plot_df$cell_region == "tissue", , drop = FALSE]
      } else if (region_mode == "core_peri") {
        fg <- plot_df[plot_df$cell_region %in% c("core", "peri"), , drop = FALSE]
        bg <- plot_df[plot_df$cell_region == "tissue", , drop = FALSE]
      } else {
        # Core only
        fg <- plot_df[plot_df$cell_region == "core", , drop = FALSE]
        bg <- plot_df[plot_df$cell_region != "core", , drop = FALSE]
      }

      p <- ggplot2::ggplot()

      color_bg <- isTRUE(input$color_background)

      # Background layer
      if (nrow(bg) > 0) {
        if (color_bg && !use_leiden) {
          # Color background cells by phenotype (dimmed)
          bg_phenos <- sort(unique(bg$phenotype))
          bg_pal <- stats::setNames(unname(palette()[bg_phenos]), bg_phenos)
          bg_pal[is.na(bg_pal)] <- "#CCCCCC"
          p <- p + ggplot2::geom_point(
            data = bg,
            ggplot2::aes(x = X_centroid, y = Y_centroid, color = phenotype),
            size = 0.15, alpha = 0.25,
            inherit.aes = FALSE
          )
        } else {
          p <- p + ggplot2::geom_point(
            data = bg,
            ggplot2::aes(x = X_centroid, y = Y_centroid),
            color = "#d9d9d9", size = 0.15, alpha = 0.3,
            inherit.aes = FALSE
          )
        }
      }

      # Foreground layer: colored by phenotype or leiden
      if (nrow(fg) > 0) {
        if (use_leiden && "cluster" %in% colnames(fg)) {
          cluster_levels <- leiden_levels(fg$cluster[fg$cluster != "tissue"])
          fg$cluster <- factor(fg$cluster, levels = c(cluster_levels, "tissue"))
          n_clusters <- length(cluster_levels)
          pal_colors <- rep_len(leiden_palette, n_clusters)
          pal <- setNames(pal_colors, cluster_levels)
          pal["tissue"] <- "#d9d9d9"

          p <- p + ggplot2::geom_point(
            data = fg,
            ggplot2::aes(x = X_centroid, y = Y_centroid, color = cluster),
            size = 0.4, alpha = 0.6,
            inherit.aes = FALSE
          ) +
          ggplot2::scale_color_manual(values = pal, name = "Cluster", na.value = "#d9d9d9",
                                        guide = ggplot2::guide_legend(ncol = 6, byrow = TRUE,
                                                                     override.aes = list(size = 4)))
        } else {
          # Phenotype coloring — legend in canonical palette order
          canonical_order <- names(palette())
          if (!is.null(selected_phenos) && length(selected_phenos) > 0) {
            present <- intersect(selected_phenos,
                                 unique(c(fg$phenotype, if (color_bg) bg$phenotype)))
          } else {
            present <- unique(c(fg$phenotype, if (color_bg) bg$phenotype))
          }
          legend_phenos <- intersect(canonical_order, present)
          legend_phenos <- c(legend_phenos, setdiff(present, legend_phenos))
          legend_phenos <- legend_phenos[order(match(phenotype_group(legend_phenos), PHENOTYPE_GROUP_ORDER))]
          pal <- stats::setNames(unname(palette()[legend_phenos]), legend_phenos)
          pal[is.na(pal)] <- "#CCCCCC"

          p <- p + ggplot2::geom_point(
            data = fg,
            ggplot2::aes(x = X_centroid, y = Y_centroid, color = phenotype),
            size = 0.4, alpha = 0.6,
            inherit.aes = FALSE
          ) +
          ggplot2::scale_color_manual(values = pal, breaks = legend_phenos, labels = spatial_tissue_legend_labels,
                                        name = "Phenotype", na.value = "#CCCCCC",
                                        guide = ggplot2::guide_legend(ncol = 2, byrow = TRUE,
                                                                     override.aes = list(size = 4)))
        }
      }

      ds <- donor_status()

      # Zoom via coord_cartesian (clips display without dropping data)
      zoomed <- !is.null(scatter_zoom$xmin)
      if (zoomed) {
        # Brush coords are in display space (y already reversed by scale_y_reverse).
        # coord_cartesian ylim is in *data* space (before reversal), so swap y.
        p <- p + ggplot2::scale_y_reverse() +
          ggplot2::coord_cartesian(
            xlim = c(scatter_zoom$xmin, scatter_zoom$xmax),
            ylim = sort(c(scatter_zoom$ymin, scatter_zoom$ymax))
          )
      } else {
        p <- p + ggplot2::coord_fixed() +
          ggplot2::scale_y_reverse()
      }

      donor_id <- selected_donor()
      req(donor_id)
      p + ggplot2::labs(
          title = paste0("Donor ", donor_id, " (", ds, ")"),
          subtitle = paste0(nrow(fg), " foreground cells | ",
                            nrow(bg), " background cells"),
          x = expression(paste("X centroid (", mu, "m)")),
          y = expression(paste("Y centroid (", mu, "m)"))
        ) +
        theme_islet() +
        ggplot2::theme(
          # Keep the rightmost x-axis tick inside the narrow tissue PNG.
          plot.margin = ggplot2::margin_part(r = 24),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title.position = "top",
          legend.key.size = ggplot2::unit(0.7, "cm")
        )
    })

    # Cached across sessions (Plan A Level 2, B3): the 170 k-cell PNG of a donor and drawing options is drawn once
    # per process. The key is everything spatial_tissue_scatter_plot() reads; the output size is part of it.
    tissue_scatter_height <- reactive(spatial_tissue_plot_height(spatial_tissue_scatter_plot()))
    output$tissue_scatter <- bindCache(
      shiny::renderPlot({ spatial_tissue_scatter_plot() }, height = function() tissue_scatter_height(), bg = "transparent"),
      selected_donor(), donor_status(), input$region_filter, input$color_by, input$pheno_filter,
      input$color_background, palette(), has_leiden(), macrophage_markers(), tissue_scatter_height(),
      if (identical(input$color_by, "leiden") && has_leiden()) list(input$leiden_res, islet_leiden_map()) else NULL,
      scatter_zoom$xmin, scatter_zoom$xmax, scatter_zoom$ymin, scatter_zoom$ymax, islet_theme_mode(),
      # bindCache supplies its own render dimensions; the auto-height container
      # initially reports zero. Keep width bucketing, but use the calculated
      # height instead of feeding that transient client height back into the PNG.
      sizePolicy = function(size) c(shiny::sizeGrowthRatio()(c(size[1], 800))[1], tissue_scatter_height()))
    output$dl_tissue_scatter <- shiny::downloadHandler(
      filename = function() {
        donor_id <- selected_donor()
        islet_png_filename(paste0("spatial_tissue_donor_", donor_id %||% "all"))
      },
      content = function(file) {
        p <- spatial_tissue_scatter_plot(); if (is.null(p)) return()
        islet_ggsave(file, plot = p,
                        width = 12, height = 10,
                        dpi = ISLET_PNG_DIMS$dpi, units = ISLET_PNG_DIMS$units)
      }
    )

    # ==== Card 3: Leiden UMAP (plotly; the objects in the embedding) ====
    output$leiden_not_available <- renderUI({
      if (!has_leiden()) {
        return(div(style = "color: var(--islet-fg-subtle); font-style: italic; margin-bottom: 10px;",
          "Leiden clustering not available in current H5AD. ",
          "Rebuild with build_h5ad_for_app.py after running Leiden clustering."))
      }
      NULL
    })

    # One islet-level UMAP (points from islet_umap_plotly(), interleaved so no group hides another): the two
    # Spatial UMAPs share their marker size and full-width layout. Legends and numeric color bars sit below
    # the axes so either coloring mode can use the card's full width for the embedding.
    # Both are linked plots (trajectory_link_plot(), keyed by `.key`).
    spatial_umap_plotly <- function(plot_df, colour_col, pal, legend_title, title, filename_root,
                                    colorbar_title = NULL) {
      fonts <- plotly_islet_fonts()
      if (is.null(plot_df$.hover)) plot_df$.hover <- paste0(legend_title, ": ", plot_df[[colour_col]], "<br>",
                               if (colour_col != "donor_status") paste0("Status: ", plot_df$donor_status, "<br>") else "",
                               plot_df$islet_key)
      # These titles previously wrapped beside a vertical bar. A horizontal
      # bar needs one line so the second line does not cross its gradient.
      if (!is.null(colorbar_title))
        colorbar_title <- gsub("<br\\s*/?>", " ", colorbar_title, ignore.case = TRUE, perl = TRUE)
      p <- islet_umap_plotly(plot_df, "umap1", "umap2", colour_col, pal, ".hover", key = ".key",
                             colorbar_title = colorbar_title)
      # Numeric colors use the same lower strip as categorical legends.
      if (is.numeric(plot_df[[colour_col]])) {
        p <- plotly::plotly_build(p)
        p$x$data <- lapply(p$x$data, function(tr) {
          if (!is.null(tr$marker$colorbar))
            tr$marker$colorbar <- utils::modifyList(tr$marker$colorbar, list(
              orientation = "h", x = 0.5, xanchor = "center", y = -0.3, yanchor = "bottom",
              len = 0.8, thickness = 14, title = list(side = "top"),
              tickfont = list(size = ISLET_FONT_SIZES$dense_min)))
          tr
        })
      }
      p <- p %>%
        layout(
          font = fonts$global,
          title = list(text = title, font = fonts$title),
          margin = list(l = 30, r = 10, t = 50, b = 110, autoexpand = FALSE),
          xaxis = list(title = list(text = "UMAP 1", font = fonts$axis_title), domain = c(0, 1),
                       zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, ticks = "",
                       showline = FALSE, constrain = "domain"),
          yaxis = list(title = list(text = "UMAP 2", font = fonts$axis_title),
                       zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, ticks = "",
                       showline = FALSE, scaleanchor = "x", scaleratio = 1,
                       constrain = "domain"),
          legend = list(title = list(text = legend_title, font = fonts$legend_title, side = "top"),
                        orientation = "h", font = fonts$legend_text,
                        x = 0.5, xanchor = "center", y = -0.3, yanchor = "bottom",
                        itemsizing = "constant", tracegroupgap = 2)
        ) %>%
        plotly_islet_config(filename_root)
      p <- trajectory_link_plot(p, ns)
      # Plotly otherwise expands each plot's lower margin independently. Measure
      # the rendered guides (including wrapped donor legends) and reserve only
      # the space needed by the taller one, equally on both plots. Afterplot also
      # runs on resize and coloring changes; equal margins need no further draw.
      htmlwidgets::onRender(p, "function(el, x, ids) {
        var prior = el._isletGuideSync;
        if (prior) {
          el.removeListener('plotly_afterplot', prior.schedule);
          if (prior.frame) window.cancelAnimationFrame(prior.frame);
        }
        var state = {frame: null};
        state.schedule = function() {
          if (state.frame) return;
          state.frame = window.requestAnimationFrame(function() {
            state.frame = null;
            var plots = ids.map(function(id) { return document.getElementById(id); });
            if (plots.some(function(plot) { return !plot || !plot._fullLayout; })) return;
            var tallest = 0;
            var clearance = 40;
            plots.forEach(function(plot) {
              plot.querySelectorAll('g.colorbar').forEach(function(guide) {
                var title = guide.querySelector('g.cbtitle > text');
                var fill = guide.querySelector('.cbfill');
                if (!title || !fill) return;
                // Reset first: repeated draws must not accumulate an offset.
                title.setAttribute('dy', '0');
                var gap = fill.getBoundingClientRect().top - title.getBoundingClientRect().bottom;
                title.setAttribute('dy', String(Math.min(0, gap - 6)));
              });
              plot.querySelectorAll('g.legend, g.colorbar').forEach(function(guide) {
                tallest = Math.max(tallest, guide.getBoundingClientRect().height);
              });
              var axisTitle = plot.querySelector('.g-xtitle');
              if (axisTitle) {
                var layout = plot._fullLayout;
                var paperBottom = plot.getBoundingClientRect().top + layout._size.t + layout._size.h;
                // Keep the guide clear of the actual axis title too; font and
                // title offsets differ between the supported Plotly versions.
                clearance = Math.max(clearance,
                  axisTitle.getBoundingClientRect().bottom - paperBottom + 8 + 0.02 * layout.height);
              }
            });
            var bottom = Math.ceil(Math.max(110, tallest + clearance));
            plots.forEach(function(plot) {
              var layout = plot._fullLayout;
              // Both supported Plotly versions default to paper coordinates;
              // Plotly.js 2.11 (the pinned image) has no yref attribute.
              // Use paper coordinates to put each guide bottom 2% above the
              // widget bottom, using the height AFTER the shared margin.
              var height = layout.height - layout.margin.t - bottom;
              if (height <= 0) return;
              var y = (0.02 * layout.height - bottom) / height;
              var update = {};
              if (Math.abs(layout.margin.b - bottom) > 0.5) update['margin.b'] = bottom;
              if (layout.legend && Math.abs(layout.legend.y - y) > 0.000001)
                update['legend.y'] = y;
              var bars = [];
              plot.data.forEach(function(trace, index) {
                if (trace.marker && trace.marker.colorbar &&
                    Math.abs(trace.marker.colorbar.y - y) > 0.000001) bars.push(index);
              });
              if (bars.length)
                Plotly.update(plot, {'marker.colorbar.y': y}, update, bars);
              else if (Object.keys(update).length) Plotly.relayout(plot, update);
            });
          });
        };
        el._isletGuideSync = state;
        el.on('plotly_afterplot', state.schedule);
        state.schedule();
      }", data = c(ns("leiden_umap"), ns("spatial_umap_immune")))
    }

    current_leiden_col <- function(comp) {
      input$leiden_res %||% {
        nbr <- get_nbr_columns(colnames(comp))
        nbr$leiden[2] %||% nbr$leiden[1]
      }
    }

    # The objects in the embedding (UMAP coordinates), for the checked donor statuses.
    spatial_umap_rows <- reactive({
      req(has_leiden())
      comp <- nbr_comp()
      req("leiden_umap_1" %in% colnames(comp), "leiden_umap_2" %in% colnames(comp))
      umap1 <- suppressWarnings(as.numeric(comp$leiden_umap_1))
      umap2 <- suppressWarnings(as.numeric(comp$leiden_umap_2))
      keep <- is.finite(umap1) & is.finite(umap2)
      comp <- comp[keep, , drop = FALSE]
      comp$.umap1 <- umap1[keep]
      comp$.umap2 <- umap2[keep]
      comp$.key <- trajectory_link_key(comp$`Case ID`, comp$islet_key)
      comp
    })

    # Offer only the colourings the served root has (the UI lists them all).
    observe({
      choices <- spatial_umap_colour_choices(colnames(prepared()$comp))
      if (setequal(unlist(choices), unlist(SPATIAL_UMAP_COLOR_CHOICES))) return()
      for (id in c("umap_color_left", "umap_color_right")) {
        current <- isolate(input[[id]])
        fallback <- if (id == "umap_color_left") "leiden" else
          intersect(c("immune_frac_peri", "donor_status"), unlist(choices))[1]
        updateSelectInput(session, id, choices = choices,
                          selected = if (!is.null(current) && current %in% unlist(choices)) current else fallback)
      }
    })

    # One UMAP figure colored by `key` (a SPATIAL_UMAP_COLOR_CHOICES value).
    spatial_umap_figure <- function(key, filename_root) {
      comp <- spatial_umap_rows()
      allowed <- unlist(unname(spatial_umap_colour_choices(colnames(comp))))  # names are the labels
      if (!islet_valid_scalar_choice(key, allowed, "spatial_umap_colour_invalid")) key <- "leiden"
      if (nrow(comp) == 0) return(trajectory_link_plot(plotly_islet_empty("No UMAP data"), ns))
      label <- names(allowed)[match(key, allowed)]
      status <- as.character(comp$`Donor Status`)
      plot_df <- data.frame(umap1 = comp$.umap1, umap2 = comp$.umap2, donor_status = status,
                            islet_key = as.character(comp$islet_key), .key = comp$.key, stringsAsFactors = FALSE)
      donor <- paste0(comp$`Case ID`, " (", status, ")")
      object_line <- paste0("Donor ", donor, "<br>", plot_df$islet_key)
      title_n <- function(t) paste0(t, " (", format(nrow(plot_df), big.mark = ","), " objects)")
      discrete <- function(values, levels, pal, legend_title, title = title_n(label)) {
        plot_df$colour <- factor(values, levels = levels)
        plot_df$.hover <- paste0(legend_title, ": ", ifelse(is.na(values), "not recorded", values), "<br>",
                                 object_line)
        spatial_umap_plotly(plot_df, "colour", pal, legend_title, title, filename_root)
      }
      continuous <- function(values, legend_title, colorbar_title, cap = FALSE, log_scale = FALSE, digits = 3,
                             missing_text = "not recorded") {
        values <- suppressWarnings(as.numeric(values))
        values[!is.finite(values)] <- NA
        colour <- values
        if (cap && any(!is.na(values))) {
          # Skewed (most bands hold no immune cell): saturate at the 95th percentile; hover shows the true value.
          top <- stats::quantile(values, 0.95, names = FALSE, na.rm = TRUE)
          colour <- pmin(values, if (top > 0) top else max(values, na.rm = TRUE))
        }
        if (log_scale) colour <- log10(pmax(values, 1))
        plot_df$colour <- colour
        plot_df$.hover <- paste0(legend_title, ": ", ifelse(is.na(values), missing_text, signif(values, digits)),
                                 "<br>", object_line)
        spatial_umap_plotly(plot_df, "colour", NULL, legend_title, title_n(legend_title), filename_root,
                            colorbar_title = colorbar_title)
      }

      switch(key,
        leiden = {
          leiden_col <- current_leiden_col(comp)
          req(leiden_col %in% colnames(comp))
          cluster <- as.character(comp[[leiden_col]])
          lv <- leiden_levels(cluster)
          discrete(cluster, lv, setNames(rep_len(leiden_palette, length(lv)), lv), "Cluster",
                   title_n(paste("Leiden", gsub("^leiden_", "", leiden_col))))
        },
        donor_status = {
          lv <- intersect(c("ND", "Aab+", "T1D"), status)
          discrete(status, lv, donor_colors_reactive()[lv], "Donor status")
        },
        donor = {
          ids <- unique(data.frame(id = comp$`Case ID`, status = status, label = donor))
          ids <- ids[order(match(ids$status, c("ND", "Aab+", "T1D")), suppressWarnings(as.numeric(ids$id))), ]
          discrete(donor, ids$label, setNames(rep_len(leiden_palette, nrow(ids)), ids$label), "Donor")
        },
        n_aab = {
          n <- suppressWarnings(as.integer(round(as.numeric(comp$n_aab))))
          lv <- as.character(sort(unique(n[!is.na(n)])))
          discrete(as.character(n), lv, setNames(viridisLite::viridis(length(lv), begin = 0.15, end = 0.95), lv),
                   "Autoantibodies")
        },
        age_group = {
          g <- as.character(comp$age_group)
          lv <- intersect(c("<13", ">=13"), g)
          lv <- c(lv, setdiff(sort(unique(g[!is.na(g)])), lv))
          discrete(g, lv, setNames(rep_len(c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"), length(lv)), lv),
                   "Age group")
        },
        race_ethnicity = , recovery_type = , cause_of_death = , ABO_group = {
          g <- as.character(comp[[key]])
          g[!nzchar(g)] <- NA
          lv <- sort(unique(g[!is.na(g)]))
          discrete(g, lv, setNames(rep_len(c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2",
                                             "#D55E00"), length(lv)), lv), label)
        },
        gender = {
          g <- as.character(comp$gender)
          g <- ifelse(g %in% c("F", "Female"), "Female", ifelse(g %in% c("M", "Male"), "Male", NA))
          discrete(g, intersect(c("Female", "Male"), g), c(Female = "#CC79A7", Male = "#009E73"), "Sex")
        },
        HbA1c_percent = {
          v <- suppressWarnings(as.numeric(comp$HbA1c_percent))
          v[!is.na(v) & v <= 0] <- NA  # two donors' HbA1c is recorded as 0: not measured
          continuous(v, label, "HbA1c (%)")
        },
        immune_frac_peri = continuous(comp$immune_frac_peri, label, "Immune<br>proportion", cap = TRUE),
        immune_frac_core = continuous(comp$immune_frac_core, label, "Immune<br>proportion", cap = TRUE),
        total_cells_core = continuous(comp$total_cells_core, label, "Core cells<br>(log10)", log_scale = TRUE),
        islet_diam_um = continuous(comp$islet_diam_um, "Diameter (\u00b5m)", "Diameter<br>(\u00b5m)"),
        age = continuous(comp$age, label, "Age<br>(years)"),
        BMI = continuous(comp$BMI, label, "BMI"),
        C_peptide_ng_mL = continuous(comp$C_peptide_ng_mL, label, "C-peptide<br>(ng/mL)"),
        # ND and Aab+ donors have no diabetes, hence no duration (gray).
        t1d_duration_years = continuous(comp$t1d_duration_years, label, "T1D duration<br>(years)",
                                        missing_text = ifelse(status == "T1D", "not recorded", "no T1D")),
        continuous(comp[[key]], label, sub(" (", "<br>(", label, fixed = TRUE))
      )
    }

    # Output ids keep their earlier names (the browser check and downloads name them); each now has a Color by.
    output$leiden_umap <- renderPlotly({
      spatial_umap_figure(input$umap_color_left %||% "leiden", "spatial_umap_left")
    })
    output$spatial_umap_immune <- renderPlotly({
      spatial_umap_figure(input$umap_color_right %||% "immune_frac_peri", "spatial_umap_right")
    })

    # ---- Linked UMAPs: one validated selection route (hover stays in the browser, www/trajectory-links.js) ----
    spatial_link_selection_key <- function() {
      info <- selected_islet()
      if (is.null(info)) return(NULL)
      key <- trajectory_link_key(info$case_id, info$islet_key)
      if (length(key) != 1L || is.na(key)) NULL else key
    }
    send_spatial_link_state <- function() {
      session$sendCustomMessage(ns("linked_state"), list(
        key = spatial_link_selection_key(), active = identical(active_tab(), "Spatial")))
    }
    observeEvent(input$linked_click, {
      req(identical(active_tab(), "Spatial"))
      event <- input$linked_click
      key <- event$key
      rows <- tryCatch(spatial_umap_rows(), error = function(e) NULL)
      if (!is.character(key) || length(key) != 1L || is.na(key) ||
          !is.character(event$source) || length(event$source) != 1L ||
          !(event$source %in% ns(c("leiden_umap", "spatial_umap_immune"))) || is.null(rows) ||
          !(key %in% rows$.key)) {
        islet_log_validation_failure("spatial_link_invalid")
        send_spatial_link_state()
        return()
      }
      row <- rows[match(key, rows$.key), , drop = FALSE]
      select_islet_from_key(as.character(row$`Case ID`), as.character(row$islet_key), prepared, selected_islet)
      send_spatial_link_state()
    }, ignoreInit = TRUE)
    observeEvent(input$linked_ready, { send_spatial_link_state() }, ignoreInit = TRUE)
    observe({ send_spatial_link_state() })
    observeEvent(input$clear_link_selection, {
      selected_islet(NULL)
      send_spatial_link_state()
    })

    # The segmentation map of the selected object, under the UMAPs. Only while this tab is active: its outputs are
    # root-level and shared with the Plot and Trajectory panels.
    output$segmentation_viewer_panel <- renderUI({
      if (!identical(active_tab(), "Spatial")) return(NULL)
      info <- selected_islet()
      if (is.null(info)) return(NULL)
      segmentation_viewer_panel_ui(info, ns("clear_segmentation"), layout = "wide")
    })
    observeEvent(input$clear_segmentation, { selected_islet(NULL) })

    # ==== Cluster card, Composition view: each Leiden cluster's mean core phenotype proportions ====
    # The owner kept it (2026-10-20: "There still needs to be a cluster composition plot"): it is how a cluster is
    # read, although its endocrine bars restate what the clusters are built from. The Membership views sit beside it.
    cluster_composition_plot <- function() {
      req(has_leiden())
      pd <- prepared()
      comp <- pd$comp

      leiden_col <- input$leiden_res %||% {
        nbr <- get_nbr_columns(colnames(comp))
        nbr$leiden[2] %||% nbr$leiden[1]
      }
      req(leiden_col %in% colnames(comp))

      # Filter by donor status
      plot_comp <- comp
      if (!is.null(input$groups) && "Donor Status" %in% colnames(plot_comp)) {
        plot_comp <- plot_comp[plot_comp$`Donor Status` %in% input$groups, , drop = FALSE]
      }

      # Get phenotype proportion columns (prop_*), filtered by selected phenotypes
      prop_cols <- grep("^prop_", colnames(plot_comp), value = TRUE)
      prop_labels <- sub("^prop_", "", prop_cols)
      prop_cols <- prop_cols[prop_labels %in% phenotype_partition_labels(prop_labels)]
      if (length(prop_cols) == 0) {
        return(plotly_islet_empty("No phenotype proportion data"))
      }
      if (!is.null(input$pheno_filter) && length(input$pheno_filter) > 0) {
        keep_cols <- paste0("prop_", input$pheno_filter)
        prop_cols <- intersect(prop_cols, keep_cols)
        if (length(prop_cols) == 0) {
          return(plotly_islet_empty("No matching phenotypes selected"))
        }
      }

      # Objects outside the embedding (under ten core cells on a composition root) have no cluster.
      plot_comp <- plot_comp[!is.na(plot_comp[[leiden_col]]), , drop = FALSE]
      if (nrow(plot_comp) == 0) return(plotly_islet_empty("No clustered objects"))
      cluster <- as.character(plot_comp[[leiden_col]])

      # Compute mean phenotype proportions per cluster
      cluster_levels <- leiden_levels(cluster)
      mat <- matrix(NA_real_, nrow = length(prop_cols), ncol = length(cluster_levels),
                    dimnames = list(prop_cols, cluster_levels))
      for (cl in cluster_levels) {
        sub <- plot_comp[cluster == cl, prop_cols, drop = FALSE]
        if (nrow(sub) > 0) {
          mat[, cl] <- colMeans(sub, na.rm = TRUE)
        }
      }

      # Build stacked bar data
      bar_rows <- list()
      for (i in seq_along(prop_cols)) {
        for (j in seq_along(cluster_levels)) {
          bar_rows[[length(bar_rows) + 1]] <- data.frame(
            phenotype = gsub("^prop_", "", prop_cols[i]),
            cluster = cluster_levels[j],
            proportion = mat[i, j],
            stringsAsFactors = FALSE
          )
        }
      }
      bar_df <- do.call(rbind, bar_rows)

      # Remove phenotypes with zero proportion across all clusters
      max_per_pheno <- tapply(bar_df$proportion, bar_df$phenotype, max, na.rm = TRUE)
      keep_phenos <- names(max_per_pheno[max_per_pheno > 0])
      bar_df <- bar_df[bar_df$phenotype %in% keep_phenos, , drop = FALSE]

      # Use active phenotype palette where available
      pheno_present <- unique(bar_df$phenotype)
      pal <- stats::setNames(unname(palette()[pheno_present]), pheno_present)
      pal[is.na(pal)] <- "#CCCCCC"

      canonical_order <- names(palette())
      ordered_phenos <- intersect(canonical_order, keep_phenos)
      ordered_phenos <- c(ordered_phenos, setdiff(keep_phenos, canonical_order))
      # One trace per phenotype, under its group in the legend (owner, 2026-10-27): a click on a group's title shows
      # or hides the whole group.
      ordered_phenos <- ordered_phenos[order(match(phenotype_group(ordered_phenos), PHENOTYPE_GROUP_ORDER))]
      bar_df$cluster <- factor(bar_df$cluster, levels = cluster_levels)

      fonts <- plotly_islet_fonts()
      p <- plot_ly()
      titled <- character(0)
      for (ph in ordered_phenos) {
        d <- bar_df[bar_df$phenotype == ph, , drop = FALSE]
        g <- phenotype_group(ph)
        p <- add_trace(p, x = d$cluster, y = d$proportion, type = "bar", name = phenotype_text(ph),
                       marker = list(color = unname(pal[ph])), legendgroup = g,
                       legendgrouptitle = if (!(g %in% titled)) list(text = g),
                       hovertemplate = paste0(phenotype_text(ph), "<br>Cluster %{x}: %{y:.3f}<extra></extra>"))
        titled <- union(titled, g)
      }
      p %>%
        layout(
          font = fonts$global,
          barmode = "stack",
          title = list(text = "Mean Phenotype Composition by Cluster", font = fonts$title),
          xaxis = list(title = list(text = "Cluster", font = fonts$axis_title),
                       tickfont = fonts$axis_tick),
          yaxis = list(title = list(text = "Mean Proportion", font = fonts$axis_title),
                       tickfont = fonts$axis_tick, range = c(0, 1)),
          # Up to 12 reviewer-facing phenotype entries — use the dense floor directly rather than
          # a max() formula that already evaluates to it.
          legend = list(font = list(size = ISLET_FONT_SIZES$dense_min),
                        title = list(text = "", font = fonts$legend_title),
                        groupclick = "toggleitem", tracegroupgap = 6)
        ) %>%
        plotly_islet_config("cluster_composition")
    }

    # ==== Cluster membership (owner review, 2026-10-20) ====
    # Replaced a mean-composition-by-cluster chart that plotted back the endocrine shares the clusters are built
    # from. Each cluster's objects by donor status (or by donor): which composition states are disease-associated,
    # and whether any belongs to a single donor.
    cluster_membership_data <- reactive({
      req(has_leiden())
      comp <- nbr_comp()
      leiden_col <- input$leiden_res %||% {
        nbr <- get_nbr_columns(colnames(comp))
        nbr$leiden[2] %||% nbr$leiden[1]
      }
      req(leiden_col %in% colnames(comp))
      m <- cluster_membership(comp, leiden_col)
      req(m)
      m$leiden_col <- leiden_col
      m
    })

    output$cluster_composition <- renderPlotly({
      by <- input$cluster_view %||% "composition"
      if (identical(by, "composition")) return(cluster_composition_plot())
      m <- cluster_membership_data()
      fonts <- plotly_islet_fonts()
      donors <- m$donors
      hover_cluster <- paste0("Cluster ", donors$cluster, ": ", donors$objects, " objects, ", donors$donors,
                              " donors; largest donor ", donors$top_donor, " holds ",
                              sprintf("%.0f%%", 100 * donors$top_donor_share))
      p <- plot_ly()
      if (identical(by, "donor")) {
        comp <- nbr_comp()
        cl <- as.character(comp[[m$leiden_col]]); ok <- !is.na(cl)
        d <- as.data.frame(table(cluster = cl[ok], donor = as.character(comp$`Case ID`[ok])), stringsAsFactors = FALSE)
        d <- d[d$Freq > 0, , drop = FALSE]
        d$share <- d$Freq / donors$objects[match(d$cluster, as.character(donors$cluster))]
        status_of <- tapply(as.character(comp$`Donor Status`), as.character(comp$`Case ID`), `[`, 1)
        donor_ids <- names(status_of)[order(match(status_of, STATUS_ORDER), suppressWarnings(as.numeric(names(status_of))))]
        pal <- rep_len(leiden_palette, length(donor_ids))
        for (i in seq_along(donor_ids)) {
          s <- d[d$donor == donor_ids[i], , drop = FALSE]
          if (!nrow(s)) next
          p <- p %>% add_trace(x = s$cluster, y = s$share, type = "bar",
                               name = paste0(donor_ids[i], " (", status_of[donor_ids[i]], ")"),
                               marker = list(color = pal[i]),
                               hovertemplate = paste0("Donor ", donor_ids[i], " (", status_of[donor_ids[i]], ")<br>",
                                                      "Cluster %{x}: %{y:.0%} of its objects<extra></extra>"))
        }
        legend_title <- "Donor"
      } else {
        bs <- m$by_status
        dcols <- donor_colors_reactive()
        for (s in levels(bs$status)) {
          b <- bs[bs$status == s, , drop = FALSE]
          p <- p %>% add_trace(x = as.character(b$cluster), y = b$share, type = "bar", name = s,
                               marker = list(color = dcols[s]),
                               text = hover_cluster[match(as.character(b$cluster), as.character(donors$cluster))],
                               customdata = b$objects,
                               hovertemplate = paste0(s, ": %{customdata} objects (%{y:.0%})<br>%{text}<extra></extra>"),
                               textposition = "none")
        }
        legend_title <- "Status"
      }
      # Each cluster's donor count above its bar.
      p %>% add_annotations(x = as.character(donors$cluster), y = 1, text = paste0(donors$donors, " d"),
                            showarrow = FALSE, yanchor = "bottom", font = list(size = ISLET_FONT_SIZES$dense_min)) %>%
        layout(
          font = fonts$global, barmode = "stack",
          title = list(text = paste0("Cluster membership by ", if (identical(by, "donor")) "donor" else "donor status",
                                     " (Leiden ", sub("^leiden_", "", m$leiden_col), ")"),
                       font = fonts$title),
          margin = list(t = 70),
          xaxis = list(title = list(text = "Cluster", font = fonts$axis_title), tickfont = fonts$axis_tick,
                       type = "category", categoryorder = "array", categoryarray = levels(donors$cluster)),
          yaxis = list(title = list(text = "Proportion of the cluster's objects", font = fonts$axis_title),
                       tickfont = fonts$axis_tick, range = c(0, 1.08), tickformat = ".0%"),
          legend = list(title = list(text = legend_title, font = fonts$legend_title),
                        font = list(size = ISLET_FONT_SIZES$dense_min))
        ) %>%
        plotly_islet_config("cluster_membership")
    })

    # ---- Download handler ----
    #
    # This CSV used to correspond to NO plot on the page. It re-derived its own
    # subset from prepared()$comp instead of the shared `nbr_comp()` reactive
    # every on-screen card uses, so it:
    #   * ignored `nbr_min_cells`   -> exported islets the cards had dropped
    #   * applied `total_cells_peri > 0`, which the cards do NOT -> dropped
    #     islets the cards keep
    #   * ignored `pheno_filter`    -> dumped every phenotype's columns
    #   * carried no outlier flags and no provenance header at all
    # It now consumes nbr_comp() and records what produced it.
    output$download_spatial <- downloadHandler(
      filename = function() {
        paste0("spatial_neighborhood_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        df <- tryCatch(nbr_comp(), error = function(e) NULL)
        sel_phenos <- phenotype_partition_labels(input$pheno_filter)

        nbr <- get_nbr_columns(colnames(df %||% data.frame()))
        metric_cols <- unlist(nbr, use.names = FALSE)

        # Restrict the per-phenotype metric families to the selected phenotypes,
        # matching the tissue scatter and Card B composition. Columns that carry
        # no phenotype (immune_frac_*, leiden_*, totals) always pass through.
        metric_cols <- metric_cols[vapply(metric_cols, function(cl) {
          ph <- extract_phenotype_from_column(cl)
          is.null(ph) || identical(ph, cl) || ph %in% sel_phenos
        }, logical(1))]

        keep_cols <- c("Case ID", "Donor Status", "islet_key", "islet_diam_um",
                       "total_cells_core", "total_cells_peri", metric_cols)
        keep_cols <- intersect(unique(keep_cols), colnames(df %||% data.frame()))
        out <- if (is.null(df)) NULL else phenotype_filter_conflict_columns(df[, keep_cols, drop = FALSE])

        write_islet_csv(file, out, islet_csv_header(
          "Spatial Neighborhood — per-islet metrics",
          paste("Donor groups:", paste(input$groups %||% "(all)", collapse = ", ")),
          paste("Endocrine objects:", paste(sort(unique(as.character(prepared()$comp$object_class %||% "Islet"))), collapse = ", ")),
          paste("Phenotypes:", if (is.null(sel_phenos) || length(sel_phenos) == 0)
                                 "(none selected — no per-phenotype columns)"
                               else paste(phenotype_text(sel_phenos), collapse = ", ")),
          paste("Rows:", if (is.null(out)) 0 else nrow(out), "islets"),
          "Scope: all donors in the selected groups (the tissue scatter above shows one donor at a time)."
        ))
      }
    )

    # ======================================================================
    # NEIGHBORHOOD ANALYSIS CARDS (A: Infiltration, B: Enrichment, C: Proximity)
    # ======================================================================

    # ---- Shared reactive: filtered comp data for neighborhood cards ----
    nbr_comp <- reactive({
      pd <- prepared()
      req(pd$comp)
      comp <- phenotype_filter_conflict_columns(pd$comp)
      if (!is.null(input$groups) && "Donor Status" %in% colnames(comp))
        comp <- comp[comp$`Donor Status` %in% input$groups, , drop = FALSE]
      # No min-cells filter since 2026-10-17: the endocrine-object classes (prepared()) select the objects.
      comp
    })

    # Phenotype options derived dynamically from prepared()$comp columns.
    # Each entry maps a phenotype label to its associated metric columns.
    # Replaces the previous hand-curated 7-immune-type list so users can
    # interrogate any phenotype the pipeline computed metrics for.
    phenotype_options <- reactive({
      comp <- nbr_comp()
      if (is.null(comp)) return(list())
      cn <- colnames(comp)
      ez_cols <- grep("^enrich_z_", cn, value = TRUE)
      if (length(ez_cols) == 0) return(list())
      build <- function(ez) {
        safe <- sub("^enrich_z_", "", ez)
        # Map safe-name back to the original phenotype label by checking
        # possible un-sanitised forms of the name against existing prop_* cols.
        candidates <- unique(c(
          gsub("plus", "+", gsub("_", " ", safe), fixed = TRUE),
          gsub("_", " ", safe),
          gsub("plus", "+", safe, fixed = TRUE),
          safe
        ))
        prop_core_col <- NA_character_
        # Fallback label: the shared un-sanitiser, NOT the raw safe-name. If the
        # prop_* probe below fails (Excel fallback path, a reticulate name
        # mangle, a phenotype with no core column) the label previously degraded
        # to "CD4_Tcell" / "ECADplus" and the tab silently spoke a different
        # dialect from the rest of the app.
        label <- extract_phenotype_from_column(ez) %||% safe
        for (cand in candidates) {
          pc <- paste0("prop_", cand)
          if (pc %in% cn) { prop_core_col <- pc; label <- cand; break }
        }
        list(
          label      = label,
          safe       = safe,
          enrich     = ez,
          prop_core  = prop_core_col,
          peri_prop  = {p <- paste0("peri_prop_", safe);  if (p %in% cn) p else NA_character_},
          peri_count = {p <- paste0("peri_count_", safe); if (p %in% cn) p else NA_character_},
          min_dist   = {p <- paste0("min_dist_", safe);   if (p %in% cn) p else NA_character_}
        )
      }
      recs <- lapply(ez_cols, build)
      # The shared retired-phenotype gate: every neighborhood selector, the
      # composition filter, both composition plots and the CSV exports read
      # phenotype_options(), so filtering here keeps them all on the 12 labels.
      recs <- Filter(function(r) phenotype_is_kept(r$label), recs)
      if (length(recs) == 0) return(list())
      recs <- recs[order(sapply(recs, function(r) r$label))]
      names(recs) <- sapply(recs, function(r) r$label)
      recs
    })

    # Choice vector for selectInputs — names are user-facing labels, values are safe-names.
    phenotype_choices <- reactive({
      opts <- phenotype_options()
      if (length(opts) == 0) return(character(0))
      vec <- vapply(opts, function(o) o$safe, character(1))
      stats::setNames(vec, phenotype_text(vapply(opts, function(o) o$label, character(1))))
    })

    # ---- Card B: donor-level values behind the bar chart and the donor heatmap ----
    # One row per donor and phenotype: the mean over the donor's objects (avoids pseudoreplication; an object-level
    # median is 0 for rare immune types because most objects have none of them). Immune phenotypes only: the endocrine
    # and exocrine classes are left out (owner, 2026-10-20).
    # The band enrichment z uses the cohort's share of each phenotype (enrich_z_*) or, with Baseline "Own tissue",
    # the donor's own tissue share (own_tissue_enrich_z(); owner review, 2026-10-20).
    enrich_donor_values <- reactive({
      comp <- nbr_comp()
      req(comp, "Donor Status" %in% colnames(comp), "Case ID" %in% colnames(comp))
      clip <- isTRUE(input$enrich_clip)
      region <- input$enrich_region %||% "peri"
      own <- identical(region, "peri") && identical(input$enrich_baseline %||% "cohort", "own")
      profile <- if (own) donor_distance_profile() else NULL
      if (own) req(profile)
      sel <- phenotype_partition_labels(input$enrich_phenos)
      rows <- list()
      for (itm in phenotype_options()) {
        if (!(itm$label %in% sel)) next
        if (itm$label %in% NON_IMMUNE_LABELS) next
        col <- switch(region, peri = itm$enrich, core = itm$prop_core, itm$peri_prop)
        if (is.null(col) || is.na(col) || !(col %in% colnames(comp))) next
        vals <- if (own) own_tissue_enrich_z(comp, itm$safe, itm$label, profile) else comp[[col]]
        keep <- is.finite(vals)
        if (region == "peri" && clip) keep <- keep & vals >= -5 & vals <= 5
        if (!any(keep)) next
        m <- tapply(vals[keep], as.character(comp$`Case ID`[keep]), mean)
        st <- tapply(as.character(comp$`Donor Status`[keep]), as.character(comp$`Case ID`[keep]), `[`, 1)
        rows[[length(rows) + 1]] <- data.frame(
          col = if (own) paste0("own_tissue_", col) else col, cell_type = itm$label, case_id = names(m),
          donor_status = as.character(st[names(m)]), value = as.numeric(m), stringsAsFactors = FALSE)
      }
      if (length(rows) == 0) return(NULL)
      do.call(rbind, rows)
    })

    # ---- Card B: summary across donors (median and IQR, or mean and SEM) ----
    enrich_summary <- reactive({
      dv <- enrich_donor_values()
      req(dv)
      stat_mode <- input$enrich_stat %||% "Median"
      groups <- split(dv, list(dv$cell_type, dv$donor_status), drop = TRUE)
      out <- do.call(rbind, lapply(groups, function(g) {
        v <- g$value; n <- length(v)
        if (stat_mode == "Median") {
          q <- stats::quantile(v, c(0.25, 0.75), names = FALSE); s <- stats::median(v); lo <- q[1]; hi <- q[2]
        } else {
          s <- mean(v); se <- if (n > 1) stats::sd(v) / sqrt(n) else 0; lo <- s - se; hi <- s + se
        }
        data.frame(col = g$col[1], cell_type = g$cell_type[1], donor_status = g$donor_status[1],
                   z_summary = s, z_lo = lo, z_hi = hi, n = n, stringsAsFactors = FALSE)
      }))
      # Keep the phenotype order of phenotype_options().
      out[order(match(out$cell_type, names(phenotype_options())), match(out$donor_status, STATUS_ORDER)), ,
          drop = FALSE]
    })

    enrich_measure_label <- function() {
      switch(input$enrich_region %||% "peri",
             peri = if (identical(input$enrich_baseline %||% "cohort", "own")) "Enrichment z (own tissue)"
                    else "Enrichment z (cohort)",
             core = "Proportion (core)", "Proportion (peri-islet)")
    }

    # Card B CSV: the per-donor values and the summary the two plots render.
    output$dl_enrichment_csv <- downloadHandler(
      filename = function() paste0("spatial_enrichment_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        dv <- tryCatch(enrich_donor_values(), error = function(e) NULL)
        es <- tryCatch(enrich_summary(), error = function(e) NULL)
        stat_mode <- input$enrich_stat %||% "Median"
        out <- NULL
        if (!is.null(dv) && nrow(dv) > 0) {
          out <- rbind(
            data.frame(level = "donor", phenotype = dv$cell_type, donor_status = dv$donor_status,
                       case_id = dv$case_id, source_column = dv$col, value = dv$value, lo = NA_real_,
                       hi = NA_real_, n_donors = NA_integer_, stringsAsFactors = FALSE),
            data.frame(level = paste0(tolower(stat_mode), "_of_donor_means"), phenotype = es$cell_type,
                       donor_status = es$donor_status, case_id = NA_character_, source_column = es$col,
                       value = es$z_summary, lo = es$z_lo, hi = es$z_hi, n_donors = es$n,
                       stringsAsFactors = FALSE))
        }
        write_islet_csv(file, out, islet_csv_header(
          "Spatial — phenotype composition & enrichment",
          paste("Measure:", enrich_measure_label()),
          paste("Interval:", if (stat_mode == "Median") "lo/hi = 25th/75th percentile of donor means"
                             else "lo/hi = mean -/+ SEM of donor means"),
          paste("Donor groups:", paste(input$groups %||% "(all)", collapse = ", ")),
          paste("Endocrine objects:", paste(sort(unique(as.character(prepared()$comp$object_class %||% "Islet"))), collapse = ", ")),
          paste("Phenotypes shown:", paste(phenotype_text(intersect(phenotype_partition_labels(input$enrich_phenos),
                                                                  names(phenotype_options()))) %||% "(none)", collapse = ", ")),
          paste("Object z outside [-5, 5] dropped before the donor mean:", if (isTRUE(input$enrich_clip)) "ON" else "OFF"),
          "Rows: level=donor is each donor's mean over its objects; the other level summarizes those donor means."
        ))
      }
    )

    # The card's phenotypes: the immune ones (NON_IMMUNE_LABELS left out; owner, 2026-10-20).
    enrich_labels <- function() sort(phenotype_partition_labels(setdiff(names(phenotype_options()), NON_IMMUNE_LABELS)))

    # ---- Card B phenotype filter (checkbox group shared by bar chart + heatmap) ----
    # Mirrors the tissue-scatter output$phenotype_filter. Choices are the
    # phenotype_options() labels (so they match enrich_summary()'s itm$label).
    # Default hides the same bulk classes as the tissue scatter so the shared
    # proportion scale is not dominated by Acinar's exocrine abundance.
    output$enrich_pheno_filter <- renderUI({
      labels <- enrich_labels()
      if (length(labels) == 0) return(NULL)
      prev <- isolate(input$enrich_phenos)
      selected <- if (!is.null(prev)) remap_macrophage_selection(prev, labels) else labels
      # Grouped, with All | None per group and for the list (owner, 2026-10-27).
      tagList(
        tags$label("Phenotypes shown", style = "font-weight: 600; font-size: 13px; margin: 0 0 2px;"),
        phenotype_group_checkboxes(ns("enrich_phenos"), labels, selected = selected,
                                   marker_input_id = if (root_has_marker_profile_labels()) ns("mac_markers_enrich"),
                                   markers = isolate(macrophage_markers())))
    })

    # ---- Neighborhood cards UI (conditional on data availability) ----
    # Layout (Oct 2026 streamlining): the Neighborhood Analysis control bar sits
    # directly above the three islet-level cards it drives (point size and
    # transparency reach the scatter and the distance box plot; the outlier
    # toggle reaches the enrichment bars and the distance box plot). The
    # composition card and the distance-profile card follow without their own
    # section headings. Each panel and card's explanatory prose lives in its own
    # "?" popover (card_help()), not inline.
    output$neighborhood_cards <- renderUI({
      if (!has_neighborhood()) return(NULL)

      card_style <- "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%; box-sizing: border-box;"
      card_title <- function(txt, ...) div(style = "display: flex; align-items: center; gap: 4px; margin: 0 0 8px 0; align-self: center;",
        h5(txt, style = "font-size: 15px; font-weight: 600; margin: 0;"),
        card_help(paste(txt, "help"), ...))
      ctrl_row   <- "display: flex; gap: 6px 18px; align-items: flex-end; flex-wrap: wrap; margin-bottom: 4px;"
      dl_csv <- function(id) div(style = "margin-top: 6px;",
        downloadButton(ns(id), "Download CSV", style = "font-size: 12px; padding: 3px 10px;"))

      tagList(
        # ==== Section heading ====
        # The section's own Data / View panels were removed with the outlier box plots (owner review, 2026-10-20):
        # every remaining control chooses what one card plots, so it sits on that card. The donor status and
        # endocrine-object selection come from the Data controls at the top of the tab.
        fluidRow(column(12,
          div(style = "display: flex; align-items: baseline; gap: 12px; margin: 18px 0 10px;",
            h4("Neighborhood Analysis", style = "margin: 0; font-weight: 600;"),
            uiOutput(ns("nbr_islet_count"), inline = TRUE))
        )),

        # ==== Row 1: insulitis, insulitis by beta-cell content, immune gradient ====
        fluidRow(style = "display: flex; flex-wrap: wrap;",
          column(4,
            div(class = "card", style = card_style,
              card_title("Insulitis prevalence",
                tags$p("Each donor's proportion of the selected endocrine objects with at least the chosen number of",
                       " immune cells in the core and the 20 µm band. Immune cells are every T Cell, B Cell,",
                       " Macrophage and other leukocyte, subtyped or not. Bars are donor means ± SEM and each",
                       " point is a donor, so every donor counts once; the p value is a Kruskal-Wallis test on the",
                       " donor proportions."),
                tags$p("The threshold follows nPOD's criterion of 15 CD45+ cells around an islet. It is a cell",
                       " count, so small objects rarely reach it: with Clusters selected the proportions fall. Uncheck",
                       " Cluster in the Data controls to count Islets only.")),
              div(style = ctrl_row,
                sliderInput(ns("insulitis_threshold"), "Immune cells per object (at least)", min = 1, max = 50,
                            value = 15, step = 1, width = "100%")),
              plotlyOutput(ns("insulitis_bars"), height = "420px"),
              dl_csv("dl_insulitis_csv")
            )
          ),
          column(4,
            div(class = "card", style = card_style,
              card_title("Immune cells and beta-cell content",
                tags$p("For each donor, the mean immune proportion of the 20 µm band around its Insulin-containing",
                       " Islets and around its Insulin-deficient Islets, joined by a line. An object is insulin-deficient",
                       " when beta cells are under the chosen proportion of its core's Endocrine Beta, Alpha, Delta and CHGA⁺",
                       " hormone-negative cells; the two groups are the endocrine objects of the Data controls. A donor",
                       " with only one kind of object is a single point."),
                tags$p("In T1D, insulitis is expected mostly around islets that still hold beta cells.")),
              div(style = ctrl_row,
                sliderInput(ns("beta_cut"), "Insulin-deficient below this beta-cell proportion", min = 0, max = 0.5,
                            value = 0.1, step = 0.05, width = "100%")),
              plotlyOutput(ns("beta_immune"), height = "420px"),
              dl_csv("dl_beta_immune_csv")
            )
          ),
          column(4,
            div(class = "card", style = card_style,
              card_title("Immune gradient around islets",
                tags$p("For each donor and each distance from the nearest endocrine object's border, the chosen",
                       " cells' proportion of all cells there, divided by the same proportion in the donor's own tissue (its",
                       " cells outside every endocrine object and its band). 1, the dashed line, is the donor's",
                       " tissue level; above 1 the cells gather there. Lines are the median of donors with their",
                       " IQR. Core holds the cells inside objects."),
                tags$p("The dotted line is the donor chosen in the Data controls (the tissue scatter's donor), so",
                       " changing the donor shows where that donor sits against its status. Hover gives its cell",
                       " counts; a point resting on fewer than 20 of the chosen cells is faded, as one donor's B Cell or",
                       " T Cell CD4⁺ cells in a bin are often only a few."),
                tags$p("Dividing by each donor's own tissue separates gathering at islets from a donor having more",
                       " immune cells everywhere. It reads per-cell data, so the Endocrine objects selector does",
                       " not apply.")),
              div(style = ctrl_row,
                selectInput(ns("gradient_phenotype"), "Cells", choices = isolate(gradient_choices()), selected = "all",
                            width = "100%")),
              plotlyOutput(ns("immune_gradient"), height = "480px"),
              dl_csv("dl_gradient_csv")
            )
          )
        ),

        # ==== Row 2: phenotype composition (one card, two plots) ====
        # The phenotype filter and the Region / Summary radios feed
        # enrich_summary(), which BOTH plots consume.
        fluidRow(column(12,
          div(class = "card", style = card_style,
            div(style = ctrl_row,
              card_title("Phenotype composition & enrichment",
                tags$p("Every immune phenotype at once (and, from v15, the cells named only by their top marker), in the core or the 20 \u00b5m band. The left chart summarizes",
                       " each donor status across its donors; the heatmap shows each donor, grouped by status, so the",
                       " spread behind each bar is visible. Each donor's value is the mean over its objects."),
                tags$p(tags$strong("Enrichment z"),
                       " (peri-islet mode) compares each phenotype's count in an object's band with the count expected",
                       " from a reference proportion: z = (observed \u2212 expected) / \u221aexpected. z > 0 = more than",
                       " expected. Baseline Cohort uses the phenotype's proportion of all analyzed cells in the cohort, so",
                       " a donor with more of it everywhere scores high: it measures how many there are near islets.",
                       " Own tissue uses the donor's own proportion outside every endocrine object and its band: it",
                       " measures whether they gather at islets. Proportion modes show the raw fraction (0\u20131)",
                       " of each compartment's cells."),
                tags$p("The endocrine and exocrine classes are left out.",
                       " Drop object z beyond \u00b15 removes those objects' values before each donor mean. The",
                       " phenotype checkboxes and the Region / Summary choices apply to both plots.")),
              radioButtons(ns("enrich_region"), "Region",
                           c("Peri-islet (enrichment z)" = "peri",
                             "Core (proportion)" = "core",
                             "Peri-islet (proportion)" = "peri_prop"),
                           selected = "peri", inline = TRUE),
              conditionalPanel(
                condition = paste0("input['", ns("enrich_region"), "'] == 'peri'"),
                radioButtons(ns("enrich_baseline"), "Baseline",
                             c("Cohort" = "cohort", "Own tissue" = "own"), selected = "cohort", inline = TRUE)
              ),
              radioButtons(ns("enrich_stat"), "Summary", c("Median", "Mean"),
                           selected = "Median", inline = TRUE),
              conditionalPanel(
                condition = paste0("input['", ns("enrich_region"), "'] == 'peri'"),
                checkboxInput(ns("enrich_clip"), "Drop object z beyond \u00b15", value = TRUE)
              )
            ),
            uiOutput(ns("enrich_pheno_filter")),
            fluidRow(
              column(6, plotlyOutput(ns("enrichment_bars"), height = "500px")),
              column(6,
                plotlyOutput(ns("enrichment_heatmap"), height = "500px"),
                dl_csv("dl_enrichment_csv"))
            )
          )
        )),

        # ==== Row 3: tissue-wide distance profile ====
        fluidRow(column(12,
          div(class = "card", style = card_style,
            # The Senior-ported per-cell tissue-wide density / composition. CODEX phenotypes only.
            div(style = ctrl_row,
              card_title("Distance-to-islet profile",
                tags$p("Per-cell signed distance-to-islet on a symlog axis \u2014 a density split by disease status",
                       " (per phenotype or all cells) or a phenotype composition (stacked bars, one facet per status,",
                       " zero-anchored at the islet border), spanning islet interior to bulk tissue. Negative = inside",
                       " islet, positive = outside; the dashed line marks the border. Pooled per cell (descriptive)",
                       " \u2014 every cell counts once, so donors with more cells weigh more. For a donor-fair view see",
                       " Immune gradient around islets."),
                tags$p(tags$strong("Group by"),
                       " Phenotype lists the immune phenotypes (the endocrine and exocrine classes are left out);",
                       " Marker+ lists every panel marker, all 31 being measured in every donor, a cell counting as",
                       if (root_has_background_positivity()) paste0(" positive when its log value exceeds the image's",
                         " own background peak by three background widths; E-cadherin when at the epithelial edge or",
                         " above (the lower edge of the image's E-cadherin-expressing population) and panCK above a",
                         " two-component mixture gate on its log values, both floored at the image's empty-slide null.")
                       else " positive at 0.10 or above on the cohort marker scale."),
                tags$p(tags$strong("Y-axis Proportion."),
                       " In composition bars it is the proportion of all cells at that distance that are of the phenotype",
                       " or positive for the marker, so a flat profile means an even mix and a peak means the cells",
                       " gather there. Each phenotype or marker has its own row and y-axis, so a rare one is not",
                       " dwarfed by a common one; the row's ND, Aab+ and T1D columns use that one scale. In density curves it is the proportion of the status's cells of that class",
                       " per \u00b5m: where those cells lie, which also follows how much tissue lies at each distance."),
                tags$p(tags$strong("Distance to."),
                       " Nearest Islet (the default) measures each cell to the nearest Islet's border; inside an",
                       " Islet it is negative. Nearest Islet or Cluster also counts the Clusters, which are small and",
                       " scattered through the tissue, so most cells beyond 20 \u00b5m are then measured to a Cluster",
                       " and the profile around Islets is blurred."),
                tags$p("Both views run to 200 \u00b5m outside the border. Composition bars are 5 \u00b5m wide;",
                       " density curves use 160 bins on the symlog axis. The dashed line is the border and the",
                       " dotted line the outer edge of the 20 \u00b5m peri-islet band."))
            ),
            tagList(
              div(style = ctrl_row,
                # Distance to the nearest Islet, or to the nearest Islet or Cluster (owner, 2026-10-20). Roots
                # before v9 hold only the second.
                if (cell_distance_has_islet_distance())
                  radioButtons(ns("dist_to"), "Distance to",
                               choices = c("Nearest Islet" = "signed_dist_islet",
                                           "Nearest Islet or Cluster" = "signed_dist"),
                               selected = "signed_dist_islet", inline = TRUE),
                radioButtons(ns("dist_mode"), "Group by",
                             choices = c("Phenotype" = "phenotype", "Marker+" = "marker"),
                             selected = "phenotype", inline = TRUE),
                radioButtons(ns("dist_norm"), "Y-axis",
                             choices = c("Proportion" = "prop", "Raw count" = "count"),
                             selected = "prop", inline = TRUE),
                radioButtons(ns("dist_ptype"), "Plot",
                             choices = c("Density curves" = "density",
                                         "Composition bars" = "composition"),
                             selected = "density", inline = TRUE),
                checkboxGroupInput(ns("dist_status"), "Disease status",
                                   choices = c("ND", "Aab+", "T1D"),
                                   selected = c("ND", "Aab+", "T1D"), inline = TRUE),
                # dist_group / dist_comp_items choices are populated server-side by
                # observeEvent(input$dist_mode) (phenotypes vs markers).
                conditionalPanel(
                  condition = "input.dist_ptype == 'density'", ns = ns,
                  selectInput(ns("dist_group"), "Show", choices = NULL, width = "260px")
                )
              ),
              conditionalPanel(
                condition = "input.dist_ptype == 'composition'", ns = ns,
                # Flow the marker checkboxes into as many ~100px columns as fit, instead of one tall column; the
                # phenotypes are grouped (phenotype_group_checkboxes()). output$dist_comp_items_ui draws either.
                tags$style(HTML(paste0(
                  ".cd-comp-cols,.cd-comp-cols .shiny-input-checkboxgroup,",
                  ".cd-comp-cols .shiny-options-group{width:100%;max-width:none;",
                  "display:block;box-sizing:border-box;}",
                  ".cd-comp-cols .shiny-options-group{column-width:100px;column-gap:16px;}",
                  ".cd-comp-cols .checkbox{break-inside:avoid;margin:2px 0;}"))),
                tags$div(class = "cd-comp-cols", uiOutput(ns("dist_comp_items_ui")))
              )
            ),
            uiOutput(ns("distance_profile_plot"))
          )
        ))
      )
    })

    # ---- Islet count display ----
    output$nbr_islet_count <- renderUI({
      comp <- nbr_comp()
      # Count islets with neighborhood metrics. `comp` also carries QuPath
      # composition rows for islets without a paired core+peri cell census,
      # which no neighborhood card can plot.
      n <- if (is.null(comp)) 0 else if ("total_cells_peri" %in% colnames(comp))
             sum(is.finite(comp$total_cells_peri)) else nrow(comp)
      span(style = "color: var(--islet-fg-muted); font-size: 13px; white-space: nowrap;",
           paste0(formatC(n, big.mark = ","), " objects"))
    })

    # Card C-Right "cellwise": the phenotype/marker pickers follow the group-by mode. The density view's single
    # choice is grouped (option groups); the composition view's checkboxes are grouped with All | None for phenotypes
    # (owner, 2026-10-27) and keep All | None as one list for markers.
    dist_items <- reactive({
      # Phenotype mode lists the immune phenotypes (owner, 2026-10-20: no endocrine or exocrine classes); marker
      # mode every marker the data root flags (all 31 panel markers from v8).
      if (identical(input$dist_mode %||% "phenotype", "marker")) cell_distance_markers()
      else phenotype_partition_labels(macrophage_collapse(setdiff(cell_distance_phenotypes(), NON_IMMUNE_LABELS), macrophage_markers()))
    })
    observeEvent(list(input$dist_mode, macrophage_markers()), {
      mode <- input$dist_mode %||% "phenotype"
      items <- dist_items()
      ch <- if (identical(mode, "marker")) c("All cells" = "__all__", stats::setNames(items, items))
            else c(list("All cells" = "__all__"), phenotype_grouped_choices(items))
      cur_g <- isolate(input$dist_group)
      sel_g <- if (!is.null(cur_g) && cur_g %in% c("__all__", items)) cur_g
               else if (!is.null(cur_g) && is_macrophage_label(cur_g) && any(is_macrophage_label(items)))
                 items[is_macrophage_label(items)][1]
               else "__all__"
      updateSelectInput(session, "dist_group", choices = ch, selected = sel_g)
    }, ignoreNULL = FALSE)
    output$dist_comp_items_ui <- renderUI({
      mode <- input$dist_mode %||% "phenotype"
      items <- dist_items()
      def <- if (identical(mode, "marker")) intersect(c("INS", "GCG", "CD3e", "CD8a", "CD68", "CD31"), items) else items
      cur_c <- isolate(input$dist_comp_items)
      sel_c <- if (!is.null(cur_c) && length(cur_c) > 0) remap_macrophage_selection(cur_c, items) else def
      if (!length(sel_c)) sel_c <- def
      if (identical(mode, "marker")) {
        return(tagList(
          div(style = "display: flex; align-items: center; gap: 6px;",
            tags$label("Show", style = "margin: 0;"),
            tags$span(class = "pheno-select-links",
              actionLink(ns("dist_markers_all"), "All", style = "color: var(--islet-link);"),
              tags$span(class = "pheno-select-sep", "|"),
              actionLink(ns("dist_markers_none"), "None", style = "color: var(--islet-link);"))),
          checkboxGroupInput(ns("dist_comp_items"), NULL, choices = items, selected = sel_c)))
      }
      phenotype_group_checkboxes(ns("dist_comp_items"), items, selected = sel_c, label = "Show",
                                 marker_input_id = if (root_has_marker_profile_labels()) ns("mac_markers_dist"),
                                 markers = isolate(macrophage_markers()))
    })
    observeEvent(input$dist_markers_all, updateCheckboxGroupInput(session, "dist_comp_items", selected = dist_items()))
    observeEvent(input$dist_markers_none, updateCheckboxGroupInput(session, "dist_comp_items", selected = character(0)))

    # ==== Neighborhood row 1 (owner review, 2026-10-20): insulitis, insulitis by beta content, immune gradient ====
    # They replaced Peri enrichment by status (one column of the composition card), Peri vs core proportion (immune
    # cells almost never in the core, endocrine cells almost never in the band) and Distance to nearest cell (taken
    # from the core's center, so it tracked object radius). Helpers: R/spatial_islet_immune.R. Every value is
    # donor-level.
    donor_points_bars <- function(d, y_title, title, filename, fmt = "%.3f", percent = FALSE) {
      dcols <- donor_colors_reactive()
      fonts <- plotly_islet_fonts()
      kw <- donor_kw_p(d$value, d$status)
      p <- plot_ly()
      for (s in levels(d$status)) {
        g <- d[d$status == s, , drop = FALSE]
        if (!nrow(g)) next
        m <- mean(g$value); se <- if (nrow(g) > 1) stats::sd(g$value) / sqrt(nrow(g)) else 0
        p <- p %>%
          add_trace(x = s, y = m, type = "bar", showlegend = FALSE,
                    marker = list(color = paste0(dcols[s], "66"), line = list(color = dcols[s], width = 1.5)),
                    error_y = list(type = "data", array = se, color = dcols[s], thickness = 1.5, width = 6),
                    hovertext = sprintf(paste0("%s<br>Mean: ", fmt, "<br>SEM: ", fmt, "<br>N donors: %d"),
                                        s, m, se, nrow(g)), hoverinfo = "text") %>%
          add_trace(x = s, y = g$value, type = "scatter", mode = "markers", showlegend = FALSE,
                    marker = list(size = 9, opacity = 0.8, color = dcols[s]),
                    hovertext = paste0("Donor ", g$case_id, ": ", sprintf(fmt, g$value), g$detail),
                    hoverinfo = "text")
      }
      p %>% layout(
        font = fonts$global, barmode = "overlay", showlegend = FALSE, margin = list(t = 110),
        title = list(text = paste0(title, "<br><sup>Donor means ± SEM; points are donors",
                                   if (is.finite(kw)) sprintf("<br>Kruskal-Wallis p = %.3g", kw) else "", "</sup>"),
                     font = fonts$title),
        xaxis = list(title = list(text = ""), tickfont = fonts$axis_tick, categoryorder = "array",
                     categoryarray = STATUS_ORDER),
        yaxis = list(title = list(text = y_title, font = fonts$axis_title), tickfont = fonts$axis_tick,
                     tickformat = if (percent) ".0%" else NULL, rangemode = "tozero")
      ) %>% plotly_islet_config(filename)
    }

    insulitis_data <- reactive({
      d <- insulitis_prevalence(nbr_comp(), input$insulitis_threshold %||% 15)
      req(d, nrow(d) > 0)
      d$detail <- paste0(" (", d$insulitic, " of ", d$objects, " objects)")
      d
    })

    output$insulitis_bars <- renderPlotly({
      d <- insulitis_data()
      thr <- input$insulitis_threshold %||% 15
      donor_points_bars(d, paste0("Proportion of objects with ≥", thr, " immune cells"),
                        "Insulitis prevalence", "spatial_insulitis_prevalence", fmt = "%.3f", percent = TRUE)
    })

    output$dl_insulitis_csv <- downloadHandler(
      filename = function() paste0("spatial_insulitis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        d <- tryCatch(insulitis_data(), error = function(e) NULL)
        out <- if (is.null(d)) NULL else
          data.frame(case_id = d$case_id, donor_status = as.character(d$status), objects = d$objects,
                     insulitic_objects = d$insulitic, proportion = d$value, stringsAsFactors = FALSE)
        write_islet_csv(file, out, islet_csv_header(
          "Spatial — insulitis prevalence per donor",
          paste("Insulitic object: at least", input$insulitis_threshold %||% 15,
                "immune cells (every immune label) in the core and the 20 um band"),
          paste("Donor groups:", paste(input$groups %||% "(all)", collapse = ", ")),
          paste("Endocrine objects:", paste(sort(unique(as.character(prepared()$comp$object_class %||% "Islet"))), collapse = ", ")),
          paste("Kruskal-Wallis p (donor proportions):",
                if (is.null(d)) "N/A" else formatC(donor_kw_p(d$value, d$status), format = "g", digits = 4))))
      }
    )

    beta_immune_data <- reactive({
      d <- immune_by_beta_content(nbr_comp(), input$beta_cut %||% 0.1)
      req(d, nrow(d) > 0)
      d
    })

    # One line per donor joining its two kinds of object; donors with only one kind are a single point.
    output$beta_immune <- renderPlotly({
      d <- beta_immune_data()
      dcols <- donor_colors_reactive()
      fonts <- plotly_islet_fonts()
      kinds <- levels(d$kind)
      p <- plot_ly()
      shown <- character(0)
      for (id in unique(d$case_id)) {
        g <- d[d$case_id == id, , drop = FALSE]
        g <- g[order(g$kind), , drop = FALSE]
        s <- as.character(g$status[1])
        p <- p %>% add_trace(x = as.character(g$kind), y = g$value, type = "scatter",
                             mode = if (nrow(g) > 1) "lines+markers" else "markers",
                             line = list(color = dcols[s], width = 1.5),
                             marker = list(color = dcols[s], size = 8, opacity = 0.85),
                             name = s, legendgroup = s, showlegend = !(s %in% shown),
                             hovertext = sprintf("Donor %s (%s)<br>%s: %.3f over %d objects", id, s,
                                                 as.character(g$kind), g$value, g$objects),
                             hoverinfo = "text")
        shown <- union(shown, s)
      }
      both <- table(d$case_id)
      p %>% layout(
        font = fonts$global, margin = list(t = 110, b = 100),
        title = list(text = paste0("Immune cells and beta-cell content<br><sup>Mean band immune proportion per donor<br>",
                                   sum(both == 2), " donors have both kinds of object</sup>"),
                     font = fonts$title),
        # The two group names on two lines each, upright (they are long enough to be drawn at an angle otherwise).
        xaxis = list(title = list(text = ""), tickfont = fonts$axis_tick, categoryorder = "array",
                     categoryarray = kinds, tickvals = kinds, ticktext = sub(" Islets$", "<br>Islets", kinds),
                     tickangle = 0),
        yaxis = list(title = list(text = "Immune proportion of the 20 µm band", font = fonts$axis_title),
                     tickfont = fonts$axis_tick, rangemode = "tozero"),
        legend = list(title = list(text = "Status", font = fonts$legend_title, side = "top"),
                      orientation = "h", x = 0.5, xanchor = "center", y = -0.25, yanchor = "top",
                      font = fonts$legend_text)
      ) %>% plotly_islet_config("spatial_immune_by_beta_content")
    })

    output$dl_beta_immune_csv <- downloadHandler(
      filename = function() paste0("spatial_immune_by_beta_content_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        d <- tryCatch(beta_immune_data(), error = function(e) NULL)
        out <- if (is.null(d)) NULL else
          data.frame(case_id = d$case_id, donor_status = as.character(d$status), object_kind = as.character(d$kind),
                     objects = d$objects, mean_band_immune_proportion = d$value, stringsAsFactors = FALSE)
        write_islet_csv(file, out, islet_csv_header(
          "Spatial — band immune proportion by beta-cell content, per donor",
          paste0("Insulin-deficient Islets: beta cells under ", 100 * (input$beta_cut %||% 0.1),
                 "% of the core's Endocrine Beta, Alpha, Delta and CHGA⁺ hormone-negative cells"),
          paste("Donor groups:", paste(input$groups %||% "(all)", collapse = ", ")),
          paste("Endocrine objects:", paste(sort(unique(as.character(prepared()$comp$object_class %||% "Islet"))), collapse = ", "))))
      }
    )

    GRADIENT_MIN_CELLS <- 20
    # The immune phenotypes in their groups (owner, 2026-10-27), each group's aggregate first where it has two or
    # more phenotypes ("All T cells" counts the T cells without a subtype too). The macrophages are one label on the
    # frozen root and twelve marker-profile labels from v14.
    # The macrophages merged for the checked markers (owner, 2026-10-27).
    gradient_phenotypes <- reactive({
      root <- root_phenotype_labels()
      listed <- c("CD8a Tcell", "CD4 Tcell", "B cell",
                  unique(macrophage_collapse(if (length(root)) intersect(MACROPHAGE_LABELS, root) else character(0),
                                             macrophage_markers())))
      if (!any(is_macrophage_label(listed))) listed <- c(listed, "Macrophage")
      if (root_cd45_is_rule_label()) listed <- c(listed, CD45_RULE_LABEL)     # the CD45 rule's class (v16)
      if ("Hematopoietic CD45+" %in% root) listed <- c(listed, "Hematopoietic CD45+")
      listed <- phenotype_partition_labels(listed)
      c(phenotype_aggregates_for(c(listed, "T cell")), listed)
    })
    gradient_choices <- reactive(c(list("All immune cells" = "all"), phenotype_grouped_choices(gradient_phenotypes())))
    gradient_labels <- reactive(c("all" = "All immune cells",
                                  stats::setNames(phenotype_text(gradient_phenotypes()), gradient_phenotypes())))
    observeEvent(gradient_choices(), {
      cur <- isolate(input$gradient_phenotype) %||% "all"
      vals <- unlist(gradient_choices(), use.names = FALSE)
      sel <- if (cur %in% vals) cur
             else if (is_macrophage_label(cur) && any(is_macrophage_label(vals))) vals[is_macrophage_label(vals)][1]
             else "all"
      updateSelectInput(session, "gradient_phenotype", choices = gradient_choices(), selected = sel)
    }, ignoreInit = TRUE)
    gradient_data <- reactive({
      choice <- input$gradient_phenotype %||% "all"
      if (!(choice %in% unlist(gradient_choices(), use.names = FALSE))) choice <- "all"
      labels <- if (identical(choice, "all")) IMMUNE_LABELS
                else if (is_phenotype_aggregate(choice)) phenotype_aggregate_members(choice, IMMUNE_LABELS)
                else choice
      g <- immune_gradient(donor_distance_profile(), labels)
      req(g, nrow(g) > 0)
      g <- g[as.character(g$status) %in% (input$groups %||% STATUS_ORDER), , drop = FALSE]
      req(nrow(g) > 0)
      g$status <- droplevels(g$status)
      g
    })

    output$immune_gradient <- renderPlotly({
      g <- gradient_data()
      label <- unname(gradient_labels()[input$gradient_phenotype %||% "all"])
      dcols <- donor_colors_reactive()
      fonts <- plotly_islet_fonts()
      p <- plot_ly()
      for (s in levels(g$status)) {
        d <- g[g$status == s, , drop = FALSE]
        q <- do.call(rbind, lapply(split(d$value, d$bin, drop = TRUE), function(v)
          data.frame(med = stats::median(v), lo = stats::quantile(v, 0.25, names = FALSE),
                     hi = stats::quantile(v, 0.75, names = FALSE), n = length(v))))
        q$bin <- rownames(q)
        p <- p %>% add_trace(x = q$bin, y = q$med, type = "scatter", mode = "lines+markers", name = s,
                             line = list(color = dcols[s], width = 2.5), marker = list(color = dcols[s], size = 8),
                             error_y = list(type = "data", symmetric = FALSE, array = q$hi - q$med,
                                            arrayminus = q$med - q$lo, color = dcols[s], thickness = 1),
                             hovertext = sprintf("%s, %s µm<br>Median %.2f× own tissue (IQR %.2f-%.2f), %d donors",
                                                 s, q$bin, q$med, q$lo, q$hi, q$n), hoverinfo = "text")
      }
      # The donor chosen in the Data controls (the tissue scatter's donor), drawn over the status medians. A point
      # resting on fewer than GRADIENT_MIN_CELLS of the chosen cells is faded: one donor's B or CD4 T cells in a bin
      # are often a handful.
      sel <- selected_donor()
      dd <- if (is.null(sel)) g[0, ] else g[g$case_id == as.character(as.integer(sel)), , drop = FALSE]
      if (nrow(dd)) {
        dd <- dd[order(dd$bin), , drop = FALSE]
        s <- as.character(dd$status[1])
        p <- p %>% add_trace(x = as.character(dd$bin), y = dd$value, type = "scatter", mode = "lines+markers",
                             name = paste0("Donor ", dd$case_id[1], " (", s, ")"),
                             line = list(color = islet_ink(), width = 2, dash = "dot"),
                             marker = list(color = dcols[s], size = 11, line = list(color = islet_ink(), width = 1.5),
                                           opacity = ifelse(dd$count >= GRADIENT_MIN_CELLS, 1, 0.3)),
                             hovertext = sprintf("Donor %s (%s), %s \u00b5m<br>%.2f\u00d7 own tissue<br>%d of %s cells%s",
                                                 dd$case_id, s, as.character(dd$bin), dd$value, as.integer(dd$count),
                                                 format(dd$cells, big.mark = ","),
                                                 ifelse(dd$count >= GRADIENT_MIN_CELLS, "",
                                                        paste0(" (under ", GRADIENT_MIN_CELLS, ": faded)"))),
                             hoverinfo = "text")
      }
      p %>% add_trace(x = DISTANCE_BIN_LABELS, y = rep(1, length(DISTANCE_BIN_LABELS)), type = "scatter",
                      mode = "lines", line = list(color = islet_ink(), dash = "dash", width = 1),
                      showlegend = FALSE, hoverinfo = "none") %>%
        layout(
          font = fonts$global, margin = list(t = 110),
          title = list(text = paste0(label, " around endocrine objects<br><sup>Proportion at each distance ÷ the ",
                                     "donor's own tissue proportion<br>Median and IQR of donors</sup>"),
                       font = fonts$title),
          xaxis = list(title = list(text = "Distance from the object's border (µm)", font = fonts$axis_title),
                       tickfont = fonts$axis_tick, tickangle = -45, categoryorder = "array",
                       categoryarray = DISTANCE_BIN_LABELS),
          yaxis = list(title = list(text = "× own tissue proportion", font = fonts$axis_title),
                       tickfont = fonts$axis_tick, rangemode = "tozero"),
          # Below the plot: the donor's entry is long, and a side legend would narrow this card's plot.
          legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.42, yanchor = "top",
                        font = fonts$legend_text)
        ) %>% plotly_islet_config("spatial_immune_gradient")
    })

    output$dl_gradient_csv <- downloadHandler(
      filename = function() paste0("spatial_immune_gradient_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        g <- tryCatch(gradient_data(), error = function(e) NULL)
        out <- if (is.null(g)) NULL else
          data.frame(case_id = g$case_id, donor_status = as.character(g$status), distance_bin_um = as.character(g$bin),
                     cells = g$cells, chosen_cells = g$count, proportion = g$share, donor_tissue_proportion = g$tissue_share,
                     ratio_to_own_tissue = g$value, stringsAsFactors = FALSE)
        write_islet_csv(file, out, islet_csv_header(
          "Spatial — immune gradient around endocrine objects, per donor",
          paste("Cells:", unname(gradient_labels()[input$gradient_phenotype %||% "all"])),
          "Distance: from the nearest endocrine object's border (negative = inside, binned as Core)",
          "Own tissue: the donor's cells outside every endocrine object and its 20 um band",
          paste("Donor groups:", paste(input$groups %||% "(all)", collapse = ", "))))
      }
    )

    # ==== Card B-Left: Enrichment/Proportion Grouped Bar Chart ====
    output$enrichment_bars <- renderPlotly({
      es <- enrich_summary()
      req(es)

      status_order <- c("ND", "Aab+", "T1D")
      es$donor_status <- factor(es$donor_status, levels = intersect(status_order, unique(es$donor_status)))

      dcols <- donor_colors_reactive()

      stat_mode <- input$enrich_stat %||% "Median"
      region <- input$enrich_region %||% "peri"
      error_label <- if (stat_mode == "Median") "IQR" else "SEM"

      # Axis label and title depend on region
      y_label <- enrich_measure_label()
      title_prefix <- paste(stat_mode, y_label)

      p <- plot_ly()
      for (s in levels(es$donor_status)) {
        sub <- es[es$donor_status == s, , drop = FALSE]
        if (nrow(sub) == 0) next
        p <- p %>% add_trace(
          x = phenotype_text(sub$cell_type), y = sub$z_summary, name = s,
          type = "bar",
          marker = list(color = dcols[s]),
          error_y = list(
            type = "data",
            symmetric = FALSE,
            array = sub$z_hi - sub$z_summary,
            arrayminus = sub$z_summary - sub$z_lo,
            color = "#666", thickness = 1
          ),
          hovertemplate = paste0("%{x}<br>", s, "<br>", y_label, " = %{y:.3f}<extra></extra>")
        )
      }

      # Reference line at 0 for enrichment z-scores only
      shapes <- if (region == "peri") {
        list(list(type = "line", x0 = -0.5, x1 = length(unique(es$cell_type)) - 0.5, y0 = 0, y1 = 0,
                  line = list(color = "#999", dash = "dash", width = 1)))
      } else list()

      fonts <- plotly_islet_fonts()
      p %>% layout(
        font = fonts$global,
        barmode = "group",
        title = list(text = paste0(title_prefix, "<br><sup>Donor-level means, ", error_label, " (N=donors)</sup>"),
                     font = fonts$title),
        xaxis = list(title = list(text = "", font = fonts$axis_title),
                     tickfont = fonts$axis_tick, tickangle = -30, automargin = TRUE,
                     categoryorder = "array",
                     categoryarray = phenotype_text(unique(es$cell_type))),
        yaxis = list(title = list(text = y_label, font = fonts$axis_title),
                     tickfont = fonts$axis_tick, automargin = TRUE),
        legend = list(title = list(text = "Status", font = fonts$legend_title),
                      font = fonts$legend_text),
        shapes = shapes
      ) %>%
        plotly_islet_config("spatial_enrichment_bars")
    })

    # ==== Card B-Right: donor heatmap (owner review, 2026-10-20) ====
    # One row per donor, grouped by status, so the spread between donors behind each bar of the left chart is
    # visible. It replaced a status heatmap that drew the same numbers as the bar chart.
    output$enrichment_heatmap <- renderPlotly({
      dv <- enrich_donor_values()
      req(dv)
      region <- input$enrich_region %||% "peri"
      cell_types <- intersect(names(phenotype_options()), unique(dv$cell_type))
      donors <- unique(dv[, c("case_id", "donor_status")])
      donors <- donors[order(match(donors$donor_status, STATUS_ORDER), suppressWarnings(as.numeric(donors$case_id))), ]
      row_lab <- paste0(donors$case_id, " (", donors$donor_status, ")")
      mat <- matrix(NA_real_, nrow = nrow(donors), ncol = length(cell_types), dimnames = list(row_lab, cell_types))
      mat[cbind(match(dv$case_id, donors$case_id), match(dv$cell_type, cell_types))] <- dv$value
      # No values in the cells (owner, 2026-10-27: "no values"); hover gives each one.
      if (region == "peri") {
        z_abs_max <- max(abs(mat), na.rm = TRUE)
        if (!is.finite(z_abs_max) || z_abs_max == 0) z_abs_max <- 1
        colorscale <- list(c(0, "#2166AC"), c(0.5, "#FFFFFF"), c(1, "#B2182B")); zmin <- -z_abs_max; zmax <- z_abs_max
      } else {
        zmax <- max(mat, na.rm = TRUE); if (!is.finite(zmax) || zmax == 0) zmax <- 0.1
        colorscale <- list(c(0, "#FFFFFF"), c(0.5, "#FDB863"), c(1, "#B2182B")); zmin <- 0
      }
      label <- enrich_measure_label()
      # Lines between the status groups.
      breaks <- which(diff(match(donors$donor_status, STATUS_ORDER)) != 0)
      shapes <- lapply(breaks, function(b) list(type = "line", xref = "paper", x0 = 0, x1 = 1, yref = "y",
                                                 y0 = b - 0.5, y1 = b - 0.5, line = list(color = islet_ink(), width = 2)))
      fonts <- plotly_islet_fonts()
      plot_ly(x = phenotype_text(colnames(mat)), y = rownames(mat), z = mat, type = "heatmap", colorscale = colorscale,
              zmin = zmin, zmax = zmax,
              hovertemplate = paste0("%{x}<br>Donor %{y}<br>", label, " = %{z:.4f}<extra></extra>"),
              colorbar = list(title = list(text = if (region == "peri") "z" else "Proportion",
                                           font = fonts$colorbar_title), tickfont = fonts$axis_tick)) %>%
        layout(font = fonts$global, title = list(text = paste0(label, " by donor"), font = fonts$title),
               shapes = shapes,
               xaxis = list(title = list(text = ""), tickfont = fonts$axis_tick, tickangle = -30, automargin = TRUE),
               yaxis = list(title = list(text = ""), tickfont = list(size = ISLET_FONT_SIZES$dense_min),
                            autorange = "reversed", automargin = TRUE, type = "category")) %>%
        plotly_islet_config("spatial_enrichment_donor_heatmap")
    })

    # The profile's container grows with the composition rows (one per phenotype or marker).
    cd_composition_height <- function(n_items) max(520, 150 * n_items + 110)
    output$distance_profile_plot <- renderUI({
      n <- if (identical(input$dist_ptype %||% "density", "composition")) length(input$dist_comp_items) else 0
      plotlyOutput(ns("immune_distance_kde"), height = paste0(if (n) cd_composition_height(n) else 520, "px"))
    })

    # ==== Distance-to-islet profile: per-cell density and composition across the tissue ====
    # (release Parquet/DuckDB; development CSV, symlog axis; ported from Senior.) Its donor-averaged KDE view was
    # removed in the owner review of 2026-10-20: that table holds only core and band cells, so it stopped at 20 um,
    # and the Immune gradient card now shows the donor-fair profile.
    output$immune_distance_kde <- renderPlotly({
      {
        ptype <- input$dist_ptype %||% "density"
        # The distance used: the reviewer's choice when the root has the Islet-only column (v9), else the only one.
        dist_col <- if (cell_distance_has_islet_distance() &&
                        identical(input$dist_to %||% "signed_dist_islet", "signed_dist_islet")) "signed_dist_islet"
                    else "signed_dist"
        dist_label <- if (identical(dist_col, "signed_dist_islet")) "nearest Islet" else "nearest Islet or Cluster"
        mode  <- input$dist_mode %||% "phenotype"
        norm  <- input$dist_norm %||% "prop"
        status_order <- c("ND", "Aab+", "T1D")
        statuses <- intersect(status_order, input$dist_status %||% status_order)
        req(length(statuses) > 0)
        dcols <- donor_colors_reactive()
        fonts <- plotly_islet_fonts()
        item_colors <- if (identical(mode, "marker")) CELL_DIST_MARKER_COLORS else PHENOTYPE_COLORS

        if (identical(ptype, "composition")) {
          # ---- Composition: one row per phenotype or marker, one column per donor status ----
          # Owner, 2026-10-20: markers drawn side by side hid small proportions next to large ones, so each item has
          # its own row and y-axis, shared by the three statuses of that row (the same y scaling for each group).
          items <- intersect(input$dist_comp_items, dist_items())
          if (!identical(mode, "marker")) items <- phenotype_partition_labels(items)
          req(length(items) > 0)
          res <- query_cell_distance_composition(mode, items, statuses, norm, dist_col)
          req(res, nrow(res$comp) > 0)
          comp <- res$comp
          present <- intersect(status_order, unique(comp$donor_status))
          items <- items[items %in% comp$item]
          req(length(items) > 0)
          pcols <- item_colors[items]; pcols[is.na(pcols)] <- "#999999"; names(pcols) <- items
          y_lab <- if (identical(norm, "prop")) "proportion of the cells at that distance" else "cells"
          height <- cd_composition_height(length(items))
          gap <- 30 / height
          panels <- list()
          for (m in items) {
            row_max <- max(comp$val[comp$item == m], na.rm = TRUE)
            y_range <- c(0, if (is.finite(row_max) && row_max > 0) row_max * 1.08 else 1)
            for (s in present) {
              sub <- comp[comp$item == m & comp$donor_status == s, , drop = FALSE]
              sub <- sub[order(sub$binx), , drop = FALSE]
              pp <- plot_ly() %>%
                add_trace(x = sub$x, y = sub$val, type = "bar", marker = list(color = unname(pcols[m])),
                          showlegend = FALSE, hoverinfo = "text",
                          hovertext = paste0(s, " · ", phenotype_text(m), "<br>dist ≈ ", round(sub$x), " µm<br>",
                                             formatC(sub$n, big.mark = ","), " of ",
                                             formatC(sub$total, big.mark = ","), " cells (",
                                             sprintf("%.1f%%", 100 * sub$n / sub$total), ")")) %>%
                add_trace(x = c(0, 0), y = y_range, type = "scatter", mode = "lines",
                          line = list(color = islet_ink(), width = 1.2, dash = "dash"),
                          showlegend = FALSE, hoverinfo = "none") %>%
                # The outer edge of the 20 um peri-islet band.
                add_trace(x = c(CD_BAND_UM, CD_BAND_UM), y = y_range, type = "scatter", mode = "lines",
                          line = list(color = islet_ink(), width = 1.1, dash = "dot"),
                          showlegend = FALSE, hoverinfo = "text",
                          hovertext = paste0(CD_BAND_UM, " \u00b5m: band edge")) %>%
                layout(bargap = 0.1,
                       xaxis = list(title = "", tickfont = fonts$axis_tick, zeroline = FALSE,
                                    range = c(res$lo, CD_MAX_OUTSIDE_UM), dtick = 40),
                       yaxis = list(title = list(text = if (identical(s, present[1])) phenotype_text(m) else "",
                                                 font = fonts$legend_title),
                                    tickfont = list(size = ISLET_FONT_SIZES$dense_min), range = y_range,
                                    tickformat = if (identical(norm, "prop")) ".1%" else NULL))
              panels[[length(panels) + 1L]] <- pp
            }
          }
          fig <- plotly::subplot(panels, nrows = length(items), shareX = TRUE, titleY = TRUE,
                                 margin = c(0.025, 0.01, gap / 2, gap / 2))
          heads <- lapply(seq_along(present), function(j) list(
            text = paste0("<b>", present[j], "</b>"), x = (j - 0.5) / length(present), y = 1, xref = "paper",
            yref = "paper", xanchor = "center", yanchor = "bottom", showarrow = FALSE, font = fonts$legend_title))
          fig %>% layout(
            font = fonts$global, showlegend = FALSE,
            margin = list(t = 40, b = 70, l = 80),
            annotations = c(heads, list(list(
              text = paste0("Signed distance to the ", dist_label, " (µm) — inside ← 0 → outside   ·   y = ", y_lab),
              x = 0.5, y = -48 / height, xref = "paper", yref = "paper",
              xanchor = "center", yanchor = "top", showarrow = FALSE,
              font = fonts$axis_title)))
          ) %>% plotly_islet_config("spatial_cell_distance_composition")

        } else {
          # ---- Density: symlog filled curves per status (pooled per cell) ----
          group <- input$dist_group %||% "__all__"
          choices <- if (identical(mode, "marker")) dist_items() else
            unlist(phenotype_grouped_choices(dist_items()), use.names = FALSE)
          if (!(group %in% c("__all__", choices))) group <- "__all__"
          res <- query_cell_distance_hist(mode, group, statuses, norm, dist_col = dist_col)
          req(res, nrow(res$hist) > 0)
          hist <- res$hist
          present <- intersect(status_order, unique(hist$donor_status))
          y_max <- max(hist$dens)
          type_label <- if (identical(group, "__all__")) "All cells"
                        else if (identical(mode, "marker")) paste0(group, "+") else phenotype_text(group)
          y_title <- if (identical(norm, "prop")) "Fraction of class cells per µm"
                     else "Cells per µm of distance"
          leg_title <- if (identical(norm, "prop")) "Status (proportion of class):"
                       else "Status (pooled counts):"
          p <- plot_ly()
          for (s in present) {
            hs <- hist[hist$donor_status == s, , drop = FALSE]
            hs <- hs[order(hs$x), , drop = FALSE]
            n_cells <- res$n_by_status[[s]]
            p <- p %>% add_trace(
              x = hs$x, y = hs$dens, type = "scatter", mode = "lines",
              name = paste0(s, " (n=", formatC(n_cells, big.mark = ","), ")"),
              line = list(color = dcols[s], width = 2.5),
              fill = "tozeroy", fillcolor = paste0(dcols[s], "22"),
              hoverinfo = "text",
              text = paste0(s, "<br>dist = ", round(hs$real_x, 1), " µm<br>",
                            y_title, " = ", signif(hs$dens, 3)))
          }
          p <- p %>% add_trace(x = c(0, 0), y = c(0, y_max * 1.05),
            type = "scatter", mode = "lines",
            line = list(color = islet_ink(), width = 1.5, dash = "dash"),
            showlegend = FALSE, hoverinfo = "none") %>%
            # The outer edge of the 20 um peri-islet band (the axis is symlog, so place it on that scale).
            add_trace(x = rep(cd_symlog(CD_BAND_UM), 2), y = c(0, y_max * 1.05), type = "scatter", mode = "lines",
                      line = list(color = islet_ink(), width = 1.2, dash = "dot"),
                      showlegend = FALSE, hoverinfo = "text", hovertext = paste0(CD_BAND_UM, " \u00b5m: band edge"))
          ticks <- cd_symlog_ticks(res$lo, res$hi)
          p %>% layout(
            font = fonts$global,
            margin = list(t = 80),
            title = list(text = paste0(type_label, " — distance from islet"),
                         font = fonts$title, y = 0.97, yanchor = "top"),
            xaxis = list(title = list(text = paste0("Signed distance to the ", dist_label,
                                                    " (µm, symlog) — inside ← 0 → outside"),
                                      font = fonts$axis_title),
                         tickmode = "array", tickvals = ticks$vals, ticktext = ticks$text,
                         tickfont = fonts$axis_tick, zeroline = FALSE,
                         # A little past 200 um so the last tick shows.
                         range = c(cd_symlog(res$lo), cd_symlog(CD_MAX_OUTSIDE_UM) + 0.02)),
            yaxis = list(title = list(text = y_title, font = fonts$axis_title),
                         tickfont = fonts$axis_tick),
            legend = list(orientation = "h", x = 0.5, xanchor = "center",
                          y = 1.0, yanchor = "bottom",
                          title = list(text = leg_title, font = fonts$legend_title, side = "top"),
                          font = fonts$legend_text)
          ) %>% plotly_islet_config(paste0("spatial_cell_distance_density_",
                                           gsub("[^A-Za-z0-9_]", "_", type_label)))
        }
      }
    })

  })
}
