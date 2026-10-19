# ---------- Spatial Tab module server ----------
# Exports: spatial_server(id, prepared)
#
# Dependencies:
#   prepared()$comp — must contain peri_prop_*, immune_*, enrich_z_*, leiden_* columns
#   prepared()$neighborhood — raw neighborhood data with leiden_umap_1/2
#   PHENOTYPE_COLORS — from drilldown_helpers.R
#   donor_tissue_available(), get_available_donors(), load_donor_tissue() — from spatial_helpers.R

spatial_server <- function(id, prepared, palette = reactive(PHENOTYPE_COLORS),
                           donor_colors_reactive = reactive(DONOR_COLORS_BRIGHT),
                           remove_outliers = reactive(FALSE),
                           outlier_threshold = reactive(3.0)) {

  moduleServer(id, function(input, output, session) {
    output$embedding_method <- renderText({
      # Since 2026-10-19 (owner) a uniform root serves the endocrine composition of the objects with ten or
      # more core cells (scripts/uniform/composition_embedding.py); older roots serve the protein embedding.
      comp_method <- prepared()$composition_embedding
      if (!is.null(comp_method)) {
        return(sprintf(paste0(
          "Built from each object's endocrine composition: the Beta, Alpha, Delta and CHGA+ hormone-negative",
          " shares of its core endocrine cells, for the %s objects with %s or more core cells. UMAP and Leiden",
          " use one neighbour graph. No donor correction: the shares are counted within each object, so staining",
          " differences between donors do not enter them. Smaller objects have too few cells to place, so they",
          " have no UMAP point or cluster and are drawn grey, as tissue, when the tissue scatter is coloured by",
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
    # Rendered inside the closed help popover, so keep it current while hidden.
    outputOptions(output, "embedding_method", suspendWhenHidden = FALSE)
    ns <- session$ns

    # ---- Outlier-table reactiveVals (parity with Plot/Trajectory) ----
    # Populated inside the Card A and Card C renderPlotly bodies, consumed by
    # the namespaced "Show excluded-outlier table" toggles below each card.
    infiltration_outliers <- reactiveVal(NULL)
    distance_outliers <- reactiveVal(NULL)

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
      load_donor_tissue(donor_id)
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
      phenos <- canonical_phenotype_order(observed)
      counts <- table(factor(observed, levels = phenos))
      list(phenos = phenos, counts = counts)
    })

    output$phenotype_filter <- renderUI({
      info <- pheno_choice_info()
      phenos <- info$phenos
      if (length(phenos) == 0) return(NULL)
      counts <- info$counts

      # Labels carry the donor's count so a zero-cell phenotype is visibly
      # present-but-empty rather than silently absent. Values stay canonical.
      choice_names <- lapply(phenos, function(p) {
        n <- suppressWarnings(as.integer(counts[p]))
        if (!is.finite(n)) n <- 0L
        if (n == 0) {
          tags$span(p, tags$span(" (0)"), style = "color: var(--islet-fg-muted);")
        } else {
          tags$span(p)
        }
      })

      # First render applies the hide-by-default set; afterwards the user's
      # selection is preserved verbatim. The old code did
      # `intersect(prev, phenos)` against a per-donor choice list, which left any
      # newly-appearing phenotype permanently unchecked. With a donor-independent
      # list there is nothing to intersect away.
      prev <- isolate(input$pheno_filter)
      selected_phenos <- if (!is.null(prev)) prev else setdiff(phenos, PHENO_HIDE_DEFAULT)

      tagList(
        div(style = "display: flex; align-items: center; gap: 6px; margin-bottom: 4px;",
          tags$label("Phenotypes", style = "font-weight: 600; font-size: 13px; margin: 0;"),
          actionLink(ns("pheno_all"), "All", style = "font-size: 12px; color: var(--islet-link);"),
          span("|", style = "color: var(--islet-fg-muted);"),
          actionLink(ns("pheno_none"), "None", style = "font-size: 12px; color: var(--islet-link);")
        ),
        div(style = "column-count: 2; column-gap: 8px;",
          checkboxGroupInput(ns("pheno_filter"), NULL,
                             choiceNames = choice_names,
                             choiceValues = as.list(phenos),
                             selected = selected_phenos,
                             inline = FALSE)
        )
      )
    })

    observeEvent(input$pheno_all, {
      updateCheckboxGroupInput(session, "pheno_filter",
                               selected = pheno_choice_info()$phenos)
    })

    observeEvent(input$pheno_none, {
      updateCheckboxGroupInput(session, "pheno_filter", selected = character(0))
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
      selected_phenos <- input$pheno_filter
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
          bg_pal <- palette()[bg_phenos]
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
                                        guide = ggplot2::guide_legend(override.aes = list(size = 4)))
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
          pal <- palette()[legend_phenos]
          pal[is.na(pal)] <- "#CCCCCC"

          p <- p + ggplot2::geom_point(
            data = fg,
            ggplot2::aes(x = X_centroid, y = Y_centroid, color = phenotype),
            size = 0.4, alpha = 0.6,
            inherit.aes = FALSE
          ) +
          ggplot2::scale_color_manual(values = pal, breaks = legend_phenos,
                                        name = "Phenotype", na.value = "#CCCCCC",
                                        guide = ggplot2::guide_legend(override.aes = list(size = 4)))
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
          legend.position = "right",
          legend.key.size = ggplot2::unit(0.7, "cm")
        )
    })

    output$tissue_scatter <- shiny::renderPlot({ spatial_tissue_scatter_plot() }, height = 800, bg = "transparent")
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

    # One islet-level UMAP (points from islet_umap_plotly(), interleaved so no group hides another): the Leiden and
    # Donor Status UMAPs share this layout so they draw at one size --
    # the same marker size, the plot area in the left 78% of the widget and the legend in the rest, wide enough for
    # either legend (a legend that does not fit makes plotly widen the right margin and shrink that plot alone).
    spatial_umap_plotly <- function(plot_df, colour_col, pal, legend_title, title, filename_root) {
      fonts <- plotly_islet_fonts()
      plot_df$.hover <- paste0(legend_title, ": ", plot_df[[colour_col]], "<br>",
                               if (colour_col != "donor_status") paste0("Status: ", plot_df$donor_status, "<br>") else "",
                               plot_df$islet_key)
      islet_umap_plotly(plot_df, "umap1", "umap2", colour_col, pal, ".hover") %>%
        layout(
          font = fonts$global,
          title = list(text = title, font = fonts$title),
          margin = list(l = 30, r = 10, t = 50, b = 30),
          xaxis = list(title = list(text = "UMAP 1", font = fonts$axis_title), domain = c(0, 0.78),
                       zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, ticks = "",
                       showline = FALSE, constrain = "domain"),
          yaxis = list(title = list(text = "UMAP 2", font = fonts$axis_title),
                       zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, ticks = "",
                       showline = FALSE, scaleanchor = "x", scaleratio = 1,
                       constrain = "domain"),
          legend = list(title = list(text = legend_title, font = fonts$legend_title),
                        font = fonts$legend_text, x = 0.8, xanchor = "left", y = 1, yanchor = "top",
                        itemsizing = "constant", tracegroupgap = 2)
        ) %>%
        plotly_islet_config(filename_root)
    }

    output$leiden_umap <- renderPlotly({
      req(has_leiden())
      pd <- prepared()
      comp <- pd$comp
      req("leiden_umap_1" %in% colnames(comp), "leiden_umap_2" %in% colnames(comp))

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

      umap1 <- suppressWarnings(as.numeric(plot_comp$leiden_umap_1))
      umap2 <- suppressWarnings(as.numeric(plot_comp$leiden_umap_2))
      cluster <- as.character(plot_comp[[leiden_col]])

      plot_df <- data.frame(
        umap1 = umap1, umap2 = umap2, cluster = cluster,
        donor_status = if ("Donor Status" %in% colnames(plot_comp)) as.character(plot_comp$`Donor Status`) else "",
        islet_key = if ("islet_key" %in% colnames(plot_comp)) as.character(plot_comp$islet_key) else "",
        stringsAsFactors = FALSE
      )
      plot_df <- plot_df[is.finite(plot_df$umap1) & is.finite(plot_df$umap2), , drop = FALSE]
      if (nrow(plot_df) == 0) return(plotly_islet_empty("No UMAP data"))

      # Sort cluster levels numerically
      cluster_levels <- leiden_levels(plot_df$cluster)
      plot_df$cluster <- factor(plot_df$cluster, levels = cluster_levels)
      n_clusters <- length(cluster_levels)
      pal_colors <- rep_len(leiden_palette, n_clusters)
      pal <- setNames(pal_colors, cluster_levels)

      res_label <- gsub("^leiden_", "", leiden_col)

      spatial_umap_plotly(plot_df, "cluster", pal, "Cluster",
                          paste0("Leiden ", res_label, " (", nrow(plot_df), " objects)"),
                          paste0("leiden_umap_", res_label))
    })

    # ==== Donor Status UMAP: the same layout as the Leiden UMAP (spatial_umap_plotly()) ====
    output$spatial_umap_donor <- renderPlotly({
      req(has_leiden())
      pd <- prepared()
      comp <- pd$comp
      req("leiden_umap_1" %in% colnames(comp), "leiden_umap_2" %in% colnames(comp))

      plot_comp <- comp
      if (!is.null(input$groups) && "Donor Status" %in% colnames(plot_comp)) {
        plot_comp <- plot_comp[plot_comp$`Donor Status` %in% input$groups, , drop = FALSE]
      }

      plot_df <- data.frame(
        umap1 = suppressWarnings(as.numeric(plot_comp$leiden_umap_1)),
        umap2 = suppressWarnings(as.numeric(plot_comp$leiden_umap_2)),
        donor_status = if ("Donor Status" %in% colnames(plot_comp)) as.character(plot_comp$`Donor Status`) else "?",
        islet_key = if ("islet_key" %in% colnames(plot_comp)) as.character(plot_comp$islet_key) else "",
        stringsAsFactors = FALSE
      )
      plot_df <- plot_df[is.finite(plot_df$umap1) & is.finite(plot_df$umap2), , drop = FALSE]
      if (nrow(plot_df) == 0) return(plotly_islet_empty("No UMAP data"))
      plot_df$donor_status <- factor(plot_df$donor_status, levels = c("ND", "Aab+", "T1D"))
      plot_df$donor_status <- droplevels(plot_df$donor_status)

      spatial_umap_plotly(plot_df, "donor_status", donor_colors_reactive(), "Status",
                          paste0("Donor status (", nrow(plot_df), " objects)"), "spatial_umap_donor")
    })

    # ==== Card 3 (bottom): Cluster Composition ====
    output$cluster_composition <- renderPlotly({
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
      prop_cols <- prop_cols[phenotype_is_kept(sub("^prop_", "", prop_cols))]
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
      pal <- palette()[pheno_present]
      pal[is.na(pal)] <- "#CCCCCC"

      canonical_order <- names(palette())
      ordered_phenos <- intersect(canonical_order, keep_phenos)
      ordered_phenos <- c(ordered_phenos, setdiff(keep_phenos, canonical_order))
      bar_df$phenotype <- factor(bar_df$phenotype, levels = rev(ordered_phenos))
      bar_df$cluster <- factor(bar_df$cluster, levels = cluster_levels)

      fonts <- plotly_islet_fonts()
      plot_ly(bar_df, x = ~cluster, y = ~proportion, color = ~phenotype,
              colors = pal,
              type = "bar") %>%
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
                        title = list(text = "", font = fonts$legend_title))
        ) %>%
        plotly_islet_config("cluster_composition")
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
        sel_phenos <- input$pheno_filter

        nbr <- get_nbr_columns(colnames(df %||% data.frame()))
        metric_cols <- unlist(nbr, use.names = FALSE)

        # Restrict the per-phenotype metric families to the selected phenotypes,
        # matching the tissue scatter and Card B composition. Columns that carry
        # no phenotype (immune_frac_*, leiden_*, totals) always pass through.
        if (!is.null(sel_phenos) && length(sel_phenos) > 0) {
          metric_cols <- metric_cols[vapply(metric_cols, function(cl) {
            ph <- extract_phenotype_from_column(cl)
            is.null(ph) || identical(ph, cl) || ph %in% sel_phenos
          }, logical(1))]
        }

        keep_cols <- c("Case ID", "Donor Status", "islet_key", "islet_diam_um",
                       "total_cells_core", "total_cells_peri", metric_cols)
        keep_cols <- intersect(unique(keep_cols), colnames(df %||% data.frame()))
        out <- if (is.null(df)) NULL else df[, keep_cols, drop = FALSE]

        write_islet_csv(file, out, islet_csv_header(
          "Spatial Neighborhood — per-islet metrics",
          paste("Donor groups:", paste(input$groups %||% "(all)", collapse = ", ")),
          paste("Endocrine objects:", paste(sort(unique(as.character(prepared()$comp$object_class %||% "Islet"))), collapse = ", ")),
          paste("Phenotypes:", if (is.null(sel_phenos) || length(sel_phenos) == 0)
                                 "(none selected — no per-phenotype columns)"
                               else paste(sel_phenos, collapse = ", ")),
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
      comp <- pd$comp
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
      stats::setNames(vec, vapply(opts, function(o) o$label, character(1)))
    })

    # ---- Intermediate reactive for Card B: enrichment/proportion summary ----
    enrich_summary <- reactive({
      comp <- nbr_comp()
      req(comp, "Donor Status" %in% colnames(comp), "Case ID" %in% colnames(comp))

      clip <- isTRUE(input$enrich_clip)
      stat_mode <- input$enrich_stat %||% "Median"
      region <- input$enrich_region %||% "peri"

      # Card B phenotype filter — controls which phenotypes appear in BOTH the
      # bar chart and the heatmap (both consume enrich_summary()). The default
      # hidden set (Acinar etc.) is applied via the checkbox's initial `selected`
      # in output$enrich_pheno_filter; here we just honor the current selection.
      # A NULL selection ("None" clicked, or the filter not yet mounted) means
      # "show nothing" — every phenotype is skipped below.
      sel <- input$enrich_phenos

      statuses <- unique(comp$`Donor Status`)
      comp_cols <- colnames(comp)
      rows <- list()

      for (itm in phenotype_options()) {
        if (!(itm$label %in% sel)) next
        # Select column based on region
        if (region == "peri") {
          col <- itm$enrich      # enrichment z-score (peri)
        } else if (region == "core") {
          col <- itm$prop_core   # core proportion
        } else {
          col <- itm$peri_prop   # peri proportion
        }
        if (is.null(col) || is.na(col) || !(col %in% comp_cols)) next

        for (ds in statuses) {
          mask <- comp$`Donor Status` == ds
          vals <- comp[[col]][mask]
          donors <- comp$`Case ID`[mask]
          finite_mask <- is.finite(vals)
          vals <- vals[finite_mask]
          donors <- donors[finite_mask]
          if (region == "peri" && clip) {
            clip_mask <- vals >= -5 & vals <= 5
            vals <- vals[clip_mask]
            donors <- donors[clip_mask]
          }
          if (length(vals) == 0) next

          # Aggregate to donor-level means first (avoids pseudoreplication;
          # islet-level median is 0 for rare immune types because >50% of islets
          # have no cells of that type)
          donor_means <- tapply(vals, donors, mean, na.rm = TRUE)
          donor_means <- donor_means[is.finite(donor_means)]
          n <- length(donor_means)
          if (n == 0) next

          if (stat_mode == "Median") {
            z_summary <- median(donor_means)
            q <- quantile(donor_means, c(0.25, 0.75))
            z_lo <- q[1]
            z_hi <- q[2]
          } else {
            z_summary <- mean(donor_means)
            se <- if (n > 1) sd(donor_means) / sqrt(n) else 0
            z_lo <- z_summary - se
            z_hi <- z_summary + se
          }
          rows[[length(rows) + 1]] <- data.frame(
            col = col, cell_type = itm$label, donor_status = ds,
            z_summary = z_summary, z_lo = z_lo, z_hi = z_hi, n = n,
            stringsAsFactors = FALSE
          )
        }
      }
      if (length(rows) == 0) return(NULL)
      do.call(rbind, rows)
    })

    # Card B CSV: exactly the donor-aggregated matrix both the bar chart and the
    # heatmap render (enrich_summary()), not the pre-aggregation islet rows.
    output$dl_enrichment_csv <- downloadHandler(
      filename = function() paste0("spatial_enrichment_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        es <- tryCatch(enrich_summary(), error = function(e) NULL)
        region <- input$enrich_region %||% "peri"
        stat_mode <- input$enrich_stat %||% "Median"
        out <- NULL
        if (!is.null(es) && nrow(es) > 0) {
          out <- data.frame(
            phenotype    = es$cell_type,
            donor_status = es$donor_status,
            source_column = es$col,
            summary      = es$z_summary,
            lo           = es$z_lo,
            hi           = es$z_hi,
            n_donors     = es$n,
            stringsAsFactors = FALSE)
          # Name the interval columns after the statistic that produced them.
          names(out)[names(out) == "summary"] <-
            if (stat_mode == "Median") "median_of_donor_means" else "mean_of_donor_means"
          names(out)[names(out) == "lo"] <- if (stat_mode == "Median") "q25" else "mean_minus_sem"
          names(out)[names(out) == "hi"] <- if (stat_mode == "Median") "q75" else "mean_plus_sem"
        }
        write_islet_csv(file, out, islet_csv_header(
          "Spatial — phenotype composition & enrichment by donor status",
          paste("Region:", switch(region, peri = "peri-islet enrichment z-score",
                                  core = "core proportion", "peri-islet proportion")),
          paste("Statistic:", stat_mode, "of per-donor means"),
          paste("Donor groups:", paste(input$groups %||% "(all)", collapse = ", ")),
          paste("Endocrine objects:", paste(sort(unique(as.character(prepared()$comp$object_class %||% "Islet"))), collapse = ", ")),
          paste("Phenotypes shown:", paste(input$enrich_phenos %||% "(none)", collapse = ", ")),
          paste("Display clip to z in [-5, 5]:", if (isTRUE(input$enrich_clip)) "ON" else "OFF"),
          "Values are aggregated to per-donor means FIRST, then summarised across donors (avoids pseudoreplication)."
        ))
      }
    )

    # ---- Card B phenotype filter (checkbox group shared by bar chart + heatmap) ----
    # Mirrors the tissue-scatter output$phenotype_filter. Choices are the
    # phenotype_options() labels (so they match enrich_summary()'s itm$label).
    # Default hides the same bulk classes as the tissue scatter so the shared
    # proportion scale is not dominated by Acinar's exocrine abundance.
    output$enrich_pheno_filter <- renderUI({
      labels <- sort(names(phenotype_options()))
      if (length(labels) == 0) return(NULL)
      hide_default <- c("Acinar", "Ductal", "Endocrine", "T cell", "Unknown")
      prev <- isolate(input$enrich_phenos)
      selected <- if (!is.null(prev)) intersect(prev, labels) else setdiff(labels, hide_default)
      tagList(
        div(style = "display: flex; align-items: center; gap: 6px; margin-bottom: 4px;",
          tags$label("Phenotypes shown", style = "font-weight: 600; font-size: 13px; margin: 0;"),
          actionLink(ns("enrich_pheno_all"), "All", style = "font-size: 12px; color: var(--islet-link);"),
          span("|", style = "color: var(--islet-fg-muted);"),
          actionLink(ns("enrich_pheno_none"), "None", style = "font-size: 12px; color: var(--islet-link);")
        ),
        # Flow the checkboxes across the full-width card (same pattern as
        # .cd-comp-cols for the distance-profile composition picker).
        tags$style(HTML(paste0(
          ".enrich-pheno-cols,.enrich-pheno-cols .shiny-input-checkboxgroup,",
          ".enrich-pheno-cols .shiny-options-group{width:100%;max-width:none;",
          "display:block;box-sizing:border-box;}",
          ".enrich-pheno-cols .shiny-options-group{column-width:150px;column-gap:12px;}",
          ".enrich-pheno-cols .checkbox{break-inside:avoid;margin:2px 0;}"))),
        div(class = "enrich-pheno-cols",
          checkboxGroupInput(ns("enrich_phenos"), NULL,
                             choices = labels, selected = selected, inline = FALSE)
        )
      )
    })

    observeEvent(input$enrich_pheno_all, {
      updateCheckboxGroupInput(session, "enrich_phenos",
                               selected = sort(names(phenotype_options())))
    })
    observeEvent(input$enrich_pheno_none, {
      updateCheckboxGroupInput(session, "enrich_phenos", selected = character(0))
    })

    # ---- Neighborhood cards UI (conditional on data availability) ----
    # Layout (Oct 2026 streamlining): the Neighborhood Analysis control bar sits
    # directly above the three islet-level cards it drives (point size and
    # transparency reach the scatter and the distance box plot; the outlier
    # toggle reaches the enrichment bars and the distance box plot). The
    # composition card and the distance-profile card follow without their own
    # section headings. All explanatory prose lives in the tab help popover
    # (`tab_help()` in mod_spatial_ui.R), not inline.
    output$neighborhood_cards <- renderUI({
      if (!has_neighborhood()) return(NULL)

      card_style <- "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%; box-sizing: border-box;"
      card_title <- function(txt) h5(txt, style = "font-size: 15px; font-weight: 600; margin: 0 0 8px 0; align-self: center;")
      ctrl_row   <- "display: flex; gap: 6px 18px; align-items: flex-end; flex-wrap: wrap; margin-bottom: 4px;"
      dl_csv <- function(id) div(style = "margin-top: 6px;",
        downloadButton(ns(id), "Download CSV", style = "font-size: 12px; padding: 3px 10px;"))

      tagList(
        # ==== Section controls ====
        # The "Exclude outliers" checkbox and z slider are the non-namespaced
        # mirrors of the Plot/Trajectory controls — app.R synchronisation
        # observers keep all three tabs in step.
        # The section's controls use the same Data / View panels as every tab (control_panel()).
        fluidRow(column(12,
          h4("Neighborhood Analysis", style = "margin: 18px 0 10px; font-weight: 600;"),
          control_panel("Data controls", open = TRUE,
            fluidRow(
              column(3, uiOutput(ns("nbr_islet_count"))),
              column(3,
                tags$div(title = GLOBAL_OUTLIER_TOOLTIP,
                  checkboxInput("spatial_remove_outliers", "Exclude outliers", value = FALSE))
              ),
              column(3,
                sliderInput("spatial_outlier_z", "Outlier z-threshold",
                            min = 0.5, max = 10, value = 3, step = 0.5)
              )
            )
          ),
          control_panel("View controls",
            fluidRow(
              column(3,
                sliderInput(ns("nbr_pt_size"), "Point size", min = 2, max = 12, value = 5, step = 1)
              ),
              column(3,
                sliderInput(ns("nbr_pt_alpha"), "Point transparency", min = 0, max = 0.95, value = 0.6,
                            step = 0.05)
              )
            )
          )
        )),

        # ==== Row 1: the islet-level cards the control bar drives ====
        fluidRow(style = "display: flex; flex-wrap: wrap;",
          column(4,
            div(class = "card", style = card_style,
              card_title("Peri enrichment by status"),
              div(style = ctrl_row,
                div(style = "flex: 1; min-width: 180px;",
                  nbr_select("infiltration_phenotype", "100%")
                ),
                checkboxInput(ns("infiltration_show_outlier_table"),
                              "Show outlier table", value = FALSE)
              ),
              plotlyOutput(ns("infiltration_bars"), height = "420px"),
              dl_csv("dl_infiltration_csv"),
              uiOutput(ns("infiltration_outlier_info"))
            )
          ),
          column(4,
            div(class = "card", style = card_style,
              card_title("Peri vs core proportion"),
              div(style = ctrl_row,
                div(style = "flex: 1; min-width: 180px;",
                  nbr_select("scatter_phenotype", "100%")
                ),
                checkboxInput(ns("scatter_sqrt"), "Sqrt scale", value = TRUE),
                checkboxInput(ns("scatter_trend"), "Trend lines", value = FALSE)
              ),
              plotlyOutput(ns("infiltration_scatter"), height = "420px")
            )
          ),
          column(4,
            div(class = "card", style = card_style,
              card_title("Distance to nearest cell"),
              div(style = ctrl_row,
                div(style = "flex: 1; min-width: 180px;",
                  # Every retained phenotype with min_dist_*.
                  nbr_select("distance_metric", "100%")
                ),
                checkboxInput(ns("distance_show_outlier_table"),
                              "Show outlier table", value = FALSE)
              ),
              plotlyOutput(ns("distance_boxplot"), height = "420px"),
              dl_csv("dl_distance_csv"),
              uiOutput(ns("distance_outlier_info"))
            )
          )
        ),

        # ==== Row 2: phenotype composition (one card, two plots) ====
        # The phenotype filter and the Region / Summary radios feed
        # enrich_summary(), which BOTH plots consume.
        fluidRow(column(12,
          div(class = "card", style = card_style,
            div(style = ctrl_row,
              card_title("Phenotype composition & enrichment"),
              radioButtons(ns("enrich_region"), "Region",
                           c("Peri-islet (enrichment z)" = "peri",
                             "Core (proportion)" = "core",
                             "Peri-islet (proportion)" = "peri_prop"),
                           selected = "peri", inline = TRUE),
              radioButtons(ns("enrich_stat"), "Summary", c("Median", "Mean"),
                           selected = "Median", inline = TRUE),
              conditionalPanel(
                condition = paste0("input['", ns("enrich_region"), "'] == 'peri'"),
                checkboxInput(ns("enrich_clip"), "Clip extreme z > 5", value = TRUE)
              ),
              checkboxInput(ns("enrich_show_values"), "Show heatmap values", value = TRUE)
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
            # Two views: the donor-averaged KDE (default) or the Senior-ported
            # per-cell tissue-wide density/composition. CODEX phenotypes only.
            div(style = ctrl_row,
              card_title("Distance-to-islet profile"),
              radioButtons(ns("dist_view"), "View",
                           choices = c("Donor-averaged KDE" = "kde",
                                       "Per-cell distance" = "cellwise"),
                           selected = "kde", inline = TRUE)
            ),
            conditionalPanel(
              condition = "input.dist_view == 'kde'", ns = ns,
              div(style = ctrl_row,
                div(style = "min-width: 240px;",
                  nbr_select("kde_immune_type", "260px")
                ),
                radioButtons(ns("kde_metric"), "Y-axis", inline = TRUE, selected = "shape",
                             choices = c("Distribution shape" = "shape",
                                         "Cell abundance" = "abundance"))
              )
            ),
            conditionalPanel(
              condition = "input.dist_view == 'cellwise'", ns = ns,
              div(style = ctrl_row,
                radioButtons(ns("dist_mode"), "Group by",
                             choices = c("Phenotype" = "phenotype", "Marker+" = "marker"),
                             selected = "phenotype", inline = TRUE),
                radioButtons(ns("dist_norm"), "Y-axis",
                             choices = c("Proportion of class" = "prop", "Raw count" = "count"),
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
                # Flow the (up to ~20) checkboxes into as many ~100px columns
                # as fit, instead of one tall column. Targets the options-group
                # so it survives the observeEvent(dist_mode) rebuild.
                tags$style(HTML(paste0(
                  ".cd-comp-cols,.cd-comp-cols .shiny-input-checkboxgroup,",
                  ".cd-comp-cols .shiny-options-group{width:100%;max-width:none;",
                  "display:block;box-sizing:border-box;}",
                  ".cd-comp-cols .shiny-options-group{column-width:100px;column-gap:16px;}",
                  ".cd-comp-cols .checkbox{break-inside:avoid;margin:2px 0;}"))),
                tags$div(class = "cd-comp-cols",
                         checkboxGroupInput(ns("dist_comp_items"), "Show",
                                            choices = NULL))
              )
            ),
            plotlyOutput(ns("immune_distance_kde"), height = "520px")
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

    # ---- Card A: excluded-outlier table ----
    # Mirrors the Plot tab's plot_outlier_info pattern: yellow warning banner
    # with header counts, body explanation, and a scrollable table. Three
    # display states — table, no-outliers note, or filter-OFF flagged note.
    output$infiltration_outlier_info <- renderUI({
      if (!isTRUE(input$infiltration_show_outlier_table)) return(NULL)
      tbl <- infiltration_outliers()
      filter_on <- isTRUE(remove_outliers())

      if (is.null(tbl) || nrow(tbl) == 0) {
        return(tags$div(
          style = "margin-top: 10px; padding: 8px 10px; background-color: var(--islet-info-bg); border: 1px solid var(--islet-info-border); border-radius: 5px; color: var(--islet-info-fg); font-size: 12px;",
          "No outliers detected for the current phenotype."
        ))
      }

      header <- if (filter_on) {
        sprintf("⚠️ %d Outlier%s Excluded (Tukey 1.5×IQR per donor status, zero-aware)",
                nrow(tbl), ifelse(nrow(tbl) > 1, "s", ""))
      } else {
        sprintf("ⓘ %d Outlier%s Flagged (global filter OFF — still in summary)",
                nrow(tbl), ifelse(nrow(tbl) > 1, "s", ""))
      }
      body_text <- if (filter_on) {
        "The following islets were excluded from the donor-level enrichment summary because their value fell beyond the Tukey fences (Q1 − 1.5·IQR, Q3 + 1.5·IQR) of the non-zero distribution for their donor status:"
      } else {
        "These points would be excluded if the global outlier filter were turned on; they currently contribute to the donor-level summary:"
      }

      tags$div(
        style = "margin-top: 10px; padding: 10px; background-color: var(--islet-warn-bg); border: 1px solid var(--islet-warn-border); border-radius: 5px; color: var(--islet-warn-fg);",
        tags$h6(style = "margin-top: 0; color: var(--islet-warn-fg);", header),
        tags$p(style = "color: var(--islet-warn-fg); font-size: 12px; margin-bottom: 10px;", body_text),
        tags$div(style = "max-height: 220px; overflow-y: auto;",
          renderTable(tbl, striped = TRUE, hover = TRUE, bordered = TRUE,
                      spacing = "xs", width = "100%")
        )
      )
    })

    # ---- Card C: excluded-outlier table ----
    output$distance_outlier_info <- renderUI({
      if (!isTRUE(input$distance_show_outlier_table)) return(NULL)
      tbl <- distance_outliers()
      filter_on <- isTRUE(remove_outliers())

      if (is.null(tbl) || nrow(tbl) == 0) {
        return(tags$div(
          style = "margin-top: 10px; padding: 8px 10px; background-color: var(--islet-info-bg); border: 1px solid var(--islet-info-border); border-radius: 5px; color: var(--islet-info-fg); font-size: 12px;",
          "No outliers detected for the current distance metric."
        ))
      }

      thr <- suppressWarnings(as.numeric(outlier_threshold()))
      if (!is.finite(thr)) thr <- 3
      header <- if (filter_on) {
        sprintf("⚠️ %d Outlier%s Excluded (|z| > %g per donor status on log(1+distance))",
                nrow(tbl), ifelse(nrow(tbl) > 1, "s", ""), thr)
      } else {
        sprintf("ⓘ %d Outlier%s Flagged (global filter OFF — still in box plot)",
                nrow(tbl), ifelse(nrow(tbl) > 1, "s", ""))
      }
      body_text <- if (filter_on) {
        sprintf("The following islets were excluded from the box plot because their robust |z| on log(1 + distance) exceeded %g within their donor status (both high- and low-side tails). They still render on the plot as faint grey × markers:", thr)
      } else {
        "These points would be excluded if the global outlier filter were turned on; they currently sit in the box plot's distribution:"
      }

      tags$div(
        style = "margin-top: 10px; padding: 10px; background-color: var(--islet-warn-bg); border: 1px solid var(--islet-warn-border); border-radius: 5px; color: var(--islet-warn-fg);",
        tags$h6(style = "margin-top: 0; color: var(--islet-warn-fg);", header),
        tags$p(style = "color: var(--islet-warn-fg); font-size: 12px; margin-bottom: 10px;", body_text),
        tags$div(style = "max-height: 220px; overflow-y: auto;",
          renderTable(tbl, striped = TRUE, hover = TRUE, bordered = TRUE,
                      spacing = "xs", width = "100%")
        )
      )
    })

    # ---- Populate dynamic phenotype dropdowns once data is available ----
    # Replaces the previous hand-curated 7-immune-type lists. Choices are
    # built from prepared()$comp columns (enrich_z_*, prop_*, peri_prop_*,
    # min_dist_*) so every phenotype the pipeline computed metrics for is
    # surfaced uniformly. Defaults pick a sensible immune phenotype when
    # available (Macrophage), else first.
    pick_default <- function(choices, prefer = c("Macrophage", "Immune", "T_cell", "CD8a_Tcell")) {
      if (length(choices) == 0) return(character(0))
      for (p in prefer) if (p %in% choices) return(p)
      choices[[1]]
    }

    # The choices are built into the selectInputs by output$neighborhood_cards itself (2026-10-19). They used to be
    # sent by an observer through updateSelectInput, but outputs on a hidden tab are suspended while observers are
    # not: the observer ran at startup, before the Spatial tab had mounted the inputs, its update was lost, and the
    # four dropdowns opened empty.
    nbr_select_choices <- reactive({
      ch <- phenotype_choices()
      opts <- phenotype_options()
      # Card A-Right: peri vs core proportion (only phenotypes with both prop_ + peri_prop_)
      paired <- vapply(opts, function(o) !is.na(o$prop_core) && !is.na(o$peri_prop), logical(1))
      # Card C-Left: distance metric — min_dist_<phenotype> for retained labels.
      # The "Immune (all)" aggregate was removed: it pools the retired residual
      # labels (Immune, T cell, APCs) with the retained immune subtypes.
      dist_opts <- Filter(function(o) !is.na(o$min_dist), opts)
      dist_ch <- stats::setNames(
                   vapply(dist_opts, function(o) o$min_dist, character(1)),
                   vapply(dist_opts, function(o) o$label, character(1)))
      # Card C-Right "kde" view: the same universe, limited to phenotypes in the KDE data. The "cellwise" view
      # sets its own choices at UI-build time from cell_distance_phenotypes().
      ch_kde <- character(0)
      kde_df <- load_immune_kde_data()
      if (!is.null(kde_df) && "phenotype" %in% colnames(kde_df)) {
        present <- unique(as.character(kde_df$phenotype))
        kopts <- Filter(function(o) o$safe %in% present, opts)
        ch_kde <- stats::setNames(
                    vapply(kopts, function(o) o$safe,  character(1)),
                    vapply(kopts, function(o) o$label, character(1)))
      }
      list(infiltration_phenotype = ch, scatter_phenotype = ch[paired], distance_metric = dist_ch,
           kde_immune_type = ch_kde)
    })
    # Keep the current choice when it is still offered, else a sensible immune default.
    nbr_select_selected <- function(id) {
      choices <- nbr_select_choices()[[id]]
      current <- isolate(input[[id]])
      if (!is.null(current) && current %in% choices) return(current)
      if (identical(id, "distance_metric"))
        return(pick_default(choices, prefer = paste0("min_dist_", c("Macrophage", "CD8a_Tcell"))))
      pick_default(choices)
    }
    nbr_select <- function(id, width) {
      selectInput(ns(id), "Phenotype", choices = nbr_select_choices()[[id]],
                  selected = nbr_select_selected(id), width = width)
    }

    # Card C-Right "cellwise": repopulate the phenotype/marker pickers when the
    # group-by mode flips. Depends on input$dist_mode so it fires post-mount.
    observeEvent(list(input$dist_view, input$dist_mode), {
      if (!identical(input$dist_view %||% "kde", "cellwise")) return()
      mode <- input$dist_mode %||% "phenotype"
      items <- if (identical(mode, "marker")) cell_distance_markers() else cell_distance_phenotypes()
      ch <- c("All cells" = "__all__", stats::setNames(items, items))
      cur_g <- isolate(input$dist_group)
      sel_g <- if (!is.null(cur_g) && cur_g %in% ch) cur_g else "__all__"
      updateSelectInput(session, "dist_group", choices = ch, selected = sel_g)
      def <- if (identical(mode, "marker"))
               intersect(c("INS", "GCG", "CD3e", "CD8a", "CD68", "CD31"), items)
             else
               intersect(phenotype_keep(c("Beta cell", "Alpha cell", "Ductal",
                                         "Macrophage", "CD8a Tcell", "B cell")), items)
      cur_c <- isolate(input$dist_comp_items)
      sel_c <- if (!is.null(cur_c) && length(cur_c) > 0 && all(cur_c %in% items)) cur_c else def
      updateCheckboxGroupInput(session, "dist_comp_items", choices = items, selected = sel_c)
    }, ignoreNULL = FALSE)

    # Selected safe-name or min_dist column -> display label. Marker-rule hints
    # are not shown inline; the definitions live in the Spatial help popover.
    .lookup_pheno_label <- function(safe_or_metric) {
      if (is.null(safe_or_metric) || !nzchar(safe_or_metric)) return(NULL)
      opts <- phenotype_options()
      for (o in opts) if (identical(o$safe, safe_or_metric)) return(o$label)
      for (o in opts) if (identical(o$min_dist, safe_or_metric)) return(o$label)
      safe_or_metric
    }

    # ==== Card A-Left: Per-phenotype peri-zone enrichment z-score (donor-level) ====
    # Card A data, extracted from the renderer so `dl_infiltration_csv` exports
    # exactly what is plotted (same outlier flags, same donor aggregation)
    # rather than re-deriving a parallel version that can drift.
    infiltration_data <- reactive({
      comp <- nbr_comp()
      req(comp, "Donor Status" %in% colnames(comp), "Case ID" %in% colnames(comp))
      opts <- phenotype_options()
      sel_safe <- input$infiltration_phenotype
      req(sel_safe)
      itm <- Filter(function(o) o$safe == sel_safe, opts)
      req(length(itm) > 0)
      itm <- itm[[1]]
      metric <- itm$enrich
      req(metric %in% colnames(comp))
      metric_label <- paste0("Peri Enrichment z-score — ", itm$label)

      df <- data.frame(
        value = comp[[metric]],
        status = comp$`Donor Status`,
        case_id = comp$`Case ID`,
        islet_key = if ("islet_key" %in% colnames(comp)) comp$islet_key else NA_character_,
        stringsAsFactors = FALSE
      )
      df <- df[is.finite(df$value), , drop = FALSE]
      if (nrow(df) == 0) return(list(empty = "No data", label = metric_label, itm = itm))

      # Outlier handling via shared helper. Tukey IQR (not the modified z) because
      # enrichment z is already a standardised, roughly symmetric quantity.
      #
      # This previously passed `exclude_zeros = TRUE` on the belief that
      # enrichment z-scores are zero-inflated. They are not: enrich_z is a Poisson
      # z, (observed - expected) / sqrt(expected) (compute_neighborhood_metrics.py),
      # so a phenotype absent from an islet's peri zone yields z = -sqrt(expected),
      # a negative continuum — never exactly 0. Measured across all 5,214 islets,
      # every enrich_z_* column contains zero exact zeros, so the flag was a no-op.
      # (The real zero-inflation lives in peri_prop_*, which this card doesn't plot.)
      df <- flag_outliers(df, value_col = "value", group_col = "status",
                          method = "iqr", threshold = 1.5,
                          enabled = isTRUE(remove_outliers()))
      filter_summary <- summarize_outlier_filter(df)
      df_full <- df

      flagged <- df_full[isTRUE(attr(df_full, "outlier_enabled")) & df_full$is_outlier, , drop = FALSE]
      if (!isTRUE(attr(df_full, "outlier_enabled"))) {
        flagged <- df_full[df_full$is_outlier, , drop = FALSE]
      }

      if (isTRUE(remove_outliers())) {
        df <- df[!df$is_outlier | is.na(df$is_outlier), , drop = FALSE]
      }
      if (nrow(df) == 0) {
        return(list(empty = "No data after outlier removal", label = metric_label,
                    itm = itm, islets = df_full, flagged = flagged))
      }

      # Aggregate to donor-level means
      donor_df <- df %>%
        dplyr::group_by(case_id, status) %>%
        dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

      status_order <- c("ND", "Aab+", "T1D")
      donor_df$status <- factor(donor_df$status, levels = intersect(status_order, unique(donor_df$status)))

      # Donor-level summary: mean +/- SEM per group
      summary_df <- donor_df %>%
        dplyr::group_by(status) %>%
        dplyr::summarise(
          mean_val = mean(value, na.rm = TRUE),
          sd_val = sd(value, na.rm = TRUE),
          n = dplyr::n(),
          sem = sd(value, na.rm = TRUE) / sqrt(dplyr::n()),
          .groups = "drop"
        )

      # Kruskal-Wallis on donor-level means (N=15)
      kw_p <- tryCatch({
        kt <- kruskal.test(value ~ status, data = donor_df)
        kt$p.value
      }, error = function(e) NA_real_)

      list(itm = itm, label = metric_label, islets = df_full, kept = df,
           flagged = flagged, donor_df = donor_df, summary_df = summary_df,
           kw_p = kw_p, filter_summary = filter_summary)
    })

    output$infiltration_bars <- renderPlotly({
      d <- infiltration_data()
      req(d)

      # Populate the per-card outlier table reactiveVal (consumed by the
      # `infiltration_outlier_info` renderUI below).
      if (!is.null(d$flagged) && nrow(d$flagged) > 0) {
        infiltration_outliers(data.frame(
          Phenotype          = d$itm$label,
          Case_ID            = d$flagged$case_id,
          Donor_Status       = d$flagged$status,
          Value              = round(d$flagged$value, 3),
          Fence_Distance_IQR = round(d$flagged$outlier_score, 2),
          stringsAsFactors   = FALSE
        ))
      } else {
        infiltration_outliers(NULL)
      }
      if (!is.null(d$empty)) return(plotly_islet_empty(d$empty))

      itm <- d$itm; metric_label <- d$label
      donor_df <- d$donor_df; summary_df <- d$summary_df
      filter_summary <- d$filter_summary
      dcols <- donor_colors_reactive()
      p_label <- if (is.finite(d$kw_p))
        paste0("KW p = ", formatC(d$kw_p, format = "g", digits = 3)) else ""

      p <- plot_ly()
      for (s in levels(summary_df$status)) {
        row <- summary_df[summary_df$status == s, , drop = FALSE]
        if (nrow(row) == 0) next
        # Individual donor points
        donor_pts <- donor_df[donor_df$status == s, , drop = FALSE]
        p <- p %>% add_trace(
          x = s, y = donor_pts$value,
          type = "scatter", mode = "markers",
          marker = list(size = 8, opacity = 0.6, color = dcols[s]),
          name = s, showlegend = FALSE,
          hoverinfo = "text",
          text = paste0("Donor ", donor_pts$case_id, ": ", round(donor_pts$value, 4))
        )
        # Bar with SEM error bars
        p <- p %>% add_trace(
          x = s, y = row$mean_val,
          type = "bar",
          marker = list(color = paste0(dcols[s], "66"),
                        line = list(color = dcols[s], width = 1.5)),
          error_y = list(type = "data", array = row$sem, visible = TRUE,
                         color = dcols[s], thickness = 1.5, width = 6),
          name = s, showlegend = FALSE,
          hoverinfo = "text",
          text = paste0(s, "<br>Mean: ", round(row$mean_val, 4),
                       "<br>SEM: ", round(row$sem, 4),
                       "<br>N donors: ", row$n)
        )
      }
      excl_text <- if (filter_summary$enabled && filter_summary$n_excluded > 0) {
        sprintf(" | %d excluded (%s)", filter_summary$n_excluded, filter_summary$label)
      } else if (!filter_summary$enabled) {
        " | outlier filter OFF"
      } else {
        ""
      }
      fonts <- plotly_islet_fonts()
      p %>% layout(
        font = fonts$global,
        title = list(text = paste0(metric_label,
                                   "<br><sup>Donor-level means \u00b1 SEM (N=donors)<br>",
                                   sub("^ \\| ", "", paste0(p_label, excl_text)), "</sup>"),
                     font = fonts$title),
        margin = list(t = 110),
        xaxis = list(title = list(text = "", font = fonts$axis_title),
                     tickfont = fonts$axis_tick,
                     categoryorder = "array",
                     categoryarray = c("ND", "Aab+", "T1D")),
        yaxis = list(title = list(text = metric_label, font = fonts$axis_title),
                     tickfont = fonts$axis_tick),
        barmode = "overlay",
        showlegend = FALSE
      ) %>%
        plotly_islet_config(paste0("spatial_infiltration_", itm$safe))
    })

    # Card A CSV: the per-islet values behind the bars, carrying the outlier
    # flags, plus the donor-level means the bars actually plot.
    output$dl_infiltration_csv <- downloadHandler(
      filename = function() paste0("spatial_infiltration_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        d <- tryCatch(infiltration_data(), error = function(e) NULL)
        out <- NULL
        if (!is.null(d) && !is.null(d$islets)) {
          isl <- d$islets
          out <- data.frame(
            level        = "islet",
            phenotype    = d$itm$label,
            case_id      = isl$case_id,
            donor_status = isl$status,
            islet_key    = isl$islet_key,
            enrich_z     = isl$value,
            is_outlier   = isl$is_outlier,
            fence_dist_iqr = isl$outlier_score,
            stringsAsFactors = FALSE)
          if (!is.null(d$donor_df) && nrow(d$donor_df) > 0) {
            out <- rbind(out, data.frame(
              level        = "donor_mean",
              phenotype    = d$itm$label,
              case_id      = d$donor_df$case_id,
              donor_status = as.character(d$donor_df$status),
              islet_key    = NA_character_,
              enrich_z     = d$donor_df$value,
              is_outlier   = NA,
              fence_dist_iqr = NA_real_,
              stringsAsFactors = FALSE))
          }
        }
        fs <- if (is.null(d)) NULL else d$filter_summary
        write_islet_csv(file, out, islet_csv_header(
          "Spatial — peri enrichment by donor status",
          paste("Phenotype:", if (is.null(d)) "(none)" else d$itm$label),
          paste("Donor groups:", paste(input$groups %||% "(all)", collapse = ", ")),
          paste("Endocrine objects:", paste(sort(unique(as.character(prepared()$comp$object_class %||% "Islet"))), collapse = ", ")),
          paste("Outlier rule:", if (is.null(fs)) "n/a"
                                 else if (isTRUE(fs$enabled))
                                   paste0(fs$label, "; ", fs$n_excluded, " of ", fs$n_total, " excluded")
                                 else "OFF"),
          paste("Kruskal-Wallis p (donor-level means):",
                if (is.null(d) || !is.finite(d$kw_p)) "N/A"
                else formatC(d$kw_p, format = "g", digits = 4)),
          "Rows: level=islet are the per-islet values (is_outlier marks the ghost x markers); level=donor_mean are the values the bars aggregate."
        ))
      }
    )

    # ==== Card A-Right: Peri vs Core Proportion Scatter ====
    output$infiltration_scatter <- renderPlotly({
      comp <- nbr_comp()
      req(comp, "Donor Status" %in% colnames(comp))

      opts <- phenotype_options()
      sel_safe <- input$scatter_phenotype
      req(sel_safe)
      itm <- Filter(function(o) o$safe == sel_safe, opts)
      req(length(itm) > 0)
      itm <- itm[[1]]
      cols <- list(peri = itm$peri_prop, core = itm$prop_core,
                   label = paste0(itm$label, " Proportion"))
      req(!is.na(cols$peri), !is.na(cols$core))
      req(cols$peri %in% colnames(comp), cols$core %in% colnames(comp))

      use_sqrt <- isTRUE(input$scatter_sqrt)

      df <- data.frame(
        peri = comp[[cols$peri]],
        core = comp[[cols$core]],
        status = comp$`Donor Status`,
        islet_key = if ("islet_key" %in% colnames(comp)) comp$islet_key else "",
        stringsAsFactors = FALSE
      )
      df <- df[is.finite(df$peri) & is.finite(df$core), , drop = FALSE]
      if (nrow(df) == 0) return(plotly_islet_empty("No data"))

      # Apply sqrt transform if requested
      if (use_sqrt) {
        df$peri_plot <- sqrt(df$peri)
        df$core_plot <- sqrt(df$core)
      } else {
        df$peri_plot <- df$peri
        df$core_plot <- df$core
      }

      status_order <- c("ND", "Aab+", "T1D")
      df$status <- factor(df$status, levels = intersect(status_order, unique(df$status)))
      dcols <- donor_colors_reactive()

      max_val <- max(c(df$peri_plot, df$core_plot), na.rm = TRUE) * 1.05

      pt_sz <- input$nbr_pt_size %||% 5
      pt_al <- 1 - (input$nbr_pt_alpha %||% 0.6)

      p <- plot_ly()
      for (s in levels(df$status)) {
        sub <- df[df$status == s, , drop = FALSE]
        if (nrow(sub) == 0) next
        p <- p %>% add_trace(
          data = sub, x = ~peri_plot, y = ~core_plot,
          text = ~paste0(islet_key, "<br>Status: ", status,
                        "<br>Peri: ", round(peri, 4),
                        "<br>Core: ", round(core, 4)),
          hoverinfo = "text", name = s,
          type = "scatter", mode = "markers",
          marker = list(size = pt_sz, opacity = pt_al, color = dcols[s])
        )
      }

      # Optional linear trend lines per status group
      if (isTRUE(input$scatter_trend)) {
        for (s in levels(df$status)) {
          sub <- df[df$status == s, , drop = FALSE]
          if (nrow(sub) < 3) next
          fit <- lm(core_plot ~ peri_plot, data = sub)
          pred_x <- seq(min(sub$peri_plot, na.rm = TRUE),
                        max(sub$peri_plot, na.rm = TRUE), length.out = 50)
          pred_y <- predict(fit, newdata = data.frame(peri_plot = pred_x))
          p <- p %>% add_trace(
            x = pred_x, y = pred_y,
            type = "scatter", mode = "lines",
            line = list(color = dcols[s], width = 2),
            name = paste0(s, " trend"),
            showlegend = FALSE,
            hoverinfo = "none"
          )
        }
      }

      # Axis labels
      scale_note <- if (use_sqrt) "\u221a" else ""
      x_title <- paste0(scale_note, cols$label, " (Peri)")
      y_title <- paste0(scale_note, cols$label, " (Core)")

      # Diagonal y=x reference
      fonts <- plotly_islet_fonts()
      p %>% add_trace(
        x = c(0, max_val), y = c(0, max_val),
        type = "scatter", mode = "lines",
        line = list(dash = "dash", color = "#999", width = 1),
        showlegend = FALSE, hoverinfo = "none"
      ) %>% layout(
        font = fonts$global,
        title = list(text = paste0("Peri vs Core: ", cols$label), font = fonts$title),
        xaxis = list(title = list(text = x_title, font = fonts$axis_title),
                     tickfont = fonts$axis_tick, range = c(0, max_val)),
        yaxis = list(title = list(text = y_title, font = fonts$axis_title),
                     tickfont = fonts$axis_tick, range = c(0, max_val)),
        legend = list(title = list(text = "Status", font = fonts$legend_title),
                      font = fonts$legend_text)
      ) %>%
        plotly_islet_config(paste0("spatial_infiltration_scatter_", itm$safe))
    })

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
      if (region == "peri") {
        y_label <- "Enrichment z-score"
        title_prefix <- paste0(stat_mode, " Enrichment z-score")
      } else if (region == "core") {
        y_label <- "Proportion (core)"
        title_prefix <- paste0(stat_mode, " Core Proportion")
      } else {
        y_label <- "Proportion (peri-islet)"
        title_prefix <- paste0(stat_mode, " Peri-islet Proportion")
      }

      p <- plot_ly()
      for (s in levels(es$donor_status)) {
        sub <- es[es$donor_status == s, , drop = FALSE]
        if (nrow(sub) == 0) next
        p <- p %>% add_trace(
          x = sub$cell_type, y = sub$z_summary, name = s,
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
                     categoryarray = unique(es$cell_type)),
        yaxis = list(title = list(text = y_label, font = fonts$axis_title),
                     tickfont = fonts$axis_tick, automargin = TRUE),
        legend = list(title = list(text = "Status", font = fonts$legend_title),
                      font = fonts$legend_text),
        shapes = shapes
      ) %>%
        plotly_islet_config("spatial_enrichment_bars")
    })

    # ==== Card B-Right: Enrichment/Proportion Heatmap ====
    output$enrichment_heatmap <- renderPlotly({
      es <- enrich_summary()
      req(es)

      region <- input$enrich_region %||% "peri"
      # Per-cell numeric overlay toggle (default on). NULL-safe so the overlay
      # is present on the first render before the checkbox registers.
      show_vals <- isTRUE(input$enrich_show_values %||% TRUE)

      # Pivot to matrix: cell_type x donor_status
      status_order <- c("ND", "Aab+", "T1D")
      cell_types <- unique(es$cell_type)
      statuses <- intersect(status_order, unique(es$donor_status))

      mat <- matrix(NA_real_, nrow = length(statuses), ncol = length(cell_types),
                    dimnames = list(statuses, cell_types))
      for (i in seq_len(nrow(es))) {
        s <- es$donor_status[i]
        ct <- es$cell_type[i]
        if (as.character(s) %in% rownames(mat) && ct %in% colnames(mat))
          mat[as.character(s), ct] <- es$z_summary[i]
      }

      # Annotation text: more decimals for proportions (small numbers)
      fmt <- if (region == "peri") "%.2f" else "%.4f"
      text_mat <- matrix(sprintf(fmt, mat), nrow = nrow(mat), ncol = ncol(mat))
      text_mat[is.na(mat)] <- ""

      if (region == "peri") {
        # Diverging for z-scores (symmetric around 0)
        z_abs_max <- max(abs(mat), na.rm = TRUE)
        if (!is.finite(z_abs_max) || z_abs_max == 0) z_abs_max <- 1
        colorscale <- list(c(0, "#2166AC"), c(0.5, "#FFFFFF"), c(1, "#B2182B"))
        zmin <- -z_abs_max; zmax <- z_abs_max
        cb_title <- "z-score"
        title_text <- paste0(input$enrich_stat %||% "Median", " Enrichment z-score")
      } else {
        # Sequential for proportions (0 to max)
        zmax <- max(mat, na.rm = TRUE)
        if (!is.finite(zmax) || zmax == 0) zmax <- 0.1
        colorscale <- list(c(0, "#FFFFFF"), c(0.5, "#FDB863"), c(1, "#B2182B"))
        zmin <- 0
        cb_title <- "Proportion"
        region_label <- if (region == "core") "Core" else "Peri-islet"
        title_text <- paste0(input$enrich_stat %||% "Median", " ", region_label, " Proportion")
      }

      fonts <- plotly_islet_fonts()
      plot_ly(
        x = colnames(mat), y = rownames(mat), z = mat,
        type = "heatmap",
        colorscale = colorscale,
        zmin = zmin, zmax = zmax,
        text = if (show_vals) text_mat else NULL,
        texttemplate = if (show_vals) "%{text}" else NULL,
        textfont = if (show_vals) list(size = ISLET_FONT_SIZES$dense_min) else NULL,
        hovertemplate = paste0("%{x}<br>%{y}<br>", cb_title, " = %{z:.4f}<extra></extra>"),
        showscale = TRUE,
        colorbar = list(title = list(text = cb_title, font = fonts$colorbar_title),
                        tickfont = fonts$axis_tick)
      ) %>% layout(
        font = fonts$global,
        title = list(text = title_text, font = fonts$title),
        # automargin: phenotype names are long and angled at -30 deg, so without
        # it plotly clips them against a fixed margin — which reads as "the font
        # is too small" when the real problem is the label is cut off.
        xaxis = list(title = list(text = "", font = fonts$axis_title),
                     tickfont = fonts$axis_tick, tickangle = -30, automargin = TRUE),
        yaxis = list(title = list(text = "", font = fonts$axis_title),
                     tickfont = fonts$axis_tick, autorange = "reversed",
                     automargin = TRUE)
      ) %>%
        plotly_islet_config("spatial_enrichment_heatmap")
    })

    # ==== Card C-Left: Distance Box Plots ====
    # Card C-Left data, extracted so `dl_distance_csv` exports the same per-islet
    # rows and outlier flags the box plot draws.
    distance_data <- reactive({
      comp <- nbr_comp()
      req(comp, "Donor Status" %in% colnames(comp))
      metric <- input$distance_metric
      req(metric, metric %in% colnames(comp))
      opts <- phenotype_options()
      match <- Filter(function(o) !is.na(o$min_dist) && o$min_dist == metric, opts)
      req(length(match) > 0)  # retained phenotypes only
      metric_label <- paste0("Min Distance to ", match[[1]]$label)

      df <- data.frame(
        value = comp[[metric]],
        status = comp$`Donor Status`,
        case_id = if ("Case ID" %in% colnames(comp)) comp$`Case ID` else NA,
        islet_key = if ("islet_key" %in% colnames(comp)) comp$islet_key else NA_character_,
        stringsAsFactors = FALSE
      )

      # Count non-NA per group before filtering
      status_order <- c("ND", "Aab+", "T1D")
      n_total <- tapply(df$value, df$status, length)
      n_valid <- tapply(df$value, df$status, function(x) sum(is.finite(x)))
      na_rate <- 1 - sum(is.finite(df$value)) / nrow(df)

      df <- df[is.finite(df$value), , drop = FALSE]
      if (nrow(df) == 0) {
        return(list(empty = paste0(metric_label, " \u2014 all NA"),
                    label = metric_label, metric = metric))
      }

      # Outlier handling via shared helper, applied to log1p(distance) to
      # tame the long right tail; ghost \u00d7 markers preserve visibility.
      df <- flag_outliers(df, value_col = "value", group_col = "status",
                          method = "zscore", threshold = outlier_threshold(),
                          transform = "log1p",
                          enabled = isTRUE(remove_outliers()))
      filter_summary <- summarize_outlier_filter(df)
      filter_summary$label <- if (filter_summary$enabled) {
        sprintf("|z| > %g per donor status on log(1+distance)",
                suppressWarnings(as.numeric(outlier_threshold())))
      } else filter_summary$label

      list(label = metric_label, metric = metric, islets = df,
           filter_summary = filter_summary, status_order = status_order,
           n_total = n_total, n_valid = n_valid, na_rate = na_rate)
    })

    output$distance_boxplot <- renderPlotly({
      d <- distance_data()
      req(d)
      metric_label <- d$label

      # Populate the per-card outlier table reactiveVal (consumed by the
      # `distance_outlier_info` renderUI below).
      flagged <- if (is.null(d$islets)) NULL else d$islets[d$islets$is_outlier, , drop = FALSE]
      if (!is.null(flagged) && nrow(flagged) > 0) {
        distance_outliers(data.frame(
          Phenotype     = metric_label,
          Case_ID       = flagged$case_id,
          Islet         = flagged$islet_key,
          Donor_Status  = flagged$status,
          Distance_um   = round(flagged$value, 1),
          Z_Score_log1p = round(flagged$outlier_score, 2),
          stringsAsFactors = FALSE
        ))
      } else {
        distance_outliers(NULL)
      }
      if (!is.null(d$empty)) return(plotly_islet_empty(d$empty))

      df <- d$islets; filter_summary <- d$filter_summary
      status_order <- d$status_order
      n_total <- d$n_total; n_valid <- d$n_valid; na_rate <- d$na_rate

      df_ghost <- df[df$is_outlier, , drop = FALSE]
      df_box   <- if (isTRUE(remove_outliers())) df[!df$is_outlier, , drop = FALSE] else df

      df_box$status <- factor(df_box$status, levels = intersect(status_order, unique(df_box$status)))
      if (nrow(df_ghost) > 0) {
        df_ghost$status <- factor(df_ghost$status, levels = levels(df_box$status))
      }
      dcols <- donor_colors_reactive()
      pt_sz <- input$nbr_pt_size %||% 5
      pt_al <- 1 - (input$nbr_pt_alpha %||% 0.6)

      # Subtitle with N per group
      n_labels <- sapply(levels(df_box$status), function(s) {
        nv <- n_valid[s]
        nt <- n_total[s]
        if (is.na(nv)) nv <- 0
        if (is.na(nt)) nt <- 0
        paste0(s, ": ", nv, "/", nt)
      })
      subtitle <- paste0("Non-NA: ", paste(n_labels, collapse = ", "),
                          "<br>", round(na_rate * 100, 0), "% NA overall")
      excl_suffix <- if (filter_summary$enabled && filter_summary$n_excluded > 0) {
        sprintf(" | %d excluded (%s)", filter_summary$n_excluded, filter_summary$label)
      } else if (!filter_summary$enabled) {
        " | outlier filter OFF"
      } else {
        ""
      }

      p <- plot_ly()
      for (s in levels(df_box$status)) {
        sub <- df_box[df_box$status == s, , drop = FALSE]
        if (nrow(sub) == 0) next
        p <- p %>% add_trace(
          y = sub$value, x = s, name = s,
          type = "box",
          boxpoints = "all", jitter = 0.3, pointpos = -1.5,
          marker = list(size = pt_sz, opacity = pt_al, color = dcols[s]),
          line = list(color = dcols[s]),
          fillcolor = paste0(dcols[s], "44"),
          hoverinfo = "y"
        )
      }
      # Ghost layer: excluded outliers as faint grey \u00d7 markers per group
      if (nrow(df_ghost) > 0) {
        for (s in levels(df_box$status)) {
          gsub_ <- df_ghost[df_ghost$status == s, , drop = FALSE]
          if (nrow(gsub_) == 0) next
          p <- p %>% add_trace(
            y = gsub_$value, x = rep(s, nrow(gsub_)),
            type = "scatter", mode = "markers",
            marker = list(symbol = "x", size = pt_sz * 1.2,
                          color = "rgba(80,80,80,0.55)",
                          line = list(color = "rgba(50,50,50,0.7)", width = 0.8)),
            showlegend = FALSE,
            hoverinfo = "text",
            text = paste0("Excluded outlier \u00b7 z = ",
                          sprintf("%.2f", gsub_$outlier_score),
                          " \u00b7 distance = ", round(gsub_$value, 1), " \u00b5m"),
            name = paste0(s, " (excluded)")
          )
        }
      }
      fonts <- plotly_islet_fonts()
      p %>% layout(
        font = fonts$global,
        title = list(text = paste0(metric_label, "<br><sup>", subtitle, excl_suffix, "</sup>"),
                     font = fonts$title),
        margin = list(t = 110),
        xaxis = list(title = list(text = "", font = fonts$axis_title),
                     tickfont = fonts$axis_tick, categoryorder = "array",
                     categoryarray = c("ND", "Aab+", "T1D")),
        yaxis = list(title = list(text = paste0("Distance (\u00b5m)"), font = fonts$axis_title),
                     tickfont = fonts$axis_tick),
        showlegend = FALSE
      ) %>%
        plotly_islet_config(paste0("spatial_distance_box_", gsub("[^A-Za-z0-9_]", "_", d$metric)))
    })

    # Card C-Left CSV: per-islet distances with the log1p outlier flags.
    output$dl_distance_csv <- downloadHandler(
      filename = function() paste0("spatial_distance_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        d <- tryCatch(distance_data(), error = function(e) NULL)
        out <- NULL
        if (!is.null(d) && !is.null(d$islets)) {
          isl <- d$islets
          out <- data.frame(
            metric        = d$label,
            case_id       = isl$case_id,
            islet_key     = isl$islet_key,
            donor_status  = isl$status,
            distance_um   = isl$value,
            is_outlier    = isl$is_outlier,
            z_log1p       = isl$outlier_score,
            stringsAsFactors = FALSE)
        }
        fs <- if (is.null(d)) NULL else d$filter_summary
        write_islet_csv(file, out, islet_csv_header(
          "Spatial — distance to nearest cell of type, per islet",
          paste("Metric:", if (is.null(d)) "(none)" else d$label),
          paste("Donor groups:", paste(input$groups %||% "(all)", collapse = ", ")),
          paste("Endocrine objects:", paste(sort(unique(as.character(prepared()$comp$object_class %||% "Islet"))), collapse = ", ")),
          paste("Outlier rule:", if (is.null(fs)) "n/a"
                                 else if (isTRUE(fs$enabled))
                                   paste0(fs$label, "; ", fs$n_excluded, " of ", fs$n_total, " excluded")
                                 else "OFF"),
          "is_outlier marks the ghost x markers; rows are NOT dropped from this file so the full distribution is exportable."
        ))
      }
    )

    # ==== Card C-Right: two views selected by input$dist_view ====
    #  "kde"      = donor-averaged signed-distance KDE over islet-associated
    #               cells (release Parquet; development CSV) — the DEFAULT.
    #  "cellwise" = per-cell tissue-wide distance density + phenotype composition
    #               (release Parquet/DuckDB; development CSV, symlog axis) — ported from Senior.
    load_immune_kde_data <- function() {
      load_signed_distance_kde()
    }

    output$immune_distance_kde <- renderPlotly({
      view <- input$dist_view %||% "kde"

      if (identical(view, "kde")) {
        kde_df <- load_immune_kde_data()
        req(kde_df)

        # Safety net: a stale CSV without case_id collapses to one donor per
        # status (degrades to pooled, never crashes). The regenerated CSV has it.
        if (!"case_id" %in% colnames(kde_df)) kde_df$case_id <- kde_df$donor_status

        cell_type <- input$kde_immune_type
        req(cell_type, phenotype_is_kept(.lookup_pheno_label(cell_type)))
        metric <- input$kde_metric %||% "shape"
        df <- kde_df[kde_df$phenotype == cell_type, , drop = FALSE]
        type_label <- .lookup_pheno_label(cell_type)

        # Filter to selected donor status groups (sidebar)
        if (!is.null(input$groups)) {
          df <- df[df$donor_status %in% input$groups, , drop = FALSE]
        }
        req(nrow(df) > 0)

        # Clip to 1st-99th percentile
        q_lo <- quantile(df$signed_dist, 0.01, na.rm = TRUE)
        q_hi <- quantile(df$signed_dist, 0.99, na.rm = TRUE)
        df <- df[df$signed_dist >= q_lo & df$signed_dist <= q_hi, , drop = FALSE]
        req(nrow(df) > 0)

        status_order <- c("ND", "Aab+", "T1D")
        df$donor_status <- factor(df$donor_status, levels = intersect(status_order, unique(df$donor_status)))
        dcols <- donor_colors_reactive()

        # Per-donor densities on a shared grid, averaged with EQUAL donor weight
        # (removes pseudoreplication, e.g. ND donor 6356 = 71% of ND's APC cells).
        #   shape     -> MEAN of per-donor densities that each integrate to 1
        #   abundance -> MEDIAN of per-donor densities scaled by that donor's count
        grid_x <- seq(q_lo, q_hi, length.out = 512)
        bw <- tryCatch(stats::bw.nrd0(df$signed_dist[is.finite(df$signed_dist)]),
                       error = function(e) NULL)
        dens_one <- function(v) {
          v <- v[is.finite(v)]
          if (length(v) < 2) return(NULL)
          d <- if (is.null(bw)) density(v, n = 512, from = q_lo, to = q_hi)
               else            density(v, bw = bw, n = 512, from = q_lo, to = q_hi)
          if (identical(metric, "abundance")) d$y * length(v) else d$y
        }
        agg_donors <- function(mat) {
          if (identical(metric, "abundance")) apply(mat, 1, median) else rowMeans(mat)
        }

        kde_list <- list(); n_donors <- list(); n_cells <- list(); y_max <- 0
        for (s in levels(df$donor_status)) {
          sdf <- df[df$donor_status == s, , drop = FALSE]
          cols <- list()
          for (dn in unique(sdf$case_id)) {
            y <- dens_one(sdf$signed_dist[sdf$case_id == dn])
            if (!is.null(y)) cols[[length(cols) + 1L]] <- y
          }
          if (length(cols) < 1) next
          ybar <- agg_donors(do.call(cbind, cols))
          kde_list[[s]] <- ybar
          n_donors[[s]] <- length(cols)
          n_cells[[s]]  <- nrow(sdf)
          y_max <- max(y_max, max(ybar))
        }
        req(length(kde_list) > 0)

        y_title <- if (identical(metric, "abundance"))
                     "Median cells per µm per donor" else "Per-donor mean density (area = 1)"

        fonts <- plotly_islet_fonts()
        p <- plot_ly()
        for (s in names(kde_list)) {
          ybar <- kde_list[[s]]
          lbl <- paste0(s, " (", n_donors[[s]], " donors, ",
                        formatC(n_cells[[s]], big.mark = ","), " cells)")
          p <- p %>% add_trace(
            x = grid_x, y = ybar, type = "scatter", mode = "lines", name = lbl,
            line = list(color = dcols[s], width = 2.5),
            fill = "tozeroy", fillcolor = paste0(dcols[s], "22"),
            hoverinfo = "text",
            text = paste0(lbl, "<br>dist = ", round(grid_x, 1), " µm<br>",
                          y_title, " = ", signif(ybar, 3)))
        }
        p <- p %>% add_trace(x = c(0, 0), y = c(0, y_max * 1.05),
          type = "scatter", mode = "lines",
          line = list(color = islet_ink(), width = 1.5, dash = "dash"),
          showlegend = FALSE, hoverinfo = "none")
        p %>% layout(
          font = fonts$global, margin = list(t = 96),
          title = list(text = paste0(type_label, " — Distance from Islet"),
                       font = fonts$title, y = 0.97, yanchor = "top"),
          xaxis = list(title = list(text = "Signed distance (µm) — inside islet ← 0 → outside",
                                    font = fonts$axis_title),
                       tickfont = fonts$axis_tick, zeroline = FALSE),
          yaxis = list(title = list(text = y_title, font = fonts$axis_title),
                       tickfont = fonts$axis_tick),
          legend = list(orientation = "h", x = 0.5, xanchor = "center",
                        y = 1.0, yanchor = "bottom",
                        title = list(text = "Status (donor-averaged):",
                                     font = fonts$legend_title, side = "top"),
                        font = fonts$legend_text)
        ) %>%
          plotly_islet_config(paste0("spatial_distance_kde_", gsub("[^A-Za-z0-9_]", "_", type_label)))

      } else {
        ptype <- input$dist_ptype %||% "density"
        mode  <- input$dist_mode %||% "phenotype"
        norm  <- input$dist_norm %||% "prop"
        status_order <- c("ND", "Aab+", "T1D")
        statuses <- intersect(status_order, input$dist_status %||% status_order)
        req(length(statuses) > 0)
        dcols <- donor_colors_reactive()
        fonts <- plotly_islet_fonts()
        item_colors <- if (identical(mode, "marker")) CELL_DIST_MARKER_COLORS else PHENOTYPE_COLORS

        if (identical(ptype, "composition")) {
          # ---- Composition: stacked bars (phenotype or marker), one facet per status ----
          items <- input$dist_comp_items
          req(length(items) > 0)
          res <- query_cell_distance_composition(mode, items, statuses, norm)
          req(res, nrow(res$comp) > 0)
          comp <- res$comp
          present <- intersect(status_order, unique(comp$donor_status))
          pcols <- item_colors[items]; pcols[is.na(pcols)] <- "#999999"; names(pcols) <- items
          y_lab <- if (identical(norm, "prop")) "fraction of class cells" else "cells"
          panels <- lapply(present, function(s) {
            ss <- comp[comp$donor_status == s, , drop = FALSE]
            pp <- plot_ly()
            for (m in items) {
              sub <- ss[ss$item == m, , drop = FALSE]
              sub <- sub[order(sub$binx), , drop = FALSE]
              if (nrow(sub) == 0) next
              pp <- pp %>% add_trace(
                x = sub$x, y = sub$val, type = "bar",
                name = m, legendgroup = m, showlegend = identical(s, present[1]),
                marker = list(color = unname(pcols[m])),
                hoverinfo = "text",
                hovertext = paste0(s, " · ", m, "<br>dist ≈ ", round(sub$x), " µm<br>",
                                   formatC(sub$n, big.mark = ","), " cells"))
            }
            tot <- stats::aggregate(val ~ binx, ss, sum)
            ymax <- if (nrow(tot)) max(tot$val) else 1
            pp <- pp %>% add_trace(x = c(0, 0), y = c(0, ymax * 1.02),
              type = "scatter", mode = "lines",
              line = list(color = islet_ink(), width = 1.4, dash = "dash"),
              showlegend = FALSE, hoverinfo = "none")
            pp %>% layout(barmode = "stack",
              xaxis = list(title = "", tickfont = fonts$axis_tick, zeroline = FALSE),
              yaxis = list(title = list(text = s, font = fonts$legend_title),
                           tickfont = fonts$axis_tick))
          })
          fig <- if (length(panels) == 1) panels[[1]]
                 else plotly::subplot(panels, nrows = length(panels), shareX = TRUE,
                                      titleY = TRUE, margin = 0.04)
          fig %>% layout(
            font = fonts$global,
            margin = list(t = 58, b = 62),
            legend = list(orientation = "h", x = 0.5, xanchor = "center",
                          y = 1.0, yanchor = "bottom", font = fonts$legend_text),
            annotations = list(list(
              text = paste0("Signed distance (µm) — inside ← 0 → outside   ·   y = ", y_lab),
              x = 0.5, y = -0.09, xref = "paper", yref = "paper",
              xanchor = "center", yanchor = "top", showarrow = FALSE,
              font = fonts$axis_title))
          ) %>% plotly_islet_config("spatial_cell_distance_composition")

        } else {
          # ---- Density: symlog filled curves per status (pooled per cell) ----
          group <- input$dist_group %||% "__all__"
          res <- query_cell_distance_hist(mode, group, statuses, norm)
          req(res, nrow(res$hist) > 0)
          hist <- res$hist
          present <- intersect(status_order, unique(hist$donor_status))
          y_max <- max(hist$dens)
          type_label <- if (identical(group, "__all__")) "All cells"
                        else if (identical(mode, "marker")) paste0(group, "+") else group
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
            showlegend = FALSE, hoverinfo = "none")
          ticks <- cd_symlog_ticks(res$lo, res$hi)
          p %>% layout(
            font = fonts$global,
            margin = list(t = 80),
            title = list(text = paste0(type_label, " — distance from islet"),
                         font = fonts$title, y = 0.97, yanchor = "top"),
            xaxis = list(title = list(text = "Signed distance (µm, symlog) — inside ← 0 → outside",
                                      font = fonts$axis_title),
                         tickmode = "array", tickvals = ticks$vals, ticktext = ticks$text,
                         tickfont = fonts$axis_tick, zeroline = FALSE),
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
