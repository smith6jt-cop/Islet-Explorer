# ---------- Spatial Tab module server ----------
# Exports: spatial_server(id, prepared)
#
# Dependencies:
#   prepared()$comp — must contain peri_prop_*, immune_*, enrich_z_*, leiden_* columns
#   prepared()$neighborhood — raw neighborhood data with leiden_umap_1/2
#   PHENOTYPE_COLORS — from drilldown_helpers.R
#   donor_tissue_available(), get_available_donors(), load_donor_tissue() — from spatial_helpers.R

spatial_server <- function(id, prepared, palette = reactive(PHENOTYPE_COLORS),
                           donor_colors_reactive = reactive(DONOR_COLORS)) {

  moduleServer(id, function(input, output, session) {
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
        return(div(style = "color: #888; font-style: italic; padding-top: 8px;",
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
      req(input$donor)
      load_donor_tissue(input$donor)
    })

    # ---- Donor status for selected donor ----
    donor_status <- reactive({
      req(input$donor)
      pd <- prepared()
      if (!is.null(pd$comp) && "Donor Status" %in% colnames(pd$comp) &&
          "Case ID" %in% colnames(pd$comp)) {
        ds <- pd$comp$`Donor Status`[pd$comp$`Case ID` == as.integer(input$donor)]
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
      donor_id <- input$donor
      req(donor_id)

      comp <- pd$comp
      donor_mask <- comp$`Case ID` == as.integer(donor_id)
      if (!any(donor_mask)) return(NULL)

      sub <- comp[donor_mask, c("islet_key", leiden_col), drop = FALSE]
      setNames(as.character(sub[[leiden_col]]), as.character(sub$islet_key))
    })

    # ---- Phenotype filter ----
    output$phenotype_filter <- renderUI({
      cells <- tryCatch(donor_cells(), error = function(e) NULL)
      if (is.null(cells)) return(NULL)
      phenos <- sort(unique(cells$phenotype))
      if (length(phenos) == 0) return(NULL)
      hide_default <- c("Acinar", "Ductal", "Endocrine", "T cell", "Unknown")
      selected_phenos <- setdiff(phenos, hide_default)
      tagList(
        div(style = "display: flex; align-items: center; gap: 6px; margin-bottom: 4px;",
          tags$label("Phenotypes", style = "font-weight: 600; font-size: 13px; margin: 0;"),
          actionLink(ns("pheno_all"), "All", style = "font-size: 11px;"),
          span("|", style = "color: #ccc;"),
          actionLink(ns("pheno_none"), "None", style = "font-size: 11px;")
        ),
        div(style = "column-count: 2; column-gap: 8px;",
          checkboxGroupInput(ns("pheno_filter"), NULL,
                             choices = phenos, selected = selected_phenos,
                             inline = FALSE)
        )
      )
    })

    observeEvent(input$pheno_all, {
      cells <- tryCatch(donor_cells(), error = function(e) NULL)
      if (!is.null(cells)) {
        phenos <- sort(unique(cells$phenotype))
        updateCheckboxGroupInput(session, "pheno_filter", selected = phenos)
      }
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
    output$tissue_scatter <- renderPlot({
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
          cluster_levels <- sort(unique(fg$cluster[fg$cluster != "tissue"]))
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
          # Phenotype coloring — legend strictly limited to checked phenotypes
          if (!is.null(selected_phenos) && length(selected_phenos) > 0) {
            legend_phenos <- sort(intersect(selected_phenos,
                                            unique(c(fg$phenotype, if (color_bg) bg$phenotype))))
          } else {
            legend_phenos <- sort(unique(c(fg$phenotype, if (color_bg) bg$phenotype)))
          }
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

      p + ggplot2::labs(
          title = paste0("Donor ", input$donor, " (", ds, ")"),
          subtitle = paste0(nrow(fg), " foreground cells | ",
                            nrow(bg), " background cells"),
          x = expression(paste("X centroid (", mu, "m)")),
          y = expression(paste("Y centroid (", mu, "m)"))
        ) +
        ggplot2::theme_minimal(base_size = 18) +
        ggplot2::theme(
          legend.position = "right",
          legend.key.size = ggplot2::unit(0.7, "cm"),
          legend.title = ggplot2::element_text(size = 18, face = "bold"),
          legend.text = ggplot2::element_text(size = 15),
          plot.title = ggplot2::element_text(size = 22, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = 16, color = "#555")
        )
    }, height = 800)

    # ==== Card 3: Leiden UMAP (plotly, 1015 islets) ====
    output$leiden_not_available <- renderUI({
      if (!has_leiden()) {
        return(div(style = "color: #888; font-style: italic; margin-bottom: 10px;",
          "Leiden clustering not available in current H5AD. ",
          "Rebuild with build_h5ad_for_app.py after running Leiden clustering."))
      }
      NULL
    })

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
      if (nrow(plot_df) == 0) return(plotly_empty() %>% layout(title = "No UMAP data"))

      # Sort cluster levels numerically
      cluster_levels <- sort(unique(plot_df$cluster))
      plot_df$cluster <- factor(plot_df$cluster, levels = cluster_levels)
      n_clusters <- length(cluster_levels)
      pal_colors <- rep_len(leiden_palette, n_clusters)
      pal <- setNames(pal_colors, cluster_levels)

      res_label <- gsub("^leiden_", "", leiden_col)

      plot_ly(plot_df, x = ~umap1, y = ~umap2, color = ~cluster,
              colors = pal,
              text = ~paste0("Cluster: ", cluster, "<br>",
                             "Status: ", donor_status, "<br>",
                             islet_key),
              hoverinfo = "text",
              type = "scatter", mode = "markers",
              marker = list(size = 3, opacity = 0.7)) %>%
        layout(
          title = list(text = paste0("Leiden ", res_label, " (", nrow(plot_df), " islets)"),
                       font = list(size = 15)),
          xaxis = list(title = list(text = "UMAP 1", font = list(size = 14)),
                       zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE,
                       showline = FALSE, constrain = "domain"),
          yaxis = list(title = list(text = "UMAP 2", font = list(size = 14)),
                       zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE,
                       showline = FALSE, scaleanchor = "x", scaleratio = 1,
                       constrain = "domain"),
          legend = list(title = list(text = "Cluster", font = list(size = 14)),
                        font = list(size = 13),
                        itemsizing = "constant")
        )
    })

    # ==== Donor Status UMAP (static ggplot, mirrors trajectory tab) ====
    output$spatial_umap_donor <- renderPlot({
      req(has_leiden())
      pd <- prepared()
      comp <- pd$comp
      req("leiden_umap_1" %in% colnames(comp), "leiden_umap_2" %in% colnames(comp))

      plot_comp <- comp
      if (!is.null(input$groups) && "Donor Status" %in% colnames(plot_comp)) {
        plot_comp <- plot_comp[plot_comp$`Donor Status` %in% input$groups, , drop = FALSE]
      }

      umap1 <- suppressWarnings(as.numeric(plot_comp$leiden_umap_1))
      umap2 <- suppressWarnings(as.numeric(plot_comp$leiden_umap_2))
      ds <- if ("Donor Status" %in% colnames(plot_comp)) as.character(plot_comp$`Donor Status`) else rep("?", nrow(plot_comp))

      plot_df <- data.frame(umap1 = umap1, umap2 = umap2, donor_status = ds, stringsAsFactors = FALSE)
      plot_df <- plot_df[is.finite(plot_df$umap1) & is.finite(plot_df$umap2), , drop = FALSE]
      if (nrow(plot_df) == 0) return(NULL)

      plot_df$donor_status <- factor(plot_df$donor_status, levels = c("ND", "Aab+", "T1D"))

      ggplot(plot_df, aes(x = umap1, y = umap2, color = donor_status)) +
        geom_point(alpha = 0.6, size = 1.0) +
        scale_color_manual(values = donor_colors_reactive()) +
        scale_x_continuous(expand = expansion(mult = 0.02)) +
        scale_y_continuous(expand = expansion(mult = 0.02)) +
        coord_fixed() +
        labs(x = "UMAP 1", y = "UMAP 2", color = "Status") +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        theme_minimal(base_size = 15) +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid = element_blank(), axis.line = element_blank(),
              axis.title = element_text(size = 14),
              legend.text = element_text(size = 13),
              legend.title = element_text(size = 14))
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
      if (length(prop_cols) == 0) {
        return(plotly_empty() %>% layout(title = "No phenotype proportion data"))
      }
      if (!is.null(input$pheno_filter) && length(input$pheno_filter) > 0) {
        keep_cols <- paste0("prop_", input$pheno_filter)
        prop_cols <- intersect(prop_cols, keep_cols)
        if (length(prop_cols) == 0) {
          return(plotly_empty() %>% layout(title = "No matching phenotypes selected"))
        }
      }

      cluster <- as.character(plot_comp[[leiden_col]])

      # Compute mean phenotype proportions per cluster
      cluster_levels <- sort(unique(cluster))
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

      # Only show top phenotypes (>1% in any cluster) to reduce clutter
      max_per_pheno <- tapply(bar_df$proportion, bar_df$phenotype, max, na.rm = TRUE)
      keep_phenos <- names(max_per_pheno[max_per_pheno > 0.01])
      bar_df <- bar_df[bar_df$phenotype %in% keep_phenos, , drop = FALSE]

      # Use active phenotype palette where available
      pheno_present <- unique(bar_df$phenotype)
      pal <- palette()[pheno_present]
      pal[is.na(pal)] <- "#CCCCCC"

      bar_df$phenotype <- factor(bar_df$phenotype, levels = names(sort(max_per_pheno[keep_phenos], decreasing = TRUE)))

      plot_ly(bar_df, x = ~cluster, y = ~proportion, color = ~phenotype,
              colors = pal,
              type = "bar") %>%
        layout(
          barmode = "stack",
          title = list(text = "Mean Phenotype Composition by Cluster", font = list(size = 13)),
          xaxis = list(title = "Cluster"),
          yaxis = list(title = "Mean Proportion", range = c(0, 1)),
          legend = list(font = list(size = 9), title = list(text = ""))
        )
    })

    # ---- Download handler ----
    output$download_spatial <- downloadHandler(
      filename = function() {
        paste0("spatial_neighborhood_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        pd <- prepared()
        df <- pd$comp
        if (is.null(df) || nrow(df) == 0) { df <- data.frame() }
        else {
          if (!is.null(input$groups) && "Donor Status" %in% colnames(df))
            df <- df[df$`Donor Status` %in% input$groups, , drop = FALSE]
          if ("total_cells_peri" %in% colnames(df))
            df <- df[!is.na(df$total_cells_peri) & df$total_cells_peri > 0, , drop = FALSE]
        }
        nbr <- get_nbr_columns(colnames(df))
        keep_cols <- c("Case ID", "Donor Status", "islet_key", "islet_diam_um",
                       unlist(nbr, use.names = FALSE))
        keep_cols <- intersect(keep_cols, colnames(df))
        write.csv(df[, keep_cols, drop = FALSE], file, row.names = FALSE)
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
      # Global min cells filter
      min_cells <- input$nbr_min_cells %||% 1
      if (min_cells > 1) {
        tc_core <- if ("total_cells_core" %in% colnames(comp)) comp$total_cells_core else NA
        tc_peri <- if ("total_cells_peri" %in% colnames(comp)) comp$total_cells_peri else NA
        tc <- ifelse(is.finite(tc_core), tc_core, 0) + ifelse(is.finite(tc_peri), tc_peri, 0)
        comp <- comp[is.finite(tc) & tc >= min_cells, , drop = FALSE]
      }
      comp
    })

    # ---- Friendly labels for enrichment columns ----
    enrich_col_labels <- c(
      "enrich_z_CD8a_Tcell" = "CD8+ T-cell",
      "enrich_z_CD4_Tcell"  = "CD4+ T-cell",
      "enrich_z_T_cell"     = "T cell",
      "enrich_z_B_cell"     = "B cell",
      "enrich_z_Macrophage" = "Macrophage",
      "enrich_z_APCs"       = "APCs",
      "enrich_z_Immune"     = "Immune (all)"
    )

    # Mapping from enrichment z-score columns to core prop_* and peri_prop_* columns
    immune_type_map <- list(
      list(label = "CD8+ T-cell",  enrich = "enrich_z_CD8a_Tcell", core = "prop_CD8a Tcell",  peri = "peri_prop_CD8a_Tcell"),
      list(label = "CD4+ T-cell",  enrich = "enrich_z_CD4_Tcell",  core = "prop_CD4 Tcell",   peri = "peri_prop_CD4_Tcell"),
      list(label = "T cell",       enrich = "enrich_z_T_cell",     core = "prop_T cell",      peri = "peri_prop_T_cell"),
      list(label = "B cell",       enrich = "enrich_z_B_cell",     core = "prop_B cell",      peri = "peri_prop_B_cell"),
      list(label = "Macrophage",   enrich = "enrich_z_Macrophage", core = "prop_Macrophage",  peri = "peri_prop_Macrophage"),
      list(label = "APCs",         enrich = "enrich_z_APCs",       core = "prop_APCs",        peri = "peri_prop_APCs"),
      list(label = "Immune (all)", enrich = "enrich_z_Immune",     core = "prop_Immune",      peri = "peri_prop_Immune")
    )

    # ---- Intermediate reactive for Card B: enrichment/proportion summary ----
    enrich_summary <- reactive({
      comp <- nbr_comp()
      req(comp, "Donor Status" %in% colnames(comp))

      clip <- isTRUE(input$enrich_clip)
      stat_mode <- input$enrich_stat %||% "Median"
      region <- input$enrich_region %||% "peri"

      statuses <- unique(comp$`Donor Status`)
      comp_cols <- colnames(comp)
      rows <- list()

      for (itm in immune_type_map) {
        # Select column based on region
        if (region == "peri") {
          col <- itm$enrich  # enrichment z-score (peri)
        } else if (region == "core") {
          col <- itm$core    # core proportion
        } else {
          col <- itm$peri    # peri proportion (core+peri mode shows raw peri prop)
        }
        if (is.null(col) || !(col %in% comp_cols)) next

        for (ds in statuses) {
          vals <- comp[[col]][comp$`Donor Status` == ds]
          vals <- vals[is.finite(vals)]
          if (region == "peri" && clip) vals <- vals[vals >= -5 & vals <= 5]
          n <- length(vals)
          if (n == 0) next
          if (stat_mode == "Median") {
            z_summary <- median(vals, na.rm = TRUE)
            q <- quantile(vals, c(0.25, 0.75), na.rm = TRUE)
            z_lo <- q[1]
            z_hi <- q[2]
          } else {
            z_summary <- mean(vals, na.rm = TRUE)
            se <- sd(vals, na.rm = TRUE) / sqrt(n)
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

    # ---- Neighborhood cards UI (conditional on data availability) ----
    output$neighborhood_cards <- renderUI({
      if (!has_neighborhood()) return(NULL)

      # section_heading helper (inline — cannot share with UI function directly)
      sec_heading <- function(step, title, subtitle) {
        div(style = "margin-bottom: 14px; margin-top: 22px; padding-bottom: 8px; border-bottom: 2px solid #d0e0f0;",
          div(style = "display: flex; align-items: baseline; gap: 10px;",
            span(step,
                 style = paste0("display: inline-block; background: linear-gradient(135deg, #4477AA, #5599CC);",
                                " color: white; font-weight: 700; font-size: 14px; padding: 2px 10px;",
                                " border-radius: 12px; min-width: 28px; text-align: center;")),
            span(title, style = "font-weight: 700; font-size: 18px;")
          ),
          tags$small(subtitle, style = "color: #777; display: block; margin-top: 4px;")
        )
      }

      tagList(
        # ==== Global filter for neighborhood cards ====
        fluidRow(column(12,
          div(class = "card", style = "padding: 12px 15px; margin-bottom: 10px; margin-top: 10px; display: flex; align-items: center; gap: 20px; overflow: visible;",
            span(style = "font-weight: 600; font-size: 15px; white-space: nowrap;", "Neighborhood Analysis"),
            numericInput(ns("nbr_min_cells"), "Min cells/islet", value = 1, min = 1, max = 100, step = 1, width = "130px"),
            sliderInput(ns("nbr_pt_size"), "Point size", min = 2, max = 12, value = 5, step = 1, width = "150px"),
            sliderInput(ns("nbr_pt_alpha"), "Opacity", min = 0.1, max = 1.0, value = 0.4, step = 0.05, width = "150px"),
            uiOutput(ns("nbr_islet_count"))
          )
        )),

        # ==== Card A: Immune Infiltration Overview ====
        fluidRow(column(12, sec_heading(
          "A", "Immune Infiltration",
          "How does immune infiltration vary across disease stages? Compare peri-islet and core immune fractions."
        ))),
        fluidRow(style = "display: flex; flex-wrap: wrap;",
          column(6,
            div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%; box-sizing: border-box;",
              div(style = "display: flex; gap: 15px; align-items: flex-end; flex-wrap: wrap; margin-bottom: 8px;",
                div(style = "flex: 1; min-width: 180px;",
                  selectInput(ns("infiltration_metric"), "Metric",
                              choices = c("Immune fraction (peri)" = "immune_frac_peri",
                                          "Immune fraction (core)" = "immune_frac_core",
                                          "T-cell density (peri)" = "tcell_density_peri",
                                          "CD8/Macrophage ratio" = "cd8_to_macro_ratio",
                                          "Peri/core immune ratio" = "immune_ratio"),
                              selected = "immune_frac_peri", width = "100%")
                ),
                checkboxInput(ns("infiltration_outliers"), "Remove outliers", value = FALSE)
              ),
              plotlyOutput(ns("infiltration_bars"), height = "400px")
            )
          ),
          column(6,
            div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%; box-sizing: border-box;",
              h5("Peri vs Core Immune Fraction", style = "font-size: 15px; margin-top: 0;"),
              plotlyOutput(ns("infiltration_scatter"), height = "400px")
            )
          )
        ),
        fluidRow(column(12,
          div(style = "background: #f0f6ff; padding: 10px 15px; border-radius: 6px; margin-bottom: 20px; font-size: 13px; color: #555;",
            tags$em("Immune infiltration quantifies the proportion of immune cells among all cells in the peri-islet zone (20\u00b5m expansion) and islet core. Higher fractions indicate increased immune surveillance or active infiltration.")
          )
        )),

        # ==== Card B: Immune Cell Enrichment by Type ====
        fluidRow(column(12, sec_heading(
          "B", "Immune Cell Composition & Enrichment",
          "Which immune cell types are present in islet core vs peri-islet zones, and which are enriched vs tissue-wide background?"
        ))),
        fluidRow(style = "display: flex; flex-wrap: wrap;",
          column(6,
            div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%; box-sizing: border-box;",
              div(style = "display: flex; gap: 15px; align-items: center; flex-wrap: wrap; margin-bottom: 8px;",
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
                )
              ),
              plotlyOutput(ns("enrichment_bars"), height = "420px")
            )
          ),
          column(6,
            div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%; box-sizing: border-box;",
              h5("Heatmap", style = "font-size: 15px; margin-top: 0;"),
              plotlyOutput(ns("enrichment_heatmap"), height = "420px")
            )
          )
        ),
        fluidRow(column(12,
          div(style = "background: #f0f6ff; padding: 10px 15px; border-radius: 6px; margin-bottom: 20px; font-size: 13px; color: #555;",
            tags$em("Enrichment z-scores (peri-islet mode) compare observed immune cell counts in the peri-islet zone vs expected counts based on tissue-wide proportions (Poisson model). z > 0 = enriched; z < 0 = depleted. Core and peri-islet proportion modes show raw cell type fractions (0\u20131) within each compartment.")
          )
        )),

        # ==== Card C: Immune Proximity to Islets ====
        fluidRow(column(12, sec_heading(
          "C", "Immune Cell Proximity",
          "How close are immune cells to islet cores? Shorter distances may indicate active immune targeting."
        ))),
        fluidRow(style = "display: flex; flex-wrap: wrap;",
          column(6,
            div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%; box-sizing: border-box;",
              selectInput(ns("distance_metric"), "Distance to nearest",
                          choices = c("Any immune cell" = "min_dist_immune_mean",
                                      "Macrophage" = "min_dist_Macrophage",
                                      "CD8+ T-cell (sparse)" = "min_dist_CD8a_Tcell"),
                          selected = "min_dist_immune_mean", width = "100%"),
              checkboxInput(ns("distance_clip"), "Clip top 1% outliers", value = TRUE),
              plotlyOutput(ns("distance_boxplot"), height = "420px")
            )
          ),
          column(6,
            div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%; box-sizing: border-box;",
              selectInput(ns("kde_immune_type"), "Immune cell type",
                          choices = c("All immune" = "all",
                                      "Macrophage" = "Macrophage",
                                      "APCs" = "APCs",
                                      "CD8+ T-cell" = "CD8a Tcell",
                                      "T cell" = "T cell",
                                      "B cell" = "B cell",
                                      "CD4+ T-cell" = "CD4 Tcell"),
                          selected = "all", width = "100%"),
              plotlyOutput(ns("immune_distance_kde"), height = "420px")
            )
          )
        ),
        fluidRow(column(12,
          div(style = "background: #f0f6ff; padding: 10px 15px; border-radius: 6px; margin-bottom: 20px; font-size: 13px; color: #555;",
            tags$em("Left: Distance metrics measure minimum Euclidean distance (\u00b5m) from islet core centroid to nearest peri-islet immune cells. NAs indicate no cells of that type in the peri-islet zone. Right: KDE of signed distance from islet boundary for individual immune cells. Negative = inside islet (core), positive = outside (peri-islet zone). The dashed line at zero marks the islet boundary.")
          )
        ))
      )
    })

    # ---- Islet count display ----
    output$nbr_islet_count <- renderUI({
      comp <- nbr_comp()
      n <- if (!is.null(comp)) nrow(comp) else 0
      span(style = "color: #555; font-size: 13px; white-space: nowrap;",
           paste0(formatC(n, big.mark = ","), " islets"))
    })

    # ==== Card A-Left: Immune Metric Bar Chart (donor-level) ====
    output$infiltration_bars <- renderPlotly({
      comp <- nbr_comp()
      req(comp, "Donor Status" %in% colnames(comp), "Case ID" %in% colnames(comp))
      metric <- input$infiltration_metric %||% "immune_frac_peri"
      req(metric %in% colnames(comp))

      # Metric label for title
      metric_labels <- c(
        "immune_frac_peri" = "Immune Fraction (Peri)",
        "immune_frac_core" = "Immune Fraction (Core)",
        "tcell_density_peri" = "T-cell Density (Peri)",
        "cd8_to_macro_ratio" = "CD8/Macrophage Ratio",
        "immune_ratio" = "Peri/Core Immune Ratio"
      )
      metric_label <- metric_labels[metric]
      if (is.na(metric_label)) metric_label <- metric

      df <- data.frame(
        value = comp[[metric]],
        status = comp$`Donor Status`,
        case_id = comp$`Case ID`,
        stringsAsFactors = FALSE
      )
      df <- df[is.finite(df$value), , drop = FALSE]
      if (nrow(df) == 0) return(plotly_empty() %>% layout(title = "No data"))

      # Outlier removal (per-group IQR on non-zero values to handle zero-inflated metrics)
      if (isTRUE(input$infiltration_outliers)) {
        keep <- logical(nrow(df))
        for (grp in unique(df$status)) {
          idx <- df$status == grp
          vals <- df$value[idx]
          nz <- vals[vals != 0]
          if (length(nz) < 4) { keep[idx] <- TRUE; next }
          q <- quantile(nz, c(0.25, 0.75), na.rm = TRUE)
          iqr <- q[2] - q[1]
          lo <- q[1] - 1.5 * iqr
          hi <- q[2] + 1.5 * iqr
          keep[idx] <- vals == 0 | (vals >= lo & vals <= hi)
        }
        df <- df[keep, , drop = FALSE]
        if (nrow(df) == 0) return(plotly_empty() %>% layout(title = "No data after outlier removal"))
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

      dcols <- donor_colors_reactive()

      # Kruskal-Wallis on donor-level means (N=15)
      kw_p <- tryCatch({
        kt <- kruskal.test(value ~ status, data = donor_df)
        kt$p.value
      }, error = function(e) NA_real_)
      p_label <- if (is.finite(kw_p)) paste0("KW p = ", formatC(kw_p, format = "g", digits = 3)) else ""

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
      p %>% layout(
        title = list(text = paste0(metric_label, "<br><sup>Donor-level means \u00b1 SEM (N=donors) | ", p_label, "</sup>"),
                     font = list(size = 14)),
        xaxis = list(title = "", categoryorder = "array",
                     categoryarray = c("ND", "Aab+", "T1D")),
        yaxis = list(title = metric_label),
        barmode = "overlay",
        showlegend = FALSE
      )
    })

    # ==== Card A-Right: Peri vs Core Scatter ====
    output$infiltration_scatter <- renderPlotly({
      comp <- nbr_comp()
      req(comp, "Donor Status" %in% colnames(comp))
      req("immune_frac_peri" %in% colnames(comp), "immune_frac_core" %in% colnames(comp))

      df <- data.frame(
        peri = comp$immune_frac_peri,
        core = comp$immune_frac_core,
        status = comp$`Donor Status`,
        islet_key = if ("islet_key" %in% colnames(comp)) comp$islet_key else "",
        stringsAsFactors = FALSE
      )
      df <- df[is.finite(df$peri) & is.finite(df$core), , drop = FALSE]
      if (nrow(df) == 0) return(plotly_empty() %>% layout(title = "No data"))

      status_order <- c("ND", "Aab+", "T1D")
      df$status <- factor(df$status, levels = intersect(status_order, unique(df$status)))
      dcols <- donor_colors_reactive()

      max_val <- max(c(df$peri, df$core), na.rm = TRUE) * 1.05

      pt_sz <- input$nbr_pt_size %||% 5
      pt_al <- input$nbr_pt_alpha %||% 0.4

      p <- plot_ly()
      for (s in levels(df$status)) {
        sub <- df[df$status == s, , drop = FALSE]
        if (nrow(sub) == 0) next
        p <- p %>% add_trace(
          data = sub, x = ~peri, y = ~core,
          text = ~paste0(islet_key, "<br>Status: ", status,
                        "<br>Peri: ", round(peri, 4),
                        "<br>Core: ", round(core, 4)),
          hoverinfo = "text", name = s,
          type = "scatter", mode = "markers",
          marker = list(size = pt_sz, opacity = pt_al, color = dcols[s])
        )
      }
      # Diagonal y=x reference
      p %>% add_trace(
        x = c(0, max_val), y = c(0, max_val),
        type = "scatter", mode = "lines",
        line = list(dash = "dash", color = "#999", width = 1),
        showlegend = FALSE, hoverinfo = "none"
      ) %>% layout(
        title = list(text = "Peri vs Core Immune Fraction", font = list(size = 14)),
        xaxis = list(title = "Immune Fraction (Peri)", range = c(0, max_val)),
        yaxis = list(title = "Immune Fraction (Core)", range = c(0, max_val)),
        legend = list(title = list(text = "Status"))
      )
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
        list(list(type = "line", x0 = -0.5, x1 = 10, y0 = 0, y1 = 0,
                  line = list(color = "#999", dash = "dash", width = 1)))
      } else list()

      p %>% layout(
        barmode = "group",
        title = list(text = paste0(title_prefix, " (error: ", error_label, ")"),
                     font = list(size = 14)),
        xaxis = list(title = "", tickangle = -30,
                     categoryorder = "array",
                     categoryarray = unique(es$cell_type)),
        yaxis = list(title = y_label),
        legend = list(title = list(text = "Status")),
        shapes = shapes
      )
    })

    # ==== Card B-Right: Enrichment/Proportion Heatmap ====
    output$enrichment_heatmap <- renderPlotly({
      es <- enrich_summary()
      req(es)

      region <- input$enrich_region %||% "peri"

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

      plot_ly(
        x = colnames(mat), y = rownames(mat), z = mat,
        type = "heatmap",
        colorscale = colorscale,
        zmin = zmin, zmax = zmax,
        text = text_mat, texttemplate = "%{text}",
        hovertemplate = paste0("%{x}<br>%{y}<br>", cb_title, " = %{z:.4f}<extra></extra>"),
        showscale = TRUE,
        colorbar = list(title = cb_title)
      ) %>% layout(
        title = list(text = title_text, font = list(size = 14)),
        xaxis = list(title = "", tickangle = -30),
        yaxis = list(title = "", autorange = "reversed")
      )
    })

    # ==== Card C-Left: Distance Box Plots ====
    output$distance_boxplot <- renderPlotly({
      comp <- nbr_comp()
      req(comp, "Donor Status" %in% colnames(comp))
      metric <- input$distance_metric %||% "min_dist_immune_mean"
      req(metric %in% colnames(comp))

      metric_labels <- c(
        "min_dist_immune_mean" = "Min Distance to Any Immune Cell",
        "min_dist_Macrophage" = "Min Distance to Macrophage",
        "min_dist_CD8a_Tcell" = "Min Distance to CD8+ T-cell"
      )
      metric_label <- metric_labels[metric]
      if (is.na(metric_label)) metric_label <- metric

      df <- data.frame(
        value = comp[[metric]],
        status = comp$`Donor Status`,
        stringsAsFactors = FALSE
      )

      # Count non-NA per group before filtering
      status_order <- c("ND", "Aab+", "T1D")
      n_total <- tapply(df$value, df$status, length)
      n_valid <- tapply(df$value, df$status, function(x) sum(is.finite(x)))
      na_rate <- 1 - sum(is.finite(df$value)) / nrow(df)

      df <- df[is.finite(df$value), , drop = FALSE]
      if (nrow(df) == 0) return(plotly_empty() %>% layout(title = paste0(metric_label, " \u2014 all NA")))

      # Clip top 1% outliers if requested
      if (isTRUE(input$distance_clip) && nrow(df) > 0) {
        q99 <- quantile(df$value, 0.99, na.rm = TRUE)
        df <- df[df$value <= q99, , drop = FALSE]
      }

      df$status <- factor(df$status, levels = intersect(status_order, unique(df$status)))
      dcols <- donor_colors_reactive()
      pt_sz <- input$nbr_pt_size %||% 5
      pt_al <- input$nbr_pt_alpha %||% 0.4

      # Subtitle with N per group
      n_labels <- sapply(levels(df$status), function(s) {
        nv <- n_valid[s]
        nt <- n_total[s]
        if (is.na(nv)) nv <- 0
        if (is.na(nt)) nt <- 0
        paste0(s, ": ", nv, "/", nt)
      })
      subtitle <- paste0("Non-NA: ", paste(n_labels, collapse = ", "),
                          " (", round(na_rate * 100, 0), "% NA overall)")

      p <- plot_ly()
      for (s in levels(df$status)) {
        sub <- df[df$status == s, , drop = FALSE]
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
      p %>% layout(
        title = list(text = paste0(metric_label, "<br><sup>", subtitle, "</sup>"),
                     font = list(size = 14)),
        xaxis = list(title = "", categoryorder = "array",
                     categoryarray = c("ND", "Aab+", "T1D")),
        yaxis = list(title = paste0("Distance (\u00b5m)")),
        showlegend = FALSE
      )
    })

    # ==== Card C-Right: Immune Cell Distance KDE (signed distance from islet boundary) ====
    # Lazy-load immune distance CSV (cached in module scope)
    .kde_data_env <- new.env(parent = emptyenv())
    load_immune_kde_data <- function() {
      if (!is.null(.kde_data_env$data)) return(.kde_data_env$data)
      kde_path <- file.path("..", "..", "data", "immune_distance_kde.csv")
      if (!file.exists(kde_path)) return(NULL)
      .kde_data_env$data <- read.csv(kde_path, stringsAsFactors = FALSE)
      .kde_data_env$data
    }

    output$immune_distance_kde <- renderPlotly({
      kde_df <- load_immune_kde_data()
      req(kde_df)

      cell_type <- input$kde_immune_type %||% "all"
      # Filter to selected immune type
      if (cell_type == "all") {
        df <- kde_df
        type_label <- "All Immune Cells"
      } else {
        df <- kde_df[kde_df$phenotype == cell_type, , drop = FALSE]
        type_label <- cell_type
      }

      # Filter to selected donor status groups
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

      # Compute KDE per disease stage
      kde_list <- list()
      y_max <- 0
      for (s in levels(df$donor_status)) {
        vals <- df$signed_dist[df$donor_status == s]
        if (length(vals) < 2) next
        d <- density(vals, n = 512, from = q_lo, to = q_hi)
        kde_list[[s]] <- d
        y_max <- max(y_max, max(d$y))
      }
      req(length(kde_list) > 0)

      p <- plot_ly()
      for (s in names(kde_list)) {
        d <- kde_list[[s]]
        p <- p %>% add_trace(
          x = d$x, y = d$y,
          type = "scatter", mode = "lines",
          name = s,
          line = list(color = dcols[s], width = 2.5),
          fill = "tozeroy",
          fillcolor = paste0(dcols[s], "22"),
          hoverinfo = "text",
          text = paste0(s, "<br>dist = ", round(d$x, 1), " \u00b5m",
                       "<br>density = ", round(d$y, 4))
        )
      }
      # Vertical line at zero (islet boundary)
      p <- p %>% add_trace(
        x = c(0, 0), y = c(0, y_max * 1.05),
        type = "scatter", mode = "lines",
        line = list(color = "#333", width = 1.5, dash = "dash"),
        showlegend = FALSE, hoverinfo = "none"
      )
      p %>% layout(
        title = list(text = paste0(type_label, " Distance from Islet"),
                     font = list(size = 14)),
        xaxis = list(title = "Signed distance (\u00b5m) \u2014 inside islet \u2190 0 \u2192 outside",
                     zeroline = FALSE),
        yaxis = list(title = "Density"),
        legend = list(title = list(text = "Status"))
      )
    })

  })
}
