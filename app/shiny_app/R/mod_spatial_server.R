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

    # 14-color qualitative palette for Leiden clusters
    leiden_palette <- c(
      "#1f77b4", "#ff7f0e", "#66c2a5", "#d62728", "#8da0cb",
      "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
      "#aec7e8", "#ffbb78", "#98df8a", "#ff9896"
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
          pal <- setNames(leiden_palette[seq_len(min(n_clusters, length(leiden_palette)))],
                          cluster_levels)
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
          # Phenotype coloring — include bg phenotypes if colored
          all_phenos <- sort(unique(fg$phenotype))
          if (color_bg && nrow(bg) > 0) {
            all_phenos <- sort(unique(c(all_phenos, bg$phenotype)))
          }
          pal <- palette()[all_phenos]
          pal[is.na(pal)] <- "#CCCCCC"

          p <- p + ggplot2::geom_point(
            data = fg,
            ggplot2::aes(x = X_centroid, y = Y_centroid, color = phenotype),
            size = 0.4, alpha = 0.6,
            inherit.aes = FALSE
          ) +
          ggplot2::scale_color_manual(values = pal, name = "Phenotype", na.value = "#CCCCCC",
                                        guide = ggplot2::guide_legend(override.aes = list(size = 4)))
        }
      }

      ds <- donor_status()
      p + ggplot2::coord_fixed() +
        ggplot2::scale_y_reverse() +
        ggplot2::labs(
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
      pal <- setNames(leiden_palette[seq_len(min(n_clusters, length(leiden_palette)))],
                      cluster_levels)

      res_label <- gsub("^leiden_", "", leiden_col)

      plot_ly(plot_df, x = ~umap1, y = ~umap2, color = ~cluster,
              colors = pal,
              text = ~paste0("Cluster: ", cluster, "<br>",
                             "Status: ", donor_status, "<br>",
                             islet_key),
              hoverinfo = "text",
              type = "scatter", mode = "markers",
              marker = list(size = 7, opacity = 0.7)) %>%
        layout(
          title = list(text = paste0("Leiden ", res_label, " (", nrow(plot_df), " islets)"),
                       font = list(size = 14)),
          xaxis = list(title = "UMAP 1", zeroline = FALSE),
          yaxis = list(title = "UMAP 2", zeroline = FALSE),
          legend = list(title = list(text = "Cluster"))
        )
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

      # Get phenotype proportion columns (prop_*)
      prop_cols <- grep("^prop_", colnames(plot_comp), value = TRUE)
      if (length(prop_cols) == 0) {
        return(plotly_empty() %>% layout(title = "No phenotype proportion data"))
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
  })
}
