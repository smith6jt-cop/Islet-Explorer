# ---------- Spatial Neighborhood module server ----------
# Exports: spatial_server(id, prepared)
#
# Dependencies:
#   prepared()$comp — must contain peri_prop_*, immune_*, enrich_z_* columns
#   cohens_d, eta_squared, pairwise_wilcox — from utils_stats.R

spatial_server <- function(id, prepared) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Donor status color palette
    donor_colors <- c("ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd")

    # ---- Helpers: column identification ----
    get_nbr_columns <- function(comp_cols) {
      list(
        peri_prop  = grep("^peri_prop_", comp_cols, value = TRUE),
        peri_count = grep("^peri_count_", comp_cols, value = TRUE),
        immune     = intersect(c("immune_frac_peri", "immune_frac_core", "immune_ratio",
                                 "cd8_to_macro_ratio", "tcell_density_peri"), comp_cols),
        enrich     = grep("^enrich_z_", comp_cols, value = TRUE),
        distance   = grep("^min_dist_", comp_cols, value = TRUE)
      )
    }

    has_neighborhood <- reactive({
      pd <- prepared()
      if (is.null(pd$comp)) return(FALSE)
      nbr <- get_nbr_columns(colnames(pd$comp))
      length(nbr$peri_prop) > 0
    })

    # ---- Filtered data reactive ----
    filtered_df <- reactive({
      req(has_neighborhood())
      pd <- prepared()
      df <- pd$comp
      if (is.null(df) || nrow(df) == 0) return(NULL)

      # Filter by donor status
      if (!is.null(input$groups) && "Donor Status" %in% colnames(df)) {
        df <- df[df$`Donor Status` %in% input$groups, , drop = FALSE]
      }

      # Filter by diameter range
      if (!is.null(input$diam_range) && "islet_diam_um" %in% colnames(df)) {
        d <- suppressWarnings(as.numeric(df$islet_diam_um))
        df <- df[!is.na(d) & d >= input$diam_range[1] & d <= input$diam_range[2], , drop = FALSE]
      }

      # Only keep rows with peri-islet data
      if ("total_cells_peri" %in% colnames(df)) {
        df <- df[!is.na(df$total_cells_peri) & df$total_cells_peri > 0, , drop = FALSE]
      }

      df
    })

    # ---- Dynamic feature selector ----
    output$feature_selector <- renderUI({
      req(has_neighborhood())
      comp_cols <- colnames(prepared()$comp)
      nbr <- get_nbr_columns(comp_cols)

      choices <- switch(input$metric_category %||% "Peri-Islet Composition",
        "Peri-Islet Composition" = {
          labs <- gsub("^peri_prop_", "", nbr$peri_prop)
          setNames(nbr$peri_prop, labs)
        },
        "Immune Infiltration" = {
          labs <- c(
            immune_frac_peri = "Immune fraction (peri)",
            immune_frac_core = "Immune fraction (core)",
            immune_ratio = "Peri/core immune ratio",
            cd8_to_macro_ratio = "CD8/macrophage ratio",
            tcell_density_peri = "T-cell density (peri)"
          )
          setNames(nbr$immune, labs[nbr$immune])
        },
        "Enrichment Z-scores" = {
          labs <- gsub("^enrich_z_", "", nbr$enrich)
          setNames(nbr$enrich, labs)
        },
        "Distance Metrics" = {
          labs <- gsub("^min_dist_", "", nbr$distance)
          setNames(nbr$distance, labs)
        }
      )

      if (length(choices) == 0) choices <- c("(none available)" = "")

      selectInput(ns("feature"), "Feature", choices = choices)
    })

    # ---- Card 1: Overview Banner ----
    output$overview_banner <- renderUI({
      if (!has_neighborhood()) {
        return(div(style = "color: #888; font-style: italic;",
          "Neighborhood metrics not available. Run compute_neighborhood_metrics.py and rebuild H5AD."))
      }
      df <- filtered_df()
      if (is.null(df) || nrow(df) == 0) {
        return(div("No data matching current filters."))
      }
      n_total <- nrow(df)
      n_by_group <- if ("Donor Status" %in% colnames(df)) {
        paste(sapply(c("ND", "Aab+", "T1D"), function(g) {
          n <- sum(df$`Donor Status` == g, na.rm = TRUE)
          paste0(g, "=", n)
        }), collapse = ", ")
      } else ""

      tags$p(style = "margin: 0; font-size: 14px;",
        strong(n_total), " islets with peri-islet data",
        if (nzchar(n_by_group)) paste0(" (", n_by_group, ")")
      )
    })

    # ---- Card 2: Feature Distribution (violin/box by disease stage) ----
    output$distribution_plot <- renderPlotly({
      df <- filtered_df()
      feat <- input$feature
      req(df, feat, feat %in% colnames(df), "Donor Status" %in% colnames(df))

      vals <- suppressWarnings(as.numeric(df[[feat]]))
      ds <- as.character(df$`Donor Status`)
      plot_df <- data.frame(value = vals, group = ds, stringsAsFactors = FALSE)
      plot_df <- plot_df[is.finite(plot_df$value), , drop = FALSE]
      if (nrow(plot_df) == 0) return(plotly_empty() %>% layout(title = "No data"))

      plot_df$group <- factor(plot_df$group, levels = c("ND", "Aab+", "T1D"))
      feat_label <- gsub("^peri_prop_|^enrich_z_|^min_dist_", "", feat)
      feat_label <- gsub("_", " ", feat_label)

      p <- plot_ly(plot_df, y = ~value, x = ~group, color = ~group,
                   colors = donor_colors[levels(plot_df$group)],
                   type = "violin", box = list(visible = TRUE),
                   meanline = list(visible = TRUE), points = "outliers") %>%
        layout(
          title = list(text = feat_label, font = list(size = 14)),
          xaxis = list(title = ""),
          yaxis = list(title = feat_label),
          showlegend = FALSE
        )
      p
    })

    # ---- Card 3: Core vs Peri Comparison ----
    output$core_peri_plot <- renderPlotly({
      df <- filtered_df()
      req(df, "Donor Status" %in% colnames(df))

      # Find matching core/peri immune fraction columns
      if (!all(c("immune_frac_peri", "immune_frac_core") %in% colnames(df))) {
        return(plotly_empty() %>% layout(title = "Core/peri immune data not available"))
      }

      plot_df <- data.frame(
        core = suppressWarnings(as.numeric(df$immune_frac_core)),
        peri = suppressWarnings(as.numeric(df$immune_frac_peri)),
        group = as.character(df$`Donor Status`),
        stringsAsFactors = FALSE
      )
      plot_df <- plot_df[is.finite(plot_df$core) & is.finite(plot_df$peri), , drop = FALSE]
      if (nrow(plot_df) == 0) return(plotly_empty() %>% layout(title = "No data"))

      # Reshape to long format for grouped box plot
      long_df <- rbind(
        data.frame(region = "Core", value = plot_df$core, group = plot_df$group),
        data.frame(region = "Peri-Islet", value = plot_df$peri, group = plot_df$group)
      )
      long_df$group <- factor(long_df$group, levels = c("ND", "Aab+", "T1D"))

      p <- plot_ly(long_df, x = ~group, y = ~value, color = ~region,
                   type = "box") %>%
        layout(
          title = list(text = "Immune Fraction: Core vs Peri-Islet", font = list(size = 14)),
          xaxis = list(title = ""),
          yaxis = list(title = "Immune cell fraction"),
          boxmode = "group"
        )
      p
    })

    # ---- Card 4: Peri-Islet Phenotype Heatmap ----
    output$phenotype_heatmap <- renderPlotly({
      df <- filtered_df()
      req(df, "Donor Status" %in% colnames(df))

      comp_cols <- colnames(df)
      peri_prop_cols <- grep("^peri_prop_", comp_cols, value = TRUE)
      if (length(peri_prop_cols) == 0) {
        return(plotly_empty() %>% layout(title = "No peri-islet proportion data"))
      }

      # Compute mean proportion per disease stage
      groups <- c("ND", "Aab+", "T1D")
      mat <- matrix(NA_real_, nrow = length(peri_prop_cols), ncol = length(groups),
                    dimnames = list(peri_prop_cols, groups))
      for (g in groups) {
        sub <- df[df$`Donor Status` == g, peri_prop_cols, drop = FALSE]
        if (nrow(sub) > 0) {
          mat[, g] <- colMeans(sub, na.rm = TRUE)
        }
      }

      # Clean labels
      row_labels <- gsub("^peri_prop_", "", rownames(mat))
      row_labels <- gsub("_", " ", row_labels)

      plot_ly(
        z = mat, x = groups, y = row_labels,
        type = "heatmap",
        colorscale = list(c(0, "#f7fbff"), c(0.5, "#6baed6"), c(1, "#08306b")),
        colorbar = list(title = "Mean\nProportion")
      ) %>%
        layout(
          title = list(text = "Mean Peri-Islet Phenotype Proportions", font = list(size = 14)),
          xaxis = list(title = ""),
          yaxis = list(title = "", tickfont = list(size = 10))
        )
    })

    # ---- Card 5: Immune Enrichment ----
    output$enrichment_plot <- renderPlotly({
      df <- filtered_df()
      req(df, "Donor Status" %in% colnames(df))

      enrich_cols <- grep("^enrich_z_", colnames(df), value = TRUE)
      if (length(enrich_cols) == 0) {
        return(plotly_empty() %>% layout(title = "No enrichment data"))
      }

      groups <- c("ND", "Aab+", "T1D")
      plot_rows <- list()
      for (g in groups) {
        sub <- df[df$`Donor Status` == g, , drop = FALSE]
        if (nrow(sub) == 0) next
        for (ec in enrich_cols) {
          vals <- suppressWarnings(as.numeric(sub[[ec]]))
          vals <- vals[is.finite(vals)]
          if (length(vals) > 0) {
            plot_rows[[length(plot_rows) + 1]] <- data.frame(
              cell_type = gsub("^enrich_z_", "", ec),
              group = g,
              mean_z = mean(vals),
              se = sd(vals) / sqrt(length(vals)),
              stringsAsFactors = FALSE
            )
          }
        }
      }

      if (length(plot_rows) == 0) {
        return(plotly_empty() %>% layout(title = "No enrichment data"))
      }

      bar_df <- do.call(rbind, plot_rows)
      bar_df$cell_type <- gsub("_", " ", bar_df$cell_type)
      bar_df$group <- factor(bar_df$group, levels = groups)

      plot_ly(bar_df, x = ~cell_type, y = ~mean_z, color = ~group,
              colors = donor_colors[groups],
              error_y = list(type = "data", array = ~se),
              type = "bar") %>%
        layout(
          title = list(text = "Enrichment Z-scores (Peri-Islet vs Tissue-wide)", font = list(size = 14)),
          xaxis = list(title = "", tickangle = 45),
          yaxis = list(title = "Mean enrichment z-score"),
          barmode = "group",
          legend = list(title = list(text = ""))
        )
    })

    # ---- Card 6: Pseudotime Correlation ----
    output$pseudotime_plot <- renderPlotly({
      df <- filtered_df()
      feat <- input$feature
      req(df, feat, feat %in% colnames(df))

      # Check for pseudotime (DPT) — it's in the prepared trajectory data
      if (!"dpt_pseudotime" %in% colnames(df)) {
        return(plotly_empty() %>% layout(
          title = "Pseudotime not available in composition data",
          annotations = list(list(
            text = "Pseudotime correlation requires DPT in comp. See Trajectory tab.",
            xref = "paper", yref = "paper", x = 0.5, y = 0.5, showarrow = FALSE
          ))
        ))
      }

      pt <- suppressWarnings(as.numeric(df$dpt_pseudotime))
      vals <- suppressWarnings(as.numeric(df[[feat]]))
      ds <- as.character(df$`Donor Status`)
      plot_df <- data.frame(pt = pt, value = vals, group = ds, stringsAsFactors = FALSE)
      plot_df <- plot_df[is.finite(plot_df$pt) & is.finite(plot_df$value), , drop = FALSE]

      if (nrow(plot_df) == 0) {
        return(plotly_empty() %>% layout(title = "No data for pseudotime correlation"))
      }

      feat_label <- gsub("^peri_prop_|^enrich_z_|^min_dist_", "", feat)
      feat_label <- gsub("_", " ", feat_label)
      plot_df$group <- factor(plot_df$group, levels = c("ND", "Aab+", "T1D"))

      cor_test <- tryCatch(cor.test(plot_df$pt, plot_df$value, method = "kendall"),
                           error = function(e) NULL)
      subtitle <- if (!is.null(cor_test)) {
        sprintf("Kendall \u03c4 = %.3f, p = %.2e", cor_test$estimate, cor_test$p.value)
      } else ""

      plot_ly(plot_df, x = ~pt, y = ~value, color = ~group,
              colors = donor_colors[levels(plot_df$group)],
              type = "scatter", mode = "markers",
              marker = list(size = 5, opacity = 0.6)) %>%
        layout(
          title = list(text = paste0(feat_label, " vs Pseudotime<br><sub>", subtitle, "</sub>"),
                       font = list(size = 14)),
          xaxis = list(title = "DPT Pseudotime"),
          yaxis = list(title = feat_label),
          legend = list(title = list(text = ""))
        )
    })

    # ---- Card 7: Statistical Tests ----
    output$test_table <- renderTable({
      df <- filtered_df()
      feat <- input$feature
      req(df, feat, feat %in% colnames(df), "Donor Status" %in% colnames(df))

      vals <- suppressWarnings(as.numeric(df[[feat]]))
      groups <- as.character(df$`Donor Status`)
      test_df <- data.frame(value = vals, group = groups, stringsAsFactors = FALSE)
      test_df <- test_df[is.finite(test_df$value), , drop = FALSE]

      if (nrow(test_df) < 3) return(data.frame(Note = "Insufficient data"))

      # Global test (Kruskal-Wallis)
      kw <- tryCatch(kruskal.test(value ~ group, data = test_df), error = function(e) NULL)

      # Pairwise tests
      pairs <- pairwise_wilcox(test_df, "group", "value")

      # Build results table
      rows <- list()
      if (!is.null(kw)) {
        rows[[1]] <- data.frame(
          Test = "Kruskal-Wallis (global)", Comparison = "All groups",
          Statistic = round(kw$statistic, 2), `p-value` = formatC(kw$p.value, format = "e", digits = 2),
          `Effect Size` = "", check.names = FALSE, stringsAsFactors = FALSE
        )
      }
      if (!is.null(pairs) && nrow(pairs) > 0) {
        for (i in seq_len(nrow(pairs))) {
          g1 <- as.character(pairs$group1[i])
          g2 <- as.character(pairs$group2[i])
          cd <- tryCatch({
            v1 <- test_df$value[test_df$group == g1]
            v2 <- test_df$value[test_df$group == g2]
            cohens_d(v1, v2)
          }, error = function(e) list(d = NA, ci_lower = NA, ci_upper = NA))

          rows[[length(rows) + 1]] <- data.frame(
            Test = "Wilcoxon (pairwise)", Comparison = paste(g1, "vs", g2),
            Statistic = NA_real_,
            `p-value` = formatC(pairs$p_value[i], format = "e", digits = 2),
            `Effect Size` = if (!is.na(cd$d)) sprintf("d=%.2f [%.2f,%.2f]", cd$d, cd$ci_lower, cd$ci_upper) else "",
            check.names = FALSE, stringsAsFactors = FALSE
          )
        }
      }

      if (length(rows) == 0) return(data.frame(Note = "Tests could not be computed"))
      do.call(rbind, rows)
    }, striped = TRUE, hover = TRUE, spacing = "s")

    output$methods_text <- renderUI({
      tags$div(style = "font-size: 12px; color: #555;",
        tags$p("Neighborhood metrics quantify the cellular microenvironment in the 20\u00b5m ",
               "expansion zone around each islet (peri-islet region)."),
        tags$ul(
          tags$li("Peri-islet composition: proportion of each phenotype among peri-islet cells"),
          tags$li("Immune infiltration: fraction of immune cells in peri vs core regions"),
          tags$li("Enrichment z-scores: Poisson model comparing peri-islet proportion to tissue-wide baseline"),
          tags$li("Distance metrics: minimum distance from islet centroid to nearest immune cell subtypes")
        ),
        tags$p("Statistical tests: Kruskal-Wallis for global comparisons, ",
               "pairwise Wilcoxon with BH correction. Effect sizes: Cohen's d with 95% CI.")
      )
    })

    # ---- Download handler ----
    output$download_spatial <- downloadHandler(
      filename = function() {
        paste0("spatial_neighborhood_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        df <- filtered_df()
        if (is.null(df)) df <- data.frame()
        nbr <- get_nbr_columns(colnames(df))
        keep_cols <- c("Case ID", "Donor Status", "islet_key", "islet_diam_um",
                       unlist(nbr, use.names = FALSE))
        keep_cols <- intersect(keep_cols, colnames(df))
        write.csv(df[, keep_cols, drop = FALSE], file, row.names = FALSE)
      }
    )
  })
}
