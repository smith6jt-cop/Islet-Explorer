statistics_server <- function(id, raw_df, summary_df, get_selection_description) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    stats_run <- eventReactive(input$run_tests, {
      rdf <- raw_df()
      if (is.null(rdf) || !nrow(rdf)) return(NULL)

      if (isTRUE(input$stats_remove_outliers)) {
        value_mean <- mean(rdf$value, na.rm = TRUE)
        value_sd <- sd(rdf$value, na.rm = TRUE)
        outlier_threshold <- 3
        rdf$is_outlier <- abs(rdf$value - value_mean) > (outlier_threshold * value_sd)
        n_outliers <- sum(rdf$is_outlier, na.rm = TRUE)
        rdf <- rdf[!rdf$is_outlier | is.na(rdf$is_outlier), ]
        cat(sprintf("[STATISTICS] Removed %d outliers (>3 SD)\n", n_outliers))
      }

      rdf$donor_status <- factor(rdf$donor_status, levels = c("ND","Aab+","T1D"))
      fit <- tryCatch(lm(value ~ donor_status + islet_diam_um, data = rdf), error = function(e) NULL)
      p_global <- NA_real_
      if (!is.null(fit)) {
        at <- tryCatch(anova(fit), error = function(e) NULL)
        if (!is.null(at)) p_global <- suppressWarnings(as.numeric(at[["Pr(>F)"]][1]))
      }
      res <- tryCatch({
        fit_res <- lm(value ~ islet_diam_um, data = rdf)
        r <- resid(fit_res)
        pt <- pairwise.t.test(r, rdf$donor_status, p.adjust.method = "BH")
        mat <- as.data.frame(as.table(pt$p.value))
        colnames(mat) <- c("group1","group2","p_value")
        mat <- mat[!is.na(mat$p_value), , drop = FALSE]
        mat
      }, error = function(e) NULL)
      global <- data.frame(type = "global", contrast = "donor_status", p_value = p_global, stringsAsFactors = FALSE)
      if (!is.null(res) && nrow(res)) {
        pairs <- data.frame(type = "pairwise",
                            contrast = paste(res$group1, "vs", res$group2),
                            p_value = res$p_value,
                            stringsAsFactors = FALSE)
        out <- rbind(global, pairs)
      } else {
        out <- global
      }
      out
    }, ignoreInit = TRUE)

    stats_data <- reactive({
      st <- stats_run()
      if (is.null(st) || nrow(st) == 0) return(NULL)
      alpha_num <- as.numeric(input$alpha)
      st$p_adj <- NA_real_
      if (any(st$type == "pairwise")) {
        idx <- which(st$type == "pairwise")
        st$p_adj[idx] <- p.adjust(st$p_value[idx], method = "BH")
      }
      st$sig <- ifelse(!is.na(st$p_adj), st$p_adj <= alpha_num, st$p_value <= alpha_num)
      st
    })

    output$stats_explanation <- renderUI({
      st <- stats_data()
      if (is.null(st) || nrow(st) == 0) {
        return(tags$div(
          style = "padding: 15px; color: #666;",
          tags$p("Click 'Run Statistics' to perform statistical analysis.")
        ))
      }
      alpha_num <- as.numeric(input$alpha)
      n_total <- nrow(st)
      n_sig <- sum(st$sig, na.rm = TRUE)
      tags$div(
        style = "padding: 15px;",
        tags$h5("Analysis Method:", style = "margin-top: 0;"),
        tags$p(
          style = "margin-bottom: 10px;",
          tags$strong("Global ANOVA:"), " Tests if there are significant differences in the selected metric across the three donor status groups (ND, Aab+, T1D) while controlling for islet diameter.",
          tags$br(),
          tags$em("Model: value ~ donor_status + islet_diam_um")
        ),
        tags$p(
          style = "margin-bottom: 10px;",
          tags$strong("Pairwise t-tests:"), " Compares each pair of donor status groups (ND vs Aab+, ND vs T1D, Aab+ vs T1D) on residuals from the diameter-adjusted model.",
          tags$br(),
          tags$em("Multiple testing correction: Benjamini-Hochberg (FDR) method")
        ),
        tags$hr(),
        tags$h5("Results Summary:"),
        tags$p(
          sprintf("Significance threshold: \u03b1 = %s", input$alpha),
          tags$br(),
          sprintf("Significant tests: %d / %d (%.1f%%)", n_sig, n_total, 100 * n_sig / n_total)
        ),
        tags$p(
          style = "font-size: 90%; color: #666;",
          tags$strong("Interpretation:"),
          " A significant result (p < \u03b1) suggests the selected metric differs between donor status groups after accounting for islet size. ",
          "Non-significant results may indicate no biological difference or insufficient statistical power."
        )
      )
    })

    output$stats_tbl <- renderTable({
      st <- stats_data()
      if (is.null(st) || nrow(st) == 0) return(NULL)
      fmt <- function(x) ifelse(is.na(x), NA_character_, formatC(x, format = "e", digits = 2))
      display <- st
      display$p_value <- fmt(display$p_value)
      display$p_adj <- fmt(display$p_adj)
      display
    }, striped = TRUE, bordered = TRUE, spacing = "s")

    output$stats_plot <- renderPlotly({
      st <- stats_data()
      req(st, nrow(st) > 0)
      alpha_num <- as.numeric(input$alpha)
      req(alpha_num)
      df <- st %>%
        mutate(neglog = ifelse(!is.na(p_value) & p_value > 0, -log10(p_value), NA_real_)) %>%
        filter(!is.na(neglog))
      if (nrow(df) == 0) return(NULL)
      thresh <- -log10(alpha_num)
      g <- ggplot(df, aes(x = reorder(contrast, neglog), y = neglog, fill = type)) +
        geom_col(show.legend = TRUE) +
        geom_hline(yintercept = thresh, linetype = "dashed", color = "#444444") +
        coord_flip() +
        labs(x = "Contrast", y = expression(-log[10](p)), fill = "Test type",
             title = "Test significance overview") +
        theme_minimal(base_size = 14)
      ggplotly(g)
    })

    output$pairwise_plot <- renderPlotly({
      st <- stats_data()
      req(st, nrow(st) > 0)
      alpha_num <- as.numeric(input$alpha)
      req(alpha_num)
      df <- st %>% filter(type == "pairwise" & !is.na(p_value))
      if (nrow(df) == 0) return(NULL)
      df <- df %>% mutate(sig = ifelse(is.na(p_adj), p_value <= alpha_num, p_adj <= alpha_num))
      g <- ggplot(df, aes(x = reorder(contrast, p_value), y = p_value, color = sig)) +
        geom_point(size = 3) +
        geom_hline(yintercept = alpha_num, linetype = "dashed", color = "#444444") +
        scale_y_log10() +
        scale_color_manual(values = c(`FALSE` = "#1f77b4", `TRUE` = "#d62728"),
                           labels = c(`FALSE` = "NS", `TRUE` = "Significant")) +
        coord_flip() +
        labs(x = "Pairwise contrast", y = "p-value (log scale)", color = "", title = "Pairwise comparison p-values") +
        theme_minimal(base_size = 14)
      ggplotly(g)
    })

    output$auc_plot <- renderPlotly({
      sdf <- summary_df()
      if (is.null(sdf) || nrow(sdf) == 0) return(NULL)
      auc_by_group <- tryCatch({
        sdf %>%
          dplyr::group_by(donor_status) %>%
          dplyr::arrange(diam_mid) %>%
          dplyr::summarise(
            auc = if (n() < 2) 0 else sum(diff(diam_mid) * (head(y, -1) + tail(y, -1))) / 2,
            .groups = "drop"
          )
      }, error = function(e) NULL)
      if (is.null(auc_by_group) || nrow(auc_by_group) == 0) return(NULL)
      g <- ggplot(auc_by_group, aes(x = donor_status, y = auc, fill = donor_status)) +
        geom_col(width = 0.7, color = "black", alpha = 0.8) +
        scale_fill_manual(values = c("ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd")) +
        labs(x = "Donor Status", y = "Area Under Curve (AUC)", title = "Integrated Area by Donor Group") +
        theme_minimal(base_size = 14) +
        theme(legend.position = "none")
      ggplotly(g)
    })

    output$auc_table <- renderTable({
      sdf <- summary_df()
      if (is.null(sdf) || nrow(sdf) == 0) return(NULL)
      auc_by_group <- tryCatch({
        result <- sdf %>%
          dplyr::group_by(donor_status) %>%
          dplyr::arrange(diam_mid) %>%
          dplyr::summarise(
            auc = if (n() < 2) 0 else sum(diff(diam_mid) * (head(y, -1) + tail(y, -1))) / 2,
            n_bins = n(),
            .groups = "drop"
          )
        if (nrow(result) > 0) result else NULL
      }, error = function(e) NULL)
      if (is.null(auc_by_group)) {
        return(data.frame(Message = "Insufficient data for AUC calculation"))
      }
      nd_auc <- auc_by_group$auc[auc_by_group$donor_status == "ND"]
      t1d_auc <- auc_by_group$auc[auc_by_group$donor_status == "T1D"]
      result <- data.frame(
        "Donor Group" = auc_by_group$donor_status,
        "AUC" = sprintf("%.2f", auc_by_group$auc),
        "N Bins" = auc_by_group$n_bins,
        check.names = FALSE
      )
      if (length(nd_auc) > 0 && length(t1d_auc) > 0) {
        pct_change_nd_t1d <- ((t1d_auc - nd_auc) / nd_auc) * 100
        result <- rbind(result, data.frame(
          "Donor Group" = "T1D vs ND",
          "AUC" = sprintf("%.1f%%", pct_change_nd_t1d),
          "N Bins" = "",
          check.names = FALSE
        ))
      }
      result
    }, striped = TRUE, bordered = TRUE)

    output$download_stats <- downloadHandler(
      filename = function() {
        paste0("statistics_", gsub("[^0-9A-Za-z]+","_", Sys.time()), ".csv")
      },
      content = function(file) {
        st <- stats_data()
        if (is.null(st) || nrow(st) == 0) st <- data.frame()
        descriptions <- get_selection_description()
        writeLines(paste("#", descriptions), file)
        writeLines("", file, sep = "\n")
        write.table(st, file, append = TRUE, sep = ",", row.names = FALSE)
      }
    )
  })
}
