# ---------- Statistics module server ----------
# Exports: statistics_server(id, raw_df, summary_df, get_selection_description)
#
# Phase 16: Donor-level aggregation to fix pseudoreplication.
# Tests use donor-level means (N~15) instead of islet-level (N~5,214).
# Mixed-effects model (lmerTest) as sensitivity analysis.
#
# Dependencies (auto-sourced):
#   summary_stats, per_bin_anova, per_bin_kendall   -- utils_stats.R
#   aggregate_to_donor, lmer_test_donor              -- utils_stats.R
#   per_bin_donor_anova, per_bin_donor_kendall       -- utils_stats.R
#   normality_tests                                  -- utils_stats.R
#   cohens_d, eta_squared, pairwise_wilcox           -- utils_stats.R
#   bin_islet_sizes                                  -- data_loading.R

statistics_server <- function(id, raw_df, summary_df, get_selection_description,
                              donor_colors = reactive(DONOR_COLORS)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ===========================================================================
    # Core computation: eventReactive triggered by Run button
    # ===========================================================================

    stats_run <- eventReactive(input$run_tests, {
      rdf <- raw_df()
      if (is.null(rdf) || !nrow(rdf)) return(NULL)

      alpha_num <- as.numeric(input$alpha)
      is_parametric <- identical(input$test_type, "Parametric")

      # ---- Outlier removal (per donor group, >3 SD) ----
      if (isTRUE(input$stats_remove_outliers)) {
        keep <- rep(TRUE, nrow(rdf))
        for (grp in unique(rdf$donor_status)) {
          idx <- which(rdf$donor_status == grp)
          if (length(idx) < 3) next
          grp_mean <- mean(rdf$value[idx], na.rm = TRUE)
          grp_sd   <- sd(rdf$value[idx], na.rm = TRUE)
          if (!is.finite(grp_sd) || grp_sd == 0) next
          keep[idx] <- is.na(rdf$value[idx]) | abs(rdf$value[idx] - grp_mean) <= (3 * grp_sd)
        }
        rdf <- rdf[keep, , drop = FALSE]
      }

      # ---- Min-cells filter (like spatial tab) ----
      min_cells <- input$stats_min_cells %||% 1
      if (min_cells > 1 &&
          "total_cells_core" %in% colnames(rdf) &&
          "total_cells_peri" %in% colnames(rdf)) {
        tc_core <- suppressWarnings(as.numeric(rdf$total_cells_core))
        tc_peri <- suppressWarnings(as.numeric(rdf$total_cells_peri))
        tc_core <- ifelse(is.finite(tc_core), tc_core, 0)
        tc_peri <- ifelse(is.finite(tc_peri), tc_peri, 0)
        rdf <- rdf[is.finite(tc_core + tc_peri) & (tc_core + tc_peri) >= min_cells, , drop = FALSE]
      }

      # ---- Diameter range filter (from Stats inline controls) ----
      if (!is.null(input$stats_diam_range) && length(input$stats_diam_range) == 2) {
        dmin <- as.numeric(input$stats_diam_range[1])
        dmax <- as.numeric(input$stats_diam_range[2])
        if (is.finite(dmin) && is.finite(dmax)) {
          rdf <- rdf[is.finite(rdf$islet_diam_um) &
                      rdf$islet_diam_um >= dmin &
                      rdf$islet_diam_um <= dmax, , drop = FALSE]
        }
      }

      if (nrow(rdf) == 0) return(NULL)

      rdf$donor_status <- factor(rdf$donor_status, levels = c("ND", "Aab+", "T1D"))
      n_islets_filtered <- nrow(rdf)

      # Detect feature name from raw_df attributes
      feature_name <- attr(raw_df(), "selection_used") %||% "Unknown"

      # ---- Donor-level aggregation ----
      ddf <- aggregate_to_donor(rdf, agg_fn = "mean")
      ddf$donor_status <- factor(ddf$donor_status, levels = c("ND", "Aab+", "T1D"))
      n_donors <- nrow(ddf)

      # ---- Normality test on donor-level data ----
      norm_result <- normality_tests(ddf)

      # ---- Mixed-effects sensitivity (islet-level data with donor random effect) ----
      lmer_result <- tryCatch(lmer_test_donor(rdf), error = function(e) NULL)
      lmer_p <- if (!is.null(lmer_result)) lmer_result$p_value else NA_real_
      icc <- if (!is.null(lmer_result)) lmer_result$icc else NA_real_

      # ---- Global test on DONOR-LEVEL means ----
      global_test <- if (is_parametric) "ANOVA" else "Kruskal-Wallis"
      global_p <- NA_real_
      eta_sq <- NA_real_

      if (is_parametric) {
        fit <- tryCatch(lm(value ~ donor_status, data = ddf), error = function(e) NULL)
        if (!is.null(fit)) {
          at <- tryCatch(anova(fit), error = function(e) NULL)
          if (!is.null(at)) global_p <- suppressWarnings(as.numeric(at[["Pr(>F)"]][1]))
          eta_sq <- eta_squared(fit)
        }
      } else {
        kt <- tryCatch(kruskal.test(value ~ donor_status, data = ddf), error = function(e) NULL)
        if (!is.null(kt)) global_p <- kt$p.value
        # Eta-squared approximation for K-W: H / (N-1)
        if (!is.null(kt) && n_donors > 1) eta_sq <- kt$statistic / (n_donors - 1)
      }

      # ---- Pairwise comparisons with effect sizes (donor-level) ----
      groups <- levels(ddf$donor_status)
      pairs <- combn(groups, 2, simplify = FALSE)

      pairwise_rows <- lapply(pairs, function(pair) {
        g1_vals <- ddf$value[ddf$donor_status == pair[1]]
        g2_vals <- ddf$value[ddf$donor_status == pair[2]]
        contrast_label <- paste(pair[1], "vs", pair[2])

        if (is_parametric) {
          tt <- tryCatch(t.test(g1_vals, g2_vals), error = function(e) NULL)
          p_val <- if (!is.null(tt)) tt$p.value else NA_real_
        } else {
          wt <- tryCatch(wilcox.test(g1_vals, g2_vals, exact = FALSE), error = function(e) NULL)
          p_val <- if (!is.null(wt)) wt$p.value else NA_real_
        }

        cd <- cohens_d(g1_vals, g2_vals)

        data.frame(
          contrast = contrast_label,
          p_value = p_val,
          cohens_d = cd$d,
          d_ci_lo = cd$ci_lo,
          d_ci_hi = cd$ci_hi,
          n1 = length(g1_vals[is.finite(g1_vals)]),
          n2 = length(g2_vals[is.finite(g2_vals)]),
          stringsAsFactors = FALSE
        )
      })

      pairwise_df <- do.call(rbind, pairwise_rows)
      if (!is.null(pairwise_df) && nrow(pairwise_df) > 0) {
        pairwise_df$p_adj <- p.adjust(pairwise_df$p_value, method = "BH")
        pairwise_df$sig <- ifelse(!is.na(pairwise_df$p_adj),
                                  pairwise_df$p_adj < alpha_num,
                                  FALSE)
      }

      # ---- Per-bin analysis (conditional on no-binning checkbox) ----
      bin_anova <- NULL
      bin_kendall <- NULL

      if (!isTRUE(input$stats_no_binning)) {
        bw <- suppressWarnings(as.numeric(input$bin_width))
        if (!is.finite(bw) || bw <= 0) bw <- 50

        rdf_binned <- bin_islet_sizes(rdf, "islet_diam_um", bw)

        bin_anova <- per_bin_donor_anova(rdf_binned, "diam_bin", "donor_status", "value")
        bin_kendall <- per_bin_donor_kendall(rdf_binned, "diam_bin", "donor_status", "value")

        # BH correction across bins to control false discovery rate
        if (!is.null(bin_anova) && nrow(bin_anova) > 0) {
          bin_anova$p_adj <- p.adjust(bin_anova$p_anova, method = "BH")
        }
        if (!is.null(bin_kendall) && nrow(bin_kendall) > 0) {
          bin_kendall$p_adj <- p.adjust(bin_kendall$p_kendall, method = "BH")
        }
      }

      # ---- Demographics (NULL if columns absent) ----
      demo_summary <- NULL
      age_corr <- NULL
      age_model_p <- NULL
      gender_strat <- NULL

      has_age <- "age" %in% colnames(ddf)
      has_gender <- "gender" %in% colnames(ddf)

      if (has_age || has_gender) {
        # Demographic summary per donor group (already at donor level)
        demo_rows <- lapply(levels(ddf$donor_status), function(g) {
          sub <- ddf[ddf$donor_status == g, , drop = FALSE]
          n_d <- nrow(sub)
          islets_total <- sum(sub$n_islets, na.rm = TRUE)

          age_med <- NA_real_; age_range_str <- ""
          if (has_age) {
            ages <- as.numeric(sub$age)
            ages <- ages[is.finite(ages)]
            if (length(ages) > 0) {
              age_med <- median(ages)
              age_range_str <- paste0(min(ages), "-", max(ages))
            }
          }

          pct_m <- NA_real_; pct_f <- NA_real_
          if (has_gender) {
            genders <- as.character(sub$gender)
            genders <- genders[!is.na(genders) & nzchar(genders)]
            if (length(genders) > 0) {
              pct_m <- round(100 * sum(genders == "M") / length(genders), 1)
              pct_f <- round(100 * sum(genders == "F") / length(genders), 1)
            }
          }

          data.frame(
            donor_status = g, n = n_d, n_islets = islets_total,
            age_median = age_med, age_range = age_range_str,
            pct_male = pct_m, pct_female = pct_f,
            stringsAsFactors = FALSE
          )
        })
        demo_summary <- do.call(rbind, demo_rows)

        # Age-feature correlation across ALL donors (not per-group with n=5)
        if (has_age) {
          a <- as.numeric(ddf$age)
          v <- ddf$value
          keep <- is.finite(a) & is.finite(v)
          age_corr <- if (sum(keep) >= 3) {
            ct <- tryCatch(cor.test(a[keep], v[keep], method = "pearson"), error = function(e) NULL)
            if (!is.null(ct)) {
              data.frame(cor = unname(ct$estimate), p_value = ct$p.value,
                         n = sum(keep), stringsAsFactors = FALSE)
            }
          }
        }

        # Covariate-adjusted model on DONOR-level data
        if (has_age && has_gender) {
          ddf_cov <- ddf
          ddf_cov$age_num <- as.numeric(ddf_cov$age)
          ddf_cov$gender_f <- factor(ddf_cov$gender)
          cov_fit <- tryCatch(lm(value ~ donor_status + age_num + gender_f, data = ddf_cov), error = function(e) NULL)
          if (!is.null(cov_fit)) {
            at <- tryCatch(anova(cov_fit), error = function(e) NULL)
            if (!is.null(at) && "donor_status" %in% rownames(at)) {
              age_model_p <- as.numeric(at["donor_status", "Pr(>F)"])
            }
          }
        }

        # Gender-stratified tests on donor-level data
        if (has_gender) {
          genders_present <- unique(na.omit(as.character(ddf$gender)))
          gender_rows <- lapply(genders_present, function(gen) {
            sub <- ddf[as.character(ddf$gender) == gen, , drop = FALSE]
            if (n_distinct(sub$donor_status) < 2) return(NULL)
            n_per_grp <- table(sub$donor_status)
            low_power <- any(n_per_grp[n_per_grp > 0] < 3)
            pw <- if (is_parametric) {
              tryCatch({
                pt <- pairwise.t.test(sub$value, sub$donor_status, p.adjust.method = "BH")
                mat <- as.data.frame(as.table(pt$p.value))
                colnames(mat) <- c("group1", "group2", "p_value")
                mat[!is.na(mat$p_value), , drop = FALSE]
              }, error = function(e) NULL)
            } else {
              pairwise_wilcox(sub, "donor_status", "value")
            }
            if (is.null(pw) || nrow(pw) == 0) return(NULL)
            pw$gender <- gen
            pw$contrast <- paste(pw$group1, "vs", pw$group2)
            pw$low_power <- low_power
            pw
          })
          gender_strat <- do.call(rbind, gender_rows)
        }
      }

      # ---- AUC: trapezoidal integration of binned summary curve per donor group ----
      sdf <- summary_df()
      auc_by_group <- tryCatch({
        if (is.null(sdf) || nrow(sdf) == 0) NULL
        else {
          sdf %>%
            dplyr::group_by(donor_status) %>%
            dplyr::arrange(diam_mid) %>%
            dplyr::summarise(
              auc = if (n() < 2) 0 else sum(diff(diam_mid) * (head(y, -1) + tail(y, -1))) / 2,
              n_bins = n(),
              .groups = "drop"
            )
        }
      }, error = function(e) NULL)

      list(
        n_islets         = n_islets_filtered,
        n_donors         = n_donors,
        donor_df         = ddf,
        feature_name     = feature_name,
        global_test      = global_test,
        global_p         = global_p,
        eta_sq           = eta_sq,
        alpha            = alpha_num,
        is_parametric    = is_parametric,
        pairwise_df      = pairwise_df,
        bin_anova        = bin_anova,
        bin_kendall      = bin_kendall,
        auc_by_group     = auc_by_group,
        demo_summary     = demo_summary,
        age_corr         = age_corr,
        age_model_p      = age_model_p,
        gender_strat     = gender_strat,
        lmer_p           = lmer_p,
        icc              = icc,
        norm_result      = norm_result,
        rdf              = rdf
      )
    }, ignoreInit = TRUE)

    # ===========================================================================
    # Card 1: Overview Banner
    # ===========================================================================

    output$overview_banner <- renderUI({
      st <- stats_run()
      if (is.null(st)) {
        return(tags$p(style = "color: #666; margin: 0;",
                      "Select a feature in the sidebar, then click 'Run Statistics'."))
      }

      p_col <- if (is.na(st$global_p)) "#999" else if (st$global_p < st$alpha) "#d62728" else "#2e7d32"
      eta_label <- if (is.na(st$eta_sq)) "N/A" else sprintf("%.3f", st$eta_sq)
      icc_label <- if (is.na(st$icc)) "" else sprintf(" | ICC = %.2f", st$icc)

      tags$div(style = "display: flex; gap: 15px; align-items: center; flex-wrap: wrap;",
        tags$span(style = "background: #e3f2fd; padding: 4px 10px; border-radius: 4px; font-weight: 600;",
                  st$feature_name),
        tags$span(style = "padding: 4px 10px;",
                  paste0("N = ", st$n_donors, " donors (",
                         formatC(st$n_islets, format = "d", big.mark = ","), " islets)",
                         icc_label)),
        tags$span(style = paste0("background: ", p_col, "; color: white; padding: 4px 10px; border-radius: 4px;"),
                  paste0(st$global_test, " p = ", if (is.na(st$global_p)) "N/A" else formatC(st$global_p, format = "e", digits = 2))),
        tags$span(style = "padding: 4px 10px;",
                  paste0("\u03b7\u00b2 = ", eta_label))
      )
    })

    # ===========================================================================
    # Normality test results
    # ===========================================================================

    output$normality_result <- renderUI({
      st <- stats_run()
      if (is.null(st) || is.null(st$norm_result)) return(NULL)

      nr <- st$norm_result
      items <- lapply(seq_len(nrow(nr)), function(i) {
        g <- nr$donor_status[i]
        w_val <- nr$W[i]
        p_val <- nr$p_value[i]
        n_val <- nr$n[i]
        is_norm <- nr$is_normal[i]
        if (is.na(is_norm)) {
          col <- "#999"
          label <- paste0(g, ": n=", n_val, " (too few)")
        } else if (is_norm) {
          col <- "#2e7d32"
          label <- paste0(g, ": W=", sprintf("%.3f", w_val),
                          ", p=", sprintf("%.3f", p_val), " (normal)")
        } else {
          col <- "#e65100"
          label <- paste0(g, ": W=", sprintf("%.3f", w_val),
                          ", p=", sprintf("%.3f", p_val), " (non-normal)")
        }
        tags$span(style = paste0("color: ", col, "; margin-right: 12px;"), label)
      })

      all_normal <- all(nr$is_normal, na.rm = TRUE)
      any_non_normal <- any(!nr$is_normal, na.rm = TRUE)
      suggestion <- if (any_non_normal) {
        tags$span(style = "color: #e65100; font-weight: 600;",
                  "Consider non-parametric tests.")
      } else if (all_normal) {
        tags$span(style = "color: #2e7d32;",
                  "Parametric tests may be appropriate.")
      } else {
        NULL
      }

      tags$div(
        style = "background: #f5f5f5; padding: 8px 12px; border-radius: 5px; font-size: 12px; line-height: 1.8;",
        tags$strong("Shapiro-Wilk normality (donor means):"),
        tags$br(),
        tagList(items),
        if (!is.null(suggestion)) tagList(tags$br(), suggestion),
        tags$br(),
        tags$span(style = "color: #999; font-size: 11px;",
                  "Note: with n=5/group, Shapiro-Wilk has limited power to detect departures from normality.")
      )
    })

    # ===========================================================================
    # Pseudoreplication info banner
    # ===========================================================================

    output$pseudorep_banner <- renderUI({
      st <- stats_run()
      if (is.null(st)) return(NULL)

      icc_str <- if (is.na(st$icc)) "N/A" else sprintf("%.2f", st$icc)
      lmer_str <- if (is.na(st$lmer_p)) "did not converge" else
        paste0("p = ", formatC(st$lmer_p, format = "e", digits = 2))

      tags$div(
        style = "background: #e8f4fd; padding: 10px 14px; border-radius: 5px; margin-bottom: 12px; font-size: 13px; line-height: 1.5;",
        tags$strong("Pseudoreplication correction: "),
        paste0("Tests use donor-level means (N = ", st$n_donors,
               ") to avoid treating correlated islets as independent. "),
        "Islets within a donor share biology and tissue processing. ",
        paste0("Mixed-effects sensitivity: ", lmer_str, ". "),
        paste0("ICC = ", icc_str,
               if (!is.na(st$icc) && st$icc > 0.1) " (substantial donor-level clustering)."
               else if (!is.na(st$icc)) " (low donor-level clustering)."
               else ".")
      )
    })

    # ===========================================================================
    # Card 2: Hypothesis Testing — table + forest plot
    # ===========================================================================

    output$test_results_table <- renderTable({
      st <- stats_run()
      if (is.null(st)) return(NULL)

      fmt_p <- function(x) ifelse(is.na(x), "N/A", formatC(x, format = "e", digits = 2))
      fmt_d <- function(x) ifelse(is.na(x), "N/A", sprintf("%.2f", x))

      # Global row (donor-level)
      global_row <- data.frame(
        Contrast = paste0("Global (", st$global_test, ", N=", st$n_donors, ")"),
        `p-value` = fmt_p(st$global_p),
        `p (adj)` = "-",
        `Cohen's d` = "-",
        `95% CI` = "-",
        Sig = if (!is.na(st$global_p) && st$global_p < st$alpha) "*" else "",
        check.names = FALSE, stringsAsFactors = FALSE
      )

      pw <- st$pairwise_df
      rows <- global_row
      if (!is.null(pw) && nrow(pw) > 0) {
        test_label <- if (st$is_parametric) "t-test" else "Wilcoxon"
        pw_rows <- data.frame(
          Contrast = paste0(pw$contrast, " (", test_label, ")"),
          `p-value` = fmt_p(pw$p_value),
          `p (adj)` = fmt_p(pw$p_adj),
          `Cohen's d` = fmt_d(pw$cohens_d),
          `95% CI` = ifelse(is.na(pw$d_ci_lo), "N/A",
                            paste0("[", sprintf("%.2f", pw$d_ci_lo), ", ", sprintf("%.2f", pw$d_ci_hi), "]")),
          Sig = ifelse(!is.na(pw$p_adj) & pw$p_adj < st$alpha, "*", ""),
          check.names = FALSE, stringsAsFactors = FALSE
        )
        rows <- rbind(rows, pw_rows)
      }

      # Mixed-effects row
      lmer_row <- data.frame(
        Contrast = "Mixed-effects (lmer)",
        `p-value` = fmt_p(st$lmer_p),
        `p (adj)` = "-",
        `Cohen's d` = "-",
        `95% CI` = "-",
        Sig = if (!is.na(st$lmer_p) && st$lmer_p < st$alpha) "*" else "",
        check.names = FALSE, stringsAsFactors = FALSE
      )
      rbind(rows, lmer_row)
    }, striped = TRUE, bordered = TRUE, spacing = "s", align = "lccccc")

    output$forest_plot <- renderPlotly({
      st <- stats_run()
      if (is.null(st)) return(NULL)
      pw <- st$pairwise_df
      if (is.null(pw) || nrow(pw) == 0 || all(is.na(pw$cohens_d))) return(NULL)

      pw <- pw[!is.na(pw$cohens_d), , drop = FALSE]
      pw$contrast <- factor(pw$contrast, levels = rev(pw$contrast))
      pw$sig_label <- ifelse(pw$sig, "Significant", "NS")

      g <- ggplot(pw, aes(x = cohens_d, y = contrast, color = sig_label)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "#999") +
        geom_point(size = 3) +
        geom_errorbarh(aes(xmin = d_ci_lo, xmax = d_ci_hi), height = 0.2) +
        scale_color_manual(values = c("Significant" = "#d62728", "NS" = "#1f77b4"),
                           name = "") +
        labs(x = "Cohen's d (donor-level)", y = NULL, title = "Pairwise Effect Sizes") +
        theme_minimal(base_size = 13) +
        theme(legend.position = "none")

      ggplotly(g, tooltip = c("x", "y", "colour")) %>%
        layout(margin = list(b = 60),
               legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.25))
    })

    # ===========================================================================
    # Card 3: Per-Bin Significance Heatmap
    # ===========================================================================

    output$bin_heatmap <- renderPlotly({
      st <- stats_run()
      if (is.null(st)) return(NULL)

      ba <- st$bin_anova
      bk <- st$bin_kendall

      if ((is.null(ba) || nrow(ba) == 0) && (is.null(bk) || nrow(bk) == 0)) return(NULL)

      # Collect all bin midpoints
      all_mids <- sort(unique(c(ba$mid, bk$mid)))
      if (length(all_mids) == 0) return(NULL)

      # Build z matrix using BH-corrected p-values
      anova_p <- rep(NA_real_, length(all_mids))
      kendall_p <- rep(NA_real_, length(all_mids))

      if (!is.null(ba) && nrow(ba) > 0) {
        idx <- match(ba$mid, all_mids)
        anova_p[idx[!is.na(idx)]] <- ba$p_adj[!is.na(idx)]
      }
      if (!is.null(bk) && nrow(bk) > 0) {
        idx <- match(bk$mid, all_mids)
        kendall_p[idx[!is.na(idx)]] <- bk$p_adj[!is.na(idx)]
      }

      z_anova <- ifelse(!is.na(anova_p) & anova_p > 0, -log10(anova_p), NA_real_)
      z_kendall <- ifelse(!is.na(kendall_p) & kendall_p > 0, -log10(kendall_p), NA_real_)

      z_mat <- rbind(z_anova, z_kendall)

      # Star annotations (based on corrected p-values)
      star_fn <- function(p) {
        ifelse(is.na(p), "",
               ifelse(p < 0.001, "***",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.05, "*", ""))))
      }
      star_anova <- star_fn(anova_p)
      star_kendall <- star_fn(kendall_p)
      text_mat <- rbind(star_anova, star_kendall)

      mid_labels <- as.character(round(all_mids))

      plot_ly(
        x = mid_labels, y = c("ANOVA", "Kendall"), z = z_mat,
        type = "heatmap",
        colorscale = list(c(0, "#f7fbff"), c(0.5, "#6baed6"), c(1, "#08306b")),
        colorbar = list(title = "-log10(q)"),
        text = text_mat, texttemplate = "%{text}",
        hovertemplate = "Bin: %{x} \u00b5m<br>Test: %{y}<br>-log10(q): %{z:.2f}<extra></extra>"
      ) %>%
        layout(
          xaxis = list(title = "Diameter bin midpoint (\u00b5m)", tickangle = -45),
          yaxis = list(title = ""),
          title = "Stratified Significance (BH-corrected, donor-level)"
        )
    })

    # ===========================================================================
    # Card 4: Trend Analysis
    # ===========================================================================

    output$trend_plot <- renderPlotly({
      st <- stats_run()
      if (is.null(st)) return(NULL)
      bk <- st$bin_kendall
      if (is.null(bk) || nrow(bk) == 0) return(NULL)

      bk$sig_label <- ifelse(!is.na(bk$p_adj) & bk$p_adj < st$alpha, "Significant", "NS")

      # Overall Kendall tau on donor-level means
      ddf <- st$donor_df
      code_group <- function(x) {
        x <- as.character(x)
        ifelse(x == "ND", 0, ifelse(x == "Aab+", 1, ifelse(x == "T1D", 2, NA_real_)))
      }
      x_all <- code_group(ddf$donor_status)
      y_all <- ddf$value
      keep <- is.finite(x_all) & is.finite(y_all)
      overall_tau <- NA_real_
      if (sum(keep) >= 3) {
        ct <- tryCatch(cor.test(x_all[keep], y_all[keep], method = "kendall", exact = FALSE), error = function(e) NULL)
        if (!is.null(ct)) overall_tau <- unname(ct$estimate)
      }

      direction_label <- if (is.na(overall_tau)) "No clear trend" else
        if (overall_tau > 0) "\u2191 Increases with disease" else "\u2193 Decreases with disease"

      g <- ggplot(bk, aes(x = mid, y = tau, color = sig_label)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "#999") +
        geom_line(color = "#333", alpha = 0.4) +
        geom_point(size = 3) +
        scale_color_manual(values = c("Significant" = "#d62728", "NS" = "#aaaaaa"), name = "") +
        labs(x = "Diameter bin midpoint (\u00b5m)", y = "Kendall \u03c4",
             title = paste("Disease Progression Trend:", direction_label)) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "none")

      ggplotly(g, tooltip = c("x", "y", "colour")) %>%
        layout(margin = list(b = 60),
               legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.25))
    })

    # ===========================================================================
    # Card 5: Demographics (conditional — NULL when no age/gender)
    # ===========================================================================

    output$demographics_card <- renderUI({
      st <- stats_run()
      # Show placeholder before first run
      if (is.null(st)) {
        rdf <- raw_df()
        if (is.null(rdf) || !("age" %in% colnames(rdf) || "gender" %in% colnames(rdf))) return(NULL)
        return(div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
                   h5("Demographics Analysis"),
                   tags$p(style = "color: #666;", "Run statistics to see demographic analysis.")))
      }
      if (is.null(st$demo_summary)) return(NULL)

      div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
        h5("Demographics Analysis"),
        h6("Donor Summary", style = "color: #555;"),
        tableOutput(ns("demo_summary_table")),
        hr(style = "margin: 8px 0;"),
        h6("Age vs Feature (donor-level)", style = "color: #555;"),
        plotlyOutput(ns("age_scatter"), height = "280px"),
        hr(style = "margin: 8px 0;"),
        h6("Gender-Stratified Tests", style = "color: #555;"),
        tableOutput(ns("gender_strat_table")),
        uiOutput(ns("covariate_results"))
      )
    })

    output$demo_summary_table <- renderTable({
      st <- stats_run()
      if (is.null(st) || is.null(st$demo_summary)) return(NULL)
      ds <- st$demo_summary
      display <- data.frame(
        `Donor Status` = ds$donor_status,
        `N Donors` = ds$n,
        `N Islets` = ifelse(is.na(ds$n_islets), "-", formatC(ds$n_islets, format = "d", big.mark = ",")),
        `Age Median` = ifelse(is.na(ds$age_median), "-", sprintf("%.0f", ds$age_median)),
        `Age Range` = ifelse(nzchar(ds$age_range), ds$age_range, "-"),
        `% Male` = ifelse(is.na(ds$pct_male), "-", sprintf("%.0f%%", ds$pct_male)),
        `% Female` = ifelse(is.na(ds$pct_female), "-", sprintf("%.0f%%", ds$pct_female)),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      display
    }, striped = TRUE, bordered = TRUE, spacing = "s")

    output$age_scatter <- renderPlotly({
      st <- stats_run()
      if (is.null(st) || is.null(st$age_corr)) return(NULL)
      ddf <- st$donor_df
      if (!("age" %in% colnames(ddf))) return(NULL)

      ddf$age_num <- as.numeric(ddf$age)
      ddf <- ddf[is.finite(ddf$age_num) & is.finite(ddf$value), , drop = FALSE]
      if (nrow(ddf) == 0) return(NULL)

      ddf$donor_status <- factor(ddf$donor_status, levels = c("ND", "Aab+", "T1D"))

      # Annotation for overall correlation
      corr_label <- sprintf("r = %.2f, p = %.3f (n=%d)",
                            st$age_corr$cor, st$age_corr$p_value, st$age_corr$n)

      g <- ggplot(ddf, aes(x = age_num, y = value, color = donor_status)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, aes(group = 1), color = "#333") +
        scale_color_manual(values = donor_colors()) +
        labs(x = "Donor Age (years)", y = st$feature_name, color = NULL,
             title = "Age-Feature Relationship (donor means)",
             subtitle = corr_label) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")

      ggplotly(g, tooltip = c("x", "y", "colour")) %>%
        layout(legend = list(orientation = "h", x = 0.2, y = -0.2))
    })

    output$gender_strat_table <- renderTable({
      st <- stats_run()
      if (is.null(st) || is.null(st$gender_strat) || nrow(st$gender_strat) == 0) return(NULL)
      gs <- st$gender_strat
      result <- data.frame(
        Gender = gs$gender,
        Contrast = gs$contrast,
        `p-value` = formatC(gs$p_value, format = "e", digits = 2),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      if ("low_power" %in% colnames(gs)) {
        result$Note <- ifelse(gs$low_power, "Low power (n<3/group)", "")
      }
      result
    }, striped = TRUE, bordered = TRUE, spacing = "s")

    output$covariate_results <- renderUI({
      st <- stats_run()
      if (is.null(st) || is.null(st$age_model_p)) return(NULL)
      p_cov <- st$age_model_p
      survived <- !is.na(p_cov) && p_cov < st$alpha
      tags$div(style = "margin-top: 10px; padding: 10px; background: #f8f9fa; border-radius: 5px;",
        tags$strong("Covariate-Adjusted Model (donor-level)"),
        tags$p(style = "margin: 5px 0 0 0;",
          paste0("lm(value ~ donor_status + age + gender, data = donor_means): donor_status p = ",
                 if (is.na(p_cov)) "N/A" else formatC(p_cov, format = "e", digits = 2)),
          tags$br(),
          if (survived) {
            tags$span(style = "color: #d62728; font-weight: 600;",
                      "Donor status effect SURVIVES adjustment for age and gender.")
          } else {
            tags$span(style = "color: #666;",
                      "Donor status effect does NOT survive adjustment for age and gender.")
          }
        )
      )
    })

    # ===========================================================================
    # Card 6: AUC Analysis
    # ===========================================================================

    output$auc_plot <- renderPlotly({
      st <- stats_run()
      if (is.null(st) || is.null(st$auc_by_group) || nrow(st$auc_by_group) == 0) return(NULL)

      auc_df <- st$auc_by_group
      auc_df$donor_status <- factor(auc_df$donor_status, levels = c("ND", "Aab+", "T1D"))

      g <- ggplot(auc_df, aes(x = donor_status, y = auc, fill = donor_status)) +
        geom_col(width = 0.7, color = "black", alpha = 0.8) +
        scale_fill_manual(values = donor_colors()) +
        labs(x = "Donor Status", y = "Area Under Curve (AUC)",
             title = "Integrated Area by Donor Group") +
        theme_minimal(base_size = 14) +
        theme(legend.position = "none")

      ggplotly(g, tooltip = c("x", "y"))
    })

    output$auc_table <- renderTable({
      st <- stats_run()
      if (is.null(st) || is.null(st$auc_by_group) || nrow(st$auc_by_group) == 0) return(NULL)

      auc_df <- st$auc_by_group

      result <- data.frame(
        "Donor Group" = as.character(auc_df$donor_status),
        "AUC" = sprintf("%.2f", auc_df$auc),
        "N Bins" = as.character(auc_df$n_bins),
        check.names = FALSE, stringsAsFactors = FALSE
      )

      # Add T1D vs ND percentage change row (with division-by-zero guard)
      nd_auc <- auc_df$auc[auc_df$donor_status == "ND"]
      t1d_auc <- auc_df$auc[auc_df$donor_status == "T1D"]
      if (length(nd_auc) > 0 && length(t1d_auc) > 0 && is.finite(nd_auc) && nd_auc != 0) {
        pct_change <- ((t1d_auc - nd_auc) / nd_auc) * 100
        result <- rbind(result, data.frame(
          "Donor Group" = "T1D vs ND",
          "AUC" = sprintf("%.1f%%", pct_change),
          "N Bins" = "",
          check.names = FALSE, stringsAsFactors = FALSE
        ))
      }

      result
    }, striped = TRUE, bordered = TRUE, spacing = "s")

    output$auc_interpretation <- renderUI({
      st <- stats_run()
      if (is.null(st) || is.null(st$auc_by_group) || nrow(st$auc_by_group) == 0) return(NULL)

      auc_df <- st$auc_by_group
      nd_auc <- auc_df$auc[auc_df$donor_status == "ND"]
      t1d_auc <- auc_df$auc[auc_df$donor_status == "T1D"]

      if (length(nd_auc) > 0 && length(t1d_auc) > 0 && is.finite(nd_auc) && nd_auc != 0) {
        pct <- ((t1d_auc - nd_auc) / nd_auc) * 100
        direction <- if (pct < 0) "decrease" else "increase"
        tags$p(style = "margin-top: 10px; font-size: 90%; color: #555;",
          sprintf("The AUC integrates the binned mean feature value across islet diameters (trapezoidal rule). T1D shows a %.1f%% %s relative to ND.",
                  abs(pct), direction))
      } else {
        NULL
      }
    })

    # ===========================================================================
    # Card 7: Methods & Interpretation
    # ===========================================================================

    output$stats_explanation <- renderUI({
      st <- stats_run()
      if (is.null(st)) {
        return(tags$div(style = "padding: 10px; color: #666;",
          tags$p("Click 'Run Statistics' to perform analysis. ",
                 "The sidebar controls which feature and subset of data are analyzed.")))
      }

      alpha_str <- as.character(st$alpha)
      test_desc <- if (st$is_parametric) {
        "Parametric tests: one-way ANOVA (global), pairwise t-tests (BH-corrected)."
      } else {
        "Non-parametric tests: Kruskal-Wallis (global), pairwise Wilcoxon rank-sum (BH-corrected)."
      }

      tags$div(style = "padding: 10px;",
        tags$p(
          tags$strong("Pseudoreplication correction: "),
          paste0("Islets from the same donor are correlated (shared biology, tissue processing, imaging). ",
                 "Testing ", formatC(st$n_islets, big.mark = ","),
                 " islets as independent observations inflates significance (true df ~12, not ~",
                 formatC(st$n_islets - 3, big.mark = ","),
                 "). All tests use donor-level means (N = ", st$n_donors, ")."),
          tags$br(), tags$br(),
          tags$strong("Tests performed: "), test_desc,
          " All computed on donor-level aggregated means.",
          tags$br(),
          tags$strong("Significance threshold: "), paste0("\u03b1 = ", alpha_str),
          tags$br(),
          tags$strong("Mixed-effects model: "),
          "lmer(value ~ donor_status + (1 | donor)) on islet-level data provides a sensitivity analysis. ",
          "The ICC (intra-class correlation) quantifies how much variance is between donors vs within donors.",
          tags$br(),
          tags$strong("Effect sizes: "), "Cohen's d with 95% CI computed on donor-level means (N = 5 per group). ",
          "CIs are wider than islet-level, reflecting honest uncertainty about group differences.",
          tags$br(),
          tags$strong("Normality: "), "Shapiro-Wilk test on donor-level means per group. ",
          "With n = 5, this test has limited power, so failure to reject normality does not guarantee it.",
          tags$br(),
          tags$strong("Stratified analysis: "), "ANOVA and Kendall \u03c4 on donor-level means within each diameter bin. ",
          "P-values are BH-corrected across bins. Bins with <2 donors per group are untestable (greyed out). ",
          "Wider bin widths are recommended for donor-level analysis.",
          tags$br(),
          tags$strong("Multiple comparisons: "), "Benjamini-Hochberg (FDR) correction applied to pairwise p-values.",
          tags$br(),
          tags$strong("AUC: "), "Area under curve computed by trapezoidal integration of binned mean values across islet diameters, per donor group."
        ),
        if (!is.null(st$demo_summary)) {
          tags$p(style = "margin-top: 8px;",
            tags$strong("Demographic covariates: "),
            "Age-feature correlation (Pearson r) computed across all donors. ",
            "Gender-stratified pairwise tests check whether effects hold within each sex ",
            "(low power with ~2-3 donors per gender per group). ",
            if (!is.null(st$age_model_p)) {
              "A linear model on donor means adjusting for age and gender tests whether donor_status retains significance."
            } else {
              ""
            }
          )
        },
        tags$p(style = "font-size: 90%; color: #666; margin-top: 8px;",
          tags$strong("Interpretation: "),
          "A significant result (p < \u03b1) suggests the feature differs between donor groups. ",
          "With N = 5 per group, power is limited; non-significant results should not be over-interpreted. ",
          "Effect sizes (Cohen's d) provide magnitude information independent of sample size."
        )
      )
    })

    # ===========================================================================
    # Download handler
    # ===========================================================================

    output$download_stats <- downloadHandler(
      filename = function() {
        paste0("statistics_", gsub("[^0-9A-Za-z]+", "_", Sys.time()), ".csv")
      },
      content = function(file) {
        st <- stats_run()
        descriptions <- get_selection_description()

        lines <- c(paste("#", descriptions))
        lines <- c(lines, "")

        if (!is.null(st)) {
          # Methodology notes
          lines <- c(lines, "# METHODOLOGY")
          lines <- c(lines, paste0("# Tests use donor-level means (N = ", st$n_donors,
                                   ") to avoid pseudoreplication"))
          lines <- c(lines, paste0("# Total islets after filtering: ",
                                   formatC(st$n_islets, big.mark = ",")))
          lines <- c(lines, paste0("# ICC (intra-class correlation): ",
                                   if (is.na(st$icc)) "N/A" else sprintf("%.4f", st$icc)))
          lines <- c(lines, paste0("# Mixed-effects p-value: ",
                                   if (is.na(st$lmer_p)) "N/A" else formatC(st$lmer_p, format = "e", digits = 4)))
          lines <- c(lines, "")

          # Global test
          lines <- c(lines, "# Global Test (donor-level)")
          lines <- c(lines, paste0("# ", st$global_test, ", p = ",
                                   formatC(st$global_p, format = "e", digits = 4),
                                   ", eta_sq = ", round(st$eta_sq, 4)))
          lines <- c(lines, "")

          # Write header lines
          writeLines(lines, file)

          # Normality tests
          if (!is.null(st$norm_result)) {
            cat("# Shapiro-Wilk Normality Tests (donor-level)\n", file = file, append = TRUE)
            write.table(st$norm_result, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          # Donor-level data
          if (!is.null(st$donor_df)) {
            cat("\n# Donor-Level Means\n", file = file, append = TRUE)
            write.table(st$donor_df, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          # Pairwise
          cat("\n# Pairwise Comparisons (donor-level)\n", file = file, append = TRUE)
          pw <- st$pairwise_df
          if (!is.null(pw) && nrow(pw) > 0) {
            write.table(pw, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          # Per-bin ANOVA
          if (!is.null(st$bin_anova) && nrow(st$bin_anova) > 0) {
            cat("\n# Per-Bin ANOVA (donor-level)\n", file = file, append = TRUE)
            write.table(st$bin_anova, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          # Per-bin Kendall
          if (!is.null(st$bin_kendall) && nrow(st$bin_kendall) > 0) {
            cat("\n# Per-Bin Kendall (donor-level)\n", file = file, append = TRUE)
            write.table(st$bin_kendall, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          # AUC
          if (!is.null(st$auc_by_group) && nrow(st$auc_by_group) > 0) {
            cat("\n# AUC by Donor Group\n", file = file, append = TRUE)
            write.table(st$auc_by_group, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          # Demographics
          if (!is.null(st$demo_summary)) {
            cat("\n# Demographics Summary\n", file = file, append = TRUE)
            write.table(st$demo_summary, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          if (!is.null(st$age_corr)) {
            cat("\n# Age-Feature Correlation (all donors)\n", file = file, append = TRUE)
            write.table(st$age_corr, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }
        } else {
          writeLines(c(lines, "# No results computed"), file)
        }
      }
    )
  })
}
