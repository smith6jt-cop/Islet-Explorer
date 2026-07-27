# ---------- Inflammation tab module server ----------
# Exports: inflammation_server(id, prepared, donor_colors)
#
# Two cohort-level analyses adapted from Steenbuck et al. 2025 (Bodenmiller lab,
# bioRxiv 2025.03.05.641526), plus the nPOD donor table they need for context.
#
# Dependencies (auto-sourced):
#   load_cell_distance()                   -- spatial_helpers.R (cached, 2.65M cells)
#   PHENOTYPE_COLORS                       -- drilldown_helpers.R
#   theme_islet, plotly_islet_*, ISLET_*   -- 00_globals.R
#   write_islet_csv, islet_csv_header      -- 00_globals.R
#   format_p_value                         -- utils_stats.R
#
# NOTE: validate/need are called as shiny::validate(shiny::need(...)). jsonlite
# masks `validate` on this app's search path, so the bare form resolves to the
# wrong function and dies with "is.character(txt) is not TRUE". app.R:420 uses
# the namespaced form for the same reason.
#
# Both heavy computations are plain `reactive()`s so they run once per session
# per setting and are shared by the plot, the table and the CSV download — the
# rule learned the hard way in the Jul 2026 CSV pass: an export must consume the
# same reactive the plot does, never re-derive its own subset.

# Immune phenotypes pre-selected in the infiltration card. Others remain
# available in the picker; these are the ones the infiltration score exists to
# interrogate.
INFLAM_DEFAULT_PHENOTYPES <- c("CD8a Tcell", "CD4 Tcell", "B cell",
                               "T cell", "Macrophage", "APCs", "Immune")

inflammation_server <- function(id, prepared, donor_colors = reactive(DONOR_COLORS)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    status_levels <- c("ND", "Aab+", "T1D")

    # =========================================================================
    # Shared data
    # =========================================================================

    # Tissue-wide single-cell table with the signed distance to the islet edge.
    cells <- reactive({
      df <- load_cell_distance()
      shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0,
                    "data/cell_distance.csv not found — run scripts/compute_neighborhood_metrics.py"))
      df$donor_status <- factor(as.character(df$donor_status), levels = status_levels)
      df[!is.na(df$donor_status), , drop = FALSE]
    })

    # One row per donor from the nPOD metadata joined into comp by prep_data().
    cohort <- reactive({
      pd <- prepared()
      shiny::validate(shiny::need(!is.null(pd$comp) && nrow(pd$comp) > 0, "No donor data loaded."))
      cn <- colnames(pd$comp)
      want <- c("Case ID", "Donor Status", "age", "age_years_npod", "age_group",
                "t1d_duration_years", "gender", "BMI", "HbA1c_percent",
                "C_peptide_ng_mL", "C_peptide_censored", "n_aab",
                "AAb_GADA", "AAb_IA2A", "AAb_mIAA", "AAb_ZnT8A")
      have <- intersect(want, cn)
      d <- unique(pd$comp[, have, drop = FALSE])
      d <- d[!duplicated(d$`Case ID`), , drop = FALSE]
      d$`Donor Status` <- factor(as.character(d$`Donor Status`), levels = status_levels)
      d[order(d$`Donor Status`, suppressWarnings(as.numeric(d$age))), , drop = FALSE]
    })

    npod_available <- reactive("t1d_duration_years" %in% colnames(cohort()))

    # =========================================================================
    # Card A: donor cohort table
    # =========================================================================

    output$cohort_note <- renderUI({
      d <- cohort()
      if (!npod_available()) {
        return(tags$div(style = paste0("background: #fdecea; border: 1px solid #f5b7b1;",
                                       " border-radius: 4px; padding: 8px 12px; margin-bottom: 10px;",
                                       " font-size: 12px; color: #8a1c12;"),
          tags$strong("nPOD clinical metadata not loaded. "),
          "Only age and sex are available. Run ",
          tags$code("python scripts/build_donor_metadata.py"),
          " to generate ", tags$code("data/donor_metadata.csv"), "."))
      }
      dur <- suppressWarnings(as.numeric(d$t1d_duration_years))
      dur <- dur[is.finite(dur)]
      dur_txt <- if (length(dur)) {
        sprintf("All %d T1D donors are recent-onset (duration %s year%s).",
                length(dur),
                if (diff(range(dur)) == 0) sprintf("%.2g", dur[1])
                else sprintf("%.2g–%.2g", min(dur), max(dur)),
                if (length(dur) == 1) "" else "s")
      } else ""
      # The one place the age discrepancy is visible rather than buried.
      mism <- if (all(c("age", "age_years_npod") %in% colnames(d))) {
        a <- suppressWarnings(as.numeric(d$age)); b <- suppressWarnings(as.numeric(d$age_years_npod))
        d$`Case ID`[is.finite(a) & is.finite(b) & abs(a - b) > 0.01]
      } else character(0)
      tagList(
        tags$p(style = "font-size: 12px; color: #555; margin: 0 0 8px 0; line-height: 1.5;",
          tags$strong("Scope. "), dur_txt,
          " This cohort therefore corresponds to the Steenbuck et al. “Onset” group; ",
          "it contains no long-duration T1D donors. Disease duration has too little ",
          "variance to serve as a covariate — it is reported here as a descriptor."),
        if (length(mism)) {
          tags$p(style = "font-size: 11px; color: #8a6d00; margin: 0 0 8px 0;",
            tags$strong("Note: "),
            sprintf("donor %s has a different age in the H5AD (%s) than in the nPOD record (%s). ",
                    paste(mism, collapse = ", "),
                    paste(sprintf("%.2f", suppressWarnings(as.numeric(d$age))[d$`Case ID` %in% mism]), collapse = ", "),
                    paste(sprintf("%.2f", suppressWarnings(as.numeric(d$age_years_npod))[d$`Case ID` %in% mism]), collapse = ", ")),
            "The H5AD value is used for all modelling so existing results do not shift.")
        } else NULL
      )
    })

    cohort_display <- reactive({
      d <- cohort()
      fmt <- function(x, digits = 1) ifelse(is.na(x), "—", sprintf(paste0("%.", digits, "f"), as.numeric(x)))
      aab_cols <- intersect(c("AAb_GADA", "AAb_IA2A", "AAb_mIAA", "AAb_ZnT8A"), colnames(d))
      profile <- if (length(aab_cols)) {
        vapply(seq_len(nrow(d)), function(i) {
          pos <- sub("^AAb_", "", aab_cols[vapply(aab_cols, function(cc) isTRUE(as.logical(d[[cc]][i])), logical(1))])
          if (!length(pos)) "—" else paste(pos, collapse = ", ")
        }, character(1))
      } else rep("—", nrow(d))

      out <- data.frame(
        `Case ID` = as.character(d$`Case ID`),
        Status = as.character(d$`Donor Status`),
        `Age (y)` = fmt(d$age, 2),
        Sex = if ("gender" %in% colnames(d)) as.character(d$gender) else "—",
        check.names = FALSE, stringsAsFactors = FALSE
      )
      if (npod_available()) {
        cp <- suppressWarnings(as.numeric(d$C_peptide_ng_mL))
        cens <- if ("C_peptide_censored" %in% colnames(d)) as.logical(d$C_peptide_censored) else rep(FALSE, nrow(d))
        out$`Age group` <- as.character(d$age_group)
        out$`T1D duration (y)` <- ifelse(is.na(d$t1d_duration_years), "n/a",
                                         fmt(d$t1d_duration_years, 2))
        out$BMI <- fmt(d$BMI, 1)
        out$`HbA1c (%)` <- fmt(d$HbA1c_percent, 1)
        # The censored value keeps its "<" so the table never implies a precision
        # the assay did not deliver.
        out$`C-peptide (ng/mL)` <- ifelse(cens %in% TRUE,
                                          paste0("<", fmt(cp, 2)), fmt(cp, 2))
        out$`AAb+` <- profile
        out$`n AAb` <- as.character(d$n_aab)
      }
      out
    })

    output$cohort_table <- renderTable({ cohort_display() },
      striped = TRUE, bordered = TRUE, spacing = "s", align = "l")

    output$dl_cohort <- downloadHandler(
      filename = function() paste0("inflammation_donor_cohort_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        write_islet_csv(file, cohort(), islet_csv_header(
          "Donor cohort — nPOD clinical metadata",
          "Source: scripts/build_donor_metadata.py -> data/donor_metadata.csv",
          "age = H5AD value (canonical for modelling); age_years_npod = nPOD record",
          "t1d_duration_years is NA for non-diabetic donors (not applicable, not zero)"))
      }
    )

    # =========================================================================
    # Card B: edgeR differential abundance
    # =========================================================================

    # Per-donor x cell-type count matrix. Counts, not proportions: that is the
    # entire point of this card.
    abundance <- reactive({
      df <- cells()
      region <- input$edger_region %||% "core_peri"
      keep <- switch(region,
                     core_peri = df$cell_region %in% c("core", "peri"),
                     core      = df$cell_region == "core",
                     rep(TRUE, nrow(df)))
      sub <- df[keep, , drop = FALSE]
      shiny::validate(shiny::need(nrow(sub) > 0, "No cells in the selected region."))
      tab <- table(sub$phenotype, sub$case_id)
      mat <- matrix(as.integer(tab), nrow = nrow(tab),
                    dimnames = list(rownames(tab), colnames(tab)))
      meta <- unique(sub[, c("case_id", "donor_status")])
      meta <- meta[match(colnames(mat), as.character(meta$case_id)), , drop = FALSE]
      list(counts = mat,
           group = factor(as.character(meta$donor_status), levels = status_levels),
           region = region, n_cells = nrow(sub))
    })

    edger_fit <- reactive({
      shiny::validate(shiny::need(requireNamespace("edgeR", quietly = TRUE),
                    "The edgeR package is not installed."))
      ab <- abundance()
      grp <- ab$group
      shiny::validate(shiny::need(nlevels(droplevels(grp)) >= 2,
                    "Need at least 2 donor-status groups."))

      y <- edgeR::DGEList(counts = ab$counts, group = grp)
      keep <- edgeR::filterByExpr(y, group = grp)
      y <- y[keep, , keep.lib.sizes = FALSE]
      y <- edgeR::calcNormFactors(y)
      design <- stats::model.matrix(~ grp)
      y <- edgeR::estimateDisp(y, design)
      fit <- edgeR::glmQLFit(y, design)

      # `~ grp` with levels ND / Aab+ / T1D gives coef 2 = Aab+ vs ND and
      # coef 3 = T1D vs ND; T1D vs Aab+ needs an explicit contrast.
      cn <- colnames(design)
      spec <- switch(input$edger_contrast %||% "t1d_nd",
        aab_nd  = list(coef = 2, label = "Aab+ vs ND"),
        t1d_nd  = list(coef = 3, label = "T1D vs ND"),
        t1d_aab = list(contrast = c(0, -1, 1), label = "T1D vs Aab+"))
      qlf <- if (!is.null(spec$contrast)) {
        shiny::validate(shiny::need(length(spec$contrast) == ncol(design),
                      "Contrast does not match the design (a donor group may be missing)."))
        edgeR::glmQLFTest(fit, contrast = spec$contrast)
      } else {
        shiny::validate(shiny::need(spec$coef <= ncol(design),
                      "Contrast unavailable — a donor group is missing from this subset."))
        edgeR::glmQLFTest(fit, coef = spec$coef)
      }
      res <- edgeR::topTags(qlf, n = Inf, adjust.method = "BH")$table
      res$cell_type <- rownames(res)
      rownames(res) <- NULL
      list(res = res[order(res$FDR, res$PValue), , drop = FALSE],
           label = spec$label,
           dispersion = y$common.dispersion,
           n_tested = nrow(y), n_dropped = sum(!keep),
           dropped = rownames(ab$counts)[!keep],
           n_cells = ab$n_cells, group = grp)
    })

    output$edger_note <- renderUI({
      e <- edger_fit()
      tags$p(style = "font-size: 12px; color: #666; margin: 4px 0 10px 0; line-height: 1.5;",
        tags$strong(e$label), " — ",
        sprintf("%s cells from %d donors; %d cell type%s tested",
                formatC(e$n_cells, format = "d", big.mark = ","),
                length(e$group), e$n_tested, if (e$n_tested == 1) "" else "s"),
        if (e$n_dropped > 0) sprintf(", %d dropped by filterByExpr (%s)",
                                     e$n_dropped, paste(e$dropped, collapse = ", ")) else "",
        sprintf(". Common dispersion %.3f (BCV %.2f).", e$dispersion, sqrt(e$dispersion)),
        tags$br(),
        tags$span(style = "color: #888;",
          "logFC is on a log2 scale after library-size normalisation, so it is a change in ",
          "abundance rather than in share of the tissue. FDR is Benjamini-Hochberg across ",
          "the cell types tested in this contrast."))
    })

    output$edger_table <- renderTable({
      e <- edger_fit()
      r <- e$res
      data.frame(
        `Cell type` = r$cell_type,
        `log2FC` = sprintf("%.2f", r$logFC),
        `logCPM` = sprintf("%.1f", r$logCPM),
        `p` = format_p_value(r$PValue),
        `FDR` = format_p_value(r$FDR),
        Sig = ifelse(r$FDR < 0.05, "*", ""),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    }, striped = TRUE, bordered = TRUE, spacing = "s", align = "lccccc")

    output$edger_volcano <- renderPlotly({
      e <- edger_fit()
      r <- e$res
      r$neglog10 <- -log10(pmax(r$FDR, .Machine$double.xmin))
      r$sig <- ifelse(r$FDR < 0.05,
                      ifelse(r$logFC > 0, "Up (FDR<0.05)", "Down (FDR<0.05)"), "NS")
      r$hover <- sprintf("%s<br>log2FC: %.2f<br>p: %s<br>FDR: %s",
                         r$cell_type, r$logFC,
                         format_p_value(r$PValue), format_p_value(r$FDR))
      g <- ggplot(r, aes(x = logFC, y = neglog10, colour = sig, text = hover)) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#999") +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "#999") +
        geom_point(size = 3, alpha = 0.85) +
        scale_colour_manual(values = c("Up (FDR<0.05)" = "#d62728",
                                       "Down (FDR<0.05)" = "#1f77b4",
                                       "NS" = "#bbbbbb"), name = NULL) +
        labs(x = paste0("log2 fold change — ", e$label),
             y = "-log10(FDR)") +
        theme_islet() +
        theme(legend.position = "bottom")
      # No in-plot point labels: they collide with each other and with the y-axis
      # ("3Beta cell") at every width tried, and the results table sits directly
      # beside this plot listing every cell type ordered by FDR. Identity comes
      # from the hover instead.

      fonts <- plotly_islet_fonts()
      ggplotly(g, tooltip = "text") %>%
        layout(font = fonts$global,
               # b covers the x-axis title AND the horizontal legend below it;
               # measured, not guessed (legend sat 16 px off-canvas at b = 90).
               margin = list(t = 20, b = 115, l = 70, r = 20),
               xaxis = list(title = list(font = fonts$axis_title),
                            tickfont = fonts$axis_tick, automargin = TRUE),
               yaxis = list(title = list(font = fonts$axis_title),
                            tickfont = fonts$axis_tick, automargin = TRUE),
               # y is a fraction of the PLOT AREA, which shrinks as the bottom
               # margin grows -- so enlarging the margin alone never fixed this.
               legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.14,
                             font = fonts$legend_text)) %>%
        plotly_islet_config("inflammation_edger_volcano")
    })

    output$dl_edger <- downloadHandler(
      filename = function() paste0("inflammation_edger_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        e <- tryCatch(edger_fit(), error = function(err) NULL)
        if (is.null(e)) return(write_islet_csv(file, NULL, "edgeR results unavailable"))
        write_islet_csv(file, e$res, islet_csv_header(
          paste0("edgeR differential abundance — ", e$label),
          "Negative-binomial GLM, quasi-likelihood F-test, library-size normalised",
          paste("Cells included:", e$n_cells, "| region:", input$edger_region %||% "core_peri"),
          paste("Cell types tested:", e$n_tested,
                "| dropped by filterByExpr:", e$n_dropped),
          sprintf("Common dispersion: %.4f (BCV %.3f)", e$dispersion, sqrt(e$dispersion)),
          "FDR = Benjamini-Hochberg across the cell types in this contrast"))
      }
    )

    # =========================================================================
    # Card C: infiltration score
    # =========================================================================

    output$infil_pheno_ui <- renderUI({
      df <- tryCatch(cells(), error = function(e) NULL)
      avail <- if (is.null(df)) character(0) else sort(unique(as.character(df$phenotype)))
      sel <- intersect(INFLAM_DEFAULT_PHENOTYPES, avail)
      if (!length(sel)) sel <- head(avail, 5)
      checkboxGroupInput(ns("infil_phenos"), "Phenotypes",
                         choices = avail, selected = sel, inline = TRUE)
    })

    # Steenbuck et al.: score 1 inside the islet, 1/d outside, normalised by the
    # total cell count of the image. We have whole-section images rather than
    # ROIs, so the denominator is the donor section.
    infiltration <- reactive({
      df <- cells()
      d <- suppressWarnings(as.numeric(df$signed_dist))
      # 1/d explodes as d -> 0; floor at 1 µm so a cell one micron outside the
      # edge cannot outweigh a cell inside it.
      score <- ifelse(!is.finite(d), NA_real_,
                      ifelse(d <= 0, 1, 1 / pmax(d, 1)))
      keep <- is.finite(score)
      agg <- stats::aggregate(list(score = score[keep]),
                              by = list(case_id = df$case_id[keep],
                                        donor_status = df$donor_status[keep],
                                        phenotype = as.character(df$phenotype)[keep]),
                              FUN = sum)
      tot <- stats::aggregate(list(n_cells = rep(1L, sum(keep))),
                              by = list(case_id = df$case_id[keep]), FUN = sum)
      out <- merge(agg, tot, by = "case_id")
      out$infiltration_score <- out$score / out$n_cells
      out$donor_status <- factor(as.character(out$donor_status), levels = status_levels)
      out
    })

    infil_selected <- reactive({
      sel <- input$infil_phenos
      shiny::validate(shiny::need(length(sel) > 0, "Select at least one phenotype."))
      s <- infiltration()
      s[s$phenotype %in% sel, , drop = FALSE]
    })

    # Donor-level tests, the paper's way: Kruskal-Wallis across the three stages
    # on donor scores, then pairwise Wilcoxon.
    #
    # Two deliberate choices:
    #
    #  * Correction runs across ALL phenotypes present, not just the ones ticked
    #    for display. Correcting across the selection would make every q-value a
    #    function of which checkboxes are on, so unticking phenotypes would
    #    shrink the correction and manufacture significance. The picker is a
    #    display filter only. (Same lesson as the Spatial Card B phenotype
    #    filter: filter the summary, never the test.)
    #  * Pairwise contrasts are gated on the omnibus being COMPUTABLE, not on it
    #    being significant. A pairwise post-hoc is never read without an omnibus
    #    behind it, which is the rule the Statistics per-bin heatmap follows.
    #    Gating on significance instead would hide the raw signal entirely at
    #    this sample size, which is a reporting choice, not a statistical one.
    infil_tests_all <- reactive({
      s <- infiltration()
      phenos <- sort(unique(s$phenotype))
      rows <- lapply(phenos, function(ph) {
        sub <- s[s$phenotype == ph, , drop = FALSE]
        g <- split(sub$infiltration_score, droplevels(sub$donor_status))
        med <- vapply(status_levels, function(l)
          if (!is.null(g[[l]]) && length(g[[l]])) stats::median(g[[l]]) else NA_real_, numeric(1))
        # Omnibus is computable only with >= 2 groups present.
        kw <- if (length(g) >= 2) {
          tryCatch(stats::kruskal.test(sub$infiltration_score,
                                       droplevels(sub$donor_status))$p.value,
                   error = function(e) NA_real_)
        } else NA_real_
        pw <- function(a, b) {
          if (!is.finite(kw)) return(NA_real_)          # no omnibus -> no post-hoc
          if (is.null(g[[a]]) || is.null(g[[b]]) || length(g[[a]]) < 2 || length(g[[b]]) < 2)
            return(NA_real_)
          tryCatch(suppressWarnings(stats::wilcox.test(g[[a]], g[[b]], exact = FALSE)$p.value),
                   error = function(e) NA_real_)
        }
        data.frame(phenotype = ph,
                   ND = med[["ND"]], `Aab+` = med[["Aab+"]], T1D = med[["T1D"]],
                   n_donors = nrow(sub),
                   p_kw = kw,
                   p_nd_aab = pw("ND", "Aab+"),
                   p_nd_t1d = pw("ND", "T1D"),
                   p_aab_t1d = pw("Aab+", "T1D"),
                   check.names = FALSE, stringsAsFactors = FALSE)
      })
      out <- do.call(rbind, rows)
      out$q_kw <- stats::p.adjust(out$p_kw, method = "BH")
      # Pairwise BH across phenotypes x the 3 contrasts, as one family.
      pw_cols <- c("p_nd_aab", "p_nd_t1d", "p_aab_t1d")
      flat <- unlist(out[pw_cols], use.names = FALSE)
      adj <- stats::p.adjust(flat, method = "BH")
      adj <- matrix(adj, ncol = length(pw_cols))
      for (i in seq_along(pw_cols)) out[[sub("^p_", "q_", pw_cols[i])]] <- adj[, i]
      out[order(out$p_kw, out$phenotype), , drop = FALSE]
    })

    # Display slice: the ticked phenotypes only. q-values come from the full
    # family above and do not change with the selection.
    infil_tests <- reactive({
      sel <- input$infil_phenos
      shiny::validate(shiny::need(length(sel) > 0, "Select at least one phenotype."))
      t <- infil_tests_all()
      t[t$phenotype %in% sel, , drop = FALSE]
    })

    output$infil_table <- renderTable({
      t <- infil_tests()
      fmt <- function(x) ifelse(is.na(x), "—", formatC(x, format = "g", digits = 3))
      data.frame(
        Phenotype = t$phenotype,
        ND = fmt(t$ND), `Aab+` = fmt(t$`Aab+`), T1D = fmt(t$T1D),
        `KW p` = format_p_value(t$p_kw),
        `KW q` = format_p_value(t$q_kw),
        `ND vs T1D p` = format_p_value(t$p_nd_t1d),
        `ND vs T1D q` = format_p_value(t$q_nd_t1d),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    }, striped = TRUE, bordered = TRUE, spacing = "s", align = "lccccccc")

    output$infil_note <- renderUI({
      t_all <- infil_tests_all()
      t <- infil_tests()
      n_all <- nrow(t_all)
      n_sig <- sum(!is.na(t_all$q_kw) & t_all$q_kw < 0.05)
      tags$p(style = "font-size: 11px; color: #666; margin: 8px 0 0 0; line-height: 1.5;",
        tags$strong("Testing. "),
        "Kruskal-Wallis across the three disease stages on donor-level scores, then ",
        "pairwise Wilcoxon. Benjamini-Hochberg runs across ",
        tags$strong(sprintf("all %d phenotypes present", n_all)),
        " (and, for the pairwise columns, across phenotypes x 3 contrasts) \u2014 ",
        tags$em("not"), " across the ", nrow(t), " you have ticked. The picker is a display ",
        "filter, so unticking phenotypes cannot shrink the correction and manufacture ",
        "significance. Pairwise contrasts are reported wherever the omnibus is computable, ",
        "so none is read without an omnibus behind it.",
        tags$br(),
        if (n_sig == 0) {
          tags$span(style = "color: #8a6d00;",
            sprintf("No phenotype survives FDR across %d cell types at n = 5 donors per stage. ", n_all),
            "The uncorrected KW p column shows where the raw signal is; treat it as ",
            "hypothesis-generating.")
        } else {
          tags$span(sprintf("%d of %d phenotypes survive FDR.", n_sig, n_all))
        })
    })

    infil_plot_obj <- reactive({
      s <- infil_selected()
      s$phenotype <- factor(s$phenotype, levels = sort(unique(s$phenotype)))
      s$hover <- sprintf("%s<br>%s — donor %s<br>score: %.3g",
                         s$phenotype, s$donor_status, s$case_id, s$infiltration_score)
      pal <- donor_colors()
      g <- ggplot(s, aes(x = phenotype, y = infiltration_score,
                         colour = donor_status, text = hover)) +
        geom_point(position = position_dodge(width = 0.65), size = 2.8, alpha = 0.85) +
        stat_summary(aes(group = donor_status), fun = stats::median, geom = "crossbar",
                     position = position_dodge(width = 0.65), width = 0.5,
                     linewidth = 0.4, show.legend = FALSE) +
        scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
        labs(x = NULL, y = "Infiltration score (per donor)") +
        theme_islet() +
        theme(legend.position = "bottom",
              axis.text.x = element_text(angle = 30, hjust = 1))
      if (isTRUE(input$infil_log)) {
        # Scores span orders of magnitude across phenotypes; pseudo_log keeps any
        # exact zeros visible instead of dropping them at -Inf.
        g <- g + scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1e-6, base = 10))
      }
      g
    })

    output$infil_plot <- renderPlotly({
      fonts <- plotly_islet_fonts()
      ggplotly(infil_plot_obj(), tooltip = "text") %>%
        layout(font = fonts$global,
               margin = list(t = 20, b = 145, l = 80, r = 20),
               xaxis = list(title = list(font = fonts$axis_title),
                            tickfont = fonts$axis_tick, automargin = TRUE),
               yaxis = list(title = list(font = fonts$axis_title),
                            tickfont = fonts$axis_tick, automargin = TRUE),
               legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.28,
                             font = fonts$legend_text)) %>%
        plotly_islet_config("inflammation_infiltration")
    })

    output$dl_infil <- downloadHandler(
      filename = function() paste0("inflammation_infiltration_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        s <- tryCatch(infil_selected(), error = function(e) NULL)
        t <- tryCatch(infil_tests_all(), error = function(e) NULL)
        if (is.null(s)) return(write_islet_csv(file, NULL, "No phenotypes selected"))
        hdr <- islet_csv_header(
          "Infiltration score (Steenbuck et al. 2025 definition)",
          "score = 1 inside the islet, 1/d outside (d = um from the islet edge, floored at 1)",
          "normalised by the donor's total cell count (paper normalises per ROI; this cohort has whole sections)",
          "source: data/cell_distance.csv",
          "tests: Kruskal-Wallis across stages; BH across ALL phenotypes present (not the display selection)",
          "pairwise Wilcoxon reported wherever the omnibus is computable; BH across phenotypes x 3 contrasts")
        # Two blocks in one file. Written through a single open connection:
        # write.csv/write.table with append = TRUE on a PATH re-opens the file and
        # emits R's "appending column names to file" warning (the reason
        # write_islet_csv() exists). Same pattern here for the second block.
        con <- base::file(file, open = "w", encoding = "UTF-8")
        on.exit(close(con), add = TRUE)
        writeLines(paste("#", hdr), con)
        writeLines("", con)
        writeLines("# Per-donor infiltration scores", con)
        utils::write.csv(s[, c("case_id", "donor_status", "phenotype",
                               "score", "n_cells", "infiltration_score")],
                         con, row.names = FALSE)
        if (!is.null(t)) {
          writeLines("", con)
          writeLines("# Donor-level tests (BH across all phenotypes present)", con)
          utils::write.csv(t, con, row.names = FALSE)
        }
      }
    )

    # =========================================================================
    # Methods
    # =========================================================================

    output$inflam_methods <- renderUI({
      tags$div(style = "padding: 10px; font-size: 13px; line-height: 1.6;",
        tags$p(
          tags$strong("Source. "),
          "Both analyses follow Steenbuck N., Damond N., et al., ",
          tags$em("Imaging mass cytometry reveals early β-cell dysfunction and changes in ",
                  "immune signatures during type 1 diabetes progression in human pancreata"),
          ", bioRxiv 2025.03.05.641526 (Bodenmiller and Atkinson labs). Their cohort is ",
          "88 nPOD donors across five stages with 75 imaging ROIs each; ours is 15 donors ",
          "across three stages with whole-section CODEX images. Where that difference ",
          "forces an adaptation it is stated below."),
        tags$p(
          tags$strong("Differential abundance. "),
          "Per-donor cell-type counts are modelled with edgeR: ",
          tags$code("filterByExpr"), " → ", tags$code("calcNormFactors"), " → ",
          tags$code("estimateDisp"), " → ", tags$code("glmQLFit"), " → ",
          tags$code("glmQLFTest"), ", with Benjamini-Hochberg across cell types. ",
          "Quasi-likelihood F-tests are used because they hold their nominal error rate at ",
          "small replicate numbers, which matters at n = 5 donors per stage. ",
          tags$strong("Counts rather than proportions is the substantive choice: "),
          "proportions must sum to 1, so the loss of β cells in T1D inflates every other ",
          "cell type's proportion whether or not its absolute abundance moved. The ",
          "library-size offset removes that artefact."),
        tags$p(
          tags$strong("Infiltration score. "),
          "Each cell scores 1 if it lies inside an islet and 1/d otherwise, with d the ",
          "distance in µm from the islet edge; per-donor sums are divided by that donor's ",
          "total cell count. The score rewards accumulation in the peri-islet space and ",
          "penalises accumulation in distal tissue, which a plain density cannot ",
          "distinguish. ", tags$strong("Adaptations: "),
          "the paper normalises per imaging ROI and we normalise per donor section (this ",
          "cohort has no ROIs); and 1/d is floored at d = 1 µm so a cell one micron ",
          "outside the edge cannot outweigh a cell inside it."),
        tags$p(
          tags$strong("Multiple comparisons. "),
          "edgeR FDR is BH across the cell types in the displayed contrast. Infiltration ",
          "Kruskal-Wallis p-values are BH across the displayed phenotypes. Neither is ",
          "corrected across the two cards, nor across contrasts — fix a contrast before ",
          "reading the FDR column. Pairwise Wilcoxon contrasts in the infiltration table are ",
          tags$strong("protected"), ": they are reported only for phenotypes whose omnibus ",
          "test survives correction, so no pairwise result is read without an omnibus ",
          "behind it."),
        tags$p(
          tags$strong("What this tab does not do. "),
          "It does not adjust for age, sex or BMI — those adjustments live in the ",
          "Statistics tab's companion models on a per-islet feature. It does not use ",
          "disease duration, which has no usable variance in this cohort. And it does not ",
          "adjust for HbA1c or C-peptide: both are consequences of β-cell loss rather ",
          "than confounders of it, so conditioning on them would bias the disease effect ",
          "toward the null."))
    })
  })
}
