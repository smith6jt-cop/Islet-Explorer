# ---------------------------------------------------------------------------
# Islet Explorer — Modernist shell
#
# Runnable as-is on synthetic data so you can see the shell and the wiring
# before swapping in the real donor table. Replace `demo_islets()` with your
# loader and the four renderPlots with your existing plot code.
#
# Deps: shiny, bslib, ggplot2 (+ sysfonts/showtext for Archivo in plots)
# ---------------------------------------------------------------------------

library(shiny)
library(bslib)
library(ggplot2)

for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
register_modernist_fonts()

MARKERS <- c("INS", "GCG", "SST", "PPY", "CD3e", "CD8a", "CD68", "SMA", "CD31")

demo_islets <- function(n = 5214, seed = 7) {
  set.seed(seed)
  donors <- data.frame(
    donor_id  = sprintf("nPOD-%04d", 6001:6015),
    status    = factor(rep(GROUP_LEVELS, each = 5), levels = GROUP_LEVELS),
    donor_age = c(24, 41, 55, 12, 33, 19, 28, 47, 8, 36, 15, 22, 31, 44, 61),
    donor_sex = c("F","M","F","M","F","M","F","M","F","M","M","F","M","F","M"),
    stringsAsFactors = FALSE
  )
  i <- sample(nrow(donors), n, replace = TRUE)
  d <- donors[i, ]
  shift <- c(ND = 0, `Aab+` = -0.5, T1D = -1.6)[as.character(d$status)]
  d$umap1 <- rnorm(n, shift * 1.4, 1.1)
  d$umap2 <- rnorm(n, shift * -0.6, 1.1)
  d$leiden <- factor(sample(0:5, n, TRUE))
  for (m in MARKERS) d[[m]] <- pmax(0, rnorm(n, 6 + shift * ifelse(m == "INS", 2.2, 0), 1.4))
  rownames(d) <- NULL
  d
}

PIPELINE <- list(
  list(key = "plot",       num = "01", label = "Plot"),
  list(key = "trajectory", num = "02", label = "Trajectory"),
  list(key = "statistics", num = "03", label = "Statistics"),
  list(key = "spatial",    num = "04", label = "Spatial")
)
REFERENCE <- list(
  list(key = "viewer", label = "Viewer", sub = "Image review · no data edits")
)

# Point Avivator at your OME-TIFF/Zarr source. Same-origin self-host if you
# ever want to inject CSS; otherwise this is a sealed box by design.
AVIVATOR_URL <- "https://avivator.gehlenborglab.org/"

panel_card <- function(title, meta = NULL, ...) {
  div(
    class = "ie-card",
    div(class = "ie-card-head",
        span(class = "ie-card-title", title),
        if (!is.null(meta)) span(class = "ie-card-meta", meta)),
    ...
  )
}

ui <- page_fillable(
  padding = 0, gap = 0,
  tags$head(
    tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    tags$link(rel = "stylesheet",
              href = "https://fonts.googleapis.com/css2?family=Archivo:wght@400;600;700;800&display=swap"),
    tags$link(rel = "stylesheet", href = "styles.css"),
    tags$link(rel = "stylesheet", href = "modernist-shiny.css")
  ),
  div(
    class = "ie-shell",

    railUI(
      "rail", PIPELINE, REFERENCE,
      meta = "5,214 islets · 15 donors<br>nPOD CODEX · 31 markers",
      footer = tagList(
        rail_palette_note("Phenotype palette", "High Contrast"),
        rail_palette_note("Donor palette", "Bright")
      )
    ),

    div(
      class = "ie-main",

      # one cohort bar for the whole pipeline; hidden on the Viewer
      conditionalPanel(
        "input['rail-tab'] != 'viewer'",
        cohortFilterUI("cohort")
      ),

      div(
        class = "ie-panels",
        navset_hidden(
          id = "panels",

          nav_panel_hidden(
            "plot",
            panel_card("Embedding", textOutput("plot_meta", inline = TRUE),
                       div(class = "ie-field", style = "margin-bottom:16px",
                           span(class = "ie-field-label", "Feature"),
                           selectInput("feature", NULL, MARKERS, "INS", width = "150px")),
                       plotOutput("umap", height = "460px"))
          ),

          nav_panel_hidden(
            "trajectory",
            panel_card("Pseudotime", textOutput("traj_meta", inline = TRUE),
                       plotOutput("traj", height = "440px"))
          ),

          nav_panel_hidden(
            "statistics",
            panel_card(
              "Mixed-effects contrasts",
              textOutput("stat_meta", inline = TRUE),
              div(class = "ie-cohortbar", style = "padding:0 0 16px;border-bottom:none",
                  ie_field("Significance threshold",
                           sliderInput("alpha", NULL, 0.001, 0.1, 0.05, 0.001,
                                       ticks = FALSE, width = "220px")),
                  ie_field("Multiple-test correction",
                           ie_segmented("corr", c(`BH FDR` = "bh",
                                                  Bonferroni = "bonf",
                                                  None = "none"), "bh"))),
              plotOutput("stats", height = "360px")
            )
          ),

          nav_panel_hidden(
            "spatial",
            panel_card("Tissue phenotypes", textOutput("spatial_meta", inline = TRUE),
                       plotOutput("spatial", height = "460px"))
          ),

          # Avivator keeps its own UI. We supply the frame and the read-only
          # contract; nothing from the cohort bar reaches it.
          nav_panel_hidden(
            "viewer",
            div(class = "ie-viewer-head",
                span(class = "ie-card-title", "Viewer"),
                span(class = "ie-badge", "Reference · read-only"),
                span(class = "ie-viewer-note", "Display only · no data edits")),
            div(class = "ie-viewer-frame",
                tags$iframe(src = AVIVATOR_URL, allow = "fullscreen",
                            title = "Avivator image viewer"))
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  tab    <- railServer("rail", "panels")
  cohort <- cohortFilterServer("cohort", demo_islets(), feature_cols = MARKERS)

  meta_txt <- reactive({
    d <- cohort$islets()
    sprintf("%s cells · %s", format(nrow(d), big.mark = ","),
            switch(cohort$norm(), robust = "robust z", z = "std z", "unnormalized"))
  })
  output$plot_meta     <- renderText(meta_txt())
  output$traj_meta     <- renderText(meta_txt())
  output$spatial_meta  <- renderText(meta_txt())

  output$umap <- renderPlot({
    d <- cohort$islets(); req(nrow(d) > 0)
    ggplot(d, aes(umap1, umap2, colour = status)) +
      geom_point(size = cohort$point_size(), alpha = cohort$point_alpha()) +
      scale_colour_donor() +
      labs(x = "UMAP 1", y = "UMAP 2") +
      theme_modernist(grid = "none")
  }, res = 96)

  output$traj <- renderPlot({
    d <- cohort$islets(); req(nrow(d) > 0)
    d$pt <- with(d, (umap1 - min(umap1)) / diff(range(umap1)))
    ggplot(d, aes(pt, .data[[input$feature]], colour = status)) +
      geom_point(size = cohort$point_size(), alpha = cohort$point_alpha() * 0.6) +
      geom_smooth(se = FALSE, linewidth = 1, method = "loess", formula = y ~ x) +
      scale_colour_donor() +
      labs(x = "Pseudotime", y = input$feature) +
      theme_modernist()
  }, res = 96)

  contrasts <- reactive({
    d <- cohort$islets(); req(nrow(d) > 0)
    combos <- utils::combn(levels(droplevels(d$status)), 2, simplify = FALSE)
    req(length(combos) > 0)
    factor_n <- switch(input$corr, bh = 1.4, bonf = length(combos), 1)
    do.call(rbind, lapply(combos, function(cc) {
      a <- d[[input$feature]][d$status == cc[1]]
      b <- d[[input$feature]][d$status == cc[2]]
      tt <- stats::t.test(a, b)
      data.frame(
        contrast = paste(cc, collapse = " vs "),
        est = unname(diff(rev(tt$estimate))),
        lo = tt$conf.int[1], hi = tt$conf.int[2],
        p_adj = min(1, tt$p.value * factor_n)
      )
    }))
  })

  output$stat_meta <- renderText({
    ct <- contrasts()
    sprintf("%d of %d contrasts significant at α = %s",
            sum(ct$p_adj < input$alpha), nrow(ct), format(input$alpha))
  })

  output$stats <- renderPlot({
    ct <- contrasts()
    ct$sig <- factor(ct$p_adj < input$alpha, levels = c("TRUE", "FALSE"))
    ggplot(ct, aes(est, contrast, colour = sig)) +
      geom_vline(xintercept = 0, colour = IE_INK, linewidth = 0.7) +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 2.5) +
      geom_point(size = 3.2) +
      scale_colour_significance() +
      labs(x = paste0("Difference in ", input$feature), y = NULL,
           caption = "Fixed effect: disease status · random intercept per donor") +
      theme_modernist(grid = "x")
  }, res = 96)

  output$spatial <- renderPlot({
    d <- cohort$islets(); req(nrow(d) > 0)
    d$pheno <- factor(sample(names(PAL_PHENO)[1:8], nrow(d), TRUE),
                      levels = names(PAL_PHENO))
    ggplot(d, aes(umap1, umap2, colour = pheno)) +
      geom_point(size = cohort$point_size(), alpha = cohort$point_alpha()) +
      scale_colour_pheno() +
      coord_fixed() +
      labs(x = "x (µm)", y = "y (µm)") +
      theme_modernist(grid = "none")
  }, res = 96)
}

shinyApp(ui, server)
