library(shiny)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(plotly)
library(broom)

master_path <- file.path("..", "data", "master_results.xlsx")

# ---------- Data loading and wrangling ----------

safe_read_sheet <- function(path, sheet) {
  # Increase guess_max to reduce type misguesses on sparse columns
  tryCatch(readxl::read_excel(path, sheet = sheet, guess_max = 100000), error = function(e) NULL)
}

load_master <- function(path = master_path) {
  sheets <- readxl::excel_sheets(path)
  list(
    markers = safe_read_sheet(path, "Islet_Markers"),
    targets = safe_read_sheet(path, "Islet_Targets"),
    comp    = safe_read_sheet(path, "Islet_Composition"),
    lgals3  = safe_read_sheet(path, "LGALS3")
  )
}

add_islet_key <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if ("region" %in% names(df)) {
    df <- df %>% mutate(
      islet_key = str_extract(region, "Islet_\\d+")
    )
  }
  if (!"islet_key" %in% names(df) && all(c("name", "islet_id") %in% names(df))) {
    df <- df %>% mutate(islet_key = ifelse(str_detect(name, "^Islet_\\\
\\d+$"), name, paste0("Islet_", islet_id)))
  }
  df
}

compute_diameter_um <- function(area_um2) {
  area_um2 <- suppressWarnings(as.numeric(area_um2))
  ifelse(is.finite(area_um2) & area_um2 > 0, 2 * sqrt(area_um2 / pi), NA_real_)
}

prep_data <- function(master) {
  # Islet size proxy per islet (prefer core; fallback to union) for diameter
  targets <- master$targets %>% add_islet_key()
  core_area <- targets %>%
    filter(tolower(type) == "islet_core") %>%
    select(`Case ID`, `Donor Status`, islet_key, core_region_um2 = region_um2) %>%
    distinct()
  union_area <- targets %>%
    filter(tolower(type) == "islet_union") %>%
    select(`Case ID`, `Donor Status`, islet_key, union_region_um2 = region_um2) %>%
    distinct()
  size_area <- full_join(core_area, union_area, by = c("Case ID", "Donor Status", "islet_key")) %>%
    mutate(size_um2 = dplyr::coalesce(core_region_um2, union_region_um2),
           islet_diam_um = compute_diameter_um(size_um2)) %>%
    select(`Case ID`, `Donor Status`, islet_key, islet_diam_um)

  # Targets: keep all region types and later filter by user selection
  targets_all <- targets %>%
    select(`Case ID`, `Donor Status`, islet_key, type, class, area_um2, region_um2, area_density, count) %>%
    left_join(size_area, by = c("Case ID", "Donor Status", "islet_key"))

  # Markers with fraction positive / mean intensity (include LGALS3 sheet)
  markers <- master$markers %>% add_islet_key()
  markers_all <- markers %>%
    select(`Case ID`, `Donor Status`, islet_key, region_type, marker, n_cells, mean, SD, threshold, pos_count, pos_frac)
  # Add LGALS3 rows if available
  if (!is.null(master$lgals3) && nrow(master$lgals3) > 0) {
    g3 <- master$lgals3 %>% add_islet_key() %>%
      mutate(marker = as.character(marker)) %>%
      select(`Case ID`, `Donor Status`, islet_key, region_type, marker, n_cells, mean, SD, threshold, pos_count, pos_frac)
    if (!is.null(g3) && nrow(g3) > 0) {
      markers_all <- bind_rows(markers_all, g3)
    }
  }
  markers_all <- markers_all %>% left_join(size_area, by = c("Case ID", "Donor Status", "islet_key"))

  # Composition by islet
  comp <- master$comp %>% add_islet_key()
  comp <- comp %>% select(`Case ID`, `Donor Status`, islet_key, cells_total, Ins_single, Glu_single, Stt_single,
                          Multi_Pos, Triple_Neg, Ins_any, Glu_any, Stt_any)
  comp <- comp %>% left_join(size_area, by = c("Case ID", "Donor Status", "islet_key"))

  list(core_area = size_area, targets_all = targets_all, markers_all = markers_all, comp = comp)
}

bin_islet_sizes <- function(df, diam_col, width) {
  x <- suppressWarnings(as.numeric(df[[diam_col]]))
  max_x <- max(x, na.rm = TRUE)
  max_x <- ifelse(is.finite(max_x), max_x, 0)
  # Define numeric bins: [lo, hi)
  bin_lo <- floor(x / width) * width
  bin_hi <- bin_lo + width
  diam_mid <- bin_lo + width/2
  # Build a human-readable label
  diam_bin <- paste0("[", bin_lo, ", ", bin_hi, ")")
  df$diam_bin <- factor(diam_bin, levels = unique(diam_bin[order(bin_lo)]), ordered = TRUE)
  df$diam_mid <- diam_mid
  df
}

# Compute ND -> Aab+ -> T1D trajectory pseudotime using islet size within each group
compute_traj_pseudotime <- function(df, group_col = "donor_status", size_col = "islet_diam_um",
                                    traj_order = c("ND", "Aab+", "T1D")) {
  if (is.null(df) || nrow(df) == 0) return(numeric(0))
  g <- as.character(df[[group_col]])
  x <- suppressWarnings(as.numeric(df[[size_col]]))
  # Only keep trajectory groups present in data, in specified order
  present <- traj_order[traj_order %in% unique(g)]
  K <- length(present)
  if (K == 0) {
    # Fallback: size continuum
    lo <- suppressWarnings(min(x, na.rm = TRUE)); hi <- suppressWarnings(max(x, na.rm = TRUE)); rng <- hi - lo
    return(if (!is.finite(rng) || rng <= 0) rep(0, length(x)) else pmin(1, pmax(0, (x - lo) / rng)))
  }
  pt <- rep(NA_real_, length(x))
  for (i in seq_along(present)) {
    grp <- present[i]
    idx <- which(g == grp)
    if (length(idx) == 0) next
    xi <- x[idx]
    # Fractional rank within group; if only one islet, set to 0
    n_i <- sum(is.finite(xi))
    frac <- rep(0, length(idx))
    if (n_i > 1) {
      r <- rank(xi, ties.method = "average", na.last = "keep")
      frac <- ifelse(is.finite(r), (r - 1) / (n_i - 1), 0)
    }
    pt[idx] <- (i - 1 + frac) / K
  }
  # Clamp and replace any remaining NAs with segment starts
  pt[!is.finite(pt)] <- 0
  pmin(1, pmax(0, pt))
}

# Unbiased pseudotime from features (no group labels). Tries DiffusionMap -> principal curve -> PC1.
compute_unbiased_pseudotime <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(numeric(0))
  # Candidate numeric feature columns if present
  cand <- c("value", "islet_diam_um", "pos_count", "n_cells", "pos_frac", "mean", "SD",
            "threshold", "count", "area_density", "cells_total", "Ins_any", "Glu_any", "Stt_any")
  present <- intersect(cand, colnames(df))
  if (length(present) == 0) {
    v <- suppressWarnings(as.numeric(df$islet_diam_um))
    v[!is.finite(v)] <- 0
    r <- range(v, na.rm = TRUE); d <- r[2] - r[1]
    return(if (!is.finite(d) || d <= 0) rep(0, length(v)) else (v - r[1]) / d)
  }
  X <- as.data.frame(lapply(present, function(nm) suppressWarnings(as.numeric(df[[nm]]))))
  colnames(X) <- present
  # Drop columns with all NA or zero variance
  good <- vapply(X, function(x) {
    x <- as.numeric(x)
    any(is.finite(x)) && (suppressWarnings(stats::sd(x, na.rm = TRUE)) > 0)
  }, logical(1))
  X <- X[, good, drop = FALSE]
  if (ncol(X) == 0) {
    v <- suppressWarnings(as.numeric(df$islet_diam_um))
    v[!is.finite(v)] <- 0
    r <- range(v, na.rm = TRUE); d <- r[2] - r[1]
    return(if (!is.finite(d) || d <= 0) rep(0, length(v)) else (v - r[1]) / d)
  }
  # Simple impute NAs with column medians, then scale
  for (j in seq_len(ncol(X))) {
    x <- as.numeric(X[[j]])
    med <- suppressWarnings(stats::median(x, na.rm = TRUE))
    if (!is.finite(med)) med <- 0
    x[!is.finite(x)] <- med
    X[[j]] <- x
  }
  Xs <- scale(as.matrix(X))
  pt <- NULL
  # Try Diffusion Map if available
  if (requireNamespace("destiny", quietly = TRUE)) {
    pt <- tryCatch({
      dm <- destiny::DiffusionMap(Xs)
      ev <- destiny::eigenvectors(dm)
      as.numeric(ev[, 1])
    }, error = function(e) NULL)
  }
  # Try principal curve on top-2 PCs
  if (is.null(pt) && requireNamespace("princurve", quietly = TRUE)) {
    pc <- tryCatch(stats::prcomp(Xs, center = TRUE, scale. = TRUE), error = function(e) NULL)
    if (!is.null(pc) && ncol(pc$x) >= 2) {
      pt <- tryCatch({
        fit <- princurve::principal_curve(as.matrix(pc$x[, 1:2, drop = FALSE]))
        as.numeric(fit$lambda)
      }, error = function(e) NULL)
    }
  }
  # Fallback: PC1
  if (is.null(pt)) {
    pc <- tryCatch(stats::prcomp(Xs, center = TRUE, scale. = TRUE), error = function(e) NULL)
    if (!is.null(pc)) pt <- as.numeric(pc$x[, 1])
  }
  if (is.null(pt)) pt <- rep(0, nrow(df))
  # Normalize to [0,1]
  r <- range(pt, na.rm = TRUE); d <- r[2] - r[1]
  if (!is.finite(d) || d <= 0) return(rep(0, length(pt)))
  pmin(1, pmax(0, (pt - r[1]) / d))
}

summary_stats <- function(df, group_cols, value_col, stat = c("mean_se","mean_sd","median_iqr")) {
  stat <- match.arg(stat)
  df %>% group_by(across(all_of(group_cols))) %>%
    summarise(
      n = sum(!is.na(.data[[value_col]])),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd = sd(.data[[value_col]], na.rm = TRUE),
      se = ifelse(n > 0, sd / sqrt(pmax(n, 1)), NA_real_),
      median = median(.data[[value_col]], na.rm = TRUE),
      q1 = quantile(.data[[value_col]], 0.25, na.rm = TRUE, type = 7),
      q3 = quantile(.data[[value_col]], 0.75, na.rm = TRUE, type = 7),
      .groups = "drop"
    ) %>%
    mutate(
      y = dplyr::case_when(
        stat == "mean_se" ~ mean,
        stat == "mean_sd" ~ mean,
        stat == "median_iqr" ~ median,
        TRUE ~ mean
      ),
      ymin = dplyr::case_when(
        stat == "mean_se" ~ mean - se,
        stat == "mean_sd" ~ mean - sd,
        stat == "median_iqr" ~ q1,
        TRUE ~ mean - se
      ),
      ymax = dplyr::case_when(
        stat == "mean_se" ~ mean + se,
        stat == "mean_sd" ~ mean + sd,
        stat == "median_iqr" ~ q3,
        TRUE ~ mean + se
      )
    )
}

per_bin_anova <- function(df, bin_col, group_col, value_col, mid_col = "diam_mid") {
  if (nrow(df) == 0) return(tibble(bin = character(), mid = numeric(), p_anova = numeric()))
  # capture unique bins with their numeric mid
  bmeta <- df %>% filter(!is.na(.data[[bin_col]])) %>%
    group_by(.data[[bin_col]]) %>%
    summarise(mid = suppressWarnings(as.numeric(first(na.omit(.data[[mid_col]])))), .groups = "drop")
  res <- lapply(seq_len(nrow(bmeta)), function(i) {
    b <- bmeta[[bin_col]][i]
    sub <- df %>% filter(.data[[bin_col]] == b)
    # Require at least 2 donor groups present and at least some finite values
    if (n_distinct(na.omit(sub[[group_col]])) < 2 || sum(is.finite(sub[[value_col]])) < 2) return(NULL)
    fit <- tryCatch(aov(reformulate(group_col, response = value_col), data = sub), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    # Extract one-way ANOVA p-value from ANOVA table (first term)
    pval <- tryCatch({
      at <- anova(fit)
      as.numeric(at[["Pr(>F)"]][1])
    }, error = function(e) NA_real_)
    tibble(
      bin = as.character(b),
      mid = bmeta$mid[i],
      p_anova = pval
    )
  })
  out <- bind_rows(res)
  if (nrow(out) == 0) return(tibble(bin = character(), mid = numeric(), p_anova = numeric()))
  arrange(out, mid)
}

per_bin_kendall <- function(df, bin_col, group_col, value_col, mid_col = "diam_mid") {
  if (nrow(df) == 0) return(tibble(bin = character(), mid = numeric(), p_kendall = numeric(), tau = numeric()))
  # Map group to ordered numeric along pseudotime ND < Aab+ < T1D
  code_group <- function(x) {
    x <- as.character(x)
    ifelse(x == "ND", 0, ifelse(x == "Aab+", 1, ifelse(x == "T1D", 2, NA_real_)))
  }
  bmeta <- df %>% filter(!is.na(.data[[bin_col]])) %>% group_by(.data[[bin_col]]) %>% summarise(mid = suppressWarnings(as.numeric(first(na.omit(.data[[mid_col]])))), .groups = "drop")
  res <- lapply(seq_len(nrow(bmeta)), function(i) {
    b <- bmeta[[bin_col]][i]
    sub <- df %>% filter(.data[[bin_col]] == b)
    x <- code_group(sub[[group_col]])
    y <- suppressWarnings(as.numeric(sub[[value_col]]))
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]; y <- y[keep]
    if (length(unique(x)) < 2 || length(y) < 3) return(NULL)
    ct <- tryCatch(cor.test(x, y, method = "kendall", exact = FALSE), error = function(e) NULL)
    if (is.null(ct)) return(NULL)
    tibble(bin = as.character(b), mid = bmeta$mid[i], p_kendall = unname(ct$p.value), tau = unname(ct$estimate))
  })
  out <- bind_rows(res)
  if (nrow(out) == 0) return(out)
  arrange(out, mid)
}

# ---------- UI ----------

ui <- fluidPage(
  titlePanel("Islet Area Distributions"),
  uiOutput("theme_css"),
  sidebarLayout(
    sidebarPanel(
      h4("Data & Filters"),
      width = 3,
      selectInput("mode", "Focus", choices = c("Markers", "Targets", "Composition"), selected = "Composition"),
      uiOutput("region_selector"),
      uiOutput("dynamic_selector"),
      uiOutput("metric_selector"),
      uiOutput("color_selector"),
      hr(),
      h5("Autoantibody filter (Aab+ donors only)"),
      checkboxGroupInput("aab_flags", NULL,
                         choices = c("GADA" = "AAb_GADA",
                                     "IA2A" = "AAb_IA2A",
                                     "ZnT8A" = "AAb_ZnT8A",
                                     "IAA" = "AAb_IAA",
                                     "mIAA" = "AAb_mIAA"),
                         selected = character(0)),
      radioButtons("aab_logic", "Match (within Aab+)", choices = c("Any", "All"), selected = "Any", inline = TRUE),
      helpText("Default: all Aab+ donors are included. Selecting one or more autoantibodies restricts only the Aab+ group to donors matching the selection (Any = at least one selected AAb, All = all selected AAbs). ND and T1D groups are unaffected."),
      checkboxGroupInput("groups", "Donor Status", choices = c("ND", "Aab+", "T1D"), selected = c("ND", "Aab+", "T1D")),
      sliderInput("binwidth", "Diameter bin width (µm)", min = 10, max = 100, value = 50, step = 5),
      radioButtons("stat", "Statistic",
                   choices = c("Mean±SE" = "mean_se",
                               "Mean±SD" = "mean_sd",
                               "Median + IQR" = "median_iqr"),
                   selected = "mean_se"),
      selectInput("alpha", "Significance level (alpha)", choices = c("0.05","0.01","0.001"), selected = "0.05"),
      checkboxInput("show_points", "Show individual points", value = FALSE),
      sliderInput("pt_size", "Point size", min = 0.3, max = 4.0, value = 0.8, step = 0.1),
      sliderInput("pt_alpha", "Point transparency", min = 0.05, max = 1.0, value = 0.25, step = 0.05),
      radioButtons("add_smooth", "Trend line", choices = c("None", "LOESS"), selected = "None", inline = TRUE),
      radioButtons("theme_bg", "Background", choices = c("Light","Dark"), selected = "Light", inline = TRUE),
      hr(),
      h5("Export"),
      downloadButton("dl_summary", "Download summary CSV"),
      downloadButton("dl_stats", "Download stats CSV"),
      actionButton("run_tests", "Run statistics")
    ),
    mainPanel(
      tabsetPanel(id = "tabs",
        tabPanel("Plot", plotlyOutput("plt", height = 650), br(), uiOutput("pseudo_ui"), plotlyOutput("pseudo", height = 400)),
        tabPanel("Statistics", tableOutput("stats_tbl")),
        tabPanel("Data Preview", tableOutput("preview"))
      )
    )
  )
)

# ---------- Server ----------

server <- function(input, output, session) {
  validate_file <- reactive({
    validate(need(file.exists(master_path), paste("Not found:", master_path)))
    master_path
  })

  master <- reactive({
    req(validate_file())
    load_master(master_path)
  })

  prepared <- reactive({
    prep_data(master())
  })

  # Raw per-islet dataset (no binning), aligned with current selections
  raw_df <- reactive({
    pd <- prepared()
    groups <- input$groups
    w <- input$which
    # Resolve composition default early
    if (identical(input$mode, "Composition")) {
      if (is.null(w) || !(w %in% c("Ins_any","Glu_any","Stt_any"))) w <- "Ins_any"
    }
    if (identical(input$mode, "Targets")) {
      req(input$region, !is.null(w) && nzchar(w))
      region_tag <- paste0("islet_", tolower(input$region))
      df <- pd$targets_all %>% dplyr::filter(`Donor Status` %in% groups, class == w, tolower(type) == region_tag)
      df <- df %>% dplyr::filter(is.finite(islet_diam_um))
      df <- df %>% dplyr::mutate(value = as.numeric(if (!is.null(input$target_metric) && input$target_metric == "Counts") count else area_density))
    } else if (identical(input$mode, "Markers")) {
      req(input$region, !is.null(w) && nzchar(w))
      region_tag <- paste0("islet_", tolower(input$region))
      df <- pd$markers_all %>% dplyr::filter(`Donor Status` %in% groups, marker == w, tolower(region_type) == region_tag)
      df <- df %>% dplyr::filter(is.finite(islet_diam_um))
      if (!is.null(input$marker_metric) && input$marker_metric == "Counts") {
        df <- df %>% dplyr::mutate(value = suppressWarnings(as.numeric(pos_count)))
      } else {
        num <- suppressWarnings(as.numeric(df$pos_count))
        den <- suppressWarnings(as.numeric(df$n_cells))
        df <- df %>% dplyr::mutate(value = ifelse(is.finite(num) & is.finite(den) & den > 0, 100.0 * num / den, NA_real_))
      }
    } else {
      df <- pd$comp %>% dplyr::filter(`Donor Status` %in% groups) %>% dplyr::filter(is.finite(islet_diam_um))
      num <- suppressWarnings(as.numeric(df[[w]]))
      den <- suppressWarnings(as.numeric(df$cells_total))
      df <- df %>% dplyr::mutate(value = ifelse(is.finite(num) & is.finite(den) & den > 0, 100.0 * num / den, NA_real_))
    }
    out <- df %>% dplyr::mutate(donor_status = `Donor Status`) %>% dplyr::filter(!is.na(value))
    # Apply AAb filters (only within the Aab+ donor group) if selected
    flags <- input$aab_flags
    if (!is.null(flags) && length(flags) > 0 && all(flags %in% colnames(out))) {
      # Keep ND and T1D unchanged; restrict Aab+ by selected AAbs
      others <- out %>% dplyr::filter(donor_status != "Aab+")
      aabp   <- out %>% dplyr::filter(donor_status == "Aab+")
      if (nrow(aabp) > 0) {
        mat <- as.data.frame(aabp[, flags, drop = FALSE])
        for (cc in colnames(mat)) mat[[cc]] <- as.logical(mat[[cc]])
        hits <- rowSums(mat, na.rm = TRUE)
        keep <- if (identical(input$aab_logic, "All")) hits >= length(flags) else hits >= 1
        aabp <- aabp[keep, , drop = FALSE]
      }
      out <- dplyr::bind_rows(others, aabp)
    }
    out
  })

  # App-wide CSS for theme background
  output$theme_css <- renderUI({
    if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
      tags$style(HTML("body { background-color: #000000; color: #e6e6e6; } .well { background-color: #111111; }"))
    } else {
      tags$style(HTML("body { background-color: #ffffff; color: #111111; }"))
    }
  })

  output$dynamic_selector <- renderUI({
    if (identical(input$mode, "Targets")) {
      req(input$region)
      region_tag <- paste0("islet_", tolower(input$region))
      classes <- prepared()$targets_all %>%
        filter(tolower(type) == region_tag) %>%
        pull(class) %>% unique() %>% na.omit() %>% sort()
      if (length(classes) == 0) classes <- character(0)
      selectInput("which", "Target class", choices = classes, selected = if (length(classes)>0) classes[1] else NULL)
    } else if (identical(input$mode, "Markers")) {
      req(input$region)
      region_tag <- paste0("islet_", tolower(input$region))
      markers <- prepared()$markers_all %>%
        filter(tolower(region_type) == region_tag) %>%
        pull(marker) %>% unique() %>% na.omit() %>% sort()
      if (length(markers) == 0) markers <- character(0)
      selectInput("which", "Marker", choices = markers, selected = if (length(markers)>0) markers[1] else NULL)
    } else {
      selectInput("which", "Composition measure", choices = c("Ins_frac" = "Ins_any", "Glu_frac" = "Glu_any", "Stt_frac" = "Stt_any"), selected = "Ins_any")
    }
  })

  output$region_selector <- renderUI({
    if (identical(input$mode, "Composition")) return(NULL)
    # Display labels while preserving underlying values used in code
    selectInput(
      "region", "Region",
      choices = c("Islet" = "core", "Peri-Islet" = "band", "Islet+20um" = "union"),
      selected = "band"
    )
  })

  output$metric_selector <- renderUI({
    if (identical(input$mode, "Targets")) {
      radioButtons("target_metric", "Targets metric", choices = c("Counts", "Density"), selected = "Density", inline = TRUE)
    } else if (identical(input$mode, "Markers")) {
      radioButtons("marker_metric", "Markers metric", choices = c("Counts", "% positive"), selected = "% positive", inline = TRUE)
    } else {
      NULL
    }
  })

  output$color_selector <- renderUI({
    # Only relevant when counts are selected (for pseudotime scatter)
    show <- (identical(input$mode, "Targets") && !is.null(input$target_metric) && input$target_metric == "Counts") ||
            (identical(input$mode, "Markers") && !is.null(input$marker_metric) && input$marker_metric == "Counts")
    if (!show) return(NULL)
    if (identical(input$mode, "Targets")) {
      selectInput("color_by", "Color by", choices = c("Donor Status" = "donor_status", "Target class" = "class", "Case ID" = "Case ID"), selected = "donor_status")
    } else if (identical(input$mode, "Markers")) {
      selectInput("color_by", "Color by", choices = c("Donor Status" = "donor_status", "Marker" = "marker", "Case ID" = "Case ID"), selected = "donor_status")
    } else {
      NULL
    }
  })

  # Dataset for plotting given mode
  plot_df <- reactive({
    pd <- prepared()
    groups <- input$groups
    w <- input$which
    bw <- input$binwidth
    # Ensure a selection is made for modes that require a variable
    if (identical(input$mode, "Targets") || identical(input$mode, "Markers")) {
      req(!is.null(w), length(w) > 0)
    }
    # Provide a safe default for Composition before input$which initializes
    if (identical(input$mode, "Composition")) {
      if (is.null(w) || !(w %in% c("Ins_any","Glu_any","Stt_any"))) {
        w <- "Ins_any"
      }
    }

    if (identical(input$mode, "Targets")) {
      req(input$region)
      region_tag <- paste0("islet_", tolower(input$region))
      df <- pd$targets_all %>% filter(`Donor Status` %in% groups, class == w, tolower(type) == region_tag)
      df <- df %>% filter(is.finite(islet_diam_um))
      df <- bin_islet_sizes(df, "islet_diam_um", bw)
      if (!is.null(input$target_metric) && input$target_metric == "Counts") {
        df <- df %>% mutate(value = as.numeric(count))
      } else {
        df <- df %>% mutate(value = as.numeric(area_density))
      }
    } else if (identical(input$mode, "Markers")) {
      req(input$region)
      region_tag <- paste0("islet_", tolower(input$region))
      df <- pd$markers_all %>% filter(`Donor Status` %in% groups, marker == w, tolower(region_type) == region_tag)
      df <- df %>% filter(is.finite(islet_diam_um))
      df <- bin_islet_sizes(df, "islet_diam_um", bw)
      # Markers metric: counts or % positive (n positive / total * 100)
      if (!is.null(input$marker_metric) && input$marker_metric == "Counts") {
        df <- df %>% mutate(value = suppressWarnings(as.numeric(pos_count)))
      } else {
        num <- suppressWarnings(as.numeric(df$pos_count))
        den <- suppressWarnings(as.numeric(df$n_cells))
        df <- df %>% mutate(value = ifelse(is.finite(num) & is.finite(den) & den > 0, 100.0 * num / den, NA_real_))
      }
    } else {
      df <- pd$comp %>% filter(`Donor Status` %in% groups)
      df <- df %>% filter(is.finite(islet_diam_um))
      df <- bin_islet_sizes(df, "islet_diam_um", bw)
      # Composition fraction as percent of total cells
      num <- suppressWarnings(as.numeric(df[[w]]))
      den <- suppressWarnings(as.numeric(df$cells_total))
      df <- df %>% mutate(value = ifelse(is.finite(num) & is.finite(den) & den > 0, 100.0 * num / den, NA_real_))
    }

    # Apply AAb filters (only within the Aab+ donor group) if requested
    out <- df %>% mutate(donor_status = `Donor Status`) %>% filter(!is.na(value))
    flags <- input$aab_flags
    if (!is.null(flags) && length(flags) > 0 && all(flags %in% colnames(out))) {
      others <- out %>% dplyr::filter(donor_status != "Aab+")
      aabp   <- out %>% dplyr::filter(donor_status == "Aab+")
      if (nrow(aabp) > 0) {
        mat <- as.data.frame(aabp[, flags, drop = FALSE])
        for (cc in colnames(mat)) mat[[cc]] <- as.logical(mat[[cc]])
        hits <- rowSums(mat, na.rm = TRUE)
        keep <- if (identical(input$aab_logic, "All")) hits >= length(flags) else hits >= 1
        aabp <- aabp[keep, , drop = FALSE]
      }
      out <- dplyr::bind_rows(others, aabp)
    }
    out
  })

  # Summary with mean and SE per bin per group
  summary_df <- reactive({
    df <- plot_df()
    req(nrow(df) > 0)
    summary_stats(df, group_cols = c("donor_status", "diam_bin", "diam_mid"), value_col = "value", stat = input$stat)
  })

  output$plt <- renderPlotly({
    sm <- summary_df()
    req(nrow(sm) > 0)

    # Order groups ND, Aab+, T1D and set colors
    grp_levels <- c("ND","Aab+","T1D")
    sm$donor_status <- factor(sm$donor_status, levels = grp_levels)
    color_map <- c("ND" = "#1f77b4", "Aab+" = "#ff7f0e", "T1D" = "#d62728")

    p <- ggplot(sm, aes(x = diam_mid, y = y, color = donor_status, group = donor_status)) +
      geom_line(alpha = ifelse(!is.null(input$add_smooth) && input$add_smooth == "LOESS", 0, 1)) +
      geom_point() +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) +
      labs(x = "Islet diameter (µm)",
           y = if (identical(input$mode, "Targets")) {
                 if (!is.null(input$target_metric) && input$target_metric == "Counts") "Target count" else "Target density (per µm²)"
               } else if (identical(input$mode, "Markers")) {
                 if (!is.null(input$marker_metric) && input$marker_metric == "Counts") "Positive cell count" else "n positive / total (%)"
               } else {
                 "% composition"
               },
           color = "Donor Status",
           title = paste0(if (identical(input$mode, "Targets")) input$which else if (identical(input$mode, "Markers")) input$which else input$which,
                          " vs islet size")) +
      scale_color_manual(values = color_map, breaks = grp_levels, drop = FALSE) +
      theme_minimal(base_size = 14)

    if (isTRUE(input$show_points)) {
      raw <- plot_df()
      raw$donor_status <- factor(raw$donor_status, levels = grp_levels)
      jw <- max(1, as.numeric(input$binwidth) * 0.35)
      p <- p +
        geom_point(data = raw, aes(x = diam_mid, y = value, color = donor_status),
                   position = position_jitter(width = jw, height = 0), size = input$pt_size, alpha = input$pt_alpha, inherit.aes = FALSE)
    }

    # Axis breaks and minor ticks
    xmax <- suppressWarnings(max(sm$diam_mid, na.rm = TRUE))
    if (!is.finite(xmax)) xmax <- 300
    major_breaks <- seq(0, ceiling(xmax/50)*50, by = 50)
    minor_breaks <- seq(0, ceiling(xmax/10)*10, by = 10)
    p <- p + scale_x_continuous(breaks = major_breaks, minor_breaks = minor_breaks)
    p <- p + theme(panel.grid.minor = element_line(size = 0.2, colour = if (!is.null(input$theme_bg) && input$theme_bg == "Dark") "#222222" else "#eeeeee"))

    # Highlight significant bins if stats were run
    st <- NULL
    if (!is.null(st) && nrow(st) > 0) {
      alpha_num <- as.numeric(input$alpha)
      st <- st %>% mutate(p_adj = p.adjust(p_value, method = "BH"), sig = p_adj <= alpha_num)
      sig <- st %>% filter(sig)
      if (nrow(sig) > 0) {
        ytop <- suppressWarnings(max(sm$ymax, na.rm = TRUE))
        ytop <- ifelse(is.finite(ytop), ytop, suppressWarnings(max(sm$y, na.rm = TRUE)))
        sig$y <- ytop * 1.05
        p <- p + geom_text(data = sig, aes(x = mid, y = y), label = ifelse(unique(st$test) == 'kendall', '816', '*'), inherit.aes = FALSE, color = "#444444")
      }
    }

    # Optional smoothing overlay (kept subtle)
    if (!is.null(input$add_smooth) && input$add_smooth == "LOESS") {
      p <- p + geom_smooth(se = FALSE, method = "loess", span = 0.6, size = 0.9)
    }

    # Apply dark theme inside the plot if selected
    if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
      p <- p + theme(
        plot.background = element_rect(fill = "#000000", colour = NA),
        panel.background = element_rect(fill = "#000000", colour = NA),
        panel.grid.major = element_line(color = "#333333"),
        axis.text = element_text(color = "#e6e6e6"),
        axis.title = element_text(color = "#f0f0f0"),
        plot.title = element_text(color = "#f0f0f0"),
        legend.text = element_text(color = "#e6e6e6"),
        legend.title = element_text(color = "#f0f0f0")
      )
    }

    gg <- ggplotly(p)
    gg <- gg %>% layout(legend = list(orientation = "h", x = 0, y = -0.15))
    gg
  })

  # Stats: continuous tests (no binning). Global LM + post-hoc pairwise on residuals
  stats_run <- eventReactive(input$run_tests, {
    rdf <- raw_df()
    if (is.null(rdf) || !nrow(rdf)) return(NULL)
    # Ensure donor_status factor ordering
    rdf$donor_status <- factor(rdf$donor_status, levels = c("ND","Aab+","T1D"))
    # Global: value ~ donor_status + islet_diam_um (Type I with donor_status first)
    fit <- tryCatch(lm(value ~ donor_status + islet_diam_um, data = rdf), error = function(e) NULL)
    p_global <- NA_real_
    if (!is.null(fit)) {
      at <- tryCatch(anova(fit), error = function(e) NULL)
      if (!is.null(at)) p_global <- suppressWarnings(as.numeric(at[["Pr(>F)"]][1]))
    }
    # Post-hoc: residualize value ~ islet_diam_um, then pairwise t-tests across groups
    res <- tryCatch({
      fit_res <- lm(value ~ islet_diam_um, data = rdf)
      r <- resid(fit_res)
      pt <- pairwise.t.test(r, rdf$donor_status, p.adjust.method = "BH")
      mat <- as.data.frame(as.table(pt$p.value))
      colnames(mat) <- c("group1","group2","p_value")
      mat <- mat[!is.na(mat$p_value), , drop = FALSE]
      mat
    }, error = function(e) NULL)
    # Assemble results
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

  output$stats_tbl <- renderTable({
    st <- stats_run()
    if (is.null(st) || nrow(st) == 0) return(NULL)
    alpha_num <- as.numeric(input$alpha)
    # Adjust only pairwise; keep global unadjusted
    st$p_adj <- NA_real_
    if (any(st$type == "pairwise")) {
      idx <- which(st$type == "pairwise")
      st$p_adj[idx] <- p.adjust(st$p_value[idx], method = "BH")
    }
    st$sig <- ifelse(!is.na(st$p_adj), st$p_adj <= alpha_num, st$p_value <= alpha_num)
    st
  })

  output$preview <- renderTable({
    head(plot_df(), 15)
  })

  # Pseudotime scatter when counts are selected (scaled counts vs pseudotime)
  output$pseudo_ui <- renderUI({
    show <- (identical(input$mode, "Targets") && !is.null(input$target_metric) && input$target_metric == "Counts") ||
            (identical(input$mode, "Markers") && !is.null(input$marker_metric) && input$marker_metric == "Counts")
    if (!show) return(NULL)
    tagList(
      tags$h4("Pseudotime scatter (scaled counts)"),
      helpText("Unbiased pseudotime is inferred from available islet features (no group labels), using diffusion map or principal curve when available, else PC1."),
      fluidRow(
        column(4,
          checkboxInput("pseudo_smooth", "Add LOESS smooth", value = TRUE)
        ),
        column(4,
          checkboxInput("pseudo_norm_robust", "Robust per-donor normalization", value = TRUE)
        ),
        column(4,
          checkboxInput("pseudo_clip", "Clip color range (2–98%)", value = TRUE)
        )
      ),
      fluidRow(
        column(6,
          selectInput("pseudo_corr", "Correlation vs pseudotime", choices = c("Kendall", "Spearman", "None"), selected = "Kendall")
        )
      ),
      tableOutput("pseudo_stats")
    )
  })

  output$pseudo <- renderPlotly({
    show <- (identical(input$mode, "Targets") && !is.null(input$target_metric) && input$target_metric == "Counts") ||
            (identical(input$mode, "Markers") && !is.null(input$marker_metric) && input$marker_metric == "Counts")
    if (!show) return(NULL)
    rdf <- raw_df()
    req(nrow(rdf) > 0)
    # Unbiased pseudotime from features (no group labels)
    rdf$pseudotime <- compute_unbiased_pseudotime(rdf)
    # Scale values: robust per-donor normalization (median/MAD) to reduce single-donor outlier impact
    # Use the same scaling as the plot for correlation
    v_raw <- suppressWarnings(as.numeric(rdf$value))
    if (isTRUE(input$pseudo_norm_robust) && all(c("Case ID") %in% colnames(rdf))) {
      rdf$scaled <- v_raw
      rdf <- rdf %>% dplyr::group_by(`Case ID`) %>% dplyr::mutate(
        .med = suppressWarnings(stats::median(scaled, na.rm = TRUE)),
        .mad = suppressWarnings(stats::mad(scaled, center = .med, constant = 1, na.rm = TRUE)),
        .r_sd = ifelse(is.finite(.mad) & .mad > 0, .mad * 1.4826, NA_real_),
        scaled = ifelse(is.finite(.r_sd) & .r_sd > 0, (scaled - .med) / .r_sd, scaled)
      ) %>% dplyr::ungroup() %>% dplyr::select(-.med, -.mad, -`.r_sd`)
      v <- rdf$scaled
    } else {
      mu <- mean(v_raw, na.rm = TRUE)
      sdv <- sd(v_raw, na.rm = TRUE)
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      v <- (v_raw - mu) / sdv
    }
    if (isTRUE(input$pseudo_norm_robust) && all(c("Case ID") %in% colnames(rdf))) {
      rdf$scaled <- v
      # compute per-donor robust z
      rdf <- rdf %>% dplyr::group_by(`Case ID`) %>% dplyr::mutate(
        .med = suppressWarnings(stats::median(scaled, na.rm = TRUE)),
        .mad = suppressWarnings(stats::mad(scaled, center = .med, constant = 1, na.rm = TRUE)),
        .r_sd = ifelse(is.finite(.mad) & .mad > 0, .mad * 1.4826, NA_real_),
        scaled = ifelse(is.finite(.r_sd) & .r_sd > 0, (scaled - .med) / .r_sd, scaled)
      ) %>% dplyr::ungroup() %>% dplyr::select(-.med, -.mad, -`.r_sd`)
    } else {
      mu <- mean(v, na.rm = TRUE)
      sdv <- sd(v, na.rm = TRUE)
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      rdf$scaled <- (v - mu) / sdv
    }
    # color mapping
    col_by <- input$color_by
    if (is.null(col_by) || !(col_by %in% colnames(rdf))) col_by <- "donor_status"
    # When 'Marker' or 'Target class' is chosen, color by count magnitude (continuous)
    if (identical(col_by, "marker") || identical(col_by, "class")) {
      plt <- ggplot(rdf, aes(x = pseudotime, y = scaled, color = value)) +
        geom_point(alpha = input$pt_alpha, size = input$pt_size, position = position_jitter(width = 0.01, height = 0)) +
        scale_x_continuous(limits = c(0,1)) +
        labs(x = "Pseudotime", y = if (isTRUE(input$pseudo_norm_robust)) "Scaled (robust z)" else "Scaled (z)", color = "Count") +
        theme_minimal(base_size = 14) +
        {
          # Apply percentile clipping if requested
          if (isTRUE(input$pseudo_clip)) {
            lim <- stats::quantile(rdf$value, probs = c(0.02, 0.98), na.rm = TRUE)
            ggplot2::scale_color_viridis_c(option = "viridis", end = 0.95, limits = lim, oob = scales::squish)
          } else {
            ggplot2::scale_color_viridis_c(option = "viridis", end = 0.95)
          }
        } +
        guides(color = guide_colorbar(title = "Count"))
    } else {
      plt <- ggplot(rdf, aes(x = pseudotime, y = scaled, color = .data[[col_by]])) +
        geom_point(alpha = input$pt_alpha, size = input$pt_size, position = position_jitter(width = 0.01, height = 0)) +
        scale_x_continuous(limits = c(0,1)) +
        labs(x = "Pseudotime", y = "Scaled counts (z)", color = "Color by") +
        theme_minimal(base_size = 14)
      # Optional: consistent donor_status palette
      if (identical(col_by, "donor_status")) {
        plt <- plt + scale_color_manual(values = c("ND" = "#1f77b4", "Aab+" = "#ff7f0e", "T1D" = "#d62728"))
      }
    }
    # Keep a simple 0..1 axis for unbiased pseudotime
    # optional LOESS smoothing overlay (single global trend)
    if (isTRUE(input$pseudo_smooth)) {
      plt <- plt + geom_smooth(se = FALSE, method = "loess", span = 0.6, color = "#444444")
    }
    if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
      plt <- plt + theme(
        plot.background = element_rect(fill = "#000000", colour = NA),
        panel.background = element_rect(fill = "#000000", colour = NA),
        panel.grid.major = element_line(color = "#333333"),
        panel.grid.minor = element_line(color = "#222222"),
        axis.text = element_text(color = "#e6e6e6"),
        axis.title = element_text(color = "#f0f0f0"),
        plot.title = element_text(color = "#f0f0f0"),
        legend.text = element_text(color = "#e6e6e6"),
        legend.title = element_text(color = "#f0f0f0")
      )
    }
    ggplotly(plt)
  })

  # Correlation summary vs pseudotime for the current subset
  output$pseudo_stats <- renderTable({
    show <- (identical(input$mode, "Targets") && !is.null(input$target_metric) && input$target_metric == "Counts") ||
            (identical(input$mode, "Markers") && !is.null(input$marker_metric) && input$marker_metric == "Counts")
    if (!show) return(NULL)
    method <- input$pseudo_corr
    if (is.null(method) || identical(method, "None")) return(NULL)
    rdf <- raw_df()
    if (is.null(rdf) || !nrow(rdf)) return(NULL)
    # compute unbiased pseudotime same as in the plot
    pt <- compute_unbiased_pseudotime(rdf)
    v <- suppressWarnings(as.numeric(rdf$value))
    # choose method
    m <- if (identical(method, "Spearman")) "spearman" else "kendall"
    ct <- tryCatch(cor.test(pt, v, method = m, exact = FALSE), error = function(e) NULL)
    if (is.null(ct)) return(NULL)
    stat_name <- if (identical(method, "Spearman")) "rho" else "tau"
    data.frame(
      method = method,
      estimate = unname(ct$estimate),
      statistic = stat_name,
      p_value = unname(ct$p.value),
      stringsAsFactors = FALSE
    )
  })

  output$dl_summary <- downloadHandler(
    filename = function() {
      paste0("summary_", gsub("[^0-9A-Za-z]+","_", Sys.time()), ".csv")
    },
    content = function(file) {
      df <- summary_df()
      if (is.null(df) || nrow(df) == 0) df <- data.frame()
      write.csv(df, file, row.names = FALSE)
    }
  )

  output$dl_stats <- downloadHandler(
    filename = function() {
      paste0("stats_", gsub("[^0-9A-Za-z]+","_", Sys.time()), ".csv")
    },
    content = function(file) {
      st <- stats_run()
      if (is.null(st) || nrow(st) == 0) st <- data.frame()
      write.csv(st, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
