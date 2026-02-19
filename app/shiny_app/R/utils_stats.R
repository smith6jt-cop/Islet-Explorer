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
  bmeta <- df %>% filter(!is.na(.data[[bin_col]])) %>%
    group_by(.data[[bin_col]]) %>%
    summarise(mid = suppressWarnings(as.numeric(first(na.omit(.data[[mid_col]])))), .groups = "drop")
  res <- lapply(seq_len(nrow(bmeta)), function(i) {
    b <- bmeta[[bin_col]][i]
    sub <- df %>% filter(.data[[bin_col]] == b)
    if (n_distinct(na.omit(sub[[group_col]])) < 2 || sum(is.finite(sub[[value_col]])) < 2) return(NULL)
    fit <- tryCatch(aov(reformulate(group_col, response = value_col), data = sub), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    pval <- tryCatch({
      at <- anova(fit)
      as.numeric(at[["Pr(>F)"]][1])
    }, error = function(e) NA_real_)
    tibble(bin = as.character(b), mid = bmeta$mid[i], p_anova = pval)
  })
  out <- bind_rows(res)
  if (nrow(out) == 0) return(tibble(bin = character(), mid = numeric(), p_anova = numeric()))
  arrange(out, mid)
}

per_bin_kendall <- function(df, bin_col, group_col, value_col, mid_col = "diam_mid") {
  if (nrow(df) == 0) return(tibble(bin = character(), mid = numeric(), p_kendall = numeric(), tau = numeric()))
  code_group <- function(x) {
    x <- as.character(x)
    ifelse(x == "ND", 0, ifelse(x == "Aab+", 1, ifelse(x == "T1D", 2, NA_real_)))
  }
  bmeta <- df %>% filter(!is.na(.data[[bin_col]])) %>% group_by(.data[[bin_col]]) %>%
    summarise(mid = suppressWarnings(as.numeric(first(na.omit(.data[[mid_col]])))), .groups = "drop")
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

# ---------- Effect size utilities (Statistics tab) ----------

#' Cohen's d for two-sample comparison with 95% CI (normal approximation)
cohens_d <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  n1 <- length(x); n2 <- length(y)
  if (n1 < 2 || n2 < 2) return(list(d = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
  pooled_sd <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
  if (!is.finite(pooled_sd) || pooled_sd == 0) return(list(d = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
  d <- (mean(x) - mean(y)) / pooled_sd
  se_d <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))
  list(d = d, ci_lo = d - 1.96 * se_d, ci_hi = d + 1.96 * se_d)
}

#' Eta-squared from anova() output: SS_effect / SS_total for the first factor
eta_squared <- function(fit) {
  at <- tryCatch(anova(fit), error = function(e) NULL)
  if (is.null(at)) return(NA_real_)
  ss <- at[["Sum Sq"]]
  if (length(ss) < 2) return(NA_real_)
  ss[1] / sum(ss)
}

#' Non-parametric pairwise Wilcoxon tests with BH correction
#' Returns data.frame with group1, group2, p_value columns
pairwise_wilcox <- function(df, group_col, value_col) {
  pw <- tryCatch(
    pairwise.wilcox.test(df[[value_col]], df[[group_col]], p.adjust.method = "BH", exact = FALSE),
    error = function(e) NULL
  )
  if (is.null(pw)) return(data.frame(group1 = character(), group2 = character(), p_value = numeric(), stringsAsFactors = FALSE))
  mat <- as.data.frame(as.table(pw$p.value))
  colnames(mat) <- c("group1", "group2", "p_value")
  mat <- mat[!is.na(mat$p_value), , drop = FALSE]
  mat$group1 <- as.character(mat$group1)
  mat$group2 <- as.character(mat$group2)
  mat
}
