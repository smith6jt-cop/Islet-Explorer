suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr); library(tidyr); library(broom)
})

master_path <- file.path("data","master_results.xlsx")

add_islet_key <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if ("region" %in% names(df)) {
    df <- df %>% mutate(islet_key = stringr::str_extract(region, "Islet_\\\
d+"))
  }
  df
}

compute_diameter_um <- function(area_um2) {
  area_um2 <- suppressWarnings(as.numeric(area_um2))
  ifelse(is.finite(area_um2) & area_um2 > 0, 2 * sqrt(area_um2 / pi), NA_real_)
}

bin_islet_sizes <- function(df, diam_col, width) {
  x <- suppressWarnings(as.numeric(df[[diam_col]]))
  bin_lo <- floor(x / width) * width
  bin_hi <- bin_lo + width
  df$diam_mid <- bin_lo + width/2
  df$diam_bin <- factor(paste0("[", bin_lo, ", ", bin_hi, ")"), levels = unique(paste0("[", bin_lo, ", ", bin_hi, ")")[order(bin_lo)]), ordered = TRUE)
  df
}

per_bin_anova <- function(df, bin_col, group_col, value_col, mid_col = "diam_mid") {
  if (nrow(df) == 0) return(tibble(bin = character(), mid = numeric(), p_anova = numeric()))
  bmeta <- df %>% filter(!is.na(.data[[bin_col]])) %>% group_by(.data[[bin_col]]) %>% summarise(mid = suppressWarnings(as.numeric(first(na.omit(.data[[mid_col]])))), .groups = "drop")
  res <- lapply(seq_len(nrow(bmeta)), function(i) {
    b <- bmeta[[bin_col]][i]
    sub <- df %>% filter(.data[[bin_col]] == b)
    if (dplyr::n_distinct(na.omit(sub[[group_col]])) < 2 || sum(is.finite(sub[[value_col]])) < 2) return(NULL)
    fit <- tryCatch(aov(reformulate(group_col, response = value_col), data = sub), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    pval <- tryCatch({ at <- anova(fit); as.numeric(at[["Pr(>F)"]][1]) }, error=function(e) NA_real_)
    tibble(bin = as.character(b), mid = bmeta$mid[i], p_anova = pval)
  })
  out <- bind_rows(res)
  if (nrow(out) == 0) return(out)
  arrange(out, mid)
}

mk <- read_excel(master_path, sheet = "Islet_Markers")
tg <- read_excel(master_path, sheet = "Islet_Targets")
g3 <- read_excel(master_path, sheet = "LGALS3")

core_area <- tg %>% add_islet_key() %>% filter(tolower(type)=="islet_core") %>% select(`Case ID`,`Donor Status`,islet_key, core_region_um2=region_um2) %>% distinct() %>% mutate(islet_diam_um = compute_diameter_um(core_region_um2))

band_targets <- tg %>% add_islet_key() %>% filter(tolower(type)=="islet_band") %>% select(`Case ID`,`Donor Status`,islet_key,class, area_density) %>% left_join(core_area, by=c("Case ID","Donor Status","islet_key"))

band_markers <- mk %>% add_islet_key() %>% filter(tolower(region_type)=="islet_band") %>% select(`Case ID`,`Donor Status`,islet_key, marker, pos_frac) 
if (!is.null(g3) && nrow(g3)>0) {
  g3b <- g3 %>% add_islet_key() %>% filter(tolower(region_type)=="islet_band") %>% select(`Case ID`,`Donor Status`,islet_key, marker, pos_frac)
  band_markers <- bind_rows(band_markers, g3b)
}
band_markers <- band_markers %>% left_join(core_area, by=c("Case ID","Donor Status","islet_key"))

classes <- sort(unique(band_targets$class))
markers <- sort(unique(band_markers$marker))
cat("classes:", length(classes), ":", paste(head(classes,5), collapse=", "), "...\n")
cat("markers:", length(markers), ":", paste(head(markers,5), collapse=", "), "...\n")

df <- band_markers %>% filter(marker == head(markers,1)) %>% filter(is.finite(islet_diam_um)) %>% bin_islet_sizes("islet_diam_um", 50) %>% mutate(value = as.numeric(pos_frac)*100, donor_status = `Donor Status`)
cat("plot rows:", nrow(df), "\n")
sm <- df %>% group_by(donor_status, diam_bin, diam_mid) %>% summarise(n=n(), mean=mean(value, na.rm=TRUE), .groups='drop')
print(head(sm))

st <- per_bin_anova(df, "diam_bin", "donor_status", "value")
print(head(st))
