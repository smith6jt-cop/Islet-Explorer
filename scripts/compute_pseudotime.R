#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr)
  library(SingleCellExperiment); library(SummarizedExperiment)
  library(slingshot); library(destiny); library(cluster)
  library(jsonlite)
})

# Load master workbook
master_path <- file.path("data","master_results.xlsx")
safe_read_sheet <- function(path, sheet) tryCatch(readxl::read_excel(path, sheet = sheet, guess_max = 100000), error = function(e) NULL)
add_islet_key <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (!"islet_key" %in% names(df)) df$islet_key <- NA_character_
  if ("region" %in% names(df)) {
    key_from_region <- stringr::str_extract(df$region, "Islet_\\d+")
    df$islet_key <- dplyr::coalesce(df$islet_key, key_from_region)
  }
  if ("name" %in% names(df)) {
    key_from_name <- stringr::str_extract(df$name, "Islet_\\d+")
    only_digits  <- stringr::str_extract(df$name, "\\d+")
    fallback_name <- ifelse(!is.na(only_digits), paste0("Islet_", only_digits), NA_character_)
    df$islet_key <- dplyr::coalesce(df$islet_key, key_from_name, fallback_name)
  }
  if ("islet_id" %in% names(df)) {
    id_str <- suppressWarnings(as.character(df$islet_id))
    id_digits <- stringr::str_extract(id_str, "\\d+")
    fallback_id <- ifelse(!is.na(id_digits), paste0("Islet_", id_digits), NA_character_)
    df$islet_key <- dplyr::coalesce(df$islet_key, fallback_id)
  }
  df
}
compute_diameter_um <- function(area_um2) ifelse(is.finite(as.numeric(area_um2)) & as.numeric(area_um2) > 0, 2*sqrt(as.numeric(area_um2)/pi), NA_real_)

message("Reading master workbook: ", master_path)
markers <- safe_read_sheet(master_path, "Islet_Markers") %>% add_islet_key()
targets <- safe_read_sheet(master_path, "Islet_Targets") %>% add_islet_key()
comp    <- safe_read_sheet(master_path, "Islet_Composition") %>% add_islet_key()
lgals3  <- safe_read_sheet(master_path, "LGALS3") %>% add_islet_key()

# Compute core-derived diameter from targets core region area
core_area <- targets %>% filter(tolower(type) == "islet_core") %>% select(`Case ID`,`Donor Status`, islet_key, core_region_um2 = region_um2) %>% distinct()
size_area <- core_area %>% mutate(islet_diam_um = compute_diameter_um(core_region_um2)) %>% select(`Case ID`,`Donor Status`, islet_key, islet_diam_um)

# Build feature table per islet
base <- comp %>% filter(!is.na(islet_key)) %>% select(`Case ID`,`Donor Status`, islet_key, cells_total, Ins_any, Glu_any, Stt_any) %>% left_join(size_area, by = c("Case ID","Donor Status","islet_key"))
base <- base %>% mutate(
  comp_frac_INS = ifelse(cells_total > 0, as.numeric(Ins_any)/as.numeric(cells_total), NA_real_),
  comp_frac_GCG = ifelse(cells_total > 0, as.numeric(Glu_any)/as.numeric(cells_total), NA_real_),
  comp_frac_SST = ifelse(cells_total > 0, as.numeric(Stt_any)/as.numeric(cells_total), NA_real_)
) %>% select(-Ins_any,-Glu_any,-Stt_any)

markers_all <- bind_rows(markers, lgals3) %>% filter(!is.na(islet_key))
mk_union <- markers_all %>% mutate(region_type = tolower(ifelse(is.na(region_type) & !is.na(region), str_extract(region, "(core|band|union)$"), region_type))) %>%
  filter(region_type == "islet_union") %>% select(`Case ID`,`Donor Status`, islet_key, marker, pos_frac)
mk_w <- if (nrow(mk_union)) pivot_wider(mk_union, id_cols = c(`Case ID`,`Donor Status`, islet_key), names_from = marker, values_from = pos_frac, values_fn = list(pos_frac = mean)) else NULL
if (!is.null(mk_w)) colnames(mk_w) <- c(colnames(mk_w)[1:3], paste0("mk_frac_", make.names(colnames(mk_w)[-(1:3)], unique = TRUE)))

tg_union <- targets %>% mutate(type = tolower(type)) %>% filter(type == "islet_union") %>% select(`Case ID`,`Donor Status`, islet_key, class, area_density)
tg_w <- if (nrow(tg_union)) pivot_wider(tg_union, id_cols = c(`Case ID`,`Donor Status`, islet_key), names_from = class, values_from = area_density, values_fn = list(area_density = mean)) else NULL
if (!is.null(tg_w)) colnames(tg_w) <- c(colnames(tg_w)[1:3], paste0("tg_dens_", make.names(colnames(tg_w)[-(1:3)], unique = TRUE)))

ff <- base
if (!is.null(mk_w)) ff <- left_join(ff, mk_w, by = c("Case ID","Donor Status","islet_key"))
if (!is.null(tg_w)) ff <- left_join(ff, tg_w, by = c("Case ID","Donor Status","islet_key"))

# Attach donor demographics if present
donors_meta <- NULL
for (df in list(comp, targets, markers, lgals3)) {
  if (!is.null(df) && nrow(df) > 0) {
    cols <- intersect(colnames(df), c("Case ID","Donor Status","Age","Gender"))
    if (length(cols) >= 2) {
      dm <- df[, cols, drop = FALSE] %>% distinct()
      donors_meta <- if (is.null(donors_meta)) dm else bind_rows(donors_meta, dm)
    }
  }
}
if (!is.null(donors_meta)) donors_meta <- donors_meta %>% distinct()
if (!is.null(donors_meta)) ff <- ff %>% left_join(donors_meta, by = c("Case ID","Donor Status"))

# Fixed bins
fixed_bins <- function(x_um) {
  brks <- c(0,10,20,35,50,100,200, Inf)
  labs <- c("[0,10)", "[10,20)", "[20,35)", "[35,50)", "[50,100)", "[100,200)", ">200")
  cut(as.numeric(x_um), breaks = brks, labels = labs, include.lowest = TRUE, right = FALSE)
}
ff$diam_bin_fixed <- fixed_bins(ff$islet_diam_um)

# Silhouette-based k selection (2..6) using kmeans
select_k <- function(emb, k_min = 2, k_max = 6) {
  best_k <- k_min; best_avg <- -Inf
  for (k in k_min:k_max) {
    set.seed(1)
    km <- kmeans(emb, centers = k)
    sil <- try(silhouette(km$cluster, dist(emb)), silent = TRUE)
    if (inherits(sil, 'try-error')) next
    avg <- mean(sil[, 'sil_width'])
    if (is.finite(avg) && avg > best_avg) { best_avg <- avg; best_k <- k }
  }
  best_k
}

compute_pt <- function(dat) {
  # Build numeric matrix, exclude keys/meta
  key_cols <- c("Case ID","Donor Status","islet_key","islet_diam_um","diam_bin_fixed","cells_total")
  num_cols <- setdiff(colnames(dat)[vapply(dat, is.numeric, logical(1))], key_cols)
  X <- as.matrix(dat[, num_cols, drop = FALSE])
  # Drop all-zero columns
  keep_col <- colSums(abs(X), na.rm = TRUE) > 0
  X <- X[, keep_col, drop = FALSE]
  # Median impute
  for (j in seq_len(ncol(X))) {
    v <- X[, j]
    if (anyNA(v)) { v[is.na(v)] <- suppressWarnings(median(v, na.rm = TRUE)); X[, j] <- v }
  }
  if (nrow(X) < 10 || ncol(X) < 2) return(NULL)
  Xs <- scale(X)
  row_ok <- apply(Xs,1,function(r) all(is.finite(r)))
  Xs <- Xs[row_ok, , drop = FALSE]
  meta <- dat[row_ok, c("Case ID","Donor Status","islet_key","islet_diam_um","diam_bin_fixed"), drop = FALSE]
  # Per-islet scaled feature range (max-min across scaled features)
  meta$scaled_range <- apply(Xs, 1, function(r) max(r) - min(r))
  # Deduplicate identical rows to avoid DiffusionMap dropping duplicates internally
  # Create a hash per row to map back
  hash <- apply(Xs, 1, function(r) paste(signif(r, 10), collapse = ","))
  u_keep <- !duplicated(hash)
  Xs_u <- Xs[u_keep, , drop = FALSE]
  # Map each row to its unique representative index
  u_index <- match(hash, unique(hash))
  # Compute diffusion map on deduplicated matrix directly
  dmap <- DiffusionMap(Xs_u)
  emb <- eigenvectors(dmap)
  ncomp <- min(ncol(emb), 10)
  emb_use <- as.matrix(emb[, seq_len(ncomp), drop = FALSE])
  # Ensure embedding corresponds to unique rows; if destiny returned full-size, compress to unique
  if (nrow(emb_use) != nrow(Xs_u)) {
    # Assume emb_use corresponds already to Xs_u size if dmap recomputed above
    emb_u <- emb_use
  } else {
    emb_u <- emb_use
  }
  k <- select_k(emb_u)
  set.seed(1); km <- kmeans(emb_u, centers = k)
  # Initialize SCE with ncol equal to number of observations
  sce <- SingleCellExperiment(assays = list(dummy = matrix(0, nrow = 1, ncol = nrow(emb_u))))
  reducedDims(sce)$DM <- emb_u
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(cluster = as.factor(km$cluster))
  ss <- slingshot(sce, clusterLabels = 'cluster', reducedDim = 'DM')
  ptime_mat <- slingPseudotime(ss)
  pt_u <- apply(as.matrix(ptime_mat), 1, function(x) suppressWarnings(min(x, na.rm = TRUE)))
  # Map pseudotime back to all rows using u_index
  pt <- pt_u[u_index]
  pmin <- min(pt, na.rm = TRUE); pmax <- max(pt, na.rm = TRUE)
  if (is.finite(pmin) && is.finite(pmax) && pmax > pmin) pt <- (pt - pmin)/(pmax - pmin)
  list(meta = meta, embedding = emb_u, pseudotime = pt)
}

message("Computing global pseudotime ...")
global <- compute_pt(ff)

bins <- levels(ff$diam_bin_fixed)
bin_res <- list()
for (lb in bins) {
  message("Computing bin ", lb, " ...")
  sub <- ff %>% filter(as.character(diam_bin_fixed) == lb)
  res <- try(compute_pt(sub), silent = TRUE)
  if (!inherits(res, 'try-error') && !is.null(res)) bin_res[[lb]] <- res
}

# Associations helper
compute_associations <- function(dat_feat, res) {
  if (is.null(res)) return(NULL)
  meta <- res$meta
  key_df <- data.frame(islet_key = meta$islet_key, pt = as.numeric(res$pseudotime))
  dd <- dat_feat %>% inner_join(key_df, by = 'islet_key')
  # Numeric features (exclude obvious meta)
  num_feats <- names(dd)[vapply(dd, is.numeric, logical(1))]
  num_feats <- setdiff(num_feats, c('pt','islet_diam_um','cells_total'))
  # Categorical demos
  cat_feats <- intersect(names(dd), c('Gender','Donor Status'))
  rows <- list()
  for (nm in num_feats) {
    v <- dd[[nm]]; keep <- is.finite(v) & is.finite(dd$pt)
    if (sum(keep) >= 10) {
      ct <- suppressWarnings(cor.test(v[keep], dd$pt[keep], method = 'spearman', exact = FALSE))
      rows[[length(rows)+1]] <- data.frame(feature = nm, type = ifelse(grepl('^mk_', nm), 'marker', ifelse(grepl('^tg_', nm), 'target', 'composition/demographic')),
                                           test = 'Spearman', stat = unname(ct$estimate), p_value = unname(ct$p.value))
    }
  }
  for (nm in cat_feats) {
    g <- dd[[nm]]; keep <- !is.na(g) & is.finite(dd$pt)
    if (length(unique(g[keep])) >= 2 && sum(keep) >= 10) {
      kw <- suppressWarnings(kruskal.test(dd$pt[keep] ~ as.factor(g[keep])))
      rows[[length(rows)+1]] <- data.frame(feature = nm, type = 'demographic', test = 'Kruskal', stat = unname(kw$statistic), p_value = unname(kw$p.value))
    }
  }
  out <- bind_rows(rows)
  if (!nrow(out)) return(NULL)
  out$p_adj <- p.adjust(out$p_value, method = 'BH')
  out <- out %>% arrange(p_adj)
  out$stat <- round(as.numeric(out$stat), 4)
  out$p_value <- signif(out$p_value, 3)
  out$p_adj <- signif(out$p_adj, 3)
  out
}

assoc_global <- compute_associations(ff, global)
assoc_bins <- lapply(names(bin_res), function(lb) compute_associations(ff %>% filter(as.character(diam_bin_fixed) == lb), bin_res[[lb]]))
names(assoc_bins) <- names(bin_res)

out <- list(
  global = global,
  bins = bin_res,
  assoc = list(global = assoc_global, bins = assoc_bins)
)

saveRDS(out, file = file.path("data","pseudotime_results.rds"))
message("Saved: data/pseudotime_results.rds")
