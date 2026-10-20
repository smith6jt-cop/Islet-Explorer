# ---------- Islet immune measures for the Spatial tab (owner review, 2026-10-20) ----------
# The owner asked that each Spatial card answer its own question. The review found that four cards repeated another
# card or measured something other than their title: Peri enrichment by status was one column of the composition
# card; Peri vs core proportion had immune cells almost never in the core and endocrine cells almost never in the
# band; Distance to nearest cell tracked object radius (rank correlation 0.78-0.86); Cluster Composition plotted the
# endocrine shares the clusters are built from. These helpers compute what replaced them. All are donor-level:
# each donor counts once, whatever its number of objects.

# The immune labels of the uniform root (scripts/uniform/assemble_uniform.py IMMUNE), including the cells left
# without a subtype ("T cell", "Immune"), which every immune total counts.
IMMUNE_LABELS <- c("CD8a Tcell", "CD4 Tcell", "T cell", "B cell", "Macrophage", "Immune")
# The endocrine shares of the composition embedding (scripts/uniform/composition_embedding.py).
COMPOSITION_ENDOCRINE <- c("Beta cell", "Alpha cell", "Delta cell", "CHGA+ hormone-negative")
# Distance bins (um from the object's border; negative = inside) for the immune gradient.
DISTANCE_BIN_BREAKS <- c(-Inf, 0, 20, 40, 60, 100, 200, 400, Inf)
DISTANCE_BIN_LABELS <- c("Core", "0-20", "20-40", "40-60", "60-100", "100-200", "200-400", ">400")
STATUS_ORDER <- c("ND", "Aab+", "T1D")

.donor_status_factor <- function(x) factor(as.character(x), levels = intersect(STATUS_ORDER, unique(as.character(x))))

#' Kruskal-Wallis p across donor statuses on donor-level values (NA when it cannot be computed).
donor_kw_p <- function(value, status) {
  ok <- is.finite(value)
  if (length(unique(status[ok])) < 2) return(NA_real_)
  tryCatch(stats::kruskal.test(value[ok], factor(status[ok]))$p.value, error = function(e) NA_real_)
}

#' Insulitis prevalence: each donor's share of objects with `threshold` or more immune cells in the core and the
#' 20 um band (all immune labels, subtyped or not). A cell-count criterion in the spirit of nPOD's 15 CD45+ cells.
insulitis_prevalence <- function(comp, threshold = 15) {
  need <- c("Case ID", "Donor Status", "immune_count_core", "immune_count_peri")
  if (is.null(comp) || !all(need %in% names(comp))) return(NULL)
  n_imm <- comp$immune_count_core + comp$immune_count_peri
  ok <- is.finite(n_imm)
  if (!any(ok)) return(NULL)
  d <- data.frame(case_id = as.character(comp$`Case ID`[ok]), status = as.character(comp$`Donor Status`[ok]),
                  hit = n_imm[ok] >= threshold, stringsAsFactors = FALSE)
  out <- do.call(rbind, lapply(split(d, d$case_id), function(g)
    data.frame(case_id = g$case_id[1], status = g$status[1], objects = nrow(g), insulitic = sum(g$hit),
               stringsAsFactors = FALSE)))
  out$value <- out$insulitic / out$objects
  out$status <- .donor_status_factor(out$status)
  out[order(out$status, out$case_id), , drop = FALSE]
}

#' Band immune share of objects that still have beta cells against those that do not, per donor. An object counts
#' as insulin-deficient when beta cells are under `cut` of its core endocrine cells (Beta, Alpha, Delta,
#' CHGA+ hormone-negative). Donors with no object of a kind have no value for it.
immune_by_beta_content <- function(comp, cut = 0.1) {
  props <- paste0("prop_", COMPOSITION_ENDOCRINE)
  need <- c("Case ID", "Donor Status", "immune_frac_peri", props)
  if (is.null(comp) || !all(need %in% names(comp))) return(NULL)
  endo <- rowSums(as.matrix(comp[, props]))
  beta <- comp[["prop_Beta cell"]] / endo
  ok <- is.finite(beta) & is.finite(comp$immune_frac_peri)
  if (!any(ok)) return(NULL)
  d <- data.frame(case_id = as.character(comp$`Case ID`[ok]), status = as.character(comp$`Donor Status`[ok]),
                  kind = ifelse(beta[ok] >= cut, "Beta cells present", "Insulin-deficient"),
                  immune = comp$immune_frac_peri[ok], stringsAsFactors = FALSE)
  out <- do.call(rbind, lapply(split(d, list(d$case_id, d$kind), drop = TRUE), function(g)
    data.frame(case_id = g$case_id[1], status = g$status[1], kind = g$kind[1], value = mean(g$immune),
               objects = nrow(g), stringsAsFactors = FALSE)))
  out$kind <- factor(out$kind, levels = c("Insulin-deficient", "Beta cells present"))
  out$status <- .donor_status_factor(out$status)
  out[order(out$status, out$case_id, out$kind), , drop = FALSE]
}

#' Each Leiden cluster's objects by donor status and by donor (objects outside the embedding have no cluster).
cluster_membership <- function(comp, leiden_col) {
  if (is.null(comp) || !all(c(leiden_col, "Case ID", "Donor Status") %in% names(comp))) return(NULL)
  cl <- as.character(comp[[leiden_col]])
  ok <- !is.na(cl) & nzchar(cl)
  if (!any(ok)) return(NULL)
  d <- data.frame(cluster = cl[ok], status = as.character(comp$`Donor Status`[ok]),
                  case_id = as.character(comp$`Case ID`[ok]), stringsAsFactors = FALSE)
  by_status <- as.data.frame(table(cluster = d$cluster, status = d$status), stringsAsFactors = FALSE)
  names(by_status)[3] <- "objects"
  totals <- tapply(d$cluster, d$cluster, length)
  by_status$share <- by_status$objects / as.numeric(totals[by_status$cluster])
  donors <- do.call(rbind, lapply(split(d, d$cluster), function(g) {
    top <- sort(table(g$case_id), decreasing = TRUE)
    data.frame(cluster = g$cluster[1], objects = nrow(g), donors = length(top), top_donor = names(top)[1],
               top_donor_share = as.numeric(top[1]) / nrow(g), stringsAsFactors = FALSE)
  }))
  levels <- leiden_levels(d$cluster)
  by_status$cluster <- factor(by_status$cluster, levels = levels)
  by_status$status <- .donor_status_factor(by_status$status)
  donors$cluster <- factor(donors$cluster, levels = levels)
  list(by_status = by_status[order(by_status$cluster, by_status$status), , drop = FALSE],
       donors = donors[order(donors$cluster), , drop = FALSE])
}

# ---- Cell counts by donor, phenotype and distance from the nearest endocrine object ----
.distance_profile_cache <- new.env(parent = emptyenv())

#' Counts of every analysed cell by donor, phenotype and distance bin, plus each donor's tissue counts (cells
#' outside every endocrine object and its band). Read once per worker from the cell-distance table, without the
#' retired-label filter of load_cell_distance(): the immune total counts the cells without a subtype. NULL when the
#' table is absent.
donor_distance_profile <- function() {
  if (!is.null(.distance_profile_cache$data)) return(.distance_profile_cache$data)
  if (!cell_distance_available()) return(NULL)
  out <- tryCatch({
    if (.cd_using_parquet()) {
      con <- islet_duckdb_connection()
      rel <- islet_parquet_relation(CELL_DISTANCE_PARQUET_PATH, con)
      cases <- paste0("CASE ", paste(sprintf("WHEN signed_dist < %s THEN '%s'",
                                             DISTANCE_BIN_BREAKS[-1][-length(DISTANCE_BIN_LABELS)],
                                             DISTANCE_BIN_LABELS[-length(DISTANCE_BIN_LABELS)]), collapse = " "),
                      " ELSE '", DISTANCE_BIN_LABELS[length(DISTANCE_BIN_LABELS)], "' END")
      bins <- DBI::dbGetQuery(con, paste0(
        "SELECT CAST(case_id AS VARCHAR) AS case_id, donor_status, phenotype, ", cases, " AS bin, ",
        "count(*) AS n FROM ", rel, " WHERE signed_dist IS NOT NULL GROUP BY 1, 2, 3, 4"))
      tissue <- DBI::dbGetQuery(con, paste0(
        "SELECT CAST(case_id AS VARCHAR) AS case_id, donor_status, phenotype, count(*) AS n FROM ", rel,
        " WHERE cell_region = 'tissue' GROUP BY 1, 2, 3"))
      list(bins = bins, tissue = tissue)
    } else {
      df <- data.table::fread(CELL_DISTANCE_PATH, select = c("case_id", "donor_status", "cell_region",
                                                             "phenotype", "signed_dist"),
                              colClasses = list(character = c("case_id", "donor_status", "cell_region", "phenotype")))
      df <- df[is.finite(df$signed_dist)]
      df$bin <- as.character(cut(df$signed_dist, DISTANCE_BIN_BREAKS, labels = DISTANCE_BIN_LABELS, right = FALSE))
      bins <- as.data.frame(df[, list(n = .N), by = c("case_id", "donor_status", "phenotype", "bin")])
      tissue <- as.data.frame(df[df$cell_region == "tissue", list(n = .N),
                                 by = c("case_id", "donor_status", "phenotype")])
      list(bins = bins, tissue = tissue)
    }
  }, error = function(e) .cd_backend_failure("distance profile", e))
  if (!is.null(out)) {
    out$bins$n <- as.numeric(out$bins$n); out$tissue$n <- as.numeric(out$tissue$n)
    .distance_profile_cache$data <- out
  }
  out
}

#' Each donor's tissue share of the given labels (cells outside every endocrine object and its band).
donor_tissue_share <- function(profile, labels) {
  t <- profile$tissue
  total <- tapply(t$n, t$case_id, sum)
  hit <- tapply(t$n[t$phenotype %in% labels], t$case_id[t$phenotype %in% labels], sum)
  share <- as.numeric(hit[names(total)]) / as.numeric(total)
  share[is.na(share)] <- 0
  stats::setNames(share, names(total))
}

#' Immune gradient: per donor and distance bin, the share of `labels` among all cells in the bin divided by the
#' donor's tissue share of them (1 = the donor's own tissue level). Donors with none in their tissue are dropped.
immune_gradient <- function(profile, labels) {
  if (is.null(profile)) return(NULL)
  b <- profile$bins
  key <- paste(b$case_id, b$bin, sep = "\r")
  total <- tapply(b$n, key, sum)
  hit <- tapply(b$n[b$phenotype %in% labels], key[b$phenotype %in% labels], sum)
  parts <- do.call(rbind, strsplit(names(total), "\r", fixed = TRUE))
  d <- data.frame(case_id = parts[, 1], bin = parts[, 2], cells = as.numeric(total), stringsAsFactors = FALSE)
  d$share <- as.numeric(hit[names(total)]) / d$cells
  d$share[is.na(d$share)] <- 0
  ref <- donor_tissue_share(profile, labels)
  d$tissue_share <- as.numeric(ref[d$case_id])
  d <- d[is.finite(d$tissue_share) & d$tissue_share > 0, , drop = FALSE]
  d$value <- d$share / d$tissue_share
  status <- unique(b[, c("case_id", "donor_status")])
  d$status <- .donor_status_factor(status$donor_status[match(d$case_id, status$case_id)])
  d$bin <- factor(d$bin, levels = DISTANCE_BIN_LABELS)
  d[order(d$status, d$case_id, d$bin), , drop = FALSE]
}

#' Band enrichment z against the donor's own tissue: (k - f n) / sqrt(f n) with f the donor's tissue share of the
#' phenotype, k its band count and n the band's cells. The cohort version (enrich_z_*) uses the cohort's share.
own_tissue_enrich_z <- function(comp, safe, label, profile) {
  count_col <- paste0("peri_count_", safe)
  if (is.null(profile) || !all(c(count_col, "total_cells_peri", "Case ID") %in% names(comp))) return(NULL)
  f <- donor_tissue_share(profile, label)[as.character(comp$`Case ID`)]
  expected <- as.numeric(f) * comp$total_cells_peri
  z <- (comp[[count_col]] - expected) / sqrt(expected)
  z[!is.finite(z)] <- NA_real_
  z
}
