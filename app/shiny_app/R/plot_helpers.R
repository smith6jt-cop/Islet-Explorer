# Helpers shared by the Plot tab and its data-independent tests.

# Build Cell Populations choices using the retained phenotype vocabulary in
# every region. Legacy residual labels can remain in source data for provenance
# but must not reappear in a selector.
# With `comp` (one row per object), a measure that is nonzero in fewer than PLOT_MIN_NONZERO_FRACTION of the objects
# is left out too (owner, 2026-10-21: plotted, it is a line of points at zero; on v9 core B cell, CD4 Tcell and
# Exocrine, and peri-islet CHGA+ hormone-negative).
PLOT_MIN_NONZERO_FRACTION <- 0.01

#' Fraction of objects with a nonzero value of a prop_* / peri_prop_* column, or of a core+peri pair (`col` a
#' length-2 vector), over the objects where it is measured.
plot_nonzero_fraction <- function(comp, col) {
  vals <- lapply(col, function(k) if (k %in% names(comp)) suppressWarnings(as.numeric(comp[[k]])) else NA_real_)
  measured <- Reduce(`|`, lapply(vals, is.finite))
  nonzero <- Reduce(`|`, lapply(vals, function(v) is.finite(v) & v > 0))
  if (!any(measured)) return(0)
  sum(nonzero & measured) / sum(measured)
}

plot_composition_selection_choices <- function(comp_cols, combined_pairs = NULL, comp = NULL) {
  enough <- function(cols) if (is.null(comp)) rep(TRUE, length(cols)) else
    vapply(cols, function(k) plot_nonzero_fraction(comp, k) >= PLOT_MIN_NONZERO_FRACTION, logical(1))
  prop_cols <- grep("^prop_", comp_cols, value = TRUE)
  prop_cols <- prop_cols[enough(prop_cols)]
  peri_prop_cols <- grep("^peri_prop_", comp_cols, value = TRUE)
  peri_prop_cols <- peri_prop_cols[enough(peri_prop_cols)]
  combined_vals <- if (!is.null(combined_pairs) && nrow(combined_pairs) > 0) {
    pairs_ok <- if (is.null(comp) || !all(c("core_col", "peri_col") %in% names(combined_pairs))) {
      rep(TRUE, nrow(combined_pairs))
    } else {
      vapply(seq_len(nrow(combined_pairs)), function(k) plot_nonzero_fraction(
        comp, c(combined_pairs$core_col[k], combined_pairs$peri_col[k])) >= PLOT_MIN_NONZERO_FRACTION, logical(1))
    }
    setNames(paste0("combined_prop_", combined_pairs$stem), combined_pairs$label)[pairs_ok]
  } else character(0)

  # Each region's phenotypes in their groups (owner, 2026-10-27: phenotype choices "need to be grouped"), as option
  # groups "<region> — <group>" (phenotype_grouped_choices()); without that helper one group per region.
  grouped <- get0("phenotype_grouped_choices", mode = "function", ifnotfound = NULL)
  add_region <- function(choices, region, values, labels) {
    if (is.null(grouped)) {
      choices[[paste(region, "— Phenotype Proportions")]] <- setNames(values, phenotype_text(labels))
      return(choices)
    }
    c(choices, grouped(labels, values, prefix = paste(region, "— ")))
  }
  choices <- list()
  if (length(prop_cols)) {
    labels <- sub("^prop_", "", prop_cols)
    keep <- phenotype_is_kept(labels) & !phenotype_is_unclassified(labels)
    if (any(keep)) choices <- add_region(choices, "Core", prop_cols[keep], labels[keep])
  }
  if (length(peri_prop_cols)) {
    labels <- vapply(peri_prop_cols, function(col) {
      extract_phenotype_from_column(col) %||% col
    }, character(1), USE.NAMES = FALSE)
    keep <- phenotype_is_kept(labels) & !phenotype_is_unclassified(labels)
    if (any(keep)) choices <- add_region(choices, "Peri-Islet", peri_prop_cols[keep], labels[keep])
  }
  if (length(combined_vals)) {
    keep <- phenotype_is_kept(names(combined_vals)) & !phenotype_is_unclassified(names(combined_vals))
    if (any(keep)) choices <- add_region(choices, "Core+Peri", unname(combined_vals[keep]), names(combined_vals)[keep])
  }
  choices
}
