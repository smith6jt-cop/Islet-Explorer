# Helpers shared by the Plot tab and its data-independent tests.

# Build Cell Populations choices using the retained phenotype vocabulary in
# every region. Legacy residual labels can remain in source data for provenance
# but must not reappear in a selector.
plot_composition_selection_choices <- function(comp_cols, combined_pairs = NULL) {
  prop_cols <- grep("^prop_", comp_cols, value = TRUE)
  peri_prop_cols <- grep("^peri_prop_", comp_cols, value = TRUE)
  combined_vals <- if (!is.null(combined_pairs) && nrow(combined_pairs) > 0) {
    setNames(paste0("combined_prop_", combined_pairs$stem), combined_pairs$label)
  } else character(0)

  choices <- list()
  if (length(prop_cols)) {
    labels <- sub("^prop_", "", prop_cols)
    keep <- phenotype_is_kept(labels)
    if (any(keep)) choices[["Core — Phenotype Proportions"]] <-
      setNames(prop_cols[keep], labels[keep])
  }
  if (length(peri_prop_cols)) {
    labels <- vapply(peri_prop_cols, function(col) {
      extract_phenotype_from_column(col) %||% col
    }, character(1), USE.NAMES = FALSE)
    keep <- phenotype_is_kept(labels)
    if (any(keep)) choices[["Peri-Islet — Phenotype Proportions"]] <-
      setNames(peri_prop_cols[keep], labels[keep])
  }
  if (length(combined_vals)) {
    keep <- phenotype_is_kept(names(combined_vals))
    if (any(keep)) choices[["Core+Peri — Phenotype Proportions"]] <- combined_vals[keep]
  }
  choices
}
