# Macrophages by marker profile (owner, 2026-10-27): twelve ASCII data labels, shown with superscript plus and minus.
# The helpers are read out of R/00_globals.R by name, so the test needs none of the app's packages or data.
.label_env <- function() {
  env <- new.env(parent = globalenv())
  exprs <- parse(.find_source("R/00_globals.R"), keep.source = FALSE)
  wanted <- c("MACROPHAGE_LABELS", "MACROPHAGE_LABEL_PATTERN", "is_macrophage_label", "phenotype_text",
              "phenotype_choices_named", "extract_phenotype_from_column", "phenotype_label_html",
              "PHENOTYPE_DISPLAY_NAMES", "MACROPHAGE_MARKERS", "MACROPHAGE_KEY_PATTERN", "macrophage_signs",
              "macrophage_collapse", "phenotype_members", "macrophage_collapse_comp", "remap_macrophage_selection",
              "macrophage_marker_checkboxes", "TOP_MARKER_LABELS", "TOP_MARKER_BASE", "CD45_RULE_LABEL",
              "CD45_TOP_MARKER_LABEL", "TOP_MARKER_ORDER", "top_marker_labels",
              "PHENOTYPE_GROUP_ORDER", "PHENOTYPE_GROUP_AGGREGATES", "phenotype_group", "is_phenotype_aggregate",
              "phenotype_aggregate_members", "phenotype_aggregates_for", "phenotype_groups_split",
              "phenotype_grouped_choices", "phenotype_group_checkboxes", "add_phenotype_aggregate_columns")
  wanted <- c(wanted, "UNCLASSIFIED_KEYS", "UNCLASSIFIED_PARENT", "UNCLASSIFIED_RAW_TO_KEY",
              "unclassified_collapse", "phenotype_collapse", "unclassified_collapse_comp", "PHENOTYPE_RETIRED", "phenotype_retired", "phenotype_keep", "phenotype_is_kept",
              "phenotype_is_conflict", "phenotype_partition_labels")
  for (e in exprs) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% wanted) eval(e, env)
  }
  env$root_has_marker_profile_labels <- function() TRUE
  env$root_has_top_marker_labels <- function() TRUE
  env$root_cd45_is_rule_label <- function() FALSE        # v15: "CD45+" is a top-marker label
  env$root_phenotype_labels <- function() NULL
  env
}

test_that("the twelve macrophage labels match the uniform typing (codex_typing.macrophage_labels)", {
  env <- .label_env()
  expect_length(env$MACROPHAGE_LABELS, 12)
  expect_equal(env$MACROPHAGE_LABELS[c(1, 8, 12)],
               c("CD68+CD163+Macrophage_LGALS3+CD4+", "CD68+CD163-Macrophage_LGALS3-CD4-",
                 "CD68-CD163+Macrophage_LGALS3-CD4-"))
  expect_false(any(grepl("^CD68-CD163-", env$MACROPHAGE_LABELS)))
  expect_true(all(grepl(env$MACROPHAGE_LABEL_PATTERN, env$MACROPHAGE_LABELS)))
  expect_equal(env$is_macrophage_label(c("Macrophage", env$MACROPHAGE_LABELS[2], "CD8a Tcell", NA)),
               c(TRUE, TRUE, FALSE, FALSE))
})

test_that("additive choices contain leaves once while single-metric choices retain parent totals", {
  env <- .label_env()
  leaves <- c("Beta cell", "Alpha cell", "CD4 Tcell", env$MACROPHAGE_LABELS[1:2])
  labels <- c("All endocrine cells", "All T cells", "All macrophages", "Macrophage", leaves,
              "T cell", "Hematopoietic CD45+", "Hematopoietic CD45-", "Conflict B/T", leaves[1])
  expect_identical(env$phenotype_partition_labels(labels), c(leaves, unname(env$UNCLASSIFIED_KEYS["subtype"])))
  html <- as.character(env$phenotype_group_checkboxes("exclusive", labels))
  expect_false(grepl('value="All ', html, fixed = TRUE))
  expect_false(grepl('value="Macrophage"', html, fixed = TRUE))
  expect_false(grepl("Conflict|Hematopoietic", html))
  expect_equal(lengths(regmatches(html, gregexpr('name="exclusive"', html))), length(leaves) + 1L)
  single <- unname(unlist(env$phenotype_grouped_choices(labels)))
  expect_true(all(c("All endocrine cells", "All T cells", "All macrophages") %in% single))
  expect_false(any(c("Macrophage", "T cell", "Conflict B/T") %in% single))
})

test_that("root-bound parent queries include residual source members once without exposing them as leaves", {
  env <- .label_env()
  raw <- c("Beta cell", "Alpha cell", "Delta cell", "CHGA+ hormone-negative", "T cell", "CD4 Tcell", "CD8a Tcell",
           env$MACROPHAGE_LABELS[1:2], "Conflict T/Myeloid")
  env$root_phenotype_labels <- function() raw
  expect_setequal(env$phenotype_members("All endocrine cells"), raw[1:4])
  expect_setequal(env$phenotype_members(c("All T cells", "CD4 Tcell", "T cell")), raw[5:7])
  expect_setequal(env$phenotype_members(c("Macrophage", "All macrophages", env$MACROPHAGE_LABELS[1])), raw[8:9])
  expect_length(env$phenotype_members("Conflict T/Myeloid"), 0)
  expect_true("T cell" %in% env$phenotype_keep(raw))
  expect_false("Conflict T/Myeloid" %in% unlist(env$phenotype_grouped_choices(raw)))
})

test_that("display names put the broad name first, with every marker sign a superscript; values stay ASCII", {
  env <- .label_env()
  expect_equal(env$phenotype_text(c("CD68+CD163-Macrophage_LGALS3+CD4-", "CHGA+ hormone-negative", "Beta cell",
                                    "CD8a Tcell", "CD4 Tcell", "T cell", "B cell", "Macrophage", "Immune", "Exocrine",
                                    "All macrophages", "SMA+", NA)),
               c("Macrophage CD68\u207aCD163\u207bLGALS3\u207aCD4\u207b", "Endocrine CHGA\u207a hormone-negative",
                 "Endocrine Beta", "T Cell CD8a\u207a", "T Cell CD4\u207a",
                 "Unclassified \u2014 Lineage or subtype unresolved", "B Cell", "Macrophage",
                 "Immune", "Exocrine", "Macrophage (all)", "SMA\u207a", NA))
  ch <- env$phenotype_choices_named(c("Beta cell", "CD68-CD163+Macrophage_LGALS3-CD4+"))
  expect_equal(unname(ch), c("Beta cell", "CD68-CD163+Macrophage_LGALS3-CD4+"))
  expect_equal(names(ch), c("Endocrine Beta", "Macrophage CD68\u207bCD163\u207aLGALS3\u207bCD4\u207a"))
})

test_that("a sanitized column name maps back to the macrophage label, underscore kept", {
  env <- .label_env()
  expect_equal(env$extract_phenotype_from_column("peri_prop_CD68plusCD163-Macrophage_LGALS3plusCD4-"),
               "CD68+CD163-Macrophage_LGALS3+CD4-")
  expect_equal(env$extract_phenotype_from_column("enrich_z_CHGAplus_hormone-negative"), "CHGA+ hormone-negative")
  expect_equal(env$extract_phenotype_from_column("min_dist_B_cell"), "B cell")
})

test_that("phenotypes are grouped, each group's aggregate first; one-label groups stay option groups", {
  env <- .label_env()
  labels <- c("Beta cell", "Alpha cell", "CD8a Tcell", "CD4 Tcell", "B cell", env$MACROPHAGE_LABELS[1:2], "Exocrine")
  expect_equal(names(env$phenotype_groups_split(labels)), c("Endocrine", "T Cell", "B Cell", "Macrophage", "Other"))
  aggs <- env$phenotype_aggregates_for(c(labels, "T cell"))
  expect_setequal(aggs, c("All endocrine cells", "All T cells", "All macrophages"))
  expect_setequal(env$phenotype_aggregate_members("All T cells", c(labels, "T cell")), c("CD8a Tcell", "CD4 Tcell"))
  ch <- env$phenotype_grouped_choices(c(aggs, labels), prefix = "Core \u2014 ")
  expect_equal(names(ch)[1], "Core \u2014 Endocrine")
  expect_true(is.list(ch[["Core \u2014 B Cell"]]))
  expect_equal(unname(unlist(ch[["Core \u2014 Macrophage"]]))[1], "All macrophages")
  expect_equal(names(ch[["Core \u2014 Macrophage"]])[1], "Macrophage (all)")
  expect_equal(names(ch[["Core \u2014 Macrophage"]])[2], "Macrophage CD68\u207aCD163\u207aLGALS3\u207aCD4\u207a")
})

test_that("the grouped checkbox input is one Shiny checkbox group with All | None per group and overall", {
  env <- .label_env()
  html <- as.character(env$phenotype_group_checkboxes("x", c("Beta cell", "B cell", env$MACROPHAGE_LABELS[1]),
                                                      selected = "B cell"))
  expect_match(html, 'class="form-group shiny-input-checkboxgroup shiny-input-container pheno-group-checkboxes"', fixed = TRUE)
  expect_equal(lengths(regmatches(html, gregexpr('name="x"', html))), 3L)
  expect_equal(lengths(regmatches(html, gregexpr("pheno-select-all pheno-scope-group", html))), 3L)
  expect_equal(lengths(regmatches(html, gregexpr("pheno-scope-all", html))), 2L)
  expect_match(html, 'value="B cell" checked="checked"', fixed = TRUE)
  expect_match(html, "Macrophage CD68", fixed = TRUE)
})

test_that("aggregate columns sum the members' proportions of the same cells", {
  env <- .label_env()
  comp <- data.frame(`prop_CD8a Tcell` = c(0.1, NA), `prop_CD4 Tcell` = c(0.2, NA), `prop_T cell` = c(0, NA),
                     peri_prop_CD8a_Tcell = c(0.05, 0), peri_prop_CD4_Tcell = c(0, 0), peri_prop_T_cell = c(0, 0.1),
                     check.names = FALSE)
  out <- env$add_phenotype_aggregate_columns(comp)
  expect_equal(out[["prop_All T cells"]], c(0.3, NA))
  expect_equal(out[["peri_prop_All_T_cells"]], c(0.05, 0))
})

test_that("macrophages merge by the checked markers; a merged category expands back to its data labels", {
  env <- .label_env()
  L <- env$MACROPHAGE_LABELS
  expect_equal(unique(env$macrophage_collapse(L)), "Macrophage")
  keys <- env$macrophage_collapse(L, c("CD4", "CD163"))
  expect_setequal(unique(keys), c("Macrophage CD163+CD4+", "Macrophage CD163+CD4-", "Macrophage CD163-CD4+",
                                  "Macrophage CD163-CD4-"))
  expect_equal(length(unique(env$macrophage_collapse(L, env$MACROPHAGE_MARKERS))), 12L)
  expect_equal(env$macrophage_collapse(c("B cell", L[1]), "CD68"), c("B cell", "Macrophage CD68+"))
  # every data label in exactly one category, and a category's members are those that merge into it
  for (k in unique(keys)) expect_setequal(env$phenotype_members(k), L[keys == k])
  expect_setequal(env$phenotype_members("Macrophage"), c("Macrophage", L))
  expect_equal(env$phenotype_text("Macrophage CD163-CD4+"), "Macrophage CD163\u207bCD4\u207a")
  expect_equal(env$remap_macrophage_selection(c("B cell", "Macrophage"), c("B cell", "Macrophage CD4+", "Macrophage CD4-")),
               c("B cell", "Macrophage CD4+", "Macrophage CD4-"))
  expect_equal(env$remap_macrophage_selection("B cell", c("B cell", "Macrophage CD4+")), "B cell")
})

test_that("merged object columns sum the members and recompute the band enrichment from summed counts", {
  env <- .label_env()
  L <- env$MACROPHAGE_LABELS
  safe <- function(l) gsub("+", "plus", gsub(" ", "_", l, fixed = TRUE), fixed = TRUE)
  comp <- data.frame(total_cells_peri = c(100, 50))
  expected <- list(c(2, 1), c(3, 0.5))           # each member's expected band count, per object
  counts <- list(c(4, 0), c(1, 2))
  for (i in 1:2) {
    comp[[paste0("prop_", L[i])]] <- c(0.1, 0.2) * i
    comp[[paste0("peri_prop_", safe(L[i]))]] <- counts[[i]] / comp$total_cells_peri
    comp[[paste0("peri_count_", safe(L[i]))]] <- counts[[i]]
    comp[[paste0("enrich_z_", safe(L[i]))]] <- (counts[[i]] - expected[[i]]) / sqrt(expected[[i]])
    comp[[paste0("min_dist_", safe(L[i]))]] <- c(5, NA) * i
  }
  out <- env$macrophage_collapse_comp(comp)
  expect_false(any(grepl("Macrophage_LGALS3", names(out))))
  expect_equal(out[["prop_Macrophage"]], c(0.3, 0.6))
  expect_equal(out[["peri_count_Macrophage"]], c(5, 2))
  E <- expected[[1]] + expected[[2]]
  expect_equal(out[["enrich_z_Macrophage"]], (c(5, 2) - E) / sqrt(E))
  expect_equal(out[["min_dist_Macrophage"]], c(5, NA))
  split <- env$macrophage_collapse_comp(comp, env$MACROPHAGE_MARKERS)
  expect_equal(split[[paste0("prop_", env$macrophage_collapse(L[1], env$MACROPHAGE_MARKERS))]], c(0.1, 0.2))
})

test_that("mixed raw and profile macrophages retain every merged family and use original member counts", {
  env <- .label_env()
  profile <- env$MACROPHAGE_LABELS[1]
  env$root_phenotype_labels <- function() c("Macrophage", profile)
  families <- c("peri_prop_", "peri_count_", "enrich_z_", "min_dist_")
  safe <- function(x) gsub("+", "plus", gsub(" ", "_", x, fixed = TRUE), fixed = TRUE)
  comp <- data.frame(prop_Macrophage = c(0.1, 0.3, NA), peri_prop_Macrophage = c(0.1, 0.3, NA),
                     peri_count_Macrophage = c(10, 30, NA), enrich_z_Macrophage = c(5 / sqrt(5), 0, NA),
                     min_dist_Macrophage = c(4, NA, NA), source_row = 1:3, check.names = FALSE)
  comp[[paste0("prop_", profile)]] <- c(0.2, 0.1, NA)
  profile_values <- list(c(0.2, 0.1, NA), c(20, 10, NA), c(10 / sqrt(10), 5 / sqrt(5), NA), c(8, 6, NA))
  for (i in seq_along(families)) comp[[paste0(families[i], safe(profile))]] <- profile_values[[i]]
  original <- comp
  out <- env$macrophage_collapse_comp(comp)
  destinations <- c("prop_Macrophage", paste0(families, "Macrophage"))
  expect_true(all(destinations %in% names(out)))
  expect_equal(out$prop_Macrophage, c(0.3, 0.4, NA))
  expect_equal(out$peri_prop_Macrophage, c(0.3, 0.4, NA))
  expect_equal(out$peri_count_Macrophage, c(30, 40, NA))
  # Counts 10 + 20, expected 5 + 10: z = (30 - 15) / sqrt(15).
  expect_equal(out$enrich_z_Macrophage, c(15 / sqrt(15), 5 / sqrt(35), NA), tolerance = 1e-12)
  expect_equal(out$min_dist_Macrophage, c(4, 6, NA))
  expect_identical(out$source_row, original$source_row)
  expect_identical(comp, original)
  expect_identical(env$macrophage_collapse_comp(out), out)
  expect_false(any(grepl("Macrophage_LGALS3", names(out), fixed = TRUE)))
  split <- env$macrophage_collapse_comp(comp, "CD4")
  expect_identical(split[destinations], comp[destinations])
  key <- safe(env$macrophage_collapse(profile, "CD4"))
  for (i in seq_along(families)) expect_equal(split[[paste0(families[i], key)]], profile_values[[i]])
  expect_setequal(env$phenotype_members("All macrophages"), c("Macrophage", profile))
  expect_identical(env$phenotype_partition_labels(c("Macrophage", profile)), profile)
})

test_that("top-marker labels form their own group, last, with a superscript and no broad name", {
  env <- .label_env()
  labels <- c("Beta cell", "Exocrine", "CD45+", "VIM+", "SMA+", "CD4 Tcell")
  groups <- env$phenotype_groups_split(labels)
  expect_equal(names(groups), c("Endocrine", "T Cell", "Other", "Top marker only"))
  expect_equal(groups[["Top marker only"]], c("SMA+", "VIM+", "CD45+"))   # TOP_MARKER_LABELS order
  expect_equal(env$phenotype_text(c("CD45+", "CD8a+", "B3TUBB+")), c("CD45\u207a", "CD8a\u207a", "B3TUBB\u207a"))
  env$root_has_top_marker_labels <- function() FALSE        # the frozen root's SMA+ is smooth muscle, not a top marker
  expect_equal(env$phenotype_group("SMA+"), "Other")
})

test_that("v16: CD45+ is the CD45 rule's class in the Immune group; CD45+ (top marker) stays a top marker", {
  env <- .label_env()
  expect_equal(env$TOP_MARKER_LABELS, c(env$TOP_MARKER_BASE, "CD45+ (top marker)"))
  expect_equal(env$phenotype_group(c("CD45+", "Immune")), c("Top marker only", "Immune"))     # v15
  env$root_cd45_is_rule_label <- function() TRUE                                              # v16
  expect_equal(env$phenotype_group(c("CD45+", "CD45+ (top marker)", "VIM+")),
               c("Immune", "Top marker only", "Top marker only"))
  expect_equal(env$top_marker_labels(), c(env$TOP_MARKER_BASE, "CD45+ (top marker)"))
  groups <- env$phenotype_groups_split(c("CD45+ (top marker)", "CD45+", "B cell", "VIM+"))
  expect_equal(names(groups), c("B Cell", "Immune", "Top marker only"))
  expect_equal(groups[["Top marker only"]], c("VIM+", "CD45+ (top marker)"))
  expect_equal(env$phenotype_text(c("CD45+", "CD45+ (top marker)")), c("CD45\u207a", "CD45\u207a (top marker)"))
  expect_equal(env$extract_phenotype_from_column("enrich_z_CD45plus_(top_marker)"), "CD45+ (top marker)")
  expect_equal(env$extract_phenotype_from_column("prop_CD45+ (top marker)"), "CD45+ (top marker)")
})
