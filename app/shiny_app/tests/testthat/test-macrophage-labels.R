# Macrophages by marker profile (owner, 2026-10-27): twelve ASCII data labels, shown with superscript plus and minus.
# The helpers are read out of R/00_globals.R by name, so the test needs none of the app's packages or data.
.label_env <- function() {
  env <- new.env(parent = globalenv())
  exprs <- parse(.find_source("R/00_globals.R"), keep.source = FALSE)
  wanted <- c("MACROPHAGE_LABELS", "MACROPHAGE_LABEL_PATTERN", "is_macrophage_label", "phenotype_text",
              "phenotype_choices_named", "extract_phenotype_from_column", "phenotype_label_html",
              "PHENOTYPE_GROUP_ORDER", "PHENOTYPE_GROUP_AGGREGATES", "phenotype_group", "is_phenotype_aggregate",
              "phenotype_aggregate_members", "phenotype_aggregates_for", "phenotype_groups_split",
              "phenotype_grouped_choices", "phenotype_group_checkboxes", "add_phenotype_aggregate_columns")
  for (e in exprs) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% wanted) eval(e, env)
  }
  env$root_has_marker_profile_labels <- function() TRUE
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

test_that("only macrophage profile labels get superscripts; values stay ASCII", {
  env <- .label_env()
  expect_equal(env$phenotype_text(c("CD68+CD163-Macrophage_LGALS3+CD4-", "CHGA+ hormone-negative", "Macrophage", NA)),
               c("CD68\u207aCD163\u207bMacrophage_LGALS3\u207aCD4\u207b", "CHGA+ hormone-negative", "Macrophage", NA))
  ch <- env$phenotype_choices_named(c("Beta cell", "CD68-CD163+Macrophage_LGALS3-CD4+"))
  expect_equal(unname(ch), c("Beta cell", "CD68-CD163+Macrophage_LGALS3-CD4+"))
  expect_equal(names(ch), c("Beta cell", "CD68\u207bCD163\u207aMacrophage_LGALS3\u207bCD4\u207a"))
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
  expect_equal(names(env$phenotype_groups_split(labels)), c("Endocrine", "T cells", "B cells", "Macrophages", "Other"))
  aggs <- env$phenotype_aggregates_for(c(labels, "T cell"))
  expect_setequal(aggs, c("All endocrine cells", "All T cells", "All macrophages"))
  expect_setequal(env$phenotype_aggregate_members("All T cells", c(labels, "T cell")), c("CD8a Tcell", "CD4 Tcell", "T cell"))
  ch <- env$phenotype_grouped_choices(c(aggs, labels), prefix = "Core \u2014 ")
  expect_equal(names(ch)[1], "Core \u2014 Endocrine")
  expect_true(is.list(ch[["Core \u2014 B cells"]]))
  expect_equal(unname(unlist(ch[["Core \u2014 Macrophages"]]))[1], "All macrophages")
  expect_equal(names(ch[["Core \u2014 Macrophages"]])[2], env$phenotype_text(env$MACROPHAGE_LABELS[1]))
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
  expect_match(html, "Macrophage_<wbr>", fixed = TRUE)
})

test_that("aggregate columns sum the members' proportions of the same cells", {
  env <- .label_env()
  comp <- data.frame(`prop_CD8a Tcell` = c(0.1, NA), `prop_CD4 Tcell` = c(0.2, NA), `prop_T cell` = c(0, NA),
                     peri_prop_CD8a_Tcell = c(0.05, 0), peri_prop_CD4_Tcell = c(0, 0), peri_prop_T_cell = c(0, 0.1),
                     check.names = FALSE)
  out <- env$add_phenotype_aggregate_columns(comp)
  expect_equal(out[["prop_All T cells"]], c(0.3, NA))
  expect_equal(out[["peri_prop_All_T_cells"]], c(0.05, 0.1))
})
