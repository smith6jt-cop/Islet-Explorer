test_that("destiny fallback returns the same root distances without fitting branch assignments", {
  skip_if_not_installed("destiny")
  env <- new.env(parent = globalenv())
  sys.source(.find_source("R/barlow_pseudotime.R"), envir = env)
  set.seed(17)
  dm <- destiny::DiffusionMap(matrix(rnorm(600), ncol = 5), k = 15, n_eigs = 8)
  root <- 7L
  expected <- destiny::DPT(dm, tips = root)[root, ]
  actual <- env$barlow_root_distances(destiny::eigenvectors(dm), destiny::eigenvalues(dm), root)
  expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-12)
  expect_equal(actual[[root]], 0)
})
