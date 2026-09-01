test_that("trajectory extent keeps jitter and panel limits synchronized", {
  ext <- traj_pt_extent(c(0, 10), jitter_frac = 0.02)
  expect_equal(ext$dom, c(0, 10))
  expect_equal(ext$jitter, 0.2)
  expect_equal(ext$lim, c(-0.2, 10.2))

  fallback <- traj_pt_extent(c(NA_real_, NA_real_))
  expect_equal(fallback$dom, c(0, 1))
})

test_that("pseudotime scale helper preserves raw values and stable ranks", {
  df <- data.frame(pt = c(4, 2, 8), pt_rank = c(0.5, 0.2, 0.9))
  raw <- apply_pt_scale(df, "raw")
  ranked <- apply_pt_scale(df, "rank")

  expect_equal(raw$pt, c(4, 2, 8))
  expect_equal(raw$pt_raw, c(4, 2, 8))
  expect_equal(ranked$pt, c(0.5, 0.2, 0.9))
  expect_equal(ranked$pt_raw, c(4, 2, 8))
})
