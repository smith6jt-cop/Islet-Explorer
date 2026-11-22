source(.find_source("R/trajectory_donor_metadata.R"), local = TRUE)

test_that("trajectory joins curated donor metadata by normalized case ID without multiplying rows", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("case_id,donor_type,age_years_npod,sex,C_peptide_ng_mL,C_peptide_censored,AAb_GADA,t1d_duration_years,narrative",
               "0112,No Diabetes,3.72,Female,2.5,False,False,,ignore",
               "6534,Type 1 Diabetes,12,Male,0.02,True,True,0,ignore"), path)
  m <- trajectory_donor_metadata(path)
  expect_equal(m$case_id, c("112", "6534"))
  expect_false("narrative" %in% names(m))
  d <- data.frame(case_id = c("6534", "0112", "not-an-ID", "6534"), age = c(12, 3, NA, 12))
  joined <- trajectory_join_donor_metadata(d, m)
  expect_identical(joined[names(d)], d)
  expect_equal(joined$gender, c("Male", "Female", NA, "Male"))
  expect_equal(trajectory_donor_display(joined, "C_peptide_ng_mL"), c("<0.02","2.5","Unavailable","<0.02"))
  expect_equal(trajectory_donor_display(joined, "AAb_GADA"), c("Positive","Negative","Unavailable","Positive"))
  expect_equal(trajectory_donor_display(joined, "t1d_duration_years"), c("0","Not applicable","Unavailable","0"))
  expect_error(trajectory_join_donor_metadata(d, rbind(m, m[1, ])), "duplicate")
  expect_equal(trajectory_join_donor_metadata(d, trajectory_donor_metadata(tempfile())), d)
})

test_that("historical HbA1c sentinels and blank metadata remain unavailable", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("case_id,sex,HbA1c_percent,n_aab,t1d_duration_years",
               "1, ,0,0,0",
               "2,Female,-1,1,1",
               "3,Male,5.1,2,2"), path)
  d <- trajectory_donor_metadata(path)
  expect_identical(d$HbA1c_percent, c(NA_real_, NA_real_, 5.1))
  expect_true(is.na(d$gender[1]))
  expect_equal(trajectory_donor_display(d, "gender")[1], "Unavailable")
  expect_equal(d$n_aab[1], 0)
  expect_equal(d$t1d_duration_years[1], 0)
})
