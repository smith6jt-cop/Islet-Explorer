# Only curated donor fields already displayed by the Spatial or Inflammation tab.
# The trajectory graph uses none of these columns. Donor ID and disease status
# have their own existing selectors.
TRAJ_DONOR_FIELDS <- c(
  "Age (years)" = "age", "Age group" = "age_group", "Sex" = "gender",
  "Race/ethnicity" = "race_ethnicity", "Height (cm)" = "height_cm",
  "Weight (kg)" = "weight_kg", "BMI" = "BMI", "BMI percentile" = "BMI_percentile",
  "HbA1c (%)" = "HbA1c_percent", "C-peptide (ng/mL)" = "C_peptide_ng_mL",
  "T1D duration (years)" = "t1d_duration_years",
  "Admission glucose (mg/dL)" = "admission_glucose_mg_dL",
  "Peak glucose (mg/dL)" = "peak_glucose_mg_dL",
  "Autoantibodies (number positive)" = "n_aab", "GADA" = "AAb_GADA",
  "IA-2A" = "AAb_IA2A", "mIAA" = "AAb_mIAA", "ZnT8A" = "AAb_ZnT8A",
  "Cause of death" = "cause_of_death", "Recovery type" = "recovery_type",
  "ICU stay (days)" = "ICU_days", "Organ transit time (min)" = "transit_time_min",
  "Pancreas weight (g)" = "pancreas_weight_g", "RNA integrity (RIN)" = "RIN",
  "Blood group (ABO)" = "ABO_group"
)

trajectory_case_id <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !grepl("^[0-9]+$", x)] <- NA_character_
  sub("^0+(?=[0-9])", "", x, perl = TRUE)
}

trajectory_donor_metadata <- function(path = data_path("donor_metadata.csv")) {
  if (!file.exists(path)) return(data.frame(case_id = character(0)))
  d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  for (field in names(d)[vapply(d, is.character, logical(1))]) {
    text <- trimws(d[[field]])
    text[!nzchar(text)] <- NA_character_
    d[[field]] <- text
  }
  if (!"case_id" %in% names(d)) stop("Donor metadata has no case_id")
  d$case_id <- trajectory_case_id(d$case_id)
  if (anyNA(d$case_id) || anyDuplicated(d$case_id)) stop("Invalid or duplicate donor case_id")
  if ("age_years_npod" %in% names(d)) d$age <- d$age_years_npod
  if ("sex" %in% names(d)) d$gender <- d$sex
  if ("donor_type" %in% names(d)) {
    statuses <- c("No Diabetes" = "ND", "Autoantibody Positive" = "Aab+", "Type 1 Diabetes" = "T1D")
    d$donor_status <- unname(statuses[d$donor_type])
  }
  flags <- c(grep("^AAb_", TRAJ_DONOR_FIELDS, value = TRUE), "C_peptide_censored")
  for (field in intersect(flags, names(d))) {
    d[[field]] <- unname(c("TRUE" = TRUE, "FALSE" = FALSE, "1" = TRUE, "0" = FALSE)[
      toupper(as.character(d[[field]]))])
  }
  categories <- c("age_group", "gender", "race_ethnicity", "cause_of_death", "recovery_type", "ABO_group")
  numeric_fields <- setdiff(unname(TRAJ_DONOR_FIELDS), c(categories, flags))
  for (field in intersect(numeric_fields, names(d))) d[[field]] <- suppressWarnings(as.numeric(d[[field]]))
  # Older exports encode unmeasured HbA1c as zero. Zero duration and zero
  # antibodies are valid values and are preserved.
  if ("HbA1c_percent" %in% names(d)) {
    hb <- suppressWarnings(as.numeric(d$HbA1c_percent))
    hb[!is.finite(hb) | hb <= 0] <- NA_real_
    d$HbA1c_percent <- hb
  }
  keep <- c("case_id", "donor_status", TRAJ_DONOR_FIELDS, "C_peptide_censored")
  d[, intersect(keep, names(d)), drop = FALSE]
}

trajectory_join_donor_metadata <- function(df, metadata) {
  if (!"case_id" %in% names(df) || !"case_id" %in% names(metadata)) stop("Missing donor case_id")
  keys <- trajectory_case_id(metadata$case_id)
  if (anyNA(keys) || anyDuplicated(keys)) stop("Invalid or duplicate donor case_id")
  which_donor <- match(trajectory_case_id(df$case_id), keys)
  allowed <- c("donor_status", TRAJ_DONOR_FIELDS, "C_peptide_censored")
  for (field in setdiff(intersect(allowed, names(metadata)), names(df))) {
    df[[field]] <- metadata[[field]][which_donor]
  }
  df
}

trajectory_donor_display <- function(df, field) {
  if (length(field) != 1L || !field %in% c("case_id", "donor_status", TRAJ_DONOR_FIELDS) ||
      !field %in% names(df)) return(rep("Unavailable", nrow(df)))
  v <- df[[field]]
  label <- as.character(v)
  if (startsWith(field, "AAb_")) {
    label <- ifelse(v %in% TRUE, "Positive", "Negative")
  }
  if (identical(field, "C_peptide_ng_mL") && "C_peptide_censored" %in% names(df)) {
    flagged <- df$C_peptide_censored %in% TRUE
    label[flagged] <- paste0("<", label[flagged])
  }
  missing <- is.na(v) | !nzchar(trimws(as.character(v)))
  if (is.numeric(v)) missing <- missing | !is.finite(v)
  label[missing] <- "Unavailable"
  if (identical(field, "t1d_duration_years") && "donor_status" %in% names(df)) {
    label[missing & df$donor_status %in% c("ND", "Aab+")] <- "Not applicable"
  }
  label
}
