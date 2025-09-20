suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
})

master_path <- file.path("data", "master_results.xlsx")

safe_read_sheet <- function(path, sheet) {
  tryCatch(readxl::read_excel(path, sheet = sheet, guess_max = 100000), error = function(e) NULL)
}

load_master <- function(path) {
  list(
    markers = safe_read_sheet(path, "Islet_Markers"),
    targets = safe_read_sheet(path, "Islet_Targets"),
    comp    = safe_read_sheet(path, "Islet_Composition"),
    lgals3  = safe_read_sheet(path, "LGALS3")
  )
}

add_islet_key <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (!"islet_key" %in% names(df)) df$islet_key <- NA_character_
  # region-like columns to try
  region_cols <- intersect(c("region", "Region", "islet_region"), names(df))
  if (length(region_cols) > 0) {
    for (rc in region_cols) {
      key_from_region <- stringr::str_extract(df[[rc]], "Islet_\\d+")
      df$islet_key <- dplyr::coalesce(df$islet_key, key_from_region)
    }
  }
  # name-like columns
  name_cols <- intersect(c("name","Name","islet_name","islet"), names(df))
  if (length(name_cols) > 0) {
    for (nc in name_cols) {
      key_from_name <- stringr::str_extract(df[[nc]], "Islet_\\d+")
      only_digits  <- stringr::str_extract(as.character(df[[nc]]), "\\\\d+")
      fallback_name <- ifelse(!is.na(only_digits), paste0("Islet_", only_digits), NA_character_)
      df$islet_key <- dplyr::coalesce(df$islet_key, key_from_name, fallback_name)
    }
  }
  # id-like columns
  id_cols <- intersect(c("islet_id","Islet ID","Islet_ID","isletID"), names(df))
  if (length(id_cols) > 0) {
    for (ic in id_cols) {
      id_str <- suppressWarnings(as.character(df[[ic]]))
      id_digits <- stringr::str_extract(id_str, "\\\\d+")
      fallback_id <- ifelse(!is.na(id_digits), paste0("Islet_", id_digits), NA_character_)
      df$islet_key <- dplyr::coalesce(df$islet_key, fallback_id)
    }
  }
  df
}

inspect <- function(df, label) {
  if (is.null(df)) {
    cat("[", label, "] sheet missing or unreadable\n", sep="")
    return()
  }
  cat("[", label, "] nrows=", nrow(df), "\n", sep="")
  cat("cols: ", paste(names(df), collapse=", "), "\n", sep="")
  dd <- add_islet_key(df)
  na_rows <- sum(is.na(dd$islet_key))
  cat("islet_key NA:", na_rows, " (", ifelse(nrow(dd)>0, round(100*na_rows/nrow(dd),2), 0), "%)\n", sep="")
  if (na_rows > 0) {
    cat("Examples with NA islet_key (first 5):\n")
    ex <- dd %>% filter(is.na(islet_key)) %>% select(any_of(c("Case ID","Donor Status","region","Region","name","Name","islet_id","Islet ID","Islet_ID"))) %>% head(5)
    print(ex)
  }
  # distinct islet set size by donor
  if (all(c("Case ID","Donor Status") %in% names(dd))) {
    di <- dd %>% filter(!is.na(islet_key)) %>% distinct(`Case ID`,`Donor Status`, islet_key)
    cat("distinct islets:", nrow(di), "\n")
    per_donor <- di %>% count(`Case ID`,`Donor Status`, name = "n_islets") %>% arrange(`Donor Status`, `Case ID`)
    print(head(per_donor, 10))
  }
}

m <- load_master(master_path)
inspect(m$comp, "comp")
inspect(m$markers, "markers")
inspect(m$targets, "targets")
if (!is.null(m$lgals3)) inspect(m$lgals3, "lgals3")
