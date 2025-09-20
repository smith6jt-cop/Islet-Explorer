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
  if ("region" %in% names(df)) {
    key_from_region <- stringr::str_extract(df$region, "Islet_\\d+")
    df$islet_key <- dplyr::coalesce(df$islet_key, key_from_region)
  }
  if ("name" %in% names(df)) {
    key_from_name <- stringr::str_extract(df$name, "Islet_\\d+")
    only_digits  <- stringr::str_extract(df$name, "\\d+")
    fallback_name <- ifelse(!is.na(only_digits), paste0("Islet_", only_digits), NA_character_)
    df$islet_key <- dplyr::coalesce(df$islet_key, key_from_name, fallback_name)
  }
  if ("islet_id" %in% names(df)) {
    id_str <- suppressWarnings(as.character(df$islet_id))
    id_digits <- stringr::str_extract(id_str, "\\d+")
    fallback_id <- ifelse(!is.na(id_digits), paste0("Islet_", id_digits), NA_character_)
    df$islet_key <- dplyr::coalesce(df$islet_key, fallback_id)
  }
  df
}

compute_diameter_um <- function(area_um2) {
  area_um2 <- suppressWarnings(as.numeric(area_um2))
  ifelse(is.finite(area_um2) & area_um2 > 0, 2 * sqrt(area_um2 / pi), NA_real_)
}

prep_data <- function(master) {
  targets <- master$targets %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key))
  core_area <- targets %>%
    dplyr::filter(tolower(type) == "islet_core") %>%
    dplyr::select(`Case ID`, `Donor Status`, islet_key, core_region_um2 = region_um2) %>%
    dplyr::distinct()
  size_area <- core_area %>%
    dplyr::mutate(islet_diam_um = compute_diameter_um(core_region_um2)) %>%
    dplyr::select(`Case ID`, `Donor Status`, islet_key, islet_diam_um)

  targets_all <- targets %>%
    dplyr::select(`Case ID`, `Donor Status`, islet_key, type, class, area_um2, region_um2, area_density, count) %>%
    dplyr::mutate(
      type = dplyr::case_when(
        tolower(type) %in% c("islet_core", "core") ~ "islet_core",
        tolower(type) %in% c("islet_band", "band", "peri-islet", "peri_islet") ~ "islet_band",
        tolower(type) %in% c("islet_union", "union", "islet+20um", "islet_20um") ~ "islet_union",
        TRUE ~ tolower(type)
      )
    ) %>%
    dplyr::left_join(size_area, by = c("Case ID", "Donor Status", "islet_key"))

  if (nrow(targets_all) > 0) {
    keys <- c("Case ID", "Donor Status", "islet_key", "class")
    have_union <- targets_all %>%
      dplyr::filter(type == "islet_union") %>%
      dplyr::select(dplyr::all_of(keys)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(.has_union = TRUE)

    core_rows <- targets_all %>% dplyr::filter(type == "islet_core") %>% dplyr::select(dplyr::all_of(c(keys, "count"))) %>% dplyr::rename(count_core = count)
    band_rows <- targets_all %>% dplyr::filter(type == "islet_band") %>% dplyr::select(dplyr::all_of(c(keys, "count"))) %>% dplyr::rename(count_band = count)
    union_missing <- core_rows %>%
      dplyr::inner_join(band_rows, by = keys) %>%
      dplyr::left_join(have_union, by = keys) %>%
      dplyr::filter(is.na(.has_union))
    if (nrow(union_missing) > 0) {
      synth <- union_missing %>%
        dplyr::transmute(`Case ID`, `Donor Status`, islet_key, type = "islet_union", class,
                          area_um2 = NA_real_, region_um2 = NA_real_, area_density = NA_real_,
                          count = suppressWarnings(as.numeric(count_core)) + suppressWarnings(as.numeric(count_band))) %>%
        dplyr::left_join(size_area, by = c("Case ID", "Donor Status", "islet_key"))
      targets_all <- dplyr::bind_rows(targets_all, synth)
    }
  }

  markers <- master$markers %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key))
  markers_all <- markers %>%
    dplyr::select(`Case ID`, `Donor Status`, islet_key, region_type, marker, n_cells, pos_count, pos_frac) %>%
    dplyr::mutate(
      region_type = dplyr::case_when(
        tolower(region_type) %in% c("islet_core", "core") ~ "islet_core",
        tolower(region_type) %in% c("islet_band", "band", "peri-islet", "peri_islet") ~ "islet_band",
        tolower(region_type) %in% c("islet_union", "union", "islet+20um", "islet_20um") ~ "islet_union",
        TRUE ~ tolower(region_type)
      )
    )
  if (!is.null(master$lgals3) && nrow(master$lgals3) > 0) {
    g3 <- master$lgals3 %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key)) %>%
      dplyr::mutate(marker = as.character(marker)) %>%
      dplyr::select(`Case ID`, `Donor Status`, islet_key, region_type, marker, n_cells, pos_count, pos_frac)
    if (!is.null(g3) && nrow(g3) > 0) {
      markers_all <- dplyr::bind_rows(markers_all, g3)
    }
  }
  markers_all <- markers_all %>% dplyr::left_join(size_area, by = c("Case ID", "Donor Status", "islet_key"))

  if (nrow(markers_all) > 0) {
    mkeys <- c("Case ID", "Donor Status", "islet_key", "marker")
    have_union_m <- markers_all %>%
      dplyr::filter(region_type == "islet_union") %>%
      dplyr::select(dplyr::all_of(mkeys)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(.has_union = TRUE)

    core_m <- markers_all %>% dplyr::filter(region_type == "islet_core") %>% dplyr::select(dplyr::all_of(c(mkeys, "n_cells", "pos_count"))) %>% dplyr::rename(n_core = n_cells, pos_core = pos_count)
    band_m <- markers_all %>% dplyr::filter(region_type == "islet_band") %>% dplyr::select(dplyr::all_of(c(mkeys, "n_cells", "pos_count"))) %>% dplyr::rename(n_band = n_cells, pos_band = pos_count)
    union_m <- markers_all %>% dplyr::filter(region_type == "islet_union") %>% dplyr::select(dplyr::all_of(c(mkeys, "n_cells", "pos_count"))) %>% dplyr::rename(n_union = n_cells, pos_union = pos_count)
    union_missing_m <- core_m %>%
      dplyr::inner_join(band_m, by = mkeys) %>%
      dplyr::left_join(have_union_m, by = mkeys) %>%
      dplyr::filter(is.na(.has_union))
    if (nrow(union_missing_m) > 0) {
      synth_m <- union_missing_m %>%
        dplyr::transmute(`Case ID`, `Donor Status`, islet_key, region_type = "islet_union", marker,
                          n_cells = suppressWarnings(as.numeric(n_core)) + suppressWarnings(as.numeric(n_band)),
                          pos_count = suppressWarnings(as.numeric(pos_core)) + suppressWarnings(as.numeric(pos_band))) %>%
        dplyr::mutate(pos_frac = ifelse(is.finite(n_cells) & n_cells > 0,
                                        suppressWarnings(as.numeric(pos_count)) / suppressWarnings(as.numeric(n_cells)),
                                        NA_real_)) %>%
        dplyr::left_join(size_area, by = c("Case ID", "Donor Status", "islet_key"))
      markers_all <- dplyr::bind_rows(markers_all, synth_m)
    }

    # Backfill missing band rows when core and union exist but band is missing: band = union - core
    have_band_m <- markers_all %>%
      dplyr::filter(region_type == "islet_band") %>%
      dplyr::select(dplyr::all_of(mkeys)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(.has_band = TRUE)
    band_missing_m <- core_m %>%
      dplyr::inner_join(union_m, by = mkeys) %>%
      dplyr::left_join(have_band_m, by = mkeys) %>%
      dplyr::filter(is.na(.has_band))
    if (nrow(band_missing_m) > 0) {
      synth_band <- band_missing_m %>%
        dplyr::transmute(`Case ID`, `Donor Status`, islet_key, region_type = "islet_band", marker,
                          n_cells = suppressWarnings(as.numeric(n_union)) - suppressWarnings(as.numeric(n_core)),
                          pos_count = suppressWarnings(as.numeric(pos_union)) - suppressWarnings(as.numeric(pos_core))) %>%
        dplyr::mutate(
          n_cells = ifelse(is.finite(n_cells), n_cells, 0),
          pos_count = ifelse(is.finite(pos_count), pos_count, 0),
          n_cells = pmax(0, n_cells),
          pos_count = pmax(0, pos_count),
          pos_frac = ifelse(n_cells > 0, pos_count / n_cells, NA_real_)
        ) %>%
        dplyr::left_join(size_area, by = c("Case ID", "Donor Status", "islet_key"))
      markers_all <- dplyr::bind_rows(markers_all, synth_band)
    }
  }

  if (nrow(markers_all) > 0) {
    markers_all <- markers_all %>%
      dplyr::mutate(
        .n = suppressWarnings(as.numeric(n_cells)),
        .p = suppressWarnings(as.numeric(pos_count)),
        pos_frac = dplyr::coalesce(pos_frac, ifelse(is.finite(.n) & .n > 0 & is.finite(.p), .p/.n, NA_real_)),
        pos_pct = ifelse(is.finite(pos_frac), 100.0 * pos_frac, NA_real_)
      ) %>%
      dplyr::select(-.n, -.p)
  }

  comp <- master$comp %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key))
  comp <- comp %>% dplyr::select(`Case ID`, `Donor Status`, islet_key, cells_total, Ins_single, Glu_single, Stt_single,
                          Multi_Pos, Triple_Neg, Ins_any, Glu_any, Stt_any)
  comp <- comp %>% dplyr::left_join(size_area, by = c("Case ID", "Donor Status", "islet_key"))

  list(core_area = size_area, targets_all = targets_all, markers_all = markers_all, comp = comp)
}

audit_na <- function(df, label) {
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  na_cnt <- vapply(df, function(x) sum(is.na(x)), integer(1))
  total <- nrow(df)
  pct <- ifelse(total > 0, round(100 * na_cnt / total, 2), 0)
  cat(paste0("[NA audit] ", label, ": ", paste(names(na_cnt), paste0(na_cnt, " (", pct, "%)"), sep = "=", collapse = "; "), "\n"))
}

master <- load_master(master_path)
pd <- prep_data(master)

cat("Rows: targets_all=", nrow(pd$targets_all), ", markers_all=", nrow(pd$markers_all), ", comp=", nrow(pd$comp), "\n", sep="")
audit_na(pd$targets_all, "targets_all")
audit_na(pd$markers_all, "markers_all")
audit_na(pd$comp, "comp")

# Distinct region rows vs expected 3x comp
di <- pd$comp %>% dplyr::distinct(`Case ID`, `Donor Status`, islet_key)
mk_regions <- pd$markers_all %>% dplyr::distinct(`Case ID`, `Donor Status`, islet_key, region_type)
tg_regions <- pd$targets_all %>% dplyr::distinct(`Case ID`, `Donor Status`, islet_key, type)
cat("Distinct islets:", nrow(di), " Expected 3x:", 3 * nrow(di), "\n")
cat("Markers distinct region rows:", nrow(mk_regions), "\n")
cat("Targets distinct region rows:", nrow(tg_regions), "\n")

per_d <- di %>% dplyr::count(`Case ID`, `Donor Status`, name = "n_islets")
mr_d <- mk_regions %>% dplyr::count(`Case ID`, `Donor Status`, name = "n_mk_regions")
tg_d <- tg_regions %>% dplyr::count(`Case ID`, `Donor Status`, name = "n_tg_regions")
dd <- per_d %>% dplyr::left_join(mr_d, by = c("Case ID","Donor Status")) %>%
  dplyr::left_join(tg_d, by = c("Case ID","Donor Status")) %>%
  dplyr::mutate(exp = 3 * n_islets, ok_m = n_mk_regions == exp, ok_t = n_tg_regions == exp)
cat("Total donors:", nrow(dd), "\n")
cat("Donors with mismatch (either markers or targets):", sum(!(dd$ok_m & dd$ok_t), na.rm = TRUE), "\n")
print(dd %>% dplyr::arrange(`Donor Status`, `Case ID`))
if (any(!(dd$ok_m & dd$ok_t))) {
  cat("\nMismatched donors detail:\n")
  print(dd %>% dplyr::filter(!(ok_m & ok_t)) %>% dplyr::arrange(`Donor Status`, `Case ID`))
}
