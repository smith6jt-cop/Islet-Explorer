# ---------- Data loading and wrangling ----------
# Extracted from app.R
# Dependencies: safe_left_join, add_islet_key, compute_diameter_um (from R/utils_safe_join.R)
#               master_path, h5ad_path (from 00_globals.R)
#
# Data loading priority: H5AD (islet_explorer.h5ad) > Excel (master_results.xlsx)
# H5AD contains groovy-derived data stored in .uns by scripts/build_h5ad_for_app.py

safe_read_sheet <- function(path, sheet) {
  # Increase guess_max to reduce type misguesses on sparse columns
  tryCatch(readxl::read_excel(path, sheet = sheet, guess_max = 100000), error = function(e) NULL)
}

load_master <- function(path = master_path) {
  list(
    markers = safe_read_sheet(path, "Islet_Markers"),
    targets = safe_read_sheet(path, "Islet_Targets"),
    comp    = safe_read_sheet(path, "Islet_Composition"),
    lgals3  = safe_read_sheet(path, "LGALS3")
  )
}

# OPTIMIZATION: Consolidate donor metadata extraction
get_donor_metadata <- function(master) {
  # Priority order: composition > targets > markers
  for (sheet in c("comp", "targets", "markers")) {
    df <- master[[sheet]]
    if (!is.null(df) && nrow(df) > 0) {
      required_cols <- c("Case ID", "Donor Status")

      # Determine which AAb columns exist
      aab_candidates <- c("AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_IAA", "AAb_mIAA")
      aab_cols <- intersect(aab_candidates, names(df))

      if (all(required_cols %in% names(df))) {
        return(df %>%
          dplyr::select(dplyr::all_of(c(required_cols, aab_cols))) %>%
          dplyr::distinct())
      }
    }
  }
  NULL
}

#' The core+peri phenotype pair table of `comp` (see prep_data()); rebuilt when the macrophage columns are merged
#' (collapse_prepared_macrophages()).
build_combined_pairs <- function(comp) {
  cn <- colnames(comp)
  peri_cols <- grep("^peri_prop_", cn, value = TRUE)
  if (length(peri_cols) == 0) return(NULL)
  rows <- lapply(peri_cols, function(pc) {
    safe <- sub("^peri_prop_", "", pc)
    candidates <- unique(c(
      gsub("plus", "+", gsub("_", " ", safe), fixed = TRUE),
      gsub("_", " ", safe),
      gsub("plus", "+", safe, fixed = TRUE),
      safe
    ))
    core_col <- NA_character_; label <- safe
    for (cand in candidates) {
      guess <- paste0("prop_", cand)
      if (guess %in% cn) { core_col <- guess; label <- cand; break }
    }
    data.frame(stem = safe, core_col = core_col, peri_col = pc,
               label = label, stringsAsFactors = FALSE)
  })
  out <- dplyr::bind_rows(rows)
  out <- out[!is.na(out$core_col), , drop = FALSE]
  out <- out[order(out$label), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' prepared() with the macrophage labels merged for the checked markers (macrophage_collapse_comp()).
collapse_prepared_macrophages <- function(pd, markers = character(0)) {
  if (is.null(pd) || is.null(pd$comp) || !any(grepl(MACROPHAGE_LABEL_PATTERN, sub("^prop_", "", names(pd$comp))))) return(pd)
  pd$comp <- macrophage_collapse_comp(pd$comp, markers)
  pd$combined_pairs <- build_combined_pairs(pd$comp)
  pd
}

prep_data <- function(master) {
  # Determine which autoantibody columns are available
  aab_cols_targets <- intersect(c("AAb_GADA","AAb_IA2A","AAb_ZnT8A","AAb_IAA","AAb_mIAA"), colnames(master$targets))
  aab_cols_markers <- intersect(c("AAb_GADA","AAb_IA2A","AAb_ZnT8A","AAb_IAA","AAb_mIAA"), colnames(master$markers))
  aab_cols_comp    <- intersect(c("AAb_GADA","AAb_IA2A","AAb_ZnT8A","AAb_IAA","AAb_mIAA"), colnames(master$comp))
  aab_cols_all     <- unique(c(aab_cols_targets, aab_cols_markers, aab_cols_comp))
  # Use new consolidated function
  donors_meta <- get_donor_metadata(master)
  # Islet size proxy per islet corefor diameter
  targets <- master$targets %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key))
  core_area <- targets %>%
    dplyr::filter(tolower(type) == "islet_core") %>%
    dplyr::select(`Case ID`, `Donor Status`, islet_key, core_region_um2 = region_um2) %>%
    dplyr::distinct()
  # Peri-islet (20 µm expansion) annotation area, measured by QuPath alongside core
  peri_area <- targets %>%
    dplyr::filter(tolower(type) == "islet_band") %>%
    dplyr::select(`Case ID`, `Donor Status`, islet_key, peri_region_um2 = region_um2) %>%
    dplyr::distinct()
  # Diameter is computed ONLY from core area; if no core exists, diameter is NA
  size_area <- core_area %>%
    dplyr::mutate(islet_diam_um = compute_diameter_um(core_region_um2)) %>%
    dplyr::select(`Case ID`, `Donor Status`, islet_key, islet_diam_um)
  # Pre-built bundle of per-islet region areas + diameter for downstream joins
  region_areas <- size_area %>%
    { safe_left_join(., core_area, by = c("Case ID", "Donor Status", "islet_key"), context = "region_areas:core_area") } %>%
    { safe_left_join(., peri_area, by = c("Case ID", "Donor Status", "islet_key"), context = "region_areas:peri_area") }

  # Targets: keep all region types and later filter by user selection
  targets_all <- targets %>%
    dplyr::select(dplyr::all_of(c("Case ID", "Donor Status")), dplyr::all_of(aab_cols_targets), islet_key, type, class, area_um2, region_um2, area_density, count) %>%
    dplyr::mutate(
      type = dplyr::case_when(
        tolower(type) %in% c("islet_core", "core") ~ "islet_core",
        tolower(type) %in% c("islet_band", "band", "peri-islet", "peri_islet") ~ "islet_band",
        tolower(type) %in% c("islet_union", "union", "islet+20um", "islet_20um") ~ "islet_union",
        TRUE ~ tolower(type)
      )
      # Density is already in um2, no conversion needed
    ) %>%
    { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "targets_all:size_area") }

  # Synthesize missing union rows for counts by summing core + band (do NOT fabricate union area)
  if (nrow(targets_all) > 0) {
    keys <- c("Case ID", "Donor Status", "islet_key", "class")
    # Identify which key combos already have a union row
    have_union <- targets_all %>%
      dplyr::filter(type == "islet_union") %>%
      dplyr::select(dplyr::all_of(keys)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(.has_union = TRUE)

    core_rows <- targets_all %>% dplyr::filter(type == "islet_core") %>%
      dplyr::select(dplyr::all_of(c(keys, "count", "area_um2", "region_um2"))) %>%
      dplyr::rename(count_core = count, area_um2_core = area_um2, region_um2_core = region_um2)
    band_rows <- targets_all %>% dplyr::filter(type == "islet_band") %>%
      dplyr::select(dplyr::all_of(c(keys, "count", "area_um2", "region_um2"))) %>%
      dplyr::rename(count_band = count, area_um2_band = area_um2, region_um2_band = region_um2)
    union_missing <- core_rows %>%
      dplyr::inner_join(band_rows, by = keys) %>%
  { safe_left_join(., have_union, by = keys, context = "targets_union_missing:have_union") } %>%
      dplyr::filter(is.na(.has_union))
    if (nrow(union_missing) > 0) {
      synth <- union_missing %>%
        dplyr::transmute(`Case ID`, `Donor Status`, islet_key, type = "islet_union", class,
                          area_um2 = suppressWarnings(as.numeric(area_um2_core)) + suppressWarnings(as.numeric(area_um2_band)),
                          region_um2 = suppressWarnings(as.numeric(region_um2_core)) + suppressWarnings(as.numeric(region_um2_band)),
                          count = suppressWarnings(as.numeric(count_core)) + suppressWarnings(as.numeric(count_band))) %>%
        dplyr::mutate(area_density = dplyr::if_else(is.finite(region_um2) & region_um2 > 0,
                                                     area_um2 / region_um2, NA_real_))
      # Attach diameter via size_area
  synth <- synth %>% { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "targets_synth:size_area") }
      # Bind synthetic union rows to targets_all
      targets_all <- dplyr::bind_rows(targets_all, synth)
    }
  }
  # Ensure AAb flags are present for all rows, including synthetic ones
  if (!is.null(donors_meta)) {
    targets_all <- targets_all %>% dplyr::select(-dplyr::any_of(aab_cols_all)) %>%
      { safe_left_join(., donors_meta, by = c("Case ID","Donor Status"), context = "targets_all:donors_meta") }
  }

  # Markers with fraction positive / mean intensity (include LGALS3 sheet)
  markers <- master$markers %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key))
  markers_all <- markers %>%
    dplyr::select(dplyr::all_of(c("Case ID", "Donor Status")), dplyr::all_of(aab_cols_markers), islet_key, region_type, marker, n_cells, pos_count, pos_frac) %>%
    dplyr::mutate(
      region_type = dplyr::case_when(
        tolower(region_type) %in% c("islet_core", "core") ~ "islet_core",
        tolower(region_type) %in% c("islet_band", "band", "peri-islet", "peri_islet") ~ "islet_band",
        tolower(region_type) %in% c("islet_union", "union", "islet+20um", "islet_20um") ~ "islet_union",
        TRUE ~ tolower(region_type)
      )
    )
  # Add LGALS3 rows if available
  if (!is.null(master$lgals3) && nrow(master$lgals3) > 0) {
    g3 <- master$lgals3 %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key)) %>%
      mutate(marker = as.character(marker)) %>%
      dplyr::select(dplyr::all_of(c("Case ID", "Donor Status")), dplyr::all_of(intersect(colnames(master$lgals3), c("AAb_GADA","AAb_IA2A","AAb_ZnT8A","AAb_IAA","AAb_mIAA"))), islet_key, region_type, marker, n_cells, pos_count, pos_frac)
    if (!is.null(g3) && nrow(g3) > 0) {
      markers_all <- bind_rows(markers_all, g3)
    }
  }
  # Ensure AAb flags are present for all rows, including synthetic ones (markers)
  if (!is.null(donors_meta)) {
    markers_all <- markers_all %>% dplyr::select(-dplyr::any_of(aab_cols_all)) %>%
      { safe_left_join(., donors_meta, by = c("Case ID","Donor Status"), context = "markers_all:donors_meta") }
  }
  markers_all <- markers_all %>% { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "markers_all:size_area") }

  # Synthesize missing union rows for markers by summing core + band counts (do NOT fabricate area)
  if (nrow(markers_all) > 0) {
    mkeys <- c("Case ID", "Donor Status", "islet_key", "marker")
    have_union_m <- markers_all %>%
      dplyr::filter(region_type == "islet_union") %>%
      dplyr::select(dplyr::all_of(mkeys)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(.has_union = TRUE)

    core_m <- markers_all %>% dplyr::filter(region_type == "islet_core") %>% dplyr::select(dplyr::all_of(c(mkeys, "n_cells", "pos_count"))) %>% dplyr::rename(n_core = n_cells, pos_core = pos_count)
    band_m <- markers_all %>% dplyr::filter(region_type == "islet_band") %>%
      dplyr::select(dplyr::all_of(c(mkeys, "n_cells", "pos_count"))) %>%
      dplyr::rename(n_band = n_cells, pos_band = pos_count)
    union_m <- markers_all %>% dplyr::filter(region_type == "islet_union") %>% dplyr::select(dplyr::all_of(c(mkeys, "n_cells", "pos_count"))) %>% dplyr::rename(n_union = n_cells, pos_union = pos_count)
    union_missing_m <- core_m %>%
      dplyr::inner_join(band_m, by = mkeys) %>%
  { safe_left_join(., have_union_m, by = mkeys, context = "markers_union_missing:have_union_m") } %>%
      dplyr::filter(is.na(.has_union))
    if (nrow(union_missing_m) > 0) {
      synth_m <- union_missing_m %>%
        dplyr::transmute(`Case ID`, `Donor Status`, islet_key, region_type = "islet_union", marker,
                          n_cells = suppressWarnings(as.numeric(n_core)) + suppressWarnings(as.numeric(n_band)),
                          pos_count = suppressWarnings(as.numeric(pos_core)) + suppressWarnings(as.numeric(pos_band))) %>%
        dplyr::mutate(pos_frac = ifelse(is.finite(n_cells) & n_cells > 0,
                                        suppressWarnings(as.numeric(pos_count)) / suppressWarnings(as.numeric(n_cells)),
                                        NA_real_)) %>%
  { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "markers_synth_m:size_area") }
      markers_all <- dplyr::bind_rows(markers_all, synth_m)
    }

    # Backfill missing band rows when core and union exist but band is missing: band = union - core (counts only)
    have_band_m <- markers_all %>%
      dplyr::filter(region_type == "islet_band") %>%
      dplyr::select(dplyr::all_of(mkeys)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(.has_band = TRUE)
    band_missing_m <- core_m %>%
      dplyr::inner_join(union_m, by = mkeys) %>%
  { safe_left_join(., have_band_m, by = mkeys, context = "markers_band_missing:have_band_m") } %>%
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
  { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "markers_synth_band:size_area") }
      markers_all <- dplyr::bind_rows(markers_all, synth_band)
    }
  }

  # Fill pos_frac when n_cells and pos_count are present but pos_frac is NA; also compute pos_pct (0..100)
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

  # Composition by islet
  comp <- master$comp %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key))
  comp <- comp %>% dplyr::select(dplyr::all_of(c("Case ID", "Donor Status")), dplyr::all_of(aab_cols_comp), islet_key, cells_total, Ins_single, Glu_single, Stt_single,
                          Multi_Pos, Triple_Neg, Ins_any, Glu_any, Stt_any, dplyr::any_of("object_class"))
  # Endocrine object class (uniform root: Islet / EndocrineCluster; Singles are not rows since 2026-10-18). A
  # root without it holds
  # islets only. Carried onto the region tables so one app-wide selection filters every table.
  if (!"object_class" %in% names(comp)) comp$object_class <- "Islet"
  object_class_of <- dplyr::distinct(comp[, c("Case ID", "islet_key", "object_class")])
  for (tbl in c("targets_all", "markers_all")) {
    d <- get(tbl)
    if (nrow(d) > 0 && !"object_class" %in% names(d)) {
      d <- safe_left_join(d, object_class_of, by = c("Case ID", "islet_key"), context = paste0(tbl, ":object_class"))
      assign(tbl, d)
    }
  }
  comp <- comp %>% { safe_left_join(., region_areas, by = c("Case ID", "Donor Status", "islet_key"), context = "comp:region_areas") }

  # Merge phenotype proportions into comp (H5AD only; NULL from Excel path)
  if (!is.null(master$phenotypes) && nrow(master$phenotypes) > 0) {
    comp <- safe_left_join(comp, master$phenotypes,
                           by = c("Case ID", "islet_key"),
                           context = "comp:phenotypes")
    message("[prep-data] Merged ", ncol(master$phenotypes) - 2, " phenotype columns into comp")
  }

  # Merge neighborhood metrics into comp (H5AD only; NULL from Excel path)
  if (!is.null(master$neighborhood) && nrow(master$neighborhood) > 0) {
    comp <- safe_left_join(comp, master$neighborhood,
                           by = c("Case ID", "islet_key"),
                           context = "comp:neighborhood")
    message("[prep-data] Merged ", ncol(master$neighborhood) - 2, " neighborhood columns into comp")
  }

  # Carry the per-islet cell counts onto targets_all so the Microenvironment
  # mode can honor the "Min cells/islet" filter.
  #
  # `total_cells_core` / `total_cells_peri` arrive on `comp` via the neighborhood
  # merge above, but `targets_all` is built from the targets sheet with a fixed
  # select() that never had them. Both the Plot sidebar filter
  # (mod_plot_server.R::raw_df_base) and the Statistics tab filter
  # (mod_statistics_server.R::stats_run) are guarded on those column names, so in
  # Microenvironment mode the filter did not error — it silently did nothing,
  # while the exported CSV header still reported the min-cells value the user set.
  # Joining by (Case ID, islet_key) matches the phenotype / neighborhood joins
  # above; islet_key is unique within a donor, so this is many-to-one and cannot
  # grow targets_all.
  if (all(c("total_cells_core", "total_cells_peri") %in% colnames(comp))) {
    cell_counts <- comp %>%
      dplyr::select(`Case ID`, islet_key, total_cells_core, total_cells_peri) %>%
      dplyr::distinct(`Case ID`, islet_key, .keep_all = TRUE)
    n_before <- nrow(targets_all)
    targets_all <- targets_all %>%
      dplyr::select(-dplyr::any_of(c("total_cells_core", "total_cells_peri"))) %>%
      { safe_left_join(., cell_counts, by = c("Case ID", "islet_key"),
                       context = "targets_all:cell_counts") }
    if (nrow(targets_all) != n_before) {
      warning(sprintf("[prep-data] targets_all row count changed on cell-count join: %d -> %d",
                      n_before, nrow(targets_all)))
    }
    message("[prep-data] Joined per-islet cell counts onto targets_all (",
            sum(is.finite(suppressWarnings(as.numeric(targets_all$total_cells_core)))),
            " of ", nrow(targets_all), " rows matched)")
  }

  # Merge donor demographics into all dataframes (H5AD only; NULL from Excel path)
  if (!is.null(master$donor_demographics) && nrow(master$donor_demographics) > 0) {
    targets_all <- safe_left_join(targets_all, master$donor_demographics, by = "Case ID", context = "targets_all:demographics")
    markers_all <- safe_left_join(markers_all, master$donor_demographics, by = "Case ID", context = "markers_all:demographics")
    comp <- safe_left_join(comp, master$donor_demographics, by = "Case ID", context = "comp:demographics")
    message("[prep-data] Merged donor demographics (age, gender)")
  }

  # ---- nPOD clinical metadata (disease duration, BMI, HbA1c, C-peptide, AAb) ----
  # Built by scripts/build_donor_metadata.py from the nPOD "All-Download" export,
  # which is deliberately NOT in git (741 donors of free-text clinical narrative;
  # `*.xlsx` is gitignored). Only the derived 15-row CSV is committed.
  #
  # `age` comes from the H5AD and this join adds `age_years_npod` alongside it. On the frozen 15-donor data
  # donor 112's H5AD age (3.00) disagrees with nPOD (3.72); from the uniform v10 root (owner, 2026-10-22) the
  # H5AD age is nPOD's, so the two agree. Every other column of the file joins as is (v10 adds height, weight,
  # BMI percentile, glucose, ICU stay, transit time, pancreas weight, RIN, race/ethnicity, recovery type,
  # cause of death, ABO group).
  #
  # Absent file is a normal state (Excel fallback path, fresh checkout before the
  # script is run) — skip quietly, every consumer already guards on the columns.
  donor_meta_path <- data_path("donor_metadata.csv")
  if (file.exists(donor_meta_path)) {
    donor_meta <- tryCatch({
      dm <- utils::read.csv(donor_meta_path, stringsAsFactors = FALSE,
                            check.names = FALSE)
      dm$`Case ID` <- suppressWarnings(as.integer(dm$case_id))
      # `sex`/`donor_type` duplicate columns the app already carries under other
      # names (gender / Donor Status).
      dm[, setdiff(colnames(dm), c("case_id", "sex", "donor_type")), drop = FALSE]
    }, error = function(e) {
      message("event=donor_metadata_load_failed"); NULL
    })
    if (!is.null(donor_meta) && nrow(donor_meta) > 0) {
      # Drop any column the target frame ALREADY has, keeping the existing one.
      #
      # dplyr suffixes name collisions .x/.y rather than erroring, and every
      # downstream lookup in this app is by exact name — so a collision does not
      # fail loudly, it silently disables a feature. The AAb flags are the live
      # example: the H5AD already supplies AAb_GADA/IA2A/ZnT8A/mIAA, the nPOD
      # sheet supplies the same four, and joining both renamed all of them to
      # .x/.y. That turned the Plot sidebar's AAb filter into a no-op, because
      # it resolves columns via intersect(flags, colnames(df)) and matched none.
      # `n_aab` survives here because it is genuinely new.
      drop_dupes <- function(df, target) {
        dupes <- setdiff(intersect(colnames(df), colnames(target)), "Case ID")
        if (length(dupes)) {
          message("[prep-data] donor_metadata: keeping existing ",
                  paste(dupes, collapse = ", "), " (not overwriting)")
        }
        df[, setdiff(colnames(df), dupes), drop = FALSE]
      }
      n_before <- c(nrow(targets_all), nrow(markers_all), nrow(comp))
      targets_all <- safe_left_join(targets_all, drop_dupes(donor_meta, targets_all),
                                    by = "Case ID", context = "targets_all:donor_metadata")
      markers_all <- safe_left_join(markers_all, drop_dupes(donor_meta, markers_all),
                                    by = "Case ID", context = "markers_all:donor_metadata")
      comp <- safe_left_join(comp, drop_dupes(donor_meta, comp),
                             by = "Case ID", context = "comp:donor_metadata")
      n_after <- c(nrow(targets_all), nrow(markers_all), nrow(comp))
      if (!identical(n_before, n_after)) {
        warning(sprintf("[prep-data] donor_metadata join changed row counts: %s -> %s",
                        paste(n_before, collapse = "/"), paste(n_after, collapse = "/")))
      }
      message("[prep-data] Merged nPOD donor metadata (", nrow(donor_meta),
              " donors, ", ncol(donor_meta) - 1, " columns)")
    }
  }

  # Group aggregates ("All macrophages", owner 2026-10-27) as derived proportion columns, before the pair table.
  comp <- add_phenotype_aggregate_columns(comp)

  # Combined (core+peri) phenotype pair table. peri_prop_* columns sanitize
  # the phenotype name (spaces -> _, + -> plus) but prop_* columns keep the
  # original label, so a direct stem intersection misses ~half the pairs. Build
  # a lookup that resolves each peri-style safe-name back to its prop_* core
  # column by trying un-sanitization candidates (mirrors the Spatial tab's
  # phenotype_options() resolver). Consumed by mod_plot_server.R Cell
  # Populations selector and raw_df_base() Combined branch.
  combined_pairs <- build_combined_pairs(comp)
  if (!is.null(combined_pairs) && nrow(combined_pairs) > 0) {
    message("[prep-data] Built combined core+peri lookup for ", nrow(combined_pairs), " phenotypes")
  }

  message("[prep-final] size_area=", paste(class(size_area), collapse='/'),
    " targets_all=", paste(class(targets_all), collapse='/'),
    " markers_all=", paste(class(markers_all), collapse='/'),
    " comp=", paste(class(comp), collapse='/'))
  list(core_area = size_area, targets_all = targets_all, markers_all = markers_all,
       comp = comp, combined_pairs = combined_pairs,
       protein_embedding = master$protein_embedding,
       composition_embedding = master$composition_embedding)
}

# Simple NA audit
audit_na <- function(df, label) {
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  na_cnt <- vapply(df, function(x) sum(is.na(x)), integer(1))
  total <- nrow(df)
  pct <- if (total > 0) round(100 * na_cnt / total, 2) else rep(0, length(na_cnt))
  msg <- paste0("[NA audit] ", label, ": ", paste(names(na_cnt), paste0(na_cnt, " (", pct, "%)"), sep = "=", collapse = "; "))
  message(msg)
}

bin_islet_sizes <- function(df, diam_col, width) {
  x <- suppressWarnings(as.numeric(df[[diam_col]]))
  max_x <- max(x, na.rm = TRUE)
  max_x <- ifelse(is.finite(max_x), max_x, 0)
  # Define numeric bins: [lo, hi)
  bin_lo <- floor(x / width) * width
  bin_hi <- bin_lo + width
  diam_mid <- bin_lo + width/2
  # Build a human-readable label
  diam_bin <- paste0("[", bin_lo, ", ", bin_hi, ")")
  df$diam_bin <- factor(diam_bin, levels = unique(diam_bin[order(bin_lo)]), ordered = TRUE)
  df$diam_mid <- diam_mid
  df
}

# ---------- H5AD-based data loading ----------
# Reconstructs the same data frames as load_master() + prep_data() from the
# enriched H5AD built by scripts/build_h5ad_for_app.py

#' Load data from enriched H5AD file (islet_explorer.h5ad)
#' Returns the same list structure as load_master(): list(markers, targets, comp, lgals3)
#' @param path Path to the enriched H5AD file
#' @return list with targets, markers, comp, lgals3 DataFrames, or NULL on failure
load_master_h5ad <- function(path) {
  if (!requireNamespace("anndata", quietly = TRUE)) {
    message("[H5AD] anndata package not available, falling back to Excel")
    return(NULL)
  }

  tryCatch({
    message("[H5AD] Loading enriched data")
    ad <- islet_read_h5ad_shared(path)

    # Cache ad$uns as R list once (0.2s) instead of per-key bridge crossings (~14s for 62 keys)
    uns <- ad$uns

    # Reconstruct targets/markers/comp from cached uns list
    targets <- reconstruct_groovy_df_from_list(uns, "targets")
    markers <- reconstruct_groovy_df_from_list(uns, "markers")
    comp    <- reconstruct_groovy_df_from_list(uns, "composition")

    # Split LGALS3 rows from markers (they share the same groovy storage)
    lgals3 <- NULL
    if (!is.null(markers) && "marker" %in% names(markers)) {
      lgals3_mask <- markers$marker == "LGALS3"
      if (any(lgals3_mask)) {
        lgals3 <- markers[lgals3_mask, , drop = FALSE]
        markers <- markers[!lgals3_mask, , drop = FALSE]
      }
    }

    # Extract phenotype proportions from .obs
    phenotype_df <- tryCatch({
      obs <- as.data.frame(ad$obs)
      prop_cols <- grep("^prop_", colnames(obs), value = TRUE)
      if (length(prop_cols) > 0 && "imageid" %in% colnames(obs) && "base_islet_id" %in% colnames(obs)) {
        phen <- obs[, c("imageid", "base_islet_id", prop_cols), drop = FALSE]
        phen$`Case ID` <- as.integer(as.character(phen$imageid))
        phen$islet_key <- gsub("^Islet_Islet_", "Islet_", as.character(phen$base_islet_id))
        phen[, c("Case ID", "islet_key", prop_cols), drop = FALSE]
      } else NULL
    }, error = function(e) { message("event=h5ad_phenotype_extract_failed"); NULL })

    # Extract donor demographics (one row per donor)
    donor_demographics <- tryCatch({
      obs <- as.data.frame(ad$obs)
      if (all(c("imageid", "age", "gender") %in% colnames(obs))) {
        demo <- data.frame(
          `Case ID` = as.integer(as.character(obs$imageid)),
          age = as.numeric(obs$age),
          gender = as.character(obs$gender),
          check.names = FALSE, stringsAsFactors = FALSE
        )
        demo[!duplicated(demo$`Case ID`), , drop = FALSE]
      } else NULL
    }, error = function(e) { message("event=h5ad_demographics_extract_failed"); NULL })

    # Extract neighborhood metrics from .obs (added by compute_neighborhood_metrics.py → build_h5ad_for_app.py)
    neighborhood_df <- tryCatch({
      obs <- as.data.frame(ad$obs)
      nbr_cols <- grep("^peri_prop_|^peri_count_|^immune_|^cd8_|^tcell_|^enrich_z_|^min_dist_|^total_cells_peri|^total_cells_core|^immune_count_|^dpt_pseudotime|^leiden_",
                       colnames(obs), value = TRUE)
      if (length(nbr_cols) > 0 && "imageid" %in% colnames(obs) && "base_islet_id" %in% colnames(obs)) {
        nbr <- obs[, c("imageid", "base_islet_id", nbr_cols), drop = FALSE]
        nbr$`Case ID` <- as.integer(as.character(nbr$imageid))
        nbr$islet_key <- gsub("^Islet_Islet_", "Islet_", as.character(nbr$base_islet_id))
        nbr[, c("Case ID", "islet_key", nbr_cols), drop = FALSE]
      } else NULL
    }, error = function(e) { message("event=h5ad_neighborhood_extract_failed"); NULL })

    message("[H5AD] Loaded: targets=", if (!is.null(targets)) nrow(targets) else 0,
            " markers=", if (!is.null(markers)) nrow(markers) else 0,
            " comp=", if (!is.null(comp)) nrow(comp) else 0,
            " lgals3=", if (!is.null(lgals3)) nrow(lgals3) else 0,
            " phenotypes=", if (!is.null(phenotype_df)) ncol(phenotype_df) - 2 else 0,
            " demographics=", if (!is.null(donor_demographics)) nrow(donor_demographics) else 0,
            " neighborhood=", if (!is.null(neighborhood_df)) ncol(neighborhood_df) - 2 else 0)

    list(markers = markers, targets = targets, comp = comp, lgals3 = lgals3,
         phenotypes = phenotype_df, donor_demographics = donor_demographics,
         neighborhood = neighborhood_df, protein_embedding = uns[["protein_embedding"]],
         composition_embedding = uns[["composition_embedding"]])
  }, error = function(e) {
    message("event=h5ad_load_failed")
    NULL
  })
}

#' Reconstruct a groovy DataFrame from a pre-fetched .uns R list
#' @param uns_list R list from ad$uns (single bridge crossing)
#' @param sheet One of "targets", "markers", "composition"
#' @return data.frame or NULL
reconstruct_groovy_df_from_list <- function(uns_list, sheet) {
  prefix <- paste0("groovy_", sheet, "_")
  cols_key <- paste0("groovy_", sheet, "_columns")
  nrows_key <- paste0("groovy_", sheet, "_n_rows")

  if (is.null(uns_list[[cols_key]])) return(NULL)

  col_names <- uns_list[[cols_key]]
  n_rows <- as.integer(uns_list[[nrows_key]])

  if (n_rows == 0) return(NULL)

  df <- data.frame(row.names = seq_len(n_rows))
  for (col in col_names) {
    key <- paste0(prefix, col)
    vals <- uns_list[[key]]
    if (!is.null(vals)) {
      df[[col]] <- vals
    }
  }

  # Standardize Case ID to match Excel format (integer-like)
  if ("Case ID" %in% names(df)) {
    df[["Case ID"]] <- suppressWarnings(as.integer(df[["Case ID"]]))
  }

  df
}

#' Legacy: Reconstruct a groovy DataFrame via ad$uns (slow, ~14s for 62 keys)
#' Kept for backward compatibility but unused in normal path
reconstruct_groovy_df <- function(ad, sheet) {
  reconstruct_groovy_df_from_list(ad$uns, sheet)
}

#' Try loading from H5AD first, fall back to Excel
#' @param h5ad_path Path to enriched H5AD (or NULL to skip)
#' @param excel_path Path to master_results.xlsx
#' @return list(markers, targets, comp, lgals3) — same structure as load_master()
load_master_auto <- function(h5ad_path = NULL, excel_path = master_path) {
  # Try H5AD first
  if (!is.null(h5ad_path) && file.exists(h5ad_path)) {
    result <- load_master_h5ad(h5ad_path)
    if (!is.null(result)) {
      message("[DATA] Using H5AD source")
      return(result)
    }
  }

  # Fall back to Excel
  if (file.exists(excel_path)) {
    message("[DATA] Using Excel source")
    return(load_master(excel_path))
  }

  stop("No data source found. Checked:\n  H5AD: ", h5ad_path, "\n  Excel: ", excel_path)
}

#' Load the precomputed, immutable publication summaries.
#'
#' These files contain the output of prep_data(), so no H5AD bridge or joins run
#' per reviewer session. The manifest is validated before any table is exposed.
load_prepared_release <- function(root = data_root) {
  manifest <- islet_validate_release(root)
  out <- lapply(ISLET_SUMMARY_ARTIFACTS, function(rel) {
    islet_read_parquet(file.path(root, rel))
  })

  artifacts <- manifest$artifacts %||% list()
  for (name in names(out)) {
    expected <- artifacts[[name]]$row_count %||% NULL
    if (!is.null(expected) && !identical(as.numeric(nrow(out[[name]])), as.numeric(expected))) {
      stop("Release row-count mismatch for ", name, ": expected ", expected,
           ", read ", nrow(out[[name]]), call. = FALSE)
    }
  }
  # A schema-5 release carries the H5AD's embedding records (the Spatial tab's method text reads them, as
  # prep_data() supplies them in development); a schema-1 release has none.
  embedding <- islet_read_embedding_metadata(manifest, root)
  if (!is.null(embedding)) {
    out$composition_embedding <- embedding$composition_embedding
    out$protein_embedding <- embedding$protein_embedding
  }
  islet_log_event("release_data_loaded")
  out
}

#' Process-wide prepared data object.
#'
#' Production is fail-closed: it must have a complete manifest-backed release.
#' Development reads H5AD (or Excel) directly from its selected data root,
#' regardless of whether a manifest happens to be present there.
load_shared_prepared_data <- function(root = data_root) {
  cache_key <- paste0("prepared::", normalizePath(root, mustWork = FALSE))
  if (!is.null(.islet_backend[[cache_key]])) return(.islet_backend[[cache_key]])

  prepared <- if (islet_is_production()) {
    load_prepared_release(root)
  } else {
    message("[DATA] Using development H5AD/Excel source under ", root)
    prep_data(load_master_auto(h5ad_path = file.path(root, "islet_explorer.h5ad"),
                               excel_path = file.path(root, "master_results.xlsx")))
  }

  # Fail closed rather than filter: a root that still carries an excluded donor
  # was built before the removal and must not be served.
  served <- unlist(lapply(intersect(c("comp", "targets_all", "markers_all"), names(prepared)),
                          function(name) prepared[[name]][["Case ID"]]))
  excluded <- islet_excluded_donors_in(served)
  if (length(excluded)) {
    stop("Data root ", root, " contains excluded donor(s) ", paste(excluded, collapse = ", "),
         "; select a root built without them (see scripts/excluded_donors.csv).", call. = FALSE)
  }

  if (!islet_is_production()) {
    for (name in intersect(c("targets_all", "markers_all", "comp"), names(prepared))) {
      try(audit_na(prepared[[name]], name), silent = TRUE)
    }
  }
  .islet_backend[[cache_key]] <- prepared
  prepared
}


#' Keep only the selected endocrine-object classes in every per-object table of a prepared list.
#' `classes` NULL keeps everything; an empty selection keeps nothing.
filter_prepared_by_object_class <- function(prepared, classes) {
  if (is.null(prepared) || is.null(classes)) return(prepared)
  for (tbl in c("comp", "targets_all", "markers_all")) {
    d <- prepared[[tbl]]
    if (!is.null(d) && "object_class" %in% names(d)) {
      prepared[[tbl]] <- d[d$object_class %in% classes, , drop = FALSE]
    }
  }
  prepared
}
