library(shiny)
library(shinyjs)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(plotly)
library(broom)
library(jsonlite)

# Load anndata package for trajectory analysis (install with BiocManager::install('anndata'))
suppressPackageStartupMessages({
  if (!requireNamespace("anndata", quietly = TRUE)) {
    message("Package 'anndata' not found. Install with: BiocManager::install('anndata')")
  }
})

# ---- Runtime diagnostics (temporary) ----
APP_VERSION <- "traj-baseR-fix-v7"
try({
  message("[APP LOAD] Version=", APP_VERSION, " time=", Sys.time())
  message("[APP LOAD] getwd=", tryCatch(getwd(), error=function(e) NA))
  message("[APP LOAD] list.files(.) contains app.R? ", any(grepl('^app.R$', list.files('.'))))
  message("[APP LOAD] file.info(app.R)$mtime=", tryCatch(file.info('app.R')$mtime, error=function(e) NA))
  # Create a flag file to confirm this specific app.R executed.
  flag_path <- file.path(dirname(sys.frame(1)$ofile %||% 'app.R'), 'APP_LOADED_FLAG')
  writeLines(paste0('loaded ', Sys.time(), ' version=', APP_VERSION), flag_path)
  # Disable tracing that was causing issues with trajectory data
  # if (!isTRUE(getOption('traj_left_join_traced'))) {
  #   try(trace(dplyr::left_join,
  #             quote({
  #               message(sprintf('[TRACE left_join] classes lhs=%s rhs=%s by=%s',
  #                                paste(class(x), collapse='/'),
  #                                paste(class(y), collapse='/'),
  #                                paste(names(by %||% list()), collapse=',')))
  #             }), print = FALSE), silent = TRUE)
  #   options(traj_left_join_traced = TRUE)
  # }
}, silent = TRUE)


# Diagnostic-safe wrapper around dplyr::left_join to catch non-data.frame inputs early.
safe_left_join <- function(x, y, by, context) {
  cx <- paste(class(x), collapse = "/"); cy <- paste(class(y), collapse = "/")
  okx <- inherits(x, "data.frame"); oky <- inherits(y, "data.frame")
  if (!okx || !oky) {
    message(sprintf("[safe_left_join] %s: coercing inputs (lhs=%s rhs=%s)", context, cx, cy))
    x <- tryCatch(as.data.frame(x), error = function(e) { message("[safe_left_join] lhs coercion failed: ", e$message); x })
    y <- tryCatch(as.data.frame(y), error = function(e) { message("[safe_left_join] rhs coercion failed: ", e$message); y })
  }
  if (!inherits(x, "data.frame") || !inherits(y, "data.frame")) {
    stop(sprintf("[safe_left_join] %s: cannot join; lhs class=%s rhs class=%s", context, paste(class(x), collapse='/'), paste(class(y), collapse='/')))
  }
  tryCatch(
    dplyr::left_join(x, y, by = by),
    error = function(e) {
      message(sprintf("[safe_left_join] %s: left_join error: %s", context, e$message))
      stop(e)
    }
  )
}
## for base64 encoding channel_config payload to Avivator
## (installed by scripts/install_shiny_deps.R)
suppressPackageStartupMessages({
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop("Package 'base64enc' is required. Run scripts/install_shiny_deps.R to install dependencies.")
  }
})


master_path <- file.path("..", "..", "data", "master_results.xlsx")

project_root <- tryCatch(normalizePath(file.path("..", ".."), mustWork = FALSE), error = function(e) NULL)

# Restore basic viewer components needed for Avivator
# Support either 'Channel_names' or 'Channel_names.txt'
channel_names_path <- NULL
for (cand in c(file.path("Channel_names"), file.path("Channel_names.txt"))) {
  if (file.exists(cand)) { channel_names_path <- cand; break }
}
if (!is.null(channel_names_path)) {
  channel_names_vec <- tryCatch({
    lines <- readLines(channel_names_path, warn = FALSE)
    rx <- "^\\s*(.*?)\\s*\\(C(\\d+)\\)\\s*$"
    matched <- stringr::str_match(lines, rx)
    matched <- matched[!is.na(matched[, 1]), , drop = FALSE]
    if (nrow(matched) == 0) return(NULL)
    idx <- as.integer(matched[, 3])
    nm <- matched[, 2]
    order_df <- dplyr::arrange(data.frame(idx = idx, name = nm, stringsAsFactors = FALSE), idx)
    max_idx <- max(order_df$idx, na.rm = TRUE)
    names_vec <- rep(NA_character_, max_idx)
    names_vec[order_df$idx] <- order_df$name
    names_vec
  }, error = function(e) NULL)
} else {
  channel_names_vec <- NULL
}

# Local images setup
# Prefer static www/local_images so shiny-server can serve large OME-TIFFs with HTTP Range support.
# Only add a dynamic resource path if www/local_images does not exist.
local_images_env <- Sys.getenv("LOCAL_IMAGE_ROOT", unset = "")
local_images_root <- NULL
if (nzchar(local_images_env)) {
  local_images_root <- tryCatch(normalizePath(local_images_env, mustWork = TRUE), error = function(e) NULL)
}
if (is.null(local_images_root)) {
  candidate <- if (!is.null(project_root)) file.path(project_root, "local_images") else NULL
  if (!is.null(candidate) && dir.exists(candidate)) {
    local_images_root <- candidate
  } else {
    local_images_root <- tryCatch(normalizePath(file.path("..","..","local_images"), mustWork = FALSE), error = function(e) NULL)
  }
}

www_local_images_dir <- file.path("www", "local_images")
has_www_local_images <- dir.exists(www_local_images_dir)

# Only register /local_images if there is no static folder under www
if (!has_www_local_images && !is.null(local_images_root)) {
  try({ shiny::addResourcePath("local_images", local_images_root) }, silent = TRUE)
}

# ---- Viewer helpers/defaults (Avivator) ----
# Optional default image URL; if NULL, we'll auto-pick the first available under /local_images
default_image_url <- NULL

resolve_avivator_base <- function() {
  # Return the URL path to the local avivator build under www/, or NULL if missing
  local_index <- file.path("www", "avivator", "index.html")
  if (file.exists(local_index)) return("avivator/index.html")
  NULL
}

# Build a base64-encoded channel configuration understood by embedded Avivator.
# Expected shape (before base64):
#   {
#     channelNames: ["DAPI", "CD31", ...],
#     primaryChannels: [
#       { index: 0, name: "DAPI", color: "#7F7F7F", visible: true },
#       { index: 25, name: "INS",  color: "#E41A1C", visible: true },
#       { index: 19, name: "GCG",  color: "#377EB8", visible: true },
#       { index: 13, name: "SST",  color: "#FFCC00", visible: true }
#     ]
#   }
build_channel_config_b64 <- function(names_vec) {
  if (is.null(names_vec) || length(names_vec) == 0) return(NULL)

  # Helper: RGB ints (0..255) to hex string
  rgb_hex <- function(r, g, b) sprintf("#%02X%02X%02X", as.integer(r), as.integer(g), as.integer(b))
  # Palette for key channels
  hex_red    <- rgb_hex(228,  26,  28)  # INS
  hex_blue   <- rgb_hex( 55, 126, 184)  # GCG
  hex_yellow <- rgb_hex(255, 204,   0)  # SST
  hex_grey   <- rgb_hex(127, 127, 127)  # DAPI

  # Locate indices for the key channels (names_vec is 1-based by channel number)
  find_idx0 <- function(tag) {
    hit <- which(toupper(names_vec) == tag)
    if (length(hit) == 0) return(NA_integer_)
    (hit[[1]] - 1L) # zero-based for viewer
  }
  idx0_dapi <- find_idx0("DAPI")
  idx0_ins  <- find_idx0("INS")
  idx0_gcg  <- find_idx0("GCG")
  idx0_sst  <- find_idx0("SST")

  prim <- list()
  add_pc <- function(idx0, name, hex) {
    if (!is.na(idx0) && idx0 >= 0) prim[[length(prim)+1]] <<- list(index = idx0, name = name, color = hex, visible = TRUE)
  }
  # Order: INS (red), GCG (blue), SST (yellow), DAPI (grey)
  add_pc(idx0_ins,  "INS",  hex_red)
  add_pc(idx0_gcg,  "GCG",  hex_blue)
  add_pc(idx0_sst,  "SST",  hex_yellow)
  add_pc(idx0_dapi, "DAPI", hex_grey)

  # Build config
  cfg <- list(
    channelNames = as.list(as.character(names_vec)),
    primaryChannels = prim
  )
  js <- jsonlite::toJSON(cfg, auto_unbox = TRUE)
  base64enc::base64encode(charToRaw(js))
}

# ---------- Data loading and wrangling ----------

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

add_islet_key <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)

  # Start with an existing islet_key column if present; otherwise create NA placeholder
  if (!"islet_key" %in% names(df)) df$islet_key <- NA_character_

  # 1) Derive from 'region' if available, e.g., "Islet_200_union" -> "Islet_200"
  if ("region" %in% names(df)) {
    key_from_region <- stringr::str_extract(df$region, "Islet_\\d+")
    df$islet_key <- dplyr::coalesce(df$islet_key, key_from_region)
  }

  # 2) Derive from 'name' if available, e.g., "Islet_200" or any text containing digits
  if ("name" %in% names(df)) {
    key_from_name <- stringr::str_extract(df$name, "Islet_\\d+")
    only_digits  <- stringr::str_extract(df$name, "\\d+")
    fallback_name <- ifelse(!is.na(only_digits), paste0("Islet_", only_digits), NA_character_)
    df$islet_key <- dplyr::coalesce(df$islet_key, key_from_name, fallback_name)
  }

  # 3) Fallback from numeric/string 'islet_id' if present
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
  # Determine which autoantibody columns are available
  aab_cols_targets <- intersect(c("AAb_GADA","AAb_IA2A","AAb_ZnT8A","AAb_IAA","AAb_mIAA"), colnames(master$targets))
  aab_cols_markers <- intersect(c("AAb_GADA","AAb_IA2A","AAb_ZnT8A","AAb_IAA","AAb_mIAA"), colnames(master$markers))
  aab_cols_comp    <- intersect(c("AAb_GADA","AAb_IA2A","AAb_ZnT8A","AAb_IAA","AAb_mIAA"), colnames(master$comp))
  aab_cols_all     <- unique(c(aab_cols_targets, aab_cols_markers, aab_cols_comp))
  # Donor-level metadata table to backfill AAb flags on synthetic rows
  donors_meta <- NULL
  if (!is.null(master$comp) && nrow(master$comp) > 0) {
    donors_meta <- master$comp %>% dplyr::select(dplyr::all_of(c("Case ID","Donor Status", aab_cols_comp))) %>% dplyr::distinct()
  } else if (!is.null(master$targets) && nrow(master$targets) > 0) {
    donors_meta <- master$targets %>% dplyr::select(dplyr::all_of(c("Case ID","Donor Status", aab_cols_targets))) %>% dplyr::distinct()
  } else if (!is.null(master$markers) && nrow(master$markers) > 0) {
    donors_meta <- master$markers %>% dplyr::select(dplyr::all_of(c("Case ID","Donor Status", aab_cols_markers))) %>% dplyr::distinct()
  }
  # Islet size proxy per islet corefor diameter
  targets <- master$targets %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key))
  core_area <- targets %>%
    dplyr::filter(tolower(type) == "islet_core") %>%
    dplyr::select(`Case ID`, `Donor Status`, islet_key, core_region_um2 = region_um2) %>%
    dplyr::distinct()
  # Diameter is computed ONLY from core area; if no core exists, diameter is NA
  size_area <- core_area %>%
    dplyr::mutate(islet_diam_um = compute_diameter_um(core_region_um2)) %>%
    dplyr::select(`Case ID`, `Donor Status`, islet_key, islet_diam_um)

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

    core_rows <- targets_all %>% dplyr::filter(type == "islet_core") %>% dplyr::select(dplyr::all_of(c(keys, "count"))) %>% dplyr::rename(count_core = count)
    band_rows <- targets_all %>% dplyr::filter(type == "islet_band") %>% dplyr::select(dplyr::all_of(c(keys, "count"))) %>% dplyr::rename(count_band = count)
    union_missing <- core_rows %>%
      dplyr::inner_join(band_rows, by = keys) %>%
  { safe_left_join(., have_union, by = keys, context = "targets_union_missing:have_union") } %>%
      dplyr::filter(is.na(.has_union))
    if (nrow(union_missing) > 0) {
      synth <- union_missing %>%
        dplyr::transmute(`Case ID`, `Donor Status`, islet_key, type = "islet_union", class,
                          area_um2 = NA_real_, region_um2 = NA_real_,
                          area_density = NA_real_,
                          count = suppressWarnings(as.numeric(count_core)) + suppressWarnings(as.numeric(count_band)))
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
    band_m <- markers_all %>% dplyr::filter(region_type == "islet_band") %>% dplyr::select(dplyr::all_of(c(mkeys, "n_cells", "pos_count"))) %>% dplyr::rename(n_band = n_cells, pos_band = pos_count)
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
                          Multi_Pos, Triple_Neg, Ins_any, Glu_any, Stt_any)
  comp <- comp %>% { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "comp:size_area") }

  message("[prep-final] size_area=", paste(class(size_area), collapse='/'),
    " targets_all=", paste(class(targets_all), collapse='/'),
    " markers_all=", paste(class(markers_all), collapse='/'),
    " comp=", paste(class(comp), collapse='/'))
  list(core_area = size_area, targets_all = targets_all, markers_all = markers_all, comp = comp)
}

# Simple NA audit
audit_na <- function(df, label) {
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  na_cnt <- vapply(df, function(x) sum(is.na(x)), integer(1))
  total <- nrow(df)
  pct <- ifelse(total > 0, round(100 * na_cnt / total, 2), 0)
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

# Unbiased pseudotime from features (no group labels). Tries DiffusionMap -> principal curve -> PC1.
## pseudotime removed

summary_stats <- function(df, group_cols, value_col, stat = c("mean_se","mean_sd","median_iqr")) {
  stat <- match.arg(stat)
  df %>% group_by(across(all_of(group_cols))) %>%
    summarise(
      n = sum(!is.na(.data[[value_col]])),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd = sd(.data[[value_col]], na.rm = TRUE),
      se = ifelse(n > 0, sd / sqrt(pmax(n, 1)), NA_real_),
      median = median(.data[[value_col]], na.rm = TRUE),
      q1 = quantile(.data[[value_col]], 0.25, na.rm = TRUE, type = 7),
      q3 = quantile(.data[[value_col]], 0.75, na.rm = TRUE, type = 7),
      .groups = "drop"
    ) %>%
    mutate(
      y = dplyr::case_when(
        stat == "mean_se" ~ mean,
        stat == "mean_sd" ~ mean,
        stat == "median_iqr" ~ median,
        TRUE ~ mean
      ),
      ymin = dplyr::case_when(
        stat == "mean_se" ~ mean - se,
        stat == "mean_sd" ~ mean - sd,
        stat == "median_iqr" ~ q1,
        TRUE ~ mean - se
      ),
      ymax = dplyr::case_when(
        stat == "mean_se" ~ mean + se,
        stat == "mean_sd" ~ mean + sd,
        stat == "median_iqr" ~ q3,
        TRUE ~ mean + se
      )
    )
}

per_bin_anova <- function(df, bin_col, group_col, value_col, mid_col = "diam_mid") {
  if (nrow(df) == 0) return(tibble(bin = character(), mid = numeric(), p_anova = numeric()))
  # capture unique bins with their numeric mid
  bmeta <- df %>% filter(!is.na(.data[[bin_col]])) %>%
    group_by(.data[[bin_col]]) %>%
    summarise(mid = suppressWarnings(as.numeric(first(na.omit(.data[[mid_col]])))), .groups = "drop")
  res <- lapply(seq_len(nrow(bmeta)), function(i) {
    b <- bmeta[[bin_col]][i]
    sub <- df %>% filter(.data[[bin_col]] == b)
    # Require at least 2 donor groups present and at least some finite values
    if (n_distinct(na.omit(sub[[group_col]])) < 2 || sum(is.finite(sub[[value_col]])) < 2) return(NULL)
    fit <- tryCatch(aov(reformulate(group_col, response = value_col), data = sub), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    # Extract one-way ANOVA p-value from ANOVA table (first term)
    pval <- tryCatch({
      at <- anova(fit)
      as.numeric(at[["Pr(>F)"]][1])
    }, error = function(e) NA_real_)
    tibble(
      bin = as.character(b),
      mid = bmeta$mid[i],
      p_anova = pval
    )
  })
  out <- bind_rows(res)
  if (nrow(out) == 0) return(tibble(bin = character(), mid = numeric(), p_anova = numeric()))
  arrange(out, mid)
}

per_bin_kendall <- function(df, bin_col, group_col, value_col, mid_col = "diam_mid") {
  if (nrow(df) == 0) return(tibble(bin = character(), mid = numeric(), p_kendall = numeric(), tau = numeric()))
  # Map group to ordered numeric along pseudotime ND < Aab+ < T1D
  code_group <- function(x) {
    x <- as.character(x)
    ifelse(x == "ND", 0, ifelse(x == "Aab+", 1, ifelse(x == "T1D", 2, NA_real_)))
  }
  bmeta <- df %>% filter(!is.na(.data[[bin_col]])) %>% group_by(.data[[bin_col]]) %>% summarise(mid = suppressWarnings(as.numeric(first(na.omit(.data[[mid_col]])))), .groups = "drop")
  res <- lapply(seq_len(nrow(bmeta)), function(i) {
    b <- bmeta[[bin_col]][i]
    sub <- df %>% filter(.data[[bin_col]] == b)
    x <- code_group(sub[[group_col]])
    y <- suppressWarnings(as.numeric(sub[[value_col]]))
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]; y <- y[keep]
    if (length(unique(x)) < 2 || length(y) < 3) return(NULL)
    ct <- tryCatch(cor.test(x, y, method = "kendall", exact = FALSE), error = function(e) NULL)
    if (is.null(ct)) return(NULL)
    tibble(bin = as.character(b), mid = bmeta$mid[i], p_kendall = unname(ct$p.value), tau = unname(ct$estimate))
  })
  out <- bind_rows(res)
  if (nrow(out) == 0) return(out)
  arrange(out, mid)
}

# ---------- UI ----------

ui <- fluidPage(
  useShinyjs(),
  tags$head(tags$style(HTML("\n    body.viewer-mode div.col-sm-4,\n    body.viewer-mode div.col-sm-3,\n    body.viewer-mode div.col-lg-3 {\n      display: none !important;\n    }\n    body.viewer-mode div.col-sm-8,\n    body.viewer-mode div.col-lg-9 {\n      width: 100% !important;\n      max-width: 100% !important;\n      flex: 0 0 100%;\n    }\n    body.viewer-mode .tab-content {\n      padding-left: 0 !important;\n      padding-right: 0 !important;\n    }\n    body.trajectory-mode .container-fluid > .row > .col-sm-3 {\n      display: none !important;\n    }\n    body.trajectory-mode .container-fluid > .row > .col-sm-9 {\n      width: 100% !important;\n      max-width: 100% !important;\n      flex: 0 0 100%;\n    }\n  "))),
  titlePanel("Islet Area Distributions"),
  uiOutput("theme_css"),
  sidebarLayout(
    sidebarPanel(
      h4("Data & Filters"),
      width = 3,
      selectInput("mode", "Focus", choices = c("Markers", "Targets", "Composition"), selected = "Composition"),
      uiOutput("region_selector"),
      uiOutput("dynamic_selector"),
      uiOutput("metric_selector"),
  # color selector removed (pseudotime removed)
      hr(),
      h5("Autoantibody filter (Aab+ donors only)"),
      checkboxGroupInput("aab_flags", NULL,
                         choices = c("GADA" = "AAb_GADA",
                                     "IA2A" = "AAb_IA2A",
                                     "ZnT8A" = "AAb_ZnT8A",
                                     "IAA" = "AAb_IAA",
                                     "mIAA" = "AAb_mIAA"),
                         selected = character(0)),
      radioButtons("aab_logic", "Match (within Aab+)", choices = c("Any", "All"), selected = "Any", inline = TRUE),
    helpText("Default: all Aab+ donors are included. Selecting one or more autoantibodies restricts only the Aab+ group to donors matching the selection (Any = at least one selected AAb, All = all selected AAbs). ND and T1D groups are unaffected."),
      checkboxGroupInput("groups", "Donor Status", choices = c("ND", "Aab+", "T1D"), selected = c("ND", "Aab+", "T1D")),
      selectInput("curve_norm", "Normalization (plot)",
                  choices = c("None (raw)" = "none", "Global z-score" = "global", "Robust per-donor" = "robust"),
                  selected = "none"),
      radioButtons("stat", "Statistic",
                   choices = c("Mean±SE" = "mean_se",
                               "Mean±SD" = "mean_sd",
                               "Median + IQR" = "median_iqr"),
                   selected = "mean_se"),
  # moved plot-specific controls below the main plot
      radioButtons("add_smooth", "Trend line", choices = c("None", "LOESS"), selected = "None", inline = TRUE),
      radioButtons("theme_bg", "Background", choices = c("Light","Dark"), selected = "Light", inline = TRUE),
      hr(),
      h5("Export"),
      downloadButton("dl_summary", "Download summary CSV"),
    downloadButton("dl_stats", "Download stats CSV")
    ),
    mainPanel(
      tabsetPanel(id = "tabs",
        tabPanel("Plot", 
                 plotlyOutput("plt", height = 650),
                 br(),
                 fluidRow(
                   column(3, sliderInput("binwidth", "Diameter bin width (µm)", min = 10, max = 100, value = 50, step = 5)),
                   column(3, sliderInput("diam_max", "Max islet diameter (µm)", min = 50, max = 1000, value = 1000, step = 10)),
                   column(2, checkboxInput("exclude_zero_top", "Exclude zero values", value = FALSE)),
                   column(2, checkboxInput("show_points", "Show individual points", value = FALSE)),
                   column(1, sliderInput("pt_size", "Point size", min = 0.3, max = 4.0, value = 0.8, step = 0.1)),
                   column(1, sliderInput("pt_alpha", "Point transparency", min = 0.05, max = 1.0, value = 0.25, step = 0.05))
                 ),
                 tags$br(),
                 tags$hr(),
                 fluidRow(column(3, checkboxInput("exclude_zero_dist", "Exclude zero values (distribution)", value = FALSE))),
                 uiOutput("dist_ui"),
                 plotlyOutput("dist", height = 500)
        ),
        tabPanel("Statistics",
                 fluidRow(
                   column(4, selectInput("alpha", "Significance level (alpha)", choices = c("0.05","0.01","0.001"), selected = "0.05")),
                   column(3, br(), actionButton("run_tests", "Run statistics"))
                 ),
                 br(),
                 tableOutput("stats_tbl"),
                 br(),
                 plotlyOutput("stats_plot", height = 320),
                 br(),
                 plotlyOutput("pairwise_plot", height = 320)
        ),
        tabPanel("Trajectory",
          tagList(
            uiOutput("traj_status"),
            fluidRow(
              column(3,
                uiOutput("traj_feature_selector")
              ),
              column(3,
                selectInput("traj_color_by", "Color points by:",
                           choices = c("Donor Status" = "donor_status", 
                                     "Donor ID" = "donor_id"),
                           selected = "donor_status")
              ),
              column(3,
                selectInput("traj_point_size", "Point size by:",
                           choices = c("Uniform" = "uniform",
                                     "Islet Diameter" = "islet_diam_um"),
                           selected = "uniform")
              ),
              column(3,
                selectInput("traj_show_trend", "Trend lines:",
                           choices = c("None" = "none", 
                                     "Overall" = "overall",
                                     "By Donor Status" = "by_donor"),
                           selected = "by_donor")
              )
            ),
            fluidRow(
              column(4,
                sliderInput("traj_alpha", "Point transparency:",
                           min = 0.1, max = 1.0, value = 0.7, step = 0.1)
              ),
              column(4),
              column(4)
            ),
            plotlyOutput("traj_scatter", height = 420),
            br(),
            fluidRow(
              column(6, h5("Donor Status Progression")),
              column(6, checkboxInput("traj_show_bins", "Show binned average", value = TRUE))
            ),
            plotOutput("traj_heatmap", height = 120),
            br(),
            fluidRow(
              column(6, h5("UMAP: Donor Status")),
              column(6, h5("UMAP: Selected Feature"))
            ),
            fluidRow(
              column(6, plotOutput("traj_umap_donor", height = 350)),
              column(6, plotOutput("traj_umap_feature", height = 350))
            ),

          )
        ),
        
        tabPanel("Viewer",
          div(style = "max-width: 1200px;",
              uiOutput("local_image_picker"),
              uiOutput("vit_view")
          )
        )
      )  # Close tabsetPanel
    )    # Close mainPanel
  )      # Close sidebarLayout
)        # Close fluidPage

# ---------- Server ----------

server <- function(input, output, session) {
  validate_file <- reactive({
    shiny::validate(shiny::need(file.exists(master_path), paste("Not found:", master_path)))
    master_path
  })

  master <- reactive({
    req(validate_file())
    load_master(master_path)
  })

  prepared <- reactive({
    pd <- prep_data(master())
    # NA audits (console)
    try({
      audit_na(pd$targets_all, "targets_all")
      audit_na(pd$markers_all, "markers_all")
      audit_na(pd$comp, "comp")
    }, silent = TRUE)
    pd
  })

  # ---------- Trajectory (AnnData) ----------
  # Try multiple locations for the H5AD; allow ADATA_PATH override
  resolve_traj_path <- function() {
    cand <- c(
      file.path("..", "..", "data", "adata_ins_root.h5ad"),
      file.path("..", "..", "scripts", "adata_ins_root.h5ad"),
      if (!is.null(project_root)) file.path(project_root, "data", "adata_ins_root.h5ad") else NULL,
      if (!is.null(project_root)) file.path(project_root, "scripts", "adata_ins_root.h5ad") else NULL,
      Sys.getenv("ADATA_PATH", unset = "")
    )
    cand <- cand[nzchar(cand)]
    for (p in cand) if (file.exists(p)) return(p)
    NULL
  }
  traj_path <- resolve_traj_path()
  traj <- reactiveVal(NULL)
  observe({
    cur <- traj_path; if (is.null(cur) || !file.exists(cur)) cur <- resolve_traj_path()
    if (is.null(cur) || !file.exists(cur)) return(NULL)
    
    # Use R anndata package instead of reticulate
    if (!requireNamespace("anndata", quietly = TRUE)) {
      traj(list(error = "R package 'anndata' not installed. Install with BiocManager::install('anndata').", error_detail = NULL))
      return(NULL)
    }
    
    # Read AnnData using R package
    ad <- NULL; load_err <- NULL
    norm_path <- tryCatch(normalizePath(cur, mustWork = TRUE), error = function(e) cur)
    ad <- tryCatch(anndata::read_h5ad(norm_path), error = function(e) { load_err <<- conditionMessage(e); NULL })
    if (is.null(ad)) {
      traj(list(error = sprintf("Failed to read H5AD: %s", norm_path), error_detail = load_err))
      return(NULL)
    }
    
    # Extract obs data.frame directly (R anndata package provides this natively)
    obs <- NULL
    try({
      obs <- ad$obs
      if (!is.null(obs) && !inherits(obs, "data.frame")) {
        obs <- tryCatch(as.data.frame(obs, stringsAsFactors = FALSE), error = function(e) NULL)
      }
    }, silent = TRUE)
    
    # Get obsm and uns if available
    obsm <- tryCatch(ad$obsm, error = function(e) NULL)
    uns  <- tryCatch(ad$uns, error = function(e) NULL)
    
    if (is.null(obs)) {
      traj(list(error = "Could not extract obs data from AnnData object", error_detail = NULL))
      return(NULL)
    }
    # Bring along imageid if present (maps to Case ID)
    cols <- intersect(c("combined_islet_id","base_islet_id","donor_status","Case ID","case_id","imageid","dpt_pseudotime","age","gender","total_cells"), colnames(obs))
    obs2 <- as.data.frame(obs[, cols, drop = FALSE])

    # --- Robust parsing of combined_islet_id to recover Case ID & islet number ---
    if (!"Case ID" %in% names(obs2)) obs2$`Case ID` <- NA_character_
    # Primary mapping from imageid if present (preferred authoritative Case ID)
    if ("imageid" %in% colnames(obs)) {
      obs2$`Case ID` <- dplyr::coalesce(obs2$`Case ID`, as.character(obs$imageid))
    }
    if ("combined_islet_id" %in% names(obs2)) {
      # Format provided: ####_Islet_## (digits, underscore, 'Islet_', digits)
      cid <- stringr::str_extract(obs2$combined_islet_id, "^[0-9]{3,4}")
      islet_num <- stringr::str_extract(obs2$combined_islet_id, "(?<=_Islet_)[0-9]{1,3}")
      obs2$`Case ID` <- dplyr::coalesce(obs2$`Case ID`, cid)
      obs2$islet_key <- ifelse(!is.na(islet_num), paste0("Islet_", islet_num), NA_character_)
    }
    # Incorporate base_islet_id where available
    if ("base_islet_id" %in% names(obs2)) {
      # base_islet_id assumed numeric; build standard key
      base_key <- paste0("Islet_", as.character(obs2$base_islet_id))
      obs2$islet_key <- dplyr::coalesce(obs2$islet_key, base_key)
    }
    # Preserve a secondary key variant (was islet_key_pref) for fallback joins
    obs2$islet_key_pref <- if ("base_islet_id" %in% names(obs2)) paste0("Islet_", as.character(obs2$base_islet_id)) else obs2$islet_key
    # Ensure Case ID is character, zero-pad to width 4 if that matches prepared() pattern
    # (We don't know required width yet; keep as-is but trim whitespace.)
    obs2$`Case ID` <- trimws(obs2$`Case ID`)
    # Add padded variant for joining if length < 4
    obs2$case_id_padded <- ifelse(nchar(obs2$`Case ID`) > 0 & nchar(obs2$`Case ID`) < 4,
                                  stringr::str_pad(obs2$`Case ID`, width = 4, pad = "0"),
                                  obs2$`Case ID`)

    # Final coercion safeguard: ensure obs2 is a data.frame before any joins
    if (!inherits(obs2, "data.frame")) {
      # Try explicit pandas conversion if available
      raw_obs <- tryCatch(ad$obs, error = function(e) NULL)
      if (!is.null(raw_obs) && reticulate::py_available(initialize = FALSE)) {
        ok_pd <- FALSE
        if (!is.null(raw_obs) && reticulate::py_has_attr(raw_obs, "to_pandas")) {
          obs_pd <- tryCatch(reticulate::py_to_r(raw_obs$to_pandas()), error = function(e) NULL)
          if (!is.null(obs_pd) && inherits(obs_pd, "data.frame")) { obs2 <- obs_pd; ok_pd <- TRUE }
        }
        if (!ok_pd && reticulate::py_has_attr(raw_obs, "to_dict")) {
          od2 <- tryCatch(raw_obs$to_dict(), error = function(e) NULL)
          if (!is.null(od2)) {
            odr2 <- tryCatch(reticulate::py_to_r(od2), error = function(e) NULL)
            if (is.list(odr2)) {
              lens <- vapply(odr2, length, integer(1)); L <- max(lens)
              if (L > 0) {
                tmp <- lapply(odr2, function(v){ if (length(v) < L) v <- c(v, rep(NA, L - length(v))); v })
                obs_try <- tryCatch(as.data.frame(tmp, stringsAsFactors = FALSE), error = function(e) NULL)
                if (!is.null(obs_try)) obs2 <- obs_try
              }
            }
          }
        }
      }
      if (!inherits(obs2, "data.frame")) {
        # Hard failure: cannot proceed
        traj(list(error = "Obs coercion failed", error_detail = "Could not convert AnnData obs to data.frame"))
        return()
      }
    }

    # Store variable names for feature selection
    var_names <- tryCatch(rownames(ad$var), error = function(e) NULL)
    
    # (Removed diameter attachment via left_join to avoid class dispatch issues; diameter will be merged later in traj_data_clean.)
    traj(list(obs = obs2, obsm = obsm, uns = uns, var_names = var_names))
  })

  # Trajectory feature selector UI
  output$traj_feature_selector <- renderUI({
    tr <- traj()
    if (is.null(tr) || !is.null(tr$error)) {
      return(selectInput("traj_feature", "Select Feature:", choices = c("No features available" = ""), selected = ""))
    }
    
    # Get available features from AnnData
    available_features <- NULL
    if (!is.null(tr$var_names)) {
      available_features <- tr$var_names
    } else {
      # Fallback - try to get from the loaded trajectory
      tryCatch({
        adata_path <- resolve_traj_path()
        if (!is.null(adata_path) && file.exists(adata_path)) {
          adata <- anndata::read_h5ad(adata_path)
          available_features <- rownames(adata$var)
        }
      }, error = function(e) NULL)
    }
    
    if (is.null(available_features) || length(available_features) == 0) {
      return(selectInput("traj_feature", "Select Feature:", choices = c("No features available" = ""), selected = ""))
    }
    
    # Organize features into categories
    hormone_markers <- intersect(c("INS", "GCG", "SST"), available_features)
    immune_markers <- intersect(c("CD3e", "CD4", "CD8a", "CD68", "CD163", "CD20", "CD45", "HLADR"), available_features)
    vascular_markers <- intersect(c("CD31", "CD34", "SMA", "ColIV"), available_features)
    neural_markers <- intersect(c("B3TUBB", "GAP43", "PGP9.5"), available_features)
    other_markers <- setdiff(available_features, c(hormone_markers, immune_markers, vascular_markers, neural_markers))
    
    # Build choices list with categories
    choices <- list()
    if (length(hormone_markers) > 0) {
      choices[["Hormone Markers"]] <- setNames(hormone_markers, hormone_markers)
    }
    if (length(immune_markers) > 0) {
      choices[["Immune Markers"]] <- setNames(immune_markers, immune_markers)
    }
    if (length(vascular_markers) > 0) {
      choices[["Vascular Markers"]] <- setNames(vascular_markers, vascular_markers)
    }
    if (length(neural_markers) > 0) {
      choices[["Neural Markers"]] <- setNames(neural_markers, neural_markers)
    }
    if (length(other_markers) > 0) {
      choices[["Other Features"]] <- setNames(other_markers, other_markers)
    }
    
    # Default selection
    default_feature <- if ("INS" %in% available_features) "INS" else available_features[1]
    
    selectInput("traj_feature", "Select Feature:", 
                choices = choices, 
                selected = default_feature)
  })
  
  # Trajectory region selector  
  output$traj_region_selector <- renderUI({
    tagList(
      selectInput("traj_region", "Region:", 
                  choices = c("Islet" = "core", 
                             "Peri-Islet" = "band", 
                             "Islet+20um" = "union"), 
                  selected = "band"),
      tags$div(style = "font-size: 11px; color: #666; margin-top: -10px;",
               "Note: Single-cell trajectory data")
    )
  })
  
  # Build trajectory measurement from AnnData features
  traj_validate <- reactive({
    selected_feature <- input$traj_feature
    if (is.null(selected_feature) || !nzchar(selected_feature)) return(NULL)
    
    # Load AnnData and validate feature exists
    adata_path <- resolve_traj_path()
    if (is.null(adata_path) || !file.exists(adata_path)) return(NULL)
    
    adata <- tryCatch(anndata::read_h5ad(adata_path), error = function(e) NULL)
    if (is.null(adata)) return(NULL)
    
    var_names <- rownames(adata$var)
    if (!selected_feature %in% var_names) return(NULL)
    
    # Return success indicator with feature name
    return(list(feature = selected_feature, status = "ready"))
  })

  # Pure base R trajectory data function - no dplyr dependencies
  traj_data_clean <- reactive({
    # Validate inputs
    ms <- traj_validate()
    if (is.null(ms)) return(NULL)
    
    selected_feature <- ms$feature
    if (is.null(selected_feature)) return(NULL)
    
    # Load AnnData with error handling
    adata_path <- resolve_traj_path()
    if (is.null(adata_path) || !file.exists(adata_path)) return(NULL)
    
    adata <- tryCatch({
      anndata::read_h5ad(adata_path)
    }, error = function(e) {
      message("[traj_data_clean] AnnData loading error: ", e$message)
      return(NULL)
    })
    
    if (is.null(adata)) return(NULL)
    
    # Validate feature exists
    var_names <- rownames(adata$var)
    feature_idx <- match(selected_feature, var_names)
    if (is.na(feature_idx)) {
      message("[traj_data_clean] Feature '", selected_feature, "' not found in variables")
      return(NULL)
    }
    
    # Extract data using base R only
    X_matrix <- adata$X
    obs_df <- as.data.frame(adata$obs)
    
    # Get expression values
    expression_vals <- X_matrix[, feature_idx]
    
    # Convert donor status using base R
    donor_raw <- as.character(obs_df$donor_status)
    donor_clean <- ifelse(donor_raw == "0", "ND",
                         ifelse(donor_raw == "1", "Aab+", 
                               ifelse(donor_raw == "2", "T1D", donor_raw)))
    
    # Extract UMAP coordinates if available
    umap_coords <- NULL
    if (!is.null(adata$obsm) && "X_umap" %in% names(adata$obsm)) {
      umap_coords <- adata$obsm$X_umap
    }
    
    # Get donor ID (imageid) for coloring
    donor_ids <- as.character(obs_df$imageid)
    
    # Build result using base R data.frame
    # Note: base_islet_id already contains 'Islet_' prefix, so use as-is
    islet_keys <- as.character(obs_df$base_islet_id)
    # Remove extra 'Islet_' prefix if it was double-added
    islet_keys <- gsub("^Islet_Islet_", "Islet_", islet_keys)
    
    result_df <- data.frame(
      case_id = as.character(obs_df$imageid),
      islet_key = islet_keys,
      combined_islet_id = as.character(obs_df$combined_islet_id),
      pt = as.numeric(obs_df$dpt_pseudotime),
      donor_status = donor_clean,
      donor_id = donor_ids,
      value = as.numeric(expression_vals),
      feature_name = selected_feature,
      stringsAsFactors = FALSE
    )
    
    # Add UMAP coordinates if available
    if (!is.null(umap_coords) && nrow(umap_coords) == nrow(result_df)) {
      result_df$umap_1 <- as.numeric(umap_coords[, 1])
      result_df$umap_2 <- as.numeric(umap_coords[, 2])
    } else {
      result_df$umap_1 <- NA_real_
      result_df$umap_2 <- NA_real_
    }
    
    # Add islet diameter using base R merge with prepared data
    prep_data <- prepared()
    if (!is.null(prep_data) && !is.null(prep_data$comp)) {
      size_lookup <- prep_data$comp[c("Case ID", "Donor Status", "islet_key", "islet_diam_um")]
      # Convert case_id format to match
      result_df$`Case ID` <- sprintf("%04d", as.numeric(result_df$case_id))
      result_df$`Donor Status` <- result_df$donor_status
      
      # Merge using base R
      merged <- merge(result_df, size_lookup, by = c("Case ID", "Donor Status", "islet_key"), all.x = TRUE)
      if (nrow(merged) > 0 && "islet_diam_um" %in% colnames(merged)) {
        result_df$islet_diam_um <- merged$islet_diam_um[match(paste(result_df$`Case ID`, result_df$`Donor Status`, result_df$islet_key), 
                                                              paste(merged$`Case ID`, merged$`Donor Status`, merged$islet_key))]
      } else {
        result_df$islet_diam_um <- NA_real_
      }
    } else {
      result_df$islet_diam_um <- NA_real_
    }
    
    # Filter valid rows using base R
    valid_idx <- which(
      is.finite(result_df$pt) & 
      !is.na(result_df$value) & 
      !is.na(result_df$case_id) & 
      !is.na(result_df$islet_key)
    )
    
    if (length(valid_idx) == 0) {
      message("[traj_data_clean] No valid rows after filtering")
      return(NULL)
    }
    
    result_df <- result_df[valid_idx, , drop = FALSE]
    
    message("[traj_data_clean] Successfully processed ", nrow(result_df), " observations for ", selected_feature)
    return(result_df)
  })

  output$traj_status <- renderUI({
    cur <- traj_path; if (is.null(cur) || !file.exists(cur)) cur <- resolve_traj_path()
    if (is.null(cur) || !file.exists(cur)) return(tags$div(style = "color:#b00;", "AnnData not found. Place 'adata_ins_root.h5ad' under data/ or scripts/, or set ADATA_PATH."))
    if (!requireNamespace("anndata", quietly = TRUE)) return(tags$div(style = "color:#b00;", "R package 'anndata' not installed. Install with BiocManager::install('anndata')."))
    tr <- traj(); if (is.null(tr)) return(tags$div(style = "color:#b00;", "AnnData object not available or H5AD failed to load."))
    if (!is.null(tr$error)) {
      det <- if (!is.null(tr$error_detail) && nzchar(tr$error_detail)) paste0(" Details: ", tr$error_detail) else ""
      install_hint <- if (grepl("anndata", tr$error, ignore.case = TRUE)) " Install with: BiocManager::install('anndata')" else ""
      return(tags$div(style = "color:#b00;", sprintf("Trajectory load error: %s%s%s", tr$error, det, install_hint)))
    }
    ms <- traj_validate(); jn <- traj_data_clean()
    nm <- if (is.null(ms)) 0 else nrow(ms)
    nj <- if (is.null(jn)) 0 else nrow(jn)
    # Distinct joined islets (case_id + islet)
    ndist <- if (is.null(jn)) 0 else nrow(dplyr::distinct(jn, case_id, islet_key))
    src <- if (!is.null(traj_path) && file.exists(traj_path)) traj_path else cur
    
    # Show PAGA trajectory info
    paga_info <- if (!is.null(tr$obs) && "dpt_pseudotime" %in% colnames(tr$obs)) {
      pt_range <- range(tr$obs$dpt_pseudotime, na.rm = TRUE)
      sprintf("PAGA trajectory (INS-rooted): %.3f - %.3f pseudotime range. ", pt_range[1], pt_range[2])
    } else {
      "No PAGA pseudotime found. "
    }
    
    tags$div(style = "color:#555;", paste0(
      paga_info,
      sprintf("Loaded AnnData (%d obs). Joined %d rows (%d unique islets) of %d islet metrics. Source: %s", 
              nrow(tr$obs), nj, ndist, nm, basename(src))
    ))
  })









  output$traj_scatter <- renderPlotly({
    df <- traj_data_clean() 
    if (is.null(df) || nrow(df) == 0) {
      return(plotly_empty() %>% layout(title = "No trajectory data available"))
    }
    
    # Determine color mapping and aesthetics
    color_by <- input$traj_color_by %||% "donor_status"
    size_by <- input$traj_point_size %||% "uniform"
    trend_type <- input$traj_show_trend %||% "by_donor"
    alpha_val <- input$traj_alpha %||% 0.7
    
    # Create base aesthetic mapping
    aes_mapping <- aes(x = pt, y = value)
    
    # Add color aesthetic
    if (color_by == "donor_status") {
      df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))
      aes_mapping$colour <- as.name("donor_status")
      color_title <- "Donor Status"
    } else if (color_by == "donor_id") {
      aes_mapping$colour <- as.name("donor_id")
      color_title <- "Donor ID"
    }
    
    # Add size aesthetic if not uniform
    if (size_by != "uniform") {
      if (size_by == "islet_diam_um" && "islet_diam_um" %in% colnames(df)) {
        # Check if diameter data is available for sizing
        if (all(is.na(df$islet_diam_um))) {
          # Use uniform size if no diameter data available
          size_title <- "Uniform (Diameter N/A)"
        } else {
          aes_mapping$size <- as.name("islet_diam_um")
          size_title <- "Islet Diameter"
        }
      }
    }
    
    # Create plot with user-controlled transparency
    g <- ggplot(df, aes_mapping) +
      geom_point(alpha = alpha_val, size = if (size_by == "uniform") 1.8 else NULL) +
      theme_minimal(base_size = 14)
    
    # Apply color scales
    if (color_by == "donor_status") {
      g <- g + scale_color_manual(values = c("ND" = "#1f77b4", "Aab+" = "#ff7f0e", "T1D" = "#d62728"),
                                  name = color_title, drop = FALSE)
    } else if (color_by == "donor_id") {
      # Use a discrete color palette for donor IDs
      g <- g + scale_color_discrete(name = color_title)
    }
    
    # Add size scale if needed
    if (size_by != "uniform" && size_by == "islet_diam_um" && "islet_diam_um" %in% colnames(df)) {
      if (!all(is.na(df$islet_diam_um))) {
        g <- g + scale_size_continuous(name = "Islet Diameter (µm)", range = c(0.8, 3.5))
      }
    }
    
    # Add trend lines based on selection
    if (trend_type == "overall") {
      g <- g + geom_smooth(method = "loess", se = TRUE, alpha = 0.2, color = "black", size = 1)
    } else if (trend_type == "by_donor") {
      # Add separate trendlines for each donor status
      g <- g + geom_smooth(aes(group = donor_status, color = donor_status), 
                          method = "loess", se = TRUE, alpha = 0.15, size = 0.8, 
                          show.legend = FALSE)  # Don't show in legend since color already shows donor status
    }
    
    # Get current selection context for labels
    selected_feature <- input$traj_feature %||% "Selected feature"
    metric_label <- "Expression Level"
    
    g <- g + labs(
      x = "Pseudotime (PAGA trajectory, INS-rooted)",
      y = paste(selected_feature, metric_label),
      title = paste("Trajectory Analysis:", selected_feature, "vs Pseudotime")
    )
    
    # Apply dark theme if selected
    if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
      g <- g + theme(
        plot.background = element_rect(fill = "#000000", colour = NA),
        panel.background = element_rect(fill = "#000000", colour = NA),
        panel.grid.major = element_line(color = "#333333"),
        axis.text = element_text(color = "#e6e6e6"),
        axis.title = element_text(color = "#f0f0f0"),
        plot.title = element_text(color = "#f0f0f0"),
        legend.text = element_text(color = "#e6e6e6"),
        legend.title = element_text(color = "#f0f0f0"),
        legend.background = element_rect(fill = "#111111", color = NA)
      )
    }
    
    ggplotly(g, tooltip = c("x", "y", "colour")) %>%
      layout(showlegend = TRUE)
  })
  
  # UMAP plot colored by donor status
  output$traj_umap_donor <- renderPlot({
    df <- traj_data_clean()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    # Check if UMAP coordinates are available
    if (all(is.na(df$umap_1)) || all(is.na(df$umap_2))) {
      return(ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = "UMAP coordinates not available", size = 5) +
        theme_void())
    }
    
    # Create UMAP plot colored by donor status
    df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))
    
    g <- ggplot(df, aes(x = umap_1, y = umap_2, color = donor_status)) +
      geom_point(alpha = 0.7, size = 1.2) +
      scale_color_manual(values = c("ND" = "#1f77b4", "Aab+" = "#ff7f0e", "T1D" = "#d62728"),
                        name = "Donor Status") +
      labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP: Donor Status") +
      theme_minimal(base_size = 12) +
      theme(aspect.ratio = 1)
    
    # Apply dark theme if selected
    if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
      g <- g + theme(
        plot.background = element_rect(fill = "#000000", colour = NA),
        panel.background = element_rect(fill = "#000000", colour = NA),
        panel.grid.major = element_line(color = "#333333"),
        axis.text = element_text(color = "#e6e6e6"),
        axis.title = element_text(color = "#f0f0f0"),
        plot.title = element_text(color = "#f0f0f0"),
        legend.text = element_text(color = "#e6e6e6"),
        legend.title = element_text(color = "#f0f0f0"),
        legend.background = element_rect(fill = "#111111", color = NA)
      )
    }
    
    return(g)
  })
  
  # UMAP plot colored by selected feature
  output$traj_umap_feature <- renderPlot({
    df <- traj_data_clean()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    # Check if UMAP coordinates are available
    if (all(is.na(df$umap_1)) || all(is.na(df$umap_2))) {
      return(ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = "UMAP coordinates not available", size = 5) +
        theme_void())
    }
    
    selected_feature <- input$traj_feature %||% "Selected feature"
    
    # Create UMAP plot colored by feature expression
    g <- ggplot(df, aes(x = umap_1, y = umap_2, color = value)) +
      geom_point(alpha = 0.7, size = 1.2) +
      labs(x = "UMAP 1", y = "UMAP 2", 
           title = paste("UMAP:", selected_feature),
           color = "Expression") +
      theme_minimal(base_size = 12) +
      theme(aspect.ratio = 1)
    
    # Use viridis or gradient color scale
    if (requireNamespace("viridisLite", quietly = TRUE)) {
      g <- g + scale_color_viridis_c(option = "plasma", na.value = "#bbbbbb")
    } else {
      g <- g + scale_color_gradient(low = "#132B43", high = "#56B1F7", na.value = "#bbbbbb")
    }
    
    # Apply dark theme if selected
    if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
      g <- g + theme(
        plot.background = element_rect(fill = "#000000", colour = NA),
        panel.background = element_rect(fill = "#000000", colour = NA),
        panel.grid.major = element_line(color = "#333333"),
        axis.text = element_text(color = "#e6e6e6"),
        axis.title = element_text(color = "#f0f0f0"),
        plot.title = element_text(color = "#f0f0f0"),
        legend.text = element_text(color = "#e6e6e6"),
        legend.title = element_text(color = "#f0f0f0"),
        legend.background = element_rect(fill = "#111111", color = NA)
      )
    }
    
    return(g)
  })

  output$traj_heatmap <- renderPlot({
    df <- traj_data_clean() 
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    show_bins <- isTRUE(input$traj_show_bins)
    
    if (show_bins) {
      # Create binned analysis
      enc <- function(s) ifelse(s == "ND", 0, ifelse(s == "Aab+", 1, ifelse(s == "T1D", 2, NA_real_)))
      df$ds_code <- enc(df$donor_status)
      
      # Normalize pseudotime to 0-1 range for binning
      pt_range <- range(df$pt, na.rm = TRUE)
      df$pt_norm <- (df$pt - pt_range[1]) / (pt_range[2] - pt_range[1])
      
      nb <- 25  # Number of bins
      brks <- seq(0, 1, length.out = nb + 1)
      df$pt_bin <- cut(pmax(0, pmin(1, df$pt_norm)), breaks = brks, include.lowest = TRUE, right = FALSE)
      
      # Calculate averages per bin using base R
      bin_levels <- levels(df$pt_bin)
      hm_list <- list()
      
      for (i in seq_along(bin_levels)) {
        bin_name <- bin_levels[i]
        bin_data <- df[!is.na(df$pt_bin) & df$pt_bin == bin_name, ]
        
        if (nrow(bin_data) >= 3) {  # Only bins with at least 3 observations
          hm_list[[length(hm_list) + 1]] <- data.frame(
            pt_bin = factor(bin_name, levels = bin_levels),
            avg_donor_status = mean(bin_data$ds_code, na.rm = TRUE),
            count = nrow(bin_data),
            stringsAsFactors = FALSE
          )
        }
      }
      
      if (length(hm_list) == 0) return(NULL)
      hm <- do.call(rbind, hm_list)
      
      if (nrow(hm) == 0) return(NULL)
      
      # Calculate bin midpoints
      mids <- head(brks, -1) + diff(brks)[1]/2
      hm$x <- mids[as.integer(hm$pt_bin)]
      
      # Create heatmap
      g <- ggplot(hm, aes(x = x, y = 1, fill = avg_donor_status)) +
        geom_tile(height = 1) +
        scale_fill_gradientn(
          colors = c("#1f77b4", "#ff7f0e", "#d62728"), 
          limits = c(0, 2), 
          na.value = "#dddddd",
          name = "Average\nDonor Type",
          breaks = c(0, 1, 2),
          labels = c("ND", "Aab+", "T1D")
        ) +
        scale_x_continuous(limits = c(0, 1), expand = c(0, 0), 
                          labels = function(x) sprintf("%.2f", x * (pt_range[2] - pt_range[1]) + pt_range[1])) +
        labs(x = "Pseudotime", y = "", title = "Donor Status Progression Along Pseudotime") +
        theme_minimal() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 12, hjust = 0.5)
        )
      
      # Apply dark theme if selected
      if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
        g <- g + theme(
          plot.background = element_rect(fill = "#000000", colour = NA),
          panel.background = element_rect(fill = "#000000", colour = NA),
          axis.text = element_text(color = "#e6e6e6"),
          axis.title = element_text(color = "#f0f0f0"),
          plot.title = element_text(color = "#f0f0f0"),
          legend.text = element_text(color = "#e6e6e6"),
          legend.title = element_text(color = "#f0f0f0"),
          legend.background = element_rect(fill = "#111111", color = NA)
        )
      }
      
      return(g)
    } else {
      # Show a simple message when binning is disabled
      return(ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = "Enable 'Show binned average' to see donor status progression", 
                size = 4, color = if (!is.null(input$theme_bg) && input$theme_bg == "Dark") "#e6e6e6" else "#333333") +
        theme_void() +
        theme(
          plot.background = element_rect(fill = if (!is.null(input$theme_bg) && input$theme_bg == "Dark") "#000000" else "#ffffff", colour = NA)
        ))
    }
  })

  

viewer_info <- reactive({
    base <- resolve_avivator_base()
    info <- list(base = base, selection = NULL, mode = "picker", ok = FALSE,
                 iframe_src = NULL, image_url = NULL)
    if (is.null(base)) return(info)
    params <- list()
    # Select image: prefer user selection from www/local_images when available; else use explicit default if set
    sel_url <- NULL
    # Selection from UI (basename only)
    sel_basename <- tryCatch(input$selected_image, silent = TRUE)
    if (!is.null(sel_basename) && nzchar(sel_basename)) {
      # Build an absolute, same-origin URL for the image to avoid any relative path resolution issues inside the iframe
      appPath <- session$clientData$url_pathname
      if (is.null(appPath) || !nzchar(appPath)) appPath <- "/"
      if (!grepl("/$", appPath)) appPath <- paste0(appPath, "/")
      sel_url <- paste0(appPath, "local_images/", sel_basename)
    } else if (!is.null(default_image_url) && nzchar(default_image_url)) {
      sel_url <- default_image_url
    }
    if (!is.null(sel_url)) {
      params[["image_url"]] <- utils::URLencode(sel_url, reserved = TRUE)
      info$image_url <- sel_url
    }
    # Channel config based on Channel_names mapping (base64-encoded per viewer expectation)
    ch_b64 <- tryCatch(build_channel_config_b64(channel_names_vec), error = function(e) NULL)
    if (!is.null(ch_b64) && nzchar(ch_b64)) {
      # URL-encode to be safe; the viewer will decodeURIComponent before atob as needed
      params[["channel_config"]] <- utils::URLencode(ch_b64, reserved = TRUE)
    }
    query <- NULL
    if (length(params)) {
      parts <- vapply(names(params), function(nm) sprintf("%s=%s", nm, params[[nm]]), character(1))
      query <- paste(parts, collapse = "&")
    }
    info$iframe_src <- if (!is.null(query)) paste0(base, "?", query) else base
    info$ok <- TRUE
    info
  })

  observe({
    if (!is.null(input$tabs) && identical(input$tabs, "Viewer")) {
      shinyjs::runjs("document.body.classList.add('viewer-mode');")
    } else {
      shinyjs::runjs("document.body.classList.remove('viewer-mode');")
    }
    
    # Hide sidebar for trajectory tab
    if (!is.null(input$tabs) && identical(input$tabs, "Trajectory")) {
      shinyjs::runjs("document.body.classList.add('trajectory-mode');")
    } else {
      shinyjs::runjs("document.body.classList.remove('trajectory-mode');")
    }
  })

  # Base per-islet dataset (no binning), aligned with current selections, without zero filtering
  raw_df_base <- reactive({
    if (!is.null(input$tabs) && identical(input$tabs, "Viewer")) {
      return(prepared()$targets_all[FALSE, , drop = FALSE])
    }
    pd <- prepared()
    groups <- input$groups
    w <- input$which
    # Resolve composition default early
    if (identical(input$mode, "Composition")) {
      if (is.null(w) || !(w %in% c("Ins_any","Glu_any","Stt_any"))) w <- "Ins_any"
    }
    if (identical(input$mode, "Targets")) {
      req(input$region, !is.null(w) && nzchar(w))
      region_tag <- paste0("islet_", tolower(input$region))
      df <- pd$targets_all %>% dplyr::filter(`Donor Status` %in% groups, class == w, tolower(type) == region_tag)
      df <- df %>% dplyr::filter(is.finite(islet_diam_um))
      df <- df %>% dplyr::mutate(value = as.numeric(if (!is.null(input$target_metric) && input$target_metric == "Counts") count else area_density))
    } else if (identical(input$mode, "Markers")) {
      req(input$region, !is.null(w) && nzchar(w))
      region_tag <- paste0("islet_", tolower(input$region))
      df <- pd$markers_all %>% dplyr::filter(`Donor Status` %in% groups, marker == w, tolower(region_type) == region_tag)
      df <- df %>% dplyr::filter(is.finite(islet_diam_um))
      if (!is.null(input$marker_metric) && input$marker_metric == "Counts") {
        df <- df %>% dplyr::mutate(value = suppressWarnings(as.numeric(pos_count)))
      } else {
        num <- suppressWarnings(as.numeric(df$pos_count))
        den <- suppressWarnings(as.numeric(df$n_cells))
        df <- df %>% dplyr::mutate(value = ifelse(is.finite(num) & is.finite(den) & den > 0, 100.0 * num / den, NA_real_))
      }
    } else {
      df <- pd$comp %>% dplyr::filter(`Donor Status` %in% groups) %>% dplyr::filter(is.finite(islet_diam_um))
      num <- suppressWarnings(as.numeric(df[[w]]))
      den <- suppressWarnings(as.numeric(df$cells_total))
      df <- df %>% dplyr::mutate(value = ifelse(is.finite(num) & is.finite(den) & den > 0, 100.0 * num / den, NA_real_))
    }
    out <- df %>% dplyr::mutate(donor_status = `Donor Status`) %>% dplyr::filter(!is.na(value))
    if (!is.null(input$diam_max) && is.finite(as.numeric(input$diam_max))) {
      out <- out %>% dplyr::filter(is.finite(islet_diam_um) & islet_diam_um <= as.numeric(input$diam_max))
    }
    # Apply AAb filters (only within the Aab+ donor group) if selected
    flags <- input$aab_flags
    if (!is.null(flags) && length(flags) > 0) {
      avail <- intersect(flags, colnames(out))
      if (length(avail) > 0) {
      # Keep ND and T1D unchanged; restrict Aab+ by selected AAbs
      others <- out %>% dplyr::filter(donor_status != "Aab+")
      aabp   <- out %>% dplyr::filter(donor_status == "Aab+")
      if (nrow(aabp) > 0) {
        mat <- as.data.frame(aabp[, avail, drop = FALSE])
        for (cc in colnames(mat)) mat[[cc]] <- as.logical(mat[[cc]])
        hits <- rowSums(mat, na.rm = TRUE)
        keep <- if (identical(input$aab_logic, "All")) hits >= length(avail) else hits >= 1
        aabp <- aabp[keep, , drop = FALSE]
      }
      out <- dplyr::bind_rows(others, aabp)
      }
    }
    out
  })

  # Raw per-islet dataset, with top-plot zero filtering applied if selected
  raw_df <- reactive({
    out <- raw_df_base()
    if (isTRUE(input$exclude_zero_top)) {
      out <- out %>% dplyr::filter(value != 0)
    }
    out
  })

  # App-wide CSS for theme background
  output$theme_css <- renderUI({
    if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
      tags$style(HTML("body { background-color: #000000; color: #e6e6e6; } .well { background-color: #111111; }"))
    } else {
      tags$style(HTML("body { background-color: #ffffff; color: #111111; }"))
    }
  })

  output$dynamic_selector <- renderUI({
    if (identical(input$mode, "Targets")) {
      # Choices should be independent of Region to preserve selection when switching regions
      classes <- prepared()$targets_all %>% pull(class) %>% unique() %>% na.omit() %>% sort()
      if (length(classes) == 0) classes <- character(0)
      # Preserve current selection if still valid
      sel <- if (!is.null(input$which) && input$which %in% classes) input$which else if (length(classes)>0) classes[1] else NULL
      selectInput("which", "Target class", choices = classes, selected = sel)
    } else if (identical(input$mode, "Markers")) {
      # Choices should be independent of Region to preserve selection when switching regions
      markers <- prepared()$markers_all %>% pull(marker) %>% unique() %>% na.omit() %>% sort()
      if (length(markers) == 0) markers <- character(0)
      sel <- if (!is.null(input$which) && input$which %in% markers) input$which else if (length(markers)>0) markers[1] else NULL
      selectInput("which", "Marker", choices = markers, selected = sel)
    } else {
      selectInput("which", "Composition measure", choices = c("Ins_frac" = "Ins_any", "Glu_frac" = "Glu_any", "Stt_frac" = "Stt_any"), selected = "Ins_any")
    }
  })

  output$region_selector <- renderUI({
    if (identical(input$mode, "Composition")) return(NULL)
    # Display labels while preserving underlying values used in code
    selectInput(
      "region", "Region",
      choices = c("Islet" = "core", "Peri-Islet" = "band", "Islet+20um" = "union"),
      selected = "band"
    )
  })

  output$metric_selector <- renderUI({
    if (identical(input$mode, "Targets")) {
      radioButtons("target_metric", "Targets metric", choices = c("Counts", "Density"), selected = "Density", inline = TRUE)
    } else if (identical(input$mode, "Markers")) {
      radioButtons("marker_metric", "Markers metric", choices = c("Counts", "% positive"), selected = "% positive", inline = TRUE)
    } else {
      NULL
    }
  })

  # color selector removed

  # Dataset for plotting: build from raw_df then bin and normalize
  plot_df <- reactive({
    df <- raw_df()
    req(nrow(df) > 0)
    # attach bins on the fly
    df <- bin_islet_sizes(df, "islet_diam_um", input$binwidth)
    # Apply normalization for main plot if requested
    norm_mode <- input$curve_norm
    if (identical(norm_mode, "robust") && all(c("Case ID") %in% colnames(df))) {
      df <- df %>% dplyr::group_by(`Case ID`) %>% dplyr::mutate(
        .med = suppressWarnings(stats::median(value, na.rm = TRUE)),
        .mad = suppressWarnings(stats::mad(value, center = .med, constant = 1, na.rm = TRUE)),
        .r_sd = ifelse(is.finite(.mad) & .mad > 0, .mad * 1.4826, NA_real_),
        value = ifelse(is.finite(.r_sd) & .r_sd > 0, (value - .med) / .r_sd, value)
      ) %>% dplyr::ungroup() %>% dplyr::select(-.med, -.mad, -`.r_sd`)
    } else if (identical(norm_mode, "global")) {
      mu <- mean(df$value, na.rm = TRUE)
      sdv <- sd(df$value, na.rm = TRUE)
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      df$value <- (df$value - mu) / sdv
    }
    df
  })

  # Summary with mean and SE per bin per group
  summary_df <- reactive({
    if (!is.null(input$tabs) && identical(input$tabs, "Viewer")) {
      return(data.frame())
    }
    df <- plot_df()
    req(nrow(df) > 0)
    summary_stats(df, group_cols = c("donor_status", "diam_bin", "diam_mid"), value_col = "value", stat = input$stat)
  })

  output$plt <- renderPlotly({
    sm <- summary_df()
    req(nrow(sm) > 0)

    # Order groups ND, Aab+, T1D and set colors
    grp_levels <- c("ND","Aab+","T1D")
    sm$donor_status <- factor(sm$donor_status, levels = grp_levels)
    color_map <- c("ND" = "#1f77b4", "Aab+" = "#ff7f0e", "T1D" = "#d62728")

    # Build y-label reflecting normalization choice
    ylab_base <- if (identical(input$mode, "Targets")) {
      if (!is.null(input$target_metric) && input$target_metric == "Counts") "Target count" else "Target density (per µm²)"
    } else if (identical(input$mode, "Markers")) {
      if (!is.null(input$marker_metric) && input$marker_metric == "Counts") "Positive cell count" else "n positive / total (%)"
    } else {
      "% composition"
    }
    cnorm <- if (!is.null(input$curve_norm) && length(input$curve_norm) == 1 && nzchar(input$curve_norm)) input$curve_norm else "none"
    ylab <- switch(cnorm,
                   none = ylab_base,
                   global = paste0(ylab_base, " (scaled z)"),
                   robust = paste0(ylab_base, " (robust z)"),
                   ylab_base)

    # Build title reflecting selection; for Composition, indicate fraction (%)
    title_text <- if (identical(input$mode, "Targets")) {
      paste0(input$which, " vs islet size")
    } else if (identical(input$mode, "Markers")) {
      paste0(input$which, " vs islet size")
    } else {
      wsel <- input$which
      if (is.null(wsel) || length(wsel) != 1 || !nzchar(wsel)) wsel <- "Ins_any"
      nm <- switch(wsel,
                   Ins_any = "Insulin+ fraction",
                   Glu_any = "Glucagon+ fraction",
                   Stt_any = "Somatostatin+ fraction",
                   wsel)
      paste0(nm, " vs islet size")
    }

    p <- ggplot(sm, aes(x = diam_mid, y = y, color = donor_status, group = donor_status)) +
      geom_line(alpha = ifelse(!is.null(input$add_smooth) && input$add_smooth == "LOESS", 0, 1)) +
      geom_point() +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) +
   labs(x = "Islet diameter (µm)",
     y = ylab,
     color = "Donor Status",
     title = title_text) +
      scale_color_manual(values = color_map, breaks = grp_levels, drop = FALSE) +
      theme_minimal(base_size = 14)

    if (isTRUE(input$show_points)) {
      raw <- plot_df()
      raw$donor_status <- factor(raw$donor_status, levels = grp_levels)
      jw <- max(1, as.numeric(input$binwidth) * 0.35)
      # Add slight vertical jitter to avoid rows when plotting counts
      yr <- suppressWarnings(range(raw$value, na.rm = TRUE))
      ydiff <- if (all(is.finite(yr))) diff(yr) else 0
      is_counts <- (identical(input$mode, "Targets") && !is.null(input$target_metric) && input$target_metric == "Counts") ||
                   (identical(input$mode, "Markers") && !is.null(input$marker_metric) && input$marker_metric == "Counts")
      jh <- if (is_counts) {
        # 2% of range or at least 0.5 in data units
        mx <- max(0.02 * ydiff, 0.5)
        if (!is.finite(mx) || mx <= 0) 0.5 else mx
      } else {
        # 1% of range
        mx <- 0.01 * ydiff
        if (!is.finite(mx) || mx < 0) 0 else mx
      }
      p <- p +
        geom_point(data = raw, aes(x = diam_mid, y = value, color = donor_status),
                   position = position_jitter(width = jw, height = jh), size = input$pt_size, alpha = input$pt_alpha, inherit.aes = FALSE)
    }

    # Axis breaks and minor ticks
    xmax <- suppressWarnings(max(sm$diam_mid, na.rm = TRUE))
    if (!is.finite(xmax)) xmax <- 300
    major_breaks <- seq(0, ceiling(xmax/50)*50, by = 50)
    minor_breaks <- seq(0, ceiling(xmax/10)*10, by = 10)
    p <- p + scale_x_continuous(breaks = major_breaks, minor_breaks = minor_breaks)
    p <- p + theme(panel.grid.minor = element_line(size = 0.2, colour = if (!is.null(input$theme_bg) && input$theme_bg == "Dark") "#222222" else "#eeeeee"))

    # Highlight significant bins if stats were run
    st <- NULL
    if (!is.null(st) && nrow(st) > 0) {
      alpha_num <- as.numeric(input$alpha)
      st <- st %>% mutate(p_adj = p.adjust(p_value, method = "BH"), sig = p_adj <= alpha_num)
      sig <- st %>% filter(sig)
      if (nrow(sig) > 0) {
        ytop <- suppressWarnings(max(sm$ymax, na.rm = TRUE))
        ytop <- ifelse(is.finite(ytop), ytop, suppressWarnings(max(sm$y, na.rm = TRUE)))
        sig$y <- ytop * 1.05
        p <- p + geom_text(data = sig, aes(x = mid, y = y), label = ifelse(unique(st$test) == 'kendall', '816', '*'), inherit.aes = FALSE, color = "#444444")
      }
    }

    # Optional smoothing overlay (kept subtle)
    if (!is.null(input$add_smooth) && input$add_smooth == "LOESS") {
      p <- p + geom_smooth(se = FALSE, method = "loess", span = 0.6, size = 0.9)
    }

    # Apply dark theme inside the plot if selected
    if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
      p <- p + theme(
        plot.background = element_rect(fill = "#000000", colour = NA),
        panel.background = element_rect(fill = "#000000", colour = NA),
        panel.grid.major = element_line(color = "#333333"),
        axis.text = element_text(color = "#e6e6e6"),
        axis.title = element_text(color = "#f0f0f0"),
        plot.title = element_text(color = "#f0f0f0"),
        legend.text = element_text(color = "#e6e6e6"),
        legend.title = element_text(color = "#f0f0f0")
      )
    }

    gg <- ggplotly(p)
    gg <- gg %>% layout(legend = list(orientation = "h", x = 0, y = -0.15))
    gg
  })

  # Stats: continuous tests (no binning). Global LM + post-hoc pairwise on residuals
  stats_run <- eventReactive(input$run_tests, {
    rdf <- raw_df()
    if (is.null(rdf) || !nrow(rdf)) return(NULL)
    # Ensure donor_status factor ordering
    rdf$donor_status <- factor(rdf$donor_status, levels = c("ND","Aab+","T1D"))
    # Global: value ~ donor_status + islet_diam_um (Type I with donor_status first)
    fit <- tryCatch(lm(value ~ donor_status + islet_diam_um, data = rdf), error = function(e) NULL)
    p_global <- NA_real_
    if (!is.null(fit)) {
      at <- tryCatch(anova(fit), error = function(e) NULL)
      if (!is.null(at)) p_global <- suppressWarnings(as.numeric(at[["Pr(>F)"]][1]))
    }
    # Post-hoc: residualize value ~ islet_diam_um, then pairwise t-tests across groups
    res <- tryCatch({
      fit_res <- lm(value ~ islet_diam_um, data = rdf)
      r <- resid(fit_res)
      pt <- pairwise.t.test(r, rdf$donor_status, p.adjust.method = "BH")
      mat <- as.data.frame(as.table(pt$p.value))
      colnames(mat) <- c("group1","group2","p_value")
      mat <- mat[!is.na(mat$p_value), , drop = FALSE]
      mat
    }, error = function(e) NULL)
    # Assemble results
    global <- data.frame(type = "global", contrast = "donor_status", p_value = p_global, stringsAsFactors = FALSE)
    if (!is.null(res) && nrow(res)) {
      pairs <- data.frame(type = "pairwise",
                          contrast = paste(res$group1, "vs", res$group2),
                          p_value = res$p_value,
                          stringsAsFactors = FALSE)
      out <- rbind(global, pairs)
    } else {
      out <- global
    }
    out
  }, ignoreInit = TRUE)

  stats_data <- reactive({
    if (!is.null(input$tabs) && identical(input$tabs, "Viewer")) return(NULL)
    st <- stats_run()
    if (is.null(st) || nrow(st) == 0) return(NULL)
    alpha_num <- as.numeric(input$alpha)
    st$p_adj <- NA_real_
    if (any(st$type == "pairwise")) {
      idx <- which(st$type == "pairwise")
      st$p_adj[idx] <- p.adjust(st$p_value[idx], method = "BH")
    }
    st$sig <- ifelse(!is.na(st$p_adj), st$p_adj <= alpha_num, st$p_value <= alpha_num)
    st
  })

  output$stats_tbl <- renderTable({
    if (!is.null(input$tabs) && identical(input$tabs, "Viewer")) return(NULL)
    st <- stats_data()
    if (is.null(st) || nrow(st) == 0) return(NULL)
    fmt <- function(x) ifelse(is.na(x), NA_character_, formatC(x, format = "e", digits = 2))
    display <- st
    display$p_value <- fmt(display$p_value)
    display$p_adj <- fmt(display$p_adj)
    display
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$stats_plot <- renderPlotly({
    if (!is.null(input$tabs) && identical(input$tabs, "Viewer")) return(NULL)
    st <- stats_data()
    req(st, nrow(st) > 0)
    alpha_num <- as.numeric(input$alpha)
    df <- st %>%
      mutate(neglog = ifelse(!is.na(p_value) & p_value > 0, -log10(p_value), NA_real_)) %>%
      filter(!is.na(neglog))
    if (nrow(df) == 0) return(NULL)
    thresh <- -log10(alpha_num)
    g <- ggplot(df, aes(x = reorder(contrast, neglog), y = neglog, fill = type)) +
      geom_col(show.legend = TRUE) +
      geom_hline(yintercept = thresh, linetype = "dashed", color = "#444444") +
      coord_flip() +
      labs(x = "Contrast", y = expression(-log[10](p)), fill = "Test type",
           title = "Test significance overview") +
      theme_minimal(base_size = 14)
    ggplotly(g)
  })

  output$pairwise_plot <- renderPlotly({
    if (!is.null(input$tabs) && identical(input$tabs, "Viewer")) return(NULL)
    st <- stats_data()
    req(st, nrow(st) > 0)
    df <- st %>% filter(type == "pairwise" & !is.na(p_value))
    if (nrow(df) == 0) return(NULL)
    alpha_num <- as.numeric(input$alpha)
    df <- df %>% mutate(sig = ifelse(is.na(p_adj), p_value <= alpha_num, p_adj <= alpha_num))
    g <- ggplot(df, aes(x = reorder(contrast, p_value), y = p_value, color = sig)) +
      geom_point(size = 3) +
      geom_hline(yintercept = alpha_num, linetype = "dashed", color = "#444444") +
      scale_y_log10() +
      scale_color_manual(values = c(`FALSE` = "#1f77b4", `TRUE` = "#d62728"),
                         labels = c(`FALSE` = "NS", `TRUE` = "Significant")) +
      coord_flip() +
      labs(x = "Pairwise contrast", y = "p-value (log scale)", color = "", title = "Pairwise comparison p-values") +
      theme_minimal(base_size = 14)
    ggplotly(g)
  })

  # Vitessce embed
  # Image picker for local OME-TIFFs under www/local_images (preferred) or LOCAL_IMAGE_ROOT fallback
  output$local_image_picker <- renderUI({
    # Only render when on Viewer tab
    if (is.null(input$tabs) || input$tabs != "Viewer") {
      return(NULL)
    }
    
    vi <- viewer_info()
    if (is.null(vi$base)) {
      return(tagList(
        tags$div(style = "color:#b00;", "Local Avivator static build not found under shiny/www/avivator."),
        tags$div(style = "color:#666; font-size:90%;", "To install: run scripts/install_avivator.sh (requires Node ≥ 18) or place a prebuilt bundle under shiny/www/avivator.")
      ))
    }
    # List candidate images from www/local_images (static) first
    pat <- ".*[.](ome[.]tif{1,2}|tif{1,2})$"
    files <- character(0)
    if (has_www_local_images) {
      files <- tryCatch(list.files(www_local_images_dir, pattern = pat, ignore.case = TRUE, recursive = FALSE), error = function(e) character(0))
    }
    # Fallback to LOCAL_IMAGE_ROOT if www/local_images is empty or missing
    if (length(files) == 0 && !is.null(local_images_root) && dir.exists(local_images_root)) {
      basefiles <- tryCatch(list.files(local_images_root, pattern = pat, ignore.case = TRUE, recursive = FALSE), error = function(e) character(0))
      # show only basenames; these will be addressed via /local_images when no www folder
      files <- basefiles
    }
    if (length(files) == 0) {
      return(tagList(
        tags$div(style = "margin-bottom:6px;", "No local images found under www/local_images or LOCAL_IMAGE_ROOT."),
        tags$div(style = "color:#666; font-size:90%;", "Place .ome.tif(f) files under app/shiny_app/www/local_images on the server for best performance (HTTP Range support).")
      ))
    }
    # Preserve selection
    sel <- if (!is.null(input$selected_image) && input$selected_image %in% files) input$selected_image else files[[1]]
    # Build an absolute URL for range check from the main app page
    appPath <- session$clientData$url_pathname
    if (is.null(appPath) || !nzchar(appPath)) appPath <- "/"
    if (!grepl("/$", appPath)) appPath <- paste0(appPath, "/")
    check_url <- if (has_www_local_images) paste0(appPath, "local_images/", sel) else paste0("/local_images/", sel)

    tagList(
      div(style = "display:flex; gap:12px; align-items:flex-end; flex-wrap:wrap;",
        div(style = "min-width:280px;", selectInput("selected_image", "Select OME-TIFF", choices = files, selected = sel)),
        div(style = "color:#666; font-size:90%; padding-bottom:8px;", "Images are served from ", if (has_www_local_images) "www/local_images (static)" else "/local_images (dynamic)" )
      )
      , tags$div(id = "range-check", style = "margin-top:6px; font-family:monospace; color:#555;", "Checking Range support...")
      , tags$script(HTML(paste0(
        "(function(){\n",
        "  const el = document.getElementById('range-check');\n",
        "  if (!el) return;\n",
        "  el.textContent = 'Checking Range support...';\n",
        "  fetch('", check_url, "', { headers: { 'Range': 'bytes=0-1' } }).then(r => {\n",
        "    const cr = r.headers.get('content-range');\n",
        "    const ar = r.headers.get('accept-ranges');\n",
        "    el.textContent = 'HTTP ' + r.status + ' | Accept-Ranges=' + ar + ' | Content-Range=' + cr;\n",
        "  }).catch(e => {\n",
        "    el.textContent = 'Fetch error: ' + e;\n",
        "  });\n",
        "})();"
      )))
    )
  })

output$vit_view <- renderUI({
  # Only render when on Viewer tab
  if (is.null(input$tabs) || input$tabs != "Viewer") {
    return(NULL)
  }
  
  vi <- viewer_info()
  if (is.null(vi$base)) {
    return(tagList(
      tags$div(style = "color:#b00;", "Local Avivator static build not found under shiny/www/avivator."),
      tags$div(style = "color:#666; font-size:90%;", "To install: run scripts/install_avivator.sh (requires Node ≥ 18) or place a prebuilt bundle under shiny/www/avivator.")
    ))
  }
  tagList(
    tags$iframe(src = vi$iframe_src,
              width = "100%",
              frameBorder = 0,
              allowfullscreen = NA,
              style = "width:100%; height:calc(100vh - 140px); border:0;"),
    tags$div(style = "margin-top:6px;",
      tags$a(href = vi$iframe_src, target = "_blank", rel = "noopener", "Open viewer in a new tab")
    )
  )
})

  # Pseudotime scatter when counts are selected (scaled counts vs pseudotime)
  # Distribution plot UI (violin/box)
  output$dist_ui <- renderUI({
    tagList(
      tags$h4("Distribution by donor status"),
      fluidRow(
        column(4,
          radioButtons("dist_type", "Plot type", choices = c("Violin", "Box"), selected = "Violin", inline = TRUE)
        ),
        column(4,
          checkboxInput("dist_show_points", "Show jittered points", value = TRUE)
        ),
        column(4,
          sliderInput("dist_pt_alpha", "Point transparency", min = 0.05, max = 1.0, value = 0.25, step = 0.05)
        )
      ),
      fluidRow(
        column(4,
          sliderInput("dist_pt_size", "Point size", min = 0.3, max = 4.0, value = 0.7, step = 0.1)
        )
      )
    )
  })

  # Distribution plot: violin/box per donor group for current raw_df()
  output$dist <- renderPlotly({
    rdf <- raw_df_base()
    if (is.null(rdf) || nrow(rdf) == 0) return(NULL)
    if (isTRUE(input$exclude_zero_dist)) {
      rdf <- rdf %>% dplyr::filter(value != 0)
    }
    # Apply normalization consistent with main plot
    norm_mode <- input$curve_norm
    if (identical(norm_mode, "robust") && all(c("Case ID") %in% colnames(rdf))) {
      rdf <- rdf %>% dplyr::group_by(`Case ID`) %>% dplyr::mutate(
        .med = suppressWarnings(stats::median(value, na.rm = TRUE)),
        .mad = suppressWarnings(stats::mad(value, center = .med, constant = 1, na.rm = TRUE)),
        .r_sd = ifelse(is.finite(.mad) & .mad > 0, .mad * 1.4826, NA_real_),
        value = ifelse(is.finite(.r_sd) & .r_sd > 0, (value - .med) / .r_sd, value)
      ) %>% dplyr::ungroup() %>% dplyr::select(-.med, -.mad, -`.r_sd`)
    } else if (identical(norm_mode, "global")) {
      mu <- mean(rdf$value, na.rm = TRUE)
      sdv <- sd(rdf$value, na.rm = TRUE)
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      rdf$value <- (rdf$value - mu) / sdv
    }
    # Ensure factor order
    rdf$donor_status <- factor(rdf$donor_status, levels = c("ND","Aab+","T1D"))
    ylab_base <- if (identical(input$mode, "Targets")) {
      if (!is.null(input$target_metric) && input$target_metric == "Counts") "Target count" else "Target density (per µm²)"
    } else if (identical(input$mode, "Markers")) {
      if (!is.null(input$marker_metric) && input$marker_metric == "Counts") "Positive cell count" else "n positive / total (%)"
    } else {
      "% composition"
    }
    ylab <- switch(input$curve_norm,
                   none = ylab_base,
                   global = paste0(ylab_base, " (scaled z)"),
                   robust = paste0(ylab_base, " (robust z)"),
                   ylab_base)
    g <- ggplot(rdf, aes(x = donor_status, y = value, fill = donor_status))
    if (!is.null(input$dist_type) && input$dist_type == "Box") {
      g <- g + geom_boxplot(outlier.shape = NA, alpha = 0.8)
    } else {
      g <- g + geom_violin(trim = FALSE, alpha = 0.8)
    }
    if (isTRUE(input$dist_show_points)) {
      g <- g + geom_jitter(width = 0.18,
                           alpha = ifelse(is.null(input$dist_pt_alpha), 0.25, input$dist_pt_alpha),
                           size = ifelse(is.null(input$dist_pt_size), 0.7, input$dist_pt_size))
    }
    g <- g + scale_fill_manual(values = c("ND" = "#1f77b4", "Aab+" = "#ff7f0e", "T1D" = "#d62728"), guide = "none") +
      labs(x = "Donor Status", y = ylab, title = "Distribution across donor groups") +
      theme_minimal(base_size = 14)
    if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
      g <- g + theme(
        plot.background = element_rect(fill = "#000000", colour = NA),
        panel.background = element_rect(fill = "#000000", colour = NA),
        panel.grid.major = element_line(color = "#333333"),
        axis.text = element_text(color = "#e6e6e6"),
        axis.title = element_text(color = "#f0f0f0"),
        plot.title = element_text(color = "#f0f0f0"),
        legend.text = element_text(color = "#e6e6e6"),
        legend.title = element_text(color = "#f0f0f0")
      )
    }
    ggplotly(g)
  })

  # pseudotime stats removed

  output$dl_summary <- downloadHandler(
    filename = function() {
      paste0("summary_", gsub("[^0-9A-Za-z]+","_", Sys.time()), ".csv")
    },
    content = function(file) {
      df <- summary_df()
      if (is.null(df) || nrow(df) == 0) df <- data.frame()
      write.csv(df, file, row.names = FALSE)
    }
  )

  output$dl_stats <- downloadHandler(
    filename = function() {
      paste0("stats_", gsub("[^0-9A-Za-z]+","_", Sys.time()), ".csv")
    },
    content = function(file) {
      st <- stats_data()
      if (is.null(st) || nrow(st) == 0) st <- data.frame()
      write.csv(st, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
