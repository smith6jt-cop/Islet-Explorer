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
if (!is.null(local_images_root)) {
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
      # normalize type values to expected tags used throughout
      type = dplyr::case_when(
        tolower(type) %in% c("islet_core", "core") ~ "islet_core",
        tolower(type) %in% c("islet_band", "band", "peri-islet", "peri_islet") ~ "islet_band",
        tolower(type) %in% c("islet_union", "union", "islet+20um", "islet_20um") ~ "islet_union",
        TRUE ~ tolower(type)
      )
    ) %>%
    dplyr::left_join(size_area, by = c("Case ID", "Donor Status", "islet_key"))

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
      dplyr::left_join(have_union, by = keys) %>%
      dplyr::filter(is.na(.has_union))
    if (nrow(union_missing) > 0) {
      synth <- union_missing %>%
        dplyr::transmute(`Case ID`, `Donor Status`, islet_key, type = "islet_union", class,
                          area_um2 = NA_real_, region_um2 = NA_real_,
                          area_density = NA_real_,
                          count = suppressWarnings(as.numeric(count_core)) + suppressWarnings(as.numeric(count_band)))
      # Attach diameter via size_area
      synth <- synth %>% dplyr::left_join(size_area, by = c("Case ID", "Donor Status", "islet_key"))
      # Bind synthetic union rows to targets_all
      targets_all <- dplyr::bind_rows(targets_all, synth)
    }
  }
  # Ensure AAb flags are present for all rows, including synthetic ones
  if (!is.null(donors_meta)) {
    targets_all <- targets_all %>% dplyr::select(-dplyr::any_of(aab_cols_all)) %>%
      dplyr::left_join(donors_meta, by = c("Case ID","Donor Status"))
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
      dplyr::left_join(donors_meta, by = c("Case ID","Donor Status"))
  }
  markers_all <- markers_all %>% left_join(size_area, by = c("Case ID", "Donor Status", "islet_key"))

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

    # Backfill missing band rows when core and union exist but band is missing: band = union - core (counts only)
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
  comp <- comp %>% left_join(size_area, by = c("Case ID", "Donor Status", "islet_key"))

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
  tags$head(tags$style(HTML("\n    body.viewer-mode div.col-sm-4,\n    body.viewer-mode div.col-sm-3,\n    body.viewer-mode div.col-lg-3 {\n      display: none !important;\n    }\n    body.viewer-mode div.col-sm-8,\n    body.viewer-mode div.col-lg-9 {\n      width: 100% !important;\n      max-width: 100% !important;\n      flex: 0 0 100%;\n    }\n    body.viewer-mode .tab-content {\n      padding-left: 0 !important;\n      padding-right: 0 !important;\n    }\n  "))),
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
        tabPanel("Data QA",
                 h4("Per-donor data integrity checks"),
                 p("Each donor should have: (i) the same number of islets across comp/markers/targets, and (ii) 3× as many distinct region rows in markers and targets (core, band, union) as in composition."),
                 tableOutput("qa_table"),
                 br(),
                 h4("Mismatches only"),
                 tableOutput("qa_mismatch")
        ),
        tabPanel("Viewer",
          div(style = "max-width: 1200px;",
              uiOutput("local_image_picker"),
              uiOutput("vit_view")
          )
        )
      )
    )
  )
)

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

  # QA summary per donor: comp islets vs markers/targets distinct region rows
  qa_summary <- reactive({
    pd <- prepared()
    if (is.null(pd)) return(NULL)
    comp_islets <- pd$comp %>% dplyr::distinct(`Case ID`, `Donor Status`, islet_key)
    n_islets <- comp_islets %>% dplyr::count(`Case ID`, `Donor Status`, name = "n_islets")
    mk_regions <- pd$markers_all %>% dplyr::distinct(`Case ID`, `Donor Status`, islet_key, region_type) %>%
      dplyr::count(`Case ID`, `Donor Status`, name = "n_mk_regions")
    tg_regions <- pd$targets_all %>% dplyr::distinct(`Case ID`, `Donor Status`, islet_key, type) %>%
      dplyr::count(`Case ID`, `Donor Status`, name = "n_tg_regions")
    dd <- n_islets %>% dplyr::left_join(mk_regions, by = c("Case ID","Donor Status")) %>%
      dplyr::left_join(tg_regions, by = c("Case ID","Donor Status")) %>%
      dplyr::mutate(expected = 3 * n_islets,
                    ok_markers = (n_mk_regions == expected),
                    ok_targets = (n_tg_regions == expected)) %>%
      dplyr::arrange(`Donor Status`, `Case ID`)
    dd
  })

  output$qa_table <- renderTable({
    qs <- qa_summary()
    if (is.null(qs) || nrow(qs) == 0) return(NULL)
    qs
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$qa_mismatch <- renderTable({
    qs <- qa_summary()
    if (is.null(qs) || nrow(qs) == 0) return(NULL)
    mm <- qs %>% dplyr::filter(!(ok_markers & ok_targets))
    if (nrow(mm) == 0) return(data.frame(message = "No mismatches: all donors satisfy 3× rule for markers and targets."))
    mm
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  viewer_info <- reactive({
    base <- resolve_avivator_base()
    info <- list(base = base, selection = NULL, mode = "picker", ok = FALSE,
                 iframe_src = NULL, image_url = NULL)
    if (is.null(base)) return(info)
    params <- list()
    # Select image: Use explicit default only; otherwise let Avivator use its packaged default demos
    sel_url <- NULL
    if (!is.null(default_image_url) && nzchar(default_image_url)) {
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
    ylab <- switch(input$curve_norm,
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
      nm <- switch(input$which,
                   Ins_any = "Insulin+ fraction",
                   Glu_any = "Glucagon+ fraction",
                   Stt_any = "Somatostatin+ fraction",
                   input$which)
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
  output$local_image_picker <- renderUI({
    vi <- viewer_info()
    if (is.null(vi$base)) {
      return(tagList(
        tags$div(style = "color:#b00;", "Local Avivator static build not found under shiny/www/avivator."),
        tags$div(style = "color:#666; font-size:90%;", "To install: run scripts/install_avivator.sh (requires Node ≥ 18) or place a prebuilt bundle under shiny/www/avivator.")
      ))
    }
    # Suppress any 'Viewing:' text by default
    NULL
  })

output$vit_view <- renderUI({
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