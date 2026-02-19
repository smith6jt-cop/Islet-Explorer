trajectory_server <- function(id, prepared, selected_islet, forced_image) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---------- H5AD Path Resolution ----------
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

    # ---------- Reactive state ----------
    traj <- reactiveVal(NULL)
    traj_load_attempted <- reactiveVal(FALSE)
    traj_coord_lookup <- reactiveVal(NULL)
    traj_outliers <- reactiveVal(NULL)

    # ---------- Lazy H5AD loader ----------
    observe({
      if (traj_load_attempted()) return()

      cur <- traj_path
      if (is.null(cur) || !file.exists(cur)) cur <- resolve_traj_path()
      if (is.null(cur) || !file.exists(cur)) return()

      traj_load_attempted(TRUE)

      if (!requireNamespace("anndata", quietly = TRUE)) {
        traj(list(error = "R package 'anndata' not installed.", error_detail = NULL))
        return()
      }

      ad <- NULL
      load_err <- NULL
      norm_path <- tryCatch(normalizePath(cur, mustWork = TRUE), error = function(e) cur)

      ad <- tryCatch(
        anndata::read_h5ad(norm_path),
        error = function(e) {
          load_err <<- conditionMessage(e)
          NULL
        }
      )

      if (is.null(ad)) {
        traj(list(error = sprintf("Failed to read H5AD: %s", norm_path), error_detail = load_err))
        return()
      }

      # Extract obs, obsm, uns
      obs  <- tryCatch(ad$obs,  error = function(e) NULL)
      obsm <- tryCatch(ad$obsm, error = function(e) NULL)
      uns  <- tryCatch(ad$uns,  error = function(e) NULL)
      var_names <- tryCatch(rownames(ad$var), error = function(e) NULL)

      if (is.null(obs)) {
        traj(list(error = "Could not extract obs data from AnnData", error_detail = NULL))
        return()
      }

      # Ensure obs is a data.frame
      if (!inherits(obs, "data.frame")) {
        obs <- tryCatch(as.data.frame(obs, stringsAsFactors = FALSE), error = function(e) NULL)
        if (is.null(obs)) {
          traj(list(error = "Obs coercion failed", error_detail = NULL))
          return()
        }
      }

      # Extract required columns
      cols <- intersect(
        c("combined_islet_id", "base_islet_id", "donor_status", "Case ID",
          "case_id", "imageid", "dpt_pseudotime", "age", "gender", "total_cells"),
        colnames(obs)
      )
      obs2 <- as.data.frame(obs[, cols, drop = FALSE])

      # Parse combined_islet_id for Case ID and islet_key
      if (!"Case ID" %in% names(obs2)) obs2$`Case ID` <- NA_character_

      if ("imageid" %in% names(obs2)) {
        obs2$`Case ID` <- dplyr::coalesce(obs2$`Case ID`, as.character(obs2$imageid))
      }

      if ("combined_islet_id" %in% names(obs2)) {
        cid <- stringr::str_extract(obs2$combined_islet_id, "^[0-9]{3,4}")
        islet_num <- stringr::str_extract(obs2$combined_islet_id, "(?<=_Islet_)[0-9]{1,3}")
        obs2$`Case ID` <- dplyr::coalesce(obs2$`Case ID`, cid)
        obs2$islet_key <- ifelse(!is.na(islet_num), paste0("Islet_", islet_num), NA_character_)
      }

      # Normalize donor_status
      if ("donor_status" %in% names(obs2)) {
        obs2$donor_status <- stringr::str_trim(obs2$donor_status)
        obs2$donor_status <- dplyr::case_when(
          tolower(obs2$donor_status) %in% c("nd", "non-diabetic", "nondiabetic") ~ "ND",
          tolower(obs2$donor_status) %in% c("aab+", "aab", "autoantibody") ~ "Aab+",
          tolower(obs2$donor_status) %in% c("t1d", "type1", "diabetic") ~ "T1D",
          TRUE ~ obs2$donor_status
        )
      }

      # Create generic pseudotime column (default to PAGA dpt_pseudotime)
      if ("dpt_pseudotime" %in% names(obs2)) {
        obs2$pseudotime <- obs2$dpt_pseudotime
      } else {
        obs2$pseudotime <- NA_real_
      }

      traj(list(obs = obs2, obsm = obsm, uns = uns, var_names = var_names, adata = ad))
    })

    # ---------- Pseudotime reactive ----------
    traj_with_pseudotime <- reactive({
      tr <- traj()
      if (is.null(tr) || !is.null(tr$error)) return(tr)

      obs <- tr$obs
      ad  <- tr$adata

      # Always use PAGA pseudotime (default)
      if ("dpt_pseudotime" %in% names(obs)) {
        obs$pseudotime <- obs$dpt_pseudotime
      }

      list(obs = obs, obsm = tr$obsm, uns = tr$uns, var_names = tr$var_names, adata = ad)
    })

    # ---------- Feature selector UI ----------
    output$traj_feature_selector <- renderUI({
      tr <- traj()
      if (is.null(tr) || !is.null(tr$error)) {
        return(selectInput(ns("traj_feature"), "Select Feature:",
                           choices = c("No features available" = ""), selected = ""))
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
        return(selectInput(ns("traj_feature"), "Select Feature:",
                           choices = c("No features available" = ""), selected = ""))
      }

      # Organize features into categories
      hormone_markers  <- intersect(c("INS", "GCG", "SST"), available_features)
      immune_markers   <- intersect(c("CD3e", "CD4", "CD8a", "CD68", "CD163", "CD20", "CD45", "HLADR"), available_features)
      vascular_markers <- intersect(c("CD31", "CD34", "SMA", "ColIV"), available_features)
      neural_markers   <- intersect(c("B3TUBB", "GAP43", "PGP9.5"), available_features)
      other_markers    <- setdiff(available_features, c(hormone_markers, immune_markers, vascular_markers, neural_markers))

      # Build choices list with categories
      choices <- list()
      if (length(hormone_markers)  > 0) choices[["Hormone Markers"]]  <- setNames(hormone_markers, hormone_markers)
      if (length(immune_markers)   > 0) choices[["Immune Markers"]]   <- setNames(immune_markers, immune_markers)
      if (length(vascular_markers) > 0) choices[["Vascular Markers"]] <- setNames(vascular_markers, vascular_markers)
      if (length(neural_markers)   > 0) choices[["Neural Markers"]]   <- setNames(neural_markers, neural_markers)
      if (length(other_markers)    > 0) choices[["Other Features"]]   <- setNames(other_markers, other_markers)

      # Default selection
      default_feature <- if ("INS" %in% available_features) "INS" else available_features[1]

      selectInput(ns("traj_feature"), "Select Feature:",
                  choices = choices,
                  selected = default_feature)
    })

    # ---------- Validate selected feature ----------
    traj_validate <- reactive({
      selected_feature <- input$traj_feature
      if (is.null(selected_feature) || !nzchar(selected_feature)) return(NULL)

      adata_path <- resolve_traj_path()
      if (is.null(adata_path) || !file.exists(adata_path)) return(NULL)

      adata <- tryCatch(anndata::read_h5ad(adata_path), error = function(e) NULL)
      if (is.null(adata)) return(NULL)

      var_names <- rownames(adata$var)
      if (!selected_feature %in% var_names) return(NULL)

      list(feature = selected_feature, status = "ready")
    })

    # ---------- Main data pipeline ----------
    traj_data_clean <- reactive({
      ms <- traj_validate()
      if (is.null(ms)) return(NULL)

      selected_feature <- ms$feature
      if (is.null(selected_feature)) return(NULL)

      adata_path <- resolve_traj_path()
      if (is.null(adata_path) || !file.exists(adata_path)) return(NULL)

      tr <- traj_with_pseudotime()
      if (is.null(tr) || !is.null(tr$error)) return(NULL)

      # Use cached adata object instead of re-reading
      adata <- tr$adata
      if (is.null(adata)) {
        adata <- tryCatch(anndata::read_h5ad(adata_path), error = function(e) NULL)
        if (is.null(adata)) return(NULL)
      }

      var_names <- rownames(adata$var)
      if (!selected_feature %in% var_names) return(NULL)

      obs_df <- tr$obs
      if (is.null(obs_df)) return(NULL)

      # Extract expression values
      expression_vals <- tryCatch({
        if (selected_feature %in% var_names) {
          idx <- which(var_names == selected_feature)
          as.numeric(adata$X[, idx])
        } else {
          rep(NA_real_, nrow(obs_df))
        }
      }, error = function(e) {
        message("[traj_data_clean] Error extracting feature: ", e$message)
        rep(NA_real_, nrow(obs_df))
      })

      # Get donor_status with proper normalization
      donor_clean <- obs_df$donor_status
      if (is.null(donor_clean)) donor_clean <- rep(NA_character_, nrow(obs_df))

      # Extract UMAP coordinates
      umap_coords <- NULL
      if (!is.null(tr$obsm) && "X_umap" %in% names(tr$obsm)) {
        umap_coords <- tr$obsm$X_umap
      }

      # Get donor ID for coloring
      donor_ids <- as.character(obs_df$imageid)

      # Build result data.frame
      islet_keys <- as.character(obs_df$base_islet_id)
      islet_keys <- gsub("^Islet_Islet_", "Islet_", islet_keys)

      result_df <- data.frame(
        case_id = as.character(obs_df$imageid),
        islet_key = islet_keys,
        combined_islet_id = as.character(obs_df$combined_islet_id),
        pt = as.numeric(obs_df$pseudotime),
        donor_status = donor_clean,
        donor_id = donor_ids,
        value = as.numeric(expression_vals),
        feature_name = selected_feature,
        stringsAsFactors = FALSE
      )

      # Add UMAP coordinates
      if (!is.null(umap_coords) && nrow(umap_coords) == nrow(result_df)) {
        result_df$umap_1 <- as.numeric(umap_coords[, 1])
        result_df$umap_2 <- as.numeric(umap_coords[, 2])
      } else {
        result_df$umap_1 <- NA_real_
        result_df$umap_2 <- NA_real_
      }

      # Merge diameter from prepared() comp table
      prep_data <- prepared()
      if (!is.null(prep_data) && !is.null(prep_data$comp)) {
        size_lookup <- prep_data$comp[c("Case ID", "Donor Status", "islet_key", "islet_diam_um")]
        result_df$`Case ID` <- sprintf("%04d", as.numeric(result_df$case_id))
        result_df$`Donor Status` <- result_df$donor_status

        # Base R merge (faster for this use case)
        merged <- merge(result_df, size_lookup,
                        by = c("Case ID", "Donor Status", "islet_key"),
                        all.x = TRUE, sort = FALSE)

        if ("islet_diam_um" %in% colnames(merged)) {
          result_df$islet_diam_um <- merged$islet_diam_um[match(
            paste(result_df$`Case ID`, result_df$`Donor Status`, result_df$islet_key),
            paste(merged$`Case ID`, merged$`Donor Status`, merged$islet_key)
          )]
        }
      }

      if (!"islet_diam_um" %in% colnames(result_df)) {
        result_df$islet_diam_um <- NA_real_
      }

      # Add spatial coordinates from segmentation data (global)
      if (!is.null(segmentation_data)) {
        result_df$centroid_x_um <- NA_real_
        result_df$centroid_y_um <- NA_real_
        result_df$has_spatial <- FALSE

        for (i in seq_len(nrow(result_df))) {
          annot <- get_islet_annotations(result_df$case_id[i], result_df$islet_key[i])
          if (!is.null(annot)) {
            result_df$centroid_x_um[i] <- annot$centroid_x
            result_df$centroid_y_um[i] <- annot$centroid_y
            result_df$has_spatial[i] <- TRUE
          }
        }

        cat("[traj_data_clean] Added spatial coordinates for",
            sum(result_df$has_spatial), "islets\n")
      }

      # Filter valid rows
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

      message("[traj_data_clean] Successfully processed ", nrow(result_df),
              " observations for ", selected_feature)
      return(result_df)
    })

    # ---------- Status bar ----------
    output$traj_status <- renderUI({
      cur <- traj_path
      if (is.null(cur) || !file.exists(cur)) cur <- resolve_traj_path()
      if (is.null(cur) || !file.exists(cur)) {
        return(tags$div(style = "color:#b00;",
          "AnnData not found. Place 'adata_ins_root.h5ad' under data/ or scripts/, or set ADATA_PATH."))
      }
      if (!requireNamespace("anndata", quietly = TRUE)) {
        return(tags$div(style = "color:#b00;",
          "R package 'anndata' not installed. Install with BiocManager::install('anndata')."))
      }
      tr <- traj_with_pseudotime()
      if (is.null(tr)) {
        return(tags$div(style = "color:#b00;",
          "AnnData object not available or H5AD failed to load."))
      }
      if (!is.null(tr$error)) {
        det <- if (!is.null(tr$error_detail) && nzchar(tr$error_detail)) paste0(" Details: ", tr$error_detail) else ""
        install_hint <- if (grepl("anndata", tr$error, ignore.case = TRUE)) " Install with: BiocManager::install('anndata')" else ""
        return(tags$div(style = "color:#b00;",
          sprintf("Trajectory load error: %s%s%s", tr$error, det, install_hint)))
      }
      ms <- traj_validate()
      jn <- traj_data_clean()
      nm <- if (is.null(ms)) 0 else nrow(ms)
      nj <- if (is.null(jn)) 0 else nrow(jn)
      ndist <- if (is.null(jn)) 0 else nrow(dplyr::distinct(jn, case_id, islet_key))
      src <- if (!is.null(traj_path) && file.exists(traj_path)) traj_path else cur

      tags$div(style = "color:#555;",
        sprintf("Loaded AnnData (%d obs). Joined %d rows (%d unique islets) of %d islet metrics. Source: %s",
                nrow(tr$obs), nj, ndist, nm, basename(src))
      )
    })

    # ---------- Trajectory scatter plot ----------
    output$traj_scatter <- renderPlotly({
      df <- traj_data_clean()
      if (is.null(df) || nrow(df) == 0) {
        return(plotly_empty() %>% layout(title = "No trajectory data available"))
      }

      # Remove outliers (>3 SD from mean) and store for reporting
      df_original <- df
      value_mean <- mean(df$value, na.rm = TRUE)
      value_sd   <- sd(df$value, na.rm = TRUE)
      outlier_threshold <- 3

      df$is_outlier <- abs(df$value - value_mean) > (outlier_threshold * value_sd)
      outliers <- df[df$is_outlier & !is.na(df$is_outlier), ]

      # Store outliers in reactiveVal
      if (nrow(outliers) > 0) {
        traj_outliers(data.frame(
          Feature      = outliers$feature_name,
          Case_ID      = outliers$case_id,
          Islet        = outliers$islet_key,
          Donor_Status = outliers$donor_status,
          Pseudotime   = round(outliers$pt, 3),
          Value        = round(outliers$value, 3),
          Z_Score      = round((outliers$value - value_mean) / value_sd, 2),
          stringsAsFactors = FALSE
        ))
      } else {
        traj_outliers(NULL)
      }

      # Remove outliers from plotting data
      df <- df[!df$is_outlier | is.na(df$is_outlier), ]

      cat(sprintf("[TRAJECTORY] Removed %d outliers (>3 SD)\n", nrow(df_original) - nrow(df)))

      # Determine color mapping and aesthetics
      color_by   <- input$traj_color_by %||% "donor_status"
      size_by    <- input$traj_point_size %||% "uniform"
      trend_type <- input$traj_show_trend %||% "by_donor"
      alpha_val  <- input$traj_alpha %||% 0.6
      point_size <- input$traj_point_size_slider %||% 3.0

      # Always apply jitter
      set.seed(42)
      jitter_amount_x <- diff(range(df$pt, na.rm = TRUE)) * 0.01
      jitter_amount_y <- diff(range(df$value, na.rm = TRUE)) * 0.02
      df$pt    <- df$pt    + runif(nrow(df), -jitter_amount_x, jitter_amount_x)
      df$value <- df$value + runif(nrow(df), -jitter_amount_y, jitter_amount_y)

      # Create base aesthetic mapping
      aes_mapping <- aes(x = pt, y = value)

      # Add color aesthetic
      if (color_by == "donor_status") {
        df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))
        aes_mapping$colour <- as.name("donor_status")
        color_title <- "Donor Status"
      } else if (color_by == "donor_id") {
        df$donor_id <- factor(df$donor_id)
        aes_mapping$colour <- as.name("donor_id")
        color_title <- "Donor ID"
      }

      # Add size aesthetic if not uniform
      if (size_by != "uniform") {
        if (size_by == "islet_diam_um" && "islet_diam_um" %in% colnames(df)) {
          if (!all(is.na(df$islet_diam_um))) {
            aes_mapping$size <- as.name("islet_diam_um")
          }
        }
      }

      # Create plot
      g <- ggplot(df, aes_mapping)
      if (size_by == "uniform") {
        g <- g + geom_point(alpha = alpha_val, size = point_size, stroke = 0)
      } else {
        g <- g + geom_point(alpha = alpha_val, stroke = 0)
      }
      g <- g +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "none",
          plot.margin = unit(c(5, 5, 5, 5), "pt"),
          plot.title = element_text(margin = margin(b = 15))
        )

      # Apply color scales
      if (color_by == "donor_status") {
        g <- g + scale_color_manual(
          values = c("ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd"),
          name = color_title, drop = FALSE)
      } else if (color_by == "donor_id") {
        df$donor_id <- factor(df$donor_id)
        n_donors <- length(levels(df$donor_id))
        if (n_donors <= 12) {
          colors <- RColorBrewer::brewer.pal(min(12, max(3, n_donors)), "Paired")
          names(colors) <- levels(df$donor_id)[1:length(colors)]
        } else {
          colors <- rainbow(n_donors, s = 0.8, v = 0.8)
          names(colors) <- levels(df$donor_id)
        }
        g <- g + scale_color_manual(values = colors, name = color_title, drop = FALSE)
      }

      # Size scale
      if (size_by != "uniform" && size_by == "islet_diam_um" && "islet_diam_um" %in% colnames(df)) {
        if (!all(is.na(df$islet_diam_um))) {
          size_range <- c(0.5 * point_size, 2.5 * point_size)
          g <- g + scale_size_continuous(name = "Islet Diameter (um)", range = size_range, guide = "none")
        }
      }

      # Trend lines
      pt_range <- range(df$pt, na.rm = TRUE)

      if (trend_type == "overall") {
        g <- g + geom_smooth(method = "loess", se = TRUE, alpha = 0.2, color = "black",
                             size = 1, fullrange = FALSE, span = 0.75)
      } else if (trend_type == "by_donor") {
        if ("donor_status" %in% names(df) && !all(is.na(df$donor_status))) {
          df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))
          trend_colors <- c("ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd")

          if (color_by == "donor_id") {
            # Add separate geom_smooth for each donor_status with explicit color
            for (status in c("ND", "Aab+", "T1D")) {
              df_subset <- df[df$donor_status == status & !is.na(df$donor_status), ]
              if (nrow(df_subset) > 0) {
                g <- g + geom_smooth(data = df_subset,
                                     aes(x = pt, y = value),
                                     method = "loess", se = TRUE, alpha = 0.15,
                                     size = 0.8, show.legend = FALSE, fullrange = FALSE, span = 0.75,
                                     color = trend_colors[status])
              }
            }
          } else {
            g <- g + geom_smooth(
              method = "loess", se = TRUE, alpha = 0.15,
              size = 0.8, show.legend = FALSE, fullrange = FALSE, span = 0.75,
              inherit.aes = FALSE,
              mapping = aes(x = pt, y = value, group = donor_status, color = donor_status)) +
              scale_color_manual(values = trend_colors,
                                name = if (color_by == "donor_status") color_title else "Trend Lines",
                                drop = FALSE,
                                guide = if (color_by == "donor_status") "legend" else "none")
          }
        }
      }

      # Labels
      selected_feature <- input$traj_feature %||% "Selected feature"
      metric_label <- "Expression Level"
      x_label <- "Pseudotime (PAGA trajectory, INS-rooted)"

      g <- g + labs(
        x = x_label,
        y = paste(selected_feature, metric_label),
        title = NULL
      ) +
      coord_cartesian(xlim = c(0, 1))

      # Convert to plotly with custom data for reliable click handling
      p <- ggplotly(g, tooltip = c("x", "y", "colour"), source = ns("traj_scatter"))

      # Add customdata to plotly traces for click handling
      if ("combined_islet_id" %in% colnames(df) && "case_id" %in% colnames(df)) {
        # Store the dataframe in a reactiveVal for coordinate-based lookup
        traj_coord_lookup(data.frame(
          pt                = df$pt,
          value             = df$value,
          case_id           = df$case_id,
          islet_key         = df$islet_key,
          combined_islet_id = df$combined_islet_id,
          donor_status      = df$donor_status,
          stringsAsFactors  = FALSE
        ))

        # Match customdata to actual trace data by coordinates
        for (i in seq_along(p$x$data)) {
          if (!is.null(p$x$data[[i]]$type) && p$x$data[[i]]$type == "scatter") {
            trace_x <- p$x$data[[i]]$x
            trace_y <- p$x$data[[i]]$y

            if (length(trace_x) > 0 && length(trace_y) > 0) {
              trace_customdata <- matrix(NA, nrow = length(trace_x), ncol = 4)

              for (j in seq_along(trace_x)) {
                tolerance <- 0.001
                matches <- which(
                  abs(df$pt - trace_x[j]) < tolerance &
                  abs(df$value - trace_y[j]) < tolerance
                )
                if (length(matches) > 0) {
                  idx <- matches[1]
                  trace_customdata[j, ] <- c(
                    as.character(df$case_id[idx]),
                    as.character(df$islet_key[idx]),
                    as.character(df$combined_islet_id[idx]),
                    as.character(df$donor_status[idx])
                  )
                }
              }

              p$x$data[[i]]$customdata <- trace_customdata
              cat("[PLOT] Added customdata to trace", i, "with", length(trace_x), "points\n")
            }
          }
        }

        cat("[PLOT] Stored coordinate lookup table with", nrow(traj_coord_lookup()), "entries\n")
      }

      # Layout and register click events
      p %>%
        layout(
          showlegend = FALSE,
          title = NULL,
          margin = list(l = 60, r = 20, t = 20, b = 50, pad = 5),
          xaxis = list(fixedrange = FALSE, automargin = TRUE, range = c(0, 1)),
          yaxis = list(fixedrange = FALSE, automargin = TRUE)
        ) %>%
        event_register("plotly_click")
    })

    # ---------- Click handler for trajectory scatter ----------
    observeEvent(event_data("plotly_click", source = ns("traj_scatter")), {
      click_data <- event_data("plotly_click", source = ns("traj_scatter"))
      if (is.null(click_data)) return()

      # Check if segmentation viewer is available
      if (!SF_AVAILABLE) {
        showNotification("Segmentation viewer requires the 'sf' package. Install with: install.packages('sf')",
                         type = "warning", duration = 5)
        return()
      }

      if (is.null(islet_spatial_lookup)) {
        showNotification("Islet spatial lookup data not available", type = "warning", duration = 5)
        return()
      }

      cat("[SEGMENTATION CLICK] Received click event\n")

      click_x <- click_data$x
      click_y <- click_data$y

      # First try to get from customdata
      custom    <- click_data$customdata
      case_id   <- NULL
      islet_key <- NULL

      if (!is.null(custom) && length(custom) >= 2) {
        case_id   <- custom[[1]]
        islet_key <- custom[[2]]
        cat("[SEGMENTATION CLICK] From customdata: case_id=", case_id, ", islet_key=", islet_key, "\n")
      }

      # If customdata not available, look up from coordinate table
      if (is.null(case_id) || is.null(islet_key) || is.na(case_id) || is.na(islet_key)) {
        lookup <- traj_coord_lookup()
        if (!is.null(lookup) && nrow(lookup) > 0) {
          tolerance <- 0.01
          matches <- which(
            abs(lookup$pt - click_x) < tolerance &
            abs(lookup$value - click_y) < tolerance
          )
          if (length(matches) > 0) {
            idx <- matches[1]
            case_id   <- lookup$case_id[idx]
            islet_key <- lookup$islet_key[idx]
            cat("[SEGMENTATION CLICK] From coord lookup: case_id=", case_id, ", islet_key=", islet_key, "\n")
          }
        }
      }

      if (is.null(case_id) || is.null(islet_key) || is.na(case_id) || is.na(islet_key)) {
        showNotification("Could not identify clicked islet", type = "warning", duration = 3)
        return()
      }

      # Look up centroid coordinates from spatial lookup
      spatial_match <- islet_spatial_lookup[
        islet_spatial_lookup$case_id == case_id & islet_spatial_lookup$islet_key == islet_key,
      ]

      # Try zero-padded variant (e.g. "112" -> "0112")
      if (nrow(spatial_match) == 0) {
        case_id_num <- suppressWarnings(as.integer(case_id))
        if (!is.na(case_id_num)) {
          case_id_padded <- sprintf("%04d", case_id_num)
          spatial_match <- islet_spatial_lookup[
            islet_spatial_lookup$case_id == case_id_padded & islet_spatial_lookup$islet_key == islet_key,
          ]
          if (nrow(spatial_match) > 0) case_id <- case_id_padded
        }
      }

      if (nrow(spatial_match) == 0) {
        showNotification(paste("Spatial data not found for", islet_key, "in case", case_id),
                         type = "warning", duration = 3)
        return()
      }

      centroid_x <- spatial_match$centroid_x_um[1]
      centroid_y <- spatial_match$centroid_y_um[1]

      cat("[SEGMENTATION CLICK] Centroid: x=", centroid_x, ", y=", centroid_y, "\n")

      # Store selected islet info - triggers embedded viewer to update
      selected_islet(list(
        case_id    = case_id,
        islet_key  = islet_key,
        centroid_x = centroid_x,
        centroid_y = centroid_y
      ))
    })

    # ---------- Segmentation viewer panel ----------
    output$segmentation_viewer_panel <- renderUI({
      info <- selected_islet()
      if (is.null(info)) return(NULL)

      # Check if single-cell drill-down data exists for this islet
      has_cells <- drilldown_available() &&
        !is.null(load_islet_cells(info$case_id, info$islet_key))

      # Marker choices for color-by dropdown
      marker_choices <- c("phenotype" = "phenotype")
      if (has_cells) {
        cells <- load_islet_cells(info$case_id, info$islet_key)
        marker_cols <- setdiff(colnames(cells),
                               c("X_centroid", "Y_centroid", "phenotype", "cell_region",
                                 "Cell.Area", "Cell Area", "Nucleus.Area", "Nucleus Area"))
        if (length(marker_cols) > 0) {
          marker_choices <- c(marker_choices, setNames(marker_cols, marker_cols))
        }
      }

      div(class = "card", style = "padding: 12px; margin-bottom: 20px; border: 2px solid #0066CC;",
        # Header row: title + view mode toggle + close button
        div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px;",
          h5(paste("Islet:", info$islet_key, "(Case", info$case_id, ")"),
             style = "margin: 0; color: #0066CC; font-size: 14px;"),
          div(
            if (has_cells) {
              tags$div(style = "display: inline-block; margin-right: 8px;",
                radioButtons("drilldown_view_mode", NULL,
                             choices = c("Boundaries", "Single Cells"),
                             selected = "Single Cells", inline = TRUE)
              )
            },
            actionButton(ns("clear_segmentation"), "Close", class = "btn btn-sm btn-outline-secondary")
          )
        ),
        # Controls row (visible only in Single Cells mode)
        if (has_cells) {
          conditionalPanel(
            condition = "input.drilldown_view_mode == 'Single Cells'",
            div(style = "margin-bottom: 8px;",
              fluidRow(
                column(6,
                  selectInput("drilldown_color_by", "Color by",
                              choices = marker_choices, selected = "phenotype")
                ),
                column(6,
                  checkboxInput("drilldown_show_peri", "Show peri-islet cells", value = TRUE)
                )
              )
            )
          )
        },
        # Main content: plot + sidebar stacked vertically to fit col-5
        conditionalPanel(
          condition = if (has_cells) "input.drilldown_view_mode == 'Boundaries'" else "true",
          plotOutput("islet_segmentation_view", height = "350px")
        ),
        if (has_cells) {
          conditionalPanel(
            condition = "input.drilldown_view_mode == 'Single Cells'",
            plotOutput("islet_drilldown_view", height = "350px"),
            div(style = "padding: 8px; background-color: #f8f9fa; border-radius: 5px; margin-top: 8px;",
              h6("Cell Composition", style = "margin-top: 0;"),
              plotOutput("islet_drilldown_summary", height = "180px"),
              tableOutput("islet_drilldown_table")
            )
          )
        },
        # Compact legend (shown in Boundaries mode)
        conditionalPanel(
          condition = if (has_cells) "input.drilldown_view_mode == 'Boundaries'" else "true",
          div(style = "padding: 8px; background-color: #f8f9fa; border-radius: 5px; margin-top: 8px; font-size: 12px;",
            div(style = "display: flex; flex-wrap: wrap; gap: 10px; margin-bottom: 6px;",
              div(style = "display: flex; align-items: center;",
                div(style = "width: 20px; height: 3px; background-color: #0066CC; margin-right: 5px;"),
                span("Islet")
              ),
              div(style = "display: flex; align-items: center;",
                div(style = "width: 20px; height: 3px; background-color: #00CCCC; margin-right: 5px;"),
                span("+20\u00b5m")
              ),
              div(style = "display: flex; align-items: center;",
                div(style = "width: 20px; height: 3px; background-color: #CC00CC; margin-right: 5px;"),
                span("Nerve")
              ),
              div(style = "display: flex; align-items: center;",
                div(style = "width: 20px; height: 3px; background-color: #CC0000; margin-right: 5px;"),
                span("Capillary")
              ),
              div(style = "display: flex; align-items: center;",
                div(style = "width: 20px; height: 3px; background-color: #00AA00; margin-right: 5px;"),
                span("Lymphatic")
              ),
              div(style = "display: flex; align-items: center;",
                div(style = "width: 20px; height: 3px; background-color: #FFD700; margin-right: 5px;"),
                tags$strong("Selected")
              )
            ),
            span(paste("Centroid:", round(info$centroid_x, 0), "/", round(info$centroid_y, 0), "\u00b5m"),
                 style = "color: #666;")
          )
        )
      )
    })

    # ---------- Clear segmentation viewer ----------
    observeEvent(input$clear_segmentation, {
      selected_islet(NULL)
    })

    # Note: islet_segmentation_view renderPlot is at root level in app.R
    # (shared between Plot modal and Trajectory embedded panel)

    # ---------- Outlier table ----------
    output$traj_outlier_info <- renderUI({
      show_table <- isTRUE(input$show_outlier_table)
      outlier_data <- traj_outliers()

      if (show_table && !is.null(outlier_data) && nrow(outlier_data) > 0) {
        tagList(
          div(style = "background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 5px; padding: 10px; margin-top: 10px;",
            h5(style = "color: #856404; margin-top: 0;",
               sprintf("\u26a0\ufe0f %d Outlier%s Removed (>3 SD)",
                       nrow(outlier_data), ifelse(nrow(outlier_data) > 1, "s", ""))),
            p(style = "color: #856404; font-size: 12px; margin-bottom: 10px;",
              "The following data points were excluded from the plot because they exceed 3 standard deviations from the mean:"),
            div(style = "max-height: 200px; overflow-y: auto;",
              renderTable({
                outlier_data
              }, striped = TRUE, hover = TRUE, bordered = TRUE, spacing = "xs", width = "100%")
            )
          )
        )
      } else {
        NULL
      }
    })

    # ---------- Legend for trajectory scatter ----------
    output$traj_legend <- renderUI({
      color_by <- input$traj_color_by %||% "donor_status"

      if (color_by == "donor_status") {
        tagList(
          div(style = "display: flex; align-items: center; margin-bottom: 8px;",
            div(style = "width: 20px; height: 20px; background-color: #2ca02c; border-radius: 3px; margin-right: 8px;"),
            span("ND", style = "font-size: 13px;")
          ),
          div(style = "display: flex; align-items: center; margin-bottom: 8px;",
            div(style = "width: 20px; height: 20px; background-color: #ffcc00; border-radius: 3px; margin-right: 8px;"),
            span("Aab+", style = "font-size: 13px;")
          ),
          div(style = "display: flex; align-items: center;",
            div(style = "width: 20px; height: 20px; background-color: #9467bd; border-radius: 3px; margin-right: 8px;"),
            span("T1D", style = "font-size: 13px;")
          )
        )
      } else if (color_by == "donor_id") {
        df <- traj_data_clean()
        if (!is.null(df) && nrow(df) > 0) {
          donor_ids <- sort(unique(df$donor_id))
          n_donors <- length(donor_ids)

          # Generate colors matching the plot
          if (n_donors <= 12) {
            colors <- RColorBrewer::brewer.pal(min(12, max(3, n_donors)), "Paired")
          } else {
            colors <- rainbow(n_donors, s = 0.8, v = 0.8)
          }

          # Split into two columns
          mid_point <- ceiling(n_donors / 2)
          col1_ids <- donor_ids[1:mid_point]
          col2_ids <- if (n_donors > mid_point) donor_ids[(mid_point + 1):n_donors] else NULL

          col1_items <- lapply(seq_along(col1_ids), function(i) {
            div(style = "display: flex; align-items: center; margin-bottom: 6px;",
              div(style = sprintf("width: 18px; height: 18px; background-color: %s; border-radius: 3px; margin-right: 6px;", colors[i])),
              span(col1_ids[i], style = "font-size: 12px;")
            )
          })

          col2_items <- if (!is.null(col2_ids)) {
            lapply(seq_along(col2_ids), function(i) {
              idx <- mid_point + i
              div(style = "display: flex; align-items: center; margin-bottom: 6px;",
                div(style = sprintf("width: 18px; height: 18px; background-color: %s; border-radius: 3px; margin-right: 6px;", colors[idx])),
                span(col2_ids[i], style = "font-size: 12px;")
              )
            })
          } else {
            NULL
          }

          tagList(
            div(style = "display: flex; gap: 10px;",
              div(style = "flex: 1;", col1_items),
              if (!is.null(col2_items)) div(style = "flex: 1;", col2_items)
            )
          )
        } else {
          p("No data available", style = "font-size: 12px; color: #999;")
        }
      }
    })

    # ---------- UMAP: Donor Status ----------
    output$traj_umap_donor <- renderPlot({
      df <- traj_data_clean()
      if (is.null(df) || nrow(df) == 0) return(NULL)

      if (all(is.na(df$umap_1)) || all(is.na(df$umap_2))) {
        return(ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "UMAP coordinates not available", size = 5) +
          theme_void())
      }

      df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))

      ggplot(df, aes(x = umap_1, y = umap_2, color = donor_status)) +
        geom_point(alpha = 0.6, size = 3.0) +
        scale_color_manual(values = c("ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd")) +
        scale_x_continuous(expand = expansion(mult = 0.02)) +
        scale_y_continuous(expand = expansion(mult = 0.02)) +
        labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP: Donor Status") +
        theme_minimal(base_size = 12)
    })

    # ---------- UMAP: Selected Feature ----------
    output$traj_umap_feature <- renderPlot({
      df <- traj_data_clean()
      if (is.null(df) || nrow(df) == 0) return(NULL)

      if (all(is.na(df$umap_1)) || all(is.na(df$umap_2))) {
        return(ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "UMAP coordinates not available", size = 5) +
          theme_void())
      }

      selected_feature <- input$traj_feature %||% "Selected feature"

      g <- ggplot(df, aes(x = umap_1, y = umap_2, color = value)) +
        geom_point(alpha = 0.6, size = 3.0) +
        scale_x_continuous(expand = expansion(mult = 0.02)) +
        scale_y_continuous(expand = expansion(mult = 0.02)) +
        labs(x = "UMAP 1", y = "UMAP 2",
             title = paste("UMAP:", selected_feature),
             color = "Expression") +
        theme_minimal(base_size = 12)

      # Continuous colormap scaled to data min/max
      val_range <- range(df$value, na.rm = TRUE)
      g <- g + scale_color_viridis_c(option = "inferno", limits = val_range,
                                      na.value = "#bbbbbb")

      g
    })

    # ---------- Heatmap: Donor Status Progression ----------
    output$traj_heatmap <- renderPlot({
      df <- traj_data_clean()
      if (is.null(df) || nrow(df) == 0) return(NULL)

      # Create binned analysis
      enc <- function(s) ifelse(s == "ND", 0, ifelse(s == "Aab+", 1, ifelse(s == "T1D", 2, NA_real_)))
      df$ds_code <- enc(df$donor_status)

      # Normalize pseudotime to 0-1 range for binning
      pt_range <- range(df$pt, na.rm = TRUE)
      df$pt_norm <- (df$pt - pt_range[1]) / (pt_range[2] - pt_range[1])

      nb   <- 25
      brks <- seq(0, 1, length.out = nb + 1)
      df$pt_bin <- cut(pmax(0, pmin(1, df$pt_norm)), breaks = brks, include.lowest = TRUE, right = FALSE)

      # Calculate averages per bin using base R
      bin_levels <- levels(df$pt_bin)
      hm_list <- list()

      for (i in seq_along(bin_levels)) {
        bin_name <- bin_levels[i]
        bin_data <- df[!is.na(df$pt_bin) & df$pt_bin == bin_name, ]

        if (nrow(bin_data) >= 3) {
          hm_list[[length(hm_list) + 1]] <- data.frame(
            pt_bin            = factor(bin_name, levels = bin_levels),
            avg_donor_status  = mean(bin_data$ds_code, na.rm = TRUE),
            count             = nrow(bin_data),
            stringsAsFactors  = FALSE
          )
        }
      }

      if (length(hm_list) == 0) return(NULL)
      hm <- do.call(rbind, hm_list)
      if (nrow(hm) == 0) return(NULL)

      # Calculate bin midpoints in actual pseudotime range
      all_mids <- head(brks, -1) + diff(brks)[1] / 2
      bin_idx  <- match(as.character(hm$pt_bin), bin_levels)
      hm$x     <- all_mids[bin_idx] * (pt_range[2] - pt_range[1]) + pt_range[1]

      ggplot(hm, aes(x = x, y = 1, fill = avg_donor_status)) +
        geom_tile(height = 1) +
        scale_fill_gradientn(
          colors   = c("#2ca02c", "#ffcc00", "#9467bd"),
          limits   = c(0, 2),
          na.value = "#dddddd",
          name     = "Average\nDonor Type",
          breaks   = c(0, 1, 2),
          labels   = c("ND", "Aab+", "T1D")
        ) +
        scale_x_continuous(limits = c(pt_range[1], pt_range[2]), expand = c(0, 0)) +
        labs(x = "Pseudotime", y = "", title = "Donor Status Progression Along Pseudotime") +
        theme_minimal() +
        theme(
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid   = element_blank(),
          plot.title   = element_text(size = 12, hjust = 0.5)
        )
    })

    # ---------- Multi-Feature Heatmap ----------

    # Marker selector UI
    output$multi_heatmap_marker_selector <- renderUI({
      tr <- traj()
      if (is.null(tr) || !is.null(tr$error) || is.null(tr$var_names)) return(NULL)
      checkboxGroupInput(ns("multi_heatmap_markers"), "Select markers:",
                         choices = tr$var_names, selected = tr$var_names, inline = TRUE)
    })

    # Multi-feature data reactive
    traj_multi_feature_data <- reactive({
      markers <- input$multi_heatmap_markers
      if (is.null(markers) || length(markers) == 0) return(NULL)
      tr <- traj_with_pseudotime()
      if (is.null(tr) || !is.null(tr$error) || is.null(tr$adata)) return(NULL)
      valid <- intersect(markers, tr$var_names)
      if (length(valid) == 0) return(NULL)
      idx <- match(valid, tr$var_names)
      expr <- tryCatch(as.matrix(tr$adata$X[, idx, drop = FALSE]), error = function(e) NULL)
      if (is.null(expr)) return(NULL)
      colnames(expr) <- valid
      result <- data.frame(pt = as.numeric(tr$obs$pseudotime), stringsAsFactors = FALSE)
      result <- cbind(result, as.data.frame(expr))
      result <- result[is.finite(result$pt), , drop = FALSE]
      list(data = result, markers = valid)
    })

    # Multi-row heatmap renderPlot
    output$traj_multi_heatmap <- renderPlot({
      mfd <- traj_multi_feature_data()
      if (is.null(mfd)) return(NULL)

      df <- mfd$data
      markers <- mfd$markers
      nb <- input$multi_heatmap_nbins %||% 20

      # Bin pseudotime
      pt_range <- range(df$pt, na.rm = TRUE)
      df$pt_norm <- (df$pt - pt_range[1]) / (pt_range[2] - pt_range[1])
      brks <- seq(0, 1, length.out = nb + 1)
      df$pt_bin <- cut(pmax(0, pmin(1, df$pt_norm)), breaks = brks, include.lowest = TRUE, right = FALSE)
      bin_levels <- levels(df$pt_bin)

      # Compute per-bin mean for each marker, then z-score
      all_mids <- head(brks, -1) + diff(brks)[1] / 2
      hm_rows <- list()

      for (mk in markers) {
        bin_means <- numeric(length(bin_levels))
        bin_counts <- integer(length(bin_levels))
        for (i in seq_along(bin_levels)) {
          bin_data <- df[!is.na(df$pt_bin) & df$pt_bin == bin_levels[i], mk]
          bin_counts[i] <- length(bin_data)
          bin_means[i] <- if (length(bin_data) >= 3) mean(bin_data, na.rm = TRUE) else NA_real_
        }
        # Z-score across bins
        valid_means <- bin_means[is.finite(bin_means)]
        if (length(valid_means) < 2) next
        mu <- mean(valid_means)
        sd_val <- sd(valid_means)
        if (!is.finite(sd_val) || sd_val == 0) sd_val <- 1
        z <- (bin_means - mu) / sd_val
        z <- pmax(-2.5, pmin(2.5, z))  # clamp

        for (i in seq_along(bin_levels)) {
          if (is.finite(z[i])) {
            hm_rows[[length(hm_rows) + 1]] <- data.frame(
              marker = mk,
              x = all_mids[i] * (pt_range[2] - pt_range[1]) + pt_range[1],
              z = z[i],
              stringsAsFactors = FALSE
            )
          }
        }
      }

      if (length(hm_rows) == 0) return(NULL)
      hm <- do.call(rbind, hm_rows)

      # Order markers: hormones first, then immune, then others
      hormone_order <- intersect(c("INS", "GCG", "SST"), markers)
      immune_order <- intersect(c("CD3e", "CD4", "CD8a", "CD20", "CD45", "CD68", "CD163", "HLADR"), markers)
      other_order <- setdiff(markers, c(hormone_order, immune_order))
      marker_order <- c(hormone_order, immune_order, other_order)
      hm$marker <- factor(hm$marker, levels = rev(marker_order))

      n_markers <- length(marker_order)

      ggplot(hm, aes(x = x, y = marker, fill = z)) +
        geom_tile() +
        scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                             midpoint = 0, limits = c(-2.5, 2.5),
                             name = "Z-score") +
        scale_x_continuous(limits = c(pt_range[1], pt_range[2]), expand = c(0, 0)) +
        labs(x = "Pseudotime", y = "", title = "Multi-Feature Expression Along Pseudotime") +
        theme_minimal() +
        theme(
          panel.grid = element_blank(),
          plot.title = element_text(size = 12, hjust = 0.5),
          axis.text.y = element_text(size = 11)
        )
    }, height = function() {
      n <- length(input$multi_heatmap_markers %||% character(0))
      max(150, 40 + n * 30)
    })

  })
}
