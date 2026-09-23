# The schema-4 release owns its audited identities and states. Older releases
# retain their original retirement policy and palettes.
.typing_catalog_cache <- new.env(parent = emptyenv())
islet_typing_catalog <- function() {
  if (!identical(islet_release_schema(), "4")) return(NULL)
  path <- data_path("typing", "catalog.json")
  if (!file.exists(path)) stop("Schema 4 typing catalog is missing")
  key <- normalizePath(path, mustWork = TRUE)
  if (is.null(.typing_catalog_cache[[key]])) {
    catalog <- jsonlite::read_json(path, simplifyVector = FALSE)
    labels <- unlist(lapply(c(catalog$base_types, catalog$states), function(x) if (is.character(x)) x else x$name), use.names = FALSE)
    if (!length(labels) || anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) stop("Invalid typing catalog labels")
    .typing_catalog_cache[[key]] <- catalog
  }
  .typing_catalog_cache[[key]]
}
islet_catalog_names <- function(entries) {
  unlist(lapply(entries, function(x) if (is.character(x)) x else x$name), use.names = FALSE)
}

# Runtime display vocabulary is deliberately broader than the reviewed typing
# catalog. The catalog controls scientific identity rules; it must not hide
# measured cell types from exploratory Plot/Spatial selectors.
islet_app_phenotype_names <- function(extra = character(0)) {
  palette_names <- if (exists("PHENOTYPE_COLORS", inherits = TRUE)) names(PHENOTYPE_COLORS) else character(0)
  catalog <- tryCatch(islet_typing_catalog(), error = function(e) NULL)
  catalog_names <- if (is.null(catalog)) character(0) else
    islet_catalog_names(c(catalog$base_types, catalog$states, catalog$display_types %||% list()))
  unique(c(palette_names, catalog_names, as.character(extra)))
}

islet_current_phenotype_names <- function(extra = character(0)) {
  catalog <- tryCatch(islet_typing_catalog(), error = function(e) NULL)
  catalog_names <- if (is.null(catalog)) character(0) else
    islet_catalog_names(c(catalog$base_types, catalog$states))
  unique(c(catalog_names, as.character(extra)))
}
islet_phenotype_choices <- function(include_states = TRUE) {
  catalog <- islet_typing_catalog()
  if (is.null(catalog)) return(phenotype_keep(names(PHENOTYPE_COLORS)))
  base <- islet_catalog_names(catalog$base_types)
  out <- list("Base identities" = as.list(setNames(base, base)))
  if (include_states && length(catalog$states)) {
    states <- islet_catalog_names(catalog$states)
    out[["Finer phenotype states"]] <- as.list(setNames(states, states))
  }
  out
}
islet_endocrine_population <- function() {
  catalog <- islet_typing_catalog()
  if (is.null(catalog)) "Endocrine union" else catalog$endocrine_population %||% "Endocrine union"
}
islet_cell_color_choices <- function(markers) {
  catalog <- islet_typing_catalog()
  if (is.null(catalog)) return(c("Base identity" = "phenotype", markers))
  out <- list("Identity" = list("Base identity" = "phenotype"))
  states <- islet_catalog_names(catalog$states)
  if (length(states)) out[["Finer phenotype states"]] <- as.list(setNames(paste0("state::", states), states))
  out[["Marker intensity"]] <- as.list(setNames(markers, markers))
  out
}
islet_catalog_palette <- function(palette) {
  catalog <- islet_typing_catalog()
  if (is.null(catalog)) return(palette)
  entries <- c(catalog$base_types, catalog$states)
  labels <- islet_catalog_names(entries)
  generated <- setNames(grDevices::hcl.colors(length(labels), "Dynamic"), labels)
  for (entry in entries) {
    if (is.list(entry) && !is.null(entry$color)) generated[[entry$name]] <- entry$color
  }
  present <- intersect(names(palette), labels)
  generated[present] <- palette[present]
  generated
}
