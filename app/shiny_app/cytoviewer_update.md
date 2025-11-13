# Cytoviewer Integration Update

## Summary

This document contains the code updates needed to replace AVIVATOR with cytomapper in the Viewer tab.

## 1. Add CytoImageList Loading Functions

Insert these functions **before the server function** (around line 2211):

```r
# ===== CytoImageList Loading Functions =====

# Load OME-TIFF images as CytoImageList for cytomapper
load_cytoimagelist <- function(image_files, names = NULL) {
  if (!requireNamespace("cytomapper", quietly = TRUE)) {
    warning("cytomapper package not available")
    return(NULL)
  }

  if (length(image_files) == 0) {
    warning("No image files provided")
    return(NULL)
  }

  # Ensure files exist
  existing_files <- image_files[file.exists(image_files)]
  if (length(existing_files) == 0) {
    warning("None of the provided image files exist")
    return(NULL)
  }

  tryCatch({
    # Load images using cytomapper
    image_list <- cytomapper::loadImages(existing_files, pattern = NULL)

    # Set names if provided
    if (!is.null(names) && length(names) == length(existing_files)) {
      names(image_list) <- names
    } else {
      # Use case IDs from filenames
      names(image_list) <- gsub(".*_([0-9]{4})\\\\..*", "\\\\1", basename(existing_files))
    }

    cat(sprintf("[CYTOMAPPER] Loaded %d image(s)\\n", length(image_list)))
    image_list
  }, error = function(e) {
    warning(sprintf("Error loading images: %s", e$message))
    NULL
  })
}

# Get available cases from local_images directory
get_available_cases <- function(image_dir = file.path("www", "local_images")) {
  if (!dir.exists(image_dir)) {
    return(character(0))
  }

  image_files <- list.files(image_dir, pattern = "\\\\.ome\\\\.tiff?$",
                            ignore.case = TRUE, full.names = FALSE)

  # Extract case IDs from filenames
  case_ids <- gsub(".*_([0-9]{4})\\\\..*", "\\\\1", image_files)
  case_ids <- unique(case_ids[grepl("^[0-9]{4}$", case_ids)])

  sort(case_ids)
}

# Find image file for a specific case
find_image_file <- function(case_id, image_dir = file.path("www", "local_images")) {
  if (!dir.exists(image_dir)) {
    return(NULL)
  }

  pattern <- sprintf(".*%s.*\\\\.ome\\\\.tiff?$", case_id)
  image_files <- list.files(image_dir, pattern = pattern,
                            ignore.case = TRUE, full.names = TRUE)

  if (length(image_files) > 0) {
    return(image_files[1])
  }
  NULL
}
```

## 2. Replace Viewer Tab Outputs

Replace the `output$local_image_picker` section (starting around line 4811) with:

```r
# ===== Viewer Tab Server Logic =====

# Populate case selector for Viewer tab
observe({
  cases <- get_available_cases()
  updateSelectInput(session, "viewer_case",
                   choices = c("Select a case..." = "", cases),
                   selected = if(length(cases) > 0) cases[1] else "")
})

# Image picker and channel selector for Viewer tab
output$local_image_picker <- renderUI({
  if (is.null(input$tabs) || input$tabs != "Viewer") {
    return(NULL)
  }

  cases <- get_available_cases()

  if (length(cases) == 0) {
    return(div(
      style = "padding: 20px; background: #ffe0e0; border-left: 4px solid #ff6b6b; border-radius: 5px;",
      h4("No Images Found"),
      p("No OME-TIFF images found in www/local_images/")
    ))
  }

  # Default channels to display
  default_channels <- c("DAPI", "INS", "GCG", "SST", "CD31")
  available_defaults <- intersect(default_channels, channel_names)
  if (length(available_defaults) == 0 && length(channel_names) > 0) {
    available_defaults <- head(channel_names, 3)
  }

  wellPanel(
    h4("Image Selection"),
    selectInput("viewer_case", "Select Donor Case:",
                choices = c("Select a case..." = "", cases),
                selected = if(length(cases) > 0) cases[1] else ""),
    hr(),
    h5("Select Channels to Display"),
    p(style = "font-size: 90%; color: #666;", "(Maximum 6 channels)"),
    checkboxGroupInput(
      "viewer_channels",
      NULL,
      choices = channel_names,
      selected = available_defaults
    ),
    hr(),
    actionButton("refresh_viewer", "Refresh Image", class = "btn-primary btn-sm")
  )
})

# Render cytoviewer plot
output$vit_view <- renderPlot({
  req(input$viewer_case, input$viewer_channels)

  image_file <- find_image_file(input$viewer_case)
  if (is.null(image_file)) {
    plot.new()
    text(0.5, 0.5, sprintf("Image not found for case %s", input$viewer_case),
         cex = 1.5, col = "red")
    return()
  }

  # Load as CytoImageList
  img_list <- load_cytoimagelist(image_file, names = input$viewer_case)

  if (is.null(img_list)) {
    plot.new()
    text(0.5, 0.5, "Error loading image",
         cex = 1.5, col = "red")
    return()
  }

  # Check if cytomapper is available
  if (!requireNamespace("cytomapper", quietly = TRUE)) {
    plot.new()
    text(0.5, 0.5, "cytomapper package not installed\\nRun: BiocManager::install('cytomapper')",
         cex = 1.2, col = "red")
    return()
  }

  # Get selected channels
  selected_chans <- input$viewer_channels
  if (is.null(selected_chans) || length(selected_chans) == 0) {
    plot.new()
    text(0.5, 0.5, "Please select at least one channel to display",
         cex = 1.2, col = "gray")
    return()
  }

  # Limit to 6 channels for display
  if (length(selected_chans) > 6) {
    selected_chans <- selected_chans[1:6]
  }

  tryCatch({
    # Use plotPixels to display the image
    cytomapper::plotPixels(
      image = img_list,
      colour_by = selected_chans,
      bcg = lapply(selected_chans, function(x) c(0, 1, 1.5)),
      legend = list(colour_by.title.cex = 0.7,
                   colour_by.labels.cex = 0.6),
      return_plot = FALSE
    )
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, sprintf("Error displaying image:\\n%s", e$message),
         cex = 1, col = "red")
  })
}, height = 700, width = 1200)
```

## Implementation

Use the bash script below to apply these changes:

```bash
# Step 1: Add CytoImageList functions before server (around line 2211)
# Step 2: Replace output$local_image_picker and output$vit_view
# Step 3: Test the app
```
