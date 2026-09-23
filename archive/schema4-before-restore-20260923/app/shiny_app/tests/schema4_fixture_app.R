#!/usr/bin/env Rscript
# Isolated browser fixture using the actual trajectory UI, server and backend.
# No app services, scientific data or production manifests are modified.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: schema4_fixture_app.R REPOSITORY DATA_ROOT PORT")
repo <- normalizePath(args[1]); data_root <- normalizePath(args[2]); port <- as.integer(args[3])
stopifnot(isTRUE(jsonlite::read_json(file.path(data_root, "candidate.json"))$synthetic))
Sys.setenv(ISLET_DATA_SCHEMA = "4", ISLET_RUNTIME_MODE = "development")
suppressPackageStartupMessages({library(shiny); library(ggplot2); library(plotly); library(dplyr); library(tidyr)})
`%||%` <- function(a,b) if (is.null(a)) b else a
data_path <- function(...) file.path(data_root,...)
PHENOTYPE_RETIRED <- character()
PHENOTYPE_COLORS <- c("Beta cell" = "#237ca1")
DONOR_COLORS_BRIGHT <- c("ND"="#3377BB", "Aab+"="#EEAA33", "T1D"="#BB3355")
phenotype_keep <- identity
islet_boundary_notice <- function(...) NULL
islet_safe_path_segment <- function(x) grepl("^[a-zA-Z0-9_.-]+$", x)
for (file in c("01_runtime_data_backend.R","03_typing_catalog.R","trajectory_marker_helpers.R","mod_trajectory_ui.R","mod_trajectory_server.R")) source(file.path(repo,"app/shiny_app/R",file))
ui <- fluidPage(title="Schema 4 synthetic regression fixture", tabsetPanel(id="tabs",
  tabPanel("Trajectory", trajectory_ui("traj")), tabPanel("Viewer",p("Synthetic fixture; no scientific image data."))))
server <- function(input,output,session) {
  selected <- reactiveVal(NULL);forced <- reactiveVal(NULL)
  trajectory_server("traj",reactive(list()),selected,forced,reactive(input$tabs))
}
shiny::runApp(shinyApp(ui,server),host="127.0.0.1",port=port,launch.browser=FALSE)
