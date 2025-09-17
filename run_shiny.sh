#!/usr/bin/env bash
set -euo pipefail

Rscript -e "shiny::runApp('app/shiny_app', host = '0.0.0.0', port = as.integer(Sys.getenv('PORT', unset = '3838')))"
