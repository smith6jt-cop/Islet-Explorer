pkgs <- c(
  "shiny","readxl","dplyr","stringr","tidyr","ggplot2","plotly","broom","base64enc"
)
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
message("Installed/verified: ", paste(pkgs, collapse=", "))

# Install vitessceR (from CRAN if available; fallback to GitHub)
install_vitesscer <- function() {
  if (requireNamespace("vitessceR", quietly = TRUE)) {
    message("vitessceR already installed")
    return(invisible(TRUE))
  }
  ok <- FALSE
  # Try CRAN first
  try({
    install.packages("vitessceR", repos = "https://cloud.r-project.org")
    ok <- requireNamespace("vitessceR", quietly = TRUE)
  }, silent = TRUE)
  if (!ok) {
    # Fallback to GitHub
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    try({
      remotes::install_github("vitessce/vitessceR")
      ok <- requireNamespace("vitessceR", quietly = TRUE)
    }, silent = TRUE)
  }
  if (ok) message("Installed/verified: vitessceR") else warning("vitessceR install did not complete; see logs")
  invisible(ok)
}

install_vitesscer()
