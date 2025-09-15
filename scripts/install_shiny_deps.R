pkgs <- c(
  "shiny","readxl","dplyr","stringr","tidyr","ggplot2","plotly","broom"
)
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
message("Installed/verified: ", paste(pkgs, collapse=", "))

