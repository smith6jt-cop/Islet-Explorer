islet_boundary_notice <- function(case_id = "all") {
  path <- data_path("trajectory", "boundary_review.csv")
  if (!file.exists(path)) return(NULL)
  review <- read.csv(path, colClasses = c(case_id = "character"))
  review$approved <- tolower(as.character(review$approved)) %in% c("true", "t", "1")
  if (identical(case_id, "all")) return(div(class = "alert alert-info", role = "status",
    strong("Boundary data"),
    p(sprintf("%d donors available. The most recent delivered boundaries are used for display and analysis.", nrow(review)))))
  row <- review[review$case_id == case_id, , drop = FALSE]
  if (!nrow(row)) return(NULL)
  div(class = "alert alert-info", role = "status",
      strong(paste("Donor", case_id, "— current delivered boundaries")),
      p(sprintf("%s islet outlines available.%s", format(row$islets[1], big.mark = ","),
        if (row$pending_new_candidates[1] > 0) paste0(" ", format(row$pending_new_candidates[1], big.mark = ","), " regenerated candidates remain in the QuPath review project.") else "")))
}
