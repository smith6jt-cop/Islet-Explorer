viewer_ui <- function(id) {
  ns <- NS(id)
  div(style = "width: 100%;",
    uiOutput(ns("local_image_picker"))
  )
}
