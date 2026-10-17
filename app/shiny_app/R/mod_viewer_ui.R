viewer_ui <- function(id) {
  ns <- NS(id)
  div(style = "width: 100%;",
    tab_help("Viewer help", "show_phenotype_rules_viewer",
      tags$p(tags$strong("Images."),
             " One whole-slide CODEX image per donor, shown in the Avivator viewer at full resolution. Choose a donor",
             " from Select Image; the file name gives donor status and case number."),
      tags$p(tags$strong("Navigating."),
             " Click inside the viewer first to enable pan, zoom and mouse-wheel zoom. Avivator's own panel sets the",
             " channels, their colours and contrast.")),
    # In-place iframe src updater. The server pushes the new Avivator URL here
    # instead of letting renderUI recreate the <iframe>. Recreating the iframe
    # (a WebGL canvas) on every reactive tick causes heavy flicker under
    # Firefox WebRender on macOS (Mozilla bug #1555544). We update the src of
    # the existing node and dedup against the last value via a data-relsrc
    # attribute so spurious messages don't trigger a reload.
    tags$script(HTML("
Shiny.addCustomMessageHandler('islet_update_viewer_src', function(msg){
  if (!msg || !msg.id) return;
  var f = document.getElementById(msg.id);
  if (!f) return;
  if (f.getAttribute('data-relsrc') === msg.src) return;
  f.setAttribute('data-relsrc', msg.src);
  f.src = msg.src;
});
")),
    control_panel("Data controls", open = TRUE,
      uiOutput(ns("image_selector"))),
    uiOutput(ns("viewer_frame"))
  )
}
