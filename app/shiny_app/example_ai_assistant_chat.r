readRenviron("~/.Renviron")
Sys.getenv("KEY")  # should show non-empty (masked in consoles)

install.packages(c("shiny", "jsonlite", "httr2"))

library(shiny)
library(httr2)
library(jsonlite)

API_KEY <- Sys.getenv("KEY")
API_BASE <- Sys.getenv("BASE", "https://api.ai.it.ufl.edu/v1")
endpoint <- function(path) {
  base <- sub("/+\\z", "", API_BASE)
  paste0(base, "/", path)
}
req_base <- function() {
  request(endpoint("responses")) |>
    req_headers(
      "Authorization" = paste("Bearer", API_KEY),
      "Content-Type"  = "application/json"
    )
}

ui <- fluidPage(
  titlePanel("GPT-5 in Shiny"),
  textAreaInput("user", "Ask anything:", rows = 4, placeholder = "e.g., Summarize ND vs AAb+ vs T1D differences..."),
  actionButton("go", "Send"),
  tags$hr(),
  strong("Model:"), textOutput("model", inline = TRUE),
  tags$br(),
  strong("Response:"),
  verbatimTextOutput("out"),
  tags$hr(),
  strong("Logs / errors:"),
  verbatimTextOutput("log")
)

server <- function(input, output, session) {
  output$model <- renderText({"gpt-5"})
  log_txt <- reactiveVal("")

  observeEvent(input$go, {
    output$out <- renderText({ "" })
    if (API_KEY == "") {
      log_txt(paste0(log_txt(), "\nMissing KEY in ~/.Renviron"))
      return()
    }

    # Build payload for Responses API (non-streamed first)
    payload <- list(
      model = "gpt-5",
      input = list(
        list(
          role = "user",
          content = input$user %||% "Say hello and confirm you are responding from the Responses API."
        )
      )
    )

    # Try a non-streamed request (simple and robust)
    resp <- tryCatch({
      req_base() |>
        req_body_json(payload) |>
        req_perform()
    }, error = function(e) e)

    if (inherits(resp, "error")) {
      log_txt(paste0(log_txt(), "\nRequest error: ", resp$message))
      output$log <- renderText({ log_txt() })
      return()
    }

    if (resp_status(resp) >= 300) {
      log_txt(paste0(
        log_txt(), "\nHTTP ", resp_status(resp), ": ",
        paste0(resp_body_string(resp), collapse = "")
      ))
      output$log <- renderText({ log_txt() })
      return()
    }

    # Parse response (Responses API returns an array of output items)
    body <- resp_body_json(resp)
    # The text can be in output_text (helper field) or within output[] items.
    text <- if (!is.null(body$output_text)) {
      body$output_text
    } else if (!is.null(body$output) && length(body$output) > 0) {
      paste(vapply(body$output, function(x) x$content[[1]]$text %||% "", character(1)), collapse = "")
    } else {
      "(no text returned)"
    }
    output$out <- renderText({ text })
    output$log <- renderText({ log_txt() })
  })
}

shinyApp(ui, server)

# Notes on parsing: the Responses API returns a structured object with an output array and often a convenience output_text field; the code above checks both. See docs if your account returns slightly different fields.

# For live token streaming, set "stream" = TRUE and consume server-sent events (SSE). Shiny doesn’t natively consume SSE, so the simplest approach is to poll a server cache or use a websocket. If you do want pure R-side streaming with httr2::req_stream() and push partial text to reactiveVal() on a timer, follow the streaming section in the Responses docs and adapt to Shiny’s reactivity.

# Verification:

library(httr2); library(jsonlite)
endpoint_base <- function(path) {
  base <- sub("/+\\z", "", Sys.getenv("BASE", "https://api.ai.it.ufl.edu/v1"))
  paste0(base, "/", path)
}
req <- request(endpoint_base("responses")) |>
  req_headers(
    "Authorization" = paste("Bearer", Sys.getenv("KEY")),
    "Content-Type"  = "application/json"
  ) |>
  req_body_json(list(
    model = "gpt-5",
    input = list(list(role="user", content="Return the string: OK"))
  ))
resp <- req_perform(req)
jsonlite::fromJSON(resp_body_string(resp))$output_text
