# ============================================================================
# scRNA Explorer -- application entry point
# ----------------------------------------------------------------------------
# Run from R:                  shiny::runApp("scrna-explorer")
# Run from a shell:            R -e 'shiny::runApp("scrna-explorer", launch.browser = FALSE, port = 3838, host = "0.0.0.0")'
#
# This file is intentionally small. It:
#   1. Sources the R/ helpers and every R/modules/ file
#   2. Builds the static layout (sidebar + workspace) under a bslib theme
#   3. Initialises the shared app state
#   4. Starts every enabled module's server once
#
# Visual styling lives in `www/styles.css`; reusable presentation primitives
# in `R/ui_components.R`. Modules should compose those primitives instead of
# hand-rolling inline CSS.
# ============================================================================

library(shiny)
library(bslib)
library(htmltools)

# ---- Source helpers + modules ---------------------------------------------
# Order matters: state/registry helpers must be available before module files
# (which reference helper functions) and before app code that calls them.
local({
  here <- function(...) file.path("R", ...)
  # state.R defines `%||%` and the state factory; everything else may use it.
  source(here("state.R"), local = FALSE)
  # All other top-level helpers in R/ (dataset, dataset_helpers, plotting,
  # registry, ui_components, ui_sidebar, ui_workspace, ...). Adding a new
  # helper file here just means dropping it into R/.
  top_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
  top_files <- setdiff(top_files, here("state.R"))
  # Modules come last so their files can freely reference any helper above.
  module_files <- list.files(here("modules"), pattern = "\\.R$", full.names = TRUE)
  for (f in setdiff(top_files, module_files)) source(f, local = FALSE)
  for (f in module_files) source(f, local = FALSE)
})

# ---- Theme ----------------------------------------------------------------
# A restrained scientific theme: cool neutrals with a teal accent. The
# stylesheet in `www/styles.css` layers a richer design system on top.
app_theme <- bslib::bs_theme(
  version       = 5,
  bg            = "#f7f8fb",
  fg            = "#0f172a",
  primary       = "#1d6f6a",
  secondary     = "#5a6376",
  success       = "#2f7d32",
  info          = "#1e4ed8",
  warning       = "#b25400",
  danger        = "#b3261e",
  base_font     = bslib::font_collection(
    bslib::font_google("Inter", local = FALSE),
    "system-ui", "Segoe UI", "Roboto", "Helvetica Neue", "Arial", "sans-serif"
  ),
  code_font     = bslib::font_collection(
    bslib::font_google("JetBrains Mono", local = FALSE),
    "ui-monospace", "SFMono-Regular", "Menlo", "Consolas", "monospace"
  ),
  font_scale    = 0.95,
  "border-radius"      = "0.5rem",
  "card-border-color"  = "#e3e6ed",
  "card-cap-bg"        = "#fbfcfe",
  "body-color"         = "#0f172a",
  "body-bg"            = "#f7f8fb"
)

# ---- UI -------------------------------------------------------------------
ui <- tagList(
  htmltools::tags$head(
    htmltools::tags$meta(charset = "utf-8"),
    htmltools::tags$meta(name = "viewport",
                         content = "width=device-width, initial-scale=1"),
    htmltools::tags$link(rel = "stylesheet", href = "styles.css"),
    htmltools::tags$title("scRNA Explorer")
  ),
  bslib::page_fluid(
    theme  = app_theme,
    title  = "scRNA Explorer",
    htmltools::div(
      class = "app-shell",
      htmltools::tags$aside(class = "app-sidebar",   sidebar_ui()),
      htmltools::tags$main (class = "app-workspace",
        htmltools::div(class = "app-workspace__inner", workspace_ui())
      )
    )
  )
)

# ---- Server ---------------------------------------------------------------
server <- function(input, output, session) {
  state <- new_app_state()

  # Sidebar + workspace wiring
  sidebar_server(input, output, session, state)
  workspace_server(input, output, session, state)

  # Dataset loading (mock for now; real loaders live in R/dataset.R)
  observeEvent(input$load_mock_dataset, {
    set_active_dataset(state, mock_dataset())
  })

  # Start every enabled module's server exactly once. Module servers run for
  # the lifetime of the session even when their UI isn't visible -- this keeps
  # navigation instant and lets modules react to state changes in the
  # background.
  for (mod in module_registry()) {
    if (isTRUE(mod$enabled) && !is.null(mod$server_fn)) {
      mod$server_fn(mod$id, state)
    }
  }
}

shinyApp(ui, server)
