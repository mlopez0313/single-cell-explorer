# ============================================================================
# scRNA Explorer -- application entry point
# ----------------------------------------------------------------------------
# Run from R:                  shiny::runApp("scrna-explorer")
# Run from a shell:            R -e 'shiny::runApp("scrna-explorer", launch.browser = FALSE, port = 3838, host = "0.0.0.0")'
#
# This file is intentionally small. It:
#   1. Sources the R/ helpers and every R/modules/ file
#   2. Builds the static layout (sidebar + workspace)
#   3. Initialises the shared app state
#   4. Starts every enabled module's server once
# ============================================================================

library(shiny)

# ---- Source helpers + modules ---------------------------------------------
# Order matters: state/registry helpers must be available before module files
# (which reference helper functions) and before app code that calls them.
local({
  here <- function(...) file.path("R", ...)
  # state.R defines `%||%` and the state factory; everything else may use it.
  source(here("state.R"), local = FALSE)
  # All other top-level helpers in R/ (dataset, dataset_helpers, plotting,
  # registry, ui_sidebar, ui_workspace). Adding a new helper file here just
  # means dropping it into R/.
  top_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
  top_files <- setdiff(top_files, here("state.R"))
  # Modules come last so their files can freely reference any helper above.
  module_files <- list.files(here("modules"), pattern = "\\.R$", full.names = TRUE)
  for (f in setdiff(top_files, module_files)) source(f, local = FALSE)
  for (f in module_files) source(f, local = FALSE)
})

# ---- UI -------------------------------------------------------------------
ui <- fluidPage(
  title = "scRNA Explorer",
  tags$head(tags$style(HTML("
    body { background: #fff; }
    .app-shell { display: flex; min-height: 100vh; }
    .app-sidebar { width: 260px; padding: 16px; background: #fafafa; border-right: 1px solid #eee; }
    .app-workspace { flex: 1; padding: 24px 32px; }
    .app-sidebar h5 { margin: 0 0 8px 0; font-size: 14px; }
    .sidebar-module.active { font-weight: 600; }
  "))),
  div(class = "app-shell",
    div(class = "app-sidebar",   sidebar_ui()),
    div(class = "app-workspace", workspace_ui())
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
