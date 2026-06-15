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

# ---- One-shot preflight banner --------------------------------------------
# Console-only summary of which optional dependency tiers are complete.
# Silent when everything is OK; otherwise points the developer at the
# setup script. Never installs anything; never blocks the launch.
local({
  status <- sce_check_setup()
  any_missing <- !all(vapply(status, `[[`, logical(1), "complete"))
  if (any_missing) {
    message("\n", sce_preflight_message(status), "\n",
            "Run `Rscript scripts/setup_dev.R --full` once to enable the ",
            "complete feature surface (or `--demo` for just PBMC 8k).\n")
  }
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

  # Dataset loading.
  #
  # The "Load demo dataset" sidebar button drives a 3-tier resolution
  # chain implemented in R/demo_dataset.R:
  #
  #   1. Prepared artifact present  -> load it directly (instant).
  #   2. Artifact missing + the build packages are installed (Seurat +
  #      Bioc deps) + auto-build is not opted out by env var ->
  #      build the artifact in-process via `ensure_demo_dataset()`,
  #      wrapped in `shiny::withProgress()` for live status, then load.
  #      The first build pulls PBMC 8k counts via ExperimentHub (cached
  #      after first run) so subsequent clicks across sessions go
  #      straight to tier 1.
  #   3. Artifact missing AND auto-build cannot run (missing packages
  #      or `SCE_AUTO_BUILD_DEMO=0`) -> push a clear workspace warning
  #      naming what's missing and fall back to the synthetic
  #      `mock_dataset()` so the first-run UI experience always works.
  #
  # The input id stays `load_mock_dataset` for back-compat (tests /
  # external scripts may bind to it).
  observeEvent(input$load_mock_dataset, {
    # Tier 1: existing artifact.
    if (demo_dataset_exists()) {
      ds <- tryCatch(load_demo_dataset(),
                     error = function(e) {
                       push_message(state, sprintf(
                         "Could not load prepared demo artifact: %s",
                         conditionMessage(e)), "warning")
                       NULL
                     })
      if (!is.null(ds)) {
        set_active_dataset(state, ds)
        return()
      }
    }

    # Tier 2: auto-build.
    if (demo_auto_build_enabled() && can_build_demo_dataset()) {
      push_message(state, paste0(
        "Building prepared PBMC 8k demo artifact (one-time setup; ",
        "first run downloads ~30 MB via ExperimentHub, then ~30-90s of ",
        "Seurat preprocessing). Subsequent sessions reuse the artifact."),
        "info")
      ds <- tryCatch(
        shiny::withProgress(
          message = "Preparing PBMC 8k demo dataset",
          detail  = "first-run, one-time cost",
          value   = 0,
          {
            ensure_demo_dataset(progress = function(fraction, detail = NULL) {
              shiny::setProgress(value = fraction, detail = detail)
            })
          }),
        error = function(e) {
          push_message(state, sprintf(
            "Demo build failed: %s. Falling back to mock_dataset().",
            conditionMessage(e)), "warning")
          NULL
        })
      if (!is.null(ds)) {
        set_active_dataset(state, ds)
        return()
      }
    } else {
      # Tier 3: explain why we're not auto-building.
      reason <- if (!demo_auto_build_enabled())
                  "SCE_AUTO_BUILD_DEMO is disabled in the environment"
                else demo_auto_build_status()
      push_message(state,
                   sprintf(paste0(
                     "Prepared PBMC 8k demo artifact not found at %s and ",
                     "auto-build was skipped (%s). Falling back to ",
                     "mock_dataset(). Run `Rscript scripts/build_pbmc8k_demo.R` ",
                     "(or install the missing packages) to enable the ",
                     "PBMC 8k demo."),
                     demo_dataset_path(), reason),
                   "warning")
    }

    # Final fallback: synthetic mock dataset.
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
