# ============================================================================
# Placeholder modules  (DISABLED)
# ----------------------------------------------------------------------------
# These modules are registered so they show up in the sidebar (greyed out),
# but their bodies are intentionally empty. When a contributor starts work on
# one, they should:
#   1. Move it into its own file: R/modules/mod_<id>.R
#   2. Replace the body with a real UI + server
#   3. Flip `enabled = TRUE` in R/registry.R
#
# The shared `coming_soon_ui()` keeps the look consistent until then.
# ============================================================================

#' Standard "coming soon" panel used by every disabled module's UI.
coming_soon_ui <- function(id, title, description) {
  ns <- shiny::NS(id)
  shiny::div(
    class = "coming-soon-panel",
    style = paste(
      "padding:32px; border:1px dashed #bbb; border-radius:8px;",
      "background:#fafafa; color:#555; max-width:720px;"
    ),
    shiny::h2(title),
    shiny::tags$span(
      style = "display:inline-block; padding:2px 8px; background:#ffe8a1; border-radius:4px; font-size:12px; font-weight:600;",
      "COMING SOON"
    ),
    shiny::p(style = "margin-top:16px;", description),
    shiny::p(style = "color:#888; font-size:13px;",
             "This module is registered but not yet implemented. ",
             "See docs/ADDING_MODULES.md to contribute.")
  )
}

# Empty server used by every disabled module. Kept as a function so individual
# modules can override it without touching this file.
noop_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {})
}

# ---- Cell identity --------------------------------------------------------
# mod_marker_investigation_* and mod_annotation_* moved to their own files in
# R/modules/ when they became real modules.

# ---- Statistics -----------------------------------------------------------
# mod_differential_expression_* moved to its own file
# (R/modules/mod_differential_expression.R) when it became a real module.

# ---- Functional -----------------------------------------------------------
# mod_pathway_analysis_* moved to its own file
# (R/modules/mod_pathway_analysis.R) when it became a real module.

# ---- Advanced -------------------------------------------------------------
# mod_imputation_* moved to its own file (R/modules/mod_imputation.R)
# when it became a real module.

# mod_trajectory_* moved to its own file (R/modules/mod_trajectory.R)
# when it became a real module.

# mod_regulons_* (was `mod_regulatory_*`) moved to its own file
# (R/modules/mod_regulons.R) when it became a real module.
