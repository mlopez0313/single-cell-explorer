# ============================================================================
# Sidebar navigation
# ----------------------------------------------------------------------------
# Renders the brand block, dataset-load region, and module list (grouped by
# category). Each entry is a button that sets `state$active_module`. Disabled
# modules are shown but visually muted and do not change the active module.
#
# All visual styling lives in `www/styles.css` (`.app-sidebar`,
# `.sidebar-nav__*`). Reusable presentation helpers come from
# `R/ui_components.R`.
# ============================================================================

#' Bootstrap icon name to associate with each module category.
#' Used for tasteful section icons when `bsicons` is installed.
.SCE_CATEGORY_ICONS <- c(
  "Overview"            = "speedometer2",
  "Exploration"         = "binoculars",
  "Cell Identity"       = "person-badge",
  "Statistics"          = "bar-chart",
  "Functional Analysis" = "diagram-3",
  "Advanced Analysis"   = "stars"
)

#' Per-module icon overrides. Falls back to the category icon when missing.
.SCE_MODULE_ICONS <- c(
  dataset_overview        = "table",
  scrna_explorer          = "scatter-chart",
  marker_investigation    = "search",
  annotation              = "person-badge",
  differential_expression = "sliders",
  pathway_analysis        = "diagram-3",
  imputation              = "droplet-half",
  trajectory              = "graph-up-arrow",
  regulons                = "share"
)

.sce_module_icon <- function(mod) {
  override <- .SCE_MODULE_ICONS[[mod$id]]
  if (!is.null(override)) return(override)
  unname(.SCE_CATEGORY_ICONS[mod$category])
}

#' Static sidebar shell (populated by the server-side `sidebar_server`).
sidebar_ui <- function() {
  shiny::tagList(
    # Brand block: signals "this is a productized app", not a prototype.
    shiny::div(
      class = "app-sidebar__brand",
      shiny::div(class = "app-sidebar__brand-mark",
                 app_icon("droplet", class = "")),
      shiny::div(
        class = "app-sidebar__brand-text",
        shiny::span(class = "app-sidebar__brand-name", "scRNA Explorer"),
        shiny::span(class = "app-sidebar__brand-sub",  "single-cell workspace")
      )
    ),
    # Dataset region.
    shiny::div(
      class = "sidebar-section",
      shiny::p(class = "sidebar-section__label", "Dataset"),
      # Input id stays "load_mock_dataset" for back-compat with tests and
      # external scripts; the visible label now reflects that the button
      # loads the prepared PBMC 8k demo artifact when present (with a
      # graceful fallback to the synthetic mock_dataset()).
      shiny::actionButton(
        "load_mock_dataset", "Load demo dataset",
        class = "btn btn-primary btn-sm", width = "100%"),
      shiny::div(class = "sidebar-dataset-status",
                 shiny::textOutput("sidebar_dataset_status", inline = TRUE))
    ),
    # Module navigation.
    shiny::div(
      class = "sidebar-section",
      shiny::p(class = "sidebar-section__label", "Modules"),
      shiny::uiOutput("sidebar_modules", class = "sidebar-nav")
    )
  )
}

#' Server logic that draws the module list and handles click events.
#'
#' @param input,output,session  Shiny app-level objects
#' @param state                 shared app state
sidebar_server <- function(input, output, session, state) {

  output$sidebar_dataset_status <- shiny::renderText({
    ds <- state$active_dataset
    if (is.null(ds)) "No dataset loaded." else sprintf("Loaded: %s", ds$name)
  })

  output$sidebar_modules <- shiny::renderUI({
    grouped <- modules_by_category()
    groups <- lapply(names(grouped), function(cat) {
      mods <- grouped[[cat]]
      shiny::tagList(
        shiny::div(class = "sidebar-nav__category", cat),
        lapply(mods, function(m) module_button(m, state$active_module))
      )
    })
    shiny::tagList(groups)
  })

  # Wire up one observer per module id. Done lazily via lapply so future
  # modules added to the registry get observers automatically.
  shiny::isolate({
    for (m in module_registry()) {
      local({
        mod <- m
        if (!mod$enabled) return()
        shiny::observeEvent(input[[paste0("nav_", mod$id)]], {
          state$active_module <- mod$id
        }, ignoreInit = TRUE)
      })
    }
  })
}

#' Render one sidebar entry. Enabled -> actionButton; disabled -> muted div.
module_button <- function(mod, active_id) {
  is_active <- identical(mod$id, active_id)
  icon_name <- .sce_module_icon(mod)
  if (!mod$enabled) {
    return(shiny::div(
      class    = "sidebar-nav__item is-disabled",
      title    = mod$description,
      `aria-disabled` = "true",
      app_icon(icon_name),
      shiny::span(class = "sidebar-nav__label", mod$name),
      shiny::span(class = "sidebar-nav__hint",  "soon")
    ))
  }
  classes <- c("sidebar-nav__item", if (is_active) "is-active")
  shiny::actionButton(
    inputId = paste0("nav_", mod$id),
    label   = htmltools::tagList(
      app_icon(icon_name),
      shiny::span(class = "sidebar-nav__label", mod$name)
    ),
    class   = paste(classes, collapse = " "),
    title   = mod$description,
    `aria-current` = if (is_active) "page"
  )
}
