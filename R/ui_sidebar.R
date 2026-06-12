# ============================================================================
# Sidebar navigation
# ----------------------------------------------------------------------------
# Renders the module list, grouped by category. Each entry is a button that
# sets `state$active_module`. Disabled modules are shown but visually muted
# and do not change the active module on click.
# ============================================================================

#' Static sidebar shell (populated by the server-side `sidebar_server`).
sidebar_ui <- function() {
  shiny::tagList(
    shiny::div(class = "sidebar-section",
      shiny::h5("Dataset"),
      shiny::actionButton("load_mock_dataset", "Load mock dataset",
                          class = "btn btn-primary btn-sm", width = "100%"),
      shiny::div(style = "margin-top:8px; font-size:12px; color:#888;",
                 shiny::textOutput("sidebar_dataset_status", inline = TRUE))
    ),
    shiny::hr(),
    shiny::div(class = "sidebar-section",
      shiny::h5("Modules"),
      shiny::uiOutput("sidebar_modules")
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
      shiny::div(class = "sidebar-category",
        shiny::tags$div(
          style = "font-size:11px; text-transform:uppercase; letter-spacing:0.05em; color:#888; margin:12px 0 4px 0;",
          cat
        ),
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
  if (!mod$enabled) {
    return(shiny::div(
      class = "sidebar-module disabled",
      style = paste(
        "padding:6px 10px; margin:2px 0; border-radius:4px; color:#aaa;",
        "background:#f5f5f5; font-size:13px; cursor:not-allowed;",
        "display:flex; justify-content:space-between; align-items:center;"
      ),
      title = mod$description,
      shiny::span(mod$name),
      shiny::tags$span(style = "font-size:10px; background:#ddd; padding:1px 6px; border-radius:8px;",
                       "soon")
    ))
  }
  shiny::actionButton(
    inputId = paste0("nav_", mod$id),
    label   = mod$name,
    class   = if (is_active) "btn btn-primary btn-sm sidebar-module active"
              else            "btn btn-light    btn-sm sidebar-module",
    style   = "width:100%; text-align:left; margin:2px 0;",
    title   = mod$description
  )
}
