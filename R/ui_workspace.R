# ============================================================================
# Main workspace
# ----------------------------------------------------------------------------
# The workspace renders one of:
#   - empty state (no dataset loaded)
#   - "needs input" state (module's required inputs aren't satisfied)
#   - the active module's UI
#
# All module server functions are started once at app start in `app.R`; only
# the UI swaps on navigation. This keeps each module's reactive graph alive
# so navigating away and back is instant.
# ============================================================================

workspace_ui <- function() {
  shiny::tagList(
    shiny::uiOutput("workspace_annotation_banner"),
    shiny::uiOutput("workspace_messages"),
    shiny::uiOutput("workspace_body")
  )
}

workspace_server <- function(input, output, session, state) {

  # Persistent active-annotation banner. Always shown when a dataset is
  # loaded so users never lose track of which annotation context the rest
  # of the app is reading from. The string itself is built by
  # `annotation_set_label()` so its format is centralised.
  output$workspace_annotation_banner <- shiny::renderUI({
    if (is.null(state$active_dataset)) return(NULL)
    set <- get_active_annotation(state)
    label <- annotation_set_label(set)
    bg <- if (is.null(set)) "#f5f5f5" else
          if (isTRUE(set$is_demo)) "#fff8e1" else "#e8f4fd"
    fg <- if (is.null(set)) "#666" else
          if (isTRUE(set$is_demo)) "#5d4037" else "#084298"
    shiny::div(style = sprintf(
        "padding:6px 12px; background:%s; color:%s; border-radius:4px; font-size:12px; margin-bottom:8px;",
        bg, fg),
      shiny::tags$strong("Active annotation: "),
      label,
      shiny::tags$span(style = "float:right; color:#888;",
                       sprintf("dataset: %s",
                               state$active_dataset$name %||% "(unnamed)"))
    )
  })

  output$workspace_messages <- shiny::renderUI({
    msgs <- state$messages
    if (length(msgs) == 0) return(NULL)
    tail_msgs <- utils::tail(msgs, 3)
    shiny::tagList(lapply(tail_msgs, function(m) {
      color <- switch(m$level,
                      "success" = "#2e7d32",
                      "warning" = "#ed6c02",
                      "error"   = "#c62828",
                      "#1565c0")
      shiny::div(
        style = sprintf(
          "border-left:3px solid %s; padding:6px 10px; margin-bottom:6px; background:#f7f7f7; font-size:13px;",
          color),
        m$text
      )
    }))
  })

  output$workspace_body <- shiny::renderUI({
    if (is.null(state$active_dataset))   return(empty_state_ui())
    mod <- get_module(state$active_module)
    if (is.null(mod))                    return(unknown_module_ui(state$active_module))
    if (!mod$enabled)                    return(mod$ui_fn(mod$id))   # coming-soon panel
    if (!module_inputs_ready(mod, state)) return(needs_inputs_ui(mod))
    mod$ui_fn(mod$id)
  })
}

empty_state_ui <- function() {
  shiny::div(
    style = "text-align:center; padding:80px 24px; color:#666;",
    shiny::h2("Welcome to scRNA Explorer"),
    shiny::p("Load a dataset from the sidebar to get started."),
    shiny::p(style = "color:#999; font-size:13px;",
             "Tip: try \"Load mock dataset\" to explore the UI without real data.")
  )
}

needs_inputs_ui <- function(mod) {
  shiny::div(
    style = "padding:32px; color:#666;",
    shiny::h3(mod$name),
    shiny::p("This module needs the following selections before it can run:"),
    shiny::tags$ul(lapply(mod$required_inputs, shiny::tags$li)),
    shiny::p(style = "color:#888;",
             "Open \"Basic scRNA Explorer\" to set them.")
  )
}

unknown_module_ui <- function(id) {
  shiny::div(
    style = "padding:32px; color:#a00;",
    shiny::h3("Unknown module"),
    shiny::p(sprintf("No module is registered with id '%s'.", id))
  )
}
