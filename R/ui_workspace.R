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
#
# Visual styling lives in `www/styles.css`; reusable primitives in
# `R/ui_components.R`.
# ============================================================================

workspace_ui <- function() {
  shiny::tagList(
    shiny::uiOutput("workspace_annotation_banner"),
    shiny::uiOutput("workspace_messages", class = "workspace-messages"),
    shiny::uiOutput("workspace_body")
  )
}

workspace_server <- function(input, output, session, state) {

  # Persistent active-annotation context strip. Always shown when a dataset
  # is loaded so users never lose track of which annotation context the rest
  # of the app is reading from. The string itself is built by
  # `annotation_set_label()` so its format is centralised.
  output$workspace_annotation_banner <- shiny::renderUI({
    if (is.null(state$active_dataset)) return(NULL)
    set <- get_active_annotation(state)
    label <- annotation_set_label(set)
    tone  <- if (is.null(set))        "neutral"
             else if (isTRUE(set$is_demo)) "demo"
             else                           "active"
    meta_text <- sprintf("dataset: %s",
                         state$active_dataset$name %||% "(unnamed)")
    context_strip(label, meta = meta_text, tone = tone)
  })

  output$workspace_messages <- shiny::renderUI({
    msgs <- state$messages
    if (length(msgs) == 0) return(NULL)
    tail_msgs <- utils::tail(msgs, 3)
    shiny::tagList(lapply(tail_msgs, function(m) {
      tone <- switch(m$level,
                     "success" = "success",
                     "warning" = "warning",
                     "error"   = "error",
                     "info")
      shiny::div(
        class = paste0("workspace-message workspace-message--", tone),
        m$text)
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
  empty_state(
    title   = "Welcome to scRNA Explorer",
    lede    = paste("Load a dataset from the sidebar to start exploring.",
                    "Selections you make are shared across every module."),
    icon    = "database",
    hint    = "Tip: try \"Load mock dataset\" to explore the UI without real data."
  )
}

needs_inputs_ui <- function(mod) {
  app_card(
    title   = mod$name,
    caption = "Inputs required",
    info_banner(
      tone  = "info",
      title = "This module needs a few selections before it can run.",
      req_list(mod$required_inputs)
    ),
    shiny::p(class = "app-card__caption",
             "Open the Basic scRNA Explorer to set them.")
  )
}

unknown_module_ui <- function(id) {
  app_card(
    title = "Unknown module",
    caption = "configuration",
    info_banner(
      tone  = "danger",
      title = "No module is registered with that id.",
      shiny::tags$code(id)
    )
  )
}
