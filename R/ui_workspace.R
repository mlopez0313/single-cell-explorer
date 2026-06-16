# ============================================================================
# Main workspace
# ----------------------------------------------------------------------------
# The workspace renders one of:
#   - empty state (no dataset loaded)
#   - persistent tabset (one tab per enabled module) once a dataset is loaded
#
# All module server functions are started once at app start in `app.R`. The
# workspace builds each module's UI ONCE on first dataset load and switches
# the visible pane via `updateTabsetPanel()` -- it never destroys a module's
# DOM. This is what lets inputs (numericInputs, actionButtons, selectInputs,
# computed tables, etc.) survive navigating away and back.
#
# Why not `renderUI({ mod$ui_fn(mod$id) })`? It LOOKS lightweight but it
# destroys the previous module's DOM on every tab switch. That wipes
# hard-coded input defaults, resets actionButton click counters, and causes
# `eventReactive(input$recompute, ...)` to appear "empty" because
# `input$recompute %||% 0L > 0L` flips back to FALSE on re-render. Symptom:
# "everything in Marker Investigation resets when I leave that tab".
#
# Visual styling lives in `www/styles.css`; reusable primitives in
# `R/ui_components.R`.
# ============================================================================

workspace_ui <- function() {
  shiny::tagList(
    shiny::uiOutput("workspace_annotation_banner"),
    shiny::uiOutput("workspace_messages", class = "workspace-messages"),
    shiny::uiOutput("workspace_needs_banner"),
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

  # Compact "this module needs inputs" banner shown above the persistent
  # tabset. We keep the module's own UI rendered (so the user can interact
  # with controls that contribute to satisfying the requirements) and add
  # a banner on top instead of replacing the whole pane with a static
  # placeholder. That avoids destroying the pane and matches the
  # persistent-tabset model.
  output$workspace_needs_banner <- shiny::renderUI({
    if (is.null(state$active_dataset)) return(NULL)
    mod <- get_module(state$active_module)
    if (is.null(mod) || !isTRUE(mod$enabled)) return(NULL)
    if (module_inputs_ready(mod, state)) return(NULL)
    needs_inputs_banner(mod)
  })

  # ---- Workspace body -----------------------------------------------------
  # `bindEvent(is.null(state$active_dataset))` makes this renderUI fire
  # only when the *presence* of a dataset flips (null -> loaded, or loaded
  # -> null). It does NOT fire when:
  #   - the user navigates between modules (active_module changes)
  #   - the user switches between two real datasets (still not-null)
  # so the persistent tabset built below stays mounted across module
  # navigation. That's the whole point of this rewrite.
  output$workspace_body <- shiny::renderUI({
    if (is.null(state$active_dataset)) return(empty_state_ui())
    build_workspace_tabset(state)
  }) |>
    shiny::bindEvent(is.null(state$active_dataset),
                     ignoreNULL = FALSE,
                     ignoreInit = FALSE)

  # Module navigation: switch the visible pane without touching the DOM.
  shiny::observeEvent(state$active_module, {
    if (is.null(state$active_dataset)) return()
    target <- state$active_module
    if (is.null(target) || !nzchar(target)) return()
    shiny::updateTabsetPanel(session, "workspace_active_module",
                             selected = target)
  }, ignoreInit = FALSE)
}

# Build the hidden tabset containing every enabled module's UI. Called once
# per workspace lifetime (i.e., once per dataset-load transition). The
# `type = "hidden"` flavor of `tabsetPanel` renders without a visible tab
# bar and is switched purely via `updateTabsetPanel(selected = ...)`.
build_workspace_tabset <- function(state) {
  enabled <- Filter(function(m) isTRUE(m$enabled) && !is.null(m$ui_fn),
                    module_registry())

  # Initial selection: prefer the currently-active module if it's a known
  # enabled module, otherwise the first one. We `isolate()` because we're
  # inside a renderUI and don't want this to depend on active_module
  # changes (those are handled by the observeEvent above).
  initial <- shiny::isolate(state$active_module)
  ids     <- vapply(enabled, `[[`, character(1), "id")
  if (is.null(initial) || !nzchar(initial) || !(initial %in% ids)) {
    initial <- ids[1L]
  }

  panels <- lapply(enabled, function(mod) {
    shiny::tabPanel(
      title = mod$name,
      value = mod$id,
      mod$ui_fn(mod$id)
    )
  })

  do.call(shiny::tabsetPanel,
          c(list(id = "workspace_active_module",
                 type = "hidden",
                 selected = initial),
            panels))
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

# Compact banner shown above the active module's UI when its required
# inputs aren't all satisfied. Replaces the previous `needs_inputs_ui()`
# full-card placeholder, which used to swap out the module DOM entirely
# (defeating state persistence on tab returns).
needs_inputs_banner <- function(mod) {
  info_banner(
    tone  = "info",
    title = sprintf("%s needs a few selections before it can run.", mod$name),
    req_list(mod$required_inputs)
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
