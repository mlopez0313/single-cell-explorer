# ============================================================================
# Module: Cell Type Annotation & Marker Discovery  (ENABLED)
# ----------------------------------------------------------------------------
# Manages multiple annotation sets backed by the registry of annotation
# engines in `ANNOTATION_ENGINES()` (R/annotation_registry.R). Sets are
# stored as `annotation_result_v1` objects in `state$annotation_sets` and
# exposed to the rest of the app through `get_active_annotation(state)`.
#
# UI surfaces:
#   - Set manager (create / activate / rename / duplicate / freeze / delete)
#   - Engine picker (manual + marker_score; future engines plug in here
#     without UI changes)
#   - Per-engine controls + Run
#   - Editable per-cluster table (for the manual engine, seeded by the
#     marker_score result if one exists)
#   - Apply to dataset metadata (provenance-named column, no overwrite)
#   - Download set CSV
#
# Per-cluster edits are immediately re-projected to per-cell labels through
# the manual engine, so the canonical storage stays per-cell regardless of
# how the user interacts with the table.
# ============================================================================

mod_annotation_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    page_header(
      eyebrow = "Cell Identity",
      title   = "Cell Type Annotation & Marker Discovery",
      lede    = shiny::tagList(
        "Manage one or more annotation sets. Downstream modules read the ",
        "active set via ",
        shiny::tags$code("get_active_annotation(state)"), ".")
    ),

    # -- Set manager -----------------------------------------------------
    app_card(
      title   = "Annotation sets",
      caption = "create / activate / duplicate / freeze / delete",
      shiny::fluidRow(
        shiny::column(4, shiny::uiOutput(ns("active_set_ui"))),
        shiny::column(8,
          beside_input(
            action_row(
              shiny::actionButton(ns("new_set"),    "+ New set",
                                  class = "btn btn-default"),
              shiny::actionButton(ns("dup_set"),    "Duplicate",
                                  class = "btn btn-default"),
              shiny::actionButton(ns("rename_set"), "Rename",
                                  class = "btn btn-default"),
              shiny::actionButton(ns("freeze_set"), "Toggle freeze",
                                  class = "btn btn-default"),
              shiny::actionButton(ns("delete_set"), "Delete",
                                  class = "btn btn-danger")
            )
          )
        )
      ),
      shiny::uiOutput(ns("active_set_summary"))
    ),

    # -- Engine + controls ----------------------------------------------
    control_panel(
      title = "Engine + run",
      shiny::fluidRow(
        shiny::column(4, shiny::selectInput(
          ns("engine_id"), "Annotation engine",
          choices  = list_annotation_engines(),
          selected = "manual")),
        shiny::column(4, shiny::uiOutput(ns("cluster_field_ui"))),
        shiny::column(4, shiny::uiOutput(ns("engine_params_ui")))
      ),
      actions = shiny::tagList(
        shiny::actionButton(ns("run_engine"),          "Run engine",
                            class = "btn btn-primary"),
        shiny::actionButton(ns("apply_edits"),         "Apply table edits (manual)",
                            class = "btn btn-success"),
        shiny::actionButton(ns("apply_to_metadata"),   "Apply to dataset metadata",
                            class = "btn btn-default"),
        shiny::downloadButton(ns("export"),            "Download CSV",
                              class = "btn btn-default")
      )
    ),

    shiny::uiOutput(ns("status_banner")),
    shiny::uiOutput(ns("warning")),

    # -- Per-cluster editable table ------------------------------------
    app_card(
      title   = "Per-cluster labels",
      caption = "edits expand to per-cell labels on Apply",
      callout_legend(
        items = list("engine-suggested" = "warning",
                     "user-confirmed"   = "success"),
        note  = "Per-cluster edits expand to per-cell labels on Apply."
      ),
      shiny::uiOutput(ns("table"))
    )
  )
}

mod_annotation_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- Set manager ----------------------------------------------------

    output$active_set_ui <- shiny::renderUI({
      sets <- state$annotation_sets
      if (!length(sets)) {
        return(beside_input(
          helper_text(
            shiny::em("No annotation sets yet. Click ",
                      shiny::strong("+ New set"), "."))))
      }
      choices <- setNames(
        vapply(sets, `[[`, character(1), "set_id"),
        vapply(sets, function(s) s$name %||% s$set_id, character(1))
      )
      shiny::selectInput(ns("active_id"), "Active annotation set",
                         choices  = choices,
                         selected = state$active_annotation_id %||% choices[1])
    })

    shiny::observeEvent(input$active_id, {
      if (!is.null(input$active_id) && nzchar(input$active_id)) {
        set_active_annotation(state, input$active_id)
      }
    }, ignoreInit = TRUE)

    output$active_set_summary <- shiny::renderUI({
      set <- get_active_annotation(state)
      if (is.null(set))
        return(helper_text(shiny::em("No active annotation set.")))
      stale <- if (!is.null(set$cluster_field_used) &&
                   !is.na(set$cluster_field_used) &&
                   !is.null(state$active_dataset)) {
        cur_n <- length(unique(as.character(get_metadata(state$active_dataset,
                                                         set$cluster_field_used) %||%
                                            character())))
        !identical(set$n_clusters_at_creation, cur_n)
      } else FALSE
      shiny::tagList(
        shiny::div(class = "set-summary",
          shiny::strong("Active: "), annotation_set_label(set),
          shiny::tags$br(),
          helper_text(sprintf(
            "id=%s | engine=%s | cluster_field=%s | n_clusters_at_creation=%s | registry=%s",
            set$set_id,
            set$engine_id %||% "?",
            set$cluster_field_used %||% "NA",
            set$n_clusters_at_creation %||% "NA",
            set$marker_registry_version %||% "NA"))
        ),
        if (stale)
          info_banner(
            tone  = "warning",
            title = "Stale-set warning",
            "Cluster cardinality has changed since this set was built.")
      )
    })

    shiny::observeEvent(input$new_set, {
      if (is.null(state$active_dataset)) {
        push_message(state, "Load a dataset first.", "warning"); return()
      }
      id <- new_annotation_set_id("set")
      cf <- input$cluster_field %||% state$selected_metadata_field
      cf <- cf %||% (available_metadata_fields(state$active_dataset)[1] %||% "cluster")
      set <- run_annotation_engine(
        engine_id = "manual",
        dataset   = state$active_dataset, state = state,
        params    = list(cluster_field = cf, labels = list()),
        set_id    = id,
        set_name  = sprintf("Set %s", format(Sys.time(), "%H:%M:%S")),
        description = "")
      add_annotation_set(state, set)
      set_active_annotation(state, id)
      shiny::updateSelectInput(session, "active_id", selected = id)
      push_message(state, sprintf("Created annotation set '%s'.", set$name),
                   "success")
    })

    shiny::observeEvent(input$dup_set, {
      cur <- state$active_annotation_id; if (is.null(cur)) return()
      new_id <- duplicate_annotation_set(state, cur)
      shiny::updateSelectInput(session, "active_id", selected = new_id)
      push_message(state, "Duplicated active set.", "success")
    })

    shiny::observeEvent(input$rename_set, {
      cur <- state$active_annotation_id; if (is.null(cur)) return()
      shiny::showModal(shiny::modalDialog(
        title = "Rename annotation set",
        shiny::textInput(ns("rename_name"), "New name",
                         value = state$annotation_sets[[cur]]$name),
        footer = shiny::tagList(
          shiny::modalButton("Cancel"),
          shiny::actionButton(ns("rename_confirm"), "Rename",
                              class = "btn btn-primary")
        )))
    })
    shiny::observeEvent(input$rename_confirm, {
      cur <- state$active_annotation_id
      if (is.null(cur) || is.null(input$rename_name) ||
          !nzchar(trimws(input$rename_name))) {
        shiny::removeModal(); return()
      }
      tryCatch({
        rename_annotation_set(state, cur, trimws(input$rename_name))
        push_message(state, "Renamed.", "success")
      }, error = function(e) push_message(state, conditionMessage(e), "error"))
      shiny::removeModal()
    })

    shiny::observeEvent(input$freeze_set, {
      cur <- state$active_annotation_id; if (is.null(cur)) return()
      cur_set <- state$annotation_sets[[cur]]
      tryCatch({
        freeze_annotation_set(state, cur, !isTRUE(cur_set$is_frozen))
        push_message(state,
          sprintf("Set is now %s.",
                  if (isTRUE(state$annotation_sets[[cur]]$is_frozen))
                    "frozen" else "unfrozen"),
          "info")
      }, error = function(e) push_message(state, conditionMessage(e), "error"))
    })

    shiny::observeEvent(input$delete_set, {
      cur <- state$active_annotation_id; if (is.null(cur)) return()
      tryCatch({
        remove_annotation_set(state, cur)
        push_message(state, "Deleted active set.", "info")
      }, error = function(e) push_message(state, conditionMessage(e), "error"))
    })

    # ---- Cluster field + engine parameters -----------------------------

    output$cluster_field_ui <- shiny::renderUI({
      fields <- available_metadata_fields(state$active_dataset)
      # Hide annotation-provenance columns from the picker so users don't
      # recursively annotate their own annotations.
      fields <- fields[!grepl("^annotation__", fields)]
      active <- get_active_annotation(state)
      sel <- active$cluster_field_used %||%
             state$selected_metadata_field %||% fields[1]
      shiny::selectInput(ns("cluster_field"), "Cluster field",
                         choices  = fields,
                         selected = sel)
    })

    output$engine_params_ui <- shiny::renderUI({
      eng <- get_annotation_engine(input$engine_id)
      if (is.null(eng) || identical(eng$id, "manual")) return(NULL)
      if (identical(eng$id, "marker_score")) {
        species_choices <- c("", sort(unique(vapply(
          state$marker_registry$entries %||% list(),
          function(e) e$species %||% NA_character_, character(1)))))
        return(shiny::tagList(
          shiny::selectInput(ns("species"), "Species filter",
                             choices  = species_choices, selected = ""),
          shiny::numericInput(ns("min_score"), "Min top-score (else Unknown)",
                              value = 0, step = 0.1)
        ))
      }
      helper_text(shiny::em(sprintf("Engine '%s' has no extra params.",
                                    input$engine_id)))
    })

    # ---- Per-cluster table -------------------------------------------

    cluster_ids <- shiny::reactive({
      shiny::req(state$active_dataset, input$cluster_field)
      get_cluster_ids(state$active_dataset, input$cluster_field)
    })

    output$table <- shiny::renderUI({
      ids <- tryCatch(cluster_ids(), error = function(e) character())
      if (!length(ids))
        return(friendly_warning("No clusters in the chosen field."))

      set      <- get_active_annotation(state)
      # Extract a column from cluster_summary by name, falling back to a
      # vector of `default` (one per id) when either the set has no
      # cluster_summary or the column is absent. This is what kept the
      # UI from crashing with `attempt to set an attribute on NULL`
      # when an engine (Azimuth / CellTypist used to be the offenders)
      # produced a cluster_summary without the canonical `top_score`
      # column.
      summary_col <- function(col, default) {
        if (is.null(set) || is.null(set$cluster_summary) ||
            !col %in% names(set$cluster_summary)) {
          return(setNames(rep(default, length(ids)), ids))
        }
        cs   <- set$cluster_summary
        keys <- as.character(cs$cluster)
        vals <- cs[[col]]
        if (is.null(vals))   # belt + braces: cs$cluster present but the
          return(setNames(rep(default, length(ids)), ids))   # col is NULL
        setNames(vals, keys)
      }
      cur_map   <- summary_col("top_label", NA_character_)
      cur_score <- summary_col("top_score", NA_real_)
      n_cells <- count_cells_per_cluster(state$active_dataset,
                                         input$cluster_field)

      header <- shiny::tags$tr(
        lapply(c("Cluster", "n cells", "Current label", "Score",
                 "Edit label", "Notes"),
               shiny::tags$th))
      rows <- lapply(ids, function(cl) {
        lab <- cur_map[cl]
        is_confirmed <- !is.na(lab) && nzchar(lab) && !identical(lab, "Unknown")
        row_class <- if (is_confirmed) "is-confirmed" else "is-suggested"
        shiny::tags$tr(class = row_class,
          shiny::tags$td(shiny::strong(cl)),
          shiny::tags$td(class = "sce-tabular",
                         format(n_cells[cl] %||% NA_integer_, big.mark = ",")),
          shiny::tags$td(if (is.na(lab) || !nzchar(lab)) "-" else lab),
          shiny::tags$td(class = "sce-tabular",
                         if (is.na(cur_score[cl])) "-"
                         else sprintf("%.2f", cur_score[cl])),
          shiny::tags$td(shiny::textInput(
            ns(paste0("edit_", cl)), label = NULL,
            value = if (is.na(lab)) "" else as.character(lab),
            placeholder = "type a label...")),
          shiny::tags$td(shiny::textInput(
            ns(paste0("note_", cl)), label = NULL,
            value = "", placeholder = "optional"))
        )
      })
      shiny::tags$table(class = "cluster-table",
                        shiny::tags$thead(header),
                        shiny::tags$tbody(rows))
    })

    # ---- Run engine ---------------------------------------------------

    shiny::observeEvent(input$run_engine, {
      if (is.null(state$active_dataset)) {
        push_message(state, "Load a dataset first.", "warning"); return()
      }
      cur_id <- state$active_annotation_id
      if (is.null(cur_id)) {
        push_message(state, "Create or select an annotation set first.",
                     "warning"); return()
      }
      cur_set <- state$annotation_sets[[cur_id]]
      if (isTRUE(cur_set$is_frozen)) {
        push_message(state, "Active set is frozen.", "warning"); return()
      }
      eng <- input$engine_id %||% "manual"
      params <- switch(eng,
        "manual" = list(cluster_field = input$cluster_field, labels = .collect_table_labels(input, cluster_ids())),
        "marker_score" = list(cluster_field = input$cluster_field,
                              species       = if (nzchar(input$species %||% "")) input$species else NA_character_,
                              min_score     = input$min_score %||% 0),
        list(cluster_field = input$cluster_field)
      )
      new_set <- tryCatch(
        run_annotation_engine(eng, state$active_dataset, state, params,
                              set_id        = cur_id,
                              set_name      = cur_set$name,
                              description   = cur_set$description,
                              parent_set_id = cur_set$parent_set_id,
                              is_demo       = isTRUE(cur_set$is_demo)),
        error = function(e) e)
      if (inherits(new_set, "error")) {
        push_message(state, sprintf("Annotation engine failed: %s",
                                    conditionMessage(new_set)), "error")
        return()
      }
      # Preserve created_at; bump modified_at.
      new_set$created_at <- cur_set$created_at %||% Sys.time()
      new_set$modified_at <- Sys.time()
      add_annotation_set(state, new_set)
      push_message(state, sprintf(
        "Ran '%s' on %d clusters of '%s'.",
        eng, new_set$n_clusters_at_creation %||% 0L,
        new_set$cluster_field_used %||% "?"), "success")
    })

    shiny::observeEvent(input$apply_edits, {
      cur_id <- state$active_annotation_id; if (is.null(cur_id)) return()
      cur_set <- state$annotation_sets[[cur_id]]
      if (isTRUE(cur_set$is_frozen)) {
        push_message(state, "Active set is frozen.", "warning"); return()
      }
      params <- list(cluster_field = input$cluster_field,
                     labels = .collect_table_labels(input, cluster_ids()))
      new_set <- tryCatch(
        run_annotation_engine("manual", state$active_dataset, state, params,
                              set_id        = cur_id,
                              set_name      = cur_set$name,
                              description   = cur_set$description,
                              parent_set_id = cur_set$parent_set_id,
                              is_demo       = isTRUE(cur_set$is_demo)),
        error = function(e) e)
      if (inherits(new_set, "error")) {
        push_message(state, conditionMessage(new_set), "error"); return()
      }
      new_set$created_at <- cur_set$created_at %||% Sys.time()
      new_set$modified_at <- Sys.time()
      add_annotation_set(state, new_set)
      push_message(state, "Table edits applied (per-cluster -> per-cell).",
                   "success")
    })

    # ---- Apply to dataset metadata -----------------------------------

    shiny::observeEvent(input$apply_to_metadata, {
      set <- get_active_annotation(state)
      if (is.null(set)) {
        push_message(state, "No active annotation set.", "warning"); return()
      }
      ds  <- state$active_dataset; if (is.null(ds)) return()
      ds2 <- tryCatch(apply_annotations_to_dataset(ds, set),
                      error = function(e) e)
      if (inherits(ds2, "error")) {
        push_message(state, conditionMessage(ds2), "error"); return()
      }
      state$active_dataset <- ds2
      push_message(state, sprintf(
        "Applied set '%s' as a new metadata column. Use it in any module's color-by-metadata picker.",
        set$name %||% set$set_id), "success")
    })

    # ---- Banners / warnings ------------------------------------------

    output$warning <- shiny::renderUI({
      if (is.null(state$active_dataset))
        return(friendly_warning("No dataset loaded."))
      if (!length(state$annotation_sets))
        return(friendly_warning("No annotation sets yet. Click '+ New set' above."))
      NULL
    })

    output$status_banner <- shiny::renderUI({
      set <- get_active_annotation(state)
      if (is.null(set)) return(NULL)
      n_labelled <- sum(!is.na(set$cell_labels) & nzchar(set$cell_labels) &
                        set$cell_labels != "Unknown")
      n_total    <- length(set$cell_labels)
      status_banner(
        shiny::span(
          shiny::tags$strong("Per-cell labels: "),
          sprintf("%d / %d assigned (%.0f%%). ",
                  n_labelled, n_total,
                  if (n_total) 100 * n_labelled / n_total else 0),
          sprintf("engine=%s; registry=%s",
                  set$engine_id %||% "?",
                  set$marker_registry_version %||% "n/a")),
        tone  = "info",
        label = "Active set")
    })

    # ---- Export ------------------------------------------------------

    output$export <- shiny::downloadHandler(
      filename = function() {
        set <- get_active_annotation(state)
        sprintf("annotations_%s.csv",
                set$set_id %||% format(Sys.time(), "%Y%m%d_%H%M%S"))
      },
      content = function(file) {
        write_annotation_set_csv(get_active_annotation(state), file)
      }
    )
  })
}

# ---- Internal: collect per-cluster edits from text inputs -----------------
.collect_table_labels <- function(input, cluster_ids) {
  out <- setNames(vector("list", length(cluster_ids)), as.character(cluster_ids))
  for (cl in cluster_ids) {
    v <- input[[paste0("edit_", cl)]]
    out[[cl]] <- if (is.null(v)) NA_character_ else trimws(as.character(v))
  }
  out
}
