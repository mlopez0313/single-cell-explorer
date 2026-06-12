# ============================================================================
# Module: Regulons / Network Analysis  (ENABLED)
# ----------------------------------------------------------------------------
# Score TF regulons (TF -> target genes) per cell using AUCell and
# show:
#   - a regulon x cluster heatmap (mean AUC by group)
#   - an embedding colored by a selected regulon's AUC
#
# Result lives in `state$analysis_results$regulons`. It is NEVER
# silently injected into Explorer / DE / etc. metadata; a future
# `apply_regulon_scores_to_dataset()` would mirror
# `apply_pseudotime_to_dataset()` / `apply_annotations_to_dataset()`.
#
# Sources of regulons (mock built-in, DoRothEA, future SCENIC-derived)
# are pluggable via `REGULON_SOURCES()`; scoring engines (pure-R
# AUCell, Bioc AUCell, future SCENIC) via `REGULON_ENGINES()`. The UI
# is registry-driven; adding a new engine or source surfaces here
# automatically.
# ============================================================================

mod_regulons_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h2("Regulons / Network Analysis"),

    # -- Demo warning ---------------------------------------------------
    shiny::div(
      style = paste(
        "padding:10px 14px; background:#fff3cd; color:#664d03;",
        "border:1px solid #ffecb5; border-radius:4px; font-size:13px;",
        "margin-bottom:12px;"),
      shiny::tags$strong("Regulon activity \u2014 exploratory."),
      " Scores are per-cell enrichment of each TF's target set in the ",
      "cell's top-expressed genes (AUCell). The regulons themselves ",
      "come from a curated source (mock / DoRothEA / future SCENIC ",
      "runs), not from the data."
    ),

    # -- Controls -------------------------------------------------------
    shiny::fluidRow(
      shiny::column(3, shiny::uiOutput(ns("source_ui"))),
      shiny::column(3, shiny::uiOutput(ns("engine_ui"))),
      shiny::column(3, shiny::uiOutput(ns("group_field_ui"))),
      shiny::column(3, shiny::uiOutput(ns("top_n_ui")))
    ),
    shiny::div(style = "margin: 8px 0 16px 0;",
      shiny::actionButton(ns("run"), "Score regulons",
                          class = "btn btn-primary"),
      shiny::tags$span(style = "margin-left:12px; color:#888; font-size:12px;",
                       "Scoring runs only when this button is clicked.")
    ),
    shiny::uiOutput(ns("status_banner")),
    shiny::uiOutput(ns("input_warning")),

    shiny::hr(),

    shiny::fluidRow(
      shiny::column(7,
        shiny::h4("Regulon \u00d7 group heatmap (mean AUC)"),
        shiny::uiOutput(ns("heatmap_warning")),
        shiny::plotOutput(ns("heatmap"), height = "420px")
      ),
      shiny::column(5,
        shiny::h4(shiny::textOutput(ns("emb_title"), inline = TRUE)),
        shiny::uiOutput(ns("regulon_picker_ui")),
        shiny::uiOutput(ns("reduction_ui")),
        shiny::uiOutput(ns("emb_warning")),
        shiny::plotOutput(ns("emb_plot"), height = "360px")
      )
    )
  )
}

mod_regulons_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- Slot helpers --------------------------------------------------
    rg_slot     <- function() state$analysis_results$regulons
    set_rg_slot <- function(x) {
      ar <- state$analysis_results
      ar$regulons <- x
      state$analysis_results <- ar
    }

    # ---- Controls ------------------------------------------------------
    output$source_ui <- shiny::renderUI({
      shiny::selectInput(ns("source"), "Regulon source",
                         choices = list_regulon_sources(),
                         selected = "mock_pbmc")
    })
    output$engine_ui <- shiny::renderUI({
      shiny::selectInput(ns("engine"), "Scoring engine",
                         choices = list_regulon_engines(),
                         selected = "aucell_r")
    })
    output$group_field_ui <- shiny::renderUI({
      ds <- state$active_dataset; if (is.null(ds)) return(NULL)
      cat <- available_categorical_metadata_fields(ds)
      shiny::selectInput(ns("group_field"), "Group by (heatmap rows)",
                         choices = cat,
                         selected = if ("cluster" %in% cat) "cluster"
                                    else cat[1])
    })
    output$top_n_ui <- shiny::renderUI({
      shiny::numericInput(ns("top_n_fraction"), "Top-N (fraction)",
                          value = 0.05, min = 0.005, max = 0.5,
                          step = 0.005)
    })
    output$reduction_ui <- shiny::renderUI({
      ds <- state$active_dataset; if (is.null(ds)) return(NULL)
      shiny::selectInput(ns("reduction"), "Reduction",
                         choices  = available_reductions(ds),
                         selected = state$selected_reduction)
    })

    # ---- Input validation --------------------------------------------
    output$input_warning <- shiny::renderUI({
      if (is.null(state$active_dataset))
        return(friendly_warning("No dataset loaded."))
      src <- get_regulon_source(input$source %||% "")
      if (!is.null(src) && !has_optional(src$requires)) {
        return(friendly_warning(sprintf(
          "Source '%s' requires: %s. Click Run for the install hint.",
          src$name, paste(src$requires, collapse = ", "))))
      }
      eng <- get_regulon_engine(input$engine %||% "")
      if (!is.null(eng) && !regulon_engine_available(eng)) {
        return(friendly_warning(sprintf(
          "Engine '%s' requires: %s. Pick `aucell_r` (no deps).",
          eng$name, paste(eng$requires, collapse = ", "))))
      }
      NULL
    })

    # ---- Run ----------------------------------------------------------
    shiny::observeEvent(input$run, {
      ds <- state$active_dataset; if (is.null(ds)) return()
      src_id <- input$source %||% "mock_pbmc"
      eng_id <- input$engine %||% "aucell_r"
      params <- list(top_n_fraction = input$top_n_fraction %||% 0.05)
      set_rg_slot(list(status = "running", results = NULL,
                       params = list(source = src_id, engine = eng_id,
                                     engine_params = params),
                       error_message = NULL, timestamp = Sys.time(),
                       duration_ms = NULL,
                       annotation_stamp = make_annotation_stamp(state)))
      t0 <- proc.time()[["elapsed"]]
      out <- tryCatch(shiny::withProgress(
        message = "Scoring regulons...", value = 0.3, {
          set <- fetch_regulon_set(src_id)
          shiny::incProgress(0.3, detail = "AUCell scoring")
          run_regulon_engine(eng_id, ds, set, params = params)
        }), error = function(e) e)
      if (inherits(out, "error")) {
        set_rg_slot(list(status = "failed", results = NULL,
                         params = list(source = src_id, engine = eng_id),
                         error_message = conditionMessage(out),
                         timestamp = Sys.time(), duration_ms = NULL,
                         annotation_stamp = make_annotation_stamp(state)))
        push_message(state, sprintf("Regulons failed: %s",
                                    conditionMessage(out)), "error")
      } else {
        dur <- as.integer(round((proc.time()[["elapsed"]] - t0) * 1000))
        set_rg_slot(list(status = "completed", results = out,
                         params = list(source = src_id, engine = eng_id,
                                       engine_params = params),
                         error_message = NULL,
                         timestamp = Sys.time(), duration_ms = dur,
                         annotation_stamp = make_annotation_stamp(state)))
        push_message(state, sprintf(
          "Scored %d regulons across %d cells (%d ms).",
          length(out$regulon_ids), length(out$cell_ids), dur), "success")
      }
    })

    # ---- Status banner ------------------------------------------------
    output$status_banner <- shiny::renderUI({
      tr <- rg_slot()
      if (is.null(tr))
        return(shiny::div(style = "padding:8px 12px; background:#eee; border-radius:4px; font-size:13px;",
                          shiny::tags$strong("Status: "),
                          "Not run yet. Pick a source + engine and click ",
                          shiny::tags$em("Score regulons"), "."))
      bg <- switch(tr$status, running = "#cfe2ff", completed = "#d1e7dd",
                   failed = "#f8d7da", "#eee")
      fg <- switch(tr$status, running = "#084298", completed = "#0a3622",
                   failed = "#842029", "#333")
      shiny::div(style = sprintf(
        "padding:8px 12px; background:%s; color:%s; border-radius:4px; font-size:13px;",
        bg, fg),
        shiny::tags$strong("Status: "),
        switch(tr$status,
          running   = "Running...",
          completed = sprintf(
            "Completed. Source=%s, engine=%s, %d regulons, %d cells (%d ms).",
            tr$params$source %||% "?",  tr$params$engine %||% "?",
            length(tr$results$regulon_ids %||% character()),
            length(tr$results$cell_ids %||% character()),
            tr$duration_ms %||% 0L),
          failed    = sprintf("Failed: %s", tr$error_message %||% "?"),
          tr$status))
    })

    # ---- Heatmap ------------------------------------------------------
    output$heatmap_warning <- shiny::renderUI({
      tr <- rg_slot()
      if (is.null(tr) || !identical(tr$status, "completed"))
        return(friendly_warning("Run the engine to populate the heatmap."))
      NULL
    })
    output$heatmap <- shiny::renderPlot({
      tr <- rg_slot(); if (is.null(tr) || !identical(tr$status, "completed")) return()
      ds <- state$active_dataset; if (is.null(ds)) return()
      gf <- input$group_field %||% "cluster"
      grp <- get_metadata(ds, gf); if (is.null(grp)) return()
      m <- regulon_mean_by_group(tr$results, grp)
      plot_regulon_heatmap(m, title = sprintf("Mean AUC by %s", gf))
    })

    # ---- Per-regulon embedding ---------------------------------------
    output$regulon_picker_ui <- shiny::renderUI({
      tr <- rg_slot()
      if (is.null(tr) || !identical(tr$status, "completed")) return(NULL)
      rids <- tr$results$regulon_ids %||% character()
      shiny::selectInput(ns("regulon"), "Color by regulon",
                         choices = rids,
                         selected = if (length(rids)) rids[1] else NULL)
    })
    output$emb_title <- shiny::renderText({
      reg <- input$regulon %||% "(none)"
      sprintf("Embedding by AUC: %s", reg)
    })
    output$emb_warning <- shiny::renderUI({
      tr <- rg_slot()
      if (is.null(tr) || !identical(tr$status, "completed"))
        return(friendly_warning("Run the engine to populate this panel."))
      if (is.null(input$regulon) || !nzchar(input$regulon))
        return(friendly_warning("Pick a regulon to color by."))
      NULL
    })
    output$emb_plot <- shiny::renderPlot({
      tr <- rg_slot(); if (is.null(tr) || !identical(tr$status, "completed")) return()
      ds <- state$active_dataset; if (is.null(ds)) return()
      reg <- input$regulon; if (is.null(reg) || !nzchar(reg)) return()
      red <- input$reduction %||% state$selected_reduction
      emb <- get_embedding(ds, red); if (is.null(emb)) return()
      auc <- tr$results$auc_matrix[, reg, drop = TRUE]
      plot_embedding_continuous(emb, auc,
        title = sprintf("%s | AUC", red),
        legend_title = "AUC")
    })
  })
}
