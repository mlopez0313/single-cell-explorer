# ============================================================================
# Module: Data Smoothing / Imputation  (ENABLED)
# ----------------------------------------------------------------------------
# Compute smoothed expression vectors for visualization-only exploration.
# Raw expression stays the analytic source of truth for DE / Markers /
# Pathway -- they call `get_gene_expression()` directly. Visualization
# modules (the Basic scRNA Explorer) read through
# `get_gene_expression_for_view()`, which switches between raw and smoothed
# based on `state$display_mode_imputation`.
#
# Results land in `state$analysis_results$imputation`. Clearing wipes that
# slot AND resets the display mode to "raw" so no module is left in a
# dangling state.
# ============================================================================

mod_imputation_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h2("Data Smoothing / Imputation"),

    # -- Visualization-only warning (prominent) ---------------------------
    shiny::div(
      style = paste(
        "padding:10px 14px; background:#fff3cd; color:#664d03;",
        "border:1px solid #ffecb5; border-radius:4px; font-size:13px;",
        "margin-bottom:12px;"),
      shiny::tags$strong("Visualization only."),
      " Smoothed values are intended for exploration. ",
      "Marker Investigation, Differential Expression, and Pathway Analysis ",
      "continue to use ", shiny::tags$strong("raw"), " expression by design."
    ),

    # -- Controls ---------------------------------------------------------
    shiny::fluidRow(
      shiny::column(3, shiny::uiOutput(ns("assay_ui"))),
      shiny::column(3, shiny::selectInput(ns("method"), "Method",
                                          choices  = available_imputation_methods(),
                                          selected = "neighbor")),
      shiny::column(3, shiny::numericInput(ns("k"), "k (neighborhood)",
                                           value = 15, min = 1, step = 1)),
      shiny::column(3, shiny::uiOutput(ns("target_genes_ui")))
    ),

    shiny::div(style = "margin: 8px 0 16px 0;",
      shiny::actionButton(ns("run"), "Run Smoothing",
                          class = "btn btn-primary"),
      shiny::actionButton(ns("clear"), "Clear Smoothed Data",
                          class = "btn btn-default", style = "margin-left:8px;"),
      shiny::tags$span(style = "margin-left:12px; color:#888; font-size:12px;",
                       "Smoothing runs only when ", shiny::tags$em("Run Smoothing"), " is clicked.")
    ),

    shiny::uiOutput(ns("status_banner")),
    shiny::uiOutput(ns("input_warning")),

    shiny::hr(),

    # -- Side-by-side raw vs smoothed FeaturePlot -------------------------
    shiny::fluidRow(
      shiny::column(12,
        shiny::uiOutput(ns("gene_picker_ui"))
      )
    ),
    shiny::fluidRow(
      shiny::column(6,
        shiny::h4("Raw"),
        shiny::uiOutput(ns("raw_warning")),
        shiny::plotOutput(ns("raw_plot"), height = "400px")
      ),
      shiny::column(6,
        shiny::h4("Smoothed"),
        shiny::uiOutput(ns("smoothed_warning")),
        shiny::plotOutput(ns("smoothed_plot"), height = "400px")
      )
    )
  )
}

mod_imputation_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- Controls --------------------------------------------------------
    output$assay_ui <- shiny::renderUI({
      shiny::selectInput(ns("assay"), "Source assay",
                         choices  = available_assays(state$active_dataset),
                         selected = state$selected_assay)
    })
    output$target_genes_ui <- shiny::renderUI({
      genes <- available_genes(state$active_dataset)
      shiny::selectizeInput(ns("genes"), "Target genes",
                            choices  = genes,
                            selected = genes,
                            multiple = TRUE,
                            options  = list(placeholder = "blank = all genes"))
    })

    # ---- Slot helpers ----------------------------------------------------
    imp_slot     <- function() state$analysis_results$imputation
    set_imp_slot <- function(x) {
      ar <- state$analysis_results
      ar$imputation <- x
      state$analysis_results <- ar
    }

    # ---- Input warnings --------------------------------------------------
    output$input_warning <- shiny::renderUI({
      ds <- state$active_dataset
      if (is.null(ds)) return(NULL)  # workspace empty state handles this
      if (length(available_genes(ds)) == 0L)
        return(friendly_warning("Dataset has no genes; smoothing is not possible."))
      NULL
    })

    # ---- Run -------------------------------------------------------------
    shiny::observeEvent(input$run, {
      ds <- state$active_dataset; if (is.null(ds)) return()
      method <- input$method %||% "neighbor"
      requested <- input$genes
      if (is.null(requested) || length(requested) == 0L) requested <- available_genes(ds)

      params <- list(assay = input$assay, method = method,
                     k = input$k %||% 15L, genes = requested)
      set_imp_slot(list(status = "running", results = NULL,
                        params = params, error_message = NULL,
                        timestamp = Sys.time(), duration_ms = NULL,
                        annotation_stamp = make_annotation_stamp(state)))

      t0 <- proc.time()[["elapsed"]]
      out <- tryCatch(
        shiny::withProgress(message = sprintf("Smoothing (%s)...", method), value = 0.4, {
          compute_smoothed(ds, genes = requested,
                           method = method, k = params$k)
        }),
        error = function(e) e
      )
      if (inherits(out, "error")) {
        set_imp_slot(list(status = "failed", results = NULL,
                          params = params, error_message = conditionMessage(out),
                          timestamp = Sys.time(), duration_ms = NULL,
                          annotation_stamp = make_annotation_stamp(state)))
        push_message(state, sprintf("Smoothing failed: %s", conditionMessage(out)), "error")
      } else {
        dur <- as.integer(round((proc.time()[["elapsed"]] - t0) * 1000))
        set_imp_slot(list(status = "completed", results = out,
                          params = params, error_message = NULL,
                          timestamp = Sys.time(), duration_ms = dur,
                          annotation_stamp = make_annotation_stamp(state)))
        push_message(state, sprintf(
          "Smoothing done: %s on %d genes (k=%d, %d ms). Visualization-only.",
          out$method, length(out$genes),
          out$k %||% NA_integer_, dur), "success")
      }
    })

    # ---- Clear -----------------------------------------------------------
    shiny::observeEvent(input$clear, {
      set_imp_slot(NULL)
      # Reset display mode so the Explorer doesn't try to render smoothed
      # data that no longer exists.
      state$display_mode_imputation <- "raw"
      push_message(state, "Cleared smoothed data; Explorer reverted to raw.", "info")
    })

    # ---- Status banner ---------------------------------------------------
    output$status_banner <- shiny::renderUI({
      imp <- imp_slot()
      if (is.null(imp))
        return(shiny::div(style = "padding:8px 12px; background:#eee; border-radius:4px; font-size:13px;",
                          shiny::tags$strong("Status: "), "Not run yet. Configure controls and click ",
                          shiny::tags$em("Run Smoothing"), "."))
      bg <- switch(imp$status, running = "#cfe2ff", completed = "#d1e7dd",
                   failed = "#f8d7da", "#eee")
      fg <- switch(imp$status, running = "#084298", completed = "#0a3622",
                   failed = "#842029", "#333")
      txt <- switch(imp$status,
        running   = "Running...",
        completed = sprintf("Completed: %s on %d genes (k=%s, reduction=%s, %d ms).",
                            imp$results$method, length(imp$results$genes),
                            imp$results$k %||% "NA",
                            imp$results$reduction_used %||% "NA",
                            imp$duration_ms %||% 0L),
        failed    = sprintf("Failed: %s", imp$error_message %||% "(unknown)"),
        "")
      shiny::div(style = sprintf("padding:8px 12px; background:%s; color:%s; border-radius:4px; font-size:13px;",
                                 bg, fg),
                 shiny::tags$strong("Status: "), txt)
    })

    # ---- Gene picker (synced with shared state$selected_gene) ----------
    output$gene_picker_ui <- shiny::renderUI({
      ds <- state$active_dataset; if (is.null(ds)) return(NULL)
      ds_genes  <- available_genes(ds)
      sm_genes  <- if (has_smoothed_results(state)) imp_slot()$results$genes else character()
      # Smoothed genes first, then any others in the dataset.
      choices <- unique(c(sm_genes, ds_genes))
      current <- state$selected_gene %||% choices[1]
      if (!current %in% choices) current <- choices[1]
      shiny::selectizeInput(ns("gene"), "Gene (visualised in both panes)",
                            choices = choices, selected = current)
    })
    shiny::observeEvent(input$gene, {
      if (!is.null(input$gene) && nzchar(input$gene)) {
        state$selected_gene <- input$gene
      }
    }, ignoreInit = TRUE)

    # ---- Raw FeaturePlot -----------------------------------------------
    output$raw_warning <- shiny::renderUI({
      ds <- state$active_dataset
      g  <- state$selected_gene
      if (is.null(ds))           return(friendly_warning("No dataset loaded."))
      if (!validate_gene(ds, g)) return(friendly_warning(sprintf("Gene '%s' is not in the dataset.", g %||% "")))
      if (is.null(get_embedding(ds, state$selected_reduction %||% ds$default_reduction)))
        return(friendly_warning("No usable reduction."))
      NULL
    })
    output$raw_plot <- shiny::renderPlot({
      ds <- state$active_dataset
      shiny::req(ds, state$selected_gene)
      red <- state$selected_reduction %||% ds$default_reduction
      emb <- get_embedding(ds, red)
      expr <- get_gene_expression(ds, state$selected_gene)
      if (is.null(emb) || is.null(expr)) return(NULL)
      plot_embedding_continuous(emb, expr,
                                title = sprintf("Raw | %s | %s", red, state$selected_gene),
                                xlab  = paste0(red, "_1"), ylab = paste0(red, "_2"),
                                legend_title = state$selected_gene)
    })

    # ---- Smoothed FeaturePlot ------------------------------------------
    output$smoothed_warning <- shiny::renderUI({
      imp <- imp_slot()
      if (is.null(imp) || !identical(imp$status, "completed"))
        return(shiny::div(style = "padding:8px 12px; font-size:13px; color:#666;",
                          "Run smoothing to populate this panel."))
      g <- state$selected_gene
      if (is.null(imp$results$expression[[g]]))
        return(friendly_warning(sprintf("Gene '%s' was not included in the smoothing run.", g %||% "")))
      NULL
    })
    output$smoothed_plot <- shiny::renderPlot({
      imp <- imp_slot()
      ds  <- state$active_dataset
      shiny::req(imp, state$selected_gene)
      if (!identical(imp$status, "completed")) return(NULL)
      v <- imp$results$expression[[state$selected_gene]]
      if (is.null(v)) return(NULL)
      red <- imp$results$reduction_used %||% state$selected_reduction %||% ds$default_reduction
      emb <- get_embedding(ds, red)
      if (is.null(emb)) return(NULL)
      plot_embedding_continuous(emb, v,
                                title = sprintf("Smoothed (%s, k=%s) | %s",
                                                imp$results$method,
                                                imp$results$k %||% "NA",
                                                state$selected_gene),
                                xlab = paste0(red, "_1"), ylab = paste0(red, "_2"),
                                legend_title = state$selected_gene)
    })
  })
}
