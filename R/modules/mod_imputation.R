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
#
# UI composition uses the shared primitives in `R/ui_components.R`. Server
# logic (reactives, observers, run/clear handlers) is unchanged.
# ============================================================================

mod_imputation_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    page_header(
      eyebrow = "Advanced Analysis",
      title   = "Data Smoothing / Imputation",
      lede    = paste("Compute smoothed expression for visualization-only",
                      "exploration. The smoothed FeaturePlot in the Basic",
                      "Explorer reads from this module's result; DE,",
                      "Markers, and Pathway continue to use raw values.")
    ),

    info_banner(
      tone  = "warning",
      title = "Visualization only.",
      "Smoothed values are intended for exploration. ",
      "Marker Investigation, Differential Expression, and Pathway Analysis ",
      "continue to use ", shiny::tags$strong("raw"),
      " expression by design."
    ),

    control_panel(
      title = "Smoothing settings",
      shiny::fluidRow(
        shiny::column(3, shiny::uiOutput(ns("assay_ui"))),
        shiny::column(3, shiny::selectInput(ns("method"), "Method",
                                            choices  = available_imputation_methods(),
                                            selected = "neighbor")),
        shiny::column(3, shiny::numericInput(ns("k"), "k (neighborhood)",
                                             value = 15, min = 1, step = 1)),
        shiny::column(3, shiny::uiOutput(ns("target_genes_ui")))
      ),
      actions = shiny::tagList(
        shiny::actionButton(ns("run"),   "Run smoothing",
                            class = "btn btn-primary"),
        shiny::actionButton(ns("clear"), "Clear smoothed data",
                            class = "btn btn-default"),
        helper_text(
          "Smoothing runs only when ", shiny::tags$em("Run smoothing"),
          " is clicked. ", "Clear resets the Explorer to the raw view.")
      )
    ),

    shiny::uiOutput(ns("status_banner")),
    shiny::uiOutput(ns("input_warning")),

    app_card(
      title   = "Gene picker",
      caption = "visualised in both panes",
      shiny::uiOutput(ns("gene_picker_ui"))
    ),

    shiny::fluidRow(
      shiny::column(6,
        plot_card(
          title = "Raw expression",
          caption = "from `get_gene_expression()`",
          shiny::uiOutput(ns("raw_warning")),
          shiny::div(class = "plot-container",
            shiny::plotOutput(ns("raw_plot"), height = "400px"))
        )
      ),
      shiny::column(6,
        plot_card(
          title = "Smoothed expression",
          caption = "from this module's result",
          shiny::uiOutput(ns("smoothed_warning")),
          shiny::div(class = "plot-container",
            shiny::plotOutput(ns("smoothed_plot"), height = "400px"))
        )
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
        return(status_banner(
          shiny::span("Not run yet. Configure controls and click ",
                      shiny::tags$em("Run smoothing"), "."),
          tone = "idle"))
      tone <- switch(imp$status, running = "running",
                     completed = "success", failed = "danger", "idle")
      text <- switch(imp$status,
        running   = "Running...",
        completed = sprintf("Completed: %s on %d genes (k=%s, reduction=%s, %d ms).",
                            imp$results$method, length(imp$results$genes),
                            imp$results$k %||% "NA",
                            imp$results$reduction_used %||% "NA",
                            imp$duration_ms %||% 0L),
        failed    = sprintf("Failed: %s", imp$error_message %||% "(unknown)"),
        "")
      status_banner(text, tone = tone)
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
        return(helper_text("Run smoothing to populate this panel."))
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
