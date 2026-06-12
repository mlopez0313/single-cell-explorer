# ============================================================================
# Module: Basic scRNA Explorer  (ENABLED, v1)
# ----------------------------------------------------------------------------
# Lets the user inspect a dataset:
#   - pick assay / reduction / metadata field / gene
#   - see a metadata-colored embedding (DimPlot style)
#   - see a gene-colored embedding         (FeaturePlot style)
#   - select cells by brushing the metadata embedding
#
# All cross-module selections are written back to shared app state. Dataset
# access goes through R/dataset_helpers.R; plotting through R/plotting.R.
# ============================================================================

mod_scrna_explorer_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h2("Basic scRNA Explorer"),
    shiny::p("Explore the active dataset. Your selections are shared with other modules."),

    # -- Controls ------------------------------------------------------------
    shiny::fluidRow(
      shiny::column(3, shiny::uiOutput(ns("assay_ui"))),
      shiny::column(3, shiny::uiOutput(ns("reduction_ui"))),
      shiny::column(3, shiny::uiOutput(ns("metadata_ui"))),
      shiny::column(3, shiny::uiOutput(ns("gene_ui")))
    ),

    shiny::hr(),

    # -- Summary band --------------------------------------------------------
    shiny::uiOutput(ns("summary_band")),

    # -- Two plots side by side ---------------------------------------------
    shiny::fluidRow(
      shiny::column(6,
        shiny::h4(shiny::textOutput(ns("meta_title"), inline = TRUE)),
        shiny::uiOutput(ns("meta_warning")),
        shiny::plotOutput(
          ns("meta_plot"), height = "440px",
          brush = shiny::brushOpts(id = ns("meta_brush"), resetOnNew = TRUE)
        ),
        shiny::div(style = "font-size:12px; color:#888;",
                   "Tip: drag a box on the plot to select cells.")
      ),
      shiny::column(6,
        shiny::h4(shiny::textOutput(ns("gene_title"), inline = TRUE)),
        shiny::uiOutput(ns("display_mode_ui")),
        shiny::uiOutput(ns("gene_warning")),
        shiny::plotOutput(ns("gene_plot"), height = "440px")
      )
    ),

    shiny::hr(),

    # -- Selection footer ----------------------------------------------------
    shiny::fluidRow(
      shiny::column(8,
        shiny::strong("Selected cells: "),
        shiny::textOutput(ns("n_selected"), inline = TRUE),
        shiny::tags$span(style = "margin-left:12px; color:#888; font-size:12px;",
                         "(brush the left plot to select; click \"Clear\" to reset)")
      ),
      shiny::column(4, align = "right",
        shiny::actionButton(ns("clear_selection"), "Clear selection",
                            class = "btn btn-default btn-sm")
      )
    )
  )
}

mod_scrna_explorer_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- Controls --------------------------------------------------------
    output$assay_ui <- shiny::renderUI({
      shiny::selectInput(ns("assay"), "Assay",
                         choices  = available_assays(state$active_dataset),
                         selected = state$selected_assay)
    })
    output$reduction_ui <- shiny::renderUI({
      shiny::selectInput(ns("reduction"), "Reduction",
                         choices  = available_reductions(state$active_dataset),
                         selected = state$selected_reduction)
    })
    output$metadata_ui <- shiny::renderUI({
      shiny::selectInput(ns("metadata"), "Color by (metadata)",
                         choices  = available_metadata_fields(state$active_dataset),
                         selected = state$selected_metadata_field)
    })
    output$gene_ui <- shiny::renderUI({
      shiny::selectizeInput(ns("gene"), "Gene (FeaturePlot)",
                            choices  = available_genes(state$active_dataset),
                            selected = state$selected_gene,
                            options  = list(placeholder = "Pick a gene..."))
    })

    # Push local inputs back into shared state.
    shiny::observeEvent(input$assay,     state$selected_assay          <- input$assay,     ignoreInit = TRUE)
    shiny::observeEvent(input$reduction, state$selected_reduction      <- input$reduction, ignoreInit = TRUE)
    shiny::observeEvent(input$metadata,  state$selected_metadata_field <- input$metadata,  ignoreInit = TRUE)
    shiny::observeEvent(input$gene,      state$selected_gene           <- input$gene,      ignoreInit = TRUE)

    # ---- Summary band ----------------------------------------------------
    output$summary_band <- shiny::renderUI({
      ds <- state$active_dataset
      shiny::div(
        style = "padding:8px 12px; background:#f5f7fa; border-radius:4px; font-size:13px; margin-bottom:12px;",
        shiny::strong(ds$name), " | ",
        format(ds$n_cells, big.mark = ","), " cells | ",
        format(ds$n_genes, big.mark = ","), " genes | ",
        "assays: ",     paste(ds$assays,     collapse = ", "), " | ",
        "reductions: ", paste(ds$reductions, collapse = ", ")
      )
    })

    # ---- Embedding ------------------------------------------------------
    # Single source for the active embedding -- both plots reuse it.
    embedding <- shiny::reactive({
      get_embedding(state$active_dataset, state$selected_reduction)
    })

    # ---- Metadata plot --------------------------------------------------
    output$meta_title <- shiny::renderText({
      sprintf("Embedding colored by %s",
              state$selected_metadata_field %||% "(none)")
    })

    output$meta_warning <- shiny::renderUI({
      if (is.null(embedding()))
        return(friendly_warning(sprintf("Reduction '%s' is not available in this dataset.",
                                        state$selected_reduction %||% "")))
      if (!validate_metadata(state$active_dataset, state$selected_metadata_field))
        return(friendly_warning(sprintf("Metadata field '%s' is not available in this dataset.",
                                        state$selected_metadata_field %||% "")))
      NULL
    })

    output$meta_plot <- shiny::renderPlot({
      emb <- embedding(); if (is.null(emb)) return(NULL)
      meta <- get_metadata(state$active_dataset, state$selected_metadata_field)
      if (is.null(meta)) return(NULL)
      plot_embedding_categorical(
        emb, meta,
        title = sprintf("%s | colored by %s",
                        state$selected_reduction, state$selected_metadata_field),
        xlab  = paste0(state$selected_reduction, "_1"),
        ylab  = paste0(state$selected_reduction, "_2")
      )
    })

    # ---- Gene FeaturePlot ----------------------------------------------
    output$gene_title <- shiny::renderText({
      sprintf("FeaturePlot: %s", state$selected_gene %||% "(none)")
    })

    output$gene_warning <- shiny::renderUI({
      if (is.null(embedding()))
        return(friendly_warning(sprintf("Reduction '%s' is not available in this dataset.",
                                        state$selected_reduction %||% "")))
      if (!validate_gene(state$active_dataset, state$selected_gene))
        return(friendly_warning(sprintf("Gene '%s' is not available in this dataset.",
                                        state$selected_gene %||% "")))
      NULL
    })

    # Display-mode toggle: appears only when smoothed data exists. Switches
    # the FeaturePlot's data source between raw and smoothed via the helper
    # `get_gene_expression_for_view()`. This affects ONLY this plot --
    # Marker Investigation, DE, and Pathway never read smoothed values.
    output$display_mode_ui <- shiny::renderUI({
      if (!has_smoothed_results(state)) {
        return(shiny::div(style = "font-size:11px; color:#aaa; margin:4px 0 8px 0;",
          "Tip: run Data Smoothing to enable a smoothed view."))
      }
      shiny::radioButtons(ns("display_mode"), label = NULL,
        choices = c("Raw" = "raw", "Smoothed (visualization only)" = "smoothed"),
        selected = state$display_mode_imputation %||% "raw", inline = TRUE)
    })
    shiny::observeEvent(input$display_mode, {
      state$display_mode_imputation <- input$display_mode
    }, ignoreInit = TRUE)

    output$gene_plot <- shiny::renderPlot({
      emb <- embedding(); if (is.null(emb)) return(NULL)
      expr <- get_gene_expression_for_view(state, state$selected_gene)
      if (is.null(expr)) return(NULL)
      mode <- state$display_mode_imputation %||% "raw"
      tag  <- if (identical(mode, "smoothed") && has_smoothed_results(state)) "smoothed" else "raw"
      plot_embedding_continuous(
        emb, expr,
        title = sprintf("%s | %s expression  [%s]",
                        state$selected_reduction, state$selected_gene, tag),
        xlab  = paste0(state$selected_reduction, "_1"),
        ylab  = paste0(state$selected_reduction, "_2"),
        legend_title = state$selected_gene
      )
    })

    # ---- Brush selection -> shared state -------------------------------
    shiny::observeEvent(input$meta_brush, {
      emb <- embedding(); if (is.null(emb)) return()
      b   <- input$meta_brush
      if (is.null(b)) return()
      hit <- emb$x >= b$xmin & emb$x <= b$xmax &
             emb$y >= b$ymin & emb$y <= b$ymax
      state$selected_cells <- emb$cell[hit]
    }, ignoreInit = TRUE, ignoreNULL = FALSE)

    shiny::observeEvent(input$clear_selection, {
      state$selected_cells <- character()
    })

    output$n_selected <- shiny::renderText({
      n <- length(state$selected_cells)
      if (n == 0) "none" else format(n, big.mark = ",")
    })
  })
}

# `friendly_warning()` lives in R/ui_helpers.R (shared across modules).
