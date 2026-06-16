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
#
# UI composition uses the shared primitives in `R/ui_components.R`
# (`page_header`, `control_panel`, `plot_card`, `info_banner`, etc.).
# Server logic (reactives, observers, plot rendering) is unchanged from v1.
# ============================================================================

mod_scrna_explorer_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    page_header(
      eyebrow = "Exploration",
      title   = "Basic scRNA Explorer",
      lede    = paste("Inspect reductions and color cells by metadata or",
                      "by gene. Selections are shared with every other",
                      "module.")
    ),

    # -- Summary band --------------------------------------------------------
    shiny::uiOutput(ns("summary_band")),

    # -- Controls (visually distinct panel) ----------------------------------
    control_panel(
      title = "View settings",
      shiny::fluidRow(
        shiny::column(3, shiny::uiOutput(ns("assay_ui"))),
        shiny::column(3, shiny::uiOutput(ns("reduction_ui"))),
        shiny::column(3, shiny::uiOutput(ns("metadata_ui"))),
        shiny::column(3, shiny::uiOutput(ns("gene_ui")))
      )
    ),

    # -- Two plots side by side ---------------------------------------------
    shiny::fluidRow(
      shiny::column(6,
        plot_card(
          title    = shiny::textOutput(ns("meta_title"), inline = TRUE),
          caption  = "metadata embedding",
          footnote = "Drag a box on the plot to select cells.",
          shiny::uiOutput(ns("meta_warning")),
          shiny::div(class = "plot-container",
            shiny::plotOutput(
              ns("meta_plot"), height = "440px",
              brush = shiny::brushOpts(id = ns("meta_brush"),
                                       resetOnNew = TRUE)
            )
          )
        )
      ),
      shiny::column(6,
        plot_card(
          title   = shiny::textOutput(ns("gene_title"), inline = TRUE),
          caption = "feature embedding",
          shiny::uiOutput(ns("display_mode_ui")),
          shiny::uiOutput(ns("gene_warning")),
          shiny::div(class = "plot-container",
            shiny::plotOutput(ns("gene_plot"), height = "440px")
          )
        )
      )
    ),

    # -- Selection footer ----------------------------------------------------
    app_card(
      title   = "Selection",
      caption = "shared with downstream modules",
      shiny::fluidRow(
        shiny::column(8,
          shiny::div(
            class = "selection-meta",
            content_label("Selected cells"),
            shiny::span(class = "sce-tabular",
                        shiny::textOutput(ns("n_selected"), inline = TRUE)),
            helper_text("brush the left plot to select; click Clear to reset")
          )
        ),
        shiny::column(4, align = "right",
          shiny::actionButton(ns("clear_selection"), "Clear selection",
                              class = "btn btn-default btn-sm")
        )
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
    # Use the shared server-side gene picker so the 32k-gene PBMC 8k
    # dataset (and larger) doesn't ship its full vocabulary to the
    # browser. See `R/ui_components.R::gene_picker_input()` for the
    # rationale; identical pattern in mod_trajectory / mod_pathway_analysis
    # / mod_marker_investigation / mod_imputation.
    output$gene_ui <- shiny::renderUI({
      gene_picker_input(ns("gene"), "Gene (FeaturePlot)",
                        selected = state$selected_gene)
    })

    # Fill the picker's choices ONCE per dataset. Reading
    # `state$selected_gene` via `isolate()` so this observer does NOT
    # re-fire (and re-transmit the full 32k-gene choice list) every
    # time another module updates the shared gene. The lighter
    # selection-sync observer below handles selection drift.
    shiny::observe({
      ds <- state$active_dataset
      shiny::req(ds)
      update_gene_picker(session, "gene",
                         choices  = available_genes(ds),
                         selected = shiny::isolate(state$selected_gene))
    })

    # Lightweight selection sync. Fires when `state$selected_gene`
    # changes from another module (Marker Investigation push, etc.).
    # Two important guards:
    #   - `nzchar(g)` so a transient empty value from a selectize
    #     update never propagates here.
    #   - `identical(input$gene, g)` no-op: if the picker is already
    #     showing this gene, don't send another updateSelectizeInput.
    #     Without this the observer self-loops with the
    #     `observeEvent(input$gene, ...)` below every time the user
    #     picks a gene, which is what was causing the FeaturePlot to
    #     "reset almost instantly" -- the cycle's update messages
    #     race with the user's pick and one of them clears the
    #     selection.
    shiny::observe({
      g <- state$selected_gene
      if (is.null(g) || !nzchar(g)) return()
      cur <- shiny::isolate(input$gene)
      if (isTRUE(identical(cur, g))) return()
      shiny::updateSelectizeInput(session, "gene", selected = g)
    })

    # Push local inputs back into shared state. The `input$gene`
    # handler is guarded against empty values because
    # `updateSelectizeInput(server = TRUE)` can briefly fire a
    # `change` event with an empty selection while it swaps choices.
    # Without the guard, that transient `""` would overwrite
    # `state$selected_gene` and immediately clear the FeaturePlot.
    shiny::observeEvent(input$assay,     state$selected_assay          <- input$assay,     ignoreInit = TRUE)
    shiny::observeEvent(input$reduction, state$selected_reduction      <- input$reduction, ignoreInit = TRUE)
    shiny::observeEvent(input$metadata,  state$selected_metadata_field <- input$metadata,  ignoreInit = TRUE)
    shiny::observeEvent(input$gene, {
      if (!is.null(input$gene) && nzchar(input$gene)) {
        state$selected_gene <- input$gene
      }
    }, ignoreInit = TRUE)

    # ---- Summary band ----------------------------------------------------
    # Content-side persistent summary bar (distinct from the workspace's
    # annotation context strip). Lists the dataset + the columns each module
    # downstream is reading from. Uses the shared `summary_bar()` primitive
    # so it doesn't depend on page-header internals.
    output$summary_band <- shiny::renderUI({
      ds <- state$active_dataset
      if (is.null(ds)) return(NULL)
      summary_bar(shiny::tagList(
        summary_item("Dataset",    ds$name),
        summary_item("Cells",      format(ds$n_cells, big.mark = ",")),
        summary_item("Genes",      format(ds$n_genes, big.mark = ",")),
        summary_item("Assays",     paste(ds$assays,     collapse = ", ")),
        summary_item("Reductions", paste(ds$reductions, collapse = ", "))
      ))
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
        return(info_banner(
          tone = "neutral",
          "Run Data Smoothing to enable a smoothed view of this plot."
        ))
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
