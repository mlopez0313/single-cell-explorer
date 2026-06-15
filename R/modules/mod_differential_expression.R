# ============================================================================
# Module: Differential Expression  (ENABLED)
# ----------------------------------------------------------------------------
# Compare two specific groups of cells, produce a DE results table, and let
# the user explore individual genes (FeaturePlot + per-group expression).
#
# Design notes:
#   - Compute logic lives in R/de.R (`compute_de()`); this module only does
#     UI and reactive plumbing.
#   - Results are stored in `state$analysis_results$de` so downstream modules
#     (Pathway Analysis, etc.) can read them without re-running the test.
#   - DE differs from Marker Investigation: marker = "group vs rest", DE =
#     "group_1 vs group_2". The result schema (group_1, group_2, p_val_adj)
#     is therefore distinct.
#   - "Click a gene" works via the volcano plot's click input + nearPoints().
#     Clicking sets `state$selected_gene`, which other modules
#     (Basic scRNA Explorer) react to.
# ============================================================================

mod_differential_expression_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    page_header(
      eyebrow = "Statistics",
      title   = "Differential Expression",
      lede    = shiny::tagList(
        "Compare two groups of cells. Results live in ",
        shiny::tags$code("state$analysis_results$de"),
        " for downstream modules to read without re-running.")
    ),

    # -- Compute controls -------------------------------------------------
    control_panel(
      title = "Test design",
      shiny::fluidRow(
        shiny::column(3, shiny::uiOutput(ns("group_field_ui"))),
        shiny::column(3, shiny::uiOutput(ns("group1_ui"))),
        shiny::column(3, shiny::uiOutput(ns("group2_ui"))),
        shiny::column(3, shiny::uiOutput(ns("assay_ui")))
      ),
      shiny::fluidRow(
        shiny::column(3, shiny::selectInput(ns("layer"), "Layer / slot",
                                            choices  = c("data", "counts", "scale.data"),
                                            selected = "data")),
        shiny::column(3, shiny::numericInput(ns("min_pct"), "min % expressed",
                                             value = 0.1, min = 0, max = 1, step = 0.05)),
        shiny::column(3, shiny::numericInput(ns("log2fc_thr"), "Significance |log2FC| \u2265",
                                             value = 0.5, min = 0, step = 0.1)),
        shiny::column(3, shiny::numericInput(ns("padj_thr"), "Significance adj. p \u2264",
                                             value = 0.05, min = 0, max = 1, step = 0.01))
      ),
      shiny::fluidRow(
        shiny::column(6, shiny::uiOutput(ns("backend_ui"))),
        shiny::column(6,
          beside_input(
            helper_text(
              "Auto-select uses ", shiny::tags$code("presto::wilcoxauc()"),
              " when installed; otherwise falls back to the pure-R Wilcoxon. ",
              "Pseudobulk backends require a `sample_by` replicate column.")))
      ),
      # ---- Pseudobulk-only controls (revealed when a pseudobulk backend is picked)
      shiny::conditionalPanel(
        condition = sprintf(
          "['pseudobulk_naive','pseudobulk_edger','pseudobulk_deseq2'].indexOf(input['%s']) >= 0",
          ns("backend")),
        shiny::fluidRow(
          shiny::column(4, shiny::uiOutput(ns("sample_by_ui"))),
          shiny::column(4, shiny::numericInput(ns("min_cells_per_sample"),
                                               "Min cells per pseudobulk sample",
                                               value = 10, min = 1, step = 1)),
          shiny::column(4, shiny::numericInput(ns("min_samples_per_group"),
                                               "Min samples per group",
                                               value = 2, min = 2, step = 1))
        )
      ),
      actions = shiny::tagList(
        shiny::actionButton(ns("run"), "Run differential expression",
                            class = "btn btn-primary"),
        helper_text("Computation runs only when this button is clicked.")
      )
    ),

    # -- Status banner ---------------------------------------------------
    shiny::uiOutput(ns("status_banner")),
    shiny::uiOutput(ns("input_warning")),

    # -- Volcano + results table ----------------------------------------
    shiny::fluidRow(
      shiny::column(6,
        plot_card(
          title    = "Volcano",
          caption  = "log2FC vs -log10 adj. p",
          footnote = "Click a point to set the selected gene.",
          shiny::div(class = "plot-container",
            shiny::plotOutput(ns("volcano"), height = "420px",
                              click = shiny::clickOpts(id = ns("volcano_click"))))
        )
      ),
      shiny::column(6,
        table_card(
          title   = "DE results",
          caption = "filter, sort, inspect",
          toolbar = shiny::tagList(
            shiny::fluidRow(
              shiny::column(4, shiny::textInput(ns("gene_search"), "Gene search",
                                                value = "", placeholder = "regex or substring")),
              shiny::column(3, shiny::numericInput(ns("filter_log2fc"), "min |log2FC|",
                                                   value = 0.5, min = 0, step = 0.1)),
              shiny::column(3, shiny::numericInput(ns("filter_padj"), "max adj. p",
                                                   value = 0.05, min = 0, max = 1, step = 0.01)),
              shiny::column(2, shiny::numericInput(ns("max_rows"), "max rows",
                                                   value = 50, min = 1, step = 10))
            ),
            shiny::fluidRow(
              shiny::column(6, shiny::selectInput(ns("sort_by"), "Sort by",
                                                  choices = c("p_val_adj", "avg_log2FC", "pct.1", "pct.2", "gene"),
                                                  selected = "p_val_adj")),
              shiny::column(6, shiny::radioButtons(ns("sort_dir"), "Direction",
                                                   choices = c("descending", "ascending"),
                                                   selected = "ascending", inline = TRUE))
            )
          ),
          max_height = "bounded",
          shiny::uiOutput(ns("table_warning")),
          shiny::tableOutput(ns("de_table"))
        )
      )
    ),

    # -- Inspect panel ---------------------------------------------------
    shiny::fluidRow(
      shiny::column(4,
        app_card(
          title   = "Inspect gene",
          caption = "sync with Explorer",
          shiny::uiOutput(ns("inspect_ui")),
          action_row(
            shiny::actionButton(ns("send_to_explorer"), "Send to Explorer",
                                class = "btn btn-default")
          ),
          microcaption(shiny::textOutput(ns("inspect_status"), inline = TRUE))
        )
      ),
      shiny::column(4,
        plot_card(
          title   = shiny::textOutput(ns("feature_title"), inline = TRUE),
          caption = "feature embedding",
          shiny::uiOutput(ns("feature_warning")),
          shiny::div(class = "plot-container",
            shiny::plotOutput(ns("feature_plot"), height = "320px"))
        )
      ),
      shiny::column(4,
        plot_card(
          title   = shiny::textOutput(ns("violin_title"), inline = TRUE),
          caption = "expression by group",
          shiny::uiOutput(ns("violin_warning")),
          shiny::div(class = "plot-container",
            shiny::plotOutput(ns("violin_plot"), height = "320px"))
        )
      )
    )
  )
}

mod_differential_expression_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- Controls --------------------------------------------------------
    output$group_field_ui <- shiny::renderUI({
      fields <- available_metadata_fields(state$active_dataset)
      shiny::selectInput(ns("group_field"), "Grouping field",
                         choices  = fields,
                         selected = state$selected_metadata_field %||% fields[1])
    })
    output$assay_ui <- shiny::renderUI({
      shiny::selectInput(ns("assay"), "Assay",
                         choices  = available_assays(state$active_dataset),
                         selected = state$selected_assay)
    })

    # Backend picker -- "auto" + every registered backend with availability
    # marker. Disabled backends remain selectable so the user gets the clear
    # "install presto" error rather than silent fallback.
    output$backend_ui <- shiny::renderUI({
      backends <- de_available_backends()
      choices <- stats::setNames(
        vapply(backends, `[[`, character(1), "id"),
        vapply(backends, function(b)
          sprintf("%s%s", b$label,
                  if (!b$available) " (not installed)" else ""),
          character(1))
      )
      shiny::selectInput(ns("backend"), "DE backend",
                         choices = choices, selected = "auto")
    })
    group_levels <- shiny::reactive({
      shiny::req(input$group_field)
      v <- get_metadata(state$active_dataset, input$group_field)
      if (is.null(v)) character() else sort(unique(as.character(v)))
    })
    output$group1_ui <- shiny::renderUI({
      lv <- group_levels(); if (length(lv) == 0L) return(NULL)
      shiny::selectInput(ns("group_1"), "Group 1", choices = lv, selected = lv[1])
    })
    output$group2_ui <- shiny::renderUI({
      lv <- group_levels(); if (length(lv) == 0L) return(NULL)
      pick <- if (length(lv) >= 2) lv[2] else lv[1]
      shiny::selectInput(ns("group_2"), "Group 2", choices = lv, selected = pick)
    })

    # `sample_by` is the biological replicate column for pseudobulk
    # backends. Default to "sample" / "donor" / "patient" / "orig.ident"
    # if any of those exist; otherwise the first metadata column that is
    # not the current grouping field. Re-renders when the grouping
    # field changes so users can't accidentally pick the same column on
    # both axes.
    output$sample_by_ui <- shiny::renderUI({
      fields <- available_metadata_fields(state$active_dataset)
      if (length(fields) == 0L) return(NULL)
      candidates <- intersect(c("sample", "donor", "patient", "orig.ident"),
                              fields)
      fields_no_g <- setdiff(fields, input$group_field %||% "")
      default <- (candidates %||% fields_no_g)[1] %||% fields[1]
      shiny::selectInput(ns("sample_by"),
                         "Sample (replicate) column",
                         choices  = fields, selected = default)
    })

    # ---- Input validation banner ---------------------------------------
    output$input_warning <- shiny::renderUI({
      if (is.null(state$active_dataset)) return(NULL)  # workspace already shows empty state
      lv <- group_levels()
      if (length(lv) < 2L)
        return(friendly_warning(sprintf(
          "Field '%s' has fewer than 2 groups; pick a different grouping field.",
          input$group_field %||% "")))
      if (!is.null(input$group_1) && !is.null(input$group_2) &&
          identical(input$group_1, input$group_2))
        return(friendly_warning("Group 1 and Group 2 are identical -- pick distinct groups."))
      NULL
    })

    # Slot helpers
    de_slot     <- function() state$analysis_results$de
    set_de_slot <- function(x) {
      ar <- state$analysis_results
      ar$de <- x
      state$analysis_results <- ar
    }

    # ---- Run button ------------------------------------------------------
    shiny::observeEvent(input$run, {
      ds <- state$active_dataset
      if (is.null(ds)) return()
      if (identical(input$group_1, input$group_2)) {
        push_message(state, "Cannot run DE: Group 1 and Group 2 are identical.", "warning")
        return()
      }
      params <- list(grouping_field = input$group_field,
                     group_1 = input$group_1, group_2 = input$group_2,
                     assay   = input$assay,   layer    = input$layer,
                     min_pct = input$min_pct, test = "wilcox",
                     backend = input$backend %||% "auto",
                     sample_by = input$sample_by,
                     min_cells_per_sample  = input$min_cells_per_sample  %||% 10L,
                     min_samples_per_group = input$min_samples_per_group %||% 2L)
      set_de_slot(list(status = "running", results = NULL,
                       params = params, error_message = NULL,
                       timestamp = Sys.time(), duration_ms = NULL,
                       annotation_stamp = make_annotation_stamp(state)))
      t0 <- proc.time()[["elapsed"]]
      out <- tryCatch(
        shiny::withProgress(message = "Running DE...", value = 0.4, {
          compute_de(ds,
                     grouping_field = params$grouping_field,
                     group_1        = params$group_1,
                     group_2        = params$group_2,
                     assay          = params$assay,
                     layer          = params$layer,
                     min_pct        = params$min_pct,
                     test           = params$test,
                     backend        = params$backend,
                     sample_by             = params$sample_by,
                     min_cells_per_sample  = params$min_cells_per_sample,
                     min_samples_per_group = params$min_samples_per_group)
        }),
        error = function(e) e
      )
      if (inherits(out, "error")) {
        set_de_slot(list(status = "failed", results = NULL,
                         params = params, error_message = conditionMessage(out),
                         timestamp = Sys.time(), duration_ms = NULL,
                         annotation_stamp = make_annotation_stamp(state)))
        push_message(state, sprintf("DE failed: %s", conditionMessage(out)), "error")
      } else {
        dur <- as.integer(round((proc.time()[["elapsed"]] - t0) * 1000))
        set_de_slot(list(status = "completed", results = out,
                         params = params, error_message = NULL,
                         timestamp = Sys.time(), duration_ms = dur,
                         annotation_stamp = make_annotation_stamp(state)))
        push_message(state, sprintf(
          "DE complete: %d genes tested (%d ms).", nrow(out), dur), "success")
      }
    })

    # ---- Status banner ---------------------------------------------------
    output$status_banner <- shiny::renderUI({
      de <- de_slot()
      if (is.null(de) || identical(de$status, "not_run")) {
        return(status_banner(
          shiny::span("Not run yet. Configure the controls above and click ",
                      shiny::tags$em("Run differential expression"), "."),
          tone = "idle"))
      }
      tone <- switch(de$status,
                     "running"   = "running",
                     "completed" = "success",
                     "failed"    = "danger",
                     "idle")
      txt <- switch(de$status,
        "running"   = "Running...",
        "completed" = sprintf("Completed: %s vs %s on '%s'. %d genes tested in %d ms.",
                              de$params$group_1, de$params$group_2,
                              de$params$grouping_field,
                              nrow(de$results), de$duration_ms %||% 0L),
        "failed"    = sprintf("Failed: %s", de$error_message %||% "(unknown error)"),
        "")
      status_banner(txt, tone = tone)
    })

    # ---- Filtered + sorted view ----------------------------------------
    view <- shiny::reactive({
      de <- de_slot()
      if (is.null(de) || !identical(de$status, "completed") || is.null(de$results)) return(NULL)
      out <- filter_de_results(de$results,
                               gene_search    = input$gene_search %||% "",
                               min_abs_log2fc = input$filter_log2fc %||% 0,
                               max_padj       = input$filter_padj   %||% 1)
      sort_de_results(out,
                      sort_by    = input$sort_by %||% "p_val_adj",
                      descending = identical(input$sort_dir, "descending"))
    })

    # ---- Volcano ---------------------------------------------------------
    output$volcano <- shiny::renderPlot({
      de <- de_slot()
      if (is.null(de) || !identical(de$status, "completed")) {
        graphics::par(mar = c(4, 4, 3, 2)); graphics::plot.new()
        graphics::title(main = "Volcano (run DE to populate)")
        return(invisible(NULL))
      }
      plot_volcano(de$results,
                   log2fc_thr = input$log2fc_thr %||% 0.5,
                   padj_thr   = input$padj_thr   %||% 0.05)
    })

    # Volcano click -> selected gene
    shiny::observeEvent(input$volcano_click, {
      de <- de_slot()
      if (is.null(de) || is.null(de$results) || nrow(de$results) == 0L) return()
      df <- de$results
      df$neg_log10_padj <- -log10(pmax(df$p_val_adj, .Machine$double.eps))
      hit <- shiny::nearPoints(df, input$volcano_click,
                               xvar = "avg_log2FC", yvar = "neg_log10_padj",
                               threshold = 12, maxpoints = 1)
      if (nrow(hit) > 0L) {
        state$selected_gene <- hit$gene[1]
        push_message(state, sprintf(
          "Picked '%s' from the volcano plot.", hit$gene[1]), "info")
      }
    })

    # ---- Results table --------------------------------------------------
    output$table_warning <- shiny::renderUI({
      de <- de_slot()
      if (is.null(de) || identical(de$status, "not_run"))
        return(shiny::div(class = "p-3", helper_text("Run DE above to populate the table.")))
      if (identical(de$status, "running"))
        return(friendly_warning("Computing DE results..."))
      if (identical(de$status, "failed"))
        return(friendly_warning(sprintf("DE failed: %s", de$error_message %||% "")))
      v <- view()
      if (is.null(v) || nrow(v) == 0L)
        return(friendly_warning("No DE genes match current filters."))
      NULL
    })

    output$de_table <- shiny::renderTable({
      v <- view(); if (is.null(v) || nrow(v) == 0L) return(NULL)
      max_n <- input$max_rows %||% 50
      v <- utils::head(v, max_n)
      v$avg_log2FC <- round(v$avg_log2FC, 3)
      v$pct.1      <- round(v$pct.1, 3)
      v$pct.2      <- round(v$pct.2, 3)
      v$p_val      <- signif(v$p_val, 3)
      v$p_val_adj  <- signif(v$p_val_adj, 3)
      v
    }, striped = TRUE, hover = TRUE, rownames = FALSE, digits = 3)

    # ---- Inspect panel --------------------------------------------------
    inspect_choices <- shiny::reactive({
      v <- view()
      if (is.null(v) || nrow(v) == 0L)
        return(available_genes(state$active_dataset))
      v$gene
    })
    output$inspect_ui <- shiny::renderUI({
      choices <- inspect_choices()
      current <- state$selected_gene
      if (!current %in% choices) current <- choices[1]
      shiny::selectizeInput(ns("inspect_gene"), "Gene",
                            choices = choices, selected = current)
    })
    shiny::observeEvent(input$inspect_gene, {
      if (!is.null(input$inspect_gene) && nzchar(input$inspect_gene)) {
        state$selected_gene <- input$inspect_gene
      }
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$send_to_explorer, {
      g <- input$inspect_gene
      if (!is.null(g) && nzchar(g)) {
        state$selected_gene <- g
        push_message(state, sprintf("Sent '%s' to the Explorer FeaturePlot.", g), "success")
      }
    })
    output$inspect_status <- shiny::renderText({
      g <- state$selected_gene
      if (is.null(g) || !nzchar(g)) "" else sprintf("Selected: %s", g)
    })

    # ---- FeaturePlot (reuse Explorer's helpers) ------------------------
    output$feature_title <- shiny::renderText({
      sprintf("FeaturePlot: %s", state$selected_gene %||% "(none)")
    })
    output$feature_warning <- shiny::renderUI({
      g  <- state$selected_gene
      ds <- state$active_dataset
      if (is.null(ds))             return(friendly_warning("No dataset loaded."))
      if (!validate_gene(ds, g))   return(friendly_warning(sprintf("Gene '%s' is not available.", g %||% "")))
      if (is.null(get_embedding(ds, state$selected_reduction)))
        return(friendly_warning("No usable reduction."))
      NULL
    })
    output$feature_plot <- shiny::renderPlot({
      ds <- state$active_dataset
      shiny::req(ds, state$selected_gene, state$selected_reduction)
      emb  <- get_embedding(ds, state$selected_reduction)
      expr <- get_gene_expression(ds, state$selected_gene)
      if (is.null(emb) || is.null(expr)) return(NULL)
      plot_embedding_continuous(emb, expr,
                                title = sprintf("%s | %s", state$selected_reduction, state$selected_gene),
                                xlab  = paste0(state$selected_reduction, "_1"),
                                ylab  = paste0(state$selected_reduction, "_2"),
                                legend_title = state$selected_gene)
    })

    # ---- Violin/box for group_1 vs group_2 ----------------------------
    output$violin_title <- shiny::renderText({
      g <- state$selected_gene
      sprintf("%s expression by %s vs %s",
              g %||% "(none)",
              input$group_1 %||% "(g1)",
              input$group_2 %||% "(g2)")
    })
    output$violin_warning <- shiny::renderUI({
      g  <- state$selected_gene
      ds <- state$active_dataset
      if (is.null(ds))           return(friendly_warning("No dataset loaded."))
      if (!validate_gene(ds, g)) return(friendly_warning(sprintf("Gene '%s' is not available.", g %||% "")))
      if (is.null(input$group_field))
        return(friendly_warning("Pick a grouping field first."))
      lv <- group_levels()
      if (!(input$group_1 %in% lv) || !(input$group_2 %in% lv))
        return(friendly_warning("Pick two valid groups first."))
      NULL
    })
    output$violin_plot <- shiny::renderPlot({
      ds <- state$active_dataset
      shiny::req(ds, state$selected_gene, input$group_field, input$group_1, input$group_2)
      meta <- get_metadata(ds, input$group_field)
      expr <- get_gene_expression(ds, state$selected_gene)
      if (is.null(meta) || is.null(expr)) return(NULL)
      keep <- as.character(meta) %in% c(input$group_1, input$group_2)
      plot_expression_by_group(
        expr[keep],
        group = as.character(meta)[keep],
        title = "",
        xlab  = input$group_field,
        ylab  = sprintf("%s expression", state$selected_gene))
    })
  })
}
