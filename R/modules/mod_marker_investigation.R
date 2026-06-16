# ============================================================================
# Module: Marker Investigation  (ENABLED, v1)
# ----------------------------------------------------------------------------
# Rank genes that distinguish each level of a metadata field. Lets the user
# push a chosen gene back into shared state so the Explorer's FeaturePlot
# updates -- this is the first cross-module loop in the app.
#
# v1 uses mock expression in `dataset$expression` via R/markers.R. A future
# real implementation should keep the same output schema:
#   group | gene | avg_log2FC | pct_in | pct_out | p_value
# so nothing in this module has to change.
#
# UI composition uses the shared primitives in `R/ui_components.R`. Server
# logic (compute, observers) is unchanged.
# ============================================================================

mod_marker_investigation_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    page_header(
      eyebrow = "Cell Identity",
      title   = "Marker Investigation",
      lede    = paste("Find genes that distinguish groups within a metadata",
                      "field. Send any marker gene to the Explorer to color",
                      "its FeaturePlot.")
    ),

    control_panel(
      title = "Marker query",
      shiny::fluidRow(
        shiny::column(3, shiny::uiOutput(ns("group_field_ui"))),
        shiny::column(3, shiny::uiOutput(ns("group_ui"))),
        shiny::column(3, shiny::numericInput(ns("top_n"), "Top N per group",
                                             value = 6, min = 1, max = 50, step = 1)),
        shiny::column(3, shiny::numericInput(ns("min_log2fc"), "min |log2FC|",
                                             value = 0, min = 0, step = 0.1))
      ),
      actions = shiny::tagList(
        shiny::actionButton(ns("recompute"), "Find markers",
                            class = "btn btn-primary"),
        helper_text("Markers refresh only when ",
                    shiny::tags$em("Find markers"), " is clicked.")
      )
    ),

    table_card(
      title   = shiny::textOutput(ns("table_title"), inline = TRUE),
      caption = "ranked by avg_log2FC and pct_in/out",
      shiny::uiOutput(ns("table_warning")),
      shiny::tableOutput(ns("markers_table")),
      max_height = "420px"
    ),

    shiny::fluidRow(
      shiny::column(6,
        app_card(
          title   = "Push to Explorer",
          caption = "send a marker into the shared state",
          shiny::uiOutput(ns("highlight_ui")),
          action_row(
            shiny::actionButton(ns("send_to_explorer"), "Send to Explorer",
                                class = "btn btn-default"),
            helper = shiny::textOutput(ns("push_status"), inline = TRUE)
          )
        )
      ),
      shiny::column(6,
        plot_card(
          title   = shiny::textOutput(ns("box_title"), inline = TRUE),
          caption = "expression by group",
          shiny::uiOutput(ns("box_warning")),
          shiny::div(class = "plot-container",
            shiny::plotOutput(ns("box_plot"), height = "320px"))
        )
      )
    )
  )
}

mod_marker_investigation_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- Controls --------------------------------------------------------
    output$group_field_ui <- shiny::renderUI({
      fields <- available_metadata_fields(state$active_dataset)
      shiny::selectInput(ns("group_field"), "Group cells by",
                         choices  = fields,
                         selected = state$selected_metadata_field %||% fields[1])
    })
    output$group_ui <- shiny::renderUI({
      vals <- get_metadata(state$active_dataset, input$group_field %||% state$selected_metadata_field)
      if (is.null(vals)) {
        return(shiny::selectInput(ns("group"), "Group of interest",
                                  choices = "All groups"))
      }
      shiny::selectInput(ns("group"), "Group of interest",
                         choices = c("All groups", sort(unique(vals))),
                         selected = "All groups")
    })

    # Keep shared state's "selected metadata field" in sync with local choice
    # so the Explorer and Marker module agree on grouping.
    shiny::observeEvent(input$group_field, {
      state$selected_metadata_field <- input$group_field
    }, ignoreInit = TRUE)

    # ---- Marker computation ---------------------------------------------
    # The UI advertises "Markers refresh only when Find markers is
    # clicked." -- but the previous `reactive({ input$recompute; ... })`
    # pattern only ADDS `input$recompute` to the dependency graph; it
    # does NOT gate execution on a click. Result: the moment the user
    # navigated to this tab, `output$table_warning` consumed
    # `markers()`, which ran `compute_markers()` on the full active
    # dataset (8000 cells x 32000 genes for PBMC 8k) and froze the
    # session for seconds-to-minutes. Use `eventReactive(input$recompute,
    # ..., ignoreInit = TRUE)` so the compute is genuinely click-only.
    # Before the first click, `markers()` returns NULL and downstream
    # consumers render the "click to compute" hint.
    markers <- shiny::eventReactive(input$recompute, {
      shiny::req(state$active_dataset, input$group_field)
      grp_filter <- if (isTRUE(input$group != "All groups")) input$group else NULL
      params <- list(grouping_field = input$group_field,
                     group_filter   = grp_filter,
                     top_n          = input$top_n %||% Inf,
                     min_log2fc     = input$min_log2fc %||% 0)

      # Surface progress on stdout + per-run log file. The browser-side
      # `withProgress()` overlay is the only signal users *normally*
      # see for `compute_markers()`, which is fine when the compute is
      # snappy. On a real PBMC 8k Seurat object (~32k genes x 8
      # clusters = ~256k wilcox.tests in pure R) the call can sit at
      # 5-10 minutes; without console output that's indistinguishable
      # from a frozen UI. Wrap in `sce_run_with_log()` so:
      #   - stderr/stdout (incl. message()/warning()) are tee'd to
      #     the R console *and* a `marker_compute_YYYYMMDD_HHMMSS.log`
      #     under the existing setup log dir;
      #   - we get an explicit "[ts] start" / "[ts] DONE/ERROR"
      #     summary line whether the compute succeeded or threw.
      log_file <- tryCatch(sce_open_log("marker_compute"),
                           error = function(e) NULL)
      start_msg <- sprintf(
        "Marker compute starting: field=%s, group_filter=%s, top_n=%s, min_log2fc=%s",
        params$grouping_field,
        format(params$group_filter %||% "(all)"),
        format(params$top_n),
        format(params$min_log2fc))

      run_compute <- function() {
        message(start_msg)
        shiny::withProgress(
          message = "Computing markers",
          detail  = sprintf("group: %s", params$grouping_field),
          value   = 0.5,
          compute_markers(state$active_dataset,
                          grouping_field = params$grouping_field,
                          group_filter   = params$group_filter,
                          top_n          = params$top_n))
      }

      t0 <- Sys.time()
      df <- tryCatch(
        if (!is.null(log_file))
          sce_run_with_log(run_compute(), log_path = log_file)
        else
          run_compute(),
        error = function(e) {
          message(sprintf("Marker compute FAILED: %s", conditionMessage(e)))
          NULL
        })
      dur_ms <- as.integer(as.numeric(difftime(Sys.time(), t0, units = "secs")) * 1000)

      if (!is.null(df)) {
        df <- df[abs(df$avg_log2FC) >= params$min_log2fc, , drop = FALSE]
        if (nrow(df) == 0L) df <- NULL
      }
      message(sprintf("Marker compute %s in %d ms: %d rows",
                      if (is.null(df)) "FAILED" else "DONE",
                      dur_ms,
                      if (is.null(df)) 0L else nrow(df)))

      # Mirror to state$analysis_results$markers with the annotation stamp,
      # so downstream tools (and the active-set-aware "is_result_stale"
      # check) can see what annotation context produced this table.
      shiny::isolate({
        ar <- state$analysis_results
        ar$markers <- list(
          status           = if (is.null(df)) "failed" else "completed",
          results          = df,
          params           = params,
          error_message    = if (is.null(df)) "no markers passed filters" else NULL,
          timestamp        = Sys.time(),
          duration_ms      = dur_ms,
          annotation_stamp = make_annotation_stamp(state)
        )
        state$analysis_results <- ar
      })
      df
    }, ignoreInit = TRUE)

    # Cheap derived value: TRUE once a markers compute has produced a
    # result this session. Distinguishes "compute hasn't run yet"
    # (show hint) from "compute ran but produced no rows" (show
    # empty-result message).
    #
    # We deliberately do NOT key this off `input$recompute` alone.
    # `input$recompute` is the actionButton click counter -- if the
    # button's DOM ever gets re-created the counter resets to 0 even
    # though the computed markers are still safely in
    # `state$analysis_results$markers`. So fall back to the persisted
    # state, which survives any UI churn.
    has_run <- shiny::reactive({
      n <- input$recompute %||% 0L
      if (isTRUE(as.integer(n) > 0L)) return(TRUE)
      ar <- state$analysis_results
      isTRUE(!is.null(ar) && !is.null(ar$markers) &&
             identical(ar$markers$status, "completed"))
    })

    # The table-driving data frame. Prefers the freshly-computed result
    # from `markers()` (in-session click), and falls back to the
    # persisted `state$analysis_results$markers$results` so navigating
    # away and back, or otherwise re-mounting this module's UI, still
    # shows the last computed table without forcing a recompute.
    effective_markers <- shiny::reactive({
      if (isTRUE(as.integer(input$recompute %||% 0L) > 0L)) {
        df <- markers()
        if (!is.null(df)) return(df)
      }
      ar <- state$analysis_results
      if (!is.null(ar) && !is.null(ar$markers) &&
          identical(ar$markers$status, "completed")) {
        return(ar$markers$results)
      }
      NULL
    })

    output$table_title <- shiny::renderText({
      gf <- input$group_field %||% state$selected_metadata_field %||% "(none)"
      grp <- input$group %||% "All groups"
      if (grp == "All groups") sprintf("Top markers per %s", gf)
      else                     sprintf("Top markers in %s = %s", gf, grp)
    })

    output$table_warning <- shiny::renderUI({
      if (is.null(state$active_dataset))
        return(friendly_warning("No dataset loaded."))
      # `has_run()` distinguishes "user hasn't clicked yet" from
      # "compute ran but returned nothing". Without it the pre-click
      # state misleadingly reported "no markers pass the current
      # filters", which sounds like a filter problem.
      if (!has_run())
        return(friendly_warning(paste0(
          "No markers computed yet. Configure 'Group cells by' / 'Top N' ",
          "above, then click 'Find markers' to compute. On large datasets ",
          "(8000+ cells) the first run can take 10-60s.")))
      if (is.null(effective_markers()))
        return(friendly_warning("No markers pass the current filters."))
      NULL
    })

    output$markers_table <- shiny::renderTable({
      # Read from `effective_markers()` so navigating away and back
      # still shows the last result via state$analysis_results$markers.
      # `effective_markers()` only forces the heavy `markers()`
      # compute when the user has actually clicked the button this
      # session, so we don't accidentally recompute on tab open.
      if (!has_run()) return(NULL)
      df <- effective_markers(); if (is.null(df)) return(NULL)
      df$avg_log2FC <- round(df$avg_log2FC, 3)
      df$pct_in     <- round(df$pct_in,     3)
      df$pct_out    <- round(df$pct_out,    3)
      df$p_value    <- signif(df$p_value,   3)
      df
    }, striped = TRUE, hover = TRUE, digits = 3, rownames = FALSE)

    # ---- Highlight gene + push to Explorer -----------------------------
    # Server-side gene picker. See R/ui_components.R for rationale.
    # Choices flip between "all dataset genes" (no markers computed
    # yet) and "marker genes only" (after a run). Same `has_run()`
    # guard as above: we *want* the picker populated before any click
    # so the user can boxplot any gene, but we MUST NOT call
    # `markers()` until they've explicitly asked for it.
    output$highlight_ui <- shiny::renderUI({
      gene_picker_input(ns("highlight_gene"), "Highlight gene")
    })
    shiny::observe({
      ds <- state$active_dataset; shiny::req(ds)
      genes <- if (has_run()) {
        df <- effective_markers()
        if (is.null(df)) available_genes(ds) else unique(df$gene)
      } else {
        available_genes(ds)
      }
      update_gene_picker(session, "highlight_gene",
                         choices  = genes,
                         selected = genes[1])
    })

    shiny::observeEvent(input$send_to_explorer, {
      g <- input$highlight_gene
      if (!is.null(g) && nzchar(g)) {
        state$selected_gene <- g
        push_message(state, sprintf("Sent '%s' to the Explorer FeaturePlot.", g), "success")
      }
    })

    output$push_status <- shiny::renderText({
      if (is.null(state$selected_gene) || !nzchar(state$selected_gene)) ""
      else sprintf("Explorer FeaturePlot = %s", state$selected_gene)
    })

    # ---- Boxplot of expression by group --------------------------------
    output$box_title <- shiny::renderText({
      g  <- input$highlight_gene %||% "(none)"
      gf <- input$group_field    %||% "(none)"
      sprintf("%s expression by %s", g, gf)
    })

    output$box_warning <- shiny::renderUI({
      ds <- state$active_dataset
      g  <- input$highlight_gene
      gf <- input$group_field
      if (is.null(ds))                  return(friendly_warning("No dataset loaded."))
      if (!validate_gene(ds, g))        return(friendly_warning(sprintf("Gene '%s' is not available.",  g %||% "")))
      if (!validate_metadata(ds, gf))   return(friendly_warning(sprintf("Metadata field '%s' is not available.", gf %||% "")))
      NULL
    })

    output$box_plot <- shiny::renderPlot({
      ds <- state$active_dataset
      shiny::req(ds, input$highlight_gene, input$group_field)
      expr <- get_gene_expression(ds, input$highlight_gene)
      meta <- get_metadata(ds, input$group_field)
      if (is.null(expr) || is.null(meta)) return(NULL)
      f <- as.factor(meta)
      pal <- grDevices::hcl.colors(max(length(levels(f)), 2L), palette = "Dark 3")
      graphics::par(mar = c(5, 4, 2, 2))
      graphics::boxplot(expr ~ f,
                        col   = pal,
                        xlab  = input$group_field,
                        ylab  = sprintf("%s expression", input$highlight_gene),
                        main  = "",
                        outline = FALSE,
                        las   = if (length(levels(f)) > 4) 2 else 1)
    })
  })
}
