# ============================================================================
# Module: Pathway Analysis  (ENABLED)
# ----------------------------------------------------------------------------
# Consumes `state$analysis_results$de` and runs simple overrepresentation
# analysis (Fisher's exact, BH adjustment) against an in-memory gene-set
# library defined in R/pathway.R.
#
# Design notes:
#   - Enrichment maths lives in R/pathway.R; this module is UI plumbing only.
#   - Results are written to `state$analysis_results$pathway` with the same
#     shape as DE results (status / results / params / error_message /
#     timestamp / duration_ms).
#   - "Click a pathway" works two ways: clicking a bar in the plot, OR
#     selecting from the dropdown. Both populate the "overlapping genes"
#     panel; clicking a gene there sets `state$selected_gene` so the
#     Explorer's FeaturePlot updates.
# ============================================================================

mod_pathway_analysis_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h2("Pathway Analysis"),
    shiny::p("Overrepresentation enrichment from Differential Expression results. ",
             "Run DE first; this module reads from ",
             shiny::tags$code("state$analysis_results$de"), "."),

    # -- DE-required banner / instructions --------------------------------
    shiny::uiOutput(ns("de_status")),

    # -- Controls ---------------------------------------------------------
    shiny::fluidRow(
      shiny::column(3, shiny::selectInput(ns("direction"), "Direction",
                                          choices = c("up in group 1" = "up_in_g1",
                                                      "up in group 2" = "up_in_g2",
                                                      "both"          = "both"),
                                          selected = "up_in_g1")),
      shiny::column(3, shiny::numericInput(ns("padj_cutoff"), "adj. p \u2264",
                                           value = 0.05, min = 0, max = 1, step = 0.01)),
      shiny::column(3, shiny::numericInput(ns("log2fc_cutoff"), "min |log2FC|",
                                           value = 0.5, min = 0, step = 0.1)),
      shiny::column(3, shiny::selectInput(ns("collection"), "Gene set collection",
                                          choices  = available_pathway_collections(),
                                          selected = available_pathway_collections()[1]))
    ),
    shiny::fluidRow(
      shiny::column(3, shiny::selectInput(ns("ranking_metric"), "Ranking metric",
                                          choices  = c("avg_log2FC", "p_val_adj"),
                                          selected = "avg_log2FC")),
      shiny::column(9, shiny::div(style = "font-size:12px; color:#888; margin-top:25px;",
        "Ranking metric is a placeholder for future GSEA support; ORA is set-based and ignores ranking."))
    ),
    shiny::div(style = "margin: 8px 0 16px 0;",
      shiny::actionButton(ns("run"), "Run Pathway Analysis",
                          class = "btn btn-primary"),
      shiny::tags$span(style = "margin-left:12px; color:#888; font-size:12px;",
                       "Computation runs only when this button is clicked.")
    ),

    shiny::uiOutput(ns("status_banner")),
    shiny::uiOutput(ns("input_warning")),
    shiny::uiOutput(ns("selection_summary")),

    shiny::hr(),

    # -- Plot + table -----------------------------------------------------
    shiny::fluidRow(
      shiny::column(6,
        shiny::h4("Top pathways"),
        shiny::div(style = "font-size:12px; color:#888;",
                   "Click a bar to inspect overlapping genes."),
        shiny::plotOutput(ns("bar_plot"), height = "420px",
                          click = shiny::clickOpts(id = ns("bar_click")))
      ),
      shiny::column(6,
        shiny::h4("Enrichment results"),
        shiny::uiOutput(ns("table_warning")),
        shiny::div(style = "max-height:420px; overflow:auto;",
                   shiny::tableOutput(ns("results_table")))
      )
    ),

    shiny::hr(),

    # -- Pathway inspection ---------------------------------------------
    shiny::fluidRow(
      shiny::column(4,
        shiny::h4("Inspect pathway"),
        shiny::uiOutput(ns("pathway_picker_ui")),
        shiny::div(style = "font-size:13px; margin-top:8px;",
                   shiny::strong("Overlapping genes:")),
        shiny::uiOutput(ns("overlap_ui"))
      ),
      shiny::column(8,
        shiny::h4("Send a gene to the Explorer"),
        shiny::uiOutput(ns("gene_picker_ui")),
        shiny::actionButton(ns("send_to_explorer"), "Send to Explorer",
                            class = "btn btn-default"),
        shiny::div(style = "margin-top:10px; font-size:12px; color:#888;",
                   shiny::textOutput(ns("send_status")))
      )
    )
  )
}

mod_pathway_analysis_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- DE prerequisite --------------------------------------------------
    de_slot <- function() state$analysis_results$de
    de_ready <- shiny::reactive({
      d <- de_slot()
      !is.null(d) && identical(d$status, "completed") &&
        !is.null(d$results) && nrow(d$results) > 0L
    })

    output$de_status <- shiny::renderUI({
      if (de_ready()) return(NULL)
      shiny::div(
        style = "padding:10px 14px; background:#fff3cd; color:#664d03; border-radius:4px; font-size:13px;",
        shiny::strong("No DE results yet. "),
        "Open ", shiny::tags$em("Differential Expression"),
        " and click ", shiny::tags$em("Run Differential Expression"),
        " to populate ", shiny::tags$code("state$analysis_results$de"), "."
      )
    })

    # ---- Slot helpers (mirrors DE pattern) -------------------------------
    pw_slot     <- function() state$analysis_results$pathway
    set_pw_slot <- function(x) {
      ar <- state$analysis_results
      ar$pathway <- x
      state$analysis_results <- ar
    }

    # ---- Selection summary (live preview before Run) ---------------------
    selected_preview <- shiny::reactive({
      if (!de_ready()) return(NULL)
      d <- de_slot()$results
      sel <- select_de_genes(d, direction = input$direction %||% "up_in_g1",
                             padj_cutoff   = input$padj_cutoff   %||% 0.05,
                             log2fc_cutoff = input$log2fc_cutoff %||% 0.5)
      sel
    })
    output$selection_summary <- shiny::renderUI({
      if (!de_ready()) return(NULL)
      sel <- selected_preview()
      shiny::div(style = "padding:6px 12px; background:#eef; border-radius:4px; font-size:12px; margin:6px 0;",
        shiny::strong(length(sel), "genes selected"), " from DE results based on direction + thresholds.",
        if (length(sel) > 0L)
          shiny::tags$div(style = "margin-top:2px; color:#444; font-family: monospace;",
            paste(utils::head(sel, 12), collapse = ", "),
            if (length(sel) > 12) sprintf(" (+%d more)", length(sel) - 12L) else "")
      )
    })

    output$input_warning <- shiny::renderUI({
      if (!de_ready()) return(NULL)
      sel <- selected_preview()
      if (length(sel) == 0L)
        return(friendly_warning("No DE genes pass the current direction + thresholds. Loosen padj or log2FC, or switch direction."))
      coll <- input$collection %||% available_pathway_collections()[1]
      pw <- get_pathways(coll)
      if (is.null(pw) || length(pw) == 0L)
        return(friendly_warning(sprintf("Gene set collection '%s' is empty.", coll)))
      all_pw_genes <- unique(unlist(pw, use.names = FALSE))
      if (length(intersect(sel, all_pw_genes)) == 0L)
        return(friendly_warning("No selected genes match any pathway in the chosen collection."))
      NULL
    })

    # ---- Run button ------------------------------------------------------
    shiny::observeEvent(input$run, {
      if (!de_ready()) {
        push_message(state, "Cannot run pathway: DE results not available.", "warning")
        return()
      }
      coll <- input$collection
      pw <- get_pathways(coll)
      if (is.null(pw) || length(pw) == 0L) {
        push_message(state, "Cannot run pathway: empty gene set collection.", "warning")
        return()
      }
      params <- list(
        direction      = input$direction,
        padj_cutoff    = input$padj_cutoff,
        log2fc_cutoff  = input$log2fc_cutoff,
        collection     = coll,
        ranking_metric = input$ranking_metric,
        de_run_at      = de_slot()$timestamp
      )
      set_pw_slot(list(status = "running", results = NULL,
                       params = params, error_message = NULL,
                       timestamp = Sys.time(), duration_ms = NULL,
                       annotation_stamp = make_annotation_stamp(state)))
      t0 <- proc.time()[["elapsed"]]
      out <- tryCatch(
        shiny::withProgress(message = "Running pathway enrichment...", value = 0.4, {
          de_df <- de_slot()$results
          de_tested <- unique(de_df$gene)
          all_pw_genes <- unique(unlist(pw, use.names = FALSE))
          # Universe = union of tested genes and pathway genes. See R/pathway.R
          # for rationale + TODO on configurable universe.
          universe <- unique(c(de_tested, all_pw_genes))

          run_one <- function(dir) {
            sel <- select_de_genes(de_df, direction = dir,
                                   padj_cutoff   = params$padj_cutoff,
                                   log2fc_cutoff = params$log2fc_cutoff)
            compute_enrichment(sel, pw,
                               universe   = universe,
                               direction  = dir,
                               collection = coll)
          }

          if (identical(params$direction, "both")) {
            res <- rbind(run_one("up_in_g1"), run_one("up_in_g2"))
            # Re-adjust p-values across the combined run for honest BH.
            res$p_val_adj <- stats::p.adjust(res$p_val, method = "BH")
            res <- res[order(res$p_val_adj, res$p_val, -res$n_overlap), , drop = FALSE]
            rownames(res) <- NULL
            res
          } else {
            run_one(params$direction)
          }
        }),
        error = function(e) e
      )
      if (inherits(out, "error")) {
        set_pw_slot(list(status = "failed", results = NULL,
                         params = params, error_message = conditionMessage(out),
                         timestamp = Sys.time(), duration_ms = NULL,
                         annotation_stamp = make_annotation_stamp(state)))
        push_message(state, sprintf("Pathway analysis failed: %s", conditionMessage(out)), "error")
      } else {
        dur <- as.integer(round((proc.time()[["elapsed"]] - t0) * 1000))
        set_pw_slot(list(status = "completed", results = out,
                         params = params, error_message = NULL,
                         timestamp = Sys.time(), duration_ms = dur,
                         annotation_stamp = make_annotation_stamp(state)))
        n_enriched <- if (is.null(out)) 0L
                      else sum(!is.na(out$p_val_adj) & out$p_val_adj <= 0.05)
        push_message(state, sprintf(
          "Pathway done: %d pathways tested (%d at adj. p \u2264 0.05; %d ms).",
          if (is.null(out)) 0L else nrow(out), n_enriched, dur), "success")
      }
    })

    # ---- Status banner ---------------------------------------------------
    output$status_banner <- shiny::renderUI({
      p <- pw_slot()
      if (is.null(p)) return(shiny::div(
        style = "padding:8px 12px; background:#eee; border-radius:4px; font-size:13px;",
        shiny::tags$strong("Status: "), "Not run yet."))
      bg <- switch(p$status, running = "#cfe2ff", completed = "#d1e7dd",
                   failed = "#f8d7da", "#eee")
      fg <- switch(p$status, running = "#084298", completed = "#0a3622",
                   failed = "#842029", "#333")
      txt <- switch(p$status,
        running   = "Running...",
        completed = sprintf("Completed: %s | direction=%s | %d pathways tested in %d ms.",
                            p$params$collection, p$params$direction,
                            if (is.null(p$results)) 0L else nrow(p$results),
                            p$duration_ms %||% 0L),
        failed    = sprintf("Failed: %s", p$error_message %||% "(unknown)"),
        "")
      shiny::div(style = sprintf("padding:8px 12px; background:%s; color:%s; border-radius:4px; font-size:13px;",
                                 bg, fg),
                 shiny::tags$strong("Status: "), txt)
    })

    # ---- Bar plot / table ------------------------------------------------
    output$bar_plot <- shiny::renderPlot({
      p <- pw_slot()
      if (is.null(p) || !identical(p$status, "completed") ||
          is.null(p$results) || nrow(p$results) == 0L) {
        graphics::par(mar = c(4, 4, 3, 2)); graphics::plot.new()
        graphics::title(main = "Top pathways (run analysis to populate)")
        return(invisible(NULL))
      }
      plot_pathway_enrichment(p$results, top_n = 10L, metric = "padj")
    })

    # Map a bar-plot click to a pathway via its y position.
    bar_map <- shiny::reactiveVal(NULL)
    shiny::observe({
      p <- pw_slot()
      bar_map(NULL)
      if (is.null(p) || !identical(p$status, "completed") ||
          is.null(p$results) || nrow(p$results) == 0L) return()
      # Recreate the same ordering as plot_pathway_enrichment for click mapping
      df <- p$results[order(p$results$p_val_adj, p$results$p_val), , drop = FALSE]
      df <- utils::head(df, 10L)
      df <- df[rev(seq_len(nrow(df))), , drop = FALSE]
      bar_map(df$pathway)
    })

    shiny::observeEvent(input$bar_click, {
      pwn <- bar_map()
      if (is.null(pwn) || length(pwn) == 0L) return()
      y <- input$bar_click$y
      if (is.null(y) || !is.finite(y)) return()
      idx <- max(1L, min(length(pwn), round(y)))
      shiny::updateSelectInput(session, "inspect_pathway", selected = pwn[idx])
    })

    output$table_warning <- shiny::renderUI({
      p <- pw_slot()
      if (is.null(p) || identical(p$status, "not_run"))
        return(shiny::div(style = "padding:8px 12px; font-size:13px; color:#666;",
                          "Run pathway analysis above to populate the table."))
      if (identical(p$status, "running"))
        return(friendly_warning("Computing enrichment..."))
      if (identical(p$status, "failed"))
        return(friendly_warning(sprintf("Failed: %s", p$error_message %||% "")))
      if (is.null(p$results) || nrow(p$results) == 0L)
        return(friendly_warning("No enriched pathways found."))
      NULL
    })

    output$results_table <- shiny::renderTable({
      p <- pw_slot()
      if (is.null(p) || is.null(p$results) || nrow(p$results) == 0L) return(NULL)
      df <- p$results
      # Truncate the overlap_genes column for readability.
      df$overlap_genes <- vapply(df$overlap_genes, function(s) {
        g <- strsplit(s, ";", fixed = TRUE)[[1]]
        if (length(g) <= 4L) paste(g, collapse = ";")
        else sprintf("%s; (+%d more)", paste(g[1:4], collapse = ";"), length(g) - 4L)
      }, FUN.VALUE = character(1))
      df$odds_ratio <- signif(df$odds_ratio, 3)
      df$p_val      <- signif(df$p_val,      3)
      df$p_val_adj  <- signif(df$p_val_adj,  3)
      df
    }, striped = TRUE, hover = TRUE, rownames = FALSE)

    # ---- Inspect pathway -> overlapping genes -> picker -----------------
    output$pathway_picker_ui <- shiny::renderUI({
      p <- pw_slot()
      if (is.null(p) || is.null(p$results) || nrow(p$results) == 0L)
        return(shiny::div(style = "color:#888; font-size:13px;", "(no results yet)"))
      df <- p$results[order(p$results$p_val_adj, p$results$p_val), , drop = FALSE]
      shiny::selectInput(ns("inspect_pathway"), label = NULL,
                         choices = df$pathway, selected = df$pathway[1])
    })

    overlap_genes_of <- shiny::reactive({
      p <- pw_slot()
      pwn <- input$inspect_pathway
      if (is.null(p) || is.null(pwn)) return(character())
      row <- p$results[p$results$pathway == pwn, , drop = FALSE]
      if (nrow(row) == 0L) return(character())
      s <- row$overlap_genes[1]
      if (!nzchar(s)) return(character())
      strsplit(s, ";", fixed = TRUE)[[1]]
    })

    output$overlap_ui <- shiny::renderUI({
      g <- overlap_genes_of()
      if (length(g) == 0L)
        return(shiny::tags$div(style = "color:#888; font-size:13px;", "(no overlapping genes)"))
      shiny::tags$div(style = "font-family: monospace; font-size:13px; color:#222;",
                      paste(g, collapse = ", "))
    })

    output$gene_picker_ui <- shiny::renderUI({
      g <- overlap_genes_of()
      ds_genes <- available_genes(state$active_dataset)
      # Prefer overlap genes that exist in the dataset (so FeaturePlot works).
      pickable <- intersect(g, ds_genes)
      if (length(pickable) == 0L) pickable <- ds_genes
      shiny::selectInput(ns("pick_gene"), "Gene",
                         choices = pickable,
                         selected = utils::head(pickable, 1))
    })

    shiny::observeEvent(input$send_to_explorer, {
      g <- input$pick_gene
      if (is.null(g) || !nzchar(g)) return()
      if (!validate_gene(state$active_dataset, g)) {
        push_message(state, sprintf("Gene '%s' not in dataset; FeaturePlot unavailable.", g), "warning")
        return()
      }
      state$selected_gene <- g
      push_message(state, sprintf("Sent '%s' to the Explorer FeaturePlot.", g), "success")
    })

    output$send_status <- shiny::renderText({
      g <- state$selected_gene
      if (is.null(g) || !nzchar(g)) "" else sprintf("Selected: %s", g)
    })
  })
}
