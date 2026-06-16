# ============================================================================
# Module: Trajectory / Pseudotime  (ENABLED)
# ----------------------------------------------------------------------------
# Generate a deterministic mock pseudotime from a chosen root cluster, OR
# rescale a numeric metadata column directly. Result lives in
# `state$analysis_results$trajectory` -- it is never silently injected into
# DE / Markers / Pathway / Imputation. Those modules continue to read raw
# expression and the standard metadata columns.
#
# UI / state contract:
#   - `state$selected_gene` is synced with the local gene picker so the
#     gene-vs-pseudotime plot moves in lockstep with the Explorer/DE/etc.
#   - `state$active_dataset` and `cell_data` are NOT mutated. If a user
#     ever wants pseudotime as a metadata column for the Explorer to color
#     by, the `pseudotime_demo` column already provided by the mock
#     dataset is the supported precomputed path.
# ============================================================================

mod_trajectory_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    page_header(
      eyebrow = "Advanced Analysis",
      title   = "Trajectory / Pseudotime",
      lede    = paste("Generate a deterministic mock pseudotime from a",
                      "chosen root cluster, or rescale a numeric metadata",
                      "column. Optionally bake the result into a new",
                      "metadata column for use across modules.")
    ),

    info_banner(
      tone  = "warning",
      title = "Mock pseudotime \u2014 exploratory only.",
      "This is not lineage inference. DE, Marker Investigation, Pathway ",
      "Analysis, and Imputation continue to use raw expression and the ",
      "standard metadata columns; they do not consume pseudotime."
    ),

    control_panel(
      title = "Pseudotime source",
      shiny::fluidRow(
        shiny::column(3, shiny::uiOutput(ns("reduction_ui"))),
        shiny::column(3, shiny::uiOutput(ns("source_ui"))),
        shiny::column(6, shiny::uiOutput(ns("source_specific_ui")))
      ),
      shiny::fluidRow(
        shiny::column(3,
          shiny::checkboxInput(ns("apply_with_bins"),
                               "Also add bin column when applying",
                               value = FALSE)),
        shiny::column(3,
          shiny::numericInput(ns("apply_bins"), "Number of bins",
                              value = 10L, min = 2L, max = 100L, step = 1L))
      ),
      actions = shiny::tagList(
        shiny::actionButton(ns("run"), "Run / generate pseudotime",
                            class = "btn btn-primary"),
        shiny::actionButton(ns("apply_to_metadata"),
                            "Apply pseudotime to dataset",
                            class = "btn btn-default"),
        helper_text(
          "Pseudotime is generated only when Run is clicked. ",
          "Apply bakes the result into a new metadata column.")
      )
    ),

    shiny::uiOutput(ns("status_banner")),
    shiny::uiOutput(ns("input_warning")),

    # -- Summary + two embedding plots side-by-side ----------------------
    shiny::fluidRow(
      shiny::column(3,
        app_card(
          title   = "Summary",
          caption = "pseudotime stats",
          shiny::verbatimTextOutput(ns("summary"))
        )
      ),
      shiny::column(5,
        plot_card(
          title   = "Embedding by pseudotime",
          caption = "from the selected reduction",
          shiny::uiOutput(ns("pt_plot_warning")),
          shiny::div(class = "plot-container",
            shiny::plotOutput(ns("pt_plot"), height = "400px"))
        )
      ),
      shiny::column(4,
        plot_card(
          title   = "Root group",
          caption = "highlights the chosen root",
          shiny::uiOutput(ns("root_plot_warning")),
          shiny::div(class = "plot-container",
            shiny::plotOutput(ns("root_plot"), height = "400px"))
        )
      )
    ),

    # -- Gene trend along pseudotime -------------------------------------
    shiny::fluidRow(
      shiny::column(3,
        app_card(
          title   = "Gene",
          caption = "synced with Explorer",
          shiny::uiOutput(ns("gene_picker_ui")),
          action_row(
            shiny::actionButton(ns("send_to_explorer"), "Send to Explorer",
                                class = "btn btn-default")
          ),
          microcaption(shiny::textOutput(ns("send_status"), inline = TRUE))
        )
      ),
      shiny::column(9,
        plot_card(
          title   = shiny::textOutput(ns("trend_title"), inline = TRUE),
          caption = "binned average expression",
          shiny::uiOutput(ns("trend_warning")),
          shiny::div(class = "plot-container",
            shiny::plotOutput(ns("trend_plot"), height = "360px"))
        )
      )
    )
  )
}

mod_trajectory_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- Controls -------------------------------------------------------
    output$reduction_ui <- shiny::renderUI({
      shiny::selectInput(ns("reduction"), "Reduction",
                         choices  = available_reductions(state$active_dataset),
                         selected = state$selected_reduction)
    })

    # The method picker is driven by the registry so adding a new
    # backend in R/trajectory_*.R surfaces here automatically.
    # Unavailable methods (missing optional packages) are kept in the
    # list but labelled, so clicking them produces a clear "install X"
    # error rather than a silent fallback.
    output$source_ui <- shiny::renderUI({
      shiny::selectInput(ns("source"), "Pseudotime source",
                         choices  = trajectory_method_choices(),
                         selected = "mock")
    })

    # The right-hand control panel adapts to the chosen method:
    #   * metadata             -> numeric field picker only
    #   * methods needing root -> root_field + root_group
    output$source_specific_ui <- shiny::renderUI({
      ds <- state$active_dataset; if (is.null(ds)) return(NULL)
      method <- get_trajectory_method(input$source %||% "mock")

      if (identical(input$source, "metadata")) {
        nums <- available_numeric_metadata_fields(ds)
        if (length(nums) == 0L) {
          return(friendly_warning("No numeric metadata fields available; switch to a different source."))
        }
        return(shiny::selectInput(
          ns("metadata_field"), "Numeric metadata field",
          choices = nums, selected = nums[1]))
      }

      # Any method that needs a root (mock / slingshot / monocle3)
      # gets the same cluster + start-cluster controls.
      if (!is.null(method) && isTRUE(method$requires_root)) {
        cat_fields <- setdiff(available_metadata_fields(ds),
                              available_numeric_metadata_fields(ds))
        if (length(cat_fields) == 0L) cat_fields <- available_metadata_fields(ds)
        return(shiny::tagList(
          shiny::fluidRow(
            shiny::column(6, shiny::selectInput(
              ns("root_field"), "Root field (cluster labels)",
              choices = cat_fields,
              selected = if ("cluster" %in% cat_fields) "cluster" else cat_fields[1])),
            shiny::column(6, shiny::uiOutput(ns("root_group_ui")))
          )
        ))
      }
      NULL
    })

    output$root_group_ui <- shiny::renderUI({
      ds <- state$active_dataset
      if (is.null(ds) || is.null(input$root_field)) return(NULL)
      vals <- get_metadata(ds, input$root_field)
      if (is.null(vals)) return(NULL)
      lv <- sort(unique(as.character(vals)))
      shiny::selectInput(ns("root_group"), "Root group",
                         choices = lv, selected = lv[1])
    })

    # ---- Slot helpers ---------------------------------------------------
    tr_slot     <- function() state$analysis_results$trajectory
    set_tr_slot <- function(x) {
      ar <- state$analysis_results
      ar$trajectory <- x
      state$analysis_results <- ar
    }

    # ---- Input validation banner ---------------------------------------
    output$input_warning <- shiny::renderUI({
      if (is.null(state$active_dataset)) return(NULL)
      if (identical(input$source, "metadata")) {
        if (is.null(input$metadata_field)) return(NULL)
        v <- get_metadata(state$active_dataset, input$metadata_field)
        if (!is.null(v) && !is.numeric(v))
          return(friendly_warning(sprintf("Field '%s' is not numeric.",
                                          input$metadata_field)))
      }
      NULL
    })

    # ---- Run ------------------------------------------------------------
    shiny::observeEvent(input$run, {
      ds <- state$active_dataset; if (is.null(ds)) return()
      params <- list(source = input$source,
                     reduction = input$reduction,
                     root_field = input$root_field,
                     root_group = input$root_group,
                     metadata_field = input$metadata_field)
      set_tr_slot(list(status = "running", results = NULL,
                       params = params, error_message = NULL,
                       timestamp = Sys.time(), duration_ms = NULL,
                       annotation_stamp = make_annotation_stamp(state)))
      t0 <- proc.time()[["elapsed"]]
      out <- tryCatch(
        shiny::withProgress(message = "Generating pseudotime...", value = 0.4, {
          # All trajectory methods route through the registry; the
          # `run_trajectory` orchestrator hands us a fully-populated
          # results payload (n_lineages, method_details, etc.).
          run_trajectory(
            ds,
            source         = params$source,
            reduction      = params$reduction,
            root_field     = params$root_field,
            root_group     = params$root_group,
            metadata_field = params$metadata_field,
            cluster_field  = params$root_field)
        }),
        error = function(e) e
      )
      if (inherits(out, "error")) {
        set_tr_slot(list(status = "failed", results = NULL,
                         params = params, error_message = conditionMessage(out),
                         timestamp = Sys.time(), duration_ms = NULL,
                         annotation_stamp = make_annotation_stamp(state)))
        push_message(state, sprintf("Trajectory failed: %s", conditionMessage(out)), "error")
      } else {
        dur <- as.integer(round((proc.time()[["elapsed"]] - t0) * 1000))
        set_tr_slot(list(status = "completed", results = out,
                         params = params, error_message = NULL,
                         timestamp = Sys.time(), duration_ms = dur,
                         annotation_stamp = make_annotation_stamp(state)))
        push_message(state, sprintf(
          "Trajectory done: source=%s, %d cells (%d ms). Exploratory only.",
          out$source, length(out$pseudotime), dur), "success")
      }
    })

    # ---- Apply to dataset metadata -------------------------------------
    # Mirrors mod_annotation's "Apply to dataset metadata" pattern. The
    # canonical result stays in state$analysis_results$trajectory; this
    # button just bakes it into a dated metadata column so other modules
    # can see it without coupling.
    shiny::observeEvent(input$apply_to_metadata, {
      tr <- tr_slot()
      if (is.null(tr) || !identical(tr$status, "completed")) {
        push_message(state,
          "No completed trajectory result to apply. Run pseudotime first.",
          "warning")
        return()
      }
      ds <- state$active_dataset; if (is.null(ds)) return()
      n_bins <- if (isTRUE(input$apply_with_bins))
                  as.integer(input$apply_bins %||% 10L)
                else 0L
      ds2 <- tryCatch(
        apply_pseudotime_to_dataset(ds, tr, bins = n_bins),
        error = function(e) e)
      if (inherits(ds2, "error")) {
        push_message(state, conditionMessage(ds2), "error")
        return()
      }
      state$active_dataset <- ds2
      n_new <- length(setdiff(ds2$metadata_fields, ds$metadata_fields))
      push_message(state, sprintf(
        "Applied pseudotime ('%s') as %d new metadata column%s.",
        tr$results$source %||% tr$params$source %||% "?",
        n_new, if (n_new == 1L) "" else "s"), "success")
    })

    # ---- Status banner --------------------------------------------------
    output$status_banner <- shiny::renderUI({
      tr <- tr_slot()
      if (is.null(tr))
        return(status_banner(
          shiny::span("Not run yet. Configure controls and click ",
                      shiny::tags$em("Run / generate pseudotime"), "."),
          tone = "idle"))
      tone <- switch(tr$status, running = "running", completed = "success",
                     failed = "danger", "idle")
      txt <- switch(tr$status,
        running   = "Running...",
        completed = {
          r <- tr$results
          base <- sprintf("Completed: source=%s", r$source %||% "?")
          extra <- if (identical(r$source, "metadata")) {
            sprintf(" (%s)", r$metadata_field %||% "?")
          } else if (!is.null(r$reduction_used) && !is.na(r$reduction_used)) {
            sprintf(", root=%s='%s' on %s",
                    r$root_field %||% "?", r$root_group %||% "?",
                    r$reduction_used)
          } else ""
          lin <- if (!is.null(r$n_lineages) && r$n_lineages > 1L)
                   sprintf(", %d lineages", r$n_lineages) else ""
          sprintf("%s%s%s. %d cells, %d ms.",
                  base, extra, lin,
                  length(r$pseudotime), tr$duration_ms %||% 0L)
        },
        failed    = sprintf("Failed: %s", tr$error_message %||% "(unknown)"),
        "")
      status_banner(txt, tone = tone)
    })

    # ---- Summary --------------------------------------------------------
    output$summary <- shiny::renderPrint({
      tr <- tr_slot()
      if (is.null(tr) || is.null(tr$results)) {
        return(list(status = "no pseudotime yet"))
      }
      r <- tr$results
      pseudotime_summary(
        r$pseudotime,
        source         = r$source,
        root_field     = r$root_field,
        root_group     = r$root_group,
        metadata_field = r$metadata_field,
        reduction_used = r$reduction_used)
    })

    # ---- Embedding by pseudotime ---------------------------------------
    output$pt_plot_warning <- shiny::renderUI({
      ds <- state$active_dataset
      if (is.null(ds)) return(friendly_warning("No dataset loaded."))
      if (!has_trajectory_results(state))
        return(helper_text("Run pseudotime to populate this panel."))
      red <- tr_slot()$results$reduction_used %||% input$reduction
      if (is.null(get_embedding(ds, red)))
        return(friendly_warning(sprintf("Reduction '%s' is not available.", red %||% "")))
      NULL
    })
    output$pt_plot <- shiny::renderPlot({
      ds <- state$active_dataset
      if (is.null(ds) || !has_trajectory_results(state)) return(NULL)
      tr <- tr_slot()$results
      red <- tr$reduction_used %||% input$reduction
      emb <- get_embedding(ds, red); if (is.null(emb)) return(NULL)
      pt <- tr$pseudotime
      plot_embedding_continuous(
        emb, pt,
        title = sprintf("%s | pseudotime", red),
        xlab  = paste0(red, "_1"), ylab = paste0(red, "_2"),
        legend_title = "pt")
    })

    # ---- Embedding by root group (or metadata-source field) ------------
    output$root_plot_warning <- shiny::renderUI({
      ds <- state$active_dataset
      if (is.null(ds)) return(friendly_warning("No dataset loaded."))
      tr <- tr_slot()
      if (is.null(tr) || !identical(tr$status, "completed"))
        return(helper_text("Run pseudotime to populate this panel."))
      NULL
    })
    output$root_plot <- shiny::renderPlot({
      ds <- state$active_dataset
      if (is.null(ds) || !has_trajectory_results(state)) return(NULL)
      tr <- tr_slot()$results
      red <- (tr$reduction_used %||% input$reduction) %||% ds$default_reduction
      emb <- get_embedding(ds, red); if (is.null(emb)) return(NULL)
      if (identical(tr$source, "metadata")) {
        v <- get_metadata(ds, tr$metadata_field)
        if (is.null(v)) return(NULL)
        plot_embedding_continuous(
          emb, v,
          title = sprintf("Source = %s", tr$metadata_field),
          xlab  = paste0(red, "_1"), ylab = paste0(red, "_2"),
          legend_title = tr$metadata_field)
      } else {
        # Every other backend uses a root group within a cluster field.
        rf <- if (!is.na(tr$root_field)) get_metadata(ds, tr$root_field) else NULL
        if (is.null(rf) || is.na(tr$root_group)) return(NULL)
        is_root <- ifelse(as.character(rf) == as.character(tr$root_group),
                          "root", "other")
        plot_embedding_categorical(
          emb, is_root,
          title = sprintf("Root = %s == '%s'", tr$root_field, tr$root_group),
          xlab  = paste0(red, "_1"), ylab = paste0(red, "_2"))
      }
    })

    # ---- Gene picker (synced with shared state$selected_gene) ----------
    # Server-side gene picker. See R/ui_components.R for rationale.
    output$gene_picker_ui <- shiny::renderUI({
      gene_picker_input(ns("gene"), label = NULL,
                        selected = state$selected_gene)
    })
    shiny::observe({
      ds <- state$active_dataset; shiny::req(ds)
      genes <- available_genes(ds)
      current <- state$selected_gene %||% genes[1]
      if (!isTRUE(current %in% genes)) current <- genes[1]
      update_gene_picker(session, "gene",
                         choices = genes, selected = current)
    })
    shiny::observeEvent(input$gene, {
      if (!is.null(input$gene) && nzchar(input$gene)) {
        state$selected_gene <- input$gene
      }
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$send_to_explorer, {
      g <- input$gene
      if (is.null(g) || !nzchar(g)) return()
      state$selected_gene <- g
      push_message(state, sprintf("Sent '%s' to the Explorer FeaturePlot.", g), "success")
    })
    output$send_status <- shiny::renderText({
      g <- state$selected_gene
      if (is.null(g) || !nzchar(g)) "" else sprintf("Selected: %s", g)
    })

    # ---- Gene vs pseudotime trend --------------------------------------
    output$trend_title <- shiny::renderText({
      sprintf("%s expression along pseudotime", state$selected_gene %||% "(none)")
    })
    output$trend_warning <- shiny::renderUI({
      ds <- state$active_dataset
      if (is.null(ds)) return(friendly_warning("No dataset loaded."))
      if (!has_trajectory_results(state))
        return(helper_text("Run pseudotime first."))
      g <- state$selected_gene
      if (!validate_gene(ds, g))
        return(friendly_warning(sprintf("Gene '%s' is not available.", g %||% "")))
      NULL
    })
    output$trend_plot <- shiny::renderPlot({
      ds <- state$active_dataset
      if (is.null(ds) || !has_trajectory_results(state)) return(NULL)
      g <- state$selected_gene
      if (!validate_gene(ds, g)) return(NULL)
      pt <- tr_slot()$results$pseudotime
      # Raw expression on purpose -- trajectory module never consumes
      # smoothed values silently. (Future: optional toggle.)
      expr <- get_gene_expression(ds, g)
      if (is.null(expr) || is.null(pt) || length(pt) != length(expr)) return(NULL)
      plot_gene_vs_pseudotime(pt, expr, gene_name = g, n_bins = 25L)
    })
  })
}
