# ============================================================================
# Module: Dataset Overview  (ENABLED)
# ----------------------------------------------------------------------------
# Smallest reference module. Shows top-level summary information about the
# active dataset. Read-only with respect to app state.
#
# UI composition uses the shared primitives in `R/ui_components.R`
# (`page_header`, `metric_card`, `app_card`, `deflist`, `helper_text`).
# ============================================================================

mod_dataset_overview_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    page_header(
      eyebrow = "Overview",
      title   = "Dataset Overview",
      lede    = paste("Top-level summary of the active dataset. All other",
                      "modules consume the same fields shown below.")
    ),
    shiny::uiOutput(ns("metrics")),
    shiny::fluidRow(
      shiny::column(6,
        app_card(
          title   = "Provenance",
          caption = "loader-reported source",
          shiny::uiOutput(ns("provenance"))
        )
      ),
      shiny::column(6,
        app_card(
          title   = "Available assays",
          caption = "expression matrices",
          shiny::uiOutput(ns("assays"))
        )
      )
    ),
    shiny::fluidRow(
      shiny::column(6,
        app_card(
          title   = "Available reductions",
          caption = "2D embeddings",
          shiny::uiOutput(ns("reductions"))
        )
      ),
      shiny::column(6,
        app_card(
          title   = "Cell metadata fields",
          caption = "columns of cell_data",
          shiny::uiOutput(ns("meta_fields"))
        )
      )
    )
  )
}

# Internal: render a character vector as a compact styled list. Empty / NULL
# vectors render as a muted "(none)" so cards don't go blank.
.dso_field_list <- function(items) {
  if (is.null(items) || !length(items))
    return(helper_text("(none reported by the loader)"))
  shiny::div(
    class = "gene-list",
    paste(items, collapse = ", ")
  )
}

mod_dataset_overview_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {

    output$metrics <- shiny::renderUI({
      ds <- state$active_dataset
      if (is.null(ds)) return(NULL)
      metric_grid(
        metric_card("Cells",      format(ds$n_cells, big.mark = ",")),
        metric_card("Genes",      format(ds$n_genes, big.mark = ",")),
        metric_card("Assays",     length(ds$assays),
                    hint = if (length(ds$assays))
                             paste(ds$assays, collapse = ", ")),
        metric_card("Reductions", length(ds$reductions),
                    hint = if (length(ds$reductions))
                             paste(ds$reductions, collapse = ", "))
      )
    })

    output$provenance <- shiny::renderUI({
      ds <- state$active_dataset
      if (is.null(ds)) return(NULL)
      deflist(list(
        Name   = shiny::tags$strong(ds$name),
        Source = shiny::tags$code(ds$source %||% "(unknown)"),
        Cells  = shiny::tags$span(class = "sce-tabular",
                                  format(ds$n_cells, big.mark = ",")),
        Genes  = shiny::tags$span(class = "sce-tabular",
                                  format(ds$n_genes, big.mark = ","))
      ))
    })

    output$assays      <- shiny::renderUI({
      .dso_field_list(state$active_dataset$assays)
    })
    output$reductions  <- shiny::renderUI({
      .dso_field_list(state$active_dataset$reductions)
    })
    output$meta_fields <- shiny::renderUI({
      .dso_field_list(state$active_dataset$metadata_fields)
    })
  })
}
