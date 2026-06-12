# ============================================================================
# Module: Dataset Overview  (ENABLED)
# ----------------------------------------------------------------------------
# Smallest reference module. Shows top-level summary information about the
# active dataset. Read-only with respect to app state.
# ============================================================================

mod_dataset_overview_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h2("Dataset Overview"),
    shiny::uiOutput(ns("summary")),
    shiny::hr(),
    shiny::h4("Available assays"),
    shiny::verbatimTextOutput(ns("assays")),
    shiny::h4("Available reductions"),
    shiny::verbatimTextOutput(ns("reductions")),
    shiny::h4("Cell metadata fields"),
    shiny::verbatimTextOutput(ns("meta_fields"))
  )
}

mod_dataset_overview_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    output$summary <- shiny::renderUI({
      ds <- state$active_dataset
      shiny::tagList(
        shiny::tags$p(shiny::strong("Name: "),   ds$name),
        shiny::tags$p(shiny::strong("Source: "), ds$source),
        shiny::tags$p(shiny::strong("Cells: "),  format(ds$n_cells, big.mark = ",")),
        shiny::tags$p(shiny::strong("Genes: "),  format(ds$n_genes, big.mark = ","))
      )
    })
    output$assays      <- shiny::renderPrint(state$active_dataset$assays)
    output$reductions  <- shiny::renderPrint(state$active_dataset$reductions)
    output$meta_fields <- shiny::renderPrint(state$active_dataset$metadata_fields)
  })
}
