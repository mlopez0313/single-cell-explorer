# ============================================================================
# Sidebar navigation
# ----------------------------------------------------------------------------
# Renders the module list, grouped by category. Each entry is a button that
# sets `state$active_module`. Disabled modules are shown but visually muted
# and do not change the active module on click.
# ============================================================================

#' Render a lightweight Census study table for the sidebar.
render_census_study_table <- function(studies) {
  if (is.null(studies) || !is.data.frame(studies) || nrow(studies) == 0L) {
    return(shiny::div(style = "font-size:12px; color:#888; margin-top:6px;",
                      "No Census studies available for the current organism."))
  }
  title_col <- intersect(c("dataset_title", "title", "name"), names(studies))[1]
  collection_col <- intersect(c("collection_name", "collection"), names(studies))[1]
  cell_count_col <- intersect(c("cell_count", "dataset_total_cell_count", "n_cells"), names(studies))[1]
  rows <- lapply(seq_len(nrow(studies)), function(i) {
    row <- studies[i, , drop = FALSE]
    label <- as.character(row[[title_col]][1] %||% row[[1]][1])
    meta <- c()
    if (!is.na(collection_col) && nzchar(collection_col)) meta <- c(meta, as.character(row[[collection_col]][1]))
    if (!is.na(cell_count_col) && nzchar(cell_count_col)) meta <- c(meta, paste0("cells: ", as.character(row[[cell_count_col]][1])))
    shiny::tags$label(
      style = "display:block; font-size:12px; margin:6px 0; padding:6px; border:1px solid #eee; border-radius:4px; cursor:pointer;",
      shiny::tags$input(type = "radio", name = "census_study_row", value = as.character(i), style = "margin-right:6px;"),
      shiny::tags$strong(label),
      if (length(meta) > 0L) shiny::tags$div(style = "color:#666; margin-top:2px;", paste(meta, collapse = " • "))
    )
  })
  shiny::tags$div(style = "max-height:180px; overflow-y:auto; margin-top:6px;", rows)
}

#' Static sidebar shell (populated by the server-side `sidebar_server`).
sidebar_ui <- function() {
  shiny::tagList(
    shiny::div(class = "sidebar-section",
      shiny::h5("Dataset"),
      shiny::actionButton("load_mock_dataset", "Load mock dataset",
                          class = "btn btn-outline-secondary btn-sm",
                          width = "100%"),
      shiny::hr(style = "margin:10px 0;"),
      shiny::tags$label(`for` = "dataset_source",
                        style = "font-size:12px; margin-bottom:2px;",
                        "Source"),
      shiny::selectInput(
        "dataset_source",
        label = NULL,
        choices = c(
          "Auto-detect" = "auto",
          "Seurat (.rds)" = "seurat",
          "AnnData (.h5ad)" = "anndata",
          "10x directory" = "10x",
          "CELLxGENE Census (experimental)" = "census"
        ),
        selected = "auto",
        width = "100%"
      ),
      shiny::conditionalPanel(
        condition = "input.dataset_source != 'census'",
        shiny::tags$label(`for` = "dataset_path",
                          style = "font-size:12px; margin-bottom:2px;",
                          "Path"),
        shiny::textInput(
          "dataset_path",
          label = NULL,
          placeholder = "/path/to/dataset.rds",
          width = "100%"
        )
      ),
      shiny::conditionalPanel(
        condition = "input.dataset_source == 'census'",
        shiny::tags$div(
          style = "font-size:11px; color:#888; margin-bottom:8px;",
          "Experimental: requires optional dependencies and network access."
        ),
        shiny::tags$label(`for` = "census_organism",
                          style = "font-size:12px; margin-bottom:2px;",
                          "Organism"),
        shiny::selectInput(
          "census_organism",
          label = NULL,
          choices = stats::setNames(CENSUS_ORGANISMS, CENSUS_ORGANISMS),
          selected = CENSUS_ORGANISMS[1],
          width = "100%"
        ),
        shiny::actionButton("browse_census_studies", "Browse studies",
                            class = "btn btn-outline-secondary btn-sm",
                            width = "100%"),
        shiny::div(style = "margin-top:6px;", shiny::uiOutput("census_study_browser")),
        shiny::tags$label(`for` = "census_obs_filter",
                          style = "font-size:12px; margin:8px 0 2px 0;",
                          "Cell filter"),
        shiny::textInput(
          "census_obs_filter",
          label = NULL,
          placeholder = "Leave blank to use the selected study",
          width = "100%"
        ),
        shiny::tags$label(`for` = "census_var_filter",
                          style = "font-size:12px; margin-bottom:2px;",
                          "Gene filter (optional)"),
        shiny::textInput(
          "census_var_filter",
          label = NULL,
          placeholder = "feature_name == 'CD3D'",
          width = "100%"
        )
      ),
      shiny::uiOutput("dataset_path_help"),
      shiny::actionButton("load_local_dataset", "Load dataset",
                          class = "btn btn-primary btn-sm", width = "100%"),
      shiny::div(style = "margin-top:8px;",
                 shiny::uiOutput("sidebar_dataset_status"))
    ),
    shiny::hr(),
    shiny::div(class = "sidebar-section",
      shiny::h5("Modules"),
      shiny::uiOutput("sidebar_modules")
    )
  )
}

#' Helper text for the dataset path input, keyed off the selected source.
dataset_path_help_text <- function(source_choice) {
  switch(source_choice %||% "auto",
    "auto" = paste0(
      "Auto-detect: .rds (Seurat), .h5ad (AnnData), or a Cellranger ",
      "feature-barcode matrix directory."
    ),
    "seurat" = "File path to a Seurat object saved with saveRDS() (.rds).",
    "anndata" = "File path to an AnnData (.h5ad) file.",
    "10x" = paste0(
      "Directory path containing matrix.mtx, barcodes.tsv, and ",
      "features.tsv (or .gz variants)."
    ),
    "census" = paste0(
      "Experimental remote loader: browse Census studies, then load the ",
      "selected study or provide a custom SOMA cell filter. Requires ",
      "cellxgene.census + SeuratObject."
    ),
    "Enter a file or directory path on the server."
  )
}

#' Build the sidebar dataset status block (idle / loading / loaded / error).
sidebar_dataset_status_ui <- function(active_dataset, load_phase, load_error) {
  if (identical(load_phase, "loading")) {
    return(shiny::div(
      style = "font-size:12px; color:#1565c0;",
      "Loading dataset..."
    ))
  }

  pieces <- list()
  if (is.null(active_dataset)) {
    pieces[[length(pieces) + 1L]] <- shiny::div(
      style = "font-size:12px; color:#888;",
      "No dataset loaded."
    )
  } else {
    src <- active_dataset$source %||% "unknown"
    pieces[[length(pieces) + 1L]] <- shiny::div(
      style = "font-size:12px; color:#2e7d32;",
      sprintf("Loaded: %s (%s)", active_dataset$name, src)
    )
  }

  if (nzchar(load_error %||% "")) {
    pieces[[length(pieces) + 1L]] <- shiny::div(
      style = "font-size:12px; color:#c62828; margin-top:4px;",
      sprintf("Load failed: %s", load_error)
    )
  }

  shiny::tagList(pieces)
}

#' Server logic that draws the module list and handles click events.
#'
#' @param input,output,session  Shiny app-level objects
#' @param state                 shared app state
sidebar_server <- function(input, output, session, state) {
  load_phase <- shiny::reactiveVal("idle")
  load_error <- shiny::reactiveVal(NULL)
  census_browse_phase <- shiny::reactiveVal("idle")
  census_browse_error <- shiny::reactiveVal(NULL)
  census_studies <- shiny::reactiveVal(NULL)

  output$dataset_path_help <- shiny::renderUI({
    shiny::tags$div(
      style = "font-size:11px; color:#888; margin:4px 0 8px 0;",
      dataset_path_help_text(input$dataset_source)
    )
  })

  output$sidebar_dataset_status <- shiny::renderUI({
    sidebar_dataset_status_ui(state$active_dataset, load_phase(), load_error())
  })

  output$census_study_browser <- shiny::renderUI({
    if (!identical(input$dataset_source, "census")) return(NULL)
    if (identical(census_browse_phase(), "loading")) {
      return(shiny::div(style = "font-size:12px; color:#1565c0; margin-top:6px;",
                        "Loading Census studies..."))
    }
    if (nzchar(census_browse_error() %||% "")) {
      return(shiny::div(style = "font-size:12px; color:#c62828; margin-top:6px;",
                        census_browse_error()))
    }
    render_census_study_table(census_studies())
  })

  shiny::observeEvent(input$load_mock_dataset, {
    load_phase("idle")
    load_error(NULL)
    set_active_dataset(state, mock_dataset())
  }, ignoreInit = TRUE)

  shiny::observeEvent(input$browse_census_studies, {
    census_browse_phase("loading")
    census_browse_error(NULL)
    result <- try_load_census_studies(input$census_organism)
    if (isTRUE(result$ok)) {
      census_studies(result$studies)
      census_browse_phase("idle")
    } else {
      census_studies(NULL)
      census_browse_phase("error")
      census_browse_error(result$error)
    }
  }, ignoreInit = TRUE)

  shiny::observeEvent(input$load_local_dataset, {
    load_phase("loading")
    load_error(NULL)

    if (identical(input$dataset_source, "census")) {
      study_filter <- NULL
      studies <- census_studies()
      sel <- suppressWarnings(as.integer(input$census_study_row %||% NA_character_))
      if (!is.null(studies) && is.data.frame(studies) && nrow(studies) > 0L && !is.na(sel)) {
        study_filter <- tryCatch(census_study_filter(studies, sel), error = function(e) e)
        if (inherits(study_filter, "error")) {
          load_phase("error")
          load_error(conditionMessage(study_filter))
          push_message(state, sprintf("Dataset load failed: %s", conditionMessage(study_filter)), "error")
          return()
        }
      }
      custom_filter <- trimws(as.character(input$census_obs_filter %||% ""))
      obs_filter <- if (nzchar(custom_filter) && nzchar(study_filter %||% "")) {
        sprintf("(%s) and (%s)", study_filter, custom_filter)
      } else if (nzchar(custom_filter)) {
        custom_filter
      } else {
        study_filter
      }
      result <- app_load_dataset(
        state,
        path = "/cellxgene-census",
        source_choice = "census",
        census = list(
          organism = input$census_organism,
          obs_value_filter = obs_filter,
          var_value_filter = input$census_var_filter
        )
      )
    } else {
      result <- app_load_dataset(state, input$dataset_path, input$dataset_source)
    }

    if (isTRUE(result$ok)) {
      load_phase("idle")
      load_error(NULL)
    } else {
      load_phase("error")
      load_error(result$error)
    }
  }, ignoreInit = TRUE)

  output$sidebar_modules <- shiny::renderUI({
    grouped <- modules_by_category()
    groups <- lapply(names(grouped), function(cat) {
      mods <- grouped[[cat]]
      shiny::div(class = "sidebar-category",
        shiny::tags$div(
          style = "font-size:11px; text-transform:uppercase; letter-spacing:0.05em; color:#888; margin:12px 0 4px 0;",
          cat
        ),
        lapply(mods, function(m) module_button(m, state$active_module))
      )
    })
    shiny::tagList(groups)
  })

  # Wire up one observer per module id. Done lazily via lapply so future
  # modules added to the registry get observers automatically.
  shiny::isolate({
    for (m in module_registry()) {
      local({
        mod <- m
        if (!mod$enabled) return()
        shiny::observeEvent(input[[paste0("nav_", mod$id)]], {
          state$active_module <- mod$id
        }, ignoreInit = TRUE)
      })
    }
  })
}

#' Render one sidebar entry. Enabled -> actionButton; disabled -> muted div.
module_button <- function(mod, active_id) {
  is_active <- identical(mod$id, active_id)
  if (!mod$enabled) {
    return(shiny::div(
      class = "sidebar-module disabled",
      style = paste(
        "padding:6px 10px; margin:2px 0; border-radius:4px; color:#aaa;",
        "background:#f5f5f5; font-size:13px; cursor:not-allowed;",
        "display:flex; justify-content:space-between; align-items:center;"
      ),
      title = mod$description,
      shiny::span(mod$name),
      shiny::tags$span(style = "font-size:10px; background:#ddd; padding:1px 6px; border-radius:8px;",
                       "soon")
    ))
  }
  shiny::actionButton(
    inputId = paste0("nav_", mod$id),
    label   = mod$name,
    class   = if (is_active) "btn btn-primary btn-sm sidebar-module active"
              else            "btn btn-light    btn-sm sidebar-module",
    style   = "width:100%; text-align:left; margin:2px 0;",
    title   = mod$description
  )
}
