# ============================================================================
# App-level dataset loading orchestration
# ----------------------------------------------------------------------------
# Thin wrapper around `load_dataset()` / `load_census()` used by the Shiny
# sidebar. Kept as pure/testable helpers so server wiring stays minimal.
# ============================================================================

#' Map a UI source choice to the `source` argument for `load_dataset()`.
#'
#' @param source_choice one of "auto", "seurat", "anndata", "10x", "census"
#' @return `NULL` for auto-detect, otherwise the explicit source string
resolve_dataset_source <- function(source_choice) {
  choice <- trimws(as.character(source_choice %||% "auto"))
  if (!nzchar(choice) || identical(choice, "auto")) return(NULL)
  choice
}

#' Attempt to browse Census studies without mutating app state.
#'
#' @param organism Census organism name
#' @return list(ok, studies, error)
try_load_census_studies <- function(organism) {
  org <- trimws(as.character(organism %||% CENSUS_ORGANISMS[1]))
  if (!(org %in% CENSUS_ORGANISMS)) {
    return(list(
      ok = FALSE, studies = NULL,
      error = sprintf(
        "Unsupported Census organism '%s' (expected one of: %s).",
        org, paste(CENSUS_ORGANISMS, collapse = ", ")
      )
    ))
  }
  studies <- tryCatch(
    load_census_studies(org),
    error = function(e) structure(list(error = conditionMessage(e)),
                                  class = "dataset_load_error")
  )
  if (inherits(studies, "dataset_load_error")) {
    return(list(ok = FALSE, studies = NULL, error = studies$error))
  }
  list(ok = TRUE, studies = studies, error = NULL)
}

#' Attempt to load a Census slice without mutating app state.
#'
#' @param census list with `organism`, `obs_value_filter`, optional
#'               `var_value_filter` (see `normalize_census_query()`)
#' @return list(ok, dataset, error) — same shape as `try_load_dataset()`
try_load_census <- function(census) {
  query <- normalize_census_query(
    organism         = census$organism,
    obs_value_filter = census$obs_value_filter,
    var_value_filter = census$var_value_filter
  )
  if (!nzchar(query$obs_value_filter)) {
    return(list(
      ok = FALSE, dataset = NULL,
      error = paste0(
        "Census query requires a cell filter ",
        "(e.g. cell_type == 'B cell')."
      )
    ))
  }
  if (!(query$organism %in% CENSUS_ORGANISMS)) {
    return(list(
      ok = FALSE, dataset = NULL,
      error = sprintf(
        "Unsupported Census organism '%s' (expected one of: %s).",
        query$organism, paste(CENSUS_ORGANISMS, collapse = ", ")
      )
    ))
  }

  ds <- tryCatch(
    load_census(
      organism         = query$organism,
      obs_value_filter = query$obs_value_filter,
      var_value_filter = query$var_value_filter
    ),
    error = function(e) structure(list(error = conditionMessage(e)),
                                  class = "dataset_load_error")
  )
  if (inherits(ds, "dataset_load_error")) {
    return(list(ok = FALSE, dataset = NULL, error = ds$error))
  }
  list(ok = TRUE, dataset = ds, error = NULL)
}

#' Attempt to load a dataset from disk without mutating app state.
#'
#' @param path           file or directory path from the UI
#' @param source_choice  UI source selector ("auto" or an explicit source)
#' @param census         optional list for Census queries when
#'                       `source_choice = "census"`
#' @return list(ok = logical(1), dataset = dataset or NULL,
#'              error = character(1) or NULL)
try_load_dataset <- function(path, source_choice = "auto", census = NULL) {
  source <- resolve_dataset_source(source_choice)
  if (identical(source, "census") || is_census_path(path)) {
    return(try_load_census(census %||% list()))
  }

  path <- trimws(as.character(path %||% ""))
  if (!nzchar(path)) {
    return(list(ok = FALSE, dataset = NULL,
                error = "Please enter a file or directory path."))
  }
  ds <- tryCatch(
    load_dataset(path, source = source),
    error = function(e) structure(list(error = conditionMessage(e)),
                                  class = "dataset_load_error")
  )
  if (inherits(ds, "dataset_load_error")) {
    return(list(ok = FALSE, dataset = NULL, error = ds$error))
  }
  list(ok = TRUE, dataset = ds, error = NULL)
}

#' Load a dataset into app state, preserving the active dataset on failure.
#'
#' On success, delegates to `set_active_dataset()`. On failure, pushes an
#' error message to the shared log and leaves `state$active_dataset` unchanged.
#'
#' @param census optional Census query list when loading from Census
#' @return list(ok = logical(1), dataset = dataset or NULL,
#'              error = character(1) or NULL)
app_load_dataset <- function(state, path, source_choice = "auto", census = NULL) {
  result <- try_load_dataset(path, source_choice, census = census)
  if (!isTRUE(result$ok)) {
    push_message(state,
                 sprintf("Dataset load failed: %s", result$error),
                 "error")
    return(result)
  }
  set_active_dataset(state, result$dataset)
  result
}
