# ============================================================================
# CELLxGENE Census loader
# ----------------------------------------------------------------------------
# Fetches a filtered Census slice via cellxgene.census::get_seurat() and maps
# it into the app dataset schema through the existing Seurat converter.
#
# Virtual path sentinel: `/cellxgene-census` (see detect_source()).
# Query parameters are passed as named arguments to load_census() / load_dataset().
# ============================================================================

#' Default cell-metadata columns fetched from Census obs.
CENSUS_DEFAULT_OBS_COLUMNS <- c(
  "assay", "cell_type", "tissue", "tissue_general", "disease",
  "sex", "development_stage", "suspension_type"
)

#' Supported Census organisms for the UI and loader.
CENSUS_ORGANISMS <- c("Homo sapiens", "Mus musculus")

#' TRUE when `path` is the Census virtual sentinel.
is_census_path <- function(path) {
  path <- trimws(as.character(path %||% ""))
  identical(path, "/cellxgene-census") ||
    grepl("^/cellxgene-census/", path)
}

#' Build a short display name for a Census query.
census_dataset_name <- function(organism, obs_value_filter) {
  org_slug <- gsub("\\s+", "_", trimws(as.character(organism %||% "Homo sapiens")))
  filter_slug <- gsub("[^a-zA-Z0-9]+", "_", trimws(as.character(obs_value_filter)))
  filter_slug <- substr(filter_slug, 1, 48)
  paste0("census_", org_slug, "_", filter_slug)
}

#' Normalise optional Census query arguments from the UI or programmatic calls.
#'
#' @return named list(organism, obs_value_filter, var_value_filter)
normalize_census_query <- function(organism = "Homo sapiens",
                                   obs_value_filter = NULL,
                                   var_value_filter = NULL) {
  obs_value_filter <- trimws(as.character(obs_value_filter %||% ""))
  var_value_filter <- trimws(as.character(var_value_filter %||% ""))
  if (!nzchar(var_value_filter)) var_value_filter <- NULL
  list(
    organism         = trimws(as.character(organism %||% CENSUS_ORGANISMS[1])),
    obs_value_filter = obs_value_filter,
    var_value_filter = var_value_filter
  )
}

#' Load Census study metadata for browsing in the app.
#'
#' Returns one row per study / dataset summary where available. This is a
#' lightweight metadata query only; it does not materialise expression.
#'
#' @param organism       "Homo sapiens" or "Mus musculus"
#' @param census_version Census release, e.g. "stable"
#' @return data.frame with at least `dataset_id` and best-effort display
#'   columns for title / collection / cell counts when the installed Census API
#'   exposes them.
load_census_studies <- function(organism = "Homo sapiens",
                                census_version = "stable") {
  organism <- match.arg(trimws(as.character(organism %||% CENSUS_ORGANISMS[1])),
                        CENSUS_ORGANISMS)
  require_optional(
    c("cellxgene.census"),
    feature = "CELLxGENE Census study browsing",
    source  = "GitHub",
    repo    = c("cellxgene.census" = "chanzuckerberg/cellxgene-census")
  )

  census <- cellxgene.census::open_soma(census_version = census_version)
  on.exit(census$close(), add = TRUE)

  exp <- .census_experiment_handle(census, organism)
  studies <- .census_table_collect(exp$ms$datasets)
  .normalize_census_studies(studies)
}

#' Build a conservative `obs_value_filter` for a selected study row.
#'
#' @param studies data.frame returned by `load_census_studies()`
#' @param row_idx selected row index (1-based)
#' @return character(1) SOMA filter string
census_study_filter <- function(studies, row_idx) {
  stopifnot(is.data.frame(studies), length(row_idx) == 1L)
  row_idx <- as.integer(row_idx)
  if (is.na(row_idx) || row_idx < 1L || row_idx > nrow(studies)) {
    stop("Census study selection is out of range.", call. = FALSE)
  }
  row <- studies[row_idx, , drop = FALSE]
  dataset_id_col <- intersect(c("dataset_id", "soma_joinid", "id"), names(row))
  if (length(dataset_id_col) == 0L) {
    stop("Census study metadata did not include a usable dataset id column.",
         call. = FALSE)
  }
  dataset_id_col <- dataset_id_col[1]
  dataset_id <- as.character(row[[dataset_id_col]][1])
  if (!nzchar(dataset_id)) {
    stop("Selected Census study has an empty dataset id.", call. = FALSE)
  }
  sprintf("dataset_id == '%s'", gsub("'", "\\\\'", dataset_id))
}

#' Load a filtered slice from CZ CELLxGENE Census.
#'
#' Requires optional packages `cellxgene.census` and `SeuratObject`. The full
#' Census atlas is too large to load wholesale; callers must supply an
#' `obs_value_filter` (and optionally `var_value_filter`) to define a slice.
#'
#' @param organism          "Homo sapiens" or "Mus musculus"
#' @param obs_value_filter  SOMA value filter over obs columns (required)
#' @param var_value_filter  optional SOMA value filter over var columns
#' @param obs_column_names  obs columns to retain in cell metadata
#' @param census_version    Census release, e.g. "stable"
#' @param name              display name; auto-generated when NULL
#' @param obsm_layers       load dimensional reductions from Census obsm
#' @return a dataset list matching `dataset_schema()` with `source = "census"`
load_census <- function(organism = "Homo sapiens",
                        obs_value_filter,
                        var_value_filter = NULL,
                        obs_column_names = NULL,
                        census_version = "stable",
                        name = NULL,
                        obsm_layers = TRUE) {
  query <- normalize_census_query(organism, obs_value_filter, var_value_filter)
  organism <- match.arg(query$organism, CENSUS_ORGANISMS)
  if (!nzchar(query$obs_value_filter)) {
    stop(
      "Census loader: obs_value_filter is required. ",
      "The full Census atlas is too large to load; provide a cell filter ",
      "(e.g. \"cell_type == 'B cell'\").",
      call. = FALSE
    )
  }

  require_optional(
    c("cellxgene.census", "SeuratObject"),
    feature = "CELLxGENE Census loading",
    source  = "GitHub",
    repo    = c("cellxgene.census" = "chanzuckerberg/cellxgene-census")
  )

  obs_column_names <- obs_column_names %||% CENSUS_DEFAULT_OBS_COLUMNS
  ds_name <- name %||% census_dataset_name(organism, query$obs_value_filter)

  census <- cellxgene.census::open_soma(census_version = census_version)
  on.exit(census$close(), add = TRUE)

  obj <- cellxgene.census::get_seurat(
    census           = census,
    organism         = organism,
    obs_value_filter = query$obs_value_filter,
    var_value_filter = query$var_value_filter,
    obs_column_names = obs_column_names,
    obsm_layers      = obsm_layers
  )

  ds <- .seurat_to_dataset(obj, name = ds_name)
  ds$source <- "census"
  ds
}

.census_experiment_handle <- function(census, organism) {
  if (!is.list(census) && !is.environment(census) && !isS4(census)) {
    stop("Census handle is not queryable.", call. = FALSE)
  }
  exp <- census[[organism]]
  if (is.null(exp)) {
    stop("Census organism handle not found: ", organism, call. = FALSE)
  }
  exp
}

.census_table_collect <- function(tbl) {
  if (is.null(tbl)) {
    stop("Census metadata table was not available.", call. = FALSE)
  }
  if (is.function(tbl$read)) {
    res <- tbl$read()
    if (is.function(res$concat)) res <- res$concat()
    if (is.function(res$to_pandas)) {
      df <- reticulate::py_to_r(res$to_pandas())
      return(as.data.frame(df, stringsAsFactors = FALSE))
    }
    if (is.function(res$to_table)) {
      res <- res$to_table()
    }
    if (is.function(res$to_pandas)) {
      df <- reticulate::py_to_r(res$to_pandas())
      return(as.data.frame(df, stringsAsFactors = FALSE))
    }
  }
  if (is.data.frame(tbl)) return(tbl)
  stop("Could not materialise Census metadata table with the installed API.",
       call. = FALSE)
}

.normalize_census_studies <- function(studies) {
  studies <- as.data.frame(studies, stringsAsFactors = FALSE)
  if (nrow(studies) == 0L) return(studies)

  rename_first <- function(from, to) {
    hit <- intersect(from, names(studies))
    if (length(hit) > 0L && !(to %in% names(studies))) {
      names(studies)[match(hit[1], names(studies))] <<- to
    }
  }

  rename_first(c("dataset_id", "id"), "dataset_id")
  rename_first(c("dataset_title", "title", "name"), "dataset_title")
  rename_first(c("collection_name", "collection", "collection_title"), "collection_name")
  rename_first(c("dataset_total_cell_count", "cell_count", "n_cells"), "cell_count")
  rename_first(c("organism", "dataset_h5ad_path"), "organism")

  if (!("dataset_title" %in% names(studies))) {
    studies$dataset_title <- studies$dataset_id %||% seq_len(nrow(studies))
  }
  studies
}
