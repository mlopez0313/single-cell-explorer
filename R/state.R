# ============================================================================
# Shared app state
# ----------------------------------------------------------------------------
# A single reactiveValues object is created at app start and passed to every
# module's server function. Modules MUST treat this object as the single source
# of truth for cross-module selections (active dataset, active assay, etc.).
#
# Adding a new shared field?
#   1. Add a default below in `new_app_state()`.
#   2. Document it in docs/ADDING_MODULES.md.
#   3. Read it from modules via `state$your_field`; write with `state$your_field <- value`.
# ============================================================================

#' Create a fresh app state object.
#'
#' Notable fields:
#'   active_dataset           NULL or a dataset list (see R/dataset.R)
#'   active_module            character(1) module id currently shown in workspace
#'   selected_assay           character(1) e.g. "RNA", "SCT", "ADT"
#'   selected_reduction       character(1) e.g. "UMAP", "PCA", "tSNE"
#'   selected_metadata_field  character(1) column name from dataset metadata
#'   selected_gene            character(1) gene symbol to feature-plot
#'   selected_cells           character() vector of cell barcodes
#'
#'   --- Annotation system (versioned, multi-set) ---
#'   annotation_sets          named list of `annotation_result_v1` objects
#'                            (one per saved annotation set). Names are
#'                            set ids. Schema lives in R/annotation_schema.R.
#'   active_annotation_id     character(1) | NULL: id of the active set.
#'                            Downstream modules MUST read the active set
#'                            via `get_active_annotation(state)` -- never by
#'                            indexing `annotation_sets` directly.
#'   marker_registry          typed marker registry from R/marker_registry.R.
#'                            Loaded from `default_marker_registry()` on
#'                            dataset load. State-level so all modules /
#'                            engines query the same object.
#'
#'   display_mode_imputation  "raw" | "smoothed" -- visualization only.
#'   analysis_results         named list; one entry per analysis kind.
#'                            Each entry has shape:
#'                              list(status=..., results=..., params=...,
#'                                   error_message=..., timestamp=...,
#'                                   duration_ms=..., annotation_stamp=...)
#'                            `annotation_stamp` is added at the slot level
#'                            so compute helpers stay pure. See
#'                            `make_annotation_stamp()`.
#'                            Slots: $de, $pathway, $imputation, $trajectory,
#'                                   $markers, $regulons.
#'   messages                 list of list(level=, text=, time=) entries
#'
#' @return a `reactiveValues` object
new_app_state <- function() {
  shiny::reactiveValues(
    active_dataset           = NULL,
    active_module            = "dataset_overview",
    selected_assay           = NULL,
    selected_reduction       = NULL,
    selected_metadata_field  = NULL,
    selected_gene            = NULL,
    selected_cells           = character(),

    # Annotation system
    annotation_sets          = list(),
    active_annotation_id     = NULL,
    marker_registry          = NULL,

    # Visualisation / smoothing
    display_mode_imputation  = "raw",

    # Analysis result slots
    analysis_results         = list(),

    # App-wide messages
    messages                 = list()
  )
}

#' Default empty entry for a single analysis kind.
empty_analysis_result <- function() {
  list(
    status           = "not_run",
    results          = NULL,
    params           = NULL,
    error_message    = NULL,
    timestamp        = NULL,
    duration_ms      = NULL,
    annotation_stamp = NULL
  )
}

#' Append a message to the shared message log.
#'
#' @param state  reactiveValues from `new_app_state()`
#' @param text   character(1) message text
#' @param level  one of "info", "success", "warning", "error"
push_message <- function(state, text, level = c("info", "success", "warning", "error")) {
  level <- match.arg(level)
  entry <- list(level = level, text = text, time = Sys.time())
  state$messages <- c(state$messages, list(entry))
  invisible(entry)
}

#' Set the active dataset and reset dependent selections to sensible defaults.
#'
#' Loads `default_marker_registry()` into state and clears all annotation
#' sets, analysis results, smoothing display mode, and gene/cell selections.
#' If the (deprecated) `state$annotations` field carries content from an
#' older session, it is migrated into a single set named "Default".
set_active_dataset <- function(state, dataset) {
  state$active_dataset           <- dataset
  state$selected_assay           <- dataset$default_assay     %||% NA_character_
  state$selected_reduction       <- dataset$default_reduction %||% NA_character_
  state$selected_metadata_field  <- dataset$metadata_fields[1] %||% NA_character_
  state$selected_gene            <- dataset$genes[1]           %||% NA_character_
  state$selected_cells           <- character()

  # Reset visualisation toggles + analysis results.
  state$display_mode_imputation  <- "raw"
  state$analysis_results         <- list()

  # Reset / initialise annotation system.
  state$marker_registry          <- default_marker_registry()
  state$annotation_sets          <- list()
  state$active_annotation_id     <- NULL

  # Legacy migration: pre-multi-set sessions stored a data.frame on
  # `state$annotations`. Promote that into a Default set, then clear it.
  if (!is.null(state$annotations) && is.data.frame(state$annotations) &&
      nrow(state$annotations) > 0L &&
      all(c("cluster", "annotation") %in% names(state$annotations))) {
    legacy <- state$annotations
    cluster_field <- state$annotation_cluster_field %||% "cluster"
    labels_map <- setNames(
      lapply(seq_len(nrow(legacy)), function(i) {
        v <- legacy$annotation[i]
        if (is.null(v) || is.na(v) || !nzchar(trimws(as.character(v))))
          NA_character_
        else trimws(as.character(v))
      }),
      as.character(legacy$cluster)
    )
    legacy_set <- tryCatch(
      run_annotation_engine(
        engine_id     = "manual",
        dataset       = dataset,
        state         = state,
        params        = list(cluster_field = cluster_field, labels = labels_map),
        set_id        = "set_default_legacy",
        set_name      = "Default",
        description   = "Migrated from legacy state$annotations.",
        is_demo       = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(legacy_set)) {
      add_annotation_set(state, legacy_set)
      set_active_annotation(state, legacy_set$set_id)
      push_message(state, "Legacy annotations migrated to 'Default' set.", "info")
    }
  }
  state$annotations              <- NULL  # deprecated; cleared
  state$annotation_cluster_field <- NULL  # deprecated; cleared

  push_message(state, sprintf("Loaded dataset: %s", dataset$name), "success")
  invisible(state)
}

# ============================================================================
# Active-annotation access path
# ----------------------------------------------------------------------------
# THE ONLY STANDARD ACCESS PATH for downstream modules. Modules MUST NOT
# index `state$annotation_sets` directly. If a future schema bump changes
# shape, this helper migrates on read so callers stay stable.
# ============================================================================

#' Return the active annotation set object (an `annotation_result_v1`) or
#' NULL if there is no active set.
get_active_annotation <- function(state) {
  id <- state$active_annotation_id
  if (is.null(id)) return(NULL)
  sets <- state$annotation_sets %||% list()
  if (!id %in% names(sets)) return(NULL)
  sets[[id]]
}

# ============================================================================
# Annotation provenance stamp
# ----------------------------------------------------------------------------
# Used by EVERY analysis module to record which annotation set was active
# when the result was produced. Cheap to add now; expensive to retrofit.
#
# A NULL active set produces a stamp where every field is NA -- the slot
# exists, so consumers can rely on its presence.
# ============================================================================

#' Build an annotation-provenance stamp from current state.
#'
#' @return a list with fields:
#'   annotation_set_id_used    character(1) | NA
#'   annotation_set_hash_used  character(1) | NA  (content hash of cell labels)
#'   annotation_set_name       character(1) | NA  (display name at stamp time)
#'   annotation_engine_id      character(1) | NA
#'   annotation_set_is_demo    logical(1)
#'   stamped_at                POSIXct(1)
#'
#' This shape is part of the analysis-result contract; do not change it
#' without a coordinated migration of consumers.
make_annotation_stamp <- function(state) {
  set <- get_active_annotation(state)
  if (is.null(set)) {
    return(list(
      annotation_set_id_used   = NA_character_,
      annotation_set_hash_used = NA_character_,
      annotation_set_name      = NA_character_,
      annotation_engine_id     = NA_character_,
      annotation_set_is_demo   = FALSE,
      stamped_at               = Sys.time()
    ))
  }
  list(
    annotation_set_id_used   = set$set_id,
    annotation_set_hash_used = annotation_set_hash(set),
    annotation_set_name      = set$name %||% set$set_id,
    annotation_engine_id     = set$engine_id %||% NA_character_,
    annotation_set_is_demo   = isTRUE(set$is_demo),
    stamped_at               = Sys.time()
  )
}

#' Detect whether an analysis result's annotation stamp is stale relative to
#' the currently active set.
#'
#' @return TRUE iff result was computed under a different set OR the active
#'   set's content hash has changed since the result was produced.
is_result_stale <- function(result, state) {
  if (is.null(result) || is.null(result$annotation_stamp)) return(FALSE)
  stamp <- result$annotation_stamp
  if (is.na(stamp$annotation_set_id_used)) return(FALSE)
  active <- get_active_annotation(state)
  if (is.null(active)) return(TRUE)
  if (!identical(stamp$annotation_set_id_used, active$set_id)) return(TRUE)
  !identical(stamp$annotation_set_hash_used, annotation_set_hash(active))
}

# ---- Null-coalescing helper (used across the app) -------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
