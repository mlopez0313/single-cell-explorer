# ============================================================================
# Annotation orchestration
# ----------------------------------------------------------------------------
# Annotation logic lives across four files:
#   R/marker_registry.R       -- the typed marker registry (data, not code)
#   R/annotation_schema.R     -- `annotation_result_v1` constructor + helpers
#   R/annotation_registry.R   -- the `ANNOTATION_ENGINES()` registry + engines
#   R/annotation.R (this)     -- multi-set management, dispatcher, apply
#
# The Annotation module and downstream modules read through
# `get_active_annotation(state)` (defined in R/state.R) and never reach
# directly into `state$annotation_sets`.
#
# A note on schema versioning:
#   Annotation result objects carry `schema_version = "annotation_v1"`.
#   When the schema bumps to v2, this file gains a v1 -> v2 migration
#   helper and the registry's `result_schema` field is bumped per engine
#   independently. Saved sessions remain readable.
# ============================================================================

# ---- Set ids --------------------------------------------------------------

#' Generate a stable, sortable, unique set id.
#'
#' Format: `<prefix>_<YYYYmmdd_HHMMSS>_<4 random letters>`. Used for
#' provenance-named metadata columns (`annotation__<set_id>__<yyyy_mm_dd>`).
new_annotation_set_id <- function(prefix = "set") {
  paste(prefix,
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        paste(sample(letters, 4L, replace = TRUE), collapse = ""),
        sep = "_")
}

# ---- Multi-set management --------------------------------------------------
# All functions are pure with respect to the dataset; they mutate state
# only by writing to its reactive fields. Callers are responsible for using
# `shiny::isolate()` when invoking these outside a reactive context.

#' Insert (or overwrite by id) an annotation set in shared state.
#' Activates the new set if no set is currently active.
add_annotation_set <- function(state, set) {
  stopifnot(is_annotation_result_v1(set))
  sets <- state$annotation_sets %||% list()
  sets[[set$set_id]] <- set
  state$annotation_sets <- sets
  if (is.null(state$active_annotation_id)) {
    state$active_annotation_id <- set$set_id
  }
  invisible(set$set_id)
}

#' Mark `set_id` as the active annotation set. Pass NULL to deactivate.
set_active_annotation <- function(state, set_id) {
  if (is.null(set_id)) {
    state$active_annotation_id <- NULL
    return(invisible(NULL))
  }
  if (!set_id %in% names(state$annotation_sets %||% list())) {
    stop(sprintf("No annotation set with id '%s'.", set_id), call. = FALSE)
  }
  state$active_annotation_id <- set_id
  invisible(set_id)
}

#' Remove a set; refuses to delete frozen sets.
remove_annotation_set <- function(state, set_id) {
  sets <- state$annotation_sets %||% list()
  if (!set_id %in% names(sets)) return(invisible())
  if (isTRUE(sets[[set_id]]$is_frozen)) {
    stop(sprintf("Set '%s' is frozen; unfreeze before deleting.", set_id),
         call. = FALSE)
  }
  sets[[set_id]] <- NULL
  state$annotation_sets <- sets
  if (identical(state$active_annotation_id, set_id)) {
    state$active_annotation_id <- if (length(sets)) names(sets)[1] else NULL
  }
  invisible()
}

#' Duplicate an existing set. Returns the new set's id.
duplicate_annotation_set <- function(state, set_id, new_name = NULL) {
  sets <- state$annotation_sets %||% list()
  src  <- sets[[set_id]]
  if (is.null(src)) stop(sprintf("No annotation set '%s'.", set_id), call. = FALSE)
  new_id <- new_annotation_set_id("dup")
  dup    <- src
  dup$set_id        <- new_id
  dup$name          <- new_name %||% paste0(src$name %||% src$set_id, " (copy)")
  dup$parent_set_id <- src$set_id
  dup$created_at    <- Sys.time()
  dup$modified_at   <- Sys.time()
  dup$is_frozen     <- FALSE
  add_annotation_set(state, dup)
  new_id
}

#' Rename a set in place. Frozen sets cannot be renamed.
rename_annotation_set <- function(state, set_id, new_name) {
  sets <- state$annotation_sets %||% list()
  if (!set_id %in% names(sets)) {
    stop(sprintf("No annotation set '%s'.", set_id), call. = FALSE)
  }
  if (isTRUE(sets[[set_id]]$is_frozen)) {
    stop(sprintf("Set '%s' is frozen; rename refused.", set_id), call. = FALSE)
  }
  sets[[set_id]]$name        <- as.character(new_name)
  sets[[set_id]]$modified_at <- Sys.time()
  state$annotation_sets <- sets
  invisible()
}

#' Toggle the frozen flag on a set.
freeze_annotation_set <- function(state, set_id, frozen = TRUE) {
  sets <- state$annotation_sets %||% list()
  if (!set_id %in% names(sets)) {
    stop(sprintf("No annotation set '%s'.", set_id), call. = FALSE)
  }
  sets[[set_id]]$is_frozen   <- isTRUE(frozen)
  sets[[set_id]]$modified_at <- Sys.time()
  state$annotation_sets <- sets
  invisible()
}

# ---- Engine dispatch ------------------------------------------------------

#' Run an annotation engine and wrap the result as `annotation_result_v1`.
#'
#' This is the single dispatcher; the Annotation module calls only this
#' function and `add_annotation_set()`. Engine selection lives in the
#' `ANNOTATION_ENGINES()` registry.
#'
#' @param engine_id      registered engine id
#' @param dataset        active dataset
#' @param state          shared app state (read-only inside the engine)
#' @param params         named list of params for the engine
#' @param set_id         the id under which to store the result (existing
#'                       set is updated in place if the id already exists)
#' @param set_name       display name
#' @param description    free-text description
#' @param parent_set_id  set this is derived from (NA for a fresh set)
#' @param is_demo        forwarded onto the result
run_annotation_engine <- function(engine_id, dataset, state, params,
                                  set_id, set_name = set_id,
                                  description    = "",
                                  parent_set_id  = NA_character_,
                                  is_demo        = FALSE) {
  engine <- get_annotation_engine(engine_id)
  if (is.null(engine)) {
    stop(sprintf("Unknown annotation engine '%s'.", engine_id), call. = FALSE)
  }

  t0  <- proc.time()[["elapsed"]]
  raw <- engine$run_fn(dataset, params, state)
  dur <- as.integer(round((proc.time()[["elapsed"]] - t0) * 1000))

  # Build ontology map from the labels using the marker registry, if any.
  ontology_map <- NULL
  if (!is.null(state$marker_registry)) {
    ontology_map <- tryCatch(
      build_ontology_map(state$marker_registry, raw$cell_labels),
      error = function(e) NULL
    )
  }

  # NB: `engine_version` is the engine *implementation* version, not the
  # result-schema version. Earlier code mistakenly stamped
  # `engine$result_schema` here; that conflated two distinct provenance
  # axes. The annotation_engine_spec() constructor now requires an
  # explicit `version`, so this lookup is guaranteed to find one.
  annotation_result_v1(
    set_id                  = set_id,
    set_name                = set_name,
    engine_id               = engine_id,
    engine_version          = engine$version,
    params                  = params,
    cell                    = raw$cell,
    cell_labels             = raw$cell_labels,
    cell_scores             = raw$cell_scores,
    alt_labels              = raw$alt_labels,
    cluster_summary         = raw$cluster_summary,
    ontology_map            = ontology_map,
    reference_source        = raw$reference_source %||% NA_character_,
    marker_registry_version = if ("marker_registry" %in% engine$requires &&
                                  !is.null(state$marker_registry))
                                state$marker_registry$version
                              else NA_character_,
    parent_set_id           = parent_set_id,
    cluster_field_used      = raw$cluster_field_used     %||% NA_character_,
    n_clusters_at_creation  = raw$n_clusters_at_creation %||% NA_integer_,
    is_demo                 = is_demo,
    description             = description,
    warnings                = raw$warnings %||% character(),
    duration_ms             = dur
  )
}

# ---- Apply to dataset (name-driven column) --------------------------------

#' Sanitize a free-form annotation set name into a valid metadata column.
#'
#' Lower-cases, maps any non `[a-z0-9_]` character to an underscore,
#' collapses runs of underscores, and strips leading/trailing
#' underscores. Falls back to `"annotation"` when the input is empty
#' or sanitizes away. Prepends `x_` when the result starts with a
#' digit so the column is a valid R name.
#'
#' Pure; safe to unit-test.
#'
#' @param name character(1)
#' @return character(1)
.annotation_col_basename <- function(name) {
  s <- tolower(as.character(name %||% ""))
  s <- gsub("[^a-z0-9_]+", "_", s)
  s <- gsub("_+", "_", s)
  s <- gsub("^_+|_+$", "", s)
  if (!nzchar(s)) s <- "annotation"
  if (grepl("^[0-9]", s)) s <- paste0("x_", s)
  s
}

#' Find an existing dataset column produced by a given annotation set.
#'
#' Walks `dataset$cell_data` and returns the first column whose
#' `annotation_set_id` attribute matches `set_id`, or NULL.
#' Attribute-based (not name-based) so renaming the set doesn't lose
#' track of the previously-applied column.
.find_annotation_col_for_set <- function(dataset, set_id) {
  if (is.null(dataset$cell_data) || is.null(set_id)) return(NULL)
  for (col in names(dataset$cell_data)) {
    a <- attr(dataset$cell_data[[col]], "annotation_set_id", exact = TRUE)
    if (!is.null(a) && identical(a, set_id)) return(col)
  }
  NULL
}

#' Persist an annotation set onto the dataset's metadata table.
#'
#' The applied column name is derived from `set$name` (sanitized to a
#' valid R identifier). Re-applying the same set after a rename
#' renames the existing column accordingly -- the column is keyed by
#' the set's `set_id`, stored as the column's `annotation_set_id`
#' attribute, NOT by its surface name. This means:
#'
#'   * Renaming a set and clicking "Apply" again moves the
#'     previously-applied column to the new name.
#'   * Re-running the same engine on the same set (different label
#'     values, same set_id) updates the existing column in place.
#'   * Two unrelated sets that happen to share a base name get
#'     numeric suffixes (`name`, `name_2`, `name_3`, ...) so neither
#'     clobbers the other.
#'
#' The annotation set itself remains the primary source of truth --
#' this function exists for export workflows, downstream tools that
#' only know how to read metadata columns, and for users who want to
#' color by the label inside the Explorer.
#'
#' Provenance attributes are attached to the new column via `attr()`
#' so downstream readers can trace which set produced the values
#' regardless of the column's surface name.
#'
#' Refuses to write the literal name `"cell_type"` (a generic name
#' that downstream tools may treat specially).
#'
#' @param dataset  the dataset list
#' @param set      an `annotation_result_v1` object
#'
#' @return a new dataset list (does not mutate in place).
apply_annotations_to_dataset <- function(dataset, set) {
  if (is.null(dataset)) stop("No dataset provided.", call. = FALSE)
  if (!is_annotation_result_v1(set))
    stop("`set` must be an annotation_result_v1 object.", call. = FALSE)

  ds_cells  <- dataset$cell_data$cell
  if (is.null(ds_cells))
    stop("Dataset has no `cell_data$cell` column.", call. = FALSE)

  pos <- match(ds_cells, set$cell)
  if (any(is.na(pos))) {
    n_missing <- sum(is.na(pos))
    stop(sprintf(
      "Annotation set covers %d/%d cells; %d dataset cells have no label.",
      sum(!is.na(pos)), length(ds_cells), n_missing), call. = FALSE)
  }
  values <- set$cell_labels[pos]

  # Desired column name = sanitized set name. Fall back to set_id when
  # the name is missing / blank (shouldn't happen via the UI, which
  # always defaults to a non-empty name -- but keep it robust).
  base_name <- .annotation_col_basename(set$name %||% set$set_id)
  if (identical(base_name, "cell_type")) {
    stop("Refusing to write a generic 'cell_type' column. Rename the ",
         "annotation set to something more specific.", call. = FALSE)
  }

  # If a column for THIS set already exists (same `annotation_set_id`
  # attribute), remove it first. The new column may end up with a
  # different name; equivalent to renaming the column in place.
  prev_col <- .find_annotation_col_for_set(dataset, set$set_id)
  if (!is.null(prev_col)) {
    dataset$cell_data[[prev_col]] <- NULL
    dataset$metadata_fields <- setdiff(dataset$metadata_fields, prev_col)
  }

  # Disambiguate against unrelated columns (different set_id but same
  # desired base name).
  col_name <- base_name
  i <- 2L
  while (col_name %in% names(dataset$cell_data)) {
    col_name <- sprintf("%s_%d", base_name, i)
    i <- i + 1L
  }

  dataset$cell_data[[col_name]] <- values
  attr(dataset$cell_data[[col_name]], "annotation_set_id")         <- set$set_id
  attr(dataset$cell_data[[col_name]], "annotation_set_name")       <- set$name
  attr(dataset$cell_data[[col_name]], "annotation_engine_id")      <- set$engine_id
  attr(dataset$cell_data[[col_name]], "annotation_engine_version") <- set$engine_version
  attr(dataset$cell_data[[col_name]], "marker_registry_version")   <- set$marker_registry_version
  attr(dataset$cell_data[[col_name]], "schema_version")            <- set$schema_version
  attr(dataset$cell_data[[col_name]], "applied_at")                <- Sys.time()

  dataset$metadata_fields <- unique(c(dataset$metadata_fields, col_name))
  dataset
}

#' Names of columns in `dataset$cell_data` produced by an annotation set.
#'
#' Detected via the `annotation_set_id` attribute, so renames don't
#' break the lookup. Used by the Annotation module's cluster-field
#' picker to hide annotation columns (recursive annotation is rarely
#' what you want).
annotation_columns <- function(dataset) {
  if (is.null(dataset$cell_data)) return(character())
  out <- character()
  for (col in names(dataset$cell_data)) {
    a <- attr(dataset$cell_data[[col]], "annotation_set_id", exact = TRUE)
    if (!is.null(a)) out <- c(out, col)
  }
  out
}

# ---- CSV export (annotation result -> wide / tidy table) ------------------

#' Write an annotation set to CSV. Cluster summary if present, otherwise
#' a cell-level table.
write_annotation_set_csv <- function(set, file) {
  if (is.null(set)) {
    utils::write.csv(data.frame(), file = file, row.names = FALSE)
    return(invisible(file))
  }
  out <- if (!is.null(set$cluster_summary)) set$cluster_summary
         else data.frame(cell = set$cell, label = set$cell_labels,
                         score = set$cell_scores, stringsAsFactors = FALSE)
  utils::write.csv(out, file = file, row.names = FALSE, na = "")
  invisible(file)
}

# ---- Convenience accessors ------------------------------------------------

#' Active set's per-cell labels aligned with `dataset$cell_data$cell`, or
#' NULL if no active set / no dataset.
get_cell_labels <- function(state) {
  set <- get_active_annotation(state)
  if (is.null(set)) return(NULL)
  ds  <- state$active_dataset
  if (is.null(ds)) return(NULL)
  pos <- match(ds$cell_data$cell, set$cell)
  set$cell_labels[pos]
}

#' Active set's label for a given cluster id in the set's cluster_summary,
#' or NA if no summary / not found.
get_label_for_cluster <- function(state, cluster_id) {
  set <- get_active_annotation(state)
  if (is.null(set) || is.null(set$cluster_summary)) return(NA_character_)
  row <- set$cluster_summary[as.character(set$cluster_summary$cluster) ==
                             as.character(cluster_id), , drop = FALSE]
  if (!nrow(row)) NA_character_ else row$top_label[1]
}

# ---- Legacy helpers retained (only those still used elsewhere) ------------
# `top_markers_per_cluster()` is still used by the marker-score scoring
# preview and by the Annotation module's "suggest" affordance. Keep the
# pure function; everything else moves to the registry-driven path.

get_cluster_ids <- function(dataset, cluster_field = "cluster") {
  meta <- get_metadata(dataset, cluster_field)
  if (is.null(meta)) return(character())
  as.character(sort(unique(meta)))
}

count_cells_per_cluster <- function(dataset, cluster_field = "cluster") {
  meta <- get_metadata(dataset, cluster_field)
  if (is.null(meta)) return(integer())
  tbl <- table(as.character(meta))
  out <- as.integer(tbl); names(out) <- names(tbl); out
}

top_markers_per_cluster <- function(dataset, cluster_field = "cluster",
                                    n_per = 3) {
  df <- compute_markers(dataset, cluster_field, top_n = n_per)
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  by_grp <- split(df$gene, df$group)
  vapply(by_grp, paste, collapse = ", ", FUN.VALUE = character(1))
}
