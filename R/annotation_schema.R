# ============================================================================
# Annotation result schema (versioned)
# ----------------------------------------------------------------------------
# A versioned, engine-agnostic container for annotation outputs. Every
# annotation engine in `ANNOTATION_ENGINES()` produces objects of this shape
# regardless of whether it is manual, marker-score, SingleR, Azimuth,
# CellTypist, consensus-voting, or a future LLM-assisted engine.
#
# Per-cell labels are the canonical storage even when the engine works at
# the cluster level -- the cluster->cell expansion happens here so the rest
# of the app never has to care.
#
# Schema name: "annotation_v1"
#
# If the schema needs to change, bump to "annotation_v2", keep the v1 reader
# in this file, and add a forward migration. Result objects on disk and in
# saved state remain readable.
# ============================================================================

ANNOTATION_RESULT_SCHEMA_VERSION <- "annotation_v1"

#' Construct an annotation_result_v1 object.
#'
#' Manual annotation fills the full schema with `cell_scores = 1.0` (or 0.0
#' for unannotated cells) and `alt_labels = NULL`. There is no separate
#' "manual annotation" schema; engine differences are surfaced as field
#' contents, not different shapes.
#'
#' @param set_id                   character(1) unique set id (stable handle)
#' @param set_name                 character(1) display name (user-visible)
#' @param engine_id                character(1) which engine produced this
#' @param engine_version           character(1)
#' @param params                   list of engine params used
#' @param cell                     character(n_cells) cell barcodes
#' @param cell_labels              character(n_cells) per-cell labels (NA allowed)
#' @param cell_scores              numeric(n_cells) per-cell confidence
#'                                  (1.0 for manual confirmed, 0.0 for blanks)
#' @param alt_labels               NULL | data.frame(cell, rank, label, score)
#'                                  top-k alternative labels per cell
#' @param cluster_summary          NULL | data.frame describing per-cluster
#'                                  top label / score
#' @param ontology_map             NULL | named character: label -> ontology_id
#' @param reference_source         character(1) e.g. "SingleR/ImmGen", "user"
#' @param marker_registry_version  character(1) version stamp of the registry
#'                                  the engine consulted (NA for manual)
#' @param parent_set_id            character(1) id of the set this was duplicated
#'                                  from, or NA
#' @param cluster_field_used       character(1) metadata field at assignment time
#' @param n_clusters_at_creation   integer(1) cardinality of that field at time
#'                                  of creation (used for stale-set detection)
#' @param is_frozen                logical(1) freeze edits
#' @param is_demo                  logical(1) flag results derived from mock /
#'                                  demo flows
#' @param description              character(1)
#' @param warnings                 character() engine warnings
#' @param created_at, modified_at  POSIXct
#' @param timestamp                POSIXct of latest run
#' @param duration_ms              integer(1)
#' @param error_message            character(1) | NULL
#' @param edit_history             list of edit records (audit log)
annotation_result_v1 <- function(set_id,
                                 set_name               = set_id,
                                 engine_id              = "manual",
                                 engine_version         = "0.1.0",
                                 params                 = list(),
                                 cell                   = character(),
                                 cell_labels            = character(),
                                 cell_scores            = NULL,
                                 alt_labels             = NULL,
                                 cluster_summary        = NULL,
                                 ontology_map           = NULL,
                                 reference_source       = NA_character_,
                                 marker_registry_version = NA_character_,
                                 parent_set_id          = NA_character_,
                                 cluster_field_used     = NA_character_,
                                 n_clusters_at_creation = NA_integer_,
                                 is_frozen              = FALSE,
                                 is_demo                = FALSE,
                                 description            = "",
                                 warnings               = character(),
                                 created_at             = Sys.time(),
                                 modified_at            = Sys.time(),
                                 timestamp              = Sys.time(),
                                 duration_ms            = NA_integer_,
                                 error_message          = NULL,
                                 edit_history           = list()) {
  stopifnot(is.character(set_id), length(set_id) == 1L, nzchar(set_id))
  stopifnot(length(cell) == length(cell_labels))
  if (is.null(cell_scores)) {
    cell_scores <- ifelse(is.na(cell_labels) | !nzchar(as.character(cell_labels)),
                          0.0, 1.0)
  }
  stopifnot(length(cell_scores) == length(cell))
  list(
    schema_version          = ANNOTATION_RESULT_SCHEMA_VERSION,
    set_id                  = set_id,
    name                    = set_name,
    description             = description,
    engine_id               = engine_id,
    engine_version          = engine_version,
    params                  = params,
    cell                    = as.character(cell),
    cell_labels             = as.character(cell_labels),
    cell_scores             = as.numeric(cell_scores),
    alt_labels              = alt_labels,
    cluster_summary         = cluster_summary,
    ontology_map            = ontology_map,
    reference_source        = reference_source,
    marker_registry_version = marker_registry_version,
    parent_set_id           = parent_set_id,
    cluster_field_used      = cluster_field_used,
    n_clusters_at_creation  = n_clusters_at_creation,
    is_frozen               = isTRUE(is_frozen),
    is_demo                 = isTRUE(is_demo),
    warnings                = as.character(warnings),
    created_at              = created_at,
    modified_at             = modified_at,
    timestamp               = timestamp,
    duration_ms             = duration_ms,
    error_message           = error_message,
    edit_history            = edit_history
  )
}

#' Defensive predicate: does this object look like an annotation_v1 result?
is_annotation_result_v1 <- function(x) {
  is.list(x) &&
    identical(x$schema_version, ANNOTATION_RESULT_SCHEMA_VERSION) &&
    !is.null(x$cell) && !is.null(x$cell_labels) &&
    length(x$cell) == length(x$cell_labels)
}

# ---- Helpers --------------------------------------------------------------

#' Expand a cluster->label mapping to per-cell labels aligned with the dataset.
#'
#' @param dataset        dataset list
#' @param cluster_field  metadata column used as the grouping key
#' @param cluster_labels named character vector: cluster_id -> label
#'
#' @return character vector of length n_cells. Cells whose cluster has no
#'   assignment receive NA.
expand_cluster_to_cells <- function(dataset, cluster_field, cluster_labels) {
  vals <- get_metadata(dataset, cluster_field)
  if (is.null(vals)) stop(sprintf("Cluster field '%s' not in dataset.",
                                  cluster_field %||% ""), call. = FALSE)
  vals <- as.character(vals)
  if (is.null(cluster_labels) || !length(cluster_labels))
    return(rep(NA_character_, length(vals)))
  out <- rep(NA_character_, length(vals))
  for (cl in names(cluster_labels)) {
    out[vals == cl] <- as.character(cluster_labels[[cl]])
  }
  out
}

#' Deterministic content hash for an annotation set, used by downstream
#' stamping to detect stale results.
#'
#' Captures cell -> label pairs (sorted, joined). Independent of metadata
#' fields, engine params, and timestamps so unchanged labels produce the
#' same hash even when re-run.
annotation_set_hash <- function(set) {
  if (is.null(set) || is.null(set$cell) || is.null(set$cell_labels)) {
    return(NA_character_)
  }
  ord  <- order(set$cell)
  s    <- paste(set$cell[ord], set$cell_labels[ord], sep = "=", collapse = ";")
  if (requireNamespace("rlang", quietly = TRUE)) {
    return(as.character(rlang::hash(s)))
  }
  # Fallback (deterministic non-cryptographic). Not as collision-resistant
  # as rlang::hash but adequate for stale-detection.
  bytes <- utf8ToInt(s)
  h     <- sum(bytes * seq_along(bytes)) %% 2147483647L
  sprintf("h%08x", as.integer(h))
}

#' Compact summary string for use in banners and stamps.
annotation_set_label <- function(set) {
  if (is.null(set)) return("(no active annotation)")
  n_labelled <- sum(!is.na(set$cell_labels) & nzchar(as.character(set$cell_labels)))
  flags <- character()
  if (isTRUE(set$is_demo))   flags <- c(flags, "demo")
  if (isTRUE(set$is_frozen)) flags <- c(flags, "frozen")
  flag_str <- if (length(flags)) sprintf(" [%s]", paste(flags, collapse = ", ")) else ""
  sprintf("%s (engine=%s, labelled=%d/%d)%s",
          set$name %||% set$set_id,
          set$engine_id %||% "?",
          n_labelled, length(set$cell_labels),
          flag_str)
}
