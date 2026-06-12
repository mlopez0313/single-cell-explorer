# ============================================================================
# Regulon analysis schema
# ----------------------------------------------------------------------------
# Mirrors the annotation_result_v1 / trajectory result patterns. Three
# layered types:
#
#   * `regulon_spec()`     - one regulon: TF + target genes (+ optional
#                            weights / direction).
#   * `regulon_set()`      - a collection of regulons (e.g. "mock_pbmc",
#                            "dorothea_human_AB", ...). Source of truth
#                            consumed by scoring engines. Carries
#                            metadata used for provenance stamping.
#   * `regulon_result_v1()`- per-cell scores for every regulon. Stored
#                            under `state$analysis_results$regulons$results`.
#                            Auc matrix, regulon ids, schema version,
#                            engine id, regulon set id.
#
# The Regulon module never modifies `dataset$cell_data` directly. Like
# annotations + trajectory, regulon scores live in `state` until the
# user opts into baking them in. (apply_regulon_scores_to_dataset()
# will follow the same pattern; deferred to a follow-up.)
# ============================================================================

REGULON_RESULT_SCHEMA_VERSION <- "regulon_result_v1"
REGULON_SET_SCHEMA_VERSION    <- "regulon_set_v1"

#' Specification for one regulon.
#'
#' A regulon is a single transcription factor (TF) and the genes
#' regulated by it. Used by every regulon-scoring engine.
#'
#' @param tf       character(1) snake-case-or-symbol TF name. Does NOT
#'                 need to be present in the dataset (AUCell scoring
#'                 only ranks the targets).
#' @param targets  character() target genes. Empty is allowed (the
#'                 engine will skip / warn).
#' @param weights  numeric() same length as `targets`, default 1.0.
#'                 Engines may treat these as positive importance
#'                 weights; the AUCell engine ignores them and uses
#'                 ranks.
#' @param type     "activating" | "repressing" | "unknown".
regulon_spec <- function(tf, targets,
                         weights = rep(1.0, length(targets)),
                         type    = c("activating", "repressing", "unknown")) {
  type <- match.arg(type)
  stopifnot(is.character(tf), length(tf) == 1L, nzchar(tf))
  stopifnot(is.character(targets))
  stopifnot(is.numeric(weights), length(weights) == length(targets))
  list(tf = tf, targets = targets, weights = as.numeric(weights), type = type)
}

#' A collection of regulons + metadata.
#'
#' Source of truth for any scoring engine. `regulons` is a list of
#' `regulon_spec()` lists; their `tf` fields become the regulon ids.
#'
#' @param id       character(1) snake_case identifier.
#' @param name     character(1) display name.
#' @param species  "human" | "mouse" | "other".
#' @param source   character(1) free-text source label (e.g. "mock",
#'                 "dorothea_AB", "scenic_<run>").
#' @param version  character(1) versioning string for the regulon
#'                 collection itself (independent of the engine
#'                 version that scores it).
#' @param regulons list of `regulon_spec()` entries. TFs must be unique.
regulon_set <- function(id, name, species = "human", source = "user",
                        version = "1.0.0", regulons = list()) {
  stopifnot(is.character(id), length(id) == 1L, nzchar(id))
  stopifnot(is.list(regulons))
  # All entries must be regulon_spec()-shaped.
  for (r in regulons) {
    if (is.null(r$tf) || is.null(r$targets))
      stop("regulon_set(): each regulon needs `tf` and `targets`.",
           call. = FALSE)
  }
  tfs <- vapply(regulons, `[[`, character(1), "tf")
  if (anyDuplicated(tfs))
    stop("regulon_set(): duplicate TF in regulon list: ",
         paste(unique(tfs[duplicated(tfs)]), collapse = ", "),
         call. = FALSE)
  structure(
    list(
      id             = id,
      name           = name,
      species        = species,
      source         = source,
      version        = version,
      regulons       = regulons,
      schema_version = REGULON_SET_SCHEMA_VERSION
    ),
    class = c("regulon_set", "list")
  )
}

#' TRUE iff `x` looks like a regulon_set.
is_regulon_set <- function(x) {
  inherits(x, "regulon_set") &&
    identical(x$schema_version %||% "", REGULON_SET_SCHEMA_VERSION)
}

#' Convert a `regulon_set` into a simple TF -> targets named list.
#' Convenience for engines that don't care about weights / direction.
regulon_set_as_target_list <- function(set) {
  stopifnot(is_regulon_set(set))
  out <- lapply(set$regulons, `[[`, "targets")
  names(out) <- vapply(set$regulons, `[[`, character(1), "tf")
  out
}

#' Build a versioned per-cell regulon result.
#'
#' All consumers (UI, downstream modules, exports, future
#' apply-to-metadata helper) read through this schema.
#'
#' @param cell_ids        character(n_cells)  aligned with `auc_matrix`
#'                                            rows AND `dataset$cell_data$cell`.
#' @param regulon_ids     character(n_regulons) aligned with columns.
#' @param auc_matrix      numeric matrix [n_cells x n_regulons] of AUC
#'                                            scores in [0, 1]. May
#'                                            contain NA for cells with
#'                                            no expressed regulon
#'                                            targets.
#' @param regulon_set_id  the source `regulon_set$id` used.
#' @param engine_id       the engine id from `REGULON_ENGINES()`.
#' @param engine_version  engine version string.
#' @param warnings        character() free-form non-fatal warnings.
regulon_result_v1 <- function(cell_ids, regulon_ids, auc_matrix,
                              regulon_set_id, engine_id,
                              engine_version = "0.1.0",
                              warnings       = character()) {
  stopifnot(is.character(cell_ids))
  stopifnot(is.character(regulon_ids))
  stopifnot(is.numeric(auc_matrix), is.matrix(auc_matrix))
  if (nrow(auc_matrix) != length(cell_ids))
    stop("regulon_result_v1: nrow(auc_matrix) != length(cell_ids).",
         call. = FALSE)
  if (ncol(auc_matrix) != length(regulon_ids))
    stop("regulon_result_v1: ncol(auc_matrix) != length(regulon_ids).",
         call. = FALSE)
  rownames(auc_matrix) <- cell_ids
  colnames(auc_matrix) <- regulon_ids
  structure(
    list(
      cell_ids        = cell_ids,
      regulon_ids     = regulon_ids,
      auc_matrix      = auc_matrix,
      regulon_set_id  = regulon_set_id,
      engine_id       = engine_id,
      engine_version  = engine_version,
      schema_version  = REGULON_RESULT_SCHEMA_VERSION,
      created_at      = Sys.time(),
      warnings        = as.character(warnings)
    ),
    class = c("regulon_result_v1", "list")
  )
}

#' TRUE iff `x` looks like a regulon_result_v1.
is_regulon_result_v1 <- function(x) {
  inherits(x, "regulon_result_v1") &&
    identical(x$schema_version %||% "", REGULON_RESULT_SCHEMA_VERSION) &&
    is.matrix(x$auc_matrix)
}

#' Mean AUC per regulon per cluster (or any categorical grouping vec).
#'
#' Used by the regulon heatmap.
#'
#' @param result        regulon_result_v1
#' @param cluster_vec   character(n_cells) aligned with `result$cell_ids`
#' @return matrix [n_clusters x n_regulons] of mean AUC values.
regulon_mean_by_group <- function(result, cluster_vec) {
  stopifnot(is_regulon_result_v1(result))
  if (length(cluster_vec) != length(result$cell_ids))
    stop("regulon_mean_by_group: cluster_vec length mismatch.",
         call. = FALSE)
  cluster_vec <- as.character(cluster_vec)
  cl_ids <- sort(unique(cluster_vec))
  out <- matrix(NA_real_, nrow = length(cl_ids),
                ncol = ncol(result$auc_matrix),
                dimnames = list(cl_ids, colnames(result$auc_matrix)))
  for (cl in cl_ids) {
    rows <- which(cluster_vec == cl)
    if (length(rows) == 0L) next
    out[cl, ] <- colMeans(result$auc_matrix[rows, , drop = FALSE],
                          na.rm = TRUE)
  }
  out
}
