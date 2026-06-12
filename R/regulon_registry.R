# ============================================================================
# Regulon scoring engine registry
# ----------------------------------------------------------------------------
# Mirrors ANNOTATION_ENGINES / DE_BACKENDS / TRAJECTORY_METHODS. Every
# regulon-scoring backend is one entry; the Regulon module builds its
# engine picker from `list_regulon_engines()`. No UI specialisation per
# engine -- parameters are surfaced from the spec.
#
# Engine contract:
#
#   spec$run_fn(dataset, regulon_set, params)
#     - returns a list with:
#         cell        character(n_cells)
#         regulon_ids character(n_regulons)
#         auc_matrix  numeric matrix [n_cells x n_regulons]
#         warnings    character() optional
#     - must NOT mutate state / dataset.
#
# `run_regulon_engine()` (below) wraps the engine output into a
# `regulon_result_v1` object. Same `annotation_stamp` + provenance
# pattern as DE / pathway / trajectory.
#
# Engines registered:
#   aucell_r   -- pure-R AUCell. No deps. Default in CI / on the mock.
#   aucell     -- AUCell::AUCell_calcAUC(). Heavy Bioc dep.
#   scenic     -- (placeholder) full SCENIC pipeline. Not registered
#                 here yet -- needs the pySCENIC reticulate bridge.
# ============================================================================

#' Specification for one regulon scoring engine.
#'
#' @param id              snake_case identifier
#' @param name            display name for the UI
#' @param requires        character() optional packages required to
#'                        run; consulted by `regulon_engine_available()`
#'                        and `list_regulon_engines()`.
#' @param parameters      named list of parameter specs (same shape as
#'                        ANNOTATION_ENGINES parameters).
#' @param run_fn          function(dataset, regulon_set, params) -> list
#' @param version         character(1) engine implementation version.
#'                        Stamped onto results as `engine_version`.
#'                        DIFFERENT from `result_schema`: this tracks
#'                        the implementation; `result_schema` tracks
#'                        the data shape. Required (no fallback) so
#'                        provenance is always explicit.
#' @param result_schema   schema id that `run_fn`'s output adheres to.
#' @param enabled         logical(1) FALSE hides from the UI
#' @param description     character(1)
regulon_engine_spec <- function(id, name,
                                requires    = character(),
                                parameters  = list(),
                                run_fn,
                                version       = NULL,
                                result_schema = REGULON_RESULT_SCHEMA_VERSION,
                                enabled     = TRUE,
                                description = "") {
  stopifnot(is.character(id), length(id) == 1L, nzchar(id))
  stopifnot(is.function(run_fn))
  if (is.null(version) || !nzchar(version))
    stop(sprintf(
      "regulon_engine_spec('%s'): `version` must be a non-empty string. ",
      id),
      "Engine version is provenance, not optional. ",
      "Use e.g. '<engine_id>_v1.0.0'.",
      call. = FALSE)
  list(
    id            = id,
    name          = name,
    requires      = as.character(requires),
    parameters    = parameters,
    run_fn        = run_fn,
    version       = as.character(version),
    result_schema = result_schema,
    enabled       = isTRUE(enabled),
    description   = description
  )
}

#' The full registry. Function so it picks up changes during dev / testing.
REGULON_ENGINES <- function() {
  list(
    regulon_engine_spec(
      id          = "aucell_r",
      name        = "AUCell (pure R)",
      requires    = character(),
      parameters  = list(
        top_n_fraction = list(type = "numeric", required = FALSE,
                              default = 0.05,
                              description = paste(
                                "Top-N rank threshold as a fraction of",
                                "the gene catalogue. Default 0.05 (top",
                                "5%); higher = more permissive, scores",
                                "saturate slower."))
      ),
      run_fn      = .run_aucell_pure_r_regulons,
      version     = "aucell_r_v1.0.0",
      description = paste(
        "Pure-R implementation of AUCell. No optional dependencies.",
        "Recommended for the mock dataset and CI. Numerically",
        "consistent with the AUCell Bioconductor package on small",
        "datasets but slower on > 50k cells.")
    ),
    regulon_engine_spec(
      id          = "aucell",
      name        = "AUCell (Bioconductor)",
      requires    = c("AUCell"),
      parameters  = list(
        top_n_fraction = list(type = "numeric", required = FALSE,
                              default = 0.05,
                              description = paste(
                                "Top-N rank threshold as a fraction of",
                                "the gene catalogue. Forwarded to",
                                "`AUCell::AUCell_calcAUC`'s `aucMaxRank`."))
      ),
      run_fn      = .run_aucell_bioc_regulons,
      version     = "aucell_bioc_v1.0.0",
      description = paste(
        "Wraps `AUCell::AUCell_calcAUC()`. Faster on large sparse",
        "datasets; sparse-matrix aware.")
    )
  )
}

#' Look up a single engine by id. NULL if not registered.
get_regulon_engine <- function(id) {
  for (e in REGULON_ENGINES()) if (identical(e$id, id)) return(e)
  NULL
}

#' TRUE iff every package in `engine$requires` is installed.
regulon_engine_available <- function(engine) {
  if (is.null(engine)) return(FALSE)
  has_optional(engine$requires)
}

#' Enabled engines as a named character vector: display_name -> id.
#' Unavailable engines are kept but labelled, mirroring DE backends &
#' trajectory methods.
list_regulon_engines <- function(enabled_only = TRUE) {
  engines <- REGULON_ENGINES()
  if (enabled_only) engines <- Filter(function(e) isTRUE(e$enabled), engines)
  labels <- vapply(engines, function(e) {
    avail <- regulon_engine_available(e)
    if (avail) e$name
    else sprintf("%s  (not installed)", e$name)
  }, character(1))
  setNames(vapply(engines, `[[`, character(1), "id"), labels)
}

# ============================================================================
# Dispatcher
# ----------------------------------------------------------------------------
# `run_regulon_engine()` is the canonical entrypoint. The module calls
# it, the engine returns the contract list, the dispatcher wraps it in
# a regulon_result_v1. Mirrors `run_annotation_engine()` in
# R/annotation.R.
# ============================================================================

#' Run a regulon scoring engine and return a regulon_result_v1.
#'
#' @param engine_id    one of `names(list_regulon_engines(enabled_only = FALSE))`
#' @param dataset      the active dataset
#' @param regulon_set  a `regulon_set()`; supplies targets per TF
#' @param params       engine-specific parameters list
#' @return regulon_result_v1
run_regulon_engine <- function(engine_id, dataset, regulon_set,
                               params = list()) {
  engine <- get_regulon_engine(engine_id)
  if (is.null(engine))
    stop(sprintf("Unknown regulon engine '%s'. Have: %s",
                 engine_id,
                 paste(unname(list_regulon_engines(enabled_only = FALSE)),
                       collapse = ", ")), call. = FALSE)
  if (!is_regulon_set(regulon_set))
    stop("run_regulon_engine: `regulon_set` is not a regulon_set().",
         call. = FALSE)
  if (is.null(dataset)) stop("No dataset provided.", call. = FALSE)

  raw <- engine$run_fn(dataset, regulon_set, params)
  if (is.null(raw$auc_matrix) || !is.matrix(raw$auc_matrix))
    stop(sprintf(
      "Regulon engine '%s' did not return an `auc_matrix`.", engine_id),
      call. = FALSE)
  # `engine$version` is guaranteed by regulon_engine_spec()'s validator,
  # so the previous %||% fallback to REGULON_ENGINE_DEFAULT_VERSION is no
  # longer needed (and was hiding the bug where engines forgot to declare
  # a version).
  regulon_result_v1(
    cell_ids       = raw$cell %||% dataset$cell_data$cell,
    regulon_ids    = raw$regulon_ids %||% colnames(raw$auc_matrix),
    auc_matrix     = raw$auc_matrix,
    regulon_set_id = regulon_set$id,
    engine_id      = engine_id,
    engine_version = engine$version,
    warnings       = raw$warnings %||% character()
  )
}
