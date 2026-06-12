# ============================================================================
# Trajectory / pseudotime method registry
# ----------------------------------------------------------------------------
# Mirrors `ANNOTATION_ENGINES()`, `DE_BACKENDS`, and `PATHWAY_SOURCES()`.
# Each method is a `trajectory_method_spec()` with a pure `run_fn` whose
# return value is the `results` payload stored under
# `state$analysis_results$trajectory$results`. The Trajectory module
# never dispatches on method id directly -- it only consults the
# registry.
#
# Built-in methods
# ----------------
#   mock      Deterministic Euclidean distance from a chosen root group's
#             centroid in a 2D reduction. Always available. Clearly demo
#             and labelled as such in the UI banner.
#
#   metadata  Use an existing numeric metadata column directly, scaled
#             to [0, 1]. Always available.
#
#   slingshot `slingshot::slingshot(reducedDims, clusterLabels, start.clus)`
#             on the chosen embedding (PCA preferred; UMAP / DiffMap OK).
#             Real lineage inference. Optional Bioconductor dependency.
#
#   monocle3  `monocle3::learn_graph()` + `order_cells()`. Real lineage
#             inference. Optional Bioconductor dependency. Builds its
#             own CDS from the active dataset, so the mapping back to
#             per-cell pseudotime keeps the schema stable.
#
# All real backends gate on `require_optional()` so the app stays
# runnable when the heavy Bioc packages are not installed.
#
# Future routes (each will land as a new spec):
#   * Palantir   -- via reticulate
#   * Diffusion  -- destiny::DPT
#   * scVelo     -- requires a velocity matrix in the dataset schema
# ============================================================================

#' Specification for one trajectory / pseudotime method.
#'
#' @param id              registry id (used as the value of the module's
#'                        "source" selectInput).
#' @param label           display label.
#' @param requires        character() optional package names. Used by
#'                        `available_trajectory_methods()` to filter the
#'                        UI choices, and by `run_fn`s for actual gating.
#' @param requires_root   logical(1): does the method need a root group?
#'                        Drives whether the module reveals `root_field`
#'                        and `root_group` controls.
#' @param run_fn          function(dataset, params) -> list with
#'                        components used by the Trajectory module:
#'                        \describe{
#'                          \item{pseudotime}{numeric(n_cells) in [0, 1].}
#'                          \item{cell}{character(n_cells) barcode ids.}
#'                          \item{source}{the method id (echo).}
#'                          \item{reduction_used}{character(1) or NA.}
#'                          \item{root_field}{character(1) or NA.}
#'                          \item{root_group}{character(1) or NA.}
#'                          \item{metadata_field}{character(1) or NA.}
#'                          \item{n_lineages}{integer(1) (>=1).}
#'                          \item{method_details}{a free-form named list
#'                            (lineage-specific psts, graph node ids,
#'                            etc.) for downstream introspection.}
#'                        }
#'                        Throw an error on invalid input.
#' @param description     character(1) shown in tooltips / docs.
trajectory_method_spec <- function(id, label,
                                   requires = character(),
                                   requires_root = TRUE,
                                   run_fn,
                                   description = "") {
  stopifnot(is.character(id), length(id) == 1L, nzchar(id))
  stopifnot(is.function(run_fn))
  list(
    id            = id,
    label         = label,
    requires      = requires,
    requires_root = isTRUE(requires_root),
    run_fn        = run_fn,
    description   = description
  )
}

#' The full registry. Function so it re-evaluates per call (cheap), which
#' makes adding a method in an interactive session immediate.
TRAJECTORY_METHODS <- function() {
  list(
    trajectory_method_spec(
      id            = "mock",
      label         = "Mock (compute from root cluster)",
      requires      = character(),
      requires_root = TRUE,
      run_fn        = .run_mock_trajectory,
      description   = paste("Demo-grade pseudotime: Euclidean distance",
                            "from a root group's centroid in a chosen 2D",
                            "embedding, normalised to [0, 1].")
    ),
    trajectory_method_spec(
      id            = "metadata",
      label         = "Numeric metadata field",
      requires      = character(),
      requires_root = FALSE,
      run_fn        = .run_metadata_trajectory,
      description   = paste("Rescale an existing numeric metadata column",
                            "to [0, 1] (e.g. precomputed pseudotime).")
    ),
    trajectory_method_spec(
      id            = "slingshot",
      label         = "Slingshot",
      requires      = "slingshot",
      requires_root = TRUE,
      run_fn        = .run_slingshot_trajectory,
      description   = paste("Bioconductor `slingshot`: minimum spanning",
                            "tree over cluster centroids + per-cell",
                            "principal-curve fits. Real lineage",
                            "inference; reduction can be PCA / UMAP /",
                            "DiffMap.")
    ),
    trajectory_method_spec(
      id            = "monocle3",
      label         = "Monocle3",
      requires      = "monocle3",
      requires_root = TRUE,
      run_fn        = .run_monocle3_trajectory,
      description   = paste("Bioconductor `monocle3`: principal-graph",
                            "learning + `order_cells()`. Requires raw",
                            "counts (uses the `counts` layer when",
                            "present).")
    )
  )
}

#' Look up a method by id. NULL if not registered.
get_trajectory_method <- function(id) {
  for (m in TRAJECTORY_METHODS()) if (identical(m$id, id)) return(m)
  NULL
}

#' Every registered method id, including currently-unavailable ones.
list_trajectory_methods <- function() {
  vapply(TRAJECTORY_METHODS(), `[[`, character(1), "id")
}

#' Methods whose `requires` are installed. Always includes `mock` and
#' `metadata` (no deps); includes optional methods if their packages are
#' present.
available_trajectory_methods <- function() {
  ids <- character()
  for (m in TRAJECTORY_METHODS()) {
    ok <- length(m$requires) == 0L || has_optional(m$requires)
    if (ok) ids <- c(ids, m$id)
  }
  ids
}

#' label -> id named character vector used by the module's selectInput.
#'
#' Unavailable methods are kept in the choices but annotated `(not
#' installed)`. Selecting them yields the same clean `require_optional`
#' error users see from `compute_pseudotime()`.
trajectory_method_choices <- function() {
  out <- character()
  for (m in TRAJECTORY_METHODS()) {
    ok <- length(m$requires) == 0L || has_optional(m$requires)
    label <- if (ok) m$label else sprintf("%s (not installed)", m$label)
    out[label] <- m$id
  }
  out
}
