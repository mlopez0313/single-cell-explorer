# ============================================================================
# Data Smoothing / Imputation
# ----------------------------------------------------------------------------
# Produces *smoothed* expression vectors for visualization-only exploration.
# Raw expression remains the analytic source of truth -- DE, Marker
# Investigation, and Pathway Analysis must continue to read
# `get_gene_expression()`, not the smoothed values produced here.
#
# Output is stored under `state$analysis_results$imputation` with the same
# envelope shape as other analyses:
#   list(status=, results=, params=, error_message=, timestamp=, duration_ms=)
# where `results` is:
#   list(method=<chr>, genes=<chr>, expression=<named list gene -> numeric>,
#        reduction_used=<chr>, k=<int>)
#
# Method registry (all three are mock; UI is honest about it):
#   - "none"        : pass-through (returns raw expression unchanged)
#   - "neighbor"    : one-pass kNN average in the embedding space
#   - "alra_mock"   : kNN average + soft-threshold (sets low values to 0)
#   - "magic_mock"  : two-pass kNN average (diffusion-style)
#
# TODO (future, real implementations -- not in this module):
#   - ALRA            : low-rank SVD recovery (Rafii et al., 2018). Real
#                       integration via the `ALRA::alra()` package.
#   - MAGIC           : Markov-Affinity diffusion. Real integration via the
#                       `Rmagic` or `phateR` packages.
#   - SAVER / scImpute: Bayesian / regression-based imputation.
#   - kNN graph       : compute kNN on PCA (not the 2D embedding) for real
#                       data; cache the graph for reuse across methods.
# ============================================================================

#' Available imputation/smoothing methods (display label -> internal id).
IMPUTATION_METHODS <- c(
  "None / Raw"                = "none",
  "Simple neighbor smoothing" = "neighbor",
  "Mock ALRA"                 = "alra_mock",
  "Mock MAGIC"                = "magic_mock"
)

available_imputation_methods <- function() IMPUTATION_METHODS

#' Per-cell k-nearest-neighbor indices in 2D embedding space.
#'
#' Returns an integer matrix `n_cells x k`. O(n^2) -- fine for mock data
#' (n ~ thousands). Real loaders should call out to `RANN::nn2()` or
#' `BiocNeighbors` on the PCA embedding instead.
compute_knn_indices <- function(emb, k = 15L) {
  if (is.null(emb)) return(NULL)
  k <- max(1L, as.integer(k))
  x <- emb$x; y <- emb$y; n <- length(x)
  k <- min(k, n)
  out <- matrix(NA_integer_, nrow = n, ncol = k)
  for (i in seq_len(n)) {
    d2 <- (x - x[i])^2 + (y - y[i])^2
    out[i, ] <- order(d2)[seq_len(k)]
  }
  out
}

#' Run smoothing/imputation on a dataset.
#'
#' @param dataset   the active dataset (`dataset_schema()`)
#' @param genes     character() target genes; falls back to all available
#'                  genes when empty
#' @param method    one of `IMPUTATION_METHODS`
#' @param k         neighborhood size for kNN-based methods
#' @param reduction reduction to use as the smoothing graph (defaults to
#'                  the dataset's default reduction; future: PCA on real data)
#'
#' @return a list:
#'   list(method=, genes=, expression=<named list>, reduction_used=, k=)
compute_smoothed <- function(dataset, genes = NULL,
                             method = "neighbor", k = 15L,
                             reduction = NULL) {
  if (is.null(dataset)) stop("No dataset provided.", call. = FALSE)
  method <- as.character(method); if (!nzchar(method)) method <- "neighbor"
  if (!method %in% IMPUTATION_METHODS)
    stop("Unknown method: ", method, call. = FALSE)

  all_genes <- available_genes(dataset)
  if (is.null(genes) || length(genes) == 0L) genes <- all_genes
  genes <- intersect(genes, all_genes)
  if (length(genes) == 0L)
    stop("No target genes available in this dataset.", call. = FALSE)

  if (identical(method, "none")) {
    expr <- lapply(genes, function(g) get_gene_expression(dataset, g))
    names(expr) <- genes
    return(list(method = method, genes = genes, expression = expr,
                reduction_used = NA_character_, k = NA_integer_))
  }

  red <- reduction %||% dataset$default_reduction %||%
         (dataset$reductions %||% character())[1]
  emb <- get_embedding(dataset, red)
  if (is.null(emb))
    stop("No usable reduction for smoothing.", call. = FALSE)

  k_effective <- switch(method,
                        "neighbor"   = as.integer(k),
                        "alra_mock"  = as.integer(max(5L, round(k * 2/3))),
                        "magic_mock" = as.integer(max(10L, round(k * 2))))
  nn <- compute_knn_indices(emb, k = k_effective)

  smooth_one <- function(v) {
    out <- rowMeans(matrix(v[as.vector(nn)], nrow = nrow(nn), ncol = ncol(nn)))
    if (identical(method, "alra_mock")) {
      thr <- stats::quantile(out, 0.30, na.rm = TRUE)
      out[out < thr] <- 0
    } else if (identical(method, "magic_mock")) {
      # Second diffusion pass on the same graph.
      out <- rowMeans(matrix(out[as.vector(nn)], nrow = nrow(nn), ncol = ncol(nn)))
    }
    out
  }

  expression <- list()
  for (g in genes) {
    v <- get_gene_expression(dataset, g)
    if (is.null(v)) next
    expression[[g]] <- smooth_one(v)
  }

  list(method = method, genes = names(expression),
       expression = expression, reduction_used = red, k = k_effective)
}

# ----------------------------------------------------------------------------
# Display-mode switch for visualization modules.
#
# The Basic scRNA Explorer's FeaturePlot uses this helper -- it picks raw
# vs smoothed based on `state$display_mode_imputation`. Crucially: DE,
# Marker Investigation, and Pathway Analysis do NOT call this helper. They
# call `get_gene_expression()` directly, which always returns raw values.
# That's intentional: smoothed expression must never silently leak into
# statistical analyses.
# ----------------------------------------------------------------------------

#' Return the expression vector that should be *visualized* for `gene`.
#'
#' Returns smoothed values iff
#'   - `state$display_mode_imputation == "smoothed"`, AND
#'   - smoothed results exist for `gene`.
#' Falls back to raw `get_gene_expression()` in every other case.
get_gene_expression_for_view <- function(state, gene) {
  ds <- state$active_dataset
  if (is.null(ds) || is.null(gene) || !nzchar(gene)) return(NULL)
  mode <- state$display_mode_imputation %||% "raw"
  if (identical(mode, "smoothed")) {
    imp <- state$analysis_results$imputation
    if (!is.null(imp) && identical(imp$status, "completed")) {
      v <- imp$results$expression[[gene]]
      if (!is.null(v)) return(v)
    }
  }
  get_gene_expression(ds, gene)
}

#' TRUE iff `state` has a completed imputation result with at least one gene.
has_smoothed_results <- function(state) {
  imp <- state$analysis_results$imputation
  !is.null(imp) && identical(imp$status, "completed") &&
    !is.null(imp$results$expression) && length(imp$results$expression) > 0L
}
