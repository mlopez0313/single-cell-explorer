# ============================================================================
# Annotation engine: CellTypist
# ----------------------------------------------------------------------------
# Reference-based per-cell annotation via the Teichmann lab's
# CellTypist Python package. We bridge through `reticulate`. Pluggable
# through `ANNOTATION_ENGINES()`; the existing Annotation module UI
# renders parameters from the spec.
#
# Why CellTypist:
#   * 50+ pretrained models covering immune (`Immune_All_Low.pkl`,
#     `Immune_All_High.pkl`), lung, gut, brain, dev, etc.
#   * Cluster-aware "majority voting" mode reduces noise on small
#     clusters.
#   * Reference catalog is downloadable on demand; CellTypist caches.
#
# Layout (mirrors annotation_singler.R / annotation_azimuth.R):
#   * `.run_celltypist_annotation()`   -- the engine `run_fn`. Bridges
#                                          through `reticulate`. Heavy
#                                          deps gated via
#                                          `require_optional`.
#   * `.celltypist_to_engine_output()` -- pure schema converter. Takes
#                                          either the Python
#                                          `AnnotationResult.predicted_labels`
#                                          DataFrame *or* a hand-built
#                                          R data.frame of the same
#                                          shape. Testable without
#                                          Python + CellTypist
#                                          installed.
#   * `.dataset_to_anndata_minimal()`  -- bridge helper: turn the
#                                          active dataset into a
#                                          minimal AnnData via the
#                                          `anndata` R package.
#                                          Counts-layer matrix +
#                                          cell_data.
#
# Reference catalogue: CellTypist downloads model `.pkl` files on
# demand into `~/.celltypist/models/`. We surface the model name as a
# parameter and trust the Python side to fetch.
# ============================================================================

# ---- Engine run function ---------------------------------------------------

.run_celltypist_annotation <- function(dataset, params, state) {
  # Both reticulate AND CellTypist's R-side adapter (anndata) need to
  # be loadable. CellTypist itself lives only on PyPI so the install
  # hint stays generic ("install celltypist in your Python env").
  require_optional("reticulate",
                   feature = "CellTypist annotation engine",
                   source  = "CRAN")
  require_optional("anndata",
                   feature = "CellTypist annotation engine (AnnData bridge)",
                   source  = "CRAN")

  if (!isTRUE(reticulate::py_module_available("celltypist"))) {
    stop("CellTypist Python module not available in the active reticulate ",
         "environment. Install via `reticulate::py_install(\"celltypist\")` ",
         "or in your active venv: `pip install celltypist`.",
         call. = FALSE)
  }

  model            <- params$model            %||% "Immune_All_Low.pkl"
  majority_voting  <- isTRUE(params$majority_voting %||% FALSE)
  over_clustering  <- params$over_clustering
  min_score        <- params$min_score        %||% 0.0
  cluster_field    <- params$cluster_field

  cells <- dataset$cell_data$cell

  # 1. Build the AnnData object CellTypist expects.
  ad <- .dataset_to_anndata_minimal(dataset)

  # 2. Run CellTypist. Python: `celltypist.annotate(...)`.
  ct <- reticulate::import("celltypist", delay_load = TRUE)
  args <- list(model = model, majority_voting = majority_voting)
  if (majority_voting && !is.null(over_clustering) && nzchar(over_clustering)) {
    # CellTypist accepts an arbitrary `obs` column name for the
    # over-clustering reference. We require it be a real metadata
    # column on the dataset.
    if (!over_clustering %in% colnames(ad$obs)) {
      stop(sprintf(
        "CellTypist: over_clustering column '%s' not present in obs.",
        over_clustering), call. = FALSE)
    }
    args$over_clustering <- over_clustering
  }
  ann_result <- do.call(ct$annotate, c(list(filename = ad), args))

  # 3. Pull the per-cell predicted_labels DataFrame down to R.
  predicted_py <- ann_result$predicted_labels
  predicted_r  <- reticulate::py_to_r(predicted_py)
  # `predicted_r` is a data.frame with rownames = cells. CellTypist
  # column names vary by mode:
  #   - default       : `predicted_labels`, `conf_score`
  #   - with voting   : adds `over_clustering`, `majority_voting`
  predicted_r$cell <- rownames(predicted_r)

  cluster_vec <- if (is.null(cluster_field) || !nzchar(cluster_field)) NULL
                 else as.character(get_metadata(dataset, cluster_field))

  .celltypist_to_engine_output(
    predicted_r,
    cells              = cells,
    model              = model,
    majority_voting    = majority_voting,
    min_score          = min_score,
    cluster_field_used = cluster_field,
    cluster_vec        = cluster_vec
  )
}

# Build an AnnData object that CellTypist can ingest. CellTypist
# expects log-normalised expression (it explicitly checks `X.max() < 10`
# heuristically), so we prefer the `data` layer when present.
.dataset_to_anndata_minimal <- function(dataset) {
  be <- as_expression_backend(dataset$expression)
  layers <- backend_available_layers(be)
  layer <- if ("data" %in% layers) "data" else backend_default_layer(be)
  M <- backend_as_matrix(be, layer = layer)
  if (is.null(rownames(M)))
    rownames(M) <- backend_genes(be, layer = layer)
  if (is.null(colnames(M)))
    colnames(M) <- dataset$cell_data$cell

  # AnnData stores cells as rows; transpose the genes x cells matrix.
  X <- t(as.matrix(M))
  obs <- dataset$cell_data
  rownames(obs) <- obs$cell

  ad_fn <- get("AnnData", envir = asNamespace("anndata"))
  ad_fn(X = X, obs = obs)
}

#' Pure schema converter: CellTypist `predicted_labels` -> engine output.
#'
#' Accepts a data.frame with at minimum a `cell` column and a
#' `predicted_labels` (or `majority_voting`) column. Returns the engine
#' output contract documented in R/annotation_registry.R.
#'
#' Factored out so the schema mapping can be regression-tested without
#' Python + CellTypist installed.
#'
#' @param df                 data.frame with `cell` + `predicted_labels`
#'                           (default mode) and optionally `conf_score`,
#'                           `majority_voting`, `over_clustering`.
#' @param cells              character(n_cells) -- canonical cell order.
#' @param model              character(1) -- stamped on `reference_source`.
#' @param majority_voting    logical(1) -- if TRUE and `majority_voting`
#'                           column is present, that wins as the label.
#' @param min_score          numeric(1) -- cells with `conf_score`
#'                           below this become "Unknown" with score 0.
#'                           Ignored if `conf_score` is missing.
#' @param cluster_field_used character(1) or NULL.
#' @param cluster_vec        character(n_cells) or NULL.
.celltypist_to_engine_output <- function(df,
                                         cells,
                                         model              = NA_character_,
                                         majority_voting    = FALSE,
                                         min_score          = 0.0,
                                         cluster_field_used = NULL,
                                         cluster_vec        = NULL) {
  if (!is.data.frame(df))
    stop(".celltypist_to_engine_output: `df` must be a data.frame.",
         call. = FALSE)
  if (!"cell" %in% names(df))
    stop(".celltypist_to_engine_output: `df` needs a `cell` column.",
         call. = FALSE)

  label_col <- if (isTRUE(majority_voting) && "majority_voting" %in% names(df))
                 "majority_voting"
               else if ("predicted_labels" %in% names(df))
                 "predicted_labels"
               else
                 stop(".celltypist_to_engine_output: no `predicted_labels` ",
                      "(or `majority_voting`) column in df. Have: ",
                      paste(names(df), collapse = ", "), call. = FALSE)

  idx <- match(cells, df$cell)
  if (any(is.na(idx))) {
    n_missing <- sum(is.na(idx))
    stop(sprintf(
      ".celltypist_to_engine_output: %d cells missing from CellTypist result.",
      n_missing), call. = FALSE)
  }

  cell_labels <- as.character(df[[label_col]][idx])
  cell_scores <- if ("conf_score" %in% names(df))
                   as.numeric(df$conf_score[idx])
                 else
                   rep(NA_real_, length(cells))

  if (is.finite(min_score) && min_score > 0 && "conf_score" %in% names(df)) {
    low <- !is.na(cell_scores) & cell_scores < min_score
    cell_labels[low] <- "Unknown"
    cell_scores[low] <- 0
  }

  cluster_summary <- if (!is.null(cluster_vec) &&
                         length(cluster_vec) == length(cells)) {
    .summarise_per_cluster_celltypist(
      cluster_vec = as.character(cluster_vec),
      labels      = cell_labels,
      scores      = cell_scores)
  } else NULL

  list(
    cell                   = cells,
    cell_labels            = cell_labels,
    cell_scores            = cell_scores,
    alt_labels             = NULL,
    cluster_summary        = cluster_summary,
    cluster_field_used     = cluster_field_used,
    n_clusters_at_creation = if (is.null(cluster_vec)) 0L
                             else length(unique(as.character(cluster_vec))),
    reference_source       = sprintf(
      "CellTypist:%s%s", model,
      if (isTRUE(majority_voting)) ":majority_voting" else ""
    )
  )
}

# Internal: top-label / fraction summary per cluster.
.summarise_per_cluster_celltypist <- function(cluster_vec, labels, scores) {
  cluster_ids <- sort(unique(cluster_vec))
  out <- vector("list", length(cluster_ids))
  for (i in seq_along(cluster_ids)) {
    cl <- cluster_ids[i]
    in_cl <- which(cluster_vec == cl)
    tab <- sort(table(labels[in_cl]), decreasing = TRUE)
    top  <- if (length(tab)) names(tab)[1] else NA_character_
    frac <- if (length(tab)) as.numeric(tab[1] / length(in_cl)) else NA_real_
    out[[i]] <- data.frame(
      cluster      = cl,
      top_label    = top,
      top_fraction = frac,
      mean_score   = mean(scores[in_cl], na.rm = TRUE),
      n_cells      = length(in_cl),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}
