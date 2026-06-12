# ============================================================================
# Slingshot trajectory backend
# ----------------------------------------------------------------------------
# Bioconductor `slingshot` infers lineages by:
#   1. Building a minimum spanning tree (MST) over cluster centroids
#      in a chosen reduced-dim space (PCA / UMAP / DiffMap).
#   2. Fitting a principal curve along each branch of the MST.
#   3. Reporting per-cell pseudotime along each lineage (`pseudotime`)
#      and per-cell weights (`curveWeights`).
#
# Inputs we need:
#   * `reducedDim` matrix:  `cells x dims`. Extracted from the dataset
#     via `get_embedding()`. PCA is preferred; UMAP / DiffMap work too.
#   * `clusterLabels`:      character() per-cell cluster ids. We use
#     `params$cluster_field` (or `params$root_field` as fallback).
#   * `start.clus`:         optional starting cluster (`params$root_group`).
#
# Output we surface:
#   * `pseudotime`         : per-cell pseudotime in [0, 1]. When
#                             slingshot returns multiple lineages, we
#                             aggregate using the inverse-distance
#                             weights it ships with the result, so cells
#                             that belong to a single lineage get that
#                             lineage's value and cells at branch points
#                             get a weighted mean.
#   * `n_lineages`         : integer count of inferred lineages.
#   * `method_details`     : raw per-lineage matrix (`lineage_psts`) +
#                             curveWeights, for downstream introspection.
#
# Why we don't store the SlingshotDataSet itself: it pins a copy of the
# reduced-dim matrix and the clustering, which is a memory tax for a
# Shiny session that may re-run several methods. The pure converter
# below extracts everything we need and discards the heavy object.
# ============================================================================

.run_slingshot_trajectory <- function(dataset, params) {
  require_optional(c("slingshot", "SingleCellExperiment"),
                   feature = "trajectory backend 'slingshot'",
                   source  = "Bioconductor")

  red <- params$reduction %||% dataset$default_reduction %||%
         (dataset$reductions %||% character())[1]
  emb <- get_embedding(dataset, red)
  if (is.null(emb))
    stop(sprintf("Slingshot: reduction '%s' is not available.", red %||% ""),
         call. = FALSE)
  # `get_embedding` returns a 2-column data.frame (`x`, `y`); slingshot
  # accepts any cells x ndim matrix.
  rd <- as.matrix(emb[, c("x", "y"), drop = FALSE])
  colnames(rd) <- c(paste0(red, "_1"), paste0(red, "_2"))

  cluster_field <- params$cluster_field %||% params$root_field
  if (is.null(cluster_field) || !nzchar(cluster_field))
    stop("Slingshot needs `cluster_field` (or `root_field`) -- the column ",
         "of categorical cluster labels.", call. = FALSE)
  cl <- get_metadata(dataset, cluster_field)
  if (is.null(cl))
    stop(sprintf("Slingshot: cluster_field '%s' is not in the dataset.",
                 cluster_field), call. = FALSE)
  cl <- as.character(cl)
  if (length(unique(cl)) < 2L)
    stop("Slingshot needs at least 2 distinct clusters; got 1.",
         call. = FALSE)

  start_clus <- params$root_group
  if (!is.null(start_clus) && !is.na(start_clus) &&
      nzchar(as.character(start_clus))) {
    if (!as.character(start_clus) %in% cl)
      stop(sprintf("Slingshot: start cluster '%s' not present in '%s'.",
                   start_clus, cluster_field), call. = FALSE)
  } else {
    start_clus <- NULL
  }

  sling_fn <- get("slingshot", envir = asNamespace("slingshot"))
  sds <- sling_fn(data = rd, clusterLabels = cl,
                  start.clus = start_clus)

  .slingshot_to_pseudotime(
    sds,
    cells          = dataset$cell_data$cell,
    cluster_field  = cluster_field,
    start_clus     = start_clus,
    reduction_used = red
  )
}

#' Convert a `SlingshotDataSet` into the canonical trajectory result.
#'
#' Pure: accepts either a real `SlingshotDataSet` *or* a minimal stand-in
#' list with the same accessors so the schema mapping can be regression-
#' tested without the heavy Bioc package.
#'
#' The stand-in must provide either:
#'   * `slingPseudotime(sds)` and `slingCurveWeights(sds)`-equivalents
#'     as numeric matrices `cells x lineages`, OR
#'   * top-level list elements `pseudotime` and `curveWeights` of the
#'     same shape. The converter prefers the slingshot accessors but
#'     falls back to the list elements.
#'
#' @param sds             a slingshot result (or stand-in list).
#' @param cells           character(n_cells): per-cell barcode ids.
#' @param cluster_field   character(1): name of the cluster column we used.
#' @param start_clus      character(1) or NULL.
#' @param reduction_used  character(1) reduction id.
.slingshot_to_pseudotime <- function(sds, cells, cluster_field,
                                     start_clus = NULL,
                                     reduction_used = NA_character_) {
  pst <- .slingshot_extract_matrix(sds, "pseudotime",
                                   accessor = "slingPseudotime")
  wts <- .slingshot_extract_matrix(sds, "curveWeights",
                                   accessor = "slingCurveWeights")
  if (is.null(pst))
    stop(".slingshot_to_pseudotime: no `pseudotime` matrix on the result.",
         call. = FALSE)

  if (nrow(pst) != length(cells))
    stop(sprintf(
      ".slingshot_to_pseudotime: nrow(pseudotime) (%d) != length(cells) (%d).",
      nrow(pst), length(cells)), call. = FALSE)

  n_lineages <- ncol(pst)
  agg_pt <- if (is.null(wts) || any(dim(wts) != dim(pst))) {
    # No usable weights: average over non-NA lineage values per cell.
    apply(pst, 1, function(row) {
      r <- row[!is.na(row)]
      if (length(r) == 0L) NA_real_ else mean(r)
    })
  } else {
    # Weighted mean across lineages, ignoring NAs cell-by-cell.
    vapply(seq_len(nrow(pst)), function(i) {
      pi <- pst[i, ]; wi <- wts[i, ]
      keep <- !is.na(pi) & !is.na(wi) & wi > 0
      if (!any(keep)) return(NA_real_)
      sum(pi[keep] * wi[keep]) / sum(wi[keep])
    }, numeric(1))
  }

  list(
    pseudotime     = rescale01(agg_pt),
    cell           = cells,
    source         = "slingshot",
    reduction_used = reduction_used,
    root_field     = cluster_field,
    root_group     = if (is.null(start_clus)) NA_character_
                     else as.character(start_clus),
    metadata_field = NA_character_,
    n_lineages     = as.integer(n_lineages),
    method_details = list(
      lineage_psts   = pst,
      curve_weights  = wts,
      cluster_field  = cluster_field
    )
  )
}

# Internal: extract a `cells x lineages` matrix from either a real
# SlingshotDataSet (via the package's accessor) or a stand-in list.
.slingshot_extract_matrix <- function(sds, list_key, accessor) {
  if (has_optional("slingshot")) {
    fn <- tryCatch(get(accessor, envir = asNamespace("slingshot")),
                   error = function(e) NULL)
    if (!is.null(fn)) {
      m <- tryCatch(fn(sds), error = function(e) NULL)
      if (!is.null(m)) return(as.matrix(m))
    }
  }
  if (is.list(sds) && !is.null(sds[[list_key]])) {
    return(as.matrix(sds[[list_key]]))
  }
  NULL
}
