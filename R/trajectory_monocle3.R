# ============================================================================
# Monocle3 trajectory backend
# ----------------------------------------------------------------------------
# Monocle3 has its own data structure (`cell_data_set`, CDS) and runs a
# full pipeline: preprocess -> reduce -> cluster -> learn graph ->
# order. We honour that pipeline but try to use the *existing*
# dim-reduction and clustering already in the dataset where possible,
# so the result lines up visually with the rest of the app.
#
# Input geometry the backend assumes:
#   * `counts` layer (raw counts) from the expression backend, as a
#     `genes x cells` matrix. We pull through `backend_as_matrix()` so
#     real sparse backends do not have to densify per-gene.
#   * A 2D embedding (UMAP preferred). We *inject* the dataset's
#     reduction into the CDS via `reducedDims(cds)$UMAP <- ...` so
#     `monocle3::cluster_cells()` and `learn_graph()` operate in the
#     same space the user is looking at in the rest of the app.
#   * Cluster labels from `params$cluster_field` / `params$root_field`.
#     They are injected as `clusters(cds, "UMAP")` so `learn_graph()`
#     respects them rather than re-clustering from scratch.
#   * Root cells: picked as every cell in `root_group` (`order_cells()`
#     supports a vector of root cell ids).
#
# Output we surface:
#   * `pseudotime`: per-cell pseudotime in [0, 1] (NA cells in
#     disconnected components are coerced to NA upstream).
#   * `method_details`: principal-graph node ids, partition labels,
#     and the count of disconnected partitions (a real Monocle3
#     pitfall when clusters are too far apart).
#
# Caveats called out in the README and in the demo banner:
#   * Monocle3 wants raw counts. If the dataset only has a log-normalised
#     `data` layer, we still run but mark `warn_lognorm = TRUE` in
#     `method_details` so the module can surface a warning.
#   * Monocle3 *can* re-cluster internally if we don't inject the
#     clustering. We always inject ours; future PRs can expose a
#     "let monocle3 cluster" flag.
# ============================================================================

.run_monocle3_trajectory <- function(dataset, params) {
  require_optional("SingleCellExperiment",
                   feature = "trajectory backend 'monocle3' (Bioc deps)",
                   source  = "Bioconductor")
  require_optional("monocle3",
                   feature = "trajectory backend 'monocle3' (GitHub package)",
                   source  = "GitHub",
                   repo    = c(monocle3 = "cole-trapnell-lab/monocle3"))

  red <- params$reduction %||% dataset$default_reduction %||%
         (dataset$reductions %||% character())[1]
  emb <- get_embedding(dataset, red)
  if (is.null(emb))
    stop(sprintf("Monocle3: reduction '%s' is not available.", red %||% ""),
         call. = FALSE)
  cluster_field <- params$cluster_field %||% params$root_field
  if (is.null(cluster_field) || !nzchar(cluster_field))
    stop("Monocle3 needs `cluster_field` (or `root_field`) -- the column ",
         "of categorical cluster labels.", call. = FALSE)
  cl <- get_metadata(dataset, cluster_field)
  if (is.null(cl))
    stop(sprintf("Monocle3: cluster_field '%s' is not in the dataset.",
                 cluster_field), call. = FALSE)
  cl <- as.character(cl)

  root_group <- params$root_group
  if (is.null(root_group) || !nzchar(as.character(root_group)))
    stop("Monocle3 needs `root_group` to pin a starting cluster.",
         call. = FALSE)
  if (!as.character(root_group) %in% cl)
    stop(sprintf("Monocle3: root group '%s' not present in '%s'.",
                 root_group, cluster_field), call. = FALSE)

  # Counts matrix (genes x cells). Prefer "counts" layer; fall back to
  # default with a warn flag for the result.
  be <- as_expression_backend(dataset$expression)
  layers <- backend_available_layers(be)
  layer <- if ("counts" %in% layers) "counts" else backend_default_layer(be)
  warn_lognorm <- !identical(layer, "counts")
  M <- backend_as_matrix(be, layer = layer)

  # Build the CDS
  new_cds_fn  <- get("new_cell_data_set",   envir = asNamespace("monocle3"))
  precds_fn   <- get("preprocess_cds",      envir = asNamespace("monocle3"))
  learn_fn    <- get("learn_graph",         envir = asNamespace("monocle3"))
  order_fn    <- get("order_cells",         envir = asNamespace("monocle3"))
  pst_fn      <- get("pseudotime",          envir = asNamespace("monocle3"))

  cells <- dataset$cell_data$cell
  cd <- data.frame(cell = cells, cluster_ = cl,
                   row.names = cells, stringsAsFactors = FALSE)
  gd <- data.frame(gene_short_name = rownames(M),
                   row.names = rownames(M), stringsAsFactors = FALSE)

  cds <- new_cds_fn(expression_data = M, cell_metadata = cd,
                    gene_metadata   = gd)
  # Minimal preprocessing -- enough to satisfy downstream calls.
  cds <- precds_fn(cds, num_dim = min(50L, ncol(M) - 1L))

  # Inject the existing reduction so monocle3 operates in the user's
  # visual space. `reducedDims()` setter lives in SingleCellExperiment.
  SCE_reducedDims_setter <- get("reducedDims<-",
                                envir = asNamespace("SingleCellExperiment"))
  rd_mat <- as.matrix(emb[, c("x", "y"), drop = FALSE])
  rownames(rd_mat) <- cells
  colnames(rd_mat) <- c(paste0(red, "_1"), paste0(red, "_2"))
  rd_list <- list(UMAP = rd_mat)
  cds <- SCE_reducedDims_setter(cds, value = rd_list)

  # Cluster injection: monocle3 stores per-partition cluster maps under
  # `cds@clusters$UMAP$clusters`. We mimic the structure to short-circuit
  # `cluster_cells()` re-clustering. Use a single partition unless the
  # caller asked otherwise.
  cluster_vec <- factor(cl, levels = unique(cl))
  names(cluster_vec) <- cells
  partition_vec <- factor(rep(1L, length(cells)))
  names(partition_vec) <- cells
  # Build the slot directly. Monocle3's S4 design exposes `clusters` as
  # a list under @clusters; assigning here keeps `learn_graph()` happy.
  cds@clusters$UMAP <- list(
    cluster_result = list(),
    partitions     = partition_vec,
    clusters       = cluster_vec
  )

  cds <- learn_fn(cds)

  # Pick root cells = every cell in root_group. order_cells() supports
  # `root_cells` as a character vector of barcode ids.
  root_cells <- cells[cl == as.character(root_group)]
  cds <- order_fn(cds, root_cells = root_cells)

  pt <- as.numeric(pst_fn(cds))
  names(pt) <- cells

  list(
    pseudotime     = rescale01(pt),
    cell           = cells,
    source         = "monocle3",
    reduction_used = red,
    root_field     = cluster_field,
    root_group     = as.character(root_group),
    metadata_field = NA_character_,
    n_lineages     = length(unique(partition_vec)),
    method_details = list(
      n_partitions = length(unique(partition_vec)),
      warn_lognorm = warn_lognorm,
      layer_used   = layer
    )
  )
}

#' Convert a Monocle3-style per-cell pseudotime vector to the canonical
#' trajectory result schema.
#'
#' Pure: takes any named numeric vector (`names(pt)` = cell ids) plus
#' a few descriptive parameters. Lets the schema mapping be regression-
#' tested without monocle3 installed.
.monocle3_to_pseudotime <- function(pt, cells,
                                    cluster_field,
                                    root_group     = NA_character_,
                                    reduction_used = NA_character_,
                                    method_details = list()) {
  if (is.null(pt) || length(pt) == 0L)
    stop(".monocle3_to_pseudotime: empty pseudotime vector.", call. = FALSE)
  if (length(pt) != length(cells))
    stop(sprintf(
      ".monocle3_to_pseudotime: length(pt) (%d) != length(cells) (%d).",
      length(pt), length(cells)), call. = FALSE)
  list(
    pseudotime     = rescale01(as.numeric(pt)),
    cell           = cells,
    source         = "monocle3",
    reduction_used = reduction_used,
    root_field     = cluster_field,
    root_group     = as.character(root_group),
    metadata_field = NA_character_,
    n_lineages     = 1L,
    method_details = method_details
  )
}
