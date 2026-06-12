# ============================================================================
# Regulon engine: AUCell
# ----------------------------------------------------------------------------
# Two implementations registered as separate engines in REGULON_ENGINES():
#
#   .run_aucell_pure_r_regulons()  -- `aucell_r`. No optional deps. Used
#                                     for the mock dataset and CI.
#                                     Numerically equivalent to the
#                                     Bioconductor version on small
#                                     dense matrices; slower on > 50k
#                                     cells.
#   .run_aucell_bioc_regulons()    -- `aucell`. Wraps
#                                     `AUCell::AUCell_calcAUC()`.
#                                     Sparse-aware, fast.
#
# Both share the pure converter `.aucell_to_regulon_engine_output()`
# which lifts an n_cells x n_regulons matrix into the engine contract.
#
# AUCell math (Aibar et al, 2017):
#   For each cell c:
#     1. Rank genes by expression (1 = highest). Ties = average rank.
#     2. For each regulon r with target genes T_r and rank threshold T
#        (= top_n = top_n_fraction * n_genes_total):
#          AUC_T  = sum_{i in T_r, rank_i <= T} (T - rank_i + 1) / |T_r|
#          AUC_max = (T * |T_r| - |T_r| * (|T_r| - 1) / 2) / |T_r|
#                  if |T_r| <= T, else T * (T + 1) / (2 * |T_r|)
#          AUC_norm = AUC_T / AUC_max   (in [0, 1])
#   Returns AUC_norm as the per-(cell, regulon) score.
#
# The normalisation here matches the AUCell paper. Higher = more of the
# regulon's targets are concentrated in the top-T expressed genes of
# that cell.
# ============================================================================

# ---- Pure-R implementation ------------------------------------------------

.run_aucell_pure_r_regulons <- function(dataset, regulon_set, params) {
  top_n_fraction <- params$top_n_fraction %||% 0.05
  if (!is.numeric(top_n_fraction) || length(top_n_fraction) != 1L ||
      top_n_fraction <= 0 || top_n_fraction > 1)
    stop("aucell_r: `top_n_fraction` must be in (0, 1]. Got ",
         top_n_fraction, call. = FALSE)

  cells   <- dataset$cell_data$cell
  be      <- as_expression_backend(dataset$expression)
  layers  <- backend_available_layers(be)
  layer   <- if ("data" %in% layers) "data" else backend_default_layer(be)
  expr    <- backend_as_matrix(be, layer = layer)  # genes x cells
  if (is.null(rownames(expr)))
    rownames(expr) <- backend_genes(be, layer = layer)
  if (is.null(colnames(expr))) colnames(expr) <- cells

  targets_by_tf <- regulon_set_as_target_list(regulon_set)
  result <- .aucell_pure_r(expr_mat = expr,
                           regulons = targets_by_tf,
                           top_n_fraction = top_n_fraction)
  .aucell_to_regulon_engine_output(
    auc_matrix = result$auc_matrix,
    cells      = cells,
    regulon_ids = names(targets_by_tf),
    warnings   = result$warnings
  )
}

# Internal: dense AUCell on a genes-by-cells matrix. Pure function.
# Exposed for unit testing.
.aucell_pure_r <- function(expr_mat, regulons, top_n_fraction = 0.05) {
  stopifnot(is.matrix(expr_mat) || methods::is(expr_mat, "Matrix"))
  n_genes_total <- nrow(expr_mat)
  n_cells       <- ncol(expr_mat)
  top_n <- max(1L, as.integer(round(top_n_fraction * n_genes_total)))
  gene_names <- rownames(expr_mat)
  if (is.null(gene_names))
    stop(".aucell_pure_r: `expr_mat` needs row names (gene symbols).",
         call. = FALSE)

  # Convert to dense once for the rank step; dgCMatrix doesn't support
  # apply() neatly and AUCell pure-R is the *small-data* path. The
  # Bioconductor engine handles sparse properly.
  if (!is.matrix(expr_mat)) expr_mat <- as.matrix(expr_mat)

  # genes x cells matrix of ranks (1 = highest expression in that cell)
  ranks_mat <- apply(expr_mat, 2, function(col)
    rank(-col, ties.method = "average"))
  rownames(ranks_mat) <- gene_names

  auc_mat <- matrix(0, nrow = n_cells, ncol = length(regulons),
                    dimnames = list(colnames(expr_mat), names(regulons)))
  warnings <- character()

  for (j in seq_along(regulons)) {
    tf      <- names(regulons)[j]
    targets <- regulons[[j]]
    hits    <- intersect(targets, gene_names)
    n_in    <- length(hits)
    if (n_in == 0L) {
      warnings <- c(warnings, sprintf(
        "Regulon '%s': no target genes present in dataset (skipped).", tf))
      next  # leave column as zero
    }
    # |hits| x n_cells matrix of ranks for the target genes
    target_ranks <- ranks_mat[hits, , drop = FALSE]

    # Max AUC for this regulon size (paper formula).
    if (n_in <= top_n) {
      auc_max <- (top_n * n_in - n_in * (n_in - 1) / 2) / n_in
    } else {
      auc_max <- top_n * (top_n + 1) / (2 * n_in)
    }
    if (auc_max <= 0) {
      # Degenerate (top_n == 0 or n_in == 0 case handled above).
      next
    }

    auc_raw <- colSums(pmax(0, top_n - target_ranks + 1) *
                       (target_ranks <= top_n)) / n_in
    auc_mat[, j] <- auc_raw / auc_max
  }

  list(auc_matrix = auc_mat, warnings = warnings)
}

# ---- Bioconductor wrapper -------------------------------------------------

.run_aucell_bioc_regulons <- function(dataset, regulon_set, params) {
  require_optional("AUCell",
                   feature = "AUCell (Bioconductor) regulon scoring",
                   source  = "Bioconductor")

  top_n_fraction <- params$top_n_fraction %||% 0.05

  cells <- dataset$cell_data$cell
  be    <- as_expression_backend(dataset$expression)
  layers <- backend_available_layers(be)
  layer <- if ("data" %in% layers) "data" else backend_default_layer(be)
  expr  <- backend_as_matrix(be, layer = layer)

  build_rk <- get("AUCell_buildRankings", envir = asNamespace("AUCell"))
  calc_auc <- get("AUCell_calcAUC",       envir = asNamespace("AUCell"))

  rk <- build_rk(exprMat = expr, plotStats = FALSE, verbose = FALSE)
  aucMaxRank <- max(1L, as.integer(round(top_n_fraction * nrow(expr))))

  auc_obj <- calc_auc(geneSets   = regulon_set_as_target_list(regulon_set),
                      rankings   = rk,
                      aucMaxRank = aucMaxRank,
                      verbose    = FALSE)

  # AUCell returns an `aucellResults` S4. `getAUC()` -> matrix
  # (regulons x cells). Normalize to cells x regulons.
  get_auc <- get("getAUC", envir = asNamespace("AUCell"))
  auc <- t(as.matrix(get_auc(auc_obj)))

  .aucell_to_regulon_engine_output(
    auc_matrix  = auc,
    cells       = cells,
    regulon_ids = colnames(auc),
    warnings    = character()
  )
}

# ---- Pure schema converter ------------------------------------------------

#' Schema converter: cells x regulons AUC matrix -> engine output.
#'
#' Aligns rows to the dataset's cell order. Used by both engines so the
#' shape coming out of `run_fn` is identical regardless of backend.
#' Tested without AUCell installed.
#'
#' @param auc_matrix  numeric matrix [n_cells x n_regulons]; rownames
#'                    must be a subset of `cells`. Rows are reordered
#'                    to match `cells`.
#' @param cells       character() canonical dataset cell order.
#' @param regulon_ids character() column names (TFs).
#' @param warnings    character() optional warnings to forward.
.aucell_to_regulon_engine_output <- function(auc_matrix, cells,
                                             regulon_ids,
                                             warnings = character()) {
  if (!is.numeric(auc_matrix) || !is.matrix(auc_matrix))
    stop(".aucell_to_regulon_engine_output: `auc_matrix` must be a ",
         "numeric matrix.", call. = FALSE)
  if (is.null(rownames(auc_matrix)))
    rownames(auc_matrix) <- cells
  if (is.null(colnames(auc_matrix)))
    colnames(auc_matrix) <- regulon_ids
  # Align rows.
  idx <- match(cells, rownames(auc_matrix))
  if (any(is.na(idx))) {
    n_missing <- sum(is.na(idx))
    stop(sprintf(
      ".aucell_to_regulon_engine_output: %d cells missing from AUC matrix.",
      n_missing), call. = FALSE)
  }
  auc <- auc_matrix[idx, , drop = FALSE]
  rownames(auc) <- cells
  list(
    cell        = cells,
    regulon_ids = colnames(auc),
    auc_matrix  = auc,
    warnings    = as.character(warnings)
  )
}
