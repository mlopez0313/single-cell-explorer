# ============================================================================
# Marker computation
# ----------------------------------------------------------------------------
# For v1 this is a tiny, dependency-free implementation that ranks the genes
# in `dataset$expression` for each level of a categorical metadata field.
#
# It is intentionally a stand-in for `Seurat::FindAllMarkers()` /
# `scanpy.tl.rank_genes_groups()`. The output columns are a deliberate
# subset of the Seurat schema so a future swap to the real implementation
# is mechanical:
#
#   gene  | group | avg_log2FC | pct_in | pct_out | p_value
#
# Statistics computed (no extra packages):
#   - avg_log2FC : log2(mean(in_group)+1) - log2(mean(out_group)+1)
#   - pct_in     : fraction of in-group cells with expression > 0
#   - pct_out    : fraction of out-group cells with expression > 0
#   - p_value    : two-sample test, chosen by `test`:
#                    "wilcox" (default, Seurat-like) -> stats::wilcox.test
#                    "t"                              -> stats::t.test (Welch)
# ============================================================================

#' Compute mock markers for each group of a metadata field.
#'
#' @param dataset         a dataset list (see `dataset_schema()`)
#' @param grouping_field  metadata column to group cells by
#' @param group_filter    optional character() restricting which groups to score
#' @param top_n           keep only the top N genes per group (by avg_log2FC).
#'                        `Inf` (default) keeps all genes.
#' @param test            one of "wilcox" (default) or "t"
#'
#' @return data.frame with columns `group`, `gene`, `avg_log2FC`, `pct_in`,
#'   `pct_out`, `p_value`. Sorted by group, then descending `avg_log2FC`.
#'   Returns `NULL` if the field/genes are unavailable.
compute_markers <- function(dataset, grouping_field,
                            group_filter = NULL, top_n = Inf,
                            test = c("wilcox", "t")) {
  test <- match.arg(test)
  pval_fn <- switch(test,
    "wilcox" = function(x, y) tryCatch(stats::wilcox.test(x, y, exact = FALSE)$p.value,
                                       error = function(e) NA_real_,
                                       warning = function(w) suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value)),
    "t"      = function(x, y) tryCatch(stats::t.test(x, y)$p.value,
                                       error = function(e) NA_real_)
  )
  meta  <- get_metadata(dataset, grouping_field)
  genes <- available_genes(dataset)
  if (is.null(meta) || length(genes) == 0L) return(NULL)

  groups <- unique(meta)
  if (!is.null(group_filter)) groups <- intersect(groups, group_filter)
  if (length(groups) == 0L) return(NULL)

  # Heartbeat messages are gated on a session option so the test suite
  # (and any caller that just wants the data.frame) can suppress them
  # without redirecting stderr.
  hb <- isTRUE(getOption("sce.marker_progress", TRUE))
  emit <- if (hb) function(...) message(sprintf(...)) else function(...) NULL
  t_all_start <- Sys.time()

  # ---- Materialise the expression matrix ONCE -------------------------
  # The previous implementation called `get_gene_expression(dataset, g)`
  # per gene -- which for a `dgCMatrix`-backed Seurat object is row
  # indexing into a column-compressed matrix, an O(nnz) walk EVERY
  # call. At 18 340 genes x 8 381 cells that's effectively unbounded
  # (we saw it hang for hours).
  #
  # `backend_as_matrix()` returns the underlying matrix in
  # genes-as-rows orientation for both the sparse / Seurat backend
  # (returns the dgCMatrix itself, no copy) and the in-memory backend
  # (returns a dense base matrix). Then every per-group stat is a
  # single vectorised `Matrix::rowMeans(...)` call -- thousands of
  # times faster.
  emit("compute_markers: materialising expression matrix (%d genes x %d cells)",
       length(genes), length(meta))
  t_mat <- Sys.time()
  be <- as_expression_backend(dataset$expression)
  M  <- tryCatch(backend_as_matrix(be, layer = NULL),
                 error = function(e) NULL)
  if (is.null(M) || nrow(M) == 0L || ncol(M) == 0L) {
    emit("compute_markers: backend_as_matrix() unusable, aborting")
    return(NULL)
  }
  # Align row order with `genes` so downstream output uses the same
  # gene labels the rest of the app sees. `match()` is O(n).
  if (!is.null(rownames(M)) && !identical(rownames(M), genes)) {
    keep <- match(genes, rownames(M))
    keep <- keep[!is.na(keep)]
    if (!length(keep)) return(NULL)
    M     <- M[keep, , drop = FALSE]
    genes <- rownames(M)
  }
  emit("compute_markers: matrix ready in %.1fs (class=%s, %s)",
       as.numeric(difftime(Sys.time(), t_mat, units = "secs")),
       class(M)[1L],
       if (methods::is(M, "sparseMatrix"))
         sprintf("nnz=%d, density=%.3f%%", length(M@x),
                 100 * length(M@x) / (nrow(M) * ncol(M)))
       else "dense")

  # Use Matrix::rowMeans (S4) for both sparse and dense -- a single
  # C call across all genes.
  rm_ <- function(x) {
    if (methods::is(x, "Matrix")) Matrix::rowMeans(x) else base::rowMeans(x)
  }

  rows <- list()
  for (gi in seq_along(groups)) {
    grp <- groups[gi]
    in_idx  <- which(meta == grp)
    out_idx <- which(meta != grp)
    if (length(in_idx) < 2L || length(out_idx) < 2L) {
      emit("compute_markers: group %d/%d (%s) skipped (in=%d, out=%d -- need >=2 each)",
           gi, length(groups), as.character(grp), length(in_idx), length(out_idx))
      next
    }
    emit("compute_markers: group %d/%d (%s, in=%d, out=%d) starting",
         gi, length(groups), as.character(grp), length(in_idx), length(out_idx))
    t_grp <- Sys.time()

    # ---- Vectorised per-gene stats (across ALL genes at once) -------
    M_in   <- M[, in_idx,  drop = FALSE]
    M_out  <- M[, out_idx, drop = FALSE]
    m_in   <- rm_(M_in)
    m_out  <- rm_(M_out)
    # `> 0` on a dgCMatrix returns lgCMatrix -- rowMeans on it gives
    # the fraction of cells with non-zero expression. No dense copy.
    pct_in   <- rm_(M_in  > 0)
    pct_out  <- rm_(M_out > 0)
    avg_l2fc <- log2(m_in + 1) - log2(m_out + 1)

    emit("compute_markers: group %d/%d (%s) stats done in %.1fs",
         gi, length(groups), as.character(grp),
         as.numeric(difftime(Sys.time(), t_grp, units = "secs")))

    # ---- Sort + take top_n, THEN compute p-values --------------------
    # Old code computed wilcox.test for every (gene, group) pair --
    # 256k tests for PBMC 8k -- then threw away all but top_n. Now
    # we test only the survivors (top_n per group = ~6).
    ord <- order(-avg_l2fc)
    keep <- if (is.finite(top_n)) ord[seq_len(min(top_n, length(ord)))] else ord

    t_pval <- Sys.time()
    pvals <- vapply(keep, function(j) {
      # Row indexing into a dgCMatrix is the one slow op left, but
      # we do it `length(keep)` (e.g. 6) times per group instead of
      # n_genes times -- a 3000x reduction at top_n = 6.
      v_in  <- as.numeric(M_in[j, , drop = TRUE])
      v_out <- as.numeric(M_out[j, , drop = TRUE])
      pval_fn(v_in, v_out)
    }, numeric(1))
    emit("compute_markers: group %d/%d (%s) p-values for top %d in %.1fs",
         gi, length(groups), as.character(grp), length(keep),
         as.numeric(difftime(Sys.time(), t_pval, units = "secs")))

    rows[[length(rows) + 1L]] <- data.frame(
      group      = grp,
      gene       = genes[keep],
      avg_log2FC = avg_l2fc[keep],
      pct_in     = pct_in[keep],
      pct_out    = pct_out[keep],
      p_value    = pvals,
      stringsAsFactors = FALSE
    )

    emit("compute_markers: group %d/%d (%s) total %.1fs",
         gi, length(groups), as.character(grp),
         as.numeric(difftime(Sys.time(), t_grp, units = "secs")))
  }
  emit("compute_markers: all groups done in %.1fs (%d rows)",
       as.numeric(difftime(Sys.time(), t_all_start, units = "secs")),
       sum(vapply(rows, NROW, integer(1))))
  if (length(rows) == 0L) return(NULL)

  df <- do.call(rbind, rows)
  df <- df[order(df$group, -df$avg_log2FC), , drop = FALSE]
  # `top_n` was already enforced per-group above, but a final guard
  # keeps the schema stable if a future caller skips the early sort
  # (e.g. if test="t" picked a different ordering criterion).
  if (is.finite(top_n)) {
    df <- do.call(rbind, by(df, df$group, utils::head, n = top_n))
  }
  rownames(df) <- NULL
  df
}
