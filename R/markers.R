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

  # Pre-fetch every gene vector once.
  expr_mat <- vapply(genes, function(g) get_gene_expression(dataset, g),
                     FUN.VALUE = numeric(length(meta)))

  rows <- list()
  for (grp in groups) {
    in_idx  <- which(meta == grp)
    out_idx <- which(meta != grp)
    if (length(in_idx) < 2L || length(out_idx) < 2L) next
    for (j in seq_along(genes)) {
      v   <- expr_mat[, j]
      m_in  <- mean(v[in_idx])
      m_out <- mean(v[out_idx])
      pval  <- pval_fn(v[in_idx], v[out_idx])
      rows[[length(rows) + 1L]] <- data.frame(
        group      = grp,
        gene       = genes[j],
        avg_log2FC = log2(m_in + 1) - log2(m_out + 1),
        pct_in     = mean(v[in_idx]  > 0),
        pct_out    = mean(v[out_idx] > 0),
        p_value    = pval,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) return(NULL)

  df <- do.call(rbind, rows)
  df <- df[order(df$group, -df$avg_log2FC), , drop = FALSE]
  if (is.finite(top_n)) {
    df <- do.call(rbind, by(df, df$group, utils::head, n = top_n))
  }
  rownames(df) <- NULL
  df
}
