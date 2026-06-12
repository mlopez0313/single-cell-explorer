# ============================================================================
# Pseudobulk aggregation + DE scaffolding
# ----------------------------------------------------------------------------
# Pseudobulk DE flips the cell-level Wilcoxon (in R/de.R) on its head:
# instead of treating every cell as an independent observation, we
# collapse cells from the same biological replicate ("sample") within
# the same biological group into one aggregated count vector ("a
# pseudobulk sample"). The downstream test then compares pseudobulk
# samples across groups, which respects sample-level variance and is
# the recommended approach for production scRNA DE (Squair 2021,
# Crowell 2020, Soneson 2018).
#
# This file is the *aggregation + validation seam*. It is dependency-
# free and pure. Three DE backends in R/de.R consume it:
#
#   pseudobulk_naive    -- aggregation + log2(CPM+1) + per-gene t-test on
#                          pseudobulk samples. Pure R, always available.
#                          Useful for demos / CI / quick checks.
#
#   pseudobulk_edger    -- aggregation + edgeR::DGEList -> exactTest /
#                          glmQLFit. Gated by `require_optional("edgeR")`.
#                          Considered the standard for pseudobulk DE.
#
#   pseudobulk_deseq2   -- aggregation + DESeq2::DESeq. Gated by
#                          `require_optional("DESeq2")`. Slower, more
#                          memory-hungry, but produces shrunken LFC.
#
# All three return a frame matching `empty_de_results()` exactly. The
# only practical difference visible to downstream modules is that the
# stamped `params` entry includes `sample_by` and `agg`.
#
# Required metadata for pseudobulk:
#   * `grouping_field`: e.g. cluster / cell_type / condition. Splits
#     cells into the two groups being compared.
#   * `sample_by`     : biological replicate id. Must be CROSSED with
#     `grouping_field`: each (sample, group) combo gets one pseudobulk
#     sample.
#
# Counts vs. data layers:
#   Pseudobulk DE methods are designed for RAW COUNTS. The mock dataset
#   exposes a `counts` layer alongside `data` so this code path can be
#   exercised. Real loaders MUST surface a `counts` layer (Seurat:
#   layer = "counts"; AnnData: `X` or `raw.X`; 10x: matrix.mtx). If
#   only `data` (log-normalised) is available, aggregation still works
#   but the user should be warned (we set `warn_lognorm = TRUE` in the
#   pseudobulk result and downstream modules surface it).
# ============================================================================

# ---- Public: validation -----------------------------------------------

#' Validate inputs to a pseudobulk DE call.
#'
#' Throws a friendly error if the dataset / metadata can't support
#' pseudobulk DE for the requested groups. Returns silently on success.
#'
#' @param dataset                a dataset list (`dataset_schema()`).
#' @param grouping_field         metadata column to split cells on.
#' @param group_1,group_2        two distinct values in `grouping_field`.
#' @param sample_by              metadata column carrying biological
#'                               replicate ids.
#' @param min_cells_per_sample   drop pseudobulk samples assembled from
#'                               fewer than this many cells. Default 10.
#' @param min_samples_per_group  each group must have at least this many
#'                               surviving pseudobulk samples for the
#'                               downstream test to have power. Default 2
#'                               (the absolute minimum; 3+ is healthier).
validate_pseudobulk_inputs <- function(dataset, grouping_field,
                                       group_1, group_2, sample_by,
                                       min_cells_per_sample  = 10L,
                                       min_samples_per_group = 2L) {
  if (is.null(dataset)) stop("validate_pseudobulk_inputs: no dataset.",
                             call. = FALSE)
  if (is.null(grouping_field) || !nzchar(grouping_field))
    stop("validate_pseudobulk_inputs: `grouping_field` is required.",
         call. = FALSE)
  if (is.null(sample_by) || !nzchar(sample_by))
    stop("validate_pseudobulk_inputs: `sample_by` is required for ",
         "pseudobulk DE -- pass the metadata column carrying biological ",
         "replicate ids (e.g. 'sample', 'donor').", call. = FALSE)
  if (identical(grouping_field, sample_by))
    stop("validate_pseudobulk_inputs: `grouping_field` and `sample_by` ",
         "cannot be the same column (need crossed design).", call. = FALSE)

  meta_g <- get_metadata(dataset, grouping_field)
  meta_s <- get_metadata(dataset, sample_by)
  if (is.null(meta_g))
    stop("Metadata field '", grouping_field, "' not in dataset.", call. = FALSE)
  if (is.null(meta_s))
    stop("Metadata field '", sample_by, "' not in dataset.", call. = FALSE)

  in1 <- which(as.character(meta_g) == as.character(group_1))
  in2 <- which(as.character(meta_g) == as.character(group_2))
  if (length(in1) == 0L)
    stop("Group '", group_1, "' has no cells in '", grouping_field, "'.",
         call. = FALSE)
  if (length(in2) == 0L)
    stop("Group '", group_2, "' has no cells in '", grouping_field, "'.",
         call. = FALSE)

  # Cell counts per (sample, group) combo
  tab_1 <- table(meta_s[in1])
  tab_2 <- table(meta_s[in2])
  ok_1  <- sum(tab_1 >= min_cells_per_sample)
  ok_2  <- sum(tab_2 >= min_cells_per_sample)

  if (ok_1 < min_samples_per_group || ok_2 < min_samples_per_group) {
    stop(sprintf(
      paste("Pseudobulk DE needs at least %d sample(s) per group with",
            ">= %d cells.\n",
            "After filtering on '%s' (cells per sample):\n",
            "  group '%s' usable samples: %d (counts: %s)\n",
            "  group '%s' usable samples: %d (counts: %s)\n",
            "Try lowering `min_cells_per_sample` (currently %d) or use a",
            "different `sample_by`."),
      min_samples_per_group, min_cells_per_sample, sample_by,
      group_1, ok_1, paste(tab_1, collapse = ","),
      group_2, ok_2, paste(tab_2, collapse = ","),
      min_cells_per_sample), call. = FALSE)
  }
  invisible(TRUE)
}

# ---- Public: aggregation ----------------------------------------------

#' Aggregate cell-level expression into pseudobulk samples.
#'
#' Collapses cells into one row per (sample, group) combination by
#' summing (default) or averaging the per-gene expression vector.
#'
#' @param dataset           a dataset list.
#' @param grouping_field    metadata column to split cells on.
#' @param group_1,group_2   the two groups being compared.
#' @param sample_by         metadata column carrying biological replicate ids.
#' @param layer             expression layer. Defaults to `"counts"` when
#'                          available; falls back to the backend's default
#'                          layer (usually `"data"`) with `warn_lognorm = TRUE`.
#' @param agg               `"sum"` (default; correct for counts) or
#'                          `"mean"` (only used with log-normalised values).
#' @param min_cells_per_sample  drop pseudobulk samples with fewer than this
#'                              many cells. Default 10.
#'
#' @return list with shape `pseudobulk_result_v1`:
#'   \describe{
#'     \item{matrix}{numeric matrix [genes x pseudobulk_samples].}
#'     \item{sample_metadata}{data.frame: pb_sample, sample, group, n_cells.}
#'     \item{layer_used}{character(1).}
#'     \item{agg}{`"sum"` or `"mean"`.}
#'     \item{warn_lognorm}{TRUE if `counts` was unavailable and we fell back to
#'                         a log-normalised layer (caller should warn).}
#'     \item{provenance}{list: timestamp, dataset_name, grouping_field, sample_by.}
#'   }
#'
#' Returns NULL if no pseudobulk sample survives `min_cells_per_sample`.
aggregate_pseudobulk <- function(dataset, grouping_field,
                                 group_1, group_2, sample_by,
                                 layer = NULL,
                                 agg = c("sum", "mean"),
                                 min_cells_per_sample = 10L) {
  agg <- match.arg(agg)

  # Layer resolution: prefer "counts" if it exists, else fall back.
  avail_layers <- tryCatch(
    backend_available_layers(as_expression_backend(dataset$expression)),
    error = function(e) "data")
  if (is.null(layer)) {
    if ("counts" %in% avail_layers) {
      layer <- "counts"; warn_lognorm <- FALSE
    } else {
      layer <- tryCatch(backend_default_layer(as_expression_backend(dataset$expression)),
                        error = function(e) "data")
      warn_lognorm <- !identical(layer, "counts")
    }
  } else {
    warn_lognorm <- !identical(layer, "counts")
  }

  meta_g <- get_metadata(dataset, grouping_field)
  meta_s <- get_metadata(dataset, sample_by)

  # Per-cell mask
  in_pair <- as.character(meta_g) %in% c(as.character(group_1),
                                         as.character(group_2))
  if (!any(in_pair)) return(NULL)

  cell_group  <- ifelse(as.character(meta_g) == as.character(group_1),
                        as.character(group_1), as.character(group_2))
  cell_sample <- as.character(meta_s)
  pb_key      <- paste(cell_sample, cell_group, sep = "__")
  pb_key[!in_pair] <- NA_character_

  # All pseudobulk samples (sample x group combos in this pair)
  key_tab <- table(pb_key, useNA = "no")
  keep_keys <- names(key_tab)[as.integer(key_tab) >= min_cells_per_sample]
  if (length(keep_keys) == 0L) return(NULL)

  # Stable order: by group first, then sample, so columns group like-with-like.
  parts <- strsplit(keep_keys, "__", fixed = TRUE)
  sm_df <- data.frame(
    pb_sample = keep_keys,
    sample    = vapply(parts, `[`, character(1), 1L),
    group     = vapply(parts, `[`, character(1), 2L),
    n_cells   = as.integer(key_tab[keep_keys]),
    stringsAsFactors = FALSE
  )
  ord    <- order(factor(sm_df$group, levels = c(as.character(group_1),
                                                 as.character(group_2))),
                  sm_df$sample)
  sm_df  <- sm_df[ord, , drop = FALSE]
  rownames(sm_df) <- NULL
  keep_keys <- sm_df$pb_sample

  # Pull gene matrix from the backend
  be <- as_expression_backend(dataset$expression)
  genes <- backend_genes(be, layer = layer)
  if (length(genes) == 0L) return(NULL)

  pb_mat <- matrix(0, nrow = length(genes), ncol = length(keep_keys),
                   dimnames = list(genes, keep_keys))

  agg_fn <- if (identical(agg, "sum")) sum else mean
  # Group cell indices by pb_key once; per-gene per-pb_sample lookup is then O(1).
  cell_idx_by_key <- split(seq_along(pb_key), pb_key)
  for (g in genes) {
    v <- backend_get_gene(be, g, layer = layer)
    if (is.null(v)) next
    for (k in keep_keys) {
      idx <- cell_idx_by_key[[k]]
      pb_mat[g, k] <- agg_fn(v[idx])
    }
  }

  list(
    matrix          = pb_mat,
    sample_metadata = sm_df,
    layer_used      = layer,
    agg             = agg,
    warn_lognorm    = warn_lognorm,
    provenance      = list(
      timestamp       = Sys.time(),
      dataset_name    = dataset$name %||% NA_character_,
      grouping_field  = grouping_field,
      sample_by       = sample_by,
      group_1         = as.character(group_1),
      group_2         = as.character(group_2),
      min_cells_per_sample = as.integer(min_cells_per_sample)
    )
  )
}

# ---- Helpers consumed by DE backends ---------------------------------

#' Simple counts-per-million normalisation, log2 + 1 transform.
#'
#' Library size = column sum. Genes never seen in any pseudobulk
#' sample stay at 0. Used by `pseudobulk_naive`. edgeR / DESeq2
#' carry their own normalisation and ignore this helper.
pseudobulk_cpm_log2 <- function(pb_mat) {
  lib <- colSums(pb_mat)
  lib[lib == 0] <- 1 # avoid divide-by-zero; the column is all-zero anyway
  cpm <- sweep(pb_mat, 2, lib, FUN = "/") * 1e6
  log2(cpm + 1)
}

#' Compute per-gene pct.1 / pct.2 from cell-level counts.
#'
#' Pseudobulk DE methods don't natively output Seurat-style pct values,
#' but downstream modules (volcano / table filters) depend on them.
#' We compute pct.1/pct.2 once from the cell-level counts of the cells
#' that contributed to each pseudobulk sample. NA-safe.
.pseudobulk_pct <- function(dataset, in1, in2, layer = "counts") {
  be <- as_expression_backend(dataset$expression)
  genes <- backend_genes(be, layer = layer)
  pct1 <- numeric(length(genes)); pct2 <- numeric(length(genes))
  names(pct1) <- genes; names(pct2) <- genes
  for (g in genes) {
    v <- backend_get_gene(be, g, layer = layer)
    if (is.null(v)) next
    pct1[g] <- mean(v[in1] > 0)
    pct2[g] <- mean(v[in2] > 0)
  }
  list(pct.1 = pct1, pct.2 = pct2)
}

# ---- Pure converter: pseudobulk-naive output -> DE schema ----------------

#' Convert per-gene pseudobulk test stats into the canonical DE schema.
#'
#' Pure: takes a data.frame with `gene`, `avg_log2FC`, `p_val`, plus
#' optional `pct.1`/`pct.2`, and lifts it into the canonical DE schema
#' (BH-adjusts p-values, fills missing columns).
.pseudobulk_to_de_schema <- function(df, group_1, group_2,
                                     min_pct = 0) {
  if (is.null(df) || nrow(df) == 0L)
    return(empty_de_results(group_1, group_2))

  needed <- c("gene", "avg_log2FC", "p_val")
  miss   <- setdiff(needed, names(df))
  if (length(miss) > 0L)
    stop(".pseudobulk_to_de_schema: missing columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  if (is.null(df$pct.1)) df$pct.1 <- NA_real_
  if (is.null(df$pct.2)) df$pct.2 <- NA_real_

  out <- data.frame(
    gene       = as.character(df$gene),
    group_1    = as.character(group_1),
    group_2    = as.character(group_2),
    avg_log2FC = as.numeric(df$avg_log2FC),
    pct.1      = as.numeric(df$pct.1),
    pct.2      = as.numeric(df$pct.2),
    p_val      = as.numeric(df$p_val),
    stringsAsFactors = FALSE
  )
  if (!is.null(min_pct) && is.finite(min_pct) && min_pct > 0) {
    keep <- !is.na(out$pct.1) & !is.na(out$pct.2) &
      pmax(out$pct.1, out$pct.2) >= min_pct
    out <- out[keep, , drop = FALSE]
  }
  if (nrow(out) == 0L) return(empty_de_results(group_1, group_2))
  out$p_val_adj <- stats::p.adjust(out$p_val, method = "BH")
  rownames(out) <- NULL
  out
}
