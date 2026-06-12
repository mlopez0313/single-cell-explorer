# ============================================================================
# Differential Expression
# ----------------------------------------------------------------------------
# DE = compare two specific groups of cells. Conceptually distinct from
# marker detection (one group vs rest, see R/markers.R) so the logic lives
# in its own file. The Wilcoxon kernel is the same.
#
# Schema (Seurat-style; identical across every backend):
#   gene | group_1 | group_2 | avg_log2FC | pct.1 | pct.2 | p_val | p_val_adj
#
# Conventions:
#   avg_log2FC > 0   => higher in group_1
#   pct.1            => fraction of group_1 cells with expression > 0
#   pct.2            => fraction of group_2 cells with expression > 0
#   p_val_adj        => Benjamini-Hochberg over the tested genes
#
# Pre-test filter:
#   min_pct          => drop any gene whose max(pct.1, pct.2) < min_pct
#                       before testing. This is the only compute-time filter;
#                       log2FC / p_val_adj thresholds are display-only.
# ----------------------------------------------------------------------------
# Backends
# ----------------------------------------------------------------------------
# `compute_de()` is a thin orchestrator that dispatches to a backend:
#
#   backend = "wilcox_r"  -- pure-R per-gene loop using stats::wilcox.test
#                            (or stats::t.test if test = "t"). Always
#                            available. Slow but dependency-free.
#
#   backend = "presto"    -- presto::wilcoxauc(); orders of magnitude faster
#                            on sparse matrices. Optional dependency from
#                            the immunogenomics/presto GitHub repo; the
#                            dispatcher falls back to "wilcox_r" when
#                            presto is missing.
#
#   backend = "auto"      -- picks "presto" if installed AND test == "wilcox",
#                            else "wilcox_r".
#
#   backend = "pseudobulk_naive"   -- aggregate counts by (sample, group)
#                                     via `aggregate_pseudobulk()` then
#                                     run per-gene t-tests on log2(CPM+1).
#                                     Pure-R, always available.
#   backend = "pseudobulk_edger"   -- aggregation + edgeR::glmQLFit /
#                                     glmQLFTest. Gated by
#                                     `require_optional("edgeR")`.
#   backend = "pseudobulk_deseq2"  -- aggregation + DESeq2::DESeq.
#                                     Gated by `require_optional("DESeq2")`.
#
# Aggregation lives in R/pseudobulk.R. All pseudobulk backends require
# a `sample_by` metadata column (biological replicate id) and produce
# the canonical DE schema via `.pseudobulk_to_de_schema()`.
#
# The result frame from every backend matches `empty_de_results()` byte-for-byte
# so downstream filtering/plotting/stamping logic does not change.
#
# TODO (future plug-in points):
#   - covariate-aware models (MAST `~ group + nFeature_RNA`).
#   - per-cluster pseudobulk fan-out (Crowell muscat-style).
# ============================================================================

# ---- Result schema -------------------------------------------------------

#' Empty DE result frame with the canonical schema.
empty_de_results <- function(group_1 = NA_character_, group_2 = NA_character_) {
  data.frame(
    gene       = character(),
    group_1    = character(),
    group_2    = character(),
    avg_log2FC = numeric(),
    pct.1      = numeric(),
    pct.2      = numeric(),
    p_val      = numeric(),
    p_val_adj  = numeric(),
    stringsAsFactors = FALSE
  )
}

# ---- Backend registry ----------------------------------------------------

DE_BACKENDS <- c(
  "Auto-select"                    = "auto",
  "Wilcoxon (pure R)"              = "wilcox_r",
  "presto::wilcoxauc"              = "presto",
  "Pseudobulk (naive t-test)"      = "pseudobulk_naive",
  "Pseudobulk (edgeR)"             = "pseudobulk_edger",
  "Pseudobulk (DESeq2)"            = "pseudobulk_deseq2"
)

PSEUDOBULK_BACKENDS <- c("pseudobulk_naive",
                         "pseudobulk_edger",
                         "pseudobulk_deseq2")

#' List the DE backends available in the current R session.
#'
#' Each entry has `id`, `label`, `available` (logical), and `kind`
#' ("cell" or "pseudobulk"). `"auto"` always resolves to one of the real
#' cell-level backends, so it is reported as available. UI controls
#' can use the `kind` field to drive a "pseudobulk needs sample_by"
#' conditional input.
de_available_backends <- function() {
  list(
    list(id = "auto",              label = "Auto-select",
         kind = "cell",            available = TRUE),
    list(id = "wilcox_r",          label = "Wilcoxon (pure R)",
         kind = "cell",            available = TRUE),
    list(id = "presto",            label = "presto::wilcoxauc",
         kind = "cell",            available = has_optional("presto")),
    list(id = "pseudobulk_naive",  label = "Pseudobulk (naive t-test)",
         kind = "pseudobulk",      available = TRUE),
    list(id = "pseudobulk_edger",  label = "Pseudobulk (edgeR)",
         kind = "pseudobulk",      available = has_optional("edgeR")),
    list(id = "pseudobulk_deseq2", label = "Pseudobulk (DESeq2)",
         kind = "pseudobulk",      available = has_optional("DESeq2"))
  )
}

#' Resolve "auto" -> a concrete backend id. Exposed for testability.
.de_resolve_backend <- function(backend, test) {
  if (!backend %in% DE_BACKENDS) {
    stop("Unknown DE backend '", backend, "'. Have: ",
         paste(unname(DE_BACKENDS), collapse = ", "), call. = FALSE)
  }
  if (backend != "auto") return(backend)
  if (identical(test, "wilcox") && has_optional("presto")) "presto" else "wilcox_r"
}

#' Is this backend a pseudobulk backend?
.is_pseudobulk_backend <- function(backend) {
  backend %in% PSEUDOBULK_BACKENDS
}

# ---- Public API ----------------------------------------------------------

#' Run differential expression between two groups of cells.
#'
#' @param dataset         a dataset list (`dataset_schema()`)
#' @param grouping_field  metadata column to split cells on
#' @param group_1,group_2 two distinct values of `grouping_field`
#' @param assay           assay name (passed through; ignored by mock-backed loaders)
#' @param layer           layer / slot name (passed through to the expression
#'                        backend; defaults to the backend's default layer).
#'                        Pseudobulk backends default to `"counts"` when
#'                        present and warn if forced onto log-normalised data.
#'                        The gene universe is taken from the resolved
#'                        layer -- so passing an explicit `layer` will
#'                        also restrict DE to the genes that exist in
#'                        that layer (different layers on the same
#'                        backend can carry different gene sets).
#' @param min_pct         pre-test filter: drop genes with max(pct.1, pct.2) < min_pct
#' @param test            "wilcox" (default) or "t"
#' @param backend         "auto" (default), "wilcox_r", "presto",
#'                        "pseudobulk_naive", "pseudobulk_edger",
#'                        "pseudobulk_deseq2"
#' @param sample_by       (pseudobulk only) metadata column carrying biological
#'                        replicate ids. Required for pseudobulk backends.
#' @param min_cells_per_sample   (pseudobulk only) drop pseudobulk samples
#'                               assembled from fewer than this many cells.
#' @param min_samples_per_group  (pseudobulk only) each group must have at
#'                               least this many surviving pseudobulk samples.
#'
#' @return data.frame matching `empty_de_results()`, sorted by p_val_adj then
#'   |avg_log2FC|. Throws an informative error on invalid input.
compute_de <- function(dataset, grouping_field, group_1, group_2,
                       assay = NULL, layer = NULL,
                       min_pct = 0.1,
                       test    = c("wilcox", "t"),
                       backend = c("auto", "wilcox_r", "presto",
                                   "pseudobulk_naive",
                                   "pseudobulk_edger",
                                   "pseudobulk_deseq2"),
                       sample_by = NULL,
                       min_cells_per_sample  = 10L,
                       min_samples_per_group = 2L) {
  test    <- match.arg(test)
  backend <- match.arg(backend)

  if (is.null(dataset))            stop("No dataset provided.",          call. = FALSE)
  if (is.null(group_1) || is.null(group_2) || group_1 == group_2)
    stop("group_1 and group_2 must be distinct.", call. = FALSE)

  meta <- get_metadata(dataset, grouping_field)
  if (is.null(meta))
    stop(sprintf("Metadata field '%s' is not available.", grouping_field), call. = FALSE)
  in1 <- which(as.character(meta) == as.character(group_1))
  in2 <- which(as.character(meta) == as.character(group_2))
  if (length(in1) < 2L || length(in2) < 2L)
    stop("Each group needs at least 2 cells.", call. = FALSE)

  # NB: must pass `layer` so the gene universe respects the resolved
  # layer. With a NULL layer this is equivalent to the previous behaviour
  # (default-layer genes); with an explicit layer it now skips genes that
  # only exist in *other* layers. Previously this used the default-layer
  # genes regardless and silently produced empty rows for genes not
  # present in the requested layer.
  genes <- available_genes(dataset, layer = layer)
  if (length(genes) == 0L)
    stop(sprintf("No genes available in layer '%s'.",
                 layer %||% "<default>"), call. = FALSE)

  resolved <- .de_resolve_backend(backend, test)
  if (resolved == "presto" && test != "wilcox") {
    stop("backend = 'presto' only supports test = 'wilcox'.", call. = FALSE)
  }

  if (.is_pseudobulk_backend(resolved)) {
    validate_pseudobulk_inputs(
      dataset, grouping_field, group_1, group_2, sample_by,
      min_cells_per_sample  = min_cells_per_sample,
      min_samples_per_group = min_samples_per_group)
    pb <- aggregate_pseudobulk(
      dataset, grouping_field, group_1, group_2, sample_by,
      layer = layer, agg = "sum",
      min_cells_per_sample = min_cells_per_sample)
    if (is.null(pb)) return(empty_de_results(group_1, group_2))
    pcts <- .pseudobulk_pct(dataset, in1, in2,
                            layer = pb$layer_used %||% "counts")
    out <- switch(resolved,
      "pseudobulk_naive"  = .de_run_pseudobulk_naive(pb,  pcts, group_1, group_2,
                                                     min_pct = min_pct),
      "pseudobulk_edger"  = .de_run_pseudobulk_edger(pb,  pcts, group_1, group_2,
                                                     min_pct = min_pct),
      "pseudobulk_deseq2" = .de_run_pseudobulk_deseq2(pb, pcts, group_1, group_2,
                                                     min_pct = min_pct)
    )
  } else {
    out <- switch(resolved,
      "wilcox_r" = .de_run_wilcox_r(dataset, genes, in1, in2,
                                    group_1, group_2, layer = layer,
                                    min_pct = min_pct, test = test),
      "presto"   = .de_run_presto(dataset, genes, in1, in2,
                                  group_1, group_2, layer = layer,
                                  min_pct = min_pct)
    )
  }

  if (nrow(out) == 0L) return(empty_de_results(group_1, group_2))
  out <- out[order(out$p_val_adj, -abs(out$avg_log2FC)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# ---- Backend: pure R -----------------------------------------------------

.de_run_wilcox_r <- function(dataset, genes, in1, in2,
                             group_1, group_2, layer = NULL,
                             min_pct = 0.1, test = "wilcox") {
  pval_fn <- switch(test,
    "wilcox" = function(x, y) tryCatch(
      suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value),
      error = function(e) NA_real_),
    "t"      = function(x, y) tryCatch(stats::t.test(x, y)$p.value,
                                       error = function(e) NA_real_)
  )

  # FUN.VALUE template == per-gene vector length, taken from the first
  # gene (every gene aligns to the same n_cells via the backend).
  n_cells_total <- length(get_gene_expression(dataset, genes[1], layer = layer))
  expr_mat <- vapply(genes,
                     function(g) get_gene_expression(dataset, g, layer = layer),
                     FUN.VALUE = numeric(n_cells_total))

  rows <- list()
  for (j in seq_along(genes)) {
    v1 <- expr_mat[in1, j]; v2 <- expr_mat[in2, j]
    pct1 <- mean(v1 > 0);   pct2 <- mean(v2 > 0)
    if (max(pct1, pct2) < min_pct) next
    m1 <- mean(v1);         m2 <- mean(v2)
    rows[[length(rows) + 1L]] <- data.frame(
      gene       = genes[j],
      group_1    = as.character(group_1),
      group_2    = as.character(group_2),
      avg_log2FC = log2(m1 + 1) - log2(m2 + 1),
      pct.1      = pct1,
      pct.2      = pct2,
      p_val      = pval_fn(v1, v2),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) return(empty_de_results(group_1, group_2))
  df <- do.call(rbind, rows)
  df$p_val_adj <- stats::p.adjust(df$p_val, method = "BH")
  df
}

# ---- Backend: presto -----------------------------------------------------

.de_run_presto <- function(dataset, genes, in1, in2,
                           group_1, group_2, layer = NULL,
                           min_pct = 0.1) {
  # presto is GitHub-only (immunogenomics/presto); require_optional()
  # generates the correct `remotes::install_github(...)` hint.
  require_optional("presto",
                   feature = "compute_de(backend = 'presto')",
                   source  = "GitHub",
                   repo    = c(presto = "immunogenomics/presto"))

  X <- .de_build_genes_x_cells(dataset, genes,
                               cell_idx = c(in1, in2), layer = layer)
  y <- c(rep(as.character(group_1), length(in1)),
         rep(as.character(group_2), length(in2)))
  stopifnot(ncol(X) == length(y))

  presto_df <- presto::wilcoxauc(X, y)
  .presto_to_de_schema(presto_df, group_1 = group_1, group_2 = group_2,
                       min_pct = min_pct)
}

#' Build a dense `genes x cells` matrix for one backend call.
#'
#' Materialises just the union of cells across the two DE groups (so this
#' stays small even on 100k-cell datasets) and returns a `genes x cells`
#' matrix with `rownames(X) = genes` and column order following
#' `cell_idx`. The orientation matches `presto::wilcoxauc()`'s
#' feature-major contract; the same helper can be reused by any future
#' backend that wants the same shape.
#'
#' Extracted from `.de_run_presto()` so the orientation/labelling contract
#' is unit-testable without `presto` installed.
#'
#' @param cell_idx  integer vector of positional cell indices (1..n_cells)
#'                  matching the order in `dataset$cell_data$cell`. NOT
#'                  a vector of cell IDs -- the underlying backend
#'                  vectors are unnamed.
.de_build_genes_x_cells <- function(dataset, genes, cell_idx,
                                    layer = NULL) {
  if (!length(genes))    stop(".de_build_genes_x_cells: no genes.", call. = FALSE)
  if (!length(cell_idx)) stop(".de_build_genes_x_cells: no cells.", call. = FALSE)
  if (!is.numeric(cell_idx))
    stop(".de_build_genes_x_cells: `cell_idx` must be integer positions, not cell IDs.",
         call. = FALSE)
  n_used <- length(cell_idx)
  # vapply yields a [cells x genes] matrix with colnames = genes; we
  # then transpose so the result is [genes x cells] with rownames = genes.
  X <- vapply(genes,
              function(g) get_gene_expression(dataset, g, layer = layer)[cell_idx],
              FUN.VALUE = numeric(n_used))
  X <- t(X)
  rownames(X) <- as.character(genes)
  X
}

#' Convert a `presto::wilcoxauc()` output frame into the canonical DE schema.
#'
#' Pure / dependency-free: takes any data.frame that quacks like a presto
#' result and returns a `empty_de_results()`-shaped frame. Factored out so
#' the schema-mapping logic can be regression-tested without `presto`
#' installed.
.presto_to_de_schema <- function(presto_df, group_1, group_2, min_pct = 0) {
  if (is.null(presto_df) || nrow(presto_df) == 0L) {
    return(empty_de_results(group_1, group_2))
  }
  # presto can name columns either `pval`/`padj` or `pvalue`/`padj`.
  pcol  <- intersect(c("pval",  "pvalue"), names(presto_df))[1]
  pacol <- intersect(c("padj",  "p_val_adj"), names(presto_df))[1]
  if (is.na(pcol))  pcol  <- "pval"
  if (is.na(pacol)) pacol <- "padj"

  g1 <- presto_df[as.character(presto_df$group) == as.character(group_1), , drop = FALSE]
  if (nrow(g1) == 0L) return(empty_de_results(group_1, group_2))

  out <- data.frame(
    gene       = as.character(g1$feature),
    group_1    = as.character(group_1),
    group_2    = as.character(group_2),
    avg_log2FC = as.numeric(g1$logFC),
    pct.1      = as.numeric(g1$pct_in)  / 100,  # presto reports percent
    pct.2      = as.numeric(g1$pct_out) / 100,
    p_val      = as.numeric(g1[[pcol]]),
    p_val_adj  = as.numeric(g1[[pacol]]),
    stringsAsFactors = FALSE
  )
  if (!is.null(min_pct) && is.finite(min_pct) && min_pct > 0) {
    out <- out[pmax(out$pct.1, out$pct.2) >= min_pct, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

# ---- Backend: pseudobulk_naive -------------------------------------------
# Aggregation already happened in aggregate_pseudobulk(). Here we
# library-size normalise (CPM), log2-transform, and run a per-gene
# two-sample t-test across pseudobulk samples. Pure-R, dependency-free.
# Reasonable smoke-grade DE; for production users should pick edgeR or
# DESeq2.

.de_run_pseudobulk_naive <- function(pb, pcts, group_1, group_2,
                                     min_pct = 0) {
  M  <- pb$matrix
  sm <- pb$sample_metadata
  is1 <- sm$group == as.character(group_1)
  is2 <- sm$group == as.character(group_2)
  if (sum(is1) < 2L || sum(is2) < 2L) {
    return(empty_de_results(group_1, group_2))
  }
  Mn <- pseudobulk_cpm_log2(M)
  rows <- vector("list", nrow(Mn))
  for (i in seq_len(nrow(Mn))) {
    v1 <- Mn[i, is1]; v2 <- Mn[i, is2]
    pv <- tryCatch(stats::t.test(v1, v2)$p.value,
                   error = function(e) NA_real_)
    rows[[i]] <- data.frame(
      gene       = rownames(Mn)[i],
      avg_log2FC = mean(v1) - mean(v2),
      pct.1      = pcts$pct.1[[rownames(Mn)[i]]] %||% NA_real_,
      pct.2      = pcts$pct.2[[rownames(Mn)[i]]] %||% NA_real_,
      p_val      = pv,
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, rows)
  .pseudobulk_to_de_schema(df, group_1 = group_1, group_2 = group_2,
                           min_pct = min_pct)
}

# ---- Backend: pseudobulk_edger -------------------------------------------
# Counts-matrix + group vector -> edgeR::glmQLFit / glmQLFTest. Standard
# pseudobulk DE. Gated by `require_optional("edgeR")`; the function body
# below is only reached when edgeR is installed.

.de_run_pseudobulk_edger <- function(pb, pcts, group_1, group_2,
                                     min_pct = 0) {
  require_optional("edgeR",
                   feature = "compute_de(backend = 'pseudobulk_edger')",
                   source  = "Bioconductor")
  M  <- pb$matrix
  sm <- pb$sample_metadata
  groups <- factor(sm$group, levels = c(as.character(group_1),
                                        as.character(group_2)))
  if (any(table(groups) < 2L))
    return(empty_de_results(group_1, group_2))

  dge_fn   <- get("DGEList",       envir = asNamespace("edgeR"))
  norm_fn  <- get("calcNormFactors", envir = asNamespace("edgeR"))
  est_fn   <- get("estimateDisp",  envir = asNamespace("edgeR"))
  fit_fn   <- get("glmQLFit",      envir = asNamespace("edgeR"))
  test_fn  <- get("glmQLFTest",    envir = asNamespace("edgeR"))
  top_fn   <- get("topTags",       envir = asNamespace("edgeR"))

  # Design: group_2 is the reference (first level), so a positive logFC
  # means "higher in group_1" -- consistent with the cell-level kernel.
  design <- stats::model.matrix(~ stats::relevel(groups,
                                                 ref = as.character(group_2)))

  y <- dge_fn(counts = M, group = groups)
  y <- norm_fn(y)
  y <- est_fn(y, design = design)
  fit <- fit_fn(y, design = design)
  qlf <- test_fn(fit, coef = 2L)
  tt  <- top_fn(qlf, n = Inf, sort.by = "none")$table

  df <- data.frame(
    gene       = rownames(tt),
    avg_log2FC = as.numeric(tt$logFC),
    pct.1      = pcts$pct.1[rownames(tt)],
    pct.2      = pcts$pct.2[rownames(tt)],
    p_val      = as.numeric(tt$PValue),
    stringsAsFactors = FALSE
  )
  .pseudobulk_to_de_schema(df, group_1 = group_1, group_2 = group_2,
                           min_pct = min_pct)
}

# ---- Backend: pseudobulk_deseq2 ------------------------------------------
# Counts-matrix + colData -> DESeq2::DESeqDataSetFromMatrix -> DESeq.
# Gated by `require_optional("DESeq2")`. Heavier than edgeR but ships
# moderated LFC out of the box; we do NOT apply lfcShrink here -- it's
# a per-PR call whether to add it.

.de_run_pseudobulk_deseq2 <- function(pb, pcts, group_1, group_2,
                                      min_pct = 0) {
  require_optional("DESeq2",
                   feature = "compute_de(backend = 'pseudobulk_deseq2')",
                   source  = "Bioconductor")
  M  <- pb$matrix
  sm <- pb$sample_metadata
  groups <- factor(sm$group, levels = c(as.character(group_2),  # ref first
                                        as.character(group_1)))
  if (any(table(groups) < 2L))
    return(empty_de_results(group_1, group_2))

  ds_from_mat <- get("DESeqDataSetFromMatrix", envir = asNamespace("DESeq2"))
  deseq_run   <- get("DESeq",                  envir = asNamespace("DESeq2"))
  results_fn  <- get("results",                envir = asNamespace("DESeq2"))

  coldata <- data.frame(group = groups,
                        row.names = colnames(M),
                        stringsAsFactors = FALSE)
  dds <- ds_from_mat(countData = round(M), colData = coldata,
                     design = ~ group)
  dds <- deseq_run(dds, quiet = TRUE)
  res <- as.data.frame(results_fn(dds))

  df <- data.frame(
    gene       = rownames(res),
    avg_log2FC = as.numeric(res$log2FoldChange),
    pct.1      = pcts$pct.1[rownames(res)],
    pct.2      = pcts$pct.2[rownames(res)],
    p_val      = as.numeric(res$pvalue),
    stringsAsFactors = FALSE
  )
  .pseudobulk_to_de_schema(df, group_1 = group_1, group_2 = group_2,
                           min_pct = min_pct)
}

# ---- Display-side helpers (unchanged) -----------------------------------

#' Apply the live table filters (gene search + log2FC + padj) to a DE frame.
filter_de_results <- function(de_df, gene_search = "",
                              min_abs_log2fc = 0, max_padj = 1) {
  if (is.null(de_df) || nrow(de_df) == 0L) return(de_df)
  keep <- rep(TRUE, nrow(de_df))
  if (nzchar(gene_search)) {
    keep <- keep & grepl(gene_search, de_df$gene, ignore.case = TRUE, fixed = FALSE)
  }
  if (!is.null(min_abs_log2fc) && is.finite(min_abs_log2fc) && min_abs_log2fc > 0) {
    keep <- keep & abs(de_df$avg_log2FC) >= min_abs_log2fc
  }
  if (!is.null(max_padj) && is.finite(max_padj) && max_padj < 1) {
    keep <- keep & !is.na(de_df$p_val_adj) & de_df$p_val_adj <= max_padj
  }
  de_df[keep, , drop = FALSE]
}

#' Sort a DE frame by one of the standard columns.
sort_de_results <- function(de_df,
                            sort_by = c("p_val_adj", "avg_log2FC", "pct.1", "pct.2", "gene"),
                            descending = TRUE) {
  sort_by <- match.arg(sort_by)
  if (is.null(de_df) || nrow(de_df) == 0L) return(de_df)
  v <- de_df[[sort_by]]
  ord <- if (is.character(v)) order(v, decreasing = descending)
         else                  order(v, decreasing = descending, na.last = TRUE)
  de_df[ord, , drop = FALSE]
}
