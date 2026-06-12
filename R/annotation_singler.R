# ============================================================================
# Annotation engine: SingleR
# ----------------------------------------------------------------------------
# Reference-based per-cell label assignment using `SingleR::SingleR()` and
# a `celldex` reference. Pluggable through `ANNOTATION_ENGINES()` -- there
# is no UI specialisation for SingleR; the existing Annotation module
# renders engine parameters from the spec.
#
# Layout:
#   * `.run_singler_annotation()`      -- the engine `run_fn`. Heavy-deps
#                                          gated through `require_optional`.
#                                          Builds the test matrix, fetches
#                                          a `celldex` reference, runs
#                                          SingleR in per-cell or per-
#                                          cluster mode, and routes the
#                                          output through the pure
#                                          converter.
#   * `.fetch_celldex_reference()`     -- thin wrapper that maps "hpca" /
#                                          "blueprint_encode" / "monaco_immune"
#                                          to the right celldex function,
#                                          with a cache shared across runs.
#   * `.singler_to_engine_output()`    -- pure schema converter (per-cell).
#                                          No SingleR dependency; testable
#                                          with hand-built input frames.
#   * `.singler_cluster_to_cells()`    -- expands a per-cluster SingleR
#                                          DataFrame into per-cell labels
#                                          + scores via the dataset's
#                                          cluster_field.
#
# Caching: celldex downloads are non-trivial (10-100MB Bioconductor cache).
# We memoise per-process via `.singler_ref_cache` so repeated runs against
# the same reference don't re-hit the network. The cache is module-local;
# it is invalidated by reloading the source file.
# ============================================================================

# Process-local reference cache. NULL until first successful fetch.
.singler_ref_cache <- new.env(parent = emptyenv())

#' Map a friendly reference id to a `celldex` fetcher.
#' Returns a `SummarizedExperiment` with `assay(., "logcounts")` and the
#' relevant label column (`label.main` or `label.fine`).
.fetch_celldex_reference <- function(reference) {
  require_optional("celldex",
                   feature = "SingleR engine reference fetch",
                   source  = "Bioconductor")
  if (!is.null(.singler_ref_cache[[reference]])) {
    return(.singler_ref_cache[[reference]])
  }
  fn <- switch(reference,
    "hpca"             = celldex::HumanPrimaryCellAtlasData,
    "blueprint_encode" = celldex::BlueprintEncodeData,
    "monaco_immune"    = celldex::MonacoImmuneData,
    NULL
  )
  if (is.null(fn)) {
    stop("SingleR engine: unknown reference '", reference,
         "'. Available: hpca, blueprint_encode, monaco_immune.",
         call. = FALSE)
  }
  ref <- fn()
  .singler_ref_cache[[reference]] <- ref
  ref
}

# ---- Engine run function -------------------------------------------------

.run_singler_annotation <- function(dataset, params, state) {
  require_optional(
    c("SingleR", "celldex", "SingleCellExperiment", "SummarizedExperiment"),
    feature = "SingleR annotation engine",
    source  = "Bioconductor")

  reference     <- params$reference     %||% "hpca"
  labels_field  <- params$labels        %||% "main"
  cluster_field <- params$cluster_field
  min_delta     <- params$min_delta     %||% 0.0

  cells <- dataset$cell_data$cell

  # 1. Build the test matrix via the backend abstraction so sparse / h5ad
  #    backends pass through their native matrix-like (no premature
  #    materialisation). SingleR is happy with dgCMatrix or dense.
  be <- as_expression_backend(dataset$expression)
  test_mat <- backend_as_matrix(be, layer = backend_default_layer(be))
  if (is.null(rownames(test_mat))) {
    rownames(test_mat) <- backend_genes(be)
  }
  if (is.null(colnames(test_mat))) {
    colnames(test_mat) <- cells
  }

  # 2. Reference. We pull logcounts + the requested label granularity.
  ref <- .fetch_celldex_reference(reference)
  ref_assay_name <- if ("logcounts" %in% SummarizedExperiment::assayNames(ref))
                      "logcounts" else SummarizedExperiment::assayNames(ref)[1]
  ref_mat <- SummarizedExperiment::assay(ref, ref_assay_name)
  ref_labels_col <- paste0("label.", labels_field)
  if (!ref_labels_col %in% names(SingleCellExperiment::colData(ref))) {
    stop("SingleR engine: reference '", reference,
         "' has no label column '", ref_labels_col, "'.", call. = FALSE)
  }
  ref_labels <- SingleCellExperiment::colData(ref)[[ref_labels_col]]

  # 3. Restrict to shared genes. SingleR can also do this internally but
  #    the friendlier error message lives here.
  shared <- intersect(rownames(test_mat), rownames(ref_mat))
  if (length(shared) < 50L) {
    stop(sprintf(
      "SingleR engine: only %d shared genes between test and '%s'. Check ",
      length(shared), reference),
      "gene symbol orthology (the dataset's gene names must match the ",
      "reference's species and naming convention).", call. = FALSE)
  }
  test_mat <- test_mat[shared, , drop = FALSE]
  ref_mat  <- ref_mat[shared, , drop = FALSE]

  reference_source <- sprintf("celldex/%s:%s", reference, labels_field)

  # 4. Per-cluster vs per-cell mode.
  if (!is.null(cluster_field) && nzchar(cluster_field)) {
    cluster_vec <- as.character(get_metadata(dataset, cluster_field))
    if (is.null(cluster_vec)) {
      stop(sprintf("SingleR engine: cluster_field '%s' not in dataset.",
                   cluster_field), call. = FALSE)
    }
    pred <- SingleR::SingleR(test     = test_mat,
                             ref      = ref_mat,
                             labels   = ref_labels,
                             clusters = cluster_vec)
    .singler_cluster_to_cells(
      pred, dataset = dataset, cluster_vec = cluster_vec,
      cluster_field_used = cluster_field,
      reference_source   = reference_source,
      min_delta          = min_delta)
  } else {
    pred <- SingleR::SingleR(test = test_mat, ref = ref_mat, labels = ref_labels)
    .singler_to_engine_output(
      pred, cells = cells,
      cluster_vec        = NULL,
      cluster_field_used = NA_character_,
      reference_source   = reference_source,
      min_delta          = min_delta)
  }
}

# ---- Pure schema converters ---------------------------------------------

#' Convert a per-cell SingleR `DataFrame` into the engine output schema.
#'
#' Accepts any data.frame-ish object with at minimum a `labels` column,
#' and optionally `pruned.labels` and `delta.next`. The `scores` matrix
#' column (one row per cell, one column per ref label) is consumed if
#' present to fill `alt_labels` with the top-3 candidates per cell.
#'
#' Factored out so the schema mapping has its own regression coverage
#' without requiring SingleR/celldex to be installed.
#'
#' @param pred              SingleR-style per-cell output (n_cells rows)
#' @param cells             character(n_cells); dataset cell ids in order
#' @param cluster_vec       optional character(n_cells) for cluster_summary
#' @param cluster_field_used character(1) name of the cluster field (or NA)
#' @param reference_source  character(1) e.g. "celldex/hpca:main"
#' @param min_delta         numeric(1) below this delta.next -> "Unknown"
.singler_to_engine_output <- function(pred, cells,
                                      cluster_vec = NULL,
                                      cluster_field_used = NA_character_,
                                      reference_source = NA_character_,
                                      min_delta = 0.0) {
  if (is.null(pred) || nrow(pred) == 0L) {
    return(list(
      cell                   = cells,
      cell_labels            = rep(NA_character_, length(cells)),
      cell_scores            = rep(0.0, length(cells)),
      alt_labels             = NULL,
      cluster_summary        = NULL,
      cluster_field_used     = cluster_field_used %||% NA_character_,
      n_clusters_at_creation = NA_integer_,
      reference_source       = reference_source,
      warnings               = "SingleR returned an empty prediction frame."
    ))
  }
  if (nrow(pred) != length(cells)) {
    stop(sprintf(
      ".singler_to_engine_output: pred has %d rows but cells has %d.",
      nrow(pred), length(cells)), call. = FALSE)
  }

  pruned <- if ("pruned.labels" %in% names(pred)) pred$pruned.labels
            else pred$labels
  delta  <- if ("delta.next" %in% names(pred))    pred$delta.next
            else rep(NA_real_, length(cells))

  cell_labels <- as.character(pruned %||% pred$labels)
  cell_labels[is.na(cell_labels)] <- "Unknown"
  if (is.finite(min_delta) && min_delta > 0) {
    cell_labels[!is.na(delta) & delta < min_delta] <- "Unknown"
  }
  cell_scores <- as.numeric(delta)
  cell_scores[is.na(cell_scores) | cell_labels == "Unknown"] <- 0.0

  # alt_labels: top-3 candidates per cell, if SingleR returned a scores matrix.
  alt_labels <- NULL
  if ("scores" %in% names(pred)) {
    scores <- as.matrix(pred$scores)
    if (nrow(scores) == length(cells) && ncol(scores) > 0L) {
      ref_lab <- colnames(scores) %||% paste0("ref_", seq_len(ncol(scores)))
      rows <- list()
      k <- min(3L, ncol(scores))
      for (i in seq_len(nrow(scores))) {
        ord <- order(scores[i, ], decreasing = TRUE)[seq_len(k)]
        for (r in seq_along(ord)) {
          rows[[length(rows) + 1L]] <- data.frame(
            cell  = cells[i], rank = r,
            label = ref_lab[ord[r]],
            score = as.numeric(scores[i, ord[r]]),
            stringsAsFactors = FALSE)
        }
      }
      if (length(rows)) alt_labels <- do.call(rbind, rows)
    }
  }

  cluster_summary <- NULL
  n_clusters <- NA_integer_
  if (!is.null(cluster_vec) && length(cluster_vec) == length(cells)) {
    cluster_vec <- as.character(cluster_vec)
    cluster_ids <- sort(unique(cluster_vec))
    n_clusters  <- length(cluster_ids)
    cluster_summary <- .summarise_per_cluster(cell_labels, cell_scores, cluster_vec, cluster_ids)
  }

  list(
    cell                   = cells,
    cell_labels            = cell_labels,
    cell_scores            = cell_scores,
    alt_labels             = alt_labels,
    cluster_summary        = cluster_summary,
    cluster_field_used     = cluster_field_used %||% NA_character_,
    n_clusters_at_creation = n_clusters,
    reference_source       = reference_source,
    warnings               = character()
  )
}

#' Convert a per-cluster SingleR DataFrame into a per-cell engine output.
#'
#' SingleR returns one row per *cluster* when called with `clusters =`.
#' This helper expands those cluster-level predictions back to every cell
#' in the cluster, then defers the per-cell -> schema mapping to
#' `.singler_to_engine_output()`.
.singler_cluster_to_cells <- function(pred, dataset, cluster_vec,
                                      cluster_field_used,
                                      reference_source,
                                      min_delta = 0.0) {
  cells <- dataset$cell_data$cell
  cluster_vec <- as.character(cluster_vec)

  # Expand to per-cell by row-indexing the cluster pred frame with each
  # cell's cluster id. SingleR's `clusters =` mode stamps the cluster ids
  # onto rownames(pred); we require at least one cluster_vec value to
  # match (otherwise the rownames don't carry cluster ids at all, which
  # almost certainly means the caller built `pred` manually with default
  # numeric rownames).
  pred_rn <- rownames(pred)
  if (is.null(pred_rn) || !any(cluster_vec %in% pred_rn)) {
    stop(".singler_cluster_to_cells: per-cluster pred has no cluster-id ",
         "rownames; cannot map clusters back to cells. (SingleR populates ",
         "rownames(pred) when called with `clusters =`.)", call. = FALSE)
  }
  idx <- match(cluster_vec, pred_rn)
  per_cell <- pred[idx, , drop = FALSE]
  rownames(per_cell) <- cells

  .singler_to_engine_output(
    per_cell, cells = cells,
    cluster_vec        = cluster_vec,
    cluster_field_used = cluster_field_used,
    reference_source   = reference_source,
    min_delta          = min_delta)
}

# Internal: build a {cluster, top_label, top_score, n_cells} summary frame
# from per-cell label + score vectors. Used by both per-cell and per-cluster
# SingleR paths so the cluster_summary shape stays consistent with the
# manual / marker_score engines.
.summarise_per_cluster <- function(cell_labels, cell_scores,
                                   cluster_vec, cluster_ids) {
  top_label <- character(length(cluster_ids))
  top_score <- numeric(length(cluster_ids))
  n_cells   <- integer(length(cluster_ids))
  for (i in seq_along(cluster_ids)) {
    in_cl <- cluster_vec == cluster_ids[i]
    n_cells[i] <- sum(in_cl)
    lab <- cell_labels[in_cl]
    if (!length(lab) || all(is.na(lab) | lab == "Unknown")) {
      top_label[i] <- "Unknown"; top_score[i] <- 0.0; next
    }
    tab <- sort(table(lab[!is.na(lab) & lab != "Unknown"]), decreasing = TRUE)
    if (!length(tab)) {
      top_label[i] <- "Unknown"; top_score[i] <- 0.0; next
    }
    top_label[i] <- names(tab)[1]
    top_score[i] <- mean(cell_scores[in_cl & cell_labels == top_label[i]], na.rm = TRUE)
  }
  data.frame(
    cluster   = cluster_ids,
    top_label = top_label,
    top_score = top_score,
    n_cells   = n_cells,
    stringsAsFactors = FALSE
  )
}
