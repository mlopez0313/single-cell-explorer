# ============================================================================
# Expression backend: lazy HDF5 / AnnData
# ----------------------------------------------------------------------------
# `expression_backend_h5ad()` reads expression on demand from an `.h5ad`
# file (HDF5). The other two backends materialise the full matrix in
# memory; this one keeps a file path and only pulls the gene(s) modules
# actually ask for. Crucial for 100k+ cell datasets where a dense
# n_genes x n_cells materialisation is several GB.
#
# Storage formats handled (per AnnData v0.8+ spec). The conceptual
# matrix is always `n_obs x n_var` (cells x genes); CSR and CSC differ
# only in which axis they walk first:
#
#   /X (or /layers/<name>)
#     * Sparse `csr_matrix`  -- row-major over cells (the typical layout
#         scanpy writes; per-gene reads walk every row)
#         /X
#           encoding-type = "csr_matrix"
#           shape = [n_obs, n_var]
#           /X/data    (length nnz)
#           /X/indices (length nnz; entries are 0-based gene indices,
#                       i.e. column indices into the conceptual matrix)
#           /X/indptr  (length n_obs + 1)
#     * Sparse `csc_matrix`  -- column-major over genes (per-gene reads
#         are a single contiguous slice)
#         /X
#           encoding-type = "csc_matrix"
#           shape = [n_obs, n_var]
#           /X/data    (length nnz)
#           /X/indices (length nnz; entries are 0-based cell indices,
#                       i.e. row indices into the conceptual matrix)
#           /X/indptr  (length n_var + 1)
#     * Dense 2D dataset, HDF5 shape [n_obs, n_var] (rare for scRNA-seq;
#         materialises in memory either via h5read default dim order or
#         the rhdf5 `native = TRUE` transpose, both handled below).
#
# The backend caches indptr / indices / data lazily on first access in
# a per-layer env. Per-gene read cost:
#   * CSC : O(nnz_per_gene) -- one contiguous slice.
#   * CSR : O(n_obs + nnz_per_row_segment_searches) -- we have to walk
#           every cell row's segment to find the requested gene index.
#           Still far cheaper than materialising the full matrix, but
#           noticeably slower than CSC on large datasets.
#
# Layer shape (passed to the constructor):
#
#   layers = list(
#     data   = list(path = "/X",              encoding = "csr_matrix",
#                   shape = c(n_obs, n_var), genes = c("CD3D", ...)),
#     counts = list(path = "/layers/counts",  encoding = "csc_matrix",
#                   shape = c(n_obs, n_var), genes = c(...))
#   )
#
# Notes:
#   - Gene names live in /var/_index in the file but each layer carries
#     its own copy so callers don't need to thread them through.
#   - The file handle is not held open between calls; each rhdf5 call
#     opens + closes the file. For sparse layers that hit cost only
#     once (the first call populates the in-process triple cache); the
#     dense per-gene reader pays it on every call. This keeps the
#     backend cheap to serialise (a path + small metadata, plus an
#     in-process cache that's not persisted).
# ============================================================================

H5AD_SUPPORTED_ENCODINGS <- c("csr_matrix", "csc_matrix", "dense")

#' Build a lazy h5ad expression backend.
#'
#' Validates the file is readable and that each layer's encoding /
#' shape is recognised, but does not eagerly load any matrix data.
#'
#' @param path           character(1) path to a `.h5ad` file
#' @param layers         named list of layer specs (see file docstring)
#' @param n_cells        integer(1) total cells (n_obs)
#' @param default_layer  layer name used when callers omit `layer`
expression_backend_h5ad <- function(path, layers,
                                    n_cells       = NULL,
                                    default_layer = "data") {
  if (!is.character(path) || length(path) != 1L || !file.exists(path))
    stop("expression_backend_h5ad: `path` must point to an existing .h5ad file.",
         call. = FALSE)
  if (!is.list(layers) || length(layers) == 0L ||
      is.null(names(layers)) || any(!nzchar(names(layers))))
    stop("expression_backend_h5ad: `layers` must be a named, non-empty list.",
         call. = FALSE)

  for (lname in names(layers)) {
    l <- layers[[lname]]
    needed <- c("path", "encoding", "shape", "genes")
    miss <- setdiff(needed, names(l))
    if (length(miss))
      stop(sprintf(
        "expression_backend_h5ad: layer '%s' missing field(s): %s",
        lname, paste(miss, collapse = ", ")), call. = FALSE)
    if (!(l$encoding %in% H5AD_SUPPORTED_ENCODINGS))
      stop(sprintf(
        "expression_backend_h5ad: layer '%s' encoding '%s' not supported. Have: %s",
        lname, l$encoding,
        paste(H5AD_SUPPORTED_ENCODINGS, collapse = ", ")),
        call. = FALSE)
    if (length(l$shape) != 2L)
      stop(sprintf(
        "expression_backend_h5ad: layer '%s' shape must be c(n_obs, n_var).",
        lname), call. = FALSE)
  }

  if (is.null(n_cells)) n_cells <- as.integer(layers[[1]]$shape[1])
  n_cells <- as.integer(n_cells)
  if (!default_layer %in% names(layers)) default_layer <- names(layers)[1]

  # Per-layer in-memory cache (env semantics so writes inside
  # backend_get_gene survive copy-on-modify of the outer list).
  caches <- stats::setNames(
    lapply(names(layers), function(.) new.env(parent = emptyenv())),
    names(layers))

  structure(
    list(
      path          = path,
      layers        = layers,
      caches        = caches,
      n_cells       = n_cells,
      default_layer = default_layer
    ),
    class = c("expression_backend_h5ad", "expression_backend")
  )
}

# ---- Methods --------------------------------------------------------------

backend_n_cells.expression_backend_h5ad <- function(backend) {
  as.integer(backend$n_cells)
}

backend_available_layers.expression_backend_h5ad <- function(backend) {
  names(backend$layers)
}

backend_default_layer.expression_backend_h5ad <- function(backend) {
  backend$default_layer
}

backend_genes.expression_backend_h5ad <- function(backend, layer = NULL) {
  layer <- .resolve_layer(backend, layer)
  as.character(backend$layers[[layer]]$genes)
}

backend_n_genes.expression_backend_h5ad <- function(backend, layer = NULL) {
  length(backend_genes(backend, layer = layer))
}

backend_has_gene.expression_backend_h5ad <- function(backend, gene, layer = NULL) {
  if (is.null(gene) || !is.character(gene) || !nzchar(gene)) return(FALSE)
  gene %in% backend_genes(backend, layer = layer)
}

backend_get_gene.expression_backend_h5ad <- function(backend, gene, layer = NULL) {
  if (!backend_has_gene(backend, gene, layer = layer)) return(NULL)
  layer <- .resolve_layer(backend, layer)
  spec  <- backend$layers[[layer]]
  j     <- match(gene, spec$genes)  # 1-based
  switch(spec$encoding,
    "csr_matrix" = .h5ad_get_gene_csr(backend, layer, j),
    "csc_matrix" = .h5ad_get_gene_csc(backend, layer, j),
    "dense"      = .h5ad_get_gene_dense(backend, layer, j),
    stop("Unsupported encoding: ", spec$encoding, call. = FALSE)
  )
}

backend_as_matrix.expression_backend_h5ad <- function(backend, layer = NULL) {
  require_optional("rhdf5",
                   feature = "AnnData (.h5ad) full-matrix materialisation",
                   source  = "Bioconductor")
  layer <- .resolve_layer(backend, layer)
  spec  <- backend$layers[[layer]]
  genes <- as.character(spec$genes)
  n_obs <- as.integer(spec$shape[1])
  n_var <- as.integer(spec$shape[2])

  switch(spec$encoding,
    "csr_matrix" = .h5ad_materialise_csr(backend, layer, genes, n_obs, n_var),
    "csc_matrix" = .h5ad_materialise_csc(backend, layer, genes, n_obs, n_var),
    "dense"      = .h5ad_materialise_dense(backend, layer, genes, n_obs, n_var),
    stop("Unsupported encoding: ", spec$encoding, call. = FALSE)
  )
}

# ---- Cache helpers --------------------------------------------------------
#
# Lazy-loads the sparse triple (indptr, indices, data) for a layer on
# first call. `backend$caches[[layer]]` is an env, so writes here are
# visible to subsequent calls.

.h5ad_load_sparse_triple <- function(backend, layer) {
  cache <- backend$caches[[layer]]
  if (isTRUE(cache$loaded)) return(invisible(NULL))
  require_optional("rhdf5",
                   feature = "AnnData (.h5ad) sparse layer access",
                   source  = "Bioconductor")
  p <- backend$layers[[layer]]$path
  cache$indptr  <- as.integer(rhdf5::h5read(backend$path,
                                            paste0(p, "/indptr")))
  cache$indices <- as.integer(rhdf5::h5read(backend$path,
                                            paste0(p, "/indices")))
  cache$data    <- as.numeric(rhdf5::h5read(backend$path,
                                            paste0(p, "/data")))
  cache$loaded  <- TRUE
  invisible(NULL)
}

# ---- Per-gene readers -----------------------------------------------------

# CSR: the matrix is (n_obs x n_var) and storage is row-major over
# cells. To extract a single gene (one column of the conceptual matrix)
# we have to walk every cell row, look up the gene's 0-based column
# index in indices[indptr[i]:(indptr[i+1]-1)] and pull the matching
# data value when present. AnnData guarantees indices within a row are
# sorted ascending, so we use a binary search (findInterval) instead
# of a linear scan within the row segment.
.h5ad_get_gene_csr <- function(backend, layer, j) {
  .h5ad_load_sparse_triple(backend, layer)
  cache  <- backend$caches[[layer]]
  n_obs  <- as.integer(backend$layers[[layer]]$shape[1])
  j0     <- j - 1L  # AnnData stores 0-indexed columns
  out    <- numeric(n_obs)
  indptr <- cache$indptr
  indices <- cache$indices
  data    <- cache$data
  for (i in seq_len(n_obs)) {
    start <- indptr[i] + 1L
    end   <- indptr[i + 1L]
    if (start > end) next
    seg <- indices[start:end]
    pos <- .findInterval_exact(seg, j0)
    if (!is.na(pos)) out[i] <- data[start + pos - 1L]
  }
  out
}

# CSC: the matrix is still (n_obs x n_var) conceptually, but storage is
# column-major over genes. A single gene = a single contiguous segment
# of `data` / `indices`. The indices entries are 0-based cell row
# indices for the conceptual matrix.
.h5ad_get_gene_csc <- function(backend, layer, j) {
  .h5ad_load_sparse_triple(backend, layer)
  cache  <- backend$caches[[layer]]
  n_obs  <- as.integer(backend$layers[[layer]]$shape[1])
  start  <- cache$indptr[j] + 1L
  end    <- cache$indptr[j + 1L]
  out    <- numeric(n_obs)
  if (end >= start) {
    rows0 <- cache$indices[start:end]
    out[rows0 + 1L] <- cache$data[start:end]
  }
  out
}

# Dense (cells x genes per AnnData spec: HDF5 shape = (n_obs, n_var)).
# rhdf5 with the default `native = FALSE` returns dim order matching
# HDF5, so for a dense X we expect dim(emb) = c(n_obs, n_var) and gene j
# lives in column j. We don't trust this blindly: we sanity-check the
# slice length against the spec's `n_obs`, and fall back to a row slice
# (transposed layout) if needed. This is the only place in the backend
# that depends on rhdf5's dim convention; if a future rhdf5 release
# flips the default, the fallback keeps reads correct.
.h5ad_get_gene_dense <- function(backend, layer, j) {
  require_optional("rhdf5",
                   feature = "AnnData (.h5ad) dense layer access",
                   source  = "Bioconductor")
  spec  <- backend$layers[[layer]]
  p     <- spec$path
  n_obs <- as.integer(spec$shape[1])
  # Primary path: (n_obs, n_var) -> select column j.
  col <- rhdf5::h5read(backend$path, p,
                       index = list(NULL, as.integer(j)))
  vec <- as.numeric(col)
  if (length(vec) == n_obs) return(vec)
  # Fallback: dim order is reversed -> select row j.
  row <- rhdf5::h5read(backend$path, p,
                       index = list(as.integer(j), NULL))
  vec <- as.numeric(row)
  if (length(vec) != n_obs)
    stop(sprintf(
      "h5ad dense reader: expected %d values for one gene, got %d. ",
      n_obs, length(vec)),
      "Layer spec shape may not match the on-disk dataset.",
      call. = FALSE)
  vec
}

# Internal: ordered intersection. `vec` is the sorted indices segment
# of a CSR row (per AnnData spec rows are sorted); `target` is the
# integer column index to find. Returns 1-based position within `vec`
# or NA if absent.
.findInterval_exact <- function(vec, target) {
  if (length(vec) == 0L) return(NA_integer_)
  # findInterval gives last index where vec[idx] <= target.
  k <- findInterval(target, vec)
  if (k >= 1L && vec[k] == target) k else NA_integer_
}

# ---- Full-matrix materialisers --------------------------------------------
#
# Used by `backend_as_matrix()`. These return base R dense matrices
# with `genes` rownames -- consumers (AUCell, SingleR, presto, etc.)
# never see the HDF5 layer.

.h5ad_materialise_csr <- function(backend, layer, genes, n_obs, n_var) {
  .h5ad_load_sparse_triple(backend, layer)
  cache  <- backend$caches[[layer]]
  # Build a dense n_var x n_obs (genes x cells) matrix by walking every
  # row's segment. Cost: O(nnz) -- still feasible if you asked for the
  # whole thing.
  m <- matrix(0, nrow = n_var, ncol = n_obs,
              dimnames = list(genes, NULL))
  for (i in seq_len(n_obs)) {
    start <- cache$indptr[i] + 1L
    end   <- cache$indptr[i + 1L]
    if (start > end) next
    seg_cols <- cache$indices[start:end] + 1L   # 0- -> 1-indexed
    seg_vals <- cache$data[start:end]
    m[seg_cols, i] <- seg_vals
  }
  m
}

.h5ad_materialise_csc <- function(backend, layer, genes, n_obs, n_var) {
  .h5ad_load_sparse_triple(backend, layer)
  cache <- backend$caches[[layer]]
  m <- matrix(0, nrow = n_var, ncol = n_obs,
              dimnames = list(genes, NULL))
  for (j in seq_len(n_var)) {
    start <- cache$indptr[j] + 1L
    end   <- cache$indptr[j + 1L]
    if (start > end) next
    rows0 <- cache$indices[start:end] + 1L
    m[j, rows0] <- cache$data[start:end]
  }
  m
}

.h5ad_materialise_dense <- function(backend, layer, genes, n_obs, n_var) {
  require_optional("rhdf5",
                   feature = "AnnData (.h5ad) dense layer access",
                   source  = "Bioconductor")
  p   <- backend$layers[[layer]]$path
  raw <- rhdf5::h5read(backend$path, p)
  # AnnData dense storage has HDF5 shape (n_obs, n_var). Under
  # rhdf5 `native = FALSE` (the default), `dim(raw)` is (n_obs, n_var)
  # and we transpose for our genes x cells convention. Under
  # `native = TRUE`, dim is already (n_var, n_obs) and we keep it as-is.
  # Sanity-check the row count rather than trusting the default blindly.
  if (identical(dim(raw), c(n_obs, n_var))) {
    m <- t(raw)
  } else if (identical(dim(raw), c(n_var, n_obs))) {
    m <- raw
  } else {
    stop(sprintf(
      "h5ad dense materialise: expected dim (%d, %d) or (%d, %d), got (%s).",
      n_obs, n_var, n_var, n_obs,
      paste(dim(raw), collapse = ", ")), call. = FALSE)
  }
  rownames(m) <- genes
  m
}

# ---- Print method ---------------------------------------------------------

print.expression_backend_h5ad <- function(x, ...) {
  cat("<expression_backend_h5ad>\n")
  cat("  path        :", x$path, "\n")
  cat("  n_cells     :", x$n_cells, "\n")
  cat("  layers      :", paste(names(x$layers), collapse = ", "), "\n")
  cat("  default     :", x$default_layer, "\n")
  for (lname in names(x$layers)) {
    l <- x$layers[[lname]]
    cat(sprintf("    %-8s | %-11s | n_var = %d\n",
                lname, l$encoding, as.integer(l$shape[2])))
  }
  invisible(x)
}
