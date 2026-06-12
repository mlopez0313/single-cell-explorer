# ============================================================================
# Expression backend abstraction
# ----------------------------------------------------------------------------
# `dataset$expression` is no longer a bare named list of numeric vectors --
# it is an S3 `expression_backend` object. Modules never touch it directly;
# they go through R/dataset_helpers.R which dispatches to the backend.
#
# Why this matters:
#   * The mock dataset uses an in-memory list of numeric vectors, one per
#     gene. That's fine for ~10 genes x 2500 cells.
#   * Real datasets (Seurat .rds, AnnData .h5ad, 10x directories) cannot
#     afford to materialise an n_genes x n_cells dense matrix in memory.
#     They need lazy per-gene access against a sparse matrix or an HDF5
#     handle.
#   * Pseudobulk DE (Cursor prompt #5) will need a separate "counts" layer
#     beside the log-normalised "data" layer. The backend reserves a
#     per-layer interface now so we don't have to break the schema again.
#
# Today's backends:
#   * `expression_backend_inmemory()`  -- named list (or list-of-lists keyed
#                                         by layer) of numeric vectors of
#                                         length n_cells. Used by the mock
#                                         dataset and any small in-memory
#                                         dataset.
#
# Designed-for, not-yet-implemented (each is a separate Cursor prompt):
#   * `expression_backend_sparse()`    -- Matrix::dgCMatrix per layer
#                                         (Seurat / Read10X output)
#   * `expression_backend_h5ad()`      -- AnnData / HDF5 handle, lazy reads
#   * `expression_backend_seurat()`    -- thin proxy over a Seurat object
#
# The generics each new backend must implement:
#   backend_n_cells(backend)
#   backend_available_layers(backend)
#   backend_default_layer(backend)
#   backend_genes(backend, layer = NULL)            -> character()
#   backend_n_genes(backend, layer = NULL)          -> integer(1)
#   backend_has_gene(backend, gene, layer = NULL)   -> logical(1)
#   backend_get_gene(backend, gene, layer = NULL)   -> numeric(n_cells) | NULL
# ============================================================================

# ---- Generics ------------------------------------------------------------

#' Number of cells the backend describes.
backend_n_cells          <- function(backend) UseMethod("backend_n_cells")

#' Layer names supported by the backend.
backend_available_layers <- function(backend) UseMethod("backend_available_layers")

#' The layer name used when callers omit `layer`.
backend_default_layer    <- function(backend) UseMethod("backend_default_layer")

#' Genes present in the given layer (defaults to the backend's default).
backend_genes            <- function(backend, layer = NULL) UseMethod("backend_genes")

#' Number of genes present in the given layer.
backend_n_genes          <- function(backend, layer = NULL) UseMethod("backend_n_genes")

#' Is `gene` queryable in the given layer?
backend_has_gene         <- function(backend, gene, layer = NULL) UseMethod("backend_has_gene")

#' Per-cell numeric vector for `gene` in `layer`. NULL if the gene is absent.
#' Implementations must return a vector of length `backend_n_cells(backend)`
#' aligned to the dataset's `cell_data$cell` order.
backend_get_gene         <- function(backend, gene, layer = NULL) UseMethod("backend_get_gene")

#' Return a (genes x cells) matrix view for the requested layer.
#'
#' Used by backends that need to hand the underlying matrix to a third
#' party (SingleR, presto, future pseudobulk aggregators) without
#' funneling every gene through `backend_get_gene()` one at a time.
#'
#' Implementations should return whatever native matrix-like the layer
#' is already stored as -- dense matrix for in-memory, dgCMatrix for
#' sparse, DelayedArray for h5ad -- so the caller can choose whether to
#' materialise.
backend_as_matrix        <- function(backend, layer = NULL) UseMethod("backend_as_matrix")

# ---- Coercion ------------------------------------------------------------

#' Coerce any expression representation into an `expression_backend`.
#'
#' Accepts:
#'   * an existing `expression_backend`           -> returned unchanged
#'   * `NULL`                                     -> empty in-memory backend
#'   * an empty list                              -> empty in-memory backend
#'   * a flat named list of numeric vectors of
#'     equal length (legacy mock shape)           -> wrapped as a single-layer
#'                                                   in-memory backend
#'                                                   (layer name "data")
#'
#' Anything else raises a clear error. The point of this helper is that
#' every legacy dataset object (where `expression` is still a bare named
#' list) keeps working when fed to the new helpers.
as_expression_backend <- function(x) {
  if (inherits(x, "expression_backend")) return(x)
  if (is.null(x)) {
    return(expression_backend_inmemory(list(), n_cells = 0L))
  }
  if (is.list(x)) {
    if (length(x) == 0L) {
      return(expression_backend_inmemory(list(), n_cells = 0L))
    }
    if (all(vapply(x, is.numeric, logical(1)))) {
      return(expression_backend_inmemory(x))
    }
  }
  stop("Cannot coerce object of class '",
       paste(class(x), collapse = "/"),
       "' to an expression_backend.", call. = FALSE)
}

# ---- In-memory backend ---------------------------------------------------

#' Build an in-memory expression backend.
#'
#' @param layers Either:
#'   * a named list of numeric vectors (flat / legacy shape): wrapped as a
#'     single layer named `default_layer`, OR
#'   * a named list of layers, where each element is itself a named list
#'     of numeric vectors keyed by gene (e.g. `list(data = ..., counts = ...)`).
#' @param n_cells expected vector length per gene. Inferred from the first
#'   gene vector if not provided.
#' @param default_layer the layer used when callers don't specify one.
expression_backend_inmemory <- function(layers,
                                        n_cells = NULL,
                                        default_layer = "data") {
  if (!is.list(layers)) {
    stop("expression_backend_inmemory(layers = ...) must be a list.",
         call. = FALSE)
  }

  # Disambiguate flat-list shape vs. layered shape.
  flat_named_list <- length(layers) > 0L &&
    all(vapply(layers, is.numeric, logical(1)))

  if (flat_named_list) {
    layers <- stats::setNames(list(layers), default_layer)
  } else if (length(layers) == 0L) {
    layers <- stats::setNames(list(list()), default_layer)
  } else {
    bad <- !vapply(layers, function(l) {
      is.list(l) && (length(l) == 0L ||
                     all(vapply(l, is.numeric, logical(1))))
    }, logical(1))
    if (any(bad)) {
      stop("expression_backend_inmemory: layer(s) '",
           paste(names(layers)[bad], collapse = "', '"),
           "' must be named lists of numeric vectors.",
           call. = FALSE)
    }
    if (is.null(names(layers)) || any(!nzchar(names(layers)))) {
      stop("expression_backend_inmemory: layered shape requires named layers.",
           call. = FALSE)
    }
  }

  # Determine / validate n_cells.
  if (is.null(n_cells)) {
    seen <- 0L
    for (l in layers) {
      if (length(l) > 0L) { seen <- length(l[[1]]); break }
    }
    n_cells <- as.integer(seen)
  } else {
    n_cells <- as.integer(n_cells)
  }

  for (lname in names(layers)) {
    l <- layers[[lname]]
    if (length(l) == 0L) next
    if (is.null(names(l)) || any(!nzchar(names(l)))) {
      stop("expression_backend_inmemory: layer '", lname,
           "' must have named gene vectors.", call. = FALSE)
    }
    dups <- duplicated(names(l))
    if (any(dups)) {
      stop("expression_backend_inmemory: duplicate gene name(s) in layer '",
           lname, "': ", paste(unique(names(l)[dups]), collapse = ", "),
           call. = FALSE)
    }
    bad <- !vapply(l, function(v) length(v) == n_cells, logical(1))
    if (any(bad)) {
      stop("expression_backend_inmemory: in layer '", lname,
           "', gene(s) ", paste(names(l)[bad], collapse = ", "),
           " have lengths different from n_cells (", n_cells, ").",
           call. = FALSE)
    }
  }

  # Default-layer fallback if the requested default isn't actually present.
  if (!default_layer %in% names(layers)) {
    default_layer <- names(layers)[1] %||% "data"
  }

  structure(
    list(
      layers        = layers,
      n_cells       = n_cells,
      default_layer = default_layer
    ),
    class = c("expression_backend_inmemory", "expression_backend")
  )
}

# ---- In-memory methods ---------------------------------------------------

backend_n_cells.expression_backend_inmemory <- function(backend) {
  as.integer(backend$n_cells)
}

backend_available_layers.expression_backend_inmemory <- function(backend) {
  names(backend$layers)
}

backend_default_layer.expression_backend_inmemory <- function(backend) {
  backend$default_layer
}

.resolve_layer <- function(backend, layer) {
  if (is.null(layer) || !nzchar(layer)) {
    return(backend_default_layer(backend))
  }
  avail <- backend_available_layers(backend)
  if (!layer %in% avail) {
    stop("Layer '", layer, "' not available. Have: ",
         paste(avail, collapse = ", "), call. = FALSE)
  }
  layer
}

backend_genes.expression_backend_inmemory <- function(backend, layer = NULL) {
  layer <- .resolve_layer(backend, layer)
  out <- names(backend$layers[[layer]])
  if (is.null(out)) character() else out
}

backend_n_genes.expression_backend_inmemory <- function(backend, layer = NULL) {
  length(backend_genes(backend, layer = layer))
}

backend_has_gene.expression_backend_inmemory <- function(backend, gene, layer = NULL) {
  if (is.null(gene) || !is.character(gene) || !nzchar(gene)) return(FALSE)
  gene %in% backend_genes(backend, layer = layer)
}

backend_get_gene.expression_backend_inmemory <- function(backend, gene, layer = NULL) {
  if (!backend_has_gene(backend, gene, layer = layer)) return(NULL)
  layer <- .resolve_layer(backend, layer)
  as.numeric(backend$layers[[layer]][[gene]])
}

backend_as_matrix.expression_backend_inmemory <- function(backend, layer = NULL) {
  layer <- .resolve_layer(backend, layer)
  genes <- backend_genes(backend, layer = layer)
  n     <- backend_n_cells(backend)
  if (length(genes) == 0L) {
    return(matrix(numeric(0), nrow = 0, ncol = n,
                  dimnames = list(NULL, NULL)))
  }
  m <- do.call(rbind, backend$layers[[layer]])
  rownames(m) <- genes
  m
}

# =========================================================================
# Sparse backend
# -------------------------------------------------------------------------
# Wraps one or more genes-x-cells sparse matrices (`Matrix::dgCMatrix`,
# typically) so per-gene reads don't materialise a dense n_genes x n_cells
# matrix. Used by the real Seurat / 10x / AnnData loaders.
#
# Layered shape:
#   layers = list(data = <dgCMatrix>, counts = <dgCMatrix>)
# Different layers may have different gene sets (e.g. counts holds all
# features, data is variable features only) but must share `ncol == n_cells`
# and a consistent column order. Modules see whichever layer they ask for.
#
# Per-gene access: `mat[gene, , drop = TRUE]` on a dgCMatrix returns a
# dense numeric vector of length ncol -- exactly what `get_gene_expression`
# returns to the rest of the app.
# =========================================================================

#' Build a sparse-matrix expression backend.
#'
#' @param layers Either:
#'   * a single object with `Matrix`-compatible row/col indexing (a sparse
#'     `dgCMatrix`, dense matrix, or DelayedArray), lifted as the default
#'     layer; OR
#'   * a named list of such objects keyed by layer name.
#' @param n_cells expected number of columns. Inferred from the matrix if
#'   not provided. Used to validate every layer agrees.
#' @param default_layer the layer used when callers don't pass `layer`.
expression_backend_sparse <- function(layers,
                                      n_cells = NULL,
                                      default_layer = "data") {
  # Single matrix -> single-layer
  if (!is.list(layers) || is.matrix(layers) ||
      inherits(layers, c("Matrix", "DelayedMatrix"))) {
    layers <- stats::setNames(list(layers), default_layer)
  }
  if (length(layers) == 0L || is.null(names(layers)) ||
      any(!nzchar(names(layers)))) {
    stop("expression_backend_sparse: layered shape requires named layers.",
         call. = FALSE)
  }

  # Validate each layer is matrix-like with row + col names
  for (lname in names(layers)) {
    m <- layers[[lname]]
    if (is.null(dim(m)) || length(dim(m)) != 2L) {
      stop("expression_backend_sparse: layer '", lname,
           "' is not a 2D matrix.", call. = FALSE)
    }
    if (is.null(rownames(m))) {
      stop("expression_backend_sparse: layer '", lname,
           "' has no rownames; cannot resolve genes by name.",
           call. = FALSE)
    }
  }

  if (is.null(n_cells)) n_cells <- ncol(layers[[1]])
  n_cells <- as.integer(n_cells)

  # All layers must have the same number of columns.
  ncols <- vapply(layers, ncol, integer(1))
  if (any(ncols != n_cells)) {
    stop("expression_backend_sparse: layer(s) ",
         paste(names(layers)[ncols != n_cells], collapse = ", "),
         " have ncol != n_cells (", n_cells, ").", call. = FALSE)
  }

  if (!default_layer %in% names(layers)) {
    default_layer <- names(layers)[1]
  }

  structure(
    list(
      layers        = layers,
      n_cells       = n_cells,
      default_layer = default_layer
    ),
    class = c("expression_backend_sparse", "expression_backend")
  )
}

# ---- Sparse methods ------------------------------------------------------

backend_n_cells.expression_backend_sparse <- function(backend) {
  as.integer(backend$n_cells)
}

backend_available_layers.expression_backend_sparse <- function(backend) {
  names(backend$layers)
}

backend_default_layer.expression_backend_sparse <- function(backend) {
  backend$default_layer
}

backend_genes.expression_backend_sparse <- function(backend, layer = NULL) {
  layer <- .resolve_layer(backend, layer)
  rn <- rownames(backend$layers[[layer]])
  if (is.null(rn)) character() else rn
}

backend_n_genes.expression_backend_sparse <- function(backend, layer = NULL) {
  length(backend_genes(backend, layer = layer))
}

backend_has_gene.expression_backend_sparse <- function(backend, gene, layer = NULL) {
  if (is.null(gene) || !is.character(gene) || !nzchar(gene)) return(FALSE)
  gene %in% backend_genes(backend, layer = layer)
}

backend_get_gene.expression_backend_sparse <- function(backend, gene, layer = NULL) {
  if (!backend_has_gene(backend, gene, layer = layer)) return(NULL)
  layer <- .resolve_layer(backend, layer)
  # Row-indexing on a dgCMatrix returns a dense numeric vector. as.numeric()
  # handles DelayedArray / Matrix / base matrix uniformly.
  row <- backend$layers[[layer]][gene, , drop = TRUE]
  as.numeric(row)
}

backend_as_matrix.expression_backend_sparse <- function(backend, layer = NULL) {
  layer <- .resolve_layer(backend, layer)
  # Hand back the underlying matrix (dgCMatrix / DelayedArray / dense).
  # Callers that need a dense view can `as.matrix()` themselves.
  backend$layers[[layer]]
}

# ---- Print / format ------------------------------------------------------

print.expression_backend <- function(x, ...) {
  cat("<expression_backend>\n")
  cat("  class       :", paste(class(x), collapse = " / "), "\n")
  cat("  n_cells     :", backend_n_cells(x), "\n")
  cat("  layers      :", paste(backend_available_layers(x), collapse = ", "), "\n")
  cat("  default     :", backend_default_layer(x), "\n")
  cat("  n_genes[def]:", backend_n_genes(x), "\n")
  invisible(x)
}

format.expression_backend <- function(x, ...) {
  sprintf("<%s: %d cells, %d gene(s), %d layer(s)>",
          class(x)[1], backend_n_cells(x), backend_n_genes(x),
          length(backend_available_layers(x)))
}
