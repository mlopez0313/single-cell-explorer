# ============================================================================
# Dataset helpers
# ----------------------------------------------------------------------------
# Small, pure functions that read from a dataset list. Modules should go
# through these helpers instead of touching dataset internals -- this is the
# single seam where future Seurat / AnnData loaders need to behave correctly.
# ============================================================================

#' List the assays available in a dataset.
available_assays <- function(dataset) {
  if (is.null(dataset)) character() else dataset$assays
}

#' List the dimensional reductions available in a dataset.
available_reductions <- function(dataset) {
  if (is.null(dataset)) character() else dataset$reductions
}

#' List the metadata column names available in a dataset.
available_metadata_fields <- function(dataset) {
  if (is.null(dataset)) character() else dataset$metadata_fields
}

#' List gene symbols whose expression is queryable for the current dataset.
#'
#' For the mock dataset this is the small set of demo genes. Future loaders
#' should expose the full set of expressed features (or a sensible subset).
#'
#' Prefers `dataset$genes` when set (cheap, identical to the backend's
#' default-layer genes for any well-formed loader). Falls back to the
#' backend so any future loader that forgets to populate `dataset$genes`
#' still produces a usable answer.
#'
#' @param dataset a dataset list satisfying `dataset_schema()`
#' @param layer   optional layer name (defaults to the backend's default).
#'   Currently only meaningful for backends with multiple layers.
available_genes <- function(dataset, layer = NULL) {
  if (is.null(dataset)) return(character())
  if (is.null(layer) && !is.null(dataset$genes)) return(dataset$genes)
  be <- as_expression_backend(dataset$expression)
  backend_genes(be, layer = layer)
}

#' Return a data.frame with columns `x`, `y`, and `cell` for the requested
#' reduction. Returns NULL if the reduction isn't available.
get_embedding <- function(dataset, reduction) {
  if (is.null(dataset) || !nzchar(reduction %||% "")) return(NULL)
  if (!reduction %in% available_reductions(dataset))  return(NULL)
  cd <- dataset$cell_data
  xcol <- paste0(reduction, "_1")
  ycol <- paste0(reduction, "_2")
  if (!all(c(xcol, ycol) %in% names(cd))) return(NULL)
  data.frame(cell = cd$cell, x = cd[[xcol]], y = cd[[ycol]],
             stringsAsFactors = FALSE)
}

#' Return a metadata vector for `field`, or NULL if not available.
get_metadata <- function(dataset, field) {
  if (is.null(dataset) || !nzchar(field %||% ""))            return(NULL)
  if (!field %in% available_metadata_fields(dataset))        return(NULL)
  if (!field %in% names(dataset$cell_data))                  return(NULL)
  dataset$cell_data[[field]]
}

#' Return per-cell expression values for `gene` (numeric vector of length
#' n_cells), or NULL if the gene is not present in the dataset / layer.
#'
#' Dispatches through the `expression_backend` so future loaders can lazy-
#' read from sparse matrices or HDF5 without changing module code. Legacy
#' datasets where `dataset$expression` is still a flat named list are
#' coerced into an in-memory backend on the fly.
#'
#' @param dataset a dataset list satisfying `dataset_schema()`
#' @param gene    character(1) gene symbol
#' @param layer   optional layer name (defaults to the backend's default
#'   layer, currently "data" for the in-memory backend). Reserved for
#'   pseudobulk DE and future counts-vs-data workflows.
get_gene_expression <- function(dataset, gene, layer = NULL) {
  # Defensive NULL/NA/empty checks. `nzchar(NA)` is NA, which makes the
  # surrounding `if(...)` throw "missing value where TRUE/FALSE needed",
  # so callers that pass `state$selected_gene` while it's still
  # `NA_character_` (set by `set_active_dataset()` when a dataset
  # exposes no gene names) used to crash any renderUI calling
  # `validate_gene()`. Treat NULL/NA/empty string as "no gene".
  if (is.null(dataset)) return(NULL)
  if (is.null(gene) || length(gene) != 1L) return(NULL)
  if (is.na(gene) || !nzchar(gene)) return(NULL)
  be  <- as_expression_backend(dataset$expression)
  out <- backend_get_gene(be, gene, layer = layer)
  if (is.null(out)) return(NULL)
  as.numeric(out)
}

#' TRUE/FALSE: is `gene` present in the current dataset?
validate_gene <- function(dataset, gene, layer = NULL) {
  !is.null(get_gene_expression(dataset, gene, layer = layer))
}

#' TRUE/FALSE: is `field` present as a metadata column?
validate_metadata <- function(dataset, field) {
  !is.null(get_metadata(dataset, field))
}
