# ============================================================================
# Dataset loading (placeholders)
# ----------------------------------------------------------------------------
# The app is designed to eventually accept multiple dataset sources:
#   - Seurat objects (.rds)
#   - AnnData objects (.h5ad)
#   - 10x Genomics raw output directories
#
# For now we expose:
#   - `mock_dataset()`         : an in-memory fake dataset so the UI can be
#                                exercised without real data
#   - `load_dataset(path, ...)`: a dispatcher with stubs for each source type.
#                                Stubs raise a friendly "not implemented" error.
#
# A loaded dataset is a plain list with a stable schema so downstream modules
# do not need to know which file format it came from. See `dataset_schema()`.
# ============================================================================

#' The canonical schema all loaders must return.
#'
#' Future loaders (Seurat / AnnData / 10x) must populate every field below.
#' Keeping the schema flat keeps modules simple and decouples them from the
#' source format.
dataset_schema <- function() {
  c(
    "name",              # character(1) dataset display name
    "source",            # one of "mock", "seurat", "anndata", "10x"
    "n_cells",           # integer(1)
    "n_genes",           # integer(1)
    "assays",            # character() available assay names
    "default_assay",     # character(1)
    "reductions",        # character() e.g. c("PCA", "UMAP", "tSNE")
    "default_reduction", # character(1)
    "metadata_fields",   # character() cell-level metadata column names
    "cells",             # character() cell barcodes / ids
    "cell_data",         # data.frame; one row per cell. Required columns:
                         #   cell + every `metadata_fields` value +
                         #   <reduction>_1 / <reduction>_2 for every reduction
    "genes",             # character() gene symbols present in `expression`
    "expression"         # expression_backend object (see R/expression_backend.R).
                         # Modules MUST NOT touch this directly; access through
                         # `get_gene_expression(dataset, gene, layer = NULL)`
                         # and `available_genes(dataset)`. Legacy plain
                         # named-list expression is also accepted by the
                         # helpers (coerced transparently on read).
  )
}

# Example genes used by both the mock dataset and the Explorer module. Real
# loaders are free to expose any gene symbols; this list is just a convenient
# default for UI testing.
MOCK_GENES <- c("CD3D", "MS4A1", "LST1", "EPCAM", "COL1A1", "NKG7")

#' Build a small mock dataset useful for UI development and CI.
#'
#' The mock data is structured as four loose Gaussian "clusters" in UMAP
#' space, each loosely mapped to a cell-type label. Mock expression for the
#' genes in `MOCK_GENES` is generated so each gene is preferentially
#' expressed in one of the clusters. This is enough to exercise the
#' metadata-colored embedding and FeaturePlot-style UI without any real
#' biology.
#'
#' Returns a list satisfying `dataset_schema()`.
mock_dataset <- function(n_cells = 2500L, n_genes = 15000L,
                         name = "mock_pbmc_2.5k",
                         seed = 42L) {
  set.seed(seed)
  n_cells <- as.integer(n_cells)

  cells <- sprintf("cell_%05d", seq_len(n_cells))

  # Four cluster centers in 2D embedding space (UMAP). Cell-type mapping is
  # one-to-one with cluster id for clarity.
  cluster_ids   <- 0:3
  cluster_types <- c("T cell", "B cell", "Myeloid", "Epithelial")
  centers <- matrix(c( 6,  6,
                      -6,  6,
                      -6, -6,
                       6, -6), ncol = 2, byrow = TRUE)

  cluster_of_cell <- sample(cluster_ids, n_cells, replace = TRUE,
                            prob = c(0.35, 0.25, 0.25, 0.15))

  umap <- centers[cluster_of_cell + 1L, ] +
          matrix(stats::rnorm(2 * n_cells, sd = 1.4), ncol = 2)

  # PCA and tSNE are intentionally placeholder coordinates derived from UMAP
  # plus noise -- enough that the UI shows different layouts when switched.
  pca   <- cbind(umap[, 1] * 0.8 + stats::rnorm(n_cells, sd = 0.5),
                 umap[, 2] * 0.5 + stats::rnorm(n_cells, sd = 0.5))
  tsne  <- cbind(umap[, 1] + umap[, 2] * 0.3 + stats::rnorm(n_cells, sd = 0.7),
                 umap[, 2] - umap[, 1] * 0.3 + stats::rnorm(n_cells, sd = 0.7))

  # Precomputed "demo" pseudotime: normalised UMAP distance from cluster 0's
  # centroid. Exposed as a numeric metadata column so the Trajectory module
  # has a non-empty option under the "metadata" source out of the box.
  # See R/trajectory.R for the documented contract.
  root_cells_demo <- cluster_of_cell == 0L
  pt_centroid <- c(mean(umap[root_cells_demo, 1]),
                   mean(umap[root_cells_demo, 2]))
  pt_demo_raw <- sqrt((umap[, 1] - pt_centroid[1])^2 +
                      (umap[, 2] - pt_centroid[2])^2)
  pt_demo <- (pt_demo_raw - min(pt_demo_raw)) /
             (max(pt_demo_raw) - min(pt_demo_raw))

  cell_data <- data.frame(
    cell           = cells,
    sample         = sample(c("S1", "S2", "S3"), n_cells, replace = TRUE),
    cluster        = as.character(cluster_of_cell),
    condition      = sample(c("ctrl", "treat"),  n_cells, replace = TRUE),
    cell_type      = cluster_types[cluster_of_cell + 1L],
    pseudotime_demo = pt_demo,
    UMAP_1         = umap[, 1], UMAP_2 = umap[, 2],
    PCA_1          = pca[, 1],  PCA_2  = pca[, 2],
    tSNE_1         = tsne[, 1], tSNE_2 = tsne[, 2],
    stringsAsFactors = FALSE
  )

  # Mock expression: each gene is high in a chosen cluster, low elsewhere,
  # plus noise. The "data" layer holds roughly log-normalised values
  # (range ~0-5) used by cell-level DE / markers / pathway analysis.
  # A parallel "counts" layer holds raw count-style integers (Poisson
  # draws around the log-normalised mean) used by pseudobulk DE.
  #
  # Layer semantics (real loaders MUST honour these names):
  #   data    -- log-normalised expression. Default layer. Consumed by
  #              cell-level DE, markers, pathway, imputation, plotting.
  #   counts  -- raw counts. Consumed by `aggregate_pseudobulk()` and
  #              pseudobulk DE backends (`pseudobulk_naive`,
  #              `pseudobulk_edger`, `pseudobulk_deseq2`).
  gene_cluster <- c(CD3D = 0, MS4A1 = 1, LST1 = 2, EPCAM = 3, COL1A1 = 3, NKG7 = 0)
  expression_data <- lapply(names(gene_cluster), function(g) {
    base <- ifelse(cluster_of_cell == gene_cluster[[g]],
                   stats::rnorm(n_cells, mean = 3.2, sd = 0.6),
                   stats::rnorm(n_cells, mean = 0.3, sd = 0.4))
    pmax(base, 0)
  })
  names(expression_data) <- names(gene_cluster)

  # Per-cell library size used to scale the count draws (1k-5k UMIs).
  lib_sizes_demo <- as.integer(stats::runif(n_cells, min = 1000, max = 5000))
  # Raw count layer: Poisson draws with cluster-specific rate. The rate
  # reflects the same gene-cluster mapping but with a treatment / condition
  # nudge so pseudobulk DE has a real signal to recover across samples.
  condition_bump <- ifelse(cell_data$condition == "treat", 1.3, 1.0)
  expression_counts <- lapply(names(gene_cluster), function(g) {
    on_target <- cluster_of_cell == gene_cluster[[g]]
    # mean count per cell = base rate scaled by on/off-target, condition
    # bump, and per-cell library-size relative to a 2.5k baseline.
    mu <- ifelse(on_target, 8, 0.5) * condition_bump *
          (lib_sizes_demo / 2500)
    stats::rpois(n_cells, lambda = pmax(mu, 0.01))
  })
  names(expression_counts) <- names(gene_cluster)

  # Wrap the mock expression in an `expression_backend` exposing both
  # layers. Default layer is "data" so existing modules see no change.
  expression_backend <- expression_backend_inmemory(
    list(data = expression_data, counts = expression_counts),
    n_cells = n_cells,
    default_layer = "data"
  )

  list(
    name              = name,
    source            = "mock",
    n_cells           = n_cells,
    n_genes           = as.integer(n_genes),
    assays            = c("RNA", "SCT", "ADT"),
    default_assay     = "RNA",
    reductions        = c("PCA", "UMAP", "tSNE"),
    default_reduction = "UMAP",
    metadata_fields   = c("sample", "cluster", "condition", "cell_type", "pseudotime_demo"),
    cells             = cells,
    cell_data         = cell_data,
    genes             = names(expression_data),
    expression        = expression_backend
  )
}

#' Load a dataset from disk.
#'
#' Dispatches by file extension or explicit `source`. All real loaders are
#' currently stubs that raise an informative error -- they are wired here so a
#' future contributor only needs to fill in one branch.
#'
#' @param path    file path or directory
#' @param source  optional, force a source type ("seurat", "anndata", "10x")
#' @return a dataset list matching `dataset_schema()`
load_dataset <- function(path, source = NULL, ...) {
  if (is.null(source)) source <- detect_source(path)
  switch(source,
    "seurat"  = load_seurat(path, ...),
    "anndata" = load_anndata(path, ...),
    "10x"     = load_10x(path, ...),
    stop("Unknown dataset source: ", source, call. = FALSE)
  )
}

#' Infer the dataset source from a path's shape.
#'
#' Rules:
#'   * `.rds`  -> "seurat"
#'   * `.h5ad` -> "anndata"
#'   * a directory containing a plausible Cellranger v2/v3 layout
#'     (matrix.mtx{,.gz} + barcodes.tsv{,.gz} +
#'      features.tsv{,.gz} or genes.tsv{,.gz}) -> "10x"
#'   * anything else -> error
#'
#' Previously *any* directory was classified as "10x", which produced
#' misleading downstream "matrix.mtx missing" errors. With this tighter
#' check, an arbitrary directory now fails the `detect_source()` step
#' itself ("cannot infer source"), and users who really mean it can
#' still force the loader by calling `load_dataset(path, source = "10x")`.
detect_source <- function(path) {
  if (dir.exists(path)) {
    if (.looks_like_10x_dir(path)) return("10x")
    stop("Cannot infer dataset source: '", path,
         "' is a directory but does not look like a Cellranger ",
         "feature-barcode matrix (expected matrix.mtx[.gz], ",
         "barcodes.tsv[.gz], features.tsv[.gz] or genes.tsv[.gz]). ",
         "Pass `source = \"10x\"` explicitly to force the 10x loader.",
         call. = FALSE)
  }
  ext <- tolower(tools::file_ext(path))
  switch(ext,
    "rds"  = "seurat",
    "h5ad" = "anndata",
    stop("Cannot infer dataset source from path: ", path, call. = FALSE)
  )
}

#' TRUE iff `dir` plausibly contains a Cellranger v2/v3 feature-barcode
#' matrix. Mirrors `.resolve_10x_paths()`'s file-pickling rules but
#' returns a logical instead of raising. Kept private to dataset.R so
#' the inference & loader stay in lockstep.
.looks_like_10x_dir <- function(dir) {
  has_any <- function(stems) {
    for (stem in stems) {
      if (file.exists(file.path(dir, stem)) ||
          file.exists(file.path(dir, paste0(stem, ".gz")))) {
        return(TRUE)
      }
    }
    FALSE
  }
  has_any("matrix.mtx") &&
    has_any("barcodes.tsv") &&
    has_any(c("features.tsv", "genes.tsv"))
}

# ===========================================================================
# Real loaders
# ---------------------------------------------------------------------------
# Each returns a list satisfying `dataset_schema()`. The expression field is
# always an `expression_backend` (in-memory for the mock; sparse for real
# inputs) so modules don't need to know which loader ran.
#
# Required packages are gated through `require_optional()` (see
# R/optional_deps.R) -- missing deps produce a clear error with the install
# command, never a NULL-pointer-style crash.
# ===========================================================================

# ---- Seurat .rds ----------------------------------------------------------

#' Load a Seurat object saved with `saveRDS()`.
#'
#' Produces an `expression_backend_sparse` over the default assay's
#' `data` layer (and `counts` when present). All standard reductions are
#' surfaced as `<RED>_1` / `<RED>_2` columns in `cell_data`.
load_seurat <- function(path, default_assay = NULL, default_reduction = NULL) {
  require_optional(c("SeuratObject"), feature = "Seurat .rds loading")
  if (!file.exists(path)) {
    stop("Seurat loader: file does not exist: ", path, call. = FALSE)
  }

  obj <- tryCatch(readRDS(path),
                  error = function(e)
                    stop("Seurat loader: failed to readRDS('", path, "'): ",
                         conditionMessage(e), call. = FALSE))
  if (!inherits(obj, "Seurat")) {
    stop("Seurat loader: object in '", path,
         "' is class '", paste(class(obj), collapse = "/"),
         "' (expected 'Seurat'). Did you mean load_anndata() or load_10x()?",
         call. = FALSE)
  }
  .seurat_to_dataset(obj, name = .basename_noext(path),
                     default_assay     = default_assay,
                     default_reduction = default_reduction)
}

#' Convert an in-memory Seurat object into the app's dataset schema.
#' Exposed so future "Run analysis on this Seurat object" flows can reuse
#' the same source-to-schema mapping without writing/reading an .rds.
.seurat_to_dataset <- function(obj, name = "seurat",
                               default_assay = NULL,
                               default_reduction = NULL) {
  assays_avail <- SeuratObject::Assays(obj)
  if (length(assays_avail) == 0L) {
    stop("Seurat loader: object has no assays.", call. = FALSE)
  }
  da <- default_assay %||% SeuratObject::DefaultAssay(obj)
  if (!da %in% assays_avail) {
    stop("Seurat loader: requested default_assay '", da,
         "' not in object assays (", paste(assays_avail, collapse = ", "),
         ").", call. = FALSE)
  }

  reds_avail <- SeuratObject::Reductions(obj)
  dr <- default_reduction %||% (if (length(reds_avail) > 0L) reds_avail[1] else NA_character_)

  cells <- colnames(obj)
  meta  <- obj@meta.data
  if (!"cell" %in% names(meta)) meta$cell <- cells

  # Build per-reduction embedding columns. We only surface the first 2
  # dims of each reduction; that's all the FeaturePlot / DimPlot views
  # need. The full embedding matrix is still on the Seurat object if a
  # future module wants it.
  for (r in reds_avail) {
    emb <- SeuratObject::Embeddings(obj, reduction = r)
    if (ncol(emb) >= 2L) {
      meta[[paste0(toupper(r), "_1")]] <- emb[, 1]
      meta[[paste0(toupper(r), "_2")]] <- emb[, 2]
    }
  }

  meta_fields <- setdiff(names(meta),
                         c("cell", unlist(lapply(reds_avail, function(r)
                           paste0(toupper(r), c("_1", "_2"))))))

  # Expression backend: prefer Seurat v5 LayerData; fall back to v4
  # GetAssayData. Raw (un-normalised) Seurat objects ship only the
  # `counts` layer -- we don't normalise on load; we just promote
  # `counts` to the default layer in that case.
  data_mat   <- .seurat_layer(obj, assay = da, layer = "data")
  counts_mat <- .seurat_layer(obj, assay = da, layer = "counts")

  if (is.null(data_mat) && is.null(counts_mat)) {
    stop("Seurat loader: assay '", da,
         "' has neither 'data' nor 'counts' layer.", call. = FALSE)
  }
  layers <- list()
  if (!is.null(data_mat))   layers$data   <- data_mat
  if (!is.null(counts_mat)) layers$counts <- counts_mat
  default_layer <- if (!is.null(data_mat)) "data" else "counts"
  primary_mat   <- layers[[default_layer]]

  expr_be <- expression_backend_sparse(layers, n_cells = length(cells),
                                       default_layer = default_layer)

  list(
    name              = name,
    source            = "seurat",
    n_cells           = as.integer(length(cells)),
    n_genes           = as.integer(nrow(primary_mat)),
    assays            = assays_avail,
    default_assay     = da,
    reductions        = toupper(reds_avail),
    default_reduction = if (is.na(dr)) NA_character_ else toupper(dr),
    metadata_fields   = meta_fields,
    cells             = cells,
    cell_data         = meta,
    genes             = backend_genes(expr_be),
    expression        = expr_be
  )
}

.seurat_layer <- function(obj, assay, layer) {
  fetch <- function(fn, ...) {
    out <- tryCatch(suppressWarnings(fn(obj, ...)),
                    error = function(e) NULL)
    if (is.null(out)) return(NULL)
    if (!is.null(dim(out)) && nrow(out) == 0L) return(NULL)
    out
  }
  # Seurat v5 path -- present iff SeuratObject >= 5
  if (exists("LayerData", envir = asNamespace("SeuratObject"), inherits = FALSE)) {
    out <- fetch(SeuratObject::LayerData, assay = assay, layer = layer)
    if (!is.null(out)) return(out)
  }
  # Fallback v4 path
  fetch(SeuratObject::GetAssayData, assay = assay, slot = layer)
}

# ---- 10x Genomics directory -----------------------------------------------

#' Load a 10x Genomics filtered/raw feature-barcode matrix directory.
#'
#' Looks for the canonical Cellranger layout:
#'   <dir>/matrix.mtx{,.gz}
#'   <dir>/barcodes.tsv{,.gz}
#'   <dir>/features.tsv{,.gz}    (modern, columns: id, symbol, feature_type)
#'   <dir>/genes.tsv{,.gz}       (legacy, columns:  id, symbol)
#'
#' Returns a dataset with `reductions = character()` and
#' `default_reduction = NA_character_` -- raw 10x output carries no
#' dimensional reductions. Modules with reduction as a required input show
#' their empty-state. Run PCA + UMAP downstream (future), or load through
#' a Seurat / AnnData pipeline that already has embeddings.
load_10x <- function(path, gene_column = 2L, unique_features = TRUE) {
  require_optional("Matrix", feature = "10x Genomics .mtx loading")
  if (!dir.exists(path)) {
    stop("10x loader: directory does not exist: ", path, call. = FALSE)
  }
  paths <- .resolve_10x_paths(path)

  mat <- Matrix::readMM(paths$matrix)
  # readMM returns a dgTMatrix; convert to the canonical dgCMatrix so
  # row-indexing for per-gene reads is fast and stable.
  mat <- methods::as(mat, "CsparseMatrix")

  barcodes <- .read_tsv_lines(paths$barcodes)
  features <- .read_tsv_table(paths$features)
  if (ncol(features) < gene_column) {
    stop("10x loader: features table has ", ncol(features),
         " column(s) but gene_column = ", gene_column, ".", call. = FALSE)
  }
  gene_names <- as.character(features[[gene_column]])
  if (unique_features) gene_names <- make.unique(gene_names)

  if (nrow(mat) != length(gene_names)) {
    stop("10x loader: matrix has ", nrow(mat), " rows but features file lists ",
         length(gene_names), " genes.", call. = FALSE)
  }
  if (ncol(mat) != length(barcodes)) {
    stop("10x loader: matrix has ", ncol(mat),
         " columns but barcodes file lists ", length(barcodes), " cells.",
         call. = FALSE)
  }
  rownames(mat) <- gene_names
  colnames(mat) <- barcodes

  expr_be <- expression_backend_sparse(list(counts = mat),
                                       n_cells = length(barcodes),
                                       default_layer = "counts")

  cell_data <- data.frame(cell = barcodes,
                          n_counts = Matrix::colSums(mat),
                          n_features = Matrix::colSums(mat > 0),
                          stringsAsFactors = FALSE)

  list(
    name              = .basename_noext(path),
    source            = "10x",
    n_cells           = as.integer(length(barcodes)),
    n_genes           = as.integer(length(gene_names)),
    assays            = "RNA",
    default_assay     = "RNA",
    reductions        = character(),
    default_reduction = NA_character_,
    metadata_fields   = c("n_counts", "n_features"),
    cells             = barcodes,
    cell_data         = cell_data,
    genes             = gene_names,
    expression        = expr_be
  )
}

.resolve_10x_paths <- function(dir) {
  pick <- function(stems) {
    for (stem in stems) {
      cand <- c(file.path(dir, stem), file.path(dir, paste0(stem, ".gz")))
      hit  <- cand[file.exists(cand)]
      if (length(hit) > 0L) return(hit[1])
    }
    NULL
  }
  mat <- pick("matrix.mtx")
  bc  <- pick("barcodes.tsv")
  ft  <- pick(c("features.tsv", "genes.tsv"))
  missing <- c(matrix = is.null(mat),
               barcodes = is.null(bc),
               features = is.null(ft))
  if (any(missing)) {
    stop("10x loader: required file(s) missing from ", dir, ": ",
         paste(names(missing)[missing], collapse = ", "),
         ". Expected the Cellranger filtered/raw_feature_bc_matrix layout.",
         call. = FALSE)
  }
  list(matrix = mat, barcodes = bc, features = ft)
}

.read_tsv_lines <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path) else file(path)
  on.exit(close(con), add = TRUE)
  readLines(con)
}

.read_tsv_table <- function(path) {
  utils::read.delim(path, header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "", comment.char = "")
}

# ---- AnnData .h5ad -------------------------------------------------------

#' Load an AnnData (.h5ad) file as a dataset.
#'
#' Three implementations, picked in this order (first available wins):
#'
#'   1. `lazy = TRUE` (default) and `rhdf5` available -> preferred path.
#'      Reads `obs` / `var` / `obsm` eagerly but builds a lazy
#'      `expression_backend_h5ad()` so expression is only pulled
#'      gene-by-gene. Best for large files (100k+ cells); avoids
#'      materialising the full expression matrix.
#'   2. `zellkonverter` + `SingleCellExperiment` available -> eager
#'      fallback. Slurps the file via `readH5AD()` into an SCE; the
#'      expression matrix is materialised in memory.
#'   3. `anndata` available -> eager fallback. Uses the CRAN `anndata`
#'      package via `read_h5ad()`; also materialises the matrix.
#'
#' `obs` becomes `cell_data`, `var_names` becomes `genes`, and `obsm`
#' reductions (with `X_` prefixes stripped) become embedding columns.
#'
#' @param path  character(1) path to `.h5ad`
#' @param lazy  logical(1); if TRUE (default) prefer the lazy `rhdf5`
#'              reader when available. Pass `FALSE` to force one of
#'              the eager `zellkonverter` / `anndata` paths.
load_anndata <- function(path, lazy = TRUE) {
  if (!file.exists(path)) {
    stop("AnnData loader: file does not exist: ", path, call. = FALSE)
  }
  if (isTRUE(lazy) && has_optional("rhdf5")) {
    return(.load_h5ad_lazy(path))
  }
  if (has_optional(c("zellkonverter", "SingleCellExperiment"))) {
    return(.load_h5ad_via_zellkonverter(path))
  }
  if (has_optional("anndata")) {
    return(.load_h5ad_via_anndata(path))
  }
  stop(paste0(
    "AnnData (.h5ad) loading needs one of:\n",
    "  Bioconductor: BiocManager::install('rhdf5')          # recommended; lazy reads\n",
    "  Bioconductor: BiocManager::install(c('zellkonverter','SingleCellExperiment'))\n",
    "  CRAN:         install.packages('anndata')\n",
    "None of them were found on the library path."), call. = FALSE)
}

# Lazy reader via rhdf5. Builds an `expression_backend_h5ad()` so
# per-gene reads stay O(nnz) without ever materialising the full
# expression matrix. obs / var / obsm are still read eagerly (they're
# the metadata layer; small).
.load_h5ad_lazy <- function(path) {
  require_optional("rhdf5",
                   feature = "AnnData (.h5ad) lazy loading",
                   source  = "Bioconductor")

  cells <- as.character(rhdf5::h5read(path, "/obs/_index"))
  genes <- as.character(rhdf5::h5read(path, "/var/_index"))

  cell_data <- .h5ad_read_obs(path, cells)
  cell_data <- .h5ad_attach_obsm(cell_data, path)
  red_pretty <- .h5ad_obsm_reduction_names(path)

  layers <- list()
  x_spec <- .h5ad_layer_spec(path, "/X", genes,
                             n_obs = length(cells),
                             n_var = length(genes))
  layers$data <- x_spec
  if (.h5ad_group_exists(path, "/layers")) {
    sub <- .h5ad_list_layers(path)
    for (nm in sub) {
      layers[[nm]] <- .h5ad_layer_spec(path, paste0("/layers/", nm),
                                       genes,
                                       n_obs = length(cells),
                                       n_var = length(genes))
    }
  }
  expr_be <- expression_backend_h5ad(
    path          = path,
    layers        = layers,
    n_cells       = length(cells),
    default_layer = "data")

  meta_fields <- setdiff(names(cell_data),
                         c("cell",
                           unlist(lapply(red_pretty, function(r)
                             paste0(r, c("_1", "_2"))))))
  list(
    name              = .basename_noext(path),
    source            = "anndata",
    n_cells           = as.integer(length(cells)),
    n_genes           = as.integer(length(genes)),
    assays            = "RNA",
    default_assay     = "RNA",
    reductions        = unname(red_pretty),
    default_reduction = if (length(red_pretty) > 0L) red_pretty[1] else NA_character_,
    metadata_fields   = meta_fields,
    cells             = cells,
    cell_data         = cell_data,
    genes             = genes,
    expression        = expr_be
  )
}

# Read /obs as a data.frame. AnnData v0.8+ stores each column as either
# a dataset (numeric / string / bool) or as a "categorical" group with
# /codes (integer) + /categories (string). Categoricals are decoded
# back to character.
.h5ad_read_obs <- function(path, cells) {
  cols <- .h5ad_ls_group(path, "/obs")
  # _index isn't user data
  cols <- setdiff(cols, "_index")
  out  <- data.frame(cell = cells, stringsAsFactors = FALSE)
  for (col in cols) {
    val <- .h5ad_read_obs_column(path, paste0("/obs/", col))
    if (is.null(val)) next
    if (length(val) != length(cells)) next
    out[[col]] <- val
  }
  out
}

# One obs/var column. Returns NULL if shape doesn't match the expected
# 1D vector (rare; e.g. weird custom structs we don't decode here).
.h5ad_read_obs_column <- function(path, dataset_path) {
  if (.h5ad_is_categorical(path, dataset_path)) {
    codes      <- as.integer(rhdf5::h5read(path,
                                           paste0(dataset_path, "/codes")))
    categories <- as.character(rhdf5::h5read(path,
                                             paste0(dataset_path, "/categories")))
    # AnnData uses -1 for missing. Codes are 0-based.
    out <- rep(NA_character_, length(codes))
    ok  <- codes >= 0L
    out[ok] <- categories[codes[ok] + 1L]
    return(out)
  }
  v <- tryCatch(rhdf5::h5read(path, dataset_path),
                error = function(e) NULL)
  if (is.null(v)) return(NULL)
  drop_attrs <- function(x) { attributes(x) <- NULL; x }
  drop_attrs(v)
}

# Read /obsm/X_* and append 2D embeddings as <name>_1 / <name>_2.
.h5ad_attach_obsm <- function(cell_data, path) {
  if (!.h5ad_group_exists(path, "/obsm")) return(cell_data)
  names <- .h5ad_ls_group(path, "/obsm")
  for (nm in names) {
    emb <- tryCatch(rhdf5::h5read(path, paste0("/obsm/", nm)),
                    error = function(e) NULL)
    if (is.null(emb) || length(dim(emb)) != 2L) next
    pretty <- .strip_X_prefix(nm)
    # AnnData stores obsm/<name> with HDF5 shape (n_obs, d). rhdf5's
    # default `native = FALSE` mode preserves HDF5 dim order, so we
    # *usually* see (n_obs, d) here -- but that depends on rhdf5
    # build / setting, and we don't control either at the call site.
    # The heuristic below handles whichever layout rhdf5 returns by
    # locking onto the dim whose length matches n_obs. This is
    # robust regardless of the rhdf5 `native` flag.
    if (nrow(emb) == nrow(cell_data)) {
      # already (n_obs, d)
    } else if (ncol(emb) == nrow(cell_data)) {
      emb <- t(emb)
    } else next
    # Edge case: when d == n_obs the first branch fires (because of the
    # check order) and assumes (n_obs, d); this is consistent with the
    # AnnData spec for the rhdf5-default layout, but pathological for
    # square embeddings under `native = TRUE`. Such embeddings are
    # vanishingly rare for scRNA-seq (n_obs is typically 1e3-1e6, d is
    # 2-50) and we accept this asymmetry rather than introduce a flaky
    # spec-attribute probe.
    if (ncol(emb) < 2L) next
    cell_data[[paste0(pretty, "_1")]] <- emb[, 1]
    cell_data[[paste0(pretty, "_2")]] <- emb[, 2]
  }
  cell_data
}

.h5ad_obsm_reduction_names <- function(path) {
  if (!.h5ad_group_exists(path, "/obsm")) return(character())
  names <- .h5ad_ls_group(path, "/obsm")
  unname(vapply(names, .strip_X_prefix, character(1)))
}

# Build a layer spec list ({path, encoding, shape, genes}) by reading
# the encoding-type attribute. AnnData < v0.8 may omit it; we
# heuristically default to "dense" in that case.
.h5ad_layer_spec <- function(path, ds_path, genes, n_obs, n_var) {
  attrs <- tryCatch(rhdf5::h5readAttributes(path, ds_path),
                    error = function(e) list())
  enc <- attrs[["encoding-type"]]
  if (is.null(enc)) {
    enc <- "dense"
  }
  if (enc == "array") enc <- "dense"
  if (!(enc %in% H5AD_SUPPORTED_ENCODINGS)) {
    stop(sprintf(
      "AnnData lazy loader: layer '%s' has encoding-type '%s' which is not in {%s}",
      ds_path, enc,
      paste(H5AD_SUPPORTED_ENCODINGS, collapse = ", ")),
      call. = FALSE)
  }
  list(path     = ds_path,
       encoding = enc,
       shape    = c(n_obs, n_var),
       genes    = genes)
}

# ---- small rhdf5 helpers --------------------------------------------------
.h5ad_group_exists <- function(path, group) {
  tryCatch({
    ls <- rhdf5::h5ls(path, recursive = FALSE)
    base <- sub("^/", "", group)
    base %in% ls$name
  }, error = function(e) FALSE)
}

.h5ad_ls_group <- function(path, group) {
  ls <- rhdf5::h5ls(path, recursive = FALSE, datasetinfo = FALSE)
  ls$name[ls$group == group]
}

.h5ad_list_layers <- function(path) .h5ad_ls_group(path, "/layers")

.h5ad_is_categorical <- function(path, ds_path) {
  attrs <- tryCatch(rhdf5::h5readAttributes(path, ds_path),
                    error = function(e) list())
  identical(attrs[["encoding-type"]], "categorical")
}

.load_h5ad_via_zellkonverter <- function(path) {
  require_optional(c("zellkonverter", "SingleCellExperiment"),
                   feature = "AnnData (.h5ad) loading",
                   source = "Bioconductor")
  sce <- zellkonverter::readH5AD(path)
  .sce_to_dataset(sce, name = .basename_noext(path))
}

.sce_to_dataset <- function(sce, name = "anndata") {
  require_optional("SingleCellExperiment",
                   feature = "AnnData -> dataset mapping",
                   source = "Bioconductor")
  cells     <- colnames(sce)
  genes     <- rownames(sce)
  cell_data <- as.data.frame(SingleCellExperiment::colData(sce),
                             stringsAsFactors = FALSE)
  if (!"cell" %in% names(cell_data)) cell_data$cell <- cells

  red_names <- SingleCellExperiment::reducedDimNames(sce)
  red_pretty <- vapply(red_names, .strip_X_prefix, character(1))
  for (i in seq_along(red_names)) {
    emb <- SingleCellExperiment::reducedDim(sce, red_names[i])
    if (ncol(emb) >= 2L) {
      cell_data[[paste0(red_pretty[i], "_1")]] <- emb[, 1]
      cell_data[[paste0(red_pretty[i], "_2")]] <- emb[, 2]
    }
  }
  meta_fields <- setdiff(names(cell_data),
                         c("cell",
                           unlist(lapply(red_pretty, function(r)
                             paste0(r, c("_1", "_2"))))))

  # zellkonverter exposes "X" as the primary assay; preserve any additional
  # layers (e.g. "counts").
  assay_names <- SummarizedExperiment::assayNames(sce)
  primary <- if ("X" %in% assay_names) "X" else assay_names[1]
  layers  <- stats::setNames(
    lapply(assay_names, function(an) {
      m <- SummarizedExperiment::assay(sce, an)
      if (is.null(rownames(m))) rownames(m) <- genes
      if (is.null(colnames(m))) colnames(m) <- cells
      m
    }),
    ifelse(assay_names == "X", "data", assay_names))
  default_layer <- if ("X" %in% assay_names) "data" else assay_names[1]
  expr_be <- expression_backend_sparse(layers, n_cells = length(cells),
                                       default_layer = default_layer)

  list(
    name              = name,
    source            = "anndata",
    n_cells           = as.integer(length(cells)),
    n_genes           = as.integer(length(genes)),
    assays            = "RNA",
    default_assay     = "RNA",
    reductions        = unname(red_pretty),
    default_reduction = if (length(red_pretty) > 0L) red_pretty[1] else NA_character_,
    metadata_fields   = meta_fields,
    cells             = cells,
    cell_data         = cell_data,
    genes             = genes,
    expression        = expr_be
  )
}

.load_h5ad_via_anndata <- function(path) {
  require_optional("anndata", feature = "AnnData (.h5ad) loading")
  ad        <- anndata::read_h5ad(path)
  cells     <- as.character(ad$obs_names)
  genes     <- as.character(ad$var_names)
  cell_data <- as.data.frame(ad$obs, stringsAsFactors = FALSE)
  if (!"cell" %in% names(cell_data)) cell_data$cell <- cells

  obsm <- ad$obsm
  red_names <- if (is.null(obsm)) character() else names(obsm)
  red_pretty <- vapply(red_names, .strip_X_prefix, character(1))
  for (i in seq_along(red_names)) {
    emb <- obsm[[red_names[i]]]
    if (ncol(emb) >= 2L) {
      cell_data[[paste0(red_pretty[i], "_1")]] <- emb[, 1]
      cell_data[[paste0(red_pretty[i], "_2")]] <- emb[, 2]
    }
  }
  meta_fields <- setdiff(names(cell_data),
                         c("cell",
                           unlist(lapply(red_pretty, function(r)
                             paste0(r, c("_1", "_2"))))))

  X <- ad$X
  # anndata$X is cells x genes; transpose so rows are genes (canonical
  # orientation for the rest of the app).
  X <- Matrix::t(X)
  rownames(X) <- genes
  colnames(X) <- cells
  expr_be <- expression_backend_sparse(list(data = X), n_cells = length(cells))

  list(
    name              = .basename_noext(path),
    source            = "anndata",
    n_cells           = as.integer(length(cells)),
    n_genes           = as.integer(length(genes)),
    assays            = "RNA",
    default_assay     = "RNA",
    reductions        = unname(red_pretty),
    default_reduction = if (length(red_pretty) > 0L) red_pretty[1] else NA_character_,
    metadata_fields   = meta_fields,
    cells             = cells,
    cell_data         = cell_data,
    genes             = genes,
    expression        = expr_be
  )
}

# ---- Tiny shared helpers -------------------------------------------------

.basename_noext <- function(path) {
  bn <- basename(path)
  sub("\\.(rds|h5ad|tar\\.gz)$", "", bn, ignore.case = TRUE)
}

.strip_X_prefix <- function(x) sub("^X_", "", x, ignore.case = FALSE)

# ---- Back-compat aliases -------------------------------------------------
# Earlier versions exposed `*_stub()` names. Keep the names so any external
# scripts or tests that still call them get the real loaders.
load_seurat_stub  <- load_seurat
load_anndata_stub <- load_anndata
load_10x_stub     <- load_10x
