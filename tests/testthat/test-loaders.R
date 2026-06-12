# =========================================================================
# Real dataset loaders: detect_source, Seurat, 10x, AnnData
#
# Strategy:
#   * 10x and Seurat loaders run end-to-end against tiny synthetic inputs
#     created in tempdir() (so we don't ship test fixtures).
#   * AnnData is tested only for its graceful-missing-dep error path
#     unless both zellkonverter/anndata are installed.
#   * Every test that touches an optional package uses `skip_if_not_installed`
#     so the suite stays green on minimal R installs.
# =========================================================================

# ---- detect_source ------------------------------------------------------

test_that("detect_source infers source from file extensions", {
  expect_identical(detect_source("/tmp/foo.rds"),  "seurat")
  expect_identical(detect_source("/tmp/foo.h5ad"), "anndata")
  expect_error(detect_source("/tmp/foo.unknown"), "Cannot infer dataset source")
})

# P4: directory inference must be plausible-10x, not "any dir is 10x".
# Prior behaviour returned "10x" for *any* directory and then surfaced
# misleading low-level "matrix.mtx missing" errors deep inside the
# loader. Now `detect_source()` itself refuses with a clear message,
# and only directories that *look* like a Cellranger feature-barcode
# matrix come back as "10x".

test_that("detect_source returns '10x' for plausible Cellranger directories", {
  td <- file.path(tempdir(), "tenx_detect_valid")
  dir.create(td, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  file.create(file.path(td, "matrix.mtx"))
  file.create(file.path(td, "barcodes.tsv"))
  file.create(file.path(td, "features.tsv"))
  expect_identical(detect_source(td), "10x")
})

test_that("detect_source recognises the gzipped Cellranger v3 layout", {
  td <- file.path(tempdir(), "tenx_detect_v3")
  dir.create(td, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  file.create(file.path(td, "matrix.mtx.gz"))
  file.create(file.path(td, "barcodes.tsv.gz"))
  file.create(file.path(td, "features.tsv.gz"))
  expect_identical(detect_source(td), "10x")
})

test_that("detect_source recognises the legacy v2 'genes.tsv' name", {
  td <- file.path(tempdir(), "tenx_detect_v2")
  dir.create(td, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  file.create(file.path(td, "matrix.mtx"))
  file.create(file.path(td, "barcodes.tsv"))
  file.create(file.path(td, "genes.tsv"))   # v2 name
  expect_identical(detect_source(td), "10x")
})

test_that("detect_source refuses an arbitrary directory with a clear error", {
  td <- file.path(tempdir(), "tenx_detect_empty")
  dir.create(td, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  expect_error(detect_source(td),
               regexp = "does not look like a Cellranger")
  # Still useless when only a subset of the triad is present.
  file.create(file.path(td, "matrix.mtx"))
  expect_error(detect_source(td),
               regexp = "does not look like a Cellranger")
})

test_that("load_dataset(source = '10x', ...) bypasses detect_source", {
  td <- file.path(tempdir(), "tenx_explicit_force")
  dir.create(td, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  # An empty dir is not plausibly 10x; explicit `source = "10x"` must
  # still route to load_10x(), which then errors with its own (clearer)
  # "required file(s) missing" message rather than the detection error.
  expect_error(load_dataset(td, source = "10x"),
               regexp = "required file\\(s\\) missing")
})

# ---- Optional-deps helpers ---------------------------------------------

test_that("has_optional / require_optional behave as documented", {
  expect_true(has_optional("base"))
  expect_false(has_optional("definitelyNotARealPackage_xyz"))
  expect_error(
    require_optional("definitelyNotARealPackage_xyz", feature = "X"),
    "requires CRAN package"
  )
  expect_error(
    require_optional("definitelyNotARealPackage_xyz", feature = "X",
                     source = "Bioconductor"),
    "BiocManager::install"
  )
})

# P7: GitHub source generates a `remotes::install_github()` hint. When
# the caller supplies an explicit `repo` mapping, the placeholder
# `<owner>` is replaced; otherwise the placeholder is shown.

test_that("require_optional(source='GitHub') without `repo` shows a <owner> placeholder", {
  err <- tryCatch(
    require_optional("definitelyNotARealPackage_xyz", feature = "X",
                     source = "GitHub"),
    error = conditionMessage)
  expect_match(err, "remotes::install_github", fixed = TRUE)
  expect_match(err, "<owner>",                fixed = TRUE)
})

test_that("require_optional(source='GitHub', repo=...) uses the explicit owner/repo", {
  err <- tryCatch(
    require_optional("presto", feature = "X",
                     source = "GitHub",
                     repo   = c(presto = "immunogenomics/presto")),
    error = conditionMessage)
  # Note: this test runs unconditionally because `presto` is almost
  # certainly not installed in this env; if it ever is, we'd skip.
  skip_if(has_optional("presto"),
          "presto installed; missing-dep path is not reachable.")
  expect_match(err, "immunogenomics/presto", fixed = TRUE)
  expect_false(grepl("<owner>", err, fixed = TRUE))
})

# ---- 10x loader (Matrix-only path) -------------------------------------

# Write a tiny `matrix.mtx`, `barcodes.tsv`, `features.tsv` triad to a
# temp dir, then exercise `load_10x()` end-to-end. Mtx is plain text so
# we don't depend on Cellranger output.
.write_tiny_10x <- function(dir, compress = FALSE) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  # 4 genes x 3 cells with a handful of non-zero counts
  mtx_lines <- c(
    "%%MatrixMarket matrix coordinate integer general",
    "4 3 5",
    "1 1 5",   # gene 1, cell 1: 5
    "2 1 2",   # gene 2, cell 1: 2
    "3 2 7",   # gene 3, cell 2: 7
    "1 3 1",   # gene 1, cell 3: 1
    "4 3 9"    # gene 4, cell 3: 9
  )
  features <- c(
    "ENSG001\tCD3D\tGene Expression",
    "ENSG002\tMS4A1\tGene Expression",
    "ENSG003\tNKG7\tGene Expression",
    "ENSG004\tCD3D\tGene Expression"        # duplicate symbol on purpose
  )
  barcodes <- c("AAA-1", "CCC-1", "GGG-1")

  write_one <- function(path, lines, ext) {
    if (compress && ext != "mtx") {
      con <- gzfile(paste0(path, ".gz"), "wt"); on.exit(close(con), add = TRUE)
      writeLines(lines, con); return(paste0(path, ".gz"))
    }
    writeLines(lines, path); path
  }
  write_one(file.path(dir, "matrix.mtx"),   mtx_lines, "mtx")
  write_one(file.path(dir, "features.tsv"), features,  "tsv")
  write_one(file.path(dir, "barcodes.tsv"), barcodes,  "tsv")
  invisible(dir)
}

test_that("load_10x ingests a tiny synthetic Cellranger directory", {
  skip_if_not_installed("Matrix")
  td <- file.path(tempdir(), "tenx_basic")
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  .write_tiny_10x(td)

  ds <- load_10x(td)
  expect_identical(ds$source, "10x")
  expect_identical(ds$n_cells, 3L)
  expect_identical(ds$n_genes, 4L)
  expect_identical(ds$assays, "RNA")
  expect_identical(ds$default_assay, "RNA")
  expect_identical(ds$reductions, character())
  expect_true(is.na(ds$default_reduction))

  # Duplicate symbol resolution (CD3D + CD3D -> CD3D + CD3D.1)
  expect_setequal(ds$genes, c("CD3D", "MS4A1", "NKG7", "CD3D.1"))

  # cell_data carries the documented Cellranger-derived QC columns
  expect_true(all(c("cell", "n_counts", "n_features") %in% names(ds$cell_data)))
  expect_identical(ds$cell_data$cell, c("AAA-1", "CCC-1", "GGG-1"))

  # expression backend is sparse + responds to the helper API
  expect_s3_class(ds$expression, "expression_backend_sparse")
  expect_identical(backend_default_layer(ds$expression), "counts")
  expect_identical(get_gene_expression(ds, "CD3D"),  c(5, 0, 1))
  expect_identical(get_gene_expression(ds, "MS4A1"), c(2, 0, 0))
  expect_identical(get_gene_expression(ds, "NKG7"),  c(0, 7, 0))
})

test_that("load_10x detects gzipped feature/barcode files", {
  skip_if_not_installed("Matrix")
  td <- file.path(tempdir(), "tenx_gz")
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  .write_tiny_10x(td, compress = TRUE)

  expect_identical(load_10x(td)$n_cells, 3L)
})

test_that("load_10x errors clearly when files are missing", {
  td <- file.path(tempdir(), "tenx_empty")
  dir.create(td, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  expect_error(load_10x(td), "required file\\(s\\) missing")
})

test_that("load_10x errors when path is not a directory", {
  expect_error(load_10x(file.path(tempdir(), "_no_such_dir_")),
               "directory does not exist")
})

# ---- Seurat loader ------------------------------------------------------

test_that("load_seurat ingests a synthesized SeuratObject end-to-end", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  cells <- paste0("c", 1:6)
  genes <- c("CD3D", "MS4A1", "NKG7", "EPCAM")
  counts <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 3, 4, 4),
    j = c(1, 2, 3, 4, 5, 6),
    x = c(5L, 3L, 7L, 2L, 8L, 4L),
    dims = c(4, 6),
    dimnames = list(genes, cells))
  obj <- SeuratObject::CreateSeuratObject(counts = counts,
                                          project = "test_seurat")
  # Add a fake UMAP reduction
  emb <- matrix(stats::rnorm(12), ncol = 2,
                dimnames = list(cells, c("UMAP_1", "UMAP_2")))
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = emb, key = "UMAP_", assay = SeuratObject::DefaultAssay(obj))
  obj$cluster <- as.character(c(0, 0, 1, 1, 2, 2))

  rds <- tempfile(fileext = ".rds")
  saveRDS(obj, rds); on.exit(unlink(rds), add = TRUE)

  ds <- load_seurat(rds)
  expect_identical(ds$source, "seurat")
  expect_identical(ds$n_cells, 6L)
  expect_identical(ds$n_genes, 4L)
  expect_true("UMAP" %in% ds$reductions)
  expect_identical(ds$default_reduction, "UMAP")
  expect_true(all(c("UMAP_1", "UMAP_2") %in% names(ds$cell_data)))
  expect_true("cluster" %in% ds$metadata_fields)

  expect_s3_class(ds$expression, "expression_backend_sparse")
  expect_setequal(ds$genes, genes)
  expect_identical(get_gene_expression(ds, "CD3D"),
                   c(5, 3, 0, 0, 0, 0))
  expect_identical(get_gene_expression(ds, "NKG7"),
                   c(0, 0, 0, 2, 0, 0))
})

test_that("load_seurat rejects non-Seurat .rds payloads", {
  rds <- tempfile(fileext = ".rds")
  on.exit(unlink(rds), add = TRUE)
  saveRDS(list(not = "a seurat object"), rds)
  expect_error(load_seurat(rds), "expected 'Seurat'")
})

test_that("load_seurat surfaces a clear error when file is missing", {
  expect_error(load_seurat(file.path(tempdir(), "no_such_file.rds")),
               "file does not exist")
})

# ---- AnnData loader -----------------------------------------------------

test_that("load_anndata errors with install instructions when no backend is available", {
  if (has_optional("rhdf5") ||
      has_optional(c("zellkonverter", "SingleCellExperiment")) ||
      has_optional("anndata")) {
    skip("AnnData backend(s) are installed; missing-dep path is not reachable.")
  }
  # Write an empty file just so the file-existence check passes and we hit
  # the no-backend branch.
  tf <- tempfile(fileext = ".h5ad")
  file.create(tf); on.exit(unlink(tf), add = TRUE)
  err <- tryCatch(load_anndata(tf), error = function(e) conditionMessage(e))
  expect_match(err, "zellkonverter|anndata")
  expect_match(err, "install")
})

test_that("load_anndata round-trips when zellkonverter is installed", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")

  # Build a minimal SCE, write to .h5ad, read back via load_anndata.
  cells <- paste0("c", 1:4); genes <- c("CD3D", "NKG7")
  m <- Matrix::sparseMatrix(i = c(1, 2, 1), j = c(1, 2, 3),
                            x = c(5, 7, 2), dims = c(2, 4),
                            dimnames = list(genes, cells))
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(X = m),
    colData = S4Vectors::DataFrame(cluster = c("0", "0", "1", "1"))
  )
  emb <- matrix(stats::rnorm(8), ncol = 2, dimnames = list(cells, NULL))
  SingleCellExperiment::reducedDim(sce, "X_umap") <- emb

  h5 <- tempfile(fileext = ".h5ad")
  on.exit(unlink(h5), add = TRUE)
  zellkonverter::writeH5AD(sce, h5)

  ds <- load_anndata(h5)
  expect_identical(ds$source, "anndata")
  expect_identical(ds$n_cells, 4L)
  expect_true("umap" %in% ds$reductions || "UMAP" %in% ds$reductions)
})
