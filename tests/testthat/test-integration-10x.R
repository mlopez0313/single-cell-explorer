# ============================================================================
# End-to-end integration test for load_10x() against a real Cellranger-
# layout directory built from scratch (Matrix Market + gzipped tsv).
#
# Skipped when `Matrix` isn't installed (it's an Imports, so this skip
# should effectively never fire -- but we keep the guard so the test
# stays defensive on stripped-down CI runners).
# ============================================================================

skip_if_no_matrix <- function() {
  testthat::skip_if_not_installed("Matrix")
}

.write_10x_fixture <- function(dir, n_cells = 80L, n_genes = 30L,
                               seed = 11L, gzip = TRUE) {
  set.seed(seed)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  cells <- sprintf("AAACGAA%03d-1", seq_len(n_cells))
  symbols <- c("CD3D", "MS4A1", "LST1", "EPCAM", "COL1A1", "NKG7",
               sprintf("GENE%03d", seq_len(n_genes - 6L)))
  ids <- sprintf("ENSG%011d", seq_along(symbols))
  ct <- sample(c("T", "B", "Myeloid"), n_cells, replace = TRUE)
  counts <- matrix(0L, nrow = n_genes, ncol = n_cells,
                   dimnames = list(symbols, cells))
  prof <- list(T = c(CD3D = 5, NKG7 = 2),
               B = c(MS4A1 = 5),
               Myeloid = c(LST1 = 5))
  for (j in seq_len(n_cells)) {
    p <- prof[[ct[j]]]
    for (g in names(p)) counts[g, j] <- rpois(1, p[g])
    counts[7:n_genes, j] <- rpois(n_genes - 6L, 0.1)
  }
  sp <- methods::as(counts, "CsparseMatrix")

  matrix_path   <- file.path(dir, "matrix.mtx")
  barcodes_path <- file.path(dir, "barcodes.tsv")
  features_path <- file.path(dir, "features.tsv")
  Matrix::writeMM(sp, matrix_path)
  writeLines(cells, barcodes_path)
  writeLines(paste(ids, symbols, "Gene Expression", sep = "\t"),
             features_path)
  if (isTRUE(gzip)) {
    .gzfile_inplace <- function(path) {
      gz <- paste0(path, ".gz")
      con_in <- file(path, "rb"); on.exit(close(con_in), add = TRUE)
      con_out <- gzfile(gz, "wb"); on.exit(close(con_out), add = TRUE)
      writeBin(readBin(con_in, raw(), file.info(path)$size), con_out)
      unlink(path)
    }
    .gzfile_inplace(matrix_path)
    .gzfile_inplace(barcodes_path)
    .gzfile_inplace(features_path)
  }
  list(dir = dir, cells = cells, symbols = symbols, true_ct = ct)
}

test_that("load_10x parses a gzipped Cellranger v3 directory correctly", {
  skip_if_no_matrix()
  d  <- file.path(tempfile("tenx_"), "filtered_feature_bc_matrix")
  on.exit(unlink(dirname(d), recursive = TRUE), add = TRUE)
  fx <- .write_10x_fixture(d, gzip = TRUE)

  ds <- load_10x(d)
  expect_equal(ds$source, "10x")
  expect_equal(ds$n_cells, length(fx$cells))
  expect_equal(ds$n_genes, length(fx$symbols))
  expect_setequal(ds$genes, fx$symbols)
  expect_setequal(ds$cells, fx$cells)
  expect_equal(ds$reductions, character(0))
  expect_true(is.na(ds$default_reduction))
  expect_true(all(c("n_counts", "n_features") %in% ds$metadata_fields))
  # CD3D should be higher in T than in B on the round-tripped matrix
  e <- get_gene_expression(ds, "CD3D")
  expect_gt(mean(e[fx$true_ct == "T"]),
            mean(e[fx$true_ct == "B"]))
})

test_that("load_10x also parses the un-gzipped legacy layout", {
  skip_if_no_matrix()
  d  <- file.path(tempfile("tenx_legacy_"), "raw_feature_bc_matrix")
  on.exit(unlink(dirname(d), recursive = TRUE), add = TRUE)
  fx <- .write_10x_fixture(d, gzip = FALSE)

  ds <- load_10x(d)
  expect_equal(ds$n_cells, length(fx$cells))
  expect_setequal(ds$genes, fx$symbols)
})

test_that("load_10x errors clearly when the matrix / barcode / feature files are missing", {
  skip_if_no_matrix()
  empty <- tempfile("tenx_empty_"); dir.create(empty)
  on.exit(unlink(empty, recursive = TRUE), add = TRUE)
  expect_error(load_10x(empty), regexp = "matrix.*barcodes.*features")
})
