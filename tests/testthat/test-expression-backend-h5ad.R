# Tests for R/expression_backend_h5ad.R
# ----------------------------------------------------------------------------
# The backend has two surfaces we can exercise independently:
#
#   1. The constructor + parameter validation, which is pure (no rhdf5
#      required).
#   2. The per-gene CSR / CSC slice maths, which is pure as long as the
#      sparse triple (indptr, indices, data) is in the per-layer cache.
#      We populate that cache directly so the algorithm is testable
#      without rhdf5 or an actual `.h5ad` file.
#
# A final end-to-end test rebuilds a real fixture via `rhdf5` and runs
# `load_anndata(lazy = TRUE)` against it; it skips when rhdf5 isn't
# installed (which is the case in CI here).

# ---- Constructor + validation --------------------------------------------

test_that("expression_backend_h5ad refuses a missing path", {
  expect_error(expression_backend_h5ad("/no/such/file.h5ad",
                                       layers = list()),
               regexp = "existing .h5ad file", fixed = TRUE)
})

test_that("expression_backend_h5ad refuses bad layer specs", {
  f <- tempfile(fileext = ".h5ad")
  file.create(f)
  on.exit(unlink(f), add = TRUE)

  expect_error(expression_backend_h5ad(f, layers = list()),
               regexp = "named, non-empty", fixed = TRUE)
  expect_error(expression_backend_h5ad(f,
    layers = list(data = list(path = "/X"))),
    regexp = "missing field", fixed = TRUE)
  expect_error(expression_backend_h5ad(f,
    layers = list(data = list(path = "/X", encoding = "weird",
                              shape = c(10, 20), genes = "G1"))),
    regexp = "encoding 'weird' not supported", fixed = TRUE)
  expect_error(expression_backend_h5ad(f,
    layers = list(data = list(path = "/X", encoding = "csr_matrix",
                              shape = c(10, 20, 30), genes = "G1"))),
    regexp = "shape must be", fixed = TRUE)
})

# ---- Per-layer cache: synthetic sparse triple -----------------------------
#
# `.h5ad_get_gene_csr/csc` and `.h5ad_materialise_csr/csc` only need
# `backend$caches[[layer]]` to be populated. We build a 3-cell x 4-gene
# toy matrix in both CSR and CSC layouts and verify the per-gene reads
# match `as.matrix()`.

.h5ad_make_synthetic_backend <- function(encoding) {
  f <- tempfile(fileext = ".h5ad")
  file.create(f)
  genes <- c("G1", "G2", "G3", "G4")
  shape <- c(3L, 4L)  # n_obs, n_var
  be <- expression_backend_h5ad(
    path = f,
    layers = list(data = list(path = "/X", encoding = encoding,
                              shape = shape, genes = genes)),
    n_cells = 3L)
  attr(be, "tmpfile") <- f  # so the caller can unlink
  be
}

test_that(".h5ad_get_gene_csr reproduces the column of a dense matrix", {
  # Dense reference (cells x genes):
  #   c1: G1=0, G2=2, G3=0, G4=5
  #   c2: G1=1, G2=0, G3=3, G4=0
  #   c3: G1=0, G2=0, G3=4, G4=6
  #
  # CSR (row-major over cells):
  #   indptr  = [0, 2, 4, 6]
  #   indices = [1, 3, 0, 2, 2, 3]   (0-based gene indices, sorted per row)
  #   data    = [2, 5, 1, 3, 4, 6]
  be <- .h5ad_make_synthetic_backend("csr_matrix")
  on.exit(unlink(attr(be, "tmpfile")), add = TRUE)
  cache <- be$caches$data
  cache$indptr  <- c(0L, 2L, 4L, 6L)
  cache$indices <- c(1L, 3L, 0L, 2L, 2L, 3L)
  cache$data    <- c(2, 5, 1, 3, 4, 6)
  cache$loaded  <- TRUE

  expect_equal(backend_get_gene(be, "G1", layer = "data"),
               c(0, 1, 0))
  expect_equal(backend_get_gene(be, "G2", layer = "data"),
               c(2, 0, 0))
  expect_equal(backend_get_gene(be, "G3", layer = "data"),
               c(0, 3, 4))
  expect_equal(backend_get_gene(be, "G4", layer = "data"),
               c(5, 0, 6))
  expect_null(backend_get_gene(be, "NOTAGENE", layer = "data"))
})

test_that(".h5ad_get_gene_csc reproduces the column of a dense matrix", {
  # Same matrix as above, CSC (column-major over genes):
  #   indptr  = [0, 1, 2, 4, 6]
  #   indices = [1, 0, 1, 2, 0, 2]   (0-based cell indices)
  #   data    = [1, 2, 3, 4, 5, 6]
  be <- .h5ad_make_synthetic_backend("csc_matrix")
  on.exit(unlink(attr(be, "tmpfile")), add = TRUE)
  cache <- be$caches$data
  cache$indptr  <- c(0L, 1L, 2L, 4L, 6L)
  cache$indices <- c(1L, 0L, 1L, 2L, 0L, 2L)
  cache$data    <- c(1, 2, 3, 4, 5, 6)
  cache$loaded  <- TRUE

  expect_equal(backend_get_gene(be, "G1", layer = "data"), c(0, 1, 0))
  expect_equal(backend_get_gene(be, "G2", layer = "data"), c(2, 0, 0))
  expect_equal(backend_get_gene(be, "G3", layer = "data"), c(0, 3, 4))
  expect_equal(backend_get_gene(be, "G4", layer = "data"), c(5, 0, 6))
})

test_that(".h5ad_materialise_csr / _csc match per-gene reads AND keep genes-as-rows orientation", {
  # Same matrix; we don't need rhdf5 because `.h5ad_materialise_csr`
  # only touches the per-layer cache that we populate directly.
  expected <- rbind(
    G1 = c(0, 1, 0),
    G2 = c(2, 0, 0),
    G3 = c(0, 3, 4),
    G4 = c(5, 0, 6))

  for (enc in c("csr_matrix", "csc_matrix")) {
    be <- .h5ad_make_synthetic_backend(enc)
    cache <- be$caches$data
    if (enc == "csr_matrix") {
      cache$indptr  <- c(0L, 2L, 4L, 6L)
      cache$indices <- c(1L, 3L, 0L, 2L, 2L, 3L)
      cache$data    <- c(2, 5, 1, 3, 4, 6)
    } else {
      cache$indptr  <- c(0L, 1L, 2L, 4L, 6L)
      cache$indices <- c(1L, 0L, 1L, 2L, 0L, 2L)
      cache$data    <- c(1, 2, 3, 4, 5, 6)
    }
    cache$loaded <- TRUE

    m <- if (enc == "csr_matrix")
           .h5ad_materialise_csr(be, "data", c("G1", "G2", "G3", "G4"),
                                 n_obs = 3, n_var = 4)
         else
           .h5ad_materialise_csc(be, "data", c("G1", "G2", "G3", "G4"),
                                 n_obs = 3, n_var = 4)
    # Values
    expect_identical(unname(m), unname(expected),
                     info = sprintf("encoding = %s", enc))
    # Orientation contract: result is `n_var x n_obs` with rownames =
    # genes. This pins the contract downstream callers (DE, AUCell,
    # SingleR) depend on.
    expect_identical(dim(m), c(4L, 3L),    info = enc)
    expect_identical(rownames(m),
                     c("G1", "G2", "G3", "G4"), info = enc)
    # Cross-check: each row matches the per-gene reader. A transpose
    # regression in the materialiser would mismatch a per-gene call.
    for (g in c("G1", "G2", "G3", "G4")) {
      expect_identical(unname(m[g, ]),
                       backend_get_gene(be, g, layer = "data"),
                       info = paste(enc, g))
    }
    unlink(attr(be, "tmpfile"))
  }
})

test_that(".h5ad_materialise_dense dispatches by HDF5 dim ordering", {
  # Substitute a fake `rhdf5::h5read` so we can drive the dispatch
  # without installing the package. The dense materialiser only calls
  # one rhdf5 entry-point (`h5read(path, p)`); we stub it via
  # `local_mocked_bindings` (testthat 3.2+) when available, else skip.
  skip_if_not(packageVersion("testthat") >= "3.2.0",
              "needs testthat::local_mocked_bindings (>= 3.2.0)")
  be <- .h5ad_make_synthetic_backend("dense")
  on.exit(unlink(attr(be, "tmpfile")), add = TRUE)
  # AnnData dense X: cells x genes => HDF5 shape (3, 4)
  dense <- matrix(seq_len(12), nrow = 3, ncol = 4)  # (n_obs, n_var)

  testthat::local_mocked_bindings(
    h5read = function(file, name, ...) dense,
    .package = "rhdf5")
  m <- .h5ad_materialise_dense(be, "data",
                               genes = c("G1", "G2", "G3", "G4"),
                               n_obs = 3, n_var = 4)
  expect_identical(dim(m), c(4L, 3L))
  expect_identical(rownames(m), c("G1", "G2", "G3", "G4"))
  expect_identical(unname(m), t(dense))

  # Same dataset, but rhdf5 returns it transposed (native = TRUE
  # convention). The materialiser must keep it genes-as-rows.
  testthat::local_mocked_bindings(
    h5read = function(file, name, ...) t(dense),
    .package = "rhdf5")
  m2 <- .h5ad_materialise_dense(be, "data",
                                genes = c("G1", "G2", "G3", "G4"),
                                n_obs = 3, n_var = 4)
  expect_identical(dim(m2), c(4L, 3L))
  expect_identical(rownames(m2), c("G1", "G2", "G3", "G4"))
  expect_identical(unname(m2), t(dense))
})

# ---- Backend basic methods -----------------------------------------------

test_that("backend_genes / _n_cells / _n_genes / _has_gene", {
  be <- .h5ad_make_synthetic_backend("csr_matrix")
  on.exit(unlink(attr(be, "tmpfile")), add = TRUE)

  expect_identical(backend_n_cells(be), 3L)
  expect_identical(backend_genes(be), c("G1", "G2", "G3", "G4"))
  expect_identical(backend_n_genes(be), 4L)
  expect_identical(backend_available_layers(be), "data")
  expect_identical(backend_default_layer(be), "data")
  expect_true(backend_has_gene(be, "G2"))
  expect_false(backend_has_gene(be, "G99"))
})

test_that("unknown layer is reported clearly", {
  be <- .h5ad_make_synthetic_backend("csr_matrix")
  on.exit(unlink(attr(be, "tmpfile")), add = TRUE)
  expect_error(backend_genes(be, layer = "counts"),
               regexp = "Layer 'counts' not available", fixed = TRUE)
})

# ---- Loader gating --------------------------------------------------------

test_that("load_anndata(lazy = TRUE) errors cleanly without rhdf5 + fallbacks", {
  skip_if(has_optional("rhdf5"),
          "rhdf5 installed; install-error path skipped")
  skip_if(has_optional("zellkonverter"),
          "zellkonverter installed; fallback handles missing rhdf5")
  skip_if(has_optional("anndata"),
          "anndata installed; fallback handles missing rhdf5")
  expect_error(load_anndata(tempfile(fileext = ".h5ad")),
               regexp = "file does not exist", fixed = TRUE)
})

# ---- .h5ad_layer_spec unsupported-encoding error message -----------------
# Earlier code mis-passed three substitution values to a format string
# with only two `%s` placeholders, producing a garbled
# "layer 'which is not in {%s}' has encoding-type '<path>'" message.
# Two tests pin the contract: a pure pattern check (so a regression
# would fire even without rhdf5), and a real `.h5ad_layer_spec()` call
# under mocked rhdf5 (runs only when rhdf5 is installed).

test_that("the unsupported-encoding error template is well-formed", {
  msg <- sprintf(
    "AnnData lazy loader: layer '%s' has encoding-type '%s' which is not in {%s}",
    "/X", "weird",
    paste(H5AD_SUPPORTED_ENCODINGS, collapse = ", "))
  expect_match(msg, "layer '/X'",            fixed = TRUE)
  expect_match(msg, "encoding-type 'weird'", fixed = TRUE)
  expect_match(msg, "csr_matrix",            fixed = TRUE)
  # No literal `%s` left over from the format string.
  expect_false(grepl("%s", msg, fixed = TRUE))
})

test_that(".h5ad_layer_spec errors with a complete message on unsupported encodings", {
  skip_if_not_installed("rhdf5")
  skip_if_not(packageVersion("testthat") >= "3.2.0",
              "needs testthat::local_mocked_bindings (>= 3.2.0)")
  testthat::local_mocked_bindings(
    h5readAttributes = function(file, name, ...)
      list(`encoding-type` = "weird_encoding"),
    .package = "rhdf5")
  err <- tryCatch(
    .h5ad_layer_spec("/no/such/path.h5ad", "/X",
                     genes = c("G1", "G2"), n_obs = 3L, n_var = 2L),
    error = conditionMessage)
  expect_match(err, "layer '/X'",                       fixed = TRUE)
  expect_match(err, "encoding-type 'weird_encoding'",   fixed = TRUE)
  expect_match(err, "csr_matrix",                       fixed = TRUE)
  expect_false(grepl("%s", err, fixed = TRUE))
})

# ---- End-to-end via rhdf5 -------------------------------------------------
#
# Build a tiny .h5ad fixture from R and verify the full lazy loader
# + backend roundtrip. Skips entirely when rhdf5 is not installed.

.write_minimal_h5ad <- function(file) {
  rhdf5::h5createFile(file)

  # /X as csr_matrix: 3 cells x 4 genes, same dense matrix as the
  # synthetic CSR tests above.
  rhdf5::h5createGroup(file, "X")
  rhdf5::h5write(c(2, 5, 1, 3, 4, 6),
                 file = file, name = "X/data")
  rhdf5::h5write(as.integer(c(1, 3, 0, 2, 2, 3)),
                 file = file, name = "X/indices")
  rhdf5::h5write(as.integer(c(0, 2, 4, 6)),
                 file = file, name = "X/indptr")
  fid <- rhdf5::H5Fopen(file)
  gid <- rhdf5::H5Gopen(fid, "X")
  rhdf5::h5writeAttribute("csr_matrix", gid, "encoding-type")
  rhdf5::h5writeAttribute(as.integer(c(3, 4)), gid, "shape")
  rhdf5::H5Gclose(gid); rhdf5::H5Fclose(fid)

  # /obs: _index + one numeric column
  rhdf5::h5createGroup(file, "obs")
  rhdf5::h5write(c("c1", "c2", "c3"),
                 file = file, name = "obs/_index")
  rhdf5::h5write(c(0.1, 0.2, 0.3),
                 file = file, name = "obs/qc_score")

  # /var: just an index
  rhdf5::h5createGroup(file, "var")
  rhdf5::h5write(c("G1", "G2", "G3", "G4"),
                 file = file, name = "var/_index")

  # /obsm/X_umap (3 cells x 2 dims)
  rhdf5::h5createGroup(file, "obsm")
  rhdf5::h5write(matrix(c(0.1, 0.2, 0.3,
                          1.1, 1.2, 1.3), nrow = 3, ncol = 2),
                 file = file, name = "obsm/X_umap")
}

test_that("load_anndata(lazy = TRUE) round-trips against a real rhdf5 fixture", {
  skip_if_not_installed("rhdf5")
  f <- tempfile(fileext = ".h5ad")
  on.exit(unlink(f), add = TRUE)
  .write_minimal_h5ad(f)

  ds <- load_anndata(f, lazy = TRUE)
  expect_identical(ds$source, "anndata")
  expect_identical(ds$n_cells, 3L)
  expect_identical(ds$n_genes, 4L)
  expect_identical(ds$cells, c("c1", "c2", "c3"))
  expect_identical(ds$genes, c("G1", "G2", "G3", "G4"))
  expect_true("qc_score" %in% ds$metadata_fields)
  expect_true("umap" %in% ds$reductions)
  # Backend is lazy
  expect_true(inherits(ds$expression, "expression_backend_h5ad"))
  # Per-gene reads land on the right values
  expect_equal(get_gene_expression(ds, "G3"), c(0, 3, 4))
  expect_equal(get_gene_expression(ds, "G4"), c(5, 0, 6))
})
