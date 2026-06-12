# =========================================================================
# Expression backend abstraction
# Mock-only backend exists today, but the contract has to hold across every
# future backend (sparse, AnnData, Seurat-proxy). These tests pin down the
# generic surface so each new backend can be added with confidence.
# =========================================================================

# ---- Flat / legacy constructor shape ------------------------------------

test_that("flat named list is lifted into a single 'data' layer", {
  expr <- list(A = c(0, 1, 2), B = c(3, 4, 5))
  be   <- expression_backend_inmemory(expr)
  expect_s3_class(be, "expression_backend")
  expect_s3_class(be, "expression_backend_inmemory")
  expect_identical(backend_n_cells(be),          3L)
  expect_identical(backend_available_layers(be), "data")
  expect_identical(backend_default_layer(be),    "data")
  expect_identical(backend_genes(be),            c("A", "B"))
  expect_identical(backend_n_genes(be),          2L)
  expect_true(backend_has_gene(be, "A"))
  expect_false(backend_has_gene(be, "Z"))
  expect_identical(backend_get_gene(be, "A"), c(0, 1, 2))
  expect_null(backend_get_gene(be, "Z"))
})

# ---- Layered shape ------------------------------------------------------

test_that("layered shape is preserved and layer-aware accessors dispatch", {
  expr <- list(
    data   = list(A = c(0.1, 0.2), B = c(1.0, 2.0)),
    counts = list(A = c(1L,   2L), B = c(10L,  20L))
  )
  be <- expression_backend_inmemory(expr)
  expect_setequal(backend_available_layers(be), c("data", "counts"))
  expect_identical(backend_default_layer(be), "data")
  expect_identical(backend_n_cells(be), 2L)

  expect_identical(backend_get_gene(be, "A"),                 c(0.1, 0.2))
  expect_identical(backend_get_gene(be, "A", layer = "data"), c(0.1, 0.2))
  expect_identical(backend_get_gene(be, "A", layer = "counts"), c(1, 2))
  expect_identical(backend_n_genes(be, layer = "counts"),     2L)
})

test_that("default_layer falls back to the first layer when the requested name is absent", {
  expr <- list(counts = list(A = c(1, 2)))
  be   <- expression_backend_inmemory(expr, default_layer = "data")
  # "data" requested but absent; ctor silently falls back to the first
  # available layer so callers without a layer arg still get something.
  expect_identical(backend_default_layer(be), "counts")
  expect_identical(backend_get_gene(be, "A"), c(1, 2))
})

# ---- Validation ---------------------------------------------------------

test_that("constructor rejects vectors of inconsistent length", {
  expect_error(
    expression_backend_inmemory(list(A = c(0, 1, 2), B = c(0, 1))),
    "different from n_cells"
  )
})

test_that("constructor rejects duplicate gene names within a layer", {
  bad <- list(A = c(0, 1), A = c(2, 3))
  expect_error(expression_backend_inmemory(bad), "duplicate gene")
})

test_that("constructor rejects unnamed gene vectors", {
  expect_error(
    expression_backend_inmemory(list(data = list(c(0, 1, 2)))),
    "named gene vectors"
  )
})

test_that("constructor rejects unnamed layers", {
  expect_error(
    expression_backend_inmemory(list(list(A = c(0, 1)))),
    "named layers"
  )
})

test_that("layered shape with non-list / non-numeric elements is rejected", {
  expect_error(
    expression_backend_inmemory(list(data = "not a gene list")),
    "named lists of numeric vectors"
  )
})

test_that("constructor accepts an explicit empty backend", {
  be <- expression_backend_inmemory(list(), n_cells = 100)
  expect_identical(backend_n_cells(be), 100L)
  expect_identical(backend_n_genes(be), 0L)
  expect_identical(backend_genes(be),   character())
  expect_false(backend_has_gene(be, "anything"))
})

# ---- Layer resolution ---------------------------------------------------

test_that("requesting an unknown layer errors at the accessor", {
  be <- expression_backend_inmemory(list(A = c(0, 1)))
  expect_error(backend_get_gene(be, "A", layer = "nope"), "not available")
  expect_error(backend_genes(be, layer = "nope"),         "not available")
})

# ---- has_gene edge cases ------------------------------------------------

test_that("has_gene gracefully handles NA / NULL / empty input", {
  be <- expression_backend_inmemory(list(A = c(0, 1)))
  expect_false(backend_has_gene(be, NULL))
  expect_false(backend_has_gene(be, ""))
  expect_false(backend_has_gene(be, NA_character_))
})

# ---- as_expression_backend ---------------------------------------------

test_that("as_expression_backend is identity on a backend", {
  be <- expression_backend_inmemory(list(A = c(0, 1)))
  expect_identical(as_expression_backend(be), be)
})

test_that("as_expression_backend wraps NULL / empty / flat legacy shapes", {
  e1 <- as_expression_backend(NULL)
  e2 <- as_expression_backend(list())
  e3 <- as_expression_backend(list(A = c(0, 1, 2)))
  expect_s3_class(e1, "expression_backend")
  expect_s3_class(e2, "expression_backend")
  expect_s3_class(e3, "expression_backend")
  expect_identical(backend_n_cells(e1), 0L)
  expect_identical(backend_n_cells(e2), 0L)
  expect_identical(backend_n_cells(e3), 3L)
  expect_identical(backend_get_gene(e3, "A"), c(0, 1, 2))
})

test_that("as_expression_backend rejects unhandled shapes", {
  expect_error(as_expression_backend("not a list"), "Cannot coerce")
  expect_error(as_expression_backend(list(A = "x")), "Cannot coerce")
})

# ---- Print / format don't crash -----------------------------------------

test_that("print and format produce single-line summaries", {
  be  <- expression_backend_inmemory(list(A = c(0, 1, 2)))
  fmt <- format(be)
  expect_match(fmt, "3 cells")
  expect_match(fmt, "1 gene")
  expect_match(fmt, "1 layer")
  expect_output(print(be), "expression_backend")
})

# ---- mock_dataset integration ------------------------------------------

test_that("mock_dataset now stores an expression_backend, not a bare list", {
  ds <- mock_dataset(n_cells = 30)
  expect_s3_class(ds$expression, "expression_backend")
  expect_identical(backend_n_cells(ds$expression), 30L)
  # Mock dataset now exposes both layers: 'data' (log-normalised, default)
  # for cell-level analyses, and 'counts' (raw counts) for pseudobulk DE.
  expect_setequal(backend_available_layers(ds$expression), c("data", "counts"))
  expect_identical(backend_default_layer(ds$expression), "data")
  expect_setequal(backend_genes(ds$expression),
                  c("CD3D", "MS4A1", "LST1", "EPCAM", "COL1A1", "NKG7"))
  expect_setequal(backend_genes(ds$expression, layer = "counts"),
                  c("CD3D", "MS4A1", "LST1", "EPCAM", "COL1A1", "NKG7"))
})

test_that("mock_dataset counts layer carries plausible integer-style counts", {
  ds <- mock_dataset(n_cells = 200, seed = 21)
  cd3d_counts <- get_gene_expression(ds, "CD3D", layer = "counts")
  cd3d_data   <- get_gene_expression(ds, "CD3D", layer = "data")
  expect_identical(length(cd3d_counts), 200L)
  # Counts should be non-negative and integer-valued (Poisson draws).
  expect_true(all(cd3d_counts >= 0))
  expect_true(all(cd3d_counts == round(cd3d_counts)))
  # The two layers are independently generated -- they must not be
  # identical (would mean the counts layer was accidentally aliased).
  expect_false(isTRUE(all.equal(cd3d_counts, cd3d_data)))
})

test_that("get_gene_expression dispatches through the backend for mock_dataset", {
  ds <- mock_dataset(n_cells = 40)
  v  <- get_gene_expression(ds, "CD3D")
  expect_type(v, "double")
  expect_identical(length(v), 40L)
  expect_null(get_gene_expression(ds, "DOESNOTEXIST"))
})

# ---- Sparse backend (Matrix-backed) -------------------------------------

test_that("expression_backend_sparse exposes the same generic surface as in-memory", {
  skip_if_not_installed("Matrix")
  m <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 1, 2),
    j = c(1, 1, 2, 3, 3),
    x = c(5, 2, 7, 1, 9),
    dims = c(3, 3),
    dimnames = list(c("A", "B", "C"), c("c1", "c2", "c3"))
  )
  be <- expression_backend_sparse(m)
  expect_s3_class(be, "expression_backend")
  expect_s3_class(be, "expression_backend_sparse")
  expect_identical(backend_n_cells(be), 3L)
  expect_identical(backend_available_layers(be), "data")
  expect_identical(backend_default_layer(be), "data")
  expect_identical(backend_genes(be), c("A", "B", "C"))
  expect_identical(backend_n_genes(be), 3L)
  expect_true(backend_has_gene(be, "A"))
  expect_false(backend_has_gene(be, "Z"))
  # Row A is (1,1)=5 and (1,3)=1 -> c(5, 0, 1)
  # Row B is (2,1)=2 and (2,3)=9 -> c(2, 0, 9)
  # Row C is (3,2)=7              -> c(0, 7, 0)
  expect_identical(backend_get_gene(be, "A"), c(5, 0, 1))
  expect_identical(backend_get_gene(be, "B"), c(2, 0, 9))
  expect_identical(backend_get_gene(be, "C"), c(0, 7, 0))
  expect_null(backend_get_gene(be, "Z"))
})

test_that("expression_backend_sparse supports multiple layers (e.g. data + counts)", {
  skip_if_not_installed("Matrix")
  data_m <- Matrix::sparseMatrix(
    i = 1:2, j = 1:2, x = c(0.5, 1.5),
    dims = c(2, 2),
    dimnames = list(c("A", "B"), c("c1", "c2")))
  counts_m <- Matrix::sparseMatrix(
    i = 1:2, j = 1:2, x = c(1L, 3L),
    dims = c(2, 2),
    dimnames = list(c("A", "B"), c("c1", "c2")))
  be <- expression_backend_sparse(list(data = data_m, counts = counts_m))

  expect_setequal(backend_available_layers(be), c("data", "counts"))
  expect_identical(backend_default_layer(be), "data")
  expect_identical(backend_get_gene(be, "A"),                  c(0.5, 0))
  expect_identical(backend_get_gene(be, "A", layer = "counts"), c(1, 0))
})

test_that("expression_backend_sparse rejects rownames-less matrices and ncol mismatch", {
  skip_if_not_installed("Matrix")
  m <- Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(2, 2))  # no dimnames
  expect_error(expression_backend_sparse(m), "rownames")

  m1 <- Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(2, 2),
                             dimnames = list(c("A", "B"), c("c1", "c2")))
  m2 <- Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(2, 3),
                             dimnames = list(c("A", "B"), c("c1", "c2", "c3")))
  expect_error(expression_backend_sparse(list(data = m1, counts = m2)),
               "ncol != n_cells")
})

test_that("backend_as_matrix returns a (genes x cells) view for both backends", {
  # In-memory backend: reconstructs a dense matrix from the named-list layer
  be_mem <- expression_backend_inmemory(list(
    A = c(0, 1, 2, 3),
    B = c(10, 20, 30, 40)))
  m_mem <- backend_as_matrix(be_mem)
  expect_true(is.matrix(m_mem) || inherits(m_mem, "Matrix"))
  expect_identical(dim(m_mem), c(2L, 4L))
  expect_identical(rownames(m_mem), c("A", "B"))
  expect_identical(unname(m_mem["A", ]), c(0, 1, 2, 3))
  expect_identical(unname(m_mem["B", ]), c(10, 20, 30, 40))

  # Sparse backend: returns the underlying matrix unchanged (no copy).
  skip_if_not_installed("Matrix")
  sp <- Matrix::sparseMatrix(
    i = c(1, 2), j = c(1, 3), x = c(5, 7),
    dims = c(2, 3),
    dimnames = list(c("A", "B"), c("c1", "c2", "c3")))
  be_sp <- expression_backend_sparse(sp)
  m_sp  <- backend_as_matrix(be_sp)
  expect_identical(dim(m_sp), c(2L, 3L))
  expect_identical(rownames(m_sp), c("A", "B"))
  expect_identical(as.numeric(m_sp["A", ]), c(5, 0, 0))
})

test_that("backend_as_matrix handles an empty backend", {
  be <- expression_backend_inmemory(list(), n_cells = 5)
  m  <- backend_as_matrix(be)
  expect_identical(dim(m), c(0L, 5L))
})

test_that("backend_as_matrix respects the requested layer", {
  be <- expression_backend_inmemory(list(
    data   = list(A = c(0, 1), B = c(2, 3)),
    counts = list(A = c(10, 20), B = c(30, 40))))
  expect_identical(unname(backend_as_matrix(be, layer = "data")["A", ]),   c(0, 1))
  expect_identical(unname(backend_as_matrix(be, layer = "counts")["A", ]), c(10, 20))
})

test_that("get_gene_expression routes through the sparse backend transparently", {
  skip_if_not_installed("Matrix")
  m <- Matrix::sparseMatrix(
    i = c(1, 1, 2),
    j = c(1, 2, 3),
    x = c(0.7, 0.3, 5.0),
    dims = c(2, 3),
    dimnames = list(c("CD3D", "MS4A1"), c("c1", "c2", "c3")))
  ds <- list(name = "sparse_demo", n_cells = 3L, cells = colnames(m),
             genes = rownames(m),
             expression = expression_backend_sparse(m))
  expect_identical(get_gene_expression(ds, "CD3D"),  c(0.7, 0.3, 0))
  expect_identical(get_gene_expression(ds, "MS4A1"), c(0,   0,   5))
  expect_null(get_gene_expression(ds, "MISSING"))
  expect_identical(available_genes(ds), c("CD3D", "MS4A1"))
})

test_that("legacy bare-list `dataset$expression` keeps working through helpers", {
  # Simulates an older dataset object built before #2 landed: the
  # helpers must transparently coerce it to a backend on read.
  legacy <- list(
    name       = "legacy",
    cells      = paste0("c", 1:5),
    n_cells    = 5L,
    genes      = c("A", "B"),
    expression = list(A = c(0, 1, 2, 3, 4), B = c(5, 4, 3, 2, 1))
  )
  expect_identical(available_genes(legacy), c("A", "B"))
  expect_identical(get_gene_expression(legacy, "A"), c(0, 1, 2, 3, 4))
  expect_identical(get_gene_expression(legacy, "B"), c(5, 4, 3, 2, 1))
  expect_null(get_gene_expression(legacy, "Z"))
  expect_true(validate_gene(legacy, "A"))
  expect_false(validate_gene(legacy, "Z"))
})
