test_that("mock_dataset adheres to the documented schema", {
  ds <- mock_dataset(n_cells = 120, n_genes = 500, name = "ds_schema")

  schema_keys <- c(
    "name", "source", "n_cells", "n_genes", "assays", "default_assay",
    "reductions", "default_reduction", "metadata_fields", "cells",
    "cell_data", "genes", "expression"
  )
  expect_true(all(schema_keys %in% names(ds)))
  expect_identical(ds$source, "mock")
  expect_identical(ds$n_cells, 120L)
  expect_identical(ds$n_genes, 500L)
  expect_true(ds$default_assay %in% ds$assays)
  expect_true(ds$default_reduction %in% ds$reductions)
  expect_identical(length(ds$cells), 120L)
})

test_that("mock_dataset cell_data carries embeddings + metadata + pseudotime_demo", {
  ds <- mock_dataset(n_cells = 150)
  cd <- ds$cell_data
  expect_true(all(c("cell", "sample", "cluster", "condition", "cell_type",
                    "pseudotime_demo",
                    "UMAP_1", "UMAP_2", "PCA_1", "PCA_2", "tSNE_1", "tSNE_2")
                  %in% names(cd)))
  expect_identical(nrow(cd), 150L)
  expect_type(cd$pseudotime_demo, "double")
  expect_true(all(cd$pseudotime_demo >= 0 - 1e-9))
  expect_true(all(cd$pseudotime_demo <= 1 + 1e-9))
  expect_true("pseudotime_demo" %in% ds$metadata_fields)
})

test_that("mock_dataset gene expression vectors align with cells (via helper)", {
  ds <- mock_dataset(n_cells = 200)
  expect_true(all(c("CD3D", "MS4A1", "LST1", "EPCAM", "COL1A1", "NKG7") %in% ds$genes))
  # `dataset$expression` is now an expression_backend object -- modules
  # must read it via get_gene_expression() rather than `$expression[[g]]`.
  expect_s3_class(ds$expression, "expression_backend")
  for (g in c("CD3D", "MS4A1", "LST1")) {
    e <- get_gene_expression(ds, g)
    expect_type(e, "double")
    expect_identical(length(e), 200L)
  }
})

test_that("mock_dataset is reproducible with seed", {
  a <- mock_dataset(n_cells = 80, seed = 7)
  b <- mock_dataset(n_cells = 80, seed = 7)
  expect_identical(a$cell_data$cluster, b$cell_data$cluster)
  expect_identical(get_gene_expression(a, "CD3D"),
                   get_gene_expression(b, "CD3D"))
})

# ---- Dataset helpers ------------------------------------------------------

test_that("available_* helpers return what mock_dataset advertises", {
  ds <- mock_dataset(n_cells = 60)
  expect_identical(available_assays(ds), ds$assays)
  expect_identical(available_reductions(ds), ds$reductions)
  expect_identical(available_metadata_fields(ds), ds$metadata_fields)
  expect_identical(available_genes(ds), ds$genes)
})

test_that("get_embedding returns a data.frame(cell, x, y) for known reductions", {
  ds <- mock_dataset(n_cells = 75)
  emb <- get_embedding(ds, "UMAP")
  expect_s3_class(emb, "data.frame")
  expect_identical(names(emb), c("cell", "x", "y"))
  expect_identical(nrow(emb), 75L)

  # Unknown reduction -> NULL
  expect_null(get_embedding(ds, "DOESNOTEXIST"))
  expect_null(get_embedding(ds, NULL))
})

test_that("get_metadata returns aligned vectors and NULL for missing fields", {
  ds <- mock_dataset(n_cells = 90)
  cl <- get_metadata(ds, "cluster")
  expect_identical(length(cl), 90L)
  expect_null(get_metadata(ds, "not_a_field"))
})

test_that("get_gene_expression returns numeric or NULL for unknown genes", {
  ds <- mock_dataset(n_cells = 50)
  expect_type(get_gene_expression(ds, "CD3D"), "double")
  expect_identical(length(get_gene_expression(ds, "CD3D")), 50L)
  expect_null(get_gene_expression(ds, "NOTAGENE"))
})

test_that("validate_gene / validate_metadata flag the right things", {
  ds <- mock_dataset(n_cells = 50)
  expect_true(validate_gene(ds, "CD3D"))
  expect_false(validate_gene(ds, "NOTAGENE"))
  expect_true(validate_metadata(ds, "cluster"))
  expect_false(validate_metadata(ds, "not_a_field"))
})
