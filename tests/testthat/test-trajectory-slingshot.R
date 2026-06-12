# Tests for R/trajectory_slingshot.R
# ----------------------------------------------------------------------------
# The pure converter `.slingshot_to_pseudotime()` accepts either a real
# SlingshotDataSet or a list stand-in shaped like one, so we can
# regression-test schema mapping without the Bioconductor package.

# ---- Registry + missing-dep gating --------------------------------------

test_that("slingshot is registered with `requires = 'slingshot'`", {
  spec <- get_trajectory_method("slingshot")
  expect_identical(spec$id,       "slingshot")
  expect_identical(spec$requires, "slingshot")
  expect_true(isTRUE(spec$requires_root))
})

test_that("run_trajectory('slingshot') errors cleanly without the package", {
  skip_if(has_optional("slingshot"),
          "slingshot installed; install-error path skipped")
  ds <- mock_dataset(n_cells = 80, seed = 5)
  expect_error(
    run_trajectory(ds, "slingshot", reduction = "UMAP",
                   root_field = "cluster", root_group = "0",
                   cluster_field = "cluster"),
    regexp = "slingshot", fixed = FALSE
  )
})

# ---- .slingshot_to_pseudotime: pure converter --------------------------

test_that(".slingshot_to_pseudotime extracts per-cell pst + weighted aggregation", {
  cells <- sprintf("c_%02d", 1:6)
  # Two lineages: cells 1-3 belong to lineage 1, cells 4-6 to lineage 2.
  pst <- matrix(c(
    0,   NA,
    0.5, NA,
    1.0, NA,
    NA,  0,
    NA,  0.5,
    NA,  1.0
  ), nrow = 6, byrow = TRUE, dimnames = list(cells, c("L1", "L2")))
  # cells 1-3 belong to lineage 1 (weight 1 on L1, 0 on L2);
  # cells 4-6 belong to lineage 2 (weight 0 on L1, 1 on L2).
  wts <- matrix(c(
    1, 0,   # cell 1
    1, 0,   # cell 2
    1, 0,   # cell 3
    0, 1,   # cell 4
    0, 1,   # cell 5
    0, 1    # cell 6
  ), nrow = 6, byrow = TRUE,
     dimnames = list(cells, c("L1", "L2")))
  sds <- list(pseudotime = pst, curveWeights = wts)

  out <- .slingshot_to_pseudotime(sds, cells = cells,
                                  cluster_field = "cluster",
                                  start_clus     = "0",
                                  reduction_used = "UMAP")

  expect_identical(out$source,         "slingshot")
  expect_identical(out$reduction_used, "UMAP")
  expect_identical(out$root_field,     "cluster")
  expect_identical(out$root_group,     "0")
  expect_identical(out$n_lineages,     2L)
  expect_identical(length(out$pseudotime), 6L)
  expect_true(all(out$pseudotime >= 0 - 1e-9 & out$pseudotime <= 1 + 1e-9))
  # method_details exposes the raw matrices
  expect_identical(out$method_details$lineage_psts,  pst)
  expect_identical(out$method_details$curve_weights, wts)
})

test_that(".slingshot_to_pseudotime averages across lineages when weights are absent", {
  cells <- c("a", "b", "c")
  pst <- matrix(c(0,   0.0,
                  0.5, 0.5,
                  1.0, 1.0), nrow = 3, byrow = TRUE)
  sds <- list(pseudotime = pst) # no curveWeights
  out <- .slingshot_to_pseudotime(sds, cells = cells,
                                  cluster_field = "cluster",
                                  start_clus = NULL,
                                  reduction_used = "PCA")
  # Equal-valued lineages -> rescaled to [0, 1]
  expect_equal(out$pseudotime[1], 0, tolerance = 1e-9)
  expect_equal(out$pseudotime[3], 1, tolerance = 1e-9)
})

test_that(".slingshot_to_pseudotime errors when pseudotime is missing", {
  sds <- list() # no pseudotime
  expect_error(
    .slingshot_to_pseudotime(sds, cells = c("a", "b"),
                             cluster_field = "cluster"),
    regexp = "no `pseudotime` matrix", fixed = TRUE
  )
})

test_that(".slingshot_to_pseudotime errors on row-count mismatch", {
  pst <- matrix(c(0.1, 0.2), nrow = 2, ncol = 1)
  sds <- list(pseudotime = pst)
  expect_error(
    .slingshot_to_pseudotime(sds, cells = c("a", "b", "c"),
                             cluster_field = "cluster"),
    regexp = "length\\(cells\\)", fixed = FALSE
  )
})

test_that(".slingshot_to_pseudotime: NA-only cells become NA after rescale", {
  cells <- c("a", "b", "c")
  pst <- matrix(c(NA, NA, NA,
                  0,  0,  1), nrow = 3, byrow = TRUE)
  wts <- matrix(c(0, 0, 1, 1, 1, 1), nrow = 3, byrow = FALSE)
  sds <- list(pseudotime = pst, curveWeights = wts)
  out <- .slingshot_to_pseudotime(sds, cells = cells, cluster_field = "cluster")
  expect_true(is.na(out$pseudotime[1]))
})

# ---- End-to-end (skipped without slingshot) -----------------------------

test_that("run_trajectory('slingshot') round-trips on a real Slingshot run", {
  skip_if_not_installed("slingshot")
  skip_if_not_installed("SingleCellExperiment")
  ds <- mock_dataset(n_cells = 300, seed = 9)
  out <- run_trajectory(ds, "slingshot", reduction = "UMAP",
                       root_field = "cluster", root_group = "0",
                       cluster_field = "cluster")
  expect_identical(out$source, "slingshot")
  expect_identical(length(out$pseudotime), 300L)
  expect_true(all(is.finite(out$pseudotime) | is.na(out$pseudotime)))
  expect_true(out$n_lineages >= 1L)
})
