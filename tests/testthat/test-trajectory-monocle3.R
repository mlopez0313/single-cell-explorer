# Tests for R/trajectory_monocle3.R
# ----------------------------------------------------------------------------
# Monocle3 is heavy (S4-class CDS, dependencies on SingleCellExperiment).
# We test:
#   * registry entry shape
#   * missing-dep error path
#   * the pure schema converter .monocle3_to_pseudotime() against a
#     hand-built named numeric vector (no monocle3 install required)

# ---- Registry + missing-dep gating --------------------------------------

test_that("monocle3 is registered with `requires = 'monocle3'`", {
  spec <- get_trajectory_method("monocle3")
  expect_identical(spec$id,       "monocle3")
  expect_identical(spec$requires, "monocle3")
  expect_true(isTRUE(spec$requires_root))
})

test_that("run_trajectory('monocle3') errors cleanly without the package", {
  skip_if(has_optional("monocle3"),
          "monocle3 installed; install-error path skipped")
  ds <- mock_dataset(n_cells = 80, seed = 13)
  expect_error(
    run_trajectory(ds, "monocle3", reduction = "UMAP",
                   root_field = "cluster", root_group = "0"),
    regexp = "monocle3", fixed = FALSE
  )
})

# ---- .monocle3_to_pseudotime: pure converter ----------------------------

test_that(".monocle3_to_pseudotime lifts a per-cell pst vector into the schema", {
  cells <- sprintf("c_%02d", 1:5)
  pt <- stats::setNames(c(0, 0.25, 0.5, 0.75, 1.0), cells)
  out <- .monocle3_to_pseudotime(
    pt, cells = cells, cluster_field = "cluster",
    root_group = "0", reduction_used = "UMAP",
    method_details = list(n_partitions = 1L, layer_used = "counts"))
  expect_identical(out$source,         "monocle3")
  expect_identical(out$reduction_used, "UMAP")
  expect_identical(out$root_field,     "cluster")
  expect_identical(out$root_group,     "0")
  expect_identical(out$n_lineages,     1L)
  expect_identical(out$method_details$layer_used, "counts")
  expect_identical(length(out$pseudotime), 5L)
  expect_true(all(out$pseudotime >= 0 - 1e-9 & out$pseudotime <= 1 + 1e-9))
})

test_that(".monocle3_to_pseudotime errors on length mismatch", {
  expect_error(
    .monocle3_to_pseudotime(c(0, 0.5), cells = c("a", "b", "c"),
                            cluster_field = "cluster"),
    regexp = "length\\(pt\\)", fixed = FALSE
  )
})

test_that(".monocle3_to_pseudotime errors on empty input", {
  expect_error(
    .monocle3_to_pseudotime(numeric(), cells = character(),
                            cluster_field = "cluster"),
    regexp = "empty pseudotime", fixed = TRUE
  )
})

# ---- End-to-end (skipped without monocle3) ------------------------------

test_that("run_trajectory('monocle3') round-trips on a real Monocle3 run", {
  skip_if_not_installed("monocle3")
  skip_if_not_installed("SingleCellExperiment")
  ds <- mock_dataset(n_cells = 300, seed = 17)
  out <- run_trajectory(ds, "monocle3", reduction = "UMAP",
                       root_field = "cluster", root_group = "0")
  expect_identical(out$source, "monocle3")
  expect_identical(length(out$pseudotime), 300L)
  expect_true(out$n_lineages >= 1L)
})
