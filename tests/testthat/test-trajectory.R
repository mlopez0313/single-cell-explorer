test_that("available_pseudotime_sources reflects the trajectory method registry", {
  src <- available_pseudotime_sources()
  # The registry advertises every method, including optional ones.
  expect_setequal(unname(src),
                  c("mock", "metadata", "slingshot", "monocle3"))
  # Built-in methods (no optional deps) are always present.
  expect_true(all(c("mock", "metadata") %in% unname(src)))
})

test_that("available_numeric_metadata_fields finds the pre-shipped pseudotime_demo", {
  ds <- mock_dataset(n_cells = 80)
  nums <- available_numeric_metadata_fields(ds)
  expect_true("pseudotime_demo" %in% nums)
})

test_that("compute_pseudotime (mock) is deterministic and root-centred", {
  ds <- mock_dataset(n_cells = 120, seed = 11)
  pt <- compute_pseudotime(ds, source = "mock", reduction = "UMAP",
                           root_field = "cluster", root_group = "0")
  expect_identical(length(pt), 120L)
  expect_true(all(pt >= 0 - 1e-9 & pt <= 1 + 1e-9))
  # Root cells should sit lower on average than non-root
  in_root <- ds$cell_data$cluster == "0"
  expect_lt(mean(pt[in_root]), mean(pt[!in_root]))
})

test_that("compute_pseudotime (metadata) rescales a numeric column to [0, 1]", {
  ds <- mock_dataset(n_cells = 100)
  pt <- compute_pseudotime(ds, source = "metadata",
                           metadata_field = "pseudotime_demo")
  expect_identical(length(pt), 100L)
  expect_true(all(abs(pt - ds$cell_data$pseudotime_demo) < 1e-9))
})

test_that("compute_pseudotime rejects bad inputs with clear messages", {
  ds <- mock_dataset(n_cells = 60)
  expect_error(compute_pseudotime(ds, source = "metadata",
                                  metadata_field = "cluster"),
               "not numeric")
  expect_error(compute_pseudotime(ds, source = "mock", reduction = "UMAP",
                                  root_field = "cluster", root_group = "999"),
               "no cells")
  expect_error(compute_pseudotime(ds, source = "mock",
                                  reduction = "DOESNOTEXIST",
                                  root_field = "cluster", root_group = "0"))
})

test_that("bin_gene_by_pseudotime preserves total cells across bins", {
  ds <- mock_dataset(n_cells = 100)
  pt <- compute_pseudotime(ds, source = "mock", reduction = "UMAP",
                           root_field = "cluster", root_group = "0")
  expr <- get_gene_expression(ds, "CD3D")
  binned <- bin_gene_by_pseudotime(pt, expr, n_bins = 10)
  expect_s3_class(binned, "data.frame")
  expect_identical(nrow(binned), 10L)
  expect_identical(sum(binned$n), 100L)
  expect_true(all(c("bin", "pt_mid", "expr_mean", "n") %in% names(binned)))
})

test_that("has_trajectory_results checks status + non-empty pseudotime", {
  state <- new_app_state()
  ds <- mock_dataset(n_cells = 50)
  with_state(state, set_active_dataset(state, ds))
  expect_false(with_state(state, has_trajectory_results(state)))

  pt <- compute_pseudotime(ds, source = "metadata",
                           metadata_field = "pseudotime_demo")
  with_state(state, {
    state$analysis_results <- list(trajectory = list(
      status = "completed", results = list(pseudotime = pt)))
  })
  expect_true(with_state(state, has_trajectory_results(state)))
})

# ---- is_trajectory_result --------------------------------------------------

test_that("is_trajectory_result accepts a canonical payload and rejects junk", {
  ds <- mock_dataset(n_cells = 30)
  out <- run_trajectory(ds, source = "metadata",
                        metadata_field = "pseudotime_demo")
  expect_true(is_trajectory_result(out))

  expect_false(is_trajectory_result(NULL))
  expect_false(is_trajectory_result(list()))
  expect_false(is_trajectory_result(list(pseudotime = "not numeric",
                                         cell = "c", source = "x")))
  expect_false(is_trajectory_result(list(pseudotime = c(0.1, 0.2),
                                         cell = "c_01",  # wrong length
                                         source = "metadata")))
  expect_false(is_trajectory_result(list(pseudotime = c(0.1, 0.2),
                                         cell = c("c_01", "c_02"),
                                         source = c("metadata", "metadata"))))
})

# ---- apply_pseudotime_to_dataset ------------------------------------------

test_that("apply_pseudotime_to_dataset writes a dated numeric column with provenance", {
  ds  <- mock_dataset(n_cells = 60, seed = 7)
  out <- run_trajectory(ds, source = "metadata",
                        metadata_field = "pseudotime_demo")
  ts  <- as.POSIXct("2025-06-15 12:00:00", tz = "UTC")
  ds2 <- apply_pseudotime_to_dataset(ds, out, applied_at = ts)
  col <- "pseudotime__metadata__2025_06_15"
  expect_true(col %in% names(ds2$cell_data))
  expect_true(col %in% ds2$metadata_fields)
  expect_identical(length(ds2$cell_data[[col]]), 60L)
  expect_true(is.numeric(ds2$cell_data[[col]]))
  # Values must match the canonical pseudotime, aligned by cell.
  pos <- match(ds$cell_data$cell, out$cell)
  expect_equal(as.numeric(ds2$cell_data[[col]]),
               out$pseudotime[pos], tolerance = 1e-12)
  # Provenance attrs.
  expect_identical(attr(ds2$cell_data[[col]], "pseudotime_source"), "metadata")
  expect_identical(attr(ds2$cell_data[[col]], "metadata_field"),    "pseudotime_demo")
  expect_identical(attr(ds2$cell_data[[col]], "kind"),              "numeric")
  expect_identical(attr(ds2$cell_data[[col]], "applied_at"),        ts)
  # Idempotent on rerun -> refuses to overwrite.
  expect_error(apply_pseudotime_to_dataset(ds2, out, applied_at = ts),
               regexp = "already exists", fixed = TRUE)
})

test_that("apply_pseudotime_to_dataset adds a bin factor column when bins > 0", {
  ds  <- mock_dataset(n_cells = 100, seed = 8)
  out <- run_trajectory(ds, source = "mock", reduction = "UMAP",
                        root_field = "cluster", root_group = "0")
  ts  <- as.POSIXct("2025-06-15 12:00:00", tz = "UTC")
  ds2 <- apply_pseudotime_to_dataset(ds, out, bins = 5L, applied_at = ts)
  num_col <- "pseudotime__mock__2025_06_15"
  bin_col <- "pseudotime_bin__mock__2025_06_15"
  expect_true(num_col %in% names(ds2$cell_data))
  expect_true(bin_col %in% names(ds2$cell_data))
  expect_true(num_col %in% ds2$metadata_fields)
  expect_true(bin_col %in% ds2$metadata_fields)
  expect_s3_class(ds2$cell_data[[bin_col]], "factor")
  expect_identical(levels(ds2$cell_data[[bin_col]]),
                   sprintf("bin_%02d", 1:5))
  # Every non-NA bin label must align with the numeric column.
  pt_col  <- ds2$cell_data[[num_col]]
  bin_lvl <- as.integer(ds2$cell_data[[bin_col]])
  ok <- !is.na(pt_col)
  brk <- seq(min(pt_col[ok]), max(pt_col[ok]), length.out = 6L)
  exp_lvl <- as.integer(cut(pt_col[ok], breaks = brk,
                            include.lowest = TRUE, labels = FALSE))
  expect_identical(bin_lvl[ok], exp_lvl)
  expect_identical(attr(ds2$cell_data[[bin_col]], "bins"), 5L)
  expect_identical(attr(ds2$cell_data[[bin_col]], "kind"), "bin")
})

test_that("apply_pseudotime_to_dataset accepts a wrapped state slot", {
  ds  <- mock_dataset(n_cells = 40, seed = 9)
  out <- run_trajectory(ds, source = "metadata",
                        metadata_field = "pseudotime_demo")
  slot <- list(status = "completed", results = out,
               params = list(), error_message = NULL,
               timestamp = Sys.time(), duration_ms = 1L)
  ds2 <- apply_pseudotime_to_dataset(ds, slot,
                                     applied_at = as.POSIXct("2025-06-15",
                                                             tz = "UTC"))
  expect_true("pseudotime__metadata__2025_06_15" %in% names(ds2$cell_data))
})

test_that("apply_pseudotime_to_dataset rejects a non-completed slot", {
  ds  <- mock_dataset(n_cells = 20)
  out <- run_trajectory(ds, source = "metadata",
                        metadata_field = "pseudotime_demo")
  failed_slot <- list(status = "failed", results = out,
                      error_message = "boom")
  expect_error(apply_pseudotime_to_dataset(ds, failed_slot),
               regexp = "not a completed", fixed = TRUE)
  running_slot <- list(status = "running", results = NULL)
  expect_error(apply_pseudotime_to_dataset(ds, running_slot),
               regexp = "not a completed", fixed = TRUE)
})

test_that("apply_pseudotime_to_dataset errors on missing cells", {
  ds  <- mock_dataset(n_cells = 30)
  out <- run_trajectory(ds, source = "metadata",
                        metadata_field = "pseudotime_demo")
  # Drop the last cell from the trajectory result -> dataset has 1 unmatched.
  out_bad <- out
  out_bad$cell       <- out$cell[-length(out$cell)]
  out_bad$pseudotime <- out$pseudotime[-length(out$pseudotime)]
  expect_error(apply_pseudotime_to_dataset(ds, out_bad),
               regexp = "dataset cells have no value", fixed = TRUE)
})

test_that("apply_pseudotime_to_dataset rejects schema-invalid input", {
  ds <- mock_dataset(n_cells = 20)
  expect_error(apply_pseudotime_to_dataset(ds, list(foo = 1)),
               regexp = "not a completed", fixed = TRUE)
  expect_error(apply_pseudotime_to_dataset(ds, NULL),
               regexp = "not a completed", fixed = TRUE)
  expect_error(apply_pseudotime_to_dataset(NULL, list(pseudotime = 1,
                                                      cell = "a",
                                                      source = "metadata")),
               regexp = "No dataset", fixed = TRUE)
})

test_that("apply_pseudotime_to_dataset validates bins arg", {
  ds  <- mock_dataset(n_cells = 30)
  out <- run_trajectory(ds, source = "metadata",
                        metadata_field = "pseudotime_demo")
  expect_error(apply_pseudotime_to_dataset(ds, out, bins = -1L),
               regexp = "non-negative", fixed = TRUE)
  # bins = 0 (default) just emits the numeric column.
  ds2 <- apply_pseudotime_to_dataset(ds, out, bins = 0L,
                                     applied_at = as.POSIXct("2025-06-15",
                                                             tz = "UTC"))
  expect_false(any(grepl("^pseudotime_bin__", names(ds2$cell_data))))
})
