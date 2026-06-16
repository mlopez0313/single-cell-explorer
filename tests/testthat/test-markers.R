# `compute_markers()` emits per-group heartbeat messages by default so
# the Shiny module's log file shows progress on long runs. Tests don't
# care about them and would otherwise pollute testthat output -- gate
# them off for the whole file.
local({
  old <- options(sce.marker_progress = FALSE)
  withr::defer(options(old), teardown_env())
})

test_that("compute_markers defaults to wilcox", {
  expect_identical(eval(formals(compute_markers)$test)[1], "wilcox")
  ds <- mock_dataset(n_cells = 120, seed = 2)
  default_out <- compute_markers(ds, grouping_field = "cluster", top_n = 3)
  wilcox_out  <- compute_markers(ds, grouping_field = "cluster", top_n = 3,
                                 test = "wilcox")
  expect_identical(default_out, wilcox_out)
})

test_that("compute_markers (t-test) also runs and returns the same column schema", {
  ds <- mock_dataset(n_cells = 100, seed = 5)
  out <- compute_markers(ds, grouping_field = "cluster", top_n = 5,
                         test = "t")
  expect_s3_class(out, "data.frame")
  expect_true(all(c("group", "gene", "avg_log2FC", "pct_in", "pct_out", "p_value")
                  %in% names(out)))
})

test_that("compute_markers returns NULL for missing fields / empty filter", {
  ds <- mock_dataset(n_cells = 60)
  expect_null(compute_markers(ds, grouping_field = "doesnotexist"))
  expect_null(compute_markers(ds, grouping_field = "cluster",
                              group_filter = "doesnotexist"))
})

test_that("compute_markers respects group_filter and top_n", {
  ds <- mock_dataset(n_cells = 120)
  out <- compute_markers(ds, grouping_field = "cluster",
                         group_filter = "0", top_n = 2)
  expect_true(all(out$group == "0"))
  expect_true(nrow(out) <= 2L)
})
