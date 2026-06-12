test_that("IMPUTATION_METHODS exposes the registered mock methods", {
  expect_true(all(c("none", "neighbor", "alra_mock", "magic_mock")
                  %in% unname(IMPUTATION_METHODS)))
})

test_that("compute_smoothed: 'none' is the identity transform", {
  ds <- mock_dataset(n_cells = 80, seed = 3)
  out <- compute_smoothed(ds, genes = "CD3D", method = "none")
  expect_identical(out$method, "none")
  expect_identical(length(out$expression$CD3D), 80L)
  expect_identical(out$expression$CD3D, get_gene_expression(ds, "CD3D"))
})

test_that("compute_smoothed: 'neighbor' produces a different vector than raw", {
  ds <- mock_dataset(n_cells = 80, seed = 4)
  out <- compute_smoothed(ds, genes = "CD3D",
                          method = "neighbor", k = 5, reduction = "UMAP")
  raw <- get_gene_expression(ds, "CD3D")
  expect_identical(length(out$expression$CD3D), 80L)
  expect_false(identical(out$expression$CD3D, raw))
})

test_that("compute_smoothed: 'alra_mock' thresholds many values to zero", {
  ds <- mock_dataset(n_cells = 80, seed = 5)
  out <- compute_smoothed(ds, genes = "CD3D",
                          method = "alra_mock", k = 5, reduction = "UMAP")
  expect_true(any(out$expression$CD3D == 0))
  expect_true(mean(out$expression$CD3D == 0) >= 0.1)
})

test_that("get_gene_expression_for_view enforces the visualization-only contract", {
  state <- new_app_state()
  ds <- mock_dataset(n_cells = 60, seed = 8)
  with_state(state, set_active_dataset(state, ds))

  raw <- get_gene_expression(ds, "CD3D")
  view_raw <- with_state(state, get_gene_expression_for_view(state, "CD3D"))
  expect_identical(view_raw, raw)

  out <- compute_smoothed(ds, genes = "CD3D",
                          method = "neighbor", k = 5, reduction = "UMAP")
  with_state(state, {
    state$analysis_results <- list(imputation = list(
      status = "completed", results = out, params = list()))
    state$display_mode_imputation <- "raw"
  })
  expect_identical(with_state(state, get_gene_expression_for_view(state, "CD3D")),
                   raw)

  with_state(state, state$display_mode_imputation <- "smoothed")
  view_smooth <- with_state(state, get_gene_expression_for_view(state, "CD3D"))
  expect_false(identical(view_smooth, raw))
  expect_identical(view_smooth, out$expression$CD3D)
})

test_that("has_smoothed_results requires completed status + non-empty expression", {
  state <- new_app_state()
  with_state(state, set_active_dataset(state, mock_dataset(n_cells = 40)))
  expect_false(with_state(state, has_smoothed_results(state)))

  with_state(state, {
    state$analysis_results <- list(imputation = list(
      status = "completed",
      results = list(expression = list(CD3D = rep(1, 40)))))
  })
  expect_true(with_state(state, has_smoothed_results(state)))
})
