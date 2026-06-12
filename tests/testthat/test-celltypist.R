# Tests for R/annotation_celltypist.R
# ----------------------------------------------------------------------------
# Pattern mirrors test-azimuth.R / test-singler.R: registry + missing-dep
# + pure converter against hand-built CellTypist-style data.frames. The
# end-to-end test is skipped unless reticulate + the Python celltypist
# module are both available.

# ---- Registry entry --------------------------------------------------------

test_that("celltypist engine is registered with the documented spec", {
  e <- get_annotation_engine("celltypist")
  expect_false(is.null(e))
  expect_identical(e$id,       "celltypist")
  expect_identical(e$category, "reference-based")
  expect_true("dataset"    %in% e$requires)
  expect_true("py_runtime" %in% e$requires)
  expect_true(all(c("per_cell_labels", "per_cell_scores",
                    "cluster_summary") %in% e$produces))
  expect_setequal(names(e$parameters),
                  c("model", "majority_voting", "over_clustering",
                    "cluster_field", "min_score"))
  expect_identical(e$parameters$model$default,           "Immune_All_Low.pkl")
  expect_identical(e$parameters$majority_voting$default, FALSE)
})

test_that("celltypist is listed in list_annotation_engines()", {
  ids <- unname(list_annotation_engines())
  expect_true("celltypist" %in% ids)
})

# ---- Missing-dep gating ----------------------------------------------------

test_that("run_annotation_engine('celltypist') errors cleanly without reticulate", {
  skip_if(has_optional("reticulate"),
          "reticulate installed; install-error path skipped")
  res <- mock_state_with_dataset(n_cells = 60)
  state <- res$state; ds <- res$dataset
  expect_error(
    with_state(state, run_annotation_engine(
      "celltypist", ds, state,
      params = list(model = "Immune_All_Low.pkl"),
      set_id = "ct_missing_reticulate")),
    regexp = "reticulate", fixed = FALSE
  )
})

test_that(
  paste("run_annotation_engine('celltypist') errors when Python celltypist",
        "is missing"), {
  skip_if_not_installed("reticulate")
  skip_if_not_installed("anndata")
  if (isTRUE(reticulate::py_module_available("celltypist"))) {
    skip("celltypist Python module is available; missing-Python-dep path skipped")
  }
  res <- mock_state_with_dataset(n_cells = 60)
  state <- res$state; ds <- res$dataset
  expect_error(
    with_state(state, run_annotation_engine(
      "celltypist", ds, state,
      params = list(model = "Immune_All_Low.pkl"),
      set_id = "ct_missing_python")),
    regexp = "celltypist", fixed = FALSE
  )
})

# ---- Pure converter --------------------------------------------------------

test_that(".celltypist_to_engine_output: default mode (predicted_labels + conf_score)", {
  cells <- sprintf("c_%02d", 1:5)
  df <- data.frame(
    cell             = cells,
    predicted_labels = c("CD4 T", "CD8 T", "B", "NK", "Mono"),
    conf_score       = c(0.9, 0.85, 0.7, 0.6, 0.5),
    stringsAsFactors = FALSE
  )
  out <- .celltypist_to_engine_output(df, cells = cells,
                                      model = "Immune_All_Low.pkl")
  expect_identical(out$cell_labels,
                   c("CD4 T", "CD8 T", "B", "NK", "Mono"))
  expect_identical(out$cell_scores, c(0.9, 0.85, 0.7, 0.6, 0.5))
  expect_identical(out$reference_source,
                   "CellTypist:Immune_All_Low.pkl")
})

test_that(".celltypist_to_engine_output: majority_voting selects the voted label", {
  cells <- c("a", "b", "c")
  df <- data.frame(
    cell             = cells,
    predicted_labels = c("T", "T", "B"),
    majority_voting  = c("T cell", "T cell", "T cell"),
    conf_score       = c(0.9, 0.9, 0.7),
    over_clustering  = c("0", "0", "0"),
    stringsAsFactors = FALSE
  )
  out <- .celltypist_to_engine_output(df, cells = cells,
                                      model = "Immune_All_Low.pkl",
                                      majority_voting = TRUE)
  # With majority voting on, the majority_voting column wins for every cell.
  expect_identical(out$cell_labels, rep("T cell", 3))
  expect_identical(out$reference_source,
                   "CellTypist:Immune_All_Low.pkl:majority_voting")
})

test_that(".celltypist_to_engine_output: majority_voting=TRUE but no column -> falls back", {
  cells <- c("a", "b")
  df <- data.frame(
    cell             = cells,
    predicted_labels = c("X", "Y"),
    conf_score       = c(0.9, 0.8),
    stringsAsFactors = FALSE
  )
  out <- .celltypist_to_engine_output(df, cells = cells,
                                      majority_voting = TRUE)
  expect_identical(out$cell_labels, c("X", "Y"))
})

test_that(".celltypist_to_engine_output: min_score gates low-confidence cells", {
  cells <- c("a", "b", "c")
  df <- data.frame(
    cell             = cells,
    predicted_labels = c("X", "Y", "Z"),
    conf_score       = c(0.9, 0.5, 0.2),
    stringsAsFactors = FALSE
  )
  out <- .celltypist_to_engine_output(df, cells = cells, min_score = 0.6)
  expect_identical(out$cell_labels, c("X", "Unknown", "Unknown"))
  expect_identical(out$cell_scores, c(0.9, 0, 0))
})

test_that(".celltypist_to_engine_output: missing conf_score -> NA scores", {
  cells <- c("a", "b")
  df <- data.frame(
    cell             = cells,
    predicted_labels = c("X", "Y"),
    stringsAsFactors = FALSE
  )
  out <- .celltypist_to_engine_output(df, cells = cells)
  expect_true(all(is.na(out$cell_scores)))
  # min_score is silently ignored when conf_score is absent (documented).
  out2 <- .celltypist_to_engine_output(df, cells = cells, min_score = 0.9)
  expect_identical(out2$cell_labels, c("X", "Y"))
})

test_that(".celltypist_to_engine_output: row reordering via match()", {
  cells <- c("a", "b", "c")
  df <- data.frame(
    cell             = c("c", "a", "b"),
    predicted_labels = c("Z", "X", "Y"),
    conf_score       = c(0.6, 0.95, 0.8),
    stringsAsFactors = FALSE
  )
  out <- .celltypist_to_engine_output(df, cells = cells)
  expect_identical(out$cell_labels, c("X", "Y", "Z"))
  expect_identical(out$cell_scores, c(0.95, 0.8, 0.6))
})

test_that(".celltypist_to_engine_output: errors on missing cells", {
  expect_error(
    .celltypist_to_engine_output(
      data.frame(cell = "a", predicted_labels = "X",
                 stringsAsFactors = FALSE),
      cells = c("a", "b")),
    regexp = "missing from CellTypist result", fixed = TRUE
  )
})

test_that(".celltypist_to_engine_output: errors when no label column", {
  expect_error(
    .celltypist_to_engine_output(
      data.frame(cell = "a", conf_score = 0.9, stringsAsFactors = FALSE),
      cells = "a"),
    regexp = "predicted_labels", fixed = FALSE
  )
})

test_that(".celltypist_to_engine_output: cluster_summary when cluster_vec given", {
  cells <- sprintf("c_%02d", 1:6)
  df <- data.frame(
    cell             = cells,
    predicted_labels = c("T", "T", "T", "B", "B", "NK"),
    conf_score       = c(0.9, 0.85, 0.8, 0.7, 0.6, 0.5),
    stringsAsFactors = FALSE
  )
  cluster_vec <- c("0", "0", "0", "1", "1", "1")
  out <- .celltypist_to_engine_output(
    df, cells = cells, model = "Immune_All_Low.pkl",
    cluster_field_used = "cluster", cluster_vec = cluster_vec)
  expect_s3_class(out$cluster_summary, "data.frame")
  cs <- out$cluster_summary
  expect_identical(cs$top_label[cs$cluster == "0"], "T")
  expect_identical(cs$top_label[cs$cluster == "1"], "B")
  expect_identical(out$cluster_field_used, "cluster")
})

# ---- End-to-end (skipped without reticulate + celltypist) ----------------

test_that("run_annotation_engine('celltypist') round-trips against a real Python install", {
  skip_if_not_installed("reticulate")
  skip_if_not_installed("anndata")
  if (!isTRUE(reticulate::py_module_available("celltypist"))) {
    skip("celltypist Python module not available")
  }
  res <- mock_state_with_dataset(n_cells = 200)
  state <- res$state; ds <- res$dataset
  out <- with_state(state, run_annotation_engine(
    "celltypist", ds, state,
    params = list(model = "Immune_All_Low.pkl",
                  majority_voting = FALSE),
    set_id = "celltypist_smoke"))
  expect_true(is_annotation_result_v1(out))
  expect_identical(out$engine_id, "celltypist")
})
