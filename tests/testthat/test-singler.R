# =========================================================================
# SingleR annotation engine
# Pure-converter tests run on every host; the end-to-end test that
# actually calls SingleR is skipped when the Bioconductor packages are
# missing (the usual case in CI). The converter has the schema-mapping
# logic, so a green converter suite is enough to guarantee schema
# compatibility -- the live SingleR path just plumbs predictions through.
# =========================================================================

# ---- Registry ------------------------------------------------------------

test_that("singler engine is registered in ANNOTATION_ENGINES", {
  spec <- get_annotation_engine("singler")
  expect_false(is.null(spec))
  expect_identical(spec$id, "singler")
  expect_identical(spec$category, "reference-based")
  expect_true("dataset" %in% spec$requires)
  expect_setequal(spec$produces,
                  c("per_cell_labels", "per_cell_scores",
                    "alt_labels", "cluster_summary"))
  # The user-facing dropdown should list the engine.
  ids <- list_annotation_engines()
  expect_true("singler" %in% ids)
})

test_that("singler engine parameter spec covers reference / labels / cluster_field / min_delta", {
  spec <- get_annotation_engine("singler")
  expect_setequal(names(spec$parameters),
                  c("reference", "labels", "cluster_field", "min_delta"))
  ref_choices <- spec$parameters$reference$choices
  expect_setequal(ref_choices, c("hpca", "blueprint_encode", "monaco_immune"))
  lbl_choices <- spec$parameters$labels$choices
  expect_setequal(lbl_choices, c("main", "fine"))
})

# ---- Missing-dep behaviour ----------------------------------------------

test_that("running the singler engine without Bioconductor deps yields a clear install message", {
  if (has_optional(c("SingleR", "celldex", "SingleCellExperiment", "SummarizedExperiment"))) {
    skip("SingleR + celldex are installed; missing-dep path is not reachable.")
  }
  setup <- mock_state_with_dataset(n_cells = 60)
  state <- setup$state; ds <- setup$dataset
  err <- tryCatch(
    with_state(state, run_annotation_engine(
      "singler", ds, state,
      params = list(reference = "hpca"),
      set_id = "sr_missing")),
    error = function(e) conditionMessage(e))
  expect_match(err, "SingleR")
  expect_match(err, "BiocManager::install")
})

test_that(".fetch_celldex_reference errors clearly on unknown reference names", {
  if (!has_optional("celldex")) {
    # Without celldex installed the function errors at the require_optional
    # gate; we exercise the gate here so the friendly message is verified.
    err <- tryCatch(.fetch_celldex_reference("hpca"),
                    error = function(e) conditionMessage(e))
    expect_match(err, "celldex")
    return(invisible())
  }
  # With celldex available we hit the unknown-reference branch instead.
  expect_error(.fetch_celldex_reference("not_a_ref"),
               "unknown reference")
})

# ---- Pure converter: per-cell mode --------------------------------------

# Build a fake SingleR-style data.frame. `scores` is held as a list-column
# so each cell carries a numeric vector across reference labels (the same
# shape `SingleR::SingleR()` ships as a one-row-per-cell scores matrix).
.make_singler_per_cell <- function(cells,
                                   labels,
                                   pruned = labels,
                                   delta  = rep(0.5, length(cells)),
                                   scores = NULL,
                                   ref_labels = c("T_cell", "B_cell", "NK_cell")) {
  if (is.null(scores)) {
    m <- matrix(stats::runif(length(cells) * length(ref_labels)),
                nrow = length(cells), ncol = length(ref_labels),
                dimnames = list(cells, ref_labels))
  } else {
    m <- scores
  }
  data.frame(
    labels         = labels,
    pruned.labels  = pruned,
    delta.next     = delta,
    scores         = I(m),
    row.names      = cells,
    stringsAsFactors = FALSE)
}

test_that(".singler_to_engine_output (per-cell) returns the engine output schema", {
  cells <- paste0("c", 1:5)
  pred  <- .make_singler_per_cell(
    cells  = cells,
    labels = c("T_cell", "B_cell", "T_cell", "NK_cell", "B_cell"),
    delta  = c(0.8, 0.6, 0.7, 0.9, 0.5))
  out <- .singler_to_engine_output(pred, cells = cells,
                                   reference_source = "celldex/hpca:main")

  expect_setequal(names(out),
                  c("cell", "cell_labels", "cell_scores", "alt_labels",
                    "cluster_summary", "cluster_field_used",
                    "n_clusters_at_creation", "reference_source", "warnings"))
  expect_identical(out$cell, cells)
  expect_identical(out$cell_labels,
                   c("T_cell", "B_cell", "T_cell", "NK_cell", "B_cell"))
  expect_identical(out$cell_scores, c(0.8, 0.6, 0.7, 0.9, 0.5))
  expect_identical(out$reference_source, "celldex/hpca:main")
  expect_true(is.na(out$cluster_field_used))
  expect_true(is.na(out$n_clusters_at_creation))
  expect_null(out$cluster_summary)
  expect_identical(out$warnings, character())
})

test_that(".singler_to_engine_output respects pruned.labels (NA -> 'Unknown', score 0)", {
  cells <- paste0("c", 1:4)
  pred  <- .make_singler_per_cell(
    cells  = cells,
    labels = c("T_cell", "B_cell", "T_cell", "NK_cell"),
    pruned = c("T_cell", NA,        "T_cell", "NK_cell"),
    delta  = c(0.8,      0.05,      0.7,       0.9))
  out <- .singler_to_engine_output(pred, cells = cells)
  expect_identical(out$cell_labels[2], "Unknown")
  expect_identical(out$cell_scores[2], 0.0)
  # Other cells unchanged
  expect_identical(out$cell_labels[c(1, 3, 4)],
                   c("T_cell", "T_cell", "NK_cell"))
})

test_that(".singler_to_engine_output applies min_delta cutoff", {
  cells <- paste0("c", 1:4)
  pred  <- .make_singler_per_cell(
    cells = cells,
    labels = c("T_cell", "B_cell", "T_cell", "NK_cell"),
    delta  = c(0.8, 0.10, 0.7, 0.15))   # cells 2, 4 below threshold
  out <- .singler_to_engine_output(pred, cells = cells, min_delta = 0.2)
  expect_identical(out$cell_labels,
                   c("T_cell", "Unknown", "T_cell", "Unknown"))
  expect_identical(out$cell_scores, c(0.8, 0.0, 0.7, 0.0))
})

test_that(".singler_to_engine_output populates alt_labels from the scores matrix", {
  cells <- paste0("c", 1:3)
  ref_labels <- c("T_cell", "B_cell", "NK_cell")
  scores_mat <- rbind(
    c(0.90, 0.30, 0.10),    # c1: T > B > NK
    c(0.20, 0.85, 0.10),    # c2: B > T > NK
    c(0.10, 0.20, 0.75)     # c3: NK > B > T
  )
  dimnames(scores_mat) <- list(cells, ref_labels)
  pred <- .make_singler_per_cell(
    cells  = cells,
    labels = c("T_cell", "B_cell", "NK_cell"),
    scores = scores_mat,
    ref_labels = ref_labels)

  out <- .singler_to_engine_output(pred, cells = cells)
  expect_s3_class(out$alt_labels, "data.frame")
  # 3 cells x 3 ranks = 9 rows
  expect_identical(nrow(out$alt_labels), 9L)
  expect_setequal(names(out$alt_labels), c("cell", "rank", "label", "score"))
  # Top-1 per cell matches the input labels
  top1 <- out$alt_labels[out$alt_labels$rank == 1, ]
  expect_identical(top1$cell,  cells)
  expect_identical(top1$label, c("T_cell", "B_cell", "NK_cell"))
})

test_that(".singler_to_engine_output produces a cluster_summary when cluster_vec is given", {
  cells <- paste0("c", 1:6)
  pred  <- .make_singler_per_cell(
    cells  = cells,
    labels = c("T_cell", "T_cell", "B_cell", "B_cell", "NK_cell", "NK_cell"),
    delta  = c(0.8, 0.7, 0.6, 0.5, 0.9, 0.85))
  cluster_vec <- c("0", "0", "1", "1", "2", "2")
  out <- .singler_to_engine_output(pred, cells = cells,
                                   cluster_vec = cluster_vec,
                                   cluster_field_used = "cluster")

  expect_identical(out$cluster_field_used, "cluster")
  expect_identical(out$n_clusters_at_creation, 3L)
  cs <- out$cluster_summary
  expect_s3_class(cs, "data.frame")
  expect_identical(cs$cluster,   c("0", "1", "2"))
  expect_identical(cs$top_label, c("T_cell", "B_cell", "NK_cell"))
  expect_identical(cs$n_cells,   c(2L, 2L, 2L))
})

test_that(".singler_to_engine_output handles an empty prediction", {
  out <- .singler_to_engine_output(NULL, cells = paste0("c", 1:3))
  expect_identical(out$cell_labels, rep(NA_character_, 3))
  expect_identical(out$cell_scores, rep(0.0,           3))
  expect_match(out$warnings, "empty prediction")
})

test_that(".singler_to_engine_output errors on row-count mismatch", {
  pred <- .make_singler_per_cell(cells = paste0("c", 1:4),
                                 labels = rep("T_cell", 4))
  expect_error(.singler_to_engine_output(pred, cells = paste0("c", 1:3)),
               "pred has 4 rows but cells has 3")
})

# ---- Pure converter: cluster mode ---------------------------------------

test_that(".singler_cluster_to_cells expands cluster-level pred back to per-cell", {
  ds <- mock_dataset(n_cells = 12, seed = 91)
  # Force a known cluster vector for test stability
  ds$cell_data$cluster <- as.character(rep(0:2, each = 4))
  cluster_vec <- ds$cell_data$cluster

  # SingleR cluster output: one row per unique cluster id
  cluster_ids <- c("0", "1", "2")
  ref_labels  <- c("T_cell", "B_cell", "NK_cell")
  scores_mat  <- diag(3); dimnames(scores_mat) <- list(cluster_ids, ref_labels)
  pred_cluster <- data.frame(
    labels        = c("T_cell", "B_cell", "NK_cell"),
    pruned.labels = c("T_cell", "B_cell", "NK_cell"),
    delta.next    = c(0.9,      0.8,      0.7),
    scores        = I(scores_mat),
    row.names     = cluster_ids,
    stringsAsFactors = FALSE)

  out <- .singler_cluster_to_cells(pred_cluster, dataset = ds,
                                   cluster_vec = cluster_vec,
                                   cluster_field_used = "cluster",
                                   reference_source = "celldex/hpca:main")
  expect_identical(out$cell, ds$cell_data$cell)
  expect_identical(out$cell_labels,
                   rep(c("T_cell", "B_cell", "NK_cell"), each = 4))
  expect_identical(out$cluster_field_used, "cluster")
  expect_identical(out$n_clusters_at_creation, 3L)
  cs <- out$cluster_summary
  expect_identical(cs$top_label, c("T_cell", "B_cell", "NK_cell"))
})

test_that(".singler_cluster_to_cells errors when pred rownames don't carry cluster ids", {
  # data.frame() always assigns numeric rownames ("1", "2"...). If the
  # caller didn't reset them to the cluster ids, none of the cluster_vec
  # entries will match -> friendly error. Use non-numeric cluster ids
  # to guarantee no accidental overlap with default rownames.
  ds <- mock_dataset(n_cells = 4, seed = 1)
  pred <- data.frame(labels = c("A", "B"), stringsAsFactors = FALSE)
  expect_error(
    .singler_cluster_to_cells(pred, dataset = ds,
                              cluster_vec = c("alpha", "alpha", "beta", "beta"),
                              cluster_field_used = "cluster",
                              reference_source = "x"),
    "no cluster-id rownames"
  )
})

# ---- End-to-end (skipped when SingleR + celldex missing) ----------------

test_that("run_annotation_engine('singler') round-trips against celldex", {
  skip_if_not_installed("SingleR")
  skip_if_not_installed("celldex")
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if(Sys.getenv("SCRNA_EXPLORER_SKIP_CELLDEX_DOWNLOAD", "") == "1",
          "SCRNA_EXPLORER_SKIP_CELLDEX_DOWNLOAD=1; skipping reference download")

  setup <- mock_state_with_dataset(n_cells = 80, name = "singler_e2e")
  state <- setup$state; ds <- setup$dataset
  res <- with_state(state, run_annotation_engine(
    "singler", ds, state,
    params = list(reference = "hpca", labels = "main",
                  cluster_field = "cluster"),
    set_id = "sr_e2e"))
  expect_true(is_annotation_result_v1(res))
  expect_identical(res$engine_id, "singler")
  expect_identical(length(res$cell), 80L)
  expect_identical(length(res$cell_labels), 80L)
  expect_match(res$reference_source, "^celldex/hpca:main")
})
