# Tests for R/annotation_azimuth.R
# ----------------------------------------------------------------------------
# Pattern mirrors test-singler.R: registry + missing-dep + pure
# converter against hand-built Azimuth-style data.frames, with an
# end-to-end test that is skipped without Azimuth installed.

# ---- Registry entry --------------------------------------------------------

test_that("azimuth engine is registered with the documented spec", {
  e <- get_annotation_engine("azimuth")
  expect_false(is.null(e))
  expect_identical(e$id,       "azimuth")
  expect_identical(e$category, "reference-based")
  expect_true("dataset" %in% e$requires)
  expect_true(all(c("per_cell_labels", "per_cell_scores",
                    "cluster_summary") %in% e$produces))
  expect_true(is.function(e$run_fn))
  # Parameters surface the four documented controls.
  expect_setequal(names(e$parameters),
                  c("reference", "annotation_level",
                    "cluster_field", "min_mapping_score"))
  expect_identical(e$parameters$reference$default,        "pbmcref")
  expect_identical(e$parameters$annotation_level$default, "celltype.l2")
})

test_that("azimuth is listed in list_annotation_engines()", {
  ids <- unname(list_annotation_engines())
  expect_true("azimuth" %in% ids)
})

# ---- ensure_attached helper ------------------------------------------------
# General-purpose; not Azimuth-specific (Azimuth uses a separate
# runtime-setters shim).

test_that("ensure_attached() is idempotent when the package is already attached", {
  expect_true("package:methods" %in% search())
  expect_true(ensure_attached("methods"))
  expect_true(ensure_attached("methods"))
  expect_equal(sum(search() == "package:methods"), 1L)
})

test_that("ensure_attached() raises a clear error for a nonexistent package", {
  expect_error(ensure_attached("definitelynotapackage_xyz"),
               regexp = "cannot attach package")
})

test_that("ensure_attached() rejects malformed input", {
  expect_error(ensure_attached(character()))
  expect_error(ensure_attached(c("a", "b")))
  expect_error(ensure_attached(""))
})

# ---- Azimuth runtime-setters shim ------------------------------------------
# These cover the specific "could not find function `Key<-`" failure
# mode without needing SeuratObject / Azimuth installed: the helper
# behaves correctly when SeuratObject is missing (returns empty), is
# idempotent on repeat calls, and ships the expected default symbol
# list.

test_that("Azimuth shim short-circuits when SeuratObject is missing", {
  skip_if(requireNamespace("SeuratObject", quietly = TRUE),
          "SeuratObject installed; missing-pkg branch not exercised here")
  expect_identical(.azimuth_install_runtime_setters(), character())
  expect_false(.AZIMUTH_SHIM_NAME %in% search())
})

test_that("Azimuth shim symbol list covers the documented setter generics", {
  expect_true("Key<-" %in% .AZIMUTH_SHIM_SYMBOLS)
  expect_true("DefaultAssay<-" %in% .AZIMUTH_SHIM_SYMBOLS)
  expect_true("Idents<-" %in% .AZIMUTH_SHIM_SYMBOLS)
})

test_that("Azimuth shim is idempotent when already on the search path", {
  skip_if_not(requireNamespace("SeuratObject", quietly = TRUE),
              "SeuratObject not installed")
  # First call: install. Second call: no-op (returns character(0)).
  on.exit(if (.AZIMUTH_SHIM_NAME %in% search())
            base::detach(.AZIMUTH_SHIM_NAME, character.only = TRUE),
          add = TRUE)
  patched <- .azimuth_install_runtime_setters()
  expect_true(length(patched) >= 1L)
  expect_true(.AZIMUTH_SHIM_NAME %in% search())
  # Second invocation must be a no-op.
  expect_identical(.azimuth_install_runtime_setters(), character())
})

test_that("Azimuth shim actually makes the setter visible to a namespaced fn", {
  skip_if_not(requireNamespace("SeuratObject", quietly = TRUE),
              "SeuratObject not installed")
  on.exit(if (.AZIMUTH_SHIM_NAME %in% search())
            base::detach(.AZIMUTH_SHIM_NAME, character.only = TRUE),
          add = TRUE)
  .azimuth_install_runtime_setters()
  # Probe: pretend we're a function in `stats` (an arbitrary
  # package). `get(..., inherits = TRUE)` walks the same chain that
  # R uses to resolve symbols, so a TRUE result here proves the shim
  # closes the lookup gap that Azimuth hits.
  probe <- function() NULL
  environment(probe) <- asNamespace("stats")
  found <- tryCatch({
    get("Key<-", envir = environment(probe), inherits = TRUE); TRUE
  }, error = function(e) FALSE)
  expect_true(found)
})

test_that(".azimuth_diagnose_setters() returns the expected structured shape", {
  d <- .azimuth_diagnose_setters()
  expect_true(is.list(d))
  expect_true(all(c("azimuth_installed", "seuratobject_installed",
                    "seurat_installed", "shim_attached", "resolves") %in% names(d)))
  expect_true(is.logical(d$shim_attached))
  expect_true(is.logical(d$resolves))
  expect_identical(names(d$resolves), .AZIMUTH_SHIM_SYMBOLS)
})

# ---- Missing-dep gating ----------------------------------------------------

test_that("run_annotation_engine('azimuth') errors cleanly without Azimuth", {
  skip_if(has_optional("Azimuth"),
          "Azimuth installed; install-error path skipped")
  res <- mock_state_with_dataset(n_cells = 60)
  state <- res$state; ds <- res$dataset
  expect_error(
    with_state(state, run_annotation_engine(
      "azimuth", ds, state,
      params = list(reference = "pbmcref",
                    annotation_level = "celltype.l2"),
      set_id = "az_missing")),
    regexp = "Azimuth", fixed = FALSE
  )
})

# ---- Pure converter --------------------------------------------------------

test_that(".azimuth_to_engine_output: schema + label/score alignment", {
  cells <- sprintf("c_%02d", 1:6)
  df <- data.frame(
    cell                          = cells,
    predicted.celltype.l2         = c("CD4 T", "CD4 T", "CD8 T",
                                      "B naive", "NK", "Mono"),
    predicted.celltype.l2.score   = c(0.95, 0.92, 0.88, 0.7, 0.6, 0.5),
    mapping.score                 = c(0.8, 0.8, 0.7, 0.6, 0.5, 0.4),
    stringsAsFactors = FALSE
  )
  out <- .azimuth_to_engine_output(df, cells = cells,
                                   annotation_level = "celltype.l2",
                                   reference        = "pbmcref")
  expect_identical(out$cell,        cells)
  expect_identical(out$cell_labels, c("CD4 T", "CD4 T", "CD8 T",
                                      "B naive", "NK", "Mono"))
  expect_identical(out$cell_scores,
                   c(0.95, 0.92, 0.88, 0.7, 0.6, 0.5))
  expect_identical(out$reference_source, "Azimuth:pbmcref")
})

test_that(".azimuth_to_engine_output: min_mapping_score gates low-confidence", {
  cells <- c("a", "b", "c")
  df <- data.frame(
    cell                       = cells,
    predicted.celltype.l2      = c("X", "Y", "Z"),
    predicted.celltype.l2.score = c(0.9, 0.9, 0.9),
    mapping.score              = c(0.9, 0.6, 0.3),
    stringsAsFactors = FALSE
  )
  out <- .azimuth_to_engine_output(df, cells = cells,
                                   annotation_level = "celltype.l2",
                                   min_mapping_score = 0.7,
                                   reference = "pbmcref")
  expect_identical(out$cell_labels, c("X", "Unknown", "Unknown"))
  expect_identical(out$cell_scores, c(0.9, 0, 0))
})

test_that(".azimuth_to_engine_output: missing score column -> NA scores", {
  cells <- c("a", "b")
  df <- data.frame(
    cell                  = cells,
    predicted.celltype.l1 = c("T", "B"),
    stringsAsFactors = FALSE
  )
  out <- .azimuth_to_engine_output(df, cells = cells,
                                   annotation_level = "celltype.l1",
                                   reference = "pbmcref")
  expect_true(all(is.na(out$cell_scores)))
  expect_identical(out$cell_labels, c("T", "B"))
})

test_that(".azimuth_to_engine_output: row reordering via match()", {
  cells <- c("a", "b", "c")
  # df rows out of order on purpose.
  df <- data.frame(
    cell                        = c("c", "a", "b"),
    predicted.celltype.l2       = c("Z", "X", "Y"),
    predicted.celltype.l2.score = c(0.6, 0.95, 0.8),
    stringsAsFactors = FALSE
  )
  out <- .azimuth_to_engine_output(df, cells = cells,
                                   annotation_level = "celltype.l2",
                                   reference = "pbmcref")
  expect_identical(out$cell_labels, c("X", "Y", "Z"))
  expect_identical(out$cell_scores, c(0.95, 0.8, 0.6))
})

test_that(".azimuth_to_engine_output: errors on missing cells", {
  cells <- c("a", "b", "c")
  df <- data.frame(cell = c("a", "b"),
                   predicted.celltype.l2 = c("X", "Y"),
                   stringsAsFactors = FALSE)
  expect_error(.azimuth_to_engine_output(df, cells = cells,
                                         annotation_level = "celltype.l2"),
               regexp = "missing from Azimuth result", fixed = TRUE)
})

test_that(".azimuth_to_engine_output: errors on missing label column", {
  df <- data.frame(cell = "a",
                   predicted.celltype.l1 = "T",
                   stringsAsFactors = FALSE)
  expect_error(.azimuth_to_engine_output(df, cells = "a",
                                         annotation_level = "celltype.l2"),
               regexp = "predicted.celltype.l2", fixed = TRUE)
})

test_that(".azimuth_to_engine_output: errors on missing `cell` column", {
  df <- data.frame(predicted.celltype.l2 = "X", stringsAsFactors = FALSE)
  expect_error(.azimuth_to_engine_output(df, cells = "a"),
               regexp = "needs a `cell` column", fixed = TRUE)
})

test_that(".azimuth_to_engine_output: cluster_summary when cluster_vec given", {
  cells <- sprintf("c_%02d", 1:6)
  df <- data.frame(
    cell                        = cells,
    predicted.celltype.l2       = c("T", "T", "T", "B", "B", "NK"),
    predicted.celltype.l2.score = c(0.9, 0.85, 0.8, 0.7, 0.6, 0.5),
    mapping.score               = rep(0.8, 6),
    stringsAsFactors = FALSE
  )
  cluster_vec <- c("0", "0", "0", "1", "1", "1")
  out <- .azimuth_to_engine_output(df, cells = cells,
                                   annotation_level   = "celltype.l2",
                                   reference          = "pbmcref",
                                   cluster_field_used = "cluster",
                                   cluster_vec        = cluster_vec)
  expect_s3_class(out$cluster_summary, "data.frame")
  expect_identical(nrow(out$cluster_summary), 2L)
  cs <- out$cluster_summary
  # Cluster 0: top T (3/3 = 1.0), Cluster 1: top B (2/3) and NK (1/3) -> top B
  expect_identical(cs$top_label[cs$cluster == "0"], "T")
  expect_identical(cs$top_label[cs$cluster == "1"], "B")
  expect_equal(cs$top_fraction[cs$cluster == "0"], 1.0)
  expect_equal(cs$top_fraction[cs$cluster == "1"], 2 / 3, tolerance = 1e-9)
  expect_identical(out$n_clusters_at_creation, 2L)
  expect_identical(out$cluster_field_used,      "cluster")
})

# ---- End-to-end (skipped without Azimuth) -------------------------------

test_that("run_annotation_engine('azimuth') round-trips against a real reference", {
  skip_if_not_installed("Azimuth")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  # Azimuth reference data packages are heavy (>1GB) and not safe to
  # download in CI, so this test is gated on Azimuth being installed.
  res <- mock_state_with_dataset(n_cells = 200)
  state <- res$state; ds <- res$dataset
  out <- with_state(state, run_annotation_engine(
    "azimuth", ds, state,
    params = list(reference        = "pbmcref",
                  annotation_level = "celltype.l2"),
    set_id = "azimuth_smoke"))
  expect_true(is_annotation_result_v1(out))
  expect_identical(out$engine_id, "azimuth")
})
