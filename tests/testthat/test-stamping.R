# =========================================================================
# Cross-cutting regression: every analysis module stamps its results with
# `annotation_stamp = make_annotation_stamp(state)`.
#
# The shape of this contract was the most expensive piece of the annotation
# refactor to land, so we guard it two ways:
#   1. Static check: each module file contains the stamping call. Cheap, and
#      catches accidental deletions in PRs.
#   2. End-to-end shiny::testServer for the DE module, which actually runs
#      `compute_de` against the mock dataset and inspects
#      `state$analysis_results$de$annotation_stamp`.
# =========================================================================

.proj_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
                            mustWork = TRUE)
.module_path <- function(name) {
  file.path(.proj_root, "R", "modules", paste0(name, ".R"))
}

.read_module <- function(name) {
  paste(readLines(.module_path(name), warn = FALSE), collapse = "\n")
}

.stamping_modules <- c(
  "mod_differential_expression",
  "mod_pathway_analysis",
  "mod_imputation",
  "mod_trajectory",
  "mod_marker_investigation"
)

test_that("each analysis module stamps results with make_annotation_stamp", {
  for (m in .stamping_modules) {
    src <- .read_module(m)
    expect_true(
      grepl("annotation_stamp\\s*=\\s*make_annotation_stamp\\(state\\)", src),
      info = sprintf("module %s is missing annotation_stamp", m)
    )
  }
})

test_that("DE/imputation/pathway/trajectory stamp all three lifecycle slots", {
  # These four use the running / failed / completed pattern.
  for (m in c("mod_differential_expression", "mod_pathway_analysis",
              "mod_imputation", "mod_trajectory")) {
    src <- .read_module(m)
    n_stamps <- length(gregexpr("make_annotation_stamp\\(state\\)", src,
                                perl = TRUE)[[1]])
    expect_gte(n_stamps, 3L)
  }
})

test_that("empty_analysis_result reserves an annotation_stamp slot", {
  e <- empty_analysis_result()
  expect_true("annotation_stamp" %in% names(e))
  expect_null(e$annotation_stamp)
})

# ---- End-to-end: DE module actually stamps via testServer ----------------

test_that("DE module writes annotation_stamp into state$analysis_results$de", {
  skip_if_not_installed("shiny")
  state <- new_app_state()
  ds <- mock_dataset(n_cells = 150, seed = 17)
  shiny::isolate(set_active_dataset(state, ds))

  # Activate an annotation set so the stamp carries real provenance.
  set <- shiny::isolate(run_annotation_engine(
    "manual", ds, state,
    params = list(cluster_field = "cluster",
                  labels = list("0" = "T cell", "1" = "B cell")),
    set_id = "stamp_de"))
  shiny::isolate({
    add_annotation_set(state, set)
    set_active_annotation(state, "stamp_de")
  })

  shiny::testServer(
    mod_differential_expression_server,
    args = list(state = state),
    {
      session$setInputs(
        group_field = "cluster", group_1 = "0", group_2 = "1",
        assay = ds$default_assay, layer = "data",
        min_pct = 0
      )
      session$setInputs(run = 1)

      de <- state$analysis_results$de
      expect_true(!is.null(de))
      expect_true(!is.null(de$annotation_stamp))
      expect_identical(de$annotation_stamp$annotation_set_id_used, "stamp_de")
      expect_identical(de$annotation_stamp$annotation_engine_id, "manual")
      expect_identical(de$annotation_stamp$annotation_set_hash_used,
                       annotation_set_hash(get_active_annotation(state)))
    }
  )
})
