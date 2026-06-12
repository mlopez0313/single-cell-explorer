# =========================================================================
# Annotation system: schema + engines + multi-set lifecycle + apply
# Mirrors the smoke suite that drove the architectural refactor; rewritten
# as proper testthat for ongoing regression coverage.
# =========================================================================

# ---- Schema --------------------------------------------------------------

test_that("annotation_result_v1 builds a full-shape object", {
  ds <- mock_dataset(n_cells = 100)
  r <- annotation_result_v1("set_a", "A", engine_id = "manual",
                            cell = ds$cell_data$cell,
                            cell_labels = rep("T cell", 100))

  expect_true(is_annotation_result_v1(r))
  expect_identical(r$schema_version, "annotation_v1")
  expect_identical(r$set_id, "set_a")

  # Manual must NOT have a different (simpler) schema; assert full surface.
  expected <- c("schema_version", "set_id", "name", "description",
                "engine_id", "engine_version", "params",
                "cell", "cell_labels", "cell_scores",
                "alt_labels", "cluster_summary", "ontology_map",
                "reference_source", "marker_registry_version",
                "parent_set_id", "cluster_field_used",
                "n_clusters_at_creation", "is_frozen", "is_demo",
                "warnings", "created_at", "modified_at", "timestamp",
                "duration_ms", "error_message", "edit_history")
  expect_true(all(expected %in% names(r)))
})

test_that("annotation_result_v1 auto-fills cell_scores from labels", {
  cells <- paste0("c", seq_len(20))
  labs  <- c(rep(NA_character_, 5), rep("T cell", 15))
  r <- annotation_result_v1("s", "s", cell = cells, cell_labels = labs)
  expect_identical(sum(r$cell_scores == 1.0), 15L)
  expect_identical(sum(r$cell_scores == 0.0), 5L)
})

test_that("annotation_set_hash is deterministic and content-sensitive", {
  cells <- paste0("c", seq_len(30))
  r  <- annotation_result_v1("s", "s", cell = cells,
                             cell_labels = rep("T", 30))
  r2 <- r; r2$cell_labels[1] <- "NK"
  expect_identical(annotation_set_hash(r), annotation_set_hash(r))
  expect_false(identical(annotation_set_hash(r), annotation_set_hash(r2)))

  # Cell order in the hash is normalised (so reorder -> same hash).
  r3 <- r
  ord <- sample.int(30)
  r3$cell        <- r$cell[ord]
  r3$cell_labels <- r$cell_labels[ord]
  expect_identical(annotation_set_hash(r3), annotation_set_hash(r))
})

test_that("expand_cluster_to_cells handles assigned + unassigned clusters", {
  ds <- mock_dataset(n_cells = 200, seed = 1)
  labs <- expand_cluster_to_cells(ds, "cluster",
                                  c("0" = "T cell", "1" = "B cell"))
  expect_identical(length(labs), 200L)
  expect_true(all(labs[ds$cell_data$cluster == "0"] == "T cell"))
  expect_true(all(labs[ds$cell_data$cluster == "1"] == "B cell"))
  expect_true(all(is.na(labs[!ds$cell_data$cluster %in% c("0", "1")])))
})

# ---- Engines -------------------------------------------------------------

test_that("ANNOTATION_ENGINES registers every shipping engine", {
  ids <- vapply(ANNOTATION_ENGINES(), `[[`, character(1), "id")
  expect_setequal(ids, c("manual", "marker_score",
                         "singler", "azimuth", "celltypist"))

  expect_identical(get_annotation_engine("manual")$id,     "manual")
  expect_identical(get_annotation_engine("singler")$id,    "singler")
  expect_identical(get_annotation_engine("azimuth")$id,    "azimuth")
  expect_identical(get_annotation_engine("celltypist")$id, "celltypist")
  expect_null(get_annotation_engine("not_an_engine"))
})

test_that("manual engine produces per-cell labels via cluster expansion", {
  setup <- mock_state_with_dataset(n_cells = 120)
  state <- setup$state; ds <- setup$dataset
  res <- with_state(state, run_annotation_engine(
    "manual", ds, state,
    params = list(cluster_field = "cluster",
                  labels = list("0" = "CD4 T cell", "1" = "B cell")),
    set_id = "m1", set_name = "Manual 1"))

  expect_true(is_annotation_result_v1(res))
  expect_identical(res$engine_id, "manual")
  expect_true(any(res$cell_labels == "CD4 T cell"))
  expect_true(any(res$cell_labels == "B cell"))
  # Unmapped clusters yield NA labels (score 0)
  expect_true(any(is.na(res$cell_labels)))
  expect_true(any(res$cell_scores == 0.0))

  # cluster_summary present and consistent
  expect_s3_class(res$cluster_summary, "data.frame")
  expect_true(all(c("cluster", "top_label", "top_score", "n_cells")
                  %in% names(res$cluster_summary)))
})

test_that("marker_score engine populates cluster summary + alt_labels", {
  setup <- mock_state_with_dataset(n_cells = 150)
  state <- setup$state; ds <- setup$dataset
  res <- with_state(state, run_annotation_engine(
    "marker_score", ds, state,
    params = list(cluster_field = "cluster", min_score = 0),
    set_id = "ms1", set_name = "ms 1"))

  expect_identical(res$engine_id, "marker_score")
  expect_true(nrow(res$cluster_summary) > 0L)
  expect_true(all(c("top_label", "top_score") %in% names(res$cluster_summary)))
  expect_s3_class(res$alt_labels, "data.frame")
  expect_true(all(c("cluster", "rank", "label", "score") %in% names(res$alt_labels)))

  # Marker-registry version stamped on result
  reg_v <- with_state(state, state$marker_registry$version)
  expect_identical(res$marker_registry_version, reg_v)
})

test_that("run_annotation_engine rejects unknown engines", {
  setup <- mock_state_with_dataset(n_cells = 60)
  expect_error(
    with_state(setup$state,
               run_annotation_engine("not_an_engine", setup$dataset,
                                     setup$state, list(), set_id = "x")),
    "Unknown annotation engine"
  )
})

# ---- Multi-set lifecycle -------------------------------------------------

test_that("add_annotation_set activates the first set automatically", {
  setup <- mock_state_with_dataset(n_cells = 80)
  state <- setup$state; ds <- setup$dataset

  set <- with_state(state, run_annotation_engine(
    "manual", ds, state,
    params = list(cluster_field = "cluster", labels = list("0" = "T")),
    set_id = "first"))
  with_state(state, {
    add_annotation_set(state, set)
    expect_identical(state$active_annotation_id, "first")
  })
})

test_that("duplicate / rename / freeze / delete behave correctly", {
  setup <- mock_state_with_dataset(n_cells = 80)
  state <- setup$state; ds <- setup$dataset

  set <- with_state(state, run_annotation_engine(
    "manual", ds, state,
    params = list(cluster_field = "cluster", labels = list("0" = "T")),
    set_id = "src"))
  with_state(state, add_annotation_set(state, set))

  dup_id <- with_state(state, duplicate_annotation_set(state, "src", "Copy"))
  with_state(state, {
    expect_true(dup_id %in% names(state$annotation_sets))
    expect_identical(state$annotation_sets[[dup_id]]$parent_set_id, "src")
    expect_identical(state$annotation_sets[[dup_id]]$name, "Copy")
  })

  with_state(state, rename_annotation_set(state, dup_id, "Renamed"))
  expect_identical(with_state(state, state$annotation_sets[[dup_id]]$name),
                   "Renamed")

  with_state(state, freeze_annotation_set(state, dup_id, TRUE))
  expect_true(with_state(state, isTRUE(state$annotation_sets[[dup_id]]$is_frozen)))
  expect_error(with_state(state, rename_annotation_set(state, dup_id, "x")),
               "frozen")
  expect_error(with_state(state, remove_annotation_set(state, dup_id)),
               "frozen")

  with_state(state, freeze_annotation_set(state, dup_id, FALSE))
  with_state(state, remove_annotation_set(state, dup_id))
  expect_false(dup_id %in% with_state(state, names(state$annotation_sets)))
})

test_that("set_active_annotation rejects unknown ids", {
  setup <- mock_state_with_dataset(n_cells = 50)
  expect_error(
    with_state(setup$state, set_active_annotation(setup$state, "ghost")),
    "No annotation set"
  )
  expect_null(with_state(setup$state, {
    set_active_annotation(setup$state, NULL)
    get_active_annotation(setup$state)
  }))
})

# ---- Apply to dataset metadata ------------------------------------------

test_that("apply_annotations_to_dataset writes a provenance-named column", {
  setup <- mock_state_with_dataset(n_cells = 100)
  state <- setup$state; ds <- setup$dataset
  set <- with_state(state, run_annotation_engine(
    "marker_score", ds, state,
    params = list(cluster_field = "cluster"),
    set_id = "ms_apply"))

  cols_before <- names(ds$cell_data)
  ds2 <- apply_annotations_to_dataset(ds, set)
  new_cols <- setdiff(names(ds2$cell_data), cols_before)
  expect_length(new_cols, 1L)
  expect_match(new_cols, "^annotation__ms_apply__\\d{4}_\\d{2}_\\d{2}$")
  expect_true(new_cols %in% ds2$metadata_fields)

  # Generic cell_type is left untouched (mock had one already)
  expect_identical(ds$cell_data$cell_type, ds2$cell_data$cell_type)

  # Provenance attributes attached
  expect_identical(attr(ds2$cell_data[[new_cols]], "annotation_set_id"), "ms_apply")
  expect_identical(attr(ds2$cell_data[[new_cols]], "annotation_engine_id"),
                   "marker_score")
  expect_identical(attr(ds2$cell_data[[new_cols]], "schema_version"),
                   "annotation_v1")
})

test_that("apply_annotations_to_dataset refuses to overwrite an existing column", {
  setup <- mock_state_with_dataset(n_cells = 50)
  state <- setup$state; ds <- setup$dataset
  set <- with_state(state, run_annotation_engine(
    "manual", ds, state,
    params = list(cluster_field = "cluster", labels = list("0" = "T")),
    set_id = "noverwrite"))
  ds2 <- apply_annotations_to_dataset(ds, set)
  expect_error(apply_annotations_to_dataset(ds2, set), "already exists")
})

# ---- Phase-2 engines explicitly absent ----------------------------------

test_that("Remaining Phase-2 engines (consensus / gpt / etc.) are not yet registered", {
  ids <- vapply(ANNOTATION_ENGINES(), `[[`, character(1), "id")
  # `singler`, `azimuth`, and `celltypist` graduated out of Phase-2;
  # the remaining placeholders below should still be absent.
  phase2 <- c("consensus", "gpt", "scvi", "scpred", "garnett", "sctype")
  expect_true(length(intersect(ids, phase2)) == 0L)
})

# ---- Engine version metadata (P1) ---------------------------------------
# `engine_version` is the engine *implementation* version, NOT the result
# schema version. They are independent provenance axes. An earlier bug
# stamped `engine$result_schema` as `engine_version`, conflating the two.
# These tests pin the contract so the regression cannot return silently.

test_that("annotation_engine_spec requires a non-empty `version`", {
  expect_error(
    annotation_engine_spec(id = "test_engine", name = "t", category = "test",
                           run_fn = function(...) list()),
    regexp = "version.*non-empty"
  )
  expect_error(
    annotation_engine_spec(id = "test_engine", name = "t", category = "test",
                           run_fn = function(...) list(), version = ""),
    regexp = "version.*non-empty"
  )
})

test_that("every built-in annotation engine declares an explicit version", {
  for (e in ANNOTATION_ENGINES()) {
    expect_true(is.character(e$version), info = e$id)
    expect_true(nzchar(e$version), info = e$id)
    # Sanity: each version mentions its engine id, not the schema label.
    expect_false(identical(e$version, "annotation_v1"), info = e$id)
    expect_false(identical(e$version, e$result_schema), info = e$id)
  }
})

test_that("run_annotation_engine stamps engine_version from spec, not from result_schema", {
  setup <- mock_state_with_dataset(n_cells = 50)
  state <- setup$state; ds <- setup$dataset

  set <- with_state(state, run_annotation_engine(
    "marker_score", ds, state,
    params = list(cluster_field = "cluster"),
    set_id = "version_check"))

  # engine_version is the spec version (provenance), not the schema label
  expect_identical(set$engine_version, "marker_score_v1.0.0")
  expect_identical(set$schema_version, "annotation_v1")
  expect_false(identical(set$engine_version, set$schema_version))

  # manual engine uses its own distinct version
  set2 <- with_state(state, run_annotation_engine(
    "manual", ds, state,
    params = list(cluster_field = "cluster", labels = list("0" = "T")),
    set_id = "version_check_manual"))
  expect_identical(set2$engine_version, "manual_v1.0.0")
})

test_that("annotation engine_version is propagated to dataset metadata attrs", {
  setup <- mock_state_with_dataset(n_cells = 50)
  state <- setup$state; ds <- setup$dataset
  set <- with_state(state, run_annotation_engine(
    "marker_score", ds, state,
    params = list(cluster_field = "cluster"),
    set_id = "version_attr"))
  ds2 <- apply_annotations_to_dataset(ds, set)
  new_col <- grep("^annotation__version_attr__", names(ds2$cell_data),
                  value = TRUE)
  expect_identical(attr(ds2$cell_data[[new_col]], "annotation_engine_version"),
                   "marker_score_v1.0.0")
  expect_identical(attr(ds2$cell_data[[new_col]], "schema_version"),
                   "annotation_v1")
})
