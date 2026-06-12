test_that("%||% null-coalescing handles NULL and zero-length", {
  expect_identical("a" %||% "b", "a")
  expect_identical(NULL %||% "b", "b")
  expect_identical(character() %||% "b", "b")
  expect_identical(list() %||% list("b"), list("b"))
})

test_that("new_app_state() returns a reactiveValues with the documented fields", {
  state <- new_app_state()
  expect_s3_class(state, "reactivevalues")

  with_state(state, {
    expect_null(state$active_dataset)
    expect_identical(state$active_module, "dataset_overview")
    expect_null(state$selected_assay)
    expect_null(state$selected_reduction)
    expect_null(state$selected_metadata_field)
    expect_null(state$selected_gene)
    expect_identical(state$selected_cells, character())

    # New (post-refactor) multi-set annotation surface
    expect_identical(state$annotation_sets, list())
    expect_null(state$active_annotation_id)
    expect_null(state$marker_registry)

    expect_identical(state$display_mode_imputation, "raw")
    expect_identical(state$analysis_results, list())
    expect_identical(state$messages, list())
  })
})

test_that("push_message appends to state$messages with level + text + time", {
  state <- new_app_state()
  with_state(state, {
    push_message(state, "hello", "info")
    push_message(state, "danger", "error")
    expect_length(state$messages, 2L)
    expect_identical(state$messages[[1]]$text, "hello")
    expect_identical(state$messages[[1]]$level, "info")
    expect_identical(state$messages[[2]]$level, "error")
    expect_s3_class(state$messages[[1]]$time, "POSIXct")
  })
})

test_that("set_active_dataset cascades defaults and loads the marker registry", {
  state <- new_app_state()
  ds <- mock_dataset(n_cells = 100, name = "ds_a")
  with_state(state, set_active_dataset(state, ds))

  with_state(state, {
    expect_identical(state$active_dataset$name, "ds_a")
    expect_identical(state$selected_assay, ds$default_assay)
    expect_identical(state$selected_reduction, ds$default_reduction)
    expect_identical(state$selected_metadata_field, ds$metadata_fields[1])
    expect_identical(state$selected_gene, ds$genes[1])
    expect_identical(state$selected_cells, character())
    expect_false(is.null(state$marker_registry))
    expect_identical(state$marker_registry$schema_version, "marker_registry_v1")
    expect_identical(state$annotation_sets, list())
    expect_null(state$active_annotation_id)
  })
})

test_that("set_active_dataset resets analysis results and display modes", {
  state <- new_app_state()
  with_state(state, {
    set_active_dataset(state, mock_dataset(n_cells = 50))
    # Pollute state, simulating prior analyses
    state$analysis_results <- list(de = list(status = "completed"))
    state$display_mode_imputation <- "smoothed"
  })
  # Switching datasets clears everything
  with_state(state, set_active_dataset(state, mock_dataset(n_cells = 60, name = "ds_b")))
  with_state(state, {
    expect_identical(state$analysis_results, list())
    expect_identical(state$display_mode_imputation, "raw")
  })
})

test_that("set_active_dataset migrates legacy state$annotations into a Default set", {
  state <- new_app_state()
  with_state(state, {
    state$annotations <- data.frame(
      cluster = c("0", "1"),
      n_cells = c(50L, 30L),
      top_markers = c("CD3D", "MS4A1"),
      suggestion = c("T/NK-like", "B cell-like"),
      annotation = c("Legacy T", "Legacy B"),
      notes = c("", ""),
      stringsAsFactors = FALSE
    )
    state$annotation_cluster_field <- "cluster"
  })
  with_state(state, set_active_dataset(state, mock_dataset(n_cells = 100)))

  with_state(state, {
    expect_length(state$annotation_sets, 1L)
    set <- state$annotation_sets[[1]]
    expect_identical(set$name, "Default")
    expect_true(any(set$cell_labels == "Legacy T"))
    expect_true(any(set$cell_labels == "Legacy B"))
    # Legacy fields cleared after migration
    expect_null(state$annotations)
    expect_null(state$annotation_cluster_field)
  })
})

test_that("get_active_annotation returns NULL with no active set, set on activation", {
  state <- new_app_state()
  with_state(state, set_active_dataset(state, mock_dataset(n_cells = 80)))
  expect_null(with_state(state, get_active_annotation(state)))

  set <- with_state(state, run_annotation_engine(
    "manual", state$active_dataset, state,
    params = list(cluster_field = "cluster", labels = list("0" = "T")),
    set_id = "s1", set_name = "S1"))
  with_state(state, {
    add_annotation_set(state, set)
    set_active_annotation(state, "s1")
    expect_identical(get_active_annotation(state)$set_id, "s1")
    set_active_annotation(state, NULL)
    expect_null(get_active_annotation(state))
  })
})

test_that("make_annotation_stamp returns NA fields without active set", {
  state <- new_app_state()
  with_state(state, set_active_dataset(state, mock_dataset(n_cells = 80)))
  stamp <- with_state(state, make_annotation_stamp(state))
  expect_true(is.na(stamp$annotation_set_id_used))
  expect_true(is.na(stamp$annotation_set_hash_used))
  expect_true(is.na(stamp$annotation_engine_id))
  expect_false(stamp$annotation_set_is_demo)
  expect_s3_class(stamp$stamped_at, "POSIXct")
})

test_that("make_annotation_stamp encodes the active set when present", {
  state <- new_app_state()
  with_state(state, set_active_dataset(state, mock_dataset(n_cells = 80)))
  set <- with_state(state, run_annotation_engine(
    "manual", state$active_dataset, state,
    params = list(cluster_field = "cluster", labels = list("0" = "T")),
    set_id = "s_stamp"))
  with_state(state, {
    add_annotation_set(state, set)
    set_active_annotation(state, "s_stamp")
  })
  stamp <- with_state(state, make_annotation_stamp(state))
  expect_identical(stamp$annotation_set_id_used, "s_stamp")
  expect_identical(stamp$annotation_engine_id, "manual")
  expect_identical(stamp$annotation_set_hash_used,
                   annotation_set_hash(with_state(state, get_active_annotation(state))))
})

test_that("is_result_stale detects id and content drift", {
  state <- new_app_state()
  with_state(state, set_active_dataset(state, mock_dataset(n_cells = 60)))
  set <- with_state(state, run_annotation_engine(
    "manual", state$active_dataset, state,
    params = list(cluster_field = "cluster", labels = list("0" = "T")),
    set_id = "s_stale"))
  with_state(state, {
    add_annotation_set(state, set)
    set_active_annotation(state, "s_stale")
  })

  result <- list(annotation_stamp = with_state(state, make_annotation_stamp(state)))
  expect_false(with_state(state, is_result_stale(result, state)))

  # Mutate labels -> content hash changes -> stale
  with_state(state, {
    sets <- state$annotation_sets
    sets[["s_stale"]]$cell_labels[1] <- "Edited"
    state$annotation_sets <- sets
  })
  expect_true(with_state(state, is_result_stale(result, state)))

  # Deactivate -> stale
  with_state(state, set_active_annotation(state, NULL))
  expect_true(with_state(state, is_result_stale(result, state)))
})

test_that("is_result_stale is FALSE when result has no annotation stamp", {
  state <- new_app_state()
  with_state(state, set_active_dataset(state, mock_dataset(n_cells = 40)))
  expect_false(with_state(state, is_result_stale(list(), state)))
  expect_false(with_state(state, is_result_stale(NULL, state)))
})
