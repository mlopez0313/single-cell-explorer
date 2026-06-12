test_that("module_registry returns the documented set of modules", {
  reg <- module_registry()
  expect_true(length(reg) >= 8L)

  ids <- vapply(reg, `[[`, character(1), "id")
  expect_setequal(ids, c(
    "dataset_overview", "scrna_explorer", "marker_investigation",
    "annotation", "differential_expression", "pathway_analysis",
    "imputation", "trajectory", "regulons"
  ))
})

test_that("Phase-1 modules are enabled; Phase-2 modules are disabled", {
  reg <- module_registry()
  enabled <- vapply(Filter(function(m) isTRUE(m$enabled), reg),
                    `[[`, character(1), "id")
  disabled <- vapply(Filter(function(m) !isTRUE(m$enabled), reg),
                     `[[`, character(1), "id")

  # `regulons` (was `regulatory`) graduated out of Phase-2 into a real
  # module; the registry no longer has any disabled entries.
  expect_setequal(enabled, c(
    "dataset_overview", "scrna_explorer", "marker_investigation",
    "annotation", "differential_expression", "pathway_analysis",
    "imputation", "trajectory", "regulons"
  ))
  expect_identical(length(disabled), 0L)
})

test_that("each module_spec carries the required fields", {
  for (m in module_registry()) {
    expect_type(m$id, "character")
    expect_type(m$name, "character")
    expect_type(m$description, "character")
    expect_true(m$category %in% MODULE_CATEGORIES)
    expect_true(is.logical(m$enabled))
    expect_true(is.character(m$required_inputs))
    expect_true(is.function(m$ui_fn) || is.null(m$ui_fn))
    expect_true(is.function(m$server_fn) || is.null(m$server_fn))
  }
})

test_that("get_module looks up by id and returns NULL for unknowns", {
  expect_identical(get_module("trajectory")$id, "trajectory")
  expect_identical(get_module("annotation")$id, "annotation")
  expect_null(get_module("does_not_exist"))
})

test_that("module_inputs_ready respects required_inputs against state", {
  state <- new_app_state()
  mod <- get_module("scrna_explorer")
  # Empty state: no dataset -> not ready
  expect_false(with_state(state, module_inputs_ready(mod, state)))

  with_state(state, set_active_dataset(state, mock_dataset(n_cells = 50)))
  # All inputs cascade-set by set_active_dataset
  expect_true(with_state(state, module_inputs_ready(mod, state)))
})

test_that("modules_by_category returns groups in MODULE_CATEGORIES order", {
  grouped <- modules_by_category()
  expect_true(all(names(grouped) %in% MODULE_CATEGORIES))
  # Order must match the canonical order
  ordered_keys <- intersect(MODULE_CATEGORIES, names(grouped))
  expect_identical(names(grouped), ordered_keys)
})
