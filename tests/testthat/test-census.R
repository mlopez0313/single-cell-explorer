# CELLxGENE Census loader + app wiring

test_that("is_census_path recognises the virtual sentinel", {
  expect_true(is_census_path("/cellxgene-census"))
  expect_true(is_census_path("/cellxgene-census/query"))
  expect_false(is_census_path("/tmp/foo.rds"))
  expect_false(is_census_path(""))
})

test_that("detect_source maps /cellxgene-census to census", {
  expect_identical(detect_source("/cellxgene-census"), "census")
})

test_that("normalize_census_query trims filters and drops blank gene filter", {
  q <- normalize_census_query(
    organism = "Mus musculus",
    obs_value_filter = "  cell_type == 'T cell'  ",
    var_value_filter = "   "
  )
  expect_identical(q$organism, "Mus musculus")
  expect_identical(q$obs_value_filter, "cell_type == 'T cell'")
  expect_null(q$var_value_filter)
})

test_that("census_dataset_name slugifies organism and filter", {
  nm <- census_dataset_name("Homo sapiens", "cell_type == 'B cell'")
  expect_match(nm, "^census_Homo_sapiens_")
  expect_match(nm, "B_cell")
})

test_that("load_census requires obs_value_filter", {
  expect_error(
    load_census(organism = "Homo sapiens", obs_value_filter = "  "),
    "obs_value_filter is required"
  )
})

test_that("load_census errors with install instructions when deps are missing", {
  skip_if(has_optional("cellxgene.census"),
          "cellxgene.census installed; missing-dep path not reachable.")
  err <- tryCatch(
    load_census(organism = "Homo sapiens",
                obs_value_filter = "cell_type == 'B cell'"),
    error = conditionMessage
  )
  expect_match(err, "cellxgene.census")
  expect_match(err, "install_github")
})

test_that("try_load_census validates UI query before hitting the network", {
  r <- try_load_census(list(
    organism = "Homo sapiens",
    obs_value_filter = "",
    var_value_filter = NULL
  ))
  expect_false(r$ok)
  expect_match(r$error, "cell filter")

  r2 <- try_load_census(list(
    organism = "Pan troglodytes",
    obs_value_filter = "cell_type == 'B cell'"
  ))
  expect_false(r2$ok)
  expect_match(r2$error, "Unsupported Census organism")
})

test_that("try_load_census_studies validates organism before hitting the network", {
  r <- try_load_census_studies("Pan troglodytes")
  expect_false(r$ok)
  expect_match(r$error, "Unsupported Census organism")
})

# Synthetic study-metadata normalisation / selection tests (no Census dep)
test_that("normalize_census_studies standardises common columns", {
  raw <- data.frame(
    id = c("ds1", "ds2"),
    title = c("Study 1", "Study 2"),
    collection = c("Col A", "Col B"),
    dataset_total_cell_count = c(10, 20),
    stringsAsFactors = FALSE
  )
  out <- .normalize_census_studies(raw)
  expect_true(all(c("dataset_id", "dataset_title", "collection_name", "cell_count") %in% names(out)))
  expect_identical(out$dataset_id, c("ds1", "ds2"))
})

test_that("census_study_filter builds a dataset_id filter", {
  studies <- data.frame(
    dataset_id = c("abc", "def"),
    dataset_title = c("Study A", "Study B"),
    stringsAsFactors = FALSE
  )
  expect_identical(census_study_filter(studies, 2), "dataset_id == 'def'")
  expect_error(census_study_filter(studies, 3), "out of range")
})

test_that("try_load_dataset routes census source without a filesystem path", {
  r <- try_load_dataset(
    path = "/cellxgene-census",
    source_choice = "census",
    census = list(
      organism = "Homo sapiens",
      obs_value_filter = "",
      var_value_filter = NULL
    )
  )
  expect_false(r$ok)
  expect_match(r$error, "cell filter")
})

test_that("load_dataset(source = 'census') forwards query args to load_census", {
  skip_if(has_optional("cellxgene.census"),
          "cellxgene.census installed; missing-dep path not reachable.")
  expect_error(
    load_dataset("/cellxgene-census",
                 source = "census",
                 organism = "Homo sapiens",
                 obs_value_filter = "cell_type == 'B cell'"),
    "cellxgene.census"
  )
})

test_that("app_load_dataset preserves dataset when census query is invalid", {
  ms <- mock_state_with_dataset(name = "keep_me")
  state <- ms$state
  with_state(state, app_load_dataset(
    state,
    path = "/cellxgene-census",
    source_choice = "census",
    census = list(
      organism = "Homo sapiens",
      obs_value_filter = "",
      var_value_filter = NULL
    )
  ))
  with_state(state, {
    expect_identical(state$active_dataset$name, "keep_me")
    err_msgs <- Filter(function(m) identical(m$level, "error"), state$messages)
    expect_length(err_msgs, 1L)
  })
})

test_that("dataset_path_help_text documents Census browsing", {
  expect_match(dataset_path_help_text("census"), "browse Census studies|CELLxGENE Census")
  expect_match(dataset_path_help_text("census"), "custom SOMA cell filter|selected study")
})

test_that("render_census_study_table renders study labels", {
  studies <- data.frame(
    dataset_id = c("a", "b"),
    dataset_title = c("Study A", "Study B"),
    collection_name = c("Col1", "Col2"),
    cell_count = c(100, 200),
    stringsAsFactors = FALSE
  )
  html <- as.character(render_census_study_table(studies))
  expect_match(html, "Study A")
  expect_match(html, "cells: 100")
})
