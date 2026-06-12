# Sources every R/ file into the test environment so the test files can
# call any helper, engine, or module function as if the app were running.
#
# This file is auto-loaded by testthat before any test-*.R runs (any file
# named `helper-*.R` under `tests/testthat/` is sourced by `test_dir()`).
#
# We deliberately mirror the source order used by `app.R`:
#   1. state.R (defines `%||%` and the state factory)
#   2. all other top-level helpers in R/
#   3. modules under R/modules/
#
# Because this is a non-package project, we do not use pkgload::load_all.

# `testthat::test_path()` returns the tests/testthat/ directory when called
# from inside a test or helper. Robust to where the runner was invoked from.
.proj_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
                            mustWork = TRUE)

.r_dir       <- file.path(.proj_root, "R")
.modules_dir <- file.path(.r_dir, "modules")

if (!dir.exists(.r_dir)) {
  stop("helper-app.R: cannot find R/ relative to ", .proj_root)
}

.source_if_new <- function(path) {
  source(path, local = FALSE)
}

# 1. state.R first (defines `%||%` and the reactive state factory)
.source_if_new(file.path(.r_dir, "state.R"))

# 2. Other top-level R files
.top_files    <- list.files(.r_dir, pattern = "\\.R$", full.names = TRUE)
.top_files    <- setdiff(.top_files, file.path(.r_dir, "state.R"))
.module_files <- list.files(.modules_dir, pattern = "\\.R$", full.names = TRUE)

for (.f in setdiff(.top_files, .module_files)) .source_if_new(.f)
# 3. Modules
for (.f in .module_files) .source_if_new(.f)

# A tiny reactive-aware wrapper used by many tests when calling helpers
# that read from `state` (a reactiveValues object).
#
# Implementation note: the expression must be evaluated in the *caller's*
# frame, otherwise local test variables get shadowed by exported package
# symbols (e.g. testthat::setup). We capture `parent.frame()` before
# isolate enters its own evaluation context.
#
# Example:
#   with_state(state, {
#     stopifnot(get_active_annotation(state)$set_id == "x")
#   })
with_state <- function(state, expr) {
  e  <- substitute(expr)
  pf <- parent.frame()
  shiny::isolate(eval(e, envir = pf))
}

# Construct a freshly-loaded mock dataset + state for tests that need one.
# Centralised so tests don't repeat the boilerplate.
mock_state_with_dataset <- function(n_cells = 200, name = "test") {
  ds <- mock_dataset(n_cells = n_cells, name = name)
  state <- new_app_state()
  shiny::isolate(set_active_dataset(state, ds))
  list(state = state, dataset = ds)
}
