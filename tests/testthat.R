# Test runner entry point.
#
# Invocation (from the project root):
#   Rscript tests/testthat.R
#
# Or, equivalently:
#   R -e 'testthat::test_dir("tests/testthat", reporter = "summary")'
#
# This is a non-package project, so we do NOT call `test_check()`
# (which expects the package to be loaded by pkgload). Instead
# `test_dir()` runs every `test-*.R` under `tests/testthat/`, having
# first sourced any `helper-*.R`. `helper-app.R` loads the app's R/
# files into the global env so tests can call any helper directly.

suppressPackageStartupMessages({
  library(testthat)
  library(shiny)
})

# Resolve project root from the location of this script if possible, so
# that `Rscript tests/testthat.R` works regardless of the working dir.
.locate_self <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("^--file=", args)]
  if (length(file_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", file_arg)))
  }
  this <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(this) && nzchar(this)) return(normalizePath(this))
  NA_character_
}

.self <- .locate_self()
project_root <- if (is.na(.self)) getwd() else normalizePath(dirname(dirname(.self)))
test_dir_path <- file.path(project_root, "tests", "testthat")

if (!dir.exists(test_dir_path)) {
  stop("Cannot locate tests/testthat/ relative to: ", project_root)
}

reporter <- if (interactive()) "progress" else "summary"
testthat::test_dir(test_dir_path, reporter = reporter, stop_on_failure = TRUE)
