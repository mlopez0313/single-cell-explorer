# ============================================================================
# Prepared demo dataset (PBMC 8k)
# ----------------------------------------------------------------------------
# The sidebar's "Load demo dataset" button used to call `mock_dataset()` --
# a synthetic ~2.5k Gaussian-blob dataset baked into R/dataset.R. That
# remains in place for the tests (it's the no-dep, no-network fallback),
# but the user-facing demo experience is now backed by a *prepared
# PBMC 8k artifact* serialized to disk by `build_pbmc8k_demo()` (see
# R/demo_dataset_build.R + scripts/build_pbmc8k_demo.R).
#
# Runtime contract:
#
#   * The artifact is a plain `.rds` of a list satisfying `dataset_schema()`
#     with an `expression_backend_sparse` for the expression field. No
#     network access is required at app start or button click.
#   * If the artifact is present, the button loads it.
#   * If the artifact is missing, `load_demo_dataset()` raises a clear
#     error that points at the build script. The app shell catches that
#     error and falls back to `mock_dataset()` so first-run developers
#     still get a working demo experience -- with a loud message in the
#     workspace explaining how to upgrade.
#
# Why a prepared artifact instead of building at runtime:
#
#   * Building PBMC 8k from a TENxPBMCData / Seurat source requires
#     normalisation + PCA + UMAP + clustering. Doing that on demand would
#     make the button take 30s-2min and pull in Seurat as a runtime dep,
#     which neither matches the rest of the app nor scales to repeated
#     clicks during development.
#   * A prepared `.rds` is a single deserialise call; the button stays
#     instant; no surprise installs.
# ============================================================================

#' Find the project root by walking up from `start` looking for a
#' directory that contains BOTH `DESCRIPTION` and `app.R`. This pair
#' uniquely identifies the scrnaExplorer source tree -- no other repo
#' nested below the dev machine's home dir will accidentally match.
#'
#' Robust to whatever `getwd()` happens to be:
#'   * `shiny::runApp("/path/to/scrna-explorer")` -> getwd() is the
#'     project root, walk returns immediately.
#'   * `source("/abs/path/to/scrna-explorer/app.R")` -> walk finds
#'     the project root via app.R's directory.
#'   * `Rscript scripts/build_pbmc8k_demo.R` from anywhere inside the
#'     repo -> walks up to the root.
#'
#' Returns `start` unchanged if no project root can be found within 32
#' levels (i.e. the caller really is outside the project tree).
.find_project_root <- function(start = getwd()) {
  d <- tryCatch(normalizePath(start, mustWork = FALSE),
                error = function(e) start)
  for (.i in seq_len(32L)) {
    if (file.exists(file.path(d, "DESCRIPTION")) &&
        file.exists(file.path(d, "app.R"))) {
      return(d)
    }
    parent <- dirname(d)
    if (identical(parent, d)) break   # reached filesystem root
    d <- parent
  }
  tryCatch(normalizePath(start, mustWork = FALSE),
           error = function(e) start)
}

#' Absolute path to the prepared PBMC 8k demo artifact.
#'
#' Resolution order (first hit wins):
#'   1. `SCE_DEMO_DATASET` env var (full path to an `.rds`). Useful for
#'      CI / dev workflows that keep the artifact outside the repo.
#'   2. Installed-package mode: `system.file("extdata", "pbmc8k_demo.rds",
#'      package = "scrnaExplorer")` if the project has been installed as
#'      a package. Returns "" (and falls through) otherwise.
#'   3. Source-tree mode (the common case): `<project_root>/inst/extdata/
#'      pbmc8k_demo.rds`, with `<project_root>` resolved via
#'      `.find_project_root()` -- robust across machines and launch
#'      methods (`runApp`, `source`, `Rscript`).
#'   4. Legacy fallback: `<project_root>/data/pbmc8k_demo.rds`.
#'
#' The function only computes a path; it does not check whether the file
#' exists. Pair it with `demo_dataset_exists()` for that.
#'
#' @param project_root project root directory. When `NULL` (default), the
#'   resolver walks up from `getwd()` looking for the scrnaExplorer
#'   source tree. Tests pass an explicit value to pin the search.
demo_dataset_path <- function(project_root = NULL) {
  env <- Sys.getenv("SCE_DEMO_DATASET", unset = "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))

  if (is.null(project_root)) {
    # Installed-package mode (rare; project is normally run from
    # source). `system.file()` returns "" when the package is not
    # installed.
    installed <- system.file("extdata", "pbmc8k_demo.rds",
                             package = "scrnaExplorer")
    if (nzchar(installed)) return(normalizePath(installed, mustWork = FALSE))
    project_root <- .find_project_root()
  }

  candidates <- c(
    file.path(project_root, "inst", "extdata", "pbmc8k_demo.rds"),
    file.path(project_root, "data", "pbmc8k_demo.rds")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0L) return(normalizePath(hit[1], mustWork = FALSE))
  # No hit -> return the canonical path so error messages point at
  # the conventional location even before anyone has built the file.
  normalizePath(candidates[1], mustWork = FALSE)
}

#' Does a prepared demo artifact exist on disk?
#'
#' Uses the same path resolution as `demo_dataset_path()`. Pass
#' `project_root` to pin the search (mostly useful in tests).
demo_dataset_exists <- function(project_root = NULL) {
  p <- demo_dataset_path(project_root = project_root)
  file.exists(p) && !dir.exists(p)
}

#' Load the prepared PBMC 8k demo dataset from disk.
#'
#' Reads `path` with `readRDS()`, runs `.validate_demo_dataset()` to make
#' sure the artifact still matches the app's expected schema, and returns
#' the dataset list.
#'
#' Fails clearly when:
#'   * the file is missing (developer has not run the build script yet)
#'   * the file deserialises but does not match `dataset_schema()`
#'
#' This function never reaches out to the network and never installs any
#' package.
#'
#' @param path  path to the prepared `.rds`; defaults to `demo_dataset_path()`
load_demo_dataset <- function(path = demo_dataset_path()) {
  if (!file.exists(path)) {
    stop(
      "Prepared PBMC 8k demo dataset not found at:\n",
      "  ", path, "\n\n",
      "Generate it once with:\n",
      "  Rscript scripts/build_pbmc8k_demo.R\n",
      "(see R/demo_dataset_build.R for the available `source` options and ",
      "their package requirements).\n\n",
      "The build only needs to run once; the artifact is reused on every ",
      "subsequent app launch.",
      call. = FALSE)
  }
  ds <- tryCatch(
    readRDS(path),
    error = function(e)
      stop("Failed to readRDS('", path, "'): ", conditionMessage(e),
           call. = FALSE)
  )
  .validate_demo_dataset(ds, path = path)
  ds
}

#' Public alias matching the `*_dataset()` naming convention used by
#' `mock_dataset()`. Useful from scripts / tests that want to be explicit
#' about which demo dataset they're loading.
pbmc8k_demo_dataset <- function(path = demo_dataset_path()) {
  load_demo_dataset(path = path)
}

# ---- Auto-build / ensure ---------------------------------------------------

#' Packages required to auto-build the demo from the default
#' `tenx_pbmc_data` source. Single source of truth so the can/cannot
#' decision and the missing-package message stay in sync.
.DEMO_AUTO_BUILD_PKGS <- c("Seurat", "SeuratObject", "Matrix",
                           "TENxPBMCData", "SingleCellExperiment",
                           "SummarizedExperiment")

#' Can the running R session auto-build the prepared PBMC 8k artifact?
#'
#' The check is package-based only: every dependency required by the
#' default `tenx_pbmc_data` build source must be installed. The first
#' build still needs network access to fetch PBMC 8k counts via
#' ExperimentHub (cached afterwards) -- we cannot detect that ahead of
#' time, so a network-less environment will simply fail at the actual
#' `TENxPBMCData::TENxPBMCData()` call, which is caught by the
#' surrounding `tryCatch` in `app.R`.
can_build_demo_dataset <- function() {
  all(vapply(.DEMO_AUTO_BUILD_PKGS, has_optional, logical(1)))
}

#' One-line summary of which auto-build packages are missing (or "all
#' present"). Used in workspace warnings when an auto-build attempt is
#' skipped.
demo_auto_build_status <- function() {
  missing <- .DEMO_AUTO_BUILD_PKGS[!vapply(.DEMO_AUTO_BUILD_PKGS,
                                           has_optional, logical(1))]
  if (length(missing) == 0L) return("all auto-build packages present")
  sprintf("missing package(s): %s", paste(missing, collapse = ", "))
}

#' Should the runtime button attempt an auto-build on a missing artifact?
#'
#' Auto-build is enabled by default. Set `SCE_AUTO_BUILD_DEMO=0` (or
#' `"false"` / `"no"` / `"off"`) in the environment to disable -- handy
#' for CI / restricted environments where the user prefers an immediate
#' fallback to `mock_dataset()` without paying a `tryCatch` cost.
demo_auto_build_enabled <- function() {
  v <- tolower(trimws(Sys.getenv("SCE_AUTO_BUILD_DEMO", "1")))
  # Treat unset / empty / blank as enabled (the default). Only explicit
  # opt-out values disable.
  if (!nzchar(v)) return(TRUE)
  !(v %in% c("0", "false", "no", "off"))
}

#' Ensure a prepared demo dataset is on disk, building it if necessary.
#'
#' Resolution order:
#'   1. If the artifact already exists, load it.
#'   2. Else, if auto-build is enabled (`demo_auto_build_enabled()`)
#'      *and* the build packages are installed
#'      (`can_build_demo_dataset()`), run `build_pbmc8k_demo()` from the
#'      default `tenx_pbmc_data` source and then load the result.
#'   3. Else, raise the same clear error `load_demo_dataset()` raises
#'      when an artifact is missing.
#'
#' This is the helper the app shell calls; callers that just want the
#' artifact's path (without triggering a build) should use
#' `load_demo_dataset()` directly.
#'
#' @param progress optional `function(fraction, detail = NULL)` callback
#'                 passed through to `build_pbmc8k_demo()` so a Shiny
#'                 `withProgress()` can render incremental status.
#' @param force_build  if TRUE, rebuild even when an artifact already
#'                     exists. Default FALSE.
ensure_demo_dataset <- function(progress = NULL, force_build = FALSE) {
  if (!force_build && demo_dataset_exists()) {
    return(load_demo_dataset())
  }
  if (!demo_auto_build_enabled()) {
    return(load_demo_dataset())   # raises the standard "not found" error
  }
  if (!can_build_demo_dataset()) {
    stop(
      "Cannot auto-build the PBMC 8k demo artifact: ",
      demo_auto_build_status(), ".\n",
      "Either install the missing package(s) (see README's 'Demo dataset' ",
      "section), or run the build script manually:\n",
      "  Rscript scripts/build_pbmc8k_demo.R",
      call. = FALSE)
  }
  build_pbmc8k_demo(progress = progress)
  load_demo_dataset()
}

#' Validate a deserialised dataset against `dataset_schema()`.
#'
#' Internal. Catches the common corruption modes:
#'   * not a list / missing required fields
#'   * `expression` is not an `expression_backend`
#'   * cell-count / gene-count mismatches between fields
#'
#' Errors are wrapped with a hint that the file may be stale / from a
#' newer build of the app and should be rebuilt.
.validate_demo_dataset <- function(ds, path = "<demo artifact>") {
  hint <- paste0(
    "Rebuild with `Rscript scripts/build_pbmc8k_demo.R` if this artifact ",
    "was generated by an older version of the app.")
  rebuild <- function(msg) {
    stop("Prepared demo dataset at '", path, "' is invalid: ", msg, ". ",
         hint, call. = FALSE)
  }
  if (!is.list(ds))                  rebuild("not a list")
  required <- dataset_schema()
  missing  <- setdiff(required, names(ds))
  if (length(missing) > 0L) {
    rebuild(sprintf("missing schema field(s): %s",
                    paste(missing, collapse = ", ")))
  }
  if (!inherits(ds$expression, "expression_backend"))
    rebuild("`expression` is not an `expression_backend`")
  if (!is.data.frame(ds$cell_data))
    rebuild("`cell_data` is not a data.frame")
  if (length(ds$cells) != nrow(ds$cell_data))
    rebuild(sprintf("length(cells) = %d != nrow(cell_data) = %d",
                    length(ds$cells), nrow(ds$cell_data)))
  if (as.integer(ds$n_cells) != length(ds$cells))
    rebuild(sprintf("n_cells = %d != length(cells) = %d",
                    ds$n_cells, length(ds$cells)))
  # Embedding columns must match the advertised reductions.
  for (r in ds$reductions %||% character()) {
    cols <- paste0(r, c("_1", "_2"))
    miss <- setdiff(cols, names(ds$cell_data))
    if (length(miss) > 0L) {
      rebuild(sprintf("reduction '%s' missing column(s): %s", r,
                      paste(miss, collapse = ", ")))
    }
  }
  invisible(TRUE)
}
