# ============================================================================
# Optional-dependency gating
# ----------------------------------------------------------------------------
# Real dataset loaders, real annotation engines, and real DE backends each
# depend on heavy R packages (Seurat, SeuratObject, Matrix, zellkonverter,
# anndata, presto, msigdbr, SingleR, slingshot, monocle3, ...). None of them
# are required to launch the app shell or run the mock dataset, so they are
# all listed under `Suggests:` in DESCRIPTION.
#
# This file provides the two one-liners every loader/engine should use to
# probe and require those packages:
#
#   has_optional("Seurat")
#   require_optional("Seurat", feature = "Seurat .rds loading")
#
# The error message produced by `require_optional()` is intentionally
# detailed: it names the missing package(s), the feature that needed them,
# and the install command (CRAN or Bioconductor) the user should run.
# ============================================================================

#' Is an optional package installed and loadable?
#'
#' Wraps `requireNamespace(..., quietly = TRUE)` to avoid scattering quoted
#' `tryCatch` blocks through every loader. Returns `TRUE` iff every package
#' in `pkgs` is installed.
#'
#' @param pkgs character vector of package names
#' @return logical(1)
has_optional <- function(pkgs) {
  if (length(pkgs) == 0L) return(TRUE)
  all(vapply(pkgs,
             function(p) requireNamespace(p, quietly = TRUE),
             logical(1)))
}

#' Require an optional package, or raise a friendly error.
#'
#' Use at the very top of any loader / engine / backend that depends on a
#' Suggests-only package. The error includes:
#'   * the missing package name(s)
#'   * a one-line description of which feature needs them
#'   * the recommended install command
#'
#' @param pkgs    character vector of package names (all must be present).
#'                Used both for the `requireNamespace()` check AND for the
#'                printed install command.
#' @param feature short human-readable description, e.g. "Seurat .rds loading"
#' @param source  one of "CRAN", "Bioconductor", or "GitHub"; influences the
#'   printed install command. Defaults to "CRAN".
#' @param repo    optional, only consulted when `source = "GitHub"`. Named
#'                character vector mapping `pkgs` to their `owner/repo`
#'                strings, e.g. `c(presto = "immunogenomics/presto")`.
#'                When supplied, the printed `remotes::install_github(...)`
#'                command uses the explicit owner; otherwise a `<owner>`
#'                placeholder is shown.
#' @return invisible TRUE on success, errors otherwise.
require_optional <- function(pkgs, feature,
                             source = c("CRAN", "Bioconductor", "GitHub"),
                             repo   = NULL) {
  source  <- match.arg(source)
  missing <- pkgs[!vapply(pkgs,
                          function(p) requireNamespace(p, quietly = TRUE),
                          logical(1))]
  if (length(missing) == 0L) return(invisible(TRUE))

  install_cmd <- switch(source,
    "CRAN" = sprintf('install.packages(c(%s))',
                     paste(shQuote(missing), collapse = ", ")),
    "Bioconductor" = sprintf(
      'BiocManager::install(c(%s))',
      paste(shQuote(missing), collapse = ", ")),
    "GitHub" = .require_optional_github_cmd(missing, repo)
  )

  stop(sprintf(paste0(
    "%s requires %s package%s '%s' to be installed.\n",
    "  Install with:  %s\n",
    "  Or install the full optional-feature tier in one go:\n",
    "    Rscript scripts/setup_dev.R --full"),
    feature,
    source,
    if (length(missing) > 1L) "s" else "",
    paste(missing, collapse = "', '"),
    install_cmd
  ), call. = FALSE)
}

#' Attach a loaded package's namespace to the global search path.
#'
#' `requireNamespace()` loads a package's namespace but does NOT attach
#' it, which is normally what we want -- it keeps the search path
#' uncluttered. A few external packages (Azimuth, in particular) reach
#' for S4 setters like `SeuratObject::"Key<-"` through the search path
#' instead of via their NAMESPACE imports, so they fail with
#' `could not find function "Key<-"` when their dependency is
#' loaded-but-not-attached. This helper closes that gap on a
#' per-engine basis.
#'
#' Idempotent: returns invisibly if the package is already attached.
#' Raises a clear error if the package isn't installed or the attach
#' itself fails. Use sparingly -- attaching modifies global state for
#' the rest of the session.
#'
#' @param pkg character(1) package name. Must already be installed
#'   (call `require_optional()` first for a friendlier missing-pkg
#'   message).
ensure_attached <- function(pkg) {
  stopifnot(is.character(pkg), length(pkg) == 1L, nzchar(pkg))
  if (paste0("package:", pkg) %in% search()) return(invisible(TRUE))
  tryCatch({
    attachNamespace(pkg)
    invisible(TRUE)
  }, error = function(e) {
    stop(sprintf(
      "ensure_attached('%s'): cannot attach package. %s",
      pkg, conditionMessage(e)), call. = FALSE)
  })
}

.require_optional_github_cmd <- function(missing, repo = NULL) {
  # Per-package owner/repo lookup. For each missing package we prefer the
  # caller-supplied `repo[[pkg]]` mapping; otherwise we fall back to a
  # "<owner>" placeholder. Multiple missing packages get a single
  # `remotes::install_github(c(...))` call.
  paths <- vapply(missing, function(pkg) {
    if (!is.null(repo) && !is.null(repo[[pkg]]) && nzchar(repo[[pkg]])) {
      as.character(repo[[pkg]])
    } else {
      sprintf("<owner>/%s", pkg)
    }
  }, character(1))
  has_placeholder <- any(grepl("^<owner>/", paths))
  cmd <- sprintf("remotes::install_github(c(%s))",
                 paste(shQuote(paths), collapse = ", "))
  if (has_placeholder) {
    cmd <- paste0(cmd, "  # replace <owner> with the upstream repo")
  }
  cmd
}
