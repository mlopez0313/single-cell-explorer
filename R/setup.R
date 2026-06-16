# ============================================================================
# Developer-facing dependency setup
# ----------------------------------------------------------------------------
# This file does NOT install anything at app startup. It defines:
#
#   * `.SCE_PKG_TIERS`        : the single source of truth for which
#                               packages belong to which tier.
#   * `sce_check_setup()`     : returns a structured per-tier status
#                               (used by the launch-time preflight and by
#                               the user-facing CLI).
#   * `sce_preflight_message()`: a multi-line summary string suitable for
#                               `message()` / `cat()`.
#   * `sce_setup(tier, auto)` : installs missing packages for the chosen
#                               tier. Interactive by default; pass
#                               `auto = TRUE` to skip the prompt.
#
# The user-facing CLI lives at scripts/setup_dev.R and just shells through
# to `sce_setup()`.
#
# Tier semantics (cumulative -- "demo" implies "core", "full" implies both):
#
#   core   : packages needed to launch the app shell and run the mock
#            dataset. Already in DESCRIPTION's `Imports`; we still check
#            so the preflight has a meaningful "core OK" line.
#   demo   : packages needed to auto-build the prepared PBMC 8k demo
#            artifact via the default `tenx_pbmc_data` source. The app
#            still launches without these; the demo button falls back to
#            `mock_dataset()`.
#   full   : packages needed for the optional full-feature modules --
#            real DE backends (presto / edgeR / DESeq2), lazy AnnData
#            reads (rhdf5), real annotation engines (SingleR + celldex),
#            real trajectory (slingshot / monocle3), real regulons
#            (AUCell), full pathway sources (msigdbr / fgsea), etc.
# ============================================================================

# ---------------------------------------------------------------------------
# Shared install-call hygiene
# ---------------------------------------------------------------------------
# Both `sce_setup()` (CLI) and `sce_install_for_demo()` (in-app) call
# `install.packages()` / `BiocManager::install()`. Conda-bundled R and
# R sessions launched from `Rscript` typically ship with no CRAN mirror
# configured (`getOption("repos")` is the sentinel `c(CRAN = "@CRAN@")`),
# which makes `install.packages()` either prompt -- in an interactive
# session -- or hard-fail with `trying to use CRAN without setting a
# mirror` under `Rscript`. We default a sensible mirror in that case
# while leaving an already-configured mirror untouched (developers may
# have a faster regional/internal mirror).
#
# When BiocManager is installed we *additionally* fold in
# `BiocManager::repositories()` (which yields CRAN + the four
# Bioconductor sub-repos: `software`, `data/annotation`,
# `data/experiment`, `workflows`). This matters for GitHub installs:
# `remotes::install_github("satijalab/azimuth")` recursively shells
# out to `install.packages()` to fetch Signac's Bioc deps
# (`GenomeInfoDb`, `Rsamtools`, `sparseMatrixStats`, ...), and those
# packages aren't on CRAN. Without Bioc repos in scope the resolver
# reports "not available", Signac fails, Azimuth fails to import,
# and the modal surfaces a generic "the package isn't importable"
# error -- the failure mode we just debugged.
#
# Returns the prior `options()` list (suitable for
# `on.exit(options(old), add = TRUE)` in the caller).
#
# `pin_check_source = TRUE` additionally disables the
# "install from source?" prompt on macOS/Windows (no-op on Linux).
# Picks a Bioconductor mirror base URL.
#
# Resolution order (first hit wins):
#   1. `getOption("BioC_mirror")` -- explicit user choice (e.g. set in
#      an .Rprofile or a corporate-mirror config).
#   2. `Sys.getenv("SCE_BIOC_MIRROR")` -- per-deployment override.
#   3. Default: Posit Public Package Manager's Bioconductor CDN.
#
# Why the Posit default (rather than the canonical
# `https://bioconductor.org`):
#
# As of Bioc 3.20+, `bioconductor.org` issues 302 redirects on every
# release-branch PACKAGES URL, pointing at the OSN mirror
# (`mghp.osn.xsede.org/bir190004-bucket01/...`). That OSN host is
# unreachable from a non-trivial slice of networks (campus firewalls,
# datacenter outbound filters, etc.), producing 60s `install.packages`
# timeouts and a cascade of "package 'X' is not available" failures.
# In the wild we observed `remotes::install_github("satijalab/azimuth")`
# fail with `dependencies 'GenomeInfoDb', 'Rsamtools',
# 'sparseMatrixStats' are not available for package 'Signac'` for
# exactly this reason. Posit's mirror is a globally cached CDN
# that exposes the same `<base>/packages/<release>/<kind>/src/contrib`
# layout BiocManager already expects -- swap in the base URL and the
# rest of the BiocManager URL construction logic just works.
.sce_bioc_mirror <- function() {
  opt <- getOption("BioC_mirror")
  if (!is.null(opt) && nzchar(opt)) return(opt)
  env <- Sys.getenv("SCE_BIOC_MIRROR", unset = "")
  if (nzchar(env)) return(env)
  "https://packagemanager.posit.co/bioconductor"
}

.sce_set_install_options <- function(pin_check_source = TRUE) {
  current <- getOption("repos")
  needs_default <- is.null(current) ||
    !"CRAN" %in% names(current) ||
    is.na(current[["CRAN"]]) ||
    !nzchar(current[["CRAN"]]) ||
    identical(unname(current[["CRAN"]]), "@CRAN@")

  # Compute Bioc repos using our preferred mirror, then restore the
  # BioC_mirror option in place -- we don't want this side effect to
  # leak into the caller's session, and we can't reliably hand it
  # back via `prior` because `options(list(BioC_mirror = NULL))`
  # does NOT remove the key (long-standing R quirk -- only the
  # `options(foo = NULL)` calling form does). The Posit URLs we
  # capture below are baked into `repos`, so subsequent
  # install.packages()/BiocManager::install() calls still route
  # through Posit even with BioC_mirror restored.
  has_bioc_mgr <- requireNamespace("BiocManager", quietly = TRUE)
  bioc_repos <- character()
  if (has_bioc_mgr) {
    prior_bioc <- getOption("BioC_mirror")
    options(BioC_mirror = .sce_bioc_mirror())
    bioc_repos <- tryCatch(BiocManager::repositories(),
                           error = function(e) character())
    if (is.null(prior_bioc)) {
      options(BioC_mirror = NULL)
    } else {
      options(BioC_mirror = prior_bioc)
    }
  }

  new_opts <- list()
  if (needs_default) {
    # Prefer BiocManager's repository list (CRAN + four Bioc
    # sub-repos) when available. Falls back to bare CRAN otherwise.
    repos <- if (length(bioc_repos)) bioc_repos
             else c(CRAN = "https://cloud.r-project.org")
    # Guarantee CRAN is present and absolute, even if BiocManager
    # returned a sentinel/empty entry.
    if (!"CRAN" %in% names(repos) ||
        !nzchar(repos[["CRAN"]]) ||
        identical(unname(repos[["CRAN"]]), "@CRAN@")) {
      repos[["CRAN"]] <- "https://cloud.r-project.org"
    }
    new_opts$repos <- repos
  } else if (length(bioc_repos)) {
    # CRAN was already configured by the user (e.g. via .Rprofile);
    # respect their choice but additively merge Bioc sub-repos in if
    # they're missing -- otherwise a user-set CRAN-only profile would
    # still trigger the Signac/Bioc-deps failure above.
    bioc_only <- bioc_repos[setdiff(names(bioc_repos), names(current))]
    if (length(bioc_only) > 0L) {
      new_opts$repos <- c(current, bioc_only)
    }
  }
  if (pin_check_source) {
    new_opts$install.packages.check.source <- "no"
  }
  if (!length(new_opts)) return(list())
  do.call(options, new_opts)
}

#' Tier definitions. Split CRAN vs Bioc so the installer dispatches each
#' to the right repository.
.SCE_PKG_TIERS <- list(
  core = list(
    cran = c("shiny", "bslib", "htmltools", "Matrix", "rlang"),
    bioc = character(),
    desc = paste("Required to launch the app shell and run the synthetic",
                 "`mock_dataset()`.")
  ),
  demo = list(
    cran = c("Seurat", "SeuratObject"),
    bioc = character(),
    desc = paste("Required to auto-build the prepared PBMC 8k demo",
                 "artifact. Downloads the matrix from the 10x Genomics",
                 "CDN directly so no Bioconductor packages are needed",
                 "for the demo build (since Bioconductor's Azure blob",
                 "storage was retired in 2026 the legacy ExperimentHub",
                 "path is broken upstream).")
  ),
  full = list(
    cran = c("presto", "msigdbr", "anndata", "reticulate",
             "BiocManager", "remotes", "withr", "digest"),
    bioc = c("zellkonverter", "rhdf5", "edgeR", "DESeq2",
             "SingleR", "celldex", "slingshot", "monocle3",
             "AUCell", "fgsea", "dorothea", "DropletUtils"),
    # GitHub-only packages, keyed by the *installed* package name so
    # `has_optional(names(github))` works. Values are
    # `remotes::install_github()` specs (`"owner/repo[@ref]"`).
    github = c(Azimuth = "satijalab/azimuth"),
    desc = paste("Optional but recommended for the full feature surface:",
                 "real DE backends (presto / edgeR / DESeq2), lazy AnnData",
                 "(rhdf5 / zellkonverter), real annotation engines",
                 "(SingleR + celldex, Azimuth from GitHub), real",
                 "trajectory (slingshot / monocle3), regulons (AUCell +",
                 "dorothea), and full pathway sources (msigdbr / fgsea).")
  )
)

#' Resolve the cumulative package list for a tier.
#'
#' `core` -> core only; `demo` -> core + demo; `full` -> core + demo + full.
#' Returns a list with `cran` / `bioc` character vectors plus a `github`
#' *named* character vector (names = installed package names, values =
#' `remotes::install_github` specs), all deduplicated.
.sce_packages_for_tier <- function(tier) {
  tier <- match.arg(tier, c("core", "demo", "full"))
  order <- c("core", "demo", "full")
  include <- order[seq_len(match(tier, order))]
  pull_named <- function(slot) {
    # GitHub slot is a named character vector. Concatenate while
    # preserving names; downstream dedup is by *name* (the package
    # name) -- two tiers can both name the same package without
    # producing a duplicate install.
    parts <- lapply(include, function(t) .SCE_PKG_TIERS[[t]][[slot]])
    parts <- parts[lengths(parts) > 0L]
    if (!length(parts)) return(character())
    out <- do.call(c, parts)
    out[!duplicated(names(out) %||% out)]
  }
  cran <- unique(unlist(lapply(include,
                               function(t) .SCE_PKG_TIERS[[t]]$cran),
                        use.names = FALSE))
  bioc <- unique(unlist(lapply(include,
                               function(t) .SCE_PKG_TIERS[[t]]$bioc),
                        use.names = FALSE))
  github <- pull_named("github")
  if (!length(github)) github <- character()
  list(cran = cran, bioc = bioc, github = github)
}

#' Per-tier installed/missing breakdown.
#'
#' Returns a named list (one element per tier) with:
#'   * `desc`     : human-readable description
#'   * `cran`     : CRAN packages declared for this tier
#'   * `bioc`     : Bioc packages declared for this tier
#'   * `github`   : named character vector (names = package names,
#'                  values = `remotes::install_github` specs)
#'   * `missing`  : packages from CRAN+Bioc+GitHub not currently installed
#'                  (names only -- the values are unambiguous)
#'   * `present`  : number installed
#'   * `total`    : total declared
#'   * `complete` : `length(missing) == 0`
sce_check_setup <- function() {
  out <- list()
  for (tier in names(.SCE_PKG_TIERS)) {
    spec    <- .SCE_PKG_TIERS[[tier]]
    gh      <- spec$github %||% character()
    pkgs    <- c(spec$cran, spec$bioc, names(gh))
    present <- vapply(pkgs, has_optional, logical(1))
    out[[tier]] <- list(
      desc     = spec$desc,
      cran     = spec$cran,
      bioc     = spec$bioc,
      github   = gh,
      missing  = pkgs[!present],
      present  = as.integer(sum(present)),
      total    = length(pkgs),
      complete = all(present)
    )
  }
  out
}

#' Format `sce_check_setup()` output as a multi-line message.
#'
#' One line per tier with `[OK] / [missing]` prefix and counts; missing
#' tiers get a follow-up line listing the packages plus the exact
#' setup command to run.
sce_preflight_message <- function(status = sce_check_setup()) {
  lines <- "scRNA Explorer dependency preflight:"
  for (tier in names(status)) {
    s <- status[[tier]]
    mark <- if (s$complete) "OK" else "missing"
    cmd  <- if (s$complete) ""
            else sprintf(
              "  (run `Rscript scripts/setup_dev.R --%s` or `sce_setup(\"%s\")`)",
              tier, tier)
    lines <- c(lines, sprintf("  [%-7s] %-4s  %d/%d packages%s",
                              mark, tier, s$present, s$total, cmd))
    if (!s$complete) {
      lines <- c(lines, sprintf("            missing: %s",
                                paste(s$missing, collapse = ", ")))
    }
  }
  paste(lines, collapse = "\n")
}

#' Install missing packages for the chosen tier.
#'
#' @param tier      one of `"core"`, `"demo"`, `"full"`. Cumulative.
#' @param auto      if `TRUE`, install without prompting. Default `FALSE`;
#'                  in a non-interactive R session this is an error
#'                  (refuses to install silently).
#' @param dry_run   if `TRUE`, print the plan but install nothing.
#' @param ask_bioc  passed to `BiocManager::install(ask = ...)`. Default
#'                  `FALSE` so the Bioc installer never blocks on its
#'                  own prompt.
#' @return invisible list with `cran` / `bioc` character vectors of the
#'   packages this call attempted to install (empty when nothing was
#'   missing).
sce_setup <- function(tier    = c("demo", "core", "full"),
                      auto    = FALSE,
                      dry_run = FALSE,
                      ask_bioc = FALSE) {
  tier <- match.arg(tier)
  spec <- .sce_packages_for_tier(tier)

  missing_cran   <- spec$cran[!vapply(spec$cran, has_optional, logical(1))]
  missing_bioc   <- spec$bioc[!vapply(spec$bioc, has_optional, logical(1))]
  missing_github <- spec$github[!vapply(names(spec$github),
                                        has_optional, logical(1))]

  cat(sprintf("\nscRNA Explorer setup -- tier '%s'\n", tier))
  cat(sprintf("  CRAN missing (%d): %s\n", length(missing_cran),
              if (length(missing_cran)) paste(missing_cran, collapse = ", ")
              else "(none)"))
  cat(sprintf("  Bioc missing (%d): %s\n", length(missing_bioc),
              if (length(missing_bioc)) paste(missing_bioc, collapse = ", ")
              else "(none)"))
  cat(sprintf("  GitHub missing (%d): %s\n\n", length(missing_github),
              if (length(missing_github))
                paste(sprintf("%s (%s)",
                              names(missing_github), unname(missing_github)),
                      collapse = ", ")
              else "(none)"))

  if (dry_run) {
    cat("(dry-run: nothing will be installed)\n")
    return(invisible(list(cran   = missing_cran,
                          bioc   = missing_bioc,
                          github = missing_github)))
  }

  if (length(missing_cran) == 0L && length(missing_bioc) == 0L &&
      length(missing_github) == 0L) {
    cat(sprintf("Everything in tier '%s' is already installed.\n", tier))
    return(invisible(list(cran   = character(),
                          bioc   = character(),
                          github = character())))
  }

  if (!auto) {
    if (!interactive()) {
      stop("sce_setup() in a non-interactive session requires `auto = TRUE`",
           " (or `--yes` on the CLI). Refusing to install silently.",
           call. = FALSE)
    }
    ans <- readline("Install missing packages now? [y/N] ")
    if (!grepl("^[yY]", ans)) {
      cat("Aborted -- no packages installed.\n")
      return(invisible(NULL))
    }
  }

  # Default a CRAN mirror if the session doesn't already have one set.
  # Without this, `Rscript scripts/setup_dev.R` on a fresh conda env
  # fails with `trying to use CRAN without setting a mirror`.
  old_opts <- .sce_set_install_options()
  if (length(old_opts)) on.exit(options(old_opts), add = TRUE)

  if (length(missing_cran) > 0L) {
    cat("Installing CRAN packages...\n")
    utils::install.packages(missing_cran)
  }
  if (length(missing_bioc) > 0L) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      cat("Installing `BiocManager` (required to install Bioconductor packages)...\n")
      utils::install.packages("BiocManager")
    }
    cat("Installing Bioconductor packages...\n")
    BiocManager::install(missing_bioc, update = FALSE, ask = ask_bioc)
  }
  if (length(missing_github) > 0L) {
    # `remotes` is a small pure-R package; install it first if absent.
    if (!requireNamespace("remotes", quietly = TRUE)) {
      cat("Installing `remotes` (required to install GitHub packages)...\n")
      utils::install.packages("remotes")
    }
    cat(sprintf("Installing GitHub packages (%s)...\n",
                paste(sprintf("%s <- %s",
                              names(missing_github), unname(missing_github)),
                      collapse = ", ")))
    # `upgrade = "never"` so we don't surprise the user by rebuilding
    # an existing tidyverse / Seurat install; Azimuth pulls in a deep
    # dep tree (Seurat, BPCells, SeuratData) and "always" can take
    # 30+ min on a fresh machine.
    remotes::install_github(unname(missing_github),
                            upgrade = "never",
                            quiet   = FALSE)
  }

  # Verify that every package we *intended* to install is actually
  # importable. `install.packages()` issues a non-fatal warning when a
  # package fails to build, but the script that wraps us has no way to
  # know that unless we re-check. Without this, the previous behaviour
  # printed "Setup complete" even when 9 of 10 packages failed.
  attempted     <- c(missing_cran, missing_bioc, names(missing_github))
  still_missing <- attempted[
    !vapply(attempted, has_optional, logical(1))]

  if (length(still_missing) > 0L) {
    cat(sprintf(
      "\nSetup for tier '%s' FAILED for %d package(s):\n  %s\n",
      tier, length(still_missing),
      paste(still_missing, collapse = ", ")))
    cat(paste0(
      "\nLikely causes:\n",
      "  * A missing system library. Scroll back to the first\n",
      "    `ERROR: configuration failed for package 'XXX'` block and\n",
      "    look for `* deb:` / `* rpm:` / configure hints there.\n",
      "  * A partial Bioconductor install (rare with conda toolchains).\n",
      "    Re-running the script is safe and will only retry the\n",
      "    missing packages.\n",
      "\nIf the failing package needs a conda system library, add it to\n",
      "`conda/environment.yml` and run `mamba env update -f`,\n",
      "then re-run this script.\n"))
    stop(sprintf("setup tier '%s' did not complete (%d package(s) missing)",
                 tier, length(still_missing)),
         call. = FALSE)
  }

  cat(sprintf("\nSetup for tier '%s' complete.\n", tier))
  cat("You can now launch the app with `shiny::runApp(\".\")`.\n")
  invisible(list(cran   = missing_cran,
                 bioc   = missing_bioc,
                 github = missing_github))
}

# ---------------------------------------------------------------------------
# Setup / build logging
# ---------------------------------------------------------------------------
# Why this exists: `install.packages()` and `BiocManager::install()` write
# build output to the controlling terminal. When the install is driven from
# the Shiny session and the user later tries to figure out why
# `TENxPBMCData` failed to install, the actual root cause (e.g. a `KEGGREST`
# source-compile error pulling in `libxml2-dev`) is long gone. These
# helpers tee that output to a persistent file so the user always has a
# log to refer to.

#' Directory where the app writes setup / build logs.
#'
#' `tools::R_user_dir("single-cell-explorer", "cache")/logs/`. Created
#' on first use; otherwise stable across sessions so earlier logs are
#' preserved. Falls back to `tempdir()` if creating the canonical
#' location fails (read-only home, etc.).
sce_log_dir <- function() {
  d <- tryCatch(
    file.path(tools::R_user_dir("single-cell-explorer", which = "cache"),
              "logs"),
    error = function(e) NULL)
  if (is.null(d)) d <- file.path(tempdir(), "single-cell-explorer-logs")
  if (!dir.exists(d)) {
    ok <- tryCatch({ dir.create(d, recursive = TRUE, showWarnings = FALSE); TRUE },
                   error = function(e) FALSE)
    if (!ok) d <- file.path(tempdir(), "single-cell-explorer-logs")
    if (!dir.exists(d))
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  d
}

#' Allocate a fresh timestamped log file.
#'
#' @param prefix file-basename prefix; sanitized to `[A-Za-z0-9_]+`.
#'   Use a short label that names the workflow (e.g. "demo_setup",
#'   "demo_build").
#' @return absolute path; the file is touched (created empty).
sce_open_log <- function(prefix = "setup") {
  prefix <- gsub("[^A-Za-z0-9_]+", "_", as.character(prefix %||% "setup"))
  prefix <- gsub("_+", "_", prefix)            # collapse runs
  prefix <- gsub("^_+|_+$", "", prefix)        # trim ends
  if (!nzchar(prefix)) prefix <- "setup"
  d <- sce_log_dir()
  path <- file.path(d, sprintf("%s_%s.log", prefix,
                               format(Sys.time(), "%Y%m%d_%H%M%S")))
  file.create(path, showWarnings = FALSE)
  path
}

#' Run `expr` with all output mirrored to `log_path` and the terminal.
#'
#' Stdout is `sink()`-ed with `split = TRUE` so the user keeps the
#' live terminal feed AND the log captures the same content. Messages
#' (R's `message()` channel, which `remotes::install_github()` uses
#' extensively for progress updates) are intercepted via a
#' `withCallingHandlers(message = ...)` handler that copies the text
#' into the log file and then lets R's default handler emit the
#' message to stderr as usual -- effectively the same "tee" behaviour
#' that `sink(split = TRUE)` provides for stdout, since R's
#' `sink(type = "message")` does NOT itself support `split = TRUE`.
#' (Earlier versions sank messages to the log only, which made
#' GitHub installs look frozen because the user saw nothing after the
#' last subprocess message.) Warnings are routed through the same
#' handler so a CRAN install warning is also visible in both places.
#'
#' If `expr` errors, the error is re-thrown after sinks close, with
#' `log_path` attached to the condition object via `e$log_path` so
#' callers can surface it.
#'
#' Returns the value of `expr` invisibly on success.
sce_run_with_log <- function(expr, log_path) {
  if (is.null(log_path) || !nzchar(log_path))
    stop("sce_run_with_log: log_path is required.", call. = FALSE)

  out_con <- file(log_path, open = "a", encoding = "UTF-8")

  ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("\n========== [%s] sce_run_with_log: start ==========\n", ts()),
      file = out_con, append = TRUE)

  sink(out_con, type = "output", split = TRUE)

  closed <- FALSE
  close_sinks <- function() {
    if (closed) return(invisible())
    closed <<- TRUE
    try(sink(type = "output"), silent = TRUE)
    try(close(out_con),        silent = TRUE)
  }
  on.exit(close_sinks(), add = TRUE)

  # Calling handlers that tee `message()` / `warning()` payloads into
  # the log file without muffling them. By NOT calling
  # `invokeRestart("muffleMessage")` / `invokeRestart("muffleWarning")`,
  # we let R's default handler still emit to stderr -- so the user
  # sees live status in the terminal AND the log captures the same
  # text.
  log_condition <- function(cnd) {
    txt <- tryCatch(conditionMessage(cnd),
                    error = function(e) "(unloggable condition)")
    # `message()` payloads already end with "\n"; warnings do not.
    if (inherits(cnd, "message")) {
      cat(txt, file = out_con, sep = "", append = TRUE)
    } else {
      cat(sprintf("Warning: %s\n", txt),
          file = out_con, append = TRUE)
    }
    flush(out_con)
  }

  tryCatch(
    withCallingHandlers({
      result <- force(expr)
      cat(sprintf("\n========== [%s] sce_run_with_log: DONE OK ==========\n",
                  ts()),
          file = out_con, append = TRUE)
      invisible(result)
    },
    message = log_condition,
    warning = log_condition),
    error = function(e) {
      cat(sprintf(paste0(
        "\n========== [%s] sce_run_with_log: ERROR ==========\n",
        "Error: %s\n"),
        ts(), conditionMessage(e)),
        file = out_con, append = TRUE)
      # Close sinks BEFORE rethrowing so the workspace warning isn't
      # also captured into the log.
      close_sinks()
      e$log_path <- log_path
      stop(e)
    })
}

#' Read the last `n` lines of a log file (best-effort).
#'
#' Returns `character()` if the log can't be read.
sce_log_tail <- function(log_path, n = 80L) {
  if (is.null(log_path) || !file.exists(log_path)) return(character())
  lines <- tryCatch(readLines(log_path, warn = FALSE),
                    error = function(e) character())
  utils::tail(lines, n)
}

#' Extract the error-shaped lines from a log file.
#'
#' Matches `Error:`, `ERROR:`, `Warning`, `had non-zero exit status`,
#' Bioc-style `* ERROR`, and bare `* removing '...'` lines (which mark
#' rollback of a failed install). Returns the LAST `max_lines`
#' matches so the most recent cascade is visible.
sce_log_summary <- function(log_path, max_lines = 10L) {
  if (is.null(log_path) || !file.exists(log_path)) return(character())
  lines <- tryCatch(readLines(log_path, warn = FALSE),
                    error = function(e) character())
  if (!length(lines)) return(character())
  pat <- paste0(
    "(^\\* ERROR\\b)",                  "|",
    "(^ERROR:)",                        "|",
    "(^Error:)",                        "|",
    "(^Warning)",                       "|",
    "(had non-zero exit status)",       "|",
    "(\\* removing '.*'$)"
  )
  hits <- grepl(pat, lines, perl = TRUE)
  out  <- lines[hits]
  if (length(out) > max_lines) out <- utils::tail(out, max_lines)
  out
}

#' Detect a stale-session lazy-load failure from an error message.
#'
#' R's lazy-load DB (`R/<pkg>.rdb`) is mmap'd into the running session
#' on first use. If the file on disk is replaced (e.g. by an in-process
#' `install.packages()` or a partially-completed source build), the
#' next access surfaces as one of a small set of messages whose only
#' reliable fix is a fresh R session. We match them here so callers
#' can render a "please restart" hint instead of a generic warning.
#'
#' Intentionally a single regex with the patterns we've observed in the
#' wild on macOS / Linux:
#'   * `lazy-load database '...' is corrupt`
#'   * `bad restore file magic number`
#'   * `unable to load shared object` (when a `.so` was rewritten mid-session)
#'   * `internal error -3 in R_decompress1`
#'
#' Returns a single logical.
.is_lazy_load_corruption <- function(msg) {
  if (length(msg) == 0L || !nzchar(msg)) return(FALSE)
  grepl(paste0(
    "lazy-load database.*is corrupt",   "|",
    "bad restore file magic number",    "|",
    "internal error -3 in R_decompress",  "|",
    "unable to load shared object.*\\.so"),
    msg, ignore.case = TRUE)
}

# Allocate a per-attempt sub-directory inside the main log dir to
# hold the per-package `<pkg>.out` files that `install.packages()`
# preserves when `keep_outputs = <dir>` is set. Sub-directory name is
# timestamped so concurrent attempts don't stomp on each other.
.sce_pkg_outputs_dir <- function() {
  d <- file.path(sce_log_dir(),
                 sprintf("pkglogs_%s",
                         format(Sys.time(), "%Y%m%d_%H%M%S")))
  if (!dir.exists(d))
    tryCatch(dir.create(d, recursive = TRUE, showWarnings = FALSE),
             error = function(e) NULL)
  d
}

# After an install run, fold each `<pkg>.out` file into the parent
# log via a clearly delimited section so the user can read a single
# file. The per-package outputs are kept on disk too (in `pkglog_dir`)
# in case the parent log isn't accessible.
#
# Best-effort: any read / write error is silently ignored so a broken
# pkglog dir can never mask the install's actual outcome.
.sce_drain_pkg_outputs <- function(pkglog_dir) {
  if (is.null(pkglog_dir) || !dir.exists(pkglog_dir)) return(invisible(NULL))
  out_files <- list.files(pkglog_dir, pattern = "\\.out$",
                          full.names = TRUE)
  if (!length(out_files)) return(invisible(NULL))

  # Order so failures land last (most useful at end of log). We can't
  # know success/failure from a `.out` filename, so just sort
  # lexicographically.
  out_files <- sort(out_files)

  for (f in out_files) {
    pkg <- sub("\\.out$", "", basename(f))
    body <- tryCatch(readLines(f, warn = FALSE),
                     error = function(e) sprintf(
                       "(could not read %s: %s)", f, conditionMessage(e)))
    cat(sprintf("\n---------- per-package install log: %s ----------\n",
                pkg),
        paste(body, collapse = "\n"), "\n",
        sep = "")
  }
  cat(sprintf("\n(per-package logs preserved on disk in %s)\n", pkglog_dir))
  invisible(NULL)
}

#' Install just the packages the PBMC 8k demo build needs.
#'
#' Used by the sidebar's "Load demo dataset" button when the user
#' confirms an in-app install in the auto-build modal. Differs from
#' `sce_setup(tier = "demo", auto = TRUE)` in two ways:
#'
#'   * it takes a `progress` callback so a Shiny `withProgress()` can
#'     render incremental status across CRAN + Bioc installs;
#'   * it always assumes consent (the modal click is the consent) and
#'     never prompts.
#'
#' Returns invisibly on success; raises a clear error if any requested
#' package still cannot be loaded after the install attempt (most
#' commonly because of a Bioconductor compile failure or a missing
#' system library).
#'
#' @param progress optional callback `function(fraction, detail = NULL)`.
sce_install_for_demo <- function(progress = NULL) {
  tick <- if (is.null(progress)) function(...) NULL
          else function(fraction, detail = NULL)
                 tryCatch(progress(fraction, detail = detail),
                          error = function(e) NULL)

  spec         <- .sce_packages_for_tier("demo")
  missing_cran <- spec$cran[!vapply(spec$cran, has_optional, logical(1))]
  missing_bioc <- spec$bioc[!vapply(spec$bioc, has_optional, logical(1))]

  if (length(missing_cran) == 0L && length(missing_bioc) == 0L) {
    tick(1.0, "All demo-build packages already installed")
    return(invisible(TRUE))
  }

  # Refuse in-place install when any target package is already loaded
  # in this R session. Reinstalling a loaded package rewrites its
  # `Meta/` and `R/*.rdb` files on disk while R still has stale file
  # handles cached, which manifests as
  #   "lazy-load database '.../Seurat.rdb' is corrupt"
  # on the very next call into that package. The only safe fix is a
  # fresh R session.
  pre_loaded <- intersect(c(missing_cran, missing_bioc), loadedNamespaces())
  if (length(pre_loaded) > 0L) {
    stop(sprintf(paste0(
      "Cannot install package(s) (%s) in this R session because %s ",
      "already loaded. Quit and relaunch the app (Ctrl+C in the ",
      "terminal, then `shiny::runApp(\".\")` again) before clicking ",
      "\"Install + build\". The first install in a fresh session ",
      "completes safely."),
      paste(pre_loaded, collapse = ", "),
      if (length(pre_loaded) == 1L) "it is" else "they are"),
      call. = FALSE)
  }

  # In a Shiny session `interactive()` is TRUE, so `install.packages()`
  # may prompt at the controlling terminal for:
  #   * a CRAN mirror,
  #   * "Do you want to install from sources the package which needs
  #     compilation? (Yes/no/cancel)" (macOS / Windows).
  # `.sce_set_install_options()` defaults a CRAN mirror if none is
  # already set, and pins `install.packages.check.source = "no"` to
  # suppress the macOS/Windows prompt. The CLI path uses the same
  # helper so the two flows can't drift.
  old_opts <- .sce_set_install_options(pin_check_source = TRUE)
  if (length(old_opts)) on.exit(options(old_opts), add = TRUE)

  # Relocate subprocess TMPDIR + parent download dir onto the roomiest
  # filesystem we can find (typically the same volume as the app and
  # the conda R library). This is also what the GitHub-install path
  # does; centralising it here avoids the "/tmp is full" silent
  # failure on machines where the system / is small but the app lives
  # on a large mount.
  scratch_state <- .sce_setup_install_scratch()
  on.exit(.sce_restore_install_scratch(scratch_state), add = TRUE)
  # Demo build pulls CRAN-only Seurat stack -- smaller than Azimuth's
  # transitive Bioc tree but still ~1 GB worth of tarballs + build
  # scratch. Use a 1.5 GB floor.
  .sce_check_tempdir_space(1500,
                           "the in-app demo dataset install",
                           scratch = scratch_state$scratch)

  # Two-phase progress band: CRAN [0, 0.4), Bioc [0.4, 0.95], verify [0.95, 1].
  #
  # `Ncpus = 1L` is important: parallel installs pipe each subprocess's
  # output to per-package files that get deleted at the end of the
  # run, so `sink()`-based logging captures only the "Warning ... had
  # non-zero exit status" summary lines and not the actual compiler
  # errors that the user needs to debug.
  #
  # `keep_outputs = .sce_pkg_outputs_dir()` preserves `<pkg>.out`
  # files (the full per-package `R CMD INSTALL` log) under the same
  # log directory the parent caller is using. After the install
  # returns -- even if it failed -- we copy / concat those files into
  # the parent log so the user has a single artefact to read.
  pkglog_dir <- .sce_pkg_outputs_dir()
  on.exit(.sce_drain_pkg_outputs(pkglog_dir), add = TRUE)

  if (length(missing_cran) > 0L) {
    tick(0.05, sprintf("Installing CRAN: %s",
                       paste(missing_cran, collapse = ", ")))
    utils::install.packages(missing_cran,
                            Ncpus        = 1L,
                            keep_outputs = pkglog_dir,
                            destdir      = scratch_state$scratch)
  }

  if (length(missing_bioc) > 0L) {
    if (!has_optional("BiocManager")) {
      tick(0.40, "Installing BiocManager (Bioc installer)")
      utils::install.packages("BiocManager",
                              Ncpus        = 1L,
                              keep_outputs = pkglog_dir,
                              destdir      = scratch_state$scratch)
      if (!has_optional("BiocManager")) {
        stop("Failed to install BiocManager. Try `install.packages('BiocManager')` ",
             "from an R console to see the underlying error.", call. = FALSE)
      }
    }
    tick(0.45, sprintf("Installing Bioconductor: %s",
                       paste(missing_bioc, collapse = ", ")))
    # `BiocManager::install()` forwards `...` to `install.packages()`.
    BiocManager::install(missing_bioc,
                         update       = FALSE,
                         ask          = FALSE,
                         Ncpus        = 1L,
                         keep_outputs = pkglog_dir,
                         destdir      = scratch_state$scratch)
  }

  tick(0.95, "Verifying installed packages")
  # Two-step verify so we surface meaningful errors:
  #
  #   1. `requireNamespace()` -> package is installed at all.
  #   2. `loadNamespace()`    -> the lazy-load DB on disk is intact
  #                              (catches partial / interrupted Bioc
  #                              source builds that leave a corrupt
  #                              `Seurat.rdb` etc.).
  installed_ok <- vapply(c(missing_cran, missing_bioc),
                         has_optional, logical(1))
  still_missing <- c(missing_cran, missing_bioc)[!installed_ok]
  if (length(still_missing) > 0L) {
    stop(sprintf(paste0(
      "Demo-build dependency install failed for package(s): %s. ",
      "Re-run `Rscript scripts/setup_dev.R --demo` from a terminal to ",
      "see the underlying CRAN/Bioconductor error (typical causes: ",
      "missing system libraries, no network, source-compile failure)."),
      paste(still_missing, collapse = ", ")),
      call. = FALSE)
  }

  loadable <- function(pkg) {
    tryCatch({ suppressMessages(suppressWarnings(loadNamespace(pkg))); TRUE },
             error = function(e) FALSE)
  }
  corrupt <- c(missing_cran, missing_bioc)[
    !vapply(c(missing_cran, missing_bioc), loadable, logical(1))]
  if (length(corrupt) > 0L) {
    stop(sprintf(paste0(
      "Installed but cannot load: %s. The package files on disk are ",
      "present but unreadable from this R session (often a corrupt ",
      "`R/<pkg>.rdb` from a partially completed source build). ",
      "Quit the app, then either reinstall via ",
      "`Rscript scripts/setup_dev.R --demo` from a fresh terminal or ",
      "reinstall the affected package directly. The dataset will ",
      "auto-build on the next launch."),
      paste(corrupt, collapse = ", ")),
      call. = FALSE)
  }

  tick(1.00, "Done")
  invisible(TRUE)
}

#' Install a single GitHub-hosted package in-process.
#'
#' Used by the annotation module's "Install + run" confirmation modal
#' (and any other in-app feature that needs to bring an optional
#' GitHub package online without restarting R). The contract mirrors
#' `sce_install_for_demo()`:
#'
#'   * Refuses to reinstall a package whose namespace is already
#'     loaded in this R session (rewriting the on-disk `R/<pkg>.rdb`
#'     would corrupt the in-memory lazy-load cache the rest of the
#'     session is still pointing at).
#'   * Bootstraps `remotes` if absent (it's pure R, ~1s install).
#'   * Uses `.sce_set_install_options()` so a missing CRAN mirror
#'     doesn't blow up dependency installs.
#'   * `progress` callback receives fractions in `[0, 1]` and detail
#'     strings; safe to pass `NULL`.
#'   * Verifies the package loads cleanly post-install (`loadNamespace`)
#'     before returning, so a broken install surfaces immediately
#'     instead of at the next call site.
#'
#' @param pkg          installed package name (e.g. `"Azimuth"`). Used
#'                     for `has_optional()` / `loadNamespace()` checks.
#' @param github_spec  `owner/repo[@ref]` string accepted by
#'                     `remotes::install_github()`.
#' @param upgrade      passed through to `remotes::install_github`.
#'                     Default `"never"` so we don't surprise the user
#'                     by rebuilding their existing Seurat/tidyverse
#'                     install.
#' @param progress     optional `function(fraction, detail = NULL)`.
#' @return invisible `TRUE` on success; errors on any failure.
# Returns free space (MB) on the filesystem hosting `path`, or NA if
# we can't determine it (e.g. `df` isn't on PATH). Uses POSIX
# `df -Pk` so output columns are stable across Linux/macOS.
.sce_path_free_mb <- function(path) {
  tryCatch({
    out <- suppressWarnings(system2("df", c("-Pk", shQuote(path)),
                                    stdout = TRUE, stderr = FALSE))
    # Expected: header line + one data line; the 4th whitespace-
    # separated field on the data line is Avail (KB).
    if (length(out) < 2L) return(NA_real_)
    fields <- strsplit(trimws(out[length(out)]), "\\s+", perl = TRUE)[[1L]]
    if (length(fields) < 4L) return(NA_real_)
    as.numeric(fields[length(fields) - 2L]) / 1024  # KB -> MB
  }, error = function(e) NA_real_)
}

# ---------------------------------------------------------------------------
# Install scratch dir
# ---------------------------------------------------------------------------
# `install.packages()` / `remotes::install_github()` use R's `tempdir()` for
# both source-tarball downloads and build/unpack scratch. On deployments
# where the OS `/tmp` (or `/`) is small but the app itself lives on a much
# larger filesystem (e.g. an NFS mount holding both the source tree and
# the conda env), it's wasteful and brittle to push multi-GB packages
# like `BSgenome.Hsapiens.UCSC.hg38` (~870 MB tarball / ~3.4 GB unpacked)
# through the small partition.
#
# Instead, ask the runtime where the app lives, pick a scratch root that
# has substantially more free space than `tempdir()`, and:
#   - set `TMPDIR` (so spawned `R CMD INSTALL` / `R CMD build` subprocesses
#     unpack into the bigger filesystem -- this is where the bulk of disk
#     pressure lives), and
#   - pass `destdir = <scratch>` to install.packages so the parent R's
#     downloads land in the bigger filesystem too.
#
# This is opportunistic and self-locating; nothing is hard-coded. If the
# user sets `SCE_SCRATCH_DIR` we honour it verbatim. Otherwise we
# consider a handful of candidates and pick whichever filesystem has the
# most free space (and at least 1.5x what tempdir() currently has, so
# we don't relocate for tiny gains).

# Candidate scratch roots, in priority order. Non-existent / empty /
# non-writable entries are filtered out by `.sce_install_scratch_dir()`.
.sce_install_scratch_candidates <- function() {
  cands <- c(
    Sys.getenv("SCE_SCRATCH_DIR", unset = ""),
    # `shiny::runApp(".")` is the canonical launch path -- getwd() is the
    # app root in that case, and that's exactly the "wherever the app
    # is installed" location the user wants us to land scratch on.
    tryCatch(getwd(), error = function(e) ""),
    # The active R library: on conda envs this typically lives on the
    # same filesystem as the env itself, which is usually the roomy one.
    if (length(.libPaths())) .libPaths()[1L] else "",
    Sys.getenv("HOME", unset = ""))
  cands <- cands[nzchar(cands)]
  unique(tryCatch(normalizePath(cands, winslash = "/", mustWork = FALSE),
                  error = function(e) cands))
}

# Returns a writable scratch dir on the best available filesystem, or
# `tempdir()` if no candidate beats it. Honours `SCE_SCRATCH_DIR` as an
# explicit override (no free-space check applied).
#
# `candidates` and `free_mb` are injectable for testing.
.sce_install_scratch_dir <- function(
    candidates       = .sce_install_scratch_candidates(),
    free_mb          = .sce_path_free_mb,
    headroom_factor  = 1.5,
    subdir           = ".sce-install-scratch") {
  # Explicit override wins outright; we honour the user's choice without
  # second-guessing free-space.
  override <- Sys.getenv("SCE_SCRATCH_DIR", unset = "")
  if (nzchar(override)) {
    dir.create(override, recursive = TRUE, showWarnings = FALSE)
    return(normalizePath(override, winslash = "/", mustWork = FALSE))
  }

  tmp <- tempdir()
  tmp_mb <- free_mb(tmp)
  best_dir <- tmp
  best_mb  <- if (is.na(tmp_mb)) 0 else tmp_mb
  for (cand in candidates) {
    if (identical(cand, tmp)) next
    if (!dir.exists(cand)) next
    if (file.access(cand, mode = 2L) != 0L) next  # writable?
    mb <- free_mb(cand)
    if (is.na(mb)) next
    if (mb > best_mb * headroom_factor) {
      best_dir <- cand
      best_mb  <- mb
    }
  }
  if (identical(best_dir, tmp)) return(tmp)
  out <- file.path(best_dir, subdir)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out)) return(tmp)
  out
}

# Wires the chosen scratch dir into the env vars / subprocesses used by
# install.packages / remotes::install_github:
#   - exports `TMPDIR` so spawned R subprocesses unpack/build there
#     (the parent R's own `tempdir()` is fixed at startup and can't be
#     redirected mid-session -- that's why the launcher in
#     scripts/run_app.sh also exports TMPDIR if not already set);
#   - returns the scratch path so callers can use it as `destdir` for
#     install.packages.
#
# Returns the scratch dir invisibly, or invisible(tempdir()) if no
# relocation was warranted. Restores TMPDIR via the caller's on.exit
# (returned `prior` arg).
.sce_setup_install_scratch <- function(scratch = .sce_install_scratch_dir()) {
  prior <- Sys.getenv("TMPDIR", unset = NA_character_)
  if (!identical(scratch, tempdir())) {
    Sys.setenv(TMPDIR = scratch)
  }
  list(scratch = scratch, prior_tmpdir = prior)
}

.sce_restore_install_scratch <- function(state) {
  if (is.null(state)) return(invisible(NULL))
  if (is.na(state$prior_tmpdir)) {
    Sys.unsetenv("TMPDIR")
  } else {
    Sys.setenv(TMPDIR = state$prior_tmpdir)
  }
  invisible(NULL)
}

# Common preflight for install paths that may pull multi-GB
# Bioconductor data packages (BSgenome.Hsapiens.UCSC.hg38 alone is
# ~870 MB tarball / ~3.4 GB unpacked). `install.packages` /
# `remotes::install_github` use R's `tempdir()` for both downloads
# and build scratch (we redirect both via `.sce_setup_install_scratch`
# when possible), and an out-of-space failure there manifests as a
# silent zero-byte per-package log and a generic "had non-zero exit
# status" warning -- exactly the failure mode we debugged.
#
# Checks `scratch` (where downloads and subprocess builds will land);
# falls back to checking `tempdir()` when scratch is unset. The
# `scratch_mb` / `free_mb` args are injectable for tests.
.sce_check_tempdir_space <- function(min_mb, context,
                                     scratch    = tempdir(),
                                     scratch_mb = .sce_path_free_mb(scratch),
                                     # legacy arg name (still accepted)
                                     free_mb    = NULL) {
  if (!is.null(free_mb)) scratch_mb <- free_mb
  if (is.na(scratch_mb) || scratch_mb >= min_mb) return(invisible(TRUE))
  # If `scratch` is tempdir() itself, we didn't manage to relocate --
  # tell the user that explicitly and point them at the launcher.
  relocated <- !identical(scratch, tempdir())
  if (relocated) {
    where_line <- sprintf("  scratch dir:   %s\n  tempdir:       %s\n",
                          scratch, tempdir())
    fix_line <- sprintf(paste0(
      "The app already auto-selected `%s` for install scratch but it ",
      "doesn't have enough free space. Either free space on that ",
      "filesystem, or override the choice with the env var\n",
      "  SCE_SCRATCH_DIR=/path/with/lots/of/free/space\n",
      "and restart the app."), scratch)
  } else {
    where_line <- sprintf("  tempdir:       %s\n", tempdir())
    fix_line <- paste0(
      "To fix without restarting on this machine, pick a directory ",
      "on a filesystem with more free space and re-launch the app ",
      "with that directory as `TMPDIR` (preferred -- redirects ",
      "everything) or `SCE_SCRATCH_DIR` (redirects subprocess builds ",
      "and downloads, but not the parent R's own scratch). The app ",
      "ships a launcher that does this for you:\n",
      "  ./scripts/run_app.sh\n",
      "which exports TMPDIR to a scratch dir under the app root if ",
      "you don't set one yourself.")
  }
  stop(sprintf(paste0(
    "Not enough free space for %s.\n%s",
    "  free:          %.0f MB\n",
    "  required:      %.0f MB (rough estimate)\n\n%s"),
    context, where_line, scratch_mb, min_mb, fix_line), call. = FALSE)
}

sce_install_github_pkg <- function(pkg, github_spec,
                                   upgrade  = "never",
                                   progress = NULL) {
  stopifnot(is.character(pkg),         length(pkg) == 1L,         nzchar(pkg))
  stopifnot(is.character(github_spec), length(github_spec) == 1L, nzchar(github_spec))

  tick <- if (is.null(progress)) function(...) NULL
          else function(fraction, detail = NULL)
                 tryCatch(progress(fraction, detail = detail),
                          error = function(e) NULL)

  if (has_optional(pkg)) {
    tick(1.0, sprintf("'%s' is already installed", pkg))
    return(invisible(TRUE))
  }

  # Same lazy-load corruption guard as `sce_install_for_demo()`.
  if (pkg %in% loadedNamespaces()) {
    stop(sprintf(paste0(
      "Cannot install '%s' in this R session because its namespace is ",
      "already loaded. Quit and relaunch the app (Ctrl+C in the ",
      "terminal, then `shiny::runApp(\".\")` again) before retrying ",
      "the install."), pkg), call. = FALSE)
  }

  # Relocate the subprocess TMPDIR and the parent's `destdir` onto the
  # roomiest filesystem we can find before checking free space, so the
  # preflight (and the install itself) benefits from the redirect.
  scratch_state <- .sce_setup_install_scratch()
  on.exit(.sce_restore_install_scratch(scratch_state), add = TRUE)

  # Be generous: Azimuth alone pulls hg38 BSgenome (~3.4 GB unpacked)
  # plus EnsDb.Hsapiens.v86 (~1 GB), and `install.packages` doesn't
  # stream-clean intermediates -- the running peak can easily exceed
  # 5 GB. Set the bar at 6 GB so we catch low-disk situations before
  # they corrupt a half-finished install tree.
  .sce_check_tempdir_space(6000,
                           sprintf("the in-app install of '%s'", pkg),
                           scratch = scratch_state$scratch)

  if (!has_optional("remotes")) {
    tick(0.03, "Installing `remotes` (required to install GitHub packages)")
    old_opts <- .sce_set_install_options(pin_check_source = TRUE)
    if (length(old_opts)) on.exit(options(old_opts), add = TRUE)
    utils::install.packages("remotes", destdir = scratch_state$scratch)
    if (!has_optional("remotes")) {
      stop("Failed to install `remotes`. Try ",
           "`install.packages('remotes')` from an R console to see ",
           "the underlying error.", call. = FALSE)
    }
  } else {
    old_opts <- .sce_set_install_options(pin_check_source = TRUE)
    if (length(old_opts)) on.exit(options(old_opts), add = TRUE)
  }

  pkglog_dir <- .sce_pkg_outputs_dir()
  on.exit(.sce_drain_pkg_outputs(pkglog_dir), add = TRUE)

  tick(0.10, sprintf("Installing %s from GitHub (%s)", pkg, github_spec))
  # `remotes::install_github` forwards `...` to the dependency
  # install path, which lets us drop dep build logs into the same
  # `pkglog_dir` so a compile failure deep in the dep tree (e.g.
  # `BPCells` needing Rust) is captured the same way as a Bioc
  # cascade. `destdir` is forwarded to install.packages for the
  # transitive Bioc/CRAN deps -- multi-GB packages like
  # BSgenome.Hsapiens.UCSC.hg38 then download into our scratch dir
  # instead of the small parent tempdir.
  remotes::install_github(github_spec,
                          upgrade      = upgrade,
                          quiet        = FALSE,
                          Ncpus        = 1L,
                          keep_outputs = pkglog_dir,
                          destdir      = scratch_state$scratch)

  tick(0.95, sprintf("Verifying '%s'", pkg))
  if (!has_optional(pkg)) {
    stop(sprintf(paste0(
      "GitHub install of '%s' (from %s) ran but the package isn't ",
      "importable. Most likely a dependency failed to build -- see ",
      "the per-package install logs in the setup log directory for ",
      "details."), pkg, github_spec), call. = FALSE)
  }
  loadable <- tryCatch({
    suppressMessages(suppressWarnings(loadNamespace(pkg)))
    TRUE
  }, error = function(e) {
    if (.is_lazy_load_corruption(conditionMessage(e))) {
      stop(sprintf(paste0(
        "Installed '%s' but loading it produced a lazy-load corruption ",
        "error (a transient symptom of installing into the running R ",
        "session). Quit the app, restart it, and the package will load ",
        "cleanly on the next launch."), pkg), call. = FALSE)
    }
    FALSE
  })
  if (isFALSE(loadable)) {
    stop(sprintf(paste0(
      "Installed '%s' but it does not load. A session restart is the ",
      "most common fix."), pkg), call. = FALSE)
  }

  tick(1.0, sprintf("'%s' installed", pkg))
  invisible(TRUE)
}
