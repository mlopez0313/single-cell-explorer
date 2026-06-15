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
    bioc = c("TENxPBMCData", "SingleCellExperiment", "SummarizedExperiment"),
    desc = paste("Required to auto-build the prepared PBMC 8k demo",
                 "artifact via the default `tenx_pbmc_data` source.")
  ),
  full = list(
    cran = c("presto", "msigdbr", "anndata", "reticulate",
             "BiocManager", "withr", "digest"),
    bioc = c("zellkonverter", "rhdf5", "edgeR", "DESeq2",
             "SingleR", "celldex", "slingshot", "monocle3",
             "AUCell", "fgsea", "dorothea", "DropletUtils"),
    desc = paste("Optional but recommended for the full feature surface:",
                 "real DE backends (presto / edgeR / DESeq2), lazy AnnData",
                 "(rhdf5 / zellkonverter), real annotation engines",
                 "(SingleR + celldex), real trajectory (slingshot /",
                 "monocle3), regulons (AUCell + dorothea), and full",
                 "pathway sources (msigdbr / fgsea).")
  )
)

#' Resolve the cumulative package list for a tier.
#'
#' `core` -> core only; `demo` -> core + demo; `full` -> core + demo + full.
#' Returns a list with `cran` / `bioc` character vectors (deduplicated).
.sce_packages_for_tier <- function(tier) {
  tier <- match.arg(tier, c("core", "demo", "full"))
  order <- c("core", "demo", "full")
  include <- order[seq_len(match(tier, order))]
  cran <- unique(unlist(lapply(include,
                               function(t) .SCE_PKG_TIERS[[t]]$cran),
                        use.names = FALSE))
  bioc <- unique(unlist(lapply(include,
                               function(t) .SCE_PKG_TIERS[[t]]$bioc),
                        use.names = FALSE))
  list(cran = cran, bioc = bioc)
}

#' Per-tier installed/missing breakdown.
#'
#' Returns a named list (one element per tier) with:
#'   * `desc`     : human-readable description
#'   * `cran`     : CRAN packages declared for this tier
#'   * `bioc`     : Bioc packages declared for this tier
#'   * `missing`  : packages from `c(cran, bioc)` not currently installed
#'   * `present`  : number installed
#'   * `total`    : total declared
#'   * `complete` : `length(missing) == 0`
sce_check_setup <- function() {
  out <- list()
  for (tier in names(.SCE_PKG_TIERS)) {
    spec    <- .SCE_PKG_TIERS[[tier]]
    pkgs    <- c(spec$cran, spec$bioc)
    present <- vapply(pkgs, has_optional, logical(1))
    out[[tier]] <- list(
      desc     = spec$desc,
      cran     = spec$cran,
      bioc     = spec$bioc,
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

  missing_cran <- spec$cran[!vapply(spec$cran, has_optional, logical(1))]
  missing_bioc <- spec$bioc[!vapply(spec$bioc, has_optional, logical(1))]

  cat(sprintf("\nscRNA Explorer setup -- tier '%s'\n", tier))
  cat(sprintf("  CRAN missing (%d): %s\n", length(missing_cran),
              if (length(missing_cran)) paste(missing_cran, collapse = ", ")
              else "(none)"))
  cat(sprintf("  Bioc missing (%d): %s\n\n", length(missing_bioc),
              if (length(missing_bioc)) paste(missing_bioc, collapse = ", ")
              else "(none)"))

  if (dry_run) {
    cat("(dry-run: nothing will be installed)\n")
    return(invisible(list(cran = missing_cran, bioc = missing_bioc)))
  }

  if (length(missing_cran) == 0L && length(missing_bioc) == 0L) {
    cat(sprintf("Everything in tier '%s' is already installed.\n", tier))
    return(invisible(list(cran = character(), bioc = character())))
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

  cat(sprintf("\nSetup for tier '%s' complete.\n", tier))
  cat("You can now launch the app with `shiny::runApp(\".\")`.\n")
  invisible(list(cran = missing_cran, bioc = missing_bioc))
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
  # Pin both options for the duration of this call so the install is
  # truly non-interactive from the user's perspective.
  old_opts <- options(
    repos = c(CRAN = "https://cloud.r-project.org"),
    install.packages.check.source = "no")
  on.exit(options(old_opts), add = TRUE)

  # Two-phase progress band: CRAN [0, 0.4), Bioc [0.4, 0.95], verify [0.95, 1].
  if (length(missing_cran) > 0L) {
    tick(0.05, sprintf("Installing CRAN: %s",
                       paste(missing_cran, collapse = ", ")))
    utils::install.packages(missing_cran)
  }

  if (length(missing_bioc) > 0L) {
    if (!has_optional("BiocManager")) {
      tick(0.40, "Installing BiocManager (Bioc installer)")
      utils::install.packages("BiocManager")
      if (!has_optional("BiocManager")) {
        stop("Failed to install BiocManager. Try `install.packages('BiocManager')` ",
             "from an R console to see the underlying error.", call. = FALSE)
      }
    }
    tick(0.45, sprintf("Installing Bioconductor: %s",
                       paste(missing_bioc, collapse = ", ")))
    BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
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
