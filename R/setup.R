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
