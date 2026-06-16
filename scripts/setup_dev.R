#!/usr/bin/env Rscript
# ============================================================================
# scRNA Explorer dev setup
# ----------------------------------------------------------------------------
# Installs the packages a developer needs for a local run.
#
# Usage:
#
#   Rscript scripts/setup_dev.R                # default tier (demo); prompts
#   Rscript scripts/setup_dev.R --core         # core deps only
#   Rscript scripts/setup_dev.R --demo         # core + PBMC 8k build deps
#   Rscript scripts/setup_dev.R --full         # demo + all optional modules
#   Rscript scripts/setup_dev.R --yes          # non-interactive; assume yes
#   Rscript scripts/setup_dev.R --dry-run      # report only; install nothing
#
# Tiers (cumulative):
#
#   core : `shiny`, `bslib`, `htmltools`, `Matrix`, `rlang`. Required to
#          launch the app shell and run the synthetic `mock_dataset()`.
#   demo : `Seurat`, `SeuratObject`, + Bioconductor: `TENxPBMCData`,
#          `SingleCellExperiment`, `SummarizedExperiment`. Required to
#          auto-build the prepared PBMC 8k demo artifact.
#   full : the optional-module packages -- `presto`, `msigdbr`, `anndata`,
#          `reticulate`, plus Bioc: `zellkonverter`, `rhdf5`, `edgeR`,
#          `DESeq2`, `SingleR`, `celldex`, `slingshot`, `monocle3`,
#          `AUCell`, `fgsea`, `dorothea`, `DropletUtils`.
#
# This script is the single recommended setup step for a fresh machine.
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(name) name %in% args

# Default tier is "demo"; explicit flags override. Avoid a multi-line
# top-level `if / else if / else` chain (the parser closes each branch
# at the end of its line in Rscript).
tier <- "demo"
if (has_flag("--core")) tier <- "core"
if (has_flag("--full")) tier <- "full"
auto    <- has_flag("--yes") || !interactive()
dry_run <- has_flag("--dry-run")

# ---- Resolve project root + bootstrap helpers ----------------------------
# Same walk-up logic as R/demo_dataset.R::.find_project_root() so this
# script works no matter where it's invoked from.
.find_proj <- function(start = getwd()) {
  d <- tryCatch(normalizePath(start, mustWork = FALSE),
                error = function(e) start)
  for (.i in seq_len(32L)) {
    if (file.exists(file.path(d, "DESCRIPTION")) &&
        file.exists(file.path(d, "app.R"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("setup_dev.R: cannot locate scrnaExplorer project root from ",
       start, ". Run this script from inside the repo.",
       call. = FALSE)
}
setwd(.find_proj())

# Source only what we actually need (optional_deps + setup). Sourcing the
# whole R/ tree pulls in Shiny/bslib namespace work that's unnecessary
# here and adds noise to the install logs.
source("R/optional_deps.R", local = FALSE)
source("R/setup.R",         local = FALSE)

# ---- Header ---------------------------------------------------------------
cat("\n=================================================================\n")
cat("  scRNA Explorer -- developer setup\n")
cat("=================================================================\n")
cat(sprintf("  tier    : %s\n", tier))
cat(sprintf("  auto    : %s\n", if (auto) "yes" else "no (will prompt)"))
cat(sprintf("  dry-run : %s\n", if (dry_run) "yes" else "no"))
cat("\n")

# ---- Run -----------------------------------------------------------------
# `sce_setup()` raises on install failure (e.g. missing system library
# blocked a package). Capture so we can still print the post-setup
# preflight summary -- the user wants to see exactly which tiers ended
# up incomplete -- and then exit non-zero so CI / pipelines notice.
setup_err <- tryCatch({ sce_setup(tier = tier, auto = auto, dry_run = dry_run); NULL },
                     error = function(e) e)

# ---- Post-setup status ---------------------------------------------------
cat("\n=================================================================\n")
cat(sce_preflight_message())
cat("\n=================================================================\n")

if (!is.null(setup_err)) {
  message(sprintf("\nsce_setup error: %s", conditionMessage(setup_err)))
  quit(status = 1, save = "no")
}
