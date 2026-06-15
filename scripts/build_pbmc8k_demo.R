#!/usr/bin/env Rscript
# ============================================================================
# Build the prepared PBMC 8k demo dataset
# ----------------------------------------------------------------------------
# One-shot offline pipeline that writes the artifact loaded at runtime
# by the sidebar's "Load demo dataset" button.
#
# Usage (run from the project root):
#
#   Rscript scripts/build_pbmc8k_demo.R
#
# Optional environment variables:
#
#   SCE_DEMO_SOURCE        one of "tenx_pbmc_data" (default), "tenx_dir",
#                          "seurat_object_rds"
#   SCE_DEMO_INPUT         path used by "tenx_dir" / "seurat_object_rds"
#   SCE_DEMO_OUT           output `.rds` path (default: inst/extdata/pbmc8k_demo.rds)
#   SCE_DEMO_SEED          integer seed (default: 8)
#
# Required packages depend on the chosen source -- see the doc on
# `build_pbmc8k_demo()` in R/demo_dataset_build.R for the exact list.
# All deps are gated through `require_optional()` and produce a clear
# install command if missing.
#
# Expect the build to take ~30-90s on a modern laptop for the TENxPBMCData
# path (first run downloads PBMC 8k counts via ExperimentHub; subsequent
# runs reuse the cache).
# ============================================================================

# --- bootstrap: source R/ helpers so we can call build_pbmc8k_demo() -------

proj_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."),
                                     ".."),
                           mustWork = FALSE)
if (!dir.exists(file.path(proj_root, "R"))) {
  # Fallback for `Rscript scripts/build_pbmc8k_demo.R` from project root.
  proj_root <- normalizePath(".", mustWork = TRUE)
}
setwd(proj_root)

# Source order mirrors app.R / tests/testthat/helper-app.R.
source(file.path("R", "state.R"), local = FALSE)
.top <- list.files("R", pattern = "\\.R$", full.names = TRUE)
.top <- setdiff(.top, file.path("R", "state.R"))
.modules <- list.files(file.path("R", "modules"),
                       pattern = "\\.R$", full.names = TRUE)
for (.f in setdiff(.top, .modules)) source(.f, local = FALSE)
for (.f in .modules)                source(.f, local = FALSE)

# --- env-driven config -----------------------------------------------------

src      <- Sys.getenv("SCE_DEMO_SOURCE", "tenx_pbmc_data")
input    <- Sys.getenv("SCE_DEMO_INPUT",   "")
out_path <- Sys.getenv("SCE_DEMO_OUT",     "")
seed_str <- Sys.getenv("SCE_DEMO_SEED",    "8")

if (!nzchar(out_path)) out_path <- demo_dataset_path()
if (!nzchar(input))    input    <- NULL

seed <- suppressWarnings(as.integer(seed_str))
if (is.na(seed)) seed <- 8L

cat("\n== Build PBMC 8k demo dataset ==\n")
cat(sprintf("  source     : %s\n", src))
cat(sprintf("  input      : %s\n", input %||% "(none)"))
cat(sprintf("  output     : %s\n", out_path))
cat(sprintf("  seed       : %d\n", seed))
cat("\n")

# --- run -------------------------------------------------------------------

build_pbmc8k_demo(out_path   = out_path,
                  source     = src,
                  input_path = input,
                  seed       = seed)

cat("\nDone. The sidebar will pick up the artifact on the next app launch.\n")
