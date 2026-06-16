# ============================================================================
# PBMC 8k demo dataset -- offline builder
# ----------------------------------------------------------------------------
# This file is NOT exercised by the running app. It defines the developer-
# facing helpers that prepare the `.rds` artifact loaded at runtime by
# R/demo_dataset.R. Sources its inputs from one of four places:
#
#   * `tenx_pbmc_10x_cdn`  -- DEFAULT. Downloads the canonical 10x
#                             Genomics PBMC 8k filtered feature-barcode
#                             matrix tarball from `cf.10xgenomics.com`,
#                             caches it under
#                             `tools::R_user_dir("single-cell-explorer",
#                             "cache")/pbmc8k/`, extracts it, and feeds
#                             the matrix dir to `Seurat::Read10X()`.
#                             Only depends on Seurat + Matrix.
#                             Replaced the legacy `tenx_pbmc_data`
#                             source after Bioconductor's
#                             `bioconductorhubs.blob.core.windows.net`
#                             Azure storage account was retired in 2026.
#   * `tenx_pbmc_data`     -- LEGACY / probably broken. Uses
#                             `TENxPBMCData::TENxPBMCData("pbmc8k")` via
#                             ExperimentHub; many resources still point
#                             at the dead Azure URL above. Kept so users
#                             can opt back in if Bioc fixes the entry.
#   * `seurat_object_rds`  -- an existing Seurat `.rds` on disk that has
#                             been preprocessed already (assays + reductions
#                             + clusters). Cheapest path; no Seurat compute.
#   * `tenx_dir`           -- a Cellranger filtered_feature_bc_matrix dir
#                             on local disk. Same code path as the CDN
#                             source after extraction.
#
# Whichever source is used, the output is a `dataset_schema()`-compliant
# list saved as an `.rds`. The runtime loader does NOT know which source
# produced the artifact.
#
# Required (build-time only) packages:
#
#   tenx_pbmc_10x_cdn    : Seurat, SeuratObject, Matrix  (no Bioconductor)
#   tenx_pbmc_data       : Seurat, SeuratObject, TENxPBMCData, SingleCellExperiment, Matrix
#   seurat_object_rds    : SeuratObject (already enforced by .seurat_to_dataset)
#   tenx_dir             : Seurat, SeuratObject, Matrix
#
# The user-facing CLI lives at scripts/build_pbmc8k_demo.R.
# ============================================================================

# 10x Genomics canonical CDN URL for the PBMC 8k filtered feature-barcode
# matrix (Chromium v2 chemistry, GRCh38). Served via Cloudflare; ~37 MB
# tar.gz; ETag-stable since 2017. Don't rely on this never moving --
# every download path checks the HTTP status before extracting.
.PBMC8K_10X_CDN_URL <- "https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_filtered_gene_bc_matrices.tar.gz"

#' Build and save the prepared PBMC 8k demo artifact.
#'
#' Drives the build pipeline end-to-end:
#'   1. Pulls counts from the chosen source.
#'   2. Runs the standard Seurat pipeline (normalisation, variable features,
#'      ScaleData, PCA, UMAP, tSNE, neighbours, Louvain clusters) when the
#'      source doesn't already carry reductions + clusters.
#'   3. Converts to the flat dataset schema via `.seurat_to_dataset()`.
#'   4. Augments the schema with demo-friendly fields: `cluster`,
#'      `cell_type`, `condition`, `pseudotime_demo`, `sample`.
#'   5. Reorders `metadata_fields` so the app's default
#'      (`dataset$metadata_fields[1]`) lands on "sample" -- same as the
#'      mock dataset.
#'   6. `saveRDS(out_path)`.
#'
#' @param out_path  destination `.rds` path. Defaults to `demo_dataset_path()`
#'                  (canonical `inst/extdata/pbmc8k_demo.rds`).
#' @param source    one of `"tenx_pbmc_data"` (default; Bioc),
#'                  `"seurat_object_rds"`, `"tenx_dir"`.
#' @param input_path required when `source != "tenx_pbmc_data"`; path to
#'                   the `.rds` Seurat object or the Cellranger dir.
#' @param dataset_name   character(1) display name baked into `ds$name`.
#' @param n_variable_features integer(1) for `FindVariableFeatures`.
#' @param n_pcs               integer(1) PCs to keep for downstream
#'                            neighbours / UMAP / tSNE.
#' @param cluster_resolution  numeric(1) for `FindClusters(resolution=)`.
#' @param seed       integer(1) controlling all the augmentation
#'                   randomness (condition, sample, pseudotime jitter).
#' @param progress   optional callback `function(fraction, detail = NULL)`.
#'                   Called at named checkpoints throughout the build so a
#'                   Shiny `withProgress()` (or any other listener) can
#'                   render incremental status. `fraction` is in `[0, 1]`.
#'                   Errors raised inside the callback are silently
#'                   ignored so progress reporting can never break a
#'                   build.
#' @return invisibly returns `out_path`.
build_pbmc8k_demo <- function(out_path        = demo_dataset_path(),
                              source          = c("tenx_pbmc_10x_cdn",
                                                  "tenx_pbmc_data",
                                                  "seurat_object_rds",
                                                  "tenx_dir"),
                              input_path      = NULL,
                              dataset_name    = "PBMC 8k (demo)",
                              n_variable_features = 2000L,
                              n_pcs               = 30L,
                              cluster_resolution  = 0.6,
                              seed                = 8L,
                              progress            = NULL) {
  source <- match.arg(source)
  tick <- .build_progress_handle(progress)

  tick(0.02, sprintf("Loading source: %s", source))
  obj <- switch(
    source,
    "tenx_pbmc_10x_cdn"  = .build_load_tenx_pbmc_10x_cdn(progress = tick),
    "tenx_pbmc_data"     = .build_load_tenx_pbmc_data(progress = tick),
    "tenx_dir"           = .build_load_tenx_dir(input_path),
    "seurat_object_rds"  = .build_load_seurat_rds(input_path)
  )

  obj <- .build_run_seurat_pipeline(
    obj,
    n_variable_features = n_variable_features,
    n_pcs               = n_pcs,
    cluster_resolution  = cluster_resolution,
    seed                = seed,
    progress            = tick)

  tick(0.92, "Converting to dataset schema")
  ds <- .seurat_to_dataset(obj, name = dataset_name)
  ds <- .augment_demo_dataset(ds, seed = seed)

  # Schema sanity check before writing. Mirrors the runtime validator
  # so a broken build fails here, not at the user's first click.
  .validate_demo_dataset(ds, path = "<in-memory build>")

  # Ensure target directory exists.
  out_dir <- dirname(out_path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  tick(0.98, "Saving artifact")
  saveRDS(ds, out_path)
  tick(1.00, "Done")
  message(sprintf(
    "Wrote prepared PBMC 8k demo dataset (%d cells x %d genes) to:\n  %s",
    ds$n_cells, ds$n_genes, out_path))
  invisible(out_path)
}

# Wrap a user-supplied progress callback so the rest of the build code
# can call `tick(frac, msg)` unconditionally. NULL becomes a no-op.
# Errors inside the callback are swallowed so a broken UI can never
# break a build.
.build_progress_handle <- function(progress) {
  if (is.null(progress)) return(function(fraction, detail = NULL) NULL)
  function(fraction, detail = NULL) {
    tryCatch(progress(fraction, detail = detail), error = function(e) NULL)
  }
}

# ---- Bioconductor cache bootstrap ----------------------------------------
#
# Pre-create the cache directories used by ExperimentHub /
# AnnotationHub / BiocFileCache so they don't block on
#
#   "<path> does not exist, create directory? (yes/no):"
#
# during the first TENxPBMCData() call from a Shiny session. Best
# effort: if any individual path can't be created (read-only FS, weird
# permissions, etc.) we swallow the error -- the downstream Hub call
# will then surface a clearer message than we could fabricate here.
#
# `getExperimentHubOption("CACHE")` is preferred over
# `tools::R_user_dir(...)` because it respects user-set
# `EXPERIMENT_HUB_CACHE` / `ANNOTATION_HUB_CACHE` env vars.
.ensure_bioc_cache_dirs <- function() {
  mkdir <- function(p) {
    if (length(p) != 1L || is.na(p) || !nzchar(p)) return(invisible(NULL))
    if (!dir.exists(p)) {
      tryCatch(dir.create(p, recursive = TRUE, showWarnings = FALSE),
               error = function(e) NULL)
    }
    invisible(NULL)
  }

  eh_cache <- tryCatch(
    if (requireNamespace("ExperimentHub", quietly = TRUE))
      ExperimentHub::getExperimentHubOption("CACHE")
    else tools::R_user_dir("ExperimentHub", which = "cache"),
    error = function(e) tools::R_user_dir("ExperimentHub", which = "cache"))
  mkdir(eh_cache)

  ah_cache <- tryCatch(
    if (requireNamespace("AnnotationHub", quietly = TRUE))
      AnnotationHub::getAnnotationHubOption("CACHE")
    else tools::R_user_dir("AnnotationHub", which = "cache"),
    error = function(e) tools::R_user_dir("AnnotationHub", which = "cache"))
  mkdir(ah_cache)

  bfc_cache <- tryCatch(tools::R_user_dir("BiocFileCache", which = "cache"),
                        error = function(e) NULL)
  mkdir(bfc_cache)

  invisible(NULL)
}

# ---- Source loaders -------------------------------------------------------

# ---- 10x CDN PBMC 8k source ----------------------------------------------
#
# Downloads the PBMC 8k filtered feature-barcode matrix tarball from
# `cf.10xgenomics.com`, caches the tarball + extracted matrix dir under
# `tools::R_user_dir("single-cell-explorer", "cache")/pbmc8k/` so
# subsequent builds reuse them, and returns a freshly-loaded Seurat
# object via the same code path as the on-disk `tenx_dir` source.
#
# Caches both:
#   * `<cache>/pbmc8k_filtered_gene_bc_matrices.tar.gz`  -- ~37 MB
#   * `<cache>/extracted/`                                -- ~30 MB
#
# Cache hits skip the network entirely; corrupt downloads (wrong size
# or wrong tar header) are removed and retried once.
.build_load_tenx_pbmc_10x_cdn <- function(progress = function(...) NULL) {
  require_optional(c("Seurat", "SeuratObject", "Matrix"),
                   feature = "PBMC 8k demo build (10x CDN source)")
  cache_dir <- .pbmc8k_cdn_cache_dir()
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  tarball <- file.path(cache_dir, "pbmc8k_filtered_gene_bc_matrices.tar.gz")
  matrix_dir <- .pbmc8k_ensure_extracted(tarball,
                                         extract_root = file.path(cache_dir, "extracted"),
                                         progress     = progress)
  progress(0.12, sprintf("Reading 10x matrix from %s",
                         normalizePath(matrix_dir, mustWork = FALSE)))
  counts <- Seurat::Read10X(data.dir = matrix_dir)
  Seurat::CreateSeuratObject(counts = counts, project = "pbmc8k_demo",
                             min.cells = 3, min.features = 200)
}

# Resolve the on-disk cache directory for the 10x CDN source. Honours
# `SCE_DEMO_CACHE_DIR` for tests / power users; otherwise falls back to
# `tools::R_user_dir(..., "cache")/pbmc8k`.
.pbmc8k_cdn_cache_dir <- function() {
  override <- Sys.getenv("SCE_DEMO_CACHE_DIR", unset = "")
  if (nzchar(override)) return(file.path(override, "pbmc8k"))
  base <- tryCatch(tools::R_user_dir("single-cell-explorer", which = "cache"),
                   error = function(e) tempdir())
  file.path(base, "pbmc8k")
}

# Download + extract idempotently. Returns the path to the directory
# containing `matrix.mtx{,.gz}` -- which is what Seurat::Read10X wants.
.pbmc8k_ensure_extracted <- function(tarball, extract_root,
                                     progress = function(...) NULL) {
  if (!.pbmc8k_tarball_ok(tarball)) {
    if (file.exists(tarball)) {
      message("Removing stale / corrupt PBMC 8k tarball at ", tarball)
      unlink(tarball)
    }
    progress(0.06, sprintf("Downloading PBMC 8k from 10x CDN (~37 MB)"))
    .pbmc8k_download(.PBMC8K_10X_CDN_URL, tarball)
    if (!.pbmc8k_tarball_ok(tarball)) {
      stop("PBMC 8k tarball at ", tarball,
           " is still invalid after download (size mismatch or unreadable). ",
           "Delete the file and retry, or pre-download manually from ",
           .PBMC8K_10X_CDN_URL, call. = FALSE)
    }
  } else {
    progress(0.06, "Reusing cached PBMC 8k tarball")
  }

  if (!dir.exists(extract_root)) {
    dir.create(extract_root, recursive = TRUE, showWarnings = FALSE)
  }
  matrix_dir <- .pbmc8k_find_matrix_dir(extract_root)
  if (is.null(matrix_dir)) {
    progress(0.10, "Extracting PBMC 8k matrix")
    rc <- utils::untar(tarball, exdir = extract_root)
    if (!identical(as.integer(rc), 0L)) {
      stop("Failed to untar ", tarball, " (utils::untar returned ", rc, ")",
           call. = FALSE)
    }
    matrix_dir <- .pbmc8k_find_matrix_dir(extract_root)
    if (is.null(matrix_dir)) {
      stop("Extracted PBMC 8k tarball but couldn't locate the matrix ",
           "directory under ", extract_root,
           ". Inspect the directory layout manually.", call. = FALSE)
    }
  }
  matrix_dir
}

# `download.file` with a few practical defaults. Honours
# `getOption("timeout")` if the user set it, else bumps to 300s for
# the 37 MB transfer (R's default 60s is too short on slow links).
.pbmc8k_download <- function(url, dest) {
  prev_timeout <- getOption("timeout", default = 60)
  options(timeout = max(prev_timeout, 300))
  on.exit(options(timeout = prev_timeout), add = TRUE)
  rc <- utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
  if (!identical(as.integer(rc), 0L)) {
    stop("download.file(", url, ") returned exit code ", rc, call. = FALSE)
  }
}

# A tarball is acceptable iff it exists and has the expected size band
# (well above any half-downloaded blob, well below a runaway). The 10x
# PBMC 8k tarball has been ~37.5 MB since 2017; allow a generous band.
.pbmc8k_tarball_ok <- function(path) {
  if (!file.exists(path)) return(FALSE)
  size <- file.info(path)$size
  isTRUE(size > 25 * 1024^2 && size < 60 * 1024^2)
}

# Walk an extracted directory looking for a dir containing the standard
# 10x matrix triplet (`matrix.mtx{,.gz}`, `barcodes.tsv{,.gz}`,
# `genes.tsv{,.gz}` or `features.tsv{,.gz}`). Returns the first match
# (PBMC 8k tarball has only one) or NULL if none is found.
.pbmc8k_find_matrix_dir <- function(root) {
  if (!dir.exists(root)) return(NULL)
  candidates <- list.dirs(root, recursive = TRUE)
  for (d in candidates) {
    files <- list.files(d)
    has_matrix <- any(grepl("^matrix\\.mtx(\\.gz)?$", files))
    has_bc     <- any(grepl("^barcodes\\.tsv(\\.gz)?$", files))
    has_genes  <- any(grepl("^(genes|features)\\.tsv(\\.gz)?$", files))
    if (has_matrix && has_bc && has_genes) return(d)
  }
  NULL
}

.build_load_tenx_pbmc_data <- function(progress = function(...) NULL) {
  require_optional(c("TENxPBMCData", "SingleCellExperiment",
                     "SummarizedExperiment", "Matrix"),
                   feature = "PBMC 8k demo build (TENxPBMCData source)",
                   source  = "Bioconductor")
  require_optional(c("Seurat", "SeuratObject"),
                   feature = "PBMC 8k demo build (Seurat preprocessing)")
  # ExperimentHub / BiocFileCache will prompt the console
  # ("...does not exist, create directory? (yes/no):") the very first
  # time their cache dirs are touched. Inside Shiny that readline goes
  # to the controlling terminal, not the browser, so the app appears
  # hung. Pre-create the dirs to short-circuit the prompt entirely.
  .ensure_bioc_cache_dirs()
  progress(0.05, "Fetching PBMC 8k counts from ExperimentHub (cached after first run)")
  sce <- TENxPBMCData::TENxPBMCData("pbmc8k")
  # TENxPBMCData stores counts as a DelayedArray. Materialise once into
  # a dgCMatrix -- that's what Seurat::CreateSeuratObject expects, and
  # PBMC 8k is small enough (~8k x ~32k) to fit in RAM.
  counts <- methods::as(SummarizedExperiment::assay(sce, 1L), "CsparseMatrix")
  # Gene symbols: TENxPBMCData uses Symbol/Symbol_TENx as the var
  # column. Use whatever column looks like a symbol; fall back to
  # rownames (Ensembl IDs).
  rd <- SingleCellExperiment::rowData(sce)
  symbol_col <- intersect(c("Symbol", "Symbol_TENx", "gene_symbols"),
                          names(rd))
  gene_names <- if (length(symbol_col) > 0L)
    as.character(rd[[symbol_col[1]]]) else rownames(sce)
  gene_names <- make.unique(gene_names)
  rownames(counts) <- gene_names
  colnames(counts) <- colnames(sce)
  Seurat::CreateSeuratObject(counts = counts, project = "pbmc8k_demo",
                             min.cells = 3, min.features = 200)
}

.build_load_tenx_dir <- function(input_path) {
  if (is.null(input_path) || !dir.exists(input_path)) {
    stop("build_pbmc8k_demo(source = 'tenx_dir') requires `input_path` ",
         "to point at a Cellranger feature-barcode matrix directory.",
         call. = FALSE)
  }
  require_optional(c("Seurat", "SeuratObject", "Matrix"),
                   feature = "PBMC 8k demo build (10x dir source)")
  counts <- Seurat::Read10X(data.dir = input_path)
  Seurat::CreateSeuratObject(counts = counts, project = "pbmc8k_demo",
                             min.cells = 3, min.features = 200)
}

.build_load_seurat_rds <- function(input_path) {
  if (is.null(input_path) || !file.exists(input_path)) {
    stop("build_pbmc8k_demo(source = 'seurat_object_rds') requires ",
         "`input_path` to point at an existing .rds Seurat object.",
         call. = FALSE)
  }
  require_optional("SeuratObject",
                   feature = "PBMC 8k demo build (.rds source)")
  obj <- readRDS(input_path)
  if (!inherits(obj, "Seurat"))
    stop("File at '", input_path, "' is not a Seurat object (class: ",
         paste(class(obj), collapse = "/"), ")", call. = FALSE)
  obj
}

# ---- Seurat pipeline ------------------------------------------------------

# Returns a Seurat object that has the assays + reductions + clusters
# the dataset schema needs. If the input already has UMAP/PCA + a
# clustering column, we skip the recompute and reuse what's there.
.build_run_seurat_pipeline <- function(obj,
                                       n_variable_features,
                                       n_pcs,
                                       cluster_resolution,
                                       seed,
                                       progress = function(...) NULL) {
  require_optional(c("Seurat", "SeuratObject"),
                   feature = "PBMC 8k demo build (Seurat preprocessing)")
  set.seed(seed)

  reds_avail <- SeuratObject::Reductions(obj)
  has_pca  <- any(toupper(reds_avail) == "PCA")
  has_umap <- any(toupper(reds_avail) == "UMAP")
  has_clusters <- any(grepl("cluster", names(obj@meta.data),
                            ignore.case = TRUE))

  if (!has_pca || !has_umap || !has_clusters) {
    message("Running Seurat preprocessing (normalise + PCA + UMAP + clusters)...")
    progress(0.18, "Normalising"           ); obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    progress(0.28, "Variable features"     ); obj <- Seurat::FindVariableFeatures(obj, nfeatures = n_variable_features, verbose = FALSE)
    progress(0.38, "Scaling"               ); obj <- Seurat::ScaleData(obj, verbose = FALSE)
    progress(0.52, "PCA"                   ); obj <- Seurat::RunPCA(obj, npcs = n_pcs, verbose = FALSE)
    progress(0.62, "Building neighbour graph"); obj <- Seurat::FindNeighbors(obj, dims = seq_len(n_pcs), verbose = FALSE)
    progress(0.70, "Finding clusters"      ); obj <- Seurat::FindClusters(obj, resolution = cluster_resolution, verbose = FALSE)
    progress(0.80, "UMAP"                  ); obj <- Seurat::RunUMAP(obj, dims = seq_len(n_pcs), verbose = FALSE)
    progress(0.88, "tSNE (optional)"       )
    # tSNE is optional; only useful for the explorer "switch reduction"
    # affordance. Skip if it fails (e.g. Rtsne not installed).
    obj <- tryCatch(
      Seurat::RunTSNE(obj, dims = seq_len(min(n_pcs, 15L)), verbose = FALSE),
      error = function(e) {
        message("  (tSNE skipped: ", conditionMessage(e), ")")
        obj
      })
  }
  obj
}

# ---- Schema augmentation -------------------------------------------------

#' Add the demo-friendly fields the rest of the app implicitly expects.
#'
#' Real Seurat-derived datasets typically don't carry every column the
#' mock dataset advertises. This step:
#'   * promotes `seurat_clusters` (or the first cluster-like column) to
#'     a column literally called `cluster`,
#'   * adds a `cell_type` placeholder column derived from cluster id so
#'     the Annotation module has a non-empty baseline,
#'   * fabricates a `condition` column ("ctrl"/"treat") -- not real biology;
#'     it gives DE / pseudobulk demo a balanced grouping out of the box,
#'   * fabricates a `sample` column ("S1"/"S2"/"S3") -- same rationale,
#'     used by pseudobulk DE's `sample_by` argument,
#'   * derives `pseudotime_demo` from the UMAP-distance to the centroid
#'     of cluster 0 (same recipe as `mock_dataset()`), so the Trajectory
#'     module's `metadata` source has a non-empty option.
#'
#' All fabrications are deterministic for a given `seed`. No remote
#' network access; no Seurat call.
#'
#' Pure on `ds`; safe to unit-test.
.augment_demo_dataset <- function(ds, seed = 8L) {
  set.seed(seed)
  n  <- ds$n_cells
  cd <- ds$cell_data

  # ---- cluster (canonical column name expected by trajectory / regulons)
  if (!"cluster" %in% names(cd)) {
    cl_col <- intersect(c("seurat_clusters", "Cluster", "louvain",
                          "leiden"),
                        names(cd))
    if (length(cl_col) == 0L) {
      cl_col <- grep("cluster", names(cd), ignore.case = TRUE, value = TRUE)
    }
    if (length(cl_col) > 0L) {
      cd$cluster <- as.character(cd[[cl_col[1]]])
    } else {
      cd$cluster <- "0"   # last-resort uniform; keeps the schema valid
    }
  } else {
    cd$cluster <- as.character(cd$cluster)
  }
  cluster_ids <- sort(unique(cd$cluster), na.last = TRUE)

  # ---- cell_type placeholder
  if (!"cell_type" %in% names(cd)) {
    cd$cell_type <- paste0("Cluster ", cd$cluster)
  } else {
    cd$cell_type <- as.character(cd$cell_type)
  }

  # ---- condition / sample (fabricated; deterministic via seed)
  if (!"condition" %in% names(cd)) {
    cd$condition <- sample(c("ctrl", "treat"), n, replace = TRUE,
                           prob = c(0.55, 0.45))
  }
  if (!"sample" %in% names(cd)) {
    cd$sample <- sample(c("S1", "S2", "S3"), n, replace = TRUE,
                        prob = c(0.4, 0.35, 0.25))
  }

  # ---- pseudotime_demo: UMAP-distance from cluster 0 centroid, scaled 0-1.
  # Falls back to a random uniform when no UMAP is available so the
  # Trajectory module's "metadata" source still has a non-empty option.
  if (!"pseudotime_demo" %in% names(cd)) {
    if (all(c("UMAP_1", "UMAP_2") %in% names(cd))) {
      root_cells <- if (length(cluster_ids) > 0L)
        cd$cluster == cluster_ids[1] else rep(TRUE, n)
      if (!any(root_cells)) root_cells <- rep(TRUE, n)
      cx <- mean(cd$UMAP_1[root_cells])
      cy <- mean(cd$UMAP_2[root_cells])
      raw <- sqrt((cd$UMAP_1 - cx)^2 + (cd$UMAP_2 - cy)^2)
      cd$pseudotime_demo <- (raw - min(raw)) / max(diff(range(raw)), 1e-9)
    } else {
      cd$pseudotime_demo <- stats::runif(n)
    }
  }

  ds$cell_data <- cd

  # ---- Reorder metadata_fields so the workspace's default first-field
  # behaviour matches the mock dataset experience. Append any extra
  # fields the loader created (n_counts, n_features, etc.) in a stable
  # order at the end.
  preferred <- c("sample", "cluster", "condition", "cell_type",
                 "pseudotime_demo")
  embed_cols <- unlist(lapply(ds$reductions %||% character(),
                              function(r) paste0(r, c("_1", "_2"))))
  meta_cols  <- setdiff(names(cd), c("cell", embed_cols))
  extras     <- setdiff(meta_cols, preferred)
  ds$metadata_fields <- intersect(c(preferred, extras), names(cd))

  # ---- Update top-level name + source provenance.
  ds$source <- "demo_pbmc8k"
  ds
}
