# ============================================================================
# Annotation engine: Azimuth
# ----------------------------------------------------------------------------
# Reference-based per-cell annotation via the Satija lab's Azimuth
# package. Pluggable through `ANNOTATION_ENGINES()`; the existing
# Annotation module UI renders parameters from the spec, so no module
# changes are needed.
#
# Layout (mirrors annotation_singler.R deliberately):
#   * `.run_azimuth_annotation()`     -- the engine `run_fn`. Heavy deps
#                                         (Azimuth + Seurat + the
#                                         tissue-specific reference
#                                         data package, e.g.
#                                         `Azimuth.pbmcref`) gated
#                                         through `require_optional`.
#                                         Wraps the dataset in a
#                                         minimal Seurat object,
#                                         calls `Azimuth::RunAzimuth`,
#                                         routes the result through
#                                         the pure converter.
#   * `.azimuth_to_engine_output()`   -- pure schema converter. Takes
#                                         any data.frame shaped like
#                                         Azimuth's `meta.data` output
#                                         (per-cell `predicted.celltype.X`
#                                         + `predicted.celltype.X.score`
#                                         columns + `mapping.score`)
#                                         and lifts it into the engine
#                                         contract. Testable without
#                                         Azimuth installed.
#   * `.dataset_to_seurat_minimal()`  -- build a minimal Seurat object
#                                         from the active dataset
#                                         (counts layer + cell_data),
#                                         enough for Azimuth to map
#                                         against a reference.
#
# Caching: Azimuth downloads multi-GB reference data the first time a
# reference is used; the package itself caches via BiocFileCache so we
# don't add another layer.
#
# References available (Azimuth ships these as separate data packages):
#   - pbmcref          (PBMC, level 1/2/3 labels)
#   - lungref          (Human Lung Cell Atlas)
#   - kidneyref        (kidney)
#   - bonemarrowref    (bone marrow)
#   - heartref         (heart)
#   - fetusref         (developmental atlas)
# Real userse install only the references they need.
# ============================================================================

# ---- Engine run function ---------------------------------------------------

.run_azimuth_annotation <- function(dataset, params, state) {
  require_optional(
    c("Seurat", "SeuratObject"),
    feature = "Azimuth annotation engine (CRAN deps)",
    source  = "CRAN")
  require_optional(
    "Azimuth",
    feature = "Azimuth annotation engine (GitHub package)",
    source  = "GitHub",
    repo    = c(Azimuth = "satijalab/azimuth"))

  # Azimuth's internals reach for S4 setters from `SeuratObject` (most
  # notably `Key<-`) via the calling-environment search path rather
  # than through its NAMESPACE imports. The `requireNamespace()` calls
  # above load both packages but do NOT attach them, so the call into
  # `RunAzimuth` fails with `could not find function "Key<-"`. Attach
  # both packages here (idempotent if already on the search path) --
  # this is the same effect as `library(SeuratObject); library(Seurat)`
  # but scoped to engines that actually need it.
  ensure_attached("SeuratObject")
  ensure_attached("Seurat")

  reference         <- params$reference         %||% "pbmcref"
  annotation_level  <- params$annotation_level  %||% "celltype.l2"
  min_mapping_score <- params$min_mapping_score %||% 0.0
  cluster_field     <- params$cluster_field

  cells <- dataset$cell_data$cell

  # Build the minimal Seurat object Azimuth needs. We use the `counts`
  # layer (the supported way to feed Azimuth) and fall back to the
  # default layer with a warning if `counts` is missing.
  seu_obj <- .dataset_to_seurat_minimal(dataset)

  # Azimuth's main entry point. We invoke through getExportedValue so
  # the symbol lookup is explicit (helps tools that can't see Suggests).
  RunAzimuth_fn <- get("RunAzimuth", envir = asNamespace("Azimuth"))
  mapped <- RunAzimuth_fn(query = seu_obj, reference = reference)

  # `mapped` is a Seurat object whose @meta.data carries the new
  # prediction columns. Pull it out as a plain data.frame so the
  # converter never depends on Seurat slots.
  md <- as.data.frame(SeuratObject::FetchData(
    mapped, vars = .azimuth_meta_cols(mapped, annotation_level)),
    stringsAsFactors = FALSE)
  md$cell <- rownames(md)

  cluster_vec <- if (is.null(cluster_field) || !nzchar(cluster_field)) NULL
                 else as.character(get_metadata(dataset, cluster_field))

  .azimuth_to_engine_output(
    md,
    cells              = cells,
    annotation_level   = annotation_level,
    min_mapping_score  = min_mapping_score,
    reference          = reference,
    cluster_field_used = cluster_field,
    cluster_vec        = cluster_vec
  )
}

# Build a minimal Seurat object from our dataset. Just enough for
# `Azimuth::RunAzimuth` to anchor / transfer. Counts layer is the
# canonical input.
.dataset_to_seurat_minimal <- function(dataset) {
  be <- as_expression_backend(dataset$expression)
  layers <- backend_available_layers(be)
  layer <- if ("counts" %in% layers) "counts" else backend_default_layer(be)
  M <- backend_as_matrix(be, layer = layer)
  if (is.null(rownames(M)))
    rownames(M) <- backend_genes(be, layer = layer)
  if (is.null(colnames(M)))
    colnames(M) <- dataset$cell_data$cell

  create_obj <- get("CreateSeuratObject", envir = asNamespace("SeuratObject"))
  meta <- dataset$cell_data
  rownames(meta) <- dataset$cell_data$cell
  seu <- create_obj(counts = M, meta.data = meta, assay = "RNA")
  seu
}

# Choose the per-cell columns we want to read back from Azimuth's
# mapped meta.data. Defensive: only request columns that actually exist
# in the result (Azimuth's column names depend on the reference).
.azimuth_meta_cols <- function(seu_mapped, annotation_level) {
  want <- c(
    paste0("predicted.", annotation_level),
    paste0("predicted.", annotation_level, ".score"),
    "mapping.score"
  )
  have <- intersect(want, colnames(seu_mapped[[]]))
  if (length(have) == 0L)
    stop(sprintf(
      "Azimuth result lacks expected columns for annotation_level = '%s'. ",
      annotation_level),
      "Available: ",
      paste(grep("^predicted\\.", colnames(seu_mapped[[]]),
                 value = TRUE), collapse = ", "),
      call. = FALSE)
  have
}

#' Pure schema converter: Azimuth-style meta.data -> engine output.
#'
#' Accepts a data.frame with the standard Azimuth columns for the chosen
#' annotation level, plus a `cell` column carrying barcodes. Returns the
#' engine output contract documented in R/annotation_registry.R.
#'
#' Factored out so the schema mapping can be regression-tested without
#' Azimuth or its reference data packages installed.
#'
#' @param df                 data.frame with `cell` + at least one of
#'                           `predicted.<level>` and (optionally)
#'                           `predicted.<level>.score`, `mapping.score`.
#' @param cells              character(n_cells) -- canonical cell order.
#'                           Used to align rows of `df` to the dataset.
#' @param annotation_level   character(1) -- which `predicted.<X>` column
#'                           to use.
#' @param min_mapping_score  numeric(1) -- cells with `mapping.score`
#'                           below this become "Unknown" with score 0.
#'                           If `mapping.score` is absent in `df`, no
#'                           filtering happens.
#' @param reference          character(1) -- stamped on `reference_source`.
#' @param cluster_field_used character(1) or NULL.
#' @param cluster_vec        character(n_cells) or NULL -- per-cell
#'                           cluster ids used to summarise per-cluster.
.azimuth_to_engine_output <- function(df,
                                      cells,
                                      annotation_level   = "celltype.l2",
                                      min_mapping_score  = 0.0,
                                      reference          = NA_character_,
                                      cluster_field_used = NULL,
                                      cluster_vec        = NULL) {
  if (!is.data.frame(df))
    stop(".azimuth_to_engine_output: `df` must be a data.frame.",
         call. = FALSE)
  if (!"cell" %in% names(df))
    stop(".azimuth_to_engine_output: `df` needs a `cell` column with ",
         "barcode ids.", call. = FALSE)

  label_col <- paste0("predicted.", annotation_level)
  score_col <- paste0("predicted.", annotation_level, ".score")
  if (!label_col %in% names(df))
    stop(sprintf(
      ".azimuth_to_engine_output: missing column '%s'. Have: %s",
      label_col, paste(names(df), collapse = ", ")), call. = FALSE)

  # Align to the canonical cell order; cells absent from `df` -> NA.
  idx <- match(cells, df$cell)
  if (any(is.na(idx))) {
    # All cells must be in df for a real Azimuth run. Surface a
    # diagnostic rather than silently inserting NA labels.
    n_missing <- sum(is.na(idx))
    stop(sprintf(
      ".azimuth_to_engine_output: %d cells missing from Azimuth result.",
      n_missing), call. = FALSE)
  }

  cell_labels <- as.character(df[[label_col]][idx])
  cell_scores <- if (score_col %in% names(df))
                   as.numeric(df[[score_col]][idx])
                 else
                   rep(NA_real_, length(cells))

  mapping_scores <- if ("mapping.score" %in% names(df))
                      as.numeric(df$mapping.score[idx])
                    else
                      rep(NA_real_, length(cells))

  # Confidence filter: cells whose mapping.score is below the cutoff
  # become "Unknown" with score 0. Comparable to SingleR's min_delta
  # gating.
  if (is.finite(min_mapping_score) && min_mapping_score > 0) {
    low <- !is.na(mapping_scores) & mapping_scores < min_mapping_score
    cell_labels[low] <- "Unknown"
    cell_scores[low] <- 0
  }

  cluster_summary <- if (!is.null(cluster_vec) &&
                         length(cluster_vec) == length(cells)) {
    .summarise_per_cluster_azimuth(
      cluster_vec = as.character(cluster_vec),
      labels      = cell_labels,
      scores      = cell_scores)
  } else NULL

  list(
    cell                   = cells,
    cell_labels            = cell_labels,
    cell_scores            = cell_scores,
    alt_labels             = NULL,
    cluster_summary        = cluster_summary,
    cluster_field_used     = cluster_field_used,
    n_clusters_at_creation = if (is.null(cluster_vec)) 0L
                             else length(unique(as.character(cluster_vec))),
    reference_source       = sprintf("Azimuth:%s", reference)
  )
}

# Internal: top-label / fraction summary per cluster.
.summarise_per_cluster_azimuth <- function(cluster_vec, labels, scores) {
  cluster_ids <- sort(unique(cluster_vec))
  out <- vector("list", length(cluster_ids))
  for (i in seq_along(cluster_ids)) {
    cl <- cluster_ids[i]
    in_cl <- which(cluster_vec == cl)
    tab <- sort(table(labels[in_cl]), decreasing = TRUE)
    top  <- if (length(tab)) names(tab)[1] else NA_character_
    frac <- if (length(tab)) as.numeric(tab[1] / length(in_cl)) else NA_real_
    out[[i]] <- data.frame(
      cluster      = cl,
      top_label    = top,
      top_fraction = frac,
      mean_score   = mean(scores[in_cl], na.rm = TRUE),
      n_cells      = length(in_cl),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}
