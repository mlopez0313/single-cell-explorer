# ============================================================================
# Annotation engine registry
# ----------------------------------------------------------------------------
# Mirrors `IMPUTATION_METHODS` and the trajectory method scaffolding. Every
# annotation route (manual, marker-score, future SingleR / Azimuth /
# CellTypist / consensus / GPT) is declared as one entry and invoked through
# the same dispatcher. The annotation UI builds controls from the spec, not
# from `if (engine == "...")` branches.
#
# Engine contract:
#
#   spec$run_fn(dataset, params, state)
#     - returns a list with at minimum:
#         cell                    character(n_cells)
#         cell_labels             character(n_cells)
#         cell_scores             numeric(n_cells)
#       and optionally:
#         alt_labels              data.frame(cell, rank, label, score)
#         cluster_summary         data.frame
#         cluster_field_used      character(1)
#         n_clusters_at_creation  integer(1)
#         scores_matrix           matrix (engine-specific, kept under params)
#         warnings                character()
#         reference_source        character(1)
#     - must NOT mutate `state` or `dataset`. State is read-only here.
#
# `run_annotation_engine()` (in R/annotation.R) wraps the engine output into
# an `annotation_result_v1` object.
#
# Phase 2 engines that will plug into this same registry without UI changes:
#   - singler, azimuth, celltypist (reference-based)
#   - sctype, ucell (marker-score family; just different scoring functions)
#   - consensus (meta-engine that reads other results)
#   - gpt (LLM-assisted; output includes free-text rationale in warnings)
# ============================================================================

#' Specification for one annotation engine.
#'
#' @param id              snake_case identifier
#' @param name            display name for the UI
#' @param category        free-text grouping label, e.g. "manual",
#'                        "marker-score", "reference-based"
#' @param requires        character() declared resources: subset of
#'                        c("dataset", "cluster_field", "marker_registry",
#'                          "reference_model", "py_runtime")
#' @param produces        character() what this engine fills, subset of
#'                        c("per_cell_labels", "per_cell_scores",
#'                          "alt_labels", "cluster_summary")
#' @param parameters      named list of parameter specs; each spec is a
#'                        list(type, default, required, choices, ...).
#'                        Consumed by the UI to build controls.
#' @param run_fn          function(dataset, params, state) -> list (see above)
#' @param version         character(1) engine implementation version. Stamped
#'                        onto results as `engine_version` (provenance).
#'                        This is the version of the *engine code/model*,
#'                        e.g. "marker_score_v1.0.0" or "azimuth_v0.5.0";
#'                        it is distinct from `result_schema`, which is
#'                        the version of the result *data shape*.
#'                        Bumped independently of `result_schema`.
#' @param result_schema   schema id that `run_fn`'s output adheres to.
#'                        DIFFERENT from `version`: this tracks the data
#'                        shape, `version` tracks the implementation.
#' @param is_demo         logical(1) results from this engine are stamped is_demo=TRUE
#' @param enabled         logical(1) FALSE hides from the UI
#' @param description     character(1)
annotation_engine_spec <- function(id, name, category,
                                   requires       = character(),
                                   produces       = c("per_cell_labels"),
                                   parameters     = list(),
                                   run_fn,
                                   version        = NULL,
                                   result_schema  = ANNOTATION_RESULT_SCHEMA_VERSION,
                                   is_demo        = FALSE,
                                   enabled        = TRUE,
                                   description    = "") {
  stopifnot(is.character(id), length(id) == 1L, nzchar(id))
  stopifnot(is.function(run_fn))
  if (is.null(version) || !nzchar(version))
    stop(sprintf(
      "annotation_engine_spec('%s'): `version` must be a non-empty string. ",
      id),
      "Engine version is provenance, not optional. ",
      "Use e.g. '<engine_id>_v1.0.0'.",
      call. = FALSE)
  list(
    id            = id,
    name          = name,
    category      = category,
    requires      = requires,
    produces      = produces,
    parameters    = parameters,
    run_fn        = run_fn,
    version       = as.character(version),
    result_schema = result_schema,
    is_demo       = isTRUE(is_demo),
    enabled       = isTRUE(enabled),
    description   = description
  )
}

#' Per-engine GitHub install map.
#'
#' Named list, keys are engine ids, values are NAMED character vectors
#' where names are installed package names and values are
#' `remotes::install_github()` specs (`owner/repo[@ref]`). Used by the
#' annotation module to decide whether to pop the "Install + run"
#' confirmation modal before invoking an engine.
#'
#' Engines whose dependencies live entirely on CRAN / Bioconductor (e.g.
#' `manual`, `marker_score`, `singler`, `celltypist`) are intentionally
#' absent. Add an entry here when introducing a new engine that needs a
#' GitHub-only package.
ANNOTATION_ENGINE_GITHUB_DEPS <- function() {
  list(
    azimuth = c(Azimuth = "satijalab/azimuth")
  )
}

#' Missing GitHub deps for a given engine id (or `character()` if none).
#'
#' Returns a named character vector keyed by package name (same shape
#' as `ANNOTATION_ENGINE_GITHUB_DEPS[[engine_id]]`) restricted to the
#' subset that is NOT currently installed. Empty result = "ready to
#' run".
engine_missing_github_deps <- function(engine_id) {
  spec <- ANNOTATION_ENGINE_GITHUB_DEPS()[[engine_id]]
  if (is.null(spec) || !length(spec)) return(character())
  installed <- vapply(names(spec), has_optional, logical(1))
  spec[!installed]
}

#' The full registry. Function so it picks up changes during dev / testing.
ANNOTATION_ENGINES <- function() {
  list(
    annotation_engine_spec(
      id          = "manual",
      name        = "Manual annotation",
      category    = "manual",
      requires    = c("dataset", "cluster_field"),
      produces    = c("per_cell_labels", "cluster_summary"),
      parameters  = list(
        cluster_field = list(type = "metadata_field", required = TRUE),
        labels        = list(type = "named_character", required = TRUE,
                             description = "Map of cluster id -> label.")
      ),
      run_fn      = .run_manual_annotation,
      version     = "manual_v1.0.0",
      description = "User-supplied per-cluster labels, expanded to per-cell."
    ),
    annotation_engine_spec(
      id          = "marker_score",
      name        = "Marker-score (registry-driven)",
      category    = "marker-score",
      requires    = c("dataset", "cluster_field", "marker_registry"),
      produces    = c("per_cell_labels", "per_cell_scores",
                      "alt_labels", "cluster_summary"),
      parameters  = list(
        cluster_field = list(type = "metadata_field",
                             required = TRUE),
        species       = list(type = "select", required = FALSE,
                             default = NA_character_),
        tissue        = list(type = "select", required = FALSE,
                             default = NA_character_),
        min_score     = list(type = "numeric", required = FALSE,
                             default = 0.0,
                             description = "Top score below this -> Unknown.")
      ),
      run_fn      = .run_marker_score_annotation,
      version     = "marker_score_v1.0.0",
      description = paste("Score cell-type marker panels from the marker",
                          "registry per cluster; pick the top label.")
    ),
    annotation_engine_spec(
      id          = "singler",
      name        = "SingleR (reference-based)",
      category    = "reference-based",
      requires    = c("dataset"),
      produces    = c("per_cell_labels", "per_cell_scores",
                      "alt_labels", "cluster_summary"),
      parameters  = list(
        reference     = list(type = "select", required = TRUE,
                             default = "hpca",
                             choices = c("hpca", "blueprint_encode",
                                         "monaco_immune"),
                             description = paste(
                               "celldex reference panel.",
                               "hpca = HumanPrimaryCellAtlas (broad);",
                               "blueprint_encode = bulk hematopoietic;",
                               "monaco_immune = sorted PBMC.")),
        labels        = list(type = "select", required = FALSE,
                             default = "main",
                             choices = c("main", "fine"),
                             description = "Granularity of reference labels."),
        cluster_field = list(type = "metadata_field", required = FALSE,
                             description = paste(
                               "If set, SingleR runs in cluster mode",
                               "(one prediction per cluster, expanded back",
                               "to per-cell). Much faster on > 50k cells.")),
        min_delta     = list(type = "numeric", required = FALSE,
                             default = 0.0,
                             description = paste(
                               "Minimum delta.next confidence margin;",
                               "below this the label becomes 'Unknown'."))
      ),
      run_fn      = .run_singler_annotation,
      version     = "singler_v1.0.0",
      description = paste(
        "Reference-based cell-type prediction via SingleR.",
        "Requires Bioconductor packages SingleR + celldex.")
    ),
    annotation_engine_spec(
      id          = "azimuth",
      name        = "Azimuth (Seurat reference mapping)",
      category    = "reference-based",
      requires    = c("dataset"),
      produces    = c("per_cell_labels", "per_cell_scores",
                      "cluster_summary"),
      parameters  = list(
        reference         = list(type = "select", required = TRUE,
                                 default = "pbmcref",
                                 choices = c("pbmcref", "lungref",
                                             "kidneyref", "bonemarrowref",
                                             "heartref", "fetusref"),
                                 description = paste(
                                   "Tissue-specific reference shipped",
                                   "as a separate `Azimuth.<ref>`",
                                   "data package on the Satija lab",
                                   "server. The first run downloads",
                                   "several GB.")),
        annotation_level  = list(type = "select", required = FALSE,
                                 default = "celltype.l2",
                                 choices = c("celltype.l1", "celltype.l2",
                                             "celltype.l3"),
                                 description = paste(
                                   "Label granularity. Not every",
                                   "reference exposes all three levels.")),
        cluster_field     = list(type = "metadata_field", required = FALSE,
                                 description = paste(
                                   "If set, used to build a per-cluster",
                                   "summary (does NOT change the",
                                   "Azimuth call -- predictions are",
                                   "always per-cell).")),
        min_mapping_score = list(type = "numeric", required = FALSE,
                                 default = 0.0,
                                 description = paste(
                                   "Cells with `mapping.score` below",
                                   "this become 'Unknown' with score 0.",
                                   "Mirrors SingleR's `min_delta` gate."))
      ),
      run_fn      = .run_azimuth_annotation,
      version     = "azimuth_v1.0.0",
      description = paste(
        "Reference-based cell-type prediction via Azimuth.",
        "Requires Seurat + Azimuth + the tissue-specific reference",
        "package (e.g. `Azimuth.pbmcref`).")
    ),
    annotation_engine_spec(
      id          = "celltypist",
      name        = "CellTypist (Python via reticulate)",
      category    = "reference-based",
      requires    = c("dataset", "py_runtime"),
      produces    = c("per_cell_labels", "per_cell_scores",
                      "cluster_summary"),
      parameters  = list(
        model            = list(type = "select", required = TRUE,
                                default = "Immune_All_Low.pkl",
                                choices = c("Immune_All_Low.pkl",
                                            "Immune_All_High.pkl",
                                            "COVID19_HumanChallenge_Blood.pkl",
                                            "Cells_Lung_Airway.pkl",
                                            "Cells_Intestinal_Tract.pkl",
                                            "Adult_Mouse_Gut.pkl",
                                            "Healthy_COVID19_PBMC.pkl"),
                                description = paste(
                                  "Pretrained CellTypist model.",
                                  "Auto-downloads to",
                                  "~/.celltypist/models/ on first use.",
                                  "Full catalogue at celltypist.org.")),
        majority_voting  = list(type = "logical", required = FALSE,
                                default = FALSE,
                                description = paste(
                                  "Aggregate predictions over a cluster",
                                  "(per CellTypist docs). Recommended",
                                  "for noisy / sparse data.")),
        over_clustering  = list(type = "metadata_field", required = FALSE,
                                description = paste(
                                  "Metadata column to use as the",
                                  "majority-voting reference. Only",
                                  "consulted when `majority_voting` is",
                                  "TRUE.")),
        cluster_field    = list(type = "metadata_field", required = FALSE,
                                description = paste(
                                  "If set, used to build a per-cluster",
                                  "summary. Independent of",
                                  "`over_clustering`.")),
        min_score        = list(type = "numeric", required = FALSE,
                                default = 0.0,
                                description = paste(
                                  "Cells with `conf_score` below this",
                                  "become 'Unknown' with score 0."))
      ),
      run_fn      = .run_celltypist_annotation,
      version     = "celltypist_v1.0.0",
      description = paste(
        "Reference-based cell-type prediction via the CellTypist",
        "Python package, bridged through reticulate. Requires",
        "reticulate + anndata (R) + celltypist (Python).")
    )
  )
}

#' Look up a single engine by id. NULL if not registered.
get_annotation_engine <- function(id) {
  for (e in ANNOTATION_ENGINES()) if (identical(e$id, id)) return(e)
  NULL
}

#' Enabled engines as a named character vector: display_name -> id.
#' Suitable for `selectInput(choices = ...)`.
list_annotation_engines <- function(enabled_only = TRUE) {
  engines <- ANNOTATION_ENGINES()
  if (enabled_only) engines <- Filter(function(e) isTRUE(e$enabled), engines)
  setNames(
    vapply(engines, `[[`, character(1), "id"),
    vapply(engines, `[[`, character(1), "name")
  )
}

# ============================================================================
# Built-in engine implementations
# ----------------------------------------------------------------------------
# Kept in the same file as the registry so they live and die together.
# Pure functions: no shiny, no state mutation. State is read-only.
# ============================================================================

# ---- manual ---------------------------------------------------------------

.run_manual_annotation <- function(dataset, params, state) {
  cluster_field <- params$cluster_field
  labels        <- params$labels %||% list()
  if (is.null(cluster_field) || !nzchar(cluster_field))
    stop("Manual engine requires a cluster_field.", call. = FALSE)
  cluster_vec <- get_metadata(dataset, cluster_field)
  if (is.null(cluster_vec))
    stop(sprintf("Cluster field '%s' not in dataset.", cluster_field), call. = FALSE)
  cluster_vec <- as.character(cluster_vec)
  cluster_ids <- sort(unique(cluster_vec))

  # Normalise labels to named character. Blank / NA values mean "no
  # assignment for this cluster" and produce NA cells.
  norm <- setNames(rep(NA_character_, length(cluster_ids)), cluster_ids)
  for (cl in intersect(names(labels), cluster_ids)) {
    v <- labels[[cl]]
    if (!is.null(v) && length(v) == 1L && !is.na(v) && nzchar(trimws(as.character(v)))) {
      norm[cl] <- trimws(as.character(v))
    }
  }

  cell_labels <- unname(norm[cluster_vec])
  cell_scores <- ifelse(is.na(cell_labels) | !nzchar(cell_labels), 0.0, 1.0)

  cluster_summary <- data.frame(
    cluster   = cluster_ids,
    top_label = unname(norm[cluster_ids]),
    top_score = ifelse(is.na(unname(norm[cluster_ids])), 0.0, 1.0),
    n_cells   = vapply(cluster_ids, function(cl) sum(cluster_vec == cl), integer(1)),
    stringsAsFactors = FALSE
  )

  list(
    cell                   = dataset$cell_data$cell,
    cell_labels            = cell_labels,
    cell_scores            = cell_scores,
    alt_labels             = NULL,
    cluster_summary        = cluster_summary,
    cluster_field_used     = cluster_field,
    n_clusters_at_creation = length(cluster_ids),
    reference_source       = "user"
  )
}

# ---- marker_score ---------------------------------------------------------
# Score each cell type in the marker registry against each cluster using:
#   score = sum(weight_pos * mean_expr_in_cluster) -
#           sum(weight_neg * mean_expr_in_cluster)
# Then pick the top cell type per cluster (Unknown if best score < min_score).
# Mean expression is the available proxy for per-cluster prevalence in the
# mock; real implementations should swap in a proper scoring backend
# (UCell / scType / AddModuleScore) registered as additional engines.
#
# TODO (real backends, slot in via ANNOTATION_ENGINES):
#   - ucell (rank-based) - lighter, robust to depth
#   - sctype (mean-expression with FC-based weighting)
#   - module_score (Seurat AddModuleScore equivalent)

.run_marker_score_annotation <- function(dataset, params, state) {
  cluster_field <- params$cluster_field
  if (is.null(cluster_field) || !nzchar(cluster_field))
    stop("marker_score engine requires a cluster_field.", call. = FALSE)
  registry <- state$marker_registry
  if (is.null(registry))
    stop("marker_score engine requires `state$marker_registry`.", call. = FALSE)

  species   <- params$species
  tissue    <- params$tissue
  min_score <- params$min_score %||% 0.0

  # Filter registry; respect NA filters as "no filter on this axis".
  species_f <- if (is.null(species) || is.na(species) || !nzchar(species)) NULL else species
  tissue_f  <- if (is.null(tissue)  || is.na(tissue)  || !nzchar(tissue))  NULL else tissue
  entries   <- marker_registry_filter(registry, species = species_f, tissue = tissue_f)
  if (!length(entries))
    stop("Marker registry returned no entries for the chosen species/tissue.", call. = FALSE)

  cluster_vec <- as.character(get_metadata(dataset, cluster_field))
  if (is.null(cluster_vec))
    stop(sprintf("Cluster field '%s' not in dataset.", cluster_field), call. = FALSE)
  cluster_ids <- sort(unique(cluster_vec))
  cell_types  <- vapply(entries, `[[`, character(1), "cell_type")
  warnings    <- character()
  if (length(cluster_ids) <= 1L) {
    # Surfaced from real-data smoke testing on a Seurat object where
    # FindClusters collapsed everything into cluster `0`. The
    # marker_score formula compares between-cluster signal -- given
    # one bucket of mixed-identity cells, the engine returns whatever
    # cell type has the cleanest single-marker panel (often the one
    # with the fewest negative markers). Make the failure mode visible.
    warnings <- c(warnings, sprintf(
      "marker_score engine: cluster_field '%s' has %d cluster(s). ",
      cluster_field, length(cluster_ids)),
      "Results compare cells against the dataset-wide mean and will ",
      "favour cell types with small, non-overlapping marker panels. ",
      "Re-cluster at a finer resolution before trusting the labels.")
  }

  scores_mat <- matrix(0,
                       nrow = length(cluster_ids),
                       ncol = length(cell_types),
                       dimnames = list(cluster_ids, cell_types))

  for (j in seq_along(entries)) {
    e <- entries[[j]]
    for (m in e$markers) {
      expr <- get_gene_expression(dataset, m$gene)
      if (is.null(expr)) {
        warnings <- c(warnings,
          sprintf("Marker gene '%s' not in dataset; skipped for '%s'.",
                  m$gene, e$cell_type))
        next
      }
      for (i in seq_along(cluster_ids)) {
        in_cl <- cluster_vec == cluster_ids[i]
        if (!any(in_cl)) next
        mu <- mean(expr[in_cl], na.rm = TRUE)
        sign <- switch(m$role,
          "positive" = +1, "negative" = -1, "specific" = +1, 0)
        scores_mat[i, j] <- scores_mat[i, j] + sign * m$weight * mu
      }
    }
  }

  # Per-cluster best label + alt_labels (top-3).
  cluster_summary <- data.frame(
    cluster      = cluster_ids,
    top_label    = NA_character_,
    top_score    = NA_real_,
    second_label = NA_character_,
    second_score = NA_real_,
    n_cells      = vapply(cluster_ids, function(cl) sum(cluster_vec == cl), integer(1)),
    stringsAsFactors = FALSE
  )
  alt_rows <- list()
  for (i in seq_along(cluster_ids)) {
    cl  <- cluster_ids[i]
    sc  <- scores_mat[i, ]
    ord <- order(sc, decreasing = TRUE)
    best <- cell_types[ord[1]]
    if (sc[ord[1]] < min_score) best <- "Unknown"
    cluster_summary$top_label[i]    <- best
    cluster_summary$top_score[i]    <- sc[ord[1]]
    if (length(ord) >= 2L) {
      cluster_summary$second_label[i] <- cell_types[ord[2]]
      cluster_summary$second_score[i] <- sc[ord[2]]
    }
    # alt_labels: top-3 candidates per cluster, replicated to every cell of
    # that cluster downstream. For now store at cluster granularity (we
    # forward-project only if/when a per-cell alt-label view is needed).
    top_k <- utils::head(ord, 3L)
    for (k in seq_along(top_k)) {
      alt_rows[[length(alt_rows) + 1L]] <- data.frame(
        cluster = cl, rank = k,
        label   = cell_types[top_k[k]],
        score   = sc[top_k[k]],
        stringsAsFactors = FALSE
      )
    }
  }
  alt_labels <- if (length(alt_rows)) do.call(rbind, alt_rows) else NULL

  cluster_to_label <- setNames(cluster_summary$top_label, cluster_summary$cluster)
  cluster_to_score <- setNames(cluster_summary$top_score, cluster_summary$cluster)
  cell_labels      <- unname(cluster_to_label[cluster_vec])
  cell_scores      <- unname(cluster_to_score[cluster_vec])
  # Cells in clusters scored below min_score got "Unknown" -> score 0.
  cell_scores[cell_labels == "Unknown" | is.na(cell_labels)] <- 0.0

  list(
    cell                   = dataset$cell_data$cell,
    cell_labels            = cell_labels,
    cell_scores            = cell_scores,
    alt_labels             = alt_labels,
    cluster_summary        = cluster_summary,
    cluster_field_used     = cluster_field,
    n_clusters_at_creation = length(cluster_ids),
    reference_source       = sprintf("marker_registry/%s",
                                     registry$version %||% "unknown"),
    scores_matrix          = scores_mat,
    warnings               = warnings
  )
}
