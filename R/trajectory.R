# ============================================================================
# Trajectory / Pseudotime
# ----------------------------------------------------------------------------
# Treats pseudotime as an exploratory analysis result. Lives ONLY in
# `state$analysis_results$trajectory` -- it is NOT silently injected
# into DE / Markers / Pathway / Imputation. Those modules read raw
# expression and standard metadata; pseudotime is opt-in per module if
# ever needed.
#
# Pluggable methods live in `TRAJECTORY_METHODS()` (see
# R/trajectory_registry.R). This file is the orchestrator + the
# `mock` / `metadata` `run_fn`s + shared helpers (binning, rescaling,
# summary). Real backends live in their own files
# (R/trajectory_slingshot.R, R/trajectory_monocle3.R) so each can be
# added or removed without touching this orchestrator.
# ============================================================================

#' Available pseudotime sources -- back-compat alias for the registry.
#'
#' Returns a named character vector `label -> id` for every registered
#' method (whether or not its optional packages are installed).
#' Unavailable methods are annotated so users can see why a backend is
#' selectable but errors at run time.
#'
#' Old code paths that read `PSEUDOTIME_SOURCES` directly continue to
#' work; the canonical accessor is `trajectory_method_choices()`.
PSEUDOTIME_SOURCES <- function() trajectory_method_choices()

available_pseudotime_sources <- function() trajectory_method_choices()

#' Metadata fields that look numeric and could be used as a pseudotime source.
#'
#' Restricted to columns listed in `dataset$metadata_fields` AND present in
#' `cell_data` AND of numeric type. Returns character() if none exist.
available_numeric_metadata_fields <- function(dataset) {
  if (is.null(dataset)) return(character())
  fields <- intersect(dataset$metadata_fields, names(dataset$cell_data))
  Filter(function(f) is.numeric(dataset$cell_data[[f]]), fields)
}

#' Categorical metadata fields: character / factor / integer with few unique
#' levels. Used as cluster-like grouping pickers (Regulons, Annotation).
available_categorical_metadata_fields <- function(dataset) {
  if (is.null(dataset)) return(character())
  fields <- intersect(dataset$metadata_fields, names(dataset$cell_data))
  Filter(function(f) {
    v <- dataset$cell_data[[f]]
    is.character(v) || is.factor(v) ||
      (is.integer(v) && length(unique(v)) <= 50L)
  }, fields)
}

#' Run a trajectory method end-to-end. Returns the full `results` payload.
#'
#' Walks the registry, calls the method's `run_fn`, and validates the
#' return shape. The output is the value that gets stored in
#' `state$analysis_results$trajectory$results`.
#'
#' @param dataset      the active dataset.
#' @param source       a registered method id (`list_trajectory_methods()`).
#' @param ...          method-specific parameters, forwarded to `run_fn`.
#'                     The orchestrator passes through `reduction`,
#'                     `root_field`, `root_group`, `metadata_field`,
#'                     `cluster_field`, and anything else the caller adds.
#'
#' @return list with at least `pseudotime`, `cell`, `source`,
#'   `reduction_used`, `root_field`, `root_group`, `metadata_field`,
#'   `n_lineages`, `method_details`. See
#'   `trajectory_method_spec()` for the contract.
run_trajectory <- function(dataset, source, ...) {
  if (is.null(dataset)) stop("No dataset provided.", call. = FALSE)
  method <- get_trajectory_method(source)
  if (is.null(method)) {
    stop(sprintf("Unknown trajectory method '%s'. Have: %s",
                 source, paste(list_trajectory_methods(), collapse = ", ")),
         call. = FALSE)
  }
  params <- list(...)
  out <- method$run_fn(dataset, params)
  if (is.null(out$pseudotime) || !is.numeric(out$pseudotime))
    stop(sprintf(
      "Trajectory method '%s' did not return a numeric `pseudotime`.",
      source), call. = FALSE)
  # Fill canonical fields if the run_fn omitted them; downstream UI
  # relies on every key being present (possibly NA).
  defaults <- list(
    cell           = dataset$cell_data$cell,
    source         = source,
    reduction_used = NA_character_,
    root_field     = NA_character_,
    root_group     = NA_character_,
    metadata_field = NA_character_,
    n_lineages     = 1L,
    method_details = list()
  )
  for (k in names(defaults)) {
    if (is.null(out[[k]])) out[[k]] <- defaults[[k]]
  }
  out
}

#' Compute a pseudotime vector of length n_cells.
#'
#' Back-compat wrapper around `run_trajectory()`. Useful when callers
#' only need the per-cell pseudotime and don't care about lineage
#' details (e.g. binning along time, gene-trend overlays).
#'
#' @param dataset         the active dataset
#' @param source          a registered method id (any of
#'                        `list_trajectory_methods()`)
#' @param reduction       reduction used by embedding-based methods
#'                        (mock / slingshot / monocle3)
#' @param root_field      categorical metadata field defining the root
#'                        group / starting cluster (mock / slingshot /
#'                        monocle3)
#' @param root_group      value within `root_field` that is the root
#' @param metadata_field  numeric metadata column (metadata source only)
#' @param cluster_field   cluster-label metadata column for methods that
#'                        need explicit cluster ids (slingshot). Defaults
#'                        to `root_field`.
#'
#' @return numeric vector aligned with `dataset$cell_data$cell`. Values
#'   in `[0, 1]`. Throws on invalid input.
compute_pseudotime <- function(dataset, source,
                               reduction = NULL,
                               root_field = NULL, root_group = NULL,
                               metadata_field = NULL,
                               cluster_field = NULL) {
  out <- run_trajectory(dataset, source,
                        reduction      = reduction,
                        root_field     = root_field,
                        root_group     = root_group,
                        metadata_field = metadata_field,
                        cluster_field  = cluster_field)
  out$pseudotime
}

# ---- Built-in run_fns ---------------------------------------------------
#
# Each `run_fn` consumes the same `params` list and returns the canonical
# result payload documented in `trajectory_method_spec()`. Keep these
# pure: no Shiny / state mutation.

.run_metadata_trajectory <- function(dataset, params) {
  metadata_field <- params$metadata_field
  if (is.null(metadata_field) || !nzchar(metadata_field))
    stop("metadata_field is required for source = 'metadata'.",
         call. = FALSE)
  v <- get_metadata(dataset, metadata_field)
  if (is.null(v))
    stop(sprintf("Metadata field '%s' is not available.", metadata_field),
         call. = FALSE)
  if (!is.numeric(v))
    stop(sprintf("Metadata field '%s' is not numeric.", metadata_field),
         call. = FALSE)
  list(
    pseudotime     = rescale01(v),
    cell           = dataset$cell_data$cell,
    source         = "metadata",
    reduction_used = NA_character_,
    root_field     = NA_character_,
    root_group     = NA_character_,
    metadata_field = metadata_field,
    n_lineages     = 1L,
    method_details = list()
  )
}

.run_mock_trajectory <- function(dataset, params) {
  red <- params$reduction %||% dataset$default_reduction %||%
         (dataset$reductions %||% character())[1]
  emb <- get_embedding(dataset, red)
  if (is.null(emb))
    stop(sprintf("No usable reduction '%s' for mock pseudotime.", red %||% ""),
         call. = FALSE)

  root_field <- params$root_field
  root_group <- params$root_group
  if (is.null(root_field) || !nzchar(root_field))
    stop("root_field is required for source = 'mock'.", call. = FALSE)
  rf <- get_metadata(dataset, root_field)
  if (is.null(rf))
    stop(sprintf("root_field '%s' is not available.", root_field),
         call. = FALSE)
  if (is.null(root_group) || !nzchar(as.character(root_group)))
    stop("root_group is required for source = 'mock'.", call. = FALSE)
  in_root <- as.character(rf) == as.character(root_group)
  if (sum(in_root, na.rm = TRUE) == 0L)
    stop(sprintf("Root group '%s' has no cells in field '%s'.",
                 root_group, root_field), call. = FALSE)

  centroid <- c(mean(emb$x[in_root]), mean(emb$y[in_root]))
  d <- sqrt((emb$x - centroid[1])^2 + (emb$y - centroid[2])^2)
  list(
    pseudotime     = rescale01(d),
    cell           = dataset$cell_data$cell,
    source         = "mock",
    reduction_used = red,
    root_field     = root_field,
    root_group     = as.character(root_group),
    metadata_field = NA_character_,
    n_lineages     = 1L,
    method_details = list(centroid = centroid)
  )
}

#' Bin a gene's expression along pseudotime for a smooth trend overlay.
#'
#' @param pt      pseudotime vector
#' @param expr    expression vector aligned with `pt`
#' @param n_bins  number of equal-width bins along the pseudotime axis
#'
#' @return data.frame with columns `bin`, `pt_mid`, `expr_mean`, `n`. Bins
#'   with no cells get `expr_mean = NA`.
bin_gene_by_pseudotime <- function(pt, expr, n_bins = 20L) {
  stopifnot(length(pt) == length(expr))
  ok <- !is.na(pt) & !is.na(expr)
  pt <- pt[ok]; expr <- expr[ok]
  if (length(pt) < 2L) return(NULL)
  n_bins <- max(2L, as.integer(n_bins))
  brk <- seq(min(pt), max(pt), length.out = n_bins + 1L)
  bin <- cut(pt, breaks = brk, include.lowest = TRUE, labels = FALSE)
  data.frame(
    bin       = seq_len(n_bins),
    pt_mid    = (brk[-1] + brk[-length(brk)]) / 2,
    expr_mean = vapply(seq_len(n_bins), function(b) {
      idx <- which(bin == b)
      if (length(idx) == 0L) NA_real_ else mean(expr[idx])
    }, FUN.VALUE = numeric(1)),
    n         = vapply(seq_len(n_bins), function(b) sum(bin == b),
                       FUN.VALUE = integer(1))
  )
}

#' Build the summary list shown by the Trajectory module.
pseudotime_summary <- function(pt, source, root_field = NA, root_group = NA,
                               metadata_field = NA, reduction_used = NA) {
  if (is.null(pt) || length(pt) == 0L) {
    return(list(n_cells = 0L, pt_min = NA_real_, pt_max = NA_real_,
                source = source, root_field = root_field,
                root_group = root_group, metadata_field = metadata_field,
                reduction_used = reduction_used))
  }
  list(
    n_cells        = length(pt),
    pt_min         = min(pt, na.rm = TRUE),
    pt_max         = max(pt, na.rm = TRUE),
    source         = source,
    root_field     = root_field,
    root_group     = root_group,
    metadata_field = metadata_field,
    reduction_used = reduction_used
  )
}

#' TRUE iff `state` has a completed trajectory result with a pseudotime vector.
has_trajectory_results <- function(state) {
  tr <- state$analysis_results$trajectory
  !is.null(tr) && identical(tr$status, "completed") &&
    !is.null(tr$results$pseudotime) && length(tr$results$pseudotime) > 0L
}

# ---- Trajectory result schema check + dataset application ------------------
#
# There is no formal `trajectory_result_v1` class (yet). The canonical
# payload is whatever `run_trajectory()` returns -- a list with at
# minimum a numeric `pseudotime` and a character `cell` of the same
# length. `is_trajectory_result()` is a structural check on that
# contract so call sites can validate inputs without depending on
# upstream Shiny state.

#' TRUE iff `x` looks like a canonical trajectory result payload.
#'
#' Validates the minimum invariants downstream consumers rely on:
#'   * numeric `pseudotime` of length n_cells
#'   * character `cell` of the same length (aligned with `pseudotime`)
#'   * character(1) `source`
#'
#' Other optional fields (`reduction_used`, `root_field`, `root_group`,
#' `metadata_field`, `n_lineages`, `method_details`, ...) are populated
#' by `run_trajectory()` itself; their absence here is tolerated.
#'
#' @param x  any R object
#' @return logical(1)
is_trajectory_result <- function(x) {
  is.list(x) &&
    !is.null(x$pseudotime) && is.numeric(x$pseudotime) &&
    !is.null(x$cell) && is.character(x$cell) &&
    length(x$pseudotime) == length(x$cell) &&
    !is.null(x$source) && is.character(x$source) && length(x$source) == 1L
}

#' Unwrap the `state$analysis_results$trajectory` slot to a bare
#' trajectory result payload, or return `x` as-is if it's already one.
#'
#' Accepts either the full slot (`list(status, results, params, ...)`)
#' or just the inner `results` list. Returns `NULL` and never throws
#' if the slot doesn't carry a usable result.
.unwrap_trajectory_result <- function(x) {
  if (is.null(x)) return(NULL)
  if (is_trajectory_result(x)) return(x)
  if (is.list(x) && !is.null(x$results) && is_trajectory_result(x$results)) {
    if (!is.null(x$status) && !identical(x$status, "completed")) return(NULL)
    return(x$results)
  }
  NULL
}

#' Bake a trajectory result into dataset metadata columns.
#'
#' Mirrors `apply_annotations_to_dataset()`: the trajectory result
#' itself remains the primary source of truth (stored under
#' `state$analysis_results$trajectory`), and this function exists for
#' export workflows, downstream tools that only know how to read
#' metadata columns, and for users who want to color by pseudotime in
#' the Explorer / use it as a grouping variable in DE.
#'
#' Always writes a numeric column
#'
#'   pseudotime__<source>__<YYYY_MM_DD>
#'
#' and, when `bins > 0`, also a categorical bin column
#'
#'   pseudotime_bin__<source>__<YYYY_MM_DD>
#'
#' Provenance is attached to each column via `attr()` so downstream
#' readers can trace which trajectory method produced the values.
#'
#' Validation:
#'   * `trajectory_result` must be a canonical payload (or a wrapped
#'     trajectory slot whose `status` is `"completed"`);
#'   * the result's `cell` vector must cover every dataset cell;
#'   * the destination columns must not already exist;
#'   * generic names like `pseudotime` / `pseudotime_bin` are refused.
#'
#' @param dataset             the dataset list
#' @param trajectory_result   canonical payload from `run_trajectory()`,
#'                            OR the wrapped slot
#'                            `state$analysis_results$trajectory`.
#' @param bins                integer >= 0. When > 0, also writes a
#'                            categorical bin column with `bins`
#'                            equal-width bins. Defaults to 0 (no bins).
#' @param applied_at          POSIXct used to date the column name.
#'                            Defaults to `Sys.time()`. Exposed so
#'                            tests can pin the date.
#'
#' @return a new dataset list (does not mutate in place).
apply_pseudotime_to_dataset <- function(dataset, trajectory_result,
                                        bins = 0L,
                                        applied_at = Sys.time()) {
  if (is.null(dataset)) stop("No dataset provided.", call. = FALSE)

  tr <- .unwrap_trajectory_result(trajectory_result)
  if (is.null(tr))
    stop("`trajectory_result` is not a completed trajectory result.",
         call. = FALSE)
  if (!is_trajectory_result(tr))
    stop("`trajectory_result` failed schema check: needs numeric ",
         "`pseudotime`, character `cell`, character(1) `source`.",
         call. = FALSE)

  ds_cells <- dataset$cell_data$cell
  if (is.null(ds_cells))
    stop("Dataset has no `cell_data$cell` column.", call. = FALSE)

  pos <- match(ds_cells, tr$cell)
  if (any(is.na(pos))) {
    n_missing <- sum(is.na(pos))
    stop(sprintf(
      "Trajectory result covers %d/%d cells; %d dataset cells have no value.",
      sum(!is.na(pos)), length(ds_cells), n_missing), call. = FALSE)
  }
  pt_aligned <- as.numeric(tr$pseudotime[pos])

  bins <- as.integer(bins %||% 0L)
  if (is.na(bins) || bins < 0L)
    stop("`bins` must be a non-negative integer.", call. = FALSE)

  date_str <- format(applied_at, "%Y_%m_%d")
  src      <- as.character(tr$source)
  pt_col   <- sprintf("pseudotime__%s__%s",     src, date_str)
  bin_col  <- sprintf("pseudotime_bin__%s__%s", src, date_str)

  if (identical(pt_col, "pseudotime") || identical(bin_col, "pseudotime_bin"))
    stop("Refusing to write a generic 'pseudotime' / 'pseudotime_bin' column.",
         call. = FALSE)
  if (pt_col %in% names(dataset$cell_data))
    stop(sprintf("Column '%s' already exists; refusing to overwrite.",
                 pt_col), call. = FALSE)
  if (bins > 0L && bin_col %in% names(dataset$cell_data))
    stop(sprintf("Column '%s' already exists; refusing to overwrite.",
                 bin_col), call. = FALSE)

  dataset$cell_data[[pt_col]] <- .stamp_pseudotime_attrs(
    pt_aligned, tr, applied_at = applied_at, kind = "numeric")

  if (bins > 0L) {
    ok <- !is.na(pt_aligned)
    bin_int <- rep(NA_integer_, length(pt_aligned))
    if (any(ok)) {
      rng <- range(pt_aligned[ok])
      if (is.finite(diff(rng)) && diff(rng) > 0) {
        brk <- seq(rng[1], rng[2], length.out = bins + 1L)
        bin_int[ok] <- as.integer(cut(pt_aligned[ok], breaks = brk,
                                      include.lowest = TRUE,
                                      labels = FALSE))
      } else {
        bin_int[ok] <- 1L  # degenerate (constant pseudotime) -> single bin
      }
    }
    bin_lbl <- ifelse(is.na(bin_int),
                      NA_character_,
                      sprintf("bin_%02d", bin_int))
    bin_vec <- factor(bin_lbl, levels = sprintf("bin_%02d", seq_len(bins)))
    dataset$cell_data[[bin_col]] <- .stamp_pseudotime_attrs(
      bin_vec, tr, applied_at = applied_at, kind = "bin", bins = bins)
  }

  new_cols <- if (bins > 0L) c(pt_col, bin_col) else pt_col
  dataset$metadata_fields <- unique(c(dataset$metadata_fields, new_cols))
  dataset
}

# Internal: attach provenance attributes to a pseudotime column. Returns
# the (possibly factor) vector so it can be reassigned in one step.
.stamp_pseudotime_attrs <- function(col, tr, applied_at, kind, bins = NA_integer_) {
  attr(col, "pseudotime_source")    <- as.character(tr$source)
  attr(col, "reduction_used")       <- as.character(tr$reduction_used %||% NA_character_)
  attr(col, "root_field")           <- as.character(tr$root_field %||% NA_character_)
  attr(col, "root_group")           <- as.character(tr$root_group %||% NA_character_)
  attr(col, "metadata_field")       <- as.character(tr$metadata_field %||% NA_character_)
  attr(col, "n_lineages")           <- as.integer(tr$n_lineages %||% NA_integer_)
  attr(col, "applied_at")           <- applied_at
  attr(col, "kind")                 <- kind
  if (!is.na(bins)) attr(col, "bins") <- as.integer(bins)
  col
}

# ---- Internals ------------------------------------------------------------

rescale01 <- function(v) {
  if (length(v) == 0L) return(numeric())
  rng <- range(v, na.rm = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) == 0) {
    out <- rep(0.5, length(v))
    out[is.na(v)] <- NA_real_
    return(out)
  }
  (v - rng[1]) / diff(rng)
}
