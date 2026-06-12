# ============================================================================
# Regulon source registry
# ----------------------------------------------------------------------------
# Mirrors PATHWAY_SOURCES. A regulon source is "what produces a
# regulon_set()". Sources are the *catalogue* layer:
#
#   mock_pbmc        -- 4 hand-built regulons whose target genes are
#                       the MOCK_GENES used by the mock dataset. No
#                       deps. Always available so CI / demos work.
#   dorothea_human   -- DoRothEA curated TF-target collections (human).
#                       Confidence levels A-E. Routed via the
#                       `dorothea` Bioconductor package.
#   dorothea_mouse   -- DoRothEA mouse equivalent.
#   scenic_grn       -- (seam) a regulon_set returned by a SCENIC run
#                       (GENIE3 + cisTarget). Not registered until the
#                       SCENIC engine lands.
#
# Each source returns a `regulon_set()`; engines (AUCell, etc.) consume
# regulon_sets. The split mirrors how pathway sources and pathway
# methods are kept independent.
# ============================================================================

#' Specification for one regulon source.
#'
#' @param id          snake_case identifier
#' @param name        display name for the UI picker
#' @param species     "human" | "mouse" | "other"
#' @param requires    character() optional packages required to call
#'                    `fetch_fn()`. Empty = always available.
#' @param parameters  named list of parameter specs (e.g. confidence
#'                    levels for DoRothEA).
#' @param fetch_fn    function(params) -> regulon_set()
#' @param description character(1)
regulon_source_spec <- function(id, name,
                                species = "human",
                                requires    = character(),
                                parameters  = list(),
                                fetch_fn,
                                description = "") {
  stopifnot(is.character(id), length(id) == 1L, nzchar(id))
  stopifnot(is.function(fetch_fn))
  list(
    id          = id,
    name        = name,
    species     = species,
    requires    = as.character(requires),
    parameters  = parameters,
    fetch_fn    = fetch_fn,
    description = description
  )
}

#' Source registry.
REGULON_SOURCES <- function() {
  list(
    regulon_source_spec(
      id         = "mock_pbmc",
      name       = "Mock PBMC (built-in)",
      species    = "human",
      requires   = character(),
      parameters = list(),
      fetch_fn   = .mock_pbmc_regulon_set,
      description = paste(
        "Four hand-built regulons (GATA3, PAX5, SPI1, KLF5) whose",
        "target genes are the canonical markers used by",
        "`mock_dataset()`. Always available. Useful for CI, demos,",
        "and validating the AUCell engine on the mock dataset.")
    ),
    regulon_source_spec(
      id          = "dorothea_human",
      name        = "DoRothEA (human, AB)",
      species     = "human",
      requires    = c("dorothea"),
      parameters  = list(
        confidence = list(type = "select", required = FALSE,
                          default = "AB",
                          choices = c("A", "AB", "ABC", "ABCD", "ABCDE"),
                          description = paste(
                            "DoRothEA confidence prefix:",
                            "A = highest confidence (literature-curated +",
                            "ChIP-seq + motif).  AB = +- 700 TFs.",
                            "Adding C/D/E broadens but adds noise."))),
      fetch_fn    = function(params) .fetch_dorothea_regulons("human", params),
      description = paste(
        "DoRothEA curated TF-target regulons (Garcia-Alonso et al,",
        "2019). Includes ~1400 human TFs with confidence levels A-E.")
    ),
    regulon_source_spec(
      id          = "dorothea_mouse",
      name        = "DoRothEA (mouse, AB)",
      species     = "mouse",
      requires    = c("dorothea"),
      parameters  = list(
        confidence = list(type = "select", required = FALSE,
                          default = "AB",
                          choices = c("A", "AB", "ABC", "ABCD", "ABCDE"))),
      fetch_fn    = function(params) .fetch_dorothea_regulons("mouse", params),
      description = "DoRothEA mouse TF-target regulons."
    )
  )
}

#' Look up a source by id.
get_regulon_source <- function(id) {
  for (s in REGULON_SOURCES()) if (identical(s$id, id)) return(s)
  NULL
}

#' Available regulon sources as named character: display_name -> id.
#' Unavailable ones are kept but labelled.
list_regulon_sources <- function() {
  srcs <- REGULON_SOURCES()
  labels <- vapply(srcs, function(s) {
    if (has_optional(s$requires)) s$name
    else sprintf("%s  (not installed)", s$name)
  }, character(1))
  setNames(vapply(srcs, `[[`, character(1), "id"), labels)
}

#' Build a regulon_set from a registered source.
#'
#' @param id      source id
#' @param params  parameter list for the source's `fetch_fn`
#' @return regulon_set
fetch_regulon_set <- function(id, params = list()) {
  s <- get_regulon_source(id)
  if (is.null(s))
    stop(sprintf("Unknown regulon source '%s'. Have: %s",
                 id, paste(unname(list_regulon_sources()), collapse = ", ")),
         call. = FALSE)
  out <- s$fetch_fn(params)
  if (!is_regulon_set(out))
    stop(sprintf(
      "Regulon source '%s' did not return a regulon_set().", id),
      call. = FALSE)
  out
}

# ============================================================================
# Built-in source implementations
# ----------------------------------------------------------------------------
# Pure functions: no shiny, no state mutation.
# ============================================================================

# Hand-built mock regulon set. TFs themselves are not in MOCK_GENES,
# which is fine -- AUCell scoring only ranks the targets. The cluster
# bias of MOCK_GENES means each regulon should peak in exactly one of
# the four mock clusters.
.mock_pbmc_regulon_set <- function(params = list()) {
  regulon_set(
    id      = "mock_pbmc",
    name    = "Mock PBMC regulons",
    species = "human",
    source  = "builtin",
    version = "0.1.0",
    regulons = list(
      regulon_spec(tf = "GATA3", targets = c("CD3D", "NKG7"),
                   type = "activating"),
      regulon_spec(tf = "PAX5",  targets = c("MS4A1"),
                   type = "activating"),
      regulon_spec(tf = "SPI1",  targets = c("LST1"),
                   type = "activating"),
      regulon_spec(tf = "KLF5",  targets = c("EPCAM", "COL1A1"),
                   type = "activating")
    )
  )
}

# DoRothEA bridge. The package ships `dorothea_hs` / `dorothea_mm`
# data.frames with columns tf, confidence, target, mor. We filter on
# the confidence prefix and turn rows into regulon_spec()s.
.fetch_dorothea_regulons <- function(species, params = list()) {
  require_optional("dorothea",
                   feature = "DoRothEA regulon source",
                   source  = "Bioconductor")
  conf <- params$confidence %||% "AB"
  conf_levels <- strsplit(conf, "")[[1]]

  data_name <- switch(species,
    "human" = "dorothea_hs",
    "mouse" = "dorothea_mm",
    stop("DoRothEA species must be 'human' or 'mouse'.", call. = FALSE))

  env <- new.env()
  utils::data(list = data_name, package = "dorothea", envir = env)
  df  <- get(data_name, envir = env)
  df  <- df[df$confidence %in% conf_levels, , drop = FALSE]
  if (nrow(df) == 0L)
    stop(sprintf("DoRothEA: no regulons at confidence '%s'.", conf),
         call. = FALSE)

  regulons <- .dorothea_df_to_regulons(df)
  regulon_set(
    id       = sprintf("dorothea_%s_%s",
                       substr(species, 1, 2), conf),
    name     = sprintf("DoRothEA %s (%s)", species, conf),
    species  = species,
    source   = "dorothea",
    version  = sprintf("dorothea/%s", utils::packageVersion("dorothea")),
    regulons = regulons
  )
}

# Pure converter: DoRothEA-shaped data.frame -> list of regulon_spec.
# Factored out so we can test it without DoRothEA installed.
.dorothea_df_to_regulons <- function(df) {
  needed <- c("tf", "target")
  missing <- setdiff(needed, names(df))
  if (length(missing))
    stop(".dorothea_df_to_regulons: missing column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  tfs <- split(df, df$tf)
  out <- lapply(names(tfs), function(tf) {
    rows <- tfs[[tf]]
    # MoR (mode of regulation): +1 activating, -1 repressing, 0 unknown.
    mor <- if ("mor" %in% names(rows)) rows$mor else rep(1, nrow(rows))
    type <- if (all(mor >= 0)) "activating"
            else if (all(mor <= 0)) "repressing"
            else "unknown"
    regulon_spec(
      tf      = tf,
      targets = as.character(rows$target),
      weights = abs(as.numeric(mor)),
      type    = type
    )
  })
  out
}
