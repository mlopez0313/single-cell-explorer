# ============================================================================
# Pathway source registry
# ----------------------------------------------------------------------------
# Mirrors the annotation-engine + DE-backend registry pattern: the
# enrichment math (Fisher's-exact ORA in R/pathway.R; future GSEA in
# `compute_gsea()`) is fully decoupled from where the gene sets come
# from. A "source" turns a collection id into a named list of
# `pathway_name -> character(gene_symbols)`.
#
# Today's sources:
#   * `builtin`  -- the hand-curated `mock_v1` collection used for demo
#                   mode and tests. No dependencies.
#   * `msigdbr`  -- MSigDB (human / mouse) via the CRAN package `msigdbr`.
#                   The collections are pinned (H, C2:CP:REACTOME,
#                   C2:CP:KEGG_LEGACY, C5:GO:BP, C7) so the available
#                   collection list is fast and deterministic.
#
# Design:
#   * `pathway_source_spec(id, label, requires, fetcher, collections_fn)`
#       declares one source.
#   * `PATHWAY_SOURCES()` returns every spec (filtered by availability
#       on the way out to `available_pathway_collections()`).
#   * Collection ids are globally unique strings; `builtin` advertises
#       `"mock_v1"` (no prefix, for back-compat), `msigdbr` advertises
#       `"msigdbr/<gs_cat>[:<gs_subcat>]"`.
#   * `.resolve_collection_source(id)` walks the registry and returns
#       the spec that owns the id (NULL on miss).
#
# Pure helpers (no heavy deps):
#   * `.msigdbr_to_pathways(df)` -- converts a msigdbr-shaped tibble or
#                                    data.frame into a named list. Lets
#                                    the schema mapping be regression-
#                                    tested without `msigdbr` installed.
#
# Future plug-in slots (each will land as a new source):
#   * reactome.db (Bioc)
#   * KEGGREST    (Bioc, REST-backed lazy fetch)
#   * org.Hs.eg.db / GO.db (Bioc)
#   * user-supplied .gmt files
# ============================================================================

# Process-local cache for msigdbr lookups -- repeated calls against the
# same (category, subcategory, species) triple are free after the first.
.msigdbr_cache <- new.env(parent = emptyenv())

#' Specification for one pathway source.
#'
#' @param id              source identifier ("builtin", "msigdbr", ...)
#' @param label           display label for the UI
#' @param requires        character() optional R package names
#' @param fetcher         function(collection_id) -> named list
#'                        (pathway_name -> gene symbols). Receives the
#'                        bare collection id (without source prefix).
#' @param collections_fn  function() -> character() of advertised
#'                        collection ids (FULLY-QUALIFIED, including any
#'                        prefix). Called only when the source is
#'                        available; otherwise the source returns an
#'                        empty advertised set.
#' @param description     character(1)
pathway_source_spec <- function(id, label,
                                requires = character(),
                                fetcher,
                                collections_fn,
                                description = "") {
  stopifnot(is.character(id), length(id) == 1L, nzchar(id))
  stopifnot(is.function(fetcher), is.function(collections_fn))
  list(
    id              = id,
    label           = label,
    requires        = requires,
    fetcher         = fetcher,
    collections_fn  = collections_fn,
    description     = description
  )
}

#' The full registry. Function so it re-evaluates per call (cheap), which
#' makes adding a source in an interactive session immediate.
PATHWAY_SOURCES <- function() {
  list(
    pathway_source_spec(
      id             = "builtin",
      label          = "Built-in demo",
      requires       = character(),
      fetcher        = .fetch_builtin_pathways,
      collections_fn = function() names(BUILTIN_PATHWAYS),
      description    = paste("Hand-curated demo collection (mock_v1) covering",
                             "the genes used by the synthetic dataset.")
    ),
    pathway_source_spec(
      id             = "msigdbr",
      label          = "MSigDB (via msigdbr)",
      requires       = "msigdbr",
      fetcher        = .fetch_msigdbr_pathways,
      collections_fn = .msigdbr_collections,
      description    = paste(
        "MSigDB gene sets for human / mouse. Pinned subcategories:",
        "Hallmark (H), Reactome (C2:CP:REACTOME),",
        "KEGG legacy (C2:CP:KEGG_LEGACY), GO BP (C5:GO:BP),",
        "Immunologic (C7).")
    )
  )
}

#' Look up a source by id. NULL if not registered.
get_pathway_source <- function(id) {
  for (s in PATHWAY_SOURCES()) if (identical(s$id, id)) return(s)
  NULL
}

#' Find which source owns a collection id. NULL if no source claims it.
.resolve_collection_source <- function(collection_id) {
  for (s in PATHWAY_SOURCES()) {
    if (length(s$requires) > 0L && !has_optional(s$requires)) next
    collections <- tryCatch(s$collections_fn(), error = function(e) character())
    if (collection_id %in% collections) return(s)
  }
  NULL
}

#' Metadata for one collection. Empty list if not found.
pathway_collection_info <- function(collection_id) {
  src <- .resolve_collection_source(collection_id)
  if (is.null(src)) return(list())
  list(
    id           = collection_id,
    source_id    = src$id,
    source_label = src$label,
    available    = TRUE
  )
}

# ---- Built-in source -----------------------------------------------------

.fetch_builtin_pathways <- function(collection_id) {
  BUILTIN_PATHWAYS[[collection_id]]
}

# ---- msigdbr source ------------------------------------------------------

# Pinned (cat, subcat, label) triples. Keep this list short and curated;
# users who need other categories can still query `msigdbr()` directly.
.MSIGDBR_PINNED <- list(
  list(id = "msigdbr/H",                 cat = "H",  subcat = NA_character_,
       label = "MSigDB Hallmark"),
  list(id = "msigdbr/C2:CP:REACTOME",    cat = "C2", subcat = "CP:REACTOME",
       label = "Reactome (C2 canonical pathways)"),
  list(id = "msigdbr/C2:CP:KEGG_LEGACY", cat = "C2", subcat = "CP:KEGG_LEGACY",
       label = "KEGG legacy (C2 canonical pathways)"),
  list(id = "msigdbr/C5:GO:BP",          cat = "C5", subcat = "GO:BP",
       label = "GO Biological Process (C5)"),
  list(id = "msigdbr/C7",                cat = "C7", subcat = NA_character_,
       label = "Immunologic signatures (C7)")
)

.msigdbr_collections <- function() {
  vapply(.MSIGDBR_PINNED, `[[`, character(1), "id")
}

.fetch_msigdbr_pathways <- function(collection_id, species = "Homo sapiens") {
  spec <- NULL
  for (s in .MSIGDBR_PINNED) {
    if (identical(s$id, collection_id)) { spec <- s; break }
  }
  if (is.null(spec)) {
    stop("MSigDB source: unknown collection '", collection_id,
         "'. Known: ", paste(.msigdbr_collections(), collapse = ", "),
         ".", call. = FALSE)
  }
  require_optional("msigdbr",
                   feature = sprintf("MSigDB collection '%s'", collection_id),
                   source  = "CRAN")
  df <- .msigdbr_fetch(spec$cat, spec$subcat, species)
  .msigdbr_to_pathways(df)
}

# Cached fetcher. Separate from .msigdbr_to_pathways so the converter can
# be tested with a hand-built data.frame without the package installed.
.msigdbr_fetch <- function(category, subcategory, species) {
  key <- paste(species, category, subcategory, sep = "|")
  if (!is.null(.msigdbr_cache[[key]])) return(.msigdbr_cache[[key]])
  call_args <- list(species = species)
  # msigdbr's API renamed `category`/`subcategory` to `collection`/`subcollection`
  # in v8. Try both for forward + backward compat.
  fn <- get("msigdbr", envir = asNamespace("msigdbr"))
  call_args$collection    <- category
  call_args$subcollection <- if (is.na(subcategory)) NULL else subcategory
  df <- tryCatch(do.call(fn, call_args),
                 error = function(e1) {
                   call_args$collection    <- NULL
                   call_args$subcollection <- NULL
                   call_args$category      <- category
                   call_args$subcategory   <- if (is.na(subcategory)) NULL else subcategory
                   do.call(fn, call_args)
                 })
  .msigdbr_cache[[key]] <- df
  df
}

#' Convert a `msigdbr::msigdbr()`-shaped frame to a named list of pathways.
#'
#' Pure: accepts any data.frame with a `gs_name` column and one of
#' (`gene_symbol`, `human_gene_symbol`, `gene_symbol_human`) for genes.
#' Pathways are returned in the same order as their first appearance in
#' the input frame; gene order within a pathway is preserved (de-duped).
#'
#' Factored out so the schema mapping is testable without the `msigdbr`
#' or `msigdbdf` packages installed.
.msigdbr_to_pathways <- function(df) {
  if (is.null(df) || NROW(df) == 0L) return(list())
  d <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!"gs_name" %in% names(d)) {
    stop(".msigdbr_to_pathways: input lacks `gs_name` column.",
         call. = FALSE)
  }
  gene_col <- intersect(c("gene_symbol", "human_gene_symbol",
                          "gene_symbol_human"), names(d))[1]
  if (is.na(gene_col)) {
    stop(".msigdbr_to_pathways: input lacks a gene-symbol column. ",
         "Need one of: gene_symbol, human_gene_symbol, gene_symbol_human.",
         call. = FALSE)
  }
  names_in_order <- unique(as.character(d$gs_name))
  out <- vector("list", length(names_in_order))
  names(out) <- names_in_order
  for (g in names_in_order) {
    syms <- unique(as.character(d[d$gs_name == g, gene_col]))
    syms <- syms[nzchar(syms) & !is.na(syms)]
    out[[g]] <- syms
  }
  out
}
