# ============================================================================
# Pathway Analysis
# ----------------------------------------------------------------------------
# Two pieces:
#   1. A pluggable source layer (R/pathway_sources.R) that turns a
#      collection id into a named list of `pathway -> character(genes)`.
#      Today: `builtin` (mock_v1) + `msigdbr` (when installed). The
#      Pathway module never sees the source -- it only calls
#      `available_pathway_collections()` and `get_pathways(id)`.
#   2. Enrichment statistics in this file:
#        - `compute_enrichment()` -- Fisher's exact ORA with BH adjustment.
#        - `compute_gsea()`       -- ranked GSEA scaffold (errors until
#                                    `fgsea` is installed).
#
# Conceptually distinct from DE: DE produces a per-gene table; pathway
# enrichment is *set-based* and operates on a chosen significant gene
# list (ORA) or a ranked gene vector (GSEA).
#
# Universe handling (ORA): the universe is the population from which
# `selected` was drawn. Calling code should pass:
#   * `universe = de_tested_genes` -- the set of genes actually tested
#     by DE. Most defensible for small / mock data.
#   * `universe = available_genes(dataset)` -- the full filtered feature
#     set. Sensible for production-scale DE.
#   * `universe = NULL` -- fallback `union(selected, unlist(pathways))`,
#     which is liberal (small denominators inflate odds ratios). Used
#     here only when no better universe is supplied.
#
# Sources / Phase-2 routes are documented at the top of
# R/pathway_sources.R.
# ============================================================================

# ---- Mock gene-set library -------------------------------------------------
# Each collection is itself a named list of pathway -> character() of gene
# symbols. The genes were chosen to (a) cover the mock dataset's 6 demo
# genes so the UI shows real overlap, and (b) include enough additional
# symbols (~12-18 per pathway) that Fisher's contingency tables are not
# degenerate.

BUILTIN_PATHWAYS <- list(
  mock_v1 = list(
    "T cell activation" = c(
      "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD4",
      "ZAP70", "LCK", "ITK", "LAT", "TRAC", "TRBC1",
      "IL2RA", "IL7R", "NKG7", "CD2"),
    "B cell receptor signaling" = c(
      "MS4A1", "CD19", "CD79A", "CD79B", "BTK", "BLNK",
      "PLCG2", "SYK", "PAX5", "POU2AF1", "TNFRSF13B",
      "TNFRSF13C", "CR2", "VPREB3", "IGHM"),
    "Myeloid inflammatory response" = c(
      "LST1", "CD14", "CD68", "FCGR3A", "FCGR1A",
      "S100A8", "S100A9", "S100A12", "IL1B", "TNF",
      "CXCL8", "CCR2", "CSF1R", "MARCO", "CD163"),
    "Epithelial program" = c(
      "EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT14",
      "CDH1", "CLDN3", "CLDN4", "OCLN", "TJP1", "MUC1",
      "ELF3", "GRHL2", "ESRP1"),
    "Extracellular matrix organization" = c(
      "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1",
      "FN1", "MMP2", "MMP9", "MMP14", "TIMP1", "LOX",
      "BGN", "DCN", "POSTN", "SPARC", "FBN1"),
    "Cytotoxicity" = c(
      "NKG7", "GZMB", "GZMA", "GZMK", "GZMH", "PRF1",
      "GNLY", "KLRD1", "KLRF1", "FGFBP2", "FCGR3A",
      "CTSW", "CST7")
  )
)

#' Names of the available pathway collections.
#'
#' Aggregates collections from every *available* registered source
#' (see `PATHWAY_SOURCES()`). Built-in collections are listed first so
#' `available_pathway_collections()[1]` remains a stable demo default.
#' MSigDB collections appear only when the `msigdbr` package is
#' installed (else they're silently filtered).
available_pathway_collections <- function() {
  out <- character()
  for (s in PATHWAY_SOURCES()) {
    if (length(s$requires) > 0L && !has_optional(s$requires)) next
    out <- c(out, tryCatch(s$collections_fn(), error = function(e) character()))
  }
  unique(out)
}

#' Fetch the gene-set list for a collection.
#'
#' Resolves the owning source from the registry. Returns NULL on a blank
#' collection id (so reactive UI code can pass `input$collection` before
#' a selection is made), but raises an informative error for non-empty
#' ids that no registered source owns -- the silent-NULL path used to
#' show up as "0 pathways" downstream and was a recurring source of
#' confusion in real-data smoke testing. Also raises a clear
#' "install X" error when the source is registered but its package
#' isn't installed.
get_pathways <- function(collection) {
  if (is.null(collection) || !nzchar(collection)) return(NULL)
  src <- .resolve_collection_source(collection)
  if (is.null(src)) {
    avail <- available_pathway_collections()
    stop(sprintf(
      "Unknown pathway collection '%s'. Available: %s.",
      collection,
      if (length(avail)) paste(avail, collapse = ", ") else "<none>"),
      call. = FALSE)
  }
  src$fetcher(collection)
}

# ---- Gene selection from DE results ----------------------------------------

#' Pick selected genes from a DE result given a direction + thresholds.
#'
#' @param de              data.frame produced by `compute_de()` (must contain
#'                        `gene`, `avg_log2FC`, `p_val_adj`)
#' @param direction       one of "up_in_g1" (avg_log2FC > 0),
#'                        "up_in_g2" (avg_log2FC < 0), or "both"
#' @param padj_cutoff     keep genes with `p_val_adj <= padj_cutoff`
#' @param log2fc_cutoff   keep genes with `|avg_log2FC| >= log2fc_cutoff`
#'
#' @return character() of unique gene symbols
select_de_genes <- function(de, direction = c("up_in_g1", "up_in_g2", "both"),
                            padj_cutoff = 0.05, log2fc_cutoff = 0.5) {
  direction <- match.arg(direction)
  if (is.null(de) || nrow(de) == 0L) return(character())
  ok <- !is.na(de$p_val_adj) & de$p_val_adj <= padj_cutoff &
        abs(de$avg_log2FC) >= log2fc_cutoff
  if (direction == "up_in_g1") ok <- ok & de$avg_log2FC > 0
  if (direction == "up_in_g2") ok <- ok & de$avg_log2FC < 0
  unique(de$gene[ok])
}

# ---- ORA via Fisher's exact -----------------------------------------------

#' Overrepresentation analysis on a gene list against a gene-set library.
#'
#' For each pathway P, builds the 2x2 contingency table
#'
#'     selected & in P     | selected & not in P
#'     not selected & in P | not selected & not in P
#'
#' and runs `stats::fisher.test(..., alternative = "greater")`.
#'
#' @param selected    character(); the chosen gene list
#' @param pathways    named list of character(); the gene-set library
#' @param universe    character(); the background. Defaults to
#'                    `union(selected, unlist(pathways))`. Real workflows
#'                    should pass the set of DE-tested genes (or whole
#'                    transcriptome) instead.
#' @param direction   character(1); copied into every output row so callers
#'                    can rbind multi-direction runs
#' @param collection  character(1); copied into every output row
#'
#' @return data.frame with columns:
#'   pathway, collection, direction, n_genes_in_pathway, n_overlap,
#'   overlap_genes (semicolon-joined), odds_ratio, p_val, p_val_adj
compute_enrichment <- function(selected, pathways,
                               universe = NULL,
                               direction = NA_character_,
                               collection = NA_character_) {
  selected <- unique(as.character(selected))
  if (is.null(pathways) || length(pathways) == 0L) return(empty_pathway_results())
  if (is.null(universe))
    universe <- unique(c(selected, unlist(pathways, use.names = FALSE)))
  universe <- unique(as.character(universe))
  selected <- intersect(selected, universe)

  rows <- list()
  for (pname in names(pathways)) {
    pgenes <- intersect(unique(pathways[[pname]]), universe)
    overlap <- intersect(selected, pgenes)
    a <- length(overlap)
    b <- length(setdiff(selected, pgenes))                  # selected, not in P
    c <- length(setdiff(pgenes,   selected))                # in P, not selected
    d <- length(universe) - a - b - c                       # neither
    if (a + c == 0L || a + b == 0L) {
      p <- 1; or <- NA_real_
    } else {
      m <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
      ft <- tryCatch(stats::fisher.test(m, alternative = "greater"),
                     error = function(e) NULL)
      p  <- if (is.null(ft)) NA_real_ else ft$p.value
      or <- if (is.null(ft)) NA_real_ else as.numeric(ft$estimate)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      pathway            = pname,
      collection         = as.character(collection),
      direction          = as.character(direction),
      n_genes_in_pathway = length(pgenes),
      n_overlap          = a,
      overlap_genes      = if (a > 0L) paste(overlap, collapse = ";") else "",
      odds_ratio         = or,
      p_val              = p,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) return(empty_pathway_results())
  df <- do.call(rbind, rows)
  df$p_val_adj <- stats::p.adjust(df$p_val, method = "BH")
  df <- df[order(df$p_val_adj, df$p_val, -df$n_overlap), , drop = FALSE]
  rownames(df) <- NULL
  df
}

empty_pathway_results <- function() {
  data.frame(
    pathway = character(), collection = character(), direction = character(),
    n_genes_in_pathway = integer(), n_overlap = integer(),
    overlap_genes = character(), odds_ratio = numeric(),
    p_val = numeric(), p_val_adj = numeric(),
    stringsAsFactors = FALSE
  )
}

# ---- Ranked GSEA scaffold (fgsea) ------------------------------------------
#
# `compute_gsea()` is the design seam for ranked-list pathway analysis
# (e.g. fgsea). It mirrors the shape of `compute_enrichment()` so callers
# can swap ORA / GSEA without restructuring downstream code. Today the
# only implementation routes through the `fgsea` package; if `fgsea` is
# not installed the call errors with a clear install hint via
# `require_optional()`. A pure converter `.fgsea_to_gsea_schema()` keeps
# the schema mapping testable without fgsea.

#' Empty GSEA result frame (canonical schema).
empty_gsea_results <- function() {
  data.frame(
    pathway            = character(),
    collection         = character(),
    n_genes_in_pathway = integer(),
    n_leading_edge     = integer(),
    leading_edge_genes = character(),
    ES                 = numeric(),
    NES                = numeric(),
    p_val              = numeric(),
    p_val_adj          = numeric(),
    stringsAsFactors   = FALSE
  )
}

#' Ranked-list GSEA via fgsea (scaffold).
#'
#' @param ranked_genes  named numeric vector. Names are gene symbols,
#'                      values are the per-gene ranking metric (e.g.
#'                      signed log2FC or -log10(p) * sign(log2FC)).
#' @param pathways      named list of character() gene sets, as returned
#'                      by `get_pathways(collection)`.
#' @param collection    character(1), stamped on output rows.
#' @param min_size,max_size  fgsea pathway-size filters.
#' @param ...           passed through to `fgsea::fgsea()`.
#'
#' @return data.frame with the schema of `empty_gsea_results()`.
compute_gsea <- function(ranked_genes, pathways,
                         collection = NA_character_,
                         min_size = 5L, max_size = 500L, ...) {
  if (length(ranked_genes) == 0L) return(empty_gsea_results())
  if (is.null(names(ranked_genes)))
    stop("compute_gsea: `ranked_genes` must be a NAMED numeric vector.",
         call. = FALSE)
  if (is.null(pathways) || length(pathways) == 0L) return(empty_gsea_results())

  require_optional("fgsea",
                   feature = "ranked GSEA via fgsea",
                   source  = "Bioconductor")

  fgsea_fn <- get("fgsea", envir = asNamespace("fgsea"))
  res <- fgsea_fn(pathways = pathways, stats = ranked_genes,
                  minSize = min_size, maxSize = max_size, ...)
  .fgsea_to_gsea_schema(res, pathways = pathways, collection = collection)
}

#' Pure converter: fgsea-style result -> canonical GSEA schema.
#' Accepts a data.frame / tibble / data.table with at least:
#'   pathway, ES, NES, pval, padj, size  (and optionally leadingEdge,
#'   a list-column of character vectors).
#' Factored out so the schema can be regression-tested without fgsea.
.fgsea_to_gsea_schema <- function(res, pathways,
                                  collection = NA_character_) {
  if (is.null(res) || NROW(res) == 0L) return(empty_gsea_results())
  d <- as.data.frame(res, stringsAsFactors = FALSE)
  need <- c("pathway", "ES", "NES", "pval", "padj")
  miss <- setdiff(need, names(d))
  if (length(miss) > 0L)
    stop(".fgsea_to_gsea_schema: missing columns: ",
         paste(miss, collapse = ", "), call. = FALSE)

  le_col <- if ("leadingEdge" %in% names(d)) d$leadingEdge else
            vector("list", nrow(d))
  leading_edge_genes <- vapply(le_col, function(le) {
    if (is.null(le) || length(le) == 0L) "" else
      paste(as.character(le), collapse = ";")
  }, character(1))
  n_leading_edge <- vapply(le_col, function(le)
    if (is.null(le)) 0L else length(le), integer(1))

  n_in_pw <- vapply(as.character(d$pathway), function(p)
    if (!is.null(pathways[[p]])) length(unique(pathways[[p]])) else NA_integer_,
    integer(1))

  out <- data.frame(
    pathway            = as.character(d$pathway),
    collection         = as.character(collection),
    n_genes_in_pathway = as.integer(n_in_pw),
    n_leading_edge     = as.integer(n_leading_edge),
    leading_edge_genes = leading_edge_genes,
    ES                 = as.numeric(d$ES),
    NES                = as.numeric(d$NES),
    p_val              = as.numeric(d$pval),
    p_val_adj          = as.numeric(d$padj),
    stringsAsFactors   = FALSE
  )
  out <- out[order(out$p_val_adj, out$p_val), , drop = FALSE]
  rownames(out) <- NULL
  out
}
