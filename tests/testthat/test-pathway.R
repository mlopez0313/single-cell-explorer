test_that("BUILTIN_PATHWAYS exposes the mock_v1 collection", {
  cols <- available_pathway_collections()
  expect_true("mock_v1" %in% cols)
  pw <- get_pathways("mock_v1")
  expect_type(pw, "list")
  expect_true(length(pw) > 0L)
  expect_true("T cell activation" %in% names(pw))
})

test_that("select_de_genes filters by direction + thresholds", {
  ds <- mock_dataset(n_cells = 150, seed = 12)
  de <- compute_de(ds, grouping_field = "cluster",
                   group_1 = "0", group_2 = "1",
                   assay = "RNA", layer = "data",
                   min_pct = 0, test = "wilcox")

  up1   <- select_de_genes(de, direction = "up_in_g1",
                           padj_cutoff = 1, log2fc_cutoff = 0)
  up2   <- select_de_genes(de, direction = "up_in_g2",
                           padj_cutoff = 1, log2fc_cutoff = 0)
  both  <- select_de_genes(de, direction = "both",
                           padj_cutoff = 1, log2fc_cutoff = 0)
  expect_true(length(both) >= length(up1))
  expect_true(length(both) >= length(up2))
  expect_identical(select_de_genes(NULL), character())
})

test_that("compute_enrichment returns the documented schema and BH-adjusts", {
  ds <- mock_dataset(n_cells = 200, seed = 14)
  de <- compute_de(ds, grouping_field = "cluster",
                   group_1 = "0", group_2 = "1",
                   assay = "RNA", layer = "data",
                   min_pct = 0, test = "wilcox")
  selected <- select_de_genes(de, direction = "both",
                              padj_cutoff = 1, log2fc_cutoff = 0)
  pw  <- get_pathways("mock_v1")
  out <- compute_enrichment(selected, pw, universe = ds$genes,
                            direction = "both", collection = "mock_v1")
  expect_s3_class(out, "data.frame")
  expect_true(all(c("pathway", "collection", "direction", "n_genes_in_pathway",
                    "n_overlap", "overlap_genes", "odds_ratio", "p_val",
                    "p_val_adj")
                  %in% names(out)))
  # BH-adjusted p-values are >= raw p-values
  expect_true(all(out$p_val_adj >= out$p_val - 1e-12))
  # Ordered by adjusted p-value ascending
  expect_true(all(diff(out$p_val_adj) >= -1e-9))
})

test_that("compute_enrichment with no selected genes -> all pathways scored, zero overlap", {
  pw <- get_pathways("mock_v1")
  out <- compute_enrichment(character(), pw, universe = letters)
  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), length(pw))
  expect_true(all(out$n_overlap == 0L))
  # No selected, no overlap -> p = 1
  expect_true(all(abs(out$p_val - 1) < 1e-9))
})

test_that("compute_enrichment with empty pathway collection returns empty frame", {
  out <- compute_enrichment(c("CD3D", "MS4A1"), list(), universe = letters)
  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), 0L)
})

# ---- GSEA scaffold --------------------------------------------------------

test_that("empty_gsea_results has the canonical schema", {
  out <- empty_gsea_results()
  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), 0L)
  expect_true(all(c("pathway", "collection", "n_genes_in_pathway",
                    "n_leading_edge", "leading_edge_genes",
                    "ES", "NES", "p_val", "p_val_adj") %in% names(out)))
})

test_that("compute_gsea errors clearly when fgsea is missing", {
  skip_if(has_optional("fgsea"), "fgsea installed; install-error path skipped")
  ranked <- c(a = 2, b = 1, c = -1)
  pw <- list(P1 = c("a", "b"))
  expect_error(compute_gsea(ranked, pw, collection = "demo"),
               regexp = "fgsea", fixed = TRUE)
})

test_that("compute_gsea requires a NAMED numeric vector", {
  # The NAMED check fires *before* require_optional("fgsea"), so this
  # test runs whether or not fgsea is installed.
  pw <- list(P = c("a", "b"))
  expect_error(compute_gsea(c(1, 2, 3), pw), regexp = "NAMED", fixed = TRUE)
})

test_that("compute_gsea short-circuits on empty input", {
  expect_identical(nrow(compute_gsea(numeric(), list(P = "a"))), 0L)
  expect_identical(nrow(compute_gsea(c(a = 1), list())), 0L)
})

test_that(".fgsea_to_gsea_schema converts an fgsea-shaped result", {
  res <- data.frame(
    pathway = c("P1", "P2"),
    ES      = c(0.6, -0.3),
    NES     = c(1.8, -1.1),
    pval    = c(0.001, 0.4),
    padj    = c(0.01,  0.4),
    size    = c(10, 5),
    stringsAsFactors = FALSE
  )
  # Hand-attach a list-column for leadingEdge -- this is how fgsea ships it.
  res$leadingEdge <- list(c("CD3D", "CD3E"), c("MS4A1"))
  pw <- list(P1 = c("CD3D", "CD3E", "CD3G"), P2 = c("MS4A1", "CD19"))

  out <- .fgsea_to_gsea_schema(res, pathways = pw, collection = "demo")
  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), 2L)
  expect_identical(sort(out$pathway), c("P1", "P2"))
  expect_true(all(out$collection == "demo"))
  # Sorted ascending by adjusted p (P1 first).
  expect_identical(out$pathway[1], "P1")
  # Schema mapping: ES, NES, p_val, p_val_adj
  p1 <- out[out$pathway == "P1", ]
  expect_equal(p1$ES, 0.6)
  expect_equal(p1$NES, 1.8)
  expect_equal(p1$p_val, 0.001)
  expect_equal(p1$p_val_adj, 0.01)
  # Leading-edge gene count + concatenation
  expect_identical(p1$n_leading_edge, 2L)
  expect_identical(p1$leading_edge_genes, "CD3D;CD3E")
  # n_genes_in_pathway reflects the gene-set size (after de-dup), not the
  # `size` column of the fgsea result, so the function is robust to
  # filtering that may have already happened upstream.
  expect_identical(p1$n_genes_in_pathway, 3L)
})

test_that(".fgsea_to_gsea_schema handles missing leadingEdge column", {
  res <- data.frame(
    pathway = "P", ES = 0.3, NES = 1.1, pval = 0.05, padj = 0.05,
    stringsAsFactors = FALSE
  )
  out <- .fgsea_to_gsea_schema(res, pathways = list(P = c("A", "B")))
  expect_identical(out$n_leading_edge, 0L)
  expect_identical(out$leading_edge_genes, "")
})

test_that(".fgsea_to_gsea_schema errors on missing required columns", {
  res <- data.frame(pathway = "P", ES = 0.1, NES = 0.5, stringsAsFactors = FALSE)
  expect_error(.fgsea_to_gsea_schema(res, pathways = list(P = "A")),
               regexp = "pval", fixed = TRUE)
})

test_that(".fgsea_to_gsea_schema returns empty frame for empty input", {
  out <- .fgsea_to_gsea_schema(NULL, pathways = list())
  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), 0L)
})
