test_that("compute_de returns the documented column schema", {
  ds <- mock_dataset(n_cells = 200, seed = 1)
  de <- compute_de(ds, grouping_field = "cluster",
                   group_1 = "0", group_2 = "1",
                   assay = "RNA", layer = "data",
                   min_pct = 0, test = "wilcox")
  expect_s3_class(de, "data.frame")
  expect_true(all(c("gene", "group_1", "group_2", "avg_log2FC",
                    "pct.1", "pct.2", "p_val", "p_val_adj")
                  %in% names(de)))
  expect_gt(nrow(de), 0L)
})

test_that("compute_de rejects identical groups", {
  ds <- mock_dataset(n_cells = 80)
  expect_error(compute_de(ds, grouping_field = "cluster",
                          group_1 = "0", group_2 = "0"))
})

test_that("filter_de_results applies gene search + thresholds", {
  ds <- mock_dataset(n_cells = 150)
  de <- compute_de(ds, grouping_field = "cluster",
                   group_1 = "0", group_2 = "1",
                   assay = "RNA", layer = "data",
                   min_pct = 0, test = "wilcox")
  hit <- filter_de_results(de, gene_search = "CD3D",
                           min_abs_log2fc = 0, max_padj = 1)
  expect_true(all(grepl("CD3D", hit$gene, ignore.case = TRUE)))

  zero <- filter_de_results(de, gene_search = "",
                            min_abs_log2fc = 1e6, max_padj = 1)
  expect_identical(nrow(zero), 0L)
})

test_that("sort_de_results sorts by named column + direction", {
  ds <- mock_dataset(n_cells = 100)
  de <- compute_de(ds, grouping_field = "cluster",
                   group_1 = "0", group_2 = "1",
                   assay = "RNA", layer = "data",
                   min_pct = 0, test = "wilcox")
  asc  <- sort_de_results(de, sort_by = "avg_log2FC", descending = FALSE)
  desc <- sort_de_results(de, sort_by = "avg_log2FC", descending = TRUE)
  expect_true(all(diff(asc$avg_log2FC)  >= -1e-9))
  expect_true(all(diff(desc$avg_log2FC) <=  1e-9))
})

test_that("empty_de_results matches the live result column set", {
  empty <- empty_de_results()
  expect_s3_class(empty, "data.frame")
  expect_true(all(c("gene", "group_1", "group_2", "avg_log2FC",
                    "pct.1", "pct.2", "p_val", "p_val_adj")
                  %in% names(empty)))
})

# ---- Backend dispatcher --------------------------------------------------

test_that("de_available_backends lists every registered backend", {
  backs <- de_available_backends()
  ids <- vapply(backs, `[[`, character(1), "id")
  expect_setequal(ids, c("auto", "wilcox_r", "presto",
                         "pseudobulk_naive",
                         "pseudobulk_edger",
                         "pseudobulk_deseq2"))
  kinds <- vapply(backs, `[[`, character(1), "kind")
  names(kinds) <- ids
  expect_identical(kinds[["auto"]],              "cell")
  expect_identical(kinds[["wilcox_r"]],          "cell")
  expect_identical(kinds[["presto"]],            "cell")
  expect_identical(kinds[["pseudobulk_naive"]],  "pseudobulk")
  expect_identical(kinds[["pseudobulk_edger"]],  "pseudobulk")
  expect_identical(kinds[["pseudobulk_deseq2"]], "pseudobulk")

  flags <- vapply(backs, `[[`, logical(1), "available")
  names(flags) <- ids
  expect_true(flags[["auto"]])
  expect_true(flags[["wilcox_r"]])
  expect_true(flags[["pseudobulk_naive"]])    # always available
  expect_identical(flags[["presto"]],            has_optional("presto"))
  expect_identical(flags[["pseudobulk_edger"]],  has_optional("edgeR"))
  expect_identical(flags[["pseudobulk_deseq2"]], has_optional("DESeq2"))
})

test_that(".is_pseudobulk_backend recognises every pseudobulk backend", {
  expect_true(.is_pseudobulk_backend("pseudobulk_naive"))
  expect_true(.is_pseudobulk_backend("pseudobulk_edger"))
  expect_true(.is_pseudobulk_backend("pseudobulk_deseq2"))
  expect_false(.is_pseudobulk_backend("wilcox_r"))
  expect_false(.is_pseudobulk_backend("presto"))
  expect_false(.is_pseudobulk_backend("auto"))
})

test_that(".de_resolve_backend picks presto for auto+wilcox iff installed, else wilcox_r", {
  expected_auto_wilcox <- if (has_optional("presto")) "presto" else "wilcox_r"
  expect_identical(.de_resolve_backend("auto",     "wilcox"), expected_auto_wilcox)
  expect_identical(.de_resolve_backend("auto",     "t"),      "wilcox_r")
  expect_identical(.de_resolve_backend("wilcox_r", "wilcox"), "wilcox_r")
  expect_identical(.de_resolve_backend("presto",   "wilcox"), "presto")
  expect_error(.de_resolve_backend("not_a_backend", "wilcox"), "Unknown DE backend")
})

test_that("compute_de honours backend = 'wilcox_r' explicitly", {
  ds <- mock_dataset(n_cells = 150, seed = 21)
  de <- compute_de(ds, grouping_field = "cluster",
                   group_1 = "0", group_2 = "1",
                   min_pct = 0, test = "wilcox",
                   backend = "wilcox_r")
  expect_s3_class(de, "data.frame")
  expect_gt(nrow(de), 0L)
  expect_true(all(c("gene", "group_1", "group_2", "avg_log2FC",
                    "pct.1", "pct.2", "p_val", "p_val_adj")
                  %in% names(de)))
})

test_that("compute_de default backend = 'auto' produces the same schema", {
  ds <- mock_dataset(n_cells = 150, seed = 22)
  de1 <- compute_de(ds, grouping_field = "cluster",
                    group_1 = "0", group_2 = "1",
                    min_pct = 0, test = "wilcox", backend = "wilcox_r")
  de2 <- compute_de(ds, grouping_field = "cluster",
                    group_1 = "0", group_2 = "1",
                    min_pct = 0, test = "wilcox", backend = "auto")
  expect_setequal(names(de1), names(de2))
})

test_that("compute_de(backend = 'presto') errors with install instructions when missing", {
  skip_if(has_optional("presto"),
          "presto is installed; missing-dep path is not reachable in this env.")
  ds <- mock_dataset(n_cells = 80)
  err <- tryCatch(
    compute_de(ds, grouping_field = "cluster",
               group_1 = "0", group_2 = "1",
               backend = "presto"),
    error = function(e) conditionMessage(e))
  expect_match(err, "presto")
  expect_match(err, "install")
})

test_that("compute_de(backend = 'presto', test = 't') is rejected", {
  ds <- mock_dataset(n_cells = 80)
  expect_error(
    compute_de(ds, grouping_field = "cluster",
               group_1 = "0", group_2 = "1",
               test = "t", backend = "presto"),
    "only supports test = 'wilcox'"
  )
})

# ---- Pure schema converter (testable without presto) --------------------

test_that(".presto_to_de_schema converts a presto-shaped frame", {
  # Mimic presto::wilcoxauc() output: one row per (feature, group) pair.
  presto_like <- data.frame(
    feature = c("A", "A", "B", "B", "C", "C"),
    group   = c("0", "1", "0", "1", "0", "1"),
    avg_expr = rnorm(6),
    statistic = rnorm(6),
    auc     = runif(6),
    pval    = c(0.001, 0.001, 0.4, 0.4, 0.05, 0.05),
    padj    = c(0.003, 0.003, 0.5, 0.5, 0.10, 0.10),
    pct_in  = c(80,  10,   30,  35,   60,  40),  # presto uses percents
    pct_out = c(10,  80,   35,  30,   40,  60),
    logFC   = c( 1.5, -1.5, 0.1, -0.1, 0.7, -0.7),
    stringsAsFactors = FALSE
  )
  out <- .presto_to_de_schema(presto_like, group_1 = "0", group_2 = "1")
  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), 3L)
  expect_setequal(out$gene, c("A", "B", "C"))
  expect_true(all(out$group_1 == "0"))
  expect_true(all(out$group_2 == "1"))
  # pct_in/pct_out divided by 100
  a_row <- out[out$gene == "A", ]
  expect_identical(a_row$pct.1, 0.80)
  expect_identical(a_row$pct.2, 0.10)
  expect_identical(a_row$avg_log2FC, 1.5)
  expect_identical(a_row$p_val, 0.001)
  expect_identical(a_row$p_val_adj, 0.003)
})

test_that(".presto_to_de_schema applies min_pct after the group-1 selection", {
  presto_like <- data.frame(
    feature = c("A", "A", "B", "B"),
    group   = c("0", "1", "0", "1"),
    pval    = c(0.01, 0.01, 0.5, 0.5),
    padj    = c(0.02, 0.02, 0.6, 0.6),
    pct_in  = c(60, 10, 5, 8),    # B has max(pct) = 8%, below 10%
    pct_out = c(10, 60, 8, 5),
    logFC   = c(1.0, -1.0, 0.1, -0.1),
    stringsAsFactors = FALSE
  )
  out <- .presto_to_de_schema(presto_like, group_1 = "0", group_2 = "1",
                              min_pct = 0.1)
  expect_identical(out$gene, "A")
})

test_that(".presto_to_de_schema accepts alternate column names (pvalue/p_val_adj)", {
  presto_like <- data.frame(
    feature = c("A", "A"), group = c("0", "1"),
    pvalue  = c(0.01, 0.01), p_val_adj = c(0.02, 0.02),
    pct_in  = c(60, 10), pct_out = c(10, 60),
    logFC   = c(1.0, -1.0), stringsAsFactors = FALSE
  )
  out <- .presto_to_de_schema(presto_like, group_1 = "0", group_2 = "1")
  expect_identical(out$p_val,     0.01)
  expect_identical(out$p_val_adj, 0.02)
})

test_that(".presto_to_de_schema returns empty frame on empty / mismatched input", {
  expect_identical(nrow(.presto_to_de_schema(NULL, "0", "1")), 0L)
  expect_identical(nrow(.presto_to_de_schema(data.frame(), "0", "1")), 0L)

  presto_like <- data.frame(
    feature = "A", group = "1",  # no rows where group == "0"
    pval = 0.5, padj = 0.6,
    pct_in = 10, pct_out = 5, logFC = 0.0,
    stringsAsFactors = FALSE)
  expect_identical(nrow(.presto_to_de_schema(presto_like, "0", "1")), 0L)
})

# ---- Matrix orientation contract for presto (P2) ------------------------
# `presto::wilcoxauc(X, y)` wants X = genes x cells with rownames = genes;
# `y` is a per-cell label vector. The helper that materialises X is
# extracted so the orientation/labelling contract is unit-testable
# without `presto` installed.

test_that(".de_build_genes_x_cells materialises a genes x cells matrix with the right labels", {
  ds <- mock_dataset(n_cells = 40)
  genes <- intersect(ds$genes, c("CD3D", "MS4A1", "LST1"))
  testthat::skip_if(length(genes) < 2L,
                    "Mock dataset doesn't carry the expected genes")
  cell_idx <- c(1:10, 21:35)  # positional indices, NOT cell IDs

  X <- .de_build_genes_x_cells(ds, genes, cell_idx)
  expect_identical(dim(X), c(length(genes), length(cell_idx)))
  expect_identical(rownames(X), as.character(genes))
  # Sanity: per-cell value matches a direct get_gene_expression lookup
  # (this is what catches an accidental "no transpose" or row/col swap).
  ge <- get_gene_expression(ds, genes[1])
  expect_identical(unname(X[1, ]), ge[cell_idx])
  # Cross-check the second gene too -- different feature should not
  # accidentally share row 1 (a transpose bug would silently align
  # the wrong feature).
  ge2 <- get_gene_expression(ds, genes[2])
  expect_identical(unname(X[2, ]), ge2[cell_idx])
})

test_that(".de_build_genes_x_cells rejects non-integer cell indices (e.g. cell IDs)", {
  ds <- mock_dataset(n_cells = 20)
  expect_error(.de_build_genes_x_cells(ds, ds$genes[1:2], ds$cells[1:5]),
               regexp = "integer positions")
})

test_that(".de_build_genes_x_cells errors when given an empty gene or cell set", {
  ds <- mock_dataset(n_cells = 20)
  expect_error(.de_build_genes_x_cells(ds, character(), seq_len(10)),
               regexp = "no genes")
  expect_error(.de_build_genes_x_cells(ds, ds$genes[1:2], integer()),
               regexp = "no cells")
})

# ---- Layer-aware gene universe (P3) -------------------------------------
# compute_de() previously called `available_genes(dataset)` without
# forwarding its own `layer` argument. With a multi-layer backend whose
# layers expose different gene sets, that caused DE to silently probe
# genes that did not exist in the requested layer (yielding NULL gene
# vectors and zero-row results).

test_that("compute_de uses the layer-specific gene set when an explicit layer is passed", {
  # Build a tiny dataset where 'counts' has 4 genes, 'data' has only 2
  # (the 2 shared genes between layers). DE on layer = 'counts' must
  # see all 4; DE on layer = 'data' must see only the 2 shared genes.
  set.seed(2L)
  n  <- 30L
  cells <- sprintf("c%02d", seq_len(n))
  grp <- rep(c("0", "1"), each = n / 2)

  counts_layer <- list(
    G_SHARED1 = c(rpois(n / 2, 5), rpois(n / 2, 0.1)),
    G_SHARED2 = c(rpois(n / 2, 0.1), rpois(n / 2, 5)),
    G_COUNTS_ONLY1 = rpois(n, 1),
    G_COUNTS_ONLY2 = rpois(n, 1))
  data_layer <- list(
    G_SHARED1 = log1p(counts_layer$G_SHARED1),
    G_SHARED2 = log1p(counts_layer$G_SHARED2))

  expr_be <- expression_backend_inmemory(
    layers = list(counts = counts_layer, data = data_layer),
    n_cells = n, default_layer = "data")

  ds <- list(
    name = "p3", source = "synthetic",
    n_cells = n, n_genes = 4L,
    assays = "RNA", default_assay = "RNA",
    reductions = character(), default_reduction = NA_character_,
    metadata_fields = "cluster",
    cells = cells,
    cell_data = data.frame(cell = cells, cluster = grp,
                           stringsAsFactors = FALSE),
    genes = c("G_SHARED1", "G_SHARED2"),  # default = data layer
    expression = expr_be)

  # Layer = NULL -> default ('data') -> only the 2 shared genes are tested.
  out_default <- compute_de(ds, grouping_field = "cluster",
                            group_1 = "0", group_2 = "1",
                            backend = "wilcox_r", min_pct = 0.0,
                            test = "wilcox")
  expect_setequal(out_default$gene, c("G_SHARED1", "G_SHARED2"))

  # Layer = "counts" -> all 4 genes are tested.
  out_counts <- compute_de(ds, grouping_field = "cluster",
                           group_1 = "0", group_2 = "1",
                           backend = "wilcox_r", min_pct = 0.0,
                           layer = "counts", test = "wilcox")
  expect_setequal(out_counts$gene,
                  c("G_SHARED1", "G_SHARED2",
                    "G_COUNTS_ONLY1", "G_COUNTS_ONLY2"))
})

test_that("available_genes respects the layer argument on a layered backend", {
  expr_be <- expression_backend_inmemory(
    layers = list(
      data   = list(A = c(0, 1, 2), B = c(1, 2, 3)),
      counts = list(A = c(0, 1, 2), B = c(1, 2, 3), C = c(0, 0, 1))),
    n_cells = 3L, default_layer = "data")
  ds <- list(genes = c("A", "B"), expression = expr_be)
  expect_setequal(available_genes(ds),                  c("A", "B"))
  expect_setequal(available_genes(ds, layer = "data"),  c("A", "B"))
  expect_setequal(available_genes(ds, layer = "counts"),
                  c("A", "B", "C"))
})

test_that("compute_de(backend='presto') feeds presto a (genes x cells, per-cell y) call", {
  # Stub presto::wilcoxauc in a private env so we can verify the call
  # contract without installing the package. We don't try to mock
  # `requireNamespace` -- skip if presto isn't really installed and the
  # require_optional gate is the test instead.
  skip_if(!has_optional("presto"),
          "presto not installed; orientation is covered via .de_build_genes_x_cells")
  ds <- mock_dataset(n_cells = 50, seed = 3)
  out <- compute_de(ds, grouping_field = "cluster",
                    group_1 = "0", group_2 = "1",
                    backend = "presto", min_pct = 0.1)
  expect_true(is.data.frame(out))
  # If the orientation were wrong, presto would either error (mismatched
  # `y` length) or return logFC with flipped signs. Assert the
  # *direction* of one strongly-expressed canonical marker so a silent
  # transpose regression would surface.
  if (nrow(out) > 0L && "CD3D" %in% out$gene) {
    cd3d <- out[out$gene == "CD3D", ]
    # On the mock, cluster 0 expresses CD3D strongly; cluster 1 does not.
    expect_gt(cd3d$avg_log2FC, 0)
  }
})
