# Tests for R/pseudobulk.R
# ----------------------------------------------------------------------------
# Aggregation correctness, validation errors, end-to-end pseudobulk_naive
# through compute_de(), and the missing-dep error path for edgeR/DESeq2.

# ---- aggregate_pseudobulk -------------------------------------------------

test_that("aggregate_pseudobulk: shape + sample_metadata + counts default", {
  ds <- mock_dataset(n_cells = 400, seed = 17)
  pb <- aggregate_pseudobulk(ds,
    grouping_field = "condition",
    group_1 = "treat", group_2 = "ctrl",
    sample_by = "sample",
    min_cells_per_sample = 5L)

  expect_type(pb, "list")
  expect_true(all(c("matrix", "sample_metadata", "layer_used", "agg",
                    "warn_lognorm", "provenance") %in% names(pb)))

  # Default layer should be "counts" (the mock now exposes it).
  expect_identical(pb$layer_used, "counts")
  expect_false(pb$warn_lognorm)
  expect_identical(pb$agg, "sum")

  # sample_metadata columns
  expect_true(all(c("pb_sample", "sample", "group", "n_cells") %in%
                    names(pb$sample_metadata)))
  expect_true(all(pb$sample_metadata$n_cells >= 5L))
  # Each (sample x condition) combination should be present at most once.
  expect_identical(anyDuplicated(pb$sample_metadata$pb_sample), 0L)
  # Both groups should be represented.
  expect_true(all(c("treat", "ctrl") %in% pb$sample_metadata$group))

  # Matrix has rows = genes, cols = pseudobulk samples (in metadata order).
  expect_identical(rownames(pb$matrix), backend_genes(ds$expression, layer = "counts"))
  expect_identical(colnames(pb$matrix), pb$sample_metadata$pb_sample)
  expect_true(all(pb$matrix >= 0))           # counts are non-negative
  expect_true(all(pb$matrix == round(pb$matrix))) # and integer-valued
})

test_that("aggregate_pseudobulk: sum aggregation matches manual per-cell sum", {
  ds <- mock_dataset(n_cells = 200, seed = 11)
  pb <- aggregate_pseudobulk(ds,
    grouping_field = "condition",
    group_1 = "treat", group_2 = "ctrl",
    sample_by = "sample",
    min_cells_per_sample = 5L)

  # Spot-check one cell: re-aggregate manually and compare.
  k    <- pb$sample_metadata$pb_sample[1]
  meta_g <- get_metadata(ds, "condition")
  meta_s <- get_metadata(ds, "sample")
  parts  <- strsplit(k, "__", fixed = TRUE)[[1]]
  cells_in_k <- which(meta_s == parts[1] & meta_g == parts[2])
  manual <- sum(get_gene_expression(ds, "CD3D", layer = "counts")[cells_in_k])
  expect_identical(pb$matrix["CD3D", k], manual)
})

test_that("aggregate_pseudobulk: mean agg honoured and warn_lognorm for non-counts", {
  ds <- mock_dataset(n_cells = 200, seed = 13)
  pb <- aggregate_pseudobulk(ds,
    grouping_field = "condition", group_1 = "treat", group_2 = "ctrl",
    sample_by = "sample", layer = "data", agg = "mean",
    min_cells_per_sample = 5L)
  expect_identical(pb$layer_used, "data")
  expect_true(pb$warn_lognorm)
  expect_identical(pb$agg, "mean")
  # Mean of non-negative values is also non-negative; sanity check.
  expect_true(all(pb$matrix >= 0))
})

test_that("aggregate_pseudobulk: drops pseudobulk samples below min_cells_per_sample", {
  ds <- mock_dataset(n_cells = 60, seed = 19)
  # Tighten min_cells so most combos drop out
  pb <- aggregate_pseudobulk(ds,
    grouping_field = "condition", group_1 = "treat", group_2 = "ctrl",
    sample_by = "sample", min_cells_per_sample = 50L)
  expect_true(is.null(pb) || all(pb$sample_metadata$n_cells >= 50L))
})

# ---- validate_pseudobulk_inputs ------------------------------------------

test_that("validate_pseudobulk_inputs: requires sample_by", {
  ds <- mock_dataset(n_cells = 100)
  expect_error(
    validate_pseudobulk_inputs(ds, "condition", "treat", "ctrl", NULL),
    regexp = "sample_by", fixed = TRUE
  )
  expect_error(
    validate_pseudobulk_inputs(ds, "condition", "treat", "ctrl", ""),
    regexp = "sample_by", fixed = TRUE
  )
})

test_that("validate_pseudobulk_inputs: rejects identical grouping_field and sample_by", {
  ds <- mock_dataset(n_cells = 100)
  expect_error(
    validate_pseudobulk_inputs(ds, "sample", "S1", "S2", "sample"),
    regexp = "cannot be the same column", fixed = TRUE
  )
})

test_that("validate_pseudobulk_inputs: errors on too few samples per group", {
  ds <- mock_dataset(n_cells = 100, seed = 23)
  # Force an impossible min by demanding hundreds of cells per sample.
  expect_error(
    validate_pseudobulk_inputs(ds, "condition", "treat", "ctrl", "sample",
                               min_cells_per_sample = 9999L,
                               min_samples_per_group = 2L),
    regexp = "Pseudobulk DE needs at least", fixed = FALSE
  )
})

test_that("validate_pseudobulk_inputs: passes silently on a healthy mock dataset", {
  ds <- mock_dataset(n_cells = 600, seed = 27)
  expect_invisible(validate_pseudobulk_inputs(
    ds, "condition", "treat", "ctrl", "sample",
    min_cells_per_sample = 5L, min_samples_per_group = 2L))
})

# ---- pseudobulk_cpm_log2 + .pseudobulk_pct ------------------------------

test_that("pseudobulk_cpm_log2 normalises columns to log2(CPM+1)", {
  M <- matrix(c(1, 1, 0,
                3, 2, 1,
                0, 4, 5), nrow = 3, byrow = TRUE,
              dimnames = list(c("A", "B", "C"), c("s1", "s2", "s3")))
  out <- pseudobulk_cpm_log2(M)
  expect_identical(dim(out), dim(M))
  # The sum of CPM per column equals 1e6, so log2(CPM+1) columns should
  # share an upper bound around log2(1e6+1) ~= 19.93.
  expect_true(all(out <= log2(1e6 + 1) + 1e-6))
  # Zero counts in column 1 row C -> log2(0+1) = 0.
  expect_equal(out["C", "s1"], 0)
})

test_that(".pseudobulk_pct returns named vectors keyed on every gene", {
  ds <- mock_dataset(n_cells = 120, seed = 31)
  meta <- get_metadata(ds, "condition")
  in1 <- which(meta == "treat"); in2 <- which(meta == "ctrl")
  out <- .pseudobulk_pct(ds, in1, in2, layer = "counts")
  expect_setequal(names(out$pct.1), backend_genes(ds$expression, layer = "counts"))
  expect_setequal(names(out$pct.2), backend_genes(ds$expression, layer = "counts"))
  expect_true(all(out$pct.1 >= 0 & out$pct.1 <= 1))
  expect_true(all(out$pct.2 >= 0 & out$pct.2 <= 1))
})

# ---- .pseudobulk_to_de_schema --------------------------------------------

test_that(".pseudobulk_to_de_schema lifts to canonical DE schema + BH-adjusts", {
  df <- data.frame(
    gene       = c("A", "B", "C"),
    avg_log2FC = c(1.5, -0.2, 0.7),
    p_val      = c(0.001, 0.6, 0.04),
    pct.1      = c(0.8, 0.1, 0.5),
    pct.2      = c(0.2, 0.1, 0.4),
    stringsAsFactors = FALSE
  )
  out <- .pseudobulk_to_de_schema(df, group_1 = "g1", group_2 = "g2")
  expect_s3_class(out, "data.frame")
  expect_setequal(names(out), c("gene", "group_1", "group_2",
                                "avg_log2FC", "pct.1", "pct.2",
                                "p_val", "p_val_adj"))
  expect_true(all(out$group_1 == "g1"))
  expect_true(all(out$p_val_adj >= out$p_val - 1e-12))
})

test_that(".pseudobulk_to_de_schema applies min_pct filter", {
  df <- data.frame(
    gene = c("A", "B"),
    avg_log2FC = c(1.0, -1.0),
    p_val = c(0.01, 0.01),
    pct.1 = c(0.6, 0.05),
    pct.2 = c(0.5, 0.02),
    stringsAsFactors = FALSE
  )
  out <- .pseudobulk_to_de_schema(df, "g1", "g2", min_pct = 0.1)
  expect_identical(nrow(out), 1L)
  expect_identical(out$gene, "A")
})

test_that(".pseudobulk_to_de_schema fills missing pct columns with NA", {
  df <- data.frame(gene = "A", avg_log2FC = 1, p_val = 0.5,
                   stringsAsFactors = FALSE)
  out <- .pseudobulk_to_de_schema(df, "g1", "g2")
  expect_identical(nrow(out), 1L)
  expect_true(is.na(out$pct.1))
  expect_true(is.na(out$pct.2))
})

test_that(".pseudobulk_to_de_schema errors on missing required columns", {
  df <- data.frame(foo = "bar", stringsAsFactors = FALSE)
  expect_error(.pseudobulk_to_de_schema(df, "g1", "g2"),
               regexp = "missing columns", fixed = TRUE)
})

# ---- End-to-end via compute_de --------------------------------------------

test_that("compute_de(backend = 'pseudobulk_naive') round-trips through the canonical schema", {
  ds <- mock_dataset(n_cells = 800, seed = 37)
  de <- compute_de(ds,
    grouping_field = "condition",
    group_1 = "treat", group_2 = "ctrl",
    backend = "pseudobulk_naive",
    sample_by = "sample",
    min_cells_per_sample = 5L,
    min_samples_per_group = 2L,
    min_pct = 0)
  expect_s3_class(de, "data.frame")
  expect_true(nrow(de) >= 1L)
  # Canonical DE schema -- same column set as cell-level backends.
  expect_setequal(names(de), c("gene", "group_1", "group_2",
                               "avg_log2FC", "pct.1", "pct.2",
                               "p_val", "p_val_adj"))
  expect_true(all(de$group_1 == "treat"))
  expect_true(all(de$group_2 == "ctrl"))
  # Sorted by p_val_adj ascending.
  expect_true(all(diff(de$p_val_adj) >= -1e-9))
})

test_that("compute_de pseudobulk requires sample_by", {
  ds <- mock_dataset(n_cells = 200, seed = 41)
  expect_error(
    compute_de(ds, grouping_field = "condition",
               group_1 = "treat", group_2 = "ctrl",
               backend = "pseudobulk_naive"),
    regexp = "sample_by", fixed = TRUE)
})

test_that("compute_de(backend = 'pseudobulk_edger') errors cleanly without edgeR", {
  skip_if(has_optional("edgeR"), "edgeR installed; install-error path skipped")
  ds <- mock_dataset(n_cells = 400, seed = 43)
  expect_error(
    compute_de(ds, grouping_field = "condition",
               group_1 = "treat", group_2 = "ctrl",
               backend = "pseudobulk_edger",
               sample_by = "sample",
               min_cells_per_sample = 5L),
    regexp = "edgeR", fixed = FALSE)
})

test_that("compute_de(backend = 'pseudobulk_deseq2') errors cleanly without DESeq2", {
  skip_if(has_optional("DESeq2"), "DESeq2 installed; install-error path skipped")
  ds <- mock_dataset(n_cells = 400, seed = 47)
  expect_error(
    compute_de(ds, grouping_field = "condition",
               group_1 = "treat", group_2 = "ctrl",
               backend = "pseudobulk_deseq2",
               sample_by = "sample",
               min_cells_per_sample = 5L),
    regexp = "DESeq2", fixed = FALSE)
})

# End-to-end against real edgeR / DESeq2 if installed -- skipped here.
test_that("pseudobulk_edger end-to-end (requires edgeR)", {
  skip_if_not_installed("edgeR")
  ds <- mock_dataset(n_cells = 800, seed = 53)
  de <- compute_de(ds, grouping_field = "condition",
                   group_1 = "treat", group_2 = "ctrl",
                   backend = "pseudobulk_edger",
                   sample_by = "sample",
                   min_cells_per_sample = 5L, min_pct = 0)
  expect_s3_class(de, "data.frame")
  expect_setequal(names(de), c("gene", "group_1", "group_2",
                               "avg_log2FC", "pct.1", "pct.2",
                               "p_val", "p_val_adj"))
})

test_that("pseudobulk_deseq2 end-to-end (requires DESeq2)", {
  skip_if_not_installed("DESeq2")
  ds <- mock_dataset(n_cells = 800, seed = 59)
  de <- compute_de(ds, grouping_field = "condition",
                   group_1 = "treat", group_2 = "ctrl",
                   backend = "pseudobulk_deseq2",
                   sample_by = "sample",
                   min_cells_per_sample = 5L, min_pct = 0)
  expect_s3_class(de, "data.frame")
  expect_setequal(names(de), c("gene", "group_1", "group_2",
                               "avg_log2FC", "pct.1", "pct.2",
                               "p_val", "p_val_adj"))
})
