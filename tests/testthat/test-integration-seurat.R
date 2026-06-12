# ============================================================================
# End-to-end integration test against a real Seurat 5 object built from
# scratch (no remote download, no optional Bioc deps beyond SeuratObject).
#
# This is the "real-data smoke test" referenced in CURSOR_PROMPTS.md: it
# exercises load_seurat() -> markers -> DE -> pathway -> annotation ->
# regulons -> trajectory -> apply_*_to_dataset on a synthesised but
# biologically structured 3-population mini-PBMC. The integration test
# in `test-integration-10x.R` covers the same path through load_10x().
#
# The test is skipped automatically when Seurat is not installed so the
# default `Rscript tests/testthat.R` stays green on minimal CI runners.
# ============================================================================

skip_if_no_seurat <- function() {
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("Matrix")
}

# ---- Build a 3-population mini-PBMC ----------------------------------------
# Six canonical marker genes (one per known PBMC population) plus
# uninformative noise. Three cell types -- T, B, Myeloid -- are simulated
# with Poisson-distributed marker expression to give Seurat's normal
# pipeline (Normalise/HVG/Scale/PCA/UMAP/FindClusters) enough signal to
# carve out multiple clusters.
.build_smoke_seurat <- function(n_cells = 600L, n_genes = 1200L, seed = 7L) {
  set.seed(seed)
  cells <- sprintf("cell_%04d", seq_len(n_cells))
  genes <- c(c("CD3D", "MS4A1", "LST1", "EPCAM", "COL1A1", "NKG7"),
             sprintf("GENE_%05d", seq_len(n_genes - 6L)))
  ct <- sample(c("T", "B", "Myeloid"), n_cells, replace = TRUE,
               prob = c(0.5, 0.3, 0.2))
  mu_by_ct <- list(
    T       = c(CD3D = 6, NKG7 = 3, MS4A1 = 0.1, LST1 = 0.1,
                EPCAM = 0.05, COL1A1 = 0.05),
    B       = c(CD3D = 0.1, NKG7 = 0.1, MS4A1 = 6, LST1 = 0.1,
                EPCAM = 0.05, COL1A1 = 0.05),
    Myeloid = c(CD3D = 0.1, NKG7 = 0.1, MS4A1 = 0.1, LST1 = 6,
                EPCAM = 0.05, COL1A1 = 0.05))
  counts <- matrix(0L, nrow = n_genes, ncol = n_cells,
                   dimnames = list(genes, cells))
  for (j in seq_len(n_cells)) {
    mu <- mu_by_ct[[ct[j]]]
    for (g in names(mu)) counts[g, j] <- rpois(1, mu[g])
    noisy <- sample(seq.int(7L, n_genes), 200L)
    counts[noisy, j] <- rpois(200L, 0.2)
  }
  list(counts = methods::as(counts, "CsparseMatrix"),
       cells = cells, genes = genes, true_ct = ct)
}

.build_smoke_seurat_object <- function() {
  fix <- .build_smoke_seurat()
  seu <- suppressWarnings(SeuratObject::CreateSeuratObject(
    counts = fix$counts, project = "smoke",
    min.cells = 0, min.features = 0))
  seu$true_ct   <- fix$true_ct
  seu$sample    <- sample(c("S1", "S2"), length(fix$cells), replace = TRUE)
  seu$condition <- sample(c("ctrl", "treat"), length(fix$cells), replace = TRUE)
  seu <- suppressWarnings(Seurat::NormalizeData(seu, verbose = FALSE))
  seu <- suppressWarnings(Seurat::FindVariableFeatures(seu, nfeatures = 200,
                                                       verbose = FALSE))
  seu <- suppressWarnings(Seurat::ScaleData(seu, verbose = FALSE))
  seu <- suppressWarnings(Seurat::RunPCA(seu, npcs = 15L, verbose = FALSE))
  seu <- suppressWarnings(Seurat::FindNeighbors(seu, dims = 1:15,
                                                verbose = FALSE))
  # Resolution 1.2 reliably gives >1 cluster on this synthetic dataset,
  # which is what the marker_score engine needs to differentiate
  # populations (see the single-cluster warning in
  # .run_marker_score_annotation()).
  seu <- suppressWarnings(Seurat::FindClusters(seu, resolution = 1.2,
                                               verbose = FALSE))
  seu <- suppressWarnings(Seurat::RunUMAP(seu, dims = 1:15, verbose = FALSE))
  list(seurat = seu, true_ct = fix$true_ct)
}

test_that("load_seurat round-trips a fully-pipelined Seurat 5 object", {
  skip_if_no_seurat()
  fx <- .build_smoke_seurat_object()
  rds <- tempfile(fileext = ".rds")
  on.exit(unlink(rds), add = TRUE)
  saveRDS(fx$seurat, rds)

  ds <- load_seurat(rds)
  expect_equal(ds$source, "seurat")
  expect_equal(ds$n_cells, ncol(fx$seurat))
  expect_true(all(c("data", "counts") %in% backend_available_layers(ds$expression)))
  expect_true(all(c("PCA", "UMAP") %in% ds$reductions))
  expect_true("true_ct" %in% ds$metadata_fields)

  e <- get_gene_expression(ds, "CD3D")
  expect_equal(length(e), ds$n_cells)
  # CD3D should be higher in T than in B on the round-tripped data
  expect_gt(mean(e[ds$cell_data$true_ct == "T"]),
            mean(e[ds$cell_data$true_ct == "B"]))
})

test_that("marker / DE / pathway analysis chains run end-to-end on real Seurat data", {
  skip_if_no_seurat()
  fx <- .build_smoke_seurat_object()
  rds <- tempfile(fileext = ".rds")
  on.exit(unlink(rds), add = TRUE)
  saveRDS(fx$seurat, rds)
  ds <- load_seurat(rds)

  # compute_markers ranks CD3D first for the T-cell partition
  m <- compute_markers(ds, grouping_field = "true_ct",
                       group_filter = "T", top_n = 50L, test = "wilcox")
  expect_true(is.data.frame(m))
  expect_true(all(m$group == "T"))
  top_by_fc <- m$gene[order(-m$avg_log2FC)]
  expect_equal(top_by_fc[1], "CD3D")

  # compute_de via wilcox_r returns a data.frame with the documented
  # column schema (gene, group_1, group_2, avg_log2FC, pct.1, pct.2,
  # p_val, p_val_adj).
  de <- compute_de(ds, grouping_field = "true_ct",
                   group_1 = "T", group_2 = "B",
                   backend = "wilcox_r", min_pct = 0.1, test = "wilcox")
  expect_true(is.data.frame(de))
  expect_true(all(c("gene", "avg_log2FC", "p_val", "p_val_adj")
                  %in% names(de)))
  sig <- subset(de, p_val_adj < 0.05 & avg_log2FC > 0)
  expect_gt(nrow(sig), 0)
  expect_true("CD3D" %in% sig$gene)

  # Pathway ORA against the default builtin collection lands at least
  # one row and resolves T-cell signal at the top.
  coll <- available_pathway_collections()[1]
  pw <- compute_enrichment(
    selected = sig$gene[order(sig$p_val_adj)][seq_len(min(50L, nrow(sig)))],
    pathways = get_pathways(coll), universe = ds$genes)
  expect_true(is.data.frame(pw))
  expect_gt(nrow(pw), 0)
  expect_true(any(grepl("T cell", pw$pathway, ignore.case = TRUE)))
})

test_that("annotation / regulons / trajectory chain runs and bakes into dataset metadata", {
  skip_if_no_seurat()
  fx <- .build_smoke_seurat_object()
  rds <- tempfile(fileext = ".rds")
  on.exit(unlink(rds), add = TRUE)
  saveRDS(fx$seurat, rds)
  ds <- load_seurat(rds)

  # marker_score annotation on the Seurat clusters
  state <- new_app_state()
  shiny::isolate(set_active_dataset(state, ds))
  ann <- shiny::isolate(run_annotation_engine(
    "marker_score", ds, state,
    params = list(cluster_field = "seurat_clusters",
                  species = NA_character_, tissue = NA_character_,
                  min_score = 0.0),
    set_id = "smoke_marker"))
  expect_true(is_annotation_result_v1(ann))
  expect_true(nrow(ann$cluster_summary) >= 2L)
  # On a multi-cluster Seurat run we expect at least one T-cell label
  # in the top picks -- the dataset's strongest signal.
  expect_true(any(grepl("T cell", ann$cluster_summary$top_label,
                        ignore.case = TRUE)))

  # AUCell pure-R regulon scoring on the mock_pbmc set picks out the
  # expected TFs per true cell type (GATA3 -> T, PAX5 -> B,
  # SPI1 -> Myeloid). This is the strictest biological check in the
  # suite; it also confirms get_gene_expression is layer-correct on
  # the Seurat-derived sparse backend.
  reg <- run_regulon_engine("aucell_r", ds,
                            fetch_regulon_set("mock_pbmc"),
                            params = list(top_n_fraction = 0.05))
  expect_true(is_regulon_result_v1(reg))
  by_ct <- regulon_mean_by_group(reg, ds$cell_data$true_ct)
  expect_identical(rownames(by_ct)[which.max(by_ct[, "GATA3"])], "T")
  expect_identical(rownames(by_ct)[which.max(by_ct[, "PAX5"])],  "B")
  expect_identical(rownames(by_ct)[which.max(by_ct[, "SPI1"])],  "Myeloid")

  # Mock trajectory + bake pseudotime + annotation back into the
  # dataset's metadata.
  tr <- run_trajectory(ds, source = "mock", reduction = "UMAP",
                       root_field = "true_ct", root_group = "T")
  expect_true(is_trajectory_result(tr))
  ds2 <- apply_pseudotime_to_dataset(ds, tr, bins = 5L)
  ds3 <- apply_annotations_to_dataset(ds2, ann)
  new_cols <- setdiff(ds3$metadata_fields, ds$metadata_fields)
  # Three new columns: pseudotime__*, pseudotime_bin__*, annotation__*
  expect_length(new_cols, 3L)
  expect_true(any(grepl("^pseudotime__",     new_cols)))
  expect_true(any(grepl("^pseudotime_bin__", new_cols)))
  expect_true(any(grepl("^annotation__",     new_cols)))
})

test_that("marker_score warns when the cluster_field has only one cluster", {
  skip_if_no_seurat()
  fx <- .build_smoke_seurat_object()
  rds <- tempfile(fileext = ".rds")
  on.exit(unlink(rds), add = TRUE)
  saveRDS(fx$seurat, rds)
  ds <- load_seurat(rds)
  # Synthesise a degenerate cluster field
  ds$cell_data$one_cluster <- factor(rep("0", ds$n_cells))
  ds$metadata_fields <- union(ds$metadata_fields, "one_cluster")

  state <- new_app_state()
  shiny::isolate(set_active_dataset(state, ds))
  ann <- shiny::isolate(run_annotation_engine(
    "marker_score", ds, state,
    params = list(cluster_field = "one_cluster",
                  species = NA_character_, tissue = NA_character_,
                  min_score = 0.0),
    set_id = "smoke_one_cluster"))
  expect_true(is_annotation_result_v1(ann))
  # Warning surfaces the degenerate-clustering caveat
  expect_true(any(grepl("only|1 cluster|cluster.*has 1",
                        ann$warnings, ignore.case = TRUE)))
})
