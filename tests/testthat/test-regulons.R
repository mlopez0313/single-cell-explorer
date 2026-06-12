# Tests for R/regulon_schema.R, R/regulon_registry.R, R/regulon_aucell.R,
# R/regulon_sources.R.
#
# Layered so each surface is exercised in isolation:
#
#   - schema (regulon_spec / regulon_set / regulon_result_v1 +
#     regulon_mean_by_group)
#   - sources (mock_pbmc builtin + DoRothEA seam + dorothea_df_to_regulons
#     pure converter)
#   - engines (registry shape, dispatcher, missing-dep, pure AUCell math)
#   - end-to-end via the mock dataset (regulons peak in expected clusters)

# ---- Schema ----------------------------------------------------------------

test_that("regulon_spec validates and stores TF + targets + weights + type", {
  r <- regulon_spec(tf = "GATA3", targets = c("CD3D", "NKG7"),
                    weights = c(0.9, 0.7), type = "activating")
  expect_identical(r$tf, "GATA3")
  expect_identical(r$targets, c("CD3D", "NKG7"))
  expect_identical(r$type, "activating")

  # Default weights = 1
  r2 <- regulon_spec("PAX5", targets = c("MS4A1"))
  expect_identical(r2$weights, 1.0)

  # Type must be one of the allowed values
  expect_error(regulon_spec("X", "Y", type = "boom"),
               regexp = "should be one of", fixed = FALSE)
  # tf must be non-empty
  expect_error(regulon_spec(tf = "", targets = "Y"),
               regexp = "is not TRUE", fixed = FALSE)
})

test_that("regulon_set rejects duplicate TFs and tags with schema version", {
  set <- regulon_set(
    id = "mini", name = "Mini",
    regulons = list(
      regulon_spec("GATA3", "CD3D"),
      regulon_spec("PAX5",  "MS4A1")
    ))
  expect_true(is_regulon_set(set))
  expect_identical(set$schema_version, REGULON_SET_SCHEMA_VERSION)
  expect_identical(length(set$regulons), 2L)

  expect_error(regulon_set(
    id = "dup", name = "Dup",
    regulons = list(regulon_spec("GATA3", "A"),
                    regulon_spec("GATA3", "B"))),
    regexp = "duplicate TF", fixed = TRUE)
})

test_that("regulon_set_as_target_list returns a named list of target vectors", {
  set <- regulon_set("mini", "Mini", regulons = list(
    regulon_spec("GATA3", c("CD3D", "NKG7")),
    regulon_spec("PAX5",  "MS4A1")
  ))
  lst <- regulon_set_as_target_list(set)
  expect_identical(names(lst), c("GATA3", "PAX5"))
  expect_identical(lst$GATA3, c("CD3D", "NKG7"))
})

test_that("regulon_result_v1 enforces matrix shape + dimnames", {
  cells   <- c("a", "b", "c")
  reg_ids <- c("GATA3", "PAX5")
  auc     <- matrix(seq(0.1, 0.6, length.out = 6L), nrow = 3, ncol = 2)
  out <- regulon_result_v1(cells, reg_ids, auc,
                           regulon_set_id = "mini",
                           engine_id      = "aucell_r")
  expect_true(is_regulon_result_v1(out))
  expect_identical(rownames(out$auc_matrix), cells)
  expect_identical(colnames(out$auc_matrix), reg_ids)
  expect_identical(out$schema_version, REGULON_RESULT_SCHEMA_VERSION)

  expect_error(regulon_result_v1(cells, reg_ids,
                                 auc[1:2, , drop = FALSE],
                                 regulon_set_id = "x", engine_id = "y"),
               regexp = "nrow", fixed = TRUE)
})

test_that("regulon_mean_by_group computes per-cluster mean AUC", {
  cells   <- sprintf("c_%02d", 1:6)
  reg_ids <- c("R1", "R2")
  auc <- matrix(c(
    0.1, 0.9,
    0.2, 0.8,
    0.3, 0.7,
    0.7, 0.1,
    0.8, 0.2,
    0.9, 0.3
  ), nrow = 6, byrow = TRUE,
     dimnames = list(cells, reg_ids))
  result <- regulon_result_v1(cells, reg_ids, auc,
                              regulon_set_id = "mini",
                              engine_id      = "aucell_r")
  groups <- c("A", "A", "A", "B", "B", "B")
  m <- regulon_mean_by_group(result, groups)
  expect_identical(rownames(m), c("A", "B"))
  expect_identical(colnames(m), c("R1", "R2"))
  expect_equal(m["A", "R1"], 0.2, tolerance = 1e-12)
  expect_equal(m["B", "R1"], 0.8, tolerance = 1e-12)
  expect_equal(m["A", "R2"], 0.8, tolerance = 1e-12)
  expect_equal(m["B", "R2"], 0.2, tolerance = 1e-12)

  expect_error(regulon_mean_by_group(result, groups[1:2]),
               regexp = "length mismatch", fixed = TRUE)
})

# ---- Sources --------------------------------------------------------------

test_that("REGULON_SOURCES advertises mock + dorothea entries", {
  ids <- vapply(REGULON_SOURCES(), `[[`, character(1), "id")
  expect_true("mock_pbmc"      %in% ids)
  expect_true("dorothea_human" %in% ids)
  expect_true("dorothea_mouse" %in% ids)
  expect_identical(get_regulon_source("mock_pbmc")$id, "mock_pbmc")
  expect_null(get_regulon_source("not_a_source"))
})

test_that("fetch_regulon_set('mock_pbmc') returns a 4-regulon set on the mock genes", {
  set <- fetch_regulon_set("mock_pbmc")
  expect_true(is_regulon_set(set))
  tfs <- vapply(set$regulons, `[[`, character(1), "tf")
  expect_setequal(tfs, c("GATA3", "PAX5", "SPI1", "KLF5"))
  # Every target must be a MOCK_GENES symbol so AUCell on the mock
  # dataset hits these directly.
  targets <- unlist(lapply(set$regulons, `[[`, "targets"))
  expect_true(all(targets %in% MOCK_GENES))
})

test_that("list_regulon_sources labels unavailable ones", {
  labels <- names(list_regulon_sources())
  if (!has_optional("dorothea")) {
    expect_true(any(grepl("DoRothEA \\(human, AB\\)\\s+\\(not installed\\)$",
                          labels)))
  }
  expect_true(any(labels == "Mock PBMC (built-in)"))
})

test_that(".dorothea_df_to_regulons converts a DoRothEA-shaped data.frame", {
  df <- data.frame(
    tf         = c("GATA3", "GATA3", "PAX5"),
    confidence = c("A",     "A",     "A"),
    target     = c("CD3D",  "NKG7",  "MS4A1"),
    mor        = c(1,       1,       -1),
    stringsAsFactors = FALSE
  )
  out <- .dorothea_df_to_regulons(df)
  expect_identical(length(out), 2L)
  tfs <- vapply(out, `[[`, character(1), "tf")
  expect_setequal(tfs, c("GATA3", "PAX5"))
  gata <- out[[which(tfs == "GATA3")]]
  expect_identical(gata$targets, c("CD3D", "NKG7"))
  expect_identical(gata$type, "activating")
  pax  <- out[[which(tfs == "PAX5")]]
  expect_identical(pax$type, "repressing")
})

test_that(".dorothea_df_to_regulons errors on missing columns", {
  expect_error(.dorothea_df_to_regulons(data.frame(tf = "X")),
               regexp = "target", fixed = TRUE)
})

test_that("fetch_regulon_set('dorothea_human') errors cleanly without dorothea", {
  skip_if(has_optional("dorothea"),
          "dorothea installed; install-error path skipped")
  expect_error(fetch_regulon_set("dorothea_human"),
               regexp = "dorothea", fixed = FALSE)
})

# ---- Registry + dispatcher -----------------------------------------------

test_that("REGULON_ENGINES registers aucell_r + aucell", {
  ids <- vapply(REGULON_ENGINES(), `[[`, character(1), "id")
  expect_setequal(ids, c("aucell_r", "aucell"))
  expect_identical(get_regulon_engine("aucell_r")$requires, character())
  expect_true("AUCell" %in% get_regulon_engine("aucell")$requires)
})

# ---- Engine version metadata (P1) ----------------------------------------
# Same provenance contract as the annotation engines: `version` is the
# engine implementation version, distinct from `result_schema` (the
# data-shape version). Earlier `run_regulon_engine()` relied on a silent
# fallback (REGULON_ENGINE_DEFAULT_VERSION) when an engine forgot to
# declare a version. Spec now requires an explicit one.

test_that("regulon_engine_spec requires a non-empty `version`", {
  expect_error(
    regulon_engine_spec(id = "t", name = "t", run_fn = function(...) list()),
    regexp = "version.*non-empty"
  )
  expect_error(
    regulon_engine_spec(id = "t", name = "t", run_fn = function(...) list(),
                        version = ""),
    regexp = "version.*non-empty"
  )
})

test_that("every built-in regulon engine declares an explicit version", {
  for (e in REGULON_ENGINES()) {
    expect_true(is.character(e$version), info = e$id)
    expect_true(nzchar(e$version), info = e$id)
    expect_false(identical(e$version, e$result_schema), info = e$id)
  }
})

test_that("run_regulon_engine stamps engine_version from spec", {
  ds <- mock_dataset(n_cells = 60)
  res <- run_regulon_engine("aucell_r", ds, fetch_regulon_set("mock_pbmc"))
  expect_identical(res$engine_id, "aucell_r")
  expect_identical(res$engine_version, "aucell_r_v1.0.0")
  # engine_version is not the schema label
  expect_false(identical(res$engine_version, res$schema_version))
})

test_that("list_regulon_engines tags unavailable backends", {
  labels <- names(list_regulon_engines())
  # Pure-R is always available
  expect_true(any(grepl("^AUCell \\(pure R\\)$", labels)))
  # Bioc is unavailable in CI
  if (!has_optional("AUCell")) {
    expect_true(any(grepl("AUCell \\(Bioconductor\\)\\s+\\(not installed\\)$",
                          labels)))
  }
})

test_that("run_regulon_engine errors on unknown engine and on wrong inputs", {
  ds <- mock_dataset(n_cells = 30)
  expect_error(run_regulon_engine("nope", ds, fetch_regulon_set("mock_pbmc")),
               regexp = "Unknown regulon engine", fixed = TRUE)
  expect_error(run_regulon_engine("aucell_r", ds, regulon_set = list()),
               regexp = "not a regulon_set", fixed = TRUE)
  expect_error(run_regulon_engine("aucell_r", NULL, fetch_regulon_set("mock_pbmc")),
               regexp = "No dataset", fixed = TRUE)
})

test_that("run_regulon_engine('aucell') errors cleanly without AUCell", {
  skip_if(has_optional("AUCell"),
          "AUCell installed; install-error path skipped")
  ds <- mock_dataset(n_cells = 30)
  expect_error(run_regulon_engine("aucell", ds,
                                  fetch_regulon_set("mock_pbmc")),
               regexp = "AUCell", fixed = FALSE)
})

# ---- Pure AUCell math ----------------------------------------------------

test_that(".aucell_pure_r is monotone in target rank and in [0, 1]", {
  # Build a small dataset where regulon target genes are at known ranks
  # in two specific cells; verify the AUC is higher where ranks are
  # more concentrated at the top.
  set.seed(42L)
  genes <- sprintf("G%03d", 1:100)
  # cell_top: targets at ranks 1, 2 -> max AUC
  # cell_mid: targets at ranks 4, 5
  # cell_low: targets at ranks 50, 60
  expr <- matrix(0, nrow = 100, ncol = 3,
                 dimnames = list(genes, c("cell_top", "cell_mid", "cell_low")))
  expr["G001", "cell_top"] <- 10
  expr["G002", "cell_top"] <- 9
  expr["G001", "cell_mid"] <- 4
  expr["G002", "cell_mid"] <- 5
  expr["G050", "cell_top"] <- 8
  expr["G060", "cell_top"] <- 7
  expr["G050", "cell_mid"] <- 3
  expr["G060", "cell_mid"] <- 2
  # for cell_low: targets at low ranks => essentially zero AUC
  expr["G001", "cell_low"] <- 0
  expr["G002", "cell_low"] <- 0
  # Fill other genes with descending baseline so ranks are well-defined
  for (c in colnames(expr)) {
    others <- setdiff(genes, c("G001", "G002"))
    expr[others, c] <- expr[others, c] + seq_along(others) * 0.001
  }

  out <- .aucell_pure_r(
    expr_mat = expr,
    regulons = list(R1 = c("G001", "G002")),
    top_n_fraction = 0.10  # top 10 of 100 -> ranks 1..10 count
  )
  auc <- out$auc_matrix[, "R1"]
  expect_identical(length(auc), 3L)
  expect_true(all(auc >= 0 - 1e-9 & auc <= 1 + 1e-9))
  # cell_top must rank R1 highest
  expect_equal(unname(auc["cell_top"]), max(auc))
})

test_that(".aucell_pure_r emits a warning when no targets are in the dataset", {
  genes <- c("A", "B", "C")
  expr <- matrix(c(1, 2, 3, 2, 1, 0), nrow = 3, ncol = 2,
                 dimnames = list(genes, c("c1", "c2")))
  out <- .aucell_pure_r(expr, list(R1 = c("X", "Y")), top_n_fraction = 0.5)
  expect_true(any(grepl("no target genes present", out$warnings)))
  expect_true(all(out$auc_matrix[, "R1"] == 0))
})

test_that(".aucell_pure_r rejects bad top_n_fraction", {
  expr <- matrix(1, nrow = 3, ncol = 3,
                 dimnames = list(c("A", "B", "C"), c("c1", "c2", "c3")))
  expect_error(.run_aucell_pure_r_regulons(
    list(expression = expression_backend_inmemory(list(data = expr)),
         cell_data  = data.frame(cell = c("c1", "c2", "c3"))),
    regulon_set("x", "x", regulons = list(regulon_spec("R", "A"))),
    params = list(top_n_fraction = 0)),
    regexp = "top_n_fraction", fixed = TRUE)
})

# ---- Pure converter ------------------------------------------------------

test_that(".aucell_to_regulon_engine_output aligns rows by cell order", {
  cells <- c("c1", "c2", "c3")
  # rows out of order
  auc <- matrix(c(0.3, 0.1, 0.2,
                  0.9, 0.7, 0.8), nrow = 3, ncol = 2,
                 dimnames = list(c("c3", "c1", "c2"),
                                 c("R1", "R2")))
  out <- .aucell_to_regulon_engine_output(
    auc_matrix = auc, cells = cells,
    regulon_ids = c("R1", "R2"))
  expect_identical(rownames(out$auc_matrix), cells)
  expect_equal(unname(out$auc_matrix[, "R1"]), c(0.1, 0.2, 0.3),
               tolerance = 1e-12)
})

test_that(".aucell_to_regulon_engine_output errors on missing cells", {
  auc <- matrix(0.5, nrow = 1, ncol = 1,
                dimnames = list("c1", "R1"))
  expect_error(.aucell_to_regulon_engine_output(
    auc_matrix = auc, cells = c("c1", "c2"),
    regulon_ids = "R1"),
    regexp = "missing from AUC matrix", fixed = TRUE)
})

# ---- End-to-end on the mock dataset --------------------------------------

test_that("aucell_r round-trips against the mock dataset + mock_pbmc regulons", {
  ds  <- mock_dataset(n_cells = 200, seed = 7)
  set <- fetch_regulon_set("mock_pbmc")
  result <- run_regulon_engine("aucell_r", ds, set,
                               params = list(top_n_fraction = 0.05))
  expect_true(is_regulon_result_v1(result))
  expect_identical(length(result$cell_ids), 200L)
  expect_setequal(result$regulon_ids, c("GATA3", "PAX5", "SPI1", "KLF5"))
  expect_true(all(result$auc_matrix >= 0 - 1e-9 &
                  result$auc_matrix <= 1 + 1e-9))

  # MOCK_GENES is structured so each regulon should peak in the
  # corresponding cluster (see .mock_pbmc_regulon_set + mock_dataset
  # gene_cluster mapping):
  #   GATA3 -> CD3D, NKG7 -> cluster 0 (T cell)
  #   PAX5  -> MS4A1      -> cluster 1 (B cell)
  #   SPI1  -> LST1       -> cluster 2 (Myeloid)
  #   KLF5  -> EPCAM/COL1A1 -> cluster 3 (Epithelial)
  by_cluster <- regulon_mean_by_group(result, ds$cell_data$cluster)
  expect_identical(rownames(by_cluster)[which.max(by_cluster[, "GATA3"])], "0")
  expect_identical(rownames(by_cluster)[which.max(by_cluster[, "PAX5"])],  "1")
  expect_identical(rownames(by_cluster)[which.max(by_cluster[, "SPI1"])],  "2")
  expect_identical(rownames(by_cluster)[which.max(by_cluster[, "KLF5"])],  "3")
})

test_that("aucell (Bioc) round-trips against the mock dataset", {
  skip_if_not_installed("AUCell")
  ds  <- mock_dataset(n_cells = 200, seed = 7)
  set <- fetch_regulon_set("mock_pbmc")
  result <- run_regulon_engine("aucell", ds, set,
                               params = list(top_n_fraction = 0.05))
  expect_true(is_regulon_result_v1(result))
  expect_identical(length(result$cell_ids), 200L)
  expect_setequal(result$regulon_ids, c("GATA3", "PAX5", "SPI1", "KLF5"))
})

# ---- Module registry hookup ----------------------------------------------

test_that("the Regulons module is registered + enabled (replaces `regulatory`)", {
  registry <- module_registry()
  ids <- vapply(registry, `[[`, character(1), "id")
  expect_true("regulons"   %in% ids)
  expect_false("regulatory" %in% ids)  # old placeholder retired
  m <- get_module("regulons")
  expect_true(m$enabled)
  expect_true(is.function(m$ui_fn))
  expect_true(is.function(m$server_fn))
})
