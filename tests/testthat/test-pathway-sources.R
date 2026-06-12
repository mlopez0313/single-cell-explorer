test_that("PATHWAY_SOURCES exposes builtin and msigdbr sources", {
  srcs <- PATHWAY_SOURCES()
  ids  <- vapply(srcs, `[[`, character(1), "id")
  expect_true("builtin" %in% ids)
  expect_true("msigdbr" %in% ids)
  # builtin must be listed first so `available_pathway_collections()[1]`
  # remains a stable demo default regardless of which optional packages
  # the user has installed.
  expect_identical(ids[1], "builtin")
})

test_that("get_pathway_source returns spec by id and NULL on miss", {
  expect_true(is.list(get_pathway_source("builtin")))
  expect_identical(get_pathway_source("builtin")$id, "builtin")
  expect_null(get_pathway_source("nope_not_a_source"))
})

test_that("available_pathway_collections still exposes mock_v1 (back-compat)", {
  cols <- available_pathway_collections()
  expect_true("mock_v1" %in% cols)
  # mock_v1 must be the first element so existing default-selection logic
  # (`available_pathway_collections()[1]`) does not regress.
  expect_identical(cols[1], "mock_v1")
})

test_that("get_pathways(mock_v1) is unchanged by the registry refactor", {
  pw <- get_pathways("mock_v1")
  expect_type(pw, "list")
  expect_true(length(pw) > 0L)
  expect_true("T cell activation" %in% names(pw))
  # Spot-check that the gene sets themselves are intact.
  expect_true("CD3D" %in% pw[["T cell activation"]])
})

test_that("get_pathways accepts blank ids (reactive UI) but errors on unknown ids", {
  # Blank / NULL keep returning NULL so reactive callers can early-out
  # before a selection has been made.
  expect_null(get_pathways(""))
  expect_null(get_pathways(NULL))
  # Non-blank id that no registered source owns is a misuse: surface
  # it loudly instead of returning an empty result silently. Real-data
  # smoke testing turned up "compute_enrichment(... pathways=NULL) ->
  # 0 results" as a recurring confusion source.
  expect_error(get_pathways("not_a_real_collection"),
               regexp = "Unknown pathway collection")
})

test_that("pathway_collection_info reports source provenance", {
  info <- pathway_collection_info("mock_v1")
  expect_identical(info$source_id, "builtin")
  expect_true(isTRUE(info$available))
  # Unknown collection -> empty list, not an error.
  expect_identical(pathway_collection_info("does_not_exist"), list())
})

test_that("msigdbr collections are listed only when msigdbr is installed", {
  cols <- available_pathway_collections()
  has_msigdbr <- has_optional("msigdbr")
  if (has_msigdbr) {
    expect_true(any(startsWith(cols, "msigdbr/")))
  } else {
    expect_false(any(startsWith(cols, "msigdbr/")))
  }
})

test_that("get_pathways(msigdbr/...) errors cleanly when msigdbr is missing", {
  skip_if(has_optional("msigdbr"), "msigdbr is installed; install-error path skipped")
  # Force resolution by asking the source's fetcher directly, since
  # `.resolve_collection_source` would skip the source when its package
  # is missing. This exercises the exact `require_optional()` call that
  # would fire from a future UI/path that bypasses the registry check.
  src <- get_pathway_source("msigdbr")
  expect_error(src$fetcher("msigdbr/H"), regexp = "msigdbr", fixed = FALSE)
})

# ----- Pure converter -------------------------------------------------------

test_that(".msigdbr_to_pathways converts a msigdbr-shaped frame", {
  df <- data.frame(
    gs_name     = c("HALLMARK_A", "HALLMARK_A", "HALLMARK_A",
                    "HALLMARK_B", "HALLMARK_B"),
    gene_symbol = c("CD3D", "CD3E", "CD3G",
                    "MS4A1", "CD19"),
    stringsAsFactors = FALSE
  )
  out <- .msigdbr_to_pathways(df)
  expect_type(out, "list")
  expect_identical(names(out), c("HALLMARK_A", "HALLMARK_B"))
  expect_identical(out[["HALLMARK_A"]], c("CD3D", "CD3E", "CD3G"))
  expect_identical(out[["HALLMARK_B"]], c("MS4A1", "CD19"))
})

test_that(".msigdbr_to_pathways de-duplicates genes within a pathway", {
  df <- data.frame(
    gs_name     = rep("P", 4),
    gene_symbol = c("A", "B", "A", "C"),
    stringsAsFactors = FALSE
  )
  expect_identical(.msigdbr_to_pathways(df)[["P"]], c("A", "B", "C"))
})

test_that(".msigdbr_to_pathways supports human_gene_symbol fallback column", {
  df <- data.frame(
    gs_name           = c("P1", "P1", "P2"),
    human_gene_symbol = c("X", "Y", "Z"),
    stringsAsFactors  = FALSE
  )
  out <- .msigdbr_to_pathways(df)
  expect_identical(out[["P1"]], c("X", "Y"))
  expect_identical(out[["P2"]], "Z")
})

test_that(".msigdbr_to_pathways errors clearly on missing required columns", {
  expect_error(.msigdbr_to_pathways(data.frame(foo = "bar")),
               regexp = "gs_name", fixed = TRUE)
  expect_error(.msigdbr_to_pathways(data.frame(gs_name = "P", foo = "X")),
               regexp = "gene-symbol", fixed = TRUE)
})

test_that(".msigdbr_to_pathways handles empty / NULL input", {
  expect_identical(.msigdbr_to_pathways(NULL), list())
  expect_identical(
    .msigdbr_to_pathways(data.frame(gs_name = character(),
                                    gene_symbol = character(),
                                    stringsAsFactors = FALSE)),
    list()
  )
})

# ----- End-to-end (skipped without msigdbr) --------------------------------

test_that("msigdbr end-to-end fetch (requires installed package)", {
  skip_if_not_installed("msigdbr")
  # `msigdbdf` is the data backend used by recent msigdbr versions.
  # Skipping when it's missing avoids a large surprise download in CI.
  skip_if_not_installed("msigdbdf")
  pw <- get_pathways("msigdbr/H")
  expect_type(pw, "list")
  expect_true(length(pw) >= 40L) # ~50 Hallmark sets in MSigDB
  expect_true(all(vapply(pw, length, integer(1)) > 0L))
})
