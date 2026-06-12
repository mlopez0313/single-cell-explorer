test_that("default_marker_registry has the documented top-level shape", {
  reg <- default_marker_registry()
  expect_identical(reg$schema_version, "marker_registry_v1")
  expect_true(nzchar(reg$version))
  expect_identical(reg$source, "builtin")
  expect_s3_class(reg$created_at, "POSIXct")
  expect_true(length(reg$entries) >= 6L)
})

test_that("every marker_entry reserves the ontology_id slot (NA allowed)", {
  reg <- default_marker_registry()
  for (e in reg$entries) {
    expect_true("ontology_id" %in% names(e))
    expect_type(e$ontology_id, "character")
  }
})

test_that("every marker_gene has role / weight / evidence", {
  reg <- default_marker_registry()
  for (e in reg$entries) {
    for (m in e$markers) {
      expect_true(all(c("gene", "role", "weight", "evidence") %in% names(m)))
      expect_true(m$role %in% c("positive", "negative", "specific"))
      expect_type(m$weight, "double")
      expect_true(m$evidence %in% c("protein_validated", "scRNA_only", "inferred"))
    }
  }
})

test_that("marker_gene rejects invalid role / evidence values", {
  expect_error(marker_gene("CD3D", role = "magic"))
  expect_error(marker_gene("CD3D", evidence = "guessing"))
})

test_that("marker_registry_cell_types lists all registered cell types", {
  reg <- default_marker_registry()
  ct <- marker_registry_cell_types(reg)
  expect_true(all(c("T cell", "B cell", "NK cell", "Myeloid cell",
                    "Epithelial cell", "Fibroblast") %in% ct))
})

test_that("marker_registry_filter honours species + cell_type + gene filters", {
  reg <- default_marker_registry()
  expect_true(length(marker_registry_filter(reg, species = "human")) >= 6L)
  # NULL species accepts everything (including NA species)
  expect_identical(length(marker_registry_filter(reg)),
                   length(reg$entries))
  cd3 <- marker_registry_filter(reg, gene = "CD3D")
  expect_true(any(vapply(cd3, function(e) e$cell_type == "T cell", logical(1))))
  bt <- marker_registry_filter(reg, cell_type = "B cell")
  expect_length(bt, 1L)
})

test_that("marker_registry_get returns NULL for unknown labels", {
  reg <- default_marker_registry()
  expect_identical(marker_registry_get(reg, "T cell")$cell_type, "T cell")
  expect_null(marker_registry_get(reg, "Banana"))
})

test_that("marker_registry_genes returns a union (or per-cell-type subset)", {
  reg <- default_marker_registry()
  all_g <- marker_registry_genes(reg)
  expect_true(all(c("CD3D", "MS4A1", "NKG7", "EPCAM") %in% all_g))
  t_g <- marker_registry_genes(reg, "T cell")
  expect_true("CD3D" %in% t_g)
  expect_false("EPCAM" %in% t_g)
})

test_that("build_ontology_map fills known IDs and NA for unknowns", {
  reg <- default_marker_registry()
  om <- build_ontology_map(reg, c("T cell", "B cell", "Banana"))
  expect_identical(unname(om[["T cell"]]), "CL:0000084")
  expect_identical(unname(om[["B cell"]]), "CL:0000236")
  expect_true(is.na(om[["Banana"]]))
})
