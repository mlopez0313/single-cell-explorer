test_that(".ensure_bioc_cache_dirs creates the configured cache directory", {
  skip_if_not_installed("withr")
  td <- tempfile("sce-user-cache-")
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  eh_path <- file.path(td, "EH-cache")
  ah_path <- file.path(td, "AH-cache")
  withr::with_envvar(c(R_USER_CACHE_DIR     = td,
                       EXPERIMENT_HUB_CACHE = eh_path,
                       ANNOTATION_HUB_CACHE = ah_path), {
    .ensure_bioc_cache_dirs()
  })
  # The helper must have created at least one of the candidate paths:
  # the env-var-supplied one (ExperimentHub installed) or the
  # R_user_dir fallback (not installed). Both are valid -- what matters
  # is that *something* now exists so a follow-up readline() prompt
  # would not fire.
  rud_eh <- withr::with_envvar(c(R_USER_CACHE_DIR = td), {
    tools::R_user_dir("ExperimentHub", which = "cache")
  })
  expect_true(dir.exists(eh_path) || dir.exists(rud_eh))
})

test_that(".ensure_bioc_cache_dirs is idempotent on a populated cache", {
  skip_if_not_installed("withr")
  td <- tempfile("sce-user-cache-")
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  withr::with_envvar(c(R_USER_CACHE_DIR = td), {
    .ensure_bioc_cache_dirs()
    # A second call must not touch the FS in a way that surfaces as a
    # warning or error.
    expect_silent(.ensure_bioc_cache_dirs())
  })
})

test_that(".ensure_bioc_cache_dirs swallows mkdir failures (best-effort)", {
  skip_if_not_installed("withr")
  # Point ExperimentHub at a path nested under a file (not a dir), so
  # `dir.create()` cannot succeed. The helper must still return cleanly.
  blocker <- tempfile("sce-blocker-")
  writeLines("not a dir", blocker)
  on.exit(unlink(blocker), add = TRUE)
  bad_path <- file.path(blocker, "EH-cache")
  withr::with_envvar(c(EXPERIMENT_HUB_CACHE = bad_path), {
    expect_silent(.ensure_bioc_cache_dirs())
  })
})

test_that("demo_dataset_path honours SCE_DEMO_DATASET when set", {
  tmp <- tempfile(fileext = ".rds")
  withr::with_envvar(c(SCE_DEMO_DATASET = tmp), {
    expect_identical(demo_dataset_path(),
                     normalizePath(tmp, mustWork = FALSE))
  })
})

test_that("demo_dataset_path falls back to canonical inst/extdata path", {
  withr::with_envvar(c(SCE_DEMO_DATASET = ""), {
    tmp <- tempfile("sce_proj_")
    dir.create(tmp)
    p <- demo_dataset_path(project_root = tmp)
    expect_match(p, "inst/extdata/pbmc8k_demo\\.rds$")
    expect_false(demo_dataset_exists(project_root = tmp))
  })
})

test_that(".find_project_root walks up to a DESCRIPTION + app.R pair", {
  # Create a fake project: <tmp>/proj/{DESCRIPTION, app.R, R/, scripts/}
  root <- tempfile("sce_root_"); dir.create(root)
  proj <- file.path(root, "proj"); dir.create(proj)
  dir.create(file.path(proj, "R"))
  dir.create(file.path(proj, "scripts"))
  writeLines("Package: x", file.path(proj, "DESCRIPTION"))
  writeLines("# app",      file.path(proj, "app.R"))

  # From the scripts/ subdir -> walks up to proj/.
  expect_identical(normalizePath(.find_project_root(file.path(proj, "scripts")),
                                 mustWork = FALSE),
                   normalizePath(proj, mustWork = FALSE))
  # From a deeper subdir.
  deep <- file.path(proj, "R", "modules", "deep"); dir.create(deep, recursive = TRUE)
  expect_identical(normalizePath(.find_project_root(deep), mustWork = FALSE),
                   normalizePath(proj, mustWork = FALSE))
  # When no project is found above `start`, returns `start` unchanged.
  outside <- tempfile("sce_outside_"); dir.create(outside)
  expect_identical(normalizePath(.find_project_root(outside), mustWork = FALSE),
                   normalizePath(outside, mustWork = FALSE))
})

test_that("demo_dataset_path walks up from cwd to find the real project root", {
  # Drop into a deeper subdir of the real project; the resolver should
  # still land on `<project_root>/inst/extdata/pbmc8k_demo.rds`.
  withr::with_envvar(c(SCE_DEMO_DATASET = ""), {
    proj_root <- normalizePath(file.path(test_path(), "..", ".."),
                               mustWork = TRUE)
    withr::with_dir(file.path(proj_root, "R"), {
      p <- demo_dataset_path()
      expect_match(p, "inst/extdata/pbmc8k_demo\\.rds$")
      expect_match(normalizePath(p, mustWork = FALSE),
                   normalizePath(proj_root, mustWork = FALSE),
                   fixed = TRUE)
    })
  })
})

test_that("demo_dataset_exists is TRUE when an artifact is on disk", {
  withr::with_envvar(c(SCE_DEMO_DATASET = ""), {
    tmp <- tempfile("sce_proj_")
    dir.create(file.path(tmp, "inst", "extdata"), recursive = TRUE)
    art <- file.path(tmp, "inst", "extdata", "pbmc8k_demo.rds")
    saveRDS(list(name = "x"), art)   # contents don't matter for this check
    expect_true(demo_dataset_exists(project_root = tmp))
    expect_identical(demo_dataset_path(project_root = tmp),
                     normalizePath(art, mustWork = FALSE))
  })
})

test_that("load_demo_dataset errors clearly when the artifact is missing", {
  err <- tryCatch(load_demo_dataset(path = tempfile(fileext = ".rds")),
                  error = conditionMessage)
  expect_match(err, "demo dataset not found", ignore.case = TRUE)
  expect_match(err, "scripts/build_pbmc8k_demo.R", fixed = TRUE)
})

test_that("load_demo_dataset round-trips a schema-valid dataset", {
  ds  <- mock_dataset(n_cells = 80, name = "rt_demo")
  out <- tempfile(fileext = ".rds")
  saveRDS(ds, out)
  loaded <- load_demo_dataset(path = out)
  expect_true(all(dataset_schema() %in% names(loaded)))
  expect_identical(loaded$n_cells, ds$n_cells)
  expect_s3_class(loaded$expression, "expression_backend")
})

test_that(".validate_demo_dataset rejects broken artifacts", {
  ds <- mock_dataset(n_cells = 40, name = "broken")

  # Missing required field.
  bad1 <- ds; bad1$expression <- NULL
  expect_error(.validate_demo_dataset(bad1),
               "missing schema field|expression_backend",
               ignore.case = TRUE)

  # cell_data row count mismatch.
  bad2 <- ds; bad2$cell_data <- bad2$cell_data[1:5, , drop = FALSE]
  expect_error(.validate_demo_dataset(bad2),
               "length\\(cells\\) = .* nrow\\(cell_data\\)")

  # Missing embedding column for an advertised reduction.
  bad3 <- ds; bad3$cell_data$UMAP_1 <- NULL
  expect_error(.validate_demo_dataset(bad3),
               "reduction 'UMAP' missing column", fixed = TRUE)
})

test_that(".augment_demo_dataset adds the demo-friendly fields", {
  ds <- mock_dataset(n_cells = 120, name = "augment_in")
  # Pretend this is a freshly-Seurat-converted dataset by stripping the
  # demo metadata columns the mock already carries.
  ds$cell_data$cluster         <- NULL
  ds$cell_data$cell_type       <- NULL
  ds$cell_data$condition       <- NULL
  ds$cell_data$pseudotime_demo <- NULL
  ds$cell_data$sample          <- NULL
  # Give it a Seurat-style cluster column so the function can promote it.
  ds$cell_data$seurat_clusters <- sample(c("0", "1", "2"), 120, replace = TRUE)
  ds$metadata_fields <- "seurat_clusters"

  out <- .augment_demo_dataset(ds, seed = 8L)
  cd  <- out$cell_data
  expect_true(all(c("cluster", "cell_type", "condition", "sample",
                    "pseudotime_demo") %in% names(cd)))
  expect_identical(cd$cluster, as.character(ds$cell_data$seurat_clusters))
  expect_true(all(cd$pseudotime_demo >= 0 - 1e-9))
  expect_true(all(cd$pseudotime_demo <= 1 + 1e-9))
  # Default metadata field lands on "sample" -- matches the mock dataset.
  expect_identical(out$metadata_fields[1], "sample")
  expect_identical(out$source, "demo_pbmc8k")
})

test_that("demo_auto_build_enabled honours SCE_AUTO_BUILD_DEMO", {
  withr::with_envvar(c(SCE_AUTO_BUILD_DEMO = ""),       expect_true(demo_auto_build_enabled()))
  withr::with_envvar(c(SCE_AUTO_BUILD_DEMO = "1"),      expect_true(demo_auto_build_enabled()))
  withr::with_envvar(c(SCE_AUTO_BUILD_DEMO = "true"),   expect_true(demo_auto_build_enabled()))
  withr::with_envvar(c(SCE_AUTO_BUILD_DEMO = "0"),      expect_false(demo_auto_build_enabled()))
  withr::with_envvar(c(SCE_AUTO_BUILD_DEMO = "false"),  expect_false(demo_auto_build_enabled()))
  withr::with_envvar(c(SCE_AUTO_BUILD_DEMO = "no"),     expect_false(demo_auto_build_enabled()))
  withr::with_envvar(c(SCE_AUTO_BUILD_DEMO = "off"),    expect_false(demo_auto_build_enabled()))
})

test_that("can_build_demo_dataset / demo_auto_build_status agree", {
  # Either every required package is present (status reports "all
  # present") or at least one is missing and the status names it.
  status <- demo_auto_build_status()
  if (can_build_demo_dataset()) {
    expect_match(status, "all auto-build packages present", fixed = TRUE)
  } else {
    expect_match(status, "missing package")
    # Whatever names the status lists must currently fail `has_optional()`.
    pkgs <- strsplit(sub("^missing package\\(s\\): ", "", status), ", ")[[1]]
    for (p in pkgs) expect_false(has_optional(p))
  }
})

test_that("ensure_demo_dataset short-circuits on an existing artifact", {
  withr::with_envvar(c(SCE_DEMO_DATASET = "", SCE_AUTO_BUILD_DEMO = "1"), {
    tmp <- tempfile(fileext = ".rds")
    saveRDS(mock_dataset(n_cells = 32, name = "ensure_rt"), tmp)
    withr::with_envvar(c(SCE_DEMO_DATASET = tmp), {
      ds <- ensure_demo_dataset()
      expect_identical(ds$n_cells, 32L)
    })
  })
})

test_that("ensure_demo_dataset errors when auto-build is off and no artifact", {
  withr::with_envvar(c(SCE_DEMO_DATASET = tempfile(fileext = ".rds"),
                       SCE_AUTO_BUILD_DEMO = "0"), {
    err <- tryCatch(ensure_demo_dataset(), error = conditionMessage)
    expect_match(err, "demo dataset not found", ignore.case = TRUE)
  })
})

test_that("ensure_demo_dataset errors with a clear hint when deps are missing", {
  skip_if(can_build_demo_dataset(),
          "All build packages installed; missing-deps path not exercised here")
  withr::with_envvar(c(SCE_DEMO_DATASET = tempfile(fileext = ".rds"),
                       SCE_AUTO_BUILD_DEMO = "1"), {
    err <- tryCatch(ensure_demo_dataset(), error = conditionMessage)
    expect_match(err, "Cannot auto-build", fixed = TRUE)
    expect_match(err, "missing package", fixed = TRUE)
    expect_match(err, "Rscript scripts/build_pbmc8k_demo.R", fixed = TRUE)
  })
})

test_that(".build_progress_handle swallows callback errors", {
  bad <- .build_progress_handle(function(...) stop("boom"))
  # Must not propagate the error.
  expect_silent(bad(0.5, "x"))
  # And the no-op handle is silent too.
  expect_silent(.build_progress_handle(NULL)(0.5))
})

test_that(".augment_demo_dataset is reproducible for a given seed", {
  ds <- mock_dataset(n_cells = 60, name = "seed_test")
  ds$cell_data$condition       <- NULL
  ds$cell_data$sample          <- NULL
  ds$cell_data$pseudotime_demo <- NULL
  ds$cell_data$cluster         <- NULL
  ds$cell_data$cell_type       <- NULL
  ds$cell_data$seurat_clusters <- rep(c("0", "1", "2"), length.out = 60)

  a <- .augment_demo_dataset(ds, seed = 1L)
  b <- .augment_demo_dataset(ds, seed = 1L)
  expect_identical(a$cell_data$condition, b$cell_data$condition)
  expect_identical(a$cell_data$sample,    b$cell_data$sample)
})
