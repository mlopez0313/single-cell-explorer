# Tests for R/setup.R: tiers, status, preflight message, dry-run install.

test_that(".SCE_PKG_TIERS has the expected shape", {
  expect_true(all(c("core", "demo", "full") %in% names(.SCE_PKG_TIERS)))
  for (tier in names(.SCE_PKG_TIERS)) {
    spec <- .SCE_PKG_TIERS[[tier]]
    expect_true(all(c("cran", "bioc", "desc") %in% names(spec)))
    expect_type(spec$cran, "character")
    expect_type(spec$bioc, "character")
    expect_type(spec$desc, "character")
  }
})

test_that(".sce_packages_for_tier is cumulative across tiers", {
  core <- .sce_packages_for_tier("core")
  demo <- .sce_packages_for_tier("demo")
  full <- .sce_packages_for_tier("full")
  # Cumulative: core ⊂ demo ⊂ full (when expressed as character sets).
  expect_true(all(c(core$cran, core$bioc) %in% c(demo$cran, demo$bioc)))
  expect_true(all(c(demo$cran, demo$bioc) %in% c(full$cran, full$bioc)))
  # Demo specifically includes the PBMC 8k build packages.
  expect_true("TENxPBMCData" %in% demo$bioc)
  expect_true("Seurat"       %in% demo$cran)
})

test_that("sce_check_setup returns one structured entry per tier", {
  st <- sce_check_setup()
  expect_identical(sort(names(st)), sort(names(.SCE_PKG_TIERS)))
  for (tier in names(st)) {
    s <- st[[tier]]
    expect_true(all(c("desc", "cran", "bioc", "missing", "present",
                      "total", "complete") %in% names(s)))
    expect_type(s$total,    "integer")
    expect_type(s$present,  "integer")
    expect_type(s$missing,  "character")
    expect_type(s$complete, "logical")
    expect_identical(s$present + length(s$missing), s$total)
  }
})

test_that("sce_preflight_message lists OK / missing per tier", {
  msg <- sce_preflight_message()
  expect_match(msg, "scRNA Explorer dependency preflight", fixed = TRUE)
  # Every tier appears in the message.
  for (tier in names(.SCE_PKG_TIERS)) expect_match(msg, tier, fixed = TRUE)
  # If anything is missing on this machine, the message points at the
  # setup CLI.
  st <- sce_check_setup()
  if (!all(vapply(st, `[[`, logical(1), "complete"))) {
    expect_match(msg, "scripts/setup_dev.R", fixed = TRUE)
  }
})

test_that("sce_setup(dry_run = TRUE) reports the missing-packages plan", {
  # We don't actually install anything; dry_run = TRUE returns the same
  # missing-packages plan that a real run would attempt.
  out <- capture.output(res <- sce_setup(tier = "demo", auto = TRUE,
                                          dry_run = TRUE))
  expect_match(paste(out, collapse = "\n"),
               "scRNA Explorer setup -- tier 'demo'", fixed = TRUE)
  expect_match(paste(out, collapse = "\n"),
               "dry-run", fixed = TRUE)
  expect_type(res, "list")
  expect_true(all(c("cran", "bioc") %in% names(res)))
  # Every name in `res$cran` / `res$bioc` should currently NOT be
  # installed (that's the meaning of a "missing-packages plan").
  for (p in c(res$cran, res$bioc)) expect_false(has_optional(p))
})

test_that("sce_setup refuses to install silently in a non-interactive session", {
  # The test runner is non-interactive (`interactive() == FALSE`). When
  # `auto = FALSE` we should refuse rather than silently `install.packages()`.
  # We only check this when there's actually something to install; if the
  # current machine has everything in tier 'core' already, sce_setup()
  # short-circuits before the interactive check.
  st <- sce_check_setup()
  if (!st$core$complete) {
    expect_error(sce_setup(tier = "core", auto = FALSE, dry_run = FALSE),
                 "non-interactive", fixed = TRUE)
  } else {
    succeed("tier 'core' already complete on this machine; skipped")
  }
})
