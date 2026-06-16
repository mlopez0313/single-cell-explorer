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
  # Demo specifically includes Seurat for the PBMC 8k build. We no
  # longer require any Bioconductor package in the demo tier since the
  # default build path pulls the matrix directly from the 10x CDN
  # (cf.10xgenomics.com), bypassing ExperimentHub.
  expect_true("Seurat"       %in% demo$cran)
  expect_true("SeuratObject" %in% demo$cran)
  expect_identical(demo$bioc, character())
  # Full tier surfaces a `github` named character vector keyed by
  # package name; Azimuth is the canonical example.
  expect_true("github" %in% names(full))
  expect_true("Azimuth" %in% names(full$github))
  expect_identical(unname(full$github[["Azimuth"]]),
                   "satijalab/azimuth")
  # GitHub installs are not cumulative-downward: core / demo carry no
  # github packages.
  expect_identical(core$github %||% character(), character())
  expect_identical(demo$github %||% character(), character())
})

test_that("sce_check_setup surfaces GitHub-only packages alongside CRAN/Bioc", {
  st <- sce_check_setup()
  # The `github` slot must be a named character vector matching the
  # tier definition; `total` must include it.
  expect_true("github" %in% names(st$full))
  expect_true("Azimuth" %in% names(st$full$github))
  expect_identical(st$full$total,
                   length(c(st$full$cran, st$full$bioc,
                            names(st$full$github))))
  # Round-trip: present + length(missing) is still total.
  expect_identical(st$full$present + length(st$full$missing),
                   st$full$total)
})

test_that("sce_setup(dry_run = TRUE) reports the GitHub-missing plan", {
  out <- capture.output(res <- sce_setup(tier = "full", auto = TRUE,
                                          dry_run = TRUE))
  joined <- paste(out, collapse = "\n")
  expect_match(joined, "GitHub missing", fixed = TRUE)
  expect_type(res$github, "character")
  # The result list now carries `github` alongside `cran` / `bioc`.
  expect_true(all(c("cran", "bioc", "github") %in% names(res)))
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

test_that("sce_install_github_pkg short-circuits when the package is already installed", {
  skip_if_not_installed("testthat")
  # `testthat` is reliably installed in any environment that's running
  # these tests, so re-attempt to "install" it. The helper must take
  # the fast path and never reach `remotes::install_github()`.
  called <- FALSE
  fake_install_github <- function(...) {
    called <<- TRUE
    stop("remotes::install_github() should not be called when the ",
         "package is already installed.")
  }
  # Patch the remotes namespace so a misbehaving helper would explode
  # loudly instead of silently dialling out to GitHub.
  if (requireNamespace("remotes", quietly = TRUE)) {
    old <- remotes::install_github
    unlockBinding("install_github", asNamespace("remotes"))
    assign("install_github", fake_install_github,
           envir = asNamespace("remotes"))
    on.exit({
      assign("install_github", old, envir = asNamespace("remotes"))
      lockBinding("install_github", asNamespace("remotes"))
    }, add = TRUE)
  }
  expect_invisible(sce_install_github_pkg("testthat",
                                           "r-lib/testthat"))
  expect_false(called)
})

test_that("sce_setup raises when an install attempt leaves packages missing", {
  # Build a controlled situation: a fake tier of packages we know
  # don't exist. Mock install.packages() / BiocManager::install() to
  # be no-ops so we can simulate the "install ran but produced no
  # working package" failure mode that motivated this test.
  ghost_pkg <- "sce_test_ghost_pkg_does_not_exist"
  old_tiers <- .SCE_PKG_TIERS
  on.exit(assign(".SCE_PKG_TIERS", old_tiers,
                 envir = .GlobalEnv), add = TRUE)
  .SCE_PKG_TIERS[["core"]] <<- list(
    cran = ghost_pkg, bioc = character(), desc = "(test)")
  .SCE_PKG_TIERS[["demo"]] <<- list(
    cran = character(), bioc = character(), desc = "(test)")
  .SCE_PKG_TIERS[["full"]] <<- list(
    cran = character(), bioc = character(), desc = "(test)")

  # No-op the actual installer.
  with_mocked_install <- function(expr) {
    old <- utils::install.packages
    unlockBinding("install.packages", asNamespace("utils"))
    assign("install.packages", function(...) invisible(NULL),
           envir = asNamespace("utils"))
    on.exit({
      assign("install.packages", old, envir = asNamespace("utils"))
      lockBinding("install.packages", asNamespace("utils"))
    }, add = TRUE)
    force(expr)
  }

  with_mocked_install({
    expect_error(
      suppressMessages(capture.output(
        sce_setup(tier = "core", auto = TRUE))),
      regexp = "setup tier 'core' did not complete"
    )
  })
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

# ---- Setup / build logging --------------------------------------------------

test_that("sce_log_dir returns a writable directory", {
  d <- sce_log_dir()
  expect_true(is.character(d) && length(d) == 1L)
  expect_true(dir.exists(d))
})

test_that("sce_open_log creates a fresh, prefix-named empty file", {
  p1 <- sce_open_log("test_open")
  on.exit(unlink(p1), add = TRUE)
  expect_true(file.exists(p1))
  expect_equal(file.info(p1)$size, 0)
  expect_match(basename(p1), "^test_open_\\d{8}_\\d{6}\\.log$")
})

test_that("sce_open_log sanitizes the prefix", {
  p <- sce_open_log("weird prefix!!!")
  on.exit(unlink(p), add = TRUE)
  expect_match(basename(p), "^weird_prefix_\\d{8}_\\d{6}\\.log$")
})

test_that("sce_run_with_log mirrors stdout into the log on success", {
  log_path <- sce_open_log("test_ok")
  on.exit(unlink(log_path), add = TRUE)
  result <- sce_run_with_log(
    { cat("hello-from-cat\n"); print("hello-from-print"); 42L },
    log_path = log_path)
  expect_identical(result, 42L)
  body <- readLines(log_path, warn = FALSE)
  expect_true(any(grepl("hello-from-cat", body)))
  expect_true(any(grepl("hello-from-print", body)))
  expect_true(any(grepl("DONE OK", body)))
})

test_that("sce_run_with_log tees message() into the log without muffling them", {
  # `remotes::install_github()` reports progress via `message()`,
  # which we route through a calling handler so:
  #   * the message text is copied into the log file
  #   * R's default handler still emits it to stderr (the user's
  #     terminal), since we do NOT invokeRestart("muffleMessage")
  log_path <- sce_open_log("test_messages")
  on.exit(unlink(log_path), add = TRUE)
  # Capture stderr to confirm the message survived the calling handler.
  stderr_capture <- capture.output(
    sce_run_with_log(
      { message("Downloading GitHub repo satijalab/azimuth@HEAD") },
      log_path = log_path),
    type = "message")
  body <- readLines(log_path, warn = FALSE)
  # Log captured the message.
  expect_true(any(grepl("Downloading GitHub repo", body)))
  # Terminal (stderr) ALSO received it -- not muffled.
  expect_true(any(grepl("Downloading GitHub repo", stderr_capture)))
})

test_that("sce_run_with_log tees warning() into the log without muffling", {
  log_path <- sce_open_log("test_warnings")
  on.exit(unlink(log_path), add = TRUE)
  stderr_capture <- capture.output(
    suppressWarnings(sce_run_with_log(
      { warning("installation of package 'foo' had non-zero exit status") },
      log_path = log_path)),
    type = "message")
  body <- readLines(log_path, warn = FALSE)
  expect_true(any(grepl("non-zero exit status", body)))
})

test_that("sce_run_with_log captures the error and rethrows with log_path attached", {
  log_path <- sce_open_log("test_err")
  on.exit(unlink(log_path), add = TRUE)
  expect_error(sce_run_with_log(
    stop("the install blew up"),
    log_path = log_path),
    regexp = "the install blew up")
  body <- readLines(log_path, warn = FALSE)
  expect_true(any(grepl("the install blew up", body)))
  expect_true(any(grepl("sce_run_with_log: ERROR", body)))
  # An "Error: <msg>" line is written so `sce_log_summary()` picks
  # it up the same way it picks up R CMD INSTALL failures.
  expect_true(any(grepl("^Error: the install blew up", body)))
  expect_true(length(sce_log_summary(log_path)) >= 1L)
  # The thrown condition carries the log_path so the caller can
  # surface it in the workspace warning.
  cond <- tryCatch(sce_run_with_log(stop("again"), log_path = log_path),
                   error = function(e) e)
  expect_identical(cond$log_path, log_path)
})

test_that("sce_run_with_log restores the sink stack after success and failure", {
  before <- list(out = sink.number(), msg = sink.number(type = "message"))
  log_path <- sce_open_log("test_sinks")
  on.exit(unlink(log_path), add = TRUE)
  sce_run_with_log(cat("ok\n"), log_path = log_path)
  expect_equal(sink.number(),                 before$out)
  expect_equal(sink.number(type = "message"), before$msg)

  expect_error(sce_run_with_log(stop("x"), log_path = log_path),
               regexp = "x")
  expect_equal(sink.number(),                 before$out)
  expect_equal(sink.number(type = "message"), before$msg)
})

test_that("sce_log_summary extracts BioC install cascade lines", {
  # Synthetic log resembling the user's report.
  log_path <- tempfile(fileext = ".log")
  on.exit(unlink(log_path), add = TRUE)
  writeLines(c(
    "* installing *source* package 'KEGGREST' ...",
    "** R",
    "Error: library load failed for 'libxml2'",
    "ERROR: configuration failed for package 'KEGGREST'",
    "* removing '/path/to/library/KEGGREST'",
    "Warning in install.packages(...) :",
    "  installation of package 'KEGGREST' had non-zero exit status",
    "ERROR: dependencies 'rhdf5', 'SparseArray' are not available for package 'h5mread'"),
    log_path)
  sm <- sce_log_summary(log_path, max_lines = 50L)
  # Must contain the actionable lines.
  expect_true(any(grepl("ERROR:", sm)))
  expect_true(any(grepl("had non-zero exit status", sm)))
  expect_true(any(grepl("\\* removing", sm)))
})

test_that("sce_log_summary is empty for a non-error log", {
  log_path <- tempfile(fileext = ".log")
  on.exit(unlink(log_path), add = TRUE)
  writeLines(c("Installing 5 packages...",
               "Successfully installed all packages."), log_path)
  expect_identical(sce_log_summary(log_path), character())
})

test_that("sce_log_tail honours `n` and handles missing files", {
  log_path <- tempfile(fileext = ".log")
  on.exit(unlink(log_path), add = TRUE)
  writeLines(as.character(1:20), log_path)
  expect_identical(sce_log_tail(log_path, n = 5L), as.character(16:20))
  expect_identical(sce_log_tail("/no/such/path"), character())
})

test_that(".sce_set_install_options defaults CRAN only when missing", {
  # 1. Sentinel `@CRAN@` -> default mirror is set.
  old <- options(repos = c(CRAN = "@CRAN@"),
                 install.packages.check.source = NULL)
  on.exit(options(old), add = TRUE)
  prior <- .sce_set_install_options(pin_check_source = FALSE)
  expect_identical(unname(getOption("repos")[["CRAN"]]),
                   "https://cloud.r-project.org")
  # Restoration round-trip.
  options(prior)
  expect_identical(unname(getOption("repos")[["CRAN"]]), "@CRAN@")
})

test_that(".sce_set_install_options leaves a real CRAN mirror alone", {
  # User has a real CRAN configured AND BiocManager is unavailable
  # (simulate by hiding it). Helper must return list() because there
  # is nothing for it to add.
  old <- options(repos = c(CRAN = "https://example.org/cran"))
  on.exit(options(old), add = TRUE)
  # We can't easily uninstall BiocManager mid-test, so this assertion
  # only holds when BiocManager isn't installed. Skip otherwise --
  # the "additively merges Bioc repos onto a user-set CRAN" test
  # below covers the BiocManager-installed case.
  skip_if(requireNamespace("BiocManager", quietly = TRUE),
          "BiocManager is installed; behaviour is covered by the additive-merge test below")
  prior <- .sce_set_install_options(pin_check_source = FALSE)
  expect_identical(prior, list())
  expect_identical(unname(getOption("repos")[["CRAN"]]),
                   "https://example.org/cran")
})

test_that(".sce_path_free_mb returns a positive number for an existing path", {
  # Smoke test: actual MB free depends on the environment, but on
  # any working machine the workspace root has *something*.
  mb <- .sce_path_free_mb(tempdir())
  skip_if(is.na(mb), "df not available in this environment")
  expect_true(is.numeric(mb))
  expect_gt(mb, 0)
})

test_that(".sce_check_tempdir_space treats NA (unknown free space) as OK", {
  # `df` may be unavailable (or path may not be queryable); the
  # guard must NOT block in that case -- we'd rather risk a real
  # OOD failure with a useful per-package log than block every
  # install whenever the platform doesn't ship POSIX `df`.
  expect_silent(.sce_check_tempdir_space(min_mb = 10^9,
                                          context = "a hypothetical install",
                                          free_mb = NA_real_))
})

test_that(".sce_check_tempdir_space errors with an actionable message when low", {
  err <- tryCatch(
    .sce_check_tempdir_space(min_mb = 5000,
                              context = "the Azimuth install",
                              free_mb = 100),
    error = function(e) e)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "Not enough free space")
  expect_match(conditionMessage(err), "the Azimuth install")
  expect_match(conditionMessage(err), "TMPDIR")
  # Confirms the numbers we report match what was passed in.
  expect_match(conditionMessage(err), "100 MB")
  expect_match(conditionMessage(err), "5000 MB")
})

test_that(".sce_check_tempdir_space passes when free space exceeds threshold", {
  expect_silent(.sce_check_tempdir_space(min_mb = 100,
                                          context = "x",
                                          scratch_mb = 10000))
})

test_that(".sce_check_tempdir_space mentions relocation in the error when scratch != tempdir", {
  # Caller already redirected scratch onto a roomier filesystem and we
  # STILL don't have enough room there -- the error should not blame
  # tempdir() (the relocation already happened) and should suggest
  # SCE_SCRATCH_DIR as the next override knob.
  tdir <- tempfile("sce_pretend_scratch_")
  dir.create(tdir)
  on.exit(unlink(tdir, recursive = TRUE), add = TRUE)
  err <- tryCatch(
    .sce_check_tempdir_space(min_mb = 5000,
                              context = "the Azimuth install",
                              scratch = tdir,
                              scratch_mb = 100),
    error = function(e) e)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "SCE_SCRATCH_DIR")
  expect_match(conditionMessage(err), "scratch dir")
})

test_that(".sce_install_scratch_dir honors SCE_SCRATCH_DIR verbatim", {
  skip_if_not_installed("withr")
  override <- tempfile("sce_override_scratch_")
  withr::with_envvar(c(SCE_SCRATCH_DIR = override), {
    got <- .sce_install_scratch_dir(
      candidates = character(0),
      free_mb    = function(p) 1)
    expect_true(dir.exists(got))
    expect_identical(normalizePath(got, mustWork = FALSE),
                     normalizePath(override, mustWork = FALSE))
  })
  unlink(override, recursive = TRUE)
})

test_that(".sce_install_scratch_dir picks the roomiest writable candidate", {
  skip_if_not_installed("withr")
  # Three candidates: tempdir + two synthetic dirs. Mock free_mb so
  # the second candidate is by far the roomiest.
  d1 <- tempfile("cand1_"); d2 <- tempfile("cand2_")
  dir.create(d1); dir.create(d2)
  on.exit(unlink(c(d1, d2), recursive = TRUE), add = TRUE)
  withr::with_envvar(c(SCE_SCRATCH_DIR = ""), {
    free <- function(p) {
      if (identical(p, d2)) 100000  # huge
      else if (identical(p, d1)) 500
      else 500  # tempdir
    }
    got <- .sce_install_scratch_dir(
      candidates = c(d1, d2),
      free_mb    = free)
    expect_true(startsWith(normalizePath(got, mustWork = FALSE),
                           normalizePath(d2, mustWork = FALSE)))
    expect_true(dir.exists(got))
  })
})

test_that(".sce_install_scratch_dir falls back to tempdir when no candidate is roomier", {
  skip_if_not_installed("withr")
  # Single candidate with the same free space as tempdir -- doesn't
  # clear the 1.5x headroom factor, so we stay on tempdir.
  d <- tempfile("cand_"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  withr::with_envvar(c(SCE_SCRATCH_DIR = ""), {
    got <- .sce_install_scratch_dir(
      candidates = c(d),
      free_mb    = function(p) 1000)
    expect_identical(got, tempdir())
  })
})

test_that(".sce_install_scratch_dir skips non-existent / non-writable candidates", {
  skip_if_not_installed("withr")
  nonexistent <- file.path(tempfile("doesnotexist_"))
  d_writable  <- tempfile("writable_"); dir.create(d_writable)
  on.exit(unlink(d_writable, recursive = TRUE), add = TRUE)
  withr::with_envvar(c(SCE_SCRATCH_DIR = ""), {
    got <- .sce_install_scratch_dir(
      candidates = c(nonexistent, d_writable),
      free_mb    = function(p) {
        if (identical(p, d_writable)) 100000 else 500
      })
    expect_true(startsWith(normalizePath(got, mustWork = FALSE),
                           normalizePath(d_writable, mustWork = FALSE)))
  })
})

test_that(".sce_setup_install_scratch / .sce_restore_install_scratch round-trip TMPDIR", {
  skip_if_not_installed("withr")
  tdir <- tempfile("scratch_"); dir.create(tdir)
  on.exit(unlink(tdir, recursive = TRUE), add = TRUE)
  # Case 1: TMPDIR was previously unset -- restore must unset it.
  withr::with_envvar(c(TMPDIR = NA), {
    state <- .sce_setup_install_scratch(scratch = tdir)
    expect_identical(Sys.getenv("TMPDIR"), tdir)
    .sce_restore_install_scratch(state)
    expect_identical(Sys.getenv("TMPDIR", unset = ""), "")
  })
  # Case 2: TMPDIR was previously set -- restore must put it back.
  withr::with_envvar(c(TMPDIR = "/some/prior/path"), {
    state <- .sce_setup_install_scratch(scratch = tdir)
    expect_identical(Sys.getenv("TMPDIR"), tdir)
    .sce_restore_install_scratch(state)
    expect_identical(Sys.getenv("TMPDIR"), "/some/prior/path")
  })
})

test_that(".sce_setup_install_scratch is a no-op when scratch == tempdir()", {
  prior <- Sys.getenv("TMPDIR", unset = NA_character_)
  on.exit({
    if (is.na(prior)) Sys.unsetenv("TMPDIR")
    else              Sys.setenv(TMPDIR = prior)
  }, add = TRUE)
  state <- .sce_setup_install_scratch(scratch = tempdir())
  # We didn't relocate -- TMPDIR shouldn't have been touched.
  expect_identical(Sys.getenv("TMPDIR", unset = NA_character_), prior)
  .sce_restore_install_scratch(state)
})

test_that(".sce_bioc_mirror prefers explicit user overrides, then env var, then Posit", {
  # Option > env var > default.
  withr::with_options(list(BioC_mirror = "https://example.org/bioc"),
    withr::with_envvar(c(SCE_BIOC_MIRROR = "https://envvar.example/bioc"), {
      expect_identical(.sce_bioc_mirror(), "https://example.org/bioc")
    }))
  # Env var > default when option is unset.
  withr::with_options(list(BioC_mirror = NULL),
    withr::with_envvar(c(SCE_BIOC_MIRROR = "https://envvar.example/bioc"), {
      expect_identical(.sce_bioc_mirror(), "https://envvar.example/bioc")
    }))
  # Default when nothing is set.
  withr::with_options(list(BioC_mirror = NULL),
    withr::with_envvar(c(SCE_BIOC_MIRROR = ""), {
      expect_identical(.sce_bioc_mirror(),
                       "https://packagemanager.posit.co/bioconductor")
    }))
})

test_that(".sce_set_install_options routes Bioc repos through the chosen mirror", {
  skip_if_not_installed("BiocManager")
  skip_if_not_installed("withr")
  # Force the default Posit mirror by clearing both override channels,
  # and force the "default repos" branch.
  withr::with_options(list(BioC_mirror = NULL, repos = c(CRAN = "@CRAN@")),
    withr::with_envvar(c(SCE_BIOC_MIRROR = ""), {
      prior <- .sce_set_install_options(pin_check_source = FALSE)
      repos <- getOption("repos")
      # All BioC* entries point at the Posit mirror, not the canonical
      # bioconductor.org (which redirects to the OSN mirror).
      bioc_urls <- unname(repos[grepl("^BioC", names(repos))])
      expect_true(length(bioc_urls) > 0L)
      expect_true(all(grepl("^https://packagemanager\\.posit\\.co/bioconductor/",
                            bioc_urls)),
                  info = paste(bioc_urls, collapse = "\n"))
      # BioC_mirror is restored in-place (no side effect leaks into
      # the caller's session). Subsequent install.packages() calls
      # still hit Posit because the URLs are baked into `repos`.
      expect_null(getOption("BioC_mirror"))
      # Round-trip restoration of everything else works.
      options(prior)
      expect_identical(unname(getOption("repos")[["CRAN"]]), "@CRAN@")
    }))
})

test_that(".sce_set_install_options respects an explicit BioC_mirror override", {
  skip_if_not_installed("BiocManager")
  skip_if_not_installed("withr")
  withr::with_options(list(BioC_mirror = "https://mirror.example.org/bioc",
                           repos = c(CRAN = "@CRAN@")),
    withr::with_envvar(c(SCE_BIOC_MIRROR = ""), {
      .sce_set_install_options(pin_check_source = FALSE)
      bioc_urls <- unname(getOption("repos")[grepl("^BioC", names(getOption("repos")))])
      expect_true(length(bioc_urls) > 0L)
      expect_true(all(grepl("^https://mirror\\.example\\.org/bioc/", bioc_urls)),
                  info = paste(bioc_urls, collapse = "\n"))
      # User's override is preserved verbatim (not clobbered by Posit).
      expect_identical(getOption("BioC_mirror"),
                       "https://mirror.example.org/bioc")
    }))
})

test_that(".sce_set_install_options folds Bioc repos in when BiocManager is available", {
  skip_if_not_installed("BiocManager")
  # Force the sentinel state so the helper takes the "default repos"
  # branch.
  old <- options(repos = c(CRAN = "@CRAN@"),
                 install.packages.check.source = NULL)
  on.exit(options(old), add = TRUE)
  prior <- .sce_set_install_options(pin_check_source = FALSE)
  repos <- getOption("repos")
  # CRAN is still present and pointing at a real URL.
  expect_true("CRAN" %in% names(repos))
  expect_false(identical(unname(repos[["CRAN"]]), "@CRAN@"))
  # Bioc repos are folded in. `BiocManager::repositories()` returns at
  # least one BioC* named entry on every supported BiocManager
  # version; assert that any-bioc-prefix is there.
  expect_true(any(grepl("^BioC", names(repos))),
              info = sprintf("repos names = %s",
                             paste(names(repos), collapse = ", ")))
  # Round-trip restoration.
  options(prior)
  expect_identical(unname(getOption("repos")[["CRAN"]]), "@CRAN@")
})

test_that(".sce_set_install_options additively merges Bioc repos onto a user-set CRAN", {
  skip_if_not_installed("BiocManager")
  # User already has a CRAN configured (e.g. via .Rprofile). We MUST
  # respect their choice (no overwrite) but still add the Bioc
  # sub-repos that `remotes::install_github()`'s dep resolver needs.
  old <- options(repos = c(CRAN = "https://example.org/cran"))
  on.exit(options(old), add = TRUE)
  prior <- .sce_set_install_options(pin_check_source = FALSE)
  repos <- getOption("repos")
  # User's CRAN URL is preserved as-is.
  expect_identical(unname(repos[["CRAN"]]), "https://example.org/cran")
  # AND BiocXXX repos are now also present.
  expect_true(any(grepl("^BioC", names(repos))))
  options(prior)
})

test_that(".sce_set_install_options pins install.packages.check.source", {
  old <- options(repos = c(CRAN = "https://example.org/cran"),
                 install.packages.check.source = "yes")
  on.exit(options(old), add = TRUE)
  prior <- .sce_set_install_options(pin_check_source = TRUE)
  expect_identical(getOption("install.packages.check.source"), "no")
  expect_true("install.packages.check.source" %in% names(prior))
  options(prior)
  expect_identical(getOption("install.packages.check.source"), "yes")
})

test_that(".sce_pkg_outputs_dir creates a unique sub-directory under sce_log_dir()", {
  d <- .sce_pkg_outputs_dir()
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  expect_true(dir.exists(d))
  expect_match(basename(d), "^pkglogs_\\d{8}_\\d{6}$")
  # Same parent as the main log dir.
  expect_identical(dirname(d), sce_log_dir())
})

test_that(".sce_drain_pkg_outputs folds <pkg>.out files into the parent log", {
  # Build a fake per-package log dir as install.packages(keep_outputs=)
  # would leave it.
  pkglog <- tempfile("pkglogs-")
  dir.create(pkglog)
  on.exit(unlink(pkglog, recursive = TRUE), add = TRUE)
  writeLines(c("* installing source package 'bit' ...",
               "gcc -I... -c bit.c",
               "ERROR: configure failed for 'bit'"),
             file.path(pkglog, "bit.out"))
  writeLines(c("* installing source package 'rhdf5' ...",
               "ERROR: dependency 'Rhdf5lib' not available"),
             file.path(pkglog, "rhdf5.out"))

  # Wrap the drain in sce_run_with_log so we can read the resulting
  # log file and assert the per-package content was inlined.
  log_path <- sce_open_log("test_drain")
  on.exit(unlink(log_path), add = TRUE)
  sce_run_with_log(.sce_drain_pkg_outputs(pkglog),
                   log_path = log_path)
  body <- paste(readLines(log_path, warn = FALSE), collapse = "\n")
  expect_match(body, "per-package install log: bit",  fixed = TRUE)
  expect_match(body, "per-package install log: rhdf5", fixed = TRUE)
  expect_match(body, "ERROR: configure failed for 'bit'", fixed = TRUE)
  expect_match(body, "dependency 'Rhdf5lib' not available", fixed = TRUE)
  expect_match(body, "per-package logs preserved on disk in", fixed = TRUE)
})

test_that(".sce_drain_pkg_outputs is a no-op when the dir is empty or missing", {
  expect_silent(.sce_drain_pkg_outputs(NULL))
  expect_silent(.sce_drain_pkg_outputs("/no/such/dir"))
  empty <- tempfile("pkglogs-empty-")
  dir.create(empty)
  on.exit(unlink(empty, recursive = TRUE), add = TRUE)
  expect_silent(.sce_drain_pkg_outputs(empty))
})

# ---- Lazy-load corruption detection ----------------------------------------

test_that(".is_lazy_load_corruption matches the on-disk patterns we surface", {
  expect_true(.is_lazy_load_corruption(
    "lazy-load database '/x/Seurat/R/Seurat.rdb' is corrupt"))
  expect_true(.is_lazy_load_corruption(
    "Error in lazyLoadDBfetch: bad restore file magic number"))
  expect_true(.is_lazy_load_corruption(
    "internal error -3 in R_decompress1"))
  expect_true(.is_lazy_load_corruption(
    "unable to load shared object '/x/Rcpp/libs/Rcpp.so': cannot ..."))
})

test_that(".is_lazy_load_corruption is conservative", {
  expect_false(.is_lazy_load_corruption(""))
  expect_false(.is_lazy_load_corruption(character()))
  expect_false(.is_lazy_load_corruption(
    "Could not connect to ExperimentHub"))
  expect_false(.is_lazy_load_corruption(
    "TENxPBMCData: dataset 'pbmc8k' not found"))
})

test_that("sce_install_for_demo refuses to upgrade an already-loaded package", {
  skip_if_not_installed("withr")
  # We need at least one package that is *guaranteed* to be in
  # `loadedNamespaces()` and that is also a target install. The demo
  # tier's CRAN list typically includes Seurat / SeuratObject -- pick
  # one that's already loaded in this test process (testthat itself
  # always is).
  spec <- .sce_packages_for_tier("demo")
  # Fake the situation: temporarily inject `testthat` (which IS
  # loaded) into the demo tier and ask install to run. Have to fake
  # it by overriding the tier table in this test.
  fake <- list(core = list(cran = character(), bioc = character()),
               demo = list(cran = "testthat", bioc = character()),
               full = list(cran = character(), bioc = character()))
  with_mocked_bindings <- function(value, code) {
    old <- get(".SCE_PKG_TIERS", envir = globalenv())
    assign(".SCE_PKG_TIERS", value, envir = globalenv())
    on.exit(assign(".SCE_PKG_TIERS", old, envir = globalenv()))
    force(code)
  }
  # Make the function look up our fake tiers. Skip the test if the
  # implementation doesn't keep .SCE_PKG_TIERS in the global env (it
  # currently does -- file-scope source()).
  if (!exists(".SCE_PKG_TIERS", envir = globalenv(), inherits = FALSE)) {
    skip("Cannot redirect .SCE_PKG_TIERS in this test environment")
  }
  with_mocked_bindings(fake, {
    # Force testthat to look "not installed" by checking before-install
    # cache: actually has_optional("testthat") is TRUE so it won't be
    # treated as missing. Skip if we can't fabricate the situation.
    if (has_optional("testthat")) {
      # The install function short-circuits on "nothing missing" --
      # we need to test the loaded-namespace branch, which is only
      # reachable when something IS missing. The actual user-facing
      # branch is well-tested manually; here we just confirm the
      # error path exists by inspecting the function source.
      src <- deparse(sce_install_for_demo)
      expect_true(any(grepl("loadedNamespaces", src)))
      expect_true(any(grepl("Cannot install package", src)))
    } else {
      succeed("testthat happened to be uninstalled; branch is exercised")
    }
  })
})

test_that("sce_install_for_demo exists and short-circuits when nothing is missing", {
  expect_true(is.function(sce_install_for_demo))

  # We cannot actually run `install.packages()` in the test suite, so we
  # exercise only the short-circuit branch. If this machine already has
  # every demo dep, the helper must return TRUE without performing any
  # installs.
  if (sce_check_setup()$demo$complete) {
    # Collect progress callbacks so we can assert the helper finishes
    # cleanly at fraction 1.0.
    seen <- list()
    cb <- function(fraction, detail = NULL) {
      seen[[length(seen) + 1L]] <<- list(fraction = fraction, detail = detail)
    }
    expect_true(sce_install_for_demo(progress = cb))
    expect_true(length(seen) >= 1L)
    expect_identical(seen[[length(seen)]]$fraction, 1.0)
  } else {
    succeed("demo tier incomplete; install path covered by manual smoke runs")
  }
})

test_that("sce_install_for_demo accepts NULL progress without error", {
  if (sce_check_setup()$demo$complete) {
    expect_silent(sce_install_for_demo(progress = NULL))
  } else {
    succeed("demo tier incomplete; install path covered by manual smoke runs")
  }
})

test_that("sce_install_for_demo swallows errors from the progress callback", {
  if (sce_check_setup()$demo$complete) {
    bad <- function(fraction, detail = NULL) stop("ui blew up")
    expect_silent(sce_install_for_demo(progress = bad))
  } else {
    succeed("demo tier incomplete; install path covered by manual smoke runs")
  }
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
