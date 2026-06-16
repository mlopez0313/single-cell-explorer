test_that(".console_tail_file returns character(0) for missing / empty path", {
  expect_identical(.console_tail_file(NULL),                  character())
  expect_identical(.console_tail_file(""),                    character())
  expect_identical(.console_tail_file(tempfile("nope_")),     character())
})

test_that(".console_tail_file returns the last N lines of a file", {
  p <- tempfile("console_tail_")
  writeLines(sprintf("line %02d", 1:50), p)
  on.exit(unlink(p), add = TRUE)

  last_5 <- .console_tail_file(p, n = 5L)
  expect_length(last_5, 5L)
  expect_identical(last_5,
                   c("line 46", "line 47", "line 48", "line 49", "line 50"))

  # n larger than file -> full file.
  all <- .console_tail_file(p, n = 500L)
  expect_length(all, 50L)
})

test_that(".console_list_log_files returns newest-first when log dir exists", {
  d <- tempfile("sce_logs_"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  # Three log files with distinct mtimes (oldest -> newest).
  f1 <- file.path(d, "a.log"); writeLines("old",    f1); Sys.setFileTime(f1, Sys.time() - 60)
  f2 <- file.path(d, "b.log"); writeLines("middle", f2); Sys.setFileTime(f2, Sys.time() - 30)
  f3 <- file.path(d, "c.log"); writeLines("newest", f3); Sys.setFileTime(f3, Sys.time())

  got <- .console_list_log_files(dir = d)
  expect_length(got, 3L)
  expect_identical(basename(got), c("c.log", "b.log", "a.log"))
})

test_that(".console_list_log_files filters non-.log files", {
  d <- tempfile("sce_logs_"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  writeLines("x", file.path(d, "ok.log"))
  writeLines("y", file.path(d, "ignore.txt"))   # not .log
  writeLines("z", file.path(d, "README"))       # no extension
  got <- .console_list_log_files(dir = d)
  expect_length(got, 1L)
  expect_identical(basename(got), "ok.log")
})

test_that(".console_list_log_files returns empty when log dir does not exist", {
  nope <- tempfile("sce_logs_nope_")  # don't create
  expect_identical(.console_list_log_files(dir = nope), character())
})

test_that("isolate_or_default returns current when valid, default otherwise", {
  valid <- c("a.log", "b.log", "c.log")
  expect_identical(isolate_or_default("b.log", "a.log", valid), "b.log")
  expect_identical(isolate_or_default(NULL,    "a.log", valid), "a.log")
  expect_identical(isolate_or_default("",      "a.log", valid), "a.log")
  expect_identical(isolate_or_default("stale.log", "a.log", valid), "a.log")
})
