# ============================================================================
# Shared UI primitive helpers
# ----------------------------------------------------------------------------
# Focused on the server-side gene picker because it has both UI and
# server obligations (the UI half must render with `choices = NULL`,
# and the server half must request `server = TRUE`). Visual primitives
# like `helper_text()`, `microcaption()`, etc. are exercised
# indirectly by module snapshot tests.
# ============================================================================

test_that("gene_picker_input renders selectize with no client-side choices", {
  out <- gene_picker_input("g", "My label", placeholder = "Pick one")
  html <- as.character(out)

  # Must be a selectize-class input (not a vanilla selectInput) so the
  # corresponding `update_gene_picker(server = TRUE)` call wires up
  # incremental delivery.
  expect_match(html, "selectize",        fixed = TRUE)
  expect_match(html, "id=\"g\"",         fixed = TRUE)
  expect_match(html, "My label",         fixed = TRUE)
  expect_match(html, "placeholder",      fixed = TRUE)
  expect_match(html, "Pick one",         fixed = TRUE)

  # The HTML <select> must contain no <option> children: that's the
  # entire point -- the gene list is delivered server-side later.
  # `<select ...></select>` with no inner <option> tags satisfies this.
  expect_false(grepl("<option[^>]*>[^<]+</option>", html))
})

test_that("gene_picker_input supports multi-select", {
  out <- gene_picker_input("g", "lbl", multiple = TRUE)
  html <- as.character(out)
  expect_match(html, "multiple", fixed = TRUE)
})

test_that("update_gene_picker drives updateSelectizeInput(server = TRUE)", {
  # Use shiny's official `MockShinySession` so the strict
  # `validate_session_object()` check inside `updateSelectizeInput`
  # passes. We intercept `sendInputMessage` to observe what the
  # server-side update actually pushes to the (mock) client.
  session <- shiny::MockShinySession$new()
  sent <- list()
  session$sendInputMessage <- function(inputId, message) {
    sent[[length(sent) + 1L]] <<- list(inputId = inputId, message = message)
  }

  update_gene_picker(session, "g",
                     choices  = c("CD3D", "CD8A", "MS4A1"),
                     selected = "CD3D")

  expect_length(sent, 1L)
  msg <- sent[[1L]]
  expect_identical(msg$inputId, "g")
  # `server = TRUE` selectize updates carry the option list inside
  # `message$options` (a JSON-formatted data.frame). The exact shape
  # isn't part of our contract -- just that the message was emitted
  # and at least one of the expected keys is present.
  expect_true(any(c("options", "value", "url") %in% names(msg$message)))
})

test_that("update_gene_picker tolerates empty / NULL choices", {
  session <- shiny::MockShinySession$new()
  sent <- list()
  session$sendInputMessage <- function(inputId, message) {
    sent[[length(sent) + 1L]] <<- list(inputId = inputId, message = message)
  }

  expect_silent(update_gene_picker(session, "g", choices = NULL))
  expect_silent(update_gene_picker(session, "g", choices = character()))
  expect_length(sent, 2L)
})
