# ============================================================================
# Shared UI primitives
# ----------------------------------------------------------------------------
# Reusable presentation helpers used by modules to keep visual styling
# consistent across the app. All visual styling lives in `www/styles.css`
# (see `.app-card`, `.metric-card`, `.info-banner`, etc.).
#
# Conventions:
#   - Every helper returns an `htmltools` tag (compatible with `shiny::tagList`).
#   - Helpers take user content via `...` so callers can compose freely.
#   - Tone arguments take a small fixed vocabulary; unknown tones fall back to
#     "neutral" silently.
#
# If you need a brand-new primitive, add it here rather than inlining a styled
# `div(...)` in a module.
# ============================================================================

# ---- internal: present/absent test ---------------------------------------
# TRUE when `x` is "renderable content": a non-NULL tag, or a non-empty
# character vector. Avoids the `nzchar()` trap when callers pass a shiny tag
# (which is a list and can't be coerced to logical of length 1).
.has_content <- function(x) {
  if (is.null(x)) return(FALSE)
  if (is.character(x)) return(any(nzchar(x)))
  TRUE
}

# ---- icon helper ----------------------------------------------------------
# bsicons is an optional dependency. Use `app_icon("info-circle")` and it will
# render a Bootstrap icon when bsicons is installed, or an empty span (with
# the same width) otherwise. Modules can rely on the API regardless.
app_icon <- function(name = NULL, class = "sidebar-nav__icon") {
  if (is.null(name) || !nzchar(name)) {
    return(shiny::span(class = class, `aria-hidden` = "true"))
  }
  if (requireNamespace("bsicons", quietly = TRUE)) {
    return(shiny::span(class = class,
                       bsicons::bs_icon(name, a11y = "deco")))
  }
  # Fallback: empty span keeps spacing predictable.
  shiny::span(class = class, `aria-hidden` = "true",
              `data-icon` = name)
}

# ---- page header ----------------------------------------------------------
#' Standard module page framing.
#'
#' @param title    character(1) main heading
#' @param eyebrow  character(1) small uppercase label rendered above the title
#' @param lede     character(1) one-line descriptive subtitle
#' @param meta     optional tagList of `page_meta_item(...)` entries
page_header <- function(title, eyebrow = NULL, lede = NULL, meta = NULL) {
  shiny::tags$header(
    class = "page-header",
    if (.has_content(eyebrow))
      shiny::div(class = "page-header__eyebrow", eyebrow),
    shiny::h1(class = "page-header__title", title),
    if (.has_content(lede))
      shiny::p(class = "page-header__lede", lede),
    if (!is.null(meta))
      shiny::div(class = "page-header__meta", meta)
  )
}

#' Contextual metadata strip item for `page_header(meta = ...)`.
#'
#' Example: `page_meta_item("Cells", "12,345")` renders `Cells 12,345`.
page_meta_item <- function(label, value) {
  shiny::span(class = "page-header__meta-item",
              shiny::tags$span(label, ":"),
              shiny::tags$strong(value))
}

# ---- app card -------------------------------------------------------------
#' Bordered surface with optional title / caption / footer.
#'
#' @param ...      child tags rendered in the card body
#' @param title    optional heading
#' @param caption  optional sub-heading (right-aligned in the header)
#' @param footer   optional tagList rendered in the footer strip
#' @param class    extra CSS classes appended to `.app-card`
#' @param body_class extra CSS classes appended to `.app-card__body`
app_card <- function(..., title = NULL, caption = NULL,
                     footer = NULL, class = NULL, body_class = NULL) {
  header <- NULL
  if (!is.null(title) || !is.null(caption)) {
    header <- shiny::div(
      class = "app-card__header",
      if (!is.null(title))
        shiny::h3(class = "app-card__title", title),
      if (!is.null(caption))
        shiny::p(class = "app-card__caption", caption)
    )
  }
  shiny::div(
    class = paste(c("app-card", class), collapse = " "),
    header,
    shiny::div(class = paste(c("app-card__body", body_class), collapse = " "),
               ...),
    if (!is.null(footer))
      shiny::div(class = "app-card__footer", footer)
  )
}

#' Wrapper around `app_card()` tuned for plot outputs.
#'
#' @param caption  optional sub-heading shown in the card header
#' @param footnote optional small footnote rendered under the plot
plot_card <- function(..., title = NULL, caption = NULL, footnote = NULL,
                      class = NULL) {
  app_card(
    ...,
    title    = title,
    caption  = caption,
    class    = paste(c("plot-card", class), collapse = " "),
    footer   = if (!is.null(footnote))
                 shiny::div(class = "plot-caption", footnote)
  )
}

# ---- metric card ----------------------------------------------------------
#' Compact KPI tile with a small uppercase label and a tabular numeric value.
metric_card <- function(label, value, hint = NULL) {
  shiny::div(
    class = "metric-card",
    shiny::p(class = "metric-card__label", label),
    shiny::p(class = "metric-card__value", value),
    if (.has_content(hint))
      shiny::p(class = "metric-card__hint", hint)
  )
}

#' Layout helper: a responsive grid of `metric_card()` tiles.
metric_grid <- function(...) {
  shiny::div(class = "metric-grid", ...)
}

# ---- status pill ----------------------------------------------------------
#' Small uppercase pill for inline status labels.
#'
#' @param tone one of "neutral","info","success","warning","danger","accent"
status_pill <- function(label, tone = c("neutral", "info", "success",
                                        "warning", "danger", "accent")) {
  tone <- match.arg(tone)
  shiny::span(class = paste0("status-pill status-pill--", tone), label)
}

# ---- empty / loading state ------------------------------------------------
#' Friendly empty-state surface with title, lede, and optional actions.
empty_state <- function(title, lede = NULL,
                        icon = "circle",
                        actions = NULL,
                        hint = NULL) {
  shiny::div(
    class = "empty-state",
    if (!is.null(icon))
      shiny::div(class = "empty-state__icon",
                 app_icon(icon, class = "empty-state__icon-glyph")),
    shiny::h2(class = "empty-state__title", title),
    if (.has_content(lede))
      shiny::p(class = "empty-state__lede", lede),
    if (!is.null(actions))
      shiny::div(class = "empty-state__actions", actions),
    if (.has_content(hint))
      shiny::p(class = "empty-state__hint", hint)
  )
}

# ---- inline banner --------------------------------------------------------
#' Inline informational / warning / success / danger banner.
#'
#' Use this for in-page callouts (missing inputs, friendly warnings, gentle
#' tips). For dataset-level transient messages use the workspace message
#' stream which is rendered separately by the workspace server.
#'
#' @param tone   one of "info","success","warning","danger","neutral"
#' @param title  optional bold first line
#' @param ...    body content
info_banner <- function(..., title = NULL,
                        tone = c("info", "success", "warning",
                                 "danger", "neutral"),
                        icon = NULL) {
  tone <- match.arg(tone)
  default_icon <- switch(tone,
                         info    = "info-circle",
                         success = "check-circle",
                         warning = "exclamation-triangle",
                         danger  = "x-circle",
                         neutral = NULL)
  glyph <- icon %||% default_icon
  shiny::div(
    class = paste0("info-banner info-banner--", tone),
    role  = if (tone %in% c("warning", "danger")) "alert" else "status",
    if (!is.null(glyph))
      shiny::span(class = "info-banner__icon",
                  app_icon(glyph, class = "info-banner__icon-glyph")),
    shiny::div(
      class = "info-banner__body",
      if (.has_content(title))
        shiny::span(class = "info-banner__title", title),
      ...
    )
  )
}

# ---- section title --------------------------------------------------------
#' Subtitle within a card body (small caps + thin rule).
section_title <- function(text, subtitle = NULL) {
  shiny::div(
    class = "section-title",
    shiny::span(text),
    if (.has_content(subtitle))
      shiny::span(class = "section-title__subtitle", subtitle),
    shiny::span(class = "section-title__rule")
  )
}

# ---- control panel --------------------------------------------------------
#' Visually distinct container for filter / selector inputs.
#'
#' @param ...      input tags (or `shiny::fluidRow(column(...))` blocks)
#' @param title    optional small uppercase header
#' @param actions  optional tagList of action buttons rendered below the
#'                 controls (e.g. a "Run" button)
control_panel <- function(..., title = NULL, actions = NULL) {
  with_actions_class <- if (!is.null(actions)) "control-panel--with-actions" else NULL
  shiny::div(
    class = paste(c("control-panel", with_actions_class), collapse = " "),
    if (.has_content(title))
      shiny::p(class = "control-panel__title", title),
    shiny::div(...),
    if (!is.null(actions))
      shiny::div(class = "control-panel__actions", actions)
  )
}

# ---- definition list ------------------------------------------------------
#' A compact term/value list (used by Dataset Overview etc.).
#'
#' @param items named list where names are the labels and values are the
#'              rendered content (any tag or character).
deflist <- function(items) {
  if (!length(items)) return(NULL)
  rows <- lapply(seq_along(items), function(i) {
    shiny::tagList(
      shiny::div(class = "app-deflist__term",  names(items)[i]),
      shiny::div(class = "app-deflist__value", items[[i]])
    )
  })
  shiny::div(class = "app-deflist", rows)
}

# ---- context strip (dataset / annotation banner) --------------------------
#' Persistent context strip rendered at the top of the workspace.
#'
#' @param label   primary text (e.g. annotation set label)
#' @param meta    optional right-aligned text (e.g. dataset name)
#' @param tone    one of "neutral","active","demo"
context_strip <- function(label, meta = NULL,
                          tone = c("neutral", "active", "demo")) {
  tone <- match.arg(tone)
  cls  <- switch(tone,
                 neutral = "context-strip",
                 active  = "context-strip context-strip--active",
                 demo    = "context-strip context-strip--demo")
  shiny::div(
    class = cls,
    shiny::span(shiny::tags$strong("Active annotation: "), label),
    if (.has_content(meta))
      shiny::span(class = "context-strip__meta", meta)
  )
}

# ---- action row -----------------------------------------------------------
#' Standard action row used at the bottom of `control_panel(...)` etc.
#'
#' Renders a flex row with one or more buttons and an optional muted helper
#' text on the right. Use this anywhere you previously hand-rolled
#' `div(style="margin: 8px 0 16px 0;", actionButton(...), span(style="margin-left:12px; color:#888;", ...))`.
#'
#' @param ...        button tags (e.g. `shiny::actionButton(...)`).
#' @param helper     optional helper text rendered to the right of buttons.
#' @param align      "left" (default) or "right".
action_row <- function(..., helper = NULL, align = c("left", "right")) {
  align <- match.arg(align)
  shiny::div(
    class = paste(c("action-row",
                    if (align == "right") "action-row--right"), collapse = " "),
    shiny::div(class = "action-row__buttons", ...),
    if (.has_content(helper))
      shiny::span(class = "action-row__helper", helper)
  )
}

# ---- summary / meta bar ---------------------------------------------------
#' Content-side summary strip (replaces ad hoc reuse of `.page-header__meta`).
#'
#' @param items  tagList of `summary_item(label, value)` entries, or any tag.
#' @param class  optional extra CSS classes
summary_bar <- function(items, class = NULL) {
  shiny::div(class = paste(c("summary-bar", class), collapse = " "), items)
}

#' One entry in a `summary_bar(...)`.
summary_item <- function(label, value) {
  shiny::span(
    class = "summary-bar__item",
    shiny::span(class = "summary-bar__label", label),
    shiny::span(class = "summary-bar__value", value)
  )
}

# ---- callout legend -------------------------------------------------------
#' Tiny legend strip ("yellow = X, green = Y, ...") used to explain status
#' colors inside tables / panels. Replaces inline `<span style="background:...">`
#' tag soup.
#'
#' @param items  named list. Names are the legend labels; values are
#'   "swatch tones" picked from `c("neutral","info","success","warning",
#'   "danger","accent")`. Example: `list("engine-suggested"="warning",
#'   "user-confirmed"="success")`.
#' @param note   optional trailing prose appended after the swatches.
callout_legend <- function(items, note = NULL) {
  if (!length(items)) return(NULL)
  pieces <- lapply(seq_along(items), function(i) {
    label <- names(items)[i]
    tone  <- items[[i]]
    shiny::span(
      class = "callout-legend__item",
      shiny::span(class = paste0("callout-legend__swatch ",
                                 "callout-legend__swatch--", tone)),
      shiny::span(label)
    )
  })
  shiny::div(
    class = "callout-legend",
    pieces,
    if (.has_content(note))
      shiny::span(class = "callout-legend__note", note)
  )
}

# ---- coming soon ----------------------------------------------------------
#' Tasteful "this feature is on the roadmap" card. Replaces the ad hoc dashed
#' panel in `mod_placeholders.R`.
coming_soon_card <- function(title, description,
                             docs_path = "docs/ADDING_MODULES.md") {
  shiny::div(
    class = "coming-soon",
    shiny::div(
      class = "coming-soon__header",
      shiny::h2(class = "coming-soon__title", title),
      status_pill("Coming soon", tone = "warning")
    ),
    shiny::p(class = "coming-soon__description", description),
    shiny::p(class = "coming-soon__hint",
             "This module is registered but not yet implemented. See ",
             shiny::tags$code(docs_path), " to contribute.")
  )
}

# ---- compute-status banner ------------------------------------------------
#' Unified status banner for compute-oriented modules (DE, pathway,
#' imputation, regulons, trajectory, ...).
#'
#' The visual surface is shared; the body text is provided by the caller so
#' modules can format their own params/errors. Replaces every hand-rolled
#' `div(style = sprintf("padding:8px 12px; background:%s; ...", bg, fg), ...)`.
#'
#' @param text  body content (string or shiny tag)
#' @param tone  one of "idle","running","success","danger","warning","info"
#' @param label small uppercase prefix (default "Status")
status_banner <- function(text,
                          tone  = c("idle", "running", "success",
                                    "danger", "warning", "info"),
                          label = "Status") {
  tone <- match.arg(tone)
  shiny::div(
    class = paste0("status-banner status-banner--", tone),
    role  = if (tone %in% c("danger", "warning")) "alert" else "status",
    if (.has_content(label))
      shiny::span(class = "status-banner__label", label),
    shiny::span(class = "status-banner__text", text)
  )
}

# ---- card toolbar ---------------------------------------------------------
#' Sub-panel of controls rendered inside an `app_card()` / `table_card()` body.
#'
#' Replaces the older pattern of nesting a `.control-panel` inside a card and
#' overriding its borders with inline CSS. The toolbar inherits the card's
#' background and is separated from the body below it by a thin rule.
#'
#' @param ...    child tags (typically `shiny::fluidRow(shiny::column(...))`).
#' @param title  optional small uppercase label.
card_toolbar <- function(..., title = NULL) {
  shiny::div(
    class = "app-card__toolbar",
    if (.has_content(title))
      shiny::p(class = "control-panel__title", title),
    ...
  )
}

# ---- table card -----------------------------------------------------------
#' Output-card wrapper for tables (or any scrollable wide output).
#'
#' @param ...        body tags (typically `shiny::tableOutput(...)` etc.).
#'                   Rendered inside a `.table-scroll` region.
#' @param title      card title.
#' @param caption    optional sub-caption.
#' @param footnote   optional small footnote.
#' @param toolbar    optional pre-scroll toolbar (typically a row of filter
#'                   inputs). Use this instead of hand-rolling a nested
#'                   `.control-panel` inside the card body.
#' @param max_height optional max height (e.g. `"420px"`) for the scroll
#'                   region. When `NULL`, no scroll bound is applied. Pass
#'                   `"bounded"` to use the standard ~420px scroll height.
table_card <- function(..., title = NULL, caption = NULL,
                       footnote = NULL, toolbar = NULL, max_height = NULL) {
  scroll_class <- "table-scroll"
  scroll_style <- NULL
  if (!is.null(max_height)) {
    if (identical(max_height, "bounded")) {
      scroll_class <- "table-scroll table-scroll--bounded"
    } else {
      scroll_style <- sprintf("max-height: %s;", max_height)
    }
  }
  body <- shiny::div(class = scroll_class, style = scroll_style, ...)
  app_card(
    if (!is.null(toolbar)) card_toolbar(toolbar) else NULL,
    body,
    title      = title,
    caption    = caption,
    class      = "table-card",
    body_class = "app-card__body--flush",
    footer     = if (!is.null(footnote))
                   shiny::div(class = "plot-caption", footnote))
}

# ---- beside-input alignment ----------------------------------------------
#' Render content vertically aligned with the input baseline of a labeled
#' `shiny::selectInput(...)` (etc.) in an adjacent `shiny::column(...)`.
#'
#' Replaces the repeated `div(style = "padding-top: 22px;", ...)` pattern that
#' appeared next to selectInputs/checkboxes in the annotation, DE, and
#' pathway modules.
beside_input <- function(...) {
  shiny::div(class = "beside-input", ...)
}

# ---- required-inputs list ------------------------------------------------
#' Bulleted list used by the workspace's "needs inputs" state.
#'
#' @param items character vector of input names.
req_list <- function(items) {
  if (!length(items)) return(NULL)
  shiny::tags$ul(class = "req-list",
                 lapply(items, shiny::tags$li))
}

# ---- small helper / caption text -----------------------------------------
#' Small muted caption span. Use anywhere you previously wrote
#' `div(style = "font-size:12px; color:#888;", ...)`.
helper_text <- function(...) {
  shiny::span(class = "helper-text", ...)
}

#' Even smaller muted line (e.g. "Selected: GENE-X") used under buttons.
microcaption <- function(...) {
  shiny::div(class = "microcaption", ...)
}

#' Content-side uppercase label suitable inside cards (not the sidebar).
#' Use instead of `sidebar-section__label` when the label lives in workspace.
content_label <- function(text) {
  shiny::span(class = "content-label", text)
}

# ---- server-side gene picker ---------------------------------------------
# Real scRNA-seq datasets routinely carry 20k-40k genes. Rendering all of
# them client-side with `selectInput()` / `selectizeInput(choices = ...)`
# is what triggers Shiny's
#
#   The select input "..." contains a large number of options;
#   consider using server-side selectize for massively improved
#   performance.
#
# warning, and it makes the picker visibly sluggish on typing. The
# canonical fix is to render the input once with `choices = NULL` and
# populate it from the server via `updateSelectizeInput(server = TRUE)`,
# which ships only the matches for what the user is typing right now.
#
# We wrap the two halves of that pattern in tiny helpers so module code
# stays readable and so we can apply a future tweak (e.g. a
# `maxOptions` bump) in one place.

#' Render a server-side selectize gene picker (UI half).
#'
#' Use this where you'd previously call
#' `selectizeInput(id, label, choices = available_genes(ds))` --
#' typically inside a `renderUI`. Render the input with no choices,
#' then call [update_gene_picker()] from a `shiny::observe` to fill
#' it in server-side.
#'
#' @param id           input id (no namespace; the caller wraps via `ns()`).
#' @param label        input label (may be `NULL`).
#' @param selected     initially-selected value (may be `NULL`).
#' @param multiple     allow multi-select (default `FALSE`).
#' @param placeholder  placeholder shown when empty.
#' @param ...          forwarded to `shiny::selectizeInput()`.
gene_picker_input <- function(id, label = "Gene",
                              selected    = NULL,
                              multiple    = FALSE,
                              placeholder = "Pick a gene...",
                              ...) {
  shiny::selectizeInput(
    inputId  = id,
    label    = label,
    choices  = NULL,
    selected = selected,
    multiple = multiple,
    ...,
    options  = list(placeholder = placeholder))
}

#' Populate a server-side selectize gene picker (server half).
#'
#' Wraps `updateSelectizeInput(server = TRUE)` so the gene list is
#' sent to the client incrementally as the user types instead of all
#' at once on first render. Safe to call with an empty `choices`
#' vector (the picker just shows the placeholder).
#'
#' @param session   the module's `session` (NOT `session$ns(...)`).
#' @param id        input id WITHOUT the namespace prefix; pass the
#'                  same id you used in [gene_picker_input()].
#' @param choices   character vector of gene names.
#' @param selected  optional initial selection (defaults to whatever
#'                  the input already holds).
update_gene_picker <- function(session, id, choices,
                               selected = NULL) {
  if (is.null(choices)) choices <- character()
  shiny::updateSelectizeInput(
    session  = session,
    inputId  = id,
    choices  = choices,
    selected = selected,
    server   = TRUE)
}
