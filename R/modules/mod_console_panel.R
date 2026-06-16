# ============================================================================
# Module: Console Panel
# ----------------------------------------------------------------------------
# Floating, collapsible right-side drawer that live-tails whichever log file
# in `sce_log_dir()` is currently most active. Surfaces the same content the
# user can `tail -f` from a terminal -- install logs, demo-build logs,
# marker-compute logs, etc. -- without leaving the browser.
#
# Why it's a Bootstrap offcanvas rather than a panel that takes layout
# space: most modules are width-greedy (Explorer's plots, Marker table,
# annotation tables) and any inline split would shrink them noticeably.
# The offcanvas is purely overlayed; closed is the default and it adds zero
# layout impact when hidden.
#
# Polling:
#   - File list refreshed every 5 s (cheap: just stat()s under the log dir).
#   - Current-file contents re-read every 1.5 s using `tail` (with a
#     readLines fallback for portability). Cap output at ~500 lines so the
#     payload over the websocket stays small even for huge install logs.
#
# Auto-scroll:
#   A tiny inline JS handler scrolls the body to the bottom every time the
#   server pushes a fresh chunk; intent matches what `tail -f` does in a
#   terminal.
# ============================================================================

# Returns the last `n` lines of a file efficiently, using POSIX `tail` when
# available. Falls back to readLines() (full read + tail in R) on platforms
# where `tail` isn't on PATH -- e.g. fresh Windows installs.
.console_tail_file <- function(path, n = 500L) {
  if (is.null(path) || !nzchar(path) || !file.exists(path))
    return(character())
  out <- tryCatch(
    suppressWarnings(system2("tail", c("-n", as.character(n), shQuote(path)),
                              stdout = TRUE, stderr = FALSE)),
    error = function(e) NULL)
  if (!is.null(out) && length(out) > 0L) return(out)
  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  utils::tail(lines, n)
}

# Sorted, newest-first list of log files in the log dir. Returns
# character(0) when the dir doesn't exist yet (no logged operation has
# run in this session). `dir` is injectable for testing.
.console_list_log_files <- function(dir = NULL) {
  d <- dir %||% tryCatch(sce_log_dir(), error = function(e) NULL)
  if (is.null(d) || !dir.exists(d)) return(character())
  fs <- list.files(d, pattern = "\\.log$", full.names = TRUE)
  if (!length(fs)) return(character())
  fs[order(file.mtime(fs), decreasing = TRUE)]
}

mod_console_panel_ui <- function(id) {
  ns <- shiny::NS(id)
  panel_id <- ns("panel")
  shiny::tagList(
    # ---- Floating toggle (vertical tab on the right edge) ---------------
    shiny::tags$button(
      class = "app-console-toggle",
      type  = "button",
      `data-bs-toggle` = "offcanvas",
      `data-bs-target` = paste0("#", panel_id),
      `aria-controls`  = panel_id,
      `aria-label`     = "Open console",
      title            = "Open console (live log tail)",
      shiny::tags$span(class = "app-console-toggle__icon",
                       app_icon("terminal", class = "")),
      shiny::tags$span(class = "app-console-toggle__label", "Console")
    ),

    # ---- Offcanvas drawer ----------------------------------------------
    shiny::div(
      id = panel_id, tabindex = "-1",
      class = "offcanvas offcanvas-end app-console",
      `aria-labelledby` = ns("title"),
      `data-bs-backdrop` = "false",   # allow interacting with app while open
      `data-bs-scroll`   = "true",
      shiny::div(class = "offcanvas-header app-console__header",
        shiny::h2(id = ns("title"),
                  class = "offcanvas-title app-console__title",
                  app_icon("terminal", class = ""), "Console"),
        shiny::div(class = "app-console__controls",
          shiny::uiOutput(ns("file_picker_ui"), inline = TRUE),
          shiny::actionButton(ns("clear"), label = NULL,
                              icon  = app_icon("eraser-fill", class = ""),
                              class = "btn btn-sm btn-outline-light app-console__btn",
                              title = "Clear (this view only -- log on disk is kept)")
        ),
        shiny::tags$button(type = "button",
                           class = "btn-close btn-close-white app-console__close",
                           `data-bs-dismiss` = "offcanvas",
                           `aria-label` = "Close console")
      ),
      shiny::div(class = "app-console__meta",
                 shiny::textOutput(ns("meta"), inline = TRUE)),
      shiny::div(class = "offcanvas-body app-console__body",
                 shiny::tags$pre(class = "app-console__pre",
                                 shiny::textOutput(ns("body"), inline = FALSE)))
    ),

    # ---- Auto-scroll: scroll body to bottom on every render -----------
    # No external JS file -- a 6-line inline handler is plenty and keeps
    # the module self-contained.
    shiny::tags$script(shiny::HTML(sprintf(
      "Shiny.addCustomMessageHandler('sce_console_scroll', function(panelId) {
         var p = document.getElementById(panelId);
         if (!p) return;
         var body = p.querySelector('.app-console__body');
         if (body) body.scrollTop = body.scrollHeight;
       });", panel_id)))
  )
}

mod_console_panel_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    panel_id <- ns("panel")

    # ---- Log-dir watcher (cheap, polled every 5 s) ---------------------
    # Uses sum(mtime) as a change signal; any add/modify/delete invalidates.
    log_files <- shiny::reactivePoll(5000, session,
      checkFunc = function() {
        fs <- .console_list_log_files()
        if (!length(fs)) return(0)
        sum(as.numeric(file.mtime(fs)))
      },
      valueFunc = .console_list_log_files)

    # State for "cleared" view: the user can click eraser to hide the
    # current contents until the log grows further. We compare current
    # tail length to the value at clear-time.
    cleared_n <- shiny::reactiveVal(NA_integer_)
    shiny::observeEvent(input$clear, {
      lines <- .console_tail_file(chosen(), n = 500L)
      cleared_n(length(lines))
    })
    # Resetting `cleared_n` when the user picks a different file makes the
    # eraser scoped to the file you erased -- behaving like `clear` in a
    # terminal that's tailing a particular path.
    shiny::observeEvent(input$which_file, {
      cleared_n(NA_integer_)
    }, ignoreInit = TRUE)

    output$file_picker_ui <- shiny::renderUI({
      fs <- log_files()
      if (!length(fs)) {
        return(shiny::tags$em(class = "app-console__hint",
                              "No logs yet."))
      }
      shiny::selectInput(
        ns("which_file"),
        label    = NULL,
        choices  = stats::setNames(fs, basename(fs)),
        selected = isolate_or_default(input$which_file, fs[1L], fs),
        width    = "260px")
    })

    chosen <- shiny::reactive({
      fs <- log_files()
      if (!length(fs)) return(NA_character_)
      pick <- input$which_file %||% fs[1L]
      if (!pick %in% fs) pick <- fs[1L]
      pick
    })

    # ---- Live tail of the chosen file (polled every 1.5 s) -------------
    body_text <- shiny::reactivePoll(1500, session,
      checkFunc = function() {
        p <- chosen()
        if (is.na(p) || !file.exists(p)) return(0)
        # Combine mtime + size so a same-second append still invalidates.
        info <- file.info(p)
        sprintf("%s|%s", info$mtime, info$size)
      },
      valueFunc = function() {
        p <- chosen()
        if (is.na(p)) return("(no logs yet -- run an install, build, or compute)")
        if (!file.exists(p)) return(sprintf("(file disappeared: %s)", p))
        lines <- .console_tail_file(p, n = 500L)
        cn <- cleared_n()
        if (!is.na(cn) && length(lines) > cn) {
          lines <- utils::tail(lines, length(lines) - cn)
        } else if (!is.na(cn)) {
          lines <- character()
        }
        if (!length(lines)) return("(no new output since last clear)")
        paste(lines, collapse = "\n")
      })

    output$body <- shiny::renderText({ body_text() })

    # Header strip: filename + age + size, so the user can tell at a
    # glance which file they're looking at and whether it's still moving.
    output$meta <- shiny::renderText({
      p <- chosen()
      if (is.na(p) || !file.exists(p)) return("")
      info <- file.info(p)
      age_s <- as.numeric(difftime(Sys.time(), info$mtime, units = "secs"))
      age_label <- if (age_s < 60)        sprintf("updated %.0fs ago",      age_s)
                   else if (age_s < 3600) sprintf("updated %.0f min ago",   age_s / 60)
                   else                    sprintf("updated %.1f h ago",     age_s / 3600)
      sprintf("%s | %s | %s",
              basename(p),
              format(structure(info$size, class = "object_size"),
                     units = "auto"),
              age_label)
    })

    # ---- Auto-scroll to bottom on each refresh -------------------------
    shiny::observe({
      body_text()  # take dep
      session$sendCustomMessage("sce_console_scroll", panel_id)
    })
  })
}

# Helper: return `current` if it's in `valid`, else `default`. Used for
# preserving the user's file selection across log-dir refreshes that
# don't drop their current file.
isolate_or_default <- function(current, default, valid) {
  if (is.null(current) || !nzchar(current %||% "")) return(default)
  if (!current %in% valid)                          return(default)
  current
}
