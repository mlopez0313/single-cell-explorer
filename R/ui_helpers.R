# ============================================================================
# Shared UI helpers
# ----------------------------------------------------------------------------
# Small Shiny tag builders reused across modules. Keep this file dependency-
# free (Shiny only). If a helper grows specific to one module, leave it in
# that module file instead.
# ============================================================================

#' Inline yellow "friendly" warning banner.
#'
#' Used by modules to surface missing data (e.g. a gene not in the dataset)
#' without crashing the workspace.
friendly_warning <- function(text) {
  shiny::div(
    style = paste(
      "padding:8px 12px; margin-bottom:8px; border-radius:4px;",
      "background:#fff3cd; color:#664d03; border:1px solid #ffecb5;",
      "font-size:13px;"),
    text
  )
}
