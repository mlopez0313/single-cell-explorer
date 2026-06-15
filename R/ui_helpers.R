# ============================================================================
# Shared UI helpers (legacy shims)
# ----------------------------------------------------------------------------
# Reusable presentation primitives now live in `R/ui_components.R`. This file
# keeps a thin compatibility shim so older modules that still call
# `friendly_warning("...")` automatically pick up the new design system.
# ============================================================================

#' Inline friendly warning banner used by modules to surface missing data.
#'
#' Delegates to `info_banner()` (see `R/ui_components.R`) so any call site
#' picks up the modern styling without changing its API.
friendly_warning <- function(text) {
  info_banner(text, tone = "warning")
}
