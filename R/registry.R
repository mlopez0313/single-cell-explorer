# ============================================================================
# Module registry
# ----------------------------------------------------------------------------
# A module is a small, self-contained Shiny module (UI + server) that renders
# inside the main workspace. Each module declares its identity and inputs once;
# the registry below is the single source of truth used by the sidebar to draw
# navigation and by the workspace to render the active module.
#
# To register a new module:
#   1. Create R/modules/mod_<id>.R that exports `<id>_ui(id)` and
#      `<id>_server(id, state)`. See R/modules/mod_dataset_overview.R for the
#      smallest working example.
#   2. Add a `module_spec(...)` entry to `module_registry()` below.
#   3. Set `enabled = TRUE` once it's ready for users.
#
# Categories are used purely for sidebar grouping; pick one of the constants
# defined in `MODULE_CATEGORIES` or add a new category there.
# ============================================================================

MODULE_CATEGORIES <- c(
  overview    = "Overview",
  exploration = "Exploration",
  identity    = "Cell Identity",
  statistics  = "Statistics",
  functional  = "Functional Analysis",
  advanced    = "Advanced Analysis"
)

#' Describe one module.
#'
#' @param id              character(1) unique module id (snake_case)
#' @param name            character(1) display name shown in the sidebar
#' @param description     character(1) one-line summary
#' @param category        one of `MODULE_CATEGORIES`
#' @param enabled         logical(1); FALSE renders a "coming soon" placeholder
#' @param required_inputs character(); fields on app state that must be set
#'                        before the module is usable (e.g. "dataset", "assay")
#' @param ui_fn           function(id) -> Shiny tagList
#' @param server_fn       function(id, state) -> Shiny moduleServer
module_spec <- function(id, name, description, category,
                        enabled = FALSE,
                        required_inputs = character(),
                        ui_fn = NULL, server_fn = NULL) {
  stopifnot(
    is.character(id), length(id) == 1L,
    category %in% MODULE_CATEGORIES,
    is.logical(enabled), length(enabled) == 1L
  )
  list(
    id              = id,
    name            = name,
    description     = description,
    category        = category,
    enabled         = enabled,
    required_inputs = required_inputs,
    ui_fn           = ui_fn,
    server_fn       = server_fn
  )
}

#' The full module registry.
#'
#' Order here controls order in the sidebar within each category.
module_registry <- function() {
  list(
    # -- Overview ------------------------------------------------------------
    module_spec(
      id = "dataset_overview",
      name = "Dataset Overview",
      description = "Summary stats and metadata of the active dataset.",
      category = MODULE_CATEGORIES[["overview"]],
      enabled = TRUE,
      required_inputs = c("dataset"),
      ui_fn = mod_dataset_overview_ui,
      server_fn = mod_dataset_overview_server
    ),

    # -- Exploration ---------------------------------------------------------
    module_spec(
      id = "scrna_explorer",
      name = "Basic scRNA Explorer",
      description = "Inspect reductions, color cells by metadata or gene.",
      category = MODULE_CATEGORIES[["exploration"]],
      enabled = TRUE,
      required_inputs = c("dataset", "assay", "reduction"),
      ui_fn = mod_scrna_explorer_ui,
      server_fn = mod_scrna_explorer_server
    ),

    # -- Cell identity (placeholders) ---------------------------------------
    module_spec(
      id = "marker_investigation",
      name = "Marker Investigation",
      description = "Rank genes per group of a metadata field; push picks to the Explorer.",
      category = MODULE_CATEGORIES[["identity"]],
      enabled = TRUE,
      required_inputs = c("dataset", "metadata"),
      ui_fn = mod_marker_investigation_ui,
      server_fn = mod_marker_investigation_server
    ),
    module_spec(
      id = "annotation",
      name = "Cell Type Annotation",
      description = "Multi-set annotations via registered engines (manual + marker-score).",
      category = MODULE_CATEGORIES[["identity"]],
      enabled = TRUE,
      required_inputs = c("dataset", "metadata"),
      ui_fn = mod_annotation_ui,
      server_fn = mod_annotation_server
    ),

    # -- Statistics (placeholders) ------------------------------------------
    module_spec(
      id = "differential_expression",
      name = "Differential Expression",
      description = "Wilcoxon DE between two cell groups; volcano + table; results stored for downstream modules.",
      category = MODULE_CATEGORIES[["statistics"]],
      enabled = TRUE,
      required_inputs = c("dataset", "metadata"),
      ui_fn = mod_differential_expression_ui,
      server_fn = mod_differential_expression_server
    ),

    # -- Functional (placeholders) ------------------------------------------
    module_spec(
      id = "pathway_analysis",
      name = "Pathway Analysis",
      description = "Fisher's-exact ORA on DE results against built-in gene-set library.",
      category = MODULE_CATEGORIES[["functional"]],
      enabled = TRUE,
      required_inputs = c("dataset"),
      ui_fn = mod_pathway_analysis_ui,
      server_fn = mod_pathway_analysis_server
    ),

    # -- Advanced (placeholders) --------------------------------------------
    module_spec(
      id = "imputation",
      name = "Data Smoothing / Imputation",
      description = "Visualization-only smoothing (mock ALRA/MAGIC/kNN); DE/markers still use raw.",
      category = MODULE_CATEGORIES[["advanced"]],
      enabled = TRUE,
      required_inputs = c("dataset"),
      ui_fn = mod_imputation_ui,
      server_fn = mod_imputation_server
    ),
    module_spec(
      id = "trajectory",
      name = "Trajectory / Pseudotime",
      description = "Deterministic mock pseudotime from a chosen root cluster; gene-vs-time trends.",
      category = MODULE_CATEGORIES[["advanced"]],
      enabled = TRUE,
      required_inputs = c("dataset", "reduction"),
      ui_fn = mod_trajectory_ui,
      server_fn = mod_trajectory_server
    ),
    module_spec(
      id = "regulons",
      name = "Regulons / Network Analysis",
      description = "TF regulon activity per cell via AUCell + group heatmap.",
      category = MODULE_CATEGORIES[["advanced"]],
      enabled = TRUE,
      required_inputs = c("dataset"),
      ui_fn = mod_regulons_ui,
      server_fn = mod_regulons_server
    )
  )
}

#' Get a single module by id, or NULL if missing.
get_module <- function(id, registry = module_registry()) {
  for (m in registry) if (identical(m$id, id)) return(m)
  NULL
}

#' Group modules by category, preserving registration order within each group.
modules_by_category <- function(registry = module_registry()) {
  split(registry, vapply(registry, `[[`, character(1), "category"))[
    intersect(MODULE_CATEGORIES, vapply(registry, `[[`, character(1), "category"))
  ]
}

#' Check whether a module's required inputs are satisfied by the current state.
#'
#' Returns TRUE/FALSE. The mapping from input name -> state field is kept here
#' so modules don't need to know about state internals.
module_inputs_ready <- function(mod, state) {
  if (length(mod$required_inputs) == 0L) return(TRUE)
  for (req in mod$required_inputs) {
    ok <- switch(req,
      "dataset"   = !is.null(state$active_dataset),
      "assay"     = !is.null(state$selected_assay)     && nzchar(state$selected_assay),
      "reduction" = !is.null(state$selected_reduction) && nzchar(state$selected_reduction),
      "metadata"  = !is.null(state$selected_metadata_field) && nzchar(state$selected_metadata_field),
      TRUE  # unknown requirement -> don't block
    )
    if (!isTRUE(ok)) return(FALSE)
  }
  TRUE
}
