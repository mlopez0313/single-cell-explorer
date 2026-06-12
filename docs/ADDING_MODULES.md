# Adding a new module

Follow these steps to turn a "coming soon" placeholder into a real module, or
to add a brand-new module to the app.

## 1. Create the module file

Each module is one R file in `R/modules/` that exports two functions:

```r
# R/modules/mod_my_module.R

my_module_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h2("My Module"),
    shiny::plotOutput(ns("plot"))
  )
}

my_module_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    output$plot <- shiny::renderPlot({
      # Read from shared state. Never mutate state from inside renderers.
      ds <- state$active_dataset
      plot(seq_len(ds$n_cells), runif(ds$n_cells),
           main = sprintf("Mock plot for %s", ds$name))
    })
  })
}
```

Conventions:

- File name:        `mod_<id>.R`
- UI function:      `<id>_ui`
- Server function:  `<id>_server`
- Always namespace inputs/outputs with `shiny::NS(id)`.
- Read shared state via `state$...`; write back only for fields owned by
  this module (e.g. `state$selected_assay <- input$assay`).

## 2. Register the module

Open `R/registry.R` and add an entry to `module_registry()`:

```r
module_spec(
  id              = "my_module",
  name            = "My Module",
  description     = "What this module does, one sentence.",
  category        = MODULE_CATEGORIES[["exploration"]],
  enabled         = TRUE,
  required_inputs = c("dataset", "assay"),
  ui_fn           = my_module_ui,
  server_fn       = my_module_server
)
```

`required_inputs` can include any of: `dataset`, `assay`, `reduction`,
`metadata`. The workspace will block the module with a "needs inputs"
screen until they're satisfied.

If you need a new requirement type, extend `module_inputs_ready()` in
`R/registry.R` to map it to a state field.

## 3. (Optional) Add a category

Categories control sidebar grouping. To add one, edit
`MODULE_CATEGORIES` at the top of `R/registry.R`.

## 4. Test

```r
shiny::runApp("scrna-explorer")
```

Then in the sidebar click **Load mock dataset** and your new module name.
The mock dataset is a plain list whose schema matches what real loaders
will produce, so any module that works against the mock will work against
real data once loaders are implemented.

## Shared state cheat sheet

| field                       | written by                            | read by                              |
| --------------------------- | ------------------------------------- | ------------------------------------ |
| `active_dataset`            | dataset loader                        | every module                         |
| `active_module`             | sidebar                               | workspace                            |
| `selected_assay`            | `scrna_explorer`                                                    | DE, markers, imputation, ...      |
| `selected_reduction`        | `scrna_explorer`                                                    | explorer, trajectory              |
| `selected_metadata_field`   | `scrna_explorer`, `marker_investigation`, `annotation`              | explorer, markers, DE, annotation |
| `selected_gene`             | `scrna_explorer`, `marker_investigation`, `differential_expression`, `pathway_analysis` | explorer, markers, pathway, DE |
| `selected_cells`            | `scrna_explorer` (brush) + future modules                           | DE, pathway, annotation           |
| `annotation_sets`           | `annotation`                                                        | every module via `get_active_annotation(state)` ONLY |
| `active_annotation_id`      | `annotation`                                                        | every module via `get_active_annotation(state)` ONLY |
| `marker_registry`           | `set_active_dataset()` (loads default)                              | annotation engines (`marker_score`, future `sctype` / `ucell` etc.) |
| `analysis_results$de`         | `differential_expression`                                         | DE, pathway                       |
| `analysis_results$pathway`    | `pathway_analysis`                                                | pathway, downstream               |
| `analysis_results$imputation` | `imputation`                                                      | explorer (via display-mode helper); **never** DE/markers/pathway |
| `display_mode_imputation`     | `scrna_explorer` (toggle), `imputation` (Clear resets to "raw")   | explorer FeaturePlot only         |
| `analysis_results$trajectory` | `trajectory`                                                      | trajectory module only; **never** DE/markers/pathway/imputation  |
| `analysis_results$markers`    | `marker_investigation`                                            | marker module (live); future tools may stamp-detect stale runs   |
| `messages`                  | `push_message(state, text, level)`    | workspace banner                     |

### Annotation system

The annotation surface is engine-driven and supports multiple coexisting
sets. Three rules every downstream module follows:

1. **Read annotations only through `get_active_annotation(state)`.**
   Never index `state$annotation_sets` directly. If the schema bumps
   to v2, the helper migrates on read so callers stay stable.
2. **Wrap every analysis result with an `annotation_stamp`.** Every
   `state$analysis_results$<slot>` entry has an `annotation_stamp` field
   built by `make_annotation_stamp(state)`. Use it in the same line as
   you write the slot's `status / results / params` fields. Detect
   staleness with `is_result_stale(result, state)`.
3. **`apply_annotations_to_dataset()` is the only sanctioned path** to
   turn an annotation set into a dataset metadata column. It creates
   `annotation__<set_id>__<YYYY_MM_DD>`, refuses generic `cell_type`,
   and never overwrites existing columns.
4. **`apply_pseudotime_to_dataset()` is the equivalent for trajectory
   results.** It creates `pseudotime__<source>__<YYYY_MM_DD>` (always)
   and `pseudotime_bin__<source>__<YYYY_MM_DD>` (when `bins > 0`).
   Same anti-clobber rules apply; never modify trajectory metadata
   columns by hand.
5. **Regulon scores stay in `state$analysis_results$regulons`.** A
   future `apply_regulon_scores_to_dataset()` will mirror the
   annotation / pseudotime apply helpers; until then no module should
   bake AUC values into `dataset$cell_data`.

See `R/marker_registry.R`, `R/annotation_schema.R`,
`R/annotation_registry.R`, and `R/annotation.R` for the full surface.
Reference-based backends live next to the registry as
`R/annotation_singler.R`, `R/annotation_azimuth.R`, and
`R/annotation_celltypist.R`; each ships a pure
`.<engine>_to_engine_output()` converter so the schema mapping is
testable without the heavy dependency installed (the `singler` pattern
is the reference; `azimuth` and `celltypist` mirror it).

Treat this table as authoritative — if a new module needs a field not here,
add it to `new_app_state()` in `R/state.R` and document it here.

## Reading dataset content

Modules should not poke at dataset internals. Use the helpers in
`R/dataset_helpers.R`:

```r
available_assays(state$active_dataset)
available_reductions(state$active_dataset)
available_metadata_fields(state$active_dataset)
available_genes(state$active_dataset)

get_embedding(state$active_dataset, state$selected_reduction)
# -> data.frame(cell, x, y) or NULL

get_metadata(state$active_dataset, state$selected_metadata_field)
# -> vector or NULL

get_gene_expression(state$active_dataset, state$selected_gene)
# -> numeric vector or NULL

# Optional `layer` argument routes through the expression backend.
# The mock dataset and every real loader expose two layers:
#   "data"   -- log-normalised expression (default; cell-level DE / markers / pathway)
#   "counts" -- raw counts (pseudobulk DE; aggregate_pseudobulk() defaults here)
get_gene_expression(state$active_dataset, "CD3D", layer = "counts")

validate_gene(state$active_dataset, "CD3D")        # TRUE/FALSE
validate_metadata(state$active_dataset, "cluster") # TRUE/FALSE
```

When a helper returns `NULL` your module should render a friendly warning
(see `friendly_warning()` in `R/modules/mod_scrna_explorer.R`) rather than
erroring.

> **Never** access `dataset$expression[[gene]]` directly. The expression
> field is an `expression_backend` object (see
> [`R/expression_backend.R`](../R/expression_backend.R)); future loaders
> will swap the in-memory backend for sparse / HDF5 / lazy backends and
> direct indexing will silently return `NULL`. Always go through
> `get_gene_expression()`.

## Plotting

Base-R helpers in `R/plotting.R` (no extra dependencies):

- `plot_embedding_categorical(emb, values, title, xlab, ylab)` -- DimPlot style
- `plot_embedding_continuous(emb, values, title, xlab, ylab, legend_title)` -- FeaturePlot style
- `plot_volcano(de_df, ...)` -- DE volcano (used by `mod_differential_expression`)
- `plot_expression_by_group(values, groups, gene, ...)` -- box + jitter
- `plot_pathway_enrichment(df, ...)` -- horizontal bars for ORA results
- `plot_gene_vs_pseudotime(pt, expr, gene_name, n_bins, ...)` -- scatter + binned trend

`emb` arguments are always the data.frame from `get_embedding()`.
