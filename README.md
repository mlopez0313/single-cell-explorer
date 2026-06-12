# scRNA Explorer

A modular R + Shiny application for interactive single-cell RNA-seq
exploration. The repository is no longer a shell with stubs: every
analysis surface listed below is implemented and tested end-to-end.

**Implemented modules** (driven by `R/registry.R` -> a single
`module_registry`; each module is independently testable via
`shiny::testServer` and pluggable through its own registry):

- **Explorer** -- gene / metadata plots over any registered reduction.
- **Differential expression** -- cell-level (pure-R Wilcoxon, optional
  presto) and pseudobulk (naive t-test, optional edgeR / DESeq2), with
  layer-aware gene universe.
- **Markers** -- per-cluster marker ranking.
- **Pathway analysis** -- ORA + GSEA scaffolding against pluggable
  pathway sources (built-in mock collection, optional `msigdbr`).
- **Annotation** -- multi-set / freeze / apply lifecycle, with engines:
  `manual`, `marker_score` (registry-driven), `singler`, `azimuth`,
  `celltypist`. Every engine stamps explicit, distinct
  `engine_version` and `schema_version` provenance.
- **Imputation** -- per-gene smoothing on KNN graphs.
- **Trajectory** -- pseudotime via a pluggable registry (`mock`,
  `metadata`, optional `slingshot`, optional `monocle3`); results can
  be baked into dataset metadata via `apply_pseudotime_to_dataset()`.
- **Regulons** -- per-cell TF activity scoring (pure-R AUCell + optional
  Bioconductor AUCell) over pluggable regulon sources (built-in
  mock, optional `dorothea`).

**Loaders**: real `Seurat`, `10x` Cellranger, and lazy `AnnData` (`.h5ad`
via `rhdf5` when present, with `zellkonverter` / `anndata` fallbacks).
Every loader feeds the same flat dataset schema; modules don't know which
loader ran.

**Optional dependencies** are all gated through `require_optional()`:
missing-dep errors print the exact install command (CRAN /
`BiocManager::install` / `remotes::install_github("owner/repo")`).

The mock dataset (`mock_dataset()`) exists so the app is fully usable
with no external data and so every test runs with no optional deps.

## Run it

```bash
R -e 'shiny::runApp("scrna-explorer", port = 3838, host = "0.0.0.0", launch.browser = FALSE)'
```

Or from R:

```r
shiny::runApp("scrna-explorer")
```

Use the sidebar **Dataset** section to load data:

- **Load mock dataset** — in-memory demo data (no files required).
- **Load dataset** — enter a server-side path and choose a source
  (`Auto-detect`, Seurat `.rds`, AnnData `.h5ad`, or a 10x Cellranger
  directory). Load errors appear in the sidebar and the message log without
  replacing the current dataset.
- **CELLxGENE Census (experimental)** — available from the dataset source
  picker in the sidebar. Requires `cellxgene.census` and `SeuratObject`
  plus network access. You can browse studies for a selected organism,
  pick one, and optionally add an extra SOMA cell filter to further narrow
  the selected study before loading. Treat this as experimental app-level support.

### Requirements

The app shell itself needs only base R + Shiny:

```r
install.packages("shiny")
```

Real dataset loaders are optional dependencies — install only the packages
needed by the format(s) you actually load:

| input format          | required packages                                            | install |
| --------------------- | ------------------------------------------------------------ | ------- |
| **Seurat `.rds`**     | `SeuratObject` (read-side); `Seurat` (for create/inspect)    | `install.packages(c("SeuratObject", "Seurat"))` |
| **10x Cellranger dir**| `Matrix` (no `Seurat` requirement)                           | `install.packages("Matrix")` |
| **AnnData `.h5ad`**   | **preferred (lazy, gene-by-gene reads):** `rhdf5` -- builds an `expression_backend_h5ad` that never materialises the full matrix | `BiocManager::install("rhdf5")` |
|                       | eager fallback #1: `zellkonverter` + `SingleCellExperiment` -- slurps the file into an SCE in memory | `BiocManager::install(c("zellkonverter", "SingleCellExperiment"))` |
|                       | eager fallback #2: `anndata` (CRAN) -- materialises the matrix via the CRAN `anndata` package | `install.packages("anndata")` |

Each loader fails with a clear "install X" message if its dependency is
missing — no NULL-pointer crashes.

### Loading data programmatically

```r
ds <- load_dataset("/path/to/object.rds")              # Seurat
ds <- load_dataset("/path/to/10x_filtered_bc_matrix/") # 10x directory
ds <- load_dataset("/path/to/data.h5ad")               # AnnData
ds <- load_dataset("/cellxgene-census", source = "census",
                   organism = "Homo sapiens",
                   obs_value_filter = "cell_type == 'B cell'")
ds <- mock_dataset()                                   # synthetic demo

# All `ds` objects satisfy `dataset_schema()` and feed every module.
```

`detect_source()` infers the source from the path (`.rds` → seurat,
`.h5ad` → anndata, directory → 10x); pass `source = "..."` explicitly to
override.

### Source → schema mapping

The three loaders normalise their sources into the same flat schema:

| schema field        | Seurat (`.rds`)                          | 10x directory                          | AnnData (`.h5ad`)                              |
| ------------------- | ---------------------------------------- | -------------------------------------- | ---------------------------------------------- |
| `assays`            | `SeuratObject::Assays(obj)`              | `"RNA"`                                | `"RNA"`                                        |
| `default_assay`     | `DefaultAssay(obj)`                      | `"RNA"`                                | `"RNA"`                                        |
| `reductions`        | uppercased `Reductions(obj)`             | `character()` (none in raw 10x)        | `obsm` / `reducedDims` (strips `X_` prefix)    |
| `cells`             | `colnames(obj)`                          | `barcodes.tsv`                         | `colnames(sce)` / `ad$obs_names`               |
| `cell_data`         | `obj@meta.data` + embedding columns      | barcodes + `n_counts` / `n_features`   | `colData(sce)` / `ad$obs` + embedding columns  |
| `genes`             | `rownames(default_assay)`                | feature symbols (column 2 by default; duplicates resolved via `make.unique`) | `rownames(sce)` / `ad$var_names` |
| `expression`        | `expression_backend_sparse({data, counts})` | `expression_backend_sparse({counts})`  | `expression_backend_sparse({data[, ...]})`     |
| `default_reduction` | uppercased first reduction               | `NA_character_`                        | first `obsm` key (after `X_` strip)            |

Per-reduction embedding columns are surfaced as `<NAME>_1` / `<NAME>_2` in
`cell_data` so every module's existing `<NAME>_1`/`<NAME>_2` lookup in
`get_embedding()` keeps working unchanged.

### Running the tests

The app ships with a `testthat`-based regression suite covering state,
the module registry, mock dataset + helpers, the marker registry, the
annotation system (schema + engines + multi-set + apply), trajectory
(method registry + mock/metadata + slingshot/monocle3 converters and
missing-dep gating), imputation, DE (cell-level + pseudobulk dispatcher
+ naive pseudobulk end-to-end), markers, pathway analysis (built-in +
msigdbr converter + GSEA scaffold), real dataset loaders (10x + Seurat
end-to-end; AnnData lazy via `rhdf5` end-to-end when installed, and the
eager `zellkonverter` round-trip when that's installed -- both skipped
otherwise), the in-memory + sparse + lazy-h5ad expression backends, and
the cross-cutting "every module stamps results with annotation
provenance" contract.

Two integration tests give the loaders true end-to-end coverage against
realistically structured data (no remote downloads):

- `tests/testthat/test-integration-seurat.R` builds a 3-population
  mini-PBMC (600 cells × 1.2 k genes), runs the full Seurat 5 pipeline
  (NormalizeData → HVG → ScaleData → PCA → FindClusters → UMAP),
  saves to `.rds`, and round-trips through `load_seurat()`. From the
  reconstructed dataset it then drives `compute_markers()`,
  `compute_de(..., backend = "wilcox_r")`, `compute_enrichment()`,
  `run_annotation_engine("marker_score")`,
  `run_regulon_engine("aucell_r")`, `run_trajectory("mock")`, and the
  `apply_pseudotime_to_dataset()` / `apply_annotations_to_dataset()`
  bakers. Skipped when `Seurat` / `SeuratObject` are not installed.

- `tests/testthat/test-integration-10x.R` writes a Cellranger v3
  layout (`matrix.mtx`, `barcodes.tsv`, `features.tsv`, both gzipped
  and plain variants) from scratch via `Matrix::writeMM()`, then
  exercises `load_10x()` against it.

```bash
# from the project root
Rscript tests/testthat.R
```

or, equivalently, from an interactive R session:

```r
testthat::test_dir("tests/testthat", reporter = "summary")
```

Tests source the app's `R/` files directly via `tests/testthat/helper-app.R`,
so there is no package install step. `testthat` is the only extra
dependency required to run them (`install.packages("testthat")`).

## Layout

```
scrna-explorer/
├── app.R                       # Entry point. Sources everything and builds the app.
├── R/
│   ├── state.R                 # new_app_state(): the shared reactiveValues object
│   ├── dataset.R               # dataset_schema(), mock_dataset(), real Seurat/AnnData/10x loaders
│   ├── dataset_helpers.R       # available_*, get_embedding, get_gene_expression, validators
│   ├── markers.R               # compute_markers(): one-vs-rest (Wilcoxon default)
│   ├── de.R                    # compute_de(): group vs group (Wilcoxon + BH); filter + sort helpers
│   ├── pathway.R               # compute_enrichment() (Fisher's + BH) + compute_gsea() scaffold
│   ├── pathway_sources.R       # PATHWAY_SOURCES registry: builtin + msigdbr
│   ├── annotation.R            # set management + dispatcher + apply_annotations_to_dataset()
│   ├── annotation_schema.R     # annotation_result_v1() + hash + cluster->cell expand
│   ├── annotation_registry.R   # ANNOTATION_ENGINES (manual + marker_score + singler + azimuth + celltypist)
│   ├── annotation_singler.R    # SingleR engine impl + pure schema converter
│   ├── annotation_azimuth.R    # Azimuth engine impl + pure schema converter
│   ├── annotation_celltypist.R # CellTypist engine impl (reticulate) + pure converter
│   ├── marker_registry.R       # typed/versioned marker registry (data, not code)
│   ├── expression_backend.R    # S3 expression_backend (in-memory + sparse)
│   ├── expression_backend_h5ad.R # Lazy HDF5 backend (csr / csc / dense)
│   ├── pseudobulk.R            # aggregate_pseudobulk + validate + naive backend helpers
│   ├── optional_deps.R         # has_optional / require_optional gating
│   ├── imputation.R            # mock smoothers + get_gene_expression_for_view()
│   ├── trajectory.R            # run_trajectory() dispatcher + mock/metadata run_fns + binning helpers
│   ├── trajectory_registry.R   # TRAJECTORY_METHODS registry + spec
│   ├── trajectory_slingshot.R  # Slingshot backend + pure converter
│   ├── trajectory_monocle3.R   # Monocle3 backend + pure converter
│   ├── regulon_schema.R        # regulon_spec / regulon_set / regulon_result_v1
│   ├── regulon_registry.R      # REGULON_ENGINES + dispatcher run_regulon_engine
│   ├── regulon_aucell.R        # Pure-R AUCell + Bioc AUCell + pure converter
│   ├── regulon_sources.R       # REGULON_SOURCES: mock_pbmc + dorothea
│   ├── plotting.R              # embedding, volcano, expression-by-group, pathway bars, gene-vs-pt
│   ├── ui_helpers.R            # friendly_warning() and other shared tag builders
│   ├── registry.R              # module_spec() + module_registry()
│   ├── ui_sidebar.R            # Sidebar navigation grouped by category
│   ├── ui_workspace.R          # Main content area; empty / coming-soon / module UI + active-set banner
│   └── modules/
│       ├── mod_dataset_overview.R           # ENABLED  (smallest reference module)
│       ├── mod_scrna_explorer.R             # ENABLED  (assay/reduction/metadata/gene + plots)
│       ├── mod_marker_investigation.R       # ENABLED  (ranks markers, stamps result)
│       ├── mod_annotation.R                 # ENABLED  (multi-set annotation, engines)
│       ├── mod_differential_expression.R    # ENABLED  (group vs group; stamps result)
│       ├── mod_pathway_analysis.R           # ENABLED  (ORA on DE; stamps result)
│       ├── mod_imputation.R                 # ENABLED  (visualization-only smoothing)
│       ├── mod_trajectory.R                 # ENABLED  (mock pseudotime + gene-vs-time)
│       ├── mod_regulons.R                   # ENABLED  (AUCell regulon scoring + heatmap)
│       └── mod_placeholders.R               # shared coming_soon_ui (no placeholders left)
├── tests/
│   ├── testthat.R                  # Runner: `Rscript tests/testthat.R`
│   └── testthat/
│       ├── helper-app.R                  # Sources R/ files into the test env
│       ├── test-state.R                  # state factory + dataset cascade + annotation stamp
│       ├── test-registry.R               # module_registry shape + enabled/disabled set
│       ├── test-dataset.R                # mock_dataset schema + dataset_helpers
│       ├── test-expression-backend.R     # in-memory + sparse backend contracts
│       ├── test-expression-backend-h5ad.R # lazy h5ad backend math + fixture round-trip
│       ├── test-pseudobulk.R              # aggregation + naive backend + missing-dep errors
│       ├── test-loaders.R                # detect_source + 10x e2e + Seurat e2e + AnnData skip-if-missing
│       ├── test-marker-registry.R        # marker registry schema + queries
│       ├── test-annotation.R             # schema + engines + multi-set + apply
│       ├── test-singler.R                # SingleR engine spec + converter + e2e skip-if-missing
│       ├── test-azimuth.R                # Azimuth engine spec + converter + e2e skip-if-missing
│       ├── test-celltypist.R             # CellTypist engine spec + converter + e2e skip-if-missing
│       ├── test-trajectory.R             # compute_pseudotime + bin (mock/metadata)
│       ├── test-trajectory-registry.R    # TRAJECTORY_METHODS registry + dispatcher
│       ├── test-trajectory-slingshot.R   # Slingshot converter + missing-dep
│       ├── test-trajectory-monocle3.R    # Monocle3 converter + missing-dep
│       ├── test-regulons.R               # schema + sources + engines + AUCell math + e2e
│       ├── test-imputation.R             # smoothers + viz-only contract
│       ├── test-de.R                     # compute_de + filter/sort
│       ├── test-markers.R                # compute_markers (Wilcoxon default)
│       ├── test-pathway.R                # compute_enrichment + compute_gsea scaffold
│       ├── test-pathway-sources.R        # PATHWAY_SOURCES registry + .msigdbr_to_pathways
│       └── test-stamping.R               # every module stamps annotation provenance
├── DESCRIPTION                    # dependency manifest for tooling/tests
├── .Rbuildignore
├── .gitignore
├── docs/
│   └── ADDING_MODULES.md
└── README.md
```

## How modules work

A **module** is a Shiny module (UI + server pair) that renders inside the
main workspace. The list of modules lives in
[`R/registry.R`](R/registry.R) and is the **single source of truth** consumed
by the sidebar and workspace.

Each entry has:

| field             | purpose                                                       |
| ----------------- | ------------------------------------------------------------- |
| `id`              | unique snake_case id, used for Shiny namespacing and routing  |
| `name`            | display name in the sidebar                                   |
| `description`     | one-line tooltip / sub-text                                   |
| `category`        | sidebar group (see `MODULE_CATEGORIES`)                       |
| `enabled`         | `FALSE` => greyed out with a "coming soon" panel              |
| `required_inputs` | strings like `"dataset"`, `"assay"`, `"reduction"`            |
| `ui_fn`           | `function(id)` returning a `tagList`                          |
| `server_fn`       | `function(id, state)` returning a `moduleServer`              |

`state` is the shared `reactiveValues` object created in
[`R/state.R`](R/state.R). It is the only object passed between modules.

## Adding a new module

See [`docs/ADDING_MODULES.md`](docs/ADDING_MODULES.md). Short version:

1. Create `R/modules/mod_<id>.R` exporting `<id>_ui(id)` and `<id>_server(id, state)`.
2. Add a `module_spec(...)` entry to `module_registry()` in `R/registry.R`.
3. Flip `enabled = TRUE` when it's ready.

## Shared app state

All cross-module state lives on a single `reactiveValues` object:

- `active_dataset`
- `active_module`
- `selected_assay`
- `selected_reduction`
- `selected_metadata_field`
- `selected_gene`
- `selected_cells`
- `messages`

Modules should never maintain their own duplicate copies of these fields —
read and write through `state`.

## Dataset structure (mock + future loaders)

Every dataset — mock today, Seurat/AnnData/10x later — must be a plain
list matching `dataset_schema()` in `R/dataset.R`. The mock dataset is
the reference implementation:

| field               | type            | notes                                                |
| ------------------- | --------------- | ---------------------------------------------------- |
| `name`              | character(1)    | display name                                         |
| `source`            | character(1)    | `"mock"` / `"seurat"` / `"anndata"` / `"10x"`        |
| `n_cells`           | integer(1)      |                                                      |
| `n_genes`           | integer(1)      |                                                      |
| `assays`            | character()     | e.g. `c("RNA", "SCT", "ADT")`                        |
| `default_assay`     | character(1)    |                                                      |
| `reductions`        | character()     | e.g. `c("PCA", "UMAP", "tSNE")`                      |
| `default_reduction` | character(1)    |                                                      |
| `metadata_fields`   | character()     | column names exposed in `cell_data`                  |
| `cells`             | character()     | cell barcodes / ids, length `n_cells`                |
| `cell_data`         | data.frame      | one row per cell. **Required columns:** `cell`, every name in `metadata_fields`, and `<RED>_1`/`<RED>_2` for every reduction in `reductions` |
| `genes`             | character()     | symbols queryable via `get_gene_expression()`        |
| `expression`        | `expression_backend` | structured backend object (see "Expression backend" below); legacy named-list shape is also accepted by the helpers |

The Basic scRNA Explorer never touches these fields directly — it goes
through `R/dataset_helpers.R` (`available_assays`, `get_embedding`,
`get_gene_expression`, etc.). Future loaders only need to produce a list
satisfying the schema above to "just work" with every module.

### Expression backend

`dataset$expression` is an S3 `expression_backend` object defined in
[`R/expression_backend.R`](R/expression_backend.R). Modules **never**
touch it directly — they call `get_gene_expression(dataset, gene)` /
`available_genes(dataset)` and the helpers dispatch to the backend.

Concrete backends today:

| backend                             | source                                 | when used                              |
| ----------------------------------- | -------------------------------------- | -------------------------------------- |
| `expression_backend_inmemory()`     | named list / list-of-lists of numeric  | mock dataset and any small in-memory data |
| `expression_backend_sparse()`       | `Matrix::dgCMatrix` (or any matrix-like with `[gene, ]` row indexing) per layer | Seurat `.rds`, 10x directories, AnnData via `zellkonverter` (`DelayedArray`s also work) |
| `expression_backend_h5ad()`         | direct HDF5 / `.h5ad` file handle via `rhdf5` | Large AnnData files where materialising n_obs \u00d7 n_var is impractical. Reads obs / var / obsm eagerly but pulls expression gene-by-gene from HDF5. |

Reserved for follow-up prompts:

| backend (planned)                | source                          | why                                  |
| -------------------------------- | ------------------------------- | ------------------------------------ |
| `expression_backend_seurat()`    | thin proxy over a Seurat object | avoid even constructing an SCE       |

#### Lazy `.h5ad` reading

`load_anndata(path, lazy = TRUE)` (the default) prefers `rhdf5` when it
is installed:

1. Reads `/obs`, `/var/_index`, and `/obsm/X_*` eagerly (metadata
   layer; small).
2. Records each expression layer's encoding (`csr_matrix`,
   `csc_matrix`, or `dense`) and shape without reading any matrix
   data.
3. Returns a dataset whose `expression` field is an
   `expression_backend_h5ad()`. Per-gene reads pull only that gene's
   sparse slice from HDF5; the first sparse read populates an
   in-process cache (`indptr`, `indices`, `data`) so subsequent reads
   on the same layer are O(nnz_per_gene) without re-hitting disk.
4. Falls back to the existing `zellkonverter` / `anndata` paths when
   `rhdf5` is absent. Pass `lazy = FALSE` to force the eager paths
   even when `rhdf5` is installed.

Supported AnnData v0.8+ encodings: `csr_matrix`, `csc_matrix`, and 2D
dense datasets. Categoricals in `/obs` (the
`encoding-type = "categorical"` group with `codes` + `categories`) are
decoded back to character columns automatically.

Limitations: the sparse triple is loaded entirely on first access for
the relevant layer, so this is best for *moderately* sparse datasets
(< nnz \u2248 200M values fits comfortably in 1 GB). For larger files
the next step is per-chunk HDF5 reads, which would extend this
backend without changing the schema.

The generic surface a new backend must implement:

```r
backend_n_cells(backend)
backend_available_layers(backend)
backend_default_layer(backend)
backend_genes(backend, layer = NULL)
backend_n_genes(backend, layer = NULL)
backend_has_gene(backend, gene, layer = NULL)
backend_get_gene(backend, gene, layer = NULL)   # -> numeric(n_cells) | NULL
```

The optional `layer` argument is already wired through
`get_gene_expression(dataset, gene, layer = NULL)` so pseudobulk DE
(Cursor #5) can later add a `"counts"` layer beside `"data"` without
touching the dataset schema or any analysis module.

Legacy compatibility: any pre-backend dataset where `expression` is still
a flat named list of numeric vectors is coerced into an in-memory backend
on read via `as_expression_backend()`. No module code change is required
for older fixtures.

## Cross-module flow (live)

Two multi-module loops are in place.

### Marker Investigation → Explorer

1. In **Marker Investigation**, pick a grouping field and click *Find markers*.
2. Pick a gene in *Highlight gene* and click *Send to Explorer*.
3. The Explorer's FeaturePlot immediately re-colors by that gene
   (it reads `state$selected_gene`, which the marker module just wrote).

The marker module also keeps `state$selected_metadata_field` in sync with
its grouping control, so the Explorer's DimPlot color also updates when
you change the grouping.

### Cell Type Annotation (multi-set, engine-driven)

The annotation system is engine-agnostic and supports multiple coexisting
sets. Files involved:

| File | Role |
| --- | --- |
| `R/marker_registry.R` | typed, versioned marker registry (data, not code) |
| `R/annotation_schema.R` | `annotation_result_v1()` + hash + cluster→cell expand |
| `R/annotation_registry.R` | `ANNOTATION_ENGINES()` + the two no-dep engines (`manual`, `marker_score`) |
| `R/annotation_singler.R` | reference-based engine: SingleR + celldex |
| `R/annotation_azimuth.R` | reference-based engine: Azimuth (Seurat) |
| `R/annotation_celltypist.R` | reference-based engine: CellTypist via reticulate |
| `R/annotation.R` | set management (add/dup/rename/freeze/delete) + dispatcher + `apply_annotations_to_dataset()` |
| `R/modules/mod_annotation.R` | the module UI (set manager + engine controls) |

#### Workflow

1. Open **Cell Type Annotation**.
2. Click **+ New set** — creates an empty annotation set named "Set HH:MM:SS".
3. Pick an **Annotation engine**:
   - **Manual annotation** — edit the per-cluster table, click *Apply
     table edits*. Per-cluster edits expand to per-cell labels.
   - **Marker-score (registry-driven)** — set species/min-score, click
     *Run engine*. Cell-type marker panels from `state$marker_registry`
     are scored against each cluster; top label per cluster is written
     to the set (per-cell labels are stored, not per-cluster).
4. **Apply to dataset metadata** writes a provenance-named column
   `annotation__<set_id>__<YYYY_MM_DD>` onto `dataset$cell_data` and
   registers it in `metadata_fields`. The Explorer's Color-by-metadata
   picker will list it. Generic `cell_type` is refused; existing columns
   are never overwritten.
5. **Download CSV** exports the active set.

#### Multi-set operations

- Top of the module: dropdown lists every set; selecting one activates it.
- **Duplicate** copies the active set under a new id (with `parent_set_id`
  pointing to the source — useful for "manuscript" vs "revision" lineage).
- **Rename** / **Toggle freeze** / **Delete** — frozen sets refuse renames
  and deletes; useful for locking the version you submitted.

#### Annotation set schema (`annotation_result_v1`)

Every set is an `annotation_result_v1` object:

```r
list(
  schema_version          = "annotation_v1",
  set_id, name, description,
  engine_id, engine_version, params,
  cell                    = character(n_cells),
  cell_labels             = character(n_cells),     # per-cell labels (NA allowed)
  cell_scores             = numeric(n_cells),       # 0 for unlabelled, 1 for manual confirmed
  alt_labels              = NULL | data.frame(cluster, rank, label, score),
  cluster_summary         = NULL | data.frame(cluster, top_label, top_score, ...),
  ontology_map            = NULL | named character: label -> Cell Ontology id,
  reference_source        = "user" | "marker_registry/builtin_v0.1.0" | ...,
  marker_registry_version = "builtin_v0.1.0" | NA,
  parent_set_id, cluster_field_used, n_clusters_at_creation,
  is_frozen, is_demo, warnings,
  created_at, modified_at, timestamp, duration_ms, error_message,
  edit_history
)
```

#### State layout

| state field | type | written by | read by |
| --- | --- | --- | --- |
| `state$annotation_sets` | named list of `annotation_result_v1` | annotation module | `get_active_annotation(state)` |
| `state$active_annotation_id` | character(1) or NULL | annotation module | every module via `get_active_annotation(state)` |
| `state$marker_registry` | typed registry list | `set_active_dataset()` | annotation engines, future SingleR/scType/CellTypist |

`get_active_annotation(state)` is the ONLY supported read path for
downstream modules. They MUST NOT index `state$annotation_sets` directly
— if the schema bumps to v2, that helper migrates on read so callers
stay stable.

#### Provenance stamping on analysis results

Every analysis-result slot in `state$analysis_results` (`$de`, `$pathway`,
`$imputation`, `$trajectory`, `$markers`) now carries an
`annotation_stamp` field:

```r
list(
  annotation_set_id_used   = "<set_id>" | NA,
  annotation_set_hash_used = "<content hash>" | NA,
  annotation_set_name      = "<display name>" | NA,
  annotation_engine_id     = "<engine>" | NA,
  annotation_set_is_demo   = logical(1),
  stamped_at               = POSIXct(1)
)
```

The stamp is built by `make_annotation_stamp(state)` at the moment a
result is wrapped — compute functions stay pure. `is_result_stale(result,
state)` returns TRUE when the active set has changed (different id or
different content hash) since the result was produced.

#### Workspace banner

`R/ui_workspace.R` renders a persistent active-annotation banner above
every module: dataset name on the right, active set summary on the left.
When no active set exists, the banner says "(no active annotation)" in
grey; demo sets show in amber.

#### Built-in engines

| engine          | category          | what it does                                                          | dependencies                                            |
| --------------- | ----------------- | --------------------------------------------------------------------- | ------------------------------------------------------- |
| `manual`        | manual            | User-supplied per-cluster labels, expanded to per-cell.                | none                                                    |
| `marker_score`  | marker-score      | Weighted sum of marker-registry genes per cluster; top label per cluster. | none                                                    |
| `singler`       | reference-based   | `SingleR::SingleR()` against a `celldex` reference (HPCA / BlueprintEncode / MonacoImmune). Per-cell or per-cluster mode. | `BiocManager::install(c("SingleR", "celldex", "SingleCellExperiment", "SummarizedExperiment"))` |
| `azimuth`       | reference-based   | `Azimuth::RunAzimuth()` against a tissue-specific Seurat reference (PBMC / lung / kidney / bone marrow / heart / fetus). Per-cell predictions at three label levels. | `remotes::install_github("satijalab/azimuth")` plus the reference data package (e.g. `Azimuth.pbmcref`) |
| `celltypist`    | reference-based   | CellTypist (Python) bridged through `reticulate` + `anndata`. 50+ pretrained models (immune, lung, gut, brain, dev...) with optional cluster-aware majority voting. | `install.packages(c("reticulate", "anndata"))` + `pip install celltypist` in the active Python env |

##### `singler` parameters

| param            | type             | notes                                                                                        |
| ---------------- | ---------------- | -------------------------------------------------------------------------------------------- |
| `reference`      | select (default `hpca`)         | One of `hpca`, `blueprint_encode`, `monaco_immune`. Fetched via `celldex` and cached per session. |
| `labels`         | select (default `main`)         | `main` (broad) or `fine` (high-resolution) reference labels.                                  |
| `cluster_field`  | metadata field (optional)       | If set, SingleR runs in cluster mode (one prediction per cluster, expanded back to per-cell). Much faster on > 50k cells. |
| `min_delta`      | numeric (default `0`)           | Minimum `delta.next` confidence margin. Cells below this become `"Unknown"` (score 0).        |

Per-cell mode is `O(n_cells × n_genes)`; cluster mode is
`O(n_clusters × n_genes)` and recommended for any production-sized
dataset. Both modes produce identical schema (`annotation_result_v1`)
with `cell_labels` and `cell_scores` aligned to `dataset$cell_data$cell`,
plus a top-3 `alt_labels` frame derived from SingleR's per-row scores
matrix.

##### `azimuth` parameters

| param                | type                            | notes                                                                                              |
| -------------------- | ------------------------------- | -------------------------------------------------------------------------------------------------- |
| `reference`          | select (default `pbmcref`)      | One of `pbmcref`, `lungref`, `kidneyref`, `bonemarrowref`, `heartref`, `fetusref`. Each ships as a separate data package on the Satija lab server and downloads several GB on first use. |
| `annotation_level`   | select (default `celltype.l2`)  | `celltype.l1`, `l2`, or `l3` — granularity depends on which the chosen reference exposes.            |
| `cluster_field`      | metadata field (optional)       | Only used to build a per-cluster summary table. Azimuth predictions are always per-cell.            |
| `min_mapping_score`  | numeric (default `0`)           | Cells with `mapping.score` below this become `"Unknown"` (score 0). Azimuth's mapping confidence is independent of the label probability. |

Internally we wrap the active dataset in a minimal Seurat object built
from the `counts` layer and call `Azimuth::RunAzimuth(query, reference)`.
The pure converter `.azimuth_to_engine_output()` accepts any data.frame
shaped like Azimuth's `meta.data` output, which means the schema layer
is regression-tested without Azimuth or its multi-GB reference data
packages installed.

##### `celltypist` parameters

| param              | type                                      | notes                                                                                         |
| ------------------ | ----------------------------------------- | --------------------------------------------------------------------------------------------- |
| `model`            | select (default `Immune_All_Low.pkl`)     | Any model name from the CellTypist catalogue (celltypist.org). Auto-downloads to `~/.celltypist/models/`. |
| `majority_voting`  | logical (default `FALSE`)                 | If TRUE, CellTypist aggregates predictions over a cluster (recommended for noisy / sparse data). |
| `over_clustering`  | metadata field (optional)                 | Column to use as the majority-voting reference. Only consulted when `majority_voting = TRUE`.   |
| `cluster_field`    | metadata field (optional)                 | Only used to build a per-cluster summary table; independent of `over_clustering`.               |
| `min_score`        | numeric (default `0`)                     | Cells with `conf_score` below this become `"Unknown"` (score 0).                                |

Internally we build an AnnData object from the `data` (log-normalised)
layer, run `celltypist.annotate(adata, model=...)` via reticulate, and
route the resulting `predicted_labels` DataFrame through the pure
converter `.celltypist_to_engine_output()`. That same converter is
covered by hand-built data.frame tests, so the schema layer remains
testable without a working Python environment.

#### Adding a new annotation engine

1. Implement a `.run_<id>_annotation(dataset, params, state)` function
   that returns a list with at least `cell`, `cell_labels`, `cell_scores`.
   Optional: `alt_labels`, `cluster_summary`, `reference_source`,
   `warnings`. Pure: no shiny, no state mutation. Heavy deps must go
   through `require_optional()` (see `R/optional_deps.R`).
2. Register it in `ANNOTATION_ENGINES()`. The annotation module's UI
   rebuilds automatically; no edits there.
3. If the engine wraps a third-party tool that ships in its own shape
   (SingleR DataFrame, CellTypist json, etc.), factor the
   shape→engine-output mapping into a pure helper so it can be
   regression-tested without the heavy dependency installed. The
   SingleR engine (`.singler_to_engine_output`) is the reference
   pattern.

Still on the Phase-2 list: scType, UCell, GPT-assisted, consensus
voting, cross-species mapping (with an orthology ingest step in the
marker registry).

### Differential Expression → Explorer

1. Open **Differential Expression**, pick a grouping field, group 1, and group 2.
2. Click **Run Differential Expression**. The module computes Wilcoxon
   p-values + BH-adjusted q-values and stores everything in
   `state$analysis_results$de` (see schema below).
3. The volcano plot is clickable — clicking a point sets
   `state$selected_gene`, the inspect panel updates, and a *Send to
   Explorer* button writes the same value so the Basic scRNA Explorer's
   FeaturePlot re-colors.
4. The table is filterable (gene search, |log2FC|, max adj p) and sortable
   on every numeric column; sorting/filtering does not re-run the test.

### `state$analysis_results`

A named list, one entry per analysis kind. Each entry has the shape:

```r
list(
  status        = "not_run" | "running" | "completed" | "failed",
  results       = <data.frame or NULL>,
  params        = <list of run parameters>,
  error_message = <character(1) or NULL>,
  timestamp     = <POSIXct or NULL>,
  duration_ms   = <integer or NULL>
)
```

Currently populated slots:

| slot                                | producer                       | `results` shape                                                                 |
| ----------------------------------- | ------------------------------ | ------------------------------------------------------------------------------- |
| `state$analysis_results$de`         | `mod_differential_expression`  | data.frame: `gene, group_1, group_2, avg_log2FC, pct.1, pct.2, p_val, p_val_adj` |
| `state$analysis_results$pathway`    | `mod_pathway_analysis`         | data.frame: `pathway, collection, direction, n_genes_in_pathway, n_overlap, overlap_genes, odds_ratio, p_val, p_val_adj` |
| `state$analysis_results$imputation` | `mod_imputation`               | list: `method, genes, expression (gene -> numeric), reduction_used, k`           |
| `state$analysis_results$trajectory` | `mod_trajectory`               | list: `pseudotime (numeric, length=n_cells), cell, source, reduction_used, root_field, root_group, metadata_field` |

Future modules should reuse this convention: add a new named slot,
populate `status` first, then `results` once the run finishes.

### DE vs. Marker Investigation

| | Marker Investigation | Differential Expression |
| --- | --- | --- |
| Comparison | one group vs *the rest* | one group vs *another group* |
| Output schema | `group, gene, avg_log2FC, pct_in, pct_out, p_value`        | `gene, group_1, group_2, avg_log2FC, pct.1, pct.2, p_val, p_val_adj` |
| Multiplicity correction | none yet | Benjamini-Hochberg on the tested gene set |
| Stored in | (recomputed on demand)                                   | `state$analysis_results$de` (persists until next run / dataset switch) |
| Run trigger | reactive (auto)                                          | explicit "Run" button                                        |

### DE backends

`compute_de(..., backend = ...)` dispatches to a pluggable backend. The
result frame is identical across every backend (same columns, same
sort, same BH-adjusted p-values) so the rest of the app — including
annotation provenance stamping — does not care which implementation
ran.

| backend              | kind       | description                                                                       | dependency               |
| -------------------- | ---------- | --------------------------------------------------------------------------------- | ------------------------ |
| `wilcox_r`           | cell       | Pure-R per-gene Wilcoxon (or t-test); always available; slow on large data        | none                     |
| `presto`             | cell       | `presto::wilcoxauc()` — sparse-aware C++ kernel, 10–100× faster                  | `remotes::install_github("immunogenomics/presto")` |
| `auto`               | cell       | Picks `presto` if installed and `test == "wilcox"`, else `wilcox_r`               | n/a                      |
| `pseudobulk_naive`   | pseudobulk | Aggregate counts → log2(CPM+1) → per-gene t-test across pseudobulk samples         | none                     |
| `pseudobulk_edger`   | pseudobulk | Aggregate counts → `edgeR::glmQLFit` / `glmQLFTest`                                | `edgeR` (Bioconductor)   |
| `pseudobulk_deseq2`  | pseudobulk | Aggregate counts → `DESeq2::DESeq` / `results`                                     | `DESeq2` (Bioconductor)  |

The DE module exposes a backend dropdown labelled "DE backend" beneath
the main controls; default is `auto`. Unavailable backends are still
selectable so picking them yields a clear `install presto` / `install edgeR`
error rather than silent fallback. `de_available_backends()` is the helper
future UI work should consult; each entry includes a `kind` field
(`"cell"` vs `"pseudobulk"`).

### Pseudobulk DE (scaffolding)

Pseudobulk DE collapses cells from the same biological replicate into
aggregated pseudobulk samples *before* testing, which respects
sample-level variance and is the recommended approach for production
scRNA-seq DE (Squair 2021, Crowell 2020, Soneson 2018).

**Aggregation** lives in [`R/pseudobulk.R`](R/pseudobulk.R):

```r
pb <- aggregate_pseudobulk(dataset,
  grouping_field       = "condition",
  group_1              = "treat", group_2 = "ctrl",
  sample_by            = "sample",        # biological replicate id
  layer                = NULL,            # auto: prefers "counts"
  agg                  = "sum",
  min_cells_per_sample = 10L)
# Returns: list(matrix, sample_metadata, layer_used, agg, warn_lognorm, provenance)
```

**Required metadata.**

| field             | meaning                                              | required by                 |
| ----------------- | ---------------------------------------------------- | --------------------------- |
| `grouping_field`  | splits cells into the two groups being compared      | every backend               |
| `sample_by`       | biological replicate id (donor / sample / patient)   | every pseudobulk backend    |

`sample_by` and `grouping_field` must be *different* columns (crossed
design). The mock dataset ships `sample ∈ {S1, S2, S3}` × `condition ∈
{ctrl, treat}` × four `cluster`s, so pseudobulk DE on
`condition ~ sample` works out of the box.

**Raw counts vs log-normalised values.**

| layer    | content                                                              | used by                                                |
| -------- | -------------------------------------------------------------------- | ------------------------------------------------------ |
| `data`   | log-normalised expression (default)                                  | cell-level DE (`wilcox_r`, `presto`), markers, pathway |
| `counts` | raw UMI counts (Poisson draws in the mock; loader-provided otherwise)| pseudobulk DE (`pseudobulk_*`)                         |

`aggregate_pseudobulk()` defaults to `layer = "counts"` when present
and falls back to the backend's default layer (usually `"data"`) with
`warn_lognorm = TRUE`. Pseudobulk on log-normalised values is *not*
recommended for `edgeR` / `DESeq2` — they expect counts — but
`pseudobulk_naive` will run regardless.

**Per-gene `pct.1` / `pct.2`** (Seurat-style fraction-of-cells-expressing
metrics) are still computed at the cell level so the canonical DE
schema is preserved across every backend.

**Validation rules** (`validate_pseudobulk_inputs()`):

- `sample_by` is required and must differ from `grouping_field`.
- Each group must have at least `min_samples_per_group` pseudobulk
  samples *after* dropping samples with fewer than
  `min_cells_per_sample` cells (defaults: 2 and 10).
- Violations raise a clear error that includes the surviving sample
  counts so the user can adjust `min_cells_per_sample` or pick a
  different `sample_by`.

**Adding a new pseudobulk backend.**

1. Implement `.de_run_pseudobulk_<name>(pb, pcts, group_1, group_2, min_pct)`
   inside `R/de.R`. The aggregation step has already happened; you
   receive the pseudobulk count matrix + cell-level pct stats.
2. Convert your engine's output to the canonical DE schema via
   `.pseudobulk_to_de_schema()`. This handles BH adjustment, schema
   conformance, and the `min_pct` filter.
3. Add the backend id to `DE_BACKENDS` / `PSEUDOBULK_BACKENDS` and
   register it in `de_available_backends()`.
4. No changes to the DE module UI or to `compute_de()` are required.

### DE: future plug-in points

- **Per-cluster pseudobulk** (Crowell muscat-style): a wrapper that
  fans the dispatcher out across every cluster and concatenates
  results.
- **Covariate-aware models** (MAST `~ group + nFeature_RNA`): another
  backend slot; same schema.
- **DESeq2 shrunken LFC** (`lfcShrink`): one flag away inside
  `.de_run_pseudobulk_deseq2()`.

### Pathway Analysis → Explorer

1. After running **Differential Expression**, open **Pathway Analysis**.
   The page shows `No DE results yet.` until DE has completed.
2. Pick direction (up in group 1 / up in group 2 / both), adj-p and
   |log2FC| cutoffs, and a gene set collection (currently only
   `mock_v1`). A live count of "N genes selected" updates above the Run
   button.
3. Click **Run Pathway Analysis** → bar plot + table populate, results
   land in `state$analysis_results$pathway`.
4. Click a bar (or pick from *Inspect pathway*) → the overlapping genes
   appear. Pick one from *Send a gene to the Explorer* and click the
   button → `state$selected_gene` is set; the Basic scRNA Explorer's
   FeaturePlot recolors.

### Pathway internals & gene-set library

- **Sources** are pluggable. The Pathway module never touches gene-set
  internals — it goes through `available_pathway_collections()` and
  `get_pathways(collection)`. Both delegate to a registry in
  `R/pathway_sources.R` (`PATHWAY_SOURCES()`), exactly mirroring
  `ANNOTATION_ENGINES()` and `DE_BACKENDS()`.
- **Built-in source.** The `builtin` source ships one collection
  (`mock_v1`) with 6 pathways (T cell activation, B cell receptor
  signaling, Myeloid inflammatory response, Epithelial program,
  Extracellular matrix organization, Cytotoxicity). Each pathway has
  13–16 genes chosen so the mock dataset's six demo genes each fall
  into at least one pathway (so the UI shows real overlap). Listed
  first, so `available_pathway_collections()[1]` is always `mock_v1`.
- **MSigDB source.** When the `msigdbr` package is installed, five
  pinned collections become available (see *Pathway sources* below).
  A process-local cache in `R/pathway_sources.R::.msigdbr_cache`
  prevents repeated MSigDB lookups within a session.
- **Statistical kernel (ORA):** per-pathway Fisher's exact test
  (one-sided "greater"), Benjamini–Hochberg adjustment across the
  tested pathways in a run. When `direction = "both"`, the module runs
  both directions independently and re-adjusts p-values jointly.
- **Universe:** the population from which `selected` was drawn.
  - `universe = de_tested_genes` — the set of genes actually tested by
    DE. Most defensible for small / mock data.
  - `universe = available_genes(dataset)` — the full filtered feature
    set. Sensible for production-scale DE.
  - `universe = NULL` — fallback `union(selected, unlist(pathways))`,
    which is liberal (small denominators inflate odds ratios).

### Pathway sources

| Source     | Collection ids                                                          | Requires    | Notes |
|------------|-------------------------------------------------------------------------|-------------|-------|
| `builtin`  | `mock_v1`                                                                | (none)      | Default, ships with the app |
| `msigdbr`  | `msigdbr/H`, `msigdbr/C2:CP:REACTOME`, `msigdbr/C2:CP:KEGG_LEGACY`, `msigdbr/C5:GO:BP`, `msigdbr/C7` | `msigdbr` (CRAN); recent `msigdbr` versions also need the `msigdbdf` data package | Pinned subcategories; defaults to `species = "Homo sapiens"` |

Resolution is single-namespace: `get_pathways(id)` walks the registry
and the first source that owns `id` wins. Built-in is checked first so
existing default-selection logic (`available_pathway_collections()[1]`)
does not regress when new sources are installed.

### Adding a new pathway source

1. Implement a fetcher `function(collection_id) -> named list` and a
   `collections_fn() -> character()` that advertises ids (already
   prefixed). Keep the schema-mapping logic in a *pure* helper so it's
   testable without the heavy data package — `.msigdbr_to_pathways()`
   is the reference pattern.
2. Register the source by adding a `pathway_source_spec(...)` entry to
   `PATHWAY_SOURCES()` in `R/pathway_sources.R`. Use `requires` to
   declare optional package dependencies — the registry skips your
   source automatically when the package isn't installed.
3. Add tests against the pure converter and a missing-dep test that
   asserts a clear "install X" error.
4. No changes to `compute_enrichment()` or the Pathway module are
   needed. The function expects `pathways` as a named list of
   character vectors.

Stubbed future sources (each will land as a new spec):
- **Reactome** via `reactome.db` or the Reactome web service
- **GO** via `org.Hs.eg.db` + `GO.db`
- **KEGG** via `KEGGREST`
- **User .gmt** — local file uploaders

### Ranked GSEA scaffold (`compute_gsea` + fgsea)

`compute_gsea(ranked_genes, pathways, ...)` is the design seam for
ranked-list enrichment. It mirrors the shape of `compute_enrichment()`
so callers can swap ORA / GSEA without restructuring downstream code.

- Backend: `fgsea::fgsea()` (Bioconductor), gated by
  `require_optional("fgsea", source = "Bioconductor")`. Until `fgsea`
  is installed, calls raise a clear install error.
- Pure converter `.fgsea_to_gsea_schema()` is regression-tested without
  fgsea, mirroring the SingleR / msigdbr pattern.
- Output schema:
  `pathway, collection, n_genes_in_pathway, n_leading_edge, leading_edge_genes, ES, NES, p_val, p_val_adj`.

The Pathway module UI is still ORA-only; wiring `compute_gsea()` into
the module is a future PR.

### Pathway limitations

- ORA-by-default. GSEA is a scaffold (`compute_gsea`) but not yet wired
  into the module UI.
- Single hypergeometric test per pathway; ignores gene-length /
  expression bias (a `goseq`-style correction would mitigate this).
- No collapsed/parent-aware multiple testing (e.g. Reactome hierarchy);
  flat BH only.

### Data Smoothing / Imputation (visualization-only)

The Imputation module produces smoothed expression vectors for
**visualization-only** exploration. **Raw expression remains the source
of truth** for DE, Marker Investigation, and Pathway Analysis — those
modules read `get_gene_expression()` directly, never the smoothed slot.

#### Where it lives

- Maths and helpers: `R/imputation.R`. Three mock methods:
  - `none` — pass-through
  - `neighbor` — one-pass kNN average in 2D embedding space (default k=15)
  - `alra_mock` — kNN average + soft-threshold (lower 30% → 0)
  - `magic_mock` — two-pass kNN diffusion (larger k)
- Module: `R/modules/mod_imputation.R`.

#### How smoothed values are stored separately from raw

```r
state$analysis_results$imputation = list(
  status        = "completed",
  results       = list(
    method         = "neighbor",
    genes          = c("CD3D", "MS4A1", ...),
    expression     = list(CD3D = numeric(n_cells), MS4A1 = ...),  # SMOOTHED
    reduction_used = "UMAP",
    k              = 15L
  ),
  params        = list(assay, method, k, genes),
  error_message = NULL,
  timestamp     = <POSIXct>,
  duration_ms   = <integer>
)
```

`dataset$expression` (the raw vectors) is never mutated. The smoothed
copy lives entirely under `state$analysis_results$imputation`.

#### How the display mode works

A single shared flag controls visualization-only switching:

```r
state$display_mode_imputation  # "raw" (default) | "smoothed"
```

- The Basic scRNA Explorer renders a radio toggle that appears **only**
  when `has_smoothed_results(state)` is `TRUE`. Switching writes to
  `state$display_mode_imputation`.
- The Explorer's FeaturePlot reads expression via
  `get_gene_expression_for_view(state, gene)`, which returns the
  smoothed vector iff the toggle is `"smoothed"` AND that gene was in
  the smoothing run; otherwise it falls back to raw.
- **No other module calls `get_gene_expression_for_view()`.** Marker
  Investigation, DE, and Pathway Analysis call `get_gene_expression()`
  directly. That's the entire mechanism: a one-helper seam keeps
  smoothed values from ever leaking into a statistical computation.

#### Why DE/markers stay raw by default

- Smoothing introduces strong correlations between neighbouring cells,
  which inflates the apparent significance of any group-wise test built
  on Wilcoxon, t-test, or Fisher's exact.
- It also collapses true biological variability that DE relies on.
- The accepted community guidance (Hou et al. 2020;
  Andrews & Hemberg 2018) is: smoothing is a visualization aid, not an
  analytic substrate. If you want to run DE on imputed values, do so
  consciously by switching the module to a "use smoothed" mode — not by
  having the app silently swap data sources.

Clearing smoothed data resets `state$display_mode_imputation` to
`"raw"`, so no module is left rendering against a deleted slot.

#### Imputation: future plug-in points

- **ALRA** — low-rank SVD recovery via `ALRA::alra()`
- **MAGIC** — diffusion via `Rmagic` / `phateR`
- **SAVER / scImpute** — Bayesian / regression-based imputation
- **Real kNN** — compute neighbours on the PCA embedding via
  `RANN::nn2()` or `BiocNeighbors`, cache the graph for reuse across
  methods. Today's mock uses 2D-UMAP-space kNN, which is the wrong
  geometry for real data.

All four would slot into `compute_smoothed()` as additional methods.
The `IMPUTATION_METHODS` registry is the only place the module touches.

### Trajectory / Pseudotime

Pseudotime lives only in `state$analysis_results$trajectory` —
DE / Markers / Pathway / Imputation never consume it automatically.

#### Where it lives

- Orchestrator + binning helpers: `R/trajectory.R`
  (`run_trajectory()` → full result payload; `compute_pseudotime()` →
  bare numeric vector).
- Method registry: `R/trajectory_registry.R` (`TRAJECTORY_METHODS()`,
  `trajectory_method_spec()`, `available_trajectory_methods()`).
- Built-in methods: `mock` + `metadata` are defined inside
  `R/trajectory.R` as `.run_mock_trajectory()` / `.run_metadata_trajectory()`.
- Optional backends: `R/trajectory_slingshot.R`, `R/trajectory_monocle3.R`.
- Plotting helper: `plot_gene_vs_pseudotime()` in `R/plotting.R`.
- Module: `R/modules/mod_trajectory.R`. UI is driven entirely by the
  registry — adding a new backend surfaces in the picker automatically.

#### Trajectory methods

| method     | kind         | requires           | needs root | description |
| ---------- | ------------ | ------------------ | ---------- | ----------- |
| `mock`     | demo         | (none)             | yes        | Euclidean distance from root cluster centroid in a 2D reduction, normalised to [0, 1]. Always available. |
| `metadata` | precomputed  | (none)             | no         | Rescale an existing numeric metadata column to [0, 1] (e.g. `pseudotime_demo`). |
| `slingshot`| real lineage | `slingshot`        | yes (start cluster) | Minimum-spanning-tree over cluster centroids + per-cell principal-curve fits. Operates on any embedding. |
| `monocle3` | real lineage | `monocle3`         | yes (root group)    | Principal-graph learning + `order_cells()`. Requires raw counts (defaults to the `counts` layer). |

`PSEUDOTIME_SOURCES()` and `available_pseudotime_sources()` are
back-compat aliases for `trajectory_method_choices()`. Unavailable
methods (missing optional packages) stay in the picker but are
labelled `(not installed)` so the user gets a clear `require_optional`
error instead of a silent fallback.

#### How pseudotime is stored

```r
state$analysis_results$trajectory = list(
  status        = "completed",
  results       = list(
    pseudotime     = numeric(n_cells),     # [0, 1]
    cell           = character(n_cells),   # cell barcodes (alignment key)
    source         = "mock" | "metadata" | "slingshot" | "monocle3",
    reduction_used = "UMAP" | NA,
    root_field     = "cluster" | NA,
    root_group     = "0" | NA,
    metadata_field = "pseudotime_demo" | NA,
    n_lineages     = 1L,                   # >1 for slingshot multi-lineage
    method_details = list()                # backend-specific (lineage psts, etc.)
  ),
  params        = list(source, reduction, root_field, root_group, metadata_field),
  error_message = NULL, timestamp = <POSIXct>, duration_ms = <integer>,
  annotation_stamp = list(...)
)
```

#### How pseudotime differs from metadata

Pseudotime is an **analysis result**, not a property of the dataset.

- Stored under `state$analysis_results$trajectory` (a Run output), not
  in `dataset$cell_data` / `dataset$metadata_fields`. Re-running gives
  a new result; the dataset is untouched.
- The Trajectory module does not silently inject pseudotime into the
  Explorer's metadata picker. To make pseudotime available as a
  metadata axis, click **Apply pseudotime to dataset** in the
  Trajectory module; this calls `apply_pseudotime_to_dataset()` to
  write a dated column:
    - `pseudotime__<source>__<YYYY_MM_DD>` (numeric, always written)
    - `pseudotime_bin__<source>__<YYYY_MM_DD>` (factor, only when the
      "+ bin column" checkbox is ticked; handy for grouping in DE or
      coloring by a categorical lineage stage)
  Each new column carries provenance attrs (`pseudotime_source`,
  `reduction_used`, `root_field`, `root_group`, `metadata_field`,
  `n_lineages`, `applied_at`, `kind`) so downstream tools can trace
  where the values came from. The pre-shipped `pseudotime_demo`
  column on the mock dataset remains a no-Run quick-look path.
- `state$selected_gene` syncs between Trajectory, Explorer, DE, Markers,
  and Pathway. The Trajectory module reads it for the gene-vs-time
  trend, and clicking *Send to Explorer* writes it back.

#### Adding a new trajectory backend

1. Implement `.run_<name>_trajectory(dataset, params)` in a new file
   `R/trajectory_<name>.R`. The function should:
   - Validate `params$reduction`, `params$cluster_field`,
     `params$root_field`, `params$root_group`, etc. as needed.
   - Call `require_optional(<your-pkg>)` for any heavy dependency.
   - Pull expression via `backend_as_matrix()` /
     `get_gene_expression()`; pull embeddings via `get_embedding()`;
     pull cluster labels via `get_metadata()`.
   - Return the canonical payload (see
     `trajectory_method_spec()` docstring). At minimum
     `list(pseudotime = numeric(n_cells), source = "<name>")` — the
     orchestrator fills in the rest.
2. Factor the schema mapping into a pure helper
   (`.<name>_to_pseudotime()`) that accepts a stand-in input shaped
   like the real engine's output. Lets you regression-test schema
   conformance without the heavy dependency.
3. Register the method by adding a `trajectory_method_spec(...)`
   entry to `TRAJECTORY_METHODS()` in `R/trajectory_registry.R`.
4. Add tests against the pure converter and a missing-dep test that
   asserts a clear "install X" error.
5. No changes to `compute_pseudotime()`, `run_trajectory()`, or the
   Trajectory module are required.

Stubbed future routes (each will land as a new spec):

- **Palantir** — via `reticulate` to the Python package
- **Diffusion pseudotime** — `destiny::DPT()`
- **scVelo / RNA velocity** — would also require a velocity matrix in
  the dataset schema; out of scope for the current shell.

#### Mock pseudotime: known limitations

- 2D-UMAP geometry is the wrong space for real biology — UMAP distorts
  distances. The mock chooses it because it's the only embedding we
  have. Real implementations (Slingshot, Monocle3) should be preferred.
- Single root, no branches. Real trajectory inference identifies
  multiple endpoints, branch points, and a partition lattice.
- No temporal direction. The "root → high-pt" direction is asserted by
  the user's root choice, not inferred from data (velocity, RNA-spliced
  ratios, etc.).
- Deterministic Euclidean distance ignores cell-cell similarities in
  gene-expression space.

### Regulons / Network Analysis

The Regulons module scores TF regulon activity (TF + its target genes)
per cell using AUCell, and renders two views: a regulon \u00d7 group
heatmap and an embedding colored by a selected regulon's AUC. The
canonical result lives in `state$analysis_results$regulons$results`;
it is never silently injected into the Explorer / DE / Markers /
Pathway / Trajectory metadata pickers.

The module is **registry-driven** on both axes:

- **Sources** (`REGULON_SOURCES()` in `R/regulon_sources.R`) supply the
  regulon catalogue. The module's "Regulon source" dropdown is
  populated from `list_regulon_sources()`; sources whose optional
  dependency is missing are kept and tagged `(not installed)`.
- **Engines** (`REGULON_ENGINES()` in `R/regulon_registry.R`) supply
  the scoring backend. The module's "Scoring engine" dropdown is
  populated from `list_regulon_engines()` and behaves the same way.

#### Regulon sources

| source             | species  | what it returns                                                                                | dependencies                                       |
| ------------------ | -------- | ---------------------------------------------------------------------------------------------- | -------------------------------------------------- |
| `mock_pbmc`        | human    | Four hand-built regulons (GATA3, PAX5, SPI1, KLF5) whose target genes are the canonical `MOCK_GENES` used by `mock_dataset()`. Always available; useful for CI, demos, and validating the engine. | none                                               |
| `dorothea_human`   | human    | DoRothEA curated TF-target regulons (Garcia-Alonso et al, 2019) at a configurable confidence prefix (A, AB, ABC, ABCD, ABCDE). | `BiocManager::install("dorothea")`                 |
| `dorothea_mouse`   | mouse    | DoRothEA mouse equivalent.                                                                       | `BiocManager::install("dorothea")`                 |

Every source returns a `regulon_set()` (see `R/regulon_schema.R`); the
TFs become regulon ids, the targets feed AUCell scoring.

#### Regulon engines

| engine     | what it does                                                          | dependencies                                |
| ---------- | --------------------------------------------------------------------- | ------------------------------------------- |
| `aucell_r` | Pure-R AUCell. Computes per-cell normalised AUC of each regulon's target ranks within the top-N expressed genes. No deps. Numerically consistent with the Bioconductor implementation on small dense matrices; slower on > 50k cells. | none                                        |
| `aucell`   | Wraps `AUCell::AUCell_calcAUC()`. Sparse-aware, fast on large datasets. | `BiocManager::install("AUCell")`            |

Both produce the same output schema (`regulon_result_v1`): an
`n_cells x n_regulons` AUC matrix plus engine / source / version
provenance. The module reads the matrix to render the heatmap
(`regulon_mean_by_group()`) and the per-regulon embedding plot.

`top_n_fraction` (default 0.05) controls the rank threshold — fraction
of the gene catalogue that counts as the cell's "top". Higher values
saturate AUC for sharper regulons and reduce sensitivity to small
ones.

#### Adding a new regulon source

1. Implement `fetch_fn(params)` returning a `regulon_set()`. Pure: no
   shiny, no state mutation. Heavy deps via `require_optional()`.
2. Register it in `REGULON_SOURCES()`. The Regulons module's source
   picker rebuilds automatically.
3. If the source wraps a third-party catalogue with its own shape
   (DoRothEA data.frame, SCENIC `regulons.gmt`, etc.), factor the
   shape \u2192 `regulon_set` mapping into a pure helper (the
   `.dorothea_df_to_regulons()` pattern) so the schema mapping is
   testable without the heavy dependency installed.

#### Adding a new regulon engine

1. Implement `.run_<name>_regulons(dataset, regulon_set, params)`
   returning a list `list(cell, regulon_ids, auc_matrix, warnings)`.
   Pure: no shiny, no state mutation.
2. Register it in `REGULON_ENGINES()`. Add a parameter spec for the
   UI to render controls.
3. Factor the schema mapping into a pure helper
   (`.<name>_to_regulon_engine_output()`) so it can be regression-
   tested without the heavy dependency installed (the
   `.aucell_to_regulon_engine_output()` pattern).

#### Regulons: known limitations

- AUCell scores regulon **enrichment**, not transcription factor
  activity in the mechanistic sense. A high AUC just means many of
  the TF's targets are highly expressed in the cell.
- The pure-R engine builds a dense rank matrix; > 50k cells will hurt.
  Use the Bioconductor engine on real-sized data (it operates on
  sparse rankings).
- DoRothEA targets are species-restricted. Use the matching
  `dorothea_human` / `dorothea_mouse` source for the species of your
  dataset.

Still on the regulons roadmap: full SCENIC pipeline (GENIE3 +
cisTarget + AUCell, plumbed through pySCENIC + reticulate),
`apply_regulon_scores_to_dataset()` (mirroring
`apply_pseudotime_to_dataset()`), and a regulon-specificity score
(RSS) view.

### How the real Seurat / AnnData loaders map onto the schema

Mapping reference for the loaders in `R/dataset.R`:

| dataset field      | Seurat source                                                   | AnnData source                              |
| ------------------ | --------------------------------------------------------------- | ------------------------------------------- |
| `assays`           | `Seurat::Assays(obj)`                                           | n/a (typically just `"RNA"`)                |
| `reductions`       | `Seurat::Reductions(obj)`                                       | `names(adata$obsm)` stripped of `X_` prefix |
| `metadata_fields`  | `colnames(obj@meta.data)`                                       | `colnames(adata$obs)`                       |
| `cells`            | `colnames(obj)`                                                 | `adata$obs_names`                           |
| `cell_data`        | `cbind(obj@meta.data, Embeddings(obj, red))` for every `red`    | `cbind(adata$obs, adata$obsm[[k]])`         |
| `genes`            | `rownames(GetAssayData(obj, assay, "data"))`                    | `adata$var_names`                           |
| `expression[[g]]`  | `as.numeric(GetAssayData(obj, assay, "data")[g, ])`             | `adata$X[, g]` or `adata$layers[[...]]`     |

The expression layer is always wrapped in an `expression_backend` so
loaders can be lazy or eager without rippling through callers. The
in-memory + sparse backends materialise everything on load; the
`expression_backend_h5ad()` reads gene-by-gene on demand via `rhdf5`.
Modules only depend on `available_genes()` / `get_gene_expression()`
and stay loader-agnostic.
