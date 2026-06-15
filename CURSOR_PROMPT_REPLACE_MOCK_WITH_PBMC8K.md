# Cursor prompt: replace the current mock dataset with PBMC 8k

You are working on an R Shiny single-cell explorer app that currently uses a built-in mock/demo dataset (approximately 2.5k PBMC cells) for the sidebar "Load mock dataset" workflow.

Your task is to replace that current mock/demo dataset with a PBMC 8k dataset, while preserving the existing app behavior and developer ergonomics as much as possible.

## Objective
Make the app load a PBMC 8k-based demo dataset instead of the current smaller mock dataset when the user clicks the sidebar dataset-load action.

The result should:
- preserve the current app flow
- preserve existing module contracts as much as possible
- remain easy for contributors to run locally
- keep the demo dataset suitable for showcasing the app

## Important constraints
- Do not redesign the UI.
- Do not change analytical module logic unless necessary for compatibility.
- Preserve the existing dataset API/shape expected by modules as much as possible.
- Keep the existing "Load mock dataset" action text unless you have a very strong reason to rename it.
- Prefer a robust, reproducible implementation over a fragile one.

## What to investigate first
Please inspect how the current demo dataset is defined and loaded.

Likely relevant files include:
- `app.R`
- `R/dataset.R`
- `R/dataset_helpers.R`
- any file defining `mock_dataset()` or equivalent
- any registry/helper code that assumes the current demo dataset size/content

Determine:
1. where the current mock dataset comes from
2. the exact object/schema expected by the rest of the app
3. whether the demo dataset is generated on the fly, loaded from disk, or bundled in code

## Recommended implementation goal
Use a real PBMC 8k-style demo dataset in a way that is reproducible for local development.

Prefer one of these strategies, in order of preference:

### Preferred option A: local bundled/prepared PBMC 8k artifact
If the app already has a data directory or an appropriate place to store a prepared serialized demo object, add a prepared PBMC 8k demo artifact and load it locally.

This is preferred if it keeps app startup and demo loading reliable.

### Option B: programmatic conversion from a standard PBMC 8k source
If bundling a prepared object is not appropriate, implement a reproducible loader that constructs the app’s expected dataset object from a PBMC 8k source.

This may involve:
- loading from SeuratData if appropriate
- loading from a local 10x matrix if that is the project convention
- converting the object into the app’s internal dataset structure

### Avoid if possible
- hidden network downloads at runtime when a user clicks the button
- brittle assumptions about packages that may not be installed
- forcing a heavyweight setup without clear fallback behavior

## Required behavior
When the user clicks the existing dataset-load button, the app should load the PBMC 8k demo dataset into the same shared state mechanism currently used by the mock dataset.

This should preserve downstream expectations for modules such as:
- overview
- basic exploration
- marker investigation
- differential expression
- pathway analysis
- annotation
- trajectory
- regulons
- imputation

## Compatibility requirements
The replacement dataset should expose the same kinds of fields the app expects, including as needed:
- dataset name
- source/provenance
- cell count
- gene count
- assays
- reductions
- metadata fields
- any expression accessors required by plotting / DE / markers / pathway / annotation / trajectory modules

If the PBMC 8k source object does not naturally contain some fields the app expects, add a clean compatibility/preparation layer rather than hacking many downstream modules.

## Strong preference
Please implement the change by adapting the dataset-preparation layer rather than editing many modules.

In other words:
- prepare/normalize the PBMC 8k dataset into the existing app dataset schema
- avoid touching module logic unless a small compatibility fix is truly necessary

## Things to consider carefully
- Does the app expect mock fields like cluster labels, reductions, pseudotime demo columns, annotation-ready metadata, regulon-ready metadata, etc.?
- If so, what is the minimal and cleanest way to derive or add those fields to the PBMC 8k demo dataset?
- Some modules may expect precomputed reductions such as PCA/UMAP.
- Some modules may work better if cluster labels are already present.
- Preserve a good out-of-the-box experience for demoing the app.

## Implementation guidance
A good solution may involve:
1. locating the current mock dataset constructor
2. introducing a `pbmc8k_demo_dataset()` (or similarly named) constructor
3. keeping `mock_dataset()` as a compatibility wrapper that now returns the PBMC 8k-backed prepared dataset
4. adding any preparation steps needed to satisfy the app’s expected schema
5. updating provenance/name text so users can tell what dataset was loaded

If you think renaming the internal function is helpful, keep a compatibility alias so existing code paths do not break.

## If external packages/data are needed
If your implementation depends on a package such as SeuratData or a locally stored PBMC 8k artifact, handle this carefully:
- fail clearly with a helpful message if the dependency is missing
- document the expected setup
- avoid surprising runtime behavior

If there is already a project convention for optional dependencies, follow that convention.

## Deliverables
Please make the code changes directly.

Also provide a concise summary covering:
1. where the old mock dataset came from
2. how the PBMC 8k dataset is now loaded/prepared
3. what files changed
4. whether any compatibility fields were added to make the PBMC 8k dataset fit the app schema
5. any local setup step required for future developers

## Non-goals
- no UI redesign
- no major refactor of analytical modules
- no unnecessary changes outside the dataset-loading layer

## Quality bar
The app should continue to feel like it has a dependable built-in demo dataset, but that demo dataset should now be PBMC 8k-based rather than the current ~2.5k mock PBMC dataset.
