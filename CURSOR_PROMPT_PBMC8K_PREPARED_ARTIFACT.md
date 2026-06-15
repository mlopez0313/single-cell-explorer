# Cursor prompt: replace the current demo dataset with a prepared PBMC 8k artifact

You are working on an R Shiny single-cell explorer app that currently uses a built-in mock/demo dataset for the sidebar dataset-load action. The current demo dataset is a smaller mock PBMC-like dataset. The goal is to replace it with a **PBMC 8k-based demo dataset**, but to do so in a way that is robust, reproducible, and friendly for local development.

## Core implementation preference
Do **not** make the demo button depend on a runtime network download.

Do **not** use SeuratData as the primary runtime data source for the app.

Instead, implement the demo dataset as a **prepared serialized local artifact** (for example an `.rds` object or equivalent), already normalized into the app’s expected dataset schema.

If SeuratData is used at all, it should be used only as an **offline/build-time source** for generating that prepared artifact, not as the live runtime dependency behind the button.

## Objective
Make the app load a prepared PBMC 8k-based demo dataset when the user clicks the existing dataset-load action.

The result should:
- preserve the current app flow
- preserve the existing shared state behavior
- minimize downstream module changes
- keep the demo experience fast and dependable
- avoid runtime install/download surprises

## Constraints
- Do not redesign the UI.
- Do not make analytical modules fetch remote data.
- Do not require a network call when the user clicks the demo dataset button.
- Prefer adapting the dataset-preparation/loading layer over modifying many modules.
- Preserve the current dataset API/schema expected throughout the app.
- Keep the existing sidebar action text unless there is a compelling reason to change it.

## What to inspect first
Please inspect the current demo dataset implementation and dataset schema.

Likely relevant files include:
- `app.R`
- `R/dataset.R`
- `R/dataset_helpers.R`
- any file defining `mock_dataset()` or equivalent
- any utility or state code that assumes the current dataset shape

Determine:
1. where the current demo/mock dataset comes from
2. the exact schema/modules expectations for a loaded dataset object
3. what metadata/reductions/fields are assumed by downstream modules
4. whether there is already a project convention for storing bundled data artifacts or local demo assets

## Required design
Implement this in two layers:

### Layer 1: runtime demo loader
At runtime, when the user clicks the existing demo-load button, the app should load a **local prepared PBMC 8k artifact**.

This runtime path should:
- be fast
- avoid network access
- avoid surprise package installs
- produce the same kind of object that the rest of the app already expects

### Layer 2: optional build/regeneration path
Separately, if useful, provide a **developer-facing script or function** to regenerate the prepared artifact from an upstream source.

That upstream source may be:
- SeuratData
- a local 10x PBMC 8k matrix
- another clearly documented PBMC 8k source

But this build path should be separate from normal app runtime.

## Strong preference for implementation approach
Please implement the change by preparing a PBMC 8k-derived dataset that already matches the app’s internal schema.

That means:
- do schema normalization once in the dataset layer
- do not scatter PBMC-specific compatibility hacks across modules
- avoid modifying modules unless a very small compatibility fix is truly necessary

## Compatibility requirements
The prepared PBMC 8k demo dataset should expose the fields the app expects, including as needed:
- dataset name
- source/provenance
- cell count
- gene count
- assays
- reductions
- metadata fields
- expression accessors used by plots / DE / markers / pathway / annotation / trajectory / regulon workflows

Please pay special attention to demo-friendly fields that the app may implicitly rely on, such as:
- cluster labels
- categorical metadata for grouping
- numeric metadata where useful
- reductions like PCA / UMAP if expected
- any pseudotime-demo-like field if the app benefits from it

If the upstream PBMC 8k source does not naturally provide everything needed, add a clean preparation step in the dataset layer.

## Recommended implementation shape
A good solution will likely involve:
1. locating the current `mock_dataset()` implementation
2. introducing a prepared-artifact loader, e.g. something like:
   - `load_demo_dataset()`
   - `pbmc8k_demo_dataset()`
   - or similar
3. keeping `mock_dataset()` as a compatibility wrapper if that avoids downstream changes
4. storing the prepared artifact in a sensible local path within the project (or a clearly documented local data path)
5. optionally adding a build script/function to regenerate the artifact from SeuratData or another source

## Preferred runtime behavior
When the artifact exists locally:
- load it directly
- validate it enough to fail clearly if corrupted or incompatible

When the artifact is missing:
- fail clearly with a helpful developer-facing message
- optionally suggest the build/regeneration command
- do not silently try to download data from the internet during normal runtime unless there is already a very explicit project convention for that

## Optional build pipeline
If you add a build script/function, it should:
- clearly document its upstream source
- clearly document required packages
- transform the upstream PBMC 8k object into the app’s internal dataset schema
- save the prepared artifact for future runtime use

If SeuratData is used here, that is fine — but only in this build/regeneration context.

## File-level expectations
At minimum, inspect and update as needed:
- `R/dataset.R`
- `R/dataset_helpers.R`
- `app.R`
- any helper defining the current demo dataset

If you add a regeneration script, place it somewhere sensible and consistent with the repo layout.

## Developer ergonomics
Please make the solution maintainable for future contributors.

This means:
- clear naming
- clear provenance metadata for the prepared dataset
- minimal runtime surprises
- concise documentation/comments where needed

## Deliverables
Please make the code changes directly.

Also provide a concise summary covering:
1. where the old demo dataset came from
2. where the new prepared PBMC 8k artifact lives
3. how the runtime loader now works
4. whether you added a regeneration/build path and where it lives
5. what compatibility fields/preparation steps were needed to fit the app schema
6. any setup step required for future developers

## Non-goals
- no UI redesign
- no major analytical refactor
- no unnecessary module rewrites
- no runtime internet dependency for the demo button

## Quality bar
The app should continue to feel like it has a dependable built-in demo dataset, but that demo dataset should now be a PBMC 8k-based prepared local artifact rather than the current smaller mock dataset.
