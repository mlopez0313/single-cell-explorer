# Cursor prompt: fix PBMC 8k artifact path handling and improve app dependency bootstrapping

You are working on an R Shiny single-cell explorer app that was recently updated to prefer a prepared PBMC 8k demo artifact. There are now two practical issues to fix:

1. The app is looking for the prepared demo artifact at an incorrect machine-specific path.
2. The app should do a better job of ensuring required dependencies are available before launch, or at least fail/help clearly in a way that makes full setup easy.

## Current observed problem
At runtime, the app reports a message like:

> Prepared PBMC 8k demo artifact not found at /Users/.../inst/extdata/pbmc8k_demo.rds and auto-build was skipped (missing package(s): TENxPBMCData). Falling back to mock_dataset(). Run `Rscript scripts/build_pbmc8k_demo.R` (or install the missing packages) to enable the PBMC 8k demo.

This suggests:
- the artifact path resolution is machine-specific or otherwise incorrect
- dependency/bootstrap behavior is not yet strong enough for a smooth first run

## Objective
Fix both of the following:

### A. PBMC 8k artifact path resolution
The app must resolve the prepared demo artifact path in a robust, machine-independent way.

### B. Dependency bootstrapping before full app launch
Make the app/repo do the best practical job of ensuring dependencies needed for a full-featured local run are installed before launching the app, or at least provide a single clear setup/bootstrapping path.

## Constraints
- Do not redesign the UI.
- Do not change analysis logic unnecessarily.
- Prefer robust repo-relative/package-relative path handling over hardcoded absolute paths.
- Do not make the normal app launch path fragile.
- Be careful about auto-install behavior: it should be explicit, safe, and developer-friendly.

## Part A: fix artifact path handling
Please inspect where the PBMC 8k artifact path is computed.

Likely relevant areas include:
- `R/dataset.R`
- `R/dataset_helpers.R`
- the PBMC 8k loader function
- any build/regeneration script
- any use of `here::here()`, `system.file()`, `normalizePath()`, `getwd()`, or user-specific hardcoded paths

### Requirements
The artifact path should resolve correctly when running the app from the project checkout on another machine.

Prefer a robust strategy such as:
- package-style path resolution with `system.file("extdata", "pbmc8k_demo.rds", package = ...)` if the project is structured/installed that way, or
- a repo-relative path strategy from the project root, or
- a small helper that can resolve dev-mode vs installed-package mode cleanly

### Important
Avoid hardcoded developer-specific absolute paths.

Please make the runtime loader and any build script agree on the artifact location.

Also verify whether the artifact is supposed to live at something like:
- `inst/extdata/pbmc8k_demo.rds`
- or another project-local data path

Use one clear convention consistently.

## Part B: improve dependency bootstrapping
I want the project to do the best practical job of getting all needed dependencies installed **before launching the app in full**, if possible.

That does **not** mean the app should silently install a huge number of packages every time it starts.

Instead, implement a sensible developer bootstrap/setup flow.

## Preferred outcome
Provide a clear setup mechanism that can ensure required packages are present for:
- core app launch
- the PBMC 8k demo artifact build path
- important optional full-feature modules where practical

This could take one or more of these forms:

### Preferred approach
Add a developer-facing setup script/function such as:
- `scripts/setup_dev.R`
- or `R/setup.R` with a callable helper

This setup path should:
- install required CRAN/Bioconductor dependencies if missing
- clearly separate required vs optional dependencies if appropriate
- include the packages needed to build/load the PBMC 8k demo artifact
- be documented as the recommended first-run step before launching the app

### Optional additional improvement
If appropriate, add a lightweight preflight check during app launch that:
- verifies core packages are installed
- gives a clear actionable message if not
- points the user to the setup script/function

This preflight should be helpful, not noisy or overly magical.

## Strong preferences for dependency handling
- Avoid silent package installation during ordinary app launch unless it is extremely targeted and clearly justified.
- Prefer explicit setup/bootstrap commands over surprising runtime installs.
- If you implement any auto-install option, gate it behind an explicit user choice or documented setup command.
- Follow any existing optional-dependency conventions already present in the repo.

## What to inspect
Please inspect:
- `DESCRIPTION`
- app startup code in `app.R`
- any optional dependency helpers already present
- the PBMC 8k build script
- any package checking/install logic already in the project

Determine:
1. what dependencies are strictly required for app launch
2. what dependencies are needed for the PBMC 8k prepared artifact build
3. what dependencies are optional but important for a full demo experience

## Implementation goals
A good implementation would likely include:

### 1. Artifact path fix
- central helper for resolving the PBMC 8k artifact path
- shared by both runtime loader and build script
- works across machines and dev environments

### 2. Setup/bootstrap path
- a documented setup script/function for installing dependencies
- support for CRAN + Bioconductor packages as needed
- inclusion of `TENxPBMCData` and any packages required for PBMC 8k preparation

### 3. Launch-time preflight (lightweight)
- check for core launch dependencies
- if missing, fail with a helpful message pointing to setup
- optionally warn when optional/full-demo packages are missing

## Documentation / messaging
Please improve error/help messages so they are concise and actionable.

For example:
- if the artifact is missing, say where the app expects it and how to generate it
- if a package is missing, say exactly which setup command to run
- avoid machine-specific confusion in the message text

## Deliverables
Please make the code changes directly.

Also provide a concise summary covering:
1. what caused the bad artifact path
2. how artifact path resolution now works
3. what setup/bootstrap mechanism was added
4. what packages are considered core vs build-time vs optional
5. what a developer should run before launching the app for a full local setup

## Non-goals
- no UI redesign
- no major module rewrites
- no brittle always-on auto-installs at app startup

## Quality bar
A fresh developer on a different machine should be able to:
1. run one clear setup/bootstrap step
2. build or access the PBMC 8k prepared artifact at the correct project-local path
3. launch the app with a much smoother full-feature experience
