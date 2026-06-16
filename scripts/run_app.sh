#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Launcher for the single-cell-explorer Shiny app.
#
# The app's install / build paths (CRAN, Bioconductor, GitHub) use R's
# `tempdir()` for both downloads and build scratch. On deployments where
# the OS `/tmp` (or `/`) is small but the app itself lives on a much
# larger filesystem (e.g. an NFS mount holding both the source tree
# AND the conda env), forcing multi-GB packages -- BSgenome.Hsapiens.UCSC.hg38
# is the worst offender at ~870 MB tarball / ~3.4 GB unpacked -- through
# the small partition causes silent install failures (out-of-space mid-
# download, zero-byte per-package log, generic "had non-zero exit status"
# warning).
#
# To avoid that, we export `TMPDIR` to a scratch dir under the app root
# BEFORE R starts. This is the only way to redirect the parent R's own
# `tempdir()`; in-process redirection (which the install helpers also
# do, via Sys.setenv) only catches subprocesses.
#
# Override knobs:
#   TMPDIR             - if already set, we leave it alone.
#   SCE_SCRATCH_DIR    - explicit scratch dir; overrides the auto-pick.
#   SHINY_HOST / PORT  - bind address / port for runApp(). Defaults
#                        match what the README documents.
# ---------------------------------------------------------------------------

set -euo pipefail

APP_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ -z "${TMPDIR:-}" ]]; then
  SCRATCH="${SCE_SCRATCH_DIR:-${APP_ROOT}/.sce-install-scratch}"
  mkdir -p "$SCRATCH"
  export TMPDIR="$SCRATCH"
  echo "run_app.sh: exporting TMPDIR=$TMPDIR" >&2
fi

HOST="${SHINY_HOST:-0.0.0.0}"
PORT="${SHINY_PORT:-3838}"

cd "$APP_ROOT"
exec R -e "shiny::runApp('.', host='${HOST}', port=${PORT}, launch.browser=FALSE)"
