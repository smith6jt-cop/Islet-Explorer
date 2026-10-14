#!/usr/bin/env bash
# DEVELOPMENT ONLY: R wrapper for deploy/islet-explorer-local.service. Not the
# production path; see deploy/README.md.
# Shiny Server deliberately forwards only HOME/LANG/PATH to R workers.
# Load the selected local release in the worker before R starts.
set -euo pipefail
islet_repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export ISLET_RUNTIME_MODE=development
export ISLET_DATA_ROOT="$islet_repo_root/data"
export RETICULATE_PYTHON=/home/smith6jt/miniconda3/envs/scvi-env/bin/python
export ISLET_ENABLE_AI=FALSE SHOW_STATISTICS_TAB=0 SHOW_INFLAMMATION_TAB=0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=2
if [[ -f "$islet_repo_root/tmp/local-app.env" ]]; then
  set -a
  source "$islet_repo_root/tmp/local-app.env"
  set +a
fi
exec /usr/local/bin/R "$@"
