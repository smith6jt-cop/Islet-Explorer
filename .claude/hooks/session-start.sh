#!/bin/bash
# SessionStart hook — ensure R + the Shiny app's test dependencies are available
# so `Rscript app/shiny_app/tests/parse_check.R` and `.../run_tests.R` work in
# Claude Code on the web sessions. Idempotent and fast on warm re-runs.
#
# It does NOT enable a full app launch: the app also needs sf / reticulate /
# anndata and the gitignored data/ artifacts, which are not present in a web
# session. This only enables the parse-check and unit tests.
#
# Uses Ubuntu's precompiled r-cran-* Debian packages, which bundle their system
# dependencies (BLAS/LAPACK, gfortran runtime, etc.) — no compilation, no P3M
# binary lookup. Always exits 0 so a failure never blocks session startup.
set -uo pipefail

# Web/remote sessions only — no-op on local machines.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

# Escalate with sudo only when needed, and only non-interactively — a bare `sudo`
# would block session startup on a password prompt, defeating the "never block"
# intent. If we're not root and passwordless sudo isn't available, warn and no-op.
SUDO=""
if [ "$(id -u)" -ne 0 ]; then
  if command -v sudo >/dev/null 2>&1 && sudo -n true >/dev/null 2>&1; then
    SUDO="sudo -n"
  else
    echo "[session-start] Not root and passwordless sudo unavailable — skipping R setup (would block on a password prompt)." >&2
    exit 0
  fi
fi

APT_UPDATED=""
apt_install() {
  export DEBIAN_FRONTEND=noninteractive
  if [ -z "$APT_UPDATED" ]; then
    $SUDO apt-get update -qq || { echo "[session-start] WARNING: apt-get update failed." >&2; return 1; }
    APT_UPDATED=1
  fi
  $SUDO apt-get install -y --no-install-recommends "$@"
}

# --- 1. Ensure R is installed ---
if ! command -v Rscript >/dev/null 2>&1; then
  echo "[session-start] R not found — installing r-base-core..."
  if ! apt_install r-base-core; then
    echo "[session-start] WARNING: could not install R; app tests unavailable this session." >&2
    exit 0
  fi
fi

# --- 2. Ensure the test packages are present (idempotent) ---
need=""
for p in testthat dplyr tidyr tibble lme4 lmerTest; do
  if ! Rscript -e "quit(status = as.integer(!requireNamespace('$p', quietly = TRUE)))" 2>/dev/null; then
    need="$need r-cran-$(printf '%s' "$p" | tr '[:upper:]' '[:lower:]')"
  fi
done

if [ -n "$need" ]; then
  echo "[session-start] Installing R test packages:$need"
  # shellcheck disable=SC2086
  apt_install $need || echo "[session-start] WARNING: some R test packages missing; a few tests may not run." >&2
else
  echo "[session-start] R test packages already present."
fi

echo "[session-start] R environment ready for app/shiny_app tests."
exit 0
