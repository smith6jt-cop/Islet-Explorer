# Shiny App Tests & CI

Automated checks for the R Shiny app (`app/shiny_app/`). Run in CI by
`.github/workflows/shiny-app-ci.yml` on every push/PR that touches
`app/shiny_app/**`.

## What runs

1. **Parse check** — `tests/parse_check.R` syntax-parses every `R/*.R` file plus
   `app.R`. Catches syntax errors without loading any package or data.
2. **Unit tests** — `tests/run_tests.R` runs the `testthat` suite in
   `tests/testthat/`. Currently covers `R/utils_stats.R` — the pure statistical
   helpers behind the Statistics tab (donor-level aggregation, mixed-effects test,
   per-bin ANOVA / Kendall τ / pairwise, Cohen's d, η², normality, outlier
   flagging). These need no app data and no Python bridge.

## Run locally

```bash
cd app/shiny_app
Rscript tests/parse_check.R      # syntax check
Rscript tests/run_tests.R        # unit tests
```

Dependencies for the tests: `testthat`, `dplyr`, `tidyr`, `tibble`, `lme4`,
`lmerTest`. Install with:

```r
install.packages(c("testthat", "dplyr", "tidyr", "tibble", "lme4", "lmerTest"))
```

## Why CI does not launch the full app

The complete app needs (a) heavy runtime deps — `sf` (GDAL/GEOS/PROJ),
`reticulate` + `anndata` (a Python bridge), `shiny`, `plotly` — and (b) the data
artifacts in `data/` (H5AD ~70 MB, per-islet/per-donor CSVs ~280 MB, GeoJSON,
`phenotype_rules.csv`). `data/` is **gitignored**, so it is absent in CI and in a
fresh clone, and `R/00_globals.R` reads some of it at startup. CI therefore
parse-checks the app files and unit-tests the data-independent utilities.

To launch the full app locally (with the data present):

```bash
cd app/shiny_app
Rscript -e 'shiny::runApp(".", port = 7777)'
```

## Adding tests

Drop `test-*.R` files in `tests/testthat/`. Dependencies and the code under test
are loaded once by `tests/testthat/helper-setup.R` (sourced automatically before
the test files). Keep new tests data-independent so they run in CI.

## Optional: auto-install R in Claude Code web sessions (SessionStart hook)

By default a fresh Claude Code on the web container has **no R**, so the agent
can't run the checks above without installing it first. To make R + the test
packages available automatically at the start of every web session, add a
`SessionStart` hook. (This is opt-in agent configuration, so it is documented
here rather than committed as an active hook — enable it yourself if you want it.)

1. Create `.claude/hooks/session-start.sh` (and `chmod +x` it):

   ```bash
   #!/bin/bash
   # Ensure R + the Shiny app's test deps are available in Claude Code web sessions.
   # Idempotent; enables tests/parse_check.R and tests/run_tests.R (not a full app
   # launch — that also needs sf/reticulate/anndata and the gitignored data/).
   set -uo pipefail
   [ "${CLAUDE_CODE_REMOTE:-}" != "true" ] && exit 0

   SUDO=""; [ "$(id -u)" -ne 0 ] && command -v sudo >/dev/null 2>&1 && SUDO="sudo"

   if ! command -v Rscript >/dev/null 2>&1; then
     export DEBIAN_FRONTEND=noninteractive
     $SUDO apt-get update -qq && \
       $SUDO apt-get install -y --no-install-recommends r-base-core || {
         echo "[session-start] WARNING: R install failed; tests unavailable." >&2; exit 0; }
   fi

   # Precompiled Debian packages pull all system libs (no gfortran/libuv yak-shave).
   need=""
   for p in testthat dplyr tidyr tibble lme4 lmerTest; do
     Rscript -e "quit(status = !requireNamespace('$p', quietly=TRUE))" 2>/dev/null || \
       need="$need r-cran-$(echo "$p" | tr '[:upper:]' '[:lower:]')"
   done
   if [ -n "$need" ]; then
     export DEBIAN_FRONTEND=noninteractive
     $SUDO apt-get install -y --no-install-recommends $need || \
       echo "[session-start] WARNING: some R test packages missing." >&2
   fi
   echo "[session-start] R ready for app/shiny_app tests."
   ```

2. Register it in `.claude/settings.json` (merge if the file exists):

   ```json
   {
     "hooks": {
       "SessionStart": [
         { "hooks": [ { "type": "command",
             "command": "$CLAUDE_PROJECT_DIR/.claude/hooks/session-start.sh" } ] }
       ]
     }
   }
   ```

Notes: the hook runs synchronously, so the session waits until R is ready (no race
where the agent runs tests before packages exist). On Ubuntu the `r-cran-*`
packages install fast and bundle their system dependencies. This was verified in a
web session by installing R this way and running the full suite (58 tests pass).
The **full app** still cannot be launched in a web session — it needs the
gitignored `data/` artifacts — so the hook only enables the automated checks.
