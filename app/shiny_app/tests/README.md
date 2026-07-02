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

## Auto-install R in Claude Code web sessions (SessionStart hook)

A fresh Claude Code on the web container has **no R**, so the checks above can't run
until R is installed. This repo ships a `SessionStart` hook that installs R + the
test packages automatically at the start of every web session:

- `.claude/hooks/session-start.sh` — installs `r-base-core` plus the precompiled
  `r-cran-*` Debian packages (testthat, dplyr, tidyr, tibble, lme4, lmerTest). It is
  idempotent (skips what's already present), remote-only (`$CLAUDE_CODE_REMOTE`
  guard, so it's a no-op on local machines), and always exits 0 so it never blocks
  session startup.
- `.claude/settings.json` — registers the hook under `hooks.SessionStart`.

The hook runs **synchronously**, so the session waits until R is ready (no race
where the agent runs tests before the packages exist). The `r-cran-*` debs install
fast on Ubuntu and bundle their own system dependencies (BLAS/LAPACK, gfortran
runtime) — no compilation. Verified by running the full suite (58 tests pass) after
a clean install.

The hook only enables the automated **checks**. The **full app** still cannot be
launched in a web session — it needs the gitignored `data/` artifacts (see above).
To disable the hook, remove the `SessionStart` entry from `.claude/settings.json`.
