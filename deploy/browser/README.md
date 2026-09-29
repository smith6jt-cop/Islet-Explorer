# Browser acceptance

Use `check-app.mjs` below for the restored 15-donor development app and for a
matching selected release. Its URL can be the shared local development URL or
an isolated loopback Shiny worker.

The schema-3 all-donor check is historical. Its prior invocation was:

```bash
node deploy/browser/check-all-donor.mjs \
  --url http://127.0.0.1:3839/islet-explorer/ \
  --output /tmp/islet-all-donor-browser.json \
  --chrome /usr/bin/google-chrome
```

It checked all 21 donor/image choices, ratio-only trajectory controls, fixed
coordinates after donor filtering, graph/ratio/heatmap downloads, raw-cell
drill-down and its donor scope, HPAP click-through, corrected-6539 TIFF Range,
all 31 heatmap defaults/exports, the visible islet-review status, absence of
research UI/downloads/routes, and axe on the four app tabs. Empty scientific results
were expected where validated inputs were unavailable. It waits for visible
reactive updates to settle before evaluating accessibility. Failure evidence
includes offending DOM fragments; keep reports in the restricted local evidence
directory. Passing the explicit unavailable-ratio state did not establish
working real ratio curves or accepted islet identities. The withdrawn 2026-09-18
corrected-boundary check is retained under
[`historical/check-reviewed-islets.mjs`](historical/check-reviewed-islets.mjs).
Neither historical harness is a gate for the restored app or phase-2 additions.

Run this test against an explicitly selected, isolated or staged four-tab app.
It exercises the Shiny session, keyboard tab navigation and skip link, Plot,
Trajectory/PAGA and Spatial outputs, four CSV downloads, the actual image
viewer canvas, a 16-byte TIFF Range request, and missing/external viewer URLs.
It runs axe on each tab and the viewer iframe. Serious or critical automated
violations, uncaught browser errors, or attempted external requests fail the run.

Use Node >=22.12 and an installed Chrome/Chromium executable. Install the pinned
browser dependencies with `npm ci --prefix deploy/browser`, then run:

```bash
node deploy/browser/check-app.mjs \
  --url http://127.0.0.1:3838/islet-explorer/ \
  --output /tmp/islet-browser-acceptance.json \
  --chrome /usr/bin/google-chrome
```

The URL must contain no credentials, query string, or fragment. Use the local
backend or an already authorized test endpoint. This harness does not establish
UFRC proxy authentication or external reviewer access. It blocks HTTP(S)
requests to other origins while testing the application's normal behavior.

The JSON evidence contains rule IDs/counts, check results, HTTP status/byte
counts, and browser/dependency versions. It does not save page text, scientific
rows, CSV bodies, image filenames, screenshots, console messages, or credentials.
The Range check cancels a non-206 response before reading its body. CSVs are
checked in memory and discarded. Output is written with mode 0600.

The workstation review identified a browser-generated `/favicon.ico` request
returning 404 as the sole console error. The report counts that missing default
icon separately from application HTTP errors; application errors still fail.
The review also found absent vendor tab-button IDs and an absent iframe
landmark. The wrapper now connects both panels to their real tab buttons and
provides a named main landmark; the harness checks those relationships.

For a workstation test, run a separate unprivileged Shiny Server with a copy of
`deploy/shiny-server.conf`: substitute the current user, an unused loopback
port, the local application/image paths, and temporary log, `socket_dir`, and
`bookmark_state_dir` directories. Keep all three explicit static image routes;
`shiny::runApp()` and the general `app_dir` route alone returned full TIFF HTTP
200 responses during testing. Do not restart the managed application.

Start that isolated process with `R_ENVIRON_USER=/dev/null`, the tested
`RETICULATE_PYTHON`, `ISLET_ENABLE_AI=FALSE`, `SHOW_STATISTICS_TAB=0`,
`SHOW_INFLAMMATION_TAB=0`, `DEBUG_CREDENTIALS=0`, and `VIEWER_DEBUG=0`; remove
`KEY` and `OPENAI_API_KEY` from its environment. Stop the isolated server after
testing. For final acceptance, point the same harness at the selected immutable
image and verified release; a workstation development-data run does not attest
to their identity or parity.

An automated pass is a subset of accessibility review. Axe's `needsReview`
entries remain in the evidence, particularly for scientific visualizations.
Screen-reader use, visualization interpretation, zoom/reflow, full keyboard
interaction, and the external Basic Auth/proxy behavior still require manual
review. The browser uses software WebGL to make headless image rendering
repeatable; this does not measure GPU performance on a reviewer's workstation.
