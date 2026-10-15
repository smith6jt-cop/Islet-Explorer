# Reviewed HPAP boundaries in the combined app candidate

The independent candidate is `data/development/all-donor-20260918/`. It replaces
all 1,234 previous HPAP-078 islet identifiers with the 1,431 owner-approved
endocrine boundaries from the six-donor QuPath review delivery. The combined
candidate contains 8,568 islets, 21 donors (seven per group), 3,592,722 raw cell
objects, and all 31 shared markers. All 201 wholly unmeasured objects are retained.

The other 20 donors' measurements and geometry are preserved. For five additions,
2,030 regenerated candidates still await visual review; the app explicitly labels
their displayed outlines as the previous boundaries. Original reference-donor
identity reviews are also pending. Boundary approval does not validate phenotype
gates. Subtype ratios remain unavailable wherever required gates are uncertified.

## Data and application behavior

HPAP's indexed geometry, spatial lookup, core/peri-islet regions, marker medians,
structural area summaries, trajectory inputs and all 1,431 per-islet cell files
use the approved geometry. The old HPAP identifiers and summaries are excluded
from the candidate, and no stale master workbook is included. Exact native
vertices and QuPath's ROI-specific containment conventions preserve all 9,569
verified native core memberships, including the former edge-rounding discrepancy.

The previous trajectory is invalidated, not transferred to the corrected islets.
Coordinates and percentiles are missing, graph edges are empty, and previous
graph/model/landmark files are absent. The status is `awaiting_reviewed_refit`;
200 donor bootstraps and leave-one-donor-out validation have not been rerun.
The technical feature eligibility count is 3,063 and is not a reviewed fit count.

The Trajectory tab offers direct islet search, corrected outline/cell inspection,
single-islet GeoJSON and cell downloads, and click-through to the correct image.
All 31 heatmap markers remain available. While pseudotime is unavailable, the
heatmap uses explicitly labeled fixed cohort-wide section-area bins. Donor and
cell-count filters apply to the heatmap, selector and all tabular downloads.
Research results remain outside the app.

## Efficient combined-data runtime

Images remain immutable hardlinks in this local candidate and use HTTP byte-range
requests; they are not loaded wholesale into R. Geometry inspection uses spatially
indexed GeoPackage bounding-box queries and a bounded region cache. Single-islet
outline downloads copy prebuilt exact-vertex files without parsing donor-wide
GeoJSON. Cell data load only for the chosen islet and use a process-wide cache
limited to 64 entries and 64 MiB. Compact Parquet trajectory tables are shared
across sessions within an R worker, with candidate-specific cache identity and
copy-on-modify behavior. Selectize searches islet choices on the server. Display
filters never fit, project or rescale a trajectory.

## Rebuild and evidence

Run `scripts/research/integrate_reviewed_boundaries.py` with explicit source and
fresh output roots, the verified QuPath project, regional root and exact geometry
root. Then run `build_cell_lookup.py --case 78`, `index_candidate_geometry.R` for
case 78, `finalize_candidate.py`, `optimize_candidate_runtime.py`, and
`validate_candidate.py`. The optimizer only works on a candidate that has not
been served. Its changed small-artifact hashes are refreshed after this candidate's
independent full hash pass; original-release passing evidence is not reused.

The candidate contains `revision.json`, `runtime_optimization.json`,
`candidate_validation.json`, `manifest.json`, `capacity.json`, and
`trajectory/boundary_review.csv`. The historical browser check is
`deploy/browser/historical/check-reviewed-islets.mjs`. Service identity, logs and exact
restart arguments are recorded in `preview/service.json`. Preserve the preceding
`all-donor-20260917` candidate and its service record for rollback.

## Verified local delivery

The preview is **http://127.0.0.1:3839/islet-explorer/**; its current PID is recorded
in `preview/service.json`. Port 3838 remains managed by its existing service.
Port 3839 was free at startup; no previous service was terminated to claim it.
The candidate's own first test process was restarted once to apply an
empty-selection fix. The final browser check passed on 2026-09-18 at 20:57:46 UTC.

- `evidence/browser-after-empty-fix.json`: all 21 donors, all 31 markers, 1,431
  corrected HPAP rows, selected-donor downloads, exact outline/cell exports,
  image click-through/canvas, subsequent donor switching, HTTP range reads,
  empty-selection/marker recovery, no Research UI/route, and no browser errors.
  All four tabs had no WCAG A/AA violations in the tested state.
- `candidate_validation.json`: 51,408 compartment/band records, marker coverage,
  missingness, CSV/Parquet parity, all 9,569 native core memberships, exact HPAP
  outline vertices and preservation of the other 20 donors' measurements passed.
- R suite: 267 assertions passed, zero failures/warnings/skips; pipeline suite:
  21 tests passed. R parse check and `git diff --check` passed.
- The first browser run and its failure are retained in `evidence/browser.json`.
  Its empty-selection observer bug was fixed and covered by a regression test;
  the later passing evidence supersedes that run.

| Local measurement | Result |
|---|---:|
| Browser app ready, including new R worker | 4.175 s |
| First 31-marker heatmap | 2.018 s |
| HPAP donor selection and heatmap | 0.300 s |
| Islet selection, outline and cells | 1.514 s |
| Single-islet outline download | 0.007 s |
| Core manifest verification | 0.325 s |
| Five trajectory tables, CSV baseline | 1.089 s |
| Five trajectory tables, Parquet first read | 0.382 s |
| Cached table reads | ≤0.001 s in 20 repetitions |
| Indexed selected-region geometry query | 0.156 s |
| Selected-islet cell file read, 118 rows | 0.002 s |

The five CSV tables total 45,999,416 bytes; their Parquet equivalents total
4,486,081 bytes. Shared table objects occupy 27,760,256 bytes in R. The tested
geometry query loads 301 nearby objects (1,454,600 bytes in R), rather than the
entire 63,098,880-byte HPAP GeoPackage. All 8,568 exact single-islet outline files
total 31,608,198 bytes. Measurements are local single-worker observations with
possible warm filesystem caches; they do not establish concurrent-user capacity.

Runtime artifacts total **164,198,139,742 bytes (152.921 GiB)**, including
163,098,139,438 image bytes. Two independent deployment copies require
328,396,279,484 bytes, before frozen releases, reports, container images, logs
and spill. Local image/source hardlinks avoid duplicating those bytes; deployment
copies still need the full allocation. Full-worker peak memory and concurrency
have not been measured for this revision. Scientific and UFRC acceptance remain false.

Restart the current preview:

```bash
python3 scripts/research/serve_preview.py \
  --data-root data/development/all-donor-20260918 \
  --app-root app/shiny_app \
  --python /home/smith6jt/miniconda3/envs/scvi-env/bin/python --restart
```

For rollback, use the previous candidate **and its preserved previous application
snapshot**, replacing only this recorded preview. Do not claim its old trajectory
is valid for the corrected boundaries:

```bash
python3 scripts/research/serve_preview.py \
  --data-root data/development/all-donor-20260917 \
  --app-root data/development/all-donor-20260917-rollback-app \
  --python /home/smith6jt/miniconda3/envs/scvi-env/bin/python \
  --restart --replace-preview data/development/all-donor-20260918
```

The preserved app snapshot was copied into the separate rollback source directory
above; the builder rejects a source located inside its destination app directory.
Logs are under the selected candidate's `preview/logs/`, with server
events in `preview/server.log`; the service JSON records the actual listener URL.

## Network cutover

On 2026-09-18, the owner requested replacement of the old network app with this
corrected candidate. The shared URL is **http://10.15.152.7:8080/islet-explorer/**.
Browser visits to **http://10.15.152.7:3838/islet-explorer/** forward to the same
host on port 8080. The nginx endpoint supplies HTTP 206 image byte ranges; the
direct Shiny static route returned full TIFF responses, so it is not the viewer
endpoint. Users need network access to this workstation; no external-host or
off-network connectivity claim is made.

The root-owned Shiny service and nginx configuration were preserved. Its existing
app symlink points to `app/shiny_app`; that app now has an ignored workstation
profile selecting the explicit corrected data root. Its `www/local_images` path
is a symlink to the candidate's 21-image directory, also used by nginx's existing
static-image alias. Six required runtime files received read ACLs for the `shiny`
service account. Credentials were neither copied into the rollback source nor
made readable. Touching the app's `restart.txt` selected a fresh worker without
restarting the separate Senior app.

The original image directory and pre-cutover application source are preserved in
`data/development/network-cutover-20260918/`; its `runtime-acls-before.txt` records
the preceding ACLs. `network_service.json` in the candidate records data/profile
identity, URLs, backup paths, service-account access and verification results.
The loopback preview on 3839 was also refreshed with the Viewer-warning removal.

`evidence/browser-network.json` passed against the network address on port 8080.
`evidence/browser-network-redirect.json` passed starting at port 3838 and resolving
to port 8080. Both checked all 21 donors, all 31 markers, 1,431 corrected HPAP
islets, raw-cell and exact-outline downloads, image canvas/range streaming,
donor switching, no Viewer boundary warning, empty-selection recovery and no
Research route. No browser errors or tested WCAG A/AA violations were reported.
Network app readiness was 3.930 s; donor filtering was 0.329 s. One live network
R worker was observed at 813,372 KiB RSS (~0.78 GiB); this is a point observation,
not a peak or concurrent-user capacity measurement.

The Viewer notice was removed because Avivator displays image channels without
islet outlines. The Trajectory tab retains boundary-review status beside its
actual outline inspection tools. Scientific and UFRC acceptance remain false.
