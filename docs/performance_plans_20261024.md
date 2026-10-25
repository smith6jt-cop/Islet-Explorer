# Making the reviewer app as fast as possible on PubApps: two plans (2026-10-24)

Companion to the [responsiveness baseline](performance_baseline_20261024.md), which
established that the hosting (VM, NFS data, UFRC proxy) adds nothing measurable and that
reviewers feel the app's own payload over their link, two one-time computations per R
process, and the single shared R process. This document turns those findings into
**Plan A**, quick wins that keep the Shiny app and the current hosting, ordered by
feasibility, and **Plan B**, a complete overhaul that removes the server and the link from
the interaction path. Neither plan changes scientific outputs; every item is verified with
the same `bench-app.mjs` scenario against the baseline tables.

## What the measurements say (inputs to both plans)

Where a Plot interaction's 1.1 s of server time goes (Rprof over six measure changes, a
drill-down, the Spatial tab and three donor changes; `prof-summary.txt` in the evidence
directory):

| Share of sampled server time | What |
|---:|---|
| 60 % | `renderPlotly` of the two Plot cards (`card_figure`: `status_figure` 22 %, `diameter_figure` 20 %) |
| 36 % | of which `ggplotly()` conversion itself (`gg2list` / `ggplotly_build`): walking the ggplot build, one R-level loop per layer and point |
| 16 % | `renderPlot` PNG drawing (`drawPlot`, grid): the tissue scatter and the drill-down map |
| 5 % | GeoPackage read for the drill-down (`CPL_read_ogr`) |
| 13 % | `data.table::fread` of the development CSV backend (absent in production, which reads Parquet) |

Nothing in that profile is data access or statistics; it is figure construction and
serialisation in R. The payload those figures produce, and what compression would do to
it (measured on staging):

| Message | Raw | Deflate (level 6) | Brotli | Composition |
|---|---:|---:|---:|---|
| First screen | 3.56 MB | 0.66 MB (5.4×) | 0.49 MB | two renders of the Plot cards (1.8 MB + 1.5 MB) |
| One Plot change | 1.70 MB | 0.32 MB (5.3×) | 0.24 MB | hover `text` 58 %, `x`/`y` at full precision 27 %, `key` 10 % |
| Spatial tab | 1.67 MB | 0.60 MB (2.8×) | 0.44 MB | base64 PNG of 170 k cells plus two UMAPs |

Fixed costs: first Trajectory open 12 s and first PAGA 5 s per R process (scanpy, cached
in-process afterwards); R process start 7 s (hidden by the unit's health check). Hosting
facts that bound any plan: 2 vCPUs, 15.7 GB, one R process per app (open-source Shiny
Server), 9 GiB cap, R resident at 2.3 GB, WebSocket through Apache `mod_proxy_wstunnel`,
Apache serving the DocumentRoot and the image aliases itself.

**Targets.** "Fast" for this app means: first screen under 2 s on campus and under 4 s on
a 25 Mbps home link; a control change under 0.5 s on campus and under 1 s at home; no
one-time cost a reviewer can meet; five simultaneous reviewers without a step waiting more
than a second on another's work.

## Plan A: quick wins (Shiny app and hosting unchanged)

### Level 1: hours each, no architecture change, no UFRC involvement

**A1. Compress the WebSocket (permessage-deflate) in the project's Shiny Server build.**
The project already builds Shiny Server from a pinned source tarball with a pinned
`package.json` and shrinkwrap (`deploy/harden-base.sh`, `deploy/shiny-server/`). SockJS
passes `options.faye_server_options` straight to `faye-websocket`, whose `extensions`
option enables permessage-deflate; Shiny Server's `lib/proxy/sockjs.js` currently passes
five options to `sockjs.createServer` and no `faye_server_options`. The change is a
patch step in `harden-base.sh` adding
`faye_server_options: { extensions: [require('permessage-deflate')] }` to that call and
`permessage-deflate@0.1.7` (one dependency, `safe-buffer`) to the pinned manifest.
Browsers negotiate the extension themselves; the proxy tunnels the handshake. Measured
effect: 1.70 MB → 0.32 MB per Plot change, 3.56 → 0.66 MB first screen, Spatial
1.67 → 0.60 MB; on the 25 Mbps profile a Plot change goes from 2.9 s to about 1.5 s and
the first screen from 8.8 s to about 4.5 s. Cost: a few tens of milliseconds of CPU per
message on the VM and about 300 KB of zlib state per connection (bounded with
`maxWindowBits`/`noContextTakeover` if memory matters). Verify on staging: the handshake
response must carry `Sec-WebSocket-Extensions: permessage-deflate` through the proxy
(`bench-app.mjs` needs a small addition to record it, since DevTools reports decompressed
frame sizes), the upstream Shiny Server tests still pass in the image build, and the
proxied step times fall under the home profile.

**A2. Stop the Plot tab's startup double render.** `output$dynamic_selector` and
`output$metric_selector` are `renderUI` outputs that create `input$which` and
`input$metric`; the cards render once on the fallback value and again when the mounted
inputs report. Call `freezeReactiveValue(input, "which")` and
`freezeReactiveValue(input, "metric")` inside those `renderUI`s so dependents wait for the
real value (the standard Shiny pattern for inputs created by `renderUI`), or build the two
selectors statically and update their choices. Effect: one full render less per page load
(1.7 MB and about 1.1 s of server time; 2–3 s on a home link). Verify with the baseline's
`valueTimeline`: each Plot output arrives once.

**A3. Warm the Trajectory computations at process start.** Schedule
`compute_scanpy_pseudotime(TRAJ_MIN_CELLS, Inf)` and the default `paga_result()` from
`app.R` right after `.shared_prepared_data` loads, in a `later::later()` call so the first
page is not delayed. The unit's health check already keeps the process alive and the
start period is 180 s, so the 17 s move to where no reviewer waits. Effect: first
Trajectory open 12 s → 0.2 s, first PAGA 5 s → 1 s for the first reviewer after each
activation. (Level 2 item B2 removes the computation from the runtime altogether.)

**A4. Round the figure numbers.** Before `ggplotly()`, round the plotted columns to the
precision the axes show (`signif(x, 4)`): the `x`/`y` arrays shrink from 0.46 MB to about
0.2 MB per interaction, and compress better. One helper applied in `diameter_figure()`,
`status_figure()`, the Trajectory scatter and the Spatial UMAP builders.

**A5. Trim the hover text.** `ggplotly()` builds each point's hover from every mapped
aesthetic (`diam_mid`, `value`, `donor_status`, `islet_click_key`, …), which is the 1 MB
of `text` per interaction. Supply a compact `text` aesthetic in R (donor, object, diameter,
value) and pass `tooltip = "text"`, keeping the `key` aesthetic the linked-plot script
reads. Effect: about 1 MB → 0.25 MB per interaction before compression; with A1 a Plot
change is then roughly 150 KB on the wire.

**A6. Serve a partial plotly.js bundle.** `plotly::partial_bundle(type = "cartesian")`
on each widget replaces the 1.09 MB gzipped `plotly-latest.min.js` with about 0.4 MB
(cartesian covers scatter, bar, box, violin, heatmap, histogram, contour, image, pie).
First visit only, since the browser caches it afterwards: about 0.3 s on 25 Mbps, 0.7 s on
8 Mbps. Confirm every trace type across the four tabs is in the bundle before switching.

**A7. Put the small data files on the VM's local disk at activation.** `geometry/`
(705 MB of GeoPackages), `cells/` (119 MB) and `spatial/` (56 MB) are read with many
small operations, and the GeoPackage window read is 0.2 s over NFS against 17 ms on local
disk. `activate-release.sh` can copy those three trees from the verified release to
`/local/maigan/releases/<name>/` and mount them in place of the NFS paths (images stay on
NFS and behind Apache). Effect: drill-down 0.75 s → about 0.55 s and steadier donor
changes. Half a day in tooling plus the existing verification.

**Level 1 status (2026-10-25: on the staging instance, commit `e10b2ee`, release
`pubapps-2026-10-uniform-v10-perf1`).** A1 (compression in the Shiny Server build), A2
(single first render), A3 (warm-up of the H5AD, pseudotime and PAGA at process start), A4
(five-significant-figure coordinates on the Plot cards, the Trajectory scatter and the
UMAPs), A5 (compact per-object hover; the violin outlines carry no hover) and A7
(`activate-release.sh --local-data`, with its recovery test) were implemented in `eba12b4`
and are deployed on `islet-explorer-dev` (`devislet-explorer.rc.ufl.edu`) since 2026-10-25
16:03 UTC, under a release name of their own beside production, which stays on
`pubapps-2026-10-uniform-v10` (`d55680a`). Development is in-house: nothing moves to
production as part of this plan until the owner decides. The first image build (`4da5ac4`)
failed the High/Critical gate on an advisory published after the previous scan
(body-parser 2.3.0, GHSA-27r3-jxmh-p828); `e10b2ee` pins 2.4.0 (`deploy/BASE_HARDENING.md`).
Measured through the staging URL (`performance_baseline_20261024.md`, "Plan A Level 1 on
staging"): `permessage-deflate` negotiated on every connection; first screen 3.5 → 1.1 MB
over the socket and sent once; first Trajectory open on a fresh process 12.7 s → 0.26 s and
first PAGA 5.2 → 1.2 s; warm revisit of the Plot tab 5.8 → 2.5 s on that instance; the
per-change server compute is unchanged, as expected (Level 2's `plot_ly()` item). A6
(partial plotly bundle) stays deferred: `plotly::partial_bundle()` fetches the bundle from
a CDN at render time, which the production container must not do, so it needs the bundle
vendored into the image first, for a first-visit-only gain. Next: Level 2 items, measured
on the same staging instance.

### Level 2: days each, app changes measured on staging

**B1. Build the Plot cards and the Trajectory scatter with `plot_ly()` directly.**
`ggplotly()` conversion is the single largest cost per interaction (36 % of server time,
close to 90 % of each Plot step). The figures are simple: markers with jitter, summary
points with error bars, LOESS lines, box or violin layers. Built with `plot_ly()` and
`add_markers()`/`add_lines()` on the data frames the module already prepares, the same
figures take 50–150 ms, produce smaller JSON (no generated hover text, no duplicated
legend traces) and keep `key` on the marker traces for the linked plots. Effect: Plot
change server time 1.1 s → about 0.2 s, which is also five times the concurrency headroom.
Two to four days including the browser checks (`check-app.mjs`, `check-ui-review`,
`check-trajectory-links`) and the download PNGs, which already use `islet_ggsave()` and
stay as they are.

**B2. Precompute the diffusion pseudotime and the PAGA layout in the release bundle.**
Both depend only on the fixed ten-cell cohort and the embedding, so they are build-time
artifacts (a Parquet table of object, pseudotime, rank, PAGA coordinates and edges) under
the next data schema. The runtime then has no scanpy call on the Trajectory path, the
first open is instant for every process, and the Python side of the R process shrinks.
One day in `scripts/build_release_data.R` and the schema validator, plus the parity check.

**B3. Cache rendered outputs across sessions.** Shiny 1.11 `bindCache()` with a shared
`cachem::cache_mem()` on the Spatial tissue scatter (keyed by donor and the drawing
options), the composition-by-donor card and the Plot figures keyed by their complete input
set. Reviewers open the same defaults, so a second viewer gets a 1.1 s figure in a few
milliseconds and the shared process is freed for the first. One day, mostly in choosing
complete cache keys (the figure reactives already enumerate their inputs for the CSV
provenance headers, which is the checklist to reuse).

**B4. Lighter tissue scatter.** The 800 px `renderPlot` PNG of 170 k cells is about 1 MB
of base64 per view and 16 % of server time. Options in increasing effort: draw with the
`ragg` device (two to three times faster than `png()`), cap the PNG at the displayed size
and a 1.5× device ratio, or pre-render the default view per donor at build time and serve
it with `renderImage` (the sidebar options then re-render on demand only). Expected:
Spatial tab 2.2 s → about 1.2 s on campus and a smaller message.

### Level 3: needs UFRC or infrastructure changes

**C1. More than one R process.** Open-source Shiny Server serves all sessions of an app
from one process, so simultaneous reviewers queue (two sessions doubled every step in the
baseline). Three routes: (a) several `location`s in `shiny-server.conf` pointing at the
same app under different paths, balanced with sticky sessions at UFRC's proxy
(`ProxyPass balancer://` with a route cookie, WebSocket-aware): a ticket to UFRC and
careful testing of the SockJS paths; (b) an in-container balancer (nginx or a small Node
proxy on 3838) in front of two or three Shiny Server workers, no proxy change, bounded by
2 vCPUs and 2.3 GB per R process inside the 9 GiB cap; (c) ask UFRC for a larger VM
(4–8 vCPUs, more memory) through the ticket, which (a) or (b) then use. Any of these only
pays after B1 reduces the per-step compute; with B1 and B3 in place, two processes on two
cores serve five reviewers comfortably.

**C2. Proxy-level HTTP tuning** (HTTP/2, static caching headers for the asset paths) is
worth at most a few hundred milliseconds on the first visit and can wait.

### What Plan A delivers

| Step | Today, campus | Today, 25/5 Mbps | After Level 1 (A1–A7) | After Levels 1–2 |
|---|---:|---:|---:|---:|
| First screen | 3.6 s | 8.8 s | about 2 s / 3 s | about 1.5 s / 2 s |
| Plot change | 1.2 s | 2.9 s | about 1.1 s / 1.4 s | about 0.3 s / 0.5 s |
| Open Spatial | 2.2 s | 3.8 s | about 2.1 s / 2.6 s | about 1.2 s / 1.5 s |
| First Trajectory open after activation | 12.7 s | 12.7 s | 0.2 s | 0.2 s |
| Two simultaneous reviewers | every step ×2 | | | queueing mostly gone with B1+B3; removed with C1 |

Level 1 is one to two working days in total, Level 2 one to two weeks; each item is
independently shippable and each is measured with `bench-app.mjs` on staging before
production (the image rebuild for A1 and A6 goes through the existing scan gate).

## Plan B: complete overhaul

Plan A reaches the targets for one reviewer on a fast link and most of them on a home
link, but three limits remain while the browser asks an R process for every change: a
round trip plus server compute on every control, a payload per interaction, and a server
whose capacity reviewers share. The overhaul removes all three by moving the interactive
work to the browser and leaving the server to serve files.

### Option B-I: thin-server Shiny (data once, render in the browser)

Keep Shiny for the page, the statistics and the drill-down images, but send each tab's
object table to the browser once per session (the release summary Parquet is 1.6 MB in
total; as Arrow IPC or compact JSON about the same compressed) and build the Plot,
Trajectory scatter and Spatial UMAP figures in JavaScript (Plotly.js is already loaded).
Filters, measure changes, palette and axis options then cost zero round trips; only
statistics, downloads and drill-down call R. Effort: three to five weeks (a JS figure
module per tab with the same hover, linking and theme contracts the R builders implement
now; `trajectory-links.js` already owns the client-side linking). Concurrency improves
because R is called rarely, but the single process and the Shiny session remain.

### Option B-II: a static reviewer app served by Apache (no R process at runtime)

The PubApps proxy already serves the DocumentRoot and byte ranges, and Avivator is
already a static bundle. A static reviewer app would be: a JavaScript application that
reads the release's Parquet tables through DuckDB-WASM (or pre-aggregated JSON per tab)
with HTTP range requests, draws with Plotly.js or deck.gl, opens drill-downs from
per-object GeoJSON and cell files cut at build time (the release already ships one cells
CSV per object), and computes the light statistics in the browser (Kruskal-Wallis,
Spearman, LOESS, medians: `simple-statistics` and `loess` ports exist and can be validated
against the R results in the parity suite). Everything heavy, including DPT, PAGA, the
composition embedding and the signed-distance KDEs, is precomputed at build time, which
the release pipeline already does for most of them.

What this buys: interaction latency independent of the server and of the link after the
first load (a few megabytes once, cached), no R process to start, warm or share, no
concurrency ceiling (Apache serves files), a far smaller attack surface for UFRC's scan,
and a deployment that is a directory under `prod/html`. What it costs: a rewrite of the
four reviewer tabs with parity to the reviewed app (figures, help text, downloads with
provenance headers, themes, keyboard navigation, the accessibility checks), roughly two
to three months for one developer, plus a new acceptance round. Shiny would stay as the
analyst's local tool and as the reference the parity checks compare against.

### Option B-III: keep R, change the serving stack

Several `shiny::runApp()` processes behind an in-container sticky balancer, Shiny
`promises` for the long computations so one session's PAGA does not block another's
click, a persistent `cachem::cache_disk()` so figures survive restarts, and the Plan A
items. This fixes concurrency and cold costs but not the per-interaction round trip or
payload; it is the fallback if B-I or B-II are not undertaken and the review round grows.

### Recommendation and sequence

1. **Do Plan A Level 1 now** (A1 compression, A2 double render, A3 warm-up, A4/A5 payload
   trims, A7 local-disk data). Together they roughly halve what an off-campus reviewer
   waits for, in days, with no change to the hosting and no scientific change.
2. **Do B1 (direct `plot_ly`) and B2 (precomputed trajectory) next.** They cut the server
   work per step by about five times and remove the last fixed cost; both are also
   prerequisites that carry over into any overhaul (the figure specifications and the
   precomputed artifacts are what the browser-side app would consume).
3. **Decide on the overhaul after the review round**, with the baseline re-measured. If
   the public release needs to serve many concurrent visitors or the reviewers' links stay
   slow, B-II (static app) is the durable answer; B-I is the intermediate step if R must
   stay in the loop for statistics. Only pursue C1 (more processes) if concurrency bites
   before that decision.

### Would more HiPerGator resources help? Would a GPU?

**A larger VM helps concurrency only, and only with more than one R process.** The
baseline measured the VM's core as fast per thread as the workstation's (sort 328 ms
against 379 ms), and every step a reviewer waits for runs single-threaded in one R
process. More vCPUs therefore leave a single reviewer's times unchanged, and with one
process the extra cores sit idle. They pay once C1 runs two or three R processes: each
needs about 2.3 GB resident, so 4 vCPUs and 24 GB would carry three processes with room
for the Podman store and the page cache. That is a line in the PubApps ticket, not a
HiPerGator allocation: the reviewer app must run on the PubApps VM, because UFRC's proxy
reaches only that VM and HiPerGator compute nodes have no public ingress. HiPerGator
could take build-time work (the 66-minute data bundle, parity checks, image scans), but
the workstation's 32 cores already do that faster than a typical node, so there is no gain
there either. Nothing on the server side changes what an off-campus reviewer feels; that
is payload over their link, which only A1, A4, A5 and the Plan B options address.

**No GPU helps on the server, and the reviewer's GPU is already in use.** The server's
time goes to R-level figure construction (`ggplotly()`), JSON serialisation, PNG
rasterisation and message transfer; none of it runs on a GPU, R's plotting stack cannot
use one, and PubApps VMs have none (this one is a two-vCPU VMware guest). Deflate
compression (A1) costs tens of milliseconds of CPU per message. Where a GPU does matter is
the browser: Avivator already renders the OME-TIFF tiles with WebGL, and in Plan B the
170 k-cell tissue scatter would move from a server-drawn PNG to a WebGL layer
(`deck.gl` or plotly `scattergl`) on the reviewer's own GPU, which removes the PNG
render from the server (16 % of its time) and makes recolouring instant; the points go
over once as binary arrays of roughly the same size as one PNG. GPUs on HiPerGator are
useful upstream, in the cell-level pipeline (segmentation, typing, embeddings over 2.6
million cells), which runs offline and does not touch the app's responsiveness; the
trajectory and PAGA computations on 5,000 objects take seconds on a CPU and belong in the
release bundle (B2), not on an accelerator.

### Verification for every item

Run `bench-app.mjs` on staging before and after each change, from campus and with the
`home` profile, and compare against the baseline tables; a change is accepted when its
predicted step improves and no other step regresses. For A1 also confirm the negotiated
extension in the handshake; for A2 the `valueTimeline`; for B1 the browser suites and the
download PNGs; for B3 cache correctness by changing every input in the key and watching
the output change. Record each accepted result in the baseline document with the release
name, and keep the two simultaneous-session runs as the concurrency check.
