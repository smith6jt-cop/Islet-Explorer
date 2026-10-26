# Responsiveness baseline for the reviewer deployment (2026-10-24)

**Question.** Off campus, the public app at `https://islet-explorer.rc.ufl.edu/islet-explorer/`
"works" but feels slower than expected to respond to control changes and to move between tabs.
Is that the app's own architecture or the PubApps hosting (VM, NFS data, UFRC proxy)? And what
is a reasonable expectation for loading and responsiveness?

**Answer in one paragraph.** For one reviewer, the hosting adds nothing measurable: through
the public URL from a campus workstation, and through an SSH tunnel straight to the VM
backend, every step of a fixed scenario took the same time as the identical release
container running on the build workstation (local NVMe data, a 32-core CPU against the
VM's two). What an off-campus reviewer feels is the app's own payload crossing their
link: the first screen sends 3.5 MB over the WebSocket (the Plot tab is rendered twice at
startup), every change of the plotted measure sends 1.7 MB (58 % of it per-point hover text),
the Spatial tab 1.5 MB, the PAGA layout 1.5 MB, and the Viewer's first image view 10.7 MB of
tile ranges. With an emulated residential link (40 ms, 25/5 Mbps) the payload-heavy steps
take about twice as long as on campus; on a poor link (120 ms, 8/2 Mbps) about three times.
Steps with small payloads (drill-down, opening the Trajectory tab) are unchanged. Two
one-time costs per R process also stand out: the first opening of the Trajectory tab
(about 12 s, the scanpy diffusion pseudotime) and the first PAGA layout (about 5 s); the
production unit's health check keeps the R process alive, so the first reviewer to open
those after each activation pays them once, and everyone after them does not.

## Method

- **Release under test:** `pubapps-2026-10-uniform-v10` (commit `d55680a`, image
  `ae95fd23…`, data schema 5). The same image and the same data bundle everywhere.
- **Vantage points.** (A) *Workstation container:* the release image started with the
  unit's mounts and environment on the build workstation (AMD Threadripper PRO 7975WX,
  data on local disk), reached on the loopback. (B) *Public URL from campus:* production
  on the VM through UFRC's proxy, TLS and Basic Auth, from workstation 10.15.152.7 with
  the reviewer account. (C) *Public URL with an emulated link:* as B, with Chrome DevTools
  network emulation added on top of the real path: `home` 40 ms latency, 25 Mbps down,
  5 Mbps up; `slow` 120 ms, 8 Mbps, 2 Mbps. (D) *Staging instance:* the first session in
  its R process. (E) *Cold R process:* the workstation container right after a restart.
  (F) *VM backend through the SSH tunnel* (workstation → HPG login node → `10.13.41.21:3838`),
  the vantage point that separates the proxy from the VM.
- **Scenario** (`deploy/browser/bench-app.mjs`, three cold runs per vantage point in a
  fresh browser context, medians reported, then two warm visits): page load until the
  Plot tab's two cards are drawn; two changes of the plotted measure; a drill-down click
  (case geometry, the object's cells, the segmentation map); the Trajectory tab and its
  PAGA layout; the Spatial tab and a donor change; the Viewer until its first tiles and
  until the Range traffic goes quiet. For each step the tool records the time until the
  server reported busy (one round trip), the time until its last output message arrived
  (server compute plus transfer), the bytes that crossed the WebSocket, and the browser's
  own draw after the last message.
- **Supporting measurements:** HTTP timings of the page, scripts and TIFF Range reads with
  curl (proxy versus workstation container); a capture of the WebSocket frames split by
  output and by plotly trace field; a server-side micro-benchmark run inside the container
  (CPU, DuckDB over the release Parquet, GeoPackage window read, cells CSVs), on the
  workstation and in the production container on the VM.
- **Evidence:** `tmp/ufrc-release-20261022/evidence/performance-baseline-20261024/`
  (bench reports, curl timings, payload analysis, micro-benchmark). Reports hold timings
  and byte counts only.

## Results

Medians of three cold runs. Steps are wall time from the action until the output is drawn.

| Step | Workstation container | VM backend via tunnel | Public URL, campus | Public URL, 25/5 Mbps + 40 ms | Public URL, 8/2 Mbps + 120 ms |
|------|----------------------:|----------------------:|-------------------:|------------------------------:|------------------------------:|
| Page: first byte | 122 ms | 118 ms | 131 ms | 167 ms | 245 ms |
| Page: Shiny connected | 268 ms | 272 ms | 324 ms | 663 ms | 1.5 s |
| Page: Plot tab drawn | 3.9 s | 3.5 s | 3.6 s | 8.8 s | 11.5 s |
| Page assets transferred (gzip) | 1.5 MB | 1.5 MB | 1.5 MB | 1.5 MB | 1.5 MB |
| WebSocket bytes to the first screen | 3.5 MB | 3.5 MB | 3.5 MB | 3.5 MB | 3.5 MB |
| Change the plotted measure | 1.8 s | 1.1 s | 1.2 s | 2.9 s | 3.7 s |
| Change it again | 1.3 s | 1.3 s | 1.3 s | 3.1 s | 3.8 s |
| Drill-down click | 513 ms | 518 ms | 753 ms | 667 ms | 753 ms |
| Open Trajectory tab | 202 ms | 182 ms | 237 ms | 214 ms | 228 ms |
| Build PAGA layout | 1.0 s | 1.1 s | 1.0 s | 2.4 s | 3.1 s |
| Open Spatial tab | 2.2 s | 2.4 s | 2.2 s | 3.8 s | 4.5 s |
| Change Spatial donor | 1.2 s | 1.3 s | 1.9 s | 1.6 s | 1.9 s |
| Viewer: first tiles | 533 ms | 579 ms | 612 ms | 1.5 s | 3.8 s |
| Viewer: Range bytes until quiet | 2.6 MB | 6.5 MB | 10.7 MB | 10.7 MB | 10.7 MB |
| Warm revisit: Plot tab drawn | 3.9 s | 3.9 s | 3.5 s | 8.1 s | 10.7 s |

Server compute plus transfer ("busy" until the last output message) and WebSocket bytes per
step, workstation container versus the public URL from campus:

| Step | Workstation | Public URL, campus | WebSocket bytes |
|------|------------:|-------------------:|----------------:|
| Change the plotted measure | 1.7 s / 1.2 s | 1.1 s / 1.2 s | 1.7 MB |
| Drill-down click | 0.49 s | 0.62 s | 0.2 MB |
| Open Trajectory tab | 0.15 s | 0.18 s | 31 KB |
| Build PAGA layout | 0.83 s | 0.86 s | 1.5 MB |
| Open Spatial tab | 2.1 s | 2.1 s | 1.5 MB |
| Change Spatial donor | 1.1 s | 1.7 s | 0.6 MB |

One-time costs per R process (first session after the process starts):

| Step | Workstation container, cold process | Staging instance on the VM, first session |
|------|------------------------------------:|------------------------------------------:|
| Page: first byte (R start and data load) | 7.1 s | 146 ms (process already running: the unit's health check requests the page every 30 s) |
| Page: Plot tab drawn | 12.8 s | 5.7 s |
| First opening of the Trajectory tab | 11.6 s | 12.7 s |
| First PAGA layout | 5.9 s | 5.2 s |
| Later sessions in the same process | 0.2 s / 1.0 s | 0.2 s / 1.0 s |

HTTP (curl, medians of five) through the proxy from campus, through the tunnel to the VM
backend, and on the workstation container: page first byte 132 / 107 / 136 ms (the R
process renders the page), the 1.09 MB gzipped plotly bundle 98 / 80 / 78 ms, an 8 MiB
mid-file TIFF Range 169 / 99 / 9 ms (about 50 MB/s through the proxy, 85 MB/s through the
tunnel), a 16-byte Range 21 ms through the proxy with a 14 ms TLS handshake. The proxy adds
about 20 ms per request over the backend.

WebSocket payload of the Plot tab (frame capture on the workstation container): the first
screen carries 3.56 MB, of which the distribution card 1.82 MB and the size scatter 1.48 MB,
each sent twice (both arrive at about 3.3 s and again at about 5.1 s after connecting: the
Population Measure and Metric selectors are `renderUI` outputs whose inputs mount after the
first render and invalidate it). One measure change carries 1.70 MB: per-point hover
`text` 0.99 MB (58 %), `x` and `y` 0.46 MB (full double precision), `key` 0.16 MB.
Compression is not available on this socket (Shiny Server's SockJS does not negotiate
permessage-deflate and the proxy cannot compress frames), so payload reduction has to
happen in the app.

The VM (`pubmaigan1`, VMware guest): 2 vCPUs of an AMD EPYC 9355, 15.7 GB of memory,
load average 0.00 with both instances up 18 hours and healthy, their R processes resident
at 2.3 GB (production) and 2.1 GB (staging) against caps of 9 and 4 GiB. `/pubapps` is
NFS v3 over TCP from the Qumulo filer (1 MiB read size); an 8 MiB read from a TIFF ran at
750–850 MB/s and the 38 MB cell-distance Parquet came from the page cache.

Server-side micro-benchmark, single thread, workstation container versus the production
container on the VM:

| Measurement | Workstation | VM |
|---|---:|---:|
| CPU: scalar loop (5 M) | 90 ms | 91 ms |
| CPU: 1000×1000 matrix product | 31 ms | 20 ms |
| CPU: sort 5 M values | 379 ms | 328 ms |
| DuckDB: summary `islets.parquet` | 15 ms | 11–13 ms |
| DuckDB: one donor's tissue table (174 k rows) | 14–20 ms | 14–16 ms |
| DuckDB: 38 MB cell-distance Parquet, grouped | 15 ms | 17 ms |
| sf: GeoPackage window read, 2,172 features | 17 ms | 201–259 ms |
| 20 cells CSVs (first touch / cached) | 15 / 14 ms | 46 / 21 ms |

The VM's core is as fast as the workstation's for R, which is why the two vantage points
match. The one data access that is slower over NFS is the GeoPackage window read (SQLite
page reads, about 0.2 s), visible as the drill-down's extra 0.2–0.3 s through the proxy.
Data access is otherwise not where the seconds go; the plotly figure construction and
serialization in R are.

## Reading

1. **Hosting versus app.** The VM backend through the tunnel and the public URL from
   campus both matched the workstation container step for step (differences of 0 to
   0.3 s, within run-to-run noise), although the workstation has sixteen times the cores
   and local data: the VM's two cores are as fast per thread, the app is single-threaded
   per session, and the data reads are small. The VM, its NFS data and the proxy are
   therefore not what the reviewer feels. The WebSocket is a real `wss://` socket
   through the proxy (no SockJS fallback), redirects and assets are served correctly, and
   the proxy's per-request overhead is about 20 ms.
2. **The off-campus experience is payload times link.** Emulating a residential link on
   the same path doubled every payload-heavy step (first screen 3.6 → 8.8 s, measure
   change 1.2 → 2.9 s, Spatial 2.2 → 3.8 s, PAGA 1.0 → 2.4 s, Viewer 0.6 → 1.5 s) and left
   the light steps alone (drill-down, Trajectory tab). A poorer link stretches them
   further (first screen 11.5 s, measure change 3.7 s, Viewer 3.8 s). The reviewer's own
   link speed therefore decides between the three right-hand columns; the hosting cannot
   change that, the app's payload can.
3. **App-side costs that would improve every link.**
   - The Plot tab's startup double render (1.7 MB and about 1.5 s per page load, about
     3 s on a home link).
   - Per-point hover text in both Plot cards (about 1 MB per interaction). Building the
     hover from `customdata` with a `hovertemplate`, or shortening it, and rounding `x`
     and `y` to the displayed precision would cut the 1.7 MB per interaction by roughly
     two thirds.
   - The first Trajectory open (12 s) and first PAGA layout (5 s) per R process. These
     could be computed when the process starts, or warmed by a request after activation,
     so no reviewer meets them.
   - The Viewer's first view pulls 10.7 MB of tile ranges through Apache (multipart
     byte ranges are honored there; the local Shiny Server static route answers a
     multi-range request with a 200 and the whole file, which only matters for local runs).
4. **Shared R process.** Shiny Server (open source) serves every session of the app from
   one R process, so simultaneous reviewers queue behind each other's compute (about
   1.1 s per Plot change, 2.1 s per Spatial open). See the concurrency note below.

## What to expect (the baseline)

For one reviewer with a fast link (campus or comparable), the reference values are the
"Public URL, campus" column: first screen in about 4 s, a Plot change in about 1.2 s,
Spatial in about 2.2 s, drill-down and Viewer under a second, Trajectory instant once the
process has computed it. On a typical home connection expect about twice those times for
the payload-heavy steps. Values well beyond the matching column point at something new:
another reviewer's concurrent compute, a changed release, or a hosting problem, which the
`vm-measure.sh` and the tunnel run would then separate.

## Concurrency note

Two simultaneous sessions against the workstation container (two benchmark processes
started together, two runs each): see the evidence directory (`local-concurrent-a.json`,
`local-concurrent-b.json`). With just two sessions active at once, each step took roughly
twice as long, and the light steps queued behind the other session's heavy ones: the
drill-down click went from 0.5 s to 4.6–5.9 s, opening the Trajectory tab from 0.2 s to
2.2–2.6 s, a Plot change from 1.3 s to 2.1–2.7 s, the Spatial tab from 2.2 s to 3.1–4.0 s.
With the five concurrent reviewers the instance request allows for, a reviewer who clicks
while others are mid-computation can wait several seconds for a step that is instant alone.
This is a property of the single R process, not of the VM or the proxy; the remedies are
the payload and startup work above (less compute per step) or more than one R process per
app, which open-source Shiny Server does not offer and would need a different serving
arrangement.

## Plan A Level 1, measured locally before deployment (2026-10-24, commit `eba12b4`)

The same scenario against the edited app on the workstation (plain `shiny::runApp` on the
development root, so the Viewer's Range step and the Parquet backend are not represented;
the Shiny Server compression is in the image build and is measured on staging):

| Measure | Before | After Level 1 |
|---|---:|---:|
| WebSocket bytes to the first screen | 3.56 MB (Plot cards sent twice) | 1.20 MB (sent once) |
| WebSocket bytes per Plot change | 1.70 MB | 1.12 MB |
| of which hover text | 0.99 MB | 0.45 MB |
| First Trajectory open in a fresh R process | 11.6–12.7 s | 0.25 s (warmed at start) |
| First PAGA layout in a fresh R process | 5.2–5.9 s | 1.3 s |
| Plot change, server compute | 1.1–1.2 s | 1.1–1.2 s (unchanged; Level 2's `plot_ly()` item) |

What remains in a Plot change is the two cards' 8,600 object points: compact hover text
0.45 MB, keys 0.16 MB, coordinates 0.18 MB, and a per-point color array on the
distribution card's single marker trace (0.09 MB). With the socket compression the
1.12 MB travels as roughly 0.2 MB. The test suite (370) and the deployment test suites
pass; the effect through the proxy is in the next section.

## Plan A Level 1 on staging (2026-10-25, commit `e10b2ee`, release `pubapps-2026-10-uniform-v10-perf1`)

The release built from the Level 1 commits (the image with the compression patch and the
body-parser pin, the data bundle of the same commit) was activated on the staging instance
`islet-explorer-dev` on 2026-10-25 16:03 UTC with `--local-data /local/maigan/releases`;
production stays on `d55680a`. The same scenario through UFRC's proxy at
`devislet-explorer.rc.ufl.edu`, five runs per column, medians. The baseline columns are
the 2026-10-24 runs: production from campus, the staging instance cold from campus (first
open of each tab in a fresh R process) and production with the emulated home link.
`socket.extensions` reported `permessage-deflate` on every connection. The browser counts
WebSocket bytes after decompression, so the bytes columns show the payload trims only; the
compression shows in the times on the emulated link.

| Step | Production, campus (baseline) | Staging cold, campus (baseline) | Staging Level 1, campus | Production, home link (baseline) | Staging Level 1, home link |
|---|---:|---:|---:|---:|---:|
| Page: Plot tab drawn | 3.6 s | 5.7 s | 3.6 s | 8.8 s | 4.1 s |
| WebSocket bytes to the first screen | 3.5 MB | 3.5 MB | 1.1 MB | 3.5 MB | 1.1 MB |
| Plot change (measure) | 1.2 s | 1.2 s | 1.5 s | 2.9 s | 1.5 s |
| Plot change again | 1.3 s | 1.3 s | 1.2 s | 3.1 s | 1.6 s |
| Drill-down click | 753 ms | 853 ms | 704 ms | 667 ms | 603 ms |
| Trajectory tab | 237 ms | 12.7 s | 264 ms | 214 ms | 192 ms |
| PAGA layout | 1.0 s | 5.2 s | 1.2 s | 2.4 s | 1.6 s |
| Spatial tab | 2.2 s | 3.1 s | 3.2 s | 3.8 s | 2.6 s |
| Spatial donor change | 1.9 s | 2.0 s | 1.9 s | 1.6 s | 1.2 s |
| Viewer first tiles | 612 ms | 658 ms | 615 ms | 1.5 s | 1.5 s |
| Warm revisit: Plot tab drawn | 3.5 s | 5.8 s | 2.5 s | 8.1 s | 3.5 s |

Reading:

- On the emulated home link every payload-heavy step roughly halves: first screen 8.8 →
  4.1 s, Plot change 2.9 → 1.5 s, Spatial 3.8 → 2.6 s, PAGA 2.4 → 1.6 s, warm revisit
  8.1 → 3.5 s. Compression and trims together bring the home-link times down to the
  baseline's campus times, which is what Level 1 set out to do.
- On campus the step times are unchanged within noise, as expected (the link was never
  the cost there), except the warm revisit (3.5 → 2.5 s, one render instead of two) and
  the first open of Trajectory and PAGA in a fresh process (12.7 s and 5.2 s → 0.26 s and
  1.2 s), which the warm-up now pays before the health check lets traffic in.
- The Plot change still costs 1.1–1.5 s of server compute. That is Level 2's direct
  `plot_ly()` item.
- The Spatial tab opens in 3.2 s on staging against 2.2 s on production from campus; the
  staging instance measured 3.1 s on the baseline day too, so this is the instance (4 GiB
  cap), not Level 1. Worth a look when Level 2 touches the tissue scatter.
- The Viewer is unchanged: its tiles come from Apache, outside the socket.

Evidence: `tmp/ufrc-release-20261022/evidence/performance-baseline-20261024/level1-staging/`
(`bench-staging-level1-campus.json`, `bench-staging-level1-home.json`, logs).

## Plan A Level 2, measured locally before deployment (2026-10-26)

The in-house service on port 3840 (development backend, no warm-up, no socket compression,
the workstation's `ragg` already in use) through nginx on port 8080, five cold runs of the
same scenario with `bench-app.mjs`, the R process restarted just before run 1 so that run 1
starts with an empty cache and the later runs show the cross-session cache. Level 1 local
numbers are the two-run measurement of 2026-10-24 on the same service.

| Step | Level 1 local | Level 2 local, run 1 (empty cache) | Level 2 local, runs 2–5 (cache warm) |
|---|---:|---:|---:|
| Plot tab drawn (first screen) | 4.1 s | 9.3 s (fresh R process, data load included) | 2.7 s |
| WebSocket bytes to first screen | 1.2 MB | 1.2 MB | 1.2 MB |
| Plot change, total | 1.2 s | 1.1 s | 0.75 s |
| Plot change, server compute + transfer | 1.1 s | 0.96 s | 0.65 s |
| Plot change again, server compute + transfer | 1.2 s | 0.77 s | 0.52 s |
| Spatial tab, server compute + transfer | 2.1 s | 4.3 s (first donor CSV of the process) | 1.7 s |
| Spatial donor change, server compute + transfer | 2.2 s | 1.0 s | 0.37 s |
| Trajectory tab, first open | 0.20 s (warm process) | 18 s (no warm-up locally) | 0.25 s |
| PAGA layout | 1.1 s | 5.9 s (first of the process) | 1.1 s |
| Warm revisit, Plot tab drawn | 2.4 s | 2.3 s | 2.2 s |

Reading:

- The Plot cards no longer go through `ggplotly()`: building a card is about 100 ms of R
  plus 120 ms of `plotly_build()` and JSON on the full data, so a cache miss on a measure
  change costs 0.55–0.92 s of server time end to end (data path, both cards, the other outputs
  of the tab) against 1.1 s before, and a cache hit 30–80 ms (`shiny::testServer`, one R
  process, 4,759 objects). In the running app the `plot_card_built` log confirms it: a two-run
  check logged six builds, all in the first five seconds of run 1 (first screen and the two
  measure changes, two cards each), and none in run 2 or the warm revisits. The 0.5–0.65 s
  the bench still counts as "server compute + transfer" on a warm Plot change is therefore
  not the figures: it is the rest of the tab's work on a measure change (the data path, the
  measure selectors' re-render and their round trips) and, mostly, the browser, which parses
  and draws the first 450 KB card before it takes delivery of the second, so the metric's
  floor is a card's draw time, not the server's.
- The first version of the direct figures doubled the payload (2.4 MB first screen, 2.3 MB a
  Plot change): a `plot_ly()` widget keeps its build inputs beside the built traces and sends
  both. Stripped, the messages are back at Level 1's 1.2 MB and 1.1 MB (`ws-payload.txt`
  against `ws-payload-before-strip.txt` in the evidence directory); the hover text is now
  58 % of a card's message, the keys 19 %, the per-point colors 11 %.
- The Spatial donor change is the clearest cache effect: 2.2 s to 0.37 s on a repeat.
- The drill-down (1.3–4.4 s here) is the development GeoJSON backend, not the app change;
  the release bundle reads indexed GeoPackages (0.6–0.7 s on staging).
- Browser suites on the same service: `check-app.mjs` 36 of 36, the UI review suite 18 of 18
  (two layout checks updated for the tab bar inside the header row), the trajectory links
  suite passed; `check-plot-load.mjs` fails on the uniform root's phenotype names (a
  historical 15-donor check, failing before this change too).

Evidence: `tmp/ufrc-release-20261022/evidence/performance-baseline-20261024/level2-local/`
(`bench-local-level2.json`, the pre-strip run, `ws-payload*.txt`, the browser reports, the
header screenshots).

## Plan A Level 2 on staging (2026-10-26, commit `0b55429`, release `pubapps-2026-10-uniform-v10-perf2`)

Same scenario and tool as the Level 1 section, through `devislet-explorer.rc.ufl.edu` with the
staging account, five cold runs from campus and five with the emulated home link, run
within minutes of the activation (the first campus run is the process's first traffic after
its warm-up). Medians of the five runs; the two right-hand columns are the Level 1 and
Level 2 candidates on the same instance.

| Step | Staging cold baseline, campus (2026-10-24) | Level 1, campus | Level 2, campus | Production baseline, home link | Level 1, home link | Level 2, home link |
|---|---:|---:|---:|---:|---:|---:|
| Plot tab drawn (first screen) | 5.7 s | 3.6 s | 2.7 s | 8.8 s | 4.1 s | 3.1 s |
| WebSocket bytes to first screen | 3.5 MB | 1.1 MB | 1.2 MB | 3.5 MB | 1.1 MB | 1.2 MB |
| Plot change | 1.2 s | 1.5 s | 0.72 s | 2.9 s | 1.5 s | 1.1 s |
| Plot change, server compute + transfer | 1.1 s | 1.5 s | 0.64 s | 2.8 s | 1.4 s | 1.0 s |
| Plot change again | 1.3 s | 1.2 s | 0.54 s | 3.1 s | 1.6 s | 0.94 s |
| Drill-down click | 853 ms | 704 ms | 654 ms | 667 ms | 603 ms | 703 ms |
| Trajectory tab (first open per browser session) | 12.7 s | 264 ms | 711 ms | 214 ms | 192 ms | 217 ms |
| PAGA layout | 5.2 s | 1.2 s | 1.2 s | 2.4 s | 1.6 s | 1.2 s |
| Spatial tab | 3.1 s | 3.2 s | 1.3 s | 3.8 s | 2.6 s | 1.6 s |
| Spatial donor change | 2.0 s | 1.9 s | 403 ms | 1.6 s | 1.2 s | 414 ms |
| Viewer first tiles | 658 ms | 615 ms | 593 ms | 1.5 s | 1.5 s | 1.5 s |
| Warm revisit, Plot tab drawn | 5.8 s | 2.5 s | 1.8 s | 8.1 s | 3.5 s | 3.1 s |

Reading:

- The Plot change halves again on campus (1.5 s → 0.72 s total, 0.64 s of it server
  compute plus transfer) and reaches 1.1 s on the home link; its message stays at 1.1 MB
  with `permessage-deflate` negotiated on every connection. The first run of the process
  built its six cards (0.93 s per change) and every later run was served from the cache:
  the staging log holds six `plot_card_built` events for ten cold runs and ten warm
  revisits across both benchmarks.
- The Spatial tab is where the cache shows most: 3.2 s → 1.3 s on campus, 2.6 s → 1.6 s at
  home, and a donor change 1.9 s → 0.40 s (server part 987 → 133 ms). The first Spatial
  open of the process still took 3.3 s, so the tissue scatter's first draw is unchanged by
  `ragg`; it is the repeat draws that are gone.
- First screen 3.6 → 2.7 s on campus (2.3 s once the process has drawn it) and 4.1 → 3.1 s
  at home; warm revisit 2.5 → 1.8 s. The PAGA message shrank from 1.45 MB to 0.97 MB
  (the build inputs no longer travel) without a change to the step time on campus.
- One regression to note: the first Trajectory open of the process took 4.2 s on the
  first campus run (later sessions 0.15 s of server time, 0.6–0.7 s drawn), although the
  warm-up had completed 26 s after the start. The warm-up computes the pseudotime and the
  PAGA layout; what it does not do is the process's first PNG draw, which with `ragg` and
  `systemfonts` in a fresh container builds the font cache once. A one-line addition to
  the warm-up (draw one small PNG at start) would move that cost before the health check;
  it is the next Level 2 follow-up, not a reason to hold the release.
- `check-app.mjs --auth-file` through the proxy: passed, every check (36), no failed
  checks, no browser errors, no accessibility rule violations.

Evidence: `tmp/ufrc-release-20261022/evidence/performance-baseline-20261024/level2-staging/`
(`bench-staging-level2-campus.json`, `bench-staging-level2-home.json`, logs) and
`evidence/browser-staging-0b55429/`.

## What to do about it

The [performance plans](performance_plans_20261024.md) turn these findings into Plan A
(quick wins by feasibility level: socket compression in the project's Shiny Server
build, the startup double render, warming the trajectory computations, payload trims,
direct `plot_ly()` figures, caching) and Plan B (a complete overhaul that moves the
interactive work to the browser), with the profile that shows `ggplotly()` conversion as
the dominant server cost.

## Re-running

`deploy/README.md` § 5 "Responsiveness baseline" has the commands. Keep the same scenario
and compare medians against the tables above; record new baselines in this file with the
release name.
