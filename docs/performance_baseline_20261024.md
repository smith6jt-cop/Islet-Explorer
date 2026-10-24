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
serialisation in R are.

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
     byte ranges are honoured there; the local Shiny Server static route answers a
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

## Re-running

`deploy/README.md` § 5 "Responsiveness baseline" has the commands. Keep the same scenario
and compare medians against the tables above; record new baselines in this file with the
release name.
