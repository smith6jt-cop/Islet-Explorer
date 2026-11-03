# Local DuckDB and display review — 2026-11-03

The owner's latest direction is **local review only**. Open
[the local preview](http://127.0.0.1:3845/islet-explorer/) on this workstation.
The final runtime and check results are recorded in the
[review receipt](evidence/local_duckdb_review_20261103.json).

The candidate contains the reviewed DuckDB source fix, accepted component-v2 data, corrected
Spatial guide placement, complete superscript marker signs, and room for the
rightmost tissue-axis tick. Its two producer identities remain explicit:

| Component | Identity |
|---|---|
| Application source | `7ac421ed7d4b5ab5425d0130e8ff7445f4e87240` |
| Immutable local image | `sha256:46b81d665cacb0f8abcd4b9fe34d7f8c786cb7916d6f9cd30e63a1223358b006` |
| Accepted data producer | `5e6580fe68af12dd18bee430988bc1ed64061526` |
| Read-only data manifest | `03038ebdcd172e82ff7afba1f2025cb905e5859bcee423b05618b04c993f7320` |
| Local service | `islet-explorer-local-review-20261103-r2.service` |

The image passed a fresh strict scan with **0 High / 0 Critical**, no DuckDB
deferral, real Parquet BROTLI round trips and installed-runtime checks. The
running decoder rejects the regression input with result `0`, error `-17`, and
verified source provenance. The native build reused 23 completed dependency
layers; it did not recompile DuckDB. Its font registration and renderer checks
ran during this build and cover four faces, ragg/Cairo PNG, 12 distinct macrophage
profiles and Cairo PDF Unicode mappings.

The running dataset contains **20 donors, 3,199,606 cells, 78,680 Ki67-positive
cells and 7,398 Unclassified cells (0.2312%)**. Every phenotype count matches the
accepted census. The [data binding](evidence/local_review_data_binding_20261103.json)
records 18 rechecked evidence hashes and 4,987 unchanged file identities.
Scientific acceptance compares actual artifacts: canonical and nongeometry
content matches, and all 20 GPKG schemas/tables/rows match except valid UTC build
timestamps. Original accepted parity retains its original producer, verifier and
composer; no new full parity execution or same-source deployable release is claimed.

All four UI scopes and all eight visual inspections passed. The local checks
include HTTP/health, exclusive loopback publication, read-only
mounts, disabled optional/debug features, exact TIFF byte ranges, 49 browser
checks, and the final UI/visual results in the receipt. The UI checks cover
Plot/Trajectory layout, taxonomy, paired Spatial modes at 1200/1600px, and tissue
default/all-27-category views with four macrophage markers. Eight actual native
raw/card images are reviewed for complete marker signs, legends and axis ticks;
real mouse brushing must change the domain and Reset must restore it exactly.

The first local preview used port 3840. Final routing review found an existing
nginx listener on port 8080 forwarding to 3840, so that preview was stopped and
the corrected candidate moved to **unproxied loopback port 3845**. Nginx was not
changed. Initial failed rendering reports and screenshots remain preserved under
`tmp/local-app-review-20261103/`; the final receipt binds them and the separate
brush diagnostic. The corrected test records the real brush event before Shiny
clears its transient input, then requires a stable changed domain and exact reset.

Preserve the active backing files while this preview is in use or may be
relaunched. The data root is
`tmp/duckdb-dev-ui-release-20261103/final-5e6580f/artifacts/release-5e6580f/data`;
the launcher, runtime files and evidence are under
`tmp/duckdb-local-review-20261103/r2/`, with image evidence under `r2-image/`.
These paths and immutable image `46b81d66…` are active preview content, not generic
cleanup candidates.

The service persists after this chat or launcher shell exits, with restart on
failure. It is a transient user service and is not enabled across reboot. Use:

```sh
# Status
systemctl --user status islet-explorer-local-review-20261103-r2.service

# Stop
systemctl --user stop islet-explorer-local-review-20261103-r2.service

# Relaunch the exact local candidate
bash /data/panc_CODEX/Islet-Explorer/tmp/duckdb-local-review-20261103/r2/start-local-service.sh
```

The independent local audit records a workstation limitation: Podman/systemd
request 4 GiB and 2048 tasks, but this cgroup-v1 user session does not enforce
those per-service limits. This is not a claim about the separately verified
PubApps limits. External authentication, public proxy acceptance, off-campus
access, reboot persistence and exhaustive accessibility were not tested here.

Neither public PubApp was replaced. Public dev remains the previously deployed
patched image `e704ac2a…`; production remains `89b3f840…`, with its separate
DuckDB deadline of 2026-11-10. Before the local-only direction, versioned tooling
had been staged and two proven unused archived image caches removed; both live
app identities were verified unchanged. No new candidate image import, data
relocation or service activation occurred on PubApps, and no host calls followed
the local-only direction. Historical deployment receipts are unchanged. No
security-ticket reply was sent.
