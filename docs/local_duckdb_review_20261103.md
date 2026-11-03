# Local DuckDB and display review — 2026-11-03

The current app at [the local preview](http://127.0.0.1:3845/islet-explorer/)
is source `26b66131668c01326d4b95bdb1efbf1d37ef7527`, image
`sha256:2931955ab95330c7804bfea7f9ca25774aa4813f7c9db645a3fd63e35f99035e`.
It includes the follow-on Trajectory and Spatial changes. The strict image
gate, read-only mounts and local runtime are verified; the browser and visual
acceptance receipts below are for the earlier `3c81075` review baseline and are
not relabeled as fresh checks of this successor. The current UI is packaged in
the fully verified [PubApps build candidate](pubapps_ui_refresh_20261103.md),
which uses accepted production data and has not been transferred or activated.

| Component | Earlier native review baseline |
|---|---|
| Application source | `3c81075307c97356cf6ded5df0ac10e810c8b9fa` |
| Immutable local image | `sha256:eb6f3c32779a22e4c1ca7a843bd9de5011de253ac42fb3fd4e574d7b34b2c920` |
| Accepted data producer | `5e6580fe68af12dd18bee430988bc1ed64061526` |
| Read-only data manifest | `03038ebdcd172e82ff7afba1f2025cb905e5859bcee423b05618b04c993f7320` |
| Persistent user service | `islet-explorer-local-review-20261103-r2.service` |

The [successor image validation](../tmp/unclassified-ui-accessibility-20261103/image-3c81075307c97356cf6ded5df0ac10e810c8b9fa/validation-status.json)
passes with **0 High / 0 Critical**, no DuckDB deferral, and the patched decoder
result `0`, error `-17`. The DuckDB binary and installed font receipt match the
previous accepted `46b81d66…` image. Prior Parquet and font-renderer checks remain
inherited evidence; they were not rerun. The intermediate `d6405d6…` /
`aceb89dc…` build was never activated.

The installed successor passed runtime, TIFF Range, read-only data and exclusive
loopback checks. Its [standard native browser run](../tmp/unclassified-ui-accessibility-20261103/runtime-3c81075307c97356cf6ded5df0ac10e810c8b9fa/local-browser.json)
passed **49/49**, with zero Spatial accessibility violations. The
[targeted successor smoke check](../tmp/unclassified-native-20261103/native-3c81075-smoke-receipt.json)
passed: the parent and four subclass choices, 4,842 plotted objects, an exact
readable five-row census, no legacy visible category names and no page errors.
Final Composition checks passed **12/12**, with three images reviewed. Final
Spatial checks passed **12/12**, covering default/dense layouts at 1200/1600px,
brush/reset, donor scope and Distance controls. Direct visual review passed
for four complete raw maps, two default browser views and four Distance views.
Two dense browser crops were excluded because sticky controls overlapped the
captures; this was capture occlusion, not application clipping. The
[card-controls receipt](evidence/card_controls_validation_20261103.json)
records these scopes and limits; the container is healthy.

The [predecessor native family receipt](../tmp/unclassified-native-20261103/native-ccda497-acceptance-receipt.json)
passed on source `ccda49712c7df0e28634d843155c98977c961eae` and image
`sha256:c7608905e097dc4573cf6566b9848bb9fc8ff827d4da5eff3e2657c29ec2431d`.
It covers exported Plot and Trajectory family arithmetic, four help tables,
eight inspected images and family selection controls. Its separate standard
browser run passed 48 of 49 checks; the sole failure was a missing accessible
label on Spatial's distance-group selector. The successor adds an off-screen
“Group by” label and changes two inline flex styles in the donor-composition
toolbar to remove an observed blank row. Classification and count code are
unchanged. The [exact source bridge](../tmp/unclassified-native-20261103/ccda497-to-3c81075-source-bridge.json)
records those three substitutions across two files; 102 other tracked app files
are identical. The predecessor proof stays bound to its original image. The
15 PAGA and 20 Plot-title checks also remain native predecessor evidence linked
by this bridge; they are not rerun claims for the successor. Two ancillary PAGA
segmentation captures with sticky-header occlusion are excluded from visual-pass
claims; their interaction checks passed separately.

Data and View controls start collapsed. The new **Unclassified (all)** parent
contains **823,947 cells** in four mutually exclusive subclasses. The **7,398**
cells with the literal source label `Unclassified` form only the reassessment
subclass. The [hierarchy and disease-status totals](unclassified_display_20261103.md)
explain each subclass and distinguish global tissue counts from filtered card
values and pooled object memberships. The accepted cohort still contains
**20 donors, 3,199,606 cells and 78,680 Ki67-positive cells**. This display change
does not retype source cells or rewrite scientific artifacts.

The [production-data integration receipt](../tmp/unclassified-local-review-20261103/actual-data-run3-production-loader/actual-data-integration.json)
records **131 passing checks on source `337513e…`**, including the real packaged
production loader, raw-to-display counts, both distance definitions and
core/peri metrics. Biological **All T cells** retains CD4, CD8 and raw unresolved
T cells; **All immune cells** remains 279,516 and excludes the 44,684
CD45-negative residual cells. The Barlow graph-input matrix is unchanged, with
no new graph fit. These source/data checks retain their original tested source
and are separate from the native successor checks above.

The immutable [earlier data binding](evidence/local_review_data_binding_20261103.json)
records 18 evidence hashes and 4,987 unchanged file identities for its original
consumer source. The new candidate mounts that same accepted `5e6580f…` data
read-only. Original scientific parity retains its original producer, verifier
and composer; no fresh full parity run or same-source deployable data release
is claimed.

The service uses a persistent unit in `~/.config/systemd/user/`, enabled for the
user's default target, with lingering enabled. This corrected the earlier
transient-unit outage exposed by the 16:04 UTC reboot. Boot persistence is
configured; no additional reboot was performed for testing.

The retained companion timer checks the stable container name every 30 seconds
through the working user bus. It follows the app service's lifecycle and does
not restart the app on an unhealthy result. Podman's automatic health interval
is disabled in the new launcher; the custom timer owns scheduling. The user
manager's private socket remains an unrepaired limitation. No manager
reexecution was performed. Two [new scheduled checks](../tmp/unclassified-ui-accessibility-20261103/runtime-3c81075307c97356cf6ded5df0ac10e810c8b9fa/periodic-health-observations/periodic-health.json)
passed on the successor with container ID, PID, start time, service and mounts
unchanged. No manual health check or restart was used for that observation.

```sh
# Status
systemctl --user status islet-explorer-local-review-20261103-r2.service

# Stop for this session
systemctl --user stop islet-explorer-local-review-20261103-r2.service

# Start again
systemctl --user start islet-explorer-local-review-20261103-r2.service

# Retire the review and disable automatic startup
systemctl --user disable --now islet-explorer-local-review-20261103-r2.service
```

Preserve the active backing content while the preview is used or may be relaunched: accepted data under
`tmp/duckdb-dev-ui-release-20261103/final-5e6580f/artifacts/release-5e6580f/data`,
launcher/runtime under
`tmp/unclassified-ui-accessibility-20261103/runtime-3c81075307c97356cf6ded5df0ac10e810c8b9fa/`,
image evidence under
`tmp/unclassified-ui-accessibility-20261103/image-3c81075307c97356cf6ded5df0ac10e810c8b9fa/`,
and immutable image `eb6f3c32…`.

Historical receipts remain unchanged: the
[1596/846 controls review](evidence/local_controls_review_20261103.json) records
its 49 browser checks and seven-panel/four-screen acceptance; the
[7ac/46b review](evidence/local_duckdb_review_20261103.json) records its earlier
49 checks, four UI scopes and eight tissue-image inspections. The
[earlier category query](evidence/category_counts_by_disease_status_20261103.json)
uses the original raw names. None of those historical checks is relabeled as
fresh acceptance of this candidate.

Port **3845 is unproxied loopback**; nginx was not changed. The workstation's
cgroup-v1 delegation does not enforce the requested 4 GiB/2048-task limits.
External authentication, off-campus access and exhaustive accessibility were
not retested. Public dev remains `e704ac2a…`; production remains `89b3f840…`.
No PubApps host calls, public replacements or security-ticket reply were made
for this local update.
