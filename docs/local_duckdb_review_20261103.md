# Local DuckDB and display review — 2026-11-03

Open [the local preview](http://127.0.0.1:3845/islet-explorer/) on this workstation.
**Data controls and View controls now start collapsed.** Click a heading to
expand it. All seven panels across Plot, Trajectory, Spatial and Viewer retain
working controls; their initial plots and images load while the panels are closed.
The current [review receipt](evidence/local_controls_review_20261103.json)
records the actual candidate and checks. This remains **local review only**.

| Component | Current identity |
|---|---|
| Application source | `1596b0336096a949a954c9793fa6600ed09c7f1a` |
| Immutable local image | `sha256:846962f791ed226fe22c50f8590eef6330fa76726fa616f78e4d9b106e6a0ddc` |
| Accepted data producer | `5e6580fe68af12dd18bee430988bc1ed64061526` |
| Read-only data manifest | `03038ebdcd172e82ff7afba1f2025cb905e5859bcee423b05618b04c993f7320` |
| Persistent user service | `islet-explorer-local-review-20261103-r2.service` |

The fresh image scan passes with **0 High / 0 Critical**, no DuckDB deferral,
and the patched decoder result `0`, error `-17`. The DuckDB binary and installed
font evidence are identical to the preceding accepted `46b81d66…` image. Its
Parquet and font-renderer validation is retained as prior evidence, not claimed
as newly rerun. The code change only removes the four Data-panel open defaults
and updates two comments; View panels already defaulted closed. No scientific
computation or data changed.

Fresh native acceptance passed all **49 browser checks**, plus a focused check
of seven initially closed panels, seven physical open/close toggles, meaningful
rendered content on all four tabs, and four inspected post-render screenshots.
Excluding T1D changes the Plot count from 4,842 to 3,202 objects; restoring it
returns 4,842. Those are plot-object counts, not the total cell census.

The accepted component-v2 data contains **20 donors, 3,199,606 cells, 78,680
Ki67-positive cells and 7,398 Unclassified cells (0.2312%)**. Its earlier
[data binding](evidence/local_review_data_binding_20261103.json) records 18
verified evidence hashes and 4,987 unchanged file identities for the previous
consumer source. Current runtime bindings separately identify `1596b03` and the
same read-only data root. Original scientific parity retains its original
producer, verifier and composer; no new full parity execution or same-source
deployable release is claimed.

The [category counts by disease status](evidence/category_counts_by_disease_status_20261103.json)
were freshly queried from the running candidate. Across all tissue cells, ND /
Aab+ / T1D respectively contain 311 / 590 / 418 Doublet suspected cells,
1,701 / 3,190 / 2,507 Unclassified cells, and 30,061 / 32,679 / 36,039 cells
without a subtype. In this dataset, the last category combines general
Hematopoietic CD45+, Hematopoietic CD45- and T cell labels, excluding conflicts.
These whole-dataset totals do not apply the composition chart's object, region
or donor filters.

The workstation reboot at 16:04 UTC exposed a lifecycle problem: the earlier
review service existed only under `/run/user/1001/systemd/transient` and was not
recreated after reboot. The replacement is a persistent unit in
`~/.config/systemd/user/`, enabled for the user's default target. User lingering
is enabled. The current service has been started and verified; an additional
workstation reboot was not performed for testing.

A persistent companion timer checks the stable container name every 30 seconds;
two scheduled checks passed with the app's container ID, PID and start time
unchanged. The timer follows the app service's lifecycle and does not restart
the app on an unhealthy result. The user manager's private socket currently
refuses connections, preventing Podman's automatic timer creation. The companion
timer uses the working user bus. The private socket failure was observed after
unit verification against the live runtime directory; do not repeat that
verification command there. No manager reexecution was attempted, since its
ability to repair this socket pathname without disruption was not established.

```sh
# Status
systemctl --user status islet-explorer-local-review-20261103-r2.service

# Stop for this session
systemctl --user stop islet-explorer-local-review-20261103-r2.service

# Start again
systemctl --user start islet-explorer-local-review-20261103-r2.service

# Stop and turn off automatic startup when this review is retired
systemctl --user disable --now islet-explorer-local-review-20261103-r2.service
```

Preserve the active backing content while the preview is used or may be
relaunched: accepted data under
`tmp/duckdb-dev-ui-release-20261103/final-5e6580f/artifacts/release-5e6580f/data`,
the current launcher/runtime under
`tmp/duckdb-local-controls-20261103/runtime-1596b03/`, image evidence under
`tmp/duckdb-local-controls-20261103/image-1596b03/`, and immutable image
`846962f7…`. These are not generic cleanup candidates.

The initial preview on port 3840 was stopped after discovering an existing
nginx port-8080 proxy to it. Port **3845 is unproxied loopback**; nginx was not
changed. The earlier [7ac/46b review receipt](evidence/local_duckdb_review_20261103.json)
preserves its original 49 browser checks, four UI scopes and eight tissue-image
inspections, as well as the earlier rendering failures. The new focused controls
check does not claim those broader visual checks were repeated.

The workstation's cgroup-v1 delegation still does not enforce the requested
per-service 4 GiB/2048-task limits. This is not a claim about PubApps limits.
External authentication, off-campus access and exhaustive accessibility were not
retested for this local update. Public dev remains `e704ac2a…`; production remains
`89b3f840…`. Neither public PubApp was replaced, no PubApps host calls were made,
and no security-ticket reply was sent.
