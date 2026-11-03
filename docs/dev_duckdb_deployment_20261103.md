# Dev DuckDB fix — live 2026-11-03

The [dev PubApp](https://devislet-explorer.rc.ufl.edu/islet-explorer/) is live with
the reviewed DuckDB/Brotli backport. Activation completed at **12:41:40 UTC**;
running-container and public verification completed at **12:42:17 UTC**.
The strict security gate passes **without a DuckDB deferral**. Production remained
healthy on its existing image throughout preparation, activation and verification.

| Item | Verified identity/result |
|---|---|
| Dev release | `pubapps-2026-11-dev-duckdb-v1` |
| Source and data producer | Clean `378843e06c0cd66d8dfff5e3b60664f8a45569ee`, branch `security/dev-duckdb-current-20261103` |
| Immutable image | `sha256:e704ac2a7f01c0fdec5905cc5b1dff43041d45fb0ba7ebfe3c821a13da4a9e9e` |
| Running DuckDB | R package 1.5.2; actual decoder result `0`, expected error `-17`; source provenance passes |
| Installed DuckDB library SHA-256 | `554342ee815993564529021db383f81038538d1b992f3a99e074894c7b9ce4f6` |
| Security | `pass-no-high-or-critical`; 0 High / 0 Critical package findings; no CVE suppression or DuckDB deferral |
| Outer release manifest SHA-256 | `c0ce0656aea0f3e164019381f0a66480c02e054290d97b6584ff0980552b0bb6` |
| Data manifest SHA-256 | `e4b2dbd3bf709513bc6a1d5f8cfdb76ac856efacec25940e13c9975f03936b41` |
| OCI archive | 2,188,275,200 bytes; SHA-256 `392c9b69d33cd124142c90713036a7898cde2faff688893d04e681d071824988` |
| Host service | `islet-explorer-dev`, healthy, 4 GiB limit, private backend `10.13.41.21:3839` |
| Public verification | Dev 15/15 access checks and 49/49 browser checks; production 15/15 access checks |

The [deployment receipt](evidence/dev_duckdb_deployment_20261103.json) binds the
actual running decoder, image scan, data acceptance, host state, public checks,
and before/after production snapshots. The
[data completion receipt](evidence/dev_release_data_20261103.json) is copied
byte-for-byte from the data acceptance run.

## Fix and validation

The reviewed [DuckDB port, PR 27791](https://github.com/duckdb/duckdb/pull/27791)
remained open and unmerged when checked on November 3. This build applies its
immutable commit `2b942af6319783a5cc3c21d76a8e560e191978b2`, patch SHA-256
`d00033af815fcb214daa9630d9da11b39a191dfab17b0f8cff367157af5806ca`, to the pinned
CRAN DuckDB 1.5.2 source. It preserves DuckDB's R API/version. The fix corrects
Brotli's invalid block-switch error handling, using the reviewed
[upstream change](https://github.com/google/brotli/commit/0ed4268a4401abcb4472096340477cf4af7d1a2b).

The installer patches the nested archive that CRAN's configure extracts, verifies
the exact two modified files and unchanged unrelated members, updates the package
checksum, and compiles afresh without prebuilt objects. The actual running
`duckdb.so` rejects the regression vector with error `-17`; its SHA-256 matches
the library tested in the immutable image. A package version or source receipt
alone was not treated as sufficient evidence.

All 35 app R files parse. R/Python dependencies, Node imports, Pillow rendering,
the installed Shiny/R/static-range smoke and palette retention checks pass.
BROTLI round trips pass on the current release's `summary/islets.parquet`
(4,853 rows, 286 columns) and `spatial/cell_distance.parquet`
(3,199,606 rows, 37 columns), preserving schemas, values, missingness, duplicate
multiplicity and aggregates. The test inputs remained unchanged.

## Data and host acceptance

The dev data bundle was rebuilt as an independent copy from the accepted current
scientific input, with the same clean `378843e` producer as the image. All
nongeometry artifacts and canonical-source hashes match accepted production
exactly. All 20 GeoPackages match in schema and every table/row, including geometry
and spatial-index data; only valid UTC `gpkg_contents.last_change` build timestamps
differ. Fresh core byte/schema/row-count verification, eight deliberate corruption
controls and a fresh inventory passed. The original accepted production parity
retains its original producer/verifier identities; this does not claim a new full
parity execution.

The running dev query matches every accepted phenotype count: **3,199,606 cells,
20 donors, 78,680 Ki67-positive cells, and 7,398 Unclassified cells (0.2312%)**.
The [independent classification review](unclassified_component_review_20261103.md)
continues to apply to these unchanged scientific contents.

Host activation verified every release byte before starting dev and verified the
local geometry/cells/spatial copies. All scientific mounts are read-only; all
three dev image aliases resolve to the new release. To fit the host, only the
unreferenced cached image `61110e9a…` was evicted after retained-archive checks,
and the independently copied, fully verified image tree from incomplete unused
v19 staging was moved into this release. Production files, container/start time,
image, service, state, manifests, mounts, aliases and global tooling stayed
unchanged. Dev uses 4 GiB alongside production's existing 9 GiB limit.

Public checks exercised authentication, anonymous denial, exact TIFF byte ranges,
sidecars, WebSocket upgrade and browser interactions. There were no browser or
application HTTP errors. These checks ran from the workstation through public
HTTPS; independent off-campus access and a host reboot were not tested. Automated
accessibility checks are incomplete, with color-contrast items still requiring
manual review. No institutional approval is claimed.

## Production remains separate

Production still runs release `pubapps-2026-11-components-v1`, source `e5689f6`,
image `89b3f840…`. Its DuckDB finding remains tracked for **2026-11-10 UTC** under
the [dated owner direction](duckdb_remediation_window_20261102.md). This dev
deployment does not patch production or extend that deadline. The
[historical production receipt](evidence/pubapps_live_20261102.json) is unchanged;
its stopped-dev observation describes its earlier verification time.

The [November 2 dev validation](dev_duckdb_fix_20261102.md) remains historical
evidence for the older `b26558a9…` image. The current release uses the fresh
`378843e` image/data build above. The security ticket reply is updated as a
[draft](container_security_20261102.md#ticket-reply-draft-only-not-sent); it has
not been sent.
