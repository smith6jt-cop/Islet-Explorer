# Islet Explorer UFRC readiness and evidence index

Updated 2026-10-22 for the uniform 20-donor release candidate (data schema 5); the 2026-10-16 update covered UFRC comment 13 and the owner-directed processing reset. This is the review index for the supplied
[`UFIT_PubApps_Instance_Request.docx`](../UFIT_PubApps_Instance_Request.docx),
including its hosting policy, project plan, and development standards. It records
implementation separately from acceptance. It is not a signed attestation or
production approval.

## Current status

UFRC support request `83810` and configuration reference `SYSCONFIG-8173` remain
the active project records. The owner confirms that the instance request form was
submitted with the original ticket on 2026-08-31. UFIT canceled individual
assessment `IRM-456246` on 2026-10-15; its closure notice says the project will be
added to the PubApps master assessment. Confirm that the original form and current
evidence are linked there, and ask UFRC whether it wants an amended form after the
final release inventory is accepted. Do not reopen the canceled assessment or
submit a duplicate form unless UFRC directs. Closure of the individual assessment
is not production approval. Trial allocation `RC83797` for group `maigan` is
recorded for 2026-08-31 through 2026-11-30.

UFRC comment 13, received 2026-10-16, says the VM will be provisioned next and
asks for the first app's DNS and proxy requirements. The project reply is drafted
in the [launch-preparation note](ufrc_launch_preparation_20261015.md); the
proxy-reachable backend address remains pending UFRC's VM assignment.

The selected reviewer profile remains Public Restricted access to a Sensitive
collection: Plot, Trajectory, Viewer, and Spatial; AI, Statistics, and Inflammation
disabled; immutable local scientific files; no managed database, upload endpoint,
or production external API. The proposed URL is
`https://islet-explorer.rc.ufl.edu/islet-explorer/`. The backend address and port
must come from UFRC; `127.0.0.1:3838` is a request, not an established proxy route.

### Current local app and release candidate (2026-10-22)

The only intended local development URL is
`http://10.15.152.7:8080/islet-explorer/`, proxied to loopback port 3840. It
was offline at the owner's direction from 2026-10-15, served a frozen copy of
the app (tag `frozen-app-20261020`) from 2026-10-20 for the owner's review of
the reprocessed data, and since 2026-10-22 serves the main tree on the uniform
root `data/development/uniform-20donor-20261022-v10/` (20 donors, 4,763
endocrine objects: 1,487 Islets and 3,276 Clusters, with 5,207 Singles counted
but not analysed). Donor 78 (HDL078 / HPAP-078) stays permanently removed from
the app and all analysis by owner directive. All 20 donors were typed and
measured by one workflow (`scripts/uniform/`, [record](../scripts/uniform/README.md)),
and every donor field comes from the nPOD export. See the
[owner directive](owner_directive_20261015.md) and the
[uniform reprocessing plan](uniform_reprocessing_plan_20261015.md).

The launch candidate is the **data-schema-5 release of that root built from
commit `d55680a`**, staged on the workstation under
`tmp/ufrc-release-20261022/` (data bundle `artifacts/release-d55680a/data/`,
built from a clean linked worktree of the commit), and since 2026-10-22 18:33 UTC
**the active release on the PubApps VM** `pubmaigan1` (Quadlet
`islet-explorer.service`, published on the VM's private address
`10.13.41.21:3838` for UFRC's proxy since 21:16 UTC). UFRC's reverse proxy with
Basic Auth to that backend, UFRC's security and availability approval and the
PI/data-owner attestation are the remaining gates before reviewers can be given
the URL. Schema 5 packages the
schema-1 summaries and spatial tables together with the uniform root's own
files and the H5AD's embedding records; `deploy/README.md` § 2 has the
commands. Its evidence is recorded in
[the candidate record below](#schema-5-release-candidate-d55680a-2026-10-22) as
each step completes. The container image is built on the workstation with
Podman, Syft 1.51.1 and Grype 0.118.0 (the pinned CI versions), so no
repository push is required for the image.

The frozen `5809be4` 15-image production artifact and its evidence are a
historical reference, superseded by this candidate. Earlier October
protein-expansion and September all-donor checks document historical roots
only. Neither this checkout nor its static image route is a PubApps release.
The current ticket, form, and IRM status are tracked in the
[launch-preparation note](ufrc_launch_preparation_20261015.md).

### Schema-5 release candidate d55680a (2026-10-22)

| Step | Result |
|---|---|
| Source commit | `d55680ae2a478ff66420c8ec57321e61501a0833` (branch `spatial-review-20261020`; the schema-5 tooling, commits `2f30ab0`, `961fdee`, `d74d2ef`, `6ffb6f2` and `d55680a`), built from the clean linked worktree `/data/panc_CODEX/Islet-Explorer-release-20261022` |
| Source root | `data/development/uniform-20donor-20261022-v10/` (v10 = v9 with donor 112's nPOD age and the 29-column nPOD donor table) |
| Data bundle | built 2026-10-22 17:20 UTC by `build_release_data.R --schema 5 --materialize copy` from the clean worktree at `d55680a` (66 min, under concurrent load): 19 artifacts; summaries 4,906 objects / 44,154 target rows / 58,872 marker rows / 12 combined pairs; 3,074,596 cell-distance rows; 300,045 signed-distance rows; 20 donor tissue tables; 20 GeoPackages with 674,717 features; 4,763 cell files; 80 image files (20 OME-TIFFs with XML and offsets sidecars, 155,515,691,675 bytes); `summary/embedding_metadata.json`. Log `tmp/ufrc-release-20261022/build-data-d55680a.log`. Three earlier bundles (`961fdee`, `d74d2ef`, `6ffb6f2`) were superseded by tooling fixes found by the checks below and removed |
| Full-hash verification (`verify_release_data.R --mode full --schema 5`) | passed 2026-10-22 17:31 UTC: every artifact hash and Parquet row count |
| Scientific parity (`verify_release_parity.R`) | passed 2026-10-22 17:58 UTC, 23 of 23 checks: canonical H5AD preparation, manifest and full hashes, source identity (clean `d55680a`), canonical source hashes, the five summary tables, passthrough files and trees, embedding metadata, signed-distance table, 20 donor tissue tables, 8 spatial histogram and composition families over both distance columns, and 60 geometry windows. Report `tmp/ufrc-release-20261022/evidence/release-parity-d55680a.{json,md}` |
| Inventory (`inventory_release_data.py`, scvi-env Python) | completed 2026-10-22 17:58 UTC, bound to the parity report: 19 release artifacts, 20 OME-TIFFs / 690 pages, names and counts only. `tmp/ufrc-release-20261022/evidence/release-inventory-d55680a.{json,md}` |
| Container image, SBOM and High/Critical gate | passed 2026-10-22 16:39 UTC on the workstation (rootless Podman 4.9.3, Syft 1.51.1, Grype 0.118.0): image `sha256:ae95fd238c72e0ec239942be0050ad72eb441bf34149218973edd41f5c7e9af9` from `d55680a` on base `rocker/shiny@sha256:dd10b638…`; upstream Shiny Server tests 315 + 8 passed; Python and R probes passed; Grype 0 High / 0 Critical (4,255 Medium, 189 Low, 196 Negligible); OCI archive `islet-explorer-pubapps-2026-10-uniform-v10-d55680a.oci.tar` (2.2 GB) with metadata and checksum. The first build (`6ffb6f2`) had failed the gate on handlebars, proxy-addr, compression, follow-redirects, pillow and certifi, all pinned since (see `deploy/BASE_HARDENING.md`). Evidence: `tmp/ufrc-release-20261022/evidence/image-scan-d55680a/` |
| Release manifest (`generate-release-manifest.sh ... 5 20`) and `verify-release.sh` | passed 2026-10-22 18:01 UTC: release `pubapps-2026-10-uniform-v10` bound to image `ae95fd23…` and commit `d55680a` (`MANIFEST_VERSION=2`, `DATA_SCHEMA=5`, `EXPECTED_IMAGES=20`, 4,819 data files / 917,879,043 bytes, 80 image files / 155,515,691,675 bytes); nine attestation files (SBOM, raw SBOM, corrections, Grype report and database status, R and Python manifests, scan attestation, published image); `MANIFEST.sha256` 4,909 entries, manifest hash `96423b92cbaba71a7d44f2bfc10f8310c76817f671c4a27e62798c2b9bb4a9fc`; full verification stamp written on the workstation. `tmp/ufrc-release-20261022/release-status.json` summarizes the candidate |
| Staging on the VM | complete 2026-10-22 18:11 UTC under `/pubapps/maigan/islet-explorer/releases/pubapps-2026-10-uniform-v10/`: the 80 image files (155.5 GB, through the HPG session at 40-107 MB/s, then matched by checksum against the bound bundle with nothing re-sent), the 4,819 data files, the nine attestations, `release.env` and `MANIFEST.sha256`; the OCI archive in `incoming/`, loaded by `load-image.sh` (archive checksum, image ID, revision label and base verified) |
| Activation on the VM | **activated 2026-10-22 18:33 UTC** (`activate-release.sh`, deploy tooling `38e87d4`): `verify-release.sh` re-hashed the whole release on the VM (manifest `96423b92…`, 20 OME-TIFFs), the in-container `verify_release_data.R --mode core --schema 5` passed, the keep-id probe mapped `maigan` to the image's `shiny` (997), the Quadlet `islet-explorer.service` started and is active and enabled under the lingering user session, and `verify-deployment.sh` passed: healthy container, backend bound to `127.0.0.1:3838` only, data and images read-only, feature flags and credentials in the reviewer profile, Python and R dependency probes, OME-TIFF Range response `206 bytes 0-15/7792217538`. Two earlier attempts failed on the VM's Podman 5.8 and were fixed in the deploy tooling without touching the release: tmpfs `uid=`/`gid=` options are rejected (now `U`, commit `2f36aaf`), and a `podman exec` does not inherit the entrypoint's cache directories, so numba could not import scanpy on the read-only root (unit environment, commit `38e87d4`); the CPU quota is omitted because the session delegates only `memory` and `pids` (`acc1ce6`). **Re-activated 2026-10-22 21:16 UTC** by the owner with `--publish-address 10.13.41.21` (tooling `9f29c4d`): the release was re-hashed on the VM, the service restarted, and `verify-deployment.sh` passed with the backend bound to `10.13.41.21:3838` only (Range `206 bytes 0-15/7792217538`). Logs: `tmp/ufrc-release-20261022/deploy-vm-d55680a*.log` |
| Browser checks against the deployed container | passed 2026-10-22 18:37-18:39 UTC from the workstation through a two-hop SSH tunnel to the VM's loopback backend: `check-app.mjs` (every check, no accessibility rule violations, no browser or HTTP errors), `check-ui-review-20261021.mjs` (15 of 15: logo, plot types in both cards, composition by donor over 20 donors, near-empty measures hidden, axis fit, frozen controls on four tabs, Spatial colourings including the nPOD fields, hover and click linking, shared selection, clear) and `check-trajectory-links.mjs` (all steps, no errors). Reports: `tmp/ufrc-release-20261022/evidence/browser-vm-d55680a/` |
| UFRC VM | provisioned 2026-10-22: `pubmaigan1.ufhpc`, service account `maigan` (home `/pubapps/maigan`, 1 TiB Qumulo share), RHEL 9.8, rootless Podman 5.8.2 (store on the 40 GB `/local`), cgroups v2, linger on, 2 vCPU / 15 GiB; reached from HPG login nodes as `ssh pubmaigan1`; project path `/pubapps/maigan/islet-explorer/` created, release tooling copied, 500 MB test transfer at about 100 MB/s |
| UFRC proxy and Basic Auth acceptance | pending (ticket `83810`). On 2026-10-22 the hostname resolves to UFRC's proxy `vpubweb.rc.ufl.edu` (wildcard `*.rc.ufl.edu` certificate) but no vhost forwards `/islet-explorer/` yet (Apache default 404 on every path). The proxy is a separate host and cannot reach the VM loopback, so the backend is the VM's private address `10.13.41.21:3838`, published since 21:16 UTC (`activate-release.sh --publish-address`, tooling `9f29c4d`; the verifier requires that address as the exclusive listener). Measured before publishing, because the application relies on the proxy's Basic Auth: an HPG login node reaches the VM on port 22 only (a test listener on `10.13.41.21:3839` received no connection, port 3838 timed out as well), the campus workstation cannot reach the VM at all, and the VM's own firewalls are inactive, so the port is filtered by UFRC's network; the owner also judged HPG's MFA login stricter than the app's Basic Auth. Still pending from UFRC: the vhost with Basic Auth for `/islet-explorer/` and confirmation that the proxy reaches the backend. SSH access list requested: `smith6jt` (VM keys received from UFRC by email on 2026-10-22) |

### Historical regenerated 21-image candidate

The September 21 schema-4 results below describe an earlier all-cohort
regeneration. Its original-donor geometry and trajectory do not satisfy the
current preservation contract. The results remain experimental evidence,
not current app counts or production acceptance.

| Historical candidate item | Recorded result and limit |
|---|---|
| Input selection and cohort size | Source mix-up resolved: 6090 matches all 36 raw channel planes; 6539 reconstructed and pixel-checked against its distinct raw acquisition. The candidate targets 20 nPOD images plus one HPAP image |
| Viewer image conversion | Corrected six TIFFs totaled 46,877,746,573 bytes. The historical `reports/image_conversion/reconciled_inventory.json` excluded the duplicate 6539 TIFF; it is not an approved production image manifest |
| Image preparation and scientific integration | The regenerated candidate retained 3,619,534 canonical raw objects and 3,592,722 archival source objects, including 248 wholly unmeasured canonical objects. It reported 11,944 islet annotations and 31 shared markers. These are not preservation-compatible expansion counts. See [local delivery](all_donor_delivery.md) |
| Cell-phenotyping provenance | Schema 4 used experimental typed-cell outputs; these do not establish equivalence to the original typing rules, gates and scVI model. Unsupported identities and marker-positive calls remained unavailable rather than becoming zeros |
| Boundary review and integration | The regenerated geometry covered all 21 donors, including 1,431 corrected HPAP-078 annotations. Its original-donor boundaries were not preserved, and scientific/deployment acceptance remained false |
| Fixed trajectory and Spatial embedding | The experimental classifier flagged 11,136 eligible islets across all 21 donors, with at least 137 per donor; typing did not revoke eligibility. Typed composition was available for 7,493 islets. The diffusion graph trained on 1,514 valid islets and projected feature-supported islets. The primary trajectory supported 7,490 coordinates, median leave-one-donor-out absolute Spearman 0.754 and 67.3% classifier-eligible coverage, below its prespecified stability/coverage checks. Spatial UMAP and four Leiden resolutions used 7,493 feature-supported rows. Orientation was arbitrary; these are not an accepted preservation-compatible 21-donor refit |
| Release count contract | Version-2 `release.env` declares `EXPECTED_IMAGES`; generation checks the staged count. Verification/activation/rollback use the declaration; legacy version-1 releases remain fixed at 15 |
| Count and recovery tests | Local synthetic image-evidence and recovery suites passed, including 21-image activation, 15-image legacy rollback, missing/duplicate/conflicting count failures, and checksum failures; these did not validate real expanded data |
| Scientific counts and capacity | The historical runtime measured 188,900,659,274 bytes (175.927 GiB); two full copies would require 377,801,318,548 bytes before frozen releases/reports/OCI/logs/spill. This is not a capacity estimate for a preservation-compatible expansion. Full-worker peak/concurrency and assigned-host tests remained pending |
| Local candidate evidence | The schema-4 candidate `all-donor-20260922-trajectory-fit-support-v1` previously used the local URL now serving the restored 15-donor app. Its validator and browser checks recorded 21 donors, 31 markers, 11,944 annotations, 11,136 classifier-eligible islets, 7,490 ranked trajectory coordinates, 7,493 Spatial UMAP points and click-through/download behavior. Its display ranks preserved 1,729 unique raw values and a 2,234-islet raw tie for audit/inference. Those checks do not certify the current app or UFRC acceptance |
| Exact production release evidence | No matching immutable image/data release, source parity, SBOM/scan, or assigned-runtime browser evidence was completed for this candidate. Production loading rejects `scientific_acceptance=false` |
| Institutional and deployed acceptance | Pending updated owner-reviewed inventory/disclosure terms, UFRC assessment follow-up, assigned-host/proxy checks, and production acceptance |

The validation tables below describe the original release and their stated dates.
Their passing results do not certify the current development UI or any expansion.

Conversion and preparation evidence is retained under the restricted, Git-ignored
`data/development/expanded-20260911/` directory. The initial duplicate audit and
mislabeled image remain as rejected evidence. The reconciled inventory recorded
the selected source paths, distinct source/pixel signatures, and TIFF/sidecar
hashes against completed pixel-validation reports. The
recovered 6539 source hash is
`1866d5e8856dc096c1829e04b3622bc796cd9b81be4225bff12b864d4df1c69e`.
Its old exclusion annotation was drawn on duplicated 6090 pixels and was preserved
separately. Corrected source paths and the raw acquisition evidence are documented
in [the expansion status](cohort_expansion_status.md).

Historical method-recovery work found that the original phenotype rules and
trained scVI model behind the first 15 donors were lost. On 2026-10-15, the
owner directed that all 20 retained donors be reprocessed through one
reproducible workflow rather than append results from a different method. The
uniform plan is proposed and its compute work has not started; D1–D5 and the
acceptance gates remain open. See the [owner directive](owner_directive_20261015.md)
and [uniform reprocessing plan](uniform_reprocessing_plan_20261015.md). This is
a scientific-provenance requirement, not an additional UFRC rule.

## Requirements mapped to evidence

The rows group the major requirements in the supplied DOCX. Detailed project
recommendations stay in the deployment/log runbooks and do not become additional
UFRC approval gates.

| Supplied requirement | Current implementation or document | Evidence still needed to close it |
|---|---|---|
| Project plan 1: purpose, audience, complete data-flow map, data categories | Instance request form submitted with ticket `83810` on 2026-08-31 (owner-confirmed); [handoff](../PUBAPPS_UFRC_HANDOFF.md); [SVG](pubapps_architecture.svg)/[PNG](pubapps_architecture.png) describe the planned profile | Confirm form linkage to the master assessment; ask UFRC whether it wants an amendment after the final inventory; verify actual PI acknowledgment and assigned configuration |
| Project plan 3–4: controls appropriate to the data and lawful data handling | Sensitive profile; whole-site proxy authentication request; read-only releases; production rejects AI credentials/debugging | Owner review of the retained asset/field inventory, authorized classification and disclosure approval, applicable IRB/DUA/MTA disposition; no Restricted asset may enter the reviewer release |
| Project plan 2: UF branding and accessibility | App skip link/landmarks and UF accessibility footer; the browser viewer URL input is disabled | Record automated and manual checks across all four tabs, keyboard/focus, labels, contrast/zoom, plots and viewer alternatives; fix defects before claiming conformance |
| Development standards: input validation and minimal, appropriate output | Allowlisted browser inputs and asset basenames; centralized fixed-field JSON events; viewer bootstrap rejects external/missing/traversal URLs before vendor import | Source CI and native production-profile negative-input/log evidence retained; final deployed runtime and UFRC log-collection checks remain |
| Development standards: lint, secure coding and static analysis | R runtime lint gate; security-boundary tests; focused deployment and viewer tests | Dated CI/lint/SAST results retained for `5809be4`; browser-vendor coverage limitations remain explicit |
| Development standards: maintained SBOM and code/package/container vulnerability review | Build/scan tooling records image identity, resolved R/Python manifests, CycloneDX and Grype results; DuckDB source/version and Rocker interpreter handling corrected | The dated `5809be4` image passed the no High/Critical gate. A new exact-release SBOM/scan and scanner/database identity remain required, as does UFRC pre-production review; lower-severity findings and incomplete browser-vendor coverage remain recorded |
| Development standards: one event per line and UFIT-monitored logging | Central timestamped JSON emitter in `00_audit.R`; [log inventory](pubapps_log_inventory.md); dedicated worker log mount | Candidate runtime capture, native/framework output review, effective permissions, and UFRC collector destination. Record applicable access/retention rules |
| UFRC team: pre-production vulnerability/availability review and owner attestation; continuing UFIT network scans | Deployment verification and rollback tooling; review request text | UFRC production acceptance and network-scan process; exact HTTPS/auth/WebSocket/download/TIFF Range tests; genuine PI/data-owner go-live attestation |
| Conditional project-plan obligations and master-assessment routing | Individual assessment `IRM-456246` canceled 2026-10-15; UFIT says project will be added to the PubApps master assessment | Confirm linkage and any remaining review/control requirements with UFRC; quarterly access reviews if UFRC applies item 6 |

Public UFRC guidance normally permits attaching the instance form to its main
assessment, with individual assessment when UFRC requires it. The original form
was submitted with ticket `83810` on 2026-08-31; confirm its master-assessment
linkage and requested update path through that ticket. [UFRC assessment guidance](https://docs.rc.ufl.edu/services/web_hosting/risk_assessment/).

## Implementation and validation record

The implementation is retained in draft [PR #11](https://github.com/smith6jt-cop/Islet-Explorer/pull/11).
It includes application logging, viewer containment, pinned dependencies, and
deployment verification/recovery. The frozen release revision is `5809be4`: source CI, the immutable container
High/Critical gate, OCI transfer checksums, full data hashes, all 22 scientific
parity checks, and the assembled release checksums passed. The matching container,
data bundle, inventory, and evidence are retained in the restricted release record. Neither source tests nor data parity establish
UFRC production approval.

| Check | Recorded result | Limit |
|---|---|---|
| `node deploy/tests/test-viewer-boundary.mjs` | 31 checks passed: accepted local prefixes, malicious/missing/duplicate URLs, XML compatibility, CSP placement, and vendor hashes | Pure Node/static checks; does not prove real TIFF rendering or proxy behavior |
| Application and logging tests | 192 assertions passed across the app suite with zero skips, including 26 logging and 12 viewer-readiness assertions | Local source tests; native worker output and UFRC collection still require exact-runtime validation |
| R source parsing and runtime lint | All 28 runtime R sources parsed; `deploy/lint-runtime-r.R` passed | Local source validation before final commit; SAST and image dependency results are separate |
| Deployment shell contracts | Image evidence, deployment verifier, and deployment recovery suites all passed | Three mocked/local contract suites; does not establish behavior on the provisioned VM |
| Avivator dependency signature scan | Retire.js 5.7.0 standard/deep scans exited 0 with no reported vulnerabilities in detected signatures; JSON and recognized-component CycloneDX retained | Recognized React DOM/react-is only; incomplete transitive coverage and no exact image claim |
| DOCX integrity | All 18 archive parts retained; 129 paragraphs and XML element/attribute topology preserved; supplied policy/project-plan/standards content unchanged | Draft fields updated; no PI attestation inferred |
| Architecture rendering | PNG regenerated from editable SVG and visually checked | Requested topology; UFRC target/collector confirmation pending |
| Authored-browser SAST | Semgrep 1.176.1: 68 rules, 6 of 6 extracted/authored targets checked, zero findings/errors/skipped rules | Retained hosted CI result at `5809be4`; vendored dependencies have a separate coverage record |
| Browser accessibility and behavior | 34 of 34 checks passed, including four tabs, keyboard navigation, four CSV downloads, actual tissue canvas, 16-byte HTTP 206, and blocked missing/external image URLs; zero axe violations across five scopes | Native production-profile rerun passed on 2026-09-09 with the verified `606fb53` data and corrected R worker launcher; Chrome 152.0.7977.82, axe 4.13.0. Manual-review items retained; no final-image or external proxy acceptance claimed |
| External HTTPS verifier | 10 synthetic tests passed without skips, including whole-site BasicAuth, bounded Range/body verification, sidecar MIME, WebSocket handshake and redirect/credential boundaries | [Acceptance procedure](../deploy/EXTERNAL_ACCEPTANCE.md) requires the assigned endpoint and real protected credential file |
| Hosted CI | Both `parse-and-test` and `deployment-quality` passed in [run 34375550621](https://github.com/smith6jt-cop/Islet-Explorer/actions/runs/34375550621) at `5809be4`, including all 22 parity assertions and 19 authentication/SBOM tests | Later runtime changes require their own CI and image acceptance |
| Later hosted CI at `000e970` | [Source run 34385986809](https://github.com/smith6jt-cop/Islet-Explorer/actions/runs/34385986809) and [image run 34385986858](https://github.com/smith6jt-cop/Islet-Explorer/actions/runs/34385986858) failed before application checks because the runner's Chrome package index had a hash mismatch | Local workflows now use Ubuntu 24.04 repositories for their Ubuntu packages. YAML/shell parsing and offline APT repository isolation pass; hosted rerun pending |
| Supported Shiny Server runtime | Checksum-locked official source build and Node 24.20.0; 315 upstream and 8 real-R tests passed without skips; stripped-runtime app/static-route smoke passed; installed-package Syft/Grype scan reports zero matches | The original base-only scan found 33 High and 3 Critical entries. The complete `5809be4` application image subsequently passed the High/Critical gate. See [base remediation](../deploy/BASE_HARDENING.md) |
| Data verification tooling | 22 parity assertions, 14 browser extraction assertions, and 4 inventory tests passed without skips | Synthetic contract tests supplement the required real source-to-release comparison |
| Canonical source preparation | Five prepared tables measured; see counts below | All five match the final `5809be4` release; counts describe different scientific tables |
| Release data | `5809be4`: full data hashes and all 22 real-data parity checks passed; all five summary tables had zero numeric difference and all 15 donor tables matched. Inventory: 5,317 files, 117,092,780,533 bytes. Assembled release checksums also passed | The data and container identify the same clean source revision. These checks do not establish deployed UFRC acceptance |
| Immutable image and SBOM/scan | `5809be4`: build, High/Critical gate, OCI export and checksums passed in [run 34375550622](https://github.com/smith6jt-cop/Islet-Explorer/actions/runs/34375550622) on 2026-09-09 | Unused system build tools were removed; one exact R VIM/editor CPE correction is recorded alongside the original SBOM. Lower-severity findings and browser-vendor coverage limits remain in the evidence. UFRC deployment acceptance is separate |
| Production proxy and operational acceptance | Pending UFRC configuration and deployed checks | VM provisioning alone does not establish reviewer readiness |

Local evidence is retained under the restricted, Git-ignored
`tmp/ufrc-release-20260908/evidence/` directory. It must also be copied to the
approved durable operating record with the final release artifacts.

Keep raw logs, contact/key details, credentials, donor inventories, agreement
approvals, and assessment exports in the restricted operating record. Store
sanitized evidence references here, not sensitive evidence contents.

### Rejected candidate and conversion correction

Candidate `537c0e3` is not acceptable for release parity. The real donor-table
probe found that CSV conversion changed 2,295,503 blank `islet_name` strings to
missing values across the 15 donor tables. Canonical files were not changed.
The builder now preserves blanks in that character field while retaining numeric
missingness. A read-only probe of all 15 corrected tables matches exactly; the
regression fixture covers blank strings, literal `NA`, numeric missingness and
partition columns. The corrected `606fb53` release passed independent full checksum verification
and all 22 real-data parity checks on 2026-09-08. Its manifest SHA-256 is
`32badc77220273f39815690dcddacfaf857c8237aae2cf483d21ec5428b1a80c`.

The raw-table comparator also now disables DuckDB's inferred Hive partition
column; separate runtime-query comparisons continue to exercise the application's
normal reader. This prevents an extra inferred `case_id` from masking the actual
stored-table comparison.

### Scientific counts and parity scope

The 2026-09-08 source preparation record reports the following distinct counts.
The trajectory H5AD's observation count is not the row count of every prepared
runtime table.

| Artifact or prepared table | Count |
|---|---:|
| Trajectory H5AD observations | 5,214 |
| `core_area` prepared rows | 5,382 |
| `comp` prepared rows | 5,382 |
| `targets_all` prepared rows | 48,438 |
| `markers_all` prepared rows | 64,584 |
| `combined_pairs` definitions | 21 |
| Donor image cohort | 15 donors |

Local source evidence: `/tmp/islet-parity-source-check.json` and
`/tmp/islet-parity-source-check.md`, created 2026-09-08T19:01:08Z. The record names
base commit `e58996168198187840636901acdaad723317fdb5` with `git_dirty=true`
and status `source_prepared_release_pending`; copy it to the restricted durable
release record. Comparison settings are absolute numeric tolerance `1e-12` with
exact categories, counts, missingness, and infinities. The final `5809be4` manifest, full-hash/parity report, inventory, and release
verification record are retained in the restricted release record; earlier
`606fb53` evidence remains historical.

## Remaining release sequence

1. Continue current UFRC coordination on ticket `83810` / `SYSCONFIG-8173`:
   confirm that the August 31 instance form is linked to the PubApps master
   assessment, and ask whether UFRC wants an amendment. Individual assessment
   `IRM-456246` was canceled on 2026-10-15; do not reopen it unless UFRC directs.
2. Keep the local app offline while the owner resolves decisions D1–D5 and the
   20-donor workflow is implemented and run. Donor 78/HDL078 remains permanently
   excluded. Follow the gates and acceptance criteria in the
   [uniform reprocessing plan](uniform_reprocessing_plan_20261015.md); the
   October 15 v3 and frozen 15-donor root are comparison history, not launch
   candidates.
3. After the owner accepts the reprocessed data, create a fresh release
   inventory and update the submitted request only in the way UFRC specifies.
   Then produce candidate-specific image/data artifacts, hashes, source parity,
   SBOM/scan, browser, capacity, and operational evidence.
4. Confirm the assigned host, account, storage quotas, SSH route, proxy/backend,
   authentication, logging, and ongoing scan process with UFRC. Retain UFRC
   production acceptance and PI/data-owner go-live approval before bringing up
   a production service or delivering reviewer credentials. Arrange allocation
   continuation before the recorded trial term ends on 2026-11-30.

## Browser dependency coverage

The checked-in Avivator bundle has partial version strings and an exact asset hash
inventory, but no verified upstream revision or complete dependency lockfile. See
[`VENDOR_PROVENANCE.md`](../app/shiny_app/www/avivator/VENDOR_PROVENANCE.md).
Do not label a container SBOM or an empty scanner report as complete transitive
JavaScript coverage.

Use a pinned Retire.js CLI as supplemental detection against the exact deployed
JavaScript files. Retain its version, local signature database and hash, JSON
findings, and verbose detected-library inventory; its supported `--path`,
`--jsrepo`, and `--outputformat` options allow an offline repeat of the scan.
Experimental `--deep` scanning may find additional bundled libraries but does
not establish completeness. [Retire.js CLI documentation](https://github.com/RetireJS/retire.js/blob/master/node/README.md).

The local 2026-09-08 scan used Retire.js **5.7.0** and signature repository
commit `7ddda6314785f7e4dbca442296aa399520e13c13` (database SHA-256
`79b3fc3efb2c880bcf9a5745d08db7067fd4b65c74de470e0f433ece823915d7`).
Standard detection found React DOM 18.3.1 and react-is 16.13.1; experimental
deep detection additionally found the React DOM prerelease string
18.3.1-next-f1338f8080-20240426. Both scans reported no vulnerability entries
and no scanner errors for those signatures. Other bundled packages remain
unaccounted for by these results.

Reports, database snapshot, scanner lockfile, source hashes, and provenance are
in `/tmp/islet-ufrc-evidence/browser-vendor/`; copy them to the restricted durable
release record. The scan ran on the implementation working tree before its final
commit and used a downloaded local signature database. No application or
scientific data were uploaded.

Review findings and unmatched/minified components explicitly with the existing
assessment. The durable way to establish the full dependency graph is to recover
the actual source revision and lockfile or produce a pinned, reviewed upstream
rebuild with scientific viewer parity. Until then, keep unidentified dependency
coverage as an unresolved evidence limitation.
