# Islet Explorer UFRC readiness and evidence index

Updated 2026-09-28 for the restored 15-donor local app; institutional status remains as
last verified on 2026-09-16. This is the review index for the supplied
[`UFIT_PubApps_Instance_Request.docx`](../UFIT_PubApps_Instance_Request.docx),
including its hosting policy, project plan, and development standards. It records
implementation separately from acceptance. It is not a signed attestation or
production approval.

## Current status

The existing records are UFRC support request `83810`, configuration reference
`SYSCONFIG-8173`, and individual assessment `IRM-456246` (submission confirmed
2026-09-01). Read-only UF Outlook search on 2026-09-16 found an assigned ISM's
request that day for a brief Teams meeting for the risk review. Scheduling and
completion are unconfirmed. Comment 8, received 2026-09-11, had reported the ISO
waiting for local ISM review. Comment 7 requested an assessment technical contact;
completion is not established. The owner previously reported a provisioned VM,
but hostname, service-account, SSH and proxy details remain unconfirmed.
Continue the existing requests. Trial allocation
`RC83797` for group `maigan` is recorded for 2026-08-31 through 2026-11-30.

The selected reviewer profile remains Public Restricted access to a Sensitive
collection: Plot, Trajectory, Viewer, and Spatial; AI, Statistics, and Inflammation
disabled; immutable local scientific files; no managed database, upload endpoint,
or production external API. The proposed URL is
`https://islet-explorer.rc.ufl.edu/islet-explorer/`. The backend address and port
must come from UFRC; `127.0.0.1:3838` is a request, not an established proxy route.

### Current local app and raw-only additions

The only intended local development URL is
`http://10.15.152.7:8080/islet-explorer/`, proxied to loopback port 3840.
The app uses the restored original 15 donors and 5,214 islets. Its image route
temporarily points to an earlier 21-image collection; the owner permits that
development-only mismatch, not an expanded release. Neither this checkout nor
its image route is a PubApps release.
Keep the dated `5809be4` 15-image production artifact and its evidence separate
from the current development app. See [the September 28 handoff](HANDOFF_2026-09-28.md).

The six donor additions—HPAP-078, 6090, 6232, 6520, recovered 6539, and
6558—are staged only as raw offline evidence. Their bundle contains 943,068
objects, including 201 wholly unmeasured objects; 852 zero-tolerance checks
preserve the original donor artifacts. Its recorded `app_ready` and
`scientific_acceptance` values are both false. No additions-only phenotype,
normalized-expression, scVI, islet-composition, or compatible trajectory
method has been accepted. HPAP-078 boundaries are owner-approved; the other
five donors still need donor-specific visual review. There is no accepted
combined islet count, optional app loader, or 21-image release. See the
[phase-2 contract](donor_expansion_phase2.md). Continue request `83810` and
assessment `IRM-456246` with an approved final inventory through the owner;
no new UFRC acceptance or classification disposition has been established.

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
| Cell-phenotyping provenance | Schema 4 used experimental typed-cell outputs; these do not establish an accepted additions-only typing method. Unsupported identities and marker-positive calls remained unavailable rather than becoming zeros |
| Boundary review and integration | The regenerated geometry covered all 21 donors, including 1,431 corrected HPAP-078 annotations. Its original-donor boundaries were not preserved, and scientific/deployment acceptance remained false |
| Fixed trajectory and Spatial embedding | The experimental classifier flagged 11,136 eligible islets across all 21 donors, with at least 137 per donor; typing did not revoke eligibility. Typed composition was available for 7,493 islets. The diffusion graph trained on 1,514 valid islets and projected feature-supported islets. The primary trajectory supported 7,490 coordinates, median leave-one-donor-out absolute Spearman 0.754 and 67.3% classifier-eligible coverage, below its prespecified stability/coverage checks. Spatial UMAP and four Leiden resolutions used 7,493 feature-supported rows. Orientation was arbitrary; none of these coordinates is an accepted projection onto the preserved 15-donor reference |
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

A second scientific blocker is recorded in the historical
`7571b6d:scripts/phenotype_evaluator.py` note: the original rules that produced
the shipped ECAD-positive/Structural labels were unavailable; v2 rules were a
reconstruction requiring biological review. The current `data/phenotype_rules.csv`
cannot produce all 21 shipped labels. The current phase-2 contract requires an
independently validated additions-only typing method that preserves the original
donor records. Completing segmentation does not establish equivalent phenotype
assignments. Raw-marker and mask analyses retain explicit limitations;
certified new-cell phenotypes still require reviewed gate evidence.
The tested transfer models did not meet their predeclared accuracy/F1/recall
screen; their design and per-class results are retained. This is a scientific
provenance issue, not a newly inferred UFRC rule.

## Requirements mapped to evidence

The rows group the major requirements in the supplied DOCX. Detailed project
recommendations stay in the deployment/log runbooks and do not become additional
UFRC approval gates.

| Supplied requirement | Current implementation or document | Evidence still needed to close it |
|---|---|---|
| Project plan 1: purpose, audience, complete data-flow map, data categories | Populated DOCX; [handoff](../PUBAPPS_UFRC_HANDOFF.md); [SVG](pubapps_architecture.svg)/[PNG](pubapps_architecture.png) describe the planned profile and reported provisioning status | Assigned configuration and reconciled provisioning status; final reviewed form with genuine dated PI acknowledgment |
| Project plan 3–4: controls appropriate to the data and lawful data handling | Sensitive profile; whole-site proxy authentication request; read-only releases; production rejects AI credentials/debugging | Owner review of the retained asset/field inventory, authorized classification and disclosure approval, applicable IRB/DUA/MTA disposition; no Restricted asset may enter the reviewer release |
| Project plan 2: UF branding and accessibility | App skip link/landmarks and UF accessibility footer; the browser viewer URL input is disabled | Record automated and manual checks across all four tabs, keyboard/focus, labels, contrast/zoom, plots and viewer alternatives; fix defects before claiming conformance |
| Development standards: input validation and minimal, appropriate output | Allowlisted browser inputs and asset basenames; centralized fixed-field JSON events; viewer bootstrap rejects external/missing/traversal URLs before vendor import | Source CI and native production-profile negative-input/log evidence retained; final deployed runtime and UFRC log-collection checks remain |
| Development standards: lint, secure coding and static analysis | R runtime lint gate; security-boundary tests; focused deployment and viewer tests | Dated CI/lint/SAST results retained for `5809be4`; browser-vendor coverage limitations remain explicit |
| Development standards: maintained SBOM and code/package/container vulnerability review | Build/scan tooling records image identity, resolved R/Python manifests, CycloneDX and Grype results; DuckDB source/version and Rocker interpreter handling corrected | The dated `5809be4` image passed the no High/Critical gate. A new exact-release SBOM/scan and scanner/database identity remain required, as does UFRC pre-production review; lower-severity findings and incomplete browser-vendor coverage remain recorded |
| Development standards: one event per line and UFIT-monitored logging | Central timestamped JSON emitter in `00_audit.R`; [log inventory](pubapps_log_inventory.md); dedicated worker log mount | Candidate runtime capture, native/framework output review, effective permissions, and UFRC collector destination. Record applicable access/retention rules |
| UFRC team: pre-production vulnerability/availability review and owner attestation; continuing UFIT network scans | Deployment verification and rollback tooling; review request text | UFRC production acceptance and network-scan process; exact HTTPS/auth/WebSocket/download/TIFF Range tests; genuine PI/data-owner go-live attestation |
| Conditional project-plan obligations and the existing individual assessment | Individual Risk Assessment already open; Sensitive access-review applicability awaiting UFRC | Assessment follow-up responses and disposition; quarterly access reviews if UFRC applies item 6; any additional controls explicitly assigned by the assessment |

Public UFRC guidance normally permits attaching the instance form to its main
assessment, with individual assessment when required. This project's already-open
individual assessment is the record to continue. [UFRC assessment guidance](https://docs.rc.ufl.edu/services/web_hosting/risk_assessment/).

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

For the 15-donor reviewer release, select and freeze the exact app revision,
then produce a matching immutable image and copied schema-1 data/image bundle
with 15 approved images. Retain new full checksums, source-to-release parity,
SBOM/scan, browser, and capacity evidence for that exact revision. The dated
`5809be4` evidence remains a reference; it does not certify the current UI.

For any later expansion, finish the five donor-specific boundary reviews and
validate additions-only typing and compatible trajectory projection while
preserving all original donor records exactly. Only after scientific and owner
acceptance should a new final image/data inventory, count, capacity estimate,
and matching release evidence be documented. Do not replace the preserved
original release or reuse its passing results as expansion evidence.

1. Reconcile provisioning and record the remaining owner/UFRC details in request
   `83810` and assessment `IRM-456246`;
   complete classification/disclosure review and genuine PI acknowledgment.
2. Apply the UFRC-confirmed backend/proxy configuration, transfer and activate
   the exact artifacts, then verify authentication, WebSocket, downloads, small
   TIFF ranges, logs, health/restart/rollback, and representative five-session use.
3. Retain UFRC acceptance and owner go-live approval before delivering reviewer
   credentials. Assign both maintainers, incident/rotation responsibilities, and
   continuation before the trial expires on 2026-11-30.

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
