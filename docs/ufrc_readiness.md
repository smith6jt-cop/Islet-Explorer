# Islet Explorer UFRC readiness and evidence index

Updated 2026-09-17 for the local candidate; institutional status remains as
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

### Expanded cohort candidate

The owner authorized six additions on 2026-09-11, targeting 21 donor images:
6090, 6232, 6539, `PAN-6558-AAb_test35_35p-LED.tif`, 6520-01A, and HPAP-078.
Exclude the alternate `CX_PAN-6558` acquisition and 6520-01B. Preserve HPAP-078's
HPAP provenance and separately review its applicable disclosure/usage terms;
the expanded cohort must not be described as 21 nPOD donors. Continue request
`83810` and assessment `IRM-456246` with the updated inventory through the owner.
No new UFRC acceptance or classification disposition has been established here.

InstanSeg 0.1.1 image preparation, fixed diffusion ordering and hypothesis
fitting are offline work. They do not add model training, GPU requirements, uploads, or
external API traffic to the PubApps runtime. Keep the frozen original `5809be4`
release and its source/evidence intact. Record shifts in original-donor
trajectories or clustering as scientific changes; source-to-release parity for
the expanded candidate compares against its own verified canonical preparation.

| Expanded-candidate item | Status |
|---|---|
| Input selection and cohort size | Source mix-up resolved: 6090 matches all 36 raw channel planes; 6539 reconstructed and pixel-checked against its distinct raw acquisition. The candidate targets 20 nPOD images plus one HPAP image |
| Viewer image conversion | Corrected six TIFFs total 46,877,746,573 bytes. Use `reports/image_conversion/reconciled_inventory.json`; the original duplicate 6539 TIFF is excluded |
| Image preparation and scientific integration | Candidate retains 3,592,722 raw objects, including 201 wholly unmeasured objects; 8,371 islets and 31 raw markers are available. Six primary estimates and four unavailable tests are recorded. See [local delivery](all_donor_delivery.md) |
| Cell-phenotyping provenance | Retirement policy supersedes reconstruction of all 21 labels. Retained assignments and complete parent chains are audited; unvalidated upstream gates leave subtype ratios and several endpoints unavailable |
| Fixed trajectory | Previous fit: 3,499 islets, 100% primary support, median LODO absolute Spearman 0.9143. Islet-identity review now required: 508 previously eligible candidates have no cells above exported hormone thresholds, including 143 in HPAP-078. Numerical stability does not validate biological identity |
| Release count contract | Version-2 `release.env` declares `EXPECTED_IMAGES`; generation checks the staged count. Verification/activation/rollback use the declaration; legacy version-1 releases remain fixed at 15 |
| Count and recovery tests | Local synthetic image-evidence and recovery suites pass, including 21-image activation, 15-image legacy rollback, missing/duplicate/conflicting count failures, and checksum failures; these do not validate real expanded data |
| Scientific counts and capacity | Revised runtime 164,147,382,976 bytes (152.874 GiB); frozen original plus candidate plus staging copy 414.799 GiB before review reports/OCI/logs/spill. New worker observed at about 0.70 GiB RSS; peak/concurrency and assigned-host tests pending |
| Local candidate evidence | Schema 3 at http://127.0.0.1:3839/islet-explorer/; all 31 heatmap markers, four app tabs and offline research reports. Fresh tests and hashes in the [2026-09-17 revision record](trajectory_revision_20260917.md). Ratios still lack certified gates; islet identity warning visible |
| Exact production release evidence | Pending clean source revision, matching immutable image/data identity, production source parity, SBOM/scan and assigned-runtime browser checks. Production loading rejects scientific_acceptance=false |
| Institutional and deployed acceptance | Pending updated owner-reviewed inventory/disclosure terms, UFRC assessment follow-up, assigned-host/proxy checks, and production acceptance |

The validation tables below describe the original release and their stated dates.
Their passing results are not expanded-candidate evidence.

Conversion and preparation evidence is retained under the restricted, Git-ignored
`data/development/expanded-20260911/` directory. The initial duplicate audit and
mislabeled image remain as rejected evidence. The reconciled inventory selects
the correct source paths, records distinct source/pixel signatures and verifies
current TIFF/sidecar hashes against completed pixel-validation reports. The
recovered 6539 source hash is
`1866d5e8856dc096c1829e04b3622bc796cd9b81be4225bff12b864d4df1c69e`.
Its old exclusion annotation was drawn on duplicated 6090 pixels and was preserved
separately. Corrected source paths and the raw acquisition evidence are documented
in [the expansion status](cohort_expansion_status.md).

A second scientific blocker is recorded in the historical
`7571b6d:scripts/phenotype_evaluator.py` note: the original rules that produced
the shipped ECAD-positive/Structural labels were unavailable; v2 rules were a
reconstruction requiring biological review. The current `data/phenotype_rules.csv`
cannot produce all 21 shipped labels. The revised scope requires validation of
retained labels and their upstream gates rather than recovery of retired labels.
Completing segmentation does not establish equivalent retained assignments.
Raw-marker and mask analyses proceed with explicit limitations; certified
new-cell phenotypes still require reviewed gate evidence.
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
| Development standards: maintained SBOM and code/package/container vulnerability review | Build/scan tooling records image identity, resolved R/Python manifests, CycloneDX and Grype results; DuckDB source/version and Rocker interpreter handling corrected | Exact-candidate SBOM/scan and scanner/database identity retained; the no High/Critical gate passed. UFRC pre-production review remains; lower-severity findings and incomplete browser-vendor coverage are recorded |
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

For the expanded candidate, use the reconciled six-source inventory, recover or
review the cell-phenotyping method, and finish biological review and the joint
scientific analysis. Record the final retained cohort/row counts and storage
requirements, then produce a
matching clean-revision image and copied immutable data bundle. Generate the
manifest with the explicitly reviewed final count (21 remains the target) and retain new full checksum,
source-to-release parity, SBOM/scan, browser, and capacity evidence. Do not replace
the preserved original release or reuse its passing results as candidate results.

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
