# Islet Explorer UFRC readiness and evidence index

Updated 2026-09-09. This is the review index for the supplied
[`UFIT_PubApps_Instance_Request.docx`](../UFIT_PubApps_Instance_Request.docx),
including its hosting policy, project plan, and development standards. It records
implementation separately from acceptance. It is not a signed attestation or
production approval.

## Current status

The existing records are UFRC support request `83810`, configuration reference
`SYSCONFIG-8173`, and individual assessment `IRM-456246` (submission confirmed
2026-09-01). Read-only UF Outlook search results supplied these references on
2026-09-08. The owner reports that the VM is provisioned; however, UFRC comment 4
on 2026-09-08 says its operations team normally waits for assessment completion
before VM setup. Reconcile that status with the assigned host details before
deployment. Hostname, service-account, SSH, and proxy details remain unknown.
Continue the existing requests. Trial allocation
`RC83797` for group `maigan` is recorded for 2026-08-31 through 2026-11-30.

The selected reviewer profile remains Public Restricted access to a Sensitive
collection: Plot, Trajectory, Viewer, and Spatial; AI, Statistics, and Inflammation
disabled; immutable local scientific files; no managed database, upload endpoint,
or production external API. The proposed URL is
`https://islet-explorer.rc.ufl.edu/islet-explorer/`. The backend address and port
must come from UFRC; `127.0.0.1:3838` is a request, not an established proxy route.

## Requirements mapped to evidence

The rows group the major requirements in the supplied DOCX. Detailed project
recommendations stay in the deployment/log runbooks and do not become additional
UFRC approval gates.

| Supplied requirement | Current implementation or document | Evidence still needed to close it |
|---|---|---|
| Project plan 1: purpose, audience, complete data-flow map, data categories | Populated DOCX; [handoff](../PUBAPPS_UFRC_HANDOFF.md); [SVG](pubapps_architecture.svg)/[PNG](pubapps_architecture.png) describe the planned profile and reported provisioning status | Assigned configuration and reconciled provisioning status; final reviewed form with genuine dated PI acknowledgment |
| Project plan 3–4: controls appropriate to the data and lawful data handling | Sensitive profile; whole-site proxy authentication request; read-only releases; production rejects AI credentials/debugging | Owner review of the retained asset/field inventory, authorized classification and disclosure approval, applicable IRB/DUA/MTA disposition; no Restricted asset may enter the reviewer release |
| Project plan 2: UF branding and accessibility | App skip link/landmarks and UF accessibility footer; the browser viewer URL input is disabled | Record automated and manual checks across all four tabs, keyboard/focus, labels, contrast/zoom, plots and viewer alternatives; fix defects before claiming conformance |
| Development standards: input validation and minimal, appropriate output | Allowlisted browser inputs and asset basenames; centralized fixed-field JSON events; viewer bootstrap rejects external/missing/traversal URLs before vendor import | Retained lint/security test results for the chosen revision; exact-runtime negative-input and log review |
| Development standards: lint, secure coding and static analysis | R runtime lint gate; security-boundary tests; focused deployment and viewer tests | Dated lint/SAST results with zero lint errors and documented residual warnings/findings, bound to the release revision |
| Development standards: maintained SBOM and code/package/container vulnerability review | Build/scan tooling records image identity, resolved R/Python manifests, CycloneDX and Grype results; DuckDB source/version and Rocker interpreter handling corrected | A successful exact-candidate build and scan, current scanner/database identity, remaining-CVE disposition, and browser-vendor coverage record. The project's image gate permits no High/Critical findings |
| Development standards: one event per line and UFIT-monitored logging | Central timestamped JSON emitter in `00_audit.R`; [log inventory](pubapps_log_inventory.md); dedicated worker log mount | Candidate runtime capture, native/framework output review, effective permissions, and UFRC collector destination. Record applicable access/retention rules |
| UFRC team: pre-production vulnerability/availability review and owner attestation; continuing UFIT network scans | Deployment verification and rollback tooling; review request text | UFRC production acceptance and network-scan process; exact HTTPS/auth/WebSocket/download/TIFF Range tests; genuine PI/data-owner go-live attestation |
| Conditional project-plan obligations and the existing individual assessment | Individual Risk Assessment already open; Sensitive access-review applicability awaiting UFRC | Assessment follow-up responses and disposition; quarterly access reviews if UFRC applies item 6; any additional controls explicitly assigned by the assessment |

Public UFRC guidance normally permits attaching the instance form to its main
assessment, with individual assessment when required. This project's already-open
individual assessment is the record to continue. [UFRC assessment guidance](https://docs.rc.ufl.edu/services/web_hosting/risk_assessment/).

## Implementation and validation record

The implementation is retained in draft [PR #11](https://github.com/smith6jt-cop/Islet-Explorer/pull/11).
It includes application logging, viewer containment, pinned dependencies, and
deployment verification/recovery. Source CI passed at `cc1c1bc`; the corrected
`606fb53` data release passed full verification and all 22 scientific parity checks.
Image acceptance and production startup validation are tracked separately in
the restricted release record. Neither source tests nor data parity establish
UFRC production approval.

| Check | Recorded result on 2026-09-08 | Limit |
|---|---|---|
| `node deploy/tests/test-viewer-boundary.mjs` | 31 checks passed: accepted local prefixes, malicious/missing/duplicate URLs, XML compatibility, CSP placement, and vendor hashes | Pure Node/static checks; does not prove real TIFF rendering or proxy behavior |
| Application and logging tests | 192 assertions passed across the app suite with zero skips, including 26 logging and 12 viewer-readiness assertions | Local source tests; native worker output and UFRC collection still require exact-runtime validation |
| R source parsing and runtime lint | All 28 runtime R sources parsed; `deploy/lint-runtime-r.R` passed | Local source validation before final commit; SAST and image dependency results are separate |
| Deployment shell contracts | Image evidence, deployment verifier, and deployment recovery suites all passed | Three mocked/local contract suites; does not establish behavior on the provisioned VM |
| Avivator dependency signature scan | Retire.js 5.7.0 standard/deep scans exited 0 with no reported vulnerabilities in detected signatures; JSON and recognized-component CycloneDX retained | Recognized React DOM/react-is only; incomplete transitive coverage and no exact image claim |
| DOCX integrity | All 18 archive parts retained; 129 paragraphs and XML element/attribute topology preserved; supplied policy/project-plan/standards content unchanged | Draft fields updated; no PI attestation inferred |
| Architecture rendering | PNG regenerated from editable SVG and visually checked | Requested topology; UFRC target/collector confirmation pending |
| Authored-browser SAST | Semgrep 1.176.1: 68 rules, 6 of 6 extracted/authored targets checked, zero findings/errors/skipped rules | Working-tree result; rerun after final browser edits. Vendored dependencies have a separate coverage record |
| Browser accessibility and behavior | 34 of 34 checks passed, including four tabs, keyboard navigation, four CSV downloads, actual tissue canvas, 16-byte HTTP 206, and blocked missing/external image URLs; zero axe violations across five scopes | Native production-profile rerun passed on 2026-09-09 with the verified `606fb53` data and corrected R worker launcher; Chrome 152.0.7977.82, axe 4.13.0. Manual-review items retained; no final-image or external proxy acceptance claimed |
| External HTTPS verifier | 10 synthetic tests passed without skips, including whole-site BasicAuth, bounded Range/body verification, sidecar MIME, WebSocket handshake and redirect/credential boundaries | [Acceptance procedure](../deploy/EXTERNAL_ACCEPTANCE.md) requires the assigned endpoint and real protected credential file |
| Hosted CI | Both `parse-and-test` and `deployment-quality` passed in [run 34355263510](https://github.com/smith6jt-cop/Islet-Explorer/actions/runs/34355263510) at `cc1c1bc`, including all 22 parity test assertions | Later runtime changes require their own CI and image acceptance |
| Supported Shiny Server runtime | Checksum-locked official source build and Node 24.20.0; 315 upstream and 8 real-R tests passed without skips; stripped-runtime app/static-route smoke passed; installed-package Syft/Grype scan reports zero matches | The original base-only scan found 33 High and 3 Critical entries. The rebuilt exact application image still must pass the image security gate. See [base remediation](../deploy/BASE_HARDENING.md) |
| Data verification tooling | 22 parity assertions, 14 browser extraction assertions, and 4 inventory tests passed without skips | Synthetic contract tests supplement the required real source-to-release comparison |
| Canonical source preparation | Five prepared tables measured; see counts below | All five match the corrected `606fb53` release; counts describe different scientific tables |
| Release data | `606fb53`: independent full checksum verification and all 22 real-data parity checks passed; all five summary tables had zero numeric difference; all 15 donor tables matched | Retained manifest and reports identify this data candidate; final image/runtime acceptance remains separate |
| Immutable image and SBOM/scan | Build and required security gate pending | Application images built successfully; scanning was blocked by an unprefixed Podman image ID. Strict ID normalization is corrected; the scan must complete before acceptance |
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
exact categories, counts, missingness, and infinities. The later `606fb53` release manifest, full-verification report, and successful
source-to-release comparisons are retained in the restricted release record.

## Remaining release sequence

1. Reconcile provisioning and record the remaining owner/UFRC details in request
   `83810` and assessment `IRM-456246`;
   complete classification/disclosure review and genuine PI acknowledgment.
2. Finish and retain local tests, accessibility checks, and security review; freeze
   one clean revision. Build its immutable image, manifests, SBOM, and scan reports;
   document browser-vendor coverage limitations and disposition.
3. Build the versioned data/image release, verify complete checksums and scientific
   parity, and retain counts and resource results. The chosen image and data must
   reference the same source revision.
4. Apply only the UFRC-confirmed backend/proxy configuration, transfer and activate
   the exact artifacts, then verify authentication, WebSocket, downloads, small
   TIFF ranges, logs, health/restart/rollback, and representative five-session use.
5. Retain UFRC acceptance and owner go-live approval before delivering reviewer
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
