# Panc Floyd / Islet Explorer UFRC PubApps Request and Reviewer-Hosting Handoff

This document is the actionable companion to
`UFIT_PubApps_Instance_Request.docx`. It provides copy-ready form and support-ticket
text, identifies information that the project owner must supply, and records the
operational decisions for the journal-review deployment. Do not submit it with
bracketed placeholders left unresolved.

## Deployment decisions

| Item | Decision |
|---|---|
| Service | One trial PubApps Instance (`PA-INSTANCE`) allocated to `maigan`; reference `RC83797` |
| Current allocation term | 2026-08-31 through 2026-11-30; renew or purchase before reviewer access extends beyond that date |
| Initial resources | 2 vCPU, 16 GiB RAM, 1 TiB `/pubapps`, standard local container storage |
| Add-ons | None initially; scale only after a measured five-session load test |
| Access during peer review | **Public Restricted**, with whole-site shared Basic Auth |
| UF data classification during peer review | **Sensitive** |
| Expected audience | 1–8 editors/reviewers total, up to 5 concurrent sessions, plus 2 named UF maintainers |
| Reviewer-visible features | Plot, Trajectory, Viewer, and Spatial |
| Disabled features | AI, Statistics, and Inflammation |
| Production URL request | `https://islet-explorer.rc.ufl.edu/islet-explorer/` |
| Production backend | Shiny listens on container `0.0.0.0:3838`; rootless Podman publishes only VM `127.0.0.1:3838`; UFRC reverse-proxies it |
| Development backend | HPG-restricted staging service on `127.0.0.1:3839`, if UFRC permits it on the same allocation |
| Database requirement | No managed database and no database network service |
| Runtime writes | Temporary files and operational logs only; scientific data are immutable and read-only |
| External AI/API traffic | None |

The base allocation is the right starting size: the current TIFF collection is
approximately 109 GiB and the broader `data/` directory is approximately 1.9 GiB.
The requested 1 TiB persistent allocation provides ample room for a curated release,
staging, one rollback release, logs, and temporary transfer files. The active data and
image payload is approximately 111 GiB. A 120 GiB working estimate assumes that unchanged
TIFFs are shared between application releases; a full duplicate image release would be
approximately 220 GiB and still fits the base allocation. UFRC's price page listed the
one-year base at $300 when this handoff was prepared, but price and availability must be
confirmed in the purchase form.

## Information the project owner must supply

Complete this table before filling the DOCX or opening the provisioning ticket.

| Required information | Value |
|---|---|
| HPG group | `maigan` |
| PI/HPG sponsor full name | Maigan Brusko, PhD |
| PI GatorLink and email | `[PI GATORLINK / EMAIL]` |
| Primary technical contact | Justin Smith, PhD — `[GATORLINK / EMAIL / PHONE]` |
| Backup technical contact | `[NAME / GATORLINK / EMAIL]` |
| Maintainer HPG public SSH keys and identifying comments | `[PRIMARY KEY / BACKUP KEY]` |
| Department/unit | `[DEPARTMENT]` |
| Fiscal contact | `[NAME / EMAIL]` |
| Purchase ChartFields | `[CHARTFIELDS]` |
| PubApps allocation or order number | `RC83797` — trial `PA-INSTANCE` allocation |
| Current allocation dates | 2026-08-31 through 2026-11-30 |
| Target staging date | `[YYYY-MM-DD]` |
| Target reviewer go-live date | `[YYYY-MM-DD]` |
| Expected end of first review round | `[YYYY-MM-DD]` |
| Journal and manuscript identifier | `[JOURNAL / MANUSCRIPT ID WHEN AVAILABLE]` |
| Data owner | `[NAME / ROLE]` |
| Applicable IRB, DUA, or MTA identifiers | `[IDENTIFIERS OR NOT APPLICABLE]` |
| UFRC support ticket | `[TICKET NUMBER WHEN ISSUED]` |
| Risk-assessment number, if required | `[NUMBER OR UFRC CONFIRMED NOT SEPARATELY REQUIRED]` |
| Confirmed reviewer/concurrency estimate | 1–8 total / up to 5 concurrent |
| Required post-publication retention period | `[JOURNAL/SPONSOR REQUIREMENT]` |
| Authoritative HPG source path/custodian | `[PATH / CUSTODIAN]` |

## External-owner checklist

These actions happen outside the repository and cannot be completed by the application
code.

### 1. Security and authority

- [x] Navigator key revoked (owner confirmation received 2026-08-31). Do not provision
  a replacement production key because AI is disabled.
- [ ] Confirm that the PI, primary maintainer, and backup maintainer have active
  HiPerGator accounts and can use GatorLink/MFA and SSH from HPG.
- [ ] Collect each maintainer's HPG public SSH key with an identifying comment for the
  instance access list; never send a private key.
- [ ] Give the PI/data owner, IRB/privacy contact where applicable, and the nPOD
  agreement owner an exact inventory of the proposed reviewer files and fields.
- [ ] Obtain written confirmation that the reviewer release contains no PHI or other
  Restricted data and that the governing IRB, DUA/MTA, consent, and nPOD terms permit
  delivery to journal editors and reviewers.
- [ ] Obtain separate written approval before making the same assets unauthenticated or
  publicly downloadable after publication.
- [ ] Ask the journal whether UF/nPOD branding, acknowledgments, the UF hostname, or
  case identifiers conflict with double-blind review.

The authorization inventory must cover the 15 donor OME-TIFF images and sidecars,
case identifiers, age/sex and clinical attributes used by the app, cell- and islet-level
derived measurements, spatial boundaries, downloadable data, proxy authentication logs,
and application logs that can contain case/islet selections or error details.
An open-source software license does not grant redistribution rights to the underlying
research material.

If an agreement owner determines that any included asset is Restricted or cannot be
shown to reviewers, stop the PubApps request. Remove the asset or ask UFRC whether the
application must move to PrivApps.

### 2. Account and resource purchase

- [x] Trial resource allocated to HPG group `maigan`: one `PA-Instance`, reference
  `RC83797`, active 2026-08-31 through 2026-11-30.
- [ ] Request any missing HiPerGator accounts through
  <https://it.ufl.edu/rc/get-started/request-hipergator-account/>.
- [ ] Have the fiscal contact submit UFRC's Internal Service Purchase Request through
  <https://it.ufl.edu/rc/get-started/purchase-request/>.
- [ ] Purchase one annual base PubApps allocation for **Islet Explorer** under the HPG
  group above: select **Basic PubApps Instance (`PA-INSTANCE`)** on the Internal
  HiPerGator Service order form. Request no GPU, additional storage, or additional compute
  initially.
- [ ] Allow for the stated five-business-day purchase change period and Business Center
  purchase-order coordination.
- [ ] Record the order/allocation number, activation date, annual renewal date, and
  fiscal owner. Create renewal reminders 90 and 30 days before expiration.
- [ ] Obtain a written UFRC quote before later adding resources; UFRC's public pages
  currently present inconsistent add-on pricing.

The trial allocation is sufficient to request and validate the instance now. Do not rely
on its current expiration for journal review: arrange continuation before 2026-11-30 and
retain service through editorial review, publication, and the agreed post-publication
period. Because the 90-day point is effectively immediate, start that continuation
decision with the setup request; set a second reminder for 2026-10-31.

### 3. Form, ticket, and approvals

- [ ] Fill `UFIT_PubApps_Instance_Request.docx` with the copy-ready content below.
- [ ] Read the supplied project plan, PubApps policy, and web application development
  standard before signing the form's receipt/compliance acknowledgment.
- [ ] Attach `docs/pubapps_architecture.svg` (or its PNG rendering) and the written
  data-classification and redistribution determination.
- [ ] Have the PI/HPG sponsor sign and date the policy acknowledgment.
- [ ] Open a UFRC support ticket at <https://support.rc.ufl.edu/> using the ticket text
  below. Attach the signed DOCX and purchase/allocation number.
- [ ] Complete an individual risk assessment only if UFRC directs the project to do so.
  Follow <https://docs.rc.ufl.edu/services/web_hosting/risk_assessment/> and return the
  assessment number in the support ticket.
- [ ] Obtain unit accessibility and UF branding reviews before production approval.
- [ ] Ask UFRC to perform its required security review and vulnerability scan, and
  retain the approval and scan results with the release record.
- [ ] Confirm the planning assumptions with the editor/PI: total reviewers, maximum
  concurrent sessions, two named maintainers, and the required retention period.

### 4. Journal and reviewer operations

- [ ] Ask the editor how credentials must be delivered. Prefer the journal's
  confidential reviewer-instructions field rather than the manuscript, repository, or
  a public support-ticket comment.
- [ ] Test the final URL and credential from a clean browser on a non-UF network.
- [ ] Give reviewers the short handoff described below only after production acceptance
  is complete.
- [ ] Review the authorized reviewer roster and maintainer SSH access quarterly while
  the Sensitive deployment is active; retain the dated review record.
- [ ] Rotate the shared password after each review round and immediately after suspected
  disclosure. Revoke it when editorial review ends unless continued access is approved.
- [ ] If UFRC cannot support shared Basic Auth for Sensitive-data access review and
  brute-force controls, stop go-live and agree on named Basic Auth accounts or another
  UFRC-supported identity mechanism.

## Copy-ready PubApps form content

### Project Name

Panc Floyd (Islet Explorer reviewer application)

### The HiPerGator group responsible for this project

`maigan`

### The primary technical contact (developer) for the project

Justin Smith, PhD — `[GATORLINK, EMAIL, PHONE]`

Backup technical contact: `[BACKUP CONTACT NAME, GATORLINK, EMAIL]`

### The principal responsible investigator for this project

Maigan Brusko, PhD — `[GATORLINK, EMAIL, DEPARTMENT]`

### A detailed description of the project

> Islet Explorer is an R Shiny research application through which confidential journal
> editors and reviewers explore de-identified and derived pancreatic-islet data and
> browser-streamed OME-TIFF multiplexed images from 15 nPOD donors. It provides marker
> and cell-population plots, composition-based trajectory/PAGA analysis, spatial
> analysis, and image/segmentation drill-down. The application is non-mutating: it
> accepts no user uploads and cannot alter its canonical scientific data. Reviewers may
> generate filtered CSV exports, which are created in session-scoped temporary storage
> and delivered through the browser. Only temporary files and operational logs are
> writable. The service has no managed or network-accessible database, submits no Slurm
> jobs, and has AI and external API integrations disabled.

### The objective of the web application

> Provide journal editors and reviewers a stable, reproducible interface for evaluating
> the figures, analysis methods, and supporting data used in the Islet Explorer
> manuscript. The exact application image and data release reviewed by the journal will
> remain frozen and recoverable through the editorial process and publication.

### The expected user community / number of users

> Approximately 1–8 journal editors and reviewers in total, with no more than five
> simultaneous sessions, plus two named UF technical maintainers. Reviewer access will
> use a shared, temporary Basic Auth credential at the UFRC reverse proxy. Maintainer
> access will use GatorLink/MFA, HiPerGator, and SSH.

Replace the estimate if the journal provides a materially different reviewer count.

### Map of external components, connections, and data

```text
Journal editors/reviewers
        | HTTPS + whole-site shared Basic Auth
        v
UFRC high-availability reverse proxy
        | TLS termination, WebSocket upgrade, HTTP byte-range passthrough
        v
PubApps VM 127.0.0.1:3838
        | rootless Podman forwarding
        v
Islet Explorer/Shiny container 0.0.0.0:3838
        |-- read-only versioned scientific release
        |     |-- Parquet queried by embedded DuckDB (no service/port)
        |     |-- indexed GeoPackage segmentation files
        |     |-- per-islet CSV drill-down files
        |     `-- H5AD scientific/provenance artifacts
        |-- read-only OME-TIFF image store
        `-- writable bounded temporary files and operational logs only

Temporary CSV/PNG response: container -> proxy -> reviewer browser
Operational logs: container/proxy -> UFRC-approved log collection and retention

UF maintainers -- GatorLink/MFA -> HiPerGator login -> SSH -> PubApps VM

Release ingress only:
Authoritative lab/HPG data -> checksum manifest -> transfer through HiPerGator
Tested OCI archive -> checksum verification -> transfer through HiPerGator -> Podman import

Optional staging, only if UFRC approves it:
HPG-restricted proxy route -> PubApps VM 127.0.0.1:3839 -> staging container
```

There is no inbound database port, no external analytics service, no user-upload path,
no email service, and no production AI/API call.

### UF data classification and data elements being processed

The deployed collection must be handled at its highest applicable classification.

| Classification | Answer | Explanation |
|---|---|---|
| Restricted | **No, conditional on written approval** | No PHI, direct identifiers, FERPA, export-controlled data, or contractually Restricted elements will be deployed. No API, database, or application credentials are present in the image/data release. The reviewer Basic Auth secret is held and processed only by the UFRC proxy. The PI/data owner and applicable agreement owners must confirm this before submission. |
| Sensitive | **Yes** | De-identified donor/patient attributes, multiplexed images, derived cell/islet measurements, unpublished research results, and minimal authentication/access-log metadata are used during peer review. |
| Open | **Yes, for selected elements only** | Source code, approved documentation, phenotype definitions, and already released descriptions may be Open. This does not lower the classification of the combined reviewer deployment. |

Example Sensitive elements:

> De-identified nPOD donor case identifiers and approved demographic/clinical fields;
> OME-TIFF pancreatic tissue images; cell, islet, trajectory, phenotype, and spatial
> measurements derived from those images; unpublished manuscript results; and minimal
> reverse-proxy/application access logs.

If the form asks, "If the app processes Restricted data, what type of restriction is in
place?", enter:

> Not applicable; none, contingent on the attached PI/data-owner and agreement-owner
> confirmation.

Example Open elements:

> MIT-licensed application source code, approved public documentation, and approved
> phenotype marker definitions. These Open elements do not lower the classification of
> the combined reviewer deployment.

Faculty PI/HPG Group Sponsor Name: Maigan Brusko, PhD

Date: 2026-08-31

## Copy-ready UFRC support ticket

**Subject:** `PubApps setup request — Panc Floyd / Islet Explorer — maigan — RC83797`

**Body:**

> Please provision and review a new Public Restricted PubApps application for the Panc
> Floyd project, deployed as Islet Explorer, under HPG group `maigan`. Trial resource:
> one `PA-Instance`, reference `RC83797`, active 2026-08-31 through 2026-11-30. The
> completed UFIT-RC Web Hosting Instance Request and data-flow/classification materials
> are attached. Desired staging date: `[DATE]`; reviewer production date: `[DATE]`.
>
> The requested starting allocation is the standard PubApps base: 2 vCPU, 16 GiB RAM,
> 1 TiB `/pubapps`, and standard local container storage. The curated deployment is
> expected to use approximately 111 GiB for its active image/data payload. The 1 TiB
> allocation also provides headroom for staging, rollback, logs, and transfer files; a
> 120 GiB working estimate assumes unchanged TIFFs are shared between releases. No GPU,
> additional storage, Slurm access,
> database service, or HPG filesystem mount is requested.
>
> The application is a non-mutating R Shiny service in a rootless Podman OCI container,
> managed by user-mode Quadlet/systemd. It creates requested CSV downloads only in
> session-scoped temporary storage; it cannot modify canonical scientific data. Please
> confirm Podman/Quadlet as the supported
> runtime, the assigned VM/service user and UID/GID, `/pubapps` project path, local
> container path, quotas, SSH route, and login-linger/reboot behavior.
>
> Requested production mapping:
> `https://islet-explorer.rc.ufl.edu/islet-explorer/` to
> `127.0.0.1:3838/islet-explorer/`. If supported on the same allocation, please also
> provide an HPG-restricted staging mapping to `127.0.0.1:3839/islet-explorer/`.
>
> Please protect the entire production site with shared Basic Auth, including the app
> HTML, Shiny WebSocket, downloads, OME-TIFF, XML, JSON, and other static assets. Please
> preserve WebSocket Upgrade connections and long-lived Shiny sessions; preserve HTTP
> Range requests and `206 Content-Range` responses; and avoid buffering or caching the
> complete multi-gigabyte TIFF files. Please provide TLS and redirect HTTP to HTTPS.
>
> The application uses immutable read-only scientific assets. It accepts no uploads and
> performs no scientific-data writes. AI, Statistics, and Inflammation are disabled. It
> makes no production AI/API calls. The tested data-free OCI archive will be transferred
> through HiPerGator, checksum-verified, and imported into rootless Podman. Production has
> no registry or other third-party runtime dependency.
>
> Database statement: no managed PostgreSQL, MySQL, or other database service is needed.
> The release may use an embedded DuckDB library to issue read-only analytical queries
> over immutable local Parquet files. This creates no listener, credentials, shared
> transaction log, or mutable database state. OME-TIFF images remain ordinary files for
> HTTP byte-range delivery, and H5AD remains a versioned scientific provenance artifact.
>
> Data are treated as Sensitive during peer review. The project has classified the
> release as not Restricted, contingent on the attached PI/data-owner and agreement
> approvals. Expected use is 1–8 reviewers total, no more than five concurrent sessions,
> and two UF maintainers.
>
> Maintainer GatorLink usernames: `[PRIMARY]`, `[BACKUP]`.
>
> Please advise whether a separate individual risk assessment is required, where
> application and proxy logs must be written and how long they must be retained, what
> backup/snapshot coverage UFRC provides, whether UFRC or the project should create and
> rotate the Basic Auth credential, and whether UFRC provides failed-login rate limiting
> or lockout for this shared credential. Please also provide the required vulnerability,
> accessibility, and production-approval steps.

Attach to the ticket:

1. Signed `UFIT_PubApps_Instance_Request.docx`.
2. The finalized architecture/data-flow map.
3. PI/data-owner classification and redistribution approval.
4. nPOD/IRB/DUA/MTA determination or identifier, as applicable.
5. Allocation/purchase confirmation.
6. Risk-assessment number if UFRC required a separate assessment.

## Database and storage position

The reviewer deployment does **not** benefit from a managed PostgreSQL service. The
workload is immutable, analytical, small-user-count, and local to one PubApps instance;
PostgreSQL would add a network service, credentials, patching, backup, and availability
obligations without providing a needed multi-writer capability.

The production release uses a hybrid, embedded layout:

- Keep OME-TIFF images as files; storing approximately 109 GiB of imagery as database
  BLOBs would break efficient HTTP Range delivery and complicate recovery.
- Keep H5AD as the canonical scientific/provenance artifact. Flatten only the fields the
  app needs into versioned runtime tables after output-parity tests.
- Store the largest tabular hot paths as immutable Parquet: all-cell distance data,
  signed-distance KDE rows, and donor tissue cells.
- Query those Parquet files with an embedded DuckDB connection created once per Shiny
  worker. Mount the files read-only, expose no database port, and use parameter binding
  plus an allowlist for selectable marker columns.
- Keep small prepared islet summaries in process-global memory rather than rewriting all
  Plot/Trajectory logic as SQL.
- Keep the 5,214 per-islet drill-down CSVs as direct-lookup files. Measured point lookup
  was faster than consolidating them into DuckDB; package them as one archive for release
  transfer, then extract them into the immutable versioned data directory.
- Package each donor's segmentation boundaries as an indexed, read-only GeoPackage and
  query only the selected islet's bounding box. Do not use SQLite as the primary
  analytical store.

Provisioning and purchasing may proceed while the release is being built, but the
manifest-backed runtime is a production prerequisite. Before promotion it must reproduce
every legacy downloadable table and plotted statistic, including row counts and identifiers,
with exact categorical values and documented numeric tolerance. If parity is not complete,
do not promote that release; keep the reviewer site closed until a validated release exists.

## Information UFRC must return or confirm

Record these answers in the ticket or release record:

- VM hostname and SSH route from HPG.
- Service-account name and UID/GID.
- Persistent and local container paths and quotas.
- Approved container runtime and Quadlet/systemd conventions.
- Production and staging URLs, backend ports, and final proxy path.
- WebSocket and TIFF HTTP Range proxy behavior.
- Basic Auth ownership, credential-change procedure, rate limiting, and lockout behavior.
- Login linger, service restart, and VM reboot behavior.
- Log path, format, collection, access, and retention requirements.
- Backup/snapshot scope and restoration process.
- Required risk-assessment number or confirmation that the form is attached to UFRC's
  main PubApps assessment.
- Security scan, SBOM, SAST/container-scan evidence, accessibility review, and production
  approval process.
- Support escalation contact and expected response route during journal review.

## Release and go-live evidence

Maintain one release record containing:

- Git release tag and commit.
- OCI image digest and prior rollback digest.
- SHA-256 data/image manifest with file counts and sizes.
- Reconciled islet, phenotype, donor, and image counts used in the manuscript and app.
- Resolved R package/version manifest, pinned Python requirements, base/image digests,
  and SBOM.
- Lint reports with zero errors and as few documented warnings as practical, plus an
  input-validation attestation for every browser-controlled parameter.
- Secret, static-analysis, dependency, and container-scan reports with no unresolved
  Critical or High findings.
- A reviewed log inventory showing one event per line, no credentials, no unnecessary
  Sensitive fields, the approved UFRC collection path, and retention.
- Automated and manual WCAG 2.1 AA review results.
- UFRC production approval and PI attestation.
- Five-concurrent-session load-test result.
- Container restart, VM reboot, health-check, log-rotation, and rollback evidence.
- An off-campus test showing Basic Auth protects the app, WebSocket, downloads, and every
  image/sidecar path.
- A TIFF test showing a small byte-range request returns `206 Partial Content` with the
  correct `Content-Range` without transferring the entire image.

Freeze the accepted Git tag, image digest, data manifest, and feature flags throughout a
review round. Do not run `git pull`, rebuild from a moving branch, or automatically update
production during review.

## Reviewer handoff template

Send the password separately through the journal's confidential mechanism.

> **Islet Explorer reviewer site**
>
> URL: `https://islet-explorer.rc.ufl.edu/islet-explorer/`
>
> Username: `[NEUTRAL REVIEWER USERNAME]`
>
> Password: supplied separately through the confidential reviewer channel.
>
> The reviewer release contains the Plot, Trajectory, Viewer, and Spatial tabs. A first
> session or trajectory recalculation may take longer than an already warmed interaction.
> Please use a current desktop version of Firefox, Chrome, Edge, or Safari. If access or
> an analysis fails, contact `[SUPPORT NAME / EMAIL]` and include the time, browser, active
> tab, and visible error message; do not include the password.

Use one neutral account name that does not reveal reviewer or author identity. Never put
the credential in the manuscript, repository, container image, data manifest, analytics,
or application logs.

## Ongoing and post-review operations

- Monitor HTTPS/readiness, container health, disk consumption, memory, and repeated 5xx
  responses. Route alerts to both maintainers.
- Apply a security fix through a newly tested immutable image digest and document any
  reviewer-visible change. Retain the previous known-good digest for rollback.
- Keep HPG as the authoritative scientific-data backup. Confirm whether UFRC snapshots
  cover only operational recovery and document a checksum-verified re-sync procedure.
- Retain the accepted reviewer release through the editorial decision and for the period
  confirmed by the journal, PI, agreement owner, and sponsor. Use 12 months after
  publication only as the planning default until that requirement is confirmed.
- At the end of peer review, rotate or revoke the shared reviewer credential and document
  the action.
- Remove Basic Auth and request a Public Open configuration only after the PI/data owner,
  nPOD/agreement owner, and UFRC approve the exact post-publication assets. Otherwise keep
  authenticated access or migrate to the service appropriate to the final classification.

## Reference links

- UFRC web hosting overview: <https://docs.rc.ufl.edu/services/web_hosting/>
- UFRC current price list and PubApps service item: <https://it.ufl.edu/rc/get-started/price-list/>
- UFRC deployment request details: <https://docs.rc.ufl.edu/services/web_hosting/deployment/>
- UFRC instance request guidance: <https://docs.rc.ufl.edu/services/web_hosting/instance_request/>
- UFRC systemd/Quadlet guidance: <https://docs.rc.ufl.edu/services/web_hosting/systemd/>
- UFRC Shiny guidance: <https://docs.rc.ufl.edu/services/web_hosting/shiny/>
- UFRC web application standards: <https://docs.rc.ufl.edu/services/web_hosting/standards/>
- UF Data Classification Policy: <https://policy.ufl.edu/policy/data-classification-policy/>
- UF accessibility policy: <https://policy.ufl.edu/policy/electronic-information-technology-and-communication-accessibility-policy/>
- UF Brand and User Experience Policy: <https://policy.ufl.edu/policy/brand-and-user-experience/>
- nPOD publication policies: <https://npod.org/publications/policies/>
