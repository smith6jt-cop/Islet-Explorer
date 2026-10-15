# Panc Floyd / Islet Explorer UFRC PubApps Request and Reviewer-Hosting Handoff

Local operating update, 2026-09-28: the development app again serves the
original 15 donors and 5,214 islets at
`http://10.15.152.7:8080/islet-explorer/` (nginx to loopback port 3840).
Its `/local_images/` route still serves an earlier 21-image collection; the
owner permits that mismatch only during local development. It is not an
approved expanded app or a production configuration. Six donor additions are
staged separately as raw-only evidence, with `app_ready: false` and
`scientific_acceptance: false`. See [the current handoff](docs/HANDOFF_2026-09-28.md)
and [phase-2 contract](docs/donor_expansion_phase2.md). Earlier regenerated
21-donor schema-2/3/4 candidates remain historical experiments. Institutional
status below remains as last verified on 2026-09-16; local checks do not
establish a new UFRC acceptance.

This document is the actionable companion to
`UFIT_PubApps_Instance_Request.docx`. It provides copy-ready form and support-ticket
text, identifies information that the project owner must supply, and records the
operational decisions for the journal-review deployment. The technical deployment
procedure is in [`deploy/README.md`](deploy/README.md), and the detailed project log
inventory is in [`docs/pubapps_log_inventory.md`](docs/pubapps_log_inventory.md).
Do not submit the request with bracketed placeholders left unresolved.

Status updated 2026-09-16: the existing records are UFRC request `83810`,
configuration reference `SYSCONFIG-8173`, and individual assessment `IRM-456246`
(submission confirmed 2026-09-01). These references were recovered through
read-only UF Outlook search. On 2026-09-16 an assigned ISM requested a brief
Microsoft Teams meeting for the risk review; scheduling and completion are not
established. UFRC comment 8, received 2026-09-11, had reported its ISO waiting
for that local review. Comment 7 requested that the UFRC-designated technical
contact be added to the assessment; completion is not established. The owner
previously reported a provisioned VM, but assigned host/account/SSH/proxy details
remain unconfirmed.
Continue the existing records; do not open duplicate requests. Neither allocation
nor an open assessment establishes production approval. The current requirement-to-evidence index is
[`docs/ufrc_readiness.md`](docs/ufrc_readiness.md).

The checklists below distinguish UFRC/policy requirements from conditional and
recurring owner actions and from project recommendations. Project hardening in the
deployment runbook remains part of this application's release process, but each such
control is not a separate UFRC approval gate.

Expansion work began on 2026-09-16 with a target of
21 donor images: the original 15 nPOD donors plus 6090, 6232, 6539, the
`PAN-6558-AAb_test35_35p-LED.tif` acquisition, 6520-01A, and HPAP-078. The alternate
`CX_PAN-6558` image and 6520-01B are excluded. HPAP-078 retains its HPAP provenance;
do not describe the target as 21 nPOD donors. Image preparation and experimental
model fitting took place offline, outside PubApps. The current phase-2 contract
requires preserving the original 15 donors exactly and adding only independently
validated new-donor measurements. Preserve the frozen `5809be4` release and its
dated 15-image evidence. **The 6090/6539 source mix-up is resolved:** 6090 matches all
36 original raw channel planes; the distinct 6539 stack was reconstructed from
its own raw acquisition and checked pixel-for-pixel. The corrected six viewer
TIFFs total 46,877,746,573 bytes, but no expanded image release is approved.
Preserve the rejected duplicate and its old annotation separately. Source
preparation and corrected paths remain in [the expansion status](docs/cohort_expansion_status.md).

The offline raw-only bundle at
`data/development/donor-additions-20260923-raw-v1/` contains 943,068 new-donor
objects: 942,867 fully measured and 201 wholly unmeasured. Its 852 zero-tolerance
checks preserve the original donor artifacts, but it has no accepted new-cell
phenotypes, normalized expression, scVI fields, islet composition, compatible
trajectory, optional loader, or combined islet count. HPAP-078 boundaries are
owner-approved; 6090, 6232, 6520, recovered 6539, and 6558 still require
donor-specific visual review. Update the existing request and assessment with
an approved final inventory and applicable HPAP disclosure/usage terms through
the authorized owner; this document does not infer an institutional disposition.

The earlier schema-2/3/4 all-donor runtimes, phenotype-retirement analyses,
and raw-compartment diffusion axes are experimental evidence, not the current
app or an accepted additions-only method. The September 16 schema-2 candidate
reported 8,371 islets and 3,499 eligible trajectory islets; these are not
preservation-compatible expansion counts. See [the historical delivery record](docs/all_donor_delivery.md) and
[September 17 revision](docs/trajectory_revision_20260917.md).

## Deployment decisions

| Item | Decision |
|---|---|
| Service | Assigned ISM requested a risk-review Teams meeting on 2026-09-16; assigned host details unconfirmed. Trial `PA-INSTANCE` allocation for `maigan`, reference `RC83797` |
| Current allocation term | 2026-08-31 through 2026-11-30; renew or purchase before reviewer access extends beyond that date |
| Initial resources | Requested base: 2 vCPU, 16 GiB RAM, 1 TiB `/pubapps`, standard local container storage; verify assigned resources |
| Add-ons | None initially; observe representative use at the expected five-session peak before requesting more capacity |
| Access during peer review | **Public Restricted**, with whole-site shared Basic Auth |
| UF data classification during peer review | **Sensitive** |
| Expected audience | 1–8 editors/reviewers total, up to 5 concurrent sessions, plus 2 named UF maintainers |
| Reviewer-visible features | Plot, Trajectory, Viewer, Spatial in the planned reviewer profile; the exact production release still needs acceptance |
| Disabled features | AI, Statistics, and Inflammation |
| Production URL request | `https://islet-explorer.rc.ufl.edu/islet-explorer/` |
| Production backend request | Shiny listens on container `0.0.0.0:3838`; the project requests rootless Podman publication on VM `127.0.0.1:3838`, subject to UFRC confirming the reachable proxy target |
| Database requirement | No managed database and no database network service |
| Runtime writes | Temporary files and operational logs only; scientific data are immutable and read-only |
| External AI/API traffic | None |

The frozen 15-image release inventory measures 117,092,780,533 bytes (109.05 GiB),
including 116,225,236,881 bytes of images and sidecars. The corrected six
experimental TIFFs total 46,877,746,573 bytes; they are not part of an
accepted expanded release. Earlier planning estimated 151.90 GiB for all 21
images and a 260.95 GiB lower bound when retaining the frozen bundle. The
regenerated candidates reported 164,130,834,281 bytes for schema 2 and
164,147,382,976 bytes for schema 3; the associated frozen-plus-candidate-staging
estimates were 414.769 and 414.799 GiB. These dated sizes are historical planning
evidence, not current expanded-release capacity acceptance.
Full deployment copies require their own bytes even where local development
used hardlinks. UFRC's price page listed the one-year base at $300 when this
handoff was prepared; price and availability must be confirmed in the purchase form.
One earlier candidate preview R worker used approximately 1 GiB RSS. Five-session memory/latency
and assigned-host measurements remain pending; local availability does not
establish capacity or deployment acceptance.

## Information the project owner must supply

Complete this table in the existing restricted ticket and operating record before
finalizing the DOCX and production review. Details have been requested from the
project owner; unknown values below are not new work requests.

| Required information | Value |
|---|---|
| HPG group | `maigan` |
| PI/HPG sponsor full name | Maigan Brusko, PhD |
| PI GatorLink and email | `[PI GATORLINK / EMAIL]` |
| Primary technical contact | Justin Smith, PhD — `[GATORLINK / EMAIL / PHONE]` |
| Backup technical contact | `[NAME / GATORLINK / EMAIL]` |
| Maintainer HPG public SSH keys and identifying comments | `[RESTRICTED TICKET REFERENCE; DO NOT PASTE KEYS HERE]` |
| Department/unit | `[DEPARTMENT]` |
| Fiscal contact | `[NAME / EMAIL]` |
| Purchase ChartFields | `[CHARTFIELDS]` |
| PubApps allocation or order number | `RC83797` — trial `PA-INSTANCE` allocation |
| Current allocation dates | 2026-08-31 through 2026-11-30 |
| Target reviewer go-live date | `[YYYY-MM-DD]` |
| Expected end of first review round | `[YYYY-MM-DD]` |
| Journal and manuscript identifier | `[JOURNAL / MANUSCRIPT ID WHEN AVAILABLE]` |
| Data owner | `[NAME / ROLE]` |
| Applicable IRB, DUA, or MTA identifiers | `[IDENTIFIERS OR NOT APPLICABLE]` |
| UFRC support ticket | `83810`, configuration reference `SYSCONFIG-8173` |
| Individual Risk Assessment | `IRM-456246`; assigned ISM requested a Teams meeting on 2026-09-16; review and requested assessment contact still unconfirmed |
| Provisioned VM, service account, and SSH route | Owner previously reported provisioning; `[ASSIGNED HOST / ACCOUNT / SSH ROUTE]` pending |
| Production URL, proxy target, and backend port | `[UFRC-CONFIRMED VALUES]` |
| `/pubapps` and local container paths/quotas | `[UFRC-CONFIRMED VALUES]` |
| Confirmed reviewer/concurrency estimate | 1–8 total / up to 5 concurrent |
| Required post-publication retention period | `[JOURNAL/SPONSOR REQUIREMENT]` |
| Authoritative HPG source path/custodian | `[PATH / CUSTODIAN]` |

## Required provisioning and release checklist

These are the external or owner-controlled requirements for the selected Sensitive,
Public Restricted deployment. A checked item should point to dated evidence in the
restricted support ticket or release record; do not put personal contact information,
SSH keys, credentials, or sensitive approvals in this repository.

### Complete the existing request and assessment

- [x] Trial `PA-INSTANCE` allocation recorded for HPG group `maigan`, reference
  `RC83797`, active 2026-08-31 through 2026-11-30.
- [x] Existing UFRC request `83810`, configuration reference `SYSCONFIG-8173`,
  and assessment `IRM-456246` identified in UF Outlook. Comment 8 of 2026-09-11
  reports local ISM review pending; this is not an approval record.
- [ ] Arrange the brief Teams meeting requested by the assigned ISM on 2026-09-16,
  complete the local risk-assessment review identified in comment 8, and
  confirm the UFRC-designated assessment technical contact requested in comment 7.
  Keep personal contact details in the restricted assessment/ticket.
- [ ] Give the PI/data owner, applicable IRB/privacy contact, and nPOD or other agreement
  owner an exact inventory of the proposed reviewer files and fields. Obtain written
  confirmation that the release contains no PHI or other Restricted data and that its
  delivery to journal editors and reviewers is permitted.
- [ ] Have the PI and technical owner verify the already populated
  `UFIT_PubApps_Instance_Request.docx`, confirm that its expected-user estimate remains
  accurate, read the supplied project plan/policy/development standard, and complete the
  named-and-dated compliance acknowledgment.
- [ ] Confirm that the PI and two named maintainers have active HiPerGator accounts.
  Supply their HPG usernames and public SSH keys through the restricted ticket; never
  send a private key or commit a public key to this repository.
- [ ] Attach the reviewed form, allocation evidence, and requested technical
  evidence to existing request `83810` and assessment `IRM-456246` through the
  authorized owner. Resolve assessment follow-up items and
  retain its disposition; opening an assessment is not its acceptance.

### Confirm the provisioned instance configuration

- [ ] Reconcile the owner-reported provisioning status with the current assessment
  status in UFRC comment 8 of 2026-09-11, then record the assigned VM/service account, storage/container
  paths/quotas, SSH route, and proxy/backend assignment.
- [ ] Have both named maintainers independently confirm that they can reach the service
  account by SSH from their own HPG accounts; retain the dated acceptance in the
  restricted ticket.

The authorization inventory must cover every donor OME-TIFF image and sidecar
(15 in the frozen release; 21 targeted in the expanded candidate),
case identifiers, age/sex and clinical attributes used by the app, cell- and islet-level
derived measurements, spatial boundaries, downloadable data, proxy authentication logs,
and application logs that can contain case/islet selections or error details.
An open-source software license does not grant redistribution rights to the underlying
research material.

If an agreement owner determines that any included asset is Restricted or cannot be
shown to reviewers, keep reviewer access closed. Remove the asset or ask UFRC whether the
application must move to PrivApps.

### Before reviewer access

- [x] Navigator key revoked (owner confirmation received 2026-08-31). Do not provision
  a replacement production key because AI is disabled.
- [ ] Document UF accessibility and branding compliance, application input validation,
  required lint results, and appropriate one-event-per-line application logging.
- [ ] Build the exact release candidate and retain its SBOM and project security results.
  The project release policy requires the Syft/Grype image gate to have no High or
  Critical findings; provide that evidence to UFRC for review.
- [ ] Obtain UFRC's production security/availability acceptance and its confirmation of
  the ongoing UFIT network-scan process. UFRC review does not replace the project scan,
  and the project scan does not replace UFRC acceptance.
- [ ] Test the confirmed HTTPS URL and neutral reviewer credential from a clean non-UF
  network. Whole-site authentication must cover the app, WebSocket, downloads, TIFFs,
  and image sidecars; TIFF byte ranges must work through the proxy.
- [ ] Obtain PI/data-owner go-live approval, agree with the editor on a confidential
  credential-delivery channel, and give reviewers the handoff below only after the other
  required checks pass.

### Conditional and recurring obligations

- [ ] Arrange continued or annual PubApps service before reviewer access extends beyond
  2026-11-30. Record the fiscal owner, activation/renewal dates, and reminders in the
  restricted operating record.
- [ ] If UFRC confirms that project-plan item 6 applies to this Sensitive deployment,
  review authorized reviewer access and maintainer SSH access quarterly while the
  deployment remains active and retain the dated review record.
- [ ] Rotate the shared password after each review round and immediately after suspected
  disclosure. Revoke it when editorial review ends unless continued access is approved.
- [ ] If an ongoing UFIT network scan reports an exposed High or Critical finding,
  remediate it on the UFIT-directed timeline and retain the rescan/clearance record.
- [ ] Obtain separate PI/data-owner, agreement-owner, and UFRC approval before making any
  scientific assets unauthenticated or publicly downloadable after publication.

Project recommendations, not separate UFRC gates: ask the journal whether branding,
the UF hostname, or case identifiers conflict with double-blind review; exercise a
representative five-session smoke test before requesting additional resources; and
obtain a current UFRC quote before adding capacity.

The trial allocation supports instance validation once setup is confirmed. Do not rely
on its current expiration for journal review: arrange continuation before 2026-11-30 and
retain service through editorial review, publication, and the agreed post-publication
period.

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
> browser-streamed OME-TIFF multiplexed images. The frozen reference release contains
> 15 nPOD donors. Six additional donors, including HPAP-078 with separate HPAP
> provenance and applicable usage terms, are staged as raw-only offline evidence;
> no expanded reviewer release is accepted. The application provides marker
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
Requested backend: PubApps VM 127.0.0.1:3838
        | actual proxy target pending UFRC confirmation
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

Date: `[ACTUAL PI ATTESTATION DATE]` — this draft was updated 2026-09-16 and does
not record a completed attestation.

## Copy-ready follow-up for the existing UFRC ticket

**Subject:** `Request 83810: Islet Explorer production review — IRM-456246`

**Body:**

> Please continue the production review for the Panc Floyd project, deployed as
> Islet Explorer under HPG group `maigan`, in existing ticket `83810`, configuration
> reference `SYSCONFIG-8173`, and individual Risk Assessment `IRM-456246`.
> The assigned ISM requested a brief risk-review meeting on 2026-09-16, following
> comment 8's report that the ISO was waiting for local review. Please confirm
> that review's status and whether the assessment technical contact requested in
> comment 7 has been added, then return the current instance setup status.
> `RC83797` records one trial `PA-INSTANCE` allocation
> active 2026-08-31 through 2026-11-30. The reviewed UFIT-RC Web Hosting Instance
> Request and current technical evidence are attached. Desired reviewer production
> date: `[DATE]`.
>
> The requested starting allocation is the standard PubApps base: 2 vCPU, 16 GiB RAM,
> 1 TiB `/pubapps`, and standard local container storage. The frozen 15-image
> runtime bundle measures 109.05 GiB. Six new donors are staged separately as
> raw-only research evidence, not as an app-ready or deployable 21-image release.
> Any future expansion requires a reviewed image/data inventory and new capacity
> measurements before deployment.
> Segmentation and model fitting are offline; no GPU,
> additional storage, Slurm access, database service, or HPG filesystem mount is
> requested.
>
> The application is a non-mutating R Shiny service in a rootless Podman OCI container,
> managed by project-supplied user-mode Quadlet/systemd configuration. It creates CSV
> downloads only in session-scoped temporary storage and cannot modify canonical
> scientific data. Please return the assigned VM/service user and UID/GID, `/pubapps`
> project path and quota, local container path/quota, SSH route, available Podman
> version, and login-linger state.
>
> The requested production URL is
> `https://islet-explorer.rc.ufl.edu/islet-explorer/` to
> `/islet-explorer/`. The project requests publication of its backend on VM loopback
> at `127.0.0.1:3838`; please confirm that the UFRC proxy can reach that target or
> return the required private backend address and agreed port. A different assigned
> target requires the project to update and re-review its deployment configuration,
> activation, and verification scripts before deployment.
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
> release as not Restricted, contingent on retained PI/data-owner and agreement-owner
> approvals. Expected use is 1–8 reviewers total, no more than five concurrent sessions,
> and two UF maintainers. Maintainer HPG usernames are `[PRIMARY]` and `[BACKUP]`; their
> public keys will be supplied through this restricted ticket.
>
> The project will provide the exact release candidate's CycloneDX SBOM, lint/static
> analysis results, dependency information, and Syft/Grype container-scan evidence. The
> project release gate permits no High or Critical image findings. Please confirm the
> evidence is sufficient for review, complete UFRC's required production
> security/availability review, and identify the ongoing UFIT network-scan and approval
> process.
>
> Please confirm the remaining items and disposition of the existing individual
> Risk Assessment, identify the
> UFRC-monitored application/proxy log destinations and retention requirements, and
> provide the procedure for creating, changing, and revoking the shared Basic Auth
> credential. Please return the production-approval reference and support escalation
> route when setup and review are complete.

Attach the completed, genuinely attested `UFIT_PubApps_Instance_Request.docx` and
the `RC83797` allocation confirmation to the existing ticket. The form contains
the architecture/data-flow map. Cross-reference the existing individual assessment
and provide additional requested owner/legal evidence through the restricted
project record. Do not send this template with unresolved placeholders.

## Database and storage position

The reviewer deployment does **not** benefit from a managed PostgreSQL service. The
workload is immutable, analytical, small-user-count, and local to one PubApps instance;
PostgreSQL would add a network service, credentials, patching, backup, and availability
obligations without providing a needed multi-writer capability.

The production release uses a hybrid, embedded layout:

- Keep OME-TIFF images as files; storing the original approximately 109 GiB of imagery as database
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
- Keep per-islet drill-down CSVs as direct-lookup files (5,214 in the original release). Measured point lookup
  was faster than consolidating them into DuckDB; package them as one archive for release
  transfer, then extract them into the immutable versioned data directory.
- Package each donor's segmentation boundaries as an indexed, read-only GeoPackage and
  query only the selected islet's bounding box. Do not use SQLite as the primary
  analytical store.

The 2026-09-08 canonical source preparation produced 5,382 `core_area` rows,
5,382 `comp` rows, 48,438 `targets_all` rows, 64,584 `markers_all` rows, and
21 `combined_pairs` definitions. The separate trajectory H5AD has 5,214
observations; the drill-down file count and trajectory count must not be used as
the count of all prepared runtime summaries. These measurements describe the
original 15-image cohort; any accepted future additions require new release counts
while preserving the original 15-donor records exactly.
See [`docs/ufrc_readiness.md`](docs/ufrc_readiness.md) for the source evidence and
the remaining release-parity check.

Configuration and service continuation may proceed while the release is being built, but the
manifest-backed runtime is a production prerequisite. Before promotion it must reproduce
the selected source's canonical prepared tables and plotted statistics, including row
counts and identifiers, with exact categorical values and documented numeric tolerance.
Any future expansion release must also pass zero-tolerance preservation checks
for original IDs, geometry, measurements, missingness and phenotype labels,
apply the verified original rules, gates and scVI model to the additions, and
validate a common 21-donor trajectory, UMAP and Leiden refit. The earlier joint
refits do not satisfy that contract. If parity or scientific acceptance is
incomplete, keep the reviewer site closed until a validated release exists.

## Information UFRC must return or confirm

Record these answers in the ticket or release record:

- VM hostname and SSH route from HPG.
- Service-account name and UID/GID.
- Persistent and local container paths and quotas.
- Available container runtime/version and login-linger state.
- Production URL, reachable backend address/port, and final proxy path.
- WebSocket and TIFF HTTP Range proxy behavior.
- Basic Auth ownership and credential creation/change/revocation procedure.
- UFRC-monitored log destination and retention/access requirements.
- Existing individual Risk Assessment number, remaining requests, and final
  disposition; any relationship UFRC records to its main PubApps assessment.
- Security-evidence submission, UFRC scan/review, and production-approval process.
- Support escalation contact and expected response route during journal review.

## Release and go-live evidence

Maintain one release record. The UFIT/UFRC compliance evidence is:

- Lint results with zero errors and as few documented warnings as practical; documented
  input validation for browser-controlled parameters.
- A maintained SBOM plus project code/dependency/container security results. The project
  image gate has no High or Critical findings; any other remaining CVE is documented.
- Appropriate industry-standard one-event-per-line application logs sent to a
  UFIT-monitored location.
- Accessibility and UF branding review/remediation results.
- UFRC production acceptance and PI/data-owner attestation.

The project's reproducibility and operations record also retains:

- The 40-character Git commit and release tag if one was created.
- Immutable local image ID (or registry digest) and prior rollback image reference.
- SHA-256 data/image manifest with file counts and sizes.
- Reconciled islet, phenotype, donor, and image counts used in the manuscript and app.
- Resolved R package/version manifest, pinned Python requirements, base/image digests,
  and SBOM.
- Health-check, restart, rollback, and representative five-session smoke-test notes.
- Project log-minimization checks and the recorded UFRC answers about log access and
  retention.
- An off-campus test showing Basic Auth protects the app, WebSocket, downloads, and every
  image/sidecar path.
- A TIFF test showing a small byte-range request returns `206 Partial Content` with the
  correct `Content-Range` without transferring the entire image.

Freeze the accepted Git commit, immutable image reference, data manifest, and feature
flags throughout a review round. Do not run `git pull`, rebuild from a moving branch,
or automatically update production during review.

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
- Apply a security fix through a newly tested immutable image reference and document any
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
