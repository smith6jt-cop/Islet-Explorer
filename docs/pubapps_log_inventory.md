# Islet Explorer PubApps Log Inventory and Handling Standard

Status: **pre-go-live security artifact; owner and UFRC confirmations remain open**
Repository snapshot reviewed: **2026-09-01**
Production profile: **journal-review deployment, Public Restricted access, Sensitive collection**

## Purpose and scope

This document inventories the log and error output that the current Islet Explorer
repository can produce when deployed on UFRC PubApps. It also defines the minimum
format, redaction, access, rotation, retention, and acceptance controls for reviewer
hosting.

The inventory covers:

- R/Shiny application output and uncaught errors;
- Shiny Server worker logs;
- rootless Podman and user-systemd service output;
- UFRC reverse-proxy, TLS, and Basic-Auth logs;
- deployment-command output and release evidence; and
- operational access to those records.

Browser downloads and immutable scientific release files are not logs. They are noted
where needed so that download content is not accidentally copied into an operational
log or support ticket.

## Current disposition

The deployment is **not log-ready for production until the unchecked controls in the
go-live checklist are closed**. In particular:

1. Current R output includes case IDs, islet keys, image IDs or filenames, selected
   image URLs, full paths, and raw `conditionMessage()` text.
2. Most application output is free-form and some credential diagnostics are multiline;
   it does not meet the one-event-per-line standard below.
3. `preserve_logs true` retains Shiny worker logs, but the repository defines no
   host-side rotation, byte cap, or deletion schedule.
4. The UFRC reverse proxy is outside this repository. Its access-log fields, query-string
   treatment, authentication-event handling, central collection, and retention are not
   yet known.
5. A shared Basic-Auth account cannot identify an individual reviewer. If UFRC or the
   data owner requires person-level access attribution, named accounts or another
   approved access method are required before go-live.

Existing controls that reduce risk are:

- `sanitize_errors true` keeps detailed R errors out of the reviewer browser while
  retaining server-side diagnostics (`deploy/shiny-server.conf:6-7`).
- Shiny worker logs are directed to `/var/log/shiny-server`, which is a dedicated
  read-write host mount (`deploy/shiny-server.conf:21` and
  `deploy/islet-explorer.container.in:23`).
- The activation script creates the runtime log directory with mode `0750`
  (`deploy/activate-release.sh:80-85`). Actual file modes and ACLs still require a
  deployment-time check.
- Production sets `DEBUG_CREDENTIALS=0` and `VIEWER_DEBUG=0`, disables AI, and rejects
  `KEY` or `OPENAI_API_KEY` if supplied (`deploy/islet-explorer.container.in:39-43` and
  `deploy/container-entrypoint.sh:15-24`).

## Log-flow inventory

| Stream | Current producer and transport | Current or expected destination | Security disposition |
|---|---|---|---|
| Shiny worker output | R `message()`, `cat()`, `warning()`, uncaught `stop()`, package output | Container `/var/log/shiny-server`; host `runtime/logs` mount | Sensitive operational record. Identifier-bearing emissions must be removed or redacted. Rotate explicitly because worker logs are preserved. |
| Shiny Server supervisor output | Shiny Server lifecycle, worker start/stop, configuration and fatal errors | Container stdout/stderr and/or Shiny Server's native log destination; exact base-image behavior must be verified | Internal operational record. Do not assume it shares the worker-log mount. Confirm with `podman logs` and the user journal. |
| Container and unit output | Entrypoint checks, health failures, Podman runtime, systemd restart events | User journal for `islet-explorer.service`; visible through `journalctl --user` and `podman logs` | Internal operational/security record. Current messages do not print secret values, but can reveal UID/GID and paths. UFRC must set journal limits and retention. |
| Reverse-proxy access | HTTPS, WebSocket upgrade, downloads, and TIFF Range requests | UFRC-managed proxy/access log and any central UFRC collector | Sensitive access metadata. Log a normalized route, never a query string, credential, cookie, or full TIFF/case filename. |
| Basic-Auth events | Successful and failed authentication at the UFRC proxy | UFRC-managed authentication/security log | Sensitive security record. Never log the `Authorization` header or password. A shared username is not individual attribution. |
| Proxy/TLS errors | TLS, upstream connection, timeout, WebSocket, and Range failures | UFRC-managed error/security log | Potentially Sensitive because native error text may include URI, backend, header, username, or client address. Redact before collection. |
| Deployment commands | `build-image.sh`, `scan-image.sh`, release verification, activation, rollback, and deployment verification | Operator terminal unless deliberately captured; some state is written under release evidence/deployments | Internal release record, not an application log. It may contain image digests, Git commits, hostnames, UIDs, absolute paths, and registry names. Never capture shell history containing credentials. |
| Host administrative access | SSH, login, service-account and privilege events | UFRC host/security logs | UFRC-owned security record. Restrict to named maintainers and UFRC administrators; confirm review and retention. |

The Quadlet does not select a Podman log driver. UFRC must confirm the effective driver,
whether `podman logs` and the user journal contain the same records, and which store is
authoritative.

## Current R application emissions

The table groups repository evidence by risk. Line numbers identify the reviewed
snapshot and may move as the app changes.

| Category | Current evidence and example content | Risk | Required production handling |
|---|---|---|---|
| Startup and capability checks | `app/shiny_app/R/00_globals.R:42-119,174-183` reports package state, app version, working directory, and Python path. `app/shiny_app/R/viewer_helpers.R:45-114` reports channel-file and image/resource paths and file counts. | Paths, usernames, runtime topology, release state, and unpublished aggregate counts may be exposed. | Keep stable capability/result event codes and release ID. Remove working directory, home path, interpreter path, and source-file paths. Log counts only where operationally necessary. |
| Release and data loading | `app/shiny_app/R/data_loading.R:235-397,425-502,555-608` reports joins, row counts, H5AD/Excel paths, manifest-backed release path, and raw read errors. `app/shiny_app/R/01_runtime_data_backend.R:22-223` and `app/shiny_app/R/utils_safe_join.R:6-21` raise detailed validation errors. | Raw errors can contain absolute paths, artifact names, schema details, values, or package text. Multiline `stop()` text breaks line discipline. | Emit an allowlisted event such as `data.release.validation_failed` with release ID and error code. Keep raw exception text only in a short-lived, access-restricted diagnostic capture after approval. |
| Plot interactions | `app/shiny_app/R/mod_plot_server.R:545-547,994,1158-1170` reports outlier state, donor count, and raw plot-click key. | A click key combines case/islet pseudonyms; it is Sensitive scientific and reviewer-behavior metadata. | Do not log clicks, selections, donor IDs, case IDs, islet keys, or phenotype values. Aggregate performance metrics only if approved. |
| Trajectory processing and interactions | `app/shiny_app/R/mod_trajectory_server.R:467,568-613,932-994`, `app/shiny_app/R/paga_layout.R:34-70`, and `app/shiny_app/R/barlow_pseudotime.R:140` report feature-processing errors, row counts, cell-range selections, case IDs, islet keys, and raw Python/R errors. | Case/islet identifiers and reviewer-selected ranges are Sensitive; Python exceptions can include data paths and values. | Suppress routine click and selection events. Use stable compute-success/failure codes, duration, and non-identifying row counts only. Never emit raw selected feature, case, islet, donor, or stack trace. |
| Drilldown and spatial data | `app/shiny_app/R/drilldown_helpers.R:152-158` reports `combined_id`; `app/shiny_app/R/spatial_helpers.R:87-181` reports `imageid`, row counts, and raw errors. | `combined_id` and `imageid` are donor-derived pseudonyms. Errors may include Parquet paths or query detail. | Replace identifiers with no value or, only if incident response requires correlation, a short-lived keyed digest controlled by UFRC. Prefer event code plus result/duration. |
| Segmentation | `app/shiny_app/R/segmentation_helpers.R:59-208` reports case IDs, GeoPackage/GeoJSON paths, feature counts, and raw geometry errors. `app/shiny_app/app.R:648` reports raw segmentation exceptions. | Case IDs, image filenames, paths, geometry errors, and counts are Sensitive in combination. | Log `segmentation.query_failed` or `segmentation.asset_missing` without case, filename, path, bounding box, or raw exception. |
| Image lookup and viewer | `app/shiny_app/R/viewer_helpers.R:218-232` reports matched filenames and missing case IDs. With `VIEWER_DEBUG=1`, `app/shiny_app/R/mod_viewer_server.R:122-127` logs the selected image and full iframe source. | The iframe source contains `image_url` and encoded `channel_config`; filenames can encode a case ID. Query strings may also reach proxy access logs. | Keep `VIEWER_DEBUG=0` in every production and staging unit. Do not log full URI/query strings, image names, case IDs, or channel configuration. Normalize proxy routes to route classes. |
| Credential and AI diagnostics | `app/shiny_app/R/ai_helpers.R:18-140` can report user/home, `.Renviron` search paths, file size/readability, key presence and length, and API base. `app/shiny_app/R/ai_helpers.R:255-564` and `app/shiny_app/R/mod_ai_assistant_server.R:203-205` can report endpoint/model and raw upstream error detail. | Although the key value is not intentionally printed, these are credential-discovery signals. Upstream errors can include submitted content or vendor detail. A future debug change could expose a secret. | AI must remain disabled, no AI credential may enter the container, and `DEBUG_CREDENTIALS` must remain `0`. Treat any credential-debug or AI request event in production as a security alert. Never enable this diagnostics path on the reviewer instance. |
| Generic caught and uncaught errors | `app/shiny_app/app.R:633-648` prints raw drilldown/segmentation exceptions; R/package warnings and Shiny worker failures can also enter worker logs. | `conditionMessage()` and stack traces are uncontrolled fields and can contain Sensitive values, paths, or injected newlines. | Map expected failures to stable codes. Sanitize unexpected messages before a single-line record; do not copy raw traces into tickets without a security review. Browser error sanitization must remain enabled. |

### Output that is not an operational log

- `app/shiny_app/R/00_globals.R` writes an app-load flag only when an operator explicitly
  sets `ISLET_APP_LOAD_FLAG` to a path beneath R's temporary directory. Production leaves
  it unset. It is an optional state marker, not a log, and must not become a durable audit
  source.
- Statistics and inflammation download handlers write report content to a Shiny-provided
  response file. Those files are reviewer downloads, not server logs. The tabs are
  disabled in the reviewer profile.
- Plot, Spatial, and other browser exports are session responses. Do not copy exported
  rows into access logs, error messages, monitoring labels, or support tickets.
- `release.env`, `MANIFEST.sha256`, `.verified`, SBOMs, scan results, and activation state
  are release evidence. Protect them as internal records, but do not mix them with
  rotating runtime logs.

## Production event standard

### One event per line

Application-owned operational events must be UTF-8 JSON Lines. Each physical line is
one complete JSON object and one event. Values must be JSON-escaped; embedded CR/LF,
tabs, control characters, and terminal escape sequences must never create another log
line. A native component that cannot emit JSON must still emit exactly one escaped
event per physical line and remain in its own stream.

Required application fields:

| Field | Rule |
|---|---|
| `ts` | UTC RFC 3339 timestamp with `Z` suffix |
| `level` | One of `DEBUG`, `INFO`, `WARN`, `ERROR`, `SECURITY`; production defaults to `INFO` |
| `event` | Stable allowlisted identifier such as `data.release.loaded` |
| `component` | `shiny`, `shiny-server`, `container`, `proxy`, or `auth` |
| `release` | Immutable `ISLET_RELEASE`; no mutable tag |
| `result` | Small allowlisted value such as `ok`, `denied`, `failed`, `timeout` |
| `error_code` | Stable application code when applicable; never raw exception text |
| Metrics | Numeric `duration_ms`, `rows`, `bytes`, or status only when needed and approved |
| Correlation | UFRC proxy request ID, if available; never a cookie, Shiny session token, Basic-Auth password, or stable reviewer fingerprint |

Example acceptable application event:

```json
{"ts":"2026-08-31T14:05:16Z","level":"INFO","event":"data.release.loaded","component":"shiny","release":"pubapps-2026-01","result":"ok","rows":5214}
```

Example acceptable proxy event after field minimization:

```json
{"ts":"2026-08-31T14:06:02Z","level":"INFO","event":"http.request","component":"proxy","route":"viewer.image","method":"GET","status":206,"bytes":1048576,"duration_ms":42,"request_id":"opaque-proxy-id"}
```

The examples define shape, not a requirement to retain the shown count or request ID.
UFRC and the data owner may choose fewer fields.

### Prohibited fields

The following must never appear in application, proxy, container, journal, monitoring,
or support-ticket logs:

- Basic-Auth password or raw `Authorization` header;
- API keys, tokens, cookies, Shiny session tokens, private keys, or `.Renviron` contents;
- full query strings, request bodies, WebSocket payloads, CSV/PNG contents, or AI prompts;
- full referrer URLs when they can contain query strings;
- case ID, islet key, `combined_id`, donor ID, `imageid`, image filename, or source-data row;
- selected phenotype/marker values, coordinates, bounding boxes, or reviewer export data;
- full filesystem paths, home directory, username, environment dump, SQL text, or stack trace;
- raw `conditionMessage()` or other untrusted multiline text; and
- the neutral reviewer password in test commands, shell history, screenshots, or ticket
  comments.

The UFRC proxy must consume Basic Auth and strip the `Authorization` header before
forwarding to the Shiny backend. It must log a normalized route such as `app`,
`websocket`, `download`, `viewer.image`, or `static`, rather than `$request_uri` or the
TIFF basename. In particular, the Avivator iframe query contains `image_url` and
`channel_config` and must not be recorded.

### Allowed Sensitive security metadata

Source IP address, failed username, neutral shared-account username, user agent, and
request time may be necessary for brute-force detection or incident response. If UFRC
collects them, they remain Sensitive, are restricted to UFRC security and named
maintainers, and follow the approved security-log retention schedule. They must not be
copied into the application worker log. With a shared account, these fields must not be
represented as proof of an individual reviewer's identity.

## Rotation, retention, access, and disposal

The following is the proposed minimum baseline. UFRC policy, an incident hold, an
agreement, or a written data-owner requirement may require longer retention. The final
values and the responsible system must be recorded before go-live.

| Record | Proposed baseline pending approval | Access and disposal |
|---|---|---|
| Shiny worker logs in `runtime/logs` | Rotate daily or at 100 MiB, whichever occurs first; compress completed logs; retain 30 days; cap the set at 2 GiB | Service account, two named maintainers, and UFRC administrators only. Delete automatically after the approved period; no copy in a public release or repository. |
| Container/user-journal output | Retain 30 days searchable, subject to UFRC journal limits; cap by bytes as well as time | Named maintainers and UFRC administrators. UFRC must identify the authoritative store and deletion mechanism. |
| Proxy access/auth/security logs | Retain through the active review period and for 90 days after reviewer access closes, or the longer UFRC-required period | UFRC security/operations and explicitly authorized maintainers. Central storage, backup, export, and deletion inherit the Sensitive classification. |
| Approved temporary diagnostic capture | Maximum 72 hours, scoped to a specific incident, with no credentials or scientific identifiers | Start only with maintainer/UFRC approval. Record the ticket and delete or sanitize when the incident closes. Never enable credential debugging. |
| Release/security evidence | Retain with the publication release record for the owner-approved publication retention period | This is not runtime telemetry. Restrict write access, preserve integrity, and store separately from reviewer-facing files. |

Additional controls:

- Keep runtime log directories at `0750` or narrower and files at `0640` or narrower.
  Verify effective UID, GID, umask, ACLs, and parent-directory traversal on the VM.
- Do not grant reviewers, journal personnel, or the shared Basic-Auth identity filesystem
  or log-collector access.
- Do not ship logs to an analytics, SaaS, or external error-reporting service unless the
  data owner and UFRC approve that destination in writing.
- Encrypt any approved transfer of logs, attach only a sanitized excerpt to a restricted
  UFRC ticket, and remove local analyst copies after the ticket closes.
- Synchronize VM and proxy clocks through UFRC-managed time service so events can be
  correlated without logging extra identifiers.
- Alert before either the log filesystem or journal reaches its cap. Range requests can
  be high volume and must not be allowed to exhaust the PubApps allocation.
- An incident, litigation, or records hold overrides automatic deletion only when issued
  by the authorized UF/UFRC office. Document the hold and subsequent release.
- Record quarterly access reviews while reviewer authentication remains active and
  document final deletion when the review deployment is retired.

## Owner confirmations required

The PI/data owner and technical owner must complete and retain these decisions:

- [ ] Confirm that de-identified case/islet/image codes, scientific selections, and
  authentication metadata are Sensitive for this combined reviewer deployment.
- [ ] Approve the allowed application event fields and the prohibition on scientific
  identifiers and raw errors.
- [ ] Name the two technical maintainers who may read runtime logs.
- [ ] Name the security/incident contact and the escalation route for suspected
  disclosure or credential exposure.
- [ ] Approve the final application, unit-journal, proxy/security, diagnostic, and
  release-evidence retention periods.
- [ ] Confirm whether shared Basic Auth supplies sufficient auditability. If individual
  attribution is required, replace it with named accounts before go-live.
- [ ] Confirm that no AI credential, external telemetry key, or external log destination
  is authorized for the reviewer instance.
- [ ] Confirm that support tickets and journal correspondence will not contain raw logs,
  reviewer credentials, case IDs, islet IDs, image filenames, or scientific rows.
- [ ] Approve the final log-deletion record when reviewer hosting ends.

## UFRC confirmations required

Ask UFRC to answer these items in the restricted PubApps ticket:

- [ ] Identify the reverse-proxy product, access/auth/error log locations, authoritative
  collector, and staff roles that can read each store.
- [ ] Provide the exact proxy log field list. Confirm that it excludes the
  `Authorization` header, password, cookies, session tokens, request/response bodies,
  query strings, full referrer, TIFF basename, and `image_url`/`channel_config`.
- [ ] Confirm that the proxy strips the Basic-Auth `Authorization` header before sending
  the request to `127.0.0.1:3838`.
- [ ] Describe successful/failed-login logging, source-IP handling, rate limiting,
  lockout/brute-force alerts, and how the shared credential is rotated and revoked.
- [ ] Confirm the effective Podman log driver and whether user-journal records duplicate
  `podman logs`; designate one authoritative service stream.
- [ ] Configure and document host-side rotation for `runtime/logs`. `preserve_logs true`
  without rotation is not acceptable.
- [ ] Provide journal and proxy byte/time caps, central retention, backup behavior,
  deletion mechanism, and disk-usage alerts.
- [ ] Confirm directory/file ownership, group membership, umask/ACLs, and that only the
  service account, named maintainers, and authorized UFRC administrators can read logs.
- [ ] Confirm clock synchronization and the proxy request-ID format, if request IDs are
  retained.
- [ ] Provide the security incident, credential-exposure, and emergency log-preservation
  process, including after-hours contact.
- [ ] Confirm whether UFRC policy requires a longer schedule than the proposed baseline;
  record the approved schedule in this document and the ticket.

## Pre-go-live leakage test

Run this test against staging using only synthetic canary strings, never a real password
or API key:

1. Record the image digest, data-release ID, test start/end time, and canary values in a
   restricted test record.
2. Start the app, exercise all four reviewer tabs, click plot and trajectory points,
   change an image, load segmentation/spatial views, export a CSV/PNG, and request a TIFF
   byte range.
3. Trigger a controlled missing-image or invalid-selection error and one failed Basic
   Auth attempt with a synthetic username.
4. Inspect `runtime/logs`, `journalctl --user -u islet-explorer.service`, `podman logs`,
   and a UFRC-provided proxy/auth log export for the test interval.
5. Verify that application-owned records are one JSON object per line and native
   component records are one escaped event per line.
6. Search for the synthetic password, `Authorization`, `Cookie`, session-token names,
   `KEY=`, `OPENAI_API_KEY`, `image_url`, `channel_config`, full query strings, source
   paths, and known test case/islet/image identifiers. The required result is zero.
7. Verify proxy events use route classes and that failed authentication is visible
   without recording a password or raw header.
8. Force a worker-log rotation and verify the running service continues, completed logs
   are compressed with restricted permissions, and expiration removes the correct file.
9. Verify an unauthorized local account and the reviewer identity cannot read any log
   store.
10. Attach only the pass/fail record and sanitized field samples to the release evidence;
    do not attach a raw log bundle.

Any prohibited-field match is a failed security gate. Keep production closed, correct
the producer or collector, repeat the full test, and retain the remediation record.

## Go-live checklist

### Application and container

- [x] Shiny browser errors are sanitized.
- [x] Reviewer production flags disable AI, Statistics, Inflammation, credential debug,
  and viewer debug.
- [x] The entrypoint rejects AI credentials in the reviewer profile without printing
  their values.
- [x] Shiny worker logs have a dedicated host-mounted directory created as `0750`.
- [ ] Identifier-bearing plot, trajectory, drilldown, spatial, segmentation, and image
  lookup messages are removed or transformed to the approved event schema.
- [ ] Raw `conditionMessage()`, full paths, multiline errors, and environment diagnostics
  are removed from routine production output.
- [ ] Application-owned events pass the one-event-per-line validation.
- [ ] Actual worker-log and temporary-file modes/ACLs are verified on the PubApps VM.
- [ ] A credential-debug or AI event in production produces an operational alert.

### UFRC proxy and host

- [ ] UFRC has supplied and approved the proxy/auth/error field inventory.
- [ ] Query strings, authentication headers, cookies, bodies, full referrers, and image
  filenames are absent from proxy logs.
- [ ] The proxy consumes Basic Auth and does not forward its header to Shiny.
- [ ] Failed-login controls, rotation/revocation, rate limiting, and alert ownership are
  documented and tested.
- [ ] Podman/user-journal destination, duplication, access, byte cap, and retention are
  documented.
- [ ] Host-side worker-log rotation, compression, byte cap, retention, and deletion are
  installed and tested.
- [ ] Disk and journal-capacity alerts have named recipients.
- [ ] Named maintainer, service-account, UFRC administrator, and unauthorized-user access
  tests pass.

### Approval and evidence

- [ ] Owner and UFRC retention decisions are written into the restricted ticket.
- [ ] Shared-account auditability is accepted in writing, or named authentication is in
  place.
- [ ] The staging leakage test passes with zero prohibited-field matches.
- [ ] A sanitized log-format sample, rotation result, access test, and deletion test are
  retained with the immutable reviewer release evidence.
- [ ] Incident contacts and the after-hours escalation path have been exercised.
- [ ] Technical owner and UFRC security/operations sign off before reviewer credentials
  are delivered.

## Final approval record

Complete this section in the restricted deployment record; do not add personal contact
details or credentials to the public repository.

| Decision | Approved value |
|---|---|
| Application log retention/cap | `[VALUE / OWNER / DATE]` |
| User-journal retention/cap | `[VALUE / OWNER / DATE]` |
| Proxy/auth log retention/cap | `[VALUE / OWNER / DATE]` |
| Release-evidence retention | `[VALUE / OWNER / DATE]` |
| Authoritative container log store | `[VALUE]` |
| UFRC proxy/collector | `[VALUE]` |
| Named maintainer group | `[RESTRICTED RECORD REFERENCE]` |
| Shared-auth auditability accepted | `[YES / NO / APPROVER / DATE]` |
| Leakage-test evidence | `[RELEASE-EVIDENCE PATH OR TICKET]` |
| UFRC approval | `[TICKET / APPROVER / DATE]` |
| Technical-owner approval | `[APPROVER / DATE]` |
