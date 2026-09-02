# Islet Explorer PubApps Log Inventory and Project Handling Guidance

Status: **application hardening recorded; exact-release and UFRC confirmations remain open**
Repository snapshot reviewed: **2026-09-02**
Production profile: **journal-review deployment, Public Restricted access, Sensitive collection**

## Purpose and scope

This document inventories the log and error output that Islet Explorer can produce on
UFRC PubApps. It distinguishes the required hosting outcome—appropriate
industry-standard one-event-per-line logging in a UFIT-monitored location—from the
project's more detailed JSONL, minimization, rotation, retention, and verification
recommendations. UFRC may satisfy collection, access, and retention centrally; the
project must record the values UFRC assigns rather than treating every recommendation
below as a separate UFRC approval gate.

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

The documented UFIT logging outcome remains open until the exact release candidate
produces appropriate one-event-per-line application records and UFRC confirms the
UFIT-monitored location. Separately, the project must resolve its acceptance questions
about minimization, access, and retention. Repository hardening has removed the reviewed identifier-bearing
Plot, Trajectory, viewer, drilldown, spatial, segmentation, and root-render error
messages and has replaced their failure paths with fixed event codes. Production now
fails closed if credential or viewer debugging is enabled.

The remaining application output is primarily one-line capability and data-preparation
text plus native R/package output. Credential and AI failures now use fixed event codes,
provider details are not logged or shown, and production rejects the AI path before any
request; the production fail-closed controls must remain in place. The exact release
output still needs a source/runtime review for line discipline and minimization. Two
external questions also remain: UFRC owns the proxy/auth collector and its retention
rules, and a shared Basic-Auth identity cannot attribute activity to an individual
reviewer. If UFRC or the data owner requires person-level attribution, use named
accounts or another approved access method.

Existing controls that reduce risk are:

- `sanitize_errors true` keeps detailed R errors out of the reviewer browser;
  `preserve_logs false` avoids indefinitely preserving per-worker files and
  `log_file_mode 0640` restricts any created Shiny log (`deploy/shiny-server.conf`).
- Shiny worker logs are directed to `/var/log/shiny-server`, which is a dedicated
  read-write host mount (`deploy/shiny-server.conf:24` and
  `deploy/islet-explorer.container.in:23`).
- The activation script creates the runtime log directory with mode `0750`
  (`deploy/activate-release.sh:88-89`). Actual file modes and ACLs still require a
  deployment-time check.
- Production sets `DEBUG_CREDENTIALS=0` and `VIEWER_DEBUG=0`, disables AI, and rejects
  either debug mode or `KEY`/`OPENAI_API_KEY` if supplied
  (`deploy/islet-explorer.container.in`, `deploy/container-entrypoint.sh`, and
  `app/shiny_app/R/00_globals.R`).
- Browser-controlled validation failures use a fixed one-line JSON event that cannot
  accept the rejected value (`app/shiny_app/R/00_audit.R`).

## Log-flow inventory

| Stream | Current producer and transport | Current or expected destination | Security disposition |
|---|---|---|---|
| Shiny worker output | R `message()`, `cat()`, `warning()`, uncaught `stop()`, package output | Container `/var/log/shiny-server`; host `runtime/logs` mount | Potentially Sensitive operational record. Reviewed identifier-bearing emissions were removed; validate the exact release. Shiny does not preserve completed worker logs, and any separately retained host copy follows the applicable UFRC or project limits. |
| Shiny Server supervisor output | Shiny Server lifecycle, worker start/stop, configuration and fatal errors | Container stdout/stderr and/or Shiny Server's native log destination; exact base-image behavior must be verified | Internal operational record. Do not assume it shares the worker-log mount. Confirm with `podman logs` and the user journal. |
| Container and unit output | Entrypoint checks, health failures, Podman runtime, systemd restart events | User journal for `islet-explorer.service`; visible through `journalctl --user` and `podman logs` | Internal operational/security record. Current messages do not print secret values, but can reveal UID/GID and paths. Record the effective UFRC host limits and retention. |
| Reverse-proxy access | HTTPS, WebSocket upgrade, downloads, and TIFF Range requests | UFRC-managed proxy/access log and any central UFRC collector | Sensitive access metadata. Log a normalized route, never a query string, credential, cookie, or full TIFF/case filename. |
| Basic-Auth events | Successful and failed authentication at the UFRC proxy | UFRC-managed authentication/security log | Sensitive security record. Never log the `Authorization` header or password. A shared username is not individual attribution. |
| Proxy/TLS errors | TLS, upstream connection, timeout, WebSocket, and Range failures | UFRC-managed error/security log | Potentially Sensitive because native error text may include URI, backend, header, username, or client address. Redact before collection. |
| Deployment commands | `build-image.sh`, `scan-image.sh`, release verification, activation, rollback, and deployment verification | Operator terminal unless deliberately captured; some state is written under release evidence/deployments | Internal release record, not an application log. It may contain image digests, Git commits, hostnames, UIDs, absolute paths, and registry names. Never capture shell history containing credentials. |
| Host administrative access | SSH, login, service-account and privilege events | UFRC host/security logs | UFRC-owned security record. Restrict to named maintainers and UFRC administrators; confirm review and retention. |

The Quadlet does not select a Podman log driver. On the provisioned host, record the
effective driver, whether `podman logs` and the user journal contain the same records,
and which store is authoritative.

## Reviewed R application emissions

The table records the disposition of the reviewed snapshot. Re-run the same source and
runtime review against the exact candidate because dependencies and future code changes
can add output.

| Category | Reviewed disposition | Remaining release check |
|---|---|---|
| Startup and capability | Working-directory, interpreter-path, viewer resource-path, image-filename, and setup-count diagnostics were removed. Some non-identifying package/version and capability messages remain. | Confirm remaining native and application records stay on one physical line and disclose no unnecessary environment detail. |
| Release and data loading | Production failure paths use fixed event codes; data-preparation code retains free-form aggregate/status messages and development fallbacks. | Exercise production startup with the immutable manifest and inspect failures without copying raw data or paths into tickets. |
| Plot, Trajectory, PAGA, and Barlow | Reviewed click keys, selected features/ranges, case/islet values, and raw compute exceptions were removed from operational output. Fixed failure codes and validation events remain. | Confirm UI values cannot enter log fields and validation JSON remains one line. |
| Drilldown, spatial, and segmentation | Reviewed combined IDs, image IDs, paths, geometry details, and raw render/read exceptions were replaced with fixed event codes. | Confirm expected errors expose no case, islet, filename, path, coordinate, or raw exception. |
| Viewer | Reviewed filename/path lookup and iframe-source output was removed. `VIEWER_DEBUG` now fails closed in both the entrypoint and R production startup. | Keep debug disabled and verify UFRC proxy logs do not retain `image_url`, `channel_config`, or the full query string. |
| Credential and AI diagnostics | Credential discovery and provider failures now emit fixed codes without paths, endpoint/model values, raw provider errors, or secret metadata. Production disables AI, rejects credentials, skips credential-file discovery, and fails closed if credential debug is enabled. | Treat any AI/credential-debug activity on the reviewer instance as a configuration failure; never place a real secret in a canary test. |
| Framework/native errors | Browser sanitization remains enabled; reviewed root render catches now emit fixed event codes. R, Shiny Server, and package failures can still produce native text. | Inspect the exact-candidate stream for credentials, unnecessary Sensitive values, uncontrolled newlines, and paths before reviewer access. |

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

## Recommended structured application format

### One event per line

UFIT requires an appropriate industry-standard one-event-per-line format. UTF-8 JSON
Lines is the project's preferred structured format for application-owned events: each
physical line is one complete JSON object and one event. Values must be JSON-escaped;
embedded CR/LF, tabs, control characters, and terminal escape sequences must never
create another log line. The exact UFIT outcome does not require every native framework
or capability message to adopt this project schema; those records should still remain
one event per physical line and in a separately identified stream.

Recommended application fields:

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

Ask UFRC to confirm that its proxy consumes Basic Auth without forwarding the
`Authorization` header to Shiny and that its logs do not retain the password, raw
header, or sensitive query values. Normalized route classes such as `app`, `websocket`,
`download`, `viewer.image`, or `static` are a project recommendation in place of a full
request URI or TIFF basename. In particular, the Avivator iframe query contains
`image_url` and `channel_config` and should not be recorded.

### Allowed Sensitive security metadata

Source IP address, failed username, neutral shared-account username, user agent, and
request time may be necessary for brute-force detection or incident response. If UFRC
collects them, they remain Sensitive, are restricted to UFRC security and named
maintainers, and follow the approved security-log retention schedule. They must not be
copied into the application worker log. With a shared account, these fields must not be
represented as proof of an individual reviewer's identity.

## Rotation, retention, access, and disposal

The following is the project's recommended baseline, not a statement of UFRC policy.
UFRC may provide host- or collector-level controls, and UFRC policy, an incident hold,
an agreement, or a written data-owner requirement may require different retention.
Before reviewer access, record the applicable destination, access rule, retention, and
responsible system; use the project values below only where the project owns the store
and no controlling requirement supersedes them.

| Record | Proposed baseline pending approval | Access and disposal |
|---|---|---|
| Shiny worker logs in `runtime/logs` | Rotate daily or at 100 MiB, whichever occurs first; compress completed logs; retain 30 days; cap the set at 2 GiB | Service account, two named maintainers, and UFRC administrators only. Delete automatically after the approved period; no copy in a public release or repository. |
| Container/user-journal output | Retain 30 days searchable, subject to UFRC journal limits; cap by bytes as well as time | Named maintainers and UFRC administrators. Record the authoritative store and deletion mechanism UFRC provides. |
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
- If UFRC confirms that project-plan item 6 applies to this Sensitive deployment,
  record quarterly access reviews while reviewer authentication remains active.
  Document final disposal when the review deployment is retired.

## Project owner acceptance questions

The PI/data owner and technical owner should resolve and retain these project decisions:

- [ ] Confirm that de-identified case/islet/image codes, scientific selections, and
  authentication metadata are Sensitive for this combined reviewer deployment.
- [ ] Approve minimization of application logs, including exclusion of credentials,
  unnecessary scientific identifiers, and uncontrolled raw errors.
- [ ] Name the two technical maintainers who may read runtime logs.
- [ ] Confirm whether shared Basic Auth supplies sufficient auditability. If individual
  attribution is required, replace it with named accounts before go-live.
- [ ] Confirm that no AI credential, external telemetry key, or external log destination
  is authorized for the reviewer instance.

## Documented UFRC logging outcome

- [ ] Application logs use an appropriate industry-standard one-event-per-line format
  and are sent to a UFIT-monitored location.

## Project acceptance questions for UFRC

Ask UFRC to answer these operational questions in the restricted PubApps ticket. They
support project acceptance but are not presented as separately documented UFRC gates:

- [ ] Identify the UFIT-monitored destinations for application and proxy/auth records,
  who may read them, and the applicable retention/deletion rules.
- [ ] Confirm that proxy/auth logging does not retain credentials or unnecessary
  Sensitive query values and that the Basic-Auth header is consumed before the request
  reaches Shiny.
- [ ] Provide the procedure for creating, changing, and revoking the reviewer credential
  and the support route for a suspected disclosure.

Recommended operating decisions, not separate UFRC gates, include selecting an
authoritative Podman/journal stream, documenting effective rotation and byte caps,
confirming clock/request-ID behavior, recording alert ownership, and exercising log
preservation and deletion. Adopt details from the tables above when the project owns
the relevant store; otherwise record UFRC's implementation.

## Recommended leakage exercise

Run this project assurance exercise against the exact release candidate through the
confirmed UFRC route before reviewer credentials are distributed. Use synthetic canary
strings, never a real password or API key:

1. Record the image digest, data-release ID, test start/end time, and canary values in a
   restricted test record.
2. Start the app, exercise all four reviewer tabs, click plot and trajectory points,
   change an image, load segmentation/spatial views, export a CSV/PNG, and request a TIFF
   byte range.
3. Trigger a controlled missing-image or invalid-selection error and one failed Basic
   Auth attempt with a synthetic username.
4. Inspect `runtime/logs`, `journalctl --user -u islet-explorer.service`, `podman logs`,
   and a UFRC-provided proxy/auth log export for the test interval.
5. Verify that application and native-component events remain one physical line. For
   records using the project JSONL format, also validate the JSON object and schema.
6. Search for the synthetic password, `Authorization`, `Cookie`, session-token names,
   `KEY=`, `OPENAI_API_KEY`, `image_url`, `channel_config`, full query strings, source
   paths, and known test case/islet/image identifiers. The required result is zero.
7. Verify proxy events use route classes and that failed authentication is visible
   without recording a password or raw header.
8. If the project retains a durable worker-log store, force a rotation and verify the
   running service continues, completed logs are compressed with restricted permissions,
   and expiration removes the correct file.
9. Verify an unauthorized local account and the reviewer identity cannot read any log
   store.
10. Attach only the pass/fail record and sanitized field samples to the release evidence;
    do not attach a raw log bundle.

This detailed exercise is recommended project evidence rather than a separately
mandated UFRC test. Any credential or prohibited-field match still fails the project's
release assurance: keep reviewer access closed, correct the producer or collector,
repeat the relevant checks, and retain the remediation record.

## Release-acceptance checklist

### Documented UFRC logging outcome

- [ ] Application logs use an appropriate industry-standard one-event-per-line format
  and are sent to a UFIT-monitored location.

### Project acceptance questions and recommended evidence

- [ ] The exact release candidate emits no credentials, uncontrolled multiline errors,
  or unnecessary Sensitive scientific fields.
- [ ] Browser errors remain sanitized and the reviewer profile fails closed if AI,
  credential debug, viewer debug, or an AI credential is enabled.
- [ ] UFRC has supplied the monitored application and proxy/auth destinations plus the
  applicable access and retention rules.
- [ ] UFRC has confirmed safe Basic-Auth handling and a credential
  creation/change/revocation procedure; the owner has accepted shared-account
  auditability or selected an approved named-account alternative.
- [ ] The technical owner has accepted the project logging controls and recorded UFRC's
  answers before reviewer credentials are distributed.

- Validate project JSONL events against the project schema and retain a sanitized sample.
- Run the synthetic-canary leakage exercise across the application, relevant host
  stream, and a UFRC-provided proxy/auth sample.
- Where the project owns a log store, test its rotation, permissions, byte cap,
  retention, and deletion; otherwise record the corresponding UFRC controls.
- Record the authoritative runtime stream, alert/incident contact, and any approved log
  preservation procedure.

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
| Recommended leakage-exercise evidence | `[RELEASE-EVIDENCE PATH OR TICKET / NOT RUN]` |
| UFRC approval | `[TICKET / APPROVER / DATE]` |
| Technical-owner approval | `[APPROVER / DATE]` |
