# Panc Floyd UFRC launch preparation

Status: **administrative preparation in progress; no production release or UFRC go-live acceptance**
Updated: **2026-10-16**

This note is the current migration-preparation reference while a uniform workflow for all retained donors is being planned. It prepares the UFRC coordination packet without treating the offline development app as a production release. The project and reviewer-facing application are both named **Panc Floyd**. The existing request lists `islet-explorer.rc.ufl.edu/islet-explorer/`; the concise reply below gives that submitted route as the DNS/proxy target. Application and data-processing changes are deferred until the reprocessing decisions are resolved and the workflow is accepted.

## Current migration status

- UFIT canceled individual assessment `IRM-456246` on 2026-10-15. The closure notice says the unit will use fully de-identified data, the PubApps packet was submitted to Research Computing, and the project will be added to the PubApps master assessment. Keep the closure notice in the restricted operating record. Cancellation closes that individual assessment; it does not establish production approval.
- Continue through existing UFRC request `83810` and configuration reference `SYSCONFIG-8173`. The owner confirms that the PubApps instance request form was submitted with the original ticket on **2026-08-31**. Treat that submission as the existing request record; do not create a duplicate form or ticket. Ask UFRC to confirm that the submitted form and evidence are linked to the master assessment, and whether it wants an amended form once the final release inventory is accepted. UFRC's [risk-assessment guidance](https://docs.rc.ufl.edu/services/web_hosting/risk_assessment/) states that the instance request form is attached to the main PubApps assessment and that UFIT-RC requests an individual assessment when it determines one is needed (checked 2026-10-15).
- The current target is 20 donors after the permanent removal of donor 78/HPAP-078, documented in the [owner directive](owner_directive_20261015.md). Uniform processing has not started; D1–D5 in the [reprocessing plan](uniform_reprocessing_plan_20261015.md) remain unresolved. The app is offline pending accepted processing and a new release candidate. Its last development root is not a production candidate. Preserve the frozen 15-donor release as a historical comparison reference; it is not the selected launch candidate.
- The collection remains treated as Sensitive for confidential peer review. The cancellation notice's de-identification statement is not itself a data-owner disclosure authorization or a change to the deployment's overall classification.
- Trial allocation `RC83797` is recorded through 2026-11-30. The owner plans to purchase long-term allocation; treat the purchase and final order reference as an owner action, not as a launch deadline.
- The current request lists `https://islet-explorer.rc.ufl.edu/islet-explorer/` for the Panc Floyd application. The backend address and port are not confirmed; `127.0.0.1:3838` is still only the project request.
- UFRC comment 13 on ticket `83810`, received 2026-10-16, says the VM will be provisioned next and asks for DNS and proxy requirements for the first application. The proposed reply below identifies the application as Panc Floyd and gives the requested hostname/path and proxy requirements; UFRC must supply the final proxy-reachable backend after assigning the VM.

## Concise DNS/proxy reply for UFRC comment 13 (not saved or sent)

The Outlook reply-draft request returned `403 Access denied`. The text below is ready to paste into the existing ticket reply and is limited to DNS and proxy requirements.

**Outlook thread:** `[Request 83810] PubApps Instance setup for the maigan group`

**Reply body:**

> Hi Alex,
>
> The first application is named Panc Floyd. Please configure `https://islet-explorer.rc.ufl.edu/islet-explorer/` on HTTPS port 443, with public DNS pointing to the UFRC-managed reverse proxy, not directly to the VM. Please preserve the `/islet-explorer/` path.
>
> The Shiny backend listens on port 3838; please provide the proxy-reachable VM address and port. The proxy must support Shiny WebSockets and pass HTTP Range requests through for OME-TIFF files. Keep the application, downloads, and image routes behind Public Restricted whole-site Basic Auth. Publication reviewers need external access to this authenticated URL during peer review.
>
> Thanks,
>
> Justin

This reply gives the requested DNS/proxy target and the minimum proxy behaviors needed for the Panc Floyd Shiny application and publication reviewers. The ticket's original instance form remains the request record; the DNS/proxy request does not authorize app activation or change the planned data classification.

## Existing request form and later amendment

The instance request form was submitted with the original ticket on 2026-08-31. No replacement form is needed for the current status update. If UFRC requests an amendment after uniform processing is accepted, update the submitted form or provide the format UFRC specifies, then append the final inventory and production details. Suggested status wording for that later amendment:

> Status as of October 15, 2026: UFRC support ticket `83810` (configuration reference `SYSCONFIG-8173`) remains the active PubApps request. UFIT canceled individual assessment `IRM-456246` on October 15 and stated that the project will be added to the PubApps master assessment. The cancellation notice states that the project will use fully de-identified data. The last development root contained 20 donors after the owner-directed removal of donor 78/HPAP-078; the app is offline while uniform processing and the final production inventory are in progress. That root is not the production release. Production deployment remains pending UFRC assessment and configuration, final owner-approved data and disclosure documentation, candidate-specific release evidence, and UFRC and PI/data-owner acceptance.

Keep the form's PI name and compliance acknowledgment fields for the PI to verify and date. Keep personal contacts, public keys, fiscal details, data-use approvals, and restricted assessment correspondence in the approved restricted channel. Do not turn the draft classification into an open-data or public-access request.

## Owner and UFRC details to close

Record personal contact information, public SSH keys, fiscal details, and approvals only in the approved restricted ticket or operating record.

| Detail | Current state / next action |
|---|---|
| Master-assessment linkage | UFIT closure notice says the project will be added; obtain confirmation in ticket `83810` |
| Individual assessment | `IRM-456246` canceled 2026-10-15; retain closure notice and do not reopen unless UFRC directs |
| Allocation | Owner plans long-term purchase; obtain current quote/order route and retain final order reference |
| PI and maintainers | Confirm active HiPerGator accounts, named primary/backup maintainers, and their public keys through the restricted channel |
| Form fields | Confirm department, fiscal contact/ChartFields, journal and review dates, data owner, applicable IRB/DUA/MTA or not applicable, retention, and actual PI attestation date |
| VM and access | UFRC comment 13 says VM provisioning is next; after assignment confirm hostname, proxy-reachable address/port, service account and UID/GID, SSH route, Podman version, linger, storage paths, quotas, and effective permissions |
| Public route | Panc Floyd reply text retains the requested `islet-explorer.rc.ufl.edu/islet-explorer/` route; DNS target and final backend remain unconfirmed |
| Reviewer availability | Confirm external reviewer access, credential issuance, setup/acceptance lead time, and the journal review window; plan activation before reviewers need the site |
| Operations | Confirm monitored log destinations, retention/access, UFIT scan process, security/availability review, escalation route, and allocation renewal/retention timing to prevent an interruption during review |

## Handoff to the processing work

After the current uniform-processing work is complete, reconcile its accepted donor list, methods, counts, limitations, and output size with the request form and support-ticket packet. Only then start production-builder/runtime changes, candidate generation, and candidate-specific release evidence. This note does not authorize a release build, deployment, credential distribution, or production launch.
