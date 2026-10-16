# Panc Floyd UFRC launch preparation

Status: **administrative preparation in progress; no production release or UFRC go-live acceptance**
Updated: **2026-10-16**

This note is the current migration-preparation reference while a uniform workflow for all retained donors is being planned. It prepares the UFRC coordination packet without treating the offline development app as a production release. The project and reviewer-facing application are both named **Panc Floyd**. The existing request lists `islet-explorer.rc.ufl.edu/islet-explorer/`; the proposed reply keeps that submitted route and asks UFRC whether the hostname should be updated to match the application name. Application and data-processing changes are deferred until the reprocessing decisions are resolved and the workflow is accepted.

## Current migration status

- UFIT canceled individual assessment `IRM-456246` on 2026-10-15. The closure notice says the unit will use fully de-identified data, the PubApps packet was submitted to Research Computing, and the project will be added to the PubApps master assessment. Keep the closure notice in the restricted operating record. Cancellation closes that individual assessment; it does not establish production approval.
- Continue through existing UFRC request `83810` and configuration reference `SYSCONFIG-8173`. The owner confirms that the PubApps instance request form was submitted with the original ticket on **2026-08-31**. Treat that submission as the existing request record; do not create a duplicate form or ticket. Ask UFRC to confirm that the submitted form and evidence are linked to the master assessment, and whether it wants an amended form once the final release inventory is accepted. UFRC's [risk-assessment guidance](https://docs.rc.ufl.edu/services/web_hosting/risk_assessment/) states that the instance request form is attached to the main PubApps assessment and that UFIT-RC requests an individual assessment when it determines one is needed (checked 2026-10-15).
- The current target is 20 donors after the permanent removal of donor 78/HPAP-078, documented in the [owner directive](owner_directive_20261015.md). Uniform processing has not started; D1–D5 in the [reprocessing plan](uniform_reprocessing_plan_20261015.md) remain unresolved. The app is offline pending accepted processing and a new release candidate. Its last development root is not a production candidate. Preserve the frozen 15-donor release as a historical comparison reference; it is not the selected launch candidate.
- The collection remains treated as Sensitive for confidential peer review. The cancellation notice's de-identification statement is not itself a data-owner disclosure authorization or a change to the deployment's overall classification.
- Trial allocation `RC83797` is recorded through 2026-11-30. The owner plans to purchase long-term allocation; treat the purchase and final order reference as an owner action, not as a launch deadline.
- The current request lists `https://islet-explorer.rc.ufl.edu/islet-explorer/` for the Panc Floyd application. The backend address and port are not confirmed; `127.0.0.1:3838` is still only the project request.
- UFRC comment 13 on ticket `83810`, received 2026-10-16, says the VM will be provisioned next and asks for DNS and proxy requirements for the first application. The proposed reply below identifies the application as Panc Floyd, asks whether the requested hostname/path should change to match, and provides the routing contract; UFRC must supply the final proxy-reachable backend after assigning the VM.

## Reply text for UFRC comment 13 (not saved or sent)

The Outlook reply-draft request returned `403 Access denied`. The text below is ready to paste into the existing ticket reply.

**Outlook thread:** `[Request 83810] PubApps Instance setup for the maigan group`

**Reply body:**

> Hi Alex,
>
> Thanks. The first application we plan to deploy for the Panc Floyd project is also named Panc Floyd. The August 31 request lists `islet-explorer.rc.ufl.edu/islet-explorer/`; please confirm whether that hostname and path can remain the public route for the Panc Floyd application or whether UFRC recommends updating them to match the application name. Our requirement is that publication reviewers can reach the authenticated application throughout the manuscript-review period. We do not yet have the journal's exact review dates, so please advise on setup and production-acceptance lead time and any steps we should schedule before review begins. The requested DNS and proxy behavior is:
>
> - **Public DNS and URL:** the currently requested route is `https://islet-explorer.rc.ufl.edu/islet-explorer/` on HTTPS port 443. Please confirm whether UFRC recommends retaining this route or changing the hostname/path to match the Panc Floyd application name, and whether a change requires an amendment to the submitted form. Point the public DNS record to the UFRC-managed high-availability reverse proxy and its TLS certificate. No public DNS record should point directly to the backend VM.
>
> - **Path and routing:** preserve the `/islet-explorer/` prefix, including its trailing slash, for the Shiny app. Its static image routes are `/islet-explorer/local_images/`, `/islet-explorer/images/`, and `/islet-explorer/img_data/`; they remain on the same host and behind the same authentication. If your proxy requires a path rewrite, please provide the exact rule so we can align the Shiny configuration.
>
> - **Backend:** the container listens on `0.0.0.0:3838`. Our current rootless Podman configuration requests publication on the VM at `127.0.0.1:3838`. Since the reverse proxy is separate from the VM, please confirm whether UFRC provides a forwarding arrangement to that loopback port or assign the private VM address and port the proxy should use. Please return the assigned VM hostname/address, reachable backend address/port, and any network restrictions; we will align and validate our deployment configuration before activation.
>
> - **Proxy behavior and access:** the requested access mode is Public Restricted, with whole-site Basic Auth for a Sensitive collection. Please require authentication for the app, WebSocket, downloads, all three image aliases, OME-TIFFs, and XML/JSON sidecars. Publication reviewers must be able to reach the authenticated site from outside UFRC for the full manuscript-review window; please confirm that access policy, credential issuance, or network restrictions will not require reviewers to have UFRC accounts or be on a campus network. Preserve the WebSocket upgrade and long-lived Shiny connections; pass HTTP Range requests through so TIFF requests return `206` with the correct `Content-Range` without buffering a complete multi-gigabyte file. Please preserve the forwarded host and HTTPS scheme, and confirm proxy timeouts, log destination/retention, and any other required headers or settings.
>
> There is one application and one canonical public hostname; no separate database or additional public service endpoint is needed. The instance request form was submitted with ticket `83810` on August 31, 2026. UFIT canceled individual assessment `IRM-456246` on October 15 and said the project would be added to the PubApps master assessment; please let us know if any assessment update is needed in this ticket.
>
> The local app is offline pending uniform reprocessing and owner acceptance of the 20-donor root; donor 78/HDL078 is permanently excluded. This request is to record DNS/proxy requirements and coordinate setup during provisioning, not to activate reviewer access. Trial allocation `RC83797` is recorded through November 30, 2026, and the owner plans a long-term purchase. Please advise the renewal/purchase lead time and any steps needed to avoid an interruption during the manuscript-review period. We will schedule activation after processing is accepted and provide the final release inventory and launch evidence before reviewers need access.
>
> Thanks,
>
> Justin

This reply answers the current DNS/proxy question without guessing the VM address, makes the Panc Floyd name explicit, and asks UFRC to resolve the hostname choice and reviewer-access timeline in the same exchange. The ticket's original instance form remains the request record; the DNS/proxy request does not authorize app activation or change the planned data classification.

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
| Public route | Panc Floyd reply text retains the requested `islet-explorer.rc.ufl.edu/islet-explorer/` route pending UFRC advice on whether the hostname/path should change to match the application name; DNS target and final backend remain unconfirmed |
| Reviewer availability | Confirm external reviewer access, credential issuance, setup/acceptance lead time, and the journal review window; plan activation before reviewers need the site |
| Operations | Confirm monitored log destinations, retention/access, UFIT scan process, security/availability review, escalation route, and allocation renewal/retention timing to prevent an interruption during review |

## Handoff to the processing work

After the current uniform-processing work is complete, reconcile its accepted donor list, methods, counts, limitations, and output size with the request form and support-ticket packet. Only then start production-builder/runtime changes, candidate generation, and candidate-specific release evidence. This note does not authorize a release build, deployment, credential distribution, or production launch.
