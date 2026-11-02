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

UFRC provisioned the VM `pubmaigan1` on 2026-10-22 and the DNS aliases on
2026-10-23: `islet-explorer.rc.ufl.edu` and `devislet-explorer.rc.ufl.edu`, the
production/development pair UFRC's deployment guide sets up for every
application; since 2026-10-23 a staging instance (`islet-explorer-dev`,
`pubmaigan1:3839`, 4 GiB cap, 6 GiB for the v20 candidate) runs beside production (9 GiB cap) for the
development name, so a candidate can be exercised through the real proxy and
Basic Auth before it replaces what reviewers use (since 2026-11-02 06:07 UTC it
runs the v20 candidate `c974946`, release `pubapps-2026-11-uniform-v20`, with a 6 GiB cap;
before that the Plan A Level 2 candidate `0b55429`, release `pubapps-2026-10-uniform-v10-perf2`,
from 2026-10-26 19:23 UTC, and the Level 1 candidate `e10b2ee` of 2026-10-25; all recorded
below). Production was promoted to perf2 on 2026-10-26 20:51 UTC by the owner ("Push what
is complete so far to the live app"); **both instances were stopped at 21:31 UTC that day,
and production has stayed stopped** (the staging instance was started again only by the v20
activation). `pubapps-2026-10-uniform-v10` (`d55680a`) stays on the VM for
`rollback-release.sh`. The
schema-5 release is active on the VM with its backend published on the VM's
private address `10.13.41.21:3838`, the host:port target UFRC's proxy needs
(`vpubweb.rc.ufl.edu` cannot reach the VM loopback). UFRC configured both
vhosts with whole-site Basic Auth on 2026-10-23: the proxy mounts the project's
`/pubapps` tree, reads the project-maintained `htpasswd` files `prod/auth` and
`dev/auth`, serves a landing page from `prod/html` and `dev/html`, and serves
the three image aliases itself from symlinks there into the active release
(`deploy/README.md` § 4). Both public URLs answer; from the campus network the
reviewer-access verifier passed all 15 checks on each hostname with UFRC's
temporary test account. The project manages the credential files itself
(bcrypt entries hashed locally; UFRC's guide calls Public Restricted Basic Auth
"credential lists provided by you"), and the guide also expects the project to
manage the service account's `authorized_keys` after initial setup.

The selected reviewer profile remains Public Restricted access to a Sensitive
collection: Plot, Trajectory, Viewer, and Spatial; AI, Statistics, and Inflammation
disabled; immutable local scientific files; no managed database, upload endpoint,
or production external API. The URL is
`https://islet-explorer.rc.ufl.edu/islet-explorer/`; the UFRC project name is
Panc Floyd.

### Current local app and release candidate (2026-10-22)

The only intended local development URL is
`http://10.15.152.7:8080/islet-explorer/`, proxied to loopback port 3840. It
was offline at the owner's direction from 2026-10-15, served a frozen copy of
the app (tag `frozen-app-20261020`) from 2026-10-20 for the owner's review of
the reprocessed data, and since 2026-10-22 serves the main tree on the uniform
root `data/development/uniform-20donor-20261022-v10/` (20 donors, 4,763
endocrine objects: 1,487 Islets and 3,276 Clusters, with 5,207 Singles counted
but not analyzed). Donor 78 (HDL078 / HPAP-078) stays permanently removed from
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
| Staging instance beside production | **activated 2026-10-23 17:35 UTC** (`activate-release.sh --instance islet-explorer-dev --port 3839 --memory 4g --publish-address 10.13.41.21`, tooling `ee6e3e3`) on the same release: own Quadlet unit and container, `runtime/islet-explorer-dev/{logs,tmp}`, `deployments/islet-explorer-dev.active.state`; `verify-deployment.sh` passed (healthy, `10.13.41.21:3839` exclusive, Range `206 bytes 0-15/7792217538`). Production was then re-activated at 17:42 UTC with the 9 GiB cap (verifier passed again); both services active and healthy, resident memory about 0.8 GiB each at idle on the 15 GiB VM. Serves UFRC's `devislet-explorer.rc.ufl.edu` once its vhost exists (backend `http://pubmaigan1:3839/islet-explorer/`, own Basic Auth entry for maintainers). Logs: `tmp/ufrc-release-20261022/deploy-vm-d55680a-staging-3839-*.log`, `deploy-vm-d55680a-production-9g-*.log` |
| Browser checks against the deployed container | passed 2026-10-22 18:37-18:39 UTC from the workstation through a two-hop SSH tunnel to the VM's loopback backend: `check-app.mjs` (every check, no accessibility rule violations, no browser or HTTP errors), `check-ui-review-20261021.mjs` (15 of 15: logo, plot types in both cards, composition by donor over 20 donors, near-empty measures hidden, axis fit, frozen controls on four tabs, Spatial colourings including the nPOD fields, hover and click linking, shared selection, clear) and `check-trajectory-links.mjs` (all steps, no errors). Reports: `tmp/ufrc-release-20261022/evidence/browser-vm-d55680a/` |
| UFRC VM | provisioned 2026-10-22: `pubmaigan1.ufhpc`, service account `maigan` (home `/pubapps/maigan`, 1 TiB Qumulo share), RHEL 9.8, rootless Podman 5.8.2 (store on the 40 GB `/local`), cgroups v2, linger on, 2 vCPU / 15 GiB; reached from HPG login nodes as `ssh pubmaigan1`; project path `/pubapps/maigan/islet-explorer/` created, release tooling copied, 500 MB test transfer at about 100 MB/s |
| UFRC proxy and Basic Auth acceptance | configured 2026-10-23 (ticket `83810`); campus acceptance passed; the external run and the reviewer credentials are pending. UFRC's two vhosts on `vpubweb.rc.ufl.edu` proxy `/islet-explorer/` (WebSocket included) to `pubmaigan1:3838` and `pubmaigan1:3839`, with Basic Auth on `<Location />` against the project-maintained files `prod/auth` and `dev/auth` (mode 640; Apache reads them as the service account), a `DocumentRoot` per hostname (`prod/html`, `dev/html`: the landing page `deploy/docroot/index.html`, favicon, logos) and the three image aliases excluded from proxying and served by Apache from symlinks into the active release's `data/images` (installed the same day; both instances were re-activated the same evening with `activate-release.sh --docroot`, staging at 21:08 UTC and production at 21:15 UTC, each re-hashing the release and passing the on-VM verifier, so their deployment state records the DocumentRoot and the aliases follow activation and rollback; tooling `5bfde65` on the VM). Measured from the campus network with UFRC's temporary test account: `deploy/verify-reviewer-access.py` 15/15 on both hostnames (anonymous 401 on the app, the Avivator entry point, the three aliases, both sidecars, the download route and the WebSocket; authenticated 200, the exact 16-byte 206 through Apache, 101 on Shiny Server's SockJS route); mid-file 8 MiB Range reads about 60 MB/s; `deploy/browser/check-app.mjs --auth-file` through the dev hostname 49/49 (Shiny session over the proxied WebSocket, four tabs, CSV downloads, the Viewer canvas on an Apache-served TIFF with the 16-byte 206, 20 images inventoried, no browser errors, no external requests, zero axe violations in five scopes; Chrome 156, 2026-10-23 19:29 UTC); the proxy's access and error logs are readable at `/pubapps/logs/maigan/`. UFRC applied both requested corrections on 2026-10-23 (about 16:35 EDT): the backend lines read `http://` and `ServerName https://...` makes Apache's own redirects name `https://` (verified: `/islet-explorer` without the slash is `301 https://.../islet-explorer/` on both names). `.htaccess` in the DocumentRoots is not read (probed twice for 90 s with valid credentials). UFRC also lowered the vhosts' error-log level; the access logs (`common` format: client IP, account name, time, request line, status, bytes) are still written, so the retention question stands. After the changes: verifier 15/15 on production as `reviewer` and on staging as `maintainer`; browser suite 49/49 through production with the reviewer account (Chrome 156, 20:50 UTC; `check-app-prod-reviewer-account.json`). History: the backend was published on `10.13.41.21:3838` on 2026-10-22 21:16 UTC after measuring that UFRC's network filters the port (an HPG login node reaches the VM on port 22 only, the campus workstation not at all, the VM's firewalls are inactive); DNS for both names was done on 2026-10-23. Evidence: `tmp/ufrc-release-20261022/evidence/proxy-acceptance-20261023/`. Credentials: the production `reviewer` account was set on 2026-10-23 about 19:50 UTC with the credential helper (bcrypt, password from a 0600 file) and UFRC's `test` entry removed from `prod/auth`; the verifier passed 15/15 on production with that account. The staging `maintainer` account was set the same evening and `test` removed from `dev/auth`, so both auth files hold only project accounts. SSH access list requested: `smith6jt` (VM keys received from UFRC by email on 2026-10-22) |

### Plan A Level 1 candidate e10b2ee on the staging instance (2026-10-25)

Development is in-house (owner, 2026-10-25): candidates go to the staging instance
first, and production stays on `d55680a` until the owner decides. This candidate
carries the Plan A Level 1 responsiveness work (`docs/performance_plans_20261024.md`,
`docs/performance_baseline_20261024.md`).

| Item | Result |
|---|---|
| Source commit | `e10b2ee6a9bbb0f526d9f1cd3f39f088a24e6b5d` (branch `spatial-review-20261020`; `eba12b4` Level 1: WebSocket compression in the Shiny Server build, single first render, trajectory warm-up, payload trims, `activate-release.sh --local-data`; `4da5ac4` docs; `e10b2ee` body-parser pin), built from the clean linked worktree `/data/panc_CODEX/Islet-Explorer-release-20261022` |
| Container image, SBOM and High/Critical gate | the build at `4da5ac4` failed the gate on 2026-10-25 14:41 UTC on one advisory published after the previous scan's database: body-parser 2.3.0 (GHSA-27r3-jxmh-p828, CVE-2026-111181, High; fixed in 2.4.0, released 2026-10-22). Pinned (`overrides` 2.4.0 with the lock bump qs 6.16.1, `deploy/BASE_HARDENING.md`) and rebuilt: **passed 2026-10-25 15:06 UTC** on the workstation (rootless Podman 4.9.3, Syft 1.51.1, Grype 0.118.0, database of 2026-10-25): image `sha256:d8c470dea33fa1d040197c1ebf1dc594265a618e594fc18cfade03e2de772627` on base `rocker/shiny@sha256:dd10b638…`; Python and R probes passed; Grype 0 High / 0 Critical (4,305 Medium, 189 Low, 196 Negligible); OCI archive `islet-explorer-pubapps-2026-10-uniform-v10-perf1-e10b2ee.oci.tar` (2.2 GB). Evidence `tmp/ufrc-release-20261022/evidence/image-scan-{4da5ac4,e10b2ee}/`, logs `build-image-{4da5ac4,e10b2ee}.log` |
| Data bundle | built 2026-10-25 14:47–15:07 UTC from the worktree at `e10b2ee` by the unchanged builder from the same v10 root (`tmp/ufrc-release-20261022/build-data.sh`; 20 min with the sources in the page cache); full hashes passed 15:19 UTC, parity passed 15:44 UTC, inventory 15:44 UTC (`verify-data-e10b2ee.log`, `evidence/release-parity-e10b2ee.{json,md}`, `evidence/release-inventory-e10b2ee.{json,md}`) |
| Release manifest and `verify-release.sh` | passed 2026-10-25 15:47 UTC: release `pubapps-2026-10-uniform-v10-perf1` (its own name, so the running release's directory is untouched; `deploy/README.md` § 2) bound to image `d8c470de…` and commit `e10b2ee` (`DATA_SCHEMA=5`, `EXPECTED_IMAGES=20`, 4,819 data files / 917,879,043 bytes, 80 image files / 155,515,691,675 bytes); `MANIFEST.sha256` 4,909 entries, manifest hash `ea6527da8d2f498b176d5cc99a4f49643b78f77945ba36ccde1a3ef7ec9ff809` |
| Staging on the VM | the unchanged 145 GB image tree was copied on the VM from the v10 release (`cp -a`) and matched by checksum against the bound bundle with nothing re-sent; data files, attestations, `release.env`, `MANIFEST.sha256` and the OCI archive uploaded (archive checksum verified on the VM); deploy tooling on the VM at `e10b2ee` (`tooling/DEPLOY_TOOLING_COMMIT`) |
| Activation (staging instance only) | **activated 2026-10-25 16:03 UTC** as `islet-explorer-dev` on `10.13.41.21:3839` (4 GiB, `--docroot dev/html`, `--local-data /local/maigan/releases`): release re-hashed on the VM, `data/{geometry,cells,spatial}` (880 MB) copied to the local disk and verified against `MANIFEST.sha256`, unit mounts the copies over the container paths, service active and healthy, `verify-deployment.sh` passed (backend bound to `10.13.41.21:3839` only, Range `206 bytes 0-15/7792217538`), `trajectory_warm_ok` logged 29 s after the R process started. Production unchanged (`islet-explorer.service` on `pubapps-2026-10-uniform-v10`, `prod/html` aliases unchanged). Log `tmp/ufrc-release-20261022/deploy-vm-e10b2ee-staging-3839-20261025T154754Z.log` |
| Responsiveness through the proxy | `deploy/browser/bench-app.mjs` against `devislet-explorer.rc.ufl.edu` from campus and with the emulated home link, 2026-10-25 16:05–16:40 UTC: `permessage-deflate` negotiated on every connection; first screen 3.5 → 1.1 MB over the socket; on the emulated home link first screen 8.8 → 4.1 s, Plot change 2.9 → 1.5 s, Spatial 3.8 → 2.6 s; cold-process Trajectory 12.7 s → 0.26 s and PAGA 5.2 → 1.2 s. Table in `docs/performance_baseline_20261024.md`; evidence `tmp/ufrc-release-20261022/evidence/performance-baseline-20261024/level1-staging/` |
| Browser acceptance suite through the proxy | `deploy/browser/check-app.mjs --auth-file` against `devislet-explorer.rc.ufl.edu` with the staging account, 2026-10-25 16:44 UTC: passed, every check (36), no failed checks, no browser or console errors, no application HTTP errors, no accessibility rule violations. Report `tmp/ufrc-release-20261022/evidence/browser-staging-e10b2ee/check-app-staging-e10b2ee.json` |
### Plan A Level 2 candidate 0b55429 on the staging instance (2026-10-26)

Second candidate for the staging instance, and since 20:51 UTC the production release. It
carries the Plan A Level 2 responsiveness work and two owner requests of the day
(`docs/performance_plans_20261024.md`, `docs/performance_baseline_20261024.md`).

| Item | Result |
|---|---|
| Source commit | `0b554293395f15ffdc7e8fe6da59339152da1fbc` (branch `spatial-review-20261020`): `5095863` Plan A Level 2 (B1 direct `plot_ly()` figures for the two Plot cards and the Trajectory scatter with the plotly build inputs stripped from every linked plot, B3 `bindCache()` across sessions on the Plot cards, the composition card and the tissue scatter PNG, B4 `ragg` in the image; the tab bar in the header row and Clear selection in the Data controls at the owner's request), `0b55429` the outlier exclusion no longer shown and every object shown by default (owner). |
| Container image, SBOM and High/Critical gate | built twice the same day with the same dependency set: at `5095863` (gate passed 2026-10-26 17:37 UTC, image `6438e48e…`, 28 min) and, after the owner's outlier request, at `0b55429`: **gate passed 2026-10-26 18:03 UTC**, image `sha256:1065d532fe7bc57ca72eb7cbecd2557b5dc645195e22c4b09d0db8c0b41d6d5a`, 0 High / 0 Critical (4,305 Medium, 189 Low, 196 Negligible, the Level 1 counts), Syft 1.51.1 / Grype 0.118.0; `ragg` 1.5.0 and `systemfonts` 1.3.1 are in the image (`evidence/image-scan-0b55429/`). |
| Data bundle | built 2026-10-26 17:43–18:26 UTC from the worktree at `0b55429` by the unchanged builder from the same v10 root (`tmp/ufrc-release-20261022/build-data.sh`, 44 min, in parallel with the image rebuild because the same dependency set had just passed the gate; 4,819 data files, 918 MB beside the 155.5 GB of images); full hashes passed 18:38 UTC, scientific parity against the v10 root passed 19:04 UTC (`evidence/release-parity-0b55429.{json,md}`), privacy-safe inventory written 19:04 UTC (`evidence/release-inventory-0b55429.{json,md}`). |
| Release manifest and `verify-release.sh` | passed 2026-10-26 19:07 UTC: release `pubapps-2026-10-uniform-v10-perf2` (its own name beside the running perf1 release) bound to image `1065d532…` and commit `0b55429` (`DATA_SCHEMA=5`, 20 OME-TIFFs), manifest `c9c02a8267ac2f650e218b60bfde5e8a91704b4f6f374710a86a23482cb4a33d`; the scan evidence is in `attestations/`. |
| Staging on the VM | the unchanged 145 GB image tree was copied on the VM from the perf1 release (`cp -a`, about 15 minutes) into `releases/pubapps-2026-10-uniform-v10-perf2/data/images` before the build finished; deploy tooling synced to `0b55429` (`tooling/DEPLOY_TOOLING_COMMIT`). Upload 2026-10-26 19:07–19:16 UTC: the image tree matched by checksum with nothing re-sent, then the data files, attestations, `release.env`, `MANIFEST.sha256` and the 2.2 GB OCI archive. |
| Activation (staging instance only) | **activated 2026-10-26 19:23 UTC** as `islet-explorer-dev` on `10.13.41.21:3839` (4 GiB, `--docroot dev/html`, `--local-data /local/maigan/releases`): the image loaded from the archive, the release re-hashed on the VM (manifest `c9c02a82…`), the local copies of `data/{geometry,cells,spatial}` verified against `MANIFEST.sha256`, the service healthy within 40 s, `verify-deployment.sh` passed (backend bound to the private address only, OME-TIFF Range response verified), `dev/html` image aliases moved to the perf2 release. Production (`islet-explorer`, `pubapps-2026-10-uniform-v10`, image `ae95fd23…`) untouched. Log `tmp/ufrc-release-20261022/deploy-vm-0b55429-staging-3839-20261026T190742Z.log`. |
| Responsiveness through the proxy | `deploy/browser/bench-app.mjs` against `devislet-explorer.rc.ufl.edu` from campus and with the emulated home link, 2026-10-26 19:24–19:26 UTC, five cold runs each: `permessage-deflate` on every connection; against the Level 1 candidate on the same instance, first screen 3.6 → 2.7 s (campus) and 4.1 → 3.1 s (home), Plot change 1.5 → 0.72 s and 1.5 → 1.1 s, Spatial tab 3.2 → 1.3 s and 2.6 → 1.6 s, Spatial donor change 1.9 → 0.40 s and 1.2 → 0.41 s, warm revisit 2.5 → 1.8 s and 3.5 → 3.1 s; the process built its Plot cards once (six `plot_card_built` events for ten cold runs). One regression: the process's first Trajectory open took 4.2 s (the first PNG draw builds `ragg`'s font cache; later sessions 0.6–0.7 s), to be folded into the warm-up. Table and reading in `docs/performance_baseline_20261024.md`, "Plan A Level 2 on staging". |
| Browser acceptance suite through the proxy | `deploy/browser/check-app.mjs --auth-file` against `devislet-explorer.rc.ufl.edu` with the staging account, 2026-10-26 19:26 UTC: passed, every check (36), no failed checks, no browser or console errors, no accessibility rule violations; report `tmp/ufrc-release-20261022/evidence/browser-staging-0b55429/check-app-staging-0b55429.json`. |
| **Production activation** | **activated 2026-10-26 20:51 UTC** by the owner ("Push what is complete so far to the live app"), run in the session with the `!` prefix because the auto-mode classifier refuses production deployments from the assistant: the release worktree on `0b55429`, then `RELEASE_NAME=pubapps-2026-10-uniform-v10-perf2 PUBLISH_ADDRESS=10.13.41.21 deploy-vm.sh --skip-stage` (`islet-explorer` on `10.13.41.21:3838`, 9 GiB, `--docroot prod/html`, `--local-data /local/maigan/releases`): `verify-release.sh` re-hashed the release on the VM, the local data copies verified against `MANIFEST.sha256`, `verify-deployment.sh` passed (health healthy, backend bound to the private address only, OME-TIFF Range response verified), `prod/html` image aliases moved to the perf2 release; app loaded 20:50:52, release data loaded 20:51:00, trajectory warm-up ok 20:51:19 UTC. The previous release `pubapps-2026-10-uniform-v10` (`d55680a`, image `ae95fd23…`) stays on the VM for `rollback-release.sh`; the staging instance is unchanged on the same release. |
| Production checks through the proxy | `check-app.mjs --auth-file` against `islet-explorer.rc.ufl.edu` with the production account, 2026-10-26 20:54 UTC: passed, every check (36), no errors, no accessibility rule violations (`evidence/browser-production-0b55429/`). `bench-app.mjs`, three cold runs from campus against the 2026-10-24 production baseline: first screen 3.6 → 2.6 s, Plot change 1.2 → 0.75 s, Spatial tab 2.2 → 1.3 s, Spatial donor change 1.9 → 0.44 s, warm revisit 3.5 → 1.8 s, socket bytes to the first screen 3.5 → 1.2 MB, `permessage-deflate` on every connection (`evidence/performance-baseline-20261024/level2-production/`). |

### Candidate cc5cf8f on the staging instance (2026-10-26: American spelling, the warm-up's first PNG draw, three owner requests)

Third candidate of the day for the staging instance, on top of the Level 2 release now in
production: the owner's spelling rule ("American English throughout") applied to every
label, help text, comment and record, the process warm-up's first PNG draw (the 4 s
font-cache cost seen on the first Trajectory open of the Level 2 process), and three owner
requests (the Spatial card's groups "Insulin-deficient Islets" / "Insulin-containing
Islets", one IAA autoantibody box covering IAA and micro-IAA, readable hover labels on the
dark theme). Commit `080714b` (spelling and warm-up only) was built and gated the same way
(image `52057049…`, gate passed 20:07 UTC, bundle verified) but superseded by `cc5cf8f`
before it was staged. This candidate goes to staging first.

| Item | Result |
|---|---|
| Source commit | `cc5cf8f` (branch `spatial-review-20261020`): `080714b` American spelling in UI text, comments and documentation (code identifiers and API names unchanged) and the warm-up's 64 px PNG draw before the health check; `cc5cf8f` the Spatial card's group names, the merged IAA box (`apply_plot_filters()` expands it to the IAA and mIAA columns; the root's nPOD `n_aab` already counts the assay once) and dark-theme hover labels as a dark box with white text (measured: the text was `#6b6b6b` on the donor blue, now white). 437 tests pass; `check-app.mjs` (36) and the UI review suite (18) pass against the in-house service on each commit. |
| Container image, SBOM and High/Critical gate | **gate passed 2026-10-26 21:28 UTC** at `cc5cf8f`: image `sha256:695dae399c0f7c50134afd22c71717ba538020299142191fa0c904d19ad6e1ef`, 0 High / 0 Critical (4,305 Medium, 189 Low, 196 Negligible, the same dependency set as the day's earlier builds), Syft 1.51.1 / Grype 0.118.0, 21 min (`evidence/image-scan-cc5cf8f/`). |
| Data bundle | built 2026-10-26 21:07–21:39 UTC from the worktree at `cc5cf8f` by the unchanged builder (32 min, in parallel with the image build); full hashes passed 21:50 UTC, scientific parity against the v10 root passed 22:18 UTC (`evidence/release-parity-cc5cf8f.{json,md}`), privacy-safe inventory written 22:18 UTC (`evidence/release-inventory-cc5cf8f.{json,md}`). |
| Release manifest and `verify-release.sh` | passed 2026-10-26 22:21 UTC: release `pubapps-2026-10-uniform-v10-perf3` bound to image `695dae39…` and commit `cc5cf8f` (`DATA_SCHEMA=5`, 20 OME-TIFFs, 4,819 data files, 918 MB), manifest `574ac1bc381b9b76df9b0f37376540e78500b1953440ce98a35c45756a4941b9`. Deploy tooling on the VM synced to `18981c5` (the unit template's `SuccessExitStatus=143`) before the activation. |
| Staging on the VM | the unchanged 145 GB image tree copied on the VM from the perf2 release (`cp -a`) into `releases/pubapps-2026-10-uniform-v10-perf3/data/images` before the build finished; deploy tooling synced to `cc5cf8f`. The upload of 2026-10-26 22:21 UTC stopped in its first step (`rsync` exit 255; `deploy-vm-cc5cf8f-staging-3839-20261026T222148Z.log`), so nothing was activated. |
| Activation (staging instance only) | never activated; superseded on 2026-11-01 by the v20 candidate below, which carries these changes. |
| Responsiveness through the proxy | not measured (never staged). |
| Browser acceptance suite through the proxy | not run (never staged). |

### Candidate c974946 with the v20 root on the staging instance (2026-11-01)

The owner asked to "use the v20 data to rebuild the data for the app, then deploy to the dev
site". The data root is `data/development/uniform-20donor-20261030-v20`, produced on HiPerGator
and returned on 2026-11-01 (`docs/hpg_processing_20261031.md`, `docs/typing_handoff_20261031.md`);
the app code is the main branch with the phenotype-label work done since `cc5cf8f`. Release name
`pubapps-2026-11-uniform-v20`, for the staging instance only. VM state as found: both instances have
been stopped since 2026-10-26 21:31 UTC (their units, rendered before `18981c5`, show the container's
exit status 143 as `failed`); production stays stopped and untouched.

| Item | Result |
|---|---|
| Source commit | `c974946` (branch `spatial-review-20261020`). Since `cc5cf8f`: the app's phenotype labels of roots v15-v17 (`82ce23e`-`270d000`: macrophage marker profiles, broad-first display names, grouped choices, top-marker and CD45 labels); the release builder and parity check for those labels (`b75fc8b` the phenotype group aggregate helpers that `prep_data()` calls, `c8d942b` GDAL's single-feature GeoJSON limit lifted for the geometry sources, `70df797` and `38793ed` the parity check's definitions and their order); the image pins `2d59d3b`; the warm-up fix `c974946`. 502 app tests pass; the runtime lint passes. |
| Pre-build checks on the v20 root | a trial bundle built from the v20 root (hard links, `tmp/release-trial-v20/data`) passed `verify_release_data.R` (core) and every scientific parity check (summaries exact, maximum absolute difference 0, 7e-15 on the signed-distance table; 3,199,606 cells in the donor tissue tables and the spatial load; 30,029 geometry features in 60 regions); the one failed check was the expected source-identity mismatch of a trial bundle. The `c8d942b` image (app code as `c974946` except the warm-up fix) served that bundle on the workstation with `deploy/run.sh`'s production flags: `check-app.mjs --expanded true` passed 49 of 49 checks with no accessibility violations, and the R process peaked at 2.25 GiB resident. Its log showed `trajectory_warm_failed`, which led to `c974946` (next row). |
| Warm-up fix | the 64 px font-cache draw added in `080714b` called `plot.new()` with R's default margins, which a 64 px figure cannot hold ("figure margins too large"); the error ended the warm-up before the scanpy pseudotime and PAGA fits, so every production start from `080714b` on would have left both to the first Trajectory visitor (run alone in the container, 10.8 s and 5.6 s). No deployed release carried it (`cc5cf8f` was never staged). `islet_warm_png()` now draws on zero margins and closes its device on any exit; `test-warm-up.R`. |
| Container image, SBOM and High/Critical gate | the first build (`c8d942b`, Grype database built 2026-11-01) **failed the gate** on advisories published after 2026-10-26, all in the Shiny Server runtime: Node 24.20.0 (CVE-2026-111497 Critical; CVE-2026-105482, -105483, -106069 High) and follow-redirects 1.16.1 (GHSA-mrwp-f3q8-3g9p High). Pinned in `2d59d3b`: Node v24.22.2 (the newest 24 LTS patch release, containing the 24.22.1 security release) and follow-redirects 1.16.2 (`deploy/BASE_HARDENING.md`, 2026-11-01). At `c974946`: upstream Shiny Server tests 315 + 8 passed on Node v24.22.2; **gate passed 2026-11-02 04:34 UTC**, image `sha256:517e969d2333744f609533b8f4a7fd217ce8a852fe30ef3551e0d805c2a0db00`, 0 High / 0 Critical (4,724 Medium, 189 Low, 196 Negligible against the 2026-11-01 database), Syft 1.51.1 / Grype 0.118.0, 21 min; the SBOM lists Node 24.22.2 and follow-redirects 1.16.2 (`evidence/image-scan-c974946/`; the failed scan is in `evidence/image-scan-c8d942b/`). OCI archive `islet-explorer-pubapps-2026-11-uniform-v20-c974946.oci.tar` (2.2 GB). |
| Data bundle | built 2026-11-02 04:20–04:46 UTC from the worktree at `c974946` from the v20 root (`SOURCE_ROOT=... build-data.sh`, 26 min, in parallel with the image build once the trial parity had passed): schema 5, 19 artifacts, 4,904 data files (1.7 GB) beside the 145 GB of images. `verify-data.sh` against the v20 root: full hashes passed 04:58 UTC; scientific parity passed all 23 checks 05:28 UTC, the release commit equal to the clean checkout (`evidence/release-parity-c974946.{json,md}`); privacy-safe inventory written 05:28 UTC (20 TIFFs, 3,199,606 cells; `evidence/release-inventory-c974946.{json,md}`). |
| Final image on the final bundle (workstation) | image `517e969d…` serving `artifacts/release-c974946/data` on loopback with `deploy/run.sh`'s production flags: `trajectory_warm_ok` 18 s after the release data loaded; `check-app.mjs --expanded true` passed 49 of 49 checks, no accessibility violations, no application HTTP errors (the one console error is the bare origin's missing `/favicon.ico`, which the proxy's DocumentRoot serves); the R process peaked at 2.27 GiB resident. |
| Release manifest and `verify-release.sh` | passed 2026-11-02 05:50 UTC: release `pubapps-2026-11-uniform-v20` bound to image `517e969d…` and commit `c974946` (`DATA_SCHEMA=5`, 20 OME-TIFFs; 4,904 data files, 1.74 GB; 80 image files with sidecars, 155.5 GB), manifest `d9a29beda767eeb3beab57f6a5eb67836229dd8d9807c27208bd646366398d7c`; the scan evidence is in `attestations/`. The VM's deploy tooling stays at `18981c5`: since then only two browser-check scripts changed under `deploy/`, and the VM runs neither. |
| Staging on the VM | the 145 GB image tree was copied on the VM (`cp`) into `releases/pubapps-2026-11-uniform-v20/data/images` on 2026-11-01 (the v20 root's 20 image links resolve to the same files as v10's). Upload 2026-11-02 05:50–06:00 UTC: the image tree matched by checksum with nothing re-sent, then the 1.74 GB of data files, attestations, `release.env`, `MANIFEST.sha256` and the 2.2 GB OCI archive. |
| Activation (staging instance only) | **activated 2026-11-02 06:07 UTC** with `RELEASE_NAME=pubapps-2026-11-uniform-v20 INSTANCE=islet-explorer-dev PORT=3839 MEMORY=6g PUBLISH_ADDRESS=10.13.41.21 deploy-vm.sh` (`--docroot dev/html`, `--local-data /local/maigan/releases`): the image loaded from the archive and verified, the release re-hashed on the VM (manifest `d9a29bed…`), the local copies of `data/{geometry,cells,spatial}` verified against `MANIFEST.sha256`, `verify-deployment.sh` passed (health healthy, backend bound to `10.13.41.21:3839` only, OME-TIFF Range `206 bytes 0-15/7792217538`), `dev/html` image aliases moved to the v20 release; release data loaded 06:07:33, `trajectory_warm_ok` 06:07:54 UTC. Cap 6 GiB instead of staging's usual 4 GiB: v20's data files are 1.74 GB against the v10 releases' 918 MB, and production's 9 GiB cap is idle while that instance is stopped; the unit's memory peaked at 3.0 GiB (page cache included) during the browser suite. Production (`islet-explorer`) stays stopped, its state file and `prod/html` aliases unchanged on the perf2 release. Log `tmp/ufrc-release-20261022/deploy-vm-c974946-staging-3839-20261102T055046Z.log`. |
| Reviewer-access verifier through the proxy | `deploy/verify-reviewer-access.py` against `devislet-explorer.rc.ufl.edu` from campus with the staging `maintainer` account, 2026-11-02 06:09 UTC: passed 15 of 15 (anonymous 401 on the app, the Avivator entry point, the three image aliases, both sidecars, the download route and the WebSocket; authenticated 200, the exact 16-byte 206 on `ND_0112.ome.tiff`, 200 on both sidecars, 101 on the WebSocket). Evidence `tmp/ufrc-release-20261022/evidence/staging-c974946/verify-reviewer-access-dev-c974946.json`. |
| Browser acceptance suite through the proxy | `deploy/browser/check-app.mjs --auth-file --expanded true` against `devislet-explorer.rc.ufl.edu` with the staging account, 2026-11-02 06:09 UTC: passed 49 of 49 checks (Shiny session over the proxied WebSocket, four tabs, CSV downloads, the 20-donor census and the five additions' tissue views, the Viewer canvas on an Apache-served TIFF with the 16-byte 206, 20 images inventoried), no browser or console errors, no application HTTP errors, no accessibility rule violations in five scopes. Report `tmp/ufrc-release-20261022/evidence/staging-c974946/check-app-staging-c974946.json`. |

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
| Browser accessibility and behavior | 34 of 34 checks passed, including four tabs, keyboard navigation, four CSV downloads, actual tissue canvas, 16-byte HTTP 206, and blocked missing/external image URLs; zero axe violations across five scopes | 2026-10-23: 49 of 49 through UFRC's proxy with Basic Auth on the dev hostname (`--auth-file`), same scopes, zero violations. | Native production-profile rerun passed on 2026-09-09 with the verified `606fb53` data and corrected R worker launcher; Chrome 152.0.7977.82, axe 4.13.0. Manual-review items retained; no final-image or external proxy acceptance claimed |
| External HTTPS verifier | 10 synthetic tests passed without skips, including whole-site BasicAuth, bounded Range/body verification, sidecar MIME, WebSocket handshake and redirect/credential boundaries | [Acceptance procedure](../deploy/EXTERNAL_ACCEPTANCE.md) requires the assigned endpoint and real protected credential file |
| Hosted CI | Both `parse-and-test` and `deployment-quality` passed in [run 34375550621](https://github.com/smith6jt-cop/Islet-Explorer/actions/runs/34375550621) at `5809be4`, including all 22 parity assertions and 19 authentication/SBOM tests | Later runtime changes require their own CI and image acceptance |
| Later hosted CI at `000e970` | [Source run 34385986809](https://github.com/smith6jt-cop/Islet-Explorer/actions/runs/34385986809) and [image run 34385986858](https://github.com/smith6jt-cop/Islet-Explorer/actions/runs/34385986858) failed before application checks because the runner's Chrome package index had a hash mismatch | Local workflows now use Ubuntu 24.04 repositories for their Ubuntu packages. YAML/shell parsing and offline APT repository isolation pass; hosted rerun pending |
| Supported Shiny Server runtime | Checksum-locked official source build and Node 24.22.2 (since 2026-11-01; 24.20.0 before); 315 upstream and 8 real-R tests passed without skips; stripped-runtime app/static-route smoke passed; installed-package Syft/Grype scan reports zero matches | The original base-only scan found 33 High and 3 Critical entries. The complete `5809be4` application image subsequently passed the High/Critical gate. See [base remediation](../deploy/BASE_HARDENING.md) |
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

1. Credentials (project side; the vhosts went live on 2026-10-23): the
   production `reviewer` account is set and UFRC's `test` entry is gone from
   `prod/auth` (2026-10-23, credential helper: bcrypt, password read from a
   0600 file, never on a command line); the staging `maintainer` account is
   set and `test` is gone from `dev/auth` too. Still to do: delivery of the
   reviewer credential through the journal's confidential channel. UFRC applied the two vhost corrections on 2026-10-23
   (`http://` backend lines, `ServerName https://...`). Still open with UFRC:
   the retention of the access logs under `/pubapps/logs/maigan/` (they
   record reviewer IP addresses) and whether the proxy forwards the
   `Authorization` header to the backend (`RequestHeader unset Authorization`
   would stop the credential crossing the private hop). UFRC's own account for
   its availability checks is a separate `htpasswd` entry, removed afterwards.
2. External acceptance from a clean non-UF network with
   `deploy/verify-reviewer-access.py` (`deploy/EXTERNAL_ACCEPTANCE.md`; passed
   from campus on 2026-10-23 with the test account on both hostnames) and the
   browser check through the proxy (`deploy/browser/check-app.mjs
   --auth-file`), then the four tabs by hand. Record the evidence in the
   candidate table above.
3. Institutional items still open: confirm that the August 31 instance form is
   linked to the main PubApps assessment and whether UFRC wants an amendment
   for the 20-donor inventory (individual assessment `IRM-456246` was canceled
   on 2026-10-15; do not reopen it unless UFRC directs); UFRC's own
   vulnerability scan and its record of the PI attestation (the repository
   copy of the form still shows the attestation date pending); the
   UFIT-monitored destination and retention for `runtime/logs` and the user
   journal ([log inventory](pubapps_log_inventory.md)); PI/data-owner go-live
   approval and the confidential credential-delivery channel.
4. Operations before reviewer access: a VM reboot test of the lingering user
   service (UFRC's Shiny guide warns that system-level restart is not
   configured in PubApps; the enabled Quadlet under linger is unproven across
   a reboot, which only UFRC can perform); allocation continuation before the
   trial term ends on 2026-11-30; a second named maintainer's key when one is
   designated.

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
