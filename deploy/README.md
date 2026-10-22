# Islet Explorer on UFRC PubApps

Production runbook for the journal-review deployment requested at:

`https://islet-explorer.rc.ufl.edu/islet-explorer/`

UFRC must confirm this hostname and path before it is given to reviewers; until then it
is a request value, not evidence that the service has been provisioned.

The UFRC access request is **Public Restricted** (whole-site Basic Auth), while
the collection is classified **Sensitive** during peer review. The app is
non-mutating: scientific inputs are read-only, browser CSV downloads are built
in session-scoped temporary storage, and only temporary files and operational
logs are writable.

This runbook defines the project's chosen production implementation and release
controls. A requirement here may be a project release invariant rather than a
separate UFRC approval gate. External ownership, provisioning, and approval are
tracked in [`PUBAPPS_UFRC_HANDOFF.md`](../PUBAPPS_UFRC_HANDOFF.md); proportional
logging requirements and additional hardening guidance are in
[`docs/pubapps_log_inventory.md`](../docs/pubapps_log_inventory.md).

The local development app at `http://10.15.152.7:8080/islet-explorer/`
(loopback port 3840) serves the main working tree on the reprocessed uniform root
`data/development/uniform-20donor-20261022-v10/` since 2026-10-22 (owner: "make the
test app the live app"), selected by `ISLET_DATA_ROOT` in the untracked
`tmp/local-app.env`. It was offline from 2026-10-15 until the uniform reprocessing
(see [the directive record](../docs/owner_directive_20261015.md#follow-up-decisions--2026-10-15)),
and served a frozen copy (tag `frozen-app-20261020`) from 2026-10-20 to 2026-10-22. Its
development-only service files are `deploy/islet-explorer-local.service`,
`deploy/shiny-server-local.conf` and `deploy/run-local-r.sh`; they are not the
production path described below. See
[the protein expansion record](../docs/protein_expansion_20261003.md) for
lineage and rollback to the frozen 15-donor root. The frozen 15-image release at
`5809be4` and its dated image, scan, and parity evidence were the reference
production release until 2026-10-22; the production candidate is now the
schema-5 release of the uniform root described in the next paragraph. See
[the September 28 handoff](../docs/HANDOFF_2026-09-28.md) and
[the readiness record](../docs/ufrc_readiness.md).

On 2026-10-03 the owner approved six donor boundary sets (HPAP-078, 6090,
6232, 6520, recovered 6539 and 6558) and the reconstruction phenotype rules,
and selected Harmony on protein PCA with the uncorrected fit retained. On
2026-10-15 the owner removed donor 78 (HDL078 / HPAP-078) from the app and all
analysis (`scripts/excluded_donors.csv`); the app refuses any data root that
still contains it. The 20 donors were then reprocessed through one workflow
(`scripts/uniform/`), reviewed by the owner against the frozen app, and made
the live local app on 2026-10-22 (root
`data/development/uniform-20donor-20261022-v10/`).

A uniform root is released under **data schema 5** (2026-10-22;
`scripts/build_release_data.R --schema 5`, `ISLET_DATA_SCHEMA=5`): the schema-1
Parquet summaries and spatial tables, plus the root's own files
(`uniform_release.json`, `endocrine_objects_by_donor.csv`,
`protein_phenotype_rules.csv`; `phenotype_rules.csv` only when the root has
one) and the H5AD's embedding records as `summary/embedding_metadata.json`,
which the production loader attaches to `prepared()`. The builder refuses to
package a uniform root as schema 1, because the app would then run with the
frozen 15-donor behaviour without saying so, and it copies the XML and offsets
sidecars beside each linked TIFF, so the release's `images/` is complete even
though the root links its images from other roots. The candidate built from
commit `d55680a` is staged under `tmp/ufrc-release-20261022/` (its data
bundle at `artifacts/release-d55680a/data/`) and has been the active release
`pubapps-2026-10-uniform-v10` on `pubmaigan1` since 2026-10-22; its
verification, parity, image and deployment evidence are recorded in
[the readiness record](../docs/ufrc_readiness.md).
Schema ids 2-4 belonged to the earlier all-donor experiments, whose runtimes and
capacity figures are historical evidence only; see
[the phase-2 contract](../docs/donor_expansion_phase2.md).

## Fixed production contract

- Rootless Podman and a user Quadlet service; Apptainer is not a supported
  production path.
- Code and dependencies baked into one immutable OCI image.
- A clean Git commit, immutable image reference, and checksummed data/image release
  are required before the service can start.
- `data/` is mounted read-only at `/srv/islet-explorer/data`.
- OME-TIFFs are mounted read-only directly at
  `/srv/islet-explorer/app/shiny_app/www/local_images`. This is deliberate:
  The explicit `site_dir` routes in `shiny-server.conf` serve `local_images/`,
  `images/`, and the legacy `img_data/` alias from that same mount with directory
  listings disabled. They are required for HTTP Range delivery: the general
  `app_dir` route forwards requests to R/httpuv and returned full-file HTTP 200
  responses in the local browser test. All three image routes remain covered
  by the UFRC proxy's application-wide authentication.
- Shiny listens on the container interface at 3838. Podman publishes one VM
  address only: the loopback `127.0.0.1:3838` by default, or, with
  `activate-release.sh --publish-address`, the private VM address UFRC's reverse
  proxy connects to (the proxy, `vpubweb.rc.ufl.edu`, is a separate host and
  cannot reach the VM loopback). The production activation uses
  `--publish-address 10.13.41.21`; `verify-deployment.sh` then requires that
  address to be the exclusive listener. The application relies on the proxy's
  Basic Auth, so the private address is acceptable only because UFRC's network
  filters it: measured on 2026-10-22, an HPG login node reaches the VM on port
  22 only (a test listener on `10.13.41.21:3839` received nothing), the campus
  network cannot reach the VM at all, and the VM's own firewalls are inactive.
- The container and Shiny workers use the same `shiny` UID/GID. Rootless
  `keep-id` maps the UFRC service account to that identity.
- AI, Statistics, and Inflammation are off. Production rejects a supplied
  `KEY` or `OPENAI_API_KEY`; no Navigator secret or AI egress is needed.
- The root filesystem, code, data, and images are read-only. Writable host
  mounts are limited to `runtime/logs` and `runtime/tmp`; `/run` and Shiny's
  ephemeral state use bounded tmpfs mounts.
- Resource guardrails are 2 CPUs, 14 GiB container memory, one numerical-library
  thread per R worker, and 2,048 PIDs, leaving host headroom on the base PubApps
  allocation. The CPU quota is rendered only where the user session has the
  `cpu` cgroup controller delegated; `pubmaigan1` delegates `memory` and `pids`
  only (2026-10-22), so there the 2-vCPU host size bounds CPU and activation says
  so.

The legacy `islet-explorer.def` is retained only as historical documentation and
must not be used for reviewer hosting.

## Deployment assets

| File | Purpose |
|---|---|
| `Containerfile` | Data-free production OCI image |
| `containerignore` | Whitelist build context; excludes data, images, Git history, credentials, and caches |
| `requirements-python.lock` | Exact tested Python/PAGA environment |
| `install-production-r-deps.R` | Mandatory, exact `destiny` and `duckdb` gates over the app's R manifest |
| `shiny-server.conf` | Sanitized `/islet-explorer/` Shiny service |
| `islet-explorer.container.in` | Hardened Quadlet template |
| `build-image.sh` / `scan-image.sh` | Clean build, dependency probes, SBOM, and High/Critical scan gate |
| `export-image.sh` / `load-image.sh` | Checksummed OCI archive transfer through HPG and identity-checked local load |
| `generate-release-manifest.sh` / `verify-release.sh` | Immutable release metadata, hashes, and full-verification stamp |
| `activate-release.sh` / `rollback-release.sh` | Atomic activation with health rollback and explicit prior-release rollback |
| `verify-deployment.sh` | Local health, UID, mount, dependency, network, and optional Range checks |
| `run.sh` | Manual Podman smoke launcher; not the persistent production service |

Core R packages are resolved by `scripts/install_shiny_deps.R` against the
base-image repository configuration. That inherited snapshot predates the tested
DuckDB 1.5.2 release, so the production installer obtains that exact CRAN source
archive and checks its pinned SHA-256 before installation. It bounds DuckDB
compilation to two jobs; `destiny` remains pinned to 3.22.0 through Bioconductor
3.21. The image build records every resolved R
package and version at `/opt/islet-explorer-build/R-packages.tsv`; the required
release scan writes ID-derived CycloneDX, Grype, R-package, Python-package, and
`.scan.txt` evidence plus `published-image.txt`. The scan attestation hashes the
four required per-image artifacts, and `published-image.txt` binds the local image
ID or registry digest to the source commit and base image. Changing the base digest
or any R dependency requires a new image and scientific parity run.

## 1. Project build, test, scan, and OCI export

Start from a clean, committed release revision; tag it if the project's release
practice calls for a tag. In particular, the PAGA R/Python files must be tracked;
the build script refuses dirty or untracked source.

Syft's `podman:` source reads images through the Podman API. In a dedicated build
shell, start a private Unix socket for the same local rootless image store before
running the build/scan commands. This does not open a TCP listener:

```bash
islet_socket_dir="$(mktemp -d "${XDG_RUNTIME_DIR:-/tmp}/islet-scan.XXXXXX")"
env -u CONTAINER_HOST podman system service --time=0 \
  "unix://${islet_socket_dir}/podman.sock" > "${islet_socket_dir}/service.log" 2>&1 &
islet_socket_pid=$!
cleanup_islet_scan_socket() {
  kill "${islet_socket_pid}" 2>/dev/null || true
  wait "${islet_socket_pid}" 2>/dev/null || true
  unset CONTAINER_HOST
  rm -f -- "${islet_socket_dir}/podman.sock" "${islet_socket_dir}/service.log"
  rmdir -- "${islet_socket_dir}"
}
trap cleanup_islet_scan_socket EXIT
for _ in {1..50}; do
  [[ -S "${islet_socket_dir}/podman.sock" ]] && break
  sleep 0.2
done
test -S "${islet_socket_dir}/podman.sock"
chmod 0600 "${islet_socket_dir}/podman.sock"
export CONTAINER_HOST="unix://${islet_socket_dir}/podman.sock"
curl --fail --silent --show-error --unix-socket "${islet_socket_dir}/podman.sock" \
  http://localhost/_ping
```

The temporary directory is private to the build user. Keep this shell open through
scan and export, then exit it to stop the API and remove the socket. If the ping
fails, inspect `service.log` before retrying; a working Podman CLI alone does not
establish scanner connectivity. The separate `PubApps Image Build and Security
Gate` GitHub workflow supplies its own isolated socket and retains build/security
evidence on failures. It exports an OCI artifact only after the security gate
passes and does not deploy or push an image.

Resolve the chosen `rocker/shiny:4.5.1` base to a registry digest and then run,
from the repository root:

```bash
BASE_IMAGE='docker.io/rocker/shiny@sha256:<verified-base-digest>'
RELEASE='pubapps-2026-01'
IMAGE_TAG="localhost/islet-explorer:${RELEASE}"
EVIDENCE_DIR="/staging/islet-explorer-evidence/${RELEASE}"
OCI_ARCHIVE="/staging/islet-explorer-${RELEASE}.oci.tar"

deploy/build-image.sh \
  --base-image "${BASE_IMAGE}" \
  --version "${RELEASE}" \
  --tag "${IMAGE_TAG}" \
  --scan-dir "${EVIDENCE_DIR}"

# Copy the immutable ID printed by build-image.sh; do not export by tag.
IMAGE_ID='sha256:<image-id-printed-by-build-image>'
deploy/export-image.sh "${IMAGE_ID}" "${OCI_ARCHIVE}"
```

Release builds require `podman`, `syft`, `grype`, and `python3`. The build performs Python
dependency consistency checks and imports, verifies exact `destiny`/`duckdb`
versions, parses all application R files, creates the SBOM and resolved dependency
manifests, and fails on a High
or Critical vulnerability. `export-image.sh` creates an OCI archive, metadata,
and a SHA-256 file covering both from the resolved immutable image ID. The
mutable build tag is never placed in release metadata.

Use the project workstation as the build host when its rootless `podman info`
check succeeds. Install Podman and its rootless UID/network helpers plus Syft and
Grype there before building, and retain the tool versions with the evidence.
The wrapper first probes the digest-selected base for x86_64, R 4.5.1, the `shiny`
account, and Shiny Server before the long dependency installation. Rscript is
resolved through the image's PATH because Rocker's R is installed under
`/usr/local`. The image is built locally; the PubApps VM receives its OCI archive.

UFIT standards require a maintained SBOM, lint/static-analysis evidence, and
pre-deployment code/package/container vulnerability testing, with remaining CVEs
documented. This project's image policy is stricter: Grype `--fail-on high` must
report no High or Critical image finding. Supply the exact-candidate evidence to
UFRC for review; UFRC's production review and ongoing network scans are a separate
requirement and do not replace this project gate.

The scan retains the original Syft SBOM, its scan copy, and a hash-bound correction
record. `correct-syft-vim-cpe.py` removes only Syft 1.51.1's exact editor CPE on
`pkg:cran/VIM@6.2.6` when the R cataloger, metadata, and DESCRIPTION location all
match. [CRAN](https://cran.r-project.org/package=VIM) and the
[package source](https://github.com/statistikat/VIM) identify this as the R
missing-value package. Every component and other field remains; Grype runs with
`--add-cpes-if-none=false` and the unchanged High/Critical gate. Genuine editor
components and other versions receive no correction.

Do not put scientific data, `.Renviron`, credentials, or the 109 GiB image tree
in the OCI archive. Inspect the external evidence directory and retain it with
the archive and publication release record; do not create build evidence inside
the clean Git checkout.

The primary production topology has no registry dependency: transfer the OCI
archive and its `.metadata` and `.sha256` sidecars from the build host to HPG,
then from HPG to the PubApps VM. A GHCR image pulled by immutable digest may be
used only if UFRC explicitly approves registry egress; `build-image.sh --push`
supports that alternative.

## 2. Build and stage an immutable data release

On the data-build host, generate the one authoritative runtime bundle. Its
`images/` subtree is part of that bundle and its `manifest.json` binds the data
schema, the source commit and the OME-TIFF count. The output directory must not
already exist. A uniform 20-donor root takes `--schema 5`; the frozen 15-donor
root under `data/` took `--schema 1`. The builder checks the root against the
schema and refuses the wrong one. Set `RETICULATE_PYTHON` to the app's
interpreter first (the H5AD is read through it).

```bash
Rscript scripts/build_release_data.R \
  --source-root data/development/uniform-20donor-20261022-v10 \
  --image-root data/development/uniform-20donor-20261022-v10/images \
  --output /staging/islet-explorer-data-pubapps-2026-01 \
  --materialize copy \
  --schema 5

Rscript scripts/verify_release_data.R \
  --root /staging/islet-explorer-data-pubapps-2026-01 \
  --mode full --schema 5

Rscript scripts/verify_release_parity.R \
  --source-root data/development/uniform-20donor-20261022-v10 \
  --image-root data/development/uniform-20donor-20261022-v10/images \
  --release-root /staging/islet-explorer-data-pubapps-2026-01 \
  --output-prefix /staging/islet-explorer-evidence/pubapps-2026-01/release-parity
```

Use `--materialize copy` for the reviewer bundle so subsequent edits to canonical
inputs cannot change the release through shared inodes. Hardlinks are only a
space-saving option for disposable local tests.
The builder and image must come from the same clean Git commit; activation
compares the OCI revision label, generated data manifest, and release metadata
and refuses a mismatch. The parity check also requires that commit to be the
checked-out, clean one, so run the builder and the checks from a clean checkout
or linked worktree of the release commit.

The VM was provisioned on 2026-10-22: `pubmaigan1.ufhpc`, service account
`maigan` (home `/pubapps/maigan` on the 1 TiB Qumulo share; RHEL 9.8; rootless
Podman 5.8.2 with its image store on the 40 GB local disk at
`/local/maigan/podman/storage`; cgroups v2; linger enabled; 2 vCPU, 15 GiB).
It is reached from an HPG login node as `ssh pubmaigan1` (HPG's `maigan@hpg` key
is authorized there), so the workstation stages through its HPG session with
`rsync -e 'ssh hpg ssh' ... pubmaigan1:...`, at about 100 MB/s. The project path
is `/pubapps/maigan/islet-explorer/`; the `deploy/` tooling of the release
commit is copied to `tooling/deploy/` beside it.

```text
/pubapps/maigan/islet-explorer/
  releases/
    pubapps-2026-01/
      data/                  # output of the release-data builder
        images/              # declared OME-TIFF count + XML/offset sidecars
      attestations/          # SBOM, scan, dependency manifests
      release.env            # generated; do not hand-edit
      MANIFEST.sha256        # generated
      .verified              # generated after full hash verification
  runtime/
    logs/
    tmp/
    deployments/             # successful active/previous state
```

Stage with resumable `rsync` through HPG; PubApps is not mounted on HPG and the
VM's SSH route is available from HPG. Transfer only the runtime whitelist, not
the repository, `.git` history, backups, caches, or upstream/raw trees.

```bash
PROJECT_ROOT='/pubapps/maigan/islet-explorer'
RELEASE='pubapps-2026-10-uniform-v10'
RELEASE_DIR="${PROJECT_ROOT}/releases/${RELEASE}"
HOP="ssh hpg ssh"   # the workstation's HPG control session, then the VM

ssh hpg ssh pubmaigan1 mkdir -p "${RELEASE_DIR}/data" "${RELEASE_DIR}/attestations" \
  "${PROJECT_ROOT}/incoming"

# From the workstation, substituting the staged source paths.
rsync -a --partial --info=progress2 -e "${HOP}" <generated-release-data>/ \
  pubmaigan1:"${RELEASE_DIR}/data/"
rsync -a -e "${HOP}" /staging/islet-explorer-evidence/"${RELEASE}"/ \
  pubmaigan1:"${RELEASE_DIR}/attestations/"
rsync -a -e "${HOP}" /staging/islet-explorer-"${RELEASE}".oci.tar{,.metadata,.sha256} \
  pubmaigan1:"${PROJECT_ROOT}/incoming/"
```

Copy this `deploy/` tooling directory to the VM separately. On the VM, generate
the release metadata with the exact immutable image reference and source commit, then
perform the initial full verification. Hashing the TIFF store is intentionally slow:
manifest generation hashes it once, this explicit verification rehashes it and writes
the stamp, and activation rehashes it immediately before service replacement. Ordinary
service restarts do not repeat the hash pass.

```bash
deploy/load-image.sh \
  "${PROJECT_ROOT}/incoming/islet-explorer-${RELEASE}.oci.tar"

IMAGE_REF='sha256:<image-id-printed-by-load-image>'
GIT_COMMIT='<40-character-release-commit>'
EXPECTED_IMAGES=20
DATA_SCHEMA=5

deploy/generate-release-manifest.sh \
  "${RELEASE_DIR}" "${RELEASE}" "${IMAGE_REF}" "${GIT_COMMIT}" "${DATA_SCHEMA}" "${EXPECTED_IMAGES}"
deploy/verify-release.sh "${RELEASE_DIR}"
```

This example is for an exact, accepted 20-image schema-5 release of the
uniform root, built from a clean revision and matched to its immutable image
(the frozen 15-image schema-1 release took `15` and `1`). Do not substitute a
development route or a rejected bundle for the approved image source. Stage
only declared runtime artifacts; do not deploy raw working directories,
preview logs, or rejected images.

Generation hashes the builder's `data/` output once, including its `images/`
subtree, and rejects symlinks that escape the release plus filenames associated
with credentials or mutable backups. Supply the independently reviewed OME-TIFF
count as the sixth argument; generation rejects a different staged count and
writes `MANIFEST_VERSION=2` and `EXPECTED_IMAGES` into the hashed `release.env`.
The preserved 15-donor runtime uses data schema `1`; the uniform 20-donor root
uses schema `5`, and `verify-release.sh` checks that a schema-5 data manifest
declares the same image count and every uniform-root artifact. Schema ids 2-4
were the retired all-donor experiments and never an accepted contract. The
explicit image declaration here binds the selected cohort; omitting the count
retains the generator's original default of 15.

Verification reads that declaration automatically, checks every manifest hash,
and writes `.verified` tied to the manifest hash. Activation, the manual launcher,
and rollback all use the same verifier; they need no image-count override.
Legacy `MANIFEST_VERSION=1` releases without `EXPECTED_IMAGES` retain their
15-image contract, including the frozen original release. For an optional
operator cross-check, `--expected-images N` must agree with the release metadata;
it cannot override it. Missing version-2 counts and duplicate, malformed, or
conflicting metadata fail verification. Even `--quick` checks the metadata hash
and actual image count, but does not rehash TIFF contents. If the bundle changes,
regenerate and fully verify it.

## 3. Activate the Quadlet service

On the provisioned host, first verify that the available Podman version accepts
`keep-id:uid=<shiny_uid>,gid=<shiny_gid>` and that login linger is enabled for
the service account. Quadlet/user-systemd is the implementation selected by this
project, not a separate UFRC approval gate. Activation performs a
bind-mount/writeability probe with the resolved image UID before replacing the
service.

```bash
loginctl show-user "${USER}" -p Linger
deploy/activate-release.sh \
  --project-root '/pubapps/<project>/islet-explorer' \
  --release 'pubapps-2026-01'
```

The command repeats the full content-hash gate immediately before replacement,
renders `~/.config/containers/systemd/islet-explorer.container`, reloads user
systemd, starts `islet-explorer.service`, and runs local verification. If the
new service fails to start or become healthy, it restores the prior unit
automatically. It records active/previous state only after a successful health
check. Recovery checks the restored service's health and immutable image identity.
If service replacement fails, exit 70 reports completed recovery and exit 71
reports that recovery could not be verified and requires operator attention.
Preflight failures can also use exit 70 before any service is replaced.

Useful operations:

```bash
systemctl --user status islet-explorer.service
journalctl --user -u islet-explorer.service -f
tail -F /pubapps/<project>/islet-explorer/runtime/logs/*.log

deploy/rollback-release.sh \
  --project-root '/pubapps/<project>/islet-explorer'
```

Never edit an active release or point the service at a mutable image tag. Build
and verify a new release, then activate it. Keep at least the active and previous
image/data pairs until the new release has survived a VM reboot test.

## 4. UFRC reverse-proxy and Basic-Auth acceptance

Request public path `/islet-explorer/` and the backend. The proxy host cannot
reach the VM loopback, so the backend is the VM's private address (`pubmaigan1`,
`10.13.41.21:3838`), published with `--publish-address 10.13.41.21` since
2026-10-22 21:16 UTC (the publish-address rule above records the measured
network filtering that makes this acceptable). Ask UFRC to preserve WebSocket upgrades,
forwarded host/protocol, long Shiny connections, and unbuffered HTTP byte ranges for
multi-gigabyte TIFFs. Basic Auth must cover the page, WebSocket endpoint, downloads,
TIFFs, XML, JSON, and other served static assets. Parquet, GeoPackage, H5AD, and
other release-data paths must not be directly web-accessible.

Run locally on the VM, substituting one TIFF relative path:

```bash
deploy/verify-deployment.sh \
  --url http://10.13.41.21:3838/islet-explorer/ \
  --container islet-explorer \
  --publish-address 10.13.41.21 \
  --range-path 'local_images/<image>.ome.tiff'
```

The verifier fails if the backend port has any listener outside its expected
binding (`127.0.0.1`, or the `--publish-address` given), including a wildcard
or any other host IP, if data or
images are writable, if feature flags or credentials violate the reviewer
profile, if dependencies are absent, or if the byte request is not exactly a
`206 Content-Range` response. Activation also checks that the running container
uses the requested immutable image. The manual `run.sh` launcher applies the
same verification when an alternate loopback `PORT` is selected for a smoke test.

Before reviewer credentials are released, use a clean non-UF browser/network to
verify the required route behavior:

1. An unauthenticated request to the site and to a TIFF returns `401`.
2. The neutral reviewer credential unlocks all four tabs and their CSV exports.
3. Plot, Trajectory/PAGA, Viewer, and Spatial match the frozen manuscript data.
4. The Viewer streams rather than downloads a complete TIFF and the Shiny
   WebSocket stays connected through long interactions.
5. Directory listing, direct release-data access, and access outside the application
   path are denied.

As project operational checks, exercise a representative five-session smoke at the
declared peak and confirm that `systemctl --user` brings the same immutable image
reference back after a controlled VM reboot. These checks inform capacity and
recovery; they are not additional UFRC policy gates and have no invented latency
threshold.

Send the shared reviewer credential only through the journal's confidential
field, rotate it after each review round, and never store it in this repository,
the image, release data, shell history, or publicly visible ticket comments.

## 5. Release evidence and operations

Retain the following standards and approval evidence per reviewer release:

- Lint/static-analysis results, documented input validation, a maintained SBOM,
  and project code/package/container vulnerability results. The project image gate
  has no High or Critical findings; any other remaining CVE is documented.
- Appropriate industry-standard one-event-per-line application logs sent to a
  UFIT-monitored location.
- Accessibility/UF-branding review or remediation record.
- UFRC security/availability approval and PI/data-owner attestation.

Also retain the project's reproducibility and recovery record:

- The 40-character source commit and release tag if one was created.
- Immutable local image ID (or UFRC-approved registry digest), plus the OCI
  archive, transfer metadata, and SHA-256 sidecar.
- Base-image digest.
- `release.env`, `MANIFEST.sha256`, and verification record.
- The ID-derived CycloneDX SBOM, Grype JSON report, resolved R package TSV,
  resolved Python package list, `.scan.txt` attestation, and
  `published-image.txt`; retain the optional Grype database-status JSON when
  generated.
- Application/unit/scientific-parity and Range results, plus representative
  concurrency, restart, and rollback notes.
- Project log-minimization checks and the recorded UFRC answers about log access and
  retention.

Monitor HTTPS health, disconnects, CPU/memory, disk usage, and log growth. Logs
must use an appropriate one-event-per-line format and must not contain exported
scientific rows, credentials, or unnecessary Sensitive content. Follow the
UFRC-provided collection, access, and retention rules; the detailed inventory adds
project recommendations for redaction, rotation, and leakage checks. If UFRC confirms
that project-plan item 6 applies to this Sensitive deployment, review Basic Auth and
SSH access quarterly while peer review remains active. Renew the PubApps allocation
before expiration. Remove Basic Auth only after written PI/data-owner/nPOD and UFRC
approval of the exact public-release bundle.
