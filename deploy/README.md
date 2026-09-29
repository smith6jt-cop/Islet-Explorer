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

The current local development app serves the restored original 15 donors and
5,214 islets at `http://10.15.152.7:8080/islet-explorer/` through loopback
port 3840. Its `/local_images/` route still points to an earlier 21-image
collection under an explicit development-only allowance; this is not a
production image/data pairing. The frozen 15-image release at `5809be4` and
its dated image, scan, and parity evidence remain the reference release, not
evidence for the current development UI. See [the September 28 handoff](../docs/HANDOFF_2026-09-28.md)
and [the readiness record](../docs/ufrc_readiness.md).

Six donor additions—HPAP-078, 6090, 6232, 6520, recovered 6539, and
6558—are staged separately in a raw-only offline bundle. Its 943,068 objects
and 852 passing original-preservation checks do not make it app-ready:
`app_ready` and `scientific_acceptance` are false. Five donors still need
boundary review. The exact original phenotype rules, numeric gates, scVI
weights and preprocessing have not been recovered and verified for their
typing; a common 21-donor trajectory/UMAP/Leiden refit remains unvalidated.
There is no approved combined islet count or 21-image production release. The
earlier schema-2/3/4 all-donor runtimes and their capacity figures are historical experimental
evidence; see [the phase-2 contract](../docs/donor_expansion_phase2.md).
Production verification rejects unaccepted expanded manifests, and the
raw-only bundle has no runtime manifest or app loader.

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
- Shiny listens on the container interface at 3838. The project requests Podman
  publication on VM loopback at `127.0.0.1:3838`, pending UFRC confirmation that
  the proxy can reach that target.
- The container and Shiny workers use the same `shiny` UID/GID. Rootless
  `keep-id` maps the UFRC service account to that identity.
- AI, Statistics, and Inflammation are off. Production rejects a supplied
  `KEY` or `OPENAI_API_KEY`; no Navigator secret or AI egress is needed.
- The root filesystem, code, data, and images are read-only. Writable host
  mounts are limited to `runtime/logs` and `runtime/tmp`; `/run` and Shiny's
  ephemeral state use bounded tmpfs mounts.
- Resource guardrails are 2 CPUs, 14 GiB container memory, one numerical-library
  thread per R worker, and 2,048 PIDs, leaving host headroom on the base PubApps
  allocation.

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
schema and source commit. The output directory must not already exist.

```bash
Rscript scripts/build_release_data.R \
  --source-root data \
  --image-root app/shiny_app/www/local_images \
  --output /staging/islet-explorer-data-pubapps-2026-01 \
  --materialize copy \
  --schema 1

Rscript scripts/verify_release_data.R \
  --root /staging/islet-explorer-data-pubapps-2026-01 \
  --mode full
```

Use `--materialize copy` for the reviewer bundle so subsequent edits to canonical
inputs cannot change the release through shared inodes. Hardlinks are only a
space-saving option for disposable local tests.
The builder and image must come from the same clean Git commit; activation
compares the OCI revision label, generated data manifest, and release metadata
and refuses a mismatch.

On the PubApps VM, use the assigned project path:

```text
/pubapps/<project>/islet-explorer/
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

Stage from HiPerGator with resumable `rsync`; PubApps is not mounted on HPG and
its VM SSH route is available from HPG. Transfer only the runtime whitelist, not
the repository, `.git` history, backups, caches, or upstream/raw trees.

```bash
PROJECT_ROOT='/pubapps/<project>/islet-explorer'
RELEASE='pubapps-2026-01'
RELEASE_DIR="${PROJECT_ROOT}/releases/${RELEASE}"

mkdir -p "${RELEASE_DIR}/data" "${RELEASE_DIR}/attestations" \
  "${PROJECT_ROOT}/incoming"

# Run these from HPG, substituting the assigned VM and staged source paths.
rsync -a --partial --info=progress2 <generated-release-data>/ \
  <service-user>@<vm>:"${RELEASE_DIR}/data/"
rsync -a /staging/islet-explorer-evidence/"${RELEASE}"/ \
  <service-user>@<vm>:"${RELEASE_DIR}/attestations/"
rsync -a /staging/islet-explorer-"${RELEASE}".oci.tar{,.metadata,.sha256} \
  <service-user>@<vm>:/pubapps/<project>/islet-explorer/incoming/
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
EXPECTED_IMAGES=15
DATA_SCHEMA=1

deploy/generate-release-manifest.sh \
  "${RELEASE_DIR}" "${RELEASE}" "${IMAGE_REF}" "${GIT_COMMIT}" "${DATA_SCHEMA}" "${EXPECTED_IMAGES}"
deploy/verify-release.sh "${RELEASE_DIR}"
```

This example is for an exact, accepted 15-image release built from a clean
revision and matched to its immutable image. Do not substitute the local
21-image development route or raw-only additions bundle for the approved
image source. Stage only declared runtime artifacts; do not deploy raw working
directories, preview logs, or rejected images.

Generation hashes the builder's `data/` output once, including its `images/`
subtree, and rejects symlinks that escape the release plus filenames associated
with credentials or mutable backups. Supply the independently reviewed OME-TIFF
count as the sixth argument; generation rejects a different staged count and
writes `MANIFEST_VERSION=2` and `EXPECTED_IMAGES` into the hashed `release.env`.
The preserved 15-donor runtime uses data schema `1`. Later schema-2/3/4
experiments do not establish an accepted production contract for the six
additions. The explicit 15-image declaration here binds the selected cohort;
omitting the count retains the generator's original default of 15.

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

Request public path `/islet-explorer/` and backend `127.0.0.1:3838`. UFRC must
confirm that its proxy can reach VM loopback or return the required private backend
address and agreed port. The current Quadlet, activation, and verification tooling
assumes `127.0.0.1:3838`; if UFRC assigns a different target, do not activate until the
project updates and re-reviews that configuration and tooling and repeats its tests.
Ask UFRC to preserve WebSocket upgrades,
forwarded host/protocol, long Shiny connections, and unbuffered HTTP byte ranges for
multi-gigabyte TIFFs. Basic Auth must cover the page, WebSocket endpoint, downloads,
TIFFs, XML, JSON, and other served static assets. Parquet, GeoPackage, H5AD, and
other release-data paths must not be directly web-accessible.

Run locally on the VM, substituting one TIFF relative path:

```bash
deploy/verify-deployment.sh \
  --url http://127.0.0.1:3838/islet-explorer/ \
  --container islet-explorer \
  --range-path 'local_images/<image>.ome.tiff'
```

The verifier fails if the backend port has any listener outside its expected
`127.0.0.1` binding, including a specific non-loopback host IP, if data or
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
