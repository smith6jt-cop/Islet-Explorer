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

## Fixed production contract

- Rootless Podman and a user Quadlet service; Apptainer is not a supported
  production path.
- Code and dependencies baked into one immutable OCI image.
- A clean Git commit, immutable image digest, and checksummed data/image release
  are required before the service can start.
- `data/` is mounted read-only at `/srv/islet-explorer/data`.
- OME-TIFFs are mounted read-only directly at
  `/srv/islet-explorer/app/shiny_app/www/local_images`. This is deliberate:
  Shiny Server's static handler must see them beneath `www/` to preserve HTTP
  Range delivery through the UFRC proxy.
- Shiny listens on the container interface at 3838; Podman publishes only
  `127.0.0.1:3838` on the VM.
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
base-image repository configuration. The image build records every resolved R
package and version at `/opt/islet-explorer-build/R-packages.tsv`; the required
release scan exports that file alongside the committed Python lock, CycloneDX
SBOM, Grype report, and hashes. The OCI digest binds that resolved environment.
Changing the base digest or any R dependency requires a new image and scientific
parity run.

## 1. Build, test, scan, and export the OCI image

Start from a clean, committed release tag. In particular, the PAGA R/Python
files must be tracked; the build script refuses dirty or untracked source.

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

deploy/export-image.sh "${IMAGE_TAG}" "${OCI_ARCHIVE}"
```

Release builds require `podman`, `syft`, and `grype`. The build performs Python
imports, verifies exact `destiny`/`duckdb` versions, parses all application R
files, creates the SBOM and resolved dependency manifests, and fails on a High
or Critical vulnerability. `export-image.sh` creates an OCI archive, metadata,
and a SHA-256 file covering both. Record the immutable local `sha256:...` image
ID; the mutable build tag is never placed in release metadata.

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
  --materialize hardlink \
  --schema 1

Rscript scripts/verify_release_data.R \
  --root /staging/islet-explorer-data-pubapps-2026-01 \
  --mode full
```

Use `--materialize copy` if the staging destination is on another filesystem.
The builder and image must come from the same clean Git commit; activation
compares the OCI revision label, generated data manifest, and release metadata
and refuses a mismatch.

On the PubApps VM, use the assigned project path:

```text
/pubapps/<project>/islet-explorer/
  releases/
    pubapps-2026-01/
      data/                  # output of the release-data builder
        images/              # 15 OME-TIFFs + XML/offset sidecars
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
the release metadata with the exact image digest and source commit, then perform
the full hash pass. Hashing the TIFF store is intentionally slow and happens
once before activation, not at every restart.

```bash
deploy/load-image.sh \
  "${PROJECT_ROOT}/incoming/islet-explorer-${RELEASE}.oci.tar"

IMAGE_REF='sha256:<image-id-printed-by-load-image>'
GIT_COMMIT='<40-character-release-commit>'

deploy/generate-release-manifest.sh \
  "${RELEASE_DIR}" "${RELEASE}" "${IMAGE_REF}" "${GIT_COMMIT}" 1
deploy/verify-release.sh "${RELEASE_DIR}"
```

Generation hashes the builder's `data/` output once, including its `images/`
subtree, and rejects symlinks that escape the release plus filenames associated
with credentials or mutable backups. Verification requires exactly 15
OME-TIFFs by default, checks every manifest hash, and writes `.verified` tied to
the manifest hash. If the bundle changes, regenerate and fully verify it.

## 3. Activate the Quadlet service

UFRC must first confirm that its Podman version accepts
`keep-id:uid=<shiny_uid>,gid=<shiny_gid>` and that login linger is enabled for
the service account. Activation itself performs a bind-mount/writeability probe
with the resolved image UID before replacing the service.

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
check.

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

Give UFRC backend `127.0.0.1:3838` and public path `/islet-explorer/`. Ask it to
preserve WebSocket upgrades, forwarded host/protocol, long Shiny connections,
and unbuffered HTTP byte ranges for multi-gigabyte TIFFs. Basic Auth must cover
the page, WebSocket endpoint, downloads, TIFFs, XML, JSON, and any Parquet,
GeoPackage, or H5AD paths.

Run locally on the VM, substituting one TIFF relative path:

```bash
deploy/verify-deployment.sh \
  --url http://127.0.0.1:3838/islet-explorer/ \
  --container islet-explorer \
  --range-path 'local_images/<image>.ome.tiff'
```

The verifier fails if host port 3838 is bound to `0.0.0.0`/IPv6-any, if data or
images are writable, if feature flags or credentials violate the reviewer
profile, if dependencies are absent, or if the byte request is not exactly a
`206 Content-Range` response.

From a clean non-UF browser/network, additionally verify:

1. An unauthenticated request to the site and to a TIFF returns `401`.
2. The neutral reviewer credential unlocks all four tabs and their CSV exports.
3. Plot, Trajectory/PAGA, Viewer, and Spatial match the frozen manuscript data.
4. The Viewer streams rather than downloads a complete TIFF and the Shiny
   WebSocket stays connected through long interactions.
5. Directory listing and direct access outside the application path are denied.
6. Five simultaneous reviewer sessions meet the documented latency/memory gate.
7. `systemctl --user` brings the same digest back after a controlled VM reboot.

Send the shared reviewer credential only through the journal's confidential
field, rotate it after each review round, and never store it in this repository,
the image, release data, shell history, or publicly visible ticket comments.

## 5. Release evidence and operations

Retain, per reviewer release:

- Signed Git tag and 40-character source commit.
- Immutable local image ID (or UFRC-approved registry digest), plus the OCI
  archive, transfer metadata, and SHA-256 sidecar.
- Base-image digest.
- `release.env`, `MANIFEST.sha256`, and verification record.
- CycloneDX SBOM, Grype JSON/attestation, resolved R package TSV, and resolved
  Python package list.
- Application/unit/scientific-parity, accessibility, concurrency, Range,
  restart, and rollback test results.
- UFRC security approval and PI/data-owner attestation.

Monitor HTTPS health, disconnects, CPU/memory, disk usage, and log rotation.
Logs must not contain exported scientific rows, passwords, or other sensitive
content. Review Basic Auth and SSH access quarterly while peer review remains
active. Renew the PubApps allocation before expiration. Remove Basic Auth only
after written PI/data-owner/nPOD approval of the exact public-release bundle.
