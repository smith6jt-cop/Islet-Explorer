# Islet Explorer on UFRC PubApps

**Current status, 2026-11-03:** production `pubapps-2026-11-components-v1`
and the separate strict-DuckDB dev release `pubapps-2026-11-dev-duckdb-v1` are
active and verified. Their exact source, image, data, checks and service state
are bound by the [production receipt](../docs/evidence/pubapps_live_20261102.json)
and [dev receipt](../docs/evidence/dev_duckdb_deployment_20261103.json).
Production retains its dated DuckDB window through November 10 UTC; dev passes
the strict decoder gate without a deferral. The older
`islet-explorer-local.service` remains inactive.
Use the [workflow guide](../docs/workflow_user_guide.md) and
[current v2 build record](../docs/component_build_v2_20261102.md). Earlier masked-service,
whole-v19 and v1 diagnostic records are historical. The [November 3 display fixes](../docs/app_display_fixes_20261103.md)
are running in a separate transient [local review](http://127.0.0.1:3845/islet-explorer/)
from source `7ac421e`, image `46b81d66…`, using the separately accepted data built
from `5e6580f`; see the [local review record](../docs/local_duckdb_review_20261103.md).
The owner's November 3 direction is **LOCAL REVIEW ONLY**: dev activation is
paused pending the owner's review and decision. Both public services remain unchanged.

Production runbook for the journal-review deployment requested at:

`https://islet-explorer.rc.ufl.edu/islet-explorer/`

UFRC provisioned this hostname/path on 2026-10-23. The production receipt above
records successful activation and checks; the address alone does not establish
scientific acceptance.

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
(loopback port 3840) previously served the main working tree on the reprocessed uniform root
`data/development/uniform-20donor-20261022-v10/` from 2026-10-22 (owner: "make the
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
frozen 15-donor behavior without saying so, and it copies the XML and offsets
sidecars beside each linked TIFF, so the release's `images/` is complete even
though the root links its images from other roots. The candidate built from
commit `d55680a` is staged under `tmp/ufrc-release-20261022/` (its data
bundle at `artifacts/release-d55680a/data/`) and was the active release
`pubapps-2026-10-uniform-v10` on `pubmaigan1` from 2026-10-22 to 2026-10-26, when
the owner promoted `pubapps-2026-10-uniform-v10-perf2` (commit `0b55429`, Plan A
Levels 1-2) from the staging instance to production; the earlier release stays on the
VM for `rollback-release.sh`. Verification, parity, image and deployment evidence of
every release are recorded in [the readiness record](../docs/ufrc_readiness.md).
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
  mounts are limited to `runtime/logs` and `runtime/tmp` (for the staging
  instance, `runtime/islet-explorer-dev/logs` and `.../tmp`); `/run` and Shiny's
  ephemeral state use bounded tmpfs mounts.
- Resource guardrails are 2 CPUs, a memory cap per instance (9 GiB production,
  4 GiB staging; `activate-release.sh --memory`), one numerical-library thread
  per R worker, and 2,048 PIDs, leaving the 15 GiB VM its own headroom when
  both instances run. The CPU quota is rendered only where the user session has the
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
| `verify-brotli.py` | Required actual decoder regression, package/library evidence, and system/Pillow linkage gate |
| `tests/check-duckdb-brotli.R` | Optional future real release Parquet BROTLI round-trip and query compatibility check with read-only input |
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
archive and checks its pinned SHA-256 before installation, with compilation
bounded to two jobs. Production retains the original decoder under the dated
November 10 remediation window below; the separate dev branch fixes its decoder
and passes strict validation, with its source provenance in the dev receipt.
`destiny` remains pinned to 3.22.0 through Bioconductor 3.21.
The working branch now contains the exact reviewed DuckDB patch from `378843e`
(cherry-pick `b5fae5c`), so new builds use the strict fixed decoder. This code
integration does not change the older image currently serving production.
The image build records every resolved R
package and version at `/opt/islet-explorer-build/R-packages.tsv`; the required
release scan writes ID-derived CycloneDX, Grype, R-package, Python-package, and
supplementary Brotli JSON evidence, plus `.scan.txt` and `published-image.txt`.
The scan attestation hashes the required evidence artifacts, and
`published-image.txt` binds the local image
ID or registry digest to the source commit and base image. Changing the base digest
or any R dependency requires a new image and scientific parity run.

### Brotli source and runtime verification

The 2026-11-02 security review found copies of the affected decoder outside the
system package inventory: Pillow's wheel bundles a private library, and DuckDB
statically includes one in `duckdb.so`. Updating Ubuntu packages and obtaining
a zero-High/Critical package scan did not fix or detect both copies. The Python
lock keeps Pillow **13.0.0** but builds it from source against the patched
Ubuntu-managed FreeType/Brotli libraries, with two compilation workers.

**Production decision, 2026-11-02:** the owner directed the then-current build to the live PubApp after its real data/release/runtime checks, within the notice's DuckDB patch window through **2026-11-10 UTC**, and directed development/testing of our DuckDB fix separately for dev. See [the dated decision](../docs/duckdb_remediation_window_20261102.md). Production retains pinned CRAN DuckDB **1.5.2** (archive SHA-256 `55f98a336513a94354f328791c8b94521b9f2e918a57d7e9fe6c0cb6b9a38db6`); the historical backport at `3a12daa` is preserved separately. The current isolated dev fix is deployed and strictly validated in the [dev receipt](../docs/evidence/dev_duckdb_deployment_20261103.json). The November 3 display candidate is local review only under the direction at the top of this runbook.

`scan-image.sh` runs this checkout's `verify-brotli.py` inside the immutable image, with no network or data mounts. The original production decoder has the known result `0`/error `-9`. `evaluate-brotli-report.py` accepts that single complete finding during the recorded window and emits `RESULT=pass-with-duckdb-remediation-window`, with the raw failure, authorization copy and deadline bound by hash. New fixed builds must return `0`/error `-17` and pass the reviewed source-provenance check without a deferral. Pillow/Ubuntu/Node failures, unknown decoder behavior, crashes, timeouts, missing evidence and other High/Critical findings still block. No general CVE suppression is added.

`generate-release-manifest.sh` verifies those hashes, the exact authorization and its date. It records the decision/deadline in `release.env`; `verify-release.sh`, including `--quick`, rejects an invalid or expired deferral. Copy the updated release tooling when staging the new release. The raw finding remains unresolved in production until the separately validated dev fix is promoted. Old strict-gate evidence stays unchanged; rescan the exact current candidate into a fresh evidence directory. No extra approval is needed for the already-authorized DuckDB window. Scientific/data acceptance and runtime checks remain required.

For a future candidate, the valid-data check below is available. Its successful
runs on the completed patched dev image are recorded in the [dev fix evidence](../docs/dev_duckdb_fix_20261102.md). Set `IMAGE_ID` to its immutable
`sha256:` ID, `REAL_PARQUET` to an absolute path for a representative nonempty
release Parquet containing text and numeric columns, and `PARQUET_REPORT` to a
fresh report path outside the release. From the repository root:

```bash
podman run --rm --pull=never --network none --read-only \
  --cap-drop=ALL --security-opt=no-new-privileges \
  --tmpfs /tmp:rw,size=1g,mode=1777 \
  --volume "${REAL_PARQUET}:/input/islets.parquet:ro" \
  --volume "${PWD}/deploy/tests/check-duckdb-brotli.R:/checks/check-duckdb-brotli.R:ro" \
  --entrypoint Rscript "${IMAGE_ID}" \
  --vanilla /checks/check-duckdb-brotli.R /input/islets.parquet \
  > "${PARQUET_REPORT}"
```

The script verifies BROTLI compression in every generated column chunk, exact
schema and ordered values, missingness, duplicate counts using bidirectional
`EXCEPT ALL`, numeric and grouped aggregates, and an unchanged input SHA-256.
Only the container's temporary filesystem receives the round-trip file; the
scientific input is mounted read-only. Check the command exit status and report
before recording success. This is separate from the malformed-stream probe
and does not replace full release scientific parity or browser checks.

**Historical early November 2 build status:** the combined backport build at `3a12daa` was stopped following the owner's
newer direction. No completed patched image, passing final scan, or patched-image
compatibility result is claimed. The offline available-fixes image `2bbb928d…`
from `23bd2fd` is built and imported under a hold: package scan 0 High/Critical
and dependency/rendering smoke pass, while DuckDB alone fails the strict gate.
It has no passing release attestation. The [dated security record](../docs/container_security_20261102.md)
binds its full image ID, export hashes, checks and unchanged service state.
Both services were stopped and masked at that time. The current receipts at
the top of this guide supersede that operational state; building or testing an
image alone does not establish scientific acceptance.

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
versions, runs the required Brotli binary/linkage gate, parses all application R
files, creates the SBOM and resolved dependency manifests, and fails on a High
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

A candidate staged beside a running release needs a release name of its own:
`releases/<name>/` is the running release's data directory, and staging it again
under the same name would rewrite the files under the live container. The name
is free text (the Plan A Level 1 candidate of 2026-10-25 used
`pubapps-2026-10-uniform-v10-perf1`); the earlier release stays in place for
rollback. When the image tree is unchanged, copy it on the VM from the earlier
release first (`cp -a .../data/images`, about 15 minutes for 145 GB on the share)
so the upload only has to match it by checksum. Copy the `deploy/` tooling of the
new commit to `tooling/deploy/` before activating with a flag the VM's copy does
not know yet, and record the commit in `tooling/DEPLOY_TOOLING_COMMIT`.

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

Since 2026-10-23 UFRC's proxy serves the three image aliases
(`/islet-explorer/local_images/`, `images/`, `img_data/`) itself from an Apache
`DocumentRoot` on the shared `/pubapps` filesystem instead of proxying them
(§ 4). `--docroot DIR` makes those aliases follow the deployment: once the new
service is verified, `DIR/islet-explorer/{local_images,images,img_data}` are
replaced atomically by symlinks to the release's `data/images`; if a later
step fails, they follow the restored release; an alias that is a real directory
is refused (exit 65) before any service is touched. The DocumentRoot is
recorded in the deployment state, so `rollback-release.sh` replays it. The
proxy's NFS lookup cache can hand out the previous target for up to a minute
after a swap.

`--local-data DIR` (2026-10-24) copies the release's `data/geometry`, `data/cells`
and `data/spatial` trees, about 880 MB whose many small reads are slow over NFS
(a drill-down's GeoPackage window took 0.2 s there against 17 ms on local disk),
to `DIR/<release>/data/` on the VM's local disk, checks every copied file
against `MANIFEST.sha256` (a mismatch stops the activation, exit 65, before any
service is touched) and mounts the copies over the same paths inside the
container; the image tree stays on `/pubapps`. The directory is recorded in the
deployment state and replayed by `rollback-release.sh`; copies of earlier
releases are kept so a restored unit still finds its data, and a retired
release's copy is removed by hand. On `pubmaigan1` the directory is
`/local/maigan/releases` (the 40 GB local disk that also holds the Podman store).

```bash
deploy/activate-release.sh \
  --project-root '/pubapps/<project>/islet-explorer' \
  --release 'pubapps-2026-01' \
  --publish-address 10.13.41.21 \
  --docroot '/pubapps/<project>/islet-explorer/prod/html' \
  --local-data /local/<service-account>/releases
```

### Staging instance

UFRC pairs every application hostname with a development one
(`devislet-explorer.rc.ufl.edu` beside `islet-explorer.rc.ufl.edu`, both on the
proxy). The project uses it for a second, independent service on the same VM:
a release candidate can be exercised through the real proxy and Basic Auth
before it replaces what reviewers are using. `--instance` names the service;
everything it owns is derived from that name, so the production service is
never touched:

```bash
deploy/activate-release.sh \
  --project-root '/pubapps/<project>/islet-explorer' \
  --release 'pubapps-2026-11' \
  --instance islet-explorer-dev --port 3839 --memory 4g \
  --publish-address 10.13.41.21 \
  --docroot '/pubapps/<project>/islet-explorer/dev/html'
```

| | production | staging |
|---|---|---|
| Quadlet unit / container | `islet-explorer.container`, `islet-explorer` | `islet-explorer-dev.container`, `islet-explorer-dev` |
| Host port (container 3838) | 3838 | 3839 |
| Memory cap | 9g | 4g (6g for the v20 release, 2026-11-02) |
| Runtime mounts | `runtime/logs`, `runtime/tmp` | `runtime/islet-explorer-dev/logs`, `.../tmp` |
| Deployment state | `deployments/active.state`, `previous.state` | `deployments/islet-explorer-dev.active.state`, `.previous.state` |
| Proxy hostname | `islet-explorer.rc.ufl.edu` → `pubmaigan1:3838` | `devislet-explorer.rc.ufl.edu` → `pubmaigan1:3839` |
| Proxy DocumentRoot and auth file | `prod/html`, `prod/auth` | `dev/html`, `dev/auth` |

Both instances mount the same read-only release directories, so staging the
active release costs no storage; a candidate release is a second directory
under `releases/`. Activation refuses a host port another unit already
publishes. The staging hostname carries the same Sensitive data and needs its
own Basic Auth entry at the proxy (maintainers only). Rollback takes the same
`--instance`.

Useful operations:

```bash
systemctl --user status islet-explorer.service islet-explorer-dev.service
journalctl --user -u islet-explorer.service -f
tail -F /pubapps/<project>/islet-explorer/runtime/logs/*.log

# Stop both instances (a problem on the site): the proxy then answers 502/503; the
# image aliases under prod/html and dev/html are static files and keep answering.
# A manual stop is never auto-restarted; `start` brings an instance back in about a
# minute (R start, release data, trajectory warm-up, then the health check).
systemctl --user stop islet-explorer.service islet-explorer-dev.service
systemctl --user start islet-explorer.service islet-explorer-dev.service
# The units are wanted by default.target and the account lingers, so a VM reboot
# starts them again; to keep one down across a reboot, mask it (unmask to undo).
systemctl --user mask islet-explorer-dev.service

deploy/rollback-release.sh \
  --project-root '/pubapps/<project>/islet-explorer'
deploy/rollback-release.sh \
  --project-root '/pubapps/<project>/islet-explorer' --instance islet-explorer-dev
```

Never edit an active release or point the service at a mutable image tag. Build
and verify a new release, then activate it. Keep at least the active and previous
image/data pairs until the new release has survived a VM reboot test.

## 4. UFRC reverse-proxy and Basic-Auth acceptance

UFRC configured the two vhosts on its proxy `vpubweb.rc.ufl.edu` on 2026-10-23
(ticket `83810`). The proxy is a separate host that mounts the same `/pubapps`
filesystem as the VM and runs each vhost as the service account (`mpm_itk`);
it cannot reach the VM loopback, so the backend is the VM's private address
(`pubmaigan1`, published with `--publish-address 10.13.41.21` since 2026-10-22
21:16 UTC; the publish-address rule above records the measured network
filtering that makes this acceptable). What the project operates, and where,
under the project root `/pubapps/maigan/islet-explorer/`:

| | production | staging |
|---|---|---|
| Public URL | `https://islet-explorer.rc.ufl.edu/islet-explorer/` | `https://devislet-explorer.rc.ufl.edu/islet-explorer/` |
| Backend (`ProxyPass`, WebSocket included) | `pubmaigan1:3838` | `pubmaigan1:3839` |
| Basic Auth (`<Location />`, the whole site) | `AuthUserFile` `prod/auth` | `dev/auth` |
| Apache `DocumentRoot` (landing page at `/`, `/favicon.ico`) | `prod/html` | `dev/html` |
| Image aliases served by Apache, not proxied (`ProxyPass ... !`) | `prod/html/islet-explorer/{local_images,images,img_data}` | the same under `dev/html` |
| Proxy access and error logs (`nobody:maigan`, readable) | `/pubapps/logs/maigan/islet-explorer-{access,error}.log` | `.../devislet-explorer-{access,error}.log` |

Because the proxy serves the three image aliases from the filesystem itself,
the Viewer's Range requests for multi-gigabyte OME-TIFFs never pass through
the container (measured 2026-10-23: exact 16-byte `206`, about 60 MB/s on
8 MiB mid-file reads from campus). `activate-release.sh --docroot` keeps each
alias a symlink into the verified release (§ 3); `shiny-server.conf` still
serves the same names inside the container for the VM-local verifier. The
auth files are bcrypt `htpasswd` files the project maintains itself (mode 640;
Apache reads them as the service account): hash locally with `htpasswd -B`,
never put a cleartext password in the ticket or on a command line, and replace
the file atomically (rsync or `mv`). The landing page is
`deploy/docroot/index.html` (its staging banner shows itself on the `dev`
hostname); install it with the favicon and header logos:

```bash
rsync -a deploy/docroot/index.html app/shiny_app/www/favicon.ico \
  app/shiny_app/www/logo-header-dark.png app/shiny_app/www/logo-header-light.png \
  '<vm>:/pubapps/<project>/islet-explorer/prod/html/'
```

`.htaccess` files in these DocumentRoots are not read by the proxy (measured
2026-10-23). Apache's own redirects (for example `/islet-explorer` without the
slash) name `http://` until UFRC sets `ServerName https://...` in the vhosts;
the load balancer bounces them back to HTTPS meanwhile. UFRC's first
configuration named `https://pubmaigan1:...` on the HTTP backend lines; the
`ws://` mapping listed first serves every request, so nothing fails, but the
correct scheme (`http://`) has been requested. Both corrections are in ticket
`83810`.

The contract the vhosts must keep: WebSocket upgrades, forwarded host, long
Shiny connections, unbuffered byte ranges for the TIFFs, and Basic Auth on the
page, WebSocket endpoint, downloads, TIFFs, XML, JSON and every other served
path. Parquet, GeoPackage, H5AD and other release-data paths must not be
web-accessible: only `data/images` is linked into the DocumentRoots.

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

Acceptance so far: on 2026-10-23, from the campus network and with UFRC's
temporary test account, `deploy/verify-reviewer-access.py` passed all 15
checks on both hostnames (anonymous `401` on the app, the Avivator entry
point, the three image aliases, both sidecars, the download route and the
WebSocket; authenticated `200`, the exact 16-byte `206`, `101` on Shiny
Server's SockJS route, which is the verifier's default WebSocket path since
that day). The run from a clean non-UF network with the reviewer credential is
still required.

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
  resolved Python package list, v2 Brotli verification JSON with the DuckDB
  binary regression and installed-library hash, `.scan.txt` attestation, and
  `published-image.txt`; retain the optional Grype database-status JSON when
  generated.
- Any future candidate's real Parquet compatibility report and input file hash, plus
  its command record identifying the immutable image, alongside the separate
  decoder regression result.
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

### Responsiveness baseline

`deploy/browser/bench-app.mjs` drives one fixed reviewer scenario against an
explicitly selected application URL and records, per step, how long the browser
waited and where the time went: the time until the server reported busy (one
round trip), the time until its last output message arrived (server compute plus
transfer), the bytes that crossed the WebSocket, and the browser's own draw. The
scenario is: page load until the Plot tab's two cards are drawn, two changes of the
plotted measure, a drill-down click, the Trajectory tab and its PAGA layout, the
Spatial tab and a donor change, and the Viewer's first tiles. Each cold run uses a
fresh browser context (empty cache); two warm visits follow. `--throttle home`
(40 ms, 25/5 Mbps) and `--throttle slow` (120 ms, 8/2 Mbps) add emulated link
conditions on top of the real path. Reports hold timings and byte counts only.

Run the same scenario from the vantage points that separate the app from the
hosting, and compare medians:

```bash
# the release image and data bundle on the build workstation (app alone)
podman run -d --name islet-bench ... localhost/islet-explorer:<release>   # same mounts and env as the unit
node deploy/browser/bench-app.mjs --url http://127.0.0.1:3843/islet-explorer/ \
  --output evidence/bench-workstation.json --runs 3 --label workstation-container
# the VM backend through the SSH tunnel (adds the VM and its NFS data)
node deploy/browser/bench-app.mjs --url http://127.0.0.1:3842/islet-explorer/ \
  --output evidence/bench-vm-direct.json --runs 3 --label vm-direct
# the public URL (adds UFRC's proxy, TLS and Basic Auth)
node deploy/browser/bench-app.mjs --url https://islet-explorer.rc.ufl.edu/islet-explorer/ \
  --auth-file ~/restricted/prod-auth.json --output evidence/bench-proxy.json --runs 3 --label proxy
node deploy/browser/bench-app.mjs --url https://islet-explorer.rc.ufl.edu/islet-explorer/ \
  --auth-file ~/restricted/prod-auth.json --output evidence/bench-proxy-home.json --runs 3 --throttle home
```

The 2026-10-24 baseline and its reading are in `docs/performance_baseline_20261024.md`:
the VM backend through the tunnel and the public URL from campus both matched the
workstation container step for step, so the hosting adds nothing measurable for one reviewer; what an off-campus reviewer feels
is the app's own payload (1.7 MB per Plot interaction, 3.6 MB for the first screen,
which is rendered twice) over their link, plus two one-time costs per R process
(the first Trajectory open, about 12 s, and the first PAGA layout, about 5 s).
