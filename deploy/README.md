# Islet Explorer on UFRC PubApps

**Split-release workflow:** use the [independent app/runtime/data procedure](#independent-app-runtime-and-data-releases) for development migration and later app releases. Dated evidence identifies the descriptor actually running.

**Recorded development status, 2026-11-05:** `islet-explorer-dev` completed split
migration, an app-only update using the same runtime/data, and both legacy and
split rollback checks. It now selects `split-dev-20261105-v1` (`app-v5`). Production
and both sites' access settings were unchanged. See the
[observed development acceptance](#observed-development-acceptance-2026-11-05)
and [dated receipt](../docs/evidence/split_release_20261105.json).

**Recorded live status, 2026-11-04:** production `pubapps-2026-11-ui-corrections-22bb011`
and the separate strict-DuckDB dev release `pubapps-2026-11-dev-duckdb-v1` were
active and healthy. Production used image
`sha256:fd6c9a990028892d43f45af53756a81b6f6ac13acaf7e10a75f61c1725dd483b`,
source `22bb011`, and a schema-5 20-donor bundle that passed 25/25 data-parity
checks. Local and host full hashes passed, the TIFF tree was reused by hard link,
and authenticated public-proxy access passed 15/15 checks. See the
[live receipt](../docs/evidence/pubapps_live_20261104_ui_corrections.json).
That production update did not change dev. The dated production DuckDB direction remains recorded in
the [remediation decision](../docs/duckdb_remediation_window_20261102.md).

Use the [workflow guide](../docs/workflow_user_guide.md) for source-to-release
steps and the [current v2 scientific build record](../docs/component_build_v2_20261102.md)
for the selected data decisions. Earlier masked-service, whole-v19, v1
diagnostic, and local-preview records are historical.

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
- New split releases use a software-only immutable OCI image plus a versioned, read-only app bundle. Existing legacy releases keep their baked-in app until migrated.
- App, runtime, and data each retain their own committed provenance and checksums. A deployment descriptor binds the exact compatible combination. Artifact builds export the selected committed revision; unrelated dirty research work is neither included nor changed. Legacy releases retain their original single-commit contract.
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
| `Containerfile` | Software-only production OCI runtime; app code and scientific data are separate read-only mounts |
| `containerignore` | Whitelist build context; excludes data, images, Git history, credentials, and caches |
| `requirements-python.lock` | Exact tested Python/PAGA environment |
| `install-production-r-deps.R` | Mandatory, exact `destiny` and `duckdb` gates over the app's R manifest |
| `verify-brotli.py` | Required actual decoder regression, package/library evidence, and system/Pillow linkage gate |
| `tests/check-duckdb-brotli.R` | Optional future real release Parquet BROTLI round-trip and query compatibility check with read-only input |
| `shiny-server.conf` | Sanitized `/islet-explorer/` Shiny service |
| `islet-explorer.container.in` | Hardened Quadlet template |
| `build-image.sh` / `scan-image.sh` | Committed runtime export, dependency probes, SBOM, and High/Critical scan gate |
| `export-image.sh` / `load-image.sh` | Checksummed OCI archive transfer through HPG and identity-checked local load |
| `generate-release-manifest.sh` / `verify-release.sh` | Immutable release metadata, hashes, and full-verification stamp |
| `activate-release.sh` / `rollback-release.sh` | Atomic activation with health rollback and explicit prior-release rollback |
| `verify-deployment.sh` | Local health, UID, mount, dependency, network, and optional Range checks |
| `run.sh` | Manual Podman smoke launcher for a legacy release, split descriptor, or app/data/runtime candidate |
| `release-artifacts.py` | Package app, register runtime/data, bind validation, assemble and verify split deployment descriptors |
| `runtime-inputs.py` | Fingerprint the runtime build inputs independently of app code |
| `prepare-local-data.py` | Verify a data-keyed local cache once and reuse it across app releases |

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

**Production decision, 2026-11-02; final UI release activated 2026-11-04:** production is within the notice's DuckDB patch window through **2026-11-10 UTC**, and the separate strict DuckDB fix remains in dev. See [the dated decision](../docs/duckdb_remediation_window_20261102.md). Production retains pinned CRAN DuckDB **1.5.2** (archive SHA-256 `55f98a336513a94354f328791c8b94521b9f2e918a57d7e9fe6c0cb6b9a38db6`); the historical backport at `3a12daa` is preserved separately. The isolated dev fix is deployed and strictly validated in the [dev receipt](../docs/evidence/dev_duckdb_deployment_20261103.json). The November 3 display and card-control changes are included in the production UI release recorded at the top of this runbook; the local preview is stopped.

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

Select a committed runtime revision with `--revision` (default `HEAD`); tag it
if the project's release practice calls for a tag. The current build script
exports only that revision, preserving unrelated dirty or untracked work. It
builds a software-only image; app code is packaged and tested separately in the
[split-release procedure](#independent-app-runtime-and-data-releases). Old
baked-app images remain valid only for their existing legacy descriptors.

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
versions, runs the required Brotli binary/linkage gate, verifies that no app
fallback is embedded, creates the SBOM and resolved dependency manifests, and fails on a High
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

The staging and hard-link safeguards in this section apply to both layouts.
The single-commit `generate-release-manifest.sh` examples describe **legacy
baked-app releases**. New split releases use the independent `seal-data` and
`register-data` procedure below; never combine a new software-only image with
these old activation assumptions.

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
For legacy baked-app releases, the builder and image share one clean Git
commit; legacy activation compares the OCI revision, data manifest, and release
metadata. Split releases retain separate app, runtime-build, and data-producer
commits. The data parity check still needs its intended committed source, so run
the builder and checks from an isolated clean checkout or worktree; do not clean
or stash unrelated work in the shared checkout.

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
rollback. When the image tree is byte-identical to a verified release on the
same Qumulo filesystem, hard-link it on the VM before transferring the other
files. First compare the complete image subtree hash/manifest entries. Hard
links share inodes, so treat both release trees as immutable and keep their
container mounts read-only. If the tree is not identical, transfer the
candidate images normally. New data always receives full content verification
after staging. Legacy activation repeats its full release hash pass; split
activation reuses the registered data receipt and checks its metadata/inventory
without a new TIFF content pass. A hard link does not change the viewer's path
or read behavior. Stage split tooling in a fresh `tooling/releases/COMMIT/deploy`
directory and invoke that exact version; preserve earlier tooling for rollback.
The historical legacy procedure used `tooling/deploy/` and recorded its commit
in `tooling/DEPLOY_TOOLING_COMMIT`.

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

# From the workstation, substituting the staged release root. Only reuse the
# previous image tree after confirming its entire hash matches the candidate.
PREVIOUS='pubapps-2026-11-components-v1'
ssh hpg ssh pubmaigan1 cp -al \
  "${PROJECT_ROOT}/releases/${PREVIOUS}/data/images" "${RELEASE_DIR}/data/images"
rsync -a --partial --info=progress2 --exclude='data/images/***' -e "${HOP}" \
  <generated-release-root>/ pubmaigan1:"${RELEASE_DIR}/"
rsync -a -e "${HOP}" /staging/islet-explorer-"${RELEASE}".oci.tar{,.metadata,.sha256} \
  pubmaigan1:"${PROJECT_ROOT}/incoming/"
```

For the legacy combined layout, copy the matching `deploy/` tooling to the VM
separately, generate release metadata with the exact image reference and shared
source commit, then perform full verification. That legacy workflow hashes TIFFs
during manifest generation, explicit verification, and activation. Ordinary
service restarts do not repeat the hash pass. For split releases, use `seal-data`
and `register-data` for new data; app-only activation does not repeat these
content passes or regenerate data metadata.

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

Legacy verification reads that declaration automatically, checks every manifest
hash, and writes `.verified` tied to the manifest hash. Legacy activation, the
manual launcher, and rollback use that verifier without an image-count override.
For a split descriptor those entrypoints instead resolve and verify its separate
app, runtime, and registered data artifacts.
Legacy `MANIFEST_VERSION=1` releases without `EXPECTED_IMAGES` retain their
15-image contract, including the frozen original release. For an optional
operator cross-check, `--expected-images N` must agree with the release metadata;
it cannot override it. Missing version-2 counts and duplicate, malformed, or
conflicting metadata fail verification. Even `--quick` checks the metadata hash
and actual image count, but does not rehash TIFF contents. If the bundle changes,
regenerate and fully verify it.

## 3. Activate the Quadlet service

The full-rehash description here is the legacy release behavior. Split
descriptors use the independent artifact checks and activation commands below.

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

For a legacy release, the command repeats the full content-hash gate before
replacement. For a split descriptor, it checks app contents, runtime evidence,
and registered data metadata/inventory. Both paths render
`~/.config/containers/systemd/islet-explorer.container`, reload user
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
step fails, their exact pre-attempt targets (including absent aliases) are
restored, whether the prior service was active or stopped. An alias that is a real
directory is refused (exit 65) before any service is touched. The DocumentRoot is
recorded in the deployment state, so `rollback-release.sh` replays it. The
proxy's NFS lookup cache can hand out the previous target for up to a minute
after a swap.

For legacy releases, `--local-data DIR` (2026-10-24) copies `data/geometry`, `data/cells`
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
Split descriptors use `DIR/datasets/DATA_ID/data/`: copies are fully verified
once, then reused across app releases after checking their data-manifest identity
and file-stat receipt. A changed published cache is rejected rather than updated
in place. Account for this initial cache's disk requirements before migration.

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

## Independent app, runtime, and data releases

This is the **development-only migration and app-release procedure**. It
describes commands and acceptance requirements; dated evidence records their
actual completion. Keep production and both sites' authentication, network
bindings, and proxy configuration unchanged. The earlier full-release
instructions remain the legacy path and rollback reference.

The ordinary app release packages all app code and assets, including cards,
plots, controls, downloads, and client-side behavior. There is no palette or UI
file allowlist. Installed software changes require a new runtime; data changes
require a new data release. Scientific calculations or loaders may change in an
app bundle without rebuilding the runtime or copying unchanged data, but require
scientific/compatibility validation before deployment.

```text
PROJECT_ROOT/
  apps/APP_ID/             # app/, scripts/, app.json, MANIFEST.sha256
  runtimes/RUNTIME_ID/     # runtime.json, image scan evidence, manifest
  datasets/DATA_ID/        # data.json verification receipt and manifest
  releases/DATA_SOURCE/    # existing immutable data/, release.env, manifest
  releases/DEPLOYMENT_ID/  # deployment.json, validation.json, manifest
  runtime/                # existing instance logs/tmp/deployment state
```

The app mounts at `/srv/islet-explorer/app`, including `app/docs/logo.png` and
its existing relative link. The app-coupled validator mounts at
`/srv/islet-explorer/scripts`. Data stays at `/srv/islet-explorer/data`, with
OME-TIFFs separately mounted at the existing nested `www/local_images` path.
All are read-only. The three static aliases and byte-range delivery are preserved.
The new image has no embedded app fallback.

Each descriptor binds the app manifest and commit, immutable runtime and build
commit, verified data receipt and producer commit, and combination-specific
validation. Runtime-input fingerprints prevent selecting a runtime built from
incompatible package/install/server inputs. A matching fingerprint does not
replace app or scientific tests. App bundles also record the normalized,
digest-pinned base image from their committed `deploy/base-image.lock`. Runtime
registration binds the published base-image evidence, and launch verifies it
against the actual image's `org.opencontainers.image.base.name` label. Changing
only the base lock therefore refuses reuse of an incompatible runtime, even when
its software-input fingerprint is unchanged. `build-image.sh --base-image` must
match the selected source revision's lock.
`change-kind` is an explicit operator statement:
`app` requires app tests and smoke checks; `science`, `data`, and `runtime` also
require passing scientific validation. `app` and `science` assembly require
`--base-release` pointing to a prior split descriptor and must reuse its exact
runtime manifest and data receipt. A runtime/data replacement cannot be labeled
app-only. First migration uses `runtime` with scientific acceptance.
Classification is based on effect, not
which filenames changed. Validation JSON records completed checks; it is not a
command that runs them or independently proves an operator's assertion.

### Observed development acceptance, 2026-11-05

The development PubApp completed initial split activation, rollback to the legacy
release, reactivation, a subsequent app-only update, and rollback to the previous
split descriptor. The final selected deployment is
`split-dev-20261105-v1`, app `app-v5` at
`b4867688c13a772a76527d8e547f1b0e11bdf4af`. The temporary component-change
acceptance app was `app-v6` at `2f123ec3e2cacf5d5e9c8d53499fa0e2cb92c783`;
its head metadata component was observed after update and absent after rollback.
Both descriptors selected runtime artifact `split-runtime-20261105-v2`, image
`sha256:e14329ba52d81dc85f65e6fd13023a8bb0ac47cbfd9ee37c04ad8681ea4068c2`,
and data artifact `split-dev-data-20261105`, produced at
`378843e06c0cd66d8dfff5e3b60664f8a45569ee`, with data manifest
`e4b2dbd3bf709513bc6a1d5f8cfdb76ac856efacec25940e13c9975f03936b41`.
Activation and rollback used release tooling
`9f682a71943b56f242f9ae372b391d0843ea0fa4`.
The [dated acceptance receipt](../docs/evidence/split_release_20261105.json)
records complete artifact identities, evidence hashes, commands, and timings.

| Development operation | Observed elapsed time |
| --- | ---: |
| First split activation | 63.05 s |
| Rollback to the legacy combined release | 252.49 s |
| Reactivation of the split release | 90.00 s |
| App-only activation | 52.70 s |
| Rollback to the previous split descriptor | 58.10 s |

App-only descriptor assembly took 2.134 s; the full bounded verification and
staging step took 5.069 s. Those timings exclude the separately recorded app tests and
browser acceptance. The app-only update reused the exact runtime image and
verified data, with no image build and no full TIFF-content hashing. Legacy
rollback retains full legacy verification, which explains its different cost.
These are observed timings for this change and VM, not an upper bound for every
app change or cold start.

Validation included 72 Python deployment/artifact/cache tests, 1,142 R unit
expectations for each candidate bundle, and the complete 32-check scientific
report. An initial Python 3.9 compatibility preflight failure occurred before
service changes; corrected tooling passed nine compatibility checks on the VM.
The full scientific run kept all data/image hashes and verified both the
historical producer's raw summaries and the candidate's five served summaries.
Each full browser acceptance stage passed 49 checks; authenticated public access
passed 15 checks per stage. Backend and development-proxy image requests returned
HTTP 206 with the expected range bytes, and unauthenticated public requests were
rejected. The first final-rollback browser attempt retained a highlight after
unhover and failed that check; its evidence is retained. An unchanged fresh-session
rerun passed all 49 checks. The harness was then changed to wait up to six seconds
for the actual rendered null highlight after each hover direction, preserving the
strict assertion. A further confirmation passed all 49 browser and 15 access
checks, with the update marker absent. This records an intermittent check failure;
it does not establish a proven application defect.

Production service/artifact state and both sites' access settings remained
unchanged. Host-level protection of verified data against mutation, off-campus
UFRC reachability, and service recovery after a VM reboot still need their
separate operational confirmation. A metadata receipt and read-only container
mounts do not themselves make host-owned files immutable.

### First migration: build and stage

On the workstation, commit only the intended deployment/app changes, leaving
unrelated local work untouched. Choose the full committed revision containing the
split tooling. Both commands export committed files in isolation; uncommitted
changes are never silently released. Use fresh IDs and directories; artifacts
cannot be overwritten.

```bash
REPO=/data/panc_CODEX/Islet-Explorer
BUILD_PROJECT="$REPO/tmp/split-release-build"
REVISION='<full-committed-split-revision>'
RUNTIME=runtime-split-v1
APP=app-split-v1
EVIDENCE="$BUILD_PROJECT/evidence/$RUNTIME"
mkdir -p "$EVIDENCE" "$BUILD_PROJECT/incoming"
BASE_IMAGE=$(sed '/^#/d; /^$/d' "$REPO/deploy/base-image.lock")

"$REPO/deploy/build-image.sh" --revision "$REVISION" \
  --tag "localhost/islet-explorer-runtime:$RUNTIME" --version "$RUNTIME" \
  --base-image "$BASE_IMAGE" --scan-dir "$EVIDENCE/scan"
IMAGE_REF=$(sed -n 's/^LOCAL_IMAGE_ID=//p' "$EVIDENCE/scan/published-image.txt")
RUNTIME_INPUTS=$(podman image inspect --format \
  '{{ index .Labels "io.islet-explorer.runtime-inputs-sha256" }}' "$IMAGE_REF")
python3 "$REPO/deploy/release-artifacts.py" register-runtime \
  --project-root "$BUILD_PROJECT" --runtime "$RUNTIME" \
  --scan-dir "$EVIDENCE/scan" --runtime-input-sha256 "$RUNTIME_INPUTS"
python3 "$REPO/deploy/release-artifacts.py" package-app \
  --project-root "$BUILD_PROJECT" --app "$APP" --revision "$REVISION" --repo "$REPO"
"$REPO/deploy/export-image.sh" "$IMAGE_REF" "$BUILD_PROJECT/incoming/$RUNTIME.oci.tar"
```

Before transferring, inspect actual VM capacity. A read-only check on November 5
found only 3.6 GB free on `/local` (the Podman store) and 29 GB on `/pubapps`.
Measure the new image's additional layer requirements; do not delete active or
rollback images/data to force it to fit. Resolve a capacity shortage before load.
Stage the committed tooling into a new versioned directory, not over the tooling
used by an existing deployment. The commands below use the established HPG hop:

```bash
VM_PROJECT=/pubapps/maigan/islet-explorer
HOP='ssh hpg ssh'
ssh hpg "ssh pubmaigan1 'df -h /local/maigan /pubapps/maigan; podman system df'"
ssh hpg "ssh pubmaigan1 'mkdir -p $VM_PROJECT/apps $VM_PROJECT/runtimes $VM_PROJECT/incoming $VM_PROJECT/tooling/releases/$REVISION'"
git -C "$REPO" archive "$REVISION" deploy | \
  ssh hpg "ssh pubmaigan1 'tar -x -C $VM_PROJECT/tooling/releases/$REVISION'"
rsync -a -e "$HOP" "$BUILD_PROJECT/apps/$APP" pubmaigan1:"$VM_PROJECT/apps/"
rsync -a -e "$HOP" "$BUILD_PROJECT/runtimes/$RUNTIME" pubmaigan1:"$VM_PROJECT/runtimes/"
rsync -a --partial -e "$HOP" "$BUILD_PROJECT/incoming/$RUNTIME.oci.tar"{,.metadata,.sha256} \
  pubmaigan1:"$VM_PROJECT/incoming/"
```

On the VM as the existing service account, set the same `REVISION`, `RUNTIME`,
and `APP`. Confirm the current dev state before selecting its data source. The
November 5 read-only inspection found the following settings; these values are
not permission to overwrite a newer deployment without checking it:

```bash
PROJECT_ROOT=/pubapps/maigan/islet-explorer
TOOLING="$PROJECT_ROOT/tooling/releases/$REVISION/deploy"
DATA=dev-data-split-v1
DATA_SOURCE=pubapps-2026-11-dev-duckdb-v1
RELEASE=dev-split-v1
EVIDENCE="$PROJECT_ROOT/incoming/$RELEASE-evidence"
mkdir "$EVIDENCE"
systemctl --user status islet-explorer-dev.service --no-pager
cat "$PROJECT_ROOT/runtime/deployments/islet-explorer-dev.active.state"
cp -p "$HOME/.config/containers/systemd/islet-explorer-dev.container" "$EVIDENCE/prior-dev.container"
cp -p "$PROJECT_ROOT/runtime/deployments/islet-explorer-dev.active.state" "$EVIDENCE/prior-dev.state"
"$TOOLING/load-image.sh" "$PROJECT_ROOT/incoming/$RUNTIME.oci.tar"
python3 "$TOOLING/release-artifacts.py" register-data \
  --project-root "$PROJECT_ROOT" --data "$DATA" --source-release "$DATA_SOURCE"
```

`register-data` performs the initial complete content check, including TIFFs,
and writes an independent data receipt. It does not rewrite the source release,
copy its images, or inherit the old image's security decision. The selected
runtime's own evidence is checked separately. Subsequent app deployments verify
the receipt, source manifest, metadata, image count, and file inventory instead
of rereading TIFF contents. Published data and hard-linked image trees must stay
immutable: metadata/stat checks are not protection against an owner deliberately
rewriting both data and receipts. Changed data requires a fresh ID and full
registration. Preserve the existing hard-link staging procedure for new data
releases, including its complete image-tree equivalence check.

### Candidate checks and descriptor assembly

Run these checks against the exact candidate artifacts before assembly. The
manual launcher also works on a workstation project containing the registered
artifacts and data; use an unused loopback port there. On the 15 GiB VM, stop
only dev while testing a 4 GiB candidate: production plus the existing dev and a
third 4 GiB container would exceed the planned memory caps. Keep the original dev
unit/state available and restore that service if any candidate check fails.

On the VM, in a dedicated shell:

```bash
set -euo pipefail
export PROJECT_ROOT APP DATA RUNTIME RELEASE EVIDENCE TOOLING
systemctl --user is-active --quiet islet-explorer-dev.service
CANDIDATE=islet-explorer-split-smoke
cleanup_candidate() {
  podman stop "$CANDIDATE" >/dev/null 2>&1 || true
  systemctl --user start islet-explorer-dev.service
}
trap cleanup_candidate EXIT
systemctl --user stop islet-explorer-dev.service
PROJECT_ROOT="$PROJECT_ROOT" APP_BUNDLE="$APP" DATA_RELEASE="$DATA" RUNTIME="$RUNTIME" \
  CONTAINER_NAME="$CANDIDATE" PORT=3839 MEMORY=4g "$TOOLING/run.sh" \
  > "$EVIDENCE/container-smoke.log" 2>&1
podman exec --workdir /srv/islet-explorer/app/shiny_app "$CANDIDATE" \
  Rscript tests/parse_check.R > "$EVIDENCE/parse.log" 2>&1
python3 "$TOOLING/release-artifacts.py" resolve-candidate --project-root "$PROJECT_ROOT" \
  --app "$APP" --data "$DATA" --runtime "$RUNTIME" > "$EVIDENCE/candidate.env"
```

Use `set -e` in a scripted run so a failed test exits through the cleanup trap.
Inspect warnings/skips as well as exit codes. The runtime deliberately excludes
`testthat`; run unit tests on the workstation's test environment against the
packaged app, with runtime support sources exported from the selected runtime
commit. These tests complement actual container parse/startup and browser checks;
they are not evidence that the production image contains the test environment.

```bash
# Workstation: use the same fresh app bundle transferred to the VM.
mkdir -p "$LOCAL_EVIDENCE/runtime-source"
RUNTIME_REVISION=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["git_commit"])' \
  "$BUILD_PROJECT/runtimes/$RUNTIME/runtime.json")
git -C "$REPO" archive "$RUNTIME_REVISION" deploy | \
  tar -x -C "$LOCAL_EVIDENCE/runtime-source"
(
  cd "$BUILD_PROJECT/apps/$APP/app/shiny_app"
  ISLET_TEST_RUNTIME_SOURCE="$LOCAL_EVIDENCE/runtime-source" \
    Rscript tests/run_tests.R
) > "$LOCAL_EVIDENCE/app-tests.log" 2>&1
rsync -a -e 'ssh hpg ssh' "$LOCAL_EVIDENCE/app-tests.log" \
  pubmaigan1:"$EVIDENCE/app-tests.log"
```

Use an absolute `LOCAL_EVIDENCE` path. From the workstation, open an SSH tunnel
using the configured HPG jump route, then run browser checks against it:

```bash
ssh -N -J hpg -L 3845:127.0.0.1:3839 maigan@pubmaigan1.ufhpc
# In another workstation shell, with a fresh local evidence directory:
node deploy/browser/check-app.mjs --url http://127.0.0.1:3845/islet-explorer/ \
  --output "$LOCAL_EVIDENCE/browser.json" --chrome /usr/bin/google-chrome --expanded true
rsync -a -e 'ssh hpg ssh' "$LOCAL_EVIDENCE/browser.json" \
  pubmaigan1:"$EVIDENCE/browser.json"
```

Run the targeted checks for affected components as well. For first migration,
scientific/loader changes, data changes, or runtime updates, run and retain the
applicable full scientific/compatibility report against this exact combination.
For full parity against an older data producer, run from a clean, isolated
checkout of the candidate app commit on the scientific validation host. Supply
both explicit commits; neither flag alone is accepted:

```bash
# In the clean, isolated candidate app checkout:
test "$(git rev-parse HEAD)" = "$APP_COMMIT"
test -z "$(git status --porcelain --untracked-files=all)"
git cat-file -e "${DATA_COMMIT}^{commit}"
Rscript scripts/verify_release_parity.R \
  --source-root "$CANONICAL_SOURCE" --image-root "$CANONICAL_IMAGES" \
  --release-root "$REGISTERED_DATA_DIR" --output-prefix "$PARITY_EVIDENCE/release-parity" \
  --app-commit "$APP_COMMIT" --data-commit "$DATA_COMMIT"
```

`APP_COMMIT` must be that clean checkout's full commit; `DATA_COMMIT` must be
the full clean producer recorded in the selected data manifest. The report keeps
app, validator, and data-producer provenance separately. Raw stored summaries
are checked against preparation by the exact producer's committed helpers,
exported into a temporary directory. Separately, the candidate's actual shared
loader is exercised on canonical source and production data, and all five served
summary tables must agree. This tests both stored scientific values and current
app behavior without hand-normalizing differences. All full hashes, source,
spatial, and geometry comparisons still run. Missing producer source/factory
fails closed. The exact `DATA_COMMIT` must exist in the validation checkout's
local Git object store: a shallow clone that omits it cannot perform this check.
Obtain the reviewed producer commit before validation rather than substituting
another revision. The temporary export contains committed app/scripts inputs,
not scientific data, and is removed on success or failure. Without the flags,
the historical same-commit/verifier-only rule is retained.

During November 5 migration preparation, the first independent-commit run
(`parity-r2`) retained a failed result for two summary checks because it compared
current preparation with older raw stored tables. The focused diagnostic showed
that the producer's raw tables and the candidate's served outputs each agree at
their own boundary. The corrected verifier now checks both boundaries, including
all five raw and all five served summaries. The subsequent complete `parity-r3`
run passed all 32 checks, including the full hashes, for the November 5 migration.
The failed report and diagnostic remain part of the dated evidence; the final
acceptance uses the completed corrected run.

Set `CHANGE_KIND=runtime` for the first migration, and set `SCIENTIFIC_REPORT` to
that passing report. For an ordinary non-scientific app change use
`CHANGE_KIND=app`; no new scientific attestation is needed. The script below
requires an explicit `SCIENTIFIC_RESULT=pass` for the other change kinds. Set it
only after inspecting the completed scientific report, never to bypass a failed
or missing check.

Back in the VM candidate shell, after the checks have actually passed:

```bash
export CHANGE_KIND SCIENTIFIC_RESULT SCIENTIFIC_REPORT
python3 - "$EVIDENCE" <<'PY'
import hashlib, json, os, pathlib, sys
root = pathlib.Path(sys.argv[1])
values = dict(line.split('=', 1) for line in (root/'candidate.env').read_text().splitlines())
kind = os.environ['CHANGE_KIND']
assert json.loads((root/'browser.json').read_text()).get('passed') is True, 'Browser acceptance failed'
checks = {'app_tests': 'pass', 'smoke': 'pass'}
reports = [root/'container-smoke.log', root/'parse.log', root/'app-tests.log', root/'browser.json']
if kind != 'app':
    assert os.environ.get('SCIENTIFIC_RESULT') == 'pass', 'Scientific acceptance is required'
    reports.append(pathlib.Path(os.environ['SCIENTIFIC_REPORT']))
    checks['scientific'] = 'pass'
validation = {
    'app_manifest_sha256': values['APP_MANIFEST_SHA256'],
    'image_ref': values['IMAGE_REF'],
    'data_receipt_sha256': values['DATA_RECEIPT_SHA256'],
    'change_kind': kind, 'checks': checks,
    'evidence': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in reports},
}
(root/'validation.json').write_text(json.dumps(validation, indent=2)+'\n')
PY
base_args=()
if [[ "$CHANGE_KIND" == app || "$CHANGE_KIND" == science ]]; then
  base_args=(--base-release "${BASE_RELEASE:?set the previous split deployment ID}")
fi
python3 "$TOOLING/release-artifacts.py" assemble --project-root "$PROJECT_ROOT" \
  --release "$RELEASE" --app "$APP" --data "$DATA" --runtime "$RUNTIME" \
  --change-kind "$CHANGE_KIND" --validation "$EVIDENCE/validation.json" "${base_args[@]}"
python3 "$TOOLING/release-artifacts.py" verify --release-dir "$PROJECT_ROOT/releases/$RELEASE"
cleanup_candidate
trap - EXIT
```

The old dev service is running again before activation so normal failure recovery
restores its previously active state. Activate only dev, using the checked VM
settings, then exercise its authenticated public endpoint and viewer Range path:

```bash
"$TOOLING/activate-release.sh" --project-root "$PROJECT_ROOT" --release "$RELEASE" \
  --instance islet-explorer-dev --port 3839 --memory 4g --publish-address 10.13.41.21 \
  --docroot "$PROJECT_ROOT/dev/html" --local-data /local/maigan/releases
# Read verified values as data; never source/eval an artifact metadata file.
IMAGE_REF=$(sed -n 's/^IMAGE_REF=//p' "$EVIDENCE/candidate.env")
APP_SHA=$(sed -n 's/^APP_MANIFEST_SHA256=//p' "$EVIDENCE/candidate.env")
DATA_DIR=$(sed -n 's/^DATA_DIR=//p' "$EVIDENCE/candidate.env")
"$TOOLING/verify-deployment.sh" --url http://10.13.41.21:3839/islet-explorer/ \
  --container islet-explorer-dev --publish-address 10.13.41.21 \
  --expected-image "$IMAGE_REF" --expected-app "$PROJECT_ROOT/apps/$APP/app" \
  --expected-app-sha256 "$APP_SHA" --expected-data "$DATA_DIR" \
  --expected-local-data "/local/maigan/releases/datasets/$DATA/data" \
  --range-path 'local_images/<selected-image>.ome.tiff'
```

Use the existing external reviewer-access procedure for
`https://devislet-explorer.rc.ufl.edu/islet-explorer/`; do not change credentials
or proxy rules. Save actual container inspection, app/data mount identities,
active state, backend/public Range results, and elapsed release timings. A new
local cache under `/local/maigan/releases/datasets/DATA_ID` is copied and checked
once; later app deployments reuse it. Missing or mutated published caches fail
instead of being repaired under a running service.

### Next app-only release and rollback

Select another committed app revision and fresh `APP`/`RELEASE` IDs. Keep `DATA`
and `RUNTIME` from the active descriptor. Package and transfer just the new app:

```bash
python3 "$REPO/deploy/release-artifacts.py" package-app --project-root "$BUILD_PROJECT" \
  --app "$APP" --revision "$REVISION" --repo "$REPO"
rsync -a -e 'ssh hpg ssh' "$BUILD_PROJECT/apps/$APP" pubmaigan1:"$VM_PROJECT/apps/"
```

Set `BASE_RELEASE` to the current split deployment ID, then repeat candidate
checks and assembly above with `CHANGE_KIND=app`, or `science` when scientific
behavior changed. The `--base-release` check pins both reused artifact identities;
then activate dev. Do not rerun runtime build,
image export/load, data registration, or data transfer for unchanged artifacts.
If runtime inputs differ, resolution refuses the combination; follow the runtime
build/scan path and validate the new combination. Never edit an existing app
bundle, data receipt, descriptor, or image tag to change a running release.

Explicit rollback selects the previously successful complete descriptor and its
recorded instance settings. It supports old baked-image releases as well as
split releases; legacy rollback retains its old full verification behavior:

```bash
"$TOOLING/rollback-release.sh" --project-root "$PROJECT_ROOT" --instance islet-explorer-dev
```

Demonstrate both same-runtime app rollback and rollback to the legacy dev release
before considering migration complete. Keep referenced app bundles, runtime
images, original data sources/receipts, descriptors, and local caches for active
and rollback releases. Do not remove shared data merely because one app release
is retired.

### New data and runtime releases

For new data, build into a fresh `releases/DATA_SOURCE/data` directory using the
existing scientific workflow, stage it with the preserved TIFF hard-link rules,
and retain full scientific acceptance. Seal its data-only metadata and register
it on the target host:

```bash
python3 "$TOOLING/release-artifacts.py" seal-data --project-root "$PROJECT_ROOT" \
  --source-release "$DATA_SOURCE" --data-schema 5 --expected-images 20
python3 "$TOOLING/release-artifacts.py" register-data --project-root "$PROJECT_ROOT" \
  --data "$DATA" --source-release "$DATA_SOURCE"
```

Use the approved schema/count for that data release, not an inferred file count.
The producer commit comes from its data manifest. Select the new data ID with the
compatible app/runtime and complete `CHANGE_KIND=data` validation. Runtime
changes similarly use a new runtime ID and build/scan evidence, followed by
`CHANGE_KIND=runtime` validation; they do not require copying unchanged data.
