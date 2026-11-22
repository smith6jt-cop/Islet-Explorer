#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: deploy/build-image.sh --tag IMAGE:TAG --version RELEASE \
         --base-image IMAGE@sha256:DIGEST --scan-dir DIR [--push] [--revision COMMIT]

Builds the software-only OCI runtime from an isolated committed Git export
(default HEAD). Local uncommitted work is never included or altered. A
digest-pinned base matching that revision's deploy/base-image.lock is required.
With --push, the immutable registry reference is printed at the end.
Syft and Grype are required; the SBOM and High/Critical vulnerability gate are
written to --scan-dir.
EOF
}

tag=""
version=""
base_image=""
push_image=0
scan_dir=""
requested_revision="HEAD"

while (($#)); do
  case "$1" in
    --tag) tag="${2:?missing value for --tag}"; shift 2 ;;
    --version) version="${2:?missing value for --version}"; shift 2 ;;
    --base-image) base_image="${2:?missing value for --base-image}"; shift 2 ;;
    --push) push_image=1; shift ;;
    --scan-dir) scan_dir="${2:?missing value for --scan-dir}"; shift 2 ;;
    --revision) requested_revision="${2:?missing value for --revision}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 64 ;;
  esac
done

[[ -n "${tag}" && -n "${version}" && -n "${base_image}" && -n "${scan_dir}" ]] || {
  usage >&2
  exit 64
}
[[ "${base_image}" =~ ^[A-Za-z0-9._:/-]+@sha256:[0-9a-f]{64}$ ]] || {
  echo "ERROR: --base-image must include an immutable @sha256: digest." >&2
  exit 64
}
[[ "${tag}" =~ ^[A-Za-z0-9._:/-]+:[A-Za-z0-9._-]+$ && "${tag}" != *@sha256:* ]] || {
  echo "ERROR: --tag must be a mutable build tag without whitespace or a digest." >&2
  exit 64
}
[[ "${version}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
  echo "ERROR: --version may contain only letters, digits, dot, underscore, and hyphen." >&2
  exit 64
}

for command_name in git podman syft grype python3 tar; do
  command -v "${command_name}" >/dev/null || {
    echo "ERROR: required command not found: ${command_name}" >&2
    exit 69
  }
done

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd -- "${script_dir}/.." && pwd -P)"

revision="$(git -C "${repo_root}" rev-parse --verify "${requested_revision}^{commit}")"
created="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
iid_file="$(mktemp)"
digest_file=""
build_context="$(mktemp -d)"
cleanup_temp_files() {
  rm -f -- "${iid_file}"
  rm -rf -- "${build_context}"
  if [[ -n "${digest_file}" ]]; then
    rm -f -- "${digest_file}"
  fi
}
trap cleanup_temp_files EXIT

# Both the build inputs and its scanner come from the selected revision. This
# also allows releases while unrelated research changes remain in the checkout.
git -C "${repo_root}" archive "${revision}" deploy scripts/install_shiny_deps.R | \
  tar -x -C "${build_context}"
build_deploy="${build_context}/deploy"
[[ -f "${build_deploy}/runtime-inputs.py" ]] &&
  grep -q 'io.islet-explorer.runtime-layout="split-v1"' "${build_deploy}/Containerfile" || {
    echo "ERROR: selected revision predates the split runtime; commit the runtime changes before building." >&2
    exit 65
  }
runtime_inputs="$(python3 "${build_deploy}/runtime-inputs.py" "${build_context}")"
# Apply today's selection guard to the exported lock, including when rebuilding
# an older split revision whose helper predates the --base-image CLI option.
required_base="$(python3 "${script_dir}/runtime-inputs.py" --base-image "${build_context}")"
[[ "${base_image}" == "${required_base}" ]] || {
  echo "ERROR: --base-image must match selected revision's deploy/base-image.lock (${required_base})." >&2
  exit 65
}
echo "Building runtime from committed Git ${revision}; runtime inputs ${runtime_inputs}."

# Check the digest-selected base before installing the full scientific stack.
# Rocker builds R under /usr/local; resolve it through the image's PATH.
podman run --rm --pull=missing --entrypoint /bin/sh "${base_image}" -ceu '
  test "$(uname -m)" = x86_64
  command -v Rscript
  command -v shiny-server
  id shiny
  Rscript -e '\''stopifnot(getRversion() == "4.5.1"); cat("Base R repository: ", getOption("repos"), "\n")'\''
'

podman build \
  --pull=missing \
  --iidfile "${iid_file}" \
  --ignorefile "${build_deploy}/containerignore" \
  --file "${build_deploy}/Containerfile" \
  --tag "${tag}" \
  --build-arg "BASE_IMAGE=${base_image}" \
  --build-arg "OCI_CREATED=${created}" \
  --build-arg "OCI_REVISION=${revision}" \
  --build-arg "OCI_VERSION=${version}" \
  --build-arg "RUNTIME_INPUTS_SHA256=${runtime_inputs}" \
  "${build_context}"

image_id="$(tr -d '[:space:]' < "${iid_file}")"
[[ "${image_id}" =~ ^sha256:[0-9a-f]{64}$ ]] || {
  echo "ERROR: Podman did not report a valid immutable build image ID." >&2
  exit 70
}

image_revision="$(podman image inspect --format '{{ index .Labels "org.opencontainers.image.revision" }}' "${image_id}")"
image_version="$(podman image inspect --format '{{ index .Labels "org.opencontainers.image.version" }}' "${image_id}")"
image_inputs="$(podman image inspect --format '{{ index .Labels "io.islet-explorer.runtime-inputs-sha256" }}' "${image_id}")"
image_layout="$(podman image inspect --format '{{ index .Labels "io.islet-explorer.runtime-layout" }}' "${image_id}")"
[[ "${image_revision}" == "${revision}" && "${image_version}" == "${version}" &&
   "${image_inputs}" == "${runtime_inputs}" && "${image_layout}" == split-v1 ]] || {
  echo "ERROR: built image labels do not match the requested Git revision and release version." >&2
  exit 70
}

podman run --rm --entrypoint /opt/islet-explorer-py/bin/python "${image_id}" -c \
  'import anndata, scanpy, igraph, leidenalg, fa2_modified; print("python dependency probe: ok")'
podman run --rm --entrypoint Rscript "${image_id}" -e \
  'stopifnot(packageVersion("destiny")=="3.22.0", packageVersion("duckdb")=="1.5.2"); cat("R dependency probe: ok\n")'
# The runtime must not carry a fallback app or a stale app-specific validator.
podman run --rm --entrypoint /bin/sh "${image_id}" -ceu \
  'test ! -e /srv/islet-explorer/app/shiny_app/app.R; test ! -e /srv/islet-explorer/scripts/verify_release_data.R'

"${build_deploy}/scan-image.sh" "${image_id}" "${scan_dir}"

echo "Built ${tag} as ${image_id} from Git ${revision}."

if ((push_image)); then
  digest_file="$(mktemp)"
  # Push the already scanned image ID, not the mutable tag a second time.
  # Supplying the tag as the destination preserves the requested registry name
  # without reopening a tag-swap window between scanning and publication.
  podman push --digestfile "${digest_file}" "${image_id}" "${tag}"
  digest="$(tr -d '[:space:]' < "${digest_file}")"
  repository="${tag%:*}"
  [[ "${repository}" != "${tag}" && "${digest}" =~ ^sha256:[0-9a-f]{64}$ ]] || {
    echo "ERROR: could not derive immutable registry reference after push." >&2
    exit 70
  }
  echo "Immutable image: ${repository}@${digest}"
  {
    printf 'REGISTRY_REF=%s@%s\n' "${repository}" "${digest}"
    printf 'IMAGE_ID=%s\n' "${image_id}"
    printf 'GIT_COMMIT=%s\n' "${revision}"
    printf 'BASE_IMAGE=%s\n' "${base_image}"
    printf 'RUNTIME_INPUTS_SHA256=%s\n' "${runtime_inputs}"
    printf 'RUNTIME_LAYOUT=split-v1\n'
  } > "${scan_dir}/published-image.txt"
else
  echo "Immutable local image: ${image_id}"
  {
    printf 'LOCAL_IMAGE_ID=%s\n' "${image_id}"
    printf 'GIT_COMMIT=%s\n' "${revision}"
    printf 'BASE_IMAGE=%s\n' "${base_image}"
    printf 'RUNTIME_INPUTS_SHA256=%s\n' "${runtime_inputs}"
    printf 'RUNTIME_LAYOUT=split-v1\n'
  } > "${scan_dir}/published-image.txt"
fi
