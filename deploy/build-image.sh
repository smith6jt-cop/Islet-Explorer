#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: deploy/build-image.sh --tag IMAGE:TAG --version RELEASE \
         --base-image IMAGE@sha256:DIGEST --scan-dir DIR [--push]

Builds the data-free OCI image from a clean Git release. A digest-pinned base is
required. With --push, the immutable registry reference is printed at the end.
Syft and Grype are required; the SBOM and High/Critical vulnerability gate are
written to --scan-dir.
EOF
}

tag=""
version=""
base_image=""
push_image=0
scan_dir=""

while (($#)); do
  case "$1" in
    --tag) tag="${2:?missing value for --tag}"; shift 2 ;;
    --version) version="${2:?missing value for --version}"; shift 2 ;;
    --base-image) base_image="${2:?missing value for --base-image}"; shift 2 ;;
    --push) push_image=1; shift ;;
    --scan-dir) scan_dir="${2:?missing value for --scan-dir}"; shift 2 ;;
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

for command_name in git podman syft grype; do
  command -v "${command_name}" >/dev/null || {
    echo "ERROR: required command not found: ${command_name}" >&2
    exit 69
  }
done

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd -- "${script_dir}/.." && pwd -P)"

if [[ -n "$(git -C "${repo_root}" status --porcelain --untracked-files=all)" ]]; then
  echo "ERROR: image builds require a clean, committed worktree." >&2
  git -C "${repo_root}" status --short >&2
  exit 65
fi

revision="$(git -C "${repo_root}" rev-parse HEAD)"
created="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
iid_file="$(mktemp)"
digest_file=""
cleanup_temp_files() {
  rm -f -- "${iid_file}"
  if [[ -n "${digest_file}" ]]; then
    rm -f -- "${digest_file}"
  fi
}
trap cleanup_temp_files EXIT

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
  --ignorefile "${script_dir}/containerignore" \
  --file "${script_dir}/Containerfile" \
  --tag "${tag}" \
  --build-arg "BASE_IMAGE=${base_image}" \
  --build-arg "OCI_CREATED=${created}" \
  --build-arg "OCI_REVISION=${revision}" \
  --build-arg "OCI_VERSION=${version}" \
  "${repo_root}"

image_id="$(tr -d '[:space:]' < "${iid_file}")"
[[ "${image_id}" =~ ^sha256:[0-9a-f]{64}$ ]] || {
  echo "ERROR: Podman did not report a valid immutable build image ID." >&2
  exit 70
}

podman run --rm --entrypoint /opt/islet-explorer-py/bin/python "${image_id}" -c \
  'import anndata, scanpy, igraph, leidenalg, fa2_modified; print("python dependency probe: ok")'
podman run --rm --entrypoint Rscript "${image_id}" -e \
  'stopifnot(packageVersion("destiny")=="3.22.0", packageVersion("duckdb")=="1.5.2"); cat("R dependency probe: ok\n")'
podman run --rm --entrypoint Rscript \
  --workdir /srv/islet-explorer/app/shiny_app "${image_id}" tests/parse_check.R

"${script_dir}/scan-image.sh" "${image_id}" "${scan_dir}"

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
  } > "${scan_dir}/published-image.txt"
else
  echo "Immutable local image: ${image_id}"
  {
    printf 'LOCAL_IMAGE_ID=%s\n' "${image_id}"
    printf 'GIT_COMMIT=%s\n' "${revision}"
    printf 'BASE_IMAGE=%s\n' "${base_image}"
  } > "${scan_dir}/published-image.txt"
fi
