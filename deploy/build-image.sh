#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: deploy/build-image.sh --tag IMAGE:TAG --version RELEASE \
         --base-image IMAGE@sha256:DIGEST [--push] [--scan-dir DIR]

Builds the data-free OCI image from a clean Git release. A digest-pinned base is
required. With --push, the immutable registry reference is printed at the end.
With --scan-dir, syft and grype are required and an SBOM plus a High/Critical
vulnerability gate are written there.
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

[[ -n "${tag}" && -n "${version}" && -n "${base_image}" ]] || {
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

for command_name in git podman; do
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

podman build \
  --pull=missing \
  --ignorefile "${script_dir}/containerignore" \
  --file "${script_dir}/Containerfile" \
  --tag "${tag}" \
  --build-arg "BASE_IMAGE=${base_image}" \
  --build-arg "OCI_CREATED=${created}" \
  --build-arg "OCI_REVISION=${revision}" \
  --build-arg "OCI_VERSION=${version}" \
  "${repo_root}"

podman run --rm --entrypoint /opt/islet-explorer-py/bin/python "${tag}" -c \
  'import anndata, scanpy, igraph, leidenalg, fa2_modified; print("python dependency probe: ok")'
podman run --rm --entrypoint /usr/bin/Rscript "${tag}" -e \
  'stopifnot(packageVersion("destiny")=="3.22.0", packageVersion("duckdb")=="1.5.2"); cat("R dependency probe: ok\n")'
podman run --rm --entrypoint /usr/bin/Rscript \
  --workdir /srv/islet-explorer/app/shiny_app "${tag}" tests/parse_check.R

if [[ -n "${scan_dir}" ]]; then
  "${script_dir}/scan-image.sh" "${tag}" "${scan_dir}"
fi

image_id="$(podman image inspect --format '{{.Id}}' "${tag}")"
echo "Built ${tag} as ${image_id} from Git ${revision}."

if ((push_image)); then
  digest_file="$(mktemp)"
  trap 'rm -f -- "${digest_file}"' EXIT
  podman push --digestfile "${digest_file}" "${tag}"
  digest="$(tr -d '[:space:]' < "${digest_file}")"
  repository="${tag%:*}"
  [[ "${repository}" != "${tag}" && "${digest}" == sha256:* ]] || {
    echo "ERROR: could not derive immutable registry reference after push." >&2
    exit 70
  }
  echo "Immutable image: ${repository}@${digest}"
  if [[ -n "${scan_dir}" ]]; then
    {
      printf 'REGISTRY_REF=%s@%s\n' "${repository}" "${digest}"
      printf 'IMAGE_ID=%s\n' "${image_id}"
      printf 'GIT_COMMIT=%s\n' "${revision}"
      printf 'BASE_IMAGE=%s\n' "${base_image}"
    } > "${scan_dir}/published-image.txt"
  fi
else
  echo "Immutable local image: ${image_id}"
  if [[ -n "${scan_dir}" ]]; then
    {
      printf 'LOCAL_IMAGE_ID=%s\n' "${image_id}"
      printf 'GIT_COMMIT=%s\n' "${revision}"
      printf 'BASE_IMAGE=%s\n' "${base_image}"
    } > "${scan_dir}/published-image.txt"
  fi
fi
