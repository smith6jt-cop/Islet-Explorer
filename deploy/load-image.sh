#!/usr/bin/env bash
set -euo pipefail

if (($# != 1)); then
  echo "Usage: deploy/load-image.sh IMAGE.oci.tar" >&2
  exit 64
fi
archive="$(realpath -e -- "$1")"
metadata="${archive}.metadata"
checksums="${archive}.sha256"

command -v podman >/dev/null || { echo "ERROR: podman is required." >&2; exit 69; }
[[ -f "${metadata}" && -f "${checksums}" ]] || {
  echo "ERROR: archive metadata/checksum sidecars are missing." >&2
  exit 66
}

archive_dir="$(dirname -- "${archive}")"
(
  cd -- "${archive_dir}"
  sha256sum --check --strict "$(basename -- "${checksums}")"
)

read_meta() {
  local key="$1"
  awk -F= -v key="${key}" '$1 == key {sub(/^[^=]*=/, ""); print; exit}' "${metadata}"
}
format="$(read_meta FORMAT)"
expected_id="$(read_meta IMAGE_ID)"
expected_commit="$(read_meta GIT_COMMIT)"
[[ "${format}" == islet-explorer-oci-transfer-v1 &&
   "${expected_id}" =~ ^sha256:[0-9a-f]{64}$ &&
   "${expected_commit}" =~ ^[0-9a-f]{40}$ ]] || {
  echo "ERROR: transfer metadata is malformed." >&2
  exit 65
}

podman load --input "${archive}"
podman image exists "${expected_id}" || {
  echo "ERROR: loaded archive did not produce expected image ${expected_id}." >&2
  exit 65
}
actual_id="$(podman image inspect --format '{{.Id}}' "${expected_id}")"
actual_commit="$(podman image inspect --format '{{ index .Labels "org.opencontainers.image.revision" }}' "${expected_id}")"
[[ "${actual_id}" == "${expected_id}" && "${actual_commit}" == "${expected_commit}" ]] || {
  echo "ERROR: loaded image identity/revision does not match transfer metadata." >&2
  exit 65
}

echo "Verified and loaded immutable image: ${expected_id}"
