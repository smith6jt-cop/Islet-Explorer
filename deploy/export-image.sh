#!/usr/bin/env bash
set -euo pipefail

if (($# != 2)); then
  echo "Usage: deploy/export-image.sh IMAGE_REF OUTPUT.oci.tar" >&2
  exit 64
fi
image_ref="$1"
output="$2"

command -v podman >/dev/null || { echo "ERROR: podman is required." >&2; exit 69; }
[[ "${output}" != *$'\n'* ]] || { echo "ERROR: invalid output path." >&2; exit 64; }
output_dir="$(realpath -e -- "$(dirname -- "${output}")")"
output="${output_dir}/$(basename -- "${output}")"
metadata="${output}.metadata"
checksums="${output}.sha256"
for target in "${output}" "${metadata}" "${checksums}"; do
  [[ ! -e "${target}" ]] || { echo "ERROR: refusing to overwrite ${target}." >&2; exit 65; }
done
podman image exists "${image_ref}" || { echo "ERROR: image not found: ${image_ref}" >&2; exit 66; }

image_id="$(podman image inspect --format '{{.Id}}' "${image_ref}")"
# Keep archive metadata canonical across Podman's inspect output formats.
image_id="sha256:${image_id#sha256:}"
[[ "${image_id}" =~ ^sha256:[0-9a-f]{64}$ ]] || {
  echo "ERROR: image lacks a valid immutable ID." >&2
  exit 65
}
# Resolve all labels and archive content through the ID captured above. The
# caller may supply a convenient tag, but a concurrent retag cannot change the
# image that is exported after this point.
revision="$(podman image inspect --format '{{ index .Labels "org.opencontainers.image.revision" }}' "${image_id}")"
base_image="$(podman image inspect --format '{{ index .Labels "org.opencontainers.image.base.name" }}' "${image_id}")"
[[ "${revision}" =~ ^[0-9a-f]{40}$ &&
   "${base_image}" =~ ^[A-Za-z0-9._:/-]+@sha256:[0-9a-f]{64}$ ]] || {
  echo "ERROR: image lacks a valid source revision or digest-pinned base-image label." >&2
  exit 65
}

archive_tmp="$(mktemp "${output_dir}/.$(basename -- "${output}").XXXXXX")"
metadata_tmp="$(mktemp "${output_dir}/.$(basename -- "${metadata}").XXXXXX")"
checksum_tmp="$(mktemp "${output_dir}/.$(basename -- "${checksums}").XXXXXX")"
trap 'rm -f -- "${archive_tmp}" "${metadata_tmp}" "${checksum_tmp}"' EXIT

podman save --format oci-archive --output "${archive_tmp}" "${image_id}"
{
  printf 'FORMAT=islet-explorer-oci-transfer-v1\n'
  printf 'IMAGE_ID=%s\n' "${image_id}"
  printf 'SOURCE_REF=%s\n' "${image_ref}"
  printf 'GIT_COMMIT=%s\n' "${revision}"
  printf 'BASE_IMAGE=%s\n' "${base_image}"
  printf 'EXPORTED_AT=%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
} > "${metadata_tmp}"

mv -- "${archive_tmp}" "${output}"
mv -- "${metadata_tmp}" "${metadata}"
(
  cd -- "${output_dir}"
  sha256sum "$(basename -- "${output}")" "$(basename -- "${metadata}")"
) > "${checksum_tmp}"
mv -- "${checksum_tmp}" "${checksums}"
trap - EXIT

echo "Exported ${image_id} to ${output}."
echo "Transfer and retain: ${output}, ${metadata}, ${checksums}"
