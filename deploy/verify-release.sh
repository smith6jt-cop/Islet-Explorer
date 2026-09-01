#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: deploy/verify-release.sh RELEASE_DIR [--quick] [--expected-images N]" >&2
}

if (($# < 1)); then usage; exit 64; fi
release_dir="$(realpath -e -- "$1")"
shift
quick=0
expected_images=15
while (($#)); do
  case "$1" in
    --quick) quick=1; shift ;;
    --expected-images) expected_images="${2:?missing value}"; shift 2 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 64 ;;
  esac
done
[[ "${expected_images}" =~ ^[0-9]+$ ]] || { echo "ERROR: expected image count must be numeric." >&2; exit 64; }

for required in data data/images data/manifest.json release.env MANIFEST.sha256; do
  [[ -e "${release_dir}/${required}" ]] || {
    echo "ERROR: release is missing ${required}." >&2
    exit 66
  }
done

read_meta() {
  local key="$1"
  awk -F= -v key="${key}" '$1 == key {sub(/^[^=]*=/, ""); print; exit}' "${release_dir}/release.env"
}

release="$(read_meta RELEASE)"
image_ref="$(read_meta IMAGE_REF)"
git_commit="$(read_meta GIT_COMMIT)"
data_schema="$(read_meta DATA_SCHEMA)"
[[ "${release}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ &&
   ( "${image_ref}" =~ ^sha256:[0-9a-f]{64}$ ||
     "${image_ref}" =~ ^[A-Za-z0-9._:/-]+@sha256:[0-9a-f]{64}$ ) &&
   "${git_commit}" =~ ^[0-9a-f]{40}$ && -n "${data_schema}" ]] || {
  echo "ERROR: release.env is incomplete or malformed." >&2
  exit 65
}

actual_images="$(find "${release_dir}/data/images" -type f \
  \( -iname '*.ome.tif' -o -iname '*.ome.tiff' \) | wc -l)"
if [[ "${actual_images}" != "${expected_images}" ]]; then
  echo "ERROR: expected ${expected_images} OME-TIFF files, found ${actual_images}." >&2
  exit 65
fi
[[ -n "$(find "${release_dir}/data" -type f -print -quit)" ]] || {
  echo "ERROR: release data directory is empty." >&2
  exit 65
}

manifest_hash="$(sha256sum "${release_dir}/MANIFEST.sha256" | awk '{print $1}')"
if ((quick)); then
  [[ -f "${release_dir}/.verified" ]] || {
    echo "ERROR: no full-verification stamp; run verify-release.sh without --quick." >&2
    exit 65
  }
  stamped_hash="$(awk -F= '$1 == "MANIFEST_SHA256" {print $2; exit}' "${release_dir}/.verified")"
  [[ "${stamped_hash}" == "${manifest_hash}" ]] || {
    echo "ERROR: release manifest changed after full verification." >&2
    exit 65
  }
  echo "Release ${release} has a current full-verification stamp."
  exit 0
fi

(
  cd -- "${release_dir}"
  sha256sum --check --strict MANIFEST.sha256
)

stamp_tmp="$(mktemp "${release_dir}/.verified.XXXXXX")"
trap 'rm -f -- "${stamp_tmp}"' EXIT
{
  printf 'MANIFEST_SHA256=%s\n' "${manifest_hash}"
  printf 'VERIFIED_AT=%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
  printf 'VERIFIED_HOST=%s\n' "$(hostname -f 2>/dev/null || hostname)"
  printf 'VERIFIED_UID=%s\n' "$(id -u)"
} > "${stamp_tmp}"
mv -f -- "${stamp_tmp}" "${release_dir}/.verified"
trap - EXIT

echo "Release ${release} verified: ${actual_images} OME-TIFFs; manifest ${manifest_hash}."
