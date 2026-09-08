#!/usr/bin/env bash
set -euo pipefail

if (($# != 2)); then
  echo "Usage: deploy/scan-image.sh IMAGE_REF OUTPUT_DIR" >&2
  exit 64
fi

readonly image_ref="$1"
readonly output_dir="$2"

for command_name in podman syft grype; do
  command -v "${command_name}" >/dev/null || {
    echo "ERROR: required scan command not found: ${command_name}" >&2
    exit 69
  }
done
podman image exists "${image_ref}" || {
  echo "ERROR: image is not present locally: ${image_ref}" >&2
  exit 66
}

if [[ -d "${output_dir}" && -n "$(find "${output_dir}" -mindepth 1 -print -quit)" ]]; then
  echo "ERROR: scan output directory must be empty: ${output_dir}" >&2
  exit 65
fi
mkdir -p -- "${output_dir}"
umask 022

image_id="$(podman image inspect --format '{{.Id}}' "${image_ref}")"
[[ "${image_id}" =~ ^sha256:[0-9a-f]{64}$ ]] || {
  echo "ERROR: image does not have a valid immutable image ID." >&2
  exit 65
}
safe_name="$(printf '%s' "${image_id}" | tr ':/' '__')"
readonly image_id safe_name
syft_version="$(syft --version 2>/dev/null)" || {
  echo "ERROR: could not determine the Syft version." >&2
  exit 69
}
grype_version="$(grype --version 2>/dev/null)" || {
  echo "ERROR: could not determine the Grype version." >&2
  exit 69
}
[[ "${syft_version}" == syft\ * && "${syft_version}" != *$'\n'* &&
   "${grype_version}" == grype\ * && "${grype_version}" != *$'\n'* ]] || {
  echo "ERROR: Syft/Grype returned an unusable version string." >&2
  exit 69
}
readonly syft_version grype_version
readonly sbom="${output_dir}/${safe_name}.cyclonedx.json"
readonly report="${output_dir}/${safe_name}.grype.json"
readonly grype_db_status="${output_dir}/${safe_name}.grype-db.json"
readonly r_packages="${output_dir}/${safe_name}.R-packages.tsv"
readonly python_packages="${output_dir}/${safe_name}.python-packages.txt"
readonly attestation="${output_dir}/${safe_name}.scan.txt"

syft_status=0
syft "podman:${image_id}" -o "cyclonedx-json=${sbom}" || syft_status=$?
if ((syft_status)); then
  echo "ERROR: SBOM generation failed. For a Podman connection error, start the private API socket and set CONTAINER_HOST as documented in deploy/README.md." >&2
  exit "${syft_status}"
fi
# --fail-on high rejects both High and Critical findings.
grype "sbom:${sbom}" --fail-on high -o json > "${report}"
if ! grype db status -o json > "${grype_db_status}" 2>/dev/null; then
  rm -f -- "${grype_db_status}"
fi
podman run --rm --entrypoint /bin/cat "${image_id}" \
  /opt/islet-explorer-build/R-packages.tsv > "${r_packages}"
podman run --rm --entrypoint /opt/islet-explorer-py/bin/python "${image_id}" \
  -m pip freeze --all > "${python_packages}"

{
  printf 'IMAGE_REF=%s\n' "${image_ref}"
  printf 'IMAGE_ID=%s\n' "${image_id}"
  printf 'SCANNED_AT=%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
  printf 'SYFT_VERSION=%s\n' "${syft_version}"
  printf 'GRYPE_VERSION=%s\n' "${grype_version}"
  printf 'RESULT=pass-no-high-or-critical\n'
  printf 'SBOM_SHA256=%s\n' "$(sha256sum "${sbom}" | awk '{print $1}')"
  printf 'REPORT_SHA256=%s\n' "$(sha256sum "${report}" | awk '{print $1}')"
  printf 'R_PACKAGES_SHA256=%s\n' "$(sha256sum "${r_packages}" | awk '{print $1}')"
  printf 'PYTHON_PACKAGES_SHA256=%s\n' "$(sha256sum "${python_packages}" | awk '{print $1}')"
  if [[ -f "${grype_db_status}" ]]; then
    printf 'GRYPE_DB_STATUS_SHA256=%s\n' "$(sha256sum "${grype_db_status}" | awk '{print $1}')"
  fi
} > "${attestation}"

echo "Security gate passed; artifacts written to ${output_dir}."
