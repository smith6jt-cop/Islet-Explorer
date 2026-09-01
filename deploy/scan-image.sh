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
safe_name="$(printf '%s' "${image_id}" | tr ':/' '__')"
readonly image_id safe_name
readonly sbom="${output_dir}/${safe_name}.cyclonedx.json"
readonly report="${output_dir}/${safe_name}.grype.json"
readonly r_packages="${output_dir}/${safe_name}.R-packages.tsv"
readonly python_packages="${output_dir}/${safe_name}.python-packages.txt"
readonly attestation="${output_dir}/${safe_name}.scan.txt"

syft "podman:${image_ref}" -o "cyclonedx-json=${sbom}"
# --fail-on high rejects both High and Critical findings.
grype "podman:${image_ref}" --fail-on high -o json > "${report}"
podman run --rm --entrypoint /bin/cat "${image_ref}" \
  /opt/islet-explorer-build/R-packages.tsv > "${r_packages}"
podman run --rm --entrypoint /opt/islet-explorer-py/bin/python "${image_ref}" \
  -m pip freeze --all > "${python_packages}"

{
  printf 'IMAGE_REF=%s\n' "${image_ref}"
  printf 'IMAGE_ID=%s\n' "${image_id}"
  printf 'SCANNED_AT=%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
  printf 'SYFT_VERSION=%s\n' "$(syft version 2>/dev/null | head -n 1)"
  printf 'GRYPE_VERSION=%s\n' "$(grype version 2>/dev/null | head -n 1)"
  printf 'RESULT=pass-no-high-or-critical\n'
  printf 'SBOM_SHA256=%s\n' "$(sha256sum "${sbom}" | awk '{print $1}')"
  printf 'REPORT_SHA256=%s\n' "$(sha256sum "${report}" | awk '{print $1}')"
  printf 'R_PACKAGES_SHA256=%s\n' "$(sha256sum "${r_packages}" | awk '{print $1}')"
  printf 'PYTHON_PACKAGES_SHA256=%s\n' "$(sha256sum "${python_packages}" | awk '{print $1}')"
} > "${attestation}"

echo "Security gate passed; artifacts written to ${output_dir}."
