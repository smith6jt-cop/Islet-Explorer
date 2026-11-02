#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: deploy/generate-release-manifest.sh RELEASE_DIR RELEASE IMAGE_REF GIT_COMMIT [DATA_SCHEMA [EXPECTED_IMAGES]]

RELEASE_DIR must contain the generated data bundle at data/, including its
data/images/ subtree, and the passing image evidence at attestations/. The
command writes release.env and MANIFEST.sha256, invalidating any prior
verification stamp. IMAGE_REF must be an immutable registry digest or local
image ID and must match the published scan evidence. EXPECTED_IMAGES is the
approved OME-TIFF count, defaults to the original cohort's 15, and must match
the staged files. Use 20 with DATA_SCHEMA 5 for the uniform 20-donor cohort;
never infer approval from the number of files present.
EOF
}

if (($# < 4 || $# > 6)); then
  usage >&2
  exit 64
fi

release_dir="$(realpath -e -- "$1")"
release="$2"
image_ref="$3"
git_commit="$4"
data_schema="${5:-1}"
expected_images="${6-15}"

[[ "${expected_images}" =~ ^[1-9][0-9]*$ ]] || {
  echo "ERROR: expected image count must be a positive integer without leading zeros." >&2
  exit 64
}

[[ "${release}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
  echo "ERROR: invalid release identifier: ${release}" >&2
  exit 64
}
[[ "${image_ref}" =~ ^sha256:[0-9a-f]{64}$ ||
   "${image_ref}" =~ ^[A-Za-z0-9._:/-]+@sha256:[0-9a-f]{64}$ ]] || {
  echo "ERROR: IMAGE_REF must be immutable (@sha256 digest or local sha256 image ID)." >&2
  exit 64
}
[[ "${git_commit}" =~ ^[0-9a-fA-F]{40}$ ]] || {
  echo "ERROR: invalid image reference or 40-character Git commit." >&2
  exit 64
}
[[ "${data_schema}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
  echo "ERROR: invalid data schema identifier." >&2
  exit 64
}
for required_dir in data data/images attestations; do
  [[ -d "${release_dir}/${required_dir}" && ! -L "${release_dir}/${required_dir}" ]] || {
    echo "ERROR: ${release_dir}/${required_dir} must be a real directory." >&2
    exit 66
  }
done
[[ -f "${release_dir}/data/manifest.json" ]] || {
  echo "ERROR: generated data bundle is missing data/manifest.json." >&2
  exit 66
}

read_evidence_value() {
  local file="$1"
  local key="$2"
  awk -F= -v key="${key}" '
    $1 == key { sub(/^[^=]*=/, ""); value = $0; count++ }
    END { if (count != 1) exit 1; print value }
  ' "${file}"
}

read_optional_evidence_value() {
  local file="$1"
  local key="$2"
  awk -F= -v key="${key}" '
    $1 == key { sub(/^[^=]*=/, ""); value = $0; count++ }
    END { if (count > 1) exit 1; if (count == 1) print value }
  ' "${file}"
}

required_evidence_value() {
  local file="$1"
  local key="$2"
  local value
  value="$(read_evidence_value "${file}" "${key}")" || {
    echo "ERROR: evidence file has missing or duplicate ${key}: ${file}" >&2
    return 65
  }
  [[ -n "${value}" ]] || {
    echo "ERROR: evidence field ${key} is empty: ${file}" >&2
    return 65
  }
  printf '%s' "${value}"
}

attestations_dir="${release_dir}/attestations"
mapfile -d '' scan_attestations < <(
  find "${attestations_dir}" -maxdepth 1 -type f -name '*.scan.txt' -print0
)
(( ${#scan_attestations[@]} == 1 )) || {
  echo "ERROR: attestations must contain exactly one regular *.scan.txt file." >&2
  exit 66
}
scan_attestation="${scan_attestations[0]}"
published_image="${attestations_dir}/published-image.txt"
[[ -f "${published_image}" && ! -L "${published_image}" ]] || {
  echo "ERROR: attestations/published-image.txt must be a regular file." >&2
  exit 66
}

scan_result="$(required_evidence_value "${scan_attestation}" RESULT)" || exit $?
scan_ref="$(required_evidence_value "${scan_attestation}" IMAGE_REF)" || exit $?
scan_image_id="$(required_evidence_value "${scan_attestation}" IMAGE_ID)" || exit $?
syft_version="$(required_evidence_value "${scan_attestation}" SYFT_VERSION)" || exit $?
grype_version="$(required_evidence_value "${scan_attestation}" GRYPE_VERSION)" || exit $?
[[ "${scan_result}" == pass-no-high-or-critical || "${scan_result}" == pass-with-duckdb-remediation-window ]] || {
  echo "ERROR: image evidence does not record a passing High/Critical gate." >&2
  exit 65
}
[[ "${scan_image_id}" =~ ^sha256:[0-9a-f]{64}$ && "${scan_ref}" == "${scan_image_id}" ]] || {
  echo "ERROR: scan evidence is not bound to its immutable image ID." >&2
  exit 65
}
[[ "${syft_version}" == syft\ * && "${grype_version}" == grype\ * ]] || {
  echo "ERROR: scan evidence does not contain usable Syft/Grype versions." >&2
  exit 65
}

safe_name="$(printf '%s' "${scan_image_id}" | tr ':/' '__')"
expected_scan_attestation="${attestations_dir}/${safe_name}.scan.txt"
[[ "${scan_attestation}" == "${expected_scan_attestation}" ]] || {
  echo "ERROR: scan-attestation filename does not match its IMAGE_ID." >&2
  exit 65
}

sbom="${attestations_dir}/${safe_name}.cyclonedx.json"
raw_sbom="${attestations_dir}/${safe_name}.syft.cyclonedx.json"
sbom_corrections="${attestations_dir}/${safe_name}.sbom-corrections.json"
report="${attestations_dir}/${safe_name}.grype.json"
r_packages="${attestations_dir}/${safe_name}.R-packages.tsv"
python_packages="${attestations_dir}/${safe_name}.python-packages.txt"
for evidence_file in "${sbom}" "${raw_sbom}" "${sbom_corrections}" "${report}" "${r_packages}" "${python_packages}"; do
  [[ -f "${evidence_file}" && ! -L "${evidence_file}" ]] || {
    echo "ERROR: required image evidence is missing or is a symlink: ${evidence_file}" >&2
    exit 66
  }
done

verify_evidence_hash() {
  local key="$1"
  local file="$2"
  local expected actual
  expected="$(required_evidence_value "${scan_attestation}" "${key}")" || return $?
  [[ "${expected}" =~ ^[0-9a-f]{64}$ ]] || {
    echo "ERROR: invalid ${key} in scan evidence." >&2
    return 65
  }
  actual="$(sha256sum "${file}" | awk '{print $1}')"
  [[ "${actual}" == "${expected}" ]] || {
    echo "ERROR: ${key} does not match ${file}." >&2
    return 65
  }
}

verify_evidence_hash SBOM_SHA256 "${sbom}" || exit $?
verify_evidence_hash RAW_SBOM_SHA256 "${raw_sbom}" || exit $?
verify_evidence_hash SBOM_CORRECTIONS_SHA256 "${sbom_corrections}" || exit $?
verify_evidence_hash REPORT_SHA256 "${report}" || exit $?
verify_evidence_hash R_PACKAGES_SHA256 "${r_packages}" || exit $?
verify_evidence_hash PYTHON_PACKAGES_SHA256 "${python_packages}" || exit $?

# Every newly generated manifest must carry the current supplementary gate.
# Existing historical manifests remain historical evidence; old scan receipts
# cannot authorize a new release without checking the actual bundled decoders.
brotli_report="${attestations_dir}/${safe_name}.brotli.json"
brotli_hash="$(read_optional_evidence_value "${scan_attestation}" BROTLI_REPORT_SHA256)" || exit 65
[[ -n "${brotli_hash}" && -f "${brotli_report}" && ! -L "${brotli_report}" ]] || {
  echo "ERROR: a current Brotli verification report is required; rescan the image." >&2
  exit 66
}
verify_evidence_hash BROTLI_REPORT_SHA256 "${brotli_report}" || exit $?
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
decision_args=("${brotli_report}")
if [[ "${scan_result}" == pass-with-duckdb-remediation-window ]]; then
  duckdb_deferral="${attestations_dir}/${safe_name}.duckdb-deferral.json"
  [[ -f "${duckdb_deferral}" && ! -L "${duckdb_deferral}" ]] || {
    echo "ERROR: DuckDB remediation-window record is required." >&2
    exit 66
  }
  verify_evidence_hash DUCKDB_DEFERRAL_SHA256 "${duckdb_deferral}" || exit $?
  cmp -s "${duckdb_deferral}" "${script_dir}/duckdb-remediation-window.json" || {
    echo "ERROR: DuckDB deferral does not match the recorded owner direction." >&2
    exit 65
  }
  [[ "$(required_evidence_value "${scan_attestation}" PACKAGE_RESULT)" == pass-no-high-or-critical &&
     "$(required_evidence_value "${scan_attestation}" DUCKDB_PATCH_BY)" == 2026-11-10 ]] || exit 65
  decision_args+=("${duckdb_deferral}")
fi
evaluated_result="$(python3 "${script_dir}/evaluate-brotli-report.py" "${decision_args[@]}")"
[[ "${evaluated_result}" == "${scan_result}" ]] || {
  echo "ERROR: scan result disagrees with the measured report and remediation window." >&2
  exit 65
}

grype_db_status="${attestations_dir}/${safe_name}.grype-db.json"
grype_db_hash="$(read_optional_evidence_value "${scan_attestation}" GRYPE_DB_STATUS_SHA256)" || {
  echo "ERROR: duplicate GRYPE_DB_STATUS_SHA256 in scan evidence." >&2
  exit 65
}
if [[ -e "${grype_db_status}" || -n "${grype_db_hash}" ]]; then
  [[ -f "${grype_db_status}" && ! -L "${grype_db_status}" &&
     "${grype_db_hash}" =~ ^[0-9a-f]{64}$ ]] || {
    echo "ERROR: Grype DB status evidence is incomplete or malformed." >&2
    exit 65
  }
  [[ "$(sha256sum "${grype_db_status}" | awk '{print $1}')" == "${grype_db_hash}" ]] || {
    echo "ERROR: GRYPE_DB_STATUS_SHA256 does not match ${grype_db_status}." >&2
    exit 65
  }
fi

published_local_id="$(read_optional_evidence_value "${published_image}" LOCAL_IMAGE_ID)" || {
  echo "ERROR: duplicate LOCAL_IMAGE_ID in published-image.txt." >&2
  exit 65
}
published_registry_ref="$(read_optional_evidence_value "${published_image}" REGISTRY_REF)" || {
  echo "ERROR: duplicate REGISTRY_REF in published-image.txt." >&2
  exit 65
}
published_image_id="$(read_optional_evidence_value "${published_image}" IMAGE_ID)" || {
  echo "ERROR: duplicate IMAGE_ID in published-image.txt." >&2
  exit 65
}
published_commit="$(required_evidence_value "${published_image}" GIT_COMMIT)" || exit $?
published_base="$(required_evidence_value "${published_image}" BASE_IMAGE)" || exit $?

if [[ -n "${published_local_id}" && -z "${published_registry_ref}" && -z "${published_image_id}" ]]; then
  published_ref="${published_local_id}"
  published_id="${published_local_id}"
elif [[ -z "${published_local_id}" && -n "${published_registry_ref}" && -n "${published_image_id}" ]]; then
  published_ref="${published_registry_ref}"
  published_id="${published_image_id}"
else
  echo "ERROR: published-image.txt must describe exactly one local or registry image." >&2
  exit 65
fi

[[ "${published_id}" == "${scan_image_id}" && "${published_ref}" == "${image_ref}" &&
   "${published_commit}" == "${git_commit,,}" &&
   "${published_base}" =~ ^[A-Za-z0-9._:/-]+@sha256:[0-9a-f]{64}$ ]] || {
  echo "ERROR: published image, scan evidence, and release identity do not match." >&2
  exit 65
}

unsafe="$(find "${release_dir}" -xdev \( \
  -name '.Renviron' -o -name '.Rhistory' -o -name '*.bak' -o \
  -name '*~' -o -name '*.pem' -o -name '*.key' \
\) -print -quit)"
if [[ -n "${unsafe}" ]]; then
  echo "ERROR: release contains a forbidden credential/backup filename: ${unsafe}" >&2
  exit 65
fi

attestation_link="$(find "${attestations_dir}" -type l -print -quit)"
if [[ -n "${attestation_link}" ]]; then
  echo "ERROR: attestation evidence may not contain symlinks: ${attestation_link}" >&2
  exit 65
fi

# Symlinks are allowed for TIFF offset aliases, but none may escape the release.
while IFS= read -r -d '' link_path; do
  resolved="$(realpath -e -- "${link_path}")" || {
    echo "ERROR: broken release symlink: ${link_path}" >&2
    exit 65
  }
  case "${resolved}" in
    "${release_dir}"/*) ;;
    *) echo "ERROR: release symlink escapes its release: ${link_path} -> ${resolved}" >&2; exit 65 ;;
  esac
done < <(find "${release_dir}/data" -type l -print0)

data_files="$(find "${release_dir}/data" -path "${release_dir}/data/images" -prune -o -type f -print | wc -l)"
actual_images="$(find "${release_dir}/data/images" -type f \
  \( -iname '*.ome.tif' -o -iname '*.ome.tiff' \) | wc -l)"
[[ "${actual_images}" == "${expected_images}" ]] || {
  echo "ERROR: expected ${expected_images} OME-TIFF files, found ${actual_images}." >&2
  exit 65
}
image_files="$(find "${release_dir}/data/images" -type f | wc -l)"
data_total_bytes="$(du -sb -- "${release_dir}/data" | awk '{print $1}')"
image_bytes="$(du -sb -- "${release_dir}/data/images" | awk '{print $1}')"
data_bytes="$((data_total_bytes - image_bytes))"

umask 022
env_tmp="$(mktemp "${release_dir}/.release.env.XXXXXX")"
manifest_tmp="$(mktemp "${release_dir}/.MANIFEST.sha256.XXXXXX")"
trap 'rm -f -- "${env_tmp}" "${manifest_tmp}"' EXIT

{
  printf 'MANIFEST_VERSION=2\n'
  printf 'RELEASE=%s\n' "${release}"
  printf 'IMAGE_REF=%s\n' "${image_ref}"
  printf 'GIT_COMMIT=%s\n' "${git_commit,,}"
  printf 'DATA_SCHEMA=%s\n' "${data_schema}"
  printf 'SECURITY_RESULT=%s\n' "${scan_result}"
  if [[ "${scan_result}" == pass-with-duckdb-remediation-window ]]; then
    printf 'DUCKDB_PATCH_BY=2026-11-10\n'
  fi
  printf 'EXPECTED_IMAGES=%s\n' "${expected_images}"
  printf 'CREATED_AT=%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
  printf 'DATA_FILES=%s\n' "${data_files}"
  printf 'DATA_BYTES=%s\n' "${data_bytes}"
  printf 'IMAGE_FILES=%s\n' "${image_files}"
  printf 'IMAGE_BYTES=%s\n' "${image_bytes}"
} > "${env_tmp}"
mv -f -- "${env_tmp}" "${release_dir}/release.env"

paths=(data release.env attestations)
(
  cd -- "${release_dir}"
  find "${paths[@]}" \( -type f -o -type l \) -print0 | \
    LC_ALL=C sort -z | xargs -0 -r sha256sum
) > "${manifest_tmp}"
mv -f -- "${manifest_tmp}" "${release_dir}/MANIFEST.sha256"
rm -f -- "${release_dir}/.verified"
trap - EXIT

echo "Generated ${release_dir}/MANIFEST.sha256"
echo "Run deploy/verify-release.sh '${release_dir}' before activation."
