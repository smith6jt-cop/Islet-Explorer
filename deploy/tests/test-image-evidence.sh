#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd -- "${script_dir}/../.." && pwd -P)"
work_dir="$(mktemp -d)"
trap 'rm -rf -- "${work_dir}"' EXIT

fail() {
  echo "FAIL: $*" >&2
  exit 1
}

expect_status() {
  local expected="$1"
  local actual
  shift
  set +e
  "$@" > "${work_dir}/last-command.log" 2>&1
  actual=$?
  set -e
  if [[ "${actual}" != "${expected}" ]]; then
    cat "${work_dir}/last-command.log" >&2
    fail "expected status ${expected}, got ${actual}: $*"
  fi
}

mock_bin="${work_dir}/bin"
mock_log="${work_dir}/mock-log"
mkdir -p -- "${mock_bin}" "${mock_log}"

cat > "${mock_bin}/podman" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
case "${1:-} ${2:-}" in
  "image exists") exit 0 ;;
  "image inspect")
    case "${4:-}" in
      "{{.Id}}") printf '%s\n' "${MOCK_INSPECT_IMAGE_ID-${MOCK_IMAGE_ID:?}}" ;;
      *org.opencontainers.image.revision*) printf '%s\n' "${MOCK_GIT_COMMIT:?}" ;;
      *org.opencontainers.image.base.name*) printf '%s\n' "${MOCK_BASE_IMAGE:?}" ;;
      *) echo "unexpected mock inspect format: ${4:-}" >&2; exit 70 ;;
    esac
    exit 0
    ;;
esac
if [[ "${1:-}" == build ]]; then
  printf '%s\n' "$*" >> "${MOCK_LOG_DIR:?}/podman-build.log"
  while (($#)); do
    if [[ "$1" == --iidfile ]]; then
      printf '%s\n' "${MOCK_IMAGE_ID:?}" > "$2"
      exit 0
    fi
    shift
  done
  echo "mock build did not receive --iidfile" >&2
  exit 64
fi
if [[ "${1:-}" == run ]]; then
  printf '%s\n' "$*" >> "${MOCK_LOG_DIR:?}/podman-run.log"
  case " $* " in
    *" --entrypoint /usr/bin/Rscript "*) echo 'Rscript must resolve through PATH' >&2; exit 70 ;;
    *" --pull=missing "*) if [[ "${MOCK_BASE_FAIL:-0}" == 1 ]]; then exit 65; fi ;;
    *" /bin/cat "*) printf 'Package\tVersion\nshiny\t1.0\n' ;;
    *" pip freeze --all "*) printf 'scanpy==1.0\n' ;;
    *) : ;;
  esac
  exit 0
fi
if [[ "${1:-}" == push ]]; then
  printf '%s\n' "$*" >> "${MOCK_LOG_DIR:?}/podman-push.log"
  while (($#)); do
    if [[ "$1" == --digestfile ]]; then
      printf '%s\n' "${MOCK_REGISTRY_DIGEST:?}" > "$2"
      exit 0
    fi
    shift
  done
  exit 64
fi
if [[ "${1:-}" == save ]]; then
  printf '%s\n' "$*" >> "${MOCK_LOG_DIR:?}/podman-save.log"
  while (($#)); do
    if [[ "$1" == --output ]]; then
      printf 'mock OCI archive\n' > "$2"
      exit 0
    fi
    shift
  done
  exit 64
fi
if [[ "${1:-}" == load ]]; then
  printf '%s\n' "$*" >> "${MOCK_LOG_DIR:?}/podman-load.log"
  exit 0
fi
echo "unexpected mock podman invocation: $*" >&2
exit 70
MOCK

cat > "${mock_bin}/git" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
case " $* " in
  *" status --porcelain --untracked-files=all "*) exit 0 ;;
  *" rev-parse HEAD "*) printf '%s\n' "${MOCK_GIT_COMMIT:?}"; exit 0 ;;
esac
echo "unexpected mock git invocation: $*" >&2
exit 70
MOCK

cat > "${mock_bin}/syft" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == --version ]]; then
  echo 'syft 1.42.3'
  exit 0
fi
printf '%s\n' "$*" >> "${MOCK_LOG_DIR:?}/syft.log"
output=""
while (($#)); do
  if [[ "$1" == -o ]]; then
    output="${2#cyclonedx-json=}"
    break
  fi
  shift
done
[[ -n "${output}" ]] || exit 64
printf '{"bomFormat":"CycloneDX","components":[]}\n' > "${output}"
MOCK

cat > "${mock_bin}/grype" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == --version ]]; then
  echo 'grype 0.110.0'
  exit 0
fi
if [[ "${1:-}" == db && "${2:-}" == status ]]; then
  printf '{"status":"valid","built":"2026-09-01T00:00:00Z"}\n'
  exit 0
fi
printf '%s\n' "$*" >> "${MOCK_LOG_DIR:?}/grype.log"
if [[ "${1:-}" != sbom:* ]]; then
  echo "expected Grype SBOM input, got: ${1:-}" >&2
  exit 64
fi
if [[ "${MOCK_GRYPE_FAIL:-0}" == 1 ]]; then
  printf '{"matches":[{"vulnerability":{"severity":"High"}}]}\n'
  exit 2
fi
printf '{"matches":[]}\n'
MOCK

chmod 0755 "${mock_bin}/podman" "${mock_bin}/git" \
  "${mock_bin}/syft" "${mock_bin}/grype"

image_hex="$(printf 'a%.0s' {1..64})"
image_id="sha256:${image_hex}"
safe_name="sha256_${image_hex}"
commit="$(printf 'b%.0s' {1..40})"
base_digest="$(printf 'c%.0s' {1..64})"
base_image="docker.io/rocker/shiny@sha256:${base_digest}"
registry_digest="sha256:$(printf 'e%.0s' {1..64})"
scan_dir="${work_dir}/scan-pass"

env PATH="${mock_bin}:${PATH}" MOCK_LOG_DIR="${mock_log}" MOCK_IMAGE_ID="${image_id}" \
  "${repo_root}/deploy/scan-image.sh" "${image_id}" "${scan_dir}"

attestation="${scan_dir}/${safe_name}.scan.txt"
sbom="${scan_dir}/${safe_name}.cyclonedx.json"
[[ -f "${attestation}" && -f "${scan_dir}/${safe_name}.grype-db.json" ]] || \
  fail "passing scan did not create its attestation and DB status"
grep -Fqx "SYFT_VERSION=syft 1.42.3" "${attestation}" || fail "Syft version was not recorded"
grep -Fqx "GRYPE_VERSION=grype 0.110.0" "${attestation}" || fail "Grype version was not recorded"
grep -Fq "podman:${image_id}" "${mock_log}/syft.log" || fail "Syft did not scan the immutable ID"
grep -Fq "sbom:${sbom}" "${mock_log}/grype.log" || fail "Grype did not scan the generated SBOM"
grep -Fq -- '--add-cpes-if-none=false' "${mock_log}/grype.log" || fail "Grype CPE regeneration was not disabled"
[[ -f "${scan_dir}/${safe_name}.syft.cyclonedx.json" && -f "${scan_dir}/${safe_name}.sbom-corrections.json" ]] || \
  fail "scan did not retain its original SBOM and correction record"
grep -Eq '^RAW_SBOM_SHA256=[0-9a-f]{64}$' "${attestation}" || fail "raw SBOM hash was not recorded"
grep -Eq '^SBOM_CORRECTIONS_SHA256=[0-9a-f]{64}$' "${attestation}" || fail "correction record hash was not recorded"

scan_fail_dir="${work_dir}/scan-fail"
expect_status 2 env PATH="${mock_bin}:${PATH}" MOCK_LOG_DIR="${mock_log}" \
  MOCK_IMAGE_ID="${image_id}" MOCK_GRYPE_FAIL=1 \
  "${repo_root}/deploy/scan-image.sh" "${image_id}" "${scan_fail_dir}"
[[ ! -e "${scan_fail_dir}/${safe_name}.scan.txt" ]] || \
  fail "failed Grype gate produced a passing attestation"
[[ -f "${scan_fail_dir}/${safe_name}.syft.cyclonedx.json" && -f "${scan_fail_dir}/${safe_name}.sbom-corrections.json" ]] || \
  fail "failed Grype gate lost its original SBOM or correction record"

expect_status 64 "${repo_root}/deploy/build-image.sh" \
  --tag local/test:release --version release \
  --base-image "${base_image}"

# A wrong or unusable base must fail before the expensive build and before any
# successful image evidence can be emitted.
expect_status 65 env PATH="${mock_bin}:${PATH}" MOCK_LOG_DIR="${mock_log}" \
  MOCK_IMAGE_ID="${image_id}" MOCK_GIT_COMMIT="${commit}" MOCK_BASE_FAIL=1 \
  "${repo_root}/deploy/build-image.sh" \
  --tag local/test:release --version release --base-image "${base_image}" \
  --scan-dir "${work_dir}/bad-base-scan"
[[ ! -e "${mock_log}/podman-build.log" && ! -e "${work_dir}/bad-base-scan" ]] || \
  fail 'failed base probe reached the build or emitted scan evidence'

build_scan_dir="${work_dir}/build-scan"
: > "${mock_log}/syft.log"
env PATH="${mock_bin}:${PATH}" MOCK_LOG_DIR="${mock_log}" \
  MOCK_IMAGE_ID="${image_id}" MOCK_GIT_COMMIT="${commit}" \
  "${repo_root}/deploy/build-image.sh" \
  --tag local/test:release --version release --base-image "${base_image}" \
  --scan-dir "${build_scan_dir}"
grep -Fq "podman:${image_id}" "${mock_log}/syft.log" || \
  fail "build wrapper did not scan the immutable image ID"
grep -Fq -- "--iidfile" "${mock_log}/podman-build.log" || \
  fail "build wrapper did not capture Podman's direct build image ID"
grep -Fq -- "--entrypoint Rscript" "${mock_log}/podman-run.log" || \
  fail "build wrapper did not resolve Rscript through the image PATH"
grep -Fqx "LOCAL_IMAGE_ID=${image_id}" "${build_scan_dir}/published-image.txt" || \
  fail "build wrapper did not publish immutable local-image evidence"

push_scan_dir="${work_dir}/build-push-scan"
push_tag="registry.example/islet-explorer:release"
env PATH="${mock_bin}:${PATH}" MOCK_LOG_DIR="${mock_log}" \
  MOCK_IMAGE_ID="${image_id}" MOCK_GIT_COMMIT="${commit}" \
  MOCK_BASE_IMAGE="${base_image}" MOCK_REGISTRY_DIGEST="${registry_digest}" \
  "${repo_root}/deploy/build-image.sh" \
  --tag "${push_tag}" --version release --base-image "${base_image}" \
  --scan-dir "${push_scan_dir}" --push
grep -Fq "${image_id} ${push_tag}" "${mock_log}/podman-push.log" || \
  fail "push did not use the scanned immutable image ID as its source"
grep -Fqx "REGISTRY_REF=registry.example/islet-explorer@${registry_digest}" \
  "${push_scan_dir}/published-image.txt" || \
  fail "push wrapper did not publish the registry digest"

export_dir="${work_dir}/export"
mkdir -p -- "${export_dir}"
archive="${export_dir}/islet-explorer.oci.tar"
env PATH="${mock_bin}:${PATH}" MOCK_LOG_DIR="${mock_log}" \
  MOCK_IMAGE_ID="${image_id}" MOCK_GIT_COMMIT="${commit}" \
  MOCK_BASE_IMAGE="${base_image}" \
  "${repo_root}/deploy/export-image.sh" local/test:release "${archive}"
grep -Fq "${image_id}" "${mock_log}/podman-save.log" || \
  fail "export did not save the resolved immutable image ID"
grep -Fqx "IMAGE_ID=${image_id}" "${archive}.metadata" || \
  fail "export metadata did not retain the resolved image ID"
env PATH="${mock_bin}:${PATH}" MOCK_LOG_DIR="${mock_log}" \
  MOCK_IMAGE_ID="${image_id}" MOCK_GIT_COMMIT="${commit}" \
  MOCK_BASE_IMAGE="${base_image}" \
  "${repo_root}/deploy/load-image.sh" "${archive}"
[[ -s "${mock_log}/podman-load.log" ]] || \
  fail "load wrapper did not import the verified OCI archive"

# Hosted Podman returns a prefixed build iidfile but bare inspect IDs. Scan,
# export, and import must retain the same canonical identity in that case.
bare_scan_dir="${work_dir}/bare-build-scan"
bare_archive="${export_dir}/bare-inspect.oci.tar"
env PATH="${mock_bin}:${PATH}" MOCK_LOG_DIR="${mock_log}" \
  MOCK_IMAGE_ID="${image_id}" MOCK_INSPECT_IMAGE_ID="${image_hex}" \
  MOCK_GIT_COMMIT="${commit}" \
  "${repo_root}/deploy/build-image.sh" \
  --tag local/test:bare --version release --base-image "${base_image}" \
  --scan-dir "${bare_scan_dir}"
grep -Fqx "IMAGE_ID=${image_id}" "${bare_scan_dir}/${safe_name}.scan.txt" || \
  fail "bare inspect ID was not canonicalized in scan evidence"
grep -Fqx "LOCAL_IMAGE_ID=${image_id}" "${bare_scan_dir}/published-image.txt" || \
  fail "bare inspect ID changed the immutable build identity"
env PATH="${mock_bin}:${PATH}" MOCK_LOG_DIR="${mock_log}" \
  MOCK_IMAGE_ID="${image_id}" MOCK_INSPECT_IMAGE_ID="${image_hex}" \
  MOCK_GIT_COMMIT="${commit}" MOCK_BASE_IMAGE="${base_image}" \
  "${repo_root}/deploy/export-image.sh" local/test:bare "${bare_archive}"
grep -Fqx "IMAGE_ID=${image_id}" "${bare_archive}.metadata" || \
  fail "bare inspect ID was not canonicalized in export metadata"
env PATH="${mock_bin}:${PATH}" MOCK_LOG_DIR="${mock_log}" \
  MOCK_IMAGE_ID="${image_id}" MOCK_INSPECT_IMAGE_ID="${image_hex}" \
  MOCK_GIT_COMMIT="${commit}" MOCK_BASE_IMAGE="${base_image}" \
  "${repo_root}/deploy/load-image.sh" "${bare_archive}"

# Normalization must not accept truncated IDs, other algorithms, whitespace,
# repeated prefixes, or arbitrary strings from a failed inspect operation.
invalid_index=0
for invalid_id in '' latest "${image_hex:1}" "sha512:${image_hex}" \
  "sha256:${image_id}" " ${image_hex}" "${image_hex}"$'\n'garbage; do
  invalid_index=$((invalid_index + 1))
  expect_status 65 env PATH="${mock_bin}:${PATH}" MOCK_LOG_DIR="${mock_log}" \
    MOCK_IMAGE_ID="${image_id}" MOCK_INSPECT_IMAGE_ID="${invalid_id}" \
    "${repo_root}/deploy/scan-image.sh" "${image_id}" "${work_dir}/invalid-scan-${invalid_index}"
  expect_status 65 env PATH="${mock_bin}:${PATH}" MOCK_LOG_DIR="${mock_log}" \
    MOCK_IMAGE_ID="${image_id}" MOCK_INSPECT_IMAGE_ID="${invalid_id}" \
    "${repo_root}/deploy/export-image.sh" "${image_id}" "${export_dir}/invalid-${invalid_index}.oci.tar"
  expect_status 65 env PATH="${mock_bin}:${PATH}" MOCK_LOG_DIR="${mock_log}" \
    MOCK_IMAGE_ID="${image_id}" MOCK_INSPECT_IMAGE_ID="${invalid_id}" \
    MOCK_GIT_COMMIT="${commit}" MOCK_BASE_IMAGE="${base_image}" \
    "${repo_root}/deploy/load-image.sh" "${archive}"
done

make_release() {
  local destination="$1"
  mkdir -p -- "${destination}/data/images" "${destination}/attestations"
  printf '{}\n' > "${destination}/data/manifest.json"
  cp -- "${scan_dir}"/* "${destination}/attestations/"
  {
    printf 'LOCAL_IMAGE_ID=%s\n' "${image_id}"
    printf 'GIT_COMMIT=%s\n' "${commit}"
    printf 'BASE_IMAGE=%s\n' "${base_image}"
  } > "${destination}/attestations/published-image.txt"
}

release_ok="${work_dir}/release-ok"
make_release "${release_ok}"
"${repo_root}/deploy/generate-release-manifest.sh" \
  "${release_ok}" release "${image_id}" "${commit}" 1
[[ -f "${release_ok}/MANIFEST.sha256" ]] || fail "valid evidence did not generate a manifest"

release_registry="${work_dir}/release-registry"
make_release "${release_registry}"
registry_ref="ghcr.io/example/islet-explorer@${registry_digest}"
{
  printf 'REGISTRY_REF=%s\n' "${registry_ref}"
  printf 'IMAGE_ID=%s\n' "${image_id}"
  printf 'GIT_COMMIT=%s\n' "${commit}"
  printf 'BASE_IMAGE=%s\n' "${base_image}"
} > "${release_registry}/attestations/published-image.txt"
"${repo_root}/deploy/generate-release-manifest.sh" \
  "${release_registry}" release "${registry_ref}" "${commit}" 1
[[ -f "${release_registry}/MANIFEST.sha256" ]] || \
  fail "valid registry evidence did not generate a manifest"

release_missing="${work_dir}/release-missing"
mkdir -p -- "${release_missing}/data/images"
printf '{}\n' > "${release_missing}/data/manifest.json"
expect_status 66 "${repo_root}/deploy/generate-release-manifest.sh" \
  "${release_missing}" release "${image_id}" "${commit}" 1

release_mismatch="${work_dir}/release-mismatch"
make_release "${release_mismatch}"
wrong_image="sha256:$(printf 'd%.0s' {1..64})"
sed -i "s/^LOCAL_IMAGE_ID=.*/LOCAL_IMAGE_ID=${wrong_image}/" \
  "${release_mismatch}/attestations/published-image.txt"
expect_status 65 "${repo_root}/deploy/generate-release-manifest.sh" \
  "${release_mismatch}" release "${image_id}" "${commit}" 1

release_failed="${work_dir}/release-failed"
make_release "${release_failed}"
sed -i 's/^RESULT=.*/RESULT=fail-high-or-critical/' \
  "${release_failed}/attestations/${safe_name}.scan.txt"
expect_status 65 "${repo_root}/deploy/generate-release-manifest.sh" \
  "${release_failed}" release "${image_id}" "${commit}" 1

release_tampered="${work_dir}/release-tampered"
make_release "${release_tampered}"
printf 'tampered\n' >> "${release_tampered}/attestations/${safe_name}.grype.json"
expect_status 65 "${repo_root}/deploy/generate-release-manifest.sh" \
  "${release_tampered}" release "${image_id}" "${commit}" 1

release_linked="${work_dir}/release-linked-attestation"
make_release "${release_linked}"
ln -s /etc/hosts "${release_linked}/attestations/external-evidence"
expect_status 65 "${repo_root}/deploy/generate-release-manifest.sh" \
  "${release_linked}" release "${image_id}" "${commit}" 1

echo "Image scan and release-evidence contract tests passed."
