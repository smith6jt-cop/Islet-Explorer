#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd -- "${script_dir}/../.." && pwd -P)"
work_dir="$(mktemp -d)"
trap 'rm -rf -- "${work_dir}"' EXIT

fail() { echo "FAIL: $*" >&2; exit 1; }
candidate_image="sha256:$(printf 'a%.0s' {1..64})"
prior_image="sha256:$(printf 'c%.0s' {1..64})"
commit="$(printf 'b%.0s' {1..40})"
mock_bin="${work_dir}/bin"
fixture_deploy="${work_dir}/deploy"
mkdir -p "${mock_bin}" "${fixture_deploy}"
cp "${repo_root}/deploy/activate-release.sh" \
  "${repo_root}/deploy/rollback-release.sh" "${repo_root}/deploy/verify-release.sh" \
  "${repo_root}/deploy/islet-explorer.container.in" "${fixture_deploy}/"

# Exercise real activation, release verification, and rollback with isolated
# user-systemd and image boundaries. Live container checks have separate tests.
cat > "${fixture_deploy}/verify-deployment.sh" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
expected=""
while (($#)); do
  case "$1" in
    --expected-image) expected="$2"; shift 2 ;;
    *) shift ;;
  esac
done
printf '%s\n' "${expected}" >> "${MOCK_CASE}/health-checks"
[[ -n "${expected}" && "$(cat "${MOCK_CASE}/running-image")" == "${expected}" ]] || exit 70
if [[ "${expected}" == "${MOCK_CANDIDATE_IMAGE}" ]]; then
  [[ "${MOCK_SCENARIO}" == success ]] || exit 70
elif [[ "${MOCK_SCENARIO}" == recovery_health_fail ]]; then
  exit 70
fi
MOCK
cat > "${mock_bin}/systemctl" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
unit="${XDG_CONFIG_HOME}/containers/systemd/islet-explorer.container"
image="$(awk -F= '$1 == "Image" {print $2; exit}' "${unit}" 2>/dev/null || true)"
printf '%s\n' "$* ${image}" >> "${MOCK_CASE}/systemctl-calls"
case "${2:-}" in
  daemon-reload)
    if [[ "${MOCK_SCENARIO}" == reload_fail && "${image}" == "${MOCK_CANDIDATE_IMAGE}" ]]; then exit 1; fi
    ;;
  restart)
    if [[ "${image}" == "${MOCK_CANDIDATE_IMAGE}" && "${MOCK_SCENARIO}" == start_fail ]]; then exit 1; fi
    if [[ "${image}" != "${MOCK_CANDIDATE_IMAGE}" && "${MOCK_SCENARIO}" == recovery_restart_fail ]]; then exit 1; fi
    printf '%s\n' "${image}" > "${MOCK_CASE}/running-image"
    ;;
  stop) rm -f "${MOCK_CASE}/running-image" ;;
  *) exit 64 ;;
esac
MOCK
cat > "${mock_bin}/podman" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
case "${1:-} ${2:-}" in
  "image exists") exit 0 ;;
  "image inspect") printf '%s\n' "${MOCK_COMMIT}"; exit 0 ;;
  "container exists") test -f "${MOCK_CASE}/running-image"; exit $? ;;
esac
if [[ "${1:-}" == inspect ]]; then
  if [[ -f "${MOCK_CASE}/running-image" ]]; then echo true; else echo false; fi
  exit 0
fi
if [[ "${1:-}" == run ]]; then
  case " $* " in
    *" /usr/bin/Rscript "*) echo 'Rscript must resolve through PATH' >&2; exit 70 ;;
    *" /usr/bin/id "*) echo 1000 ;;
    *"json.load("*) printf '1\t%s\n' "${MOCK_COMMIT}" ;;
    *) : ;;
  esac
  exit 0
fi
exit 64
MOCK
chmod 0755 "${mock_bin}/systemctl" "${mock_bin}/podman" \
  "${fixture_deploy}/verify-release.sh" "${fixture_deploy}/verify-deployment.sh"

make_release() {
  local root="$1" release="$2" image="$3" count="$4" version="$5" index
  local destination="${root}/releases/${release}"
  mkdir -p "${destination}/data/images"
  {
    printf 'MANIFEST_VERSION=%s\nRELEASE=%s\nIMAGE_REF=%s\nGIT_COMMIT=%s\nDATA_SCHEMA=1\n' \
      "${version}" "${release}" "${image}" "${commit}"
    if [[ "${version}" == 2 ]]; then printf 'EXPECTED_IMAGES=%s\n' "${count}"; fi
  } > "${destination}/release.env"
  printf '{}\n' > "${destination}/data/manifest.json"
  for ((index = 1; index <= count; index++)); do
    printf 'image fixture %s\n' "${index}" > "${destination}/data/images/${index}.ome.tiff"
  done
  (
    cd -- "${destination}"
    find data release.env -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum > MANIFEST.sha256
  )
}

for scenario in success reload_fail start_fail candidate_health_fail recovery_restart_fail recovery_health_fail first_health_fail; do
  case_dir="${work_dir}/${scenario}"
  project="${case_dir}/project"
  config="${case_dir}/config"
  unit="${config}/containers/systemd/islet-explorer.container"
  state="${project}/runtime/deployments"
  mkdir -p "${project}/releases/new/data/images" "$(dirname "${unit}")" "${state}"
  make_release "${project}" new "${candidate_image}" 21 2
  make_release "${project}" old "${prior_image}" 15 1
  if [[ "${scenario}" != first_health_fail ]]; then
    printf '[Container]\nImage=%s\n' "${prior_image}" > "${unit}"
    cp "${unit}" "${case_dir}/original-unit"
    printf 'RELEASE=old\nIMAGE_REF=%s\n' "${prior_image}" > "${state}/active.state"
    cp "${state}/active.state" "${case_dir}/original-state"
    printf '%s\n' "${prior_image}" > "${case_dir}/running-image"
  fi
  set +e
  env PATH="${mock_bin}:${PATH}" XDG_CONFIG_HOME="${config}" \
    MOCK_CASE="${case_dir}" MOCK_SCENARIO="${scenario}" \
    MOCK_CANDIDATE_IMAGE="${candidate_image}" MOCK_COMMIT="${commit}" \
    bash "${fixture_deploy}/activate-release.sh" --project-root "${project}" --release new \
    > "${case_dir}/output" 2>&1
  status=$?
  set -e
  expected=70
  case "${scenario}" in
    success) expected=0 ;;
    recovery_restart_fail|recovery_health_fail) expected=71 ;;
  esac
  if [[ "${status}" != "${expected}" ]]; then
    cat "${case_dir}/output" >&2
    fail "${scenario}: expected exit ${expected}, got ${status}"
  fi
  case "${scenario}" in
    success)
      grep -qx 'RELEASE=new' "${state}/active.state" || fail 'new release was not recorded'
      cmp -s "${state}/previous.state" "${case_dir}/original-state" || fail 'previous state lost'
      grep -Fqx "${candidate_image}" "${case_dir}/health-checks" || fail 'candidate identity was not verified'
      grep -Fq 'verified: 21 OME-TIFFs' "${case_dir}/output" || fail 'activation did not use the expanded release count'
      env PATH="${mock_bin}:${PATH}" XDG_CONFIG_HOME="${config}" \
        MOCK_CASE="${case_dir}" MOCK_SCENARIO=success \
        MOCK_CANDIDATE_IMAGE="${prior_image}" MOCK_COMMIT="${commit}" \
        bash "${fixture_deploy}/rollback-release.sh" --project-root "${project}" \
        > "${case_dir}/rollback-output" 2>&1 || {
          cat "${case_dir}/rollback-output" >&2
          fail 'rollback to the legacy 15-image release failed'
        }
      grep -Fq 'verified: 15 OME-TIFFs' "${case_dir}/rollback-output" || fail 'rollback did not use the original release count'
      grep -qx 'RELEASE=old' "${state}/active.state" || fail 'rollback did not restore the original release'
      grep -qx 'RELEASE=new' "${state}/previous.state" || fail 'rollback lost the expanded release state'
      ;;
    first_health_fail)
      [[ ! -e "${unit}" && ! -e "${state}/active.state" && ! -e "${case_dir}/running-image" ]] || \
        fail 'failed first activation left a candidate service or active state'
      ;;
    *)
      cmp -s "${unit}" "${case_dir}/original-unit" || fail "${scenario}: old unit was not restored"
      cmp -s "${state}/active.state" "${case_dir}/original-state" || fail "${scenario}: active state changed"
      [[ ! -e "${state}/previous.state" ]] || fail "${scenario}: prior state was overwritten"
      if [[ "${expected}" == 70 ]]; then
        grep -Fqx "${prior_image}" "${case_dir}/health-checks" || fail "${scenario}: recovery identity was not verified"
        grep -Fq 'restored and verified healthy' "${case_dir}/output" || fail "${scenario}: recovery outcome missing"
      else
        grep -Fq 'activation and recovery failed' "${case_dir}/output" || fail "${scenario}: recovery failure was hidden"
      fi
      ;;
  esac
done

echo 'Deployment activation and verified recovery tests passed.'
