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
# The unit acted on is the instance named by the *.service argument
# (activate-release.sh --instance); the production service is the default.
service="islet-explorer"
for arg in "$@"; do [[ "${arg}" == *.service ]] && service="${arg%.service}"; done
unit="${XDG_CONFIG_HOME}/containers/systemd/${service}.container"
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
  # The small-read trees that --local-data copies to the VM's local disk.
  local tree
  for tree in geometry cells spatial; do
    mkdir -p "${destination}/data/${tree}/nested"
    printf '%s fixture %s a\n' "${tree}" "${release}" > "${destination}/data/${tree}/a.bin"
    printf '%s fixture %s b\n' "${tree}" "${release}" > "${destination}/data/${tree}/nested/b.bin"
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
      grep -qx 'ContainerName=islet-explorer' "${unit}" || fail 'production unit lost its container name'
      grep -qx 'PublishPort=127.0.0.1:3838:3838' "${unit}" || fail 'production unit lost its loopback port'
      grep -q -- '--memory=9g ' "${unit}" || fail 'production unit lost its default memory cap'
      grep -qx "Volume=${project}/runtime/logs:/var/log/shiny-server:rw" "${unit}" || fail 'production logs moved'
      grep -qx 'PORT=3838' "${state}/active.state" || fail 'port not recorded'
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

# A second instance (the staging service beside production, 2026-10-23) gets
# its own unit, container name, host port, memory cap, runtime tree and state
# files, leaves the production unit and state untouched, refuses a port another
# unit already publishes, and rolls back on its own state.
case_dir="${work_dir}/second_instance"
project="${case_dir}/project"
config="${case_dir}/config"
unit_dir="${config}/containers/systemd"
state="${project}/runtime/deployments"
mkdir -p "${unit_dir}" "${state}"
make_release "${project}" new "${candidate_image}" 21 2
make_release "${project}" old "${prior_image}" 15 1
printf '[Container]\nImage=%s\nPublishPort=127.0.0.1:3838:3838\n' "${prior_image}" > "${unit_dir}/islet-explorer.container"
cp "${unit_dir}/islet-explorer.container" "${case_dir}/original-unit"
printf 'RELEASE=old\nIMAGE_REF=%s\n' "${prior_image}" > "${state}/active.state"
cp "${state}/active.state" "${case_dir}/original-state"
run_dev() {
  env PATH="${mock_bin}:${PATH}" XDG_CONFIG_HOME="${config}" \
    MOCK_CASE="${case_dir}" MOCK_SCENARIO=success \
    MOCK_CANDIDATE_IMAGE="${candidate_image}" MOCK_COMMIT="${commit}" "$@"
}
run_dev bash "${fixture_deploy}/activate-release.sh" --project-root "${project}" --release new \
  --instance islet-explorer-dev --port 3839 --memory 4g > "${case_dir}/output" 2>&1 || {
    cat "${case_dir}/output" >&2
    fail 'second instance activation failed'
  }
dev_unit="${unit_dir}/islet-explorer-dev.container"
[[ -f "${dev_unit}" ]] || fail 'second instance unit was not rendered'
grep -qx 'ContainerName=islet-explorer-dev' "${dev_unit}" || fail 'second instance container name'
grep -qx 'PublishPort=127.0.0.1:3839:3838' "${dev_unit}" || fail 'second instance port'
grep -q -- '--memory=4g ' "${dev_unit}" || fail 'second instance memory cap'
grep -qx "Volume=${project}/runtime/islet-explorer-dev/logs:/var/log/shiny-server:rw" "${dev_unit}" || fail 'second instance log directory'
grep -qx "Volume=${project}/runtime/islet-explorer-dev/tmp:/tmp/islet-explorer:rw" "${dev_unit}" || fail 'second instance tmp directory'
[[ -d "${project}/runtime/islet-explorer-dev/logs" && -d "${project}/runtime/islet-explorer-dev/tmp" ]] || fail 'second instance runtime tree'
cmp -s "${unit_dir}/islet-explorer.container" "${case_dir}/original-unit" || fail 'production unit changed by the second instance'
cmp -s "${state}/active.state" "${case_dir}/original-state" || fail 'production state changed by the second instance'
[[ ! -e "${state}/previous.state" ]] || fail 'production previous state written by the second instance'
grep -qx 'RELEASE=new' "${state}/islet-explorer-dev.active.state" || fail 'second instance state not recorded'
grep -qx 'INSTANCE=islet-explorer-dev' "${state}/islet-explorer-dev.active.state" || fail 'second instance name not recorded'
grep -qx 'PORT=3839' "${state}/islet-explorer-dev.active.state" || fail 'second instance port not recorded'
grep -qx 'MEMORY=4g' "${state}/islet-explorer-dev.active.state" || fail 'second instance memory not recorded'
[[ ! -e "${state}/islet-explorer-dev.previous.state" ]] || fail 'first second-instance activation wrote a previous state'
grep -q 'restart islet-explorer-dev.service' "${case_dir}/systemctl-calls" || fail 'second instance service was not restarted'
grep -q 'restart islet-explorer.service' "${case_dir}/systemctl-calls" && fail 'production service was restarted by the second instance'
cp "${dev_unit}" "${case_dir}/dev-unit-before-collision"
set +e
run_dev bash "${fixture_deploy}/activate-release.sh" --project-root "${project}" --release new \
  --instance islet-explorer-dev --port 3838 --memory 4g > "${case_dir}/collision-output" 2>&1
status=$?
set -e
[[ "${status}" == 65 ]] || { cat "${case_dir}/collision-output" >&2; fail "port collision: expected exit 65, got ${status}"; }
grep -Fq 'already published by islet-explorer.container' "${case_dir}/collision-output" || fail 'port collision not reported'
cmp -s "${dev_unit}" "${case_dir}/dev-unit-before-collision" || fail 'port collision changed the unit'
run_dev bash "${fixture_deploy}/activate-release.sh" --project-root "${project}" --release old \
  --instance islet-explorer-dev --port 3839 --memory 4g > "${case_dir}/output-old" 2>&1 || {
    cat "${case_dir}/output-old" >&2
    fail 'second instance re-activation failed'
  }
grep -qx 'RELEASE=old' "${state}/islet-explorer-dev.active.state" || fail 'second instance re-activation not recorded'
grep -qx 'RELEASE=new' "${state}/islet-explorer-dev.previous.state" || fail 'second instance previous state missing'
env PATH="${mock_bin}:${PATH}" XDG_CONFIG_HOME="${config}" \
  MOCK_CASE="${case_dir}" MOCK_SCENARIO=success \
  MOCK_CANDIDATE_IMAGE="${candidate_image}" MOCK_COMMIT="${commit}" \
  bash "${fixture_deploy}/rollback-release.sh" --project-root "${project}" --instance islet-explorer-dev \
  > "${case_dir}/rollback-output" 2>&1 || {
    cat "${case_dir}/rollback-output" >&2
    fail 'second instance rollback failed'
  }
grep -Fq 'Rolling back islet-explorer-dev to release new on 127.0.0.1:3839' "${case_dir}/rollback-output" || fail 'second instance rollback target'
grep -qx 'RELEASE=new' "${state}/islet-explorer-dev.active.state" || fail 'second instance rollback not recorded'
grep -qx 'PublishPort=127.0.0.1:3839:3838' "${dev_unit}" || fail 'second instance rollback lost its port'
grep -q -- '--memory=4g ' "${dev_unit}" || fail 'second instance rollback lost its memory cap'
cmp -s "${state}/active.state" "${case_dir}/original-state" || fail 'production state changed by the second instance rollback'

# UFRC's proxy serves the three image aliases itself from a DocumentRoot on the
# shared filesystem (vhost of 2026-10-23). With --docroot, activation points
# DOCROOT/islet-explorer/{local_images,images,img_data} at the verified
# release's images, a failed activation points them back at the restored
# release, an alias that is a real directory is refused before any service is
# touched, and rollback replays the recorded DocumentRoot.
case_dir="${work_dir}/docroot"
project="${case_dir}/project"
config="${case_dir}/config"
unit_dir="${config}/containers/systemd"
state="${project}/runtime/deployments"
docroot="${case_dir}/html"
mkdir -p "${unit_dir}" "${state}" "${docroot}"
make_release "${project}" new "${candidate_image}" 21 2
make_release "${project}" old "${prior_image}" 15 1
printf '[Container]\nImage=%s\nPublishPort=127.0.0.1:3838:3838\nVolume=%s/releases/old/data/images:/srv/islet-explorer/app/shiny_app/www/local_images:ro\n' \
  "${prior_image}" "${project}" > "${unit_dir}/islet-explorer.container"
printf 'RELEASE=old\nIMAGE_REF=%s\n' "${prior_image}" > "${state}/active.state"
printf '%s\n' "${prior_image}" > "${case_dir}/running-image"
run_docroot() {
  local scenario="$1" candidate="$2"
  shift 2
  env PATH="${mock_bin}:${PATH}" XDG_CONFIG_HOME="${config}" \
    MOCK_CASE="${case_dir}" MOCK_SCENARIO="${scenario}" \
    MOCK_CANDIDATE_IMAGE="${candidate}" MOCK_COMMIT="${commit}" "$@"
}
aliases_point_at() {
  local target="$1" alias
  for alias in local_images images img_data; do
    [[ -L "${docroot}/islet-explorer/${alias}" && "$(readlink "${docroot}/islet-explorer/${alias}")" == "${target}" ]] || return 1
  done
  [[ -z "$(find "${docroot}/islet-explorer" -name '.*.tmp' -print -quit)" ]]
}
run_docroot success "${candidate_image}" bash "${fixture_deploy}/activate-release.sh" --project-root "${project}" \
  --release new --docroot "${docroot}" > "${case_dir}/output" 2>&1 || {
    cat "${case_dir}/output" >&2
    fail 'activation with a DocumentRoot failed'
  }
aliases_point_at "${project}/releases/new/data/images" || fail 'image aliases do not point at the activated release'
grep -qx "DOCROOT=${docroot}" "${state}/active.state" || fail 'DocumentRoot not recorded'
grep -Fq "Proxy image aliases under ${docroot}/islet-explorer now point at ${project}/releases/new/data/images" "${case_dir}/output" || \
  fail 'alias update not reported'
cp "${unit_dir}/islet-explorer.container" "${case_dir}/unit-before-refusal"
cp "${case_dir}/systemctl-calls" "${case_dir}/calls-before-refusal"
rm "${docroot}/islet-explorer/images"
mkdir "${docroot}/islet-explorer/images"
set +e
run_docroot success "${candidate_image}" bash "${fixture_deploy}/activate-release.sh" --project-root "${project}" \
  --release old --docroot "${docroot}" > "${case_dir}/refusal-output" 2>&1
status=$?
set -e
[[ "${status}" == 65 ]] || { cat "${case_dir}/refusal-output" >&2; fail "real directory in an alias's place: expected exit 65, got ${status}"; }
grep -Fq "${docroot}/islet-explorer/images exists and is not a symlink" "${case_dir}/refusal-output" || fail 'real directory not reported'
[[ -d "${docroot}/islet-explorer/images" && ! -L "${docroot}/islet-explorer/images" ]] || fail 'the real directory was removed'
cmp -s "${unit_dir}/islet-explorer.container" "${case_dir}/unit-before-refusal" || fail 'refusal changed the unit'
cmp -s "${case_dir}/systemctl-calls" "${case_dir}/calls-before-refusal" || fail 'refusal touched the service'
rmdir "${docroot}/islet-explorer/images"
ln -s "${project}/releases/new/data/images" "${docroot}/islet-explorer/images"
# A failed activation of the old release (its image is the failing candidate
# here) restores the new unit; aliases that had drifted follow the restored release.
for alias in local_images images img_data; do ln -sfn "${project}/releases/old/data/images" "${docroot}/islet-explorer/${alias}"; done
set +e
run_docroot candidate_health_fail "${prior_image}" bash "${fixture_deploy}/activate-release.sh" --project-root "${project}" \
  --release old --docroot "${docroot}" > "${case_dir}/failed-output" 2>&1
status=$?
set -e
[[ "${status}" == 70 ]] || { cat "${case_dir}/failed-output" >&2; fail "failed activation with a DocumentRoot: expected exit 70, got ${status}"; }
grep -Fq 'restored and verified healthy' "${case_dir}/failed-output" || fail 'recovery outcome missing with a DocumentRoot'
aliases_point_at "${project}/releases/new/data/images" || fail 'image aliases do not follow the restored release'
grep -qx 'RELEASE=new' "${state}/active.state" || fail 'failed activation changed the active state'
run_docroot success "${candidate_image}" bash "${fixture_deploy}/activate-release.sh" --project-root "${project}" \
  --release old --docroot "${docroot}" > "${case_dir}/output-old" 2>&1 || {
    cat "${case_dir}/output-old" >&2
    fail 'activation of the old release with a DocumentRoot failed'
  }
aliases_point_at "${project}/releases/old/data/images" || fail 'image aliases do not point at the re-activated old release'
grep -qx "DOCROOT=${docroot}" "${state}/previous.state" || fail 'previous state lost the DocumentRoot'
run_docroot success "${candidate_image}" bash "${fixture_deploy}/rollback-release.sh" --project-root "${project}" \
  > "${case_dir}/rollback-output" 2>&1 || {
    cat "${case_dir}/rollback-output" >&2
    fail 'rollback with a DocumentRoot failed'
  }
grep -qx 'RELEASE=new' "${state}/active.state" || fail 'rollback did not restore the new release'
grep -qx "DOCROOT=${docroot}" "${state}/active.state" || fail 'rollback lost the DocumentRoot'
aliases_point_at "${project}/releases/new/data/images" || fail 'image aliases do not follow the rollback'

# --local-data copies the release's geometry, cells and spatial trees to the
# VM's local disk, checks them against MANIFEST.sha256 and mounts them over the
# container's data paths (2026-10-24). A copy that does not match the manifest
# stops the activation before any service is touched; copies of earlier
# releases stay in place; rollback replays the recorded directory.
case_dir="${work_dir}/local_data"
project="${case_dir}/project"
config="${case_dir}/config"
unit_dir="${config}/containers/systemd"
state="${project}/runtime/deployments"
local_data="${case_dir}/local"
unit="${unit_dir}/islet-explorer.container"
mkdir -p "${unit_dir}" "${state}" "${local_data}"
make_release "${project}" new "${candidate_image}" 21 2
make_release "${project}" old "${prior_image}" 15 1
run_local() {
  local scenario="$1" candidate="$2"
  shift 2
  env PATH="${mock_bin}:${PATH}" XDG_CONFIG_HOME="${config}" \
    MOCK_CASE="${case_dir}" MOCK_SCENARIO="${scenario}" \
    MOCK_CANDIDATE_IMAGE="${candidate}" MOCK_COMMIT="${commit}" "$@"
}
local_copies_match() {
  local release="$1" tree
  for tree in geometry cells spatial; do
    cmp -s "${project}/releases/${release}/data/${tree}/a.bin" "${local_data}/${release}/data/${tree}/a.bin" || return 1
    cmp -s "${project}/releases/${release}/data/${tree}/nested/b.bin" "${local_data}/${release}/data/${tree}/nested/b.bin" || return 1
  done
}
run_local success "${candidate_image}" bash "${fixture_deploy}/activate-release.sh" --project-root "${project}" \
  --release new --local-data "${local_data}" > "${case_dir}/output" 2>&1 || {
    cat "${case_dir}/output" >&2
    fail 'activation with local data copies failed'
  }
local_copies_match new || fail 'local data copies do not match the release'
for tree in geometry cells spatial; do
  grep -qx "Volume=${local_data}/new/data/${tree}:/srv/islet-explorer/data/${tree}:ro" "${unit}" || fail "unit does not mount the local ${tree} copy"
done
grep -qx "Volume=${project}/releases/new/data:/srv/islet-explorer/data:ro" "${unit}" || fail 'unit lost the shared data mount'
grep -q "^RequiresMountsFor=.* ${local_data}/new/data$" "${unit}" || fail 'unit does not require the local data mount'
grep -q '@@LOCAL_DATA' "${unit}" && fail 'unit kept a local-data placeholder'
grep -qx "LOCAL_DATA=${local_data}" "${state}/active.state" || fail 'local data directory not recorded'
grep -Fq "Local copies of data/{geometry,cells,spatial} under ${local_data}/new/data verified against MANIFEST.sha256" "${case_dir}/output" || \
  fail 'local data verification not reported'
# A stale local copy with a preserved timestamp survives rsync's quick check; the
# manifest check must catch it before the service is touched.
cp "${unit}" "${case_dir}/unit-before-mismatch"
cp "${case_dir}/systemctl-calls" "${case_dir}/calls-before-mismatch"
printf 'cells fixture new X\n' > "${local_data}/new/data/cells/a.bin"   # same length as the fixture
touch -r "${project}/releases/new/data/cells/a.bin" "${local_data}/new/data/cells/a.bin"
set +e
run_local success "${candidate_image}" bash "${fixture_deploy}/activate-release.sh" --project-root "${project}" \
  --release new --local-data "${local_data}" > "${case_dir}/mismatch-output" 2>&1
status=$?
set -e
[[ "${status}" == 65 ]] || { cat "${case_dir}/mismatch-output" >&2; fail "mismatched local copy: expected exit 65, got ${status}"; }
grep -Fq "do not match MANIFEST.sha256" "${case_dir}/mismatch-output" || fail 'mismatched local copy not reported'
cmp -s "${unit}" "${case_dir}/unit-before-mismatch" || fail 'mismatched local copy changed the unit'
cmp -s "${case_dir}/systemctl-calls" "${case_dir}/calls-before-mismatch" || fail 'mismatched local copy touched the service'
rm "${local_data}/new/data/cells/a.bin"
# Activating the old release adds its own copy and leaves the new release's in place.
run_local success "${candidate_image}" bash "${fixture_deploy}/activate-release.sh" --project-root "${project}" \
  --release old --local-data "${local_data}" > "${case_dir}/output-old" 2>&1 || {
    cat "${case_dir}/output-old" >&2
    fail 'activation of the old release with local data copies failed'
  }
local_copies_match old || fail 'old release local copies do not match'
[[ -f "${local_data}/new/data/geometry/a.bin" ]] || fail 'activation removed the local copy of another release'
grep -qx "Volume=${local_data}/old/data/geometry:/srv/islet-explorer/data/geometry:ro" "${unit}" || fail 'old release unit does not mount its local copy'
grep -qx "LOCAL_DATA=${local_data}" "${state}/previous.state" || fail 'previous state lost the local data directory'
run_local success "${candidate_image}" bash "${fixture_deploy}/rollback-release.sh" --project-root "${project}" \
  > "${case_dir}/rollback-output" 2>&1 || {
    cat "${case_dir}/rollback-output" >&2
    fail 'rollback with local data copies failed'
  }
grep -qx 'RELEASE=new' "${state}/active.state" || fail 'rollback did not restore the new release'
grep -qx "LOCAL_DATA=${local_data}" "${state}/active.state" || fail 'rollback lost the local data directory'
grep -qx "Volume=${local_data}/new/data/spatial:/srv/islet-explorer/data/spatial:ro" "${unit}" || fail 'rollback unit does not mount the local copy'
# Without --local-data the unit carries no local mounts and the placeholder is gone.
run_local success "${candidate_image}" bash "${fixture_deploy}/activate-release.sh" --project-root "${project}" \
  --release old > "${case_dir}/output-plain" 2>&1 || {
    cat "${case_dir}/output-plain" >&2
    fail 'activation without local data failed'
  }
grep -q "/srv/islet-explorer/data/geometry:ro" "${unit}" && fail 'unit without --local-data mounts a local copy'
grep -q '@@LOCAL_DATA' "${unit}" && fail 'unit without --local-data kept a placeholder'
grep -qx 'LOCAL_DATA=' "${state}/active.state" || fail 'empty local data directory not recorded'

echo 'Deployment activation and verified recovery tests passed.'
