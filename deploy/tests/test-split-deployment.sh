#!/usr/bin/env bash
# Integration coverage for two app releases sharing one runtime and dataset.
# Artifact validation itself is covered by test_release_artifacts.py; only that
# resolver, Podman, systemd and the external health boundary are replaced here.
set -euo pipefail
repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)"
work_dir="$(mktemp -d)"
trap 'rm -rf -- "${work_dir}"' EXIT
fixture_deploy="${work_dir}/deploy"
mock_bin="${work_dir}/bin"
project="${work_dir}/project"
config="${work_dir}/config"
state="${project}/runtime/deployments"
unit_dir="${config}/containers/systemd"
unit="${unit_dir}/islet-explorer-dev.container"
docroot="${project}/dev/html"
local_data="${work_dir}/local-data"
mkdir -p "${fixture_deploy}" "${mock_bin}" "${state}" "${unit_dir}" "${docroot}"
cp "${repo_root}/deploy/activate-release.sh" "${repo_root}/deploy/rollback-release.sh" \
   "${repo_root}/deploy/release-layout.sh" "${repo_root}/deploy/islet-explorer.container.in" \
   "${repo_root}/deploy/prepare-local-data.py" \
   "${fixture_deploy}/"
image="sha256:$(printf 'a%.0s' {1..64})"
runtime_commit="$(printf 'b%.0s' {1..40})"
data_commit="$(printf 'c%.0s' {1..40})"
inputs="$(printf 'd%.0s' {1..64})"
base_image="example@sha256:$(printf 'e%.0s' {1..64})"
data_release="${project}/releases/data-original"
mkdir -p "${data_release}/data/images"
printf '{"schema_version":5,"git_commit":"%s"}\n' "${data_commit}" > "${data_release}/data/manifest.json"
printf 'microscopy fixture must not be hashed by activation\n' > "${data_release}/data/images/test.ome.tiff"
for tree in geometry cells spatial; do
  mkdir -p "${data_release}/data/${tree}"
  printf 'unchanged %s\n' "${tree}" > "${data_release}/data/${tree}/table.bin"
done
(cd "${data_release}" && find data -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum > MANIFEST.sha256)

for release in first second broken; do
  mkdir -p "${project}/apps/${release}/app/shiny_app/www/local_images" \
    "${project}/apps/${release}/scripts" "${project}/releases/${release}"
  printf '{}\n' > "${project}/releases/${release}/deployment.json"
  printf '# %s app\n' "${release}" > "${project}/apps/${release}/app/shiny_app/app.R"
done

cat > "${fixture_deploy}/release-artifacts.py" <<'PY'
import hashlib, os, pathlib, sys
if sys.argv[1] == 'verify-app':
    app = pathlib.Path(sys.argv[3])
    assert sys.argv[4] == '--expected-manifest'
    assert sys.argv[5] == hashlib.sha256(app.name.encode()).hexdigest()
    assert not (app / 'tampered').exists()
    sys.exit(0)
assert sys.argv[1:3] == ['resolve', '--release-dir']
release_dir = pathlib.Path(sys.argv[3])
root = release_dir.parent.parent
release = release_dir.name
fields = {
    'IMAGE_REF': os.environ['MOCK_IMAGE'],
    'RUNTIME_COMMIT': os.environ['MOCK_RUNTIME_COMMIT'],
    'APP_COMMIT': hashlib.sha1(release.encode()).hexdigest(),
    'APP_DIR': root / 'apps' / release / 'app',
    'APP_SCRIPTS_DIR': root / 'apps' / release / 'scripts',
    'DATA_DIR': root / 'releases/data-original/data',
    'DATA_SCHEMA': '5',
    'DATA_COMMIT': os.environ['MOCK_DATA_COMMIT'],
    'DATA_ID': 'shared-data',
    'DATA_RELEASE_DIR': root / 'releases/data-original',
    'DATA_MANIFEST': root / 'releases/data-original/MANIFEST.sha256',
    'APP_MANIFEST_SHA256': hashlib.sha256(release.encode()).hexdigest(),
    'RUNTIME_INPUT_SHA256': os.environ['MOCK_INPUTS'],
    'RUNTIME_BASE_IMAGE': os.environ['MOCK_BASE_IMAGE'],
    'DESCRIPTOR_SHA256': hashlib.sha256(('descriptor-' + release).encode()).hexdigest(),
}
for key, value in fields.items():
    print(f'{key}={value}')
PY

cat > "${mock_bin}/podman" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" >> "${MOCK_WORK}/podman-calls"
case "${1:-} ${2:-}" in
  'image exists') exit 0 ;;
  'image inspect')
    case "$*" in
      *io.islet-explorer.runtime-layout*) echo "${MOCK_LAYOUT:-split-v1}" ;;
      *io.islet-explorer.runtime-inputs-sha256*) echo "${MOCK_INPUTS}" ;;
      *org.opencontainers.image.revision*) echo "${MOCK_RUNTIME_COMMIT}" ;;
      *org.opencontainers.image.base.name*) echo "${MOCK_ACTUAL_BASE_IMAGE:-${MOCK_BASE_IMAGE}}" ;;
      *) exit 64 ;;
    esac
    exit 0 ;;
  'container exists') test -f "${MOCK_WORK}/running-unit"; exit $? ;;
esac
if [[ "${1:-}" == inspect ]]; then
  test -f "${MOCK_WORK}/running-unit" && echo true || echo false
  exit 0
fi
if [[ "${1:-}" == run ]]; then
  case " $* " in
    *' /usr/bin/id '*) echo 1000 ;;
    *'json.load('*) printf '5\t%s\n' "${MOCK_DATA_COMMIT}" ;;
    *'verify_release_data.R'*)
      [[ "$*" == *'/srv/islet-explorer/app:ro'* && "$*" == *'/srv/islet-explorer/scripts:ro'* ]] ;;
    *) : ;;
  esac
  exit 0
fi
exit 64
MOCK

cat > "${mock_bin}/systemctl" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" >> "${MOCK_WORK}/systemctl-calls"
service='islet-explorer-dev'
for arg in "$@"; do [[ "${arg}" != *.service ]] || service="${arg%.service}"; done
[[ "${service}" == islet-explorer-dev ]] || exit 77
unit="${XDG_CONFIG_HOME}/containers/systemd/${service}.container"
case "${2:-}" in
  is-active) test -f "${MOCK_WORK}/running-unit" ;;
  daemon-reload) : ;;
  restart) cp "${unit}" "${MOCK_WORK}/running-unit" ;;
  stop) rm -f "${MOCK_WORK}/running-unit" ;;
  *) exit 64 ;;
esac
MOCK

real_rsync="$(command -v rsync)"
cat > "${mock_bin}/rsync" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" >> "${MOCK_WORK}/rsync-calls"
exec "${MOCK_REAL_RSYNC}" "$@"
MOCK

cat > "${fixture_deploy}/verify-deployment.sh" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
image=''; app=''; data=''; app_hash=''
while (($#)); do
  case "$1" in
    --expected-image) image="$2"; shift 2 ;;
    --expected-app) app="$2"; shift 2 ;;
    --expected-data) data="$2"; shift 2 ;;
    --expected-app-sha256) app_hash="$2"; shift 2 ;;
    *) shift ;;
  esac
done
[[ -n "${app}" && -n "${data}" && -n "${app_hash}" ]] || exit 70
running="${MOCK_WORK}/running-unit"
grep -Fxq "Image=${image}" "${running}"
grep -Fxq "Volume=${app}:/srv/islet-explorer/app:ro" "${running}"
grep -Fxq "Volume=${data}:/srv/islet-explorer/data:ro" "${running}"
grep -Fxq "Label=io.islet-explorer.app-manifest-sha256=${app_hash}" "${running}"
printf '%s %s\n' "${app}" "${app_hash}" >> "${MOCK_WORK}/health-checks"
[[ "${app}" != */broken/app ]]
MOCK
chmod +x "${mock_bin}/podman" "${mock_bin}/systemctl" "${mock_bin}/rsync" \
  "${fixture_deploy}/verify-deployment.sh" "${fixture_deploy}/activate-release.sh"

# Sentinel production unit/state must remain byte-identical throughout.
printf '[Container]\nContainerName=islet-explorer\nPublishPort=127.0.0.1:3838:3838\n' > "${unit_dir}/islet-explorer.container"
printf 'RELEASE=production-untouched\n' > "${state}/active.state"
cp "${unit_dir}/islet-explorer.container" "${work_dir}/production-unit"
cp "${state}/active.state" "${work_dir}/production-state"
run() {
  env PATH="${mock_bin}:${PATH}" XDG_CONFIG_HOME="${config}" MOCK_WORK="${work_dir}" \
    MOCK_REAL_RSYNC="${real_rsync}" \
    MOCK_IMAGE="${image}" MOCK_RUNTIME_COMMIT="${runtime_commit}" \
    MOCK_DATA_COMMIT="${data_commit}" MOCK_INPUTS="${inputs}" MOCK_BASE_IMAGE="${base_image}" "$@"
}
activate() {
  run bash "${fixture_deploy}/activate-release.sh" --project-root "${project}" \
    --release "$1" --instance islet-explorer-dev --port 3839 --memory 4g --docroot "${docroot}" \
    --local-data "${local_data}"
}
fail() { echo "FAIL: $*" >&2; exit 1; }
app_hash() { printf '%s' "$1" | sha256sum | awk '{print $1}'; }
assert_release() {
  local release="$1"
  grep -Fxq "RELEASE=${release}" "${state}/islet-explorer-dev.active.state" || fail 'wrong active release'
  grep -Fxq "Image=${image}" "${unit}" || fail 'runtime changed'
  grep -Fxq "Volume=${project}/apps/${release}/app:/srv/islet-explorer/app:ro" "${unit}" || fail 'wrong app mount'
  grep -Fxq "Volume=${project}/apps/${release}/scripts:/srv/islet-explorer/scripts:ro" "${unit}" || fail 'wrong scripts mount'
  grep -Fxq "Label=io.islet-explorer.app-manifest-sha256=$(app_hash "${release}")" "${unit}" || fail 'wrong app label'
  grep -Fxq "DATA_ID=shared-data" "${state}/islet-explorer-dev.active.state" || fail 'wrong data identity'
  for tree in geometry cells spatial; do
    grep -Fxq "Volume=${local_data}/datasets/shared-data/data/${tree}:/srv/islet-explorer/data/${tree}:ro" "${unit}" || fail 'local cache not keyed by data ID'
  done
  grep -Fxq "Volume=${data_release}/data/images:/srv/islet-explorer/app/shiny_app/www/local_images:ro" "${unit}" || fail 'image mount changed'
  cmp "${unit}" "${work_dir}/running-unit" || fail 'running and selected units differ'
  for alias in local_images images img_data; do
    [[ "$(readlink "${docroot}/islet-explorer/${alias}")" == "${data_release}/data/images" ]] || fail 'wrong image alias'
  done
  ! grep -q '@@' "${unit}" || fail 'unresolved template placeholder'
}

activate first > "${work_dir}/first.log" 2>&1 || { cat "${work_dir}/first.log"; fail 'first activation'; }
assert_release first
[[ "$(wc -l < "${work_dir}/rsync-calls")" == 3 ]] || fail 'initial cache should copy each tree once'
cache_file="${local_data}/datasets/shared-data/data/spatial/table.bin"
cache_inode="$(stat -c '%i' "${cache_file}")"
alias_inode="$(stat -c '%i' "${docroot}/islet-explorer/local_images")"
activate second > "${work_dir}/second.log" 2>&1 || { cat "${work_dir}/second.log"; fail 'second activation'; }
assert_release second
[[ "$(wc -l < "${work_dir}/rsync-calls")" == 3 ]] || fail 'second app recopied unchanged data'
[[ "$(stat -c '%i' "${cache_file}")" == "${cache_inode}" ]] || fail 'second app replaced cache'
[[ "$(stat -c '%i' "${docroot}/islet-explorer/local_images")" == "${alias_inode}" ]] || fail 'unchanged alias replaced'
grep -Fxq 'RELEASE=first' "${state}/islet-explorer-dev.previous.state" || fail 'previous descriptor lost'
run bash "${fixture_deploy}/rollback-release.sh" --project-root "${project}" --instance islet-explorer-dev \
  > "${work_dir}/rollback.log" 2>&1 || { cat "${work_dir}/rollback.log"; fail 'same-runtime rollback'; }
assert_release first
[[ "$(wc -l < "${work_dir}/rsync-calls")" == 3 ]] || fail 'rollback recopied unchanged data'
cp "${state}/islet-explorer-dev.active.state" "${work_dir}/active-before-failure"
cp "${state}/islet-explorer-dev.previous.state" "${work_dir}/previous-before-failure"
set +e
activate broken > "${work_dir}/failure.log" 2>&1
status=$?
set -e
[[ "${status}" == 70 ]] || { cat "${work_dir}/failure.log"; fail "recovery expected exit70, got${status}"; }
assert_release first
cmp "${state}/islet-explorer-dev.active.state" "${work_dir}/active-before-failure" || fail 'failure changed active state'
cmp "${state}/islet-explorer-dev.previous.state" "${work_dir}/previous-before-failure" || fail 'failure changed previous state'
grep -Fq "${project}/apps/first/app $(app_hash first)" "${work_dir}/health-checks" || fail 'restored app not checked'

# An image with the wrong runtime layout must fail before service replacement.
cp "${unit}" "${work_dir}/before-layout-failure"
set +e
MOCK_LAYOUT=legacy activate second > "${work_dir}/layout-failure.log" 2>&1
status=$?
set -e
[[ "${status}" == 65 ]] || { cat "${work_dir}/layout-failure.log"; fail 'invalid runtime layout accepted'; }
cmp "${unit}" "${work_dir}/before-layout-failure" || fail 'invalid layout changed unit'

# Runtime receipt provenance must agree with the actual immutable image label.
cp "${work_dir}/systemctl-calls" "${work_dir}/before-base-failure"
set +e
MOCK_ACTUAL_BASE_IMAGE="example@sha256:$(printf 'f%.0s' {1..64})" activate second > "${work_dir}/base-failure.log" 2>&1
status=$?
set -e
[[ "${status}" == 65 ]] || { cat "${work_dir}/base-failure.log"; fail 'invalid runtime base image accepted'; }
cmp "${unit}" "${work_dir}/before-layout-failure" || fail 'invalid base image changed unit'
cmp "${work_dir}/systemctl-calls" "${work_dir}/before-base-failure" || fail 'invalid base image touched service'

# A published cache is never repaired in place: mutation blocks before restart.
cp "${work_dir}/systemctl-calls" "${work_dir}/before-cache-failure"
printf 'corrupted local cache\n' > "${cache_file}"
set +e
activate second > "${work_dir}/cache-failure.log" 2>&1
status=$?
set -e
[[ "${status}" != 0 ]] || fail 'cache mutation accepted'
grep -Fq 'Published local cache changed' "${work_dir}/cache-failure.log" || { cat "${work_dir}/cache-failure.log"; fail 'cache mutation not explained'; }
cmp "${work_dir}/systemctl-calls" "${work_dir}/before-cache-failure" || fail 'cache mutation touched service'
cmp "${unit}" "${work_dir}/before-layout-failure" || fail 'cache mutation changed unit'
cmp "${unit_dir}/islet-explorer.container" "${work_dir}/production-unit" || fail 'production unit changed'
cmp "${state}/active.state" "${work_dir}/production-state" || fail 'production state changed'
[[ ! -e "${state}/previous.state" ]] || fail 'production previous state created'
! grep -Eq ' (build|pull|push) |restart islet-explorer.service' "${work_dir}/podman-calls" "${work_dir}/systemctl-calls" || fail 'unexpected build or production restart'
echo 'Split deployment, same-runtime app rollback, recovery, and dev isolation tests passed.'
