#!/usr/bin/env bash
set -euo pipefail

expected_uid="$(id -u shiny)"
expected_gid="$(id -g shiny)"
actual_uid="$(id -u)"
actual_gid="$(id -g)"
readonly expected_uid expected_gid actual_uid actual_gid
readonly runtime_mode="${ISLET_RUNTIME_MODE:-production}"

if [[ "${actual_uid}" != "${expected_uid}" || "${actual_gid}" != "${expected_gid}" ]]; then
  echo "ERROR: container must run as shiny (${expected_uid}:${expected_gid}); got ${actual_uid}:${actual_gid}." >&2
  echo "Use rootless keep-id mapping as documented in deploy/README.md." >&2
  exit 70
fi

if [[ "${runtime_mode,,}" == "production" ]]; then
  if [[ "${ISLET_ENABLE_AI:-FALSE}" != "FALSE" || "${SHOW_STATISTICS_TAB:-0}" != "0" ||
        "${SHOW_INFLAMMATION_TAB:-0}" != "0" || "${DEBUG_CREDENTIALS:-0}" != "0" ||
        "${VIEWER_DEBUG:-0}" != "0" ]]; then
    echo "ERROR: the reviewer release must keep optional features and debug modes disabled." >&2
    exit 78
  fi
  if [[ -n "${KEY:-}" || -n "${OPENAI_API_KEY:-}" ]]; then
    echo "ERROR: an AI credential was supplied to an AI-disabled production container." >&2
    exit 78
  fi
fi

readonly data_root="${ISLET_DATA_ROOT:-/srv/islet-explorer/data}"
readonly image_root="${LOCAL_IMAGE_ROOT:-/srv/islet-explorer/app/shiny_app/www/local_images}"
readonly runtime_tmp="${TMPDIR:-/tmp/islet-explorer}"
readonly runtime_r_user="${runtime_tmp}/r-user"
readonly runtime_cache="${runtime_tmp}/cache"
readonly runtime_matplotlib="${runtime_tmp}/matplotlib"
readonly runtime_numba="${runtime_tmp}/numba"

for required_dir in "${data_root}" "${image_root}" /var/log/shiny-server; do
  if [[ ! -d "${required_dir}" || ! -r "${required_dir}" ]]; then
    echo "ERROR: required runtime directory is not readable: ${required_dir}" >&2
    exit 66
  fi
done

mkdir -p "${runtime_tmp}" "${runtime_r_user}" "${runtime_cache}" \
  "${runtime_matplotlib}" "${runtime_numba}"
for writable_dir in "${runtime_tmp}" "${runtime_r_user}" "${runtime_cache}" \
  "${runtime_matplotlib}" "${runtime_numba}" /var/log/shiny-server; do
  if [[ ! -w "${writable_dir}" ]]; then
    echo "ERROR: runtime directory is not writable by ${actual_uid}:${actual_gid}: ${writable_dir}" >&2
    exit 73
  fi
done

export R_USER="${runtime_r_user}"
export XDG_CACHE_HOME="${runtime_cache}"
export MPLCONFIGDIR="${runtime_matplotlib}"
export NUMBA_CACHE_DIR="${runtime_numba}"
export ISLET_RUNTIME_MODE="${runtime_mode}"
export ISLET_DATA_ROOT="${data_root}"
export LOCAL_IMAGE_ROOT="${image_root}"
export TMPDIR="${runtime_tmp}"

# Shiny Server starts R with only HOME/LANG/PATH. The immutable wrapper restores
# this allowlisted profile before R can read user startup files.
/usr/local/bin/islet-r-worker --write-profile
export R=/usr/local/bin/islet-r-worker

exec "$@"
