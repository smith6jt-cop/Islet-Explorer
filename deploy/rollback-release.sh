#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: deploy/rollback-release.sh --project-root ABS_PATH [--instance NAME]" >&2
}

project_root=""
instance="islet-explorer"
while (($#)); do
  case "$1" in
    --project-root) project_root="${2:?missing value for --project-root}"; shift 2 ;;
    --instance) instance="${2:?missing value for --instance}"; shift 2 ;;
    *) usage; exit 64 ;;
  esac
done
[[ "${project_root}" == /* && "${project_root}" =~ ^[A-Za-z0-9._/-]+$ ]] || {
  echo "ERROR: --project-root must be an absolute path without whitespace." >&2
  exit 64
}
[[ "${instance}" =~ ^[a-z][a-z0-9-]{0,40}$ ]] || {
  echo "ERROR: --instance must be a lowercase name of letters, digits and hyphens." >&2
  exit 64
}

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
# The production instance records deployments/previous.state; every other
# instance records <instance>.previous.state (see activate-release.sh).
if [[ "${instance}" == islet-explorer ]]; then
  previous_state="${project_root}/runtime/deployments/previous.state"
else
  previous_state="${project_root}/runtime/deployments/${instance}.previous.state"
fi
[[ -f "${previous_state}" ]] || {
  echo "ERROR: no previously successful deployment is recorded for ${instance}." >&2
  exit 66
}
read_state() {
  awk -F= -v key="$1" '$1 == key {sub(/^[^=]*=/, ""); print; exit}' "${previous_state}"
}
previous_release="$(read_state RELEASE)"
[[ "${previous_release}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
  echo "ERROR: previous deployment state is malformed." >&2
  exit 65
}

# The previous deployment is restored on the address, port and memory cap it
# was published with; deployments recorded before those fields existed were
# loopback-only on 3838 and take the current default cap.
previous_address="$(read_state PUBLISH_ADDRESS)"
[[ "${previous_address}" =~ ^([0-9]{1,3}\.){3}[0-9]{1,3}$ ]] || previous_address="127.0.0.1"
previous_port="$(read_state PORT)"
[[ "${previous_port}" =~ ^[0-9]{1,5}$ ]] || previous_port="3838"
previous_memory="$(read_state MEMORY)"
[[ "${previous_memory}" =~ ^[1-9][0-9]*[mg]$ ]] || previous_memory="9g"

echo "Rolling back ${instance} to release ${previous_release} on ${previous_address}:${previous_port}."
exec "${script_dir}/activate-release.sh" \
  --project-root "${project_root}" --release "${previous_release}" \
  --instance "${instance}" --port "${previous_port}" --memory "${previous_memory}" \
  --publish-address "${previous_address}"
