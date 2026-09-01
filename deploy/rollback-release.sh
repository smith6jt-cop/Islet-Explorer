#!/usr/bin/env bash
set -euo pipefail

if (($# != 2)) || [[ "$1" != "--project-root" ]]; then
  echo "Usage: deploy/rollback-release.sh --project-root ABS_PATH" >&2
  exit 64
fi
project_root="$2"
[[ "${project_root}" == /* && "${project_root}" =~ ^[A-Za-z0-9._/-]+$ ]] || {
  echo "ERROR: --project-root must be an absolute path without whitespace." >&2
  exit 64
}

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
previous_state="${project_root}/runtime/deployments/previous.state"
[[ -f "${previous_state}" ]] || {
  echo "ERROR: no previously successful deployment is recorded." >&2
  exit 66
}
previous_release="$(awk -F= '$1 == "RELEASE" {print $2; exit}' "${previous_state}")"
[[ "${previous_release}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
  echo "ERROR: previous deployment state is malformed." >&2
  exit 65
}

echo "Rolling back to release ${previous_release}."
exec "${script_dir}/activate-release.sh" \
  --project-root "${project_root}" --release "${previous_release}"
