#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: deploy/activate-release.sh --project-root ABS_PATH --release RELEASE
         [--publish-address IPV4]

Installs/updates the user Quadlet from a fully verified release. The immutable
IMAGE_REF and DATA_SCHEMA are read from RELEASE_DIR/release.env. If startup or
health verification fails, the prior unit is restored automatically.
--publish-address publishes the backend on that VM address instead of the
loopback 127.0.0.1: use it only for the private address UFRC's reverse proxy
connects to, and only once UFRC has confirmed that the port is reachable from
the proxy alone, because the application relies on the proxy's Basic Auth.
EOF
}

project_root=""
release=""
publish_address="127.0.0.1"
while (($#)); do
  case "$1" in
    --project-root) project_root="${2:?missing value for --project-root}"; shift 2 ;;
    --release) release="${2:?missing value for --release}"; shift 2 ;;
    --publish-address) publish_address="${2:?missing value for --publish-address}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 64 ;;
  esac
done

[[ "${project_root}" == /* && "${project_root}" =~ ^[A-Za-z0-9._/-]+$ ]] || {
  echo "ERROR: --project-root must be an absolute path without whitespace." >&2
  exit 64
}
[[ "${release}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
  echo "ERROR: invalid release identifier." >&2
  exit 64
}
[[ "${publish_address}" =~ ^([0-9]{1,3}\.){3}[0-9]{1,3}$ ]] || {
  echo "ERROR: --publish-address must be an IPv4 address literal." >&2
  exit 64
}
for command_name in podman systemctl curl install sed awk; do
  command -v "${command_name}" >/dev/null || {
    echo "ERROR: required command not found: ${command_name}" >&2
    exit 69
  }
done

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
release_dir="${project_root}/releases/${release}"
# Re-hash immediately before service replacement. A stamp alone cannot prove
# that a host-side file has not changed since the earlier staging check.
"${script_dir}/verify-release.sh" "${release_dir}"

read_meta() {
  local key="$1"
  awk -F= -v key="${key}" '$1 == key {sub(/^[^=]*=/, ""); print; exit}' "${release_dir}/release.env"
}
manifest_release="$(read_meta RELEASE)"
image_ref="$(read_meta IMAGE_REF)"
data_schema="$(read_meta DATA_SCHEMA)"
git_commit="$(read_meta GIT_COMMIT)"
[[ "${manifest_release}" == "${release}" ]] || {
  echo "ERROR: requested release does not match release.env (${manifest_release})." >&2
  exit 65
}
[[ "${image_ref}" =~ ^sha256:[0-9a-f]{64}$ ||
   "${image_ref}" =~ ^[A-Za-z0-9._:/-]+@sha256:[0-9a-f]{64}$ ]] || {
  echo "ERROR: release image is not immutable: ${image_ref}" >&2
  exit 65
}

podman image exists "${image_ref}" || {
  if [[ "${image_ref}" == *@sha256:* ]]; then
    echo "Image is absent locally; pulling immutable digest ${image_ref}."
    podman pull "${image_ref}"
  else
    echo "ERROR: local image ID is absent: ${image_ref}" >&2
    exit 66
  fi
}

container_uid="$(podman run --rm --entrypoint /usr/bin/id "${image_ref}" -u shiny)"
container_gid="$(podman run --rm --entrypoint /usr/bin/id "${image_ref}" -g shiny)"
[[ "${container_uid}" =~ ^[0-9]+$ && "${container_gid}" =~ ^[0-9]+$ ]] || {
  echo "ERROR: could not determine the image's shiny UID/GID." >&2
  exit 70
}

data_dir="${release_dir}/data"
image_dir="${release_dir}/data/images"
runtime_dir="${project_root}/runtime"
log_dir="${runtime_dir}/logs"
tmp_dir="${runtime_dir}/tmp"
state_dir="${runtime_dir}/deployments"

mkdir -p -- "${log_dir}" "${tmp_dir}" "${state_dir}"
chmod 0750 "${runtime_dir}" "${log_dir}" "${tmp_dir}" "${state_dir}"
for readable_dir in "${data_dir}" "${image_dir}"; do
  [[ -r "${readable_dir}" && -x "${readable_dir}" ]] || {
    echo "ERROR: service account cannot read ${readable_dir}." >&2
    exit 73
  }
done
for writable_dir in "${log_dir}" "${tmp_dir}" "${state_dir}"; do
  [[ -w "${writable_dir}" && -x "${writable_dir}" ]] || {
    echo "ERROR: service account cannot write ${writable_dir}." >&2
    exit 73
  }
done

image_revision="$(podman image inspect --format \
  '{{ index .Labels "org.opencontainers.image.revision" }}' "${image_ref}")"
[[ "${image_revision}" == "${git_commit}" ]] || {
  echo "ERROR: image revision ${image_revision} does not match release commit ${git_commit}." >&2
  exit 65
}

# Confirm the installed Podman supports the exact rootless keep-id mapping
# before replacing a running service.
podman run --rm \
  --userns="keep-id:uid=${container_uid},gid=${container_gid}" \
  --user "${container_uid}:${container_gid}" \
  --entrypoint /bin/sh \
  --volume "${data_dir}:/probe-data:ro" \
  --volume "${log_dir}:/probe-log:rw" \
  "${image_ref}" -ceu \
  'test "$(id -u)" = "'"${container_uid}"'"; test -r /probe-data; touch /probe-log/.uid-probe; rm /probe-log/.uid-probe'

manifest_identity="$(podman run --rm \
  --userns="keep-id:uid=${container_uid},gid=${container_gid}" \
  --user "${container_uid}:${container_gid}" \
  --entrypoint /opt/islet-explorer-py/bin/python \
  --volume "${data_dir}:/srv/islet-explorer/data:ro" \
  "${image_ref}" -c \
  'import json; m=json.load(open("/srv/islet-explorer/data/manifest.json", encoding="utf-8")); print(str(m.get("schema_version", ""))+"\t"+str(m.get("git_commit", "")))')"
IFS=$'\t' read -r manifest_schema manifest_commit <<< "${manifest_identity}"
[[ "${manifest_schema}" == "${data_schema}" && "${manifest_commit}" == "${git_commit}" ]] || {
  echo "ERROR: data manifest identity does not match release.env." >&2
  echo "  schema=${manifest_schema} commit=${manifest_commit}" >&2
  exit 65
}

# Exercise the same schema, artifact, hash, and Parquet row-count validator used
# by the application before replacing the active unit. The outer release hash
# has already covered every TIFF; core mode avoids hashing that tree a third time.
podman run --rm \
  --userns="keep-id:uid=${container_uid},gid=${container_gid}" \
  --user "${container_uid}:${container_gid}" \
  --entrypoint Rscript \
  --volume "${data_dir}:/srv/islet-explorer/data:ro" \
  --env "ISLET_DATA_SCHEMA=${data_schema}" \
  "${image_ref}" /srv/islet-explorer/scripts/verify_release_data.R \
  --root /srv/islet-explorer/data --mode core

# A rootless CPU quota needs the cpu controller delegated to the user session
# (cgroups v2). PubApps VMs delegate memory and pids only (2026-10-22,
# pubmaigan1), and their 2-vCPU size bounds CPU by itself; a quota the kernel
# cannot apply would keep the container from starting at all.
cpu_limit=""
controllers_file="/sys/fs/cgroup/user.slice/user-$(id -u).slice/user@$(id -u).service/cgroup.controllers"
if [[ -r "${controllers_file}" ]] && grep -qw cpu "${controllers_file}"; then
  cpu_limit="--cpus=2 "
else
  echo "NOTE: the cpu controller is not delegated to this user session; the container keeps the host's vCPU bound."
fi

escape_sed() { printf '%s' "$1" | sed 's/[&|\\]/\\&/g'; }
rendered="$(mktemp)"
unit_dir="${XDG_CONFIG_HOME:-${HOME}/.config}/containers/systemd"
unit_path="${unit_dir}/islet-explorer.container"
backup=""
trap 'rm -f -- "${rendered}" "${backup:-}"' EXIT

sed \
  -e "s|@@RELEASE@@|$(escape_sed "${release}")|g" \
  -e "s|@@IMAGE_REF@@|$(escape_sed "${image_ref}")|g" \
  -e "s|@@DATA_SCHEMA@@|$(escape_sed "${data_schema}")|g" \
  -e "s|@@CPU_LIMIT@@|$(escape_sed "${cpu_limit}")|g" \
  -e "s|@@PUBLISH_ADDRESS@@|$(escape_sed "${publish_address}")|g" \
  -e "s|@@DATA_DIR@@|$(escape_sed "${data_dir}")|g" \
  -e "s|@@IMAGE_DIR@@|$(escape_sed "${image_dir}")|g" \
  -e "s|@@LOG_DIR@@|$(escape_sed "${log_dir}")|g" \
  -e "s|@@TMP_DIR@@|$(escape_sed "${tmp_dir}")|g" \
  -e "s|@@CONTAINER_UID@@|${container_uid}|g" \
  -e "s|@@CONTAINER_GID@@|${container_gid}|g" \
  "${script_dir}/islet-explorer.container.in" > "${rendered}"

mkdir -p -- "${unit_dir}"
if [[ -f "${unit_path}" ]]; then
  backup="$(mktemp)"
  cp -- "${unit_path}" "${backup}"
fi
install -m 0644 "${rendered}" "${unit_path}"

restore_previous_unit() {
  if [[ -n "${backup}" && -f "${backup}" ]]; then
    echo "Activation failed; restoring and checking the prior Quadlet service." >&2
    local previous_image previous_address
    previous_image="$(awk -F= '$1 == "Image" {sub(/^[^=]*=/, ""); print; exit}' "${backup}")"
    [[ -n "${previous_image}" ]] || return 1
    # The prior unit is verified on the address it publishes (PublishPort=ADDR:3838:3838).
    previous_address="$(awk -F= '$1 == "PublishPort" {sub(/^[^=]*=/, ""); split($0, a, ":"); print a[1]; exit}' "${backup}")"
    [[ "${previous_address}" =~ ^([0-9]{1,3}\.){3}[0-9]{1,3}$ ]] || previous_address="127.0.0.1"
    install -m 0644 "${backup}" "${unit_path}" || return 1
    systemctl --user daemon-reload || return 1
    systemctl --user restart islet-explorer.service || return 1
    "${script_dir}/verify-deployment.sh" --url \
      "http://${previous_address}:3838/islet-explorer/" --container islet-explorer \
      --publish-address "${previous_address}" \
      --expected-image "${previous_image}" || return 1
    echo "Prior release restored and verified healthy; release state is unchanged." >&2
  else
    echo "First activation failed; stopping and removing the candidate service." >&2
    systemctl --user stop islet-explorer.service || return 1
    rm -f -- "${unit_path}" || return 1
    systemctl --user daemon-reload || return 1
    if podman container exists islet-explorer; then
      local still_running
      still_running="$(podman inspect --format '{{.State.Running}}' islet-explorer)" || return 1
      [[ "${still_running}" == false ]] || return 1
    fi
    echo "Candidate service removed; no release was activated." >&2
  fi
}

activation_failed() {
  if ! restore_previous_unit; then
    echo "ERROR: activation and recovery failed; service availability is unverified. Inspect the user service before retrying." >&2
    exit 71
  fi
  exit 70
}

if ! systemctl --user daemon-reload; then
  activation_failed
fi
# Quadlet's [Install] section wires the generated service into default.target;
# generated units themselves are not enabled with `systemctl enable`.
if ! systemctl --user restart islet-explorer.service; then
  activation_failed
fi
if ! "${script_dir}/verify-deployment.sh" --url \
  "http://${publish_address}:3838/islet-explorer/" --container islet-explorer \
  --publish-address "${publish_address}" \
  --expected-image "${image_ref}"; then
  activation_failed
fi

if [[ -f "${state_dir}/active.state" ]]; then
  cp -- "${state_dir}/active.state" "${state_dir}/previous.state"
fi
state_tmp="$(mktemp "${state_dir}/.active.state.XXXXXX")"
{
  printf 'RELEASE=%s\n' "${release}"
  printf 'IMAGE_REF=%s\n' "${image_ref}"
  printf 'PUBLISH_ADDRESS=%s\n' "${publish_address}"
  printf 'ACTIVATED_AT=%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
} > "${state_tmp}"
mv -f -- "${state_tmp}" "${state_dir}/active.state"

echo "Activated Islet Explorer release ${release} (${image_ref}) on ${publish_address}:3838."
