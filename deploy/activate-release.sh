#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: deploy/activate-release.sh --project-root ABS_PATH --release RELEASE

Installs/updates the user Quadlet from a fully verified release. The immutable
IMAGE_REF and DATA_SCHEMA are read from RELEASE_DIR/release.env. If startup or
health verification fails, the prior unit is restored automatically.
EOF
}

project_root=""
release=""
while (($#)); do
  case "$1" in
    --project-root) project_root="${2:?missing value for --project-root}"; shift 2 ;;
    --release) release="${2:?missing value for --release}"; shift 2 ;;
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
  --entrypoint /usr/bin/Rscript \
  --volume "${data_dir}:/srv/islet-explorer/data:ro" \
  --env "ISLET_DATA_SCHEMA=${data_schema}" \
  "${image_ref}" /srv/islet-explorer/scripts/verify_release_data.R \
  --root /srv/islet-explorer/data --mode core

escape_sed() { printf '%s' "$1" | sed 's/[&|\\]/\\&/g'; }
rendered="$(mktemp)"
unit_dir="${HOME}/.config/containers/systemd"
unit_path="${unit_dir}/islet-explorer.container"
backup=""
trap 'rm -f -- "${rendered}" "${backup:-}"' EXIT

sed \
  -e "s|@@RELEASE@@|$(escape_sed "${release}")|g" \
  -e "s|@@IMAGE_REF@@|$(escape_sed "${image_ref}")|g" \
  -e "s|@@DATA_SCHEMA@@|$(escape_sed "${data_schema}")|g" \
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
  echo "Activation failed; restoring the prior Quadlet unit." >&2
  if [[ -n "${backup}" && -f "${backup}" ]]; then
    install -m 0644 "${backup}" "${unit_path}"
    systemctl --user daemon-reload || true
    systemctl --user restart islet-explorer.service || true
  else
    systemctl --user stop islet-explorer.service || true
    rm -f -- "${unit_path}"
    systemctl --user daemon-reload || true
  fi
}

if ! systemctl --user daemon-reload; then
  restore_previous_unit
  exit 70
fi
# Quadlet's [Install] section wires the generated service into default.target;
# generated units themselves are not enabled with `systemctl enable`.
if ! systemctl --user restart islet-explorer.service; then
  restore_previous_unit
  exit 70
fi
if ! "${script_dir}/verify-deployment.sh" --url \
  "http://127.0.0.1:3838/islet-explorer/" --container islet-explorer; then
  restore_previous_unit
  exit 70
fi

if [[ -f "${state_dir}/active.state" ]]; then
  cp -- "${state_dir}/active.state" "${state_dir}/previous.state"
fi
state_tmp="$(mktemp "${state_dir}/.active.state.XXXXXX")"
{
  printf 'RELEASE=%s\n' "${release}"
  printf 'IMAGE_REF=%s\n' "${image_ref}"
  printf 'ACTIVATED_AT=%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
} > "${state_tmp}"
mv -f -- "${state_tmp}" "${state_dir}/active.state"

echo "Activated Islet Explorer release ${release} (${image_ref})."
