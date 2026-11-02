#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: deploy/activate-release.sh --project-root ABS_PATH --release RELEASE
         [--publish-address IPV4] [--instance NAME] [--port PORT] [--memory SIZE]
         [--docroot ABS_PATH] [--local-data ABS_PATH]

Installs/updates the user Quadlet from a fully verified release. The immutable
IMAGE_REF and DATA_SCHEMA are read from RELEASE_DIR/release.env. If startup or
health verification fails, the prior unit and its active/inactive state are
restored automatically. A previously stopped service is never restarted by
failure recovery.
--publish-address publishes the backend on that VM address instead of the
loopback 127.0.0.1: use it only for the private address UFRC's reverse proxy
connects to, and only where UFRC's network keeps that port from every other
host, because the application relies on the proxy's Basic Auth.
--instance names an independent service beside the production one (default
islet-explorer): its own Quadlet unit, container, host port (--port, default
3838), memory cap (--memory, default 9g), deployment state and runtime
directories. The staging instance is islet-explorer-dev on port 3839 (4g); the
caps of all instances together must leave the 15 GiB VM its own headroom.
--docroot names the Apache DocumentRoot UFRC's proxy serves the image aliases
from itself (DOCROOT/islet-explorer/{local_images,images,img_data}, excluded
from proxying in the vhost; prod/html and dev/html under the project root).
Once the new service is verified, each alias is replaced atomically by a symlink
to this release's data/images; if a later step fails, the aliases follow the
restored release. An alias that is a real directory is refused before any
service is touched.
--local-data names a directory on the VM's local disk (for example
/local/maigan/releases). The release's data/geometry, data/cells and
data/spatial trees, whose many small reads are slow over NFS (a drill-down's
GeoPackage window took 0.2 s there against 17 ms on local disk), are copied to
LOCAL_DATA/RELEASE/data/, checked file by file against MANIFEST.sha256, and
mounted over the same paths inside the container; the image tree stays on the
shared filesystem. A mismatch stops the activation before any service is
touched. Copies of earlier releases are kept so a restored or rolled-back unit
still finds its data; remove a retired release's copy by hand.
EOF
}

project_root=""
release=""
publish_address="127.0.0.1"
instance="islet-explorer"
port="3838"
memory="9g"
docroot=""
local_data=""
while (($#)); do
  case "$1" in
    --project-root) project_root="${2:?missing value for --project-root}"; shift 2 ;;
    --release) release="${2:?missing value for --release}"; shift 2 ;;
    --publish-address) publish_address="${2:?missing value for --publish-address}"; shift 2 ;;
    --instance) instance="${2:?missing value for --instance}"; shift 2 ;;
    --port) port="${2:?missing value for --port}"; shift 2 ;;
    --memory) memory="${2:?missing value for --memory}"; shift 2 ;;
    --docroot) docroot="${2:?missing value for --docroot}"; shift 2 ;;
    --local-data) local_data="${2:?missing value for --local-data}"; shift 2 ;;
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
[[ "${instance}" =~ ^[a-z][a-z0-9-]{0,40}$ ]] || {
  echo "ERROR: --instance must be a lowercase name of letters, digits and hyphens." >&2
  exit 64
}
[[ "${port}" =~ ^[0-9]{1,5}$ ]] && port="$((10#${port}))" && ((port >= 1 && port <= 65535)) || {
  echo "ERROR: --port must be between 1 and 65535." >&2
  exit 64
}
[[ "${memory}" =~ ^[1-9][0-9]*[mg]$ ]] || {
  echo "ERROR: --memory must be a size such as 9g or 4096m." >&2
  exit 64
}
[[ -z "${docroot}" || ( "${docroot}" == /* && "${docroot}" =~ ^[A-Za-z0-9._/-]+$ ) ]] || {
  echo "ERROR: --docroot must be an absolute path without whitespace." >&2
  exit 64
}
[[ -z "${local_data}" || ( "${local_data}" == /* && "${local_data}" =~ ^[A-Za-z0-9._/-]+$ ) ]] || {
  echo "ERROR: --local-data must be an absolute path without whitespace." >&2
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
state_dir="${runtime_dir}/deployments"
# The production instance keeps the original layout (runtime/logs, runtime/tmp,
# deployments/active.state); every other instance has its own runtime tree and
# state files, so two services never share a writable directory.
if [[ "${instance}" == islet-explorer ]]; then
  instance_runtime="${runtime_dir}"
  state_name="active"
  previous_name="previous"
else
  instance_runtime="${runtime_dir}/${instance}"
  state_name="${instance}.active"
  previous_name="${instance}.previous"
fi
log_dir="${instance_runtime}/logs"
tmp_dir="${instance_runtime}/tmp"

mkdir -p -- "${log_dir}" "${tmp_dir}" "${state_dir}"
chmod 0750 "${runtime_dir}" "${instance_runtime}" "${log_dir}" "${tmp_dir}" "${state_dir}"
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
# The proxy's static image aliases (shiny-server.conf serves the same three
# names from the container). They must already be symlinks, or absent: nothing
# here deletes a real directory.
docroot_aliases=(local_images images img_data)
if [[ -n "${docroot}" ]]; then
  [[ -d "${docroot}" && -w "${docroot}" && -x "${docroot}" ]] || {
    echo "ERROR: --docroot ${docroot} is not a writable directory." >&2
    exit 73
  }
  for alias in "${docroot_aliases[@]}"; do
    alias_path="${docroot}/islet-explorer/${alias}"
    if [[ -e "${alias_path}" && ! -L "${alias_path}" ]]; then
      echo "ERROR: ${alias_path} exists and is not a symlink; the proxy's image aliases must be symlinks into releases/." >&2
      exit 65
    fi
  done
fi

unit_dir="${XDG_CONFIG_HOME:-${HOME}/.config}/containers/systemd"
unit_path="${unit_dir}/${instance}.container"
service_name="${instance}.service"
# Two instances cannot publish the same host port; refuse before touching any
# service rather than letting the restart fail and recovery run.
if [[ -d "${unit_dir}" ]]; then
  for other_unit in "${unit_dir}"/*.container; do
    [[ -f "${other_unit}" && "${other_unit}" != "${unit_path}" ]] || continue
    if grep -Eq "^PublishPort=[0-9.]+:${port}:3838$" "${other_unit}"; then
      echo "ERROR: host port ${port} is already published by $(basename "${other_unit}")." >&2
      exit 65
    fi
  done
fi

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
# Local-disk copies of the small-read data trees (--local-data). rsync brings
# the copy up to date (unchanged files are skipped), then every file is checked
# against the release manifest; a mismatch stops the activation here, before
# any service is touched. The copies are mounted over the same paths inside the
# container by the Volume lines rendered below.
local_data_trees=(geometry cells spatial)
local_data_volumes=""
local_data_mount=""
if [[ -n "${local_data}" ]]; then
  command -v rsync >/dev/null || { echo "ERROR: --local-data needs rsync." >&2; exit 69; }
  local_release_data="${local_data}/${release}/data"
  mkdir -p -- "${local_release_data}" || { echo "ERROR: cannot create ${local_release_data}." >&2; exit 73; }
  for tree in "${local_data_trees[@]}"; do
    [[ -d "${data_dir}/${tree}" ]] || { echo "ERROR: the release has no data/${tree} tree to copy." >&2; exit 65; }
    rsync -a --delete -- "${data_dir}/${tree}/" "${local_release_data}/${tree}/" || {
      echo "ERROR: copying data/${tree} to ${local_release_data} failed." >&2
      exit 74
    }
  done
  manifest_subset="$(mktemp)"
  grep -E "  data/(geometry|cells|spatial)/" "${release_dir}/MANIFEST.sha256" > "${manifest_subset}" || true
  [[ -s "${manifest_subset}" ]] || {
    rm -f -- "${manifest_subset}"
    echo "ERROR: MANIFEST.sha256 lists no files under data/{geometry,cells,spatial}." >&2
    exit 65
  }
  if ! (cd -- "${local_data}/${release}" && sha256sum --check --quiet --strict "${manifest_subset}"); then
    rm -f -- "${manifest_subset}"
    echo "ERROR: the local copies under ${local_release_data} do not match MANIFEST.sha256." >&2
    exit 65
  fi
  rm -f -- "${manifest_subset}"
  for tree in "${local_data_trees[@]}"; do
    local_data_volumes+="Volume=${local_release_data}/${tree}:/srv/islet-explorer/data/${tree}:ro"$'\n'
  done
  local_data_mount=" ${local_release_data}"
  echo "Local copies of data/{geometry,cells,spatial} under ${local_release_data} verified against MANIFEST.sha256."
fi

cpu_limit=""
controllers_file="/sys/fs/cgroup/user.slice/user-$(id -u).slice/user@$(id -u).service/cgroup.controllers"
if [[ -r "${controllers_file}" ]] && grep -qw cpu "${controllers_file}"; then
  cpu_limit="--cpus=2 "
else
  echo "NOTE: the cpu controller is not delegated to this user session; the container keeps the host's vCPU bound."
fi

escape_sed() { printf '%s' "$1" | sed 's/[&|\\]/\\&/g'; }
rendered="$(mktemp)"
backup=""
previous_active=false
if systemctl --user is-active --quiet "${service_name}"; then
  previous_active=true
fi
trap 'rm -f -- "${rendered}" "${backup:-}"' EXIT

sed \
  -e "s|@@RELEASE@@|$(escape_sed "${release}")|g" \
  -e "s|@@INSTANCE@@|$(escape_sed "${instance}")|g" \
  -e "s|@@IMAGE_REF@@|$(escape_sed "${image_ref}")|g" \
  -e "s|@@DATA_SCHEMA@@|$(escape_sed "${data_schema}")|g" \
  -e "s|@@CPU_LIMIT@@|$(escape_sed "${cpu_limit}")|g" \
  -e "s|@@MEMORY@@|$(escape_sed "${memory}")|g" \
  -e "s|@@PUBLISH_ADDRESS@@|$(escape_sed "${publish_address}")|g" \
  -e "s|@@PORT@@|${port}|g" \
  -e "s|@@DATA_DIR@@|$(escape_sed "${data_dir}")|g" \
  -e "s|@@IMAGE_DIR@@|$(escape_sed "${image_dir}")|g" \
  -e "s|@@LOG_DIR@@|$(escape_sed "${log_dir}")|g" \
  -e "s|@@TMP_DIR@@|$(escape_sed "${tmp_dir}")|g" \
  -e "s|@@CONTAINER_UID@@|${container_uid}|g" \
  -e "s|@@CONTAINER_GID@@|${container_gid}|g" \
  -e "s|@@LOCAL_DATA_MOUNT@@|$(escape_sed "${local_data_mount}")|g" \
  "${script_dir}/islet-explorer.container.in" > "${rendered}.sed"
# The local-data placeholder line becomes one Volume line per copied tree, or
# nothing at all when --local-data is not given.
awk -v volumes="${local_data_volumes}" \
  '$0 == "@@LOCAL_DATA_VOLUMES@@" { printf "%s", volumes; next } { print }' \
  "${rendered}.sed" > "${rendered}"
rm -f -- "${rendered}.sed"

mkdir -p -- "${unit_dir}"
if [[ -f "${unit_path}" ]]; then
  backup="$(mktemp)"
  cp -- "${unit_path}" "${backup}"
fi
install -m 0644 "${rendered}" "${unit_path}"

# Point the proxy's image aliases at one release's images. Each alias is
# replaced by rename(2), so the proxy never sees a missing or half-written link
# (its NFS lookup cache may still hand out the old target for up to a minute).
publish_docroot_links() {
  local target="$1" alias link tmp_link
  mkdir -p -- "${docroot}/islet-explorer" || return 1
  for alias in "${docroot_aliases[@]}"; do
    link="${docroot}/islet-explorer/${alias}"
    tmp_link="${docroot}/islet-explorer/.${alias}.${instance}.tmp"
    ln -sfn -- "${target}" "${tmp_link}" || return 1
    mv -Tf -- "${tmp_link}" "${link}" || { rm -f -- "${tmp_link}"; return 1; }
  done
  echo "Proxy image aliases under ${docroot}/islet-explorer now point at ${target}."
}

restore_previous_unit() {
  if [[ -n "${backup}" && -f "${backup}" ]]; then
    echo "Activation failed; restoring and checking the prior Quadlet service." >&2
    local previous_image previous_publish previous_address previous_port
    previous_image="$(awk -F= '$1 == "Image" {sub(/^[^=]*=/, ""); print; exit}' "${backup}")"
    [[ -n "${previous_image}" ]] || return 1
    # The prior unit is verified on the address and port it publishes
    # (PublishPort=ADDR:PORT:3838).
    previous_publish="$(awk -F= '$1 == "PublishPort" {sub(/^[^=]*=/, ""); print; exit}' "${backup}")"
    previous_address="${previous_publish%%:*}"
    previous_port="${previous_publish#*:}"
    previous_port="${previous_port%%:*}"
    [[ "${previous_address}" =~ ^([0-9]{1,3}\.){3}[0-9]{1,3}$ ]] || previous_address="127.0.0.1"
    [[ "${previous_port}" =~ ^[0-9]{1,5}$ ]] || previous_port="${port}"
    install -m 0644 "${backup}" "${unit_path}" || return 1
    if [[ "${previous_active}" != true ]]; then
      # A service may have been stopped for security remediation. Restoring
      # its configuration must not start that old image on candidate failure.
      systemctl --user stop "${service_name}" || return 1
      systemctl --user daemon-reload || return 1
      if podman container exists "${instance}"; then
        local still_running
        still_running="$(podman inspect --format '{{.State.Running}}' "${instance}")" || return 1
        [[ "${still_running}" == false ]] || return 1
      fi
      echo "Prior unit restored and left stopped; release state is unchanged." >&2
      return 0
    fi
    systemctl --user daemon-reload || return 1
    systemctl --user restart "${service_name}" || return 1
    "${script_dir}/verify-deployment.sh" --url \
      "http://${previous_address}:${previous_port}/islet-explorer/" --container "${instance}" \
      --publish-address "${previous_address}" \
      --expected-image "${previous_image}" || return 1
    if [[ -n "${docroot}" ]]; then
      # The restored unit mounts its release's images on the container's
      # local_images path; the proxy aliases follow it.
      local previous_images
      previous_images="$(awk '/^Volume=.*:\/srv\/islet-explorer\/app\/shiny_app\/www\/local_images:ro$/ {sub(/^Volume=/, ""); sub(/:\/srv\/islet-explorer\/app\/shiny_app\/www\/local_images:ro$/, ""); print; exit}' "${backup}")"
      if [[ -n "${previous_images}" ]]; then
        publish_docroot_links "${previous_images}" || return 1
      else
        echo "NOTE: the restored unit names no image directory; the proxy image aliases are unchanged." >&2
      fi
    fi
    echo "Prior release restored and verified healthy; release state is unchanged." >&2
  else
    echo "First activation failed; stopping and removing the candidate service." >&2
    systemctl --user stop "${service_name}" || return 1
    rm -f -- "${unit_path}" || return 1
    systemctl --user daemon-reload || return 1
    if podman container exists "${instance}"; then
      local still_running
      still_running="$(podman inspect --format '{{.State.Running}}' "${instance}")" || return 1
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
if ! systemctl --user restart "${service_name}"; then
  activation_failed
fi
if ! "${script_dir}/verify-deployment.sh" --url \
  "http://${publish_address}:${port}/islet-explorer/" --container "${instance}" \
  --publish-address "${publish_address}" \
  --expected-image "${image_ref}"; then
  activation_failed
fi
if [[ -n "${docroot}" ]] && ! publish_docroot_links "${image_dir}"; then
  echo "ERROR: could not point the proxy image aliases at ${image_dir}." >&2
  activation_failed
fi

if [[ -f "${state_dir}/${state_name}.state" ]]; then
  cp -- "${state_dir}/${state_name}.state" "${state_dir}/${previous_name}.state"
fi
state_tmp="$(mktemp "${state_dir}/.${state_name}.state.XXXXXX")"
{
  printf 'RELEASE=%s\n' "${release}"
  printf 'IMAGE_REF=%s\n' "${image_ref}"
  printf 'INSTANCE=%s\n' "${instance}"
  printf 'PUBLISH_ADDRESS=%s\n' "${publish_address}"
  printf 'PORT=%s\n' "${port}"
  printf 'MEMORY=%s\n' "${memory}"
  printf 'DOCROOT=%s\n' "${docroot}"
  printf 'LOCAL_DATA=%s\n' "${local_data}"
  printf 'ACTIVATED_AT=%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
} > "${state_tmp}"
mv -f -- "${state_tmp}" "${state_dir}/${state_name}.state"

docroot_note=""
[[ -n "${docroot}" ]] && docroot_note=", proxy image aliases in ${docroot}"
local_data_note=""
[[ -n "${local_data}" ]] && local_data_note=", local data copies in ${local_data}/${release}"
echo "Activated Islet Explorer release ${release} (${image_ref}) as ${instance} on ${publish_address}:${port} (memory ${memory})${docroot_note}${local_data_note}."
