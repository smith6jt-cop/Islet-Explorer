#!/usr/bin/env bash
# Manual rootless Podman smoke launcher. Production uses activate-release.sh and
# the generated Quadlet unit; this command intentionally has no Apptainer path.
set -euo pipefail

project_root="${PROJECT_ROOT:?set PROJECT_ROOT to /pubapps/<project>/islet-explorer}"
release="${RELEASE:?set RELEASE to a verified release identifier}"
port="${PORT:-3838}"
container_name="${CONTAINER_NAME:-islet-explorer-smoke}"
# A smoke container may run beside the production and staging instances, so
# its default cap is the staging one.
memory="${MEMORY:-4g}"
[[ "${memory}" =~ ^[1-9][0-9]*[mg]$ ]] || {
  echo "ERROR: MEMORY must be a size such as 4g or 4096m." >&2
  exit 64
}

[[ "${project_root}" == /* && "${project_root}" =~ ^[A-Za-z0-9._/-]+$ ]] || {
  echo "ERROR: PROJECT_ROOT must be an absolute path without whitespace." >&2
  exit 64
}
[[ "${release}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ && "${port}" =~ ^[0-9]{1,5}$ ]] || {
  echo "ERROR: invalid RELEASE or PORT." >&2
  exit 64
}
port="$((10#${port}))"
((port >= 1 && port <= 65535)) || {
  echo "ERROR: PORT must be between 1 and 65535." >&2
  exit 64
}

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
release_dir="${project_root}/releases/${release}"
"${script_dir}/verify-release.sh" "${release_dir}"

read_meta() {
  local key="$1"
  awk -F= -v key="${key}" '$1 == key {sub(/^[^=]*=/, ""); print; exit}' "${release_dir}/release.env"
}
image_ref="$(read_meta IMAGE_REF)"
data_schema="$(read_meta DATA_SCHEMA)"
[[ "${image_ref}" =~ ^sha256:[0-9a-f]{64}$ ||
   "${image_ref}" =~ ^[A-Za-z0-9._:/-]+@sha256:[0-9a-f]{64}$ ]] || {
  echo "ERROR: release image is not immutable." >&2
  exit 65
}
podman image exists "${image_ref}" || {
  echo "ERROR: image is not present locally; pull the immutable digest first: ${image_ref}" >&2
  exit 66
}
if podman container exists "${container_name}"; then
  echo "ERROR: container already exists: ${container_name}" >&2
  exit 65
fi

container_uid="$(podman run --rm --entrypoint /usr/bin/id "${image_ref}" -u shiny)"
container_gid="$(podman run --rm --entrypoint /usr/bin/id "${image_ref}" -g shiny)"
log_dir="${project_root}/runtime/logs-smoke"
tmp_dir="${project_root}/runtime/tmp-smoke"
mkdir -p -- "${log_dir}" "${tmp_dir}"
chmod 0750 "${log_dir}" "${tmp_dir}"

# As in activate-release.sh: the CPU quota only where the user session has the
# cpu controller delegated (PubApps VMs delegate memory and pids only).
cpu_limit=()
controllers_file="/sys/fs/cgroup/user.slice/user-$(id -u).slice/user@$(id -u).service/cgroup.controllers"
if [[ -r "${controllers_file}" ]] && grep -qw cpu "${controllers_file}"; then
  cpu_limit=(--cpus=2)
else
  echo "NOTE: the cpu controller is not delegated to this user session; the container keeps the host's vCPU bound."
fi

# OCI images do not retain the Containerfile HEALTHCHECK. Set it here, as the
# production Quadlet does, so verification can observe a healthy container.
container_id="$(podman run --detach --rm \
  --name "${container_name}" \
  --pull=never \
  --userns="keep-id:uid=${container_uid},gid=${container_gid}" \
  --user "${container_uid}:${container_gid}" \
  --publish "127.0.0.1:${port}:3838" \
  --read-only \
  --cap-drop=ALL \
  --security-opt=no-new-privileges \
  "${cpu_limit[@]}" --memory="${memory}" --pids-limit=2048 \
  --health-cmd='curl --fail --silent --show-error --max-time 8 http://127.0.0.1:3838/islet-explorer/' \
  --health-interval=30s --health-timeout=10s \
  --health-start-period=180s --health-retries=5 \
  --volume "${release_dir}/data:/srv/islet-explorer/data:ro" \
  --volume "${release_dir}/data/images:/srv/islet-explorer/app/shiny_app/www/local_images:ro" \
  --volume "${log_dir}:/var/log/shiny-server:rw" \
  --volume "${tmp_dir}:/tmp/islet-explorer:rw" \
  --tmpfs "/run:rw,size=16m,mode=0755,U" \
  --tmpfs "/var/lib/shiny-server:rw,size=64m,mode=0755,U" \
  --env "ISLET_RELEASE=${release}" \
  --env ISLET_RUNTIME_MODE=production \
  --env ISLET_DATA_ROOT=/srv/islet-explorer/data \
  --env "ISLET_DATA_SCHEMA=${data_schema}" \
  --env LOCAL_IMAGE_ROOT=/srv/islet-explorer/app/shiny_app/www/local_images \
  --env RETICULATE_PYTHON=/opt/islet-explorer-py/bin/python \
  --env ISLET_ENABLE_AI=FALSE \
  --env SHOW_STATISTICS_TAB=0 \
  --env SHOW_INFLAMMATION_TAB=0 \
  --env DEBUG_CREDENTIALS=0 \
  --env VIEWER_DEBUG=0 \
  --env TMPDIR=/tmp/islet-explorer \
  --env R_USER=/tmp/islet-explorer/r-user \
  --env XDG_CACHE_HOME=/tmp/islet-explorer/cache \
  --env MPLCONFIGDIR=/tmp/islet-explorer/matplotlib \
  --env NUMBA_CACHE_DIR=/tmp/islet-explorer/numba \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  --env NUMEXPR_NUM_THREADS=1 \
  "${image_ref}")"

echo "Started ${container_name}: ${container_id}"
if ! "${script_dir}/verify-deployment.sh" \
  --url "http://127.0.0.1:${port}/islet-explorer/" --container "${container_name}" \
  --expected-image "${image_ref}"; then
  podman stop "${container_name}" >/dev/null || true
  exit 70
fi
echo "Stop with: podman stop ${container_name}"
