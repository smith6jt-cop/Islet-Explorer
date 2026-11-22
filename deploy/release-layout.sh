#!/usr/bin/env bash
# shellcheck disable=SC2034
# Shared adapter for independently verified app/runtime/data artifacts. Metadata
# is parsed as values, never evaluated as shell code. The assigned variables
# and arrays form this sourced library's public interface to both launchers.
load_split_layout() {
  local resolved key value layout_dir
  layout_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
  resolved="$(python3 "${layout_dir}/release-artifacts.py" "$@")" || return $?
  app_dir=""; app_scripts_dir=""; app_manifest_sha256=""; descriptor_sha256=""
  runtime_base_image=""
  while IFS='=' read -r key value; do
    case "${key}" in
      IMAGE_REF) image_ref="${value}" ;;
      RUNTIME_COMMIT) runtime_commit="${value}" ;;
      APP_COMMIT) git_commit="${value}" ;;
      APP_DIR) app_dir="${value}" ;;
      APP_SCRIPTS_DIR) app_scripts_dir="${value}" ;;
      DATA_DIR) data_dir="${value}" ;;
      DATA_SCHEMA) data_schema="${value}" ;;
      DATA_COMMIT) data_commit="${value}" ;;
      DATA_ID) data_id="${value}" ;;
      DATA_RELEASE_DIR) data_release_dir="${value}" ;;
      DATA_MANIFEST) data_manifest="${value}" ;;
      APP_MANIFEST_SHA256) app_manifest_sha256="${value}" ;;
      RUNTIME_INPUT_SHA256) runtime_input_sha256="${value}" ;;
      RUNTIME_BASE_IMAGE) runtime_base_image="${value}" ;;
      DESCRIPTOR_SHA256) descriptor_sha256="${value}" ;;
    esac
  done <<< "${resolved}"
  [[ -n "${app_dir}" && -n "${app_scripts_dir}" && -n "${app_manifest_sha256}" && -n "${runtime_base_image}" ]] || return 65
  app_mounts=(--volume "${app_dir}:/srv/islet-explorer/app:ro"
              --volume "${app_scripts_dir}:/srv/islet-explorer/scripts:ro")
  app_checks=(--expected-app "${app_dir}" --expected-data "${data_dir}"
              --expected-app-sha256 "${app_manifest_sha256}")
}

check_split_runtime() {
  local layout inputs revision base_image
  layout="$(podman image inspect --format '{{ index .Labels "io.islet-explorer.runtime-layout" }}' "${image_ref}")" || return $?
  inputs="$(podman image inspect --format '{{ index .Labels "io.islet-explorer.runtime-inputs-sha256" }}' "${image_ref}")" || return $?
  revision="$(podman image inspect --format '{{ index .Labels "org.opencontainers.image.revision" }}' "${image_ref}")" || return $?
  base_image="$(podman image inspect --format '{{ index .Labels "org.opencontainers.image.base.name" }}' "${image_ref}")" || return $?
  [[ "${layout}" == split-v1 && "${inputs}" == "${runtime_input_sha256}" && "${revision}" == "${runtime_commit}" && "${base_image}" == "${runtime_base_image}" ]] || {
    echo 'ERROR: runtime layout, build inputs, base image or build revision differ from the verified runtime artifact.' >&2
    return 65
  }
}
