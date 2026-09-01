#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: deploy/generate-release-manifest.sh RELEASE_DIR RELEASE IMAGE_REF GIT_COMMIT [DATA_SCHEMA]

RELEASE_DIR must contain the generated data bundle at data/, including its
data/images/ subtree. Optional attestations/ files are also included. The
command writes release.env and MANIFEST.sha256, invalidating any prior
verification stamp. IMAGE_REF must be an immutable registry digest or local
image ID.
EOF
}

if (($# < 4 || $# > 5)); then
  usage >&2
  exit 64
fi

release_dir="$(realpath -e -- "$1")"
release="$2"
image_ref="$3"
git_commit="$4"
data_schema="${5:-1}"

[[ "${release}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
  echo "ERROR: invalid release identifier: ${release}" >&2
  exit 64
}
[[ "${image_ref}" =~ ^sha256:[0-9a-f]{64}$ ||
   "${image_ref}" =~ ^[A-Za-z0-9._:/-]+@sha256:[0-9a-f]{64}$ ]] || {
  echo "ERROR: IMAGE_REF must be immutable (@sha256 digest or local sha256 image ID)." >&2
  exit 64
}
[[ "${git_commit}" =~ ^[0-9a-fA-F]{40}$ ]] || {
  echo "ERROR: invalid image reference or 40-character Git commit." >&2
  exit 64
}
[[ "${data_schema}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
  echo "ERROR: invalid data schema identifier." >&2
  exit 64
}
for required_dir in data data/images; do
  [[ -d "${release_dir}/${required_dir}" && ! -L "${release_dir}/${required_dir}" ]] || {
    echo "ERROR: ${release_dir}/${required_dir} must be a real directory." >&2
    exit 66
  }
done
[[ -f "${release_dir}/data/manifest.json" ]] || {
  echo "ERROR: generated data bundle is missing data/manifest.json." >&2
  exit 66
}

unsafe="$(find "${release_dir}" -xdev \( \
  -name '.Renviron' -o -name '.Rhistory' -o -name '*.bak' -o \
  -name '*~' -o -name '*.pem' -o -name '*.key' \
\) -print -quit)"
if [[ -n "${unsafe}" ]]; then
  echo "ERROR: release contains a forbidden credential/backup filename: ${unsafe}" >&2
  exit 65
fi

# Symlinks are allowed for TIFF offset aliases, but none may escape the release.
while IFS= read -r -d '' link_path; do
  resolved="$(realpath -e -- "${link_path}")" || {
    echo "ERROR: broken release symlink: ${link_path}" >&2
    exit 65
  }
  case "${resolved}" in
    "${release_dir}"/*) ;;
    *) echo "ERROR: release symlink escapes its release: ${link_path} -> ${resolved}" >&2; exit 65 ;;
  esac
done < <(find "${release_dir}/data" -type l -print0)

data_files="$(find "${release_dir}/data" -path "${release_dir}/data/images" -prune -o -type f -print | wc -l)"
image_files="$(find "${release_dir}/data/images" -type f | wc -l)"
data_total_bytes="$(du -sb -- "${release_dir}/data" | awk '{print $1}')"
image_bytes="$(du -sb -- "${release_dir}/data/images" | awk '{print $1}')"
data_bytes="$((data_total_bytes - image_bytes))"

umask 022
env_tmp="$(mktemp "${release_dir}/.release.env.XXXXXX")"
manifest_tmp="$(mktemp "${release_dir}/.MANIFEST.sha256.XXXXXX")"
trap 'rm -f -- "${env_tmp}" "${manifest_tmp}"' EXIT

{
  printf 'MANIFEST_VERSION=1\n'
  printf 'RELEASE=%s\n' "${release}"
  printf 'IMAGE_REF=%s\n' "${image_ref}"
  printf 'GIT_COMMIT=%s\n' "${git_commit,,}"
  printf 'DATA_SCHEMA=%s\n' "${data_schema}"
  printf 'CREATED_AT=%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
  printf 'DATA_FILES=%s\n' "${data_files}"
  printf 'DATA_BYTES=%s\n' "${data_bytes}"
  printf 'IMAGE_FILES=%s\n' "${image_files}"
  printf 'IMAGE_BYTES=%s\n' "${image_bytes}"
} > "${env_tmp}"
mv -f -- "${env_tmp}" "${release_dir}/release.env"

paths=(data release.env)
[[ -d "${release_dir}/attestations" ]] && paths+=(attestations)
(
  cd -- "${release_dir}"
  find "${paths[@]}" \( -type f -o -type l \) -print0 | \
    LC_ALL=C sort -z | xargs -0 -r sha256sum
) > "${manifest_tmp}"
mv -f -- "${manifest_tmp}" "${release_dir}/MANIFEST.sha256"
rm -f -- "${release_dir}/.verified"
trap - EXIT

echo "Generated ${release_dir}/MANIFEST.sha256"
echo "Run deploy/verify-release.sh '${release_dir}' before activation."
