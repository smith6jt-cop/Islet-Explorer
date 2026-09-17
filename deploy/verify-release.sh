#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: deploy/verify-release.sh RELEASE_DIR [--quick] [--expected-images N]" >&2
}

if (($# < 1)); then usage; exit 64; fi
release_dir="$(realpath -e -- "$1")"
shift
quick=0
requested_images=""
while (($#)); do
  case "$1" in
    --quick) quick=1; shift ;;
    --expected-images)
      [[ $# -ge 2 && -z "${requested_images}" && "${2}" =~ ^[1-9][0-9]*$ ]] || {
        echo "ERROR: supply --expected-images once with a positive integer without leading zeros." >&2
        exit 64
      }
      requested_images="$2"; shift 2 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 64 ;;
  esac
done

for required in data data/images data/manifest.json release.env MANIFEST.sha256; do
  [[ -e "${release_dir}/${required}" ]] || {
    echo "ERROR: release is missing ${required}." >&2
    exit 66
  }
done

read_meta() {
  local key="$1"
  awk -F= -v key="${key}" '
    $1 == key {sub(/^[^=]*=/, ""); value = $0; count++}
    END {if (count > 1) exit 1; if (count == 1) print value}
  ' "${release_dir}/release.env" || {
    echo "ERROR: duplicate ${key} in release.env." >&2
    return 65
  }
}

manifest_version="$(read_meta MANIFEST_VERSION)"
expected_images="$(read_meta EXPECTED_IMAGES)"
# Version 1 predates the explicit cohort count and is the frozen 15-image
# contract. A caller cannot override it to accept an expanded legacy bundle.
case "${manifest_version}" in
  1)
    [[ -z "${expected_images}" || "${expected_images}" == 15 ]] || {
      echo "ERROR: version 1 releases must use the original 15-image cohort." >&2
      exit 65
    }
    expected_images=15
    ;;
  2)
    [[ "${expected_images}" =~ ^[1-9][0-9]*$ ]] || {
      echo "ERROR: version 2 release.env requires a positive EXPECTED_IMAGES count." >&2
      exit 65
    }
    ;;
  *) echo "ERROR: unsupported or missing MANIFEST_VERSION in release.env." >&2; exit 65 ;;
esac
[[ -z "${requested_images}" || "${requested_images}" == "${expected_images}" ]] || {
  echo "ERROR: requested count ${requested_images} conflicts with release count ${expected_images}." >&2
  exit 65
}

release="$(read_meta RELEASE)"
image_ref="$(read_meta IMAGE_REF)"
git_commit="$(read_meta GIT_COMMIT)"
data_schema="$(read_meta DATA_SCHEMA)"
[[ "${release}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ &&
   ( "${image_ref}" =~ ^sha256:[0-9a-f]{64}$ ||
     "${image_ref}" =~ ^[A-Za-z0-9._:/-]+@sha256:[0-9a-f]{64}$ ) &&
   "${git_commit}" =~ ^[0-9a-f]{40}$ && -n "${data_schema}" ]] || {
  echo "ERROR: release.env is incomplete or malformed." >&2
  exit 65
}

# Metadata determines the image-count contract, including during quick checks.
# Require one hash entry and verify it instead of trusting a mutable release.env.
metadata_hash="$(awk '
  length($1) == 64 && $1 !~ /[^0-9a-f]/ && $2 == "release.env" && NF == 2 {
    hash = $1; count++
  }
  END {if (count != 1) exit 1; print hash}
' "${release_dir}/MANIFEST.sha256")" || {
  echo "ERROR: MANIFEST.sha256 must contain exactly one release.env checksum." >&2
  exit 65
}
[[ "$(sha256sum "${release_dir}/release.env" | awk '{print $1}')" == "${metadata_hash}" ]] || {
  echo "ERROR: release.env does not match its manifest checksum." >&2
  exit 65
}

actual_images="$(find "${release_dir}/data/images" -type f \
  \( -iname '*.ome.tif' -o -iname '*.ome.tiff' \) | wc -l)"
if [[ "${actual_images}" != "${expected_images}" ]]; then
  echo "ERROR: expected ${expected_images} OME-TIFF files, found ${actual_images}." >&2
  exit 65
fi
if [[ "${data_schema}" == 2 || "${data_schema}" == 3 ]]; then
  python3 - "${release_dir}/data/manifest.json" "${expected_images}" "${data_schema}" <<'PY'
import json
import sys
from pathlib import Path
with open(sys.argv[1]) as stream:
    manifest = json.load(stream)
schema = sys.argv[3]
required = {"cohort", "trajectory", "cells"}
if schema == "2":
    required.add("research")
if (str(manifest.get("schema_version")) != schema
        or manifest.get("expected_images") != int(sys.argv[2])
        or not required.issubset(manifest.get("artifacts", {}))):
    raise SystemExit("ERROR: expanded runtime requires matching cohort, trajectory and cell artifacts.")
if schema == "3":
    root = Path(sys.argv[1]).parent
    files = ["counts.csv", "endocrine_counts.csv", "phenotype_availability.csv", "marker_panel.json", "count_manifest.json"]
    if "research" in manifest["artifacts"] or not all((root/"trajectory"/f).is_file() for f in files):
        raise SystemExit("ERROR: schema-3 count artifacts are required; research reports must be outside the runtime.")
if manifest.get("scientific_acceptance") is not True:
    raise SystemExit("ERROR: local expanded candidate has not passed scientific deployment acceptance.")
PY
fi
[[ -n "$(find "${release_dir}/data" -type f -print -quit)" ]] || {
  echo "ERROR: release data directory is empty." >&2
  exit 65
}

manifest_hash="$(sha256sum "${release_dir}/MANIFEST.sha256" | awk '{print $1}')"
if ((quick)); then
  [[ -f "${release_dir}/.verified" ]] || {
    echo "ERROR: no full-verification stamp; run verify-release.sh without --quick." >&2
    exit 65
  }
  stamped_hash="$(awk -F= '$1 == "MANIFEST_SHA256" {print $2; exit}' "${release_dir}/.verified")"
  [[ "${stamped_hash}" == "${manifest_hash}" ]] || {
    echo "ERROR: release manifest changed after full verification." >&2
    exit 65
  }
  echo "Release ${release} has a current full-verification stamp."
  exit 0
fi

(
  cd -- "${release_dir}"
  sha256sum --check --strict MANIFEST.sha256
)

stamp_tmp="$(mktemp "${release_dir}/.verified.XXXXXX")"
trap 'rm -f -- "${stamp_tmp}"' EXIT
{
  printf 'MANIFEST_SHA256=%s\n' "${manifest_hash}"
  printf 'VERIFIED_AT=%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
  printf 'VERIFIED_HOST=%s\n' "$(hostname -f 2>/dev/null || hostname)"
  printf 'VERIFIED_UID=%s\n' "$(id -u)"
} > "${stamp_tmp}"
mv -f -- "${stamp_tmp}" "${release_dir}/.verified"
trap - EXIT

echo "Release ${release} verified: ${actual_images} OME-TIFFs; manifest ${manifest_hash}."
