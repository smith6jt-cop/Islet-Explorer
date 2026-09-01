#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: deploy/verify-deployment.sh --url URL [--container NAME]
       [--range-path RELATIVE_OME_TIFF_PATH]

Checks container health, loopback-only publication, application HTTP response,
production flags, dependency availability, read-only scientific mounts, and an
optional OME-TIFF byte-range response. Run the external Basic-Auth checks from a
non-UF network as described in deploy/README.md.
EOF
}

url=""
container_name="islet-explorer"
range_path=""
while (($#)); do
  case "$1" in
    --url) url="${2:?missing value for --url}"; shift 2 ;;
    --container) container_name="${2:?missing value for --container}"; shift 2 ;;
    --range-path) range_path="${2:?missing value for --range-path}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 64 ;;
  esac
done
[[ "${url}" == http://127.0.0.1:*/* || "${url}" == http://localhost:*/* ]] || {
  echo "ERROR: this verifier expects the local loopback backend URL." >&2
  exit 64
}

for command_name in podman curl ss; do
  command -v "${command_name}" >/dev/null || {
    echo "ERROR: required command not found: ${command_name}" >&2
    exit 69
  }
done

deadline=$((SECONDS + 210))
until curl --fail --location --silent --show-error --max-time 15 "${url}" >/dev/null; do
  if ((SECONDS >= deadline)); then
    echo "ERROR: app did not become ready within 210 seconds: ${url}" >&2
    podman logs --tail 100 "${container_name}" >&2 || true
    exit 70
  fi
  sleep 3
done

while :; do
  health="$(podman inspect --format '{{if .State.Health}}{{.State.Health.Status}}{{else}}missing{{end}}' "${container_name}")"
  [[ "${health}" == healthy ]] && break
  if [[ "${health}" == unhealthy || "${health}" == missing ]] || ((SECONDS >= deadline)); then
    echo "ERROR: container did not reach healthy state (health=${health})." >&2
    podman logs --tail 100 "${container_name}" >&2 || true
    exit 70
  fi
  sleep 3
done

running="$(podman inspect --format '{{.State.Running}}' "${container_name}")"
[[ "${running}" == true && "${health}" == healthy ]] || {
  echo "ERROR: container state running=${running}, health=${health}." >&2
  exit 70
}

if ss -ltnH | awk '$4 ~ /(^|\])0\.0\.0\.0:3838$|^\*:3838$|^\[::\]:3838$|^:::3838$/ {found=1} END {exit !found}'; then
  echo "ERROR: backend port 3838 is exposed on a non-loopback host address." >&2
  exit 77
fi
ss -ltnH | awk '$4 ~ /127\.0\.0\.1:3838$/ {found=1} END {exit !found}' || {
  echo "ERROR: expected 127.0.0.1:3838 listener was not found." >&2
  exit 77
}

podman exec "${container_name}" /bin/sh -ceu '
  test "$(id -u)" = "$(id -u shiny)"
  test "$(id -g)" = "$(id -g shiny)"
  test "${ISLET_RUNTIME_MODE}" = production
  test "${ISLET_ENABLE_AI}" = FALSE
  test "${SHOW_STATISTICS_TAB}" = 0
  test "${SHOW_INFLAMMATION_TAB}" = 0
  test -z "${KEY:-}"
  test -z "${OPENAI_API_KEY:-}"
  test -r /srv/islet-explorer/data
  test ! -w /srv/islet-explorer/data
  test -r /srv/islet-explorer/app/shiny_app/www/local_images
  test ! -w /srv/islet-explorer/app/shiny_app/www/local_images
  test -w /tmp/islet-explorer
  test -w /var/log/shiny-server
'
podman exec "${container_name}" /opt/islet-explorer-py/bin/python -c \
  'import anndata, scanpy, igraph, leidenalg, fa2_modified'
podman exec "${container_name}" /usr/bin/Rscript -e \
  'stopifnot(packageVersion("destiny")=="3.22.0", packageVersion("duckdb")=="1.5.2")'

if [[ -n "${range_path}" ]]; then
  range_url="${url%/}/${range_path#/}"
  headers="$(mktemp)"
  body="$(mktemp)"
  trap 'rm -f -- "${headers}" "${body}"' EXIT
  curl --silent --show-error --max-time 30 \
    --range 0-15 --dump-header "${headers}" --output "${body}" "${range_url}"
  status="$(awk 'toupper($1) ~ /^HTTP\// {code=$2} END {print code}' "${headers}")"
  content_range="$(awk 'tolower($1) == "content-range:" {gsub("\r", ""); sub(/^[^:]+:[[:space:]]*/, ""); print; exit}' "${headers}")"
  bytes="$(wc -c < "${body}" | tr -d '[:space:]')"
  [[ "${status}" == 206 && "${content_range}" == bytes\ 0-15/* && "${bytes}" == 16 ]] || {
    echo "ERROR: invalid Range response: status=${status}, range=${content_range}, bytes=${bytes}." >&2
    exit 70
  }
  echo "OME-TIFF Range response verified: ${content_range}."
fi

echo "Deployment verification passed for ${url} (health=${health}, loopback-only backend)."
