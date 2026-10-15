#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd -- "${script_dir}/../.." && pwd -P)"
work_dir="$(mktemp -d)"
trap 'rm -rf -- "${work_dir}"' EXIT
mock_bin="${work_dir}/bin"
mkdir -p "${mock_bin}"
image_id="sha256:$(printf 'a%.0s' {1..64})"
other_image="sha256:$(printf 'b%.0s' {1..64})"

cat > "${mock_bin}/podman" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
case " $* " in
  *" /usr/bin/Rscript "*) echo 'Rscript must resolve through PATH' >&2; exit 70 ;;
  *" {{.Id}} "*) printf '%s\n' "${MOCK_IMAGE_ID}" ;;
  *" {{.Image}} "*) printf '%s\n' "${MOCK_ACTUAL_IMAGE:-${MOCK_IMAGE_ID}}" ;;
  *" .State.Health"*|*"{{if .State.Health}}"*) echo healthy ;;
  *" {{.State.Running}} "*) echo true ;;
  *" exec "*) : ;;
  *) echo "unexpected Podman invocation: $*" >&2; exit 64 ;;
esac
MOCK
cat > "${mock_bin}/ss" <<'MOCK'
#!/usr/bin/env bash
cat "${MOCK_LISTENERS}"
MOCK
cat > "${mock_bin}/curl" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
headers=""
body=""
while (($#)); do
  case "$1" in
    --dump-header) headers="$2"; shift 2 ;;
    --output) body="$2"; shift 2 ;;
    *) shift ;;
  esac
done
if [[ -n "${headers}" ]]; then
  printf 'HTTP/1.1 %s Partial Content\r\nContent-Range: bytes 0-15/1000\r\n\r\n' \
    "${MOCK_RANGE_STATUS:-206}" > "${headers}"
  printf '0123456789abcdef' > "${body}"
fi
MOCK
chmod 0755 "${mock_bin}/podman" "${mock_bin}/ss" "${mock_bin}/curl"

check() {
  local expected="$1"
  local status
  shift
  set +e
  env PATH="${mock_bin}:${PATH}" MOCK_IMAGE_ID="${image_id}" \
    MOCK_LISTENERS="${work_dir}/listeners" "$@" \
    > "${work_dir}/output" 2>&1
  status=$?
  set -e
  if [[ "${status}" != "${expected}" ]]; then
    cat "${work_dir}/output" >&2
    echo "FAIL: expected exit ${expected}, got ${status}: $*" >&2
    exit 1
  fi
}

verifier=(bash "${repo_root}/deploy/verify-deployment.sh")
url='http://127.0.0.1:3838/islet-explorer/'
printf 'LISTEN 0 128 127.0.0.1:3838 0.0.0.0:*\n' > "${work_dir}/listeners"
check 0 "${verifier[@]}" --url "${url}" --expected-image "${image_id}" --range-path local_images/test.ome.tiff
check 0 env MOCK_IMAGE_ID="${image_id#sha256:}" MOCK_ACTUAL_IMAGE="${image_id}" \
  "${verifier[@]}" --url "${url}" --expected-image "${image_id}"
check 0 env MOCK_IMAGE_ID="${image_id}" MOCK_ACTUAL_IMAGE="${image_id#sha256:}" \
  "${verifier[@]}" --url "${url}" --expected-image "${image_id}"
check 0 env MOCK_IMAGE_ID="${image_id#sha256:}" MOCK_ACTUAL_IMAGE="${image_id#sha256:}" \
  "${verifier[@]}" --url "${url}" --expected-image "${image_id}"
for invalid_id in latest "${image_id:1}" "sha256:${image_id}" " ${image_id}"; do
  check 70 env MOCK_IMAGE_ID="${invalid_id}" MOCK_ACTUAL_IMAGE="${invalid_id}" \
    "${verifier[@]}" --url "${url}" --expected-image "${image_id}"
  check 70 env MOCK_ACTUAL_IMAGE="${invalid_id}" \
    "${verifier[@]}" --url "${url}" --expected-image "${image_id}"
done
check 70 env MOCK_ACTUAL_IMAGE="${other_image}" "${verifier[@]}" --url "${url}" --expected-image "${image_id}"
check 70 env MOCK_RANGE_STATUS=200 "${verifier[@]}" --url "${url}" --range-path local_images/test.ome.tiff

# A correct loopback binding must not conceal an additional public/private IP
# listener. Test both wildcard families and specifically bound host addresses.
for address in '0.0.0.0' '*' '[::]' '::' '10.20.30.40' '[2001:db8::1]' '[::1]'; do
  printf 'LISTEN 0 128 127.0.0.1:3838 0.0.0.0:*\nLISTEN 0 128 %s:3838 *:*\n' \
    "${address}" > "${work_dir}/listeners"
  check 77 "${verifier[@]}" --url "${url}"
done
printf 'LISTEN 0 128 127.0.0.1:38380 0.0.0.0:*\n' > "${work_dir}/listeners"
check 77 "${verifier[@]}" --url "${url}"
printf 'LISTEN 0 128 127.0.0.1:3939 0.0.0.0:*\n' > "${work_dir}/listeners"
check 0 "${verifier[@]}" --url http://localhost:3939/islet-explorer/
check 64 "${verifier[@]}" --url http://10.20.30.40:3939/islet-explorer/
check 64 "${verifier[@]}" --url http://127.0.0.1:65536/islet-explorer/
check 64 "${verifier[@]}" --url http://127.0.0.1:0/islet-explorer/
check 64 "${verifier[@]}" --url http://127.0.0.1:3939/other-app/
check 64 "${verifier[@]}" --url http://127.0.0.1:3939/islet-explorer/ --expected-image latest

echo 'Deployment listener, image identity, and Range verifier tests passed.'
