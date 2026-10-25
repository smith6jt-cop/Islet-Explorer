#!/usr/bin/env bash
# Build the reviewed upstream server in an isolated directory. This script never
# installs into the host. Containerfile copies only its staged runtime output.
set -euo pipefail

output_dir=${1:?Usage: harden-base.sh OUTPUT_DIR EVIDENCE_DIR}
evidence_dir=${2:?Usage: harden-base.sh OUTPUT_DIR EVIDENCE_DIR}
[[ "$output_dir" = /* && "$evidence_dir" = /* && ! -e "$output_dir" ]] || {
  echo 'Use absolute paths and a new output directory.' >&2
  exit 1
}
[[ "$(uname -m)" = x86_64 ]] || {
  echo 'The reviewed Node archive is locked to Linux x86_64.' >&2
  exit 1
}

source_revision=ae268b685b234b3a97a33bf6caa28e6d395cb232
source_sha256=4a1d74aa4e2c266f4aaa04198132a017cb8b1e07ef766d163f7dec5d9dc17a62
source_url=https://api.github.com/repos/rstudio/shiny-server/tarball/$source_revision
# Upstream's dependency files at that revision, and the project's copies with
# four security pins (2026-10-22): handlebars 4.7.10 and compression 1.8.2 as
# exact direct dependencies, proxy-addr 2.0.8 and follow-redirects 1.16.1 as
# overrides. The pinned lock was regenerated from upstream's with the same npm
# (`npm install --package-lock-only`) and differs from it only in those four
# packages and compression's new dependency `destroy`; see BASE_HARDENING.md.
# 2026-10-24 adds one runtime dependency, permessage-deflate 0.1.7 (below).
# 2026-10-25 adds a fifth pin, body-parser 2.4.0 as an override (GHSA-27r3-jxmh-p828),
# whose dependency range moves qs to 6.16.1; see BASE_HARDENING.md.
upstream_package_sha256=191d5ee18b8b1263a3fff0109634e55f2bec78e13f3db004b751c8c0342bb714
upstream_shrinkwrap_sha256=914bc927be865cb448c2eceff03fad641b93bea962c2bbacc485a14309a6b4fd
pinned_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)/shiny-server
pinned_package_sha256=1b57c364d942417f015578676f2490c9c4bff72331a3747a5e456e4b14c36892
pinned_shrinkwrap_sha256=d3ddc237a34ab5bd141881a65e807dbe8a738accfc62041d3285209448a294a9
# 2026-10-24: the WebSocket transport offers permessage-deflate (the app's figure messages
# shrink about five times on the wire, see BASE_HARDENING.md). The pinned manifest adds
# permessage-deflate 0.1.7 (its one dependency, safe-buffer, is already locked), and
# lib/proxy/sockjs.js gets a six-line option block passed through SockJS to faye-websocket.
# Both the upstream file and the patched result are hash-checked.
upstream_sockjs_sha256=f5bc3729f08b2110f34d0196cb2e8b27a305792c4f7e255ac23e60d04795a253
patched_sockjs_sha256=9c671905f16d61cb1b93c71717b2762266d0f08cf31e2cffc5635b01f7f2af32
node_version=v24.20.0
node_sha256=2f2c0da162318f0de47665410c7c8c2ed3d36c8f3105de4bbc61176c70a7cbf2
node_url=https://nodejs.org/dist/$node_version/node-$node_version-linux-x64.tar.xz
work_dir=$(mktemp -d /tmp/islet-shiny-build.XXXXXXXX)
trap 'rm -rf "$work_dir"' EXIT
mkdir -p "$evidence_dir" "$work_dir/source"

curl --fail --location --retry 3 --proto '=https' --tlsv1.2 \
  --output "$work_dir/source.tar.gz" "$source_url"
curl --fail --location --retry 3 --proto '=https' --tlsv1.2 \
  --output "$work_dir/node.tar.xz" "$node_url"
printf '%s  %s\n' "$source_sha256" "$work_dir/source.tar.gz" \
  "$node_sha256" "$work_dir/node.tar.xz" | sha256sum --check --strict
tar -xzf "$work_dir/source.tar.gz" --strip-components=1 -C "$work_dir/source"
cd "$work_dir/source"
[[ "$(cat .nvmrc)" = "$node_version" ]]

# Offer permessage-deflate on the SockJS WebSocket transport (2026-10-24). SockJS hands
# `faye_server_options` to faye-websocket, whose `extensions` option negotiates the
# compression with browsers that request it; everything else in the file stays as upstream.
printf '%s  %s\n' "$upstream_sockjs_sha256" lib/proxy/sockjs.js | sha256sum --check --strict
python3 - lib/proxy/sockjs.js <<'PY'
import pathlib, sys
path = pathlib.Path(sys.argv[1])
source = path.read_text()
old = """    heartbeat_delay: heartbeatDelay,
    disconnect_delay: disconnectDelay
  });
"""
new = """    heartbeat_delay: heartbeatDelay,
    disconnect_delay: disconnectDelay,
    // Islet Explorer (deploy/BASE_HARDENING.md, 2026-10-24): offer permessage-deflate on the
    // WebSocket transport. faye-websocket negotiates it with browsers that ask for it, which
    // shrinks the app's figure messages about five times for reviewers on slow links.
    faye_server_options: {
      extensions: [require('permessage-deflate')]
    }
  });
"""
assert source.count(old) == 1, "sockjs.createServer block not found exactly once"
path.write_text(source.replace(old, new))
PY
printf '%s  %s\n' "$patched_sockjs_sha256" lib/proxy/sockjs.js | sha256sum --check --strict
mkdir -p ext/node
tar -xJf "$work_dir/node.tar.xz" --strip-components=1 -C ext/node
cp ext/node/bin/node ext/node/bin/shiny-server

export PATH="$PWD/ext/node/bin:$PATH"
export npm_config_cache="$work_dir/npm-cache"
export npm_config_nodedir="$PWD/ext/node"
export npm_config_update_notifier=false
# Upstream shrinkwrap locks three GitHub dependencies by commit. Keep the lock
# byte-for-byte and use HTTPS transport without requiring an SSH credential.
export GIT_CONFIG_COUNT=1
export GIT_CONFIG_KEY_0=url.https://github.com/.insteadOf
export GIT_CONFIG_VALUE_0=ssh://git@github.com/
# Replace upstream's dependency files with the pinned copies, checking both.
printf '%s  %s\n' "$upstream_package_sha256" package.json \
  "$upstream_shrinkwrap_sha256" npm-shrinkwrap.json \
  "$pinned_package_sha256" "$pinned_dir/package.json" \
  "$pinned_shrinkwrap_sha256" "$pinned_dir/npm-shrinkwrap.json" | sha256sum --check --strict
cp "$pinned_dir/package.json" package.json
cp "$pinned_dir/npm-shrinkwrap.json" npm-shrinkwrap.json
sha256sum npm-shrinkwrap.json > "$work_dir/shrinkwrap.sha256"
npm ci --no-audit --no-fund > "$evidence_dir/npm-ci.log" 2>&1
npm run build > "$evidence_dir/typescript-build.log" 2>&1
sha256sum --check --strict "$work_dir/shrinkwrap.sha256"

# Invoke precisely the upstream test directories, using JSON so skipped real-R
# tests cannot silently count as successful release evidence.
node node_modules/mocha/bin/mocha.js test test/integration --reporter json \
  --reporter-option "output=$evidence_dir/upstream-tests.json" \
  > "$evidence_dir/upstream-tests.stdout" 2> "$evidence_dir/upstream-tests.stderr"
node node_modules/mocha/bin/mocha.js test/integration-r --reporter json \
  --reporter-option "output=$evidence_dir/upstream-r-tests.json" \
  > "$evidence_dir/upstream-r-tests.stdout" 2> "$evidence_dir/upstream-r-tests.stderr"
python3 - "$evidence_dir" <<'PY'
import json, pathlib, sys
root = pathlib.Path(sys.argv[1])
for name, minimum in [('upstream-tests.json', 300), ('upstream-r-tests.json', 8)]:
    stats = json.loads((root / name).read_text())['stats']
    assert stats['passes'] >= minimum and stats['failures'] == 0 and stats['pending'] == 0, stats
    print(f"{name}: {stats['passes']} passed, no failures or skips")
PY

# The upstream CMake launcher target is this single C++ source and generated
# constant. Building it directly avoids CMake's unrelated, obsolete Pandoc 1.19
# download. The runtime uses the base image's maintained /usr/bin/pandoc.
printf '#define SHINY_SERVER_DEFAULT_BIN_PATH "/opt/shiny-server/bin/shiny-server"\n' > src/launcher.h
c++ -std=c++11 -O2 -o bin/shiny-server src/launcher.cc
# Reinstall just the locked runtime graph. npm prune rewrites shrinkwrap
# metadata; npm ci preserves it. The already-tested native addon is outside
# node_modules and is copied below, so no lifecycle rebuild is needed here.
npm ci --omit=dev --ignore-scripts --no-audit --no-fund \
  > "$evidence_dir/production-npm-ci.log" 2>&1
sha256sum --check --strict "$work_dir/shrinkwrap.sha256"

mkdir -p "$output_dir/bin" "$output_dir/build/Release" \
  "$output_dir/ext/node/bin" "$output_dir/ext/pandoc"
cp -a assets config lib node_modules python R scripts templates "$output_dir/"
cp -a bin/node bin/shiny-server "$output_dir/bin/"
cp build/Release/posix.node "$output_dir/build/Release/"
cp ext/node/bin/node ext/node/bin/shiny-server "$output_dir/ext/node/bin/"
cp ext/node/LICENSE "$output_dir/ext/node/LICENSE"
cp package.json npm-shrinkwrap.json COPYING NEWS README.md "$output_dir/"
cp NOTICE.md "$output_dir/NOTICE"
printf '1.5.24.islet.%s\n' "${source_revision:0:12}" > "$output_dir/VERSION"
printf '%s\n' "$source_revision" > "$output_dir/GIT_VERSION"
ln -s /usr/bin/pandoc "$output_dir/ext/pandoc/pandoc"
# npm's CLI, native build scratch, headers, tests and development dependencies
# are absent from the runtime. Installed package manifests and the unchanged
# upstream lock remain available for the image SBOM and source audit.
[[ ! -e "$output_dir/ext/node/lib/node_modules/npm" ]]
[[ ! -e "$output_dir/node_modules/mocha" ]]
"$output_dir/bin/shiny-server" --version > "$evidence_dir/staged-version.txt"
python3 - "$evidence_dir" "$source_revision" "$source_url" "$source_sha256" \
  "$node_version" "$node_url" "$node_sha256" \
  "$upstream_package_sha256" "$upstream_shrinkwrap_sha256" \
  "$pinned_package_sha256" "$pinned_shrinkwrap_sha256" <<'PY'
import json, pathlib, sys
(dest, revision, source_url, source_hash, node_version, node_url, node_hash,
 upstream_package, upstream_shrinkwrap, pinned_package, pinned_shrinkwrap) = sys.argv[1:]
pathlib.Path(dest, 'build-provenance.json').write_text(json.dumps({
    'upstream_revision': revision,
    'source_url': source_url, 'source_sha256': source_hash,
    'node_version': node_version, 'node_url': node_url, 'node_sha256': node_hash,
    'dependency_policy': ('upstream npm-shrinkwrap.json regenerated with five security pins '
                          '(handlebars 4.7.10, compression 1.8.2, proxy-addr 2.0.8, follow-redirects 1.16.1, '
                          'body-parser 2.4.0 with qs 6.16.1) '
                          'plus permessage-deflate 0.1.7 for WebSocket compression (lib/proxy/sockjs.js '
                          'patched, both hashes checked); npm ci; omit dev at runtime'),
    'upstream_package_json_sha256': upstream_package,
    'upstream_npm_shrinkwrap_sha256': upstream_shrinkwrap,
    'pinned_package_json_sha256': pinned_package,
    'pinned_npm_shrinkwrap_sha256': pinned_shrinkwrap,
    'status': 'staged runtime; exact application image scan and browser gate still required'
}, indent=2) + '\n')
PY
echo "Shiny Server runtime staged at $output_dir"
