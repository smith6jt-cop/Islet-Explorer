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
  "$node_version" "$node_url" "$node_sha256" <<'PY'
import json, pathlib, sys
dest, revision, source_url, source_hash, node_version, node_url, node_hash = sys.argv[1:]
pathlib.Path(dest, 'build-provenance.json').write_text(json.dumps({
    'upstream_revision': revision,
    'source_url': source_url, 'source_sha256': source_hash,
    'node_version': node_version, 'node_url': node_url, 'node_sha256': node_hash,
    'dependency_policy': 'unchanged upstream npm-shrinkwrap.json; npm ci; omit dev at runtime',
    'status': 'staged runtime; exact application image scan and browser gate still required'
}, indent=2) + '\n')
PY
echo "Shiny Server runtime staged at $output_dir"
