# Shiny Server base remediation

Review date: 2026-09-08, with the 2026-10-22 dependency pins below. This is a
tested source-build remediation; the exact application image and its
operating-system package upgrades still require the release build, SBOM,
vulnerability scan, and browser gates.

## 2026-10-22: four runtime dependency pins

The first workstation build of the schema-5 release image (commit `6ffb6f2`)
failed the High/Critical gate on advisories published after September against
the runtime graph of the upstream revision below: handlebars 4.7.9
(GHSA-8r5x-fm3f-whwj, GHSA-p8wg-vrv2-v86f, Critical), proxy-addr 2.0.7
(GHSA-jqcg-44mw-7w3h, Critical), compression 1.8.1 (GHSA-vc2v-76pw-4v95, High)
and follow-redirects 1.16.0 (GHSA-8r9p-f939-6h3c, GHSA-q3jh-wrc6-jgq8, High),
each fixed in the next patch release. Upstream `main` had no newer revision.
The same scan found pillow 12.3.0 and certifi 2026.7.22 in the application
venv; `requirements-python.lock` now pins pillow 13.0.0 and certifi
2026.10.19.

[`shiny-server/package.json`](shiny-server/package.json) and
[`shiny-server/npm-shrinkwrap.json`](shiny-server/npm-shrinkwrap.json) are
upstream's files with handlebars `4.7.10` and compression `1.8.2` as exact
direct dependencies and proxy-addr `2.0.8` and follow-redirects `1.16.1` as
`overrides` (npm refuses an override that conflicts with a direct spec, hence
the split). The lock was regenerated from upstream's with the locked Node
24.20.0's npm 11.19.0 (`npm install --package-lock-only --ignore-scripts`); it
differs from upstream's only in those four entries, compression 1.8.2's new
dependency `destroy 1.2.0`, and the root entry's `engines` field, which npm
refreshed from package.json. `harden-base.sh` checks the SHA-256 of upstream's
two files and of the pinned copies before substituting them, then builds,
tests and installs exactly as before; `build-provenance.json` records all four
hashes. Every other package stays at upstream's locked version.

## 2026-10-24: WebSocket compression (permessage-deflate)

The responsiveness baseline (`docs/performance_baseline_20261024.md`) showed that what an
off-campus reviewer waits for is the app's WebSocket payload over their link (3.5 MB for
the first screen, 1.7 MB per Plot change), and that Shiny Server's SockJS transport sends
it uncompressed. SockJS passes `faye_server_options` straight to `faye-websocket`, whose
`extensions` option negotiates permessage-deflate with browsers that offer it, so the
project build now adds one runtime dependency, `permessage-deflate 0.1.7` (Apache-2.0; its
only dependency, `safe-buffer`, was already locked at 5.1.2), and patches the six-line
option block into `lib/proxy/sockjs.js`'s `sockjs.createServer` call. `harden-base.sh`
checks the SHA-256 of the upstream file before and of the patched file after the
substitution, alongside the manifest hashes; `build-provenance.json` names the addition.
Measured on the app's own messages: the Plot cards compress 5.3 times, the Spatial tab
(a base64 PNG) 2.8 times. Per-message zlib on the VM costs tens of milliseconds; each
connection holds about 300 KB of compression state. Browsers that do not offer the
extension, and the proxy's health checks, are unaffected. The negotiated extension is
read from the handshake response by `deploy/browser/bench-app.mjs` (`socket.extensions`).

The inspected `rocker/shiny:4.5.1` base contained Shiny Server 1.5.24.1034,
Node 20.17.0, and npm 10.8.2. Its registry SBOM/Grype report recorded **33 High
and 3 Critical matches**. These were real package findings; none is suppressed.

| Location | Findings and resolution |
|---|---|
| Ubuntu packages | `libssl3t64` and `openssl` 3.0.13-0ubuntu3.6 need at least 3.0.13-0ubuntu3.11; `sudo` 1.9.15p5-3ubuntu5.24.04.1 needs .2; `gpgv` 2.4.4-2ubuntu17.3 needs .4. The final stage upgrades installed apt packages. Only the final image scan establishes closure. |
| `/opt/shiny-server/node_modules` | Runtime findings included handlebars 4.7.8, underscore 1.13.7, flatted 3.3.1, path-to-regexp 0.1.10, ip-address 9.0.5, and websocket-driver 0.7.4. The selected upstream graph contains 4.7.9, 1.13.8, 3.4.4, 8.4.2, 10.7.0, and 0.7.5 respectively. Express is upgraded to 5.2.1 with upstream's matching router changes. |
| `/opt/shiny-server/ext/node/lib/node_modules/npm` | Build/package-manager findings included tar, brace-expansion, minimatch, glob, cross-spawn, pacote, sigstore, and ip-address. The staged runtime contains no npm CLI or its dependency tree. The upstream launcher directly executes `ext/node/bin/shiny-server` (Node) with `lib/main.js`; npm is required for the build, not server startup. |

The [public installation guide](https://docs.posit.co/shiny-server/) still points
to 1.5.23.1030. The latest listed
[daily binary, 1.5.24.1035](https://dailies.rstudio.com/shiny-server/1.5/open-source/ubuntu-20.04/1-5-24-1035/),
was built 2026-04-13 and retains the vulnerable runtime dependencies. Its published
SHA-256, `bd834b039ddddee36e7d1d31afecec6e8b08ea8f4112708b31b9cc59b4ae4b8b`,
was verified before inspecting the extracted package. Its broader installed-tree
scan also fails the High/Critical gate (42 High, 4 Critical matches). These numbers
describe that extracted package, not the original registry-image scan.

The implementation instead builds the exact
[official upstream source revision](https://github.com/rstudio/shiny-server/tree/ae268b685b234b3a97a33bf6caa28e6d395cb232)
`ae268b685b234b3a97a33bf6caa28e6d395cb232`, with its `npm-shrinkwrap.json`
preserved byte-for-byte until 2026-10-22 and since then carrying the four
patch-level pins described above. This is a project-built upstream snapshot, and
is not presented as a released Posit binary. Upstream made the Node/runtime-dependency/Express changes
together with compatibility tests. Node 20 is
[end of life; Node 24 is LTS](https://nodejs.org/en/about/previous-releases).

Locked inputs in [`harden-base.sh`](harden-base.sh):

| Input | URL and SHA-256 |
|---|---|
| Upstream source | `https://api.github.com/repos/rstudio/shiny-server/tarball/ae268b685b234b3a97a33bf6caa28e6d395cb232`; `4a1d74aa4e2c266f4aaa04198132a017cb8b1e07ef766d163f7dec5d9dc17a62` |
| Node v24.20.0, Linux x64 | `https://nodejs.org/dist/v24.20.0/node-v24.20.0-linux-x64.tar.xz`; `2f2c0da162318f0de47665410c7c8c2ed3d36c8f3105de4bbc61176c70a7cbf2`, verified against the [official checksum list](https://nodejs.org/dist/v24.20.0/SHASUMS256.txt) |

The build uses the unmodified upstream lock, compiles TypeScript and the native
POSIX addon, and requires both upstream test tiers with no skipped tests. The
builder runs as `shiny`, matching the upstream tests' permission-denial expectations
and the production service's unprivileged operation. It then
uses `npm ci --omit=dev --ignore-scripts` to install the same locked runtime graph,
preserving the tested native addon. Only runtime files, installed package manifests,
the upstream lock, licenses, source identity, and build evidence enter the final
stage. The superseded Debian Shiny Server package and complete old server tree are
removed. The launcher is compiled from the upstream single-file C++ target; its
unrelated legacy Pandoc 1.19 downloader is bypassed, and the existing maintained
`/usr/bin/pandoc` serves the runtime's Pandoc path. R is asserted to remain 4.5.1;
the destiny and DuckDB production pins remain enforced by their existing installer.

Local evidence before the final project commit:

- The exact `harden-base.sh` completed: **315 upstream tests and 8 real-R tests
  passed, no failures or skips**, and the dependency lock stayed unchanged.
- [`tests/smoke-shiny-server.py`](tests/smoke-shiny-server.py) passed against the
  staged runtime with npm/development artifacts absent. It uses the actual
  `shiny-server.conf`, substituting only temporary paths, user and listener port.
  A real R app, its Shiny JavaScript, the server stylesheet, and all three image
  aliases passed; each `Range: bytes=0-15` request returned HTTP 206 and exactly
  the requested 16 bytes.
- Syft 1.51.1 with **`--override-default-catalogers image`** scanned the staged
  directory's installed packages (149 package manifests; 294 CycloneDX components
  including file records). Grype 0.118.0 reported **zero matches**, and
  `--fail-on high` exited 0. Directory-default Syft catalogers do not cover the
  installed JavaScript tree adequately; their earlier exploratory report is not
  the remediation evidence.

Reports are retained locally in `/tmp/islet-ufrc-evidence/base-remediation/`:
`staged-installed-sbom.json`, `staged-installed-grype.json`,
`daily-installed-sbom.json`, `daily-installed-grype.json`, and
`staged-build-evidence/` (test JSON, smoke result, logs, and input provenance).
The full candidate build must repeat these tests and scan the **entire application
image**, including upgraded Ubuntu packages, R/Python packages, and vendored browser
code. This staged server result does not waive that gate or establish release,
UFRC, or IRM approval.
