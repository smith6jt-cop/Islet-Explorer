# Shiny Server base remediation

**Latest owner direction, 2026-11-03: LOCAL REVIEW ONLY.** The display candidate
is running for [local review](http://127.0.0.1:3845/islet-explorer/) from source
`1596b03`, image `846962f7…`, using the separately accepted read-only data built from
`5e6580f`; see the [local review record](../docs/local_duckdb_review_20261103.md).
The persistent local user service is enabled on loopback port 3845. Data and View
controls start collapsed; click a heading to show or hide them. The
[controls review receipt](../docs/evidence/local_controls_review_20261103.json) records this successor’s checks.
Dev activation is paused pending the owner's review and decision. Both public
services remain unchanged.

**Production decision, 2026-11-02:** the owner directed the then-current build to the live production PubApp after its actual data/release/runtime checks passed. The known DuckDB finding is nonblocking through the notice's **2026-11-10 UTC** patch deadline; the owner directed development/testing of our DuckDB fix separately for dev. See [the dated release decision](../docs/duckdb_remediation_window_20261102.md). This superseded the earlier strict DuckDB hold and no-activation directions below for that release. Historical failed reports remain unchanged; data checks and other security findings still apply.

**Current outcome, 2026-11-03:** production and the separate strict-DuckDB dev
release are active and verified. The [production receipt](../docs/evidence/pubapps_live_20261102.json)
and [dev receipt](../docs/evidence/dev_duckdb_deployment_20261103.json) bind the
actual images and their gates. Production uses the dated DuckDB window; dev
passes strict validation without it. Current build source now includes the
reviewed source patch and strict receipt check, cherry-picked as `b5fae5c`;
this does not change the bytes of the already-running production image.
Neither live image includes the subsequent
[local display changes](../docs/app_display_fixes_20261103.md).

Historical review date: 2026-09-08, with the 2026-10-22 dependency pins below.
These source-build records are not a substitute for each new image’s build,
SBOM, vulnerability scan and browser gates.

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

## 2026-10-25: body-parser pin

The first build of the Plan A Level 1 image (commit `4da5ac4`) failed the High/Critical
gate on one advisory published after the 2026-10-22 scan: body-parser 2.3.0, the request
body parser of express 5.2.1 (GHSA-27r3-jxmh-p828, CVE-2026-111181, High: denial of
service through trailing bytes after a compressed request body), fixed in 2.4.0, released
2026-10-22. body-parser is not a direct dependency, so it joins proxy-addr and
follow-redirects as an `overrides` entry (`body-parser 2.4.0`), and its lock entry carries
2.4.0's dependency ranges. Of those, `qs ^6.16.1` was not met by the locked 6.16.0, so qs
moves to 6.16.1 (same dependencies and license, within the root's own `^6.16.0`). Both
lock entries were written from the registry metadata because the workstation's npm
refuses to regenerate the lock (`allow-git`/`allow-remote` none); `npm ci --dry-run`
accepted the pair and both tarballs' integrity hashes were verified before the rebuild.
Nothing else in the graph changes.

## 2026-11-01: Node 24.22.2 and follow-redirects 1.16.2

The first build of the v20 release image (commit `c8d942b`, Grype database built
2026-11-01) failed the High/Critical gate on advisories published after the 2026-10-26
scans, all in the Shiny Server runtime:

- Node 24.20.0 (`/opt/shiny-server/ext/node/bin/node`): CVE-2026-111497 (Critical,
  out-of-bounds write in the bundled Brotli decoder), CVE-2026-105482 (High, NUL-byte
  truncation in `dns.resolve*()`), CVE-2026-105483 and CVE-2026-106069 (High, TLS
  hostname verification of non-canonical IPv4 hosts), from the
  [October 2026 Node.js security releases](https://nodejs.org/en/blog/vulnerability/october-2026-security-releases),
  fixed in 24.22.1 (2026-10-20);
- follow-redirects 1.16.1 (GHSA-mrwp-f3q8-3g9p, CVE-2026-116689, High: custom
  credential headers can follow cross-origin redirects), fixed in 1.16.2 (2026-10-25).

The build now uses Node v24.22.2 (2026-10-27), the newest patch release of the 24 LTS
line, which contains the 24.22.1 security fixes and bundles npm 11.19.1. Its
`node-v24.22.2-linux-x64.tar.xz` SHA-256 was taken from the official
[checksum list](https://nodejs.org/dist/v24.22.2/SHASUMS256.txt) and matched a download
of the archive. The source revision's `.nvmrc` still names v24.20.0; `harden-base.sh` now
checks `.nvmrc` against that upstream value and requires the build's Node to be on the
same major line, instead of requiring the two to be equal. The upstream and real-R test
tiers run on the new Node as before. The follow-redirects override moves to 1.16.2,
whose registry metadata matches 1.16.1's (MIT, no dependencies, `engines` and the
optional `debug` peer unchanged), so its lock entry changes only in version, URL and
integrity; the tarball's SHA-512 and SHA-1 matched the registry, and `npm ci --dry-run`
accepted the lock. Nothing else in the graph changes.

## 2026-11-02: bundled Brotli in Pillow and DuckDB; services stopped

The UFRC alert still identified CVE-2026-111497 in the v20 image even after its
Node and Ubuntu packages had been patched. Inspection found Pillow 13.0.0's
wheel bundles a separate Brotli 1.2.0 decoder, reached through its FreeType
library. The package SBOM did not catalog that copy: rescanning the image with
the November 2 Grype database still reported zero High/Critical findings.

The Python lock now forces Pillow to build from source against Ubuntu's patched
FreeType/Brotli libraries; Pillow's Python API version remains 13.0.0. The
Containerfile installs the required development headers and bounds compilation
to two workers. Updating the system Brotli package cannot repair a wheel's
private decoder, so the image gate rejects that layout and checks Pillow's
actual shared-library linkage to the Ubuntu-managed decoder.
See [Ubuntu's package fix](https://ubuntu.com/security/CVE-2026-111497),
[Pillow's dependency update](https://github.com/python-pillow/Pillow/pull/10340),
and [the remediation record](../docs/container_security_20261102.md).

DuckDB 1.5.2 also compiles a private Brotli 1.1.0 decoder into `duckdb.so`.
The original pinned source lacks the fix, and the inspected installed library
returned the vulnerable error (-9) on the upstream regression vector, instead
of the patched -17. A Pillow-only replacement cannot pass the complete gate.

**Historical direction before the dated production window and separate dev work:** implement the
released fixes and defer DuckDB until its upstream PR merges. The
[DuckDB port proposal](https://github.com/duckdb/duckdb/pull/27791) was
**open and unmerged when reviewed on 2026-11-02**. Production build instructions
therefore retained the original CRAN R package **1.5.2**, its existing R API and
SHA-256 pin `55f98a336513a94354f328791c8b94521b9f2e918a57d7e9fe6c0cb6b9a38db6`,
without applying the proposed decoder patch. The known decoder vulnerability
failed the strict gate; released Pillow, Node, and system-library fixes did
not resolve it. The later production window and separately validated dev fix
are recorded at the top of this document.

The stopped combined-backport attempt is retained on local branch
`security/deferred-duckdb-brotli-20261102` at `3a12daa`, with its source helper, patch,
receipt checks, and the [dated remediation record](../docs/container_security_20261102.md).
Its historical port was `2b942af6319783a5cc3c21d76a8e560e191978b2`, carrying
upstream Brotli fix `0ed4268a4401abcb4472096340477cf4af7d1a2b`; patch SHA-256
was `d00033af815fcb214daa9630d9da11b39a191dfab17b0f8cff367157af5806ca`.
Those are archived review details, not active installation instructions or an
accepted image. The attempt identified an installation pitfall: CRAN configuration
extracts `duckdb/src/duckdb.tar.xz`, so editing only an extracted directory can
be undone. Its helper patched the nested archive and package `MD5`, preserving
unrelated members. The helper and patch were removed at that historical stage;
the subsequently reviewed dev implementation is restored in current build
source at `b5fae5c`, with its receipt and actual binary regression required.

The scan feeds the current checkout's `verify-brotli.py` into the exact immutable
image with networking disabled and no data mounts. The required
`islet.brotli_verification/2` report combines Ubuntu/Node patch checks, Pillow
linkage, DuckDB runtime package version and installed-library hash, and a
bounded subprocess calling the installed
DuckDB decoder. The malformed upstream stream must return decoder result 0
and error -17; timeout, crash, missing symbols, or a vulnerable decoder fail.
The earlier released-fixes-only gate did not require a patch source receipt.
Current build source does require the reviewed patch receipt in addition to
the actual decoder check; the live dev receipt records that strict result.

`scan-image.sh` then runs the SBOM/Grype High/Critical gate and binds the Brotli
report hash in `.scan.txt`. Every **new** release manifest requires the complete
v2 report and its DuckDB binary-check fields, plus a passing strict result or
the exact dated production exception validated by the release security gate. An old scan without the report, or
an earlier Pillow-only report, cannot authorize a new manifest. No CVE
suppression was added. Existing historical manifests remain historical records,
not evidence that an old image satisfies the new gate.

The valid-data check, run successfully for the later dev fix, is `tests/check-duckdb-brotli.R`: it reads a real
release Parquet, round-trips all its rows through BROTLI compression in temporary
storage, checks the actual compression metadata, schema/values/missingness,
duplicate multiplicities, numeric and grouped aggregates, and unchanged source
SHA-256. The [deployment guide](README.md#brotli-source-and-runtime-verification)
gives the read-only container invocation for a future eligible candidate. Local
script-development checks do not demonstrate a fixed image, and a successful
valid-data check cannot waive malformed-stream rejection.

**Historical early November 2 build status:** the combined `3a12daa` backport build was stopped following
the owner's direction. No completed patched image or patched-image compatibility
result is claimed. The released-fixes-only image `2bbb928d…` from `23bd2fd` was
built, scanned and imported offline under a hold. Its package scan has 0
High/Critical findings and dependency/rendering smoke checks pass, but the
original DuckDB decoder alone fails the strict gate. It has no passing release
attestation and cannot be activated. The dated security record binds its exact
image/export/staging evidence. Revisit
the DuckDB change after the upstream condition is met and its implementation
and compatibility have been reviewed.

At that point the owner had stopped dev and both PubApps services were
persistently masked. The later verified activations are recorded above.
Old releases are not eligible rollback targets until rebuilt and reverified.

## Original base review

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
| Node v24.22.2, Linux x64 (since 2026-11-01; upstream's `.nvmrc` names v24.20.0) | `https://nodejs.org/dist/v24.22.2/node-v24.22.2-linux-x64.tar.xz`; `48cc0e21e048790e54c0cf5749c2725005313d652651c40b7a21feb9697aaf75`, verified against the [official checksum list](https://nodejs.org/dist/v24.22.2/SHASUMS256.txt) |
| Node v24.20.0, Linux x64 (until 2026-11-01) | `https://nodejs.org/dist/v24.20.0/node-v24.20.0-linux-x64.tar.xz`; `2f2c0da162318f0de47665410c7c8c2ed3d36c8f3105de4bbc61176c70a7cbf2`, verified against the [official checksum list](https://nodejs.org/dist/v24.20.0/SHASUMS256.txt) |

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
