# Avivator vendor evidence

Reviewed 2026-09-08. This directory contains a checked-in static Avivator build,
not a reproducible source checkout. The vendor files in `assets/` were retained
unchanged for the PubApps hardening work. `VENDOR_SHA256SUMS` identifies the exact
reviewed bytes; it is not an upstream signature or a vulnerability clearance.

The current files entered this Shiny path in repository commit
`137b0c9c5ff4e9de323d33c68c4fe987edde8764` (2025-09-17). The legacy
`scripts/install_avivator.sh` can build an available local source tree or mirror
`https://viv.gehlenborglab.org/avivator/`; it neither pins an upstream revision nor
records a dependency lockfile. That script also names the older `shiny/www` path.
It must not be used to replace the reviewed production files.

## Evidence and limits

The application identifies `https://github.com/hms-dbmi/viv` as its upstream
repository. No upstream Avivator release tag/commit, complete package manifest,
dependency lockfile, source map, or reproducible build record accompanies this
copy. Minified literals provide partial dependency evidence:

| Component | Observed bundle string |
|---|---|
| React | `18.3.1` |
| React DOM | `18.3.1` renderer metadata and `18.3.1-next-f1338f8080-20240426` export |
| deck.gl | `9.1.11` |
| luma.gl | `9.1.9` |
| loaders.gl | `4.3.2` |
| JSS | `10.10.0` |
| pako | `2.0.3` license banner in the main bundle |

These strings are observations, not a complete or independently verified npm
inventory. Preserve embedded license notices. A container SBOM or file checksum
does not establish that every dependency bundled in minified JavaScript was
identified and scanned. Keep this provenance limitation visible in release
evidence. An upstream rebuild requires a pinned source revision and lockfile,
license/dependency review, and image/channel/rendering parity tests before it
replaces this build; do not silently substitute a current public-site mirror.

## Supplemental local scan

On 2026-09-08, Retire.js 5.7.0 standard and experimental deep scans ran against
the local `js`/`mjs` files with a downloaded, fixed signature database. Its
repository commit was `7ddda6314785f7e4dbca442296aa399520e13c13`, database
SHA-256 `79b3fc3efb2c880bcf9a5745d08db7067fd4b65c74de470e0f433ece823915d7`.
Standard detection found React DOM 18.3.1 and react-is 16.13.1; deep detection
also recognized the React DOM prerelease string already listed above. Neither
scan reported vulnerabilities or scanner errors for the detected signatures.

This result does not cover every component listed above or establish a complete
transitive npm inventory. Local reports, recognized-component CycloneDX,
scanner lockfile, database, and provenance are stored at
`/tmp/islet-ufrc-evidence/browser-vendor/` for transfer to the restricted evidence
record. The source tree was not yet at its final release commit. No application
or scientific data were uploaded. See [the readiness index](../../../../docs/ufrc_readiness.md)
for the remaining exact-image and assessment evidence.

## Local containment and compatibility

`index.html` applies an iframe-specific Content Security Policy. `bootstrap.mjs`
validates the URL before importing the vendor bundle and retains the existing
OME XML element-prefix compatibility fix. The vendor's channel configuration
handling remains unchanged. Only a same-origin `.ome.tif` or `.ome.tiff` file
directly below the current app's `local_images/` or `images/` route is accepted.
Missing, duplicate, external, credential-bearing, and traversal URLs display a
fixed message without importing the vendor's random remote-demo fallback.

The vendor URL-entry form is disabled after React mounts and on subsequent
remounts; capture-phase submission blocking is installed before importing React.
The outer Islet Explorer selector remains the supported image selection path.
The vendor's browser-local file picker is retained; it reads local files in the
browser and does not provide an upload endpoint. `connect-src 'self'` blocks
external fetch/XHR/WebSocket traffic, including the demo/CDN fallbacks still
present in vendor bytes. Same-origin assets, image ranges and sidecars, and
local/blob workers remain permitted. This boundary is not a replacement for
whole-site proxy authentication or dependency security review.

Run `node deploy/tests/test-viewer-boundary.mjs` from the repository root. Browser
acceptance also needs a real TIFF under the deployed app prefix: verify tile
rendering and ranges, channel names, disabled URL entry after menu remount,
absence of external requests, and no vendor import for rejected URLs.
