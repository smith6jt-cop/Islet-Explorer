# PubApps container security remediation — 2026-11-02

**Latest owner direction, 2026-11-03: LOCAL REVIEW ONLY.** The owner requested
“Instead of sending directly to pubapp, create a local version for review.”
This supersedes the planned next dev deployment. The new patched image and UI
candidate runs only at [the local preview](http://127.0.0.1:3845/islet-explorer/);
see its [review record](local_duckdb_review_20261103.md) and
[receipt](evidence/local_duckdb_review_20261103.json). Its image source is
`7ac421e`; the separately identified accepted dataset was produced by `5e6580f`
and is mounted read-only. No new pubapp image import, data relocation or service
activation occurred. Existing public dev `e704ac2a…` and production `89b3f840…`
remain as recorded below. The local candidate does not authorize their replacement.

**Earlier production direction, 2026-11-02:** the current build is intended for the live production PubApp once its actual data/release/runtime checks pass. The known DuckDB finding is nonblocking through the notice's **2026-11-10 UTC** patch deadline; develop/test our DuckDB fix separately for the dev PubApp. See [the dated release decision](duckdb_remediation_window_20261102.md). This supersedes the earlier strict DuckDB hold and no-activation directions below. Historical failed reports remain unchanged; data checks and other security findings still apply.

## Dev DuckDB fix live — 2026-11-03 12:42 UTC

Following the owner's instruction to continue the dev fix, the
[dev PubApp](https://devislet-explorer.rc.ufl.edu/islet-explorer/) now runs
`pubapps-2026-11-dev-duckdb-v1`, clean source `378843e`, immutable image
`sha256:e704ac2a7f01c0fdec5905cc5b1dff43041d45fb0ba7ebfe3c821a13da4a9e9e`.
The reviewed source backport is compiled into DuckDB 1.5.2. The actual running
decoder returns `0 / -17`, its source provenance passes, and its library hash
matches the scanned image. The strict security gate passes with **0 High /
0 Critical and no DuckDB deferral**. Dev is healthy with a 4 GiB limit.

Fresh data acceptance, full host release verification, the complete running cell
census, 15 public-access checks and 49 browser checks passed. Production's
container, start time, image, service, state, manifests, mounts and aliases stayed
unchanged; its 15 public-access checks passed again. See the
[deployment record](dev_duckdb_deployment_20261103.md) and
[hash-bound receipt](evidence/dev_duckdb_deployment_20261103.json).
Production's separate DuckDB finding remains due **2026-11-10 UTC**. Earlier
stopped/masked-dev observations below describe historical states.

## Production live — 2026-11-03 UTC

[The production PubApp](https://islet-explorer.rc.ufl.edu/islet-explorer/) is
running release `pubapps-2026-11-components-v1`. Production was activated under
the owner's explicit direction after full release/data checks passed; dev was
stopped and persistently masked at that time. The [deployment receipt](evidence/pubapps_live_20261102.json)
records the actual image, data, service state and public acceptance.

The deployed immutable image is
`sha256:89b3f84070b1edb89124347572ed43f278b5a8478e243053957037ccb1ceff8e`,
from clean source `e5689f6453e68bc9e94cbe83689cb9db4d3d0b83`.
Its package scan reports **0 High / 0 Critical** and all eight prior locally
identified High/Critical package findings are closed. The notice's three High
findings were not individually enumerated in its body. Released fixes include
Node 24.22.2, body-parser 2.4.0, follow-redirects 1.16.2, Ubuntu Brotli
1.1.0-2ubuntu0.1, and source-built Pillow 13.0.0. All 35 R files parse;
installed Shiny/R/static-range and actual plot-scale tests pass. The OCI archive
SHA-256 is `d68761c53f9b0d90678640e63b31cc3e6dac9d080b6c937f22a00f3859f85992`.

The data builder uses the same clean source. Scientific HPG job `49168359`
uses `67ee6c9`; the later app-only change preserves dynamic phenotype names in
plot legends. Full-byte data verification, release parity, inventory and the
host's full release verification passed. All 20 donors and 3,199,606 eligible
cells are deployed. A query inside the running production container matches
**every** accepted phenotype count: 7,398 cells are Unclassified (0.23%), no
previously named nonendocrine identities are lost, and all 110,721 selected
endocrine identities are preserved. See the
[independent row review](unclassified_component_review_20261103.md).

The original parity report is preserved with its one failed state-export check.
It treated floating-point distance values as exact keys: 213 CSV/Parquet values
differed by at most `7.1054e-15`, while all Ki67 flags matched exactly. Corrected
verifier `d3ec00e` applies the existing `1e-12` distance tolerance and keeps binary
states exact. Its complete real-data check, independent UUID audit and deliberate
corruption controls passed. The hash-bound acceptance record retains all 24 other
full parity checks and separately identifies the corrected verifier; runtime
source, image and scientific data remain unchanged.

Production is healthy with the expected immutable image and read-only scientific
mounts. The three public image aliases resolve to this release. All public-access
checks passed (zero failures), including anonymous denial, authenticated app and
static access, exact TIFF byte ranges, sidecars and WebSocket upgrade. All
**49 browser checks** passed through the public HTTPS URL. These are
workstation-to-public-URL checks; independent off-campus access and a host reboot
were not tested, and this record does not claim UFIT approval.

**DuckDB remains open:** the release decision is
`pass-with-duckdb-remediation-window`, with the raw failing decoder result
retained and a **2026-11-10 UTC** patch deadline. This is the owner's dated
deployment direction, not a patched-decoder claim or institutional exception.
The separate dev fix subsequently went live at 12:41:40 UTC, as recorded above.
No vulnerable old image was started during activation or recovery.

## Earlier completed builds

The earlier production candidate was image `41e62a97…`, source `b15f0fc`.
It passes package/dependency checks, the installed Shiny/R/static-image smoke,
and the policy at `7941df1` with the dated DuckDB exception. Its failed raw
DuckDB decoder report remains visible. It is built locally and has not been
exported or host-staged. The associated app and QuPath builds are complete;
actual endocrine consistency failures still prevent a data release. See the
[component build record](component_build_20261102.md) and its linked evidence.

The separate dev image `b26558a9…`, source `92c63ba`, passes the strict decoder
and full security gate **without a deferral**, plus real v19/v20 release-table
BROTLI compatibility checks. It is exported locally and has not been host-staged
or activated. The [dev validation record](dev_duckdb_fix_20261102.md) binds the
image, source patch, scans, tests and archive. A read-only host check at
22:40:10 UTC confirms both services remain stopped/masked with no backend ports.

The earlier build sequence below is retained as history; it does not override
the dated production direction or these completed results.

## Earlier request and build sequence

The owner requested a build with the fixes required by today's IT security email,
then clarified: **released fixes only; leave DuckDB aside until its upstream PR
merges**. This supersedes the attempted build containing the unmerged DuckDB
backport. The email for **Request85407**, received at
**2026-11-02 19:37:28 UTC**, names **2 containers, 2 with KEV**. The retrieved
message is a CC notification; it does not supply the detailed package findings
below. Those findings come from the retained local image inspection, rescans and
decoder probes, not from the notification body.

The owner pasted the underlying notice into the parallel security session. That
notice lists dev image `517e969d…` with 0 fixable Critical/High findings but KEV
`CVE-2026-111497`, due **2026-11-10**, and historical production image
`ae95fd23…` with 0 fixable Critical / 3 High, the same KEV due **2026-11-10**,
and High findings due **2026-12-02**. It says production was last seen on
2026-10-26 and was not running, and must be patched before restarting. Those
are the notice's counts and dates; the later local rescan counts below are
separate observations. The pasted notice does not enumerate the three High IDs.

The combined current-app/backport build started at **20:30:16 UTC** from clean
commit `3a12daa5670806e264f703d97911d1f0340df530`. After the clarification, it
was stopped at **20:39:51 UTC**, exit `130`. No completed candidate image, image
export or service restart resulted. Subsequent process/container checks found no
remaining build/compiler process or container, and the candidate tag did not
exist. Evidence: `tmp/security-build-current-20261102/stopped.json` and
`build.exit`. Source, logs and worktree remain preserved under the historical
branches described below.

The parallel security session **`01a19dc6-aab3-7382-8e98-7bb7e0012a20`** now owns
one **current-app, released-fixes-only offline build** from
`23bd2fdd223d090f4c5ff438167245c459dbeb16`, started at **20:42:45 UTC** in
`/data/panc_CODEX/Islet-Explorer-security-20261102` on
`security/available-fixes-20261102`. The tag is
`localhost/islet-explorer:security-available-duckdb-deferred-20261102`; evidence
is under `tmp/security-available-fixes-20261102/`. Its package scan and any
export must be marked **`remediation-incomplete`**: the unchanged DuckDB decoder
is expected to fail the normal security gate. No passing security attestation,
deployable release manifest or activation may be inferred from that offline
artifact. The local build, package scan, smoke checks and export have now
completed as recorded below. Building it does not assemble or accept a new
scientific data root.

## Built offline candidate — released fixes verified, DuckDB unresolved

| Item | Recorded result |
|---|---|
| Immutable image | `sha256:2bbb928db9760c86d7a1bb9adf9002bb1a27d88da96b02a600baca93c090d733` |
| Source | Clean `23bd2fdd223d090f4c5ff438167245c459dbeb16`; original CRAN DuckDB installer matches `5b3c535`, backport helper/patch absent |
| Package scan, 21:04:59 UTC | 0 Critical / 0 High; 4,724 Medium / 189 Low / 196 Negligible; Grype database built 2026-11-02 06:26:52 UTC |
| Earlier production rescan findings | All eight Critical/High package matches are no longer detected; exact package/version comparison in `fixed-findings.json` |
| Released components | Node 24.22.2, body-parser 2.4.0, follow-redirects 1.16.2; Ubuntu `libbrotli1` and `libbrotli-dev` 1.1.0-2ubuntu0.1; source-built Pillow 13.0.0 linked to the managed decoder, no vendored Brotli in the venv |
| Runtime checks | Python/R dependency probes and all 35 app source parses pass; exact Node package imports and Pillow PNG/JPEG/TIFF, FreeType text and Matplotlib Agg smoke checks pass against the immutable image |
| Full security gate | **FAIL**, solely DuckDB 1.5.2's bundled decoder returning `-9` instead of `-17` for CVE-2026-111497; the backport was not applied |
| OCI export | 2,188,395,008 bytes, exported 21:05:03 UTC; SHA-256 `6509df201e1fb9e3d1294d401a168a9ee7c523ad3bc564f876169f92da2a35c1` |
| Release authority | `remediation-incomplete`; no passing `.scan.txt`, `published-image.txt`, `release.env` or `MANIFEST.sha256`; no new data release or activation |

The archive is
`tmp/security-available-fixes-20261102/islet-explorer-available-duckdb-deferred-20261102.oci.tar`,
with `.metadata` and `.sha256` sidecars. The immutable image is bound to the
same source and base digest in its labels and export metadata. Its package scan
is deliberately separate from the failing full gate. `HOLD.txt` prohibits
activation; the active service masks remain required. No complete DuckDB
remediation or scientific acceptance is claimed by these results.

Offline import to `pubmaigan1.ufhpc` completed at **21:08:24 UTC**, under
`/pubapps/maigan/islet-explorer/incoming/security-available-20261102` with its
`HOLD.txt`. All 32 transferred artifact checksums passed; `load-image.sh`
verified the same immutable image, source and base identity. Before/after
records confirm unchanged service masks, both Quadlet files and all four
deployment-state files, no running containers and no listeners on 3838/3839.
Both services remain stopped/masked; no activation occurred. Host local storage
has 5,962,076,160 bytes free after import. `host-load-receipt.json` and
`host-after-staging.json` retain these checks. See the compact
[build and staging receipt](evidence/security_available_fixes_20261102.json).
After adding the host receipts/logs, all **38** final `ARTIFACTS.sha256` entries
passed on the server at **21:10:11 UTC**. The final check again confirmed the
same imported image, unchanged masks/configuration and no running containers
or backend listeners. Proof: `tmp/security-available-final-verification-20261102.log`.

## Host mitigation

The owner requested remediation of the UFRC container alert and then explicitly
directed: **stop dev**. `islet-explorer-dev.service` on `pubmaigan1.ufhpc` was
stopped and persistently masked on 2026-11-02. Production was already stopped;
`islet-explorer.service` was also masked to prevent its old image starting on
reboot. Neither service is authorized to restart as part of this remediation.

At 19:45:58 UTC both units were masked, `podman ps` showed no running containers,
and neither port 3838 nor 3839 had a listener. Production retains an earlier
`failed` status; that is historical state, not a running process. Evidence:
`tmp/security-remediation-20261102/host-mitigation.txt`. The same state was
confirmed at 19:56:43 UTC in `host-mitigation-final.txt` in that directory.
Read-only inspection at 20:38:13 UTC again found both masks, no containers and
no backend listeners: `tmp/security-build-current-20261102/host-status-readonly.txt`.

The masks are in `/pubapps/maigan/.config/systemd/user/`, pointing to `/dev/null`.
The Quadlet files, application data, and release manifests were preserved.
Activation and rollback cannot restart these masked services. A future restart
requires explicitly lifting the applicable mask, a patched image with current
security evidence, and the usual release/data acceptance. Do not unmask a service
while it still references a vulnerable image.

## Findings from the existing images

| Image | Evidence and disposition |
|---|---|
| Dev `sha256:517e969d2333744f609533b8f4a7fd217ce8a852fe30ef3551e0d805c2a0db00` | Node 24.22.2 and Ubuntu `libbrotli1` / `libbrotli-dev` 1.1.0-2ubuntu0.1 already contain the fixes. Pillow 13.0.0's wheel separately bundles Brotli 1.2.0, and DuckDB 1.5.2 statically bundles Brotli 1.1.0. The new library check rejects both copies. Dev is stopped and masked. |
| Older production image from the retained investigation, `sha256:ae95fd238c72e0ec239942be0050ad72eb441bf34149218973edd41f5c7e9af9` | Current rescan: 1 Critical / 7 High matches, including Node 24.20.0, system Brotli, body-parser 2.3.0 and follow-redirects 1.16.1. It also contains the affected Pillow wheel. This historical release must not be used for rollback or restarted. |
| Last configured production `sha256:1065d532fe7bc57ca72eb7cbecd2557b5dc645195e22c4b09d0db8c0b41d6d5a` | The Quadlet actually selects `pubapps-2026-10-uniform-v10-perf2` (`0b55429`), not the older image named in the retained alert investigation. Production is stopped and masked. Only the earlier, aborted production build preserved perf2's application code; the newly assigned offline build uses current app source with released fixes only. |

Grype's refreshed database is built `2026-11-02T06:26:52Z`. Its rescan of the dev
SBOM still returns 0 High / 0 Critical because this SBOM does not identify
Pillow's bundled Brotli as a component. A clean package scan alone therefore
does not close CVE-2026-111497. No CVE suppression was added. The retained counts
and package versions are in `tmp/security-remediation-20261102/before-summary.json`;
these are observations of the existing images, not results for the new build.

Upstream sources:

- [Ubuntu CVE-2026-111497](https://ubuntu.com/security/CVE-2026-111497): Noble's fixed package is 1.1.0-2ubuntu0.1.
- [Node's October security release](https://nodejs.org/en/blog/vulnerability/october-2026-security-releases): Node 24.22.1 contains the fix. Node 24.22.2 still reports bundled Brotli 1.2.0 despite carrying the backport.
- [Pillow's Brotli update](https://github.com/python-pillow/Pillow/pull/10340) and [source-build instructions](https://pillow.readthedocs.io/en/stable/installation/building-from-source.html): build the pinned Pillow 13.0.0 release against patched system libraries instead of its affected prebuilt wheel.
- [Google's decoder fix](https://github.com/google/brotli/commit/0ed4268a4401abcb4472096340477cf4af7d1a2b) distinguishes an invalid block switch from a request for more input and propagates the error through the decoder. The [v1.1.1 security release](https://github.com/google/brotli/releases/tag/v1.1.1) contains backport `da9806362874b9bd4b552a3af360b5d6347beb12`.
- [DuckDB issue 27660](https://github.com/duckdb/duckdb/issues/27660) and [PR 27791](https://github.com/duckdb/duckdb/pull/27791): the reviewed DuckDB port is **open and unmerged** as of 2026-11-02. Immutable port commit `2b942af6319783a5cc3c21d76a8e560e191978b2` was examined and used only in the now-stopped attempt. It is not part of the active released-fixes-only build method.

The pinned CRAN DuckDB 1.5.2 archive was downloaded and verified against the
existing SHA-256 `55f98a336513a94354f328791c8b94521b9f2e918a57d7e9fe6c0cb6b9a38db6`.
Its nested source archive identifies Brotli 1.1.0 and lacks the block-switch fix.
The actual installed `duckdb.so` was tested in an isolated, network-free process
with the upstream regression vector. It returned error **-9**, the vulnerable
behavior; the patched Ubuntu decoder returned **-17** for the same vector.
The bounded regression check is now part of the image gate. This establishes
why updating only the system package cannot repair DuckDB's statically bundled
decoder, and why a source receipt or version string alone is insufficient.

## Earlier build attempts, preserved as history

The first attempts used isolated worktrees to preserve the code of the existing
deployments:

- Dev: branch `security/brotli-20261102`, initial build commit `255ee8b`, based on deployed `c974946`.
- Production: branch `security/production-brotli-20261102`, initial build commit `a120094`, based on last deployed `0b55429`, also carrying Node 24.22.2 and follow-redirects 1.16.2 pins.

Both builds were deliberately terminated after confirming DuckDB's remaining
vulnerability. Neither produced an exported, transferred or activated replacement.
The branches also retain the later DuckDB gate; these commit IDs identify the
initial attempts, not completed security releases. Logs and scanner evidence
remain under `tmp/security-remediation-20261102/`.

The later combined attempt used current app source plus the reviewed decoder
backport at `3a12daa`, rather than either old deployed app version. It ran from
20:30:16 to 20:39:51 UTC and was stopped after the owner's released-fixes-only
clarification. The exact source is preserved at
`security/deferred-duckdb-brotli-20261102` and the existing
`security/current-app-20261102` branch, both at `3a12daa`, with worktree
`/data/panc_CODEX/Islet-Explorer-security-combined-20261102`. Its command, input
hashes, log and stop receipt remain in `tmp/security-build-current-20261102/`.
An intentionally untracked `DEFERRED_DUCKDB_BACKPORT_DO_NOT_BUILD.md` in that
worktree records the owner direction and causes the clean-worktree build guard
to reject accidental reuse. The committed historical source remains unchanged.
These are deferred technical evidence, not the active implementation or a
completed security release.

## Historical released-fixes-only scope and gate

`deploy/requirements-python.lock` now requires a source build of Pillow 13.0.0;
the Containerfile explicitly installs FreeType and zlib development headers and
bounds Pillow compilation to two workers. Released Node 24.22.2 and the fixed
Ubuntu Brotli packages remain the intended runtime dependencies. The active
source removes the unmerged DuckDB helper/patch and its installer, Containerfile,
context and CI wiring. The original pinned CRAN installer is restored; a
backport source receipt is no longer required or claimed.

`deploy/verify-brotli.py` retains **`islet.brotli_verification/2`** and the actual
installed-decoder check. It checks Ubuntu and Node patch levels, rejects vendored
Brotli libraries in the application venv, and confirms Pillow's FreeType path
reaches the Ubuntu-managed decoder. It records DuckDB's installed library SHA-256
and requires decoder result `0` with error `-17` for the upstream malformed-stream
vector. Timeout, crash, unknown symbol or a different decoder result fails the
gate. There is **no active requirement for the unmerged backport's source
provenance**.

The normal `scan-image.sh`/release-manifest path still requires a passing,
hash-bound supplementary report as well as its package scan. It is expected to
reject the released-fixes-only image because CRAN DuckDB 1.5.2 still bundles the
affected decoder. An offline package scan/export may document the fixes that are
available now, but must retain that failure and the **`remediation-incomplete`**
status. A clean package scan must not be converted into a passing normal security
attestation. DuckDB remediation is deferred under the owner's instruction; the
strict binary gate has not been waived.

## Deferred backport implementation at `3a12daa` — historical only

DuckDB remains pinned to the [CRAN 1.5.2 source archive](https://cran.r-project.org/src/contrib/Archive/duckdb/duckdb_1.5.2.tar.gz),
whose SHA-256 is
`55f98a336513a94354f328791c8b94521b9f2e918a57d7e9fe6c0cb6b9a38db6`.
The reviewed [immutable port patch](https://github.com/duckdb/duckdb/commit/2b942af6319783a5cc3c21d76a8e560e191978b2.patch)
has SHA-256
`d00033af815fcb214daa9630d9da11b39a191dfab17b0f8cff367157af5806ca`.

The preserved `deploy/patch-duckdb-brotli.py` verified both hashes before modifying
source. It applied the exact two-file decoder/header patch with no fuzz or offsets
and checked both files' expected before/after hashes. CRAN's `configure` extracts
`src/duckdb.tar.xz` during installation, so patching only an unpacked copy would
be overwritten. The helper instead repacked that nested archive, preserved member
order and metadata, updated its one package `MD5` entry, and verified every
unrelated member was unchanged. It rejected modified inputs, an already-patched
input, or existing output/receipt paths. This helper is retained in Git history
and the preserved worktree; it has been removed from active source.

That historical installer always rebuilt, unset `DUCKDB_R_PREBUILT_ARCHIVE`, used
`R CMD INSTALL --preclean` with two workers, and checked its exit status. Its
intended receipt was `/opt/islet-explorer-build/duckdb-brotli-source.json`. The
historical v2 gate also required exact source-receipt pins and file hashes,
alongside the binary regression. Those unmerged-source requirements and installer
changes were withdrawn with the backport. The source checks below remain valid
evidence of the attempted transformation; they do not authorize its use now.

## Validation recorded so far

| Check | What actually ran | What it establishes |
|---|---|---|
| Original dev decoder and system decoder | Actual shared-library probes; dev returned `-9`, patched Ubuntu returned `-17` | Confirms the original DuckDB defect and the expected rejection behavior; not a test of the new image |
| Pillow mitigation | Disposable container: source-built Pillow, system linkage, PNG/JPEG/TIFF, font and matplotlib smoke checks passed; DuckDB still failed | Validates the separate Pillow repair; not complete image remediation |
| Deferred DuckDB source backport | Actual pinned CRAN archive, exact patch application, two byte-identical repacks, all 229 package MD5 entries and R `tools::checkMD5sums` passed; altered inputs and existing outputs rejected | Historical source transformation/packaging evidence only; the owner has deferred this unmerged fix |
| Image/release evidence contract at `3a12daa` | Mocked Podman/scanner fixtures, including missing/partial/failed Brotli report rejection | Verifies that historical gate's plumbing; fixture image IDs and scan results are not real image attestations or proof of the successor gate |
| Active evidence contract at `23bd2fd` | Mocked contract suite passed after removing the deferred source-receipt requirement; actual old-dev image still fails on its bundled Pillow decoder and DuckDB error `-9` | Confirms the active gate remains strict; no claim that the new candidate passes |
| Combined `3a12daa` build | Stopped at 20:39:51 UTC, exit `130`; no completed image or export | No replacement-image validation or release acceptance |
| Released-fixes-only offline build | Image `2bbb928d…` constructed; package scan 0 High/Critical, dependency/rendering smoke pass; full gate fails only DuckDB; OCI exported | Released fixes verified; artifact remains `remediation-incomplete`, with no full-security pass or activation |

The source backport receipt and checks are in
`tmp/security-build-20261102-review/patch-receipt-a.json` and
`helper-validation.json`. The local repeated repacks produced SHA-256
`1f66d65720c93f856fdcf367619e7207c54e7de9c239983c1b04b91b752381b1`;
these are historical backport bytes, not an active build input. The historical
contract test log is `tmp/security-build-current-20261102/test-image-evidence.log`.
Actual old-image probe/scan and Pillow evidence remain under
`tmp/security-remediation-20261102/`.

At that earlier stage, both services remained stopped and persistently masked.
Successful source checks or mock tests do not permit restarting an old image. An offline artifact marked
`remediation-incomplete` is not deployable. A future deployable replacement still
needs its actual binary, valid-data and security checks, followed by separate
release/data acceptance and activation.

## Ticket reply (draft only; not sent)

> `maigan/islet-explorer-dev` on `pubmaigan1.ufhpc` has been rebuilt with our
> reviewed DuckDB/Brotli source backport and activated as
> `pubapps-2026-11-dev-duckdb-v1`, source `378843e`, immutable image
> `sha256:e704ac2a7f01c0fdec5905cc5b1dff43041d45fb0ba7ebfe3c821a13da4a9e9e`.
> Its actual running decoder passes the CVE-2026-111497 regression check;
> source provenance and strict security checks pass without a DuckDB deferral.
> The package scan reports 0 High / 0 Critical. Full data/host verification,
> runtime cell census, 15 public-access checks and 49 browser checks passed.
>
> Production `maigan/islet-explorer` has been rebuilt and activated as release
> `pubapps-2026-11-components-v1`, source `e5689f6`, immutable image
> `sha256:89b3f84070b1edb89124347572ed43f278b5a8478e243053957037ccb1ceff8e`.
> The available Node/system/Pillow and dependency updates are applied. The
> current package scan reports 0 High / 0 Critical; all eight prior locally
> identified High/Critical package findings are no longer detected. Full data,
> release, running-container and public-access checks passed.
>
> CVE-2026-111497 in production's DuckDB bundled Brotli decoder remains open and
> is tracked for the notice's November 10 patch deadline. The project owner directed this
> production deployment within that remediation window. The dev deployment did
> not change production; its health, identity and public access were verified
> again. We are not claiming production DuckDB remediation or an institutional
> exception.
