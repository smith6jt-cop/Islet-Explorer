# Dev DuckDB fix — validated 2026-11-02

**Update, 2026-11-03 12:42 UTC:** the current-app dev build is now live and
verified without a DuckDB deferral. See the
[new deployment record](dev_duckdb_deployment_20261103.md), source `378843e`, image
`e704ac2a…`. The older image and validation below are preserved as history;
production stayed unchanged during the dev deployment.

The owner directed the current component build to the live production PubApp
within the DuckDB patch window, with our DuckDB fix developed separately for
dev. The dev backport is now **built and validated**, without using that deferral.
Production's build remains separate; see the
[dated production decision](duckdb_remediation_window_20261102.md).

| Item | Verified result |
|---|---|
| Dev image | `sha256:b26558a9b3791c05f6bf158152efc7eba29fe1c17159d6a23678bd306a4532d6` |
| Source | `92c63ba914f71e7b96d45efae1f60b3eae81e024`, branch `security/dev-duckdb-fix-20261102`, application/data-builder base `b15f0fc` |
| Validation completed | 2026-11-02 22:37:40 UTC |
| Actual decoder | **PASS**: result `0`, error `-17`; the original vulnerable decoder returned `-9` |
| Installed DuckDB library | SHA-256 `a7f0d156ebd3d510f387289a698ce4cf165ed02eb5b44ba5cb82ebd2989ecbee` |
| Source provenance | **PASS**: exact verified CRAN source, reviewed two-file patch, unchanged unrelated members, installed source receipt |
| Full security gate | **PASS**, 0 High / 0 Critical package findings; no DuckDB deferral |
| Dependencies/application | R/Python probes, all 35 app R source parses, Node package imports, Pillow PNG/JPEG/TIFF/FreeType and Matplotlib Agg checks pass |
| Real-data compatibility | Read-only BROTLI round trips pass for v19 reference `190addc` (4,853 rows × 371 columns) and v20 `c974946` (4,853 × 296), including schema, values, missingness, duplicate multiplicities and aggregates; input files unchanged |
| OCI export | 2,188,034,048 bytes; SHA-256 `fb19874949a6e203a5207d18f590aaa8feac12f0757f09b38fc0cafdd8bae008` |
| Host deployment | Not staged or activated by this work; image and export are retained locally |

The source worktree is
`/data/panc_CODEX/Islet-Explorer-duckdb-dev-20261102`. The implementation preserves
DuckDB's R API/version **1.5.2** and applies the previously reviewed Brotli port
`2b942af6319783a5cc3c21d76a8e560e191978b2` (patch SHA-256
`d00033af815fcb214daa9630d9da11b39a191dfab17b0f8cff367157af5806ca`). It patches the
nested source archive that CRAN's configure actually extracts, updates the package
MD5 entry, and always recompiles without prebuilt object reuse. A fresh repack
reproduced the earlier reviewed archive byte for byte, SHA-256
`1f66d65720c93f856fdcf367619e7207c54e7de9c239983c1b04b91b752381b1`.

The dev build tag is
`localhost/islet-explorer:dev-duckdb-brotli-20261102`. Its evidence and archive are
under `tmp/security-duckdb-dev-20261102/`, with export metadata and checksum
sidecars. The [compact receipt](evidence/dev_duckdb_fix_20261102.json) binds the
real tests, scan, source receipt, image identity and archive.

This validates the dev security fix and compatibility with representative real
release tables. It does not establish acceptance of the new mixed scientific
data root. That build's actual data/release/runtime checks remain required for
live activation. Keep this patch on the separate dev branch; do not substitute
it for the production image as part of the current live-build instruction.
Coordinate host image imports with the production build agent because `/local`
space is limited and both images need distinct, verified release identities.
