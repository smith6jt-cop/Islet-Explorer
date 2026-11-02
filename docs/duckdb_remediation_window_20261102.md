# DuckDB remediation window and production/dev split — 2026-11-02

The owner clarified that DuckDB should not block deployment because the notice
allows several days to patch it. The notice's KEV patch date is **2026-11-10**.
The current component build is intended for the **live production PubApp** once
its actual data, release and runtime checks pass. A separate build will carry
our DuckDB fix for the **dev PubApp**. This supersedes the earlier instruction
to hold all releases until the upstream PR merges; the original production
candidate continues using pinned CRAN DuckDB 1.5.2.

The authorization is recorded in
[`deploy/duckdb-remediation-window.json`](../deploy/duckdb-remediation-window.json).
Its scope is the known `CVE-2026-111497` failure in DuckDB 1.5.2's bundled Brotli
only, through **2026-11-10 UTC**. It does not say that the decoder is fixed or
claim a separate institutional approval. Existing failed scan reports and the
earlier offline image's 38-file evidence bundle remain unchanged.

## Release behavior

`verify-brotli.py` continues to measure and report the actual decoder result.
`evaluate-brotli-report.py` permits only the complete version-2 report with the
single known DuckDB error, decoder result `0`, error `-9` and an installed-library
SHA-256. Missing binary evidence, crashes, timeouts, unknown decoder behavior,
Pillow/Node/system-package failures, and other High/Critical package findings
still block. A fixed decoder returning `-17` needs no deferral.

The scanner retains the raw report as `result: fail` and creates a distinct
`RESULT=pass-with-duckdb-remediation-window` attestation. It binds the raw report
and a copy of the dated authorization by SHA-256, records the patch date, and
separately records the passing package scan. It adds no global CVE suppression.
The manifest generator verifies this evidence and carries the decision/deadline
into `release.env`. Release verification, including `--quick`, rejects an
invalid or expired deadline before activation. It rehashes the security reports,
authorization and image identity on every check, independently of a previous
full-data verification stamp; changing or removing the release decision cannot
hide the raw failure. The authorization is not renewed
automatically. Running services are not automatically stopped by these scripts;
the recorded patch deadline remains an operational obligation.

The data build and its real validation remain the production prerequisites.
This change does not convert an incomplete diagnostic data root into an accepted
release. The production build agent owns the current matching image/data bundle,
validation and live activation. The security agent owns the isolated dev fix.
No additional approval is required merely because the recorded DuckDB finding
remains open during this window.

## Validation

The image/release contract suite passes with the new dated decision. It verifies
normal fixed-decoder releases, the nonblocking known finding, retained failure
evidence, unrelated High/Critical rejection, missing/tampered authorization,
deadline changes, and full/quick verification. Seven decision tests also cover
both date boundaries and rejection after the deadline, independently of the
live clock. The real old dev image remains rejected because its affected Pillow
library is outside the DuckDB-only deferral. Deployment/rollback recovery tests
also pass. An independent read-only review confirmed the full/quick evidence
checks after the initial quick-check gap was corrected.

Fresh scan evidence is retained separately under
`tmp/security-duckdb-window-20261102/`; it does not overwrite the prior strict-gate
receipts. The real released-fixes image `sha256:2bbb928db9760c86d7a1bb9adf9002bb1a27d88da96b02a600baca93c090d733`
was rescanned at **22:04:10 UTC**: package scan 0 High / 0 Critical, raw DuckDB
report still `fail` with error `-9`, release decision
`pass-with-duckdb-remediation-window`. The current component image must be rescanned by immutable image ID
using this tooling when its build finishes. Its original build/source identity
must remain unchanged.

The [compact validation receipt](evidence/duckdb_remediation_window_20261102.json)
binds the tests, real scans, dated authorization and security-tooling hashes.
