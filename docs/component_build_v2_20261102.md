# Selected-component v2 build — 2026-11-02

The owner delegated the remaining selection and production activation. The
[resolution contract](component_resolution_20261102.md) defines the implemented
rule: preserve the exact named v19 endocrine component, apply current
non-endocrine methods to the remaining cells, and isolate the DG Lymphatic
revision from other components. The [first component build](component_build_20261102.md)
is retained as historical evidence; its boundary placeholders are not the
final assignment policy.

## Build and provenance

HPG job `49168359` runs through the sealed site environment with scientific
source `67ee6c9d0b28f467eee589fea24829a6546e3d74`. The matched production image
and schema-5 release builder use `e5689f6453e68bc9e94cbe83689cb9db4d3d0b83`.
The latter adds only the phenotype-palette name fix and its test. Scientific
and release manifests must retain these distinct actual producer commits.

| Artifact | Path relative to the repository |
|---|---|
| Config | `tmp/component-build-20261102-v2/build-config.json` |
| Control | `data/development/uniform-20donor-20261102-component-control-v2` |
| Candidate | `data/development/uniform-20donor-20261102-component-candidate-v2` |
| App | `data/development/uniform-20donor-20261102-component-app-v2` |
| Reports | `output/component_candidate_build_20261102_v2` |
| Intended QuPath review | `data/development/qupath-review-20261102-components-v2` |
| Release worktree | `/data/panc_CODEX/Islet-Explorer-component-release-20261102-r2` |
| Release build/evidence | `tmp/pubapps-component-release-20261102-v2` |

Config SHA256: `cc6085cd7a0d18528b67b9890552c7ce28ba082ab75ae91f8aa52b3eb3bd2457`.
Source synchronization transferred an explicit committed file list and Git
objects. All seven protected scientific files retain their recorded hashes,
inodes and timestamps on both sites; the site seal and source package are
unchanged. No source update is permitted during the running job.

## Observed results

At 2026-11-03 00:16 UTC (November 2 local), all 20 control donors completed:

| Measure | Control result |
|---|---:|
| Eligible cells | 3,199,606 |
| Exact selected named v19 endocrine cells | 110,721 |
| Final unnamed endocrine definitions | 0 |
| Technical component-boundary placeholders | 0 |
| Rows reconsidered by the selected-component rule | 11,100 |
| Explicit Unclassified cells | 7,396 (0.2312%) |
| Source unnamed endocrine-object cells, retained as provenance | 7,151 |

Unclassified reasons are shared-marker-only (2,868), no marker above background
(1,412), structure-marker-only (1,479), and endocrine-marker-only outside the
selected membership (1,637). These counts are measured outcomes, not an
acceptance cutoff. All rows and their original evidence remain present.
Existing Tissue labels are separate and unchanged.

The all-20 DG candidate subsequently completed with 7,398 Unclassified cells
(0.23122%), 11,101 resolution targets, the same exact 110,721 selected named
endocrine cells, and zero unnamed endocrine definitions or boundary placeholders.
The two extra Unclassified outcomes occur in Aab_6090 and T1D_6551.

The all-20 Lymphatic invariance audit passes: zero unexpected context, identity
or REDSEA changes; marker positivity, canary and fold states are unchanged.
All **353 final label changes** and every natural label change involve Lymphatic and changed Lymphatic
geometry. The two extra Unclassified outcomes are therefore permitted Lymphatic
removals, not changes to another component. Summary SHA256:
`d04598815bd7bcb2c95ac8090783a29f8372006402277807de9441cd4997623a`.

Independent evidence review covers every control Unclassified row. In the
previous mixed app, these exact UUIDs were 5,170 Islet, 436 EndocrineCluster,
153 EndocrineSingle and 1,637 technical boundary-conflict labels. **None had a
named non-endocrine phenotype in that app.** All positivity/canary/fold fields
and Ki67 replay exactly from actual REDSEA and selected levels. No row meets
current structure qualification or non-endocrine epithelial residual evidence.
The nine historical v19 immune names were seven boundary conflicts and two
Islet labels in the prior mixed app; historical names are not current supported
calls. Full evidence is retained in
`output/unclassified_control_review_20261103/summary.json`.

The final-candidate review also passes across all 3,199,606 UUIDs. All 7,398
Unclassified cells previously had an anatomical label or technical boundary
conflict. Across the 11,101 reconsidered rows, 3,672 gain a supported
non-endocrine label and 31 retain a real conflict/doublet outcome. No other
prior label changes: all 716,451 Tissue calls, exact selected endocrine calls
and every eligible Ki67 flag are preserved. See the
[independent review](unclassified_component_review_20261103.md) and
[hash-bound receipt](unclassified_component_review_20261103.json).

The all-20 classifier-lineage gate passes: 10,080 supported final objects,
5,227 exact Single whole-cell associations and zero failures. It records 178
split descendants, including three with multiple supported parent candidates.
There are 13 Single whole-cell associations that differ from centroid-based
membership; this broader ownership comparison is distinct from the three
wrong-count anomalies in the old diagnostic. The other ten contain eleven
already excluded missing-nucleus detections; each still has exactly one eligible
member matching its source cell. [All thirteen are explained](endocrine_classifier_discrepancy_20261102.md#complete-single-association-reconciliation-in-v2)
from retained context and QC evidence. Both diagnostics remain explicit.
Summary SHA256:
`1282fd48f76664ac3e69433b6cc0d09680e5f52e80f45096da6f0ec40702c209`.

The complete scientific job finished at 2026-11-03 00:42:10 UTC, exit `0:0`,
in 30 minutes 27 seconds. All nine stages and three hash-bound all-20 gates
pass. The final app contains 4,848 paired analysis objects, 5,227 Singles
counted separately, 1,392 composition-embedding objects and the unchanged
31-marker panel. All 3,199,606 eligible cells, 4,848 object CSVs / 357,658
rows, exact four-type endocrine counts and all 78,680 Ki67-positive states
pass export verification. Export gate SHA256:
`d703084e8bfeba1811b495446c83467f9f18c4bca478419ed5e293ed8e6257ea`.

The fresh app return verifies 4,967 regular files / 1,195,349,441 bytes and
20 preserved image links; reports add 98 files / 18,361,975 bytes. Source
outputs remained unchanged during transfer. Return receipt
`tmp/component-build-20261102-v2/return-app.receipt.json` SHA256:
`0e3fb4302a4c9593a479a4cfea9bace6b0f50dc6e650f5cf08f2190dabf70da6`.

Immutable schema-5 packaging began at 00:45:40 UTC from the clean `e5689f6`
release worktree. It uses independent copied files and an image-tree
expectation bound to the original seal and the verified 80-file host inventory.
Full byte verification, source/release parity, inventory and production
activation remain pending.

QuPath job `49170425` completed `0:0` in 8 minutes 30 seconds. All 20 imports
pass, reconciling all 3,525,977 detections as 3,199,606 eligible and 326,371
excluded records. The complete app, reports, typing audit and QuPath return
verifies **5,280 regular files / 16,091,729,844 bytes**, plus 20 preserved
image links; source outputs did not change during transfer. All 20 image
URIs inside the native QuPath project resolve to the expected local images,
with no HPG-only image paths. The combined return summary is
`tmp/component-build-20261102-v2/artifact-return-summary.json`, SHA256
`08f18fdf57ea6260fdce1df037d65eca4261197c068f62f9f42832eb6f9c4d57`.
The [compact scientific receipt](evidence/component_build_v2_20261102.json)
binds completed gates, counts, diagnostics, source protection and returned
artifacts. Its producer flags remain unchanged: these checks establish rule
consistency and provenance, not independent biological ground truth.

## What the additional object proof establishes

The new audit binds exact classifier replay and geometry-provenance hashes.
It verifies original positive-mask support and preserved geometry ancestry;
split descendants retain all supported parent candidates when historical
construction did not record a unique parent. Each Single must match exactly
one original nucleated whole-cell outline. This does not assert that every
final corrected cell mean or split descendant is hormone positive.

The earlier corrected-cell diagnostic remains visible, including its negative
means and centroid-membership anomalies. Neither cell thresholds, object
geometry nor centroid memberships are changed to make that diagnostic pass.
Assembly and release checks revalidate the passing lineage report against the
actual candidate inputs.

The endocrine denominator contract is `named_four/1`: Beta, Alpha, Delta and
CHGA+ hormone-negative. All dependent summaries and embeddings are rebuilt.
Objects with zero named endocrine core cells retain structural records with
unavailable composition. Unclassified cells remain in overall counts,
phenotype controls, plots and downloads; palette lookups now preserve their
names so their grey legend entry is not silently lost.

## Production and delivery

The production release is being prepared as `pubapps-2026-11-components-v1`.
The parallel security session owns image scans, installed-runtime checks,
host staging, binding, activation and proxy verification. See the authoritative
[security and deployment record](container_security_20261102.md). A build or
copied artifact is not evidence that the site is live.

The separate development DuckDB fix remains an independent build. The owner's
dated production DuckDB exception expires at the recorded November 10 UTC
deadline; all other security, scientific and runtime checks still apply.

The new implementation groups are `ade7c46`, `d946fab`, `9d9adc9`, `b3be57a`,
`c93eebf`, `67ee6c9` and `e5689f6`. Deployment recovery fix `caf79f2` preserves
a previously stopped service on activation failure. A new immutable Senior /
monorepo / HPG / MultiPlex5000 supplement will carry these changes and final
evidence; prior delivered packages remain intact.
