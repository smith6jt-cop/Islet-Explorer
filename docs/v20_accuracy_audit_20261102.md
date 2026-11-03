# v20 accuracy audit — 2026-11-02

**Subsequent selected-component build, November 2:** the completed v2 app and
QuPath builds use scientific source `67ee6c9`; the matched image/release builder
uses `e5689f6`. Production `pubapps-2026-11-components-v1` is now live and
verified (2026-11-03 02:42 UTC). The separate strict-DuckDB dev release is also
live, verified at 12:42 UTC; see [its receipt](dev_duckdb_deployment_20261103.md). Read [the v2 build record](component_build_v2_20261102.md),
[resolution contract](component_resolution_20261102.md), and
[actual deployment status](container_security_20261102.md) for current work.
The version-specific source paths, failed diagnostics and prior service state
below remain historical evidence; they do not reopen the delegated selection
or production decision.

**Historical frozen-v20 status: not ready for production.** The audit found an app immune-list defect and an unresolved upstream typing issue. Production remains stopped. Neither a passing reconciliation nor a passing browser suite is scientific acceptance.

This records the frozen v20 audit. The subsequent selected-v19 component app,
QuPath and matching container builds are complete; their passing export checks
and remaining endocrine failures are recorded in the
[November 2 build record](component_build_20261102.md). Those newer results do
not rewrite this audit or establish replacement identities.

**Latest follow-up:** [the classifier-to-cell discrepancy trace](endocrine_classifier_discrepancy_20261102.md) establishes where original pixel detection diverges from corrected whole-cell typing. All 20 classifier masks replay exactly; 7,413 of the 7,800 unresolved cells contain original above-threshold pixels. All 1,468 islets contain positive cells; 910 single objects and 38 clusters have no positive eligible member. The new level-selection gate rejects this inconsistency. The earlier retype/exclude question is superseded by the owner's instruction to resolve the classifier discrepancy.

## Verified evidence

The audited candidate is `pubapps-2026-11-uniform-v20`, app commit `c974946`, image `sha256:517e969d2333744f609533b8f4a7fd217ce8a852fe30ef3551e0d805c2a0db00`. Its deployment `MANIFEST.sha256` hashes to `d9a29beda767eeb3beab57f6a5eb67836229dd8d9807c27208bd646366398d7c`; its separate `data/manifest.json` hashes to `81259bcaedb878f3039150006e18aded71abd8d9469090dae1064ad371bdec55` on both the workstation and VM. Live dev state and image identity were read from the VM; production still names v10-perf2 and is failed/stopped.

Evidence and runnable independent audit: `output/v20_accuracy_audit_20261102/` (local, not committed).

| Check | Result and limits |
|---|---|
| HPG final source identity | All 20 identity Parquet hashes, positivity levels and DG structure binding match the v20 assembly manifest. HPG is authoritative; the workstation processing directory predates DG. |
| Accepted geometry records | All 120 review records and 859 referenced render hashes verified. Every donor has Tissue, filtered augmentation, Nerve, Capillary, DG Lymphatic and ColIV records. This verifies evidence integrity; it is not a new visual review of all those classes. |
| Final reconciliation to assembly | All 22 assembled files pinned by the upstream reconciliation match their hashes. Historical reconciliation covers 3,525,977 detections: 3,199,606 eligible and 326,371 excluded by QC. |
| Independent donor check | All 3,199,606 unique per-donor cell IDs reconcile to census totals. Labels, display regions, associated objects and coordinates match the donor display CSVs for all 20 donors; donor 78 absent. |
| Independent object check | All 4,848 paired objects: recounting drill-down CSVs exactly reproduces every phenotype core/peri count, immune count, and trajectory count. Proportions and log1p endocrine-normalized features agree within absolute 1e-12. Displayed cells have positive nucleus area. This tests consistency, not whether an identity is correct. |
| Embedding eligibility | 1,398 objects have at least ten core cells; six have zero cells in the assembly's endocrine group. The 1,392 finite-composition objects are embedded. The initial audit assertion omitted the positive endocrine denominator condition; correcting the audit to the implemented documented condition resolves that assertion. |
| Donor metadata | The entire released metadata CSV is byte-identical to `scripts/uniform/donor_metadata.csv`. Age, disease status and sex agree with H5AD observations for all 20 donors; donor 112's age is 3.72. This checks the canonical nPOD-derived table, not a fresh extraction from the original workbook. |
| Endocrine positivity | All four stored endocrine positivity flags exactly reproduce from the current REDSEA measurements and per-image levels across all 3,199,606 eligible cells on HPG. |
| Live dev access | Fresh reviewer-access checks 15/15; browser checks 49/49. Run from the workstation, not a new off-campus acceptance test. Browser checks validate behavior and export availability, not every exported scientific value. |
| App correction tests | 506 test assertions pass, no failures/warnings/skips in testthat; all 35 R sources parse. R emits a package-build-version notice outside testthat. |
| Release definition checks | 53 assertions in `scripts/tests/test_release_parity.R` pass, including seven new regression assertions. A real-data invocation confirms the new check rejects v20's unresolved endocrine identities. |

The fresh full release-parity rerun **passed all 23 checks** against the unchanged release worktree at c974946: complete hashes, prepared summaries, spatial calculations, passthrough files and representative geometry windows. This is the original verifier matching that release's commit; it does not include the new endocrine-definition check. Source/release parity passes while scientific readiness fails.

## Findings and corrections

### Immune gradient omitted a current cell label

`IMMUNE_LABELS` in `spatial_islet_immune.R` retained the earlier `CD45+` naming and omitted `Hematopoietic CD45+`. The v20 upstream immune group includes 50,625 cells with that label. Consequently, the Spatial “All immune cells” gradient disagreed with the upstream immune population despite browser and source/release parity passing.

Local commit `34b6853` includes that label in the gradient, its selector and its Immune display group. A regression test includes CD45-positive, CD45-negative, conflict, T-cell and exocrine populations. A real-data check across 20 donors and 160 distance bins confirms the corrected gradient exactly equals a calculation using the manifest's immune membership. All 160 bin values change; the maximum absolute change in enrichment ratio is 1.1285. This is a computational finding, not a biological interpretation. The correction is **not deployed**.

### Endocrine typing cannot be accepted by changing aggregate labels

The user clarified during this audit: **“If the typing is done properly, there should never be unnamed endocrine definitions.”** The earlier proposed choice between aggregate definitions is superseded by this requirement to investigate typing.

Upstream identity files on HPG contain **7,800 cells** retaining object names; **none** has INS, GCG, SST or CHGA above its current positivity level:

| Retained label | Cells | Recorded reason |
|---|---:|---|
| EndocrineCluster | 1,714 | no_endocrine_marker_above_level |
| EndocrineSingle | 910 | no_endocrine_marker_above_level |
| Islet | 1,159 | no_marker_above_background |
| Islet | 2,656 | shared_marker_only |
| Islet | 1,361 | structure_marker_only |

The registry's `defined_by_object` assigns cluster/single cells to Epithelial/Endocrine before marker-based naming. `identity.py::name_lineage` then keeps the object name when none of its lineage markers is positive. Islet-resident cells without a resolved call retain `Islet` through the general fallback. These are explicit paths, not hormone ties or a lost export field.

The assembly counts the 2,624 unnamed cluster/single cells in its endocrine group, but excludes 5,176 unnamed Islet cells. The app's “All endocrine cells” aggregate uses only the four named types. Within paired cores, the disagreement includes 1,712 cells across 1,300 objects. The two additional cluster cells belong outside the paired-object inventory. Changing UI membership would conceal the unresolved typing issue; no such change was made.

Diagnosis must distinguish inappropriate object assignment, insufficient/different marker measurements, positivity-level selection and genuine non-endocrine cells within objects. No replacement cell type or threshold is justified merely by requiring a named result.

### Measurement and image diagnosis

HPG job **48965615** completed successfully, producing 180 diagnostic cells on 24 sheets across all 20 images. Every sheet was inspected. For each image and unresolved label, the sample contains the cell nearest the endocrine positivity boundary, the median and the furthest below it. This deliberate boundary sample is not an estimate of typing accuracy. Outputs and review notes: `output/v20_unnamed_review_20261102/`.

The white cell outlines often contain DAPI profiles with little bright cytoplasm, beside much brighter hormone signal. Small or narrow outlines recur at the edges of bright patches; some cases have substantial image-pixel coverage despite very low corrected means. The grayscale panels show despeckled image pixels **before REDSEA**; the numerical cell means are **after REDSEA**. Their difference alone does not prove compensation is wrong or identify a replacement cell type.

The reproducible measurement audit (`scripts/research/audit_unnamed_endocrine_measurements.py`, source through `b45a64b`) reads the authoritative HPG identity, raw QuPath measurements, REDSEA measurements, levels and earlier calls:

| Unresolved label | Cells | Previously named endocrine | Raw mean above a current corrected-value level | Zero cytoplasm pixels |
|---|---:|---:|---:|---:|
| EndocrineCluster | 1,714 | 556 | 735 | 371 |
| EndocrineSingle | 910 | 730 | 150 | 91 |
| Islet | 5,176 | 1,486 | 2,776 | 881 |
| Total | 7,800 | 2,772 | 3,661 | 1,343 |

All 7,800 are negative for the four endocrine markers under the current corrected-value levels. Their median nucleus/cell pixel fraction is 0.783 and median cell area is 37.46 µm². Comparing raw means to a corrected-value threshold is only a diagnostic: it is not a calibrated alternative positivity method. Earlier endocrine calls also used a superseded measurement/threshold method and cannot simply be restored. The evidence does not justify discarding all 7,800 cells or forcing hormone names.

Use `output/v20_unnamed_measurements_20261102_r2/` for this first measurement report. The first report's named-cell geometry comparison used the old hormone spelling; r2 corrects that comparison to all 113,116 currently named endocrine cells. The unresolved-cell totals above were unaffected. Thirty returned diagnostic artifacts/scripts match their HPG SHA-256 values (`HPG_DIAGNOSTIC_SHA256SUMS`). At that stage only the two diagnostic scripts were sent to HPG. The subsequent classifier trace used committed source updates and an explicit changed-file plus `.git` transfer; its separate provenance, native-raster measurements and 68 returned artifact hashes are recorded in the follow-up report. Existing scientific outputs were not overwritten.

### New release check

Commit `6cee485` adds `validate_release_label_groups()` to the builder and parity verifier. A root declaring final-typing label groups cannot be newly released while it carries `EndocrineCluster`, `EndocrineSingle` or unresolved `Endocrine` identities, counts `Islet`/`Tissue` as endocrine, or disagrees with the app's endocrine aggregate. Removing cluster/single cells only from the denominator still fails. This enforces the clarified requirement without assigning new identities. Historical roots without declared label groups retain their existing contract; passing this check alone is not scientific acceptance.

## Geometry review scope

The old raw augmentation statuses (`needs_review`, `quality_hold`, and `not_reviewed_superseded_before_consumption`) describe superseded raw geometry. The final reconciliation explicitly validates the filtered generation's reviews and hashes. `segmentation_quality_hold.resolution.json` resolves the hold only for that filtered generation; it does not release old raw generations.

T1D_6533 has two explicitly documented inherited baseline geometry exceptions. Its reconciliation records `primary_geometry_checks_passed: false` alongside `inherited_baseline_exception`; do not summarize this as every geometry passing without exception.

All 20 existing DG lymphatic sheets were visually rechecked in this audit (whole field plus six sampled fields each). Detailed per-image observations are in `visual-recheck.json`. No global registration shift was apparent at sheet scale. Variable background, retained puncta, and broad PDPN-positive walls remain visible; the sheets do not establish biological lymphatic specificity. Other geometry classes have hash-verified earlier reviews but have not yet received a complete new visual pass in this audit.

## Remaining work before promotion

1. Reconcile cell selection and corrected-mean positivity with the original classifier evidence, using the [completed causal trace](endocrine_classifier_discrepancy_20261102.md). The earlier general-typing/exclusion question is withdrawn: a negative corrected mean does not establish absent original endocrine signal. The new level-selection gate detects fully negative objects; the separate release gate still detects unnamed members of otherwise positive objects. No replacement identity or threshold has been established.
2. Apply only a scientifically justified correction; if typing or measurements change, run the required image review and HPG workflow into a fresh candidate, then regenerate assembly, embeddings and release data.
3. Complete the remaining geometry/measurement and UI/export accuracy checks; do not infer their completion from hash parity.
4. Obtain the owner's QuPath/dev acceptance of the corrected candidate. Then build and stage its immutable image/data release with the immune fix, run full parity and browser/access acceptance, and promote through the established deployment procedure.
5. Verify production and retain rollback. Off-campus acceptance and the UFRC reboot test remain separate outstanding operational checks.

No production activation, data deletion, Git push or hypothesis analysis occurred in this audit.
