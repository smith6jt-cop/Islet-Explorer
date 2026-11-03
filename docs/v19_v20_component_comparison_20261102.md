# V19 versus v20: component comparison — 2026-11-02

Owner: “I meant use v19 for the endocrine cells since they were done properly. compare v19 with v20 to see which has the better of each component.”

The task is to evaluate components separately, using v19 as the selected endocrine source. A wholesale rollback to v19 was the assistant's mistaken interpretation. The prepared whole-v19 release is retained as reference evidence only; no new scientific release was activated. Dev last served v20; dev and production are now stopped and masked under the owner's parallel [security direction](container_security_20261102.md). The mixed candidate is now built at `b15f0fc`; its actual endocrine consistency checks still fail. See [the completed build record](component_build_20261102.md). The main-tree exception that allowed unresolved v19 identities into a new release has been withdrawn.

**Subsequent resolution:** the owner delegated the remaining decision and live
release. The [resolution contract](component_resolution_20261102.md) implements
exact named v19 endocrine membership with current non-endocrine methods and
explicit unsupported outcomes. The [v2 build record](component_build_v2_20261102.md)
contains current validation; later status supersedes pending-selection language
in the historical comparison below.

## Comparison contract

* Compare frozen `uniform-20donor-20261029-v19` with frozen `uniform-20donor-20261030-v20` by exact cell UUID, donor, geometry and source hashes.
* V19's matching typing is archived in `uniform-20donor-20261027/superseded/typing-v19-20261030`. Final v20 processing is authoritative on HPG. The local working `v19-qc-work` directory is a stale intermediate and must not stand in for final v20.
* A label-count increase, fewer conflicts, or a passing assertion does not independently prove better cell identity. Distinguish compliance with an explicit owner rule, quantitative consistency and image-supported accuracy.
* Preserve v19 endocrine as the owner's selection; compare its actual calls and measurements rather than presuming every discrepancy is fixed. No arbitrary new labels, positivity thresholds or geometric changes follow from selecting a version.
* Any eventual combination must assign each UUID one identity, resolve cross-component overlaps explicitly, and recompute summaries, normalization, embeddings and spatial outputs from the selected identities. Mixing already aggregated tables is not a reproducible correction.

Evidence is being assembled under `output/v19_v20_component_comparison_20261102/`. The completed source audit for the frozen roots remains in `output/v19_selection_20261102/source_audit_r2.json`.

## Verified starting point

Both frozen roots contain the same 3,199,606 eligible UUIDs and 20 donors, with donor 78 absent. Exactly 509,408 cell labels differ. The two versions retain the same endocrine-object geometry. V19 contains 7,151 unresolved cells in the Islet/EndocrineCluster/EndocrineSingle labels; v20 contains 7,800. These figures are cells, not counts of negative islets. V19 has 953 fully negative populated objects (915 singles and 38 clusters), and v20 has 948 (910 singles and 38 clusters); all 1,468 islets in each version contain positive member cells. These measurements do not establish which individual changed call is correct.

## Component findings

| Component | Evidence and selection |
|---|---|
| Endocrine cells | **V19, as selected by the owner.** Its 110,721 named cells are the baseline to preserve. This is a component selection, not acceptance of every unnamed cell in the archive. All 80 endocrine thresholds changed in final v20 despite unchanged measurements and endocrine geometry. |
| Tissue, cell segmentation and QC universe | Identical between the frozen versions: 20 Tissue geometries, 3,199,606 eligible UUIDs and their coordinates/regions; neither version gains an advantage from these unchanged inputs. The prior reviewed geometry exceptions remain recorded. |
| Endocrine objects | Identical: 1,468 islets, 3,385 clusters and 5,227 singles. A change in cell label is not a changed endocrine object. |
| Nerve, Capillary, Collagen IV geometry/endpoints | Identical component geometries and 10,100 endpoint rows per class; Collagen IV distance Parquets are byte-identical. Retain these common inputs. |
| Lymphatic geometry | V20 DG implements the later explicit despeckling and PDPN-gate instructions. It is the selected method to retain, with its own regenerated context and endpoints. Smaller counts/area alone do not establish accuracy. |
| Neural and endothelial cell rules | V20 implements the explicit structure-only rule. V19 also assigned these identities through the marker-only route. Retain the corrected rule, while isolating calibration from Lymphatic changes. |
| T cells, macrophages, exocrine cells | V20 implements whole-cell positivity with membrane differentiation for CD3e, CD68, CD163 and ECAD. All 80 membrane thresholds are unchanged from v19. Retain the instructed measurement separation; the frozen v20 labels still contain unintended calibration effects and cannot simply be copied as the better component. No held-out image accuracy comparison establishes a universal winner among changed calls. |
| B cells | Same approved CD20/CD45/single-lineage rule. The net four-cell difference is not an accuracy measure; competing markers can change the result. |
| Ki67 state | V20 uses the instructed nuclear measurement. Both frozen app exports omit Ki67/proliferation fields, so neither supplies this requested state to the app. Selecting the upstream method alone does not repair that export gap. |
| Shared calibration and canary competition | Must be bound independently of a structure revision. These shared inputs couple otherwise separate components; independently selecting larger or smaller phenotype counts is invalid. |

The exact method review is in `output/v19_v20_component_comparison_20261102/methods/component_method_review.json` and `component_verdicts.csv`. Geometry evidence is in `geometry/summary.json` and `interpretation_and_constraints.json`.

The app's Rules modal also exposed an inherited documentation mismatch: identity-generation roots carry an older `protein_phenotype_rules.csv` describing cohort scaling, a twofold hormone margin and marker-based fallthrough. Commit `8be9391` renders the actual v19/v20 method descriptions from their manifest metadata instead, including their differing measurement and structural-cell rules. Real v10/v19/v20 render checks pass; v10 HTML is unchanged. This is a display-text fix, not a scientific data edit.

## Isolating the lymphatic regression

The frozen v19→v20 comparison includes deliberate typing-policy changes as well as Lymphatic changes. To isolate the latter, the audit compared **initial v20 option-1 with final v20 DG**, whose typing policy, identity, levels, canary and registry hashes are identical. The local working directory matches both completed option-1 journals; the final clone matches the DG journal. `geometry/lymphatic_generation_locator.json` records the exact sources and hashes.

Across all 3,199,606 matched eligible UUIDs, this controlled contrast found:

* 205 of 540 positivity thresholds changed, across 15 markers. INS, GCG, SST, CHGA, CD31, CD34, LGALS3 and VIM changed in every donor.
* 265 of 520 canary thresholds changed. Some immune thresholds changed even where plain positivity thresholds stayed fixed. The 80 differentiation rows, 80 display rows and Ki67 flags did not change.
* 3,341 cells changed at least one positivity/canary flag; 715 changed their cell label.
* 352 label changes involved a Lymphatic label; 363 did not. Of all label changes, **303 were outside both old and new Lymphatic centroid/positive-overlap footprints**. Forty-three involved a named endocrine label.
* REDSEA inputs, eligible UUIDs, non-Lymphatic memberships and non-Lymphatic geometries remained unchanged.

These are unintended computational dependencies, not evidence of biological change or an accuracy estimate. Evidence: `methods/lymphatic_controlled_contrast/summary.json`, exact changed UUIDs and transition tables. The 509,408 full-version label changes must not all be attributed to DG.

The cause is shared reference selection. `levels.load()` derived `ordinary` and structural reference populations from the current `typing/context`; canary calibration reused those populations. Regenerating Lymphatic context therefore recalculated unrelated markers' levels, folds and competition. Separately, the structure step protected clusters/singles but could override named endocrine cells inside islets.

## Correction and verification

`typing/calibration_context.py` explicitly snapshots the selected baseline's complete reference context, records its hashes and provenance, and validates UUIDs, coordinates, QC and endocrine membership against live inputs. Level/reference calculations read that snapshot; identity assignment still reads current structure membership. A selected structure generation without an explicit baseline now stops instead of silently choosing new calibration cells.

`typing/structure_update.py` binds a fresh Lymphatic revision to unchanged positivity, canary, differentiation, display and selection tables, the registry, REDSEA inputs and baseline identity hashes. It never chooses v19/v20 implicitly. The levels and canary CLIs refuse recalibration of a bound structure-only run. Context regeneration requires that binding. Identity assignment protects baseline named endocrine cells from structure overwrite and rejects any endocrine membership/subtype change rather than covering inconsistent evidence with copied labels.

The audit script `scripts/uniform/audit_structure_update.py` checks every marker's positivity, canary and fold fields, coordinates/QC, non-Lymphatic context and calibration, as well as changes outside the old/new Lymphatic footprint. Changed labels must have a Lymphatic label on one side. The identity CLI automatically runs this audit and fails if it fails. A passing isolation audit is necessary for the correction; it is not scientific release acceptance.

The final all-donor diagnostic replay **passed** in `methods/lymphatic_isolated_replay_v2`, with the automatic report in `methods/audit_isolated_replay_v2/summary.json`. It uses initial-v20 calibration with verified DG live context to test isolation:

* All 20 REDSEA Parquets and all nine calibration/selection/options/null CSVs are byte-identical to the baseline; 3,525,977 detections and 3,199,606 eligible UUIDs are preserved.
* Every marker positivity, canary and fold field, Ki67 state and non-Lymphatic context field is unchanged.
* All **113,134 initial-v20** named endocrine labels, compartments and lineages are preserved. This diagnostic count is distinct from the 110,721 named endocrine cells selected from v19.
* Exactly 352 label changes remain; each involves `Endothelium Lymphatic`, lies within the old/new Lymphatic footprint and has changed Lymphatic geometry. There are zero unexpected changes and zero non-Lymphatic-to-non-Lymphatic label changes.
* All 20 final identity Parquets reproduce the first isolated replay byte-for-byte. An unchanged-context control reproduces the original ND_0112 identity Parquet exactly (175,978 rows, 104 columns).
* Both calibration CLIs refuse a bound run without writing files. The historical DG contrast fails the audit as expected. No synthetic scientific cases were introduced.

**This is not the selected v19-endocrine mixed candidate, and is not an app data root.** Fresh authoritative processing and reconstruction of all derived app artifacts remain separate from this diagnostic. The compact tracked evidence receipt is `docs/evidence/lymphatic_isolation_20261102.json`; full per-cell evidence remains in the ignored output tree.

## Endocrine evidence and remaining release work

The paired review retains unchanged native image pixels and cell boundaries, with both versions' corrected means, thresholds and labels. All 14 subtype-switch cells and 39 representative gains/losses have paired panels under `visual/paired_v2/`. All 14 subtype-switch panels were inspected; six sheets also received rendering QA. Brightness cannot independently decide the winner because the hormone rule compares per-marker normalized scores. Eight switches cross a positivity boundary; six retain the same positivity flags but change the normalized winner.

The final verified REDSEA values reproduce both versions' stored endocrine positivity and float32 folds exactly. Replacing only reference membership reproduces all 80 endocrine thresholds. There is no evidence here that DG changed the underlying endocrine measurements.

Selecting v19 does not erase the inherited unresolved-object problem: its 915 singles and 38 clusters still lack a positive eligible member under its corrected-cell thresholds. No subtype is invented from a negative corrected mean. No new release may hide those cells by broadening endocrine totals or reuse inconsistent summaries. A coherent candidate must preserve the chosen v19 endocrine decisions, reconcile remaining component boundaries by UUID, carry required state fields, rebuild downstream outputs, and pass the existing object and release checks. No scientific release was activated; both services are stopped and masked under the separate security-remediation direction.

## Completed selected-component build

The fresh mixed-calibration control and DG candidate preserve all **110,721 named
v19 endocrine cells**. Their all-20 invariance audit passes: **351** Lymphatic-only
label changes within changed geometry, zero unrelated changes, all marker/state
and non-Lymphatic context fields preserved. This count belongs to the new mixed
component comparison, not the earlier 352-change initial-v20 diagnostic.

The app root `data/development/uniform-20donor-20261102-component-app-v1` and its
20-image QuPath review project are built. All UUID/label/Ki67 export checks pass
for 3,199,606 cells; summaries, scales, spatial outputs and embeddings were rebuilt.
The completed artifact still has 8,629 unnamed object-class labels and 2,472
boundary conflicts outside the exact selected v19 endocrine membership. Its
corrected-cell object audit retains the 953 fully negative populated objects and
three Single membership anomalies. The release checks reject those actual
failures. See [the build record](component_build_20261102.md) for paths, source
hashes, image/security results and verified return evidence. No new release was
activated; production activation is already authorized once actual checks pass.
