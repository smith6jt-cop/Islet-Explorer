# Islet Explorer: workflow user guide

This guide follows the data from the source images to the reviewed app and its release. It is for the owner, analysts, and operators who need to know which inputs are being used, what each stage changes, and what must pass before the app goes online. For the app controls themselves, use the [reviewer guide](reviewer_guide.md). For scientific decisions, the [typing rulings ledger](typing_rulings_ledger.md) and latest owner directions in [CLAUDE.md](../CLAUDE.md) take precedence over historical examples.

## Current state — 2026-11-02

**The mixed-component app and QuPath builds are complete; actual endocrine consistency checks still fail. Both dev and production remain stopped and persistently masked.** The owner selected **v19 for endocrine cells**, with the other components evaluated separately. This is not approval to replace the entire app with v19.

| Item | Current status |
|---|---|
| Endocrine component | The 110,721 named endocrine cells in frozen v19 are the selected baseline. Its unresolved cells and object-membership anomalies remain explicit problems. |
| Lymphatic component | The instructed v20 DG method uses despeckled PDPN detection plus the PDPN level gate. Its geometry must not recalibrate other components. |
| Other components | Keep the common geometry and the instructed v20 measurement/structure rules, with component boundaries reconciled by cell UUID. Frozen v20 labels cannot simply be copied because they include unintended calibration changes. |
| Isolation fix | The all-20-image diagnostic replay passes: 352 label changes, all involving Lymphatic identity inside the old/new Lymphatic footprint; no unexpected changes to other markers, identities, QC, or non-Lymphatic context. It uses initial-v20 calibration and preserves 113,134 initial-v20 named endocrine cells. It is **not** the selected v19-endocrine candidate. |
| Completed component build | HPG job `49141854`, source `b15f0fc`: all eight stages pass. The selected 110,721 named v19 endocrine cells are preserved; the mixed-control→DG audit finds 351 Lymphatic-only changes and zero unexpected changes. Export checks pass for all 3,199,606 cells, including Ki67. |
| Remaining data failures | 8,629 unresolved object-class cell labels, 2,472 component boundary conflicts; the strict corrected-cell object audit fails for 953 populated objects and three Single membership anomalies. These are actual failures, independent of the old diagnostic status flag. |
| Current matching image | `41e62a97…` from `b15f0fc` passes package, dependency, installed Shiny/R/static-image and dated security-policy checks. Built offline; not exported or staged. See the [build record](component_build_20261102.md). |
| Separate dev DuckDB fix | `b26558a9…` from `92c63ba` passes the strict decoder/security gate without a deferral and real release-table BROTLI round trips. Exported locally, not host-staged or activated. See the [dev validation record](dev_duckdb_fix_20261102.md). |
| Dev site | `https://devislet-explorer.rc.ufl.edu/islet-explorer/` last served `pubapps-2026-11-uniform-v20`. The owner stopped dev during container security remediation on 2026-11-02; it is now masked. |
| Production | `https://islet-explorer.rc.ufl.edu/islet-explorer/` has been stopped since 2026-10-26 and is now also masked. Its recorded release is `pubapps-2026-10-uniform-v10-perf2`. |
| Local app service | `islet-explorer-local.service`, port 3840, is inactive/disabled. |
| Offline security image | `2bbb928d…` from `23bd2fd` is built and imported under a hold on PubApps. Released fixes pass package/dependency checks. Its original strict gate failed DuckDB; fresh dated deferral evidence now makes that known finding nonblocking within the patch window. A matching validated data release is still required. |
| Whole-v19 reference release | The build at `190addc` was staged but never activated. It remains reference evidence; its historical compatibility exception is withdrawn from current source. |

The [component comparison](v19_v20_component_comparison_20261102.md) records the selections and evidence. The [compact isolation receipt](evidence/lymphatic_isolation_20261102.json) binds the verified source and reports. The [container security record](container_security_20261102.md) documents the service masks, affected images, and replacement-image work. The owner has now directed the current build to live production once its real data/release/runtime checks pass, with our DuckDB fix developed separately for dev. The known DuckDB finding is nonblocking through the notice's **2026-11-10 UTC** deadline under the [dated release decision](duckdb_remediation_window_20261102.md). Do not infer that a committed fix, successful diagnostic, or staged bundle has changed served data.

## 1. Find the correct inputs

Paths below are relative to `/data/panc_CODEX/Islet-Explorer` unless marked otherwise.

| Input or output | Location and purpose |
|---|---|
| Cohort | `scripts/uniform/cohort.csv`; 20 donors. `scripts/excluded_donors.csv` excludes donor 78 throughout the app and analysis. |
| Donor metadata | `scripts/uniform/donor_metadata.csv`, derived from the nPOD export by `scripts/build_donor_metadata.py`. Donor 112's age is 3.72. Do not carry an older app age forward. |
| Source QuPath project | `/data/QuPath_Projects/QuPath_project_IsletExplorer`; cohort records identify source entries, image files, and Tissue outlines. |
| Frozen v19 app root | `data/development/uniform-20donor-20261029-v19` |
| Matching v19 typing archive | `data/development/uniform-20donor-20261027/superseded/typing-v19-20261030/typing`; use its hashes, not a similarly named current working folder. |
| Frozen final-v20 app root | `data/development/uniform-20donor-20261030-v20` |
| Authoritative final-v20 processing run | **On HPG:** `data/development/uniform-20donor-20261029-v19-qc-work` under the HPG repository root. The workstation directory with the same relative name is an earlier option-1 run. |
| Verified final-v20 diagnostic clone | `output/v19_v20_component_comparison_20261102/endocrine/authoritative_run`; its transfer receipt verifies the copied identity/context/REDSEA/levels/object files. It is a partial diagnostic source, not a complete processing run. |
| Comparison evidence | `output/v19_v20_component_comparison_20261102/`; source-selection evidence is in `output/v19_selection_20261102/`. |
| Frozen QuPath reviews | `data/development/qupath-review-20261029-v19/project.qpproj` and `data/development/qupath-review-20261030-v20/project.qpproj` |

Both frozen app roots contain the same 3,199,606 eligible cell UUIDs from 3,525,977 detections and the same 20 donors. Exact UUID-to-label and coordinate checks, source hashes, donor tables, and image inventory are recorded in `output/v19_selection_20261102/source_audit_r2.json` and `donor_scope_audit.json`. The original report is retained beside the revision. These checks establish provenance and parity, not that every cell call is correct.

Before starting a new run, record the selected input roots, resolved symlink targets, image hashes, cohort and metadata hashes, code commit, model/configuration, and intended component change. Read `uniform_release.json` in an app root and the processing run's provenance instead of treating its folder name as a version guarantee. Keep source images, archived typing, accepted geometry, and prior app roots unchanged.

## 2. How the stages fit together

```text
Source images + donor/cohort records
    → reviewed Tissue and exclusions
    → cell/nucleus segmentation + measurements
    → endocrine classifier objects + reviewed structures
    → QC eligibility + despeckled REDSEA measurements
    → spatial context + explicit calibration reference + levels/canary
    → one identity per eligible cell + state fields + consistency checks
    → selected components reconciled by UUID
    → fresh QuPath review and complete app assembly
    → immutable schema-5 release + source/release parity + browser checks
    → dev review → accepted production activation → verified rollback path
```

There are two different kinds of change. A new segmentation or scientific typing method needs its own reviewed processing generation and downstream rebuild. A **Lymphatic-only revision** reuses the selected measurements, calibration, and unrelated identities, and changes only the reviewed Lymphatic geometry and its permitted consequences. Section 5 describes that separate path. Running the old all-marker typing driver after a structure edit is not a structure-only update.

Commands below use the repository root as the working directory. `RUN`, `BASELINE`, `ASSEMBLY`, and report paths must be set to explicitly selected inputs and fresh output locations. They are not defaults selecting a scientific version. Processing commands belong in allocated HPG jobs; see section 10. The workstation interpreter used for lightweight checks is:

```bash
PY=/home/smith6jt/miniconda3/envs/scvi-env/bin/python
```

## 3. Review images, geometry, segmentation, and QC

Every geometry that gates data must be inspected against the image on **every image** before acceptance: Tissue, cell/nucleus segmentation, endocrine and structure objects, exclusion masks, and edge regions. Render a whole-image view and native-resolution tiles with DAPI and the relevant stain; actually inspect them. Quantify coverage against an image-derived reference, part/hole counts, perimeter, cell density, and objects touching boundaries. Record outliers and their disposition. A reused file, passing script, or area-only report is not an image review.

This rule follows the discovery that earlier Tissue outlines omitted substantial tissue from some additions. The corrected outlines and prior review notes are historical inputs to verify, not permission to skip review on a new generation. The original images use 0.3774 µm pixels; the classifier planes in these runs use 0.7548 µm pixels. Keep coordinates, origins, downsampling, and units explicit in overlays.

For a new, initialized full processing generation, the existing entry points are:

```bash
"$PY" scripts/uniform/make_project.py --out "$RUN"
"$PY" scripts/uniform/run_uniform.py --out "$RUN" --images all
```

These are full-generation entry points, not instructions to rerun a frozen root. `run_uniform.py --help` lists the selectable steps. It imports Tissue, checks the project, runs InstanSeg and cell measurements, exports classifier planes/cell outlines, fits the endocrine classifier, makes and names objects, calculates signed distances and parenting, then exports final cells/context/annotations. All cohort cell measurements precede the first fit. Each step records its inputs and outcome under `RUN/provenance/<image>/`; dependencies become stale when their inputs change. Do not use `--force` to bypass a source, geometry, or acceptance discrepancy.

The recorded segmentation for these generations is InstanSeg `fluorescence_nuclei_and_cells` 0.1.1 with DAPI, CD68, CD44, VIM, ECAD, LGALS3, panCK, CD3e, and BActin. This is the selected historical configuration, not a claim about the newest downloadable model. Changing it is a new scientific generation.

`scripts/uniform/tissue_qc.py --out "$RUN" --report-dir "$FRESH_REPORT"` produces Tissue QC; its optional `--render` takes image names. Additional geometry-specific review is still required. The [uniform methods history](../scripts/uniform/README.md) describes the source project and exports; its older typing commands and thresholds are superseded where the rulings ledger says so.

QC determines the eligible analysis universe before identities and app summaries are built. Preserve every detection and its exclusion reason in the processing evidence, while keeping excluded cells out of analysis tables. Preserve raw UUIDs and coordinates across comparisons. Review acquisition-edge truncation and cell/nucleus eligibility through the existing generation-specific QC records; do not silently add an exclusion to make a count reconcile.

## 4. Measurements, calibration, and typing

The pixel classifier and the later cell typing answer different questions. The classifier identifies endocrine signal regions from INS/GCG/SST with CHGA support on blurred raw pixels. Whole-cell object construction can add cells through overlap and morphology. Typing then evaluates despeckled, REDSEA-corrected cell measurements against image-specific cell-level thresholds. A corrected cell mean of zero does not show that the original image lacked endocrine pixels.

The [classifier discrepancy report](endocrine_classifier_discrepancy_20261102.md) reproduced all 20 original classifier masks exactly. Of final v20's 7,800 unresolved cells, 7,413 contain original above-threshold pixels. Therefore neither exclusion nor an invented subtype follows from their negative corrected means. Pixel cutoffs must not be substituted for cell-mean thresholds.

The current measurement contract separates:

| Quantity | Role |
|---|---|
| Raw QuPath cell/nucleus measurements and boundaries | Preserved source evidence and coordinate/QC checks. |
| Despeckled REDSEA values | Cell measurements used by final typing; the 3×3 median and spillover correction are part of the method. |
| `typing/levels.csv` | Per-image positivity on the marker's selected measurement. The instructed rule is whole-cell positivity, with Ki67 read on the nucleus. |
| `typing/levels_differentiation.csv` | Membrane differentiation for CD3e, CD68, CD163, and ECAD when groups compete for a cell. It is distinct from their whole-cell positivity. |
| `typing/canary_levels.csv` | Competing-marker/background checks; these share calibration reference populations and must also be isolated from structure changes. |
| `typing/display_levels.csv` | Display-only calibration for markers outside the typing decision. |
| `typing/context/<image>.parquet` | Current geometry membership and overlap fields. In a bound structure revision, calibration uses a separately pinned full context snapshot. |
| `typing/identity/<image>.parquet` | One final eligible-cell record with identity, evidence, positivity/fold fields, and states. |

The protected `scripts/uniform/typing/redsea_cells.py` is the measurement implementation. For an authorized fresh full generation, its CLI is `--out RUN --images all --workers N`. Changing its source or rerunning it in an accepted root is not part of a structure revision.

For an **unbound new typing generation**, the current calibration and typing entry points are below. This is a stage map: each command requires the reviewed upstream artifacts and must stop if its preceding checks fail.

```bash
"$PY" scripts/uniform/typing/cell_context.py --out "$RUN" --images all --workers 6
"$PY" scripts/uniform/typing/levels.py --out "$RUN" --workers 20
"$PY" scripts/uniform/typing/canary.py --run "$RUN" --workers 20
"$PY" scripts/uniform/typing/area_coverage.py --run "$RUN" --workers 5
"$PY" scripts/uniform/typing/identity.py --run "$RUN" --workers 20
"$PY" scripts/uniform/typing/acceptance_checks.py --run "$RUN" --workers 20
"$PY" scripts/uniform/typing/display_levels.py --run "$RUN" --workers 10
"$PY" scripts/uniform/typing/structure_metrics.py --run "$RUN" --workers 3
```

`levels.py --report FRESH_DIRECTORY` can write a separate candidate calibration report. Its method-selection tables retain failed choices and reasons. A recorded level is not proof the method passed. The new endocrine object audit runs after level construction and causes failure if a populated endocrine object has no INS/GCG/SST/CHGA-positive eligible member, or if object membership is malformed. To audit existing inputs without changing them:

```bash
"$PY" scripts/uniform/typing/endocrine_object_audit.py \
  --run "$RUN" --out "$FRESH_OBJECT_REPORT"
```

Exit code 4 means the reported scientific consistency checks failed. The check does **not** require every cell within an islet to be endocrine-positive. The separate release-label check also detects unresolved endocrine definitions inside otherwise positive objects. Neither check repairs the inputs by relabeling cells.

## 5. Make a Lymphatic-only revision without changing other components

Use a fresh candidate run with explicitly selected baseline cell measurements, QC, identity evidence, and calibration tables. Prepare its unchanged inputs and initial context from that baseline, with an empty identity output directory. Retain the baseline's marker-area coverage with its verified provenance: it depends on cell pixels and fixed levels, not Lymphatic membership. Select the reviewed Lymphatic geometry using its existing hash-pinned `structure_generation.json` ancestry. `scripts/uniform/structure_generation.py` validates the reviewed generation and requires Nerve, Capillary, and Collagen IV geometry to remain unchanged. A binding is a generated provenance record, not a JSON file to hand-edit to force acceptance.

Writable outputs must be independent copies/directories. Do not link new context, identity or endpoint destinations into a historical root: a writer could follow the link and replace the baseline. Immutable input links are appropriate only where validation proves their content unchanged.

The DG generation and its original detection/review drivers are documented in [the typing handoff](typing_handoff_20261031.md) and [HPG processing record](hpg_processing_20261031.md). Reuse the exact verified geometry and review ancestry when it is unchanged. A new geometry requires new image review and records. Do not rerun the historical `output/typing_20261030/run_v20.py` all-marker calibration path for this purpose.

Before regenerating live context, bind the fresh run explicitly:

```bash
"$PY" scripts/uniform/typing/structure_update.py \
  --run "$RUN" --baseline "$BASELINE" --images all \
  --scientific-baseline "$BASELINE_DESCRIPTION" \
  --selection-reason "$SELECTION_REASON"
```

This validates identical calibration and REDSEA inputs, records baseline identity hashes and the current registry, and copies the complete baseline reference context to read-only snapshots. `typing/structure_update.json` records the structure-only contract; `typing/calibration_context.json` records the reference snapshots, SHA-256 hashes, and source provenance. UUIDs, coordinates, QC, endocrine membership, and normal context schema must agree. Raw-cell and REDSEA inputs are checked where available. An existing pin cannot be silently rebound.

Now regenerate only the allowed live structure consequences and identities:

```bash
"$PY" scripts/uniform/typing/cell_context.py --out "$RUN" --images all --workers 6
"$PY" scripts/uniform/typing/structure_metrics.py --run "$RUN" --workers 3
"$PY" scripts/uniform/typing/identity.py \
  --run "$RUN" --workers 20 --structure-audit "$FRESH_STRUCTURE_REPORT"
```

In this mode, `levels.py` and `canary.py` refuse recalibration. Reference readers validate and use the pinned context; structure membership remains live. Named baseline endocrine cells inside islets receive the same protection as those inside clusters/singles. Identity assignment fails if an endocrine decision changes; it does not cover inconsistent evidence by copying old labels onto new results.

The identity CLI automatically runs the all-component invariant audit. The standalone read-only check is:

```bash
"$PY" scripts/uniform/audit_structure_update.py \
  --baseline "$BASELINE" --candidate "$RUN" \
  --output "$FRESH_STRUCTURE_REPORT" --images all --workers 4
```

Choose a new report directory if identity already created one. The audit checks every marker's positivity/canary/fold fields, calibration, REDSEA, UUIDs/coordinates/QC, endocrine decisions, and non-Lymphatic context. Label changes must involve Lymphatic identity and remain within the old/new Lymphatic footprint. A failure stops downstream acceptance; a successful job exit without the audit report is insufficient.

The completed real-data check is `output/v19_v20_component_comparison_20261102/methods/audit_isolated_replay_v2/summary.json`. It verifies the code's isolation behavior. It does not select the scientific baseline for the eventual mixed app root.

## 6. Reconcile selected components and inspect them in QuPath

The mixed candidate must start from cell-level evidence. Preserve the selected v19 endocrine decisions and their explicitly selected calibration/measurement contract, apply the instructed non-endocrine methods, and resolve overlaps by UUID so each eligible cell has one identity. Record the origin of every selected component and every changed boundary decision. Recompute downstream tables from that coherent assignment; do not combine version-specific aggregate CSVs or reuse a v20 denominator with v19 labels.

The explicit diagnostic builder is now `scripts/hpg/build_component_candidate.py`,
with preparation in `scripts/hpg/stage_component_build.py` and the component
contract in `typing/component_selection.py`. Its [build record](component_build_20261102.md)
explains the exact source choices and remaining conflicts. `identity.py
--endocrine-from` remains a historical input option; it does not implement this
selection.

The builder requires a JSON config with these exact keys: `expected_commit`,
`selection_reason`, `source_package`, `endocrine_source`, `structure_source`,
`donor_source`, `donor_metadata`, `initial_source_run`, `control_run`,
`candidate_run`, `assembly_root`, and `report_root`. Every source must exist;
every output must be fresh and separate from all inputs and other outputs.
The package is the checksummed initial-v20 calibration/context/identity/coverage
capture, not the current DG run under the same historical workstation path.

Submit from the HPG repository root after explicitly syncing the clean source
commit and verifying the package:

```bash
/opt/slurm/bin/sbatch --job-name=component-build --qos=maigan-b \
  --cpus-per-task=24 --mem=192G --time=06:00:00 \
  scripts/hpg/jobs/site_job.sbatch \
  python scripts/hpg/build_component_candidate.py --config "$BUILD_CONFIG"
```

The serial job verifies and materializes the initial source, prepares a typing
control and a DG candidate, builds both identity sets, runs the all-component
Lymphatic invariance audit, assembles the app root, and validates the exported
cell states. Its `build_journal.json` records commands, host, source commit,
seal, stages and exit codes. A failed stage stops the remaining stages and
retains the failed output; select fresh destinations for another attempt.

The completed November 2 run uses
`tmp/component-build-20261102-v1/build-config-r2.json`; its copied config and
journal are under `output/component_candidate_build_20261102_v1`. The app root is
`data/development/uniform-20donor-20261102-component-app-v1`. Open
`data/development/qupath-review-20261102-components-v1/project.qpproj` for its
20-image review project. Keep these completed outputs immutable; a corrected run
needs new output names and a newly recorded config.

The control preserves the old context as a classification comparison and does
not claim a new geometry acceptance. The DG candidate uses a fresh derivation
binding that fully validates the original accepted geometry and exact shared
input bytes. Historical generation metadata remains unchanged; the new
validator proof is recorded separately. Neither a passing control nor a
completed app assembly proves biological accuracy. The component protocol
preserves inherited unresolved cells and records additional boundary conflicts.
The release builder verifies the evidence and rejects actual unresolved
assignments, conflicting definitions or broken provenance; a historical review
status alone does not add an approval requirement to the owner's live-build
authorization.

The paired evidence for review is under `output/v19_v20_component_comparison_20261102/endocrine/` and `visual/paired_v2/`. It includes exact changed UUIDs, all 14 named subtype switches, and 39 representative gains/losses with both versions' labels, thresholds, and corrected values. Eight subtype switches cross a positivity boundary; six retain the same binary flags but change the normalized winning hormone. Image brightness or a larger named-cell count alone cannot establish a better call.

Build a separate native QuPath project to inspect the candidate identities and
selected geometry, retaining unresolved labels visibly. Its creation is not
scientific acceptance:

```bash
"$PY" scripts/uniform/review_project.py \
  --run "$RUN" --project "$FRESH_QUPATH_PROJECT" \
  --assignments "$FRESH_ASSIGNMENTS" --workers 2 \
  --java-heap-gb 64 --java-processors 12
```

The project and assignments paths must not exist. The exporter verifies every detection is assigned exactly once, imports final eligible identities and exclusions, and includes the selected Nerve, Capillary, Lymphatic, and Collagen IV objects. Inspect `review_manifest.json`, `summary_by_image.csv`, and the per-image import logs. A successful import verifies consistency with the exported assignment; it is not the owner's visual scientific acceptance.

Open the new project's `project.qpproj`, inspect whole-image and zoomed cell/structure overlays, and record findings with image, object UUID, relevant channels, and source hashes. Avoid opening a sealed native input as an editable review: QuPath can update project state and thumbnails. Preserve historical projects and use fresh review copies for any edits.

## 7. Assemble and verify the app data

Assembly is a fresh derived output, separate from a processing run. The current CLI is:

```bash
"$PY" scripts/uniform/assemble_uniform.py \
  --run "$RUN" \
  --donor-source data/development/uniform-20donor-20261029-v18 \
  --donor-metadata scripts/uniform/donor_metadata.csv \
  --output "$ASSEMBLY" --report-root "$FRESH_ASSEMBLY_REPORT" \
  --workers 6 --include-collagen-iv
```

`--donor-source` supplies the existing donor/image-link source; it does not select that root's old phenotypes. The explicit nPOD metadata table supplies the authoritative current donor fields. Check the actual source coverage before running. `--include-collagen-iv` retains the instructed structure endpoints.

Assembly rebuilds the app's cell labels, per-object core/peri summaries, marker displays, spatial records, structure endpoints, donor counts, and embedding inputs. It records input provenance in `uniform_release.json`. Core membership uses the cell centroid. The peri region is a 20 µm ring within Tissue with endocrine objects removed. A cell can belong to multiple peri regions. Singles are counted and shown, but are not islet-level analysis rows or separate core/peri objects.

The composition UMAP/Leiden and trajectory input use the named Beta, Alpha, Delta, and CHGA+ hormone-negative shares for objects with at least ten core cells, with no donor correction. The Harmony protein embedding is retained as a separate record. Rebuild these derived representations when component identities change; do not interpret pseudotime as longitudinal disease progression.

Verify the exact donor inventory, UUID uniqueness, all required files and image links, matching counts across cell and summary tables, excluded-donor absence, and consistent endocrine/immune definitions between assembly and UI. The corrected exporter preserves binary `Ki67_pos` in cell types, per-object cells, donor cells and cell-distance tables, with per-image measurement/level provenance in `uniform_release.json`. It verifies every state against the selected upstream measurement. The release parity check verifies those states survive runtime conversion. Ki67 is not an identity or a new embedding input. Both frozen v19 and v20 omit this field. The local immune-gradient correction includes `Hematopoietic CD45+`, which the last dev app's older list omitted. It needs a rebuilt release before it can appear on a future served site.

The release builder and parity verifier reject unresolved endocrine definitions or disagreement with the app endocrine aggregate. Removing unnamed cells only from a denominator is not a correction. Preserve the failing report and return to the input discrepancy.

## 8. Build a release and review the app on dev

Use one clean committed release checkout for the OCI image, data builder, and parity checks. The existing linked release checkout is `/data/panc_CODEX/Islet-Explorer-release-20261022`. Its commit must match the image revision, data manifest, and release metadata. A source change after building requires rebuilding the affected release artifacts; editing metadata to make commits agree is invalid. Current security builds use separate clean worktrees described in the security record; do not move a checkout another task is using.

The project-specific helpers live in `tmp/ufrc-release-20261022/`. They are local operational files, not a portable tracked interface. Read their paths and defaults before reuse; several default to v10. From that helper directory, after selecting a clean release commit and a new unique release name, the established sequence is:

```bash
RELEASE_NAME="$RELEASE" ./build-image.sh
SOURCE_ROOT="$ASSEMBLY" ./build-data.sh
SOURCE_ROOT="$ASSEMBLY" ./verify-data.sh
RELEASE_NAME="$RELEASE" ./bind-release.sh
```

The tracked commands and complete instructions are in [deploy/README.md](../deploy/README.md). Data packaging uses `scripts/build_release_data.R --schema 5`, then `scripts/verify_release_data.R --mode full --schema 5` and `scripts/verify_release_parity.R`. The app interpreter for those R operations is selected by `RETICULATE_PYTHON=/home/smith6jt/.local/share/islet-explorer-py/bin/python`. Schema 1 belongs to the frozen 15-donor release, not this 20-donor workflow.

Use copied files for a reviewer bundle so later source edits cannot alter it through shared inodes. `--materialize hardlink` is only for disposable local tests. Keep image data and credentials outside the OCI image. Preserve the SBOM, vulnerability scan, dependency evidence, immutable image ID, data manifest, parity/inventory reports, and full logs. Current security validation also requires the image's Brotli library evidence from `deploy/verify-brotli.py`, invoked by `deploy/scan-image.sh`: patched Ubuntu/Node, Pillow linked to the system decoder, and DuckDB's actual decoder regression with its installed-library hash. New release manifests require the version-2 report; missing reports and older Pillow-only reports are rejected. A clean package scan alone missed the embedded decoders. Never rewrite the directory of an active release.

The owner's latest November 2 direction is to deploy the current validated build to the live PubApp within the DuckDB remediation window through **2026-11-10 UTC**, then develop and test our DuckDB fix separately for dev. The production candidate retains the original hash-pinned CRAN DuckDB 1.5.2. `scan-image.sh` records the real failed binary report plus a distinct dated nonblocking release decision; manifest and full/quick release checks enforce that deadline. Other package findings and scientific/data validation are not waived. The preserved backport at `3a12daa` is available for the isolated dev worktree. Follow the [dated decision](duckdb_remediation_window_20261102.md), [build commands](../deploy/README.md) and [security evidence](container_security_20261102.md). The real-Parquet test for the dev fix writes only temporary BROTLI-compressed output and compares schema, values, missingness, duplicates and aggregates. An older read-only payload can support a diagnostic smoke, but cannot supply the matching release/data acceptance.

A verified release can first be smoke-tested with the production runtime using `deploy/run.sh`. It takes environment variables `PROJECT_ROOT`, `RELEASE`, `PORT`, `CONTAINER_NAME`, and `MEMORY`; it binds loopback and verifies the release before launching. The documented local smoke port is 3843. Check the container's health and `trajectory_warm_ok` in its log, then run the expanded browser suite. Do not substitute a development-mode launch for this test.

The following is the established helper invocation for **future dev activation**, after explicit restart authorization, a verified patched image, scientific acceptance appropriate to that review, and deliberate handling of the service mask. It is not a current restart instruction:

```bash
RELEASE_NAME="$RELEASE" INSTANCE=islet-explorer-dev PORT=3839 MEMORY=6g \
  PUBLISH_ADDRESS=10.13.41.21 ./deploy-vm.sh
```

When the service is permitted to run, this stages through HPG to `pubmaigan1`, validates the bundle/image, and activates dev. It cannot restart the currently masked service. Keep its release name distinct from every existing release. The VM's project root is `/pubapps/maigan/islet-explorer`; dev uses `dev/html`, `dev/auth`, and `runtime/islet-explorer-dev/`. The image aliases are served by UFRC's proxy from release-linked files, while Shiny runs in rootless Podman.

After activation, run the actual proxy browser check with a restricted credential file:

```bash
node deploy/browser/check-app.mjs \
  --url https://devislet-explorer.rc.ufl.edu/islet-explorer/ \
  --auth-file "$AUTH_FILE" --expanded true --output "$FRESH_BROWSER_REPORT"
```

Inspect the same candidate yourself. Plot, Trajectory, Spatial, and Viewer are the reviewer tabs; Statistics, Inflammation, and AI are not part of the reviewer profile. Check cohort and phenotype selectors, linked object selection/clear-selection, cell drill-down, endocrine and immune totals, spatial gradients, exported CSVs, and Viewer image/geometry alignment. The Viewer should stream TIFF ranges with its XML/offset sidecars, rather than download an entire image. Compare representative displayed objects directly with the candidate's records and new QuPath project.

Browser checks establish runtime behavior; source/release parity establishes packaging equivalence. Neither establishes correct scientific cell identities. Record scientific review separately.

## 9. Production activation, rollback, and operations

Production activation follows acceptance of the exact dev-reviewed image/data pair and explicit authorization to lift the current stop. **Nothing in this guide authorizes unmasking or activating either site, the unaccepted v20, or the whole-v19 reference.** The owner controls the production decision. Record the selected patched release, source commit, image ID, manifest hash, scientific/security review evidence, and intended rollback before activation.

For the existing helper workflow, the documented promotion command, run from `tmp/ufrc-release-20261022`, is:

```bash
RELEASE_NAME="$RELEASE" PUBLISH_ADDRESS=10.13.41.21 ./deploy-vm.sh --skip-stage
```

Its defaults select production `islet-explorer`, port 3838, 9 GiB, and `prod/html`. `--skip-stage` requires the exact accepted release already staged and verified. The original handoff describes the owner's `!` invocation for this operation; it is not evidence that the present candidate has approval.

On the VM, the tracked equivalent is `tooling/deploy/activate-release.sh` with `--project-root /pubapps/maigan/islet-explorer`, `--release`, `--publish-address 10.13.41.21`, `--docroot /pubapps/maigan/islet-explorer/prod/html`, and `--local-data /local/maigan/releases`. Activation rechecks content hashes, starts the immutable image, verifies health, and only then records active/previous state. If startup fails, it attempts to restore and verify the previous unit. Check the exit status and recovery report rather than assuming rollback succeeded.

The tracked rollback commands, executed as the service account on `pubmaigan1`, are:

```bash
tooling/deploy/rollback-release.sh --project-root /pubapps/maigan/islet-explorer
tooling/deploy/rollback-release.sh --project-root /pubapps/maigan/islet-explorer \
  --instance islet-explorer-dev
```

Choose the intended instance. Rollback restores its recorded previous release, including the image aliases and local-data mounts; verify the destination's security and scientific acceptance first. Activation and rollback cannot restart the currently masked services. The obsolete image identified in the security alert must not be restarted or selected for rollback. Preserve the masks until an explicitly authorized patched release is ready; do not unmask a unit while it still points to a vulnerable image. The current masks also prevent automatic restart after reboot.

Operational evidence includes `systemctl --user status`, `journalctl --user -u <instance>.service`, the instance's runtime logs, deployment state under `runtime/deployments`, and proxy logs under `/pubapps/logs/maigan/`. Keep at least active and previous image/data pairs until the new release has passed the reboot test. Deleting old workstation bundles or VM releases is a separate owner decision.

The remaining UFRC work is recorded in [readiness](ufrc_readiness.md): the clean non-UF acceptance run, reviewer credential handoff, and reboot test remain pending; trial allocation RC83797 ends **2026-11-30**. UFRC's temporary `ufrc` login remains in `dev/auth` pending completion of its work. Do not publish credentials in Git, logs, command arguments, or tickets.

Use `deploy/verify-reviewer-access.py` with the selected URL, a real TIFF `--image-path`, its independently calculated `--image-prefix-sha256`, the XML and offsets `--sidecar-path` entries, `--auth-file`, and a fresh `--output`. It verifies anonymous rejection and authenticated page, range, sidecar, and WebSocket access. Run from the required outside network as well as the VM/browser checks; a campus pass is not the outside-network acceptance run.

## 10. HPG processing, reproducible records, and handoff

The authoritative HPG workspace is `/blue/maigan/smith6jt/islet-explorer-hpg/root`. Its site wrapper presents the workstation path `/data/panc_CODEX/Islet-Explorer` inside the container and the sealed QuPath inputs read-only. Use [the HPG processing record](hpg_processing_20261031.md) for the current seal and environment. Do not repair a seal failure by weakening the guards or editing a protected scientific file.

Workstation code is the editing authority. Commit a coherent change, transfer the explicit changed-file list and required Git metadata, verify remote hashes, and record the receipt before using it on HPG. Never perform a whole-tree `rsync`, `git checkout`, `reset`, or `stash` into the sealed HPG scientific tree. Administrative operations that change pinned file identities require the documented reseal procedure; they are not routine pipeline reruns.

Submit scientific work through SLURM from the HPG repository root. For example, with `HPG_RUN` and `HPG_REPORT` set to the selected run and a new report path inside that site:

```bash
/opt/slurm/bin/sbatch --job-name=endo-audit --cpus-per-task=4 --mem=16G --time=01:00:00 \
  scripts/hpg/jobs/site_job.sbatch \
  python scripts/uniform/typing/endocrine_object_audit.py \
  --run "$HPG_RUN" --out "$HPG_REPORT"
```

`site_job.sbatch` invokes `scripts/hpg/run_in_site.sh`, records host/job/timing/resource information, and uses the `maigan` account/QOS and `hpg-turin` partition by default. The owner has allowed `maigan-b` for this work; choose resources/QOS for the actual job. Do not run scientific extraction on an HPG login node to avoid a pending queue. Inspect `/opt/slurm/bin/squeue` and `/opt/slurm/bin/sacct`, and the log `/blue/maigan/smith6jt/islet-explorer-hpg/logs/<job-name>-<job-id>.out`. Record both SLURM state and application exit code.

Return only new results, with `SHA256SUMS` and source provenance, to fresh workstation paths. Verify the returned bytes before using them. Read-only analysis can run locally on a separately transferred, hash-verified diagnostic subset when it has all required inputs; that subset must not be mistaken for an updated canonical processing run.

For every meaningful run, keep the command/arguments, timestamps, commit and relevant file hashes, input manifests, environment/model, SLURM job ID, stdout/stderr, exit status, reports, and output hashes. Keep failed attempts with their reason; use a new output directory for a correction. Large Parquets, rendered panels, raw images, runtime logs, and credentials do not belong in Git. Track the scripts, decision record, and compact evidence receipt pointing to the retained outputs.

Related Senior, monorepo, HPG, and multiplex5000 workspaces/agents receive the [committed correction notice](typing_change_notice_20261102.md). The discovered workspaces are:

| Recipient context | Location at discovery |
|---|---|
| Workstation Senior monorepo | `/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior`, branch `fix/portability-blockers-20261102` at `610a6a1`, with active local changes. Shared code is in its `core/multiplex5000/`. |
| HPG Senior main | `/blue/maigan/smith6jt/Islet-Explorer-Senior-main`, `main` at `40fde07`, with local changes. |
| HPG current typing | `/blue/maigan/smith6jt/IES-typing`, `typing/per-image-method` at `a9549cc`, clean at discovery. |

The immutable handoff package locations are `/home/smith6jt/IO60panc2nd/CODEX-handoff-20261102-lymphatic-isolation` and HPG `/blue/maigan/smith6jt/CODEX-handoff-20261102-lymphatic-isolation`. Consult `docs/evidence/agent_delivery_20261102.json` for the checksum receipts and the Senior/MultiPlex5000 agent's explicit acknowledgement: it read the notice/workflow and independently verified all 17 payload files. The package is received/read; the port remains unapplied. This delivery must not edit active recipient checkouts. The Senior agent owns adapting and testing a port into `core/multiplex5000`; CODEX retains its candidate and app ownership. CODEX's own HPG source synchronization has its separate explicit changed-file receipt.

A supplemental package, `CODEX-handoff-20261102-component-build`, carries the
completed component/Ki67/build fixes, release preflight correction and outcome
evidence to the same workstation/HPG sibling locations. The
[change notice](typing_change_notice_20261102.md) lists its commit groups and
receipt. Keep it distinct from the initial 17-file isolation package.

The handoff carries the actual commit range, changed-file list, evidence receipt, guide, and blockers: separate calibration reference context from mutable structure context; bind all marker/canary calibration and baseline identity; protect named islet endocrine cells; run the all-component invariant audit; preserve the selected v19 endocrine component; do not infer that this diagnostic is a deployable root. File delivery is not an agent read acknowledgement or a completed port. Record those later outcomes separately and preserve workspace-specific changes. Do not transplant CODEX thresholds or data into another modality merely because its pipeline shares code ancestry.

## 11. What remains before the app can return online

1. Resolve the mixed candidate's 2,472 component boundary conflicts: preserve the exact v19 endocrine membership or explicitly allow additional calls under its hormone thresholds. Keep the unresolved decisions visible until the rule is defined.
2. Resolve the inherited endocrine discrepancy without inventing cutoffs, assigning unsupported labels, or excluding cells to hide it. The completed candidate has 8,629 unresolved cell labels and 953 populated objects with no positive eligible member under corrected-cell thresholds: 915 Singles and 38 Clusters. All 1,468 Islets have positive members. Two Singles are empty and one has two eligible members. Original pixel-classifier positivity is a different test; see the classifier trace before changing membership or thresholds.
3. Rebuild into fresh paths after the evidence-supported correction. The current app, QuPath, Ki67 and invariant checks are complete, but changed identities require regenerated summaries, scales, spatial outputs and embeddings. Pass actual component provenance, object consistency and release-label/aggregate checks. A positive islet does not resolve an unnamed member's typing.
4. Reconcile the final QuPath and app outputs with the required image/geometry checks. Retain reviewed exceptions and failed evidence; do not create a new permission requirement solely from an old diagnostic metadata flag.
5. Package and verify the matching schema-5 data release and image. The current image's installed-runtime smoke and dated security check pass. Bind the final payload to that image and validate the actual runtime queries, downloads and image routes. DuckDB's known raw failure is nonblocking only within the recorded deadline; the other data/security checks remain mandatory.
6. The owner has already directed the validated current build to live production. Once the actual checks pass, activate that exact release with the production mask handled deliberately, verify service/proxy/image aliases, and finish the outstanding outside-network and reboot acceptance work. Coordinate the separate patched-DuckDB dev build and its service changes with the security owner; do not require a new production permission prompt.

Hypothesis work remains a separate analysis stream. Retain its approved definitions and historical outputs under `output/app_hypotheses_20261026/`; the typing correction is not permission to rerun or reinterpret it against a new data root automatically. CODEX INS intensity is not an approved quantitative between-donor hypothesis endpoint under the owner's acquisition clarification. Report measurements and technical evidence without adding biological interpretations beyond the approved plan.

The practical stopping rule is simple: a reproducible inconsistency stays visible and blocks acceptance until its evidence-supported correction is recorded. Neither changing a label count nor successfully launching Shiny makes the scientific inputs accurate.
