# Expanded-cohort image processing

**18 September delivery:** [HPAP approval and five-donor rollout](../../docs/endocrine_boundary_rollout_20260918.md)
documents the current six-donor review project. `stage_adaptive_rollout.py` and
`promote_reviewed_islets.groovy` bind the owner's approval to exact saved HPAP
objects. `export_exact_reviewed_geometry.groovy` preserves native coordinate
precision for `rebuild_reviewed_islet_regions.py`. New-donor predictions use
`--minimum-votes 160`, followed by native tracing and the adaptive boundary
stage. Approval of boundaries and validation of phenotype gates stay separate.

**Current boundary work:** use the endocrine foreground followed by the owner's
[IO60 size-adaptive rules](../../docs/size_adaptive_islet_boundaries_20260917.md).
`build_size_adaptive_islets.groovy` replaces uniform boundary cleanup with
size-dependent erosion/dilation, split handling, enclosed-hole filling and
50 µm² fragment refinement. It consumes verified raw foreground from
`vectorize_endocrine_pixel_regions.groovy` and writes a fresh output folder.
New predictions must complete this stage before review import. The SAM R5
candidate remains withdrawn after repeated HPAP-078 errors. The processing
workflow below records the preserved historical run.

The separate endocrine-seeded boundary review added after the small HPAP-078
annotation finding is documented in
[`docs/islet_regeneration_20260917.md`](../../docs/islet_regeneration_20260917.md).
It preserves the processing objects described below and does not itself accept
new islets or rebuild downstream measurements.

This workflow stages the existing QuPath project without writing to its images, annotations, classifiers, or scripts. Outputs remain separate from the frozen application release.

The candidate lives at `data/development/expanded-20260911/`. Its `images/` directory holds viewer TIFFs and sidecars; `qupath/` holds the staged project and scripts; `exports/CASE/` holds measurements; `analysis/raw_cells/` holds unphenotyped cell data; `metadata/`, `reports/`, and `logs/` retain identity, validation, and processing evidence. These large or restricted artifacts are Git-ignored. The processing scripts and cohort selection are versioned separately.

The original 6539 stack in `Images/New` duplicates 6090. The correct 6090 was
verified against its raw planes and the distinct 6539 was rebuilt from its raw
channel folder. Corrected paths and measured preparation status are in
[`docs/cohort_expansion_status.md`](../../docs/cohort_expansion_status.md).
The initial conversion directory retains the rejected 6539 image and audit;
never copy it wholesale into a candidate release. Verify the corrected six-image
selection with:

```bash
/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/cohort_expansion/inventory_reconciled.py
```

The original four-case project has already been prepared and processed. These
commands describe its initial staging; do not recreate a project that already
contains completed results. Run one case at a time. GPU execution requires
access to the NVIDIA device and the installed CUDA runtime:

```bash
python3 scripts/cohort_expansion/prepare_qupath.py
python3 scripts/cohort_expansion/run_qupath.py --case 6232 --device gpu
python3 scripts/cohort_expansion/run_qupath.py --case 6558 --device gpu
ISLET_KRT8_PANEL_VERIFIED=1 python3 scripts/cohort_expansion/run_qupath.py --case 6520 --device gpu
ISLET_KRT8_PANEL_VERIFIED=1 python3 scripts/cohort_expansion/run_qupath.py --case 0078 --device gpu
```

The KRT8 workflow alias is documented in `data/development/expanded-20260911/reports/channel_mapping/krt8_workflow_alias_evidence.json`. It preserves the previous QuPath classifier convention for the same acquisition panel; original TIFF marker names remain unchanged.

The recovered pair has an independently staged project under `recovery/`.
`recover_6539.py` reconstructed all 36 selected planes; `prepare_recovered.py`
bound the project to the verified sources and preserved historical annotations.
Those one-time operations are complete. The runner requires their reconciliation
record before accepting either recovered case. Resume with completed checkpoints:

```bash
python3 scripts/cohort_expansion/run_qupath.py --case 6090 \
  --root data/development/expanded-20260911/recovery --device gpu
python3 scripts/cohort_expansion/run_qupath.py --case 6539 \
  --root data/development/expanded-20260911/recovery --device gpu
```

`cohort.csv` preserves the initial identity holds and fingerprints for the four
original imports. `cohort_recovered.csv` selects only the reconciled pair. The
combined `reports/image_conversion/reconciled_inventory.json` selects six distinct
images and supersedes the initial duplicate audit for source selection.

Stages are normalization, Tissue/Islet pixel classification, SAM rendering and VIT_B box refinement, InstanSeg cell segmentation, immune/structure measurements, and export. The runner accepts `--stages` for restarting a failed stage. A completed stage with a different script checksum requires explicit checkpoint review; downstream results must be invalidated if an upstream analysis input changes.

The runner holds an exclusive project lock for its whole run. QuPath saves the
shared project record for each image, so separate case jobs against one project
must be sequential. A metadata overwrite from earlier overlapping launches was
repaired for 6539 from its intact normalization checkpoint; the source URI,
dimensions, channels and calibration were checked before restoration. Evidence
is in `recovery/channel_metadata_repair_20260916.json`.

QuPath 0.7 runs with headless AWT. SAM uses the installed `samapi` environment, its local VIT_B checkpoint, the original per-image display ranges and colors for SST/GCG/INS, full-image downsample 6, 16-pixel box padding, and single-mask inference. Recovered 6539 uses display ranges calculated from its own pixels, documented in `recovery/display_ranges_6539.json`; the ranges affect SAM rendering, not raw measurements. Its old exclusion annotation belongs to duplicated 6090 pixels and remains in `recovery/historical_annotations/38/` for review. The synchronous adapter follows `ksugar/qupath-extension-sam`'s `SAMDetectionTask` and rendering helper. Masks are converted with the installed `samapi.utils.mask_to_geometry`, then transformed to source image pixels and filtered with the original minimum dimensions/area.

InstanSeg uses fluorescence_nuclei_and_cells **model 0.1.1**, QuPath extension 0.1.7, the original 22 segmentation input channels, 512-pixel tiles, 32-pixel padding, and full Tissue regions. Measurements include all original image channels. DJL 0.36.0 uses the installed PyTorch 2.7.1/CUDA 12.8 runtime with two GPU tile workers. `--device cpu` selects the installed PyTorch 2.5.1 CPU runtime with 32 tile workers. Each worker uses one PyTorch intra-op/inter-op thread to prevent nested thread oversubscription. The runner points directly to cached native libraries, so these stages do not download a runtime. Model inputs retain the model's percentile normalization and resolution handling.

Staged copies of the original Groovy scripts adapt removed QuPath measurement-write methods to `MeasurementList.put`, use an explicit equivalent X/Y comparator, repair the `Islets`/`Islet` clipping class typo, remove the duplicated endocrine block in the v2 counter, and restore the composition writer from the preceding counter script. The model path/version and execution runtime configuration are explicit. Source and staged script checksums are in the staged project's `processing_provenance.json`.

Each `exports/CASE/` directory contains `cells.tsv`, `annotations.tsv`, regional composition/marker/target TSVs, annotation GeoJSON, SAM inputs/results, and provenance. `export_manifest.json` is written last and records the full cell count, source filename, coordinate calibration, 31 required marker means, and table checksums. A partial export without this manifest is not ready for downstream analysis.

Numeric table fields use Java's round-trip double representation, recorded as
`numeric_encoding=java_double_to_string` in the manifest. QuPath's default
adaptive formatting rounded centroids enough to change some reconstructed islet
memberships. On 2026-09-16 the five earlier rounded exports/imports were archived
in `precision_revision_20260916/`, and all six cases were exported at full precision.
Segmentation, object identities and regional measurements were retained.

The final `qc` stage exports a deterministic representative islet crop at full resolution, its original DAPI/INS/GCG rendering, and actual cell/nucleus outlines. Render the review figure after this stage completes:

```bash
/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/cohort_expansion/render_qc.py --cases 6232
```

The resulting `reports/segmentation/CASE_cell_review.png` supports visual review of cell and nuclear boundaries. A single crop does not establish biological acceptance across the image.

Logs and completion records are in `logs/CASE/`. Finishing these measurement stages does not establish new release readiness: donor identity, phenotype-method compatibility, joint analysis, application validation, and release evidence remain separate requirements.

After a case finishes all export stages, import its cells with the installed analysis environment. Request only cases with a completed export:

```bash
/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/cohort_expansion/prepare_analysis.py \
  --cases 6232 6558 6520 0078 --retain-unmeasured
/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/cohort_expansion/prepare_analysis.py \
  --root data/development/expanded-20260911/recovery \
  --cohort scripts/cohort_expansion/cohort_recovered.csv \
  --cases 6090 6539 --retain-unmeasured
```

The importer checks image identity, export hashes and cell counts, donor metadata, full-precision encoding, and the original ordered 31-marker panel. It stores raw mean intensities in `X` as float32, stores `log1p` values in `layers['log']`, and preserves float64 centroids in micrometres using a round-trip decimal parser. Viewer TIFFs preserve their exact source calibration; QuPath measurements use the historical 0.3774 µm/pixel convention. HPAP-078 uses numeric donor ID `78` in analysis and `0078` in viewer and QuPath directory names. No missing marker is replaced with zero.

The default importer rejects missing marker means. `--retain-unmeasured` is an
explicit archival mode: it preserves objects with all 31 means missing as NaN
rows and sets `marker_measurements_complete=false`. It still rejects partial
missingness, negative/infinite values and wholly unmeasured inputs. It does not
choose a training exclusion or imputation policy.

These imported data have no `phenotype` column. The current rules cannot emit all 21 phenotype labels in the frozen dataset, and the historical replacement rules are documented as a reconstruction requiring biological review. Recover the original method or review a replacement before assigning new-cell phenotypes or running the joint expanded analysis. Raw cell matrices must not be concatenated with the baseline's rescaled `X` as if they were on the same scale. The unresolved method is recorded in `reports/phenotyping_method_audit.json`. The two `benchmark_phenotype_transfer.py` evaluations (raw log and rescaled features) failed their donor-disjoint accuracy/F1/recall screen; their reports are under `reports/phenotype_transfer/` and `reports/phenotype_transfer_rescaled/`. No benchmark assigns new-cell labels.

After exports, QC rendering and raw imports complete, verify each project and
join their evidence into the corrected inventory:

```bash
/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/cohort_expansion/audit_preparation.py
/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/cohort_expansion/audit_preparation.py \
  --root data/development/expanded-20260911/recovery \
  --cohort scripts/cohort_expansion/cohort_recovered.csv
/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/cohort_expansion/inventory_reconciled.py --with-preparation
```

The audits verify stage hashes, UUID ordering, all raw marker values and their
log transform, coordinates, cell/islet counts, and regional joins. They also
independently reproduce every core's cell count from exported polygons and
unrounded centroids. Exact boundary ties follow QuPath 0.7 `PolygonROI`'s winding
rule; ties and overlapping memberships are recorded explicitly. Other core
geometry types require a reviewed predicate and are rejected by this audit.
A report of completed offline preparation does not establish biological or
release acceptance.

The subsequent endocrine pixel refinement uses the original 15 donors' source
annotations as weak references. See
[the method and artifact record](../../docs/endocrine_reference_refinement_20260917.md).
`export_endocrine_training_references.groovy` reads source objects without saving
them. `refine_endocrine_pixel_classifier.py` prepares and harmonizes samples,
runs donor-held-out fitting, and predicts from endocrine channels and tissue
masks. `vectorize_endocrine_pixel_regions.groovy` traces candidate boundaries in
native source coordinates. These steps do not promote geometry into the app.
