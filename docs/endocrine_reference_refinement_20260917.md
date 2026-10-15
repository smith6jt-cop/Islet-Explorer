# Endocrine pixel refinement from the original 15 donors

The active boundary construction now follows the owner's
[IO60 size-adaptive rules](size_adaptive_islet_boundaries_20260917.md).
The uniform cleanup and delivery below are preserved comparison evidence;
the fitted foreground model is unchanged in the new boundary review.

The owner authorized using annotations from the QuPath project that generated
the original 15 donors. The source project is
`/data/QuPath_Projects/QuPath_project_IsletExplorer/project.qpproj`.
All 5,382 `Islet` annotations in those 15 images carry a `SAM Quality`
measurement. The project has tissue and islet regions, but no separate set of
manually verified endocrine/background pixel annotations and no encoded manual
review provenance. This refinement therefore uses **weak reference labels**;
it does not treat every old outline or its interior as pixel-level truth.

Source image data and native object files are read without saving changes.
Exports include original annotation UUIDs, source paths, image fingerprints,
native-file hashes, annotation hashes, channel mappings and software versions.
The frozen release, the owner's editable review projects and the withdrawn R5
candidate are preserved.

## Fixed first refinement

- References: donors 112, 6356, 6450, 6479, 6505, 6516, 6521, 6533, 6534,
  6538, 6548, 6549, 6550, 6551 and 6563. All six added donors are excluded.
- Inputs: INS, GCG and SST only; 2× downsampling, at 0.7548 µm per pixel.
  Native 2× TIFF pyramids are used where available; otherwise images are area
  averaged after conversion to native-endian float32.
- Reference foreground: at least 2 µm inside an old islet region, with at least
  one smoothed hormone above that donor's 99.5th percentile in tissue at least
  10 µm outside all old islet regions. This is a training selection criterion,
  **not a validated phenotype gate or definition of an islet**.
- Reference background: the same distant tissue without that high hormone
  signal. Dark reference interiors, uncertain edges and strong signal outside
  the reference islets remain unlabeled.
- Sampling: 6,000 foreground and 6,000 background pixels per donor. Foreground
  is balanced across dominant hormone and reference connected regions where
  possible; background is distributed across image tiles. Seed: 20260917.
- Normalization: each image's tissue median and interquartile range, independent
  of disease labels, annotation classes and held-out donors. Scaling has a floor
  of one raw intensity unit. Statistics are sampled on a fixed grid and saved.
- Features: each hormone smoothed with Gaussian sigma 1 and 4 prediction pixels
  (six features). OpenCV random forest: 160 trees, maximum depth 12, minimum
  node sample count 40, three active variables, equal donor/class sample counts.
- Evaluation: reconstruct the forest 15 times, excluding one donor each time.
  Report agreement both with clear training-reference pixels and with a separate
  sample of **all** source interiors/background, including ambiguous pixels.
  These metrics do not establish manual boundary accuracy.
- Transfer: apply the fixed full-reference forest to HPAP-078. The previously
  selected 30 comparison locations, six old negative examples and the whole
  field support review; they are not ground truth. The existing `Jesus.json`
  output and withdrawn SAM shapes are shown as comparisons.
- Boundary construction: 0.5 foreground vote fraction, tissue clipping, disk
  closing with a 2 µm radius, filling enclosed holes only up to 80 µm², and
  retaining connected regions of at least 80 µm². Smaller foci remain in a
  separate raster. No convex hull, broad gap filling, SAM inference or cell
  phenotype gate is used. QuPath traces the rasters back to source coordinates.

The original reference shapes have SAM provenance. The statement “no SAM” for
the new prediction means no SAM model or target SAM boundary is used to generate
its pixels; it does not erase the provenance of the source reference regions.
The existing tissue mask remains an input requiring its own review.

## Artifacts and commands

Root: `data/development/endocrine-pixel-refinement-20260917/`.
The authoritative training samples are `samples_v2/` and
`normalized_samples_v2/`. The first draft `samples/` and `normalized_samples/`
are explicitly invalid: OpenCV had resized big-endian uint16 arrays without
converting byte order. This was detected before any model fitting. Source
pixels and native objects were untouched. Both TIFF byte orders now have
regression coverage, and invalid sample folders are rejected by later steps.

The completed model is in `model_parallel/`: four independent processes fitted
the 15 held-out forests, with isolated seeds, followed by the full-reference
forest. The earlier serial `model/` attempt was superseded without changing
the statistical specification. All 15 original native object files were freshly
rehashed after training-reference extraction and remain byte-for-byte unchanged.

Each output directory must be fresh. Commands use explicit roots:

```bash
env JAVA_TOOL_OPTIONS='-Xmx24g -Djava.awt.headless=true' \
  /home/smith6jt/QuPath/bin/QuPath script \
  --args /absolute/source/project.qpproj --args /absolute/fresh/references \
  scripts/cohort_expansion/export_endocrine_training_references.groovy

python scripts/cohort_expansion/refine_endocrine_pixel_classifier.py prepare \
  --reference-root /absolute/references --output-root /absolute/fresh/raw_samples
python scripts/cohort_expansion/refine_endocrine_pixel_classifier.py harmonize \
  --samples-root /absolute/raw_samples --output-root /absolute/fresh/normalized_samples
python scripts/cohort_expansion/refine_endocrine_pixel_classifier.py train \
  --samples-root /absolute/normalized_samples --output-root /absolute/fresh/model
python scripts/cohort_expansion/refine_endocrine_pixel_classifier.py sources \
  --source-root /absolute/reconciled_image_metadata --output-root /absolute/fresh/target_sources
python scripts/cohort_expansion/refine_endocrine_pixel_classifier.py predict \
  --model-root /absolute/model --specification /absolute/target_sources/78/source.json \
  --output-root /absolute/fresh/prediction_78

env JAVA_TOOL_OPTIONS='-Xmx16g -Djava.awt.headless=true' \
  /home/smith6jt/QuPath/bin/QuPath script --args /absolute/prediction_78 \
  scripts/cohort_expansion/vectorize_endocrine_pixel_regions.groovy
```

The `sources` command extracts verified images and valid-tissue geometry from
the existing reconciled acquisition record. It excludes SAM prompts, old islet
geometry, cell gates and clinical data from prediction specifications. Previously
computed complete image hashes remain labeled as prior hashes; current path,
size and modification time are checked rather than misreported as new full-file
hash computations.

`export_refined_qupath_classifier.groovy` packages the shared forest with one
target image's fixed normalization and checks deserialization. Its feature
parity patch must be checked before using that image-specific classifier in
QuPath; copying it to another image would use the wrong normalization.

The new measurements and app trajectories must not use these boundaries until
identity and geometry have been reviewed and compartments/cell membership have
been rebuilt. Pixel foreground components are not accepted islet counts.

## Completed evaluation and diagnostic delivery

The model used 180,000 pixels (6,000 per class per donor). Median donor-held-out
balanced agreement on clear reference pixels was 0.99275. Median agreement with
all sampled old islet interiors was 0.85858, with a range of 0.56342–0.91417;
median background-reference agreement was 0.98050. These are weak-reference
agreement metrics, not independent biological accuracy.

The full HPAP-078 prediction took 227.7 seconds. The initial 80/160-tree cutoff
produced 73,425 raw components and 1,526 regions after bounded cleanup. It fixed
the reported empty SAM enclosure but still included weak/background patches
and broad peripheral contours. It was not promoted to replacement geometry.

`calibrate_endocrine_pixel_operating_point.py` evaluates cutoffs using only the
original donors' held-out predictions. No cutoff met the assessed mean
background-reference positive-fraction target of 0.5%. Unanimity reached
0.5178%, with median clear-foreground-reference recall 85.0% and minimum donor
recall 59.4%. The 160/160 version is therefore retained as a **diagnostic
sensitivity**, not retrospectively accepted as a passing primary cutoff.
Replaying the versioned calibration script reproduced all three sensitivity
rasters byte-for-byte.

The [HPAP-078 diagnostic QuPath project](../data/reviews/qupath-endocrine-pixel-diagnostic-20260917/project.qpproj)
contains 1,056 cyan candidate regions from that sensitivity, plus the original
grey geometry. Native serialization checks freshly hashed all **267,993 cells**
before and after import; cell/nuclear geometry, UUIDs, classifications, parents
and measurements are unchanged. All original annotation geometry is preserved.
The source and previously delivered user-editable projects were not modified.

The owner's reported failure is at `PixelIslet_78_00755` in this copy (historical
R5 `New_78_0242`). All six historical hormone-negative examples have zero overlap
with the conservative candidate layer. This does not establish whole-image
specificity. Weak staining, small-focus identity, grouping and boundary extent
still require review. Only HPAP-078 has been predicted with this model; the other
five added donors have not been regenerated with it.

The [comparison gallery](../data/reviews/endocrine-reference-refinement-20260917/index.html)
contains training examples from every source donor, 30 fixed HPAP locations,
whole-image views, CSV tables, and PNG/PDF plots. The
[operating-point page](../data/reviews/endocrine-reference-refinement-20260917/operating_point/index.html)
shows the source-reference tradeoff and four diagnostic comparisons. Results
remain in folders, outside the app.

The experimental native QuPath classifier loads successfully. On one 256×256
prediction-pixel patch its classifications agreed 99.22% with the offline
reader/features, with some differing pixels. This is neither exact equivalence
nor a whole-image check. That classifier is not installed in the diagnostic
project; the saved candidate objects are the exact offline rasters traced using
QuPath's native coordinate transform.

The current delivery status is **diagnostic review only**, with zero accepted
replacement islets. No app, trajectory, phenotype or regional analysis was
recomputed from these candidates.

At final verification these generated folders occupied 4,641,293,334 bytes
(4.32 GiB), including preserved draft samples, models, plots and the independent
native review copy; linked source TIFFs are excluded. Native QuPath jobs used
a 24 GiB heap ceiling for reference export/vectorization/import (16 GiB for the
classifier packaging check). These diagnostic artifacts do not add runtime
requirements to the app or establish UFRC deployment acceptance.
