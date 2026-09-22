# Endocrine pixel-classifier pilot, 2026-09-17

The owner reports many HPAP-078 segmentation errors, not an isolated bad SAM
outline, and proposed a simple endocrine-marker pixel classifier. The next
regeneration method will start from endocrine pixel classification. R5 remains
withdrawn; the original processing workflow below is historical evidence.

## Existing classifiers

The project's `Islet.json`, used before SAM in the six-donor run, takes all 35
channels at 3.0192 µm/pixel. A second existing classifier, `Jesus.json`, takes
only SST, GCG and INS at 0.7548516 µm/pixel. It uses 18 local features (Gaussian,
Laplacian and weighted standard deviation at two scales) and an OpenCV ANN with
18 inputs and two outputs, without a hidden layer.

The endocrine-only classifier was tested **unchanged**. Its original and staged
copies have SHA-256
`71678cd91794f5b763a79b8a099f440d68036cb9dd42edff034728e82ac4e802`.
Its training images/annotations are unknown, so this is not an independent
accuracy validation or a newly trained classifier.

## Test and interpretation

The [comparison gallery](../data/reviews/endocrine-pixel-pilot-20260917/78-bioformats/comparisons/index.html)
contains 30 HPAP examples, including 13 flagged and 17 unflagged R5 masks across
seed-count strata, the reported failure and two large reference structures.
Examples were selected deterministically before inspecting pixel-classifier
output. Their old outlines are comparison locations, not ground truth.

The native QuPath classifier was run over the **whole raw HPAP-078 image**.
Neither SAM proposals nor cell intensity thresholds constrained the prediction.
CHGA is displayed as context and was not an input to this three-channel model.
QuPath itself generated foreground coordinates and clipped exported regions to
the existing valid-tissue geometry. The full classification raster, including
predictions outside tissue, is also retained.

The reported background-overreach example and other sparse fields look more
faithful to hormone signal with this classifier. The large reference structures
remain visible. However, the unchanged model also detects weak/background
patches, produces small fragments and gives broad peripheral contours in some
fields. These observations support testing/refining this method; they do not
establish accepted islet boundaries.

There are 17,856 raw endocrine foreground components, including 16,298 below
80 µm² and 1,558 at least 80 µm². These are **pixel components, not islet counts**.
No size filtering, gap closing, hole filling, phenotype assignment or SAM
refinement was applied to this classifier comparison. The native raw export
also includes ignored/background objects; use its classification labels rather
than counting every feature as endocrine.

The six historical small negative-control outlines have foreground-overlap
fractions 0, 0.0568, 0, 0, 0.00764 and 0, respectively, in the order recorded in
`old_negative_controls.csv`. Even those examples are not all perfectly blank
under this model. No classifier-wide accuracy claim is made from them.

## Refinement and review

The useful next step is to refine an endocrine-only classifier using reviewed
foreground and background annotations. Samples should cover INS-dominant,
GCG-dominant and SST-dominant regions, mixed/weak staining, sparse foci, dense
islets, exocrine tissue, ducts, stroma and staining artifacts. CHGA can be
evaluated as an additional input rather than silently introduced as a mandatory
gate. The owner subsequently authorized using the original 15-donor project's
annotations. The [reference refinement record](endocrine_reference_refinement_20260917.md)
documents the new training, source SAM provenance and held-out evaluation.
This does not recover the unknown training history of the old `Jesus.json`.

Foreground classification and whole-islet boundary construction are separate
steps. Short gaps and enclosed interiors need controlled handling; widely
separated endocrine foci should not acquire a shared background enclosure.
The separate geometric-envelope experiment has five passing regression tests
but has not been selected as the final boundary method. Reference masks and
reviewed training material are still needed for biological validation and
transfer across donors. The 110 HPAP proximity flags are not a complete count
of segmentation errors.

QuPath supports training from small, diverse annotations and restricting the
input channels; cross-image staining variation needs explicit evaluation.
See the [QuPath pixel-classification documentation](https://qupath.readthedocs.io/en/latest/docs/tutorials/pixel_classification.html).

## Artifacts and reproduction

Pilot root: `data/reviews/endocrine-pixel-pilot-20260917/78-bioformats/`.
It contains the unmodified classifier, native pixel-class raster, raw region
GeoJSON and prediction manifest. `comparisons/` contains plots and CSV tables.
No native QuPath object file or user-edited project was changed. There are zero
accepted replacement islets and no new downstream measurements.

The default ImageJ reader attempt was stopped before delivery because it tried
to cache the complete 35-channel image repeatedly. The successful run uses the
explicit Bio-Formats reader and took 141.2 seconds for prediction and export.
The configured Java heap ceiling was 24 GiB. The stopped attempt and log remain
archived under the sibling `78/` directory.

```bash
env JAVA_TOOL_OPTIONS='-Xmx24g -Djava.awt.headless=true -Djava.util.concurrent.ForkJoinPool.common.parallelism=8' \
  /home/smith6jt/QuPath/bin/QuPath script \
  --args /absolute/candidate/78/prompts.json \
  --args /absolute/project/classifiers/pixel_classifiers/Jesus.json \
  --args /absolute/FRESH_PIXEL_PILOT \
  scripts/cohort_expansion/export_endocrine_pixel_pilot.groovy

/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/cohort_expansion/render_endocrine_pixel_pilot.py \
  --prediction-root /absolute/FRESH_PIXEL_PILOT \
  --candidate-root data/development/islet-regeneration-20260917-r5 \
  --selection data/reviews/islet-sam-failure-20260917/hpap-broader-review/selected_examples.csv \
  --output-root /absolute/FRESH_PIXEL_PILOT/comparisons
```

For a separate QuPath training project, copy `endocrine_classifier.json` to its
`classifiers/pixel_classifiers/` directory and load it via **Classify → Pixel
classification → Load pixel classifier**. Its original output label is `Islet`;
that classifier label does not confer biological acceptance on its foreground.
