# Endocrine-seeded islet boundary review, 2026-09-17

The owner confirmed repeated HPAP-078 errors and proposed endocrine-marker pixel
classification. The [pixel-classifier pilot](endocrine_pixel_classifier_pilot_20260917.md)
records the current replacement direction and a native full-image test.

**Quality status update:** R5 is withdrawn as a usable boundary result after the
owner identified SAM background overreach in HPAP-078 `New_78_0242`. The outline
contains 19 cells but only three operational hormone-positive seeds. Its source
support recall was 0.75; the pipeline flagged this but still exported the mask.
Mean intensity above the surround did not prevent extensive dark interior.
`quality_status.json` in the R5 project records the withdrawal; native geometry
and user edits are preserved. The delivery and verification history below
describes file integrity, not successful segmentation quality.

The [failure audit and comparison plots](../data/reviews/islet-sam-failure-20260917/signal-audit/index.html)
now cover all 1,162 R5 boundaries. The diagnostic smooths raw hormone pixels with
a 0.6-µm Gaussian and applies the existing donor cell-mean thresholds as
**provisional pixel cutoffs**. It flags a mask if less than 80% of its area is
within 3 µm of signal. These are review flags, not validated biological labels.

| Donor | R5 masks | Proximity flags | Flags among previously green masks |
|---|---:|---:|---:|
| HPAP-078 | 437 | 110 | 89 |
| 6090 | 221 | 1 | 1 |
| 6232 | 96 | 1 | 0 |
| 6520 | 136 | 2 | 1 |
| 6539 | 130 | 0 | 0 |
| 6558 | 142 | 1 | 0 |
| Total | 1,162 | 115 | 91 |

For the reported `New_78_0242`, 11.54% of the mask contains signal above those
cutoffs and 28.36% is within 3 µm of it. The per-source mask's mean brightness
and sampled seed coverage did not measure the spatial support of its final
outline. The screenshot geometry matches this annotation; the raw-channel
comparison reproduces its misplaced boundary. Other poor HPAP masks, including
previously green `New_78_0435` and `New_78_0408`, show the same broad-background
failure in the comparison sheets.

The cyan comparison is a separate, model-free signal-footprint pilot: 2-µm
closing, filling holes no larger than 80 µm², with every output pixel constrained
to within 3 µm of hormone signal. It retains only components within 4 µm of an
operational endocrine seed. This follows the sparse signal in the reported
example, but can fragment whole islets and retain internal holes in larger
structures. It has **not** been accepted or imported as replacement islet
geometry. Foreground cutoffs, whole-islet versus sparse-focus definitions and
gap/hole handling still need calibration against reviewed examples.

`quality_status.json` exists in both the candidate and native-project roots.
Finalization, native import and packaging now refuse withdrawn candidates,
including attempts to package one in a fresh project. Nine targeted support
and withdrawal tests passed; the actual QuPath import command was also checked
to exit at the withdrawal guard before loading or saving object data.

This is a separate review candidate addressing the small hormone-negative SAM
annotations found during the HPAP-078 review. It does not replace the app,
processing evidence, frozen release or the owner's existing QuPath review
project. No new annotation is accepted automatically as an islet.

**Historical R5 delivery, subsequently withdrawn:**
[`project.qpproj`](../data/reviews/qupath-regenerated-islets-20260917-r5/project.qpproj)
in QuPath 0.7. The [raw-channel gallery](../data/development/islet-regeneration-20260917-r5/index.html)
contains INS/GCG/SST/CHGA comparisons. Green and yellow outlines are newly
generated candidates; grey outlines are the previous segmentation.

| Donor | New boundary candidates | Green | Yellow: extra review | Separate point foci |
|---|---:|---:|---:|---:|
| HPAP-078 | 437 | 271 | 166 | 860 |
| 6090 | 221 | 136 | 85 | 405 |
| 6232 | 96 | 56 | 40 | 150 |
| 6520 | 136 | 44 | 92 | 327 |
| 6539 | 130 | 78 | 52 | 219 |
| 6558 | 142 | 54 | 88 | 336 |
| Total | 1,162 | 639 | 523 | 2,297 |

All six native read-back checks passed, preserving 943,068 original cells and
their measurements. The final package checks passed: boundaries are valid,
non-overlapping and contained in valid tissue, and every operational endocrine
seed is represented by a boundary or point annotation. All 770 SAM input crops
were reproduced exactly from the source pixels. The six checked hormone-negative
HPAP-078 circles have zero overlap with the delivered new boundaries. These
checks establish file and geometry integrity; all candidates still require
biological review. Seven targeted mask/coordinate and support-point tests passed.

The original six-donor SAM masks matched the saved QuPath geometry. The problem
was upstream candidate identity, rather than an import displacement. The six
HPAP-078 examples `Islet_130`, `206`, `1148`, `124`, `161` and `179` have no
INS/GCG/SST-threshold-positive core cells, and raw-channel review showed little
hormone signal in those small outlines. This does not establish that every
hormone-negative region is biologically non-endocrine.

## Candidate definition

`prepare_endocrine_candidates.py` uses the explicit candidate source manifest,
the six-donor native review manifest, source TIFF channel names and measured
pixel calibration. Corrected 6539 is required through the selected source
project. Original cells and source images are retained.

1. A cell is an operational endocrine seed if INS, GCG and SST measurements
   are finite and at least one exceeds its literal exported donor threshold.
   Its centroid must lie in the existing tissue mask minus exclusions.
   These thresholds are **not validated phenotype gates**.
2. Seed centroids are grouped with 20-µm connected-neighbor clustering
   (`DBSCAN`, minimum samples 1). Groups of three or more seed cells prompt
   SAM. Smaller groups remain visible as points; there is no rule deleting
   single-cell or two-cell endocrine foci.
3. INS/GCG/SST form RGB channels. Each image retains its archived channel
   display ranges, without fitting contrast separately to each crop. The box
   extends 12 µm beyond the seed group, with another 25 µm of crop context.
   Crops use full resolution unless their longest side exceeds 900 pixels.
   At most 32 spatially distributed seed locations are moved to the local
   hormone-signal maximum within 4 µm for a post-inference support check.
   Those points are not supplied to the selected box-only SAM method.
4. The existing local SAM VIT_B checkpoint supplies one box-prompted mask.
   Model confidence and local support are not biological validation.
   Connected mask components without a signal-support point
   are omitted from the candidate while remaining auditable in the prediction
   record. Contour conversion preserves disconnected components and holes.
5. Masks are transformed to full-resolution source pixels and clipped to
   valid tissue. Intersecting predictions are unioned, without dilation,
   convex hulls or synthetic bridges. Each final polygon must contain an
   operational hormone-positive cell. Empty-support pieces remain in the raw
   prediction artifacts but do not become review islet boundaries.
6. Green candidates have at least three seed cells, at least 80% coverage of
   their source masks' sampled support points, no crop-edge contact, and
   higher mean hormone RGB intensity inside each source mask than in its
   10-µm surround. Yellow candidates fail at least one of these checks.
   All unmasked seeds remain
   blue point annotations. Every seed must be represented in the export.

The count and distance rules are review-generation assumptions, not validated
islet identity or a minimum biological islet size. Hormone-low structures can
be missed; SAM can fragment a structure or include non-endocrine tissue.
Previous annotations remain available to evaluate these errors and decide
whether candidates should be merged, split, rejected or redrawn. Candidate
counts and areas must not be used for inference before that review.

The first GPU attempt stopped at an unsupported deterministic CUDA `cumsum`
operation. That attempt and its log are preserved. The successful candidate
uses fixed inputs and seed 20260917, with deterministic algorithms requested
in warning mode. Strict GPU determinism is not claimed. An initial preparation
attempt on the remaining donors also stopped safely at acquisition-prefixed
channel names; the retry uses the existing canonical channel mapping.

The 10-µm grouping pilots were not accepted for delivery. In donor 6520,
549 operationally positive cells within the clearly hormone-positive old
`Islet_2` split into 204 groups at 10 µm; all 549 form one group at 20 µm.
Changing point placement alone did not resolve that fragmentation. The selected
20-µm, box-only pilot covered 97.78% of the old `Islet_2` reference area and
96.51% of `Islet_196`, with raw-channel plots showing whole structures again.
Old outlines are comparison references, not ground-truth masks. Wider grouping
can merge nearby structures, which remains a review concern. These parameter
changes addressed observed image geometry; disease ordering was not evaluated
or used to choose them. Failed pilots `r2`, `r3` and `r4` remain archived.

## Output and interpretation

Candidate artifacts are under
`data/development/islet-regeneration-20260917-r5/CASE/`:

- `prompts.json`, `seed_cells.csv`, `crops/`, `valid_tissue.geojson`;
- `sam_predictions.geojson`, `prediction_audit.json`;
- `review_candidates.geojson`, `candidate_review.csv`, `candidate_summary.json`;
- `previous_islet_comparison.csv`;
- raw-channel PNG/PDF sheets for old positive controls, old negative controls
  and new candidates across their seed-count range.

The independently copied native project is
`data/reviews/qupath-regenerated-islets-20260917-r5/project.qpproj`.
`regenerated_review_manifest.json` now has `status=withdrawn_pending_boundary_correction`.
The original delivery record is preserved in
`data/reviews/islet-sam-failure-20260917/`.
`regeneration-verification/native-summary.json` has
`status=passed` for all six donors. The original `review_manifest.json` records the initial copies
before overlays were added; it is preserved rather than rewritten.

Native verification compares every cell UUID, cell/nuclear geometry,
classification, parent UUID and measurement before and after the save. Every
original annotation geometry and every new candidate coordinate must also
survive serialization. The selected native run hashes every saved cell record
again, comparing with the earlier pilot's cell digest only after verifying
that its original native source is byte-identical. The new candidate geometry
is independently checked; earlier mask acceptance is not reused. Original
islets and expanded regions are reclassified
as `Previous segmentation: Islet` and `Previous segmentation: IsletExpanded`
in the new copy. Cell parents, existing distances and old regional
measurements still describe the previous geometry. New candidate seed and
total-cell counts use direct spatial intersection, not the old hierarchy.

The initial HPAP-078 pilot removed overlap with the six problematic examples,
but that pilot is superseded. The selected delivery records its own counts,
comparison results and native verification; consult
`regenerated_donors.csv` and `regenerated_review_manifest.json` in the new
project. Positive-control plots do not validate every boundary or weak focus.

Fresh image checks compare every actual SAM crop with the corresponding raw
TIFF pixels and check resolved paths, dimensions, channel mappings, byte counts
and modification times against the reconciled acquisition records. Prior
whole-TIFF hashes remain linked provenance; they are not described as freshly
recomputed. A slow whole-file scan was stopped without a validation result,
as recorded in the earlier pilot. Raw marker tables and geometry inputs carry
their own verified hashes.

After human identity and boundary review, accepted geometry needs new
core/peri-islet measurements and a new trajectory reconstruction/validation.
Phenotype gates and the endocrine denominator still need their own validation.
No scientific acceptance, live-app replacement, UFRC acceptance or deployment
is implied by availability of this review project.

The native review project occupies approximately 3.55 GiB and the candidate
artifacts and plots occupy 0.30 GiB; full source images are referenced rather
than duplicated. Recorded SAM inference across the six selected runs took
178.8 seconds on the local RTX 4000 Ada GPU. This excludes preparation, plotting,
native serialization and verification. QuPath verification used a 24-GiB Java
heap limit; this is a configured ceiling, not a measured memory requirement.
Generated images, native objects and review artifacts remain excluded by
`.gitignore`; the scripts, tests and this method record are source files.

## Reproduction

Use new output directories; the scripts reject overwriting completed donor
artifacts or delivered native data. The analysis interpreter is the existing
`scvi-env`; inference and mask tests use the existing `samapi` environment.

```bash
/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/cohort_expansion/prepare_endocrine_candidates.py \
  --input-root data/development/all-donor-20260917 \
  --review-manifest data/reviews/qupath-six-donors-20260917/review_manifest.json \
  --output-root NEW_CANDIDATE_ROOT --cluster-distance-um 20 \
  --case-ids 78 6090 6232 6520 6539 6558

CUBLAS_WORKSPACE_CONFIG=:4096:8 /home/smith6jt/miniconda3/envs/samapi/bin/python \
  scripts/cohort_expansion/predict_endocrine_candidates.py \
  --root NEW_CANDIDATE_ROOT --checkpoint /home/smith6jt/.samapi/vit_b/sam_vit_b_01ec64.pth \
  --case-ids 78 6090 6232 6520 6539 6558 --prompt-mode box-only

/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/cohort_expansion/finalize_endocrine_candidates.py \
  --root NEW_CANDIDATE_ROOT --case-ids 78 6090 6232 6520 6539 6558
/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/cohort_expansion/render_endocrine_candidate_review.py \
  --root NEW_CANDIDATE_ROOT --case-ids 78 6090 6232 6520 6539 6558
/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/cohort_expansion/prepare_review_project.py \
  --output-root NEW_REVIEW_ROOT
env JAVA_TOOL_OPTIONS='-Xmx24g -Djava.awt.headless=true' /home/smith6jt/QuPath/bin/QuPath script \
  --args /absolute/NEW_REVIEW_ROOT --args /absolute/NEW_CANDIDATE_ROOT \
  scripts/cohort_expansion/import_endocrine_review.groovy
/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/cohort_expansion/package_endocrine_review.py \
  --candidate-root NEW_CANDIDATE_ROOT --project-root NEW_REVIEW_ROOT \
  --source-inventory data/development/expanded-20260911/reports/image_conversion/reconciled_inventory.json
/home/smith6jt/miniconda3/envs/samapi/bin/python -m unittest discover \
  -s scripts/tests -p test_endocrine_mask_geometry.py -v
/home/smith6jt/miniconda3/envs/scvi-env/bin/python -m unittest discover \
  -s scripts/tests -p test_endocrine_prompt_points.py -v
```
