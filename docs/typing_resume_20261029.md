# Cell typing continuation — October 29, 2026

The owner approved implementing the recovered continuation plan. Source Claude session:
`8ce4ae92-f09a-4a34-b446-fc296c96b58e`, interrupted at 20:15 UTC by its usage limit.
Recovered transcript, workflow journals, scratch outlines and the prior typing review are preserved under
`output/typing_resume_20261029/recovered_session/`.

Baseline run: `data/development/uniform-20donor-20261027`.
Baseline assembled candidate: `data/development/uniform-20donor-20261029-v18`.
New working run: `data/development/uniform-20donor-20261029-v19-work`.
Production v10 remains unchanged.

## Decisions in force

- Use the content rule on all 20 tissue images, preserving existing cell/nucleus and endocrine-object geometry.
- Detect cells only in newly included tissue; preserve existing IDs, raw measurements and endocrine objects.
- Retain edge-touching Capillary, Nerve and Lymphatic structures.
- Resolve equal-depth typing ties only with brighter signal and an interior-coverage lead of at least 0.15;
  measure the interior at least two pixels from the cell boundary.
- Include all tissue Collagen IV, including acinar mesh and vessel collagen, as structural area and distance
  measurements. Collagen IV cannot name a cell or alter lineage/context precedence.
- Deliver a fresh review project and candidate data root; do not promote a production release.

## Recovered evidence

v18 contains 3,197,326 analyzed cells and 325,546 detections excluded by cell QC. The prior 256-cell image review
completed. The latest tissue-content workflow launched 20 reviewers but returned no completed verdicts.
The saved coverage tables lack `covint__` columns; the new tie code was not applied to the data.
The Collagen IV re-review completed, but its renders precede the block-grid correction and the new all-collagen
definition. These outputs are comparison evidence, not validation of the new method.

The prior exact check found 87 cell centroids outside the proposed outlines, two newly outside compared with
the old outline. Coarse raster checks substantially overstated losses. Final validation must check full vector
cell/nucleus/endocrine geometry, not only centroids or rounded contained fractions.

## Execution and acceptance

1. Preserve evidence and create an independent candidate copy.
2. Complete exact footprint preservation and image review on all 20 tissue outlines.
3. Pilot dedicated augmentation, then detect and measure added cells; refresh parenting/context/exports.
4. Rebuild established structures without whole-object edge rejection; calibrate and review all-tissue Collagen IV.
5. Recompute REDSEA, context, levels, canary, coverage, identity and display levels with explicit freshness checks.
6. Run production acceptance checks and fresh nonoverlapping image review, assemble the candidate, and generate
   full and blinded QuPath review artifacts.

Every stage records its inputs and completion. Existing endocrine objects and raw cell measurements are preserved;
derived measurements can change with new neighbors and corrected structures. All detections must reconcile across
geometry, measurements, typing, assignments and assembly. No partial review is reported as completed.

## Implementation checkpoint (17:06 Eastern)

The first exact all-image check preserves all 3,522,872 original cells and every nonempty nucleus/endocrine footprint.
Eighty-five centroids remain outside the outline because they already lay outside their own concave geometry;
four zero-area source nucleus sentinels are retained unchanged. No Tissue island is invented for those sentinels.
Actual image review found two large lumens incorrectly included by smooth DAPI/stromal haze, so these first outlines
are not the final accepted inputs.

The common correction under final validation requires local contrast in addition to the existing absolute gates:
a 51-µm square morphological opening, DAPI residual above 0.6 times its positive-residual Otsu threshold,
and stromal residual above its positive-residual Otsu threshold. The image's brightest 1% of stromal composite
is retained to preserve broad bright walls. Existing footprints and the established small-hole/piece rules remain.
This is an agent calibration implementing the owner's content rule, not a new owner-specified parameter.
The trial removes 0.04–3.87% of the earlier content outline, while preserving every rasterized original footprint.
The two problematic large lumens open and the bright-wall control remains. The small approximately 1,304-µm²
source lumen still fills because its residual hole after smoothing is below the 1,000-µm² rule; this limitation is recorded.
Trial evidence is in `output/typing_resume_20261029/tissue_strong_rescue_trial/`; final outlines will be written
separately to `output/typing_resume_20261029/tissue_final/` and need exact checks and final image reviews.

The Aab_6450 augmentation pilot added 2,242 cells to its original 166,000. Detection, measurement and all exports
completed; final preservation verification is in progress. Two implementation issues were found and fixed during
the pilot: missing QuPath static imports in delegated context scripts, and a quadratic Groovy set-equality check.
Restart uses successful stage checkpoints and never re-detects or remeasures the original cells.

The REDSEA readiness watcher, full native review project and separate blinded review project are implemented.
Blinded review removes prediction metadata and excludes the original 256 reviewed cells; its blank verdicts are
not an accuracy result. The historical v17 comparison was repaired after the assembler removed its old label mapper.
Established structure contours are clipped to the exact reviewed Tissue boundary, retaining edge-touching objects.
All-tissue collagen calibration remains preliminary: the initial peak-width method failed on Aab_6549 and is under
investigation before any final cohort result is accepted.
