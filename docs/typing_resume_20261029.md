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
