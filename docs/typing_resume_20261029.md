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

## Historical implementation checkpoint (17:06 Eastern; superseded below)

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

## Implementation checkpoint (17:55 Eastern)

The first augmentation pilot passed exact preservation but failed image review: detecting inside narrow added-tissue
slivers clipped predictions into tiny fragments. Its 2,242 new detections are not accepted. That independent run is
archived at `data/development/uniform-20donor-20261029-v19-pilot1-work`; historical absolute checkpoint paths are not
resumable after the archive move. The fresh candidate was copied independently from the original baseline.

The corrected detector uses 30 µm of image context around added tissue, then accepts only intact cells whose exact
centroids lie in the unbuffered addition and whose cell and positive-area nucleus lie fully within final Tissue.
Existing cells win overlap conflicts. The second Aab_6450 pilot accepted 25 additions; none touch the context or
image boundary. Measurement and all exports have completed; native and export preservation verification and actual
new-cell image review are still required before the cohort run.

An OpenCV contour corner contact exposed another defect: default polygon repair could convert a real mask hole into
filled Tissue. Repair now uses the known shell/hole hierarchy, with an exhaustive pixel-center roundtrip check before
edge cropping. All 20 outlines at `output/typing_resume_20261029/tissue_final_v2/` pass that check with zero additions
or losses, and preserve all 3,522,872 original cells plus every nonempty nucleus/endocrine footprint exactly.
All 20 whole images and six native fields per image have been inspected, with additional changed-region views.
Nineteen have acceptable-with-notes records; ND_6516 remains flagged. Its large empty lumen survives the first
1,000 µm² cleanup, but footprint padding reduces the residual hole below that cutoff and a second cleanup fills it.
An isolated operation-order correction is being compared; `tissue_final_v2` is not a fully accepted final cohort.
Faint fibrous septa and smooth diffuse stromal signal remain explicit limitations in the review records.

The established Nerve, Capillary and Lymphatic redraw method completed a pilot and actual whole/six-field review.
It preserves visible bright structures and remains close to baseline (pixel Dice 0.995–0.997); weak profiles omitted
by the inherited classifiers remain a limitation. Those pilot masks use superseded Tissue and need regeneration.
All-tissue Collagen IV remains under calibration. Global intensity peaks and nuclear high percentiles failed to
retain dim visible mesh in some images. Local robust nuclear-background references are now being compared, without
changing cell-typing levels or using predicted cell identities as a reference. No cohort collagen method is accepted yet.

The 35 focused tissue, structure-metric and augmentation-review tests pass. Full and blinded review generators,
fresh coverage provenance, union-area metrics and independent Collagen IV distances are implemented; completing code
does not substitute for the remaining cohort execution and image review.

## Implementation checkpoint (18:21 Eastern)

The corrected augmentation pilot is verified and reviewed: 166,000 original cells, their 574 raw measurement
columns and all 578 endocrine objects are preserved; 25 intact additions were accepted. None of the additions
are below 5 µm² or touch the detection-context/image boundary. Four additions have zero-area nuclei and remain
subject to the existing cell-QC rules. Whole and all six native views were actually inspected. The pilot is
archived at `data/development/uniform-20donor-20261029-v19-pilot2-work`; the final candidate is a fresh independent
baseline copy because the pilot used superseded Tissue.

Final vector outlines are at `output/typing_resume_20261029/tissue_final_v3/`. Large bounded holes that survive
first content cleanup are restored after footprint padding and the second small-hole fill, except occupied
footprint pixels. All 20 outlines preserve all 3,522,872 original cells and every nonempty nucleus/endocrine
footprint, with zero added/lost mask centers during contour conversion. There are 86 pre-existing concave
centroid exceptions and four retained zero-area nucleus sentinels; no new outside centroid is introduced.
Final v3 is a subset of v2 on every image, removing 1.248 mm² in total. Faint septal strands and initially
subthreshold residual lumens remain documented limitations.

Final overlay comparison exposed a separate raster union bug: a later polygon hole could erase an earlier
island in the accumulated mask. Vector outlines are unchanged. Tissue review, structural measurement and
REDSEA empty-slide raster helpers now union each polygon independently. The REDSEA correction also reads
all Tissue features; its window, margin, rounding and seeded sampling rules are unchanged. The raw cell-label
raster contract is unchanged. Corrected review views are in `tissue_final_v3_review_union`; earlier views are
reused only with byte-identity evidence, and changed views receive a fresh visual review. The cohort launch
waits for all 20 reviewed records and exact checks. All 79 focused tests pass.

The all-tissue Collagen IV candidate is a common local nuclear-reference rule: log1p median plus one
lower-side robust scale, with a two-native-pixel nuclear erosion and 200 µm blocks. Missing/degenerate block
references borrow within the same image; no donor pooling or phenotype labels are used. A bounded two-scale
comparison lost visible mesh and walls without removing all weak fringe in T1D_6563, so one scale is retained
with explicit reference-contamination/unsupported-block diagnostics, not a specificity claim. Final masks
and all-image visual acceptance are still pending. Freshness checks cover current Tissue, nucleus geometry,
method source, outputs and imported-code stability. Structural region union areas and independent collagen
distances are required and validated before candidate assembly.

The three-worker augmentation, three-worker REDSEA refresh and single-worker established-structure redraw
watchers are prepared. Derived typing, final structure review, fresh independent diagnostic samples, assembly
and full/blinded native review projects remain required before completion.

## Cohort execution checkpoint (18:25 Eastern launch)

All 20 corrected-union Tissue reviews and exact preservation checks passed at 22:25:31 UTC. The authoritative
launch record is `output/typing_resume_20261029/tissue_final_v3_acceptance.json`; one subsequent ND_6232 wording
correction (empty footprint is a nucleus) has an archived original review and an explicit amendment linking both
hashes. Exact Tissue totals are 287.2654 mm² before and 288.8261 mm² after, comprising 6.5644 mm² added and
5.0038 mm² removed. These vector totals are distinct from coarse raster coverage diagnostics.

The three-worker augmentation is running on the fresh v19 work root. REDSEA refresh is gated per image on actual
whole/native augmentation review, matching review/render/input hashes and a root-recorded outlier-report checkpoint
in `augmentation_review_release.json`. The former verification-only idle watcher was stopped before any image
finished. Established structures use the already reviewed Tissue; their final actual image reviews are required
before typing context/metrics. Final ColIV uses only reviewed current augmented nuclei and current structure
outputs. An isolated original-nucleus ColIV precompute is diagnostic evidence and cannot pass final freshness checks.

Task-local drivers and logs reside under `output/typing_resume_20261029/`: `launch_cohort.py`,
`render_augmentation_when_ready.py`, `refresh_reviewed_cells.py`, `refresh_structures_when_ready.py`, and
`refresh_collagen_when_reviewed.py`. `geometry_review_gate.py` checks the reviewed augmentation contract.
`final_reconciliation.py` will check every final inventory and review hash; `blind_diagnostic_protocol/` freezes
the six evidence-only verdicts and prediction-key separation before the fresh difficult-case sample is viewed.
These diagnostics will not be presented as an external accuracy estimate.

## Quality hold checkpoint (19:20 Eastern)

The raw augmentation and established-structure generation continue. Downstream REDSEA, final Collagen IV,
typing, metrics and assembly are held by `output/typing_resume_20261029/segmentation_quality_hold.json`.
The first completed Aab_6505 REDSEA and Collagen IV outputs and their original review records are preserved
as provisional evidence. No final candidate or production release has been published.

A targeted Aab_6450 view found a positive-area added nucleus over a dark lumen. It was also present in the
second pilot but was not the selected review object. Production-style normalized model replay reproduces
the solid raw label before geometry conversion, so changing polygon hole filling would not fix it.
Analogous unsupported predictions occur in other images. A DAPI-only model trial rejects the first ghost
but misses a convincing bright nucleus and varies with tile context; it is unsuitable as an automatic gate.

The correction will use a separate filtered generation, preserving the entire unfiltered augmentation and
every original cell. A common native-image review policy applies to newly added positive nuclei: retain
coherent bright or dim DAPI nuclear profiles; reject polygons supported only by dark lumen/background,
diffuse haze or isolated noise/debris lacking a coherent nuclear profile; leave genuinely uncertain cases
unresolved. Small size or dim brightness alone is not a rejection criterion. Zero signal above the image's
log-DAPI Otsu threshold in the eroded interior is a review trigger only, never an automatic rejection.
All remaining additions also receive image review. Missing-nucleus additions remain in the raw inventory
and are excluded from analysis and spillover by existing cell QC. Decisions must identify every addition
and pin reviewed image evidence and immutable source hashes.

Filtering must preserve retained geometry, raw and signed-distance measurements, parenting and endocrine
objects exactly; re-export dynamic annotation counts; and recompute REDSEA and Collagen IV for affected
images. Established masks may be reused only with an explicit unchanged-input lineage certificate.
Original augmentation checkpoints and detector ID lists will not be rewritten. Holds require separate
hash-linked resolution records before any final downstream consumer can run.
