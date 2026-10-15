# Typed-cell reprocessing and app revision

**2026-09-23 status:** the local default app has been restored at `6b7884b`.
The owner requires exact preservation of original islet identities and values.
The all-21 regeneration contract recorded below is therefore pending revision
and evaluation; its results are not a compatible donor addition. The original
15-donor subset has 8,483 regenerated IDs versus 5,214 frozen IDs, with no exact
overlap. Candidate launch now requires the read-only original-parity gate
documented in `scripts/research/README.md`. This correction does not redeploy
or attest to the separately running public service.

The [phase-2 contract](donor_expansion_phase2.md) now governs additions.

Implementation started on 2026-09-20. The historical candidate was
`data/development/all-donor-20260922-trajectory-fit-support-v1` and supplied the shared
workstation app at that checkpoint. This remains an execution record, not UFRC or scientific
acceptance.

## Historical all-21 experiment contract (superseded)

The candidate contains seven ND, seven autoantibody-positive and seven T1D
donors. HPAP-078 provenance remains explicit. Original inputs and the frozen
release remain preserved. The working candidate is
`data/development/all-donor-20260922-trajectory-fit-support-v1`; scientific reports belong in
`data/reports/reprocessing-20260920`, outside the app runtime.

1. Regenerate all 21 donors' boundaries with the frozen INS/GCG/SST pixel
   classifier and the physical size-dependent rules from
   `~/IO60panc2nd/scripts/1_panc_objects.groovy`. CHGA is not a boundary feature.
   Record and exclude image artifacts from assessable tissue. Keep rejected
   candidate geometry and raw objects for inspection. Agent image audits may
   resolve clear cases; unresolved consequential choices remain explicit.
2. Use native QuPath detection UUIDs to join raw means, corrected means,
   exact centroids, new core/peri membership, cell labels and image displays.
   Recompute REDSEA correction from declared images and cell polygons.
3. Fit four donor-balanced cell-discovery ensembles: log1p with landmark
   median/IQR scaling and arcsinh with pooled training-derived cofactor and
   standard scaling, each with Euclidean neighbors at k=15/30 and Leiden
   resolution 1. Use seed 20260919 and up to 20,000 landmarks per donor.
   All 31 varying markers enter discovery and lineage refinement. The required
   nine markers must also have measured contribution/ablation evidence.
   Annotate clusters using profiles and disease-blinded image references.
   Require at least three supported, reviewed identity votes out of four.
   Separate exclusive base identities from overlapping marker-defined states.
4. Define E as the unique validated core alpha/beta/delta and CHGA-supported
   hormone-unresolved endocrine union. Unknown and zero E remain unavailable.
   Fit the islet trajectory on separate core and tissue-clipped peri-20-µm
   `log1p(base-type count/E)` features. Overlapping states are sensitivity
   features. Cell typing is an upstream dependency of this axis.
5. Keep the outcome-independent diffusion specification: equal donor
   representation, up to 200 size-stratified landmarks per donor, robust
   scaling, cosine k=15, ten components, seed 20260916, hash-oriented
   two-sweep root, fixed Nyström projection and largest-component ordering.
   Reconstruct for 21 donor omissions, 200 donor bootstraps and the prescribed
   parameter sensitivities. Report instability and failed fits; do not tune
   for disease ordering. Filtering the display cannot change coordinates.
6. Retain ratio-only trajectory scatter and all 31 heatmap markers. The
   default heatmap shows mean log1p corrected intensity within an explicit
   typed population, initially E. Marker-positive count/E remains an
   alternative that stays visibly unavailable until positivity gates exist.
   Use equal-islet means within donor and then equal-donor bin means.
7. Rebuild the ten registered donor-aware hypothesis tests and their
   result-dependent exploratory follow-ups. H4/H6/H7/H8 use all 21 donors
   with categorical disease adjustment. Disease contrasts use T1D versus ND
   primarily. Structural density means actual positive area/assessable area.
   Formal outcome–pseudotime tests require outcome-excluded, donor-held-out
   typing and trajectory reconstruction. Endocrine-module tests are
   non-estimable under this fixed endocrine denominator.
8. Build schema-4 runtime summaries directly from canonical Parquet. Preserve
   frozen-schema readers and the existing caching, spatial indexes and image
   range serving. Research results remain outside the app. Verify an isolated
   candidate before coordinated shared-service cutover, and retain rollback
   of app code, data, profile and service links together.

## Findings incorporated during execution

The older 15-donor H5AD-derived tables used numeric row-index IDs, whereas
native QuPath measurements and REDSEA outputs use UUIDs. The new common
universe contains **3,619,534 native cells**, 26,812 more than the prior
3,592,722-row candidate. The previous tables remain in
`provenance/archival_raw`. All 21 donors' exact native centroids were exported
read-only and reconciled; rounded legacy centroids are not membership inputs.
The native tables contain 248 wholly unmeasured objects and no partially
measured objects; raw and correction-specific exclusions are recorded separately.

An independent 651-comparison raw-pixel/native-measurement audit identified
reversed HLADR/PDPN labels in ten original OME images: 112, 6479, 6521, 6533,
6534, 6538, 6548, 6549, 6550 and 6551. Native acquisition metadata independently
confirms the permutation. Candidate image metadata is corrected on independent
copies; frozen/source pixels and images remain unchanged. Scalar REDSEA is
channel-independent, so relabeling its two output columns is mathematically
equivalent to rerunning correction with the corrected labels. The first
discovery attempt was stopped and preserved under
`provenance/withdrawn_typing_hladr_pdpn_mismatch`.

Five boundary refinements exposed a mixed-dimension geometry edge case:
polygon collisions included detached touching lines, and JTS refused their
subtraction. The fix subtracts only polygonal collision area; it does not
change size thresholds, erosion/dilation or fragment rules. Eleven native
regression checks pass, including the new mixed-geometry case. All donors'
adaptive boundaries are rebuilt with the same fixed implementation.

Visual review identified a multichannel registration band at the bottom of
donor 6563 that generated many small candidate regions. Its technical exclusion
is recorded from image evidence; usable tissue from this donor is retained.
Other uncertain tiny or diffuse candidates receive individual flags, not
blanket donor rejection.

The final boundary audit covers all 21 donors with 126 systematically selected
image regions and an all-donor edge screen. Sixteen candidate-specific flags
exclude uncertain or technically unassessable regions from primary analyses;
no donor is rejected. There are 11,944 regenerated cores. The raw native core
membership sets agree exactly with independent Python containment for every
donor. Analysis compartments additionally clip confirmed technical regions.
Overlap is recorded for each actual core/peri compartment and band width,
rather than applying 30-µm sensitivity exclusions to otherwise valid cores.

A separate native-cell review examined 36 examples from donors 6232, HPAP-078
and 6563 across four cell-area strata. All 18 sampled objects below 5 µm² were
subcellular dots or slivers, including some with high hormone signal. The
primary cell-quality rule therefore requires area at least 5 µm²; thresholds
of 1 and 10 µm² are sensitivity analyses. This is segmentation quality control,
not islet-size normalization. Confirmed technical image regions are also
excluded before both discovery and projection. No uncertain biological islet
flag alone excludes cells from whole-image discovery.

The complete native universe is retained. There are 20,550 objects below the
primary cell-area threshold and 402 objects in confirmed technical regions,
with 52 in both groups: 20,900 quality exclusions in total. Corrected means
are fully finite for 3,619,095 objects; 439 have unavailable corrected values,
including objects with no rasterized pixels. After both quality and measurement
checks, 3,598,630 cells are eligible for typing. These are distinct from the
248 wholly unmeasured raw objects. The final four-ensemble discovery
fit applies the combined eligibility policy before landmark sampling.

Marker-positivity mixture splits are diagnostic proposals only. They require
separate image evidence before becoming available count measurements. The
accepted historical RESTORE reference-pair work is being reconciled with fresh
21-donor measurements; neither historical acceptance nor mixture posterior
probability alone certifies a new gate. Continuous typed-population expression
and nullable marker-positive counts have independent availability.

The fresh RESTORE audit also exposed a historical numerical provenance defect:
many archived level-2 values labeled as negative-control maxima equal the
controls' 99.5th percentile instead. Exact historical-sample reconstruction
separates this from sampling or REDSEA changes. The rejected adaptive divisor
is not reinstated. Original outputs and their biological reference-review
history remain preserved; see `docs/restore_positivity_audit_20260920.md`.

All 21 native QuPath project metadata records specify 0.3774 µm per source
pixel. That calibration controls cell area, distances, compartment bands and
the editable QuPath delivery. The slightly different OME acquisition rounding
is recorded separately, with native-project and image-metadata hashes, in
`preprocessing/calibration.json`.

## Execution and evidence

The entry points are `scripts/research/reprocess_candidate.py`,
`correct_cells.py`, `reconcile_channel_identity.py`, `rebuild_spatial.py`,
`cell_typing.py`, `cell_typing_validation.py`, `composition.py`,
`trajectory.py`, `analysis_v2.py` and `build_runtime_schema4.py`.
Inputs, parameters, checkpoints, cell universes, software/code hashes and
failure reasons are recorded in the candidate. Images/models/reports/logs and
local profiles remain ignored by Git.

Use the installed analysis Python at
`/home/smith6jt/miniconda3/envs/scvi-env/bin/python` and QuPath 0.7 at
`/home/smith6jt/QuPath/bin/QuPath`. CPU work is budgeted at 56 cores and memory
at 384 GiB across simultaneous stages, with donor-level checkpoints and bounded
query batches. Finished stages are not rerun merely to keep processors busy.

Acceptance requires the new data's own scientific and runtime evidence. Local
preview availability is separate from UFRC acceptance. Commit/push remains
conditional on completing the authorized plan and auditing the final staged
diff and `.gitignore`; prior passing evidence is not a substitute.
### 2026-09-22 classifier eligibility and trajectory status

The reviewed endocrine refinement was applied to all 21 donors and merged with
the immutable broad hierarchy. Composition outputs are in the candidate's
`trajectory/` and `analysis/` directories. Classifier eligibility now has only
the owner-specified rules: no image-edge contact and at least one canonical
cell centroid in the endocrine pixel-classifier core. It yields 11,136 eligible
islets, with 137–1,408 per donor. Downstream typing and phenotype gates cannot
overrule eligibility. The fixed denominator is the census of all canonical
cells inside the classifier core, counted once.

Of the eligible islets, 7,493 have complete typed-composition features. The
diffusion graph is trained on 1,514 valid islets with at least ten measured core
cells, using 20 landmarks per donor. Its landmark graph has one component, and
Nyström projection supplies 7,490 supported primary coordinates. Median
leave-one-donor-out absolute Spearman is 0.754 and total classifier-eligible
coverage is 67.3%; all 21 donor omissions and 200 bootstrap reconstructions
completed. The ordering remains unstable by the prespecified checks. Spatial
UMAP and Leiden use the 7,493 feature-supported islets. Earlier gates and fit
counts are superseded by these artifacts.

The 7,490 raw primary coordinates contain 1,729 unique values; 2,234 observations
share the largest exact tie. The 2026-09-22 fit-support revision assigns 7,490
fixed midpoint display ranks across 0–1. Exact ties are
spread by SHA-256 of islet ID only. This order within a tie is descriptive and
has no biological meaning; supporting inferential models use raw pseudotime.
The public browser regression passed the evenly spaced rank grid, nonnegative
ratio axis, nonzero canonical beta/core values, Spatial views, all 31 heatmap
markers, donor filters, downloads and image click-through.
