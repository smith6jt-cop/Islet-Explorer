# Expanded-cohort grouping and comparability audit — 2026-10-14

> **Updated 2026-10-15:** donor 78 (HDL078 / HPAP-078) was removed from the app
> and all analysis by owner directive. The figures below describe the current
> v3 root; the 2026-10-14 audit covered the 21-donor v2 root.

The local app's expanded-v3 source has 20 donors: seven ND, seven Aab+, and
six T1D. The H5AD, composition, target, marker, neighborhood-feature,
distance, and spatial tables all contain records for the full cohort. The
7,028 islet rows in the H5AD have a valid status and non-missing core and
peri-islet phenotype proportions. Donor status is assigned from the explicit
cohort manifest for additions and retained from the frozen baseline for the
original donors; it is not inferred from cell measurements.

| Donor status | Original 15 | Added 5 |
|---|---|---|
| ND | 112, 6356, 6479, 6516, 6548 | 6232, 6539 |
| Aab+ | 6450, 6505, 6521, 6538, 6549 | 6090, 6558 |
| T1D | 6533, 6534, 6550, 6551, 6563 | 6520 |

The T1D group now has one donor fewer than the ND and Aab+ groups. Age, sex
and autoantibody results are present for all 20 donors.

## Processing that is shared

- The original 31-marker QC rule is applied to additions with the frozen
  original-cohort cutoff (total raw MFI at least 5234.796272 and at least two
  detected markers). The additions retain unmeasured objects separately.
- Regional CD3e, CD163, CD68, and LGALS3 positivity follows the same
  mean-plus-1.5-sample-SD rule within each region. Original values come from
  QuPath exports; added-donor values are recomputed from raw measurements.
- The protein-based UMAP and Leiden fits are refit over all 20 donors using
  Harmony on protein PCA. This gives the app one cohort-wide embedding and
  clustering rather than mixing the original 15-donor fit with a separate
  additions-only fit.
- App-wide donor-status controls use the same three labels for all 20 donors.

## Remaining method differences

- **Three hormone-threshold measures:** original INS/GCG/SST counts came from
  the historical QuPath exporter; the additions' counts were reconstructed in
  Python from raw intensity. They are no longer offered as Plot Population
  Measures. The source fields remain stored for provenance and legacy analysis.
- **Cell phenotypes:** the original 15 donors keep their saved historical
  labels. The five additions use the reviewed `protein_phenotype_rules.csv`
  reconstruction. Thus phenotype proportions, cell labels, neighborhood
  summaries, and the phenotype-count inputs to trajectory span two labeling
  procedures. Identical column names and complete coverage do not make those
  calls method-equivalent. The live trajectory now excludes the nine retired
  categories from both core and peri graph inputs. Previously those categories
  nearly separated the original and added donors by themselves (99.3% of raw
  DPT variance at 25–1902 cells was associated with cohort membership, measured
  on the 21-donor root). The
  corrected graph uses the same 12 retained categories as the app selectors;
  residual comparability limits remain. See the
  [trajectory input audit](../app/shiny_app/docs/trajectory.md#retired-feature-correction-2026-10-14).
- **Marker intensity scale:** the original `X` values preserve scimap's
  per-image piecewise min–max rescaling using historical gates. Additions use
  per-donor rescaling against recovered gates, including explicit handling for
  gates outside the observed range. Their marker values share a gate-oriented
  0–1 scale, but they were not produced by one joint rescaling fit; direct
  intensity comparisons and the embedding may retain donor-specific scale
  effects.
- **Region geometry:** original region summaries and coordinates remain from
  the frozen QuPath exports. Addition regions are recomputed from reviewed
  native QuPath boundaries, clipped to tissue minus exclusions, with a
  20-µm peri band. Region-based area, cell-membership, and distance outputs
  therefore come from distinct source geometries even though the app exposes
  matching fields.

These choices preserve the accepted original baseline and the accepted
addition reconstruction. The app now identifies the mixed phenotype provenance
where phenotype comparisons are selected. Claims of full processing
equivalence require a separate reprocessing decision; the present release
should be treated as a correctly grouped cohort with the method limits above.

Source details: [protein expansion record](protein_expansion_20261003.md),
[original processing trace](current_nontrajectory_processing.md), and the
expanded-v3 `protein_release.json` (20 donors, 7,028 islets, original subset
parity passed). A single reproducible workflow for all 20 donors is planned in
[the reprocessing plan](uniform_reprocessing_plan_20261015.md).
